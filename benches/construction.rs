//! Tree construction and structural-edit benchmarks.

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BatchSize, BenchmarkId, Criterion, Throughput};
use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::PhyloTree;

mod common;
use common::{first_two_leaf_edges, leaf_subset, yule, TAXA};

/// Simulating a random tree under the Yule model.
///
/// Worth watching the throughput column rather than the raw time: `yule` is
/// quadratic in the number of taxa, so elements per second falls off across the
/// sweep instead of holding flat. That decay is the thing to fix, and this
/// benchmark is what would show it fixed.
fn yule_simulation(c: &mut Criterion) {
    let mut group = c.benchmark_group("yule_simulation");
    group
        .sample_size(20)
        .measurement_time(Duration::from_secs(10));
    for &taxa in TAXA {
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, &taxa| {
            b.iter(|| black_box(PhyloTree::yule(taxa)))
        });
    }
    group.finish();
}

/// One subtree prune-and-regraft.
///
/// `spr` takes `&mut self`, so each iteration needs its own tree. The tree is
/// built once and *cloned* per iteration rather than re-simulated: cloning an
/// arena is linear, while `yule` is quadratic, so the per-iteration setup
/// criterion has to run between samples stops dominating the wall clock.
fn spr(c: &mut Criterion) {
    let mut group = c.benchmark_group("spr");
    group.sample_size(20);
    for &taxa in TAXA {
        let template = yule(taxa);
        let (edge1, edge2) = first_two_leaf_edges(&template);

        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter_batched_ref(
                || template.clone(),
                |tree| {
                    let _ = tree.spr(edge1, edge2);
                },
                BatchSize::LargeInput,
            )
        });
    }
    group.finish();
}

/// Contracting a tree down to 5% of its leaves, returning a new tree.
fn contract_tree(c: &mut Criterion) {
    let mut group = c.benchmark_group("contract_tree");
    for &taxa in TAXA {
        let tree = yule(taxa);
        let subset = leaf_subset(&tree, 0.05);

        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(tree.contract_tree(subset.as_slice()).unwrap()))
        });
    }
    group.finish();
}

/// The same contraction, against a caller-supplied LCA oracle.
///
/// `contract_tree` builds a throwaway `LcaOracle` per call, and that build is
/// O(n) while the contraction itself is proportional to the subset. A caller
/// contracting many subsets against one tree pays for the index once. Running
/// both groups is what separates the two costs — comparing this group against
/// `contract_tree` at the same taxa count gives the build's share directly.
fn contract_tree_with_oracle(c: &mut Criterion) {
    let mut group = c.benchmark_group("contract_tree_with_oracle");
    for &taxa in TAXA {
        let tree = yule(taxa);
        let subset = leaf_subset(&tree, 0.05);
        let oracle = tree.lca();

        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| {
                black_box(
                    tree.contract_tree_with_oracle(subset.as_slice(), &oracle)
                        .unwrap(),
                )
            })
        });
    }
    group.finish();
}

/// The node stream behind `contract_tree`, without building the result tree.
///
/// Timed separately so a regression can be attributed to the contraction walk
/// or to arena construction, rather than to "contraction" as a whole. Uses a
/// much larger leaf subset (75%), where the walk rather than the allocation
/// dominates.
///
/// Takes the `_from_root` form with the root resolved outside the timed
/// closure. The convenience form resolves it by `get_lca_id`, which builds a
/// throwaway oracle — that build is O(n) and would swamp the walk this group
/// exists to isolate.
fn contracted_tree_nodes(c: &mut Criterion) {
    let mut group = c.benchmark_group("contracted_tree_nodes");
    for &taxa in TAXA {
        let tree = yule(taxa);
        let subset = leaf_subset(&tree, 0.75);
        let root = tree.get_lca_id(subset.as_slice()).unwrap();

        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| {
                black_box(
                    tree.contracted_tree_nodes_from_root(root, subset.as_slice())
                        .collect_vec(),
                )
            })
        });
    }
    group.finish();
}

criterion_group!(
    construction,
    yule_simulation,
    spr,
    contract_tree,
    contract_tree_with_oracle,
    contracted_tree_nodes
);
criterion_main!(construction);
