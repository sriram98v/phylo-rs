//! Traversal and whole-tree summary benchmarks.
//!
//! Every one of these visits each node a constant number of times, so their
//! throughput figures are directly comparable: they say how many nodes per
//! second each traversal shape sustains, and a linear routine holds that number
//! flat across the sweep.

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use itertools::Itertools;
use phylo::prelude::*;

mod common;
use common::{yule, TAXA};

/// Post-order traversal, collected so the iterator is actually driven.
fn postorder_ids(c: &mut Criterion) {
    let mut group = c.benchmark_group("postorder_ids");
    for &taxa in TAXA {
        let tree = yule(taxa);
        let root = tree.get_root_id();
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(tree.postord_ids(root).unwrap().collect_vec()))
        });
    }
    group.finish();
}

/// Every bipartition, resolved to taxa on both sides.
///
/// This is the shape RF-style comparisons build on, so it is worth timing apart
/// from the metrics that consume it.
fn bipartitions(c: &mut Criterion) {
    let mut group = c.benchmark_group("bipartitions");
    // Quadratic in practice: a single iteration at 16k taxa runs for seconds,
    // so this group cannot afford criterion's default 100 samples.
    group
        .sample_size(10)
        .measurement_time(Duration::from_secs(30));
    for &taxa in TAXA {
        let tree = yule(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| {
                black_box(
                    tree.get_bipartitions_ids()
                        .map(|(c1, c2)| {
                            (
                                c1.map(|x| tree.get_node_taxa(x).unwrap()).collect_vec(),
                                c2.map(|x| tree.get_node_taxa(x).unwrap()).collect_vec(),
                            )
                        })
                        .collect_vec(),
                )
            })
        });
    }
    group.finish();
}

/// The median node — a full pass weighing each subtree.
fn median_node(c: &mut Criterion) {
    let mut group = c.benchmark_group("median_node");
    for &taxa in TAXA {
        let tree = yule(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(tree.get_median_node()))
        });
    }
    group.finish();
}

criterion_group!(traversal, postorder_ids, bipartitions, median_node);
criterion_main!(traversal);
