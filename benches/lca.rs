//! Lowest-common-ancestor benchmarks: oracle construction, oracle queries, and
//! the uncached walk.
//!
//! Splitting build from query is the point. `LcaOracle` trades an O(n) Euler
//! tour plus RMQ build for O(1) queries, so the two halves have to be judged
//! separately — a build-time regression and a query-time regression mean very
//! different things.

use std::hint::black_box;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use itertools::Itertools;
use phylo::prelude::*;

mod common;
use common::{yule, TAXA};

/// One query against a prebuilt oracle.
///
/// Deliberately carries no throughput: the query is O(1), so tree size is not a
/// quantity of work and an "elements per second" figure would only be
/// misleading. A flat time across the sweep is the result to look for.
fn lca_oracle_query(c: &mut Criterion) {
    let mut group = c.benchmark_group("lca_oracle_query");
    for &taxa in TAXA {
        // Build the tree and its oracle once: the oracle borrows the tree, so both must outlive the timed closure.
        let tree = yule(taxa);
        let oracle = tree.lca();
        let leaves = tree.get_leaf_ids().take(2).collect_vec();

        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(oracle.get_lca_id(leaves.as_slice())))
        });
    }
    group.finish();
}

/// Building the oracle — the Euler tour and RMQ precomputation.
fn lca_oracle_build(c: &mut Criterion) {
    let mut group = c.benchmark_group("lca_oracle_build");
    for &taxa in TAXA {
        let tree = yule(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(tree.lca()))
        });
    }
    group.finish();
}

/// A single LCA without an oracle — the root-ward walk this all replaces.
///
/// Throughput is reported against tree size so this can be read directly
/// against `lca_oracle_build`: one uncached query versus the whole
/// precomputation.
fn lca_uncached(c: &mut Criterion) {
    let mut group = c.benchmark_group("lca_uncached");
    for &taxa in TAXA {
        let tree = yule(taxa);
        let leaves = tree.get_leaf_ids().take(2).collect_vec();

        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(tree.get_lca_id(leaves.as_slice())))
        });
    }
    group.finish();
}

criterion_group!(lca, lca_oracle_query, lca_oracle_build, lca_uncached);
criterion_main!(lca);
