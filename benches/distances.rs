//! Tree-comparison metric benchmarks.
//!
//! These are quadratic in the number of taxa, so they run their own smaller
//! sweep and a reduced sample count. Throughput is still reported against taxa
//! count, which makes the quadratic term legible: each fourfold jump in taxa
//! should roughly quarter the elements-per-second figure.

use std::hint::black_box;
use std::time::Duration;

use criterion::{criterion_group, criterion_main, BenchmarkId, Criterion, Throughput};
use phylo::prelude::*;

mod common;
use common::{yule_pair, yule_pair_with_zeta, NORM, QUADRATIC_TAXA};

/// Robinson-Foulds distance.
fn robinson_foulds(c: &mut Criterion) {
    let mut group = c.benchmark_group("robinson_foulds");
    for &taxa in QUADRATIC_TAXA {
        let (t1, t2) = yule_pair(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(t1.rf(&t2)))
        });
    }
    group.finish();
}

/// Cluster matching distance.
fn cluster_matching(c: &mut Criterion) {
    let mut group = c.benchmark_group("cluster_matching");
    for &taxa in QUADRATIC_TAXA {
        let (t1, t2) = yule_pair(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(t1.cm(&t2)))
        });
    }
    group.finish();
}

/// Cluster affinity.
fn cluster_affinity(c: &mut Criterion) {
    let mut group = c.benchmark_group("cluster_affinity");
    for &taxa in QUADRATIC_TAXA {
        let (t1, t2) = yule_pair(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(t1.ca(&t2)))
        });
    }
    group.finish();
}

/// Cophenetic distance, serial.
///
/// Zeta is populated in setup, outside the timed closure — the metric reads it
/// but does not compute it.
fn cophenetic_distance(c: &mut Criterion) {
    let mut group = c.benchmark_group("cophenetic_distance");
    for &taxa in QUADRATIC_TAXA {
        let (t1, t2) = yule_pair_with_zeta(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(t1.cophen_dist(&t2, NORM)))
        });
    }
    group.finish();
}

/// Cophenetic distance under the `parallel` feature.
///
/// Named to sit next to `cophenetic_distance` in the report, so the speedup is
/// a direct read across the two groups at matching taxa counts.
#[cfg(feature = "parallel")]
fn cophenetic_distance_parallel(c: &mut Criterion) {
    let mut group = c.benchmark_group("cophenetic_distance_parallel");
    for &taxa in QUADRATIC_TAXA {
        let (t1, t2) = yule_pair_with_zeta(taxa);
        group.throughput(Throughput::Elements(taxa as u64));
        group.bench_with_input(BenchmarkId::from_parameter(taxa), &taxa, |b, _| {
            b.iter(|| black_box(t1.cophen_dist_par(&t2, NORM)))
        });
    }
    group.finish();
}

#[cfg(not(feature = "parallel"))]
fn cophenetic_distance_parallel(_: &mut Criterion) {}

/// Quadratic metrics need fewer, longer samples than criterion's defaults.
///
/// `cophenetic_distance` at 1000 taxa runs for seconds per iteration, so 10 —
/// criterion's minimum — is the only sample count that keeps the group to a
/// few minutes. The confidence intervals are correspondingly wide; read these
/// groups for order-of-magnitude and for the throughput decay across the sweep,
/// not for single-digit-percent regressions.
fn config() -> Criterion {
    Criterion::default()
        .sample_size(10)
        .measurement_time(Duration::from_secs(30))
}

criterion_group! {
    name = distances;
    config = config();
    targets = robinson_foulds,
              cluster_matching,
              cluster_affinity,
              cophenetic_distance,
              cophenetic_distance_parallel
}
criterion_main!(distances);
