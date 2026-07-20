use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::PhyloTree;
use rand::{seq::IteratorRandom, thread_rng};

const NORM: u32 = 1;

/// Taxa counts for benchmarks that scale near-linearly.
///
/// Geometric rather than linear: the point of sweeping at all is to expose the
/// growth rate, and a linear sweep spends most of its time re-measuring the
/// same complexity class.
const TAXA: &[usize] = &[1000, 4000, 16000];

/// Taxa counts for benchmarks that are quadratic in the number of taxa.
///
/// These compare every node against every node, so they need their own, much
/// smaller sweep to stay runnable.
const QUADRATIC_TAXA: &[usize] = &[250, 500, 1000];

fn main() {
    divan::main();
}

/// One LCA query against a precomputed index.
///
/// This previously timed `lca_map[i][j]` -- an index into a `Vec<Vec<usize>>`
/// the setup had filled in. The real `get_lca_id` calls all happened inside
/// `with_inputs`, which divan excludes from timing, so the benchmark reported
/// the cost of a vector index and never touched the RMQ. It also built an
/// n-by-n map, which is why it needed the small sweep.
///
/// It now times the query itself, which is the point: this is what an RMQ
/// change has to be judged against.
#[divan::bench(args = TAXA)]
fn benchmark_constant_time_lca(bencher: divan::Bencher, taxa_size: usize) {
    // Build the tree and its oracle once, then time only the constant-time
    // query -- the oracle borrows the tree, so both must outlive the closure.
    let tree = PhyloTree::yule(taxa_size);
    let oracle = tree.lca();
    let leaves = tree.get_leaf_ids().take(2).collect_vec();
    bencher.bench(|| oracle.get_lca_id(leaves.as_slice()));
}

#[divan::bench(args = TAXA,sample_size = 1, sample_count = 10)]
fn benchmark_lca(bencher: divan::Bencher, taxa_size: usize) {
    let tree = PhyloTree::yule(taxa_size);
    bencher.bench(|| tree.get_lca_id(vec![10, 20].as_slice()));
}

#[divan::bench(args = TAXA)]
fn benchmark_yule(bencher: divan::Bencher, taxa_size: usize) {
    bencher.bench(|| PhyloTree::yule(taxa_size))
}

#[divan::bench(args = TAXA)]
fn benchmark_precompute_rmq(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| PhyloTree::yule(taxa_size))
        .bench_refs(|tree| {
            divan::black_box(tree.lca());
        });
}

#[divan::bench(args = TAXA)]
fn benchmark_spr(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let tree = PhyloTree::yule(taxa_size);
            let leaf_edges = tree
                .get_leaf_ids()
                .map(|l_id| (tree.get_node_parent_id(l_id).unwrap(), l_id))
                .collect_vec();
            let e1 = leaf_edges[0];
            let e2 = leaf_edges[1];
            (tree, e1, e2)
        })
        .bench_refs(|(t, e1, e2)| {
            let _ = t.spr(*e1, *e2);
        });
}

#[divan::bench(args = QUADRATIC_TAXA)]
fn benchmark_cophen_dist_naive(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
                tree.depth(node_id) as f32
            }

            let mut t1 = PhyloTree::yule(taxa_size);
            let mut t2 = PhyloTree::yule(taxa_size);
            let _ = t1.set_zeta(depth);
            let _ = t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist(t2, NORM);
        });
}

#[divan::bench(args = QUADRATIC_TAXA)]
fn benchmark_rf(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            let t2 = PhyloTree::yule(taxa_size);

            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            let _ = t1.rf(t2);
        });
}

#[divan::bench(args = QUADRATIC_TAXA)]
fn benchmark_cm(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            let t2 = PhyloTree::yule(taxa_size);

            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            let _ = t1.cm(t2);
        });
}

#[divan::bench(args = TAXA)]
fn benchmark_bps(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| PhyloTree::yule(taxa_size))
        .bench_refs(|t1| {
            let _ = t1
                .get_bipartitions_ids()
                .map(|(c1, c2)| {
                    (
                        c1.map(|x| t1.get_node_taxa(x).unwrap()).collect_vec(),
                        c2.map(|x| t1.get_node_taxa(x).unwrap()).collect_vec(),
                    )
                })
                .collect_vec();
        });
}

#[divan::bench(args = TAXA)]
fn benchmark_postord_ids(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| PhyloTree::yule(taxa_size))
        .bench_refs(|t1| {
            let _ = t1.postord_ids(t1.get_root_id()).collect_vec();
        });
}

#[divan::bench(args = QUADRATIC_TAXA)]
fn benchmark_ca(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            let t2 = PhyloTree::yule(taxa_size);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            let _ = t1.ca(t2);
        });
}

// 1_000_000 taxa is unrunnable while `yule` is quadratic; the setup alone
// dominates any measurement of `contract_tree`.
#[divan::bench(args = [100_000])]
fn benchmark_contract(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let mut rng = thread_rng();
            let t1 = PhyloTree::yule(taxa_size);
            let taxa_set = (0..taxa_size).collect_vec();
            let taxa_subset = taxa_set
                .into_iter()
                .choose_multiple(&mut rng, ((taxa_size as f32) * 0.05) as usize);
            (t1, taxa_subset)
        })
        .bench_refs(|(t1, taxa_subset)| {
            t1.contract_tree(taxa_subset.as_slice()).unwrap();
        });
}

#[divan::bench(args = TAXA)]
fn new_contract_nodes(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let mut rng = thread_rng();
            let t1 = PhyloTree::yule(taxa_size);
            let taxa_set = (0..taxa_size).collect_vec();
            let taxa_subset = taxa_set
                .into_iter()
                .choose_multiple(&mut rng, 3 * taxa_size / 4);
            (t1, taxa_subset)
        })
        .bench_refs(|(t1, taxa_subset)| {
            t1.contracted_tree_nodes(taxa_subset.as_slice())
                .collect_vec();
        });
}

#[divan::bench(args = TAXA)]
fn benchmark_median_node(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| PhyloTree::yule(taxa_size))
        .bench_refs(|t1| {
            let _ = t1.get_median_node();
        });
}

#[cfg(feature = "parallel")]
#[divan::bench(args = QUADRATIC_TAXA)]
fn benchmark_cophen_dist_par(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
                tree.depth(node_id) as f32
            }

            let mut t1 = PhyloTree::yule(taxa_size);
            let mut t2 = PhyloTree::yule(taxa_size);
            let _ = t1.set_zeta(depth);
            let _ = t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist_par(t2, NORM);
        });
}
