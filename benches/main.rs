use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;
use rand::{seq::IteratorRandom, thread_rng};

const NUM_TAXA: usize = 4000;
const NORM: u32 = 1;

fn main() {
    divan::main();
}

#[divan::bench]
fn benchmark_constant_time_lca(bencher: divan::Bencher) {
    let mut tree = SimpleRootedTree::yule(NUM_TAXA);
    tree.precompute_constant_time_lca();
    bencher.bench(|| tree.get_lca_id(vec![10, 20].as_slice()));
}

#[divan::bench]
fn benchmark_lca(bencher: divan::Bencher) {
    let tree = SimpleRootedTree::yule(NUM_TAXA);
    bencher.bench(|| tree.get_lca_id(vec![10, 20].as_slice()));
}

#[divan::bench]
fn benchmark_yule(bencher: divan::Bencher) {
    bencher.bench(|| SimpleRootedTree::yule(NUM_TAXA));
}

#[divan::bench]
fn benchmark_precompute_rmq(bencher: divan::Bencher) {
    bencher
        .with_inputs(|| SimpleRootedTree::yule(NUM_TAXA))
        .bench_refs(|tree| {
            tree.precompute_constant_time_lca();
        });
}

#[divan::bench]
fn benchmark_cophen_dist_naive(bencher: divan::Bencher) {
    bencher
        .with_inputs(|| {
            fn depth(tree: &SimpleRootedTree, node_id: usize) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let mut t1 = SimpleRootedTree::yule(NUM_TAXA);
            let mut t2 = SimpleRootedTree::yule(NUM_TAXA);
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            let _ = t1.set_zeta(depth);
            let _ = t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist_naive(t2, NORM);
        });
}

#[divan::bench]
fn benchmark_contract(bencher: divan::Bencher) {
    bencher
        .with_inputs(|| {
            let mut rng = thread_rng();
            let mut t1 = SimpleRootedTree::yule(NUM_TAXA);
            let taxa_set = (0..NUM_TAXA).collect_vec();
            let taxa_subset = taxa_set.into_iter().choose_multiple(&mut rng, NUM_TAXA / 3);
            t1.precompute_constant_time_lca();
            (t1, taxa_subset)
        })
        .bench_refs(|(t1, taxa_subset)| {
            let _ = t1.contract_tree(taxa_subset.as_slice());
        });
}
