use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::PhyloTree;
use rand::{seq::IteratorRandom, thread_rng};

    
const NORM: u32 = 1;

fn main() {
    divan::main();
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000], sample_size = 1, sample_count = 10)]
fn benchmark_constant_time_lca(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let mut tree = PhyloTree::yule(taxa_size);
            tree.precompute_constant_time_lca();
            let leaves = vec![10, 20];
            let mut lca_map = vec![vec![0; taxa_size]; taxa_size];
            (0..taxa_size)
                // .map(|x| tree.get_taxa_node_id(&x.to_string()).unwrap())
                .combinations(2)
                .for_each(|x| {
                    let l_0 = tree.get_taxa_node_id(&x[0].to_string()).unwrap();
                    let l_1 = tree.get_taxa_node_id(&x[1].to_string()).unwrap();
                    lca_map[x[0]][x[1]] = tree.get_lca_id(vec![l_0, l_1].as_slice());
                    lca_map[x[1]][x[0]] = tree.get_lca_id(vec![l_0, l_1].as_slice())
                });
                // .map(|x| (x.clone(), tree.get_lca_id(x.as_slice())))
                // .collect::<HashMap<_,_>>();
            (lca_map, leaves)
        })
        .bench_refs(|(lca_map, leaves)| {
            lca_map[leaves[0]][leaves[1]];
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000],sample_size = 1, sample_count = 10)]
fn benchmark_lca(bencher: divan::Bencher, taxa_size: usize) {
    let tree = PhyloTree::yule(taxa_size);
    bencher.bench(|| tree.get_lca_id(vec![10, 20].as_slice()));
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_yule(bencher: divan::Bencher, taxa_size: usize) {
    bencher.bench(|| PhyloTree::yule(taxa_size))
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_precompute_rmq(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| PhyloTree::yule(taxa_size))
        .bench_refs(|tree| {
            tree.precompute_constant_time_lca();
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_spr(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let tree = PhyloTree::yule(taxa_size);
            let leaf_edges = tree.get_leaf_ids().map(|l_id| (tree.get_node_parent_id(l_id).unwrap(), l_id)).collect_vec();
            let e1 = leaf_edges[0];
            let e2 = leaf_edges[1];
            return (tree, e1, e2);
        })
        .bench_refs(|(t, e1, e2)| {
            let _ = t.spr(*e1, *e2);
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_cophen_dist_naive(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let mut t1 = PhyloTree::yule(taxa_size);
            let mut t2 = PhyloTree::yule(taxa_size);
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            let _ = t1.set_zeta(depth);
            let _ = t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist(t2, NORM);
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_rf(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            let t2 = PhyloTree::yule(taxa_size);

            (t1,t2)
        })
        .bench_refs(|(t1, t2)| {
            let _ = t1.rf(&t2);
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_cm(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            let t2 = PhyloTree::yule(taxa_size);

            (t1,t2)
        })
        .bench_refs(|(t1, t2)| {
            let _ = t1.cm(&t2);
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_bps(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            t1
        })
        .bench_refs(|t1| {
            let _ = t1.get_bipartitions_ids().map(|(c1,c2)| (c1.map(|x| t1.get_node_taxa(x).unwrap()).collect_vec(), c2.map(|x| t1.get_node_taxa(x).unwrap()).collect_vec())).collect_vec();
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_postord_ids(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            t1
        })
        .bench_refs(|t1| {
            let _ = t1.postord_ids(t1.get_root_id()).collect_vec();
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
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

#[divan::bench(args = [1000000])]
fn benchmark_contract(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let mut rng = thread_rng();
            let mut t1 = PhyloTree::yule(taxa_size);
            let taxa_set = (0..taxa_size).collect_vec();
            let taxa_subset = taxa_set
                .into_iter()
                .choose_multiple(&mut rng, ((taxa_size as f32)*0.05) as usize);
            t1.precompute_constant_time_lca();
            (t1, taxa_subset)
        })
        .bench_refs(|(t1, taxa_subset)| {
            t1.contract_tree(taxa_subset.as_slice()).unwrap();
        });
}

#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn new_contract_nodes(bencher: divan::Bencher, taxa_size: usize){
    bencher
        .with_inputs(|| {
            let mut rng = thread_rng();
            let mut t1 = PhyloTree::yule(taxa_size);
            let taxa_set = (0..taxa_size).collect_vec();
            let taxa_subset = taxa_set
                .into_iter()
                .choose_multiple(&mut rng, 3 * taxa_size / 4);
            t1.precompute_constant_time_lca();
            (t1, taxa_subset)
        })
        .bench_refs(|(t1, taxa_subset)| {
            t1.contracted_tree_nodes(taxa_subset.as_slice()).collect_vec();
        });
}


#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_median_node(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            let t1 = PhyloTree::yule(taxa_size);
            t1
        })
        .bench_refs(|t1| {
            let _ = t1.get_median_node();
        });
}

#[cfg(feature = "parallel")]
#[divan::bench(args = [200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600, 3800, 4000, 4200, 4400, 4600, 4800, 5000, 5200, 5400, 5600, 5800, 6000, 6200, 6400, 6600, 6800, 7000, 7200, 7400, 7600, 7800, 8000, 8200, 8400, 8600, 8800, 9000, 9200, 9400, 9600, 9800, 10000])]
fn benchmark_cophen_dist_par(bencher: divan::Bencher, taxa_size: usize) {
    bencher
        .with_inputs(|| {
            fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
                EulerWalk::get_node_depth(tree, node_id) as f32
            }

            let mut t1 = PhyloTree::yule(taxa_size);
            let mut t2 = PhyloTree::yule(taxa_size);
            t1.precompute_constant_time_lca();
            t2.precompute_constant_time_lca();
            let _ = t1.set_zeta(depth);
            let _ = t2.set_zeta(depth);
            (t1, t2)
        })
        .bench_refs(|(t1, t2)| {
            t1.cophen_dist_par(t2, NORM);
        });
}
