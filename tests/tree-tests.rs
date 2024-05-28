use std::collections::HashSet;
// use std::Arc::Arc;
use std::sync::Arc;
use std::time::{Duration, Instant};


use itertools::Itertools;
use phylo::iter::node_iter::{Ancestors, Clusters, EulerWalk, BFS, DFS};
use phylo::node::simple_rnode::{RootedTreeNode, RootedWeightedNode, RootedZetaNode};
use phylo::node::Node;
use phylo::tree::ops::{Balance, ContractTree, Subtree, SPR};
use phylo::tree::SimpleRootedTree;
use phylo::tree::simple_rtree::{RootedMetaTree, RootedTree};
use phylo::tree::distances::PathFunction;
use phylo::tree::io::*;
use phylo::tree::ops::CopheneticDistance;

#[test]
fn build_small_tree() {
    let mut tree = SimpleRootedTree::new(1);
    dbg!(&tree);
    let new_node = Node::new(2);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node = Node::new(3);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node: Node = Node::new(4);
    tree.add_child(2, new_node);
    let new_node: Node = Node::new(5);
    tree.add_child(2, new_node);
    dbg!(&tree, tree.get_node(1).unwrap().get_children().into_iter().collect_vec());
    dbg!(RootedTree::get_node_depth(&tree, 2));
    dbg!(&tree.to_newick().to_string());
    tree.clear();
    dbg!(&tree);

}

#[test]
fn tree_iter() {
    let mut tree = SimpleRootedTree::new(1);
    let new_node = Node::new(2);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node = Node::new(5);
    tree.add_child(tree.get_root_id(), new_node);
    let new_node: Node = Node::new(3);
    tree.add_child(2, new_node);
    let new_node: Node = Node::new(4);
    tree.add_child(2, new_node);
    let new_node: Node = Node::new(6);
    tree.add_child(5, new_node);
    let new_node: Node = Node::new(7);
    tree.add_child(5, new_node);
    dbg!(&tree.get_node(1).unwrap().get_children().into_iter().collect_vec());
    dbg!(&tree.dfs(tree.get_root_id()).into_iter().collect_vec());
    dbg!(&tree.bfs(tree.get_root_id()).collect_vec());
    dbg!(&tree.postord(tree.get_root_id()).collect_vec());
    dbg!(&tree.euler_walk(tree.get_root_id()).into_iter().collect_vec());
    dbg!(&tree.dfs(tree.get_root_id()).into_iter().collect_vec());
    dbg!(&tree.node_to_root(5).into_iter().collect_vec());
    dbg!(&tree.root_to_node(5).into_iter().collect_vec());
}
#[test]
fn read_small_tree() {
    let input_str = String::from("((A,B),C);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    dbg!(&tree.euler_walk(tree.get_root_id()).into_iter().collect_vec());
    let input_str = String::from("((A:0.1,B:0.2),C:0.6);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    dbg!(&tree.euler_walk(tree.get_root_id()).into_iter().collect_vec());
    dbg!(format!("{}", &tree.to_newick()));
    assert_eq!(&tree.get_taxa_space().into_iter().collect::<HashSet<String>>(), &HashSet::from(["A".to_string(), "B".to_string(), "C".to_string()]));
}
#[test]
fn tree_spr() {
    let input_str = String::from("((A,B),C);");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_nodes().into_iter().collect_vec());
    let p_tree = tree.prune(1);
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(format!("{}", &p_tree.to_newick()));
    tree.graft(p_tree, (0,4));
    tree.clean();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(&tree.get_node_parent(4));
    tree.spr((1,2), (5,4));
    dbg!(format!("{}", &tree.to_newick()));
}
#[test]
fn tree_cluster() {
    let input_str: String = String::from("((A,B),C);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    dbg!(&tree.get_cluster(0).into_iter().collect_vec());
    dbg!(&tree.get_cluster(1).into_iter().collect_vec());
    let bp = tree.get_bipartition((0, 1));
    dbg!(&bp.0.into_iter().collect_vec());
    dbg!(&bp.1.into_iter().collect_vec());
}
#[test]
fn balance_tree() {
    let input_str: String = String::from("(((A,B),C),D);");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    tree.balance_subtree();
    dbg!(format!("{}", &tree.to_newick()));
    let input_str: String = String::from("(D,(C,(A,B)));");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    tree.balance_subtree();
    dbg!(format!("{}", &tree.to_newick()));
    let input_str: String = String::from("(D,(A,(C,B)));");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    tree.balance_subtree();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_nodes().into_iter().collect_vec());
    dbg!(tree.get_root_id());
}
#[test]
fn induce_tree() {
    let input_str: String = String::from("(((A,B),C),D);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    dbg!(format!("{}", &tree.to_newick()));
    let mut x = tree.induce_tree(vec![3, 5, 6]);
    x.clean();
    dbg!(x.get_root().get_children().into_iter().collect_vec());
    dbg!(x.get_nodes().into_iter().collect_vec());
    dbg!(format!("{}", &x.to_newick()));
}
#[test]
fn median_node() {
    let input_str: String = String::from("(((A,B),C),D);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_cluster(tree.get_median_node_id()).into_iter().collect_vec());
}

#[test]
fn yule() {
    let tree1 = SimpleRootedTree::yule(20).ok().unwrap();
    dbg!(format!("{}", &tree1.to_newick()));
}

#[test]
fn uniform() {
    let tree1 = SimpleRootedTree::unif(20).ok().unwrap();
    dbg!(format!("{}", &tree1.to_newick()));
}

#[test]
fn benchmark_lca() {
    let num_taxa = 10000;
    println!("Simulating tree!");
    let mut tree = SimpleRootedTree::yule(num_taxa).unwrap();
    println!("Precomputing RMQ");
    tree.precompute_constant_time_lca();
    let num_pairs = (num_taxa as f32)*((num_taxa-1) as f32)/(2_f32);
    println!("Computing average time");
    let mean_time = tree.get_taxa_space()
        .tuple_combinations::<(_,_)>()
        .map(|(a,b)| {
            let a_id = tree.get_taxa_node_id(&a).unwrap();
            let b_id = tree.get_taxa_node_id(&b).unwrap();
            let now = Instant::now();
            let _ = tree.get_lca_id(&vec![a_id, b_id]);
            return now.elapsed();
            })
        .sum::<Duration>()/(num_pairs as u32);
    dbg!(mean_time);
}

#[test]
fn contract_tree() {
    fn depth(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID)->f32
    {
        EulerWalk::get_node_depth(tree, node_id) as f32
    }
    let mut tree = SimpleRootedTree::yule(10).unwrap();
    tree.precompute_constant_time_lca();
    // let input_str: String = String::from("(((A,B),C),D);");
    // let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes());
    // dbg!(&tree);
    // tree.precompute_constant_time_lca();
    tree.set_zeta(depth);
    println!("{}", tree.to_newick());
    let taxa_subset = vec!["1".to_string(), "4".to_string(), "3".to_string(), "7".to_string()].into_iter()
        .map(|x| tree.get_taxa_node_id(&x).unwrap())
        .collect_vec();
    let mut new_tree = tree.contract_tree(&taxa_subset);
    println!("{}", new_tree.to_newick());
    new_tree.precompute_constant_time_lca();
}


#[test]
fn benchmark_yule() {
    let num_taxa = 100;
    let num_iter = 100;
    println!("Computing average time");
    let mean_time = (0..num_iter).map(|_| {
            let now = Instant::now();
            let _ = SimpleRootedTree::yule(num_taxa).unwrap();
            return now.elapsed();
        })
        .sum::<Duration>()/(num_iter);
    dbg!(mean_time);
}

#[test]
fn benchmark_precompute_rmq() {
    let num_taxa = 100;
    let num_iter = 100;
    println!("Computing average time");
    let mean_rmq_time = (0..num_iter).map(|_| {
            let mut tree = SimpleRootedTree::yule(num_taxa).unwrap();
            let now = Instant::now();
            tree.precompute_constant_time_lca();
            return now.elapsed();
        })
        .sum::<Duration>()/(num_iter);
    dbg!(mean_rmq_time);
}

    
#[test]
fn cophenetic_dist() {
    fn depth(tree: &SimpleRootedTree, node_id: <SimpleRootedTree as RootedTree>::NodeID)->f32
    {
        tree.depth(node_id) as f32
    }
    let t1_input_str: String = String::from("((A,B),C);");
    let t2_input_str: String = String::from("(A,(B,C));");
    let mut t1 = SimpleRootedTree::from_newick(t1_input_str.as_bytes());
    let mut t2 = SimpleRootedTree::from_newick(t2_input_str.as_bytes());
    dbg!(format!("{}", &t1.to_newick()));
    dbg!(format!("{}", &t2.to_newick()));

    dbg!(&t1);
    dbg!(&t2);
    
    t1.precompute_constant_time_lca();
    t2.precompute_constant_time_lca();

    t1.set_zeta(depth);
    t2.set_zeta(depth);

    dbg!(&t1);
    dbg!(&t2);

    dbg!(format!("{}", &t1.cophen_dist_naive(&t2, 1)));
    assert_eq!(t1.cophen_dist_naive(&t2, 1), 2_f32);
}

