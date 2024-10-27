use std::collections::HashSet;

use itertools::Itertools;
use phylo::node::Node;
use phylo::prelude::*;
use phylo::tree::SimpleRootedTree;

#[test]
fn distance_matrix() {
    let input_str = String::from("((A:0.1,B:0.2):0.3,C:0.6);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    let matrix = tree.matrix();
    dbg!(&matrix);
}

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
    dbg!(
        &tree,
        tree.get_node(1).unwrap().get_children().collect_vec()
    );
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
    dbg!(&tree.get_node(1).unwrap().get_children().collect_vec());
    dbg!(&tree.dfs(tree.get_root_id()).collect_vec());
    dbg!(&tree.bfs_ids(tree.get_root_id()).collect_vec());
    dbg!(&tree.postord_ids(tree.get_root_id()).collect_vec());
    dbg!(&tree.euler_walk_ids(tree.get_root_id()).collect_vec());
    dbg!(&tree.dfs(tree.get_root_id()).collect_vec());
    dbg!(&tree.node_to_root(5).collect_vec());
    dbg!(&tree.root_to_node(5).collect_vec());
}
#[test]
fn read_small_tree() {
    let input_str = String::from("((A,B),C);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree.euler_walk_ids(tree.get_root_id()).collect_vec());
    let input_str = String::from("((A:0.1,B:0.2),C:0.6);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree.euler_walk_ids(tree.get_root_id()).collect_vec());
    dbg!(format!("{}", &tree.to_newick()));
    assert_eq!(
        &tree.get_taxa_space().collect::<HashSet<&String>>(),
        &HashSet::from([&"A".to_string(), &"B".to_string(), &"C".to_string()])
    );
    let input_str = String::from("((A:1e-3,B:2e-3),C:6e-3);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    for node in tree.postord_nodes(tree.get_root_id()) {
        dbg!(node.get_weight());
    }
}
#[test]
fn read_nexus() {
    let input_str =
        String::from("#NEXUS\n\nBEGIN TREES;\n\tTree tree=((A:1,B:1):1,(C:1,D:1):1);\nEND;");
    let tree = SimpleRootedTree::from_nexus(input_str.clone()).unwrap();
    assert_eq!(
        tree.to_newick().to_string(),
        "((A:1,B:1):1,(C:1,D:1):1);".to_string()
    );
    assert_eq!(tree.to_nexus().unwrap(), input_str);
}
#[test]
fn tree_spr() {
    let input_str = String::from("((A,B),C);");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_nodes().collect_vec());
    let p_tree = tree.prune(1).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(format!("{}", &p_tree.to_newick()));
    tree.graft(p_tree, (0, 4)).unwrap();
    tree.clean();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(&tree.get_node_parent(4));
    tree.spr((1, 2), (5, 4)).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
}
#[test]
fn tree_cluster() {
    let input_str: String = String::from("((A,B),C);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(&tree.get_cluster(0).collect_vec());
    dbg!(&tree.get_cluster(1).collect_vec());
    let bp = tree.get_bipartition((0, 1));
    dbg!(&bp.0.collect_vec());
    dbg!(&bp.1.collect_vec());
}
#[test]
fn balance_tree() {
    let input_str: String = String::from("(((A,B),C),D);");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    tree.balance_subtree().unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    let input_str: String = String::from("(D,(C,(A,B)));");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    tree.balance_subtree().unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    let input_str: String = String::from("(D,(A,(C,B)));");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    tree.balance_subtree().unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_nodes().collect_vec());
    dbg!(tree.get_root_id());
}
#[test]
fn induce_tree() {
    let input_str: String = String::from("(((A,B),C),D);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    let mut x = tree.induce_tree(vec![3, 5, 6]).unwrap();
    x.clean();
    dbg!(x.get_root().get_children().collect_vec());
    dbg!(x.get_nodes().collect_vec());
    dbg!(format!("{}", &x.to_newick()));
}
#[test]
fn median_node() {
    let input_str: String = String::from("(((A,B),C),D);");
    let tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    dbg!(format!("{}", &tree.to_newick()));
    dbg!(tree.get_cluster(tree.get_median_node_id()).collect_vec());
}

#[test]
fn yule() {
    let tree1 = SimpleRootedTree::yule(20);
    dbg!(format!("{}", &tree1.to_newick()));
}

#[test]
fn uniform() {
    let tree1 = SimpleRootedTree::unif(20);
    dbg!(format!("{}", &tree1.to_newick()));
}

#[test]
fn contract_tree() {
    fn depth(tree: &SimpleRootedTree, node_id: usize) -> f32 {
        EulerWalk::get_node_depth(tree, node_id) as f32
    }
    let mut tree = SimpleRootedTree::yule(10);
    tree.precompute_constant_time_lca();
    dbg!(&tree);
    tree.precompute_constant_time_lca();
    tree.set_zeta(depth).unwrap();
    println!("{}", tree.to_newick());
    let taxa_subset = vec![
        "1".to_string(),
        "4".to_string(),
        "3".to_string(),
        "7".to_string(),
    ]
    .into_iter()
    .map(|x| tree.get_taxa_node_id(&x).unwrap())
    .collect_vec();
    let mut new_tree = tree.contract_tree(taxa_subset.as_slice()).unwrap();
    println!("{}", new_tree.to_newick());
    new_tree.precompute_constant_time_lca();
}

#[test]
fn cophenetic_dist() {
    fn depth(tree: &SimpleRootedTree, node_id: usize) -> f32 {
        tree.depth(node_id) as f32
    }
    let t1_input_str: String = String::from("((A,B),C);");
    let t2_input_str: String = String::from("(A,(B,C));");
    let mut t1 = SimpleRootedTree::from_newick(t1_input_str.as_bytes()).unwrap();
    let mut t2 = SimpleRootedTree::from_newick(t2_input_str.as_bytes()).unwrap();

    t1.precompute_constant_time_lca();
    t2.precompute_constant_time_lca();

    t1.set_zeta(depth).unwrap();
    t2.set_zeta(depth).unwrap();

    assert_eq!(t1.cophen_dist_naive(&t2, 1), 4_f32);
}

#[test]
fn suppress_tree_node() {
    let input_str: String = String::from("(((A,B),C),D);");
    let mut tree = SimpleRootedTree::from_newick(input_str.as_bytes()).unwrap();
    tree.supress_node(2).expect("node id should be valid");
}
