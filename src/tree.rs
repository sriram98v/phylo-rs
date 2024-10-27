#![allow(clippy::needless_lifetimes)]

/// Module with traits and structs for distance computation
pub mod distances;
/// Module with traits and structs for tree encoding
pub mod io;
/// Module with traits and structs for tree operations
pub mod ops;
/// Module with traits and structs for general tree traits
pub mod simple_rtree;
/// Module with traits and structs for tree simulation
pub mod simulation;

use core::f32;
use fxhash::FxHashMap as HashMap;
use simulation::{Uniform, Yule};
use std::ops::Index;

use itertools::Itertools;
use rand::prelude::IteratorRandom;

use crate::prelude::*;

use crate::iter::{BFSIterator, DFSPostOrderIterator};
use crate::node::{Node, NodeID};
use vers_vecs::BinaryRmq;

/// Arena memory-managed tree struct
#[derive(Debug, Clone)]
pub struct SimpleRootedTree {
    /// Root NodeID
    root: NodeID,
    /// Nodes of the tree
    nodes: Vec<Option<Node>>,
    /// Index of nodes by taxa
    taxa_node_id_map: HashMap<String, NodeID>,
    /// Field to hold precomputed euler tour for constant-time LCA queries
    precomputed_euler: Option<Vec<NodeID>>,
    /// Field to hold precomputed first-appearance for constant-time LCA queries
    precomputed_fai: Option<Vec<Option<usize>>>,
    /// Field to hold precomputed depth-array for constant-time LCA queries
    precomputed_da: Option<Vec<usize>>,
    /// Field to hold precomputed range-minimum-query for constant-time LCA queries
    precomputed_rmq: Option<BinaryRmq>,
}

impl SimpleRootedTree {
    /// Creates new empty tree
    pub fn new(root_id: NodeID) -> Self {
        let root_node = Node::new(root_id);
        let mut nodes = vec![None; root_id + 1];
        nodes[root_id] = Some(root_node);
        SimpleRootedTree {
            root: root_id,
            nodes,
            precomputed_euler: None,
            taxa_node_id_map: [].into_iter().collect::<HashMap<_, _>>(),
            precomputed_fai: None,
            precomputed_da: None,
            precomputed_rmq: None,
        }
    }

    /// Returns new empty tree
    pub fn next_id(&self) -> usize {
        match &self.nodes.iter().position(|r| r.is_none()) {
            Some(x) => *x,
            None => self.nodes.len(),
        }
    }

    /// Creates new node with next NodeID
    pub fn next_node(&self) -> Node {
        Node::new(self.next_id())
    }
}

impl RootedTree for SimpleRootedTree {
    type Node = Node;

    /// Returns reference to node by ID
    fn get_node<'a>(&'a self, node_id: TreeNodeID<Self>) -> Option<&'a Node> {
        self.nodes[node_id].as_ref()
    }

    fn get_node_mut(&mut self, node_id: TreeNodeID<Self>) -> Option<&mut Node> {
        self.nodes[node_id].as_mut()
    }

    fn get_node_ids(&self) -> impl Iterator<Item = TreeNodeID<Self>> {
        (0..self.nodes.len()).filter(|x| self.nodes[*x].is_some())
    }

    fn get_nodes_mut<'a>(&'a mut self) -> impl Iterator<Item = &'a mut Self::Node> {
        self.nodes.iter_mut().filter_map(|x| x.as_mut())
    }

    fn set_node(&mut self, node: Node) {
        let node_id = node.get_id();
        match node.get_taxa() {
            None => {}
            Some(taxa) => {
                self.taxa_node_id_map
                    .insert(taxa.to_string(), node.get_id());
            }
        };
        match self.nodes.len() > node_id {
            true => self.nodes[node_id] = Some(node),
            false => {
                let new_len = node.get_id() - self.nodes.len();
                self.nodes.extend((0..new_len + 1).map(|_| None));
                self.nodes[node_id] = Some(node);
            }
        };
    }

    fn get_root_id(&self) -> TreeNodeID<Self> {
        self.root
    }

    fn set_root(&mut self, node_id: TreeNodeID<Self>) {
        self.root = node_id;
    }

    fn remove_node(&mut self, node_id: TreeNodeID<Self>) -> Option<Node> {
        if let Some(pid) = self.get_node_parent_id(node_id) {
            self.get_node_mut(pid).unwrap().remove_child(&node_id)
        }
        self.nodes[node_id].take()
    }

    fn delete_node(&mut self, node_id: TreeNodeID<Self>) {
        let _ = self.nodes[node_id].take();
    }

    fn clear(&mut self) {
        let root_node = self.get_root().clone();
        let root_node_id = root_node.get_id();
        self.nodes = vec![None; root_node_id + 1];
        self.nodes[root_node_id] = Some(root_node);
        self.taxa_node_id_map.clear();
    }

    /// Supresses all nodes of degree 2
    fn supress_unifurcations<'a>(&'a mut self) {}
}

impl RootedMetaTree for SimpleRootedTree {
    fn get_taxa_node(&self, taxa: &TreeNodeMeta<Self>) -> Option<&Self::Node> {
        self.get_node(*self.taxa_node_id_map.get(taxa)?)
    }

    fn set_node_taxa(&mut self, node_id: TreeNodeID<Self>, taxa: Option<TreeNodeMeta<Self>>) {
        self.get_node_mut(node_id).unwrap().set_taxa(taxa.clone());
        if let Some(t) = taxa {
            self.taxa_node_id_map.insert(t, node_id);
        }
    }

    fn num_taxa(&self) -> usize {
        self.taxa_node_id_map.len()
    }

    fn get_taxa_space<'a>(&'a self) -> impl ExactSizeIterator<Item = &'a TreeNodeMeta<Self>> {
        self.taxa_node_id_map.keys()
    }

    fn get_node_taxa_cloned(&self, node_id: TreeNodeID<Self>) -> Option<TreeNodeMeta<Self>> {
        self.get_node(node_id).unwrap().get_taxa().cloned()
    }
}

impl Yule for SimpleRootedTree {
    fn yule(num_taxa: usize) -> SimpleRootedTree {
        let mut tree = SimpleRootedTree::new(0);
        if num_taxa < 3 {
            return tree;
        }
        let new_node = Node::new(1);
        tree.add_child(0, new_node);
        tree.set_node_taxa(1, Some("0".to_string()));
        let new_node = Node::new(2);
        tree.add_child(0, new_node);
        tree.set_node_taxa(2, Some("1".to_string()));
        if num_taxa < 4 {
            return tree;
        }
        let mut current_leaf_ids = vec![1, 2];
        for i in 2..num_taxa {
            let rand_leaf_id = current_leaf_ids
                .iter()
                .choose(&mut rand::thread_rng())
                .unwrap();
            let rand_leaf_parent_id = tree.get_node_parent_id(*rand_leaf_id).unwrap();
            let split_node = Node::new(tree.next_id());
            let split_node_id = split_node.get_id();
            tree.split_edge((rand_leaf_parent_id, *rand_leaf_id), split_node);
            let new_leaf = Node::new(tree.next_id());
            let new_leaf_id = new_leaf.get_id();
            tree.add_child(split_node_id, new_leaf);
            tree.set_node_taxa(new_leaf_id, Some(i.to_string()));
            current_leaf_ids.push(new_leaf_id);
        }
        tree
    }
}

impl Uniform for SimpleRootedTree {
    fn unif(num_taxa: usize) -> SimpleRootedTree {
        let mut tree = SimpleRootedTree::new(0);
        if num_taxa < 3 {
            return tree;
        }
        let new_node = Node::new(1);
        tree.add_child(0, new_node);
        tree.set_node_taxa(1, Some("0".to_string()));
        let new_node = Node::new(2);
        tree.add_child(0, new_node);
        tree.set_node_taxa(2, Some("1".to_string()));
        if num_taxa < 3 {
            return tree;
        }
        let mut current_node_ids = vec![1, 2];
        for i in 1..num_taxa {
            let rand_leaf_id = *current_node_ids
                .iter()
                .choose(&mut rand::thread_rng())
                .unwrap();
            let rand_leaf_parent_id = tree.get_node_parent_id(rand_leaf_id).unwrap();
            let split_node = Node::new(tree.next_id());
            let split_node_id = split_node.get_id();
            current_node_ids.push(split_node_id);
            tree.split_edge((rand_leaf_parent_id, rand_leaf_id), split_node);
            let mut new_leaf = Node::new(tree.next_id());
            new_leaf.set_taxa(Some(i.to_string()));
            current_node_ids.push(new_leaf.get_id());
            tree.add_child(split_node_id, new_leaf)
        }
        tree
    }
}

impl RootedWeightedTree for SimpleRootedTree {
    fn unweight(&mut self) {
        self.nodes
            .iter_mut()
            .filter(|x| x.is_none())
            .for_each(|x| x.as_mut().unwrap().unweight());
    }
}

impl PathFunction for SimpleRootedTree {}

impl Ancestors for SimpleRootedTree {}

impl Subtree for SimpleRootedTree {}

impl PreOrder for SimpleRootedTree {}

impl ClusterAffinity for SimpleRootedTree {}

impl RobinsonFoulds for SimpleRootedTree {}

impl DistanceMatrix for SimpleRootedTree {
    fn matrix(&self) -> Vec<Vec<TreeNodeWeight<Self>>> {
        let mut out_mat = vec![vec![f32::INFINITY; self.nodes.len()]; self.nodes.len()];
        for node_ids in self.postord_ids(self.get_root_id()).combinations(2) {
            let n1 = node_ids[0];
            let n2 = node_ids[1];
            out_mat[n1][n1] = 0_f32;
            out_mat[n2][n2] = 0_f32;
            out_mat[n1][n2] = self.pairwise_distance(n1, n2);
            out_mat[n2][n1] = out_mat[n1][n2];
        }
        out_mat
    }

    fn pairwise_distance(
        &self,
        node_id_1: TreeNodeID<Self>,
        node_id_2: TreeNodeID<Self>,
    ) -> TreeNodeWeight<Self> {
        let lca = self.get_lca_id(vec![node_id_1, node_id_2].as_slice());
        let d1: TreeNodeWeight<Self> = self
            .node_to_root_ids(node_id_1)
            .map(|x| match x == self.get_root_id() {
                true => 0_f32,
                false => self.get_edge_weight(0, x).unwrap_or(1_f32),
            })
            .sum();

        let d2: TreeNodeWeight<Self> = self
            .node_to_root_ids(node_id_2)
            .map(|x| match x == self.get_root_id() {
                true => 0_f32,
                false => self.get_edge_weight(0, x).unwrap_or(1_f32),
            })
            .sum();

        let dlca: TreeNodeWeight<Self> = self
            .node_to_root_ids(lca)
            .map(|x| match x == self.get_root_id() {
                true => 0_f32,
                false => self.get_edge_weight(0, x).unwrap_or(1_f32),
            })
            .sum();

        d1 + d2 - 2_f32 * dlca
    }
}

impl DFS for SimpleRootedTree {
    fn postord_ids(&self, start_node: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>> {
        DFSPostOrderIterator::new(self, start_node).map(|x| x.get_id())
    }

    fn postord_nodes<'a>(
        &'a self,
        start_node: TreeNodeID<Self>,
    ) -> impl Iterator<Item = &'a Self::Node> {
        DFSPostOrderIterator::new(self, start_node)
    }
}

impl BFS for SimpleRootedTree {
    fn bfs_nodes<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl Iterator<Item = &'a Self::Node> {
        BFSIterator::new(self, start_node_id)
    }

    fn bfs_ids(&self, start_node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>> {
        BFSIterator::new(self, start_node_id).map(|x| x.get_id())
    }
}

impl ContractTree for SimpleRootedTree {
    fn contracted_tree_nodes(
        &self,
        leaf_ids: &[TreeNodeID<Self>],
    ) -> impl Iterator<Item = Self::Node> {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let node_postord_iter = self.postord_nodes(new_tree_root_id);
        let mut node_map: Vec<Option<Self::Node>> = vec![None; self.nodes.len()];
        node_map[new_tree_root_id] = Some(self.get_lca(leaf_ids).clone());
        let mut remove_list = vec![];
        node_postord_iter.for_each(|orig_node| {
            let mut node = orig_node.clone();
            match node.is_leaf() {
                true => if leaf_ids.contains(&node.get_id()) {
                    node_map[node.get_id()] = Some(node.clone());
                },
                false => {
                    let node_children_ids = node.get_children().collect_vec();
                    for child_id in &node_children_ids {
                        match node_map[*child_id].is_some() {
                            true => {}
                            false => node.remove_child(child_id),
                        }
                    }
                    let node_children_ids = node.get_children().collect_vec();
                    match node_children_ids.len() {
                        0 => {}
                        1 => {
                            // the node is a unifurcation
                            // node should be added to both node_map and remove_list
                            // if child of node is already in remove list, attach node children to node first
                            let child_node_id = node_children_ids[0];

                            if remove_list.contains(&child_node_id) {
                                node.remove_child(&child_node_id);
                                let grandchildren_ids = node_map[child_node_id]
                                    .as_mut()
                                    .unwrap()
                                    .get_children()
                                    .collect_vec();
                                for grandchild_id in grandchildren_ids {
                                    node_map[grandchild_id]
                                        .as_mut()
                                        .unwrap()
                                        .set_parent(Some(node.get_id()));
                                    node.add_child(grandchild_id);
                                }
                            }
                            let n_id = node.get_id();
                            remove_list.push(n_id);
                            node_map[n_id] = Some(node.clone());
                        }
                        _ => {
                            // node has multiple children
                            // for each child, suppress child if child is in remove list
                            node_children_ids.into_iter().for_each(|chid| {
                                if remove_list.contains(&chid) {
                                    // suppress chid
                                    // remove chid from node children
                                    // children of chid are node grandchildren
                                    // add grandchildren to node children
                                    // set grandchildren parent to node
                                    node.remove_child(&chid);
                                    let node_grandchildren = node_map[chid]
                                        .as_mut()
                                        .unwrap()
                                        .get_children()
                                        .collect_vec();
                                    for grandchild in node_grandchildren {
                                        node.add_child(grandchild);
                                        node_map[grandchild]
                                            .as_mut()
                                            .unwrap()
                                            .set_parent(Some(node.get_id()))
                                    }
                                }
                            });
                            if node.get_id() == new_tree_root_id {
                                node.set_parent(None);
                            }
                            node_map[node.get_id()] = Some(node.clone());
                        }
                    };
                }
            }
        });
        remove_list.into_iter().for_each(|x| {
            node_map[x] = None;
        });
        node_map.into_iter().flatten()
    }

    fn contract_tree(&self, leaf_ids: &[TreeNodeID<Self>]) -> Result<Self, ()> {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let new_nodes = self.contracted_tree_nodes(leaf_ids).collect_vec();
        let mut new_tree = self.clone();
        new_tree.set_root(new_tree_root_id);
        new_tree.clear();
        new_tree.set_nodes(new_nodes);
        Ok(new_tree)
    }

    fn contract_tree_from_iter(
        &self,
        leaf_ids: &[TreeNodeID<Self>],
        node_iter: impl Iterator<Item = TreeNodeID<Self>>,
    ) -> Result<Self, ()> {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let new_nodes = self
            .contracted_tree_nodes_from_iter(new_tree_root_id, leaf_ids, node_iter)
            .collect_vec();
        let mut new_tree = self.clone();
        new_tree.set_root(new_tree_root_id);
        new_tree.clear();
        new_tree.set_nodes(new_nodes);
        Ok(new_tree)
    }
}

impl EulerWalk for SimpleRootedTree {
    fn precompute_fai(&mut self) {
        // todo!()
        let max_id = self.get_node_ids().max().unwrap();
        let mut index = vec![None; max_id + 1];
        match self.get_precomputed_walk() {
            Some(walk) => {
                for node_id in self.get_node_ids() {
                    index[node_id] = Some(walk.iter().position(|x| x == &node_id).unwrap());
                }
            }
            None => {
                let walk = self.euler_walk_ids(self.get_root_id()).collect_vec();
                for node_id in self.get_node_ids() {
                    index[node_id] = Some(walk.iter().position(|x| x == &node_id).unwrap());
                }
            }
        }
        let mut fai_vec = vec![None; max_id + 1];
        for id in self.get_node_ids() {
            fai_vec[id] = index[id];
        }
        self.precomputed_fai = Some(fai_vec);
    }

    fn precompute_da(&mut self) {
        let da = self.depth_array();
        let rmq = BinaryRmq::from_vec(da.iter().map(|x| *x as u64).collect_vec());
        self.precomputed_rmq = Some(rmq);

        self.precomputed_da = Some(self.depth_array());
    }

    fn precompute_walk(&mut self) {
        self.precomputed_euler = Some(self.euler_walk_ids(self.get_root_id()).collect_vec());
    }

    fn get_precomputed_walk(
        &self,
    ) -> Option<&Vec<<<Self as RootedTree>::Node as RootedTreeNode>::NodeID>> {
        self.precomputed_euler.as_ref()
    }

    fn precompute_constant_time_lca(&mut self) {
        self.precompute_walk();
        self.precompute_da();
        self.precompute_fai();
    }

    fn get_precomputed_fai(&self) -> Option<impl Index<TreeNodeID<Self>, Output = Option<usize>>> {
        self.precomputed_fai.clone()
    }

    fn get_precomputed_da(&self) -> Option<&Vec<usize>> {
        self.precomputed_da.as_ref()
    }

    fn is_euler_precomputed(&self) -> bool {
        self.precomputed_euler.is_some()
    }

    fn first_appearance(&self) -> impl Index<TreeNodeID<Self>, Output = Option<usize>> {
        let max_id = self.get_node_ids().max().unwrap();
        let mut fa = vec![None; max_id + 1];
        match self.get_precomputed_walk() {
            Some(walk) => {
                for node_id in self.get_node_ids() {
                    fa[node_id] = Some(walk.iter().position(|x| x == &node_id).unwrap());
                }
            }
            None => {
                let walk = self.euler_walk_ids(self.get_root_id()).collect_vec();
                for node_id in self.get_node_ids() {
                    fa[node_id] = Some(walk.iter().position(|x| x == &node_id).unwrap());
                }
            }
        }

        fa
    }

    fn get_fa_index(&self, node_id: TreeNodeID<Self>) -> usize {
        match &self.precomputed_fai {
            Some(fai) => fai[node_id].unwrap(),
            None => self.first_appearance()[node_id].unwrap(),
        }
    }

    fn get_lca_id(&self, node_id_vec: &[TreeNodeID<Self>]) -> TreeNodeID<Self> {
        let min_pos = node_id_vec
            .iter()
            .map(|x| self.get_fa_index(*x))
            .min()
            .unwrap();
        let max_pos = node_id_vec
            .iter()
            .map(|x| self.get_fa_index(*x))
            .max()
            .unwrap();

        if min_pos == max_pos {
            return min_pos;
        }

        match self.precomputed_rmq.as_ref() {
            Some(dp) => self.get_euler_pos(dp.range_min(min_pos, max_pos)),
            None => {
                let da = self.depth_array();
                let rmq = BinaryRmq::from_vec(da.iter().map(|x| *x as u64).collect_vec());
                self.get_euler_pos(rmq.range_min(min_pos, max_pos))
            }
        }
    }
}

impl Clusters for SimpleRootedTree {}

impl Newick for SimpleRootedTree {
    fn from_newick(newick_str: &[u8]) -> std::io::Result<Self> {
        let mut tree = SimpleRootedTree::new(0);
        let mut stack: Vec<TreeNodeID<Self>> = Vec::new();
        let mut context: TreeNodeID<Self> = tree.get_root_id();
        let mut taxa_str = String::new();
        let mut decimal_str: String = String::new();
        let mut str_ptr: usize = 0;
        let newick_string = String::from_utf8(newick_str.to_vec())
            .unwrap()
            .chars()
            .filter(|c| !c.is_whitespace())
            .collect::<Vec<char>>();
        while str_ptr < newick_string.len() {
            match newick_string[str_ptr] {
                '(' => {
                    stack.push(context);
                    let new_node = Node::new(tree.next_id());
                    context = new_node.get_id();
                    tree.set_node(new_node);
                    str_ptr += 1;
                }
                ')' | ',' => {
                    // last context id
                    let last_context = stack.last().ok_or_else(|| {
                        std::io::Error::new(
                            std::io::ErrorKind::InvalidData,
                            NewickError::InvalidCharacter { idx: str_ptr },
                        )
                    })?;
                    // add current context as a child to last context
                    tree.set_child(*last_context, context);
                    tree.set_edge_weight(
                        (*last_context, context),
                        decimal_str.parse::<TreeNodeWeight<Self>>().ok(),
                    );
                    if !taxa_str.is_empty() {
                        tree.set_node_taxa(context, Some(taxa_str.to_string()));
                    }
                    // we clear the strings
                    taxa_str.clear();
                    decimal_str.clear();

                    match newick_string[str_ptr] {
                        ',' => {
                            let new_node = Node::new(tree.next_id());
                            context = new_node.get_id();
                            tree.set_node(new_node);
                            str_ptr += 1;
                        }
                        _ => {
                            context = stack.pop().ok_or_else(|| {
                                std::io::Error::new(
                                    std::io::ErrorKind::InvalidData,
                                    NewickError::InvalidCharacter { idx: str_ptr },
                                )
                            })?;
                            str_ptr += 1;
                        }
                    }
                }
                ';' => {
                    if !taxa_str.is_empty() {
                        tree.set_node_taxa(context, Some(taxa_str));
                    }
                    break;
                }
                ':' => {
                    // if the current context had a weight
                    if newick_string[str_ptr] == ':' {
                        str_ptr += 1;
                        while newick_string[str_ptr] != ';'
                            && newick_string[str_ptr] != ','
                            && newick_string[str_ptr] != ')'
                            && newick_string[str_ptr] != '('
                        {
                            decimal_str.push(newick_string[str_ptr]);
                            str_ptr += 1;
                        }
                    }
                }
                _ => {
                    // push taxa characters into taxa string
                    while newick_string[str_ptr] != ':'
                        && newick_string[str_ptr] != ')'
                        && newick_string[str_ptr] != ','
                        && newick_string[str_ptr] != '('
                        && newick_string[str_ptr] != ';'
                    {
                        taxa_str.push(newick_string[str_ptr]);
                        str_ptr += 1;
                    }
                }
            }
        }
        Ok(tree)
    }

    fn subtree_to_newick(&self, node_id: TreeNodeID<Self>) -> impl std::fmt::Display {
        dbg!(node_id);
        let node = self.get_node(node_id).unwrap();
        let mut tmp = String::new();
        if node.get_children().len() != 0 {
            if node.get_children().len() > 1 {
                tmp.push('(');
            }
            for child_id in node.get_children() {
                let child_str = format!("{},", self.subtree_to_newick(child_id));
                tmp.push_str(&child_str);
            }
            tmp.pop();
            if node.get_children().collect_vec().len() > 1 {
                tmp.push(')');
            }
        }
        tmp.push_str(node.get_taxa().unwrap_or(&"".to_string()));
        if let Some(w) = node.get_weight() {
            tmp.push(':');
            tmp.push_str(&w.to_string());
        }
        tmp
    }
}

impl Nexus for SimpleRootedTree {}

impl SPR for SimpleRootedTree {
    fn graft(&mut self, tree: Self, edge: (TreeNodeID<Self>, TreeNodeID<Self>)) -> Result<(), ()> {
        let new_node = self.next_node();
        let new_node_id = dbg!(new_node.get_id());
        for node in tree.dfs(tree.get_root_id()) {
            self.set_node(node.clone());
        }
        self.split_edge(edge, new_node);
        self.set_child(dbg!(new_node_id), dbg!(tree.get_root_id()));
        Ok(())
    }
    fn prune(&mut self, node_id: TreeNodeID<Self>) -> Result<Self, ()> {
        let mut pruned_tree = SimpleRootedTree::new(node_id);
        let p_id = self.get_node_parent_id(node_id).unwrap();
        self.get_node_mut(p_id).unwrap().remove_child(&node_id);
        pruned_tree
            .get_node_mut(pruned_tree.get_root_id())
            .unwrap()
            .add_children(self.get_node(node_id).unwrap().get_children());
        let dfs = self.dfs(node_id).collect_vec();
        for node in dfs {
            // self.nodes.remove(node.get_id());
            pruned_tree.set_node(node.clone());
        }
        Ok(pruned_tree)
    }
}

impl Balance for SimpleRootedTree {
    fn balance_subtree(&mut self) -> Result<(), ()> {
        assert!(
            self.get_cluster(self.get_root_id()).collect_vec().len() == 4,
            "Quartets have 4 leaves!"
        );
        assert!(self.is_binary(), "Cannot balance non-binary tree!");
        let root_children = self.get_node_children(self.get_root_id()).collect_vec();
        let (child1, child2) = (root_children[0].get_id(), root_children[1].get_id());
        let next_id = self.next_id();
        let split_id = self.next_id() + 1;
        match dbg!((
            (self.get_node(child1).unwrap().is_leaf()),
            (self.get_node(child2).unwrap().is_leaf())
        )) {
            (false, false) => {}
            (true, false) => {
                let mut leaf_node = self.remove_node(child1).unwrap();
                leaf_node.set_id(next_id);
                let other_leaf_id = &self
                    .get_node_children(child2)
                    .filter(|node| node.is_leaf())
                    .collect_vec()[0]
                    .get_id();
                self.split_edge((child2, *other_leaf_id), Node::new(split_id));
                self.add_child(dbg!(split_id), leaf_node);
            }
            (false, true) => {
                let mut leaf_node = self.remove_node(child2).unwrap();
                leaf_node.set_id(next_id);
                let other_leaf_id = &self
                    .get_node_children(child1)
                    .filter(|node| node.is_leaf())
                    .collect_vec()[0]
                    .get_id();
                self.split_edge((child1, *other_leaf_id), Node::new(split_id));
                self.add_child(split_id, leaf_node);
            }
            _ => {}
        }
        self.clean();
        Ok(())
    }
}

impl CopheneticDistance for SimpleRootedTree {}
