#![allow(clippy::needless_lifetimes)]

#[cfg(feature = "non_crypto_hash")]
use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::{HashMap, HashSet};

use std::{collections::VecDeque, ops::Index};

use itertools::Itertools;

use crate::{
    node::simple_rnode::RootedTreeNode,
    tree::simple_rtree::{RootedTree, TreeNodeID},
};

/// Trait describing depth-first iteration of nodes in a tree
pub trait DFS
where
    Self: RootedTree + Sized,
{
    /// Returns an iterator of immutable reference of nodes in a tree in postfix order
    fn postord_nodes<'a>(
        &'a self,
        start_node: TreeNodeID<Self>,
    ) -> impl Iterator<Item = &'a Self::Node>;

    /// Returns an iterator of NodeID's in a tree in postfix order
    fn postord_ids(&self, start_node: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>>;

    /// Returns a DFS iterator of immutable node references a tree
    fn dfs<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut out_vec = vec![];
        let mut visited = vec![];
        while let Some(x) = stack.pop_front() {
            if !visited.contains(&x.get_id()) {
                visited.push(x.get_id());
                out_vec.push(x);
                for child_id in x.get_children().collect_vec().into_iter().rev() {
                    stack.push_front(self.get_node(child_id).unwrap());
                }
            };
        }
        out_vec.into_iter()
    }
}

/// Trait describing breadth-first iteration of nodes in a tree
pub trait BFS
where
    Self: RootedTree + Sized,
{
    /// Returns an iterator of immutable reference of nodes in a tree in postfix order
    fn bfs_nodes<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl Iterator<Item = &'a Self::Node>;

    /// Returns an iterator of NodeID's in a tree in postfix order
    fn bfs_ids(&self, start_node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>>;
}

/// Trait describing breadth-first iteration of nodes in a tree
pub trait PreOrder
where
    Self: RootedTree + Sized,
{
    /// Returns an iterator of immutable reference of nodes in a tree in prefix order
    fn preord_nodes<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut out_vec = vec![];
        let mut visited = vec![];
        while let Some(x) = stack.pop_front() {
            if !visited.contains(&x.get_id()) {
                visited.push(x.get_id());
                out_vec.push(x);
                for child_id in x.get_children().collect_vec().into_iter().rev() {
                    stack.push_front(self.get_node(child_id).unwrap());
                }
            };
        }
        out_vec.into_iter()
    }

    /// Returns an iterator of NodeID's in a tree in prefix order
    fn preord_ids(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        let mut stack = VecDeque::from([start_node_id]);
        let mut out_vec = vec![];
        let mut visited = vec![];
        while let Some(x) = stack.pop_front() {
            if !visited.contains(&x) {
                visited.push(x);
                out_vec.push(x);
                for child_id in self
                    .get_node_children_ids(x)
                    .collect_vec()
                    .into_iter()
                    .rev()
                {
                    stack.push_front(child_id);
                }
            };
        }
        out_vec.into_iter()
    }
}

/// Trait describing iteration of nodes along a path
pub trait Ancestors
where
    Self: RootedTree + Sized,
{
    /// Returns an iterator of immutable references to nodes in a tree from root to node
    fn root_to_node<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        while let Some(x) = stack.pop_front() {
            stack.push_front(x);
            match x.get_parent() {
                Some(pid) => {
                    stack.push_front(self.get_node(pid).unwrap());
                }
                None => {
                    break;
                }
            }
        }
        stack.into_iter()
    }

    /// Returns an iterator of NodeID's in a tree from root to node
    fn root_to_node_ids(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        let mut stack = VecDeque::from([start_node_id]);
        while let Some(x) = stack.pop_front() {
            stack.push_front(x);
            match self.get_node_parent_id(x) {
                Some(pid) => {
                    stack.push_front(pid);
                }
                None => {
                    break;
                }
            }
        }
        stack.into_iter()
    }

    /// Returns an iterator of immutable references to nodes in a tree from node to root
    fn node_to_root<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        while let Some(x) = stack.pop_front() {
            stack.push_back(x);
            match x.get_parent() {
                Some(pid) => {
                    stack.push_front(self.get_node(pid).unwrap());
                }
                None => {
                    break;
                }
            }
        }
        stack.into_iter()
    }

    /// Returns an iterator of NodeID's in a tree from node to root
    fn node_to_root_ids(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        let mut stack = VecDeque::from([start_node_id]);
        while let Some(x) = stack.pop_front() {
            stack.push_back(x);
            match self.get_node_parent_id(x) {
                Some(pid) => {
                    stack.push_front(pid);
                }
                None => {
                    break;
                }
            }
        }
        stack.into_iter()
    }

    /// Returns depth of a node as number of edges in the path from node to root
    fn depth(&self, node_id: TreeNodeID<Self>) -> usize {
        self.node_to_root_ids(node_id).len() - 1
    }
}

/// Trait describing an Euler Tour of a tree
pub trait EulerWalk
where
    Self: RootedTree + Sized,
{
    /// Precomputes an Euler Tour of a tree
    fn precompute_walk(&mut self);

    /// Returns the euler tour of the tree
    fn get_precomputed_walk(&self) -> Option<&Vec<TreeNodeID<Self>>>;

    /// Precomputes the first-appearance index of nodes in an euler walk
    fn precompute_fai(&mut self);

    /// Returns the first-appearance index of nodes in an euler walk
    fn get_precomputed_fai(&self) -> Option<impl Index<TreeNodeID<Self>, Output = Option<usize>>>;

    /// Precomutes the depth-array of nodes in an euler walk
    fn precompute_da(&mut self);

    /// Returns the depth-array of nodes in an euler walk
    fn get_precomputed_da(&self) -> Option<&Vec<usize>>;

    /// Precomutes the structures required for constant-time lca queries
    fn precompute_constant_time_lca(&mut self) {
        self.precompute_walk();
        self.precompute_da();
        self.precompute_fai();
    }

    /// Returns euler tour of tree as iterator of immutable references to nodes
    fn euler_walk_nodes<'a>(
        &'a self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut visited = vec![];
        let mut out_vec = vec![];
        while let Some(node) = stack.pop_front() {
            match visited.contains(&node.get_id()) {
                true => {
                    if let Some(parent_id) = node.get_parent() {
                        out_vec.push(self.get_node(parent_id).unwrap())
                    }
                }
                false => {
                    visited.push(node.get_id());
                    out_vec.push(node);
                    stack.push_front(node);
                    for child_id in node.get_children().rev() {
                        stack.push_front(self.get_node(child_id).unwrap())
                    }
                }
            }
        }
        out_vec.into_iter()
    }

    /// Returns euler tour of tree as iterator of NodeID's
    fn euler_walk_ids(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        let mut stack = VecDeque::from([start_node_id]);
        let mut visited = vec![];
        let mut out_vec = vec![];
        while let Some(node_id) = stack.pop_front() {
            match visited.contains(&node_id) {
                true => {
                    if let Some(parent_id) = self.get_node_parent_id(node_id) {
                        out_vec.push(parent_id)
                    }
                }
                false => {
                    visited.push(node_id);
                    out_vec.push(node_id);
                    stack.push_front(node_id);
                    for child_id in self
                        .get_node_children_ids(node_id)
                        .collect_vec()
                        .iter()
                        .rev()
                    {
                        stack.push_front(*child_id)
                    }
                }
            }
        }
        out_vec.into_iter()
    }

    /// Returns true if euler tour is precomputed
    fn is_euler_precomputed(&self) -> bool {
        self.get_precomputed_walk().is_some()
    }

    /// Computes and returns the first appearance index of each node in the euler tour
    fn first_appearance(&self) -> impl Index<TreeNodeID<Self>, Output = Option<usize>>;

    /// Returns true if first-appearance index is precomputed
    fn is_fai_precomputed(&self) -> bool {
        self.get_precomputed_fai().is_some()
    }

    /// Returns first-appearance index of tree
    fn get_fa_index(&self, node_id: TreeNodeID<Self>) -> usize {
        match self.get_precomputed_fai() {
            Some(fai) => fai[node_id].unwrap(),
            None => self.first_appearance()[node_id].unwrap(),
        }
    }

    /// depth array for nodes
    fn depth_array(&self) -> Vec<usize> {
        let da = match self.get_precomputed_walk() {
            Some(walk) => walk
                .iter()
                .map(|x| RootedTree::get_node_depth(self, *x))
                .collect_vec(),
            None => self
                .euler_walk_ids(self.get_root_id())
                .map(|x| RootedTree::get_node_depth(self, x))
                .collect_vec(),
        };
        da
    }

    /// Returns true if depth-array is precomputed
    fn is_da_precomputed(&self) -> bool {
        self.get_precomputed_da().is_some()
    }

    /// Returns depth of node using depth-array
    fn get_node_depth(&self, node_id: TreeNodeID<Self>) -> usize {
        match self.get_precomputed_da() {
            Some(da) => da[self.get_fa_index(node_id)],
            None => RootedTree::get_node_depth(self, node_id),
        }
    }

    /// Returns slice of euler tour
    fn get_euler_slice(&self, start: usize, end: usize) -> Vec<TreeNodeID<Self>> {
        match self.get_precomputed_walk() {
            Some(walk) => walk[start..end].to_vec(),
            None => self.euler_walk_ids(self.get_root_id()).collect_vec()[start..end].to_vec(),
        }
    }

    /// Returns slice of depth-array
    fn get_da_slice(&self, start: usize, end: usize) -> Vec<usize> {
        match self.get_precomputed_da() {
            Some(da) => da[start..end].to_vec(),
            None => self.depth_array()[start..end].to_vec(),
        }
    }

    /// Returns NodeID at index pos of euler tour
    fn get_euler_pos(&self, pos: usize) -> TreeNodeID<Self> {
        match self.get_precomputed_walk() {
            Some(walk) => walk[pos],
            None => self.euler_walk_ids(self.get_root_id()).nth(pos).unwrap(),
        }
    }

    /// Constant time lca query of slice of nodes by NodeID
    fn get_lca_id(&self, node_id_vec: &[TreeNodeID<Self>]) -> TreeNodeID<Self> {
        if node_id_vec.len() == 1 {
            return node_id_vec[0];
        }
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

        let depth_subarray_min_value = self
            .get_da_slice(min_pos, max_pos)
            .into_iter()
            .min()
            .unwrap();
        let depth_subarray_min_pos = self
            .get_da_slice(min_pos, max_pos)
            .into_iter()
            .position(|x| x == depth_subarray_min_value)
            .unwrap();
        self.get_euler_pos(min_pos + depth_subarray_min_pos)
    }

    /// Constant time lca query of slice of nodes by immutable reference to node.
    fn get_lca<'a>(&'a self, node_id_vec: &[TreeNodeID<Self>]) -> &'a Self::Node {
        self.get_node(self.get_lca_id(node_id_vec)).unwrap()
    }
}

/// Trait describing iteration of clusters and bipartitions in a tree.
pub trait Clusters: DFS + BFS + Sized {
    /// Returns cluster of a node in a rooted tree (smallest cluster in an unrooted tree) as iterator of immutable reference to a node
    fn get_cluster<'a>(
        &'a self,
        node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        self.dfs(node_id)
            .filter(|x| x.is_leaf())
            .collect_vec()
            .into_iter()
    }

    /// Returns cluster of a node in a rooted tree (smallest cluster in an unrooted tree) as iterator of NodeID's
    fn get_cluster_ids(
        &self,
        node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        self.get_cluster(node_id).map(move |x| x.get_id())
    }

    /// Returns all clusters of a tree as iterator of NodeID's
    fn get_clusters_ids(
        &self,
    ) -> impl ExactSizeIterator<
        Item = (
            TreeNodeID<Self>,
            impl ExactSizeIterator<Item = TreeNodeID<Self>>,
        ),
    > {
        let mut clusters: HashMap<TreeNodeID<Self>, Vec<TreeNodeID<Self>>> =
            vec![].into_iter().collect();
        for n_id in self.postord_ids(self.get_root_id()) {
            match self.is_leaf(n_id) {
                true => {
                    clusters.insert(n_id, vec![n_id]);
                }
                false => {
                    let node_cluster = self
                        .get_node_children_ids(n_id)
                        .flat_map(|x| clusters.get(&x).cloned().unwrap())
                        .collect_vec();
                    clusters.insert(n_id, node_cluster);
                }
            };
        }

        clusters
            .into_iter()
            .map(|(n_id, cluster)| (n_id, cluster.into_iter()))
    }

    /// Returns size of a cluster of nodes
    fn get_cluster_size(&self, node_id: TreeNodeID<Self>) -> usize {
        self.get_cluster_ids(node_id).len()
    }

    /// Returns bipartition of an edge in a tree as iterator of immutable reference to a node
    fn get_bipartition<'a>(
        &'a self,
        edge: (TreeNodeID<Self>, TreeNodeID<Self>),
    ) -> (
        impl ExactSizeIterator<Item = &'a Self::Node>,
        impl ExactSizeIterator<Item = &'a Self::Node>,
    ) {
        let c2 = self.get_cluster(edge.1);
        let c2_ids = self.get_cluster_ids(edge.1).collect_vec();
        let c1 = self
            .get_cluster(edge.0)
            .filter(|x| !c2_ids.contains(&x.get_id()))
            .collect_vec();
        (c1.into_iter(), c2.into_iter())
    }

    /// Returns bipartition of an edge in a tree as iterator of NodeID's
    fn get_bipartition_ids(
        &self,
        edge: (TreeNodeID<Self>, TreeNodeID<Self>),
    ) -> (
        impl ExactSizeIterator<Item = TreeNodeID<Self>>,
        impl ExactSizeIterator<Item = TreeNodeID<Self>>,
    ) {
        let c2 = self.get_cluster_ids(edge.1);
        let c2_ids = self.get_cluster_ids(edge.1).collect_vec();
        let c1 = self
            .get_cluster_ids(edge.0)
            .filter(|x| !c2_ids.contains(x))
            .collect_vec();
        (c1.into_iter(), c2.into_iter())
    }

    /// Returns all bipartitions of a tree as iterator of NodeID's
    fn get_bipartitions_ids(
        &self,
    ) -> impl ExactSizeIterator<
        Item = (
            impl ExactSizeIterator<Item = TreeNodeID<Self>>,
            impl ExactSizeIterator<Item = TreeNodeID<Self>>,
        ),
    > {
        let leaf_ids: HashMap<TreeNodeID<Self>, usize> = self
            .get_leaf_ids()
            .enumerate()
            .map(|(idx, id)| (id, idx))
            .collect();
        let leaf_ids_rev: HashMap<usize, TreeNodeID<Self>> =
            leaf_ids.iter().map(|(id, idx)| (*idx, *id)).collect();
        let mut bps: HashMap<TreeNodeID<Self>, Vec<bool>> = vec![].into_iter().collect();
        let mut skip_clusters: HashSet<Vec<bool>> = vec![].into_iter().collect();
        let mut keep_nodes: Vec<TreeNodeID<Self>> = vec![];
        for n_id in self.postord_ids(self.get_root_id()) {
            match self.is_leaf(n_id) {
                true => {
                    let bp: (HashSet<TreeNodeID<Self>>, HashSet<TreeNodeID<Self>>) = (
                        vec![n_id].into_iter().collect(),
                        leaf_ids.keys().filter(|x| **x != n_id).copied().collect(),
                    );
                    let mut binary_str = vec![false; leaf_ids.len()];
                    for i in bp.0.iter() {
                        binary_str[*leaf_ids.get(i).unwrap()] = true;
                    }
                    if !skip_clusters.contains(&binary_str)
                        && !skip_clusters.contains(&binary_str.iter().map(|x| !x).collect_vec())
                    {
                        skip_clusters.insert(binary_str.iter().map(|x| !x).collect_vec());
                        keep_nodes.push(n_id);
                    }
                    bps.insert(n_id, binary_str);
                }
                false => {
                    if n_id == self.get_root_id() {
                        continue;
                    }
                    let children_clusters = self
                        .get_node_children_ids(n_id)
                        .map(|x| bps.get(&x).unwrap())
                        .collect_vec();
                    let binary_str = (0..leaf_ids.len())
                        .map(|x| {
                            let values = children_clusters
                                .iter()
                                .map(|y| y[x] as usize)
                                .collect_vec();
                            values.iter().sum::<usize>() > 0
                        })
                        .collect_vec();
                    if !skip_clusters.contains(&binary_str)
                        && !skip_clusters.contains(&binary_str.iter().map(|x| !x).collect_vec())
                    {
                        skip_clusters.insert(binary_str.iter().map(|x| !x).collect_vec());
                        keep_nodes.push(n_id);
                    }
                    bps.insert(n_id, binary_str);
                }
            };
        }

        return keep_nodes.into_iter().map(move |n_id| {
            let mut bp1 = vec![];
            let mut bp2 = vec![];
            bps.get(&n_id)
                .unwrap()
                .iter()
                .enumerate()
                .for_each(|(idx, x)| match x {
                    true => {
                        bp1.push(leaf_ids_rev.get(&idx).unwrap().to_owned());
                    }
                    false => {
                        bp2.push(leaf_ids_rev.get(&idx).unwrap().to_owned());
                    }
                });

            (bp1.into_iter(), bp2.into_iter())
        });
    }

    /// Returns median NodeID of a set of leaves in a tree.
    fn get_median_node_id_for_leaves(
        &self,
        taxa_set: impl ExactSizeIterator<Item = TreeNodeID<Self>>,
    ) -> TreeNodeID<Self> {
        let mut median_node_id: TreeNodeID<Self> = self.get_root_id();
        let leaf_ids = taxa_set.collect_vec();
        loop {
            median_node_id = self
                .get_node_children_ids(median_node_id)
                // .filter(|x| !self.is_leaf(x))
                .max_by(|x, y| {
                    let x_cluster = self
                        .get_cluster_ids(*x)
                        .filter(|nids| leaf_ids.contains(nids))
                        .collect_vec();
                    let y_cluster = self
                        .get_cluster_ids(*y)
                        .filter(|nids| leaf_ids.contains(nids))
                        .collect_vec();
                    x_cluster.len().cmp(&y_cluster.len())
                })
                .unwrap();
            let num_leaves = self
                .get_cluster_ids(median_node_id)
                .filter(|nids| leaf_ids.contains(nids))
                .collect_vec()
                .len();
            if (num_leaves) <= (leaf_ids.len() / 2) {
                break;
            }
        }
        median_node_id
    }

    /// Returns immutable reference to median node of a set of leaves in a tree.
    fn get_median_node_for_leaves<'a>(
        &'a self,
        taxa_set: impl ExactSizeIterator<Item = TreeNodeID<Self>>,
    ) -> &'a Self::Node {
        self.get_node(self.get_median_node_id_for_leaves(taxa_set))
            .unwrap()
    }

    /// Returns an immutable reference to median node of all leaves in a tree.
    fn get_median_node<'a>(&'a self) -> &'a Self::Node {
        let leaves = self.get_leaves().map(|x| x.get_id());
        self.get_median_node_for_leaves(leaves)
    }

    /// Returns median NodeID of all leaves in a tree.
    fn get_median_node_id(&self) -> TreeNodeID<Self> {
        let leaves = self.get_leaf_ids();
        self.get_median_node_id_for_leaves(leaves)
    }
}
