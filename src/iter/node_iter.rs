#[cfg(feature = "non_crypto_hash")]
use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::{HashMap, HashSet};

use itertools::Itertools;
use std::collections::VecDeque;
use vers_vecs::BitVec;

use crate::{
    iter::lca::LcaOracle,
    node::simple_rnode::RootedTreeNode,
    tree::simple_rtree::{RootedTree, TreeNodeID},
};

/// Trait describing depth-first iteration of nodes in a tree
pub trait DFS
where
    Self: RootedTree + Sized,
{
    /// Returns an iterator of immutable reference of nodes in a tree in postfix order
    fn postord_nodes(&self, start_node: TreeNodeID<Self>) -> impl Iterator<Item = &Self::Node>;

    /// Returns an iterator of NodeID's in a tree in postfix order
    fn postord_ids(&self, start_node: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>>;

    /// Returns a DFS iterator of immutable node references a tree
    fn dfs(&self, start_node_id: TreeNodeID<Self>) -> impl ExactSizeIterator<Item = &Self::Node> {
        // A tree reaches every node exactly once through the parent/child
        // structure, so no visited-set is needed to guard against revisits.
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut out_vec = vec![];
        while let Some(x) = stack.pop_front() {
            out_vec.push(x);
            for &child_id in x.get_children().iter().rev() {
                stack.push_front(self.get_node(child_id).unwrap());
            }
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
    fn bfs_nodes(&self, start_node_id: TreeNodeID<Self>) -> impl Iterator<Item = &Self::Node>;

    /// Returns an iterator of NodeID's in a tree in postfix order
    fn bfs_ids(&self, start_node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>>;
}

/// Trait describing breadth-first iteration of nodes in a tree
pub trait PreOrder
where
    Self: RootedTree + Sized,
{
    /// Returns an iterator of immutable reference of nodes in a tree in prefix order
    fn preord_nodes(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &Self::Node> {
        // A tree reaches every node exactly once, so the visited-set the
        // previous body carried was pure overhead.
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut out_vec = vec![];
        while let Some(x) = stack.pop_front() {
            out_vec.push(x);
            for &child_id in x.get_children().iter().rev() {
                stack.push_front(self.get_node(child_id).unwrap());
            }
        }
        out_vec.into_iter()
    }

    /// Returns an iterator of NodeID's in a tree in prefix order
    fn preord_ids(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        // Reversing the children slice in place avoids both the visited-set and
        // the per-node `collect_vec` the previous body allocated.
        let mut stack = VecDeque::from([start_node_id]);
        let mut out_vec = vec![];
        while let Some(x) = stack.pop_front() {
            out_vec.push(x);
            for &child_id in self.get_node(x).unwrap().get_children().iter().rev() {
                stack.push_front(child_id);
            }
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
    fn root_to_node(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &Self::Node> {
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
    fn node_to_root(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &Self::Node> {
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
        // Count edges by walking to the root; the previous body materialised
        // the whole path in a `VecDeque` just to read its length.
        RootedTree::get_node_depth(self, node_id)
    }
}

/// Trait describing an Euler Tour of a tree
pub trait EulerWalk
where
    Self: RootedTree + Sized,
{
    /// Returns euler tour of tree as iterator of immutable references to nodes
    fn euler_walk_nodes(
        &self,
        start_node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &Self::Node> {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut visited: HashSet<TreeNodeID<Self>> = HashSet::default();
        let mut out_vec = vec![];
        while let Some(node) = stack.pop_front() {
            let id = node.get_id();
            if !visited.insert(id) {
                if let Some(parent_id) = node.get_parent() {
                    out_vec.push(self.get_node(parent_id).unwrap())
                }
            } else {
                out_vec.push(node);
                stack.push_front(node);
                for &child_id in node.get_children().iter().rev() {
                    stack.push_front(self.get_node(child_id).unwrap())
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
        let mut visited: HashSet<TreeNodeID<Self>> = HashSet::default();
        let mut out_vec = vec![];
        while let Some(node_id) = stack.pop_front() {
            if !visited.insert(node_id) {
                if let Some(parent_id) = self.get_node_parent_id(node_id) {
                    out_vec.push(parent_id)
                }
            } else {
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
        out_vec.into_iter()
    }

    /// Builds a constant-time LCA oracle borrowing this tree immutably.
    ///
    /// The returned [`LcaOracle`] holds a shared borrow of `self`, so the tree
    /// cannot be mutated while it is alive. Build one, run every query against
    /// it, then drop it before mutating the tree again.
    fn lca(&self) -> LcaOracle<'_, Self> {
        LcaOracle::build(self)
    }

    /// Lowest common ancestor of a slice of nodes, by NodeID.
    ///
    /// Fallback for one-off queries: it builds a throwaway [`LcaOracle`] for a
    /// single lookup. Callers that query repeatedly should build one oracle
    /// with [`Self::lca`] and reuse it.
    ///
    /// # Panics
    ///
    /// Panics if `node_id_vec` is empty or contains an id that is not a node of
    /// this tree.
    fn get_lca_id(&self, node_id_vec: &[TreeNodeID<Self>]) -> TreeNodeID<Self> {
        self.lca().get_lca_id(node_id_vec)
    }

    /// Lowest common ancestor of a slice of nodes, by immutable reference.
    fn get_lca<'a>(&'a self, node_id_vec: &[TreeNodeID<Self>]) -> &'a Self::Node {
        self.get_node(self.get_lca_id(node_id_vec)).unwrap()
    }
}

/// Trait describing iteration of clusters and bipartitions in a tree.
pub trait Clusters: DFS + BFS + Sized {
    /// Returns cluster of a node in a rooted tree (smallest cluster in an unrooted tree) as iterator of immutable reference to a node
    fn get_cluster(&self, node_id: TreeNodeID<Self>) -> impl Iterator<Item = &Self::Node> {
        // `dfs` already materialises its walk, so filtering it lazily avoids a
        // second `Vec` that only existed to promise `ExactSizeIterator`.
        self.dfs(node_id).filter(|x| x.is_leaf())
    }

    /// Returns cluster of a node in a rooted tree (smallest cluster in an unrooted tree) as iterator of NodeID's
    fn get_cluster_ids(&self, node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>> {
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
        self.get_cluster_ids(node_id).count()
    }

    /// Returns bipartition of an edge in a tree as iterator of immutable reference to a node
    fn get_bipartition(
        &self,
        edge: (TreeNodeID<Self>, TreeNodeID<Self>),
    ) -> (
        impl Iterator<Item = &Self::Node>,
        impl Iterator<Item = &Self::Node>,
    ) {
        let c2 = self.get_cluster(edge.1);
        // Hash the opposite cluster so membership is O(1); the previous body
        // scanned `c2_ids` linearly for every element of `c1`. Both sides are
        // now lazy -- no `Vec` just to hand back an `ExactSizeIterator`.
        let c2_ids: HashSet<TreeNodeID<Self>> = self.get_cluster_ids(edge.1).collect();
        let c1 = self
            .get_cluster(edge.0)
            .filter(move |x| !c2_ids.contains(&x.get_id()));
        (c1, c2)
    }

    /// Returns bipartition of an edge in a tree as iterator of NodeID's
    fn get_bipartition_ids(
        &self,
        edge: (TreeNodeID<Self>, TreeNodeID<Self>),
    ) -> (
        impl Iterator<Item = TreeNodeID<Self>>,
        impl Iterator<Item = TreeNodeID<Self>>,
    ) {
        let c2 = self.get_cluster_ids(edge.1);
        // O(1) membership instead of a linear scan per element of `c1`.
        let c2_ids: HashSet<TreeNodeID<Self>> = self.get_cluster_ids(edge.1).collect();
        let c1 = self
            .get_cluster_ids(edge.0)
            .filter(move |x| !c2_ids.contains(x));
        (c1, c2)
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
        let leaf_ids_rev: Vec<TreeNodeID<Self>> = leaf_ids.keys().copied().collect();
        let num_leaves = leaf_ids.len();
        let mut bps: HashMap<TreeNodeID<Self>, BitVec> = vec![].into_iter().collect();
        for n_id in self.postord_ids(self.get_root_id()) {
            let mut bp = BitVec::from_zeros(num_leaves);
            match self.is_leaf(n_id) {
                true => {
                    bp.flip_bit(*leaf_ids.get(&n_id).unwrap());
                    bps.insert(n_id, bp.clone());
                }
                false => {
                    if n_id == self.get_root_id() {
                        continue;
                    }
                    self.get_node_children_ids(n_id)
                        .map(|x| bps.get(&x).unwrap())
                        .for_each(|x| {
                            let _ = bp.apply_mask_or(x);
                        });
                    if self.get_node_parent_id(n_id) != Some(self.get_root_id()) {
                        bps.insert(n_id, bp);
                    }
                }
            };
        }

        bps.into_values().map(move |bit_bp| {
            let mut bp1 = Vec::with_capacity(leaf_ids.len());
            let mut bp2 = Vec::with_capacity(leaf_ids.len());
            for (idx, bit) in leaf_ids_rev.iter().enumerate().take(bit_bp.len()) {
                match bit_bp.is_bit_set(idx).unwrap() {
                    true => {
                        bp1.push(bit.to_owned());
                    }
                    false => {
                        bp2.push(bit.to_owned());
                    }
                }
            }
            (bp1.into_iter(), bp2.into_iter())
        })
    }

    /// Returns median NodeID of a set of leaves in a tree.
    fn get_median_node_id_for_leaves(
        &self,
        taxa_set: impl Iterator<Item = TreeNodeID<Self>>,
    ) -> TreeNodeID<Self> {
        let mut cluster_sizes: HashMap<TreeNodeID<Self>, usize> = vec![].into_iter().collect();
        let mut median_node_id: TreeNodeID<Self> = self.get_root_id();
        let leaf_ids: HashSet<TreeNodeID<Self>> = taxa_set.collect();
        for n_id in self.postord_ids(self.get_root_id()) {
            if self.is_leaf(n_id) && leaf_ids.contains(&n_id) {
                cluster_sizes.insert(n_id, 1);
            } else {
                let mut cluster_size = 0;
                for c_id in self.get_node_children_ids(n_id) {
                    cluster_size += cluster_sizes.get(&c_id).unwrap();
                }
                cluster_sizes.insert(n_id, cluster_size);
            }
        }
        loop {
            median_node_id = self
                .get_node_children_ids(median_node_id)
                .max_by(|x, y| {
                    let x_cluster_size = cluster_sizes.get(x).unwrap();
                    let y_cluster_size = cluster_sizes.get(y).unwrap();
                    x_cluster_size.cmp(y_cluster_size)
                })
                .unwrap();
            if cluster_sizes.get(&median_node_id).unwrap() <= &(leaf_ids.len() / 2) {
                break;
            }
        }
        median_node_id
    }

    /// Returns immutable reference to median node of a set of leaves in a tree.
    fn get_median_node_for_leaves(
        &self,
        taxa_set: impl Iterator<Item = TreeNodeID<Self>>,
    ) -> &Self::Node {
        self.get_node(self.get_median_node_id_for_leaves(taxa_set))
            .unwrap()
    }

    /// Returns an immutable reference to median node of all leaves in a tree.
    fn get_median_node(&self) -> &Self::Node {
        let leaves = self.get_leaves().map(|x| x.get_id());
        self.get_median_node_for_leaves(leaves)
    }

    /// Returns median NodeID of all leaves in a tree.
    fn get_median_node_id(&self) -> TreeNodeID<Self> {
        let leaves = self.get_leaf_ids();
        self.get_median_node_id_for_leaves(leaves)
    }
}
