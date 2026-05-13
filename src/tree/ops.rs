use std::{
    fmt::{Debug, Display},
    hash::Hash,
};

use crate::node::simple_rnode::{NodeTaxa, RootedMetaNode};
use crate::prelude::{Clusters, EulerWalk, PreOrder, RootedMetaTree, RootedTree, DFS};
use crate::{
    iter::node_iter::Ancestors,
    node::simple_rnode::RootedTreeNode,
    tree::simple_rtree::{TreeNodeID, TreeNodeMeta},
};

#[cfg(feature = "non_crypto_hash")]
use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::{HashMap, HashSet};

/// A trait describing subtree-prune-regraft operations
pub trait SPR: RootedTree + DFS + Sized {
    /// Attaches input tree to self by spliting an edge
    fn graft(&mut self, tree: Self, edge: (TreeNodeID<Self>, TreeNodeID<Self>)) -> Result<(), ()>;

    /// Returns subtree starting at given node, while corresponding nodes from self.
    fn prune(&mut self, node_id: TreeNodeID<Self>) -> Result<Self, ()>;

    /// SPR function
    fn spr(
        &mut self,
        edge1: (TreeNodeID<Self>, TreeNodeID<Self>),
        edge2: (TreeNodeID<Self>, TreeNodeID<Self>),
    ) -> Result<(), ()> {
        let pruned_tree = SPR::prune(self, edge1.1)?;
        SPR::graft(self, pruned_tree, edge2)
    }
}

/// A trait describing Nearest Neighbour interchange operations
pub trait NNI
where
    Self: RootedTree + Sized,
{
    /// Performs an NNI operation
    fn nni(&mut self, node_id: TreeNodeID<Self>, left_ch: bool) -> Result<(), ()>;
}

/// A trait describing rerooting a tree
pub trait Reroot<'a>
where
    Self: RootedTree + Sized,
{
    /// Reroots tree at node. **Note: this changes the degree of a node**
    fn reroot_at_node(&mut self, node_id: TreeNodeID<Self>) -> Result<(), ()>;
    /// Reroots tree at a split node.
    fn reroot_at_edge(&mut self, edge: (TreeNodeID<Self>, TreeNodeID<Self>)) -> Result<(), ()>;
}

/// A trait describing balancing a binary tree
pub trait Balance: Clusters + SPR + Sized
where
    TreeNodeID<Self>: Display + Debug + Hash + Clone + Ord,
{
    /// Balances a binary tree
    fn balance_subtree(&mut self) -> Result<(), ()>;
}

/// A trait describing subtree queries of a tree
pub trait Subtree: Ancestors + DFS + Sized
where
    TreeNodeID<Self>: Display + Debug + Hash + Clone + Ord,
{
    /// Returns a subtree consisting of only provided nodes
    fn induce_tree(
        &self,
        node_id_list: impl IntoIterator<
            Item = TreeNodeID<Self>,
            IntoIter = impl ExactSizeIterator<Item = TreeNodeID<Self>>,
        >,
    ) -> Result<Self, ()> {
        let mut subtree = Self::new();
        subtree.set_root(self.get_root_id());
        subtree.set_node(self.get_root().clone());
        for node_id in node_id_list.into_iter() {
            let ancestors = self.root_to_node(node_id).cloned();
            subtree.set_nodes(ancestors);
        }
        subtree.clean();
        Ok(subtree)
    }

    /// Returns subtree starting at provided node.
    fn subtree(&self, node_id: TreeNodeID<Self>) -> Result<Self, ()> {
        let mut subtree = Self::new();
        subtree.set_root(node_id);
        let dfs = self.dfs(node_id).cloned();
        subtree.set_nodes(dfs);
        subtree.get_node_mut(node_id).unwrap().set_parent(None);
        Ok(subtree)
    }
}

/// A trait describing tree contraction operations
pub trait ContractTree: EulerWalk + DFS {
    /// Contracts tree that from post_ord node_id iterator.
    fn contracted_tree_nodes_from_iter(
        &self,
        new_tree_root_id: TreeNodeID<Self>,
        leaf_ids: &[TreeNodeID<Self>],
        node_iter: impl Iterator<Item = TreeNodeID<Self>>,
    ) -> impl Iterator<Item = Self::Node> {
        let mut node_map: HashMap<TreeNodeID<Self>, Self::Node> =
            HashMap::from_iter(vec![(new_tree_root_id, self.get_lca(leaf_ids).clone())]);
        let mut remove_list = HashSet::from_iter(vec![]);
        node_iter
            .map(|x| self.get_node(x).cloned().unwrap())
            .for_each(|mut node| {
                match node.is_leaf() {
                    true => {
                        if leaf_ids.contains(&node.get_id()) {
                            node_map.insert(node.get_id(), node);
                        }
                    }
                    false => {
                        let node_children_ids = node.get_children().to_vec();
                        for child_id in &node_children_ids {
                            match node_map.contains_key(child_id) {
                                true => {}
                                false => node.remove_child(child_id),
                            }
                        }
                        let node_children_ids = node.get_children().to_vec();
                        match node_children_ids.len() {
                            0 => {}
                            1 => {
                                // the node is a unifurcation
                                // node should be added to both node_map and remove_list
                                // if child of node is already in remove list, attach node children to node first
                                let child_node_id = node_children_ids[0];

                                if remove_list.contains(&child_node_id) {
                                    node.remove_child(&child_node_id);
                                    let grandchildren_ids = node_map
                                        .get(&child_node_id)
                                        .unwrap()
                                        .get_children()
                                        .to_vec();
                                    for grandchild_id in grandchildren_ids {
                                        node_map
                                            .get_mut(&grandchild_id)
                                            .unwrap()
                                            .set_parent(Some(node.get_id()));
                                        node.add_child(grandchild_id);
                                    }
                                }
                                remove_list.insert(node.get_id());
                                node_map.insert(node.get_id(), node);
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
                                        let node_grandchildren = node_map
                                            .get(&chid)
                                            .unwrap()
                                            .get_children()
                                            .to_vec();
                                        for grandchild in node_grandchildren {
                                            node.add_child(grandchild);
                                            node_map
                                                .get_mut(&grandchild)
                                                .unwrap()
                                                .set_parent(Some(node.get_id()))
                                        }
                                    }
                                });
                                if node.get_id() == new_tree_root_id {
                                    node.set_parent(None);
                                }
                                node_map.insert(node.get_id(), node);
                            }
                        };
                    }
                }
            });
        remove_list.into_iter().for_each(|x| {
            node_map.remove(&x);
        });
        node_map.into_values()
    }

    /// Returns a deep copy of the nodes in the contracted tree
    fn contracted_tree_nodes(
        &self,
        leaf_ids: &[TreeNodeID<Self>],
    ) -> impl Iterator<Item = Self::Node> {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let node_postord_iter = self.postord_nodes(new_tree_root_id);
        let mut node_map: HashMap<TreeNodeID<Self>, Self::Node> =
            HashMap::from_iter(vec![(new_tree_root_id, self.get_lca(leaf_ids).clone())]);
        let leaf_ids: HashSet<&TreeNodeID<Self>> = leaf_ids.iter().collect();
        let mut remove_list = vec![];
        node_postord_iter.for_each(|orig_node| {
            let mut node = orig_node.clone();
            match node.is_leaf() {
                true => {
                    if leaf_ids.contains(&node.get_id()) {
                        node_map.insert(node.get_id(), node);
                    }
                }
                false => {
                    let node_children_ids = node.get_children().to_vec();
                    for child_id in &node_children_ids {
                        match node_map.contains_key(child_id) {
                            true => {}
                            false => node.remove_child(child_id),
                        }
                    }
                    let node_children_ids = node.get_children().to_vec();
                    match node_children_ids.len() {
                        0 => {}
                        1 => {
                            // the node is a unifurcation
                            // node should be added to both node_map and remove_list
                            // if child of node is already in remove list, attach node children to node first
                            let child_node_id = node_children_ids[0];

                            if remove_list.contains(&child_node_id) {
                                node.remove_child(&child_node_id);
                                let grandchildren_ids = node_map
                                    .get(&child_node_id)
                                    .unwrap()
                                    .get_children()
                                    .to_vec();
                                for grandchild_id in grandchildren_ids {
                                    node_map
                                        .get_mut(&grandchild_id)
                                        .unwrap()
                                        .set_parent(Some(node.get_id()));
                                    node.add_child(grandchild_id);
                                }
                            }
                            remove_list.push(node.get_id());
                            node_map.insert(node.get_id(), node);
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
                                    let node_grandchildren =
                                        node_map.get(&chid).unwrap().get_children().to_vec();
                                    for grandchild in node_grandchildren {
                                        node.add_child(grandchild);
                                        node_map
                                            .get_mut(&grandchild)
                                            .unwrap()
                                            .set_parent(Some(node.get_id()))
                                    }
                                }
                            });
                            if node.get_id() == new_tree_root_id {
                                node.set_parent(None);
                            }
                            node_map.insert(node.get_id(), node.clone());
                        }
                    };
                }
            }
        });
        remove_list.into_iter().for_each(|x| {
            node_map.remove(&x);
        });
        node_map.into_values()
    }

    /// Returns a contracted tree from slice containing NodeID's
    fn contract_tree(&self, leaf_ids: &[TreeNodeID<Self>]) -> Result<Self, ()>;

    /// Returns a contracted tree from an iterator containing NodeID's
    fn contract_tree_from_iter(
        &self,
        leaf_ids: &[TreeNodeID<Self>],
        node_iter: impl Iterator<Item = TreeNodeID<Self>>,
    ) -> Result<Self, ()>;
}

/// A struct representing an Ordered Leaf Array tree
#[derive(Clone)]
pub struct OLATree<T: NodeTaxa> {
    /// Taxa labels in leaf ordering σ
    pub taxa: Vec<T>,
    /// OLA indices: non-negative values are leaf indices, negative values are internal node indices
    pub indices: Vec<i64>,
}

impl<T: NodeTaxa> Default for OLATree<T> {
    fn default() -> Self {
        OLATree {
            taxa: Vec::new(),
            indices: Vec::new(),
        }
    }
}

/// Returns the child of `ancestor` on the path toward `descendant`.
fn child_of<Tr>(tree: &Tr, ancestor: TreeNodeID<Tr>, descendant: TreeNodeID<Tr>) -> TreeNodeID<Tr>
where
    Tr: RootedTree,
    Tr::Node: RootedTreeNode,
{
    let mut current = descendant;
    loop {
        let parent = tree.get_node_parent_id(current).unwrap();
        if parent == ancestor {
            return current;
        }
        current = parent;
    }
}

/// A trait for converting trees to and from an Ordered Leaf Array representation
pub trait OLA: RootedMetaTree + EulerWalk + PreOrder
where
    Self::Node: RootedMetaNode,
{
    /// Constructs a tree from an OLATree representation
    fn from_vec(ola: OLATree<TreeNodeMeta<Self>>) -> Self;

    /// Converts the tree into an OLATree representation.
    ///
    /// Leaves are ordered by pre-order DFS traversal to establish the leaf ordering σ.
    /// Each index entry is either a leaf index (≥ 0) or an internal node index (< 0).
    fn to_vec(&self) -> OLATree<TreeNodeMeta<Self>> {
        // Step 1: collect leaves in pre-order to fix leaf ordering σ
        let leaf_ids: Vec<TreeNodeID<Self>> = self
            .preord_ids(self.get_root_id())
            .filter(|id| self.is_leaf(*id))
            .collect();

        let n = leaf_ids.len();
        if n <= 1 {
            return OLATree::default();
        }

        let mut ola_indices: Vec<i64> = Vec::with_capacity(n - 1);

        // Step 2: for each leaf l_i (i >= 1), find its sibling in the restricted tree T^i
        for i in 1..n {
            let li = leaf_ids[i];

            // The parent of l_i in T^i is the LCA(l_i, l_j) with the greatest depth over all j < i
            let p_id = (0..i)
                .map(|j| self.get_lca_id(&[li, leaf_ids[j]]))
                .max_by_key(|&lca| EulerWalk::get_node_depth(self, lca))
                .unwrap();

            // Sibling's leaves in T^{i-1}: those l_j (j < i) whose LCA with l_i is exactly p_id,
            // meaning they live on the opposite side of p_id from l_i
            let sibling_indices: Vec<usize> = (0..i)
                .filter(|&j| self.get_lca_id(&[li, leaf_ids[j]]) == p_id)
                .collect();

            let entry = if sibling_indices.len() == 1 {
                // Sibling is a single leaf: entry = its index in σ (non-negative)
                sibling_indices[0] as i64
            } else {
                // Sibling is an internal node.
                // index(v) = -max(μ(c1), μ(c2)), where μ(c) = min leaf index in child c's subtree.
                // The sibling node in T^i is the LCA of all sibling leaves in the original tree.
                // Split sibling_indices by which child of sib_id each leaf descends through.
                let sib_leaf_ids: Vec<TreeNodeID<Self>> =
                    sibling_indices.iter().map(|&j| leaf_ids[j]).collect();
                let sib_id = self.get_lca_id(&sib_leaf_ids);
                let mut child_min: HashMap<TreeNodeID<Self>, usize> = HashMap::default();
                for &j in &sibling_indices {
                    let child = child_of(self, sib_id, leaf_ids[j]);
                    let e = child_min.entry(child).or_insert(j);
                    if j < *e {
                        *e = j;
                    }
                }
                let mu_max = *child_min.values().max().unwrap();
                -(mu_max as i64)
            };

            ola_indices.push(entry);
        }

        // Step 3: collect taxa labels in leaf ordering σ
        let taxa: Vec<TreeNodeMeta<Self>> = leaf_ids
            .iter()
            .map(|&id| self.get_node_taxa_cloned(id).unwrap())
            .collect();

        OLATree { taxa, indices: ola_indices }
    }
}
