use fxhash::FxHashMap as HashMap;
use anyhow::Result;
use itertools::Itertools;
use std::{
    fmt::{Debug, Display},
    hash::Hash,
};

use super::{Clusters, EulerWalk, RootedTree, DFS};
use crate::{
    iter::node_iter::Ancestors,
    node::simple_rnode::RootedTreeNode,
    tree::simple_rtree::TreeNodeID,
};

/// A trait describing subtree-prune-regraft operations
pub trait SPR<'a>: RootedTree<'a> + DFS<'a> + Sized {
    /// Attaches input tree to self by spliting an edge
    fn graft(&mut self, tree: Self, edge: (TreeNodeID<'a, Self>, TreeNodeID<'a, Self>))->Result<()>;

    /// Returns subtree starting at given node, while corresponding nodes from self.
    fn prune(&mut self, node_id: TreeNodeID<'a, Self>) -> Result<Self>;

    /// SPR function
    fn spr(
        &mut self,
        edge1: (TreeNodeID<'a, Self>, TreeNodeID<'a, Self>),
        edge2: (TreeNodeID<'a, Self>, TreeNodeID<'a, Self>),
    )->Result<()> {
        let pruned_tree = SPR::prune(self, edge1.1)?;
        SPR::graft(self, pruned_tree, edge2)
    }
}

/// A trait describing Nearest Neighbour interchange operations
pub trait NNI<'a>
where
    Self: RootedTree<'a> + Sized,
{
    /// Performs an NNI operation
    fn nni(&mut self, parent_id: TreeNodeID<'a, Self>)->Result<()>;
}

/// A trait describing rerooting a tree
pub trait Reroot<'a>
where
    Self: RootedTree<'a> + Sized,
{
    /// Reroots tree at node. **Note: this changes the degree of a node**
    fn reroot_at_node(&mut self, node_id: TreeNodeID<'a, Self>)->Result<()>;
    /// Reroots tree at a split node.
    fn reroot_at_edge(&mut self, edge: (TreeNodeID<'a, Self>, TreeNodeID<'a, Self>))->Result<()>;
}

/// A trait describing blancing a binary tree
pub trait Balance<'a>: Clusters<'a> + SPR<'a> + Sized
where
    TreeNodeID<'a, Self>: Display + Debug + Hash + Clone + Ord,
{
    /// Balances a binary tree
    fn balance_subtree(&mut self)->Result<()>;
}

/// A trait describing subtree queries of a tree
pub trait Subtree<'a>: Ancestors<'a> + DFS<'a> + Sized
where
    TreeNodeID<'a, Self>: Display + Debug + Hash + Clone + Ord,
{
    /// Returns a subtree consisting of only provided nodes
    fn induce_tree(
        &'a self,
        node_id_list: impl IntoIterator<
            Item = TreeNodeID<'a, Self>,
            IntoIter = impl ExactSizeIterator<Item = TreeNodeID<'a, Self>>,
        >,
    ) -> Result<Self>;

    /// Returns subtree starting at provided node.
    fn subtree(&'a self, node_id: TreeNodeID<'a, Self>) -> Result<Self>;
}

/// A trait describing tree contraction operations
pub trait ContractTree<'a>: EulerWalk<'a> + DFS<'a> {
    /// Contracts tree that from post_ord node_id iterator.
    fn contracted_tree_nodes_from_iter(
        &'a self,
        new_tree_root_id: TreeNodeID<'a, Self>,
        leaf_ids: &[TreeNodeID<'a, Self>],
        node_iter: impl Iterator<Item = TreeNodeID<'a, Self>>,
    ) -> impl Iterator<Item = Self::Node> {
        let mut node_map: HashMap<TreeNodeID<'a, Self>, Self::Node> =
            HashMap::from_iter(vec![(new_tree_root_id, self.get_lca(leaf_ids).clone())]);
        let mut remove_list = vec![];
        node_iter
            .map(|x| self.get_node(x).cloned().unwrap())
            .for_each(|mut node| {
                match node.is_leaf() {
                    true => match leaf_ids.contains(&node.get_id()) {
                        true => {
                            node_map.insert(node.get_id(), node);
                        }
                        false => {}
                    },
                    false => {
                        let node_children_ids = node.get_children().collect_vec();
                        for child_id in &node_children_ids {
                            match node_map.contains_key(child_id) {
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

                                match remove_list.contains(&child_node_id) {
                                    true => {
                                        node.remove_child(&child_node_id);
                                        let grandchildren_ids = node_map
                                            .get(&child_node_id)
                                            .unwrap()
                                            .get_children()
                                            .collect_vec();
                                        for grandchild_id in grandchildren_ids {
                                            node_map
                                                .get_mut(&grandchild_id)
                                                .unwrap()
                                                .set_parent(Some(node.get_id()));
                                            node.add_child(grandchild_id);
                                        }
                                    }
                                    false => {}
                                }
                                remove_list.push(node.get_id());
                                node_map.insert(node.get_id(), node);
                            }
                            _ => {
                                // node has multiple children
                                // for each child, suppress child if child is in remove list
                                node_children_ids.into_iter().for_each(|chid| {
                                    match remove_list.contains(&chid) {
                                        true => {
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
                                                .collect_vec();
                                            for grandchild in node_grandchildren {
                                                node.add_child(grandchild);
                                                node_map
                                                    .get_mut(&grandchild)
                                                    .unwrap()
                                                    .set_parent(Some(node.get_id()))
                                            }
                                        }
                                        false => {}
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
        &'a self,
        leaf_ids: &[TreeNodeID<'a, Self>],
    ) -> impl Iterator<Item = Self::Node> {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let node_postord_iter = self.postord_nodes(new_tree_root_id);
        let mut node_map: HashMap<TreeNodeID<'a, Self>, Self::Node> =
            HashMap::from_iter(vec![(new_tree_root_id, self.get_lca(leaf_ids).clone())]);
        let mut remove_list = vec![];
        node_postord_iter.for_each(|orig_node| {
            let mut node = orig_node.clone();
            match node.is_leaf() {
                true => match leaf_ids.contains(&node.get_id()) {
                    true => {
                        node_map.insert(node.get_id(), node);
                    }
                    false => {}
                },
                false => {
                    let node_children_ids = node.get_children().collect_vec();
                    for child_id in &node_children_ids {
                        match node_map.contains_key(child_id) {
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

                            match remove_list.contains(&child_node_id) {
                                true => {
                                    node.remove_child(&child_node_id);
                                    let grandchildren_ids = node_map
                                        .get(&child_node_id)
                                        .unwrap()
                                        .get_children()
                                        .collect_vec();
                                    for grandchild_id in grandchildren_ids {
                                        node_map
                                            .get_mut(&grandchild_id)
                                            .unwrap()
                                            .set_parent(Some(node.get_id()));
                                        node.add_child(grandchild_id);
                                    }
                                }
                                false => {}
                            }
                            remove_list.push(node.get_id());
                            node_map.insert(node.get_id(), node);
                        }
                        _ => {
                            // node has multiple children
                            // for each child, suppress child if child is in remove list
                            node_children_ids.into_iter().for_each(|chid| {
                                match remove_list.contains(&chid) {
                                    true => {
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
                                            .collect_vec();
                                        for grandchild in node_grandchildren {
                                            node.add_child(grandchild);
                                            node_map
                                                .get_mut(&grandchild)
                                                .unwrap()
                                                .set_parent(Some(node.get_id()))
                                        }
                                    }
                                    false => {}
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
    fn contract_tree(&self, leaf_ids: &[TreeNodeID<'a, Self>]) -> Result<Self>;

    /// Returns a contracted tree from an iterator containing NodeID's
    fn contract_tree_from_iter(
        &self,
        leaf_ids: &[TreeNodeID<'a, Self>],
        node_iter: impl Iterator<Item = TreeNodeID<'a, Self>>,
    ) -> Result<Self>;
}
