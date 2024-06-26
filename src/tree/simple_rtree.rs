use crate::node::simple_rnode::*;
use itertools::Itertools;
use std::{fmt::Debug, hash::Hash};

/// A trait describing the behaviour of a rooted tree
pub trait RootedTree: Clone + Sync
where
    Self::Node: RootedTreeNode<NodeID = Self::NodeID> + Debug,
{
    // identifier for nodes
    type NodeID: Copy + PartialEq + Eq + Hash + Debug + Send;
    type Node;

    /// Returns reference to node by ID
    fn get_node(&self, node_id: Self::NodeID) -> Option<&Self::Node>;

    fn get_node_mut(&mut self, node_id: Self::NodeID) -> Option<&mut Self::Node>;

    fn get_node_ids(&self) -> impl Iterator<Item = Self::NodeID>;

    fn get_nodes(&self) -> impl ExactSizeIterator<Item = &Self::Node>;

    fn get_root_id(&self) -> Self::NodeID;

    fn set_root(&mut self, node_id: Self::NodeID);

    /// Returns reference to node by ID
    fn set_node(&mut self, node: Self::Node);

    fn add_child(&mut self, parent_id: Self::NodeID, child: Self::Node);

    fn remove_node(&mut self, node_id: Self::NodeID) -> Option<Self::Node>;

    fn delete_node(&mut self, node_id: Self::NodeID);

    fn contains_node(&self, node_id: Self::NodeID) -> bool;

    fn clean(&mut self);

    fn clear(&mut self) {
        let node_ids = self
            .get_node_ids()
            .filter(|x| x != &self.get_root_id())
            .collect_vec();
        for node_id in node_ids {
            self.delete_node(node_id);
        }
        self.remove_all_children(self.get_root_id());
    }

    fn delete_edge(&mut self, parent_id: Self::NodeID, child_id: Self::NodeID) {
        self.get_node_mut(parent_id)
            .unwrap()
            .remove_child(&child_id);
        self.get_node_mut(child_id).unwrap().set_parent(None);
    }

    fn set_nodes(
        &mut self,
        node_list: impl IntoIterator<
            Item = Self::Node,
            IntoIter = impl ExactSizeIterator<Item = Self::Node>,
        >,
    ) {
        for node in node_list.into_iter() {
            self.set_node(node);
        }
    }

    fn split_edge(&mut self, edge: (Self::NodeID, Self::NodeID), node: Self::Node) {
        let p_id = edge.0;
        let c_id = edge.1;
        let n_id = node.get_id();
        self.set_node(node);
        self.get_node_mut(p_id).unwrap().remove_child(&c_id);
        self.set_child(p_id, n_id);
        self.set_child(n_id, c_id);
    }

    fn add_sibling(
        &mut self,
        node_id: Self::NodeID,
        split_node: Self::Node,
        sibling_node: Self::Node,
    ) {
        let node_parent_id = self.get_node_parent_id(node_id).unwrap();
        let split_node_id = split_node.get_id();
        self.split_edge((node_parent_id, node_id), split_node);
        self.add_child(split_node_id, sibling_node);
    }

    fn get_leaves(&self) -> impl ExactSizeIterator<Item = &Self::Node> {
        self.get_nodes()
            .filter(|x| x.is_leaf())
            .collect_vec()
            .into_iter()
    }

    /// Get root node
    fn get_root(&self) -> &Self::Node {
        self.get_node(self.get_root_id()).unwrap()
    }

    fn get_root_mut(&mut self) -> &mut Self::Node {
        self.get_node_mut(self.get_root_id()).unwrap()
    }

    fn set_child(&mut self, parent_id: Self::NodeID, child_id: Self::NodeID) {
        self.get_node_mut(parent_id).unwrap().add_child(child_id);
        self.get_node_mut(child_id)
            .unwrap()
            .set_parent(Some(parent_id));
    }

    fn remove_child(&mut self, parent_id: Self::NodeID, child_id: Self::NodeID) {
        self.get_node_mut(parent_id)
            .unwrap()
            .remove_child(&child_id);
    }

    fn remove_children(&mut self, parent_id: Self::NodeID, child_ids: Vec<Self::NodeID>) {
        for child_id in child_ids {
            self.get_node_mut(parent_id)
                .unwrap()
                .remove_child(&child_id);
        }
    }

    fn remove_all_children(&mut self, node_id: Self::NodeID) {
        let node_children_ids = self.get_node_children_ids(node_id).collect_vec();
        self.remove_children(node_id, node_children_ids);
    }

    fn get_node_parent_id(&self, node_id: Self::NodeID) -> Option<Self::NodeID> {
        self.get_node(node_id).unwrap().get_parent()
    }

    fn get_node_parent(&self, node_id: Self::NodeID) -> Option<Self::Node> {
        match self.get_node_parent_id(node_id) {
            Some(id) => self.get_node(id).cloned(),
            None => None,
        }
    }

    fn get_node_children(
        &self,
        node_id: Self::NodeID,
    ) -> impl ExactSizeIterator<Item = Self::Node> {
        let node = self.get_node(node_id).unwrap();
        node.get_children()
            .map(|x| self.get_node(x).cloned().unwrap())
    }

    fn get_node_children_ids(
        &self,
        node_id: Self::NodeID,
    ) -> impl ExactSizeIterator<Item = Self::NodeID> {
        self.get_node(node_id)
            .unwrap()
            .get_children()
            .collect_vec()
            .into_iter()
    }

    fn node_degree(&self, node_id: Self::NodeID) -> usize {
        self.get_node(node_id).unwrap().degree()
    }
    fn get_node_depth(&self, node_id: Self::NodeID) -> usize {
        let mut start_id = node_id;
        let mut depth = 0;
        while let Some(parent_id) = self.get_node_parent_id(start_id) {
            depth += 1;
            start_id = parent_id;
        }
        depth
    }
    fn is_binary(&self) -> bool {
        for node in self.get_nodes() {
            if node.get_children().collect_vec().len() > 2 {
                return false;
            }
        }
        true
    }

    fn is_leaf(&self, node_id: &Self::NodeID) -> bool {
        self.get_node(*node_id).unwrap().is_leaf()
    }

    fn num_nodes(&self) -> usize {
        self.get_nodes().len()
    }

    fn get_siblings(&self, node_id: Self::NodeID) -> impl Iterator<Item = Self::Node> {
        let parent_id = self
            .get_node_parent_id(node_id)
            .expect("Node has no siblings!");
        return self
            .get_node_children(parent_id)
            .filter(move |x| x.get_id() != node_id);
    }

    fn get_sibling_ids(&self, node_id: Self::NodeID) -> impl Iterator<Item = Self::NodeID> {
        return self.get_siblings(node_id).map(|x| x.get_id());
    }

    fn supress_node(&mut self, _node_id: Self::NodeID) {
        todo!()
    }

    fn supress_unifurcations(&mut self) {
        todo!()
    }
}

pub trait RootedMetaTree: RootedTree
where
    Self::Node: RootedMetaNode<Meta = Self::Meta>,
{
    type Meta: PartialEq + Send + Sync;

    fn get_taxa_node(&self, taxa: &Self::Meta) -> Option<&Self::Node>;

    fn get_taxa_node_id(&self, taxa: &Self::Meta) -> Option<Self::NodeID> {
        self.get_taxa_node(taxa).map(|node| node.get_id())
    }
    fn num_taxa(&self) -> usize {
        self.get_nodes().filter(|f| f.is_leaf()).collect_vec().len()
    }
    fn set_node_taxa(&mut self, node_id: Self::NodeID, taxa: Option<Self::Meta>) {
        self.get_node_mut(node_id).unwrap().set_taxa(taxa)
    }
    fn get_node_taxa(&self, node_id: Self::NodeID) -> Option<Self::Meta> {
        self.get_node(node_id).unwrap().get_taxa()
    }
    fn get_taxa_space(&self) -> impl Iterator<Item = Self::Meta> {
        self.get_nodes()
            .filter(|x| x.get_taxa().is_some())
            .map(|x| x.get_taxa().unwrap())
    }
}

pub trait RootedWeightedTree: RootedTree
where
    Self::Node: RootedWeightedNode<Weight = Self::Weight>,
{
    type Weight: PartialEq;

    fn unweight(&mut self) {
        let node_ids = self.get_node_ids().collect_vec();
        for node_id in node_ids {
            self.get_node_mut(node_id).unwrap().set_weight(None);
        }
    }

    fn set_edge_weight(
        &mut self,
        edge: (Self::NodeID, Self::NodeID),
        edge_weight: Option<Self::Weight>,
    ) {
        self.get_node_mut(edge.1).unwrap().set_weight(edge_weight);
    }

    fn is_weighted(&self) -> bool {
        for node_id in self.get_node_ids() {
            if self.get_node(node_id).unwrap().get_weight().is_none() {
                return false;
            }
        }
        true
    }

    fn get_edge_weight(
        &self,
        _parent_id: Self::NodeID,
        child_id: Self::NodeID,
    ) -> Option<Self::Weight> {
        self.get_node(child_id).unwrap().get_weight()
    }
}
