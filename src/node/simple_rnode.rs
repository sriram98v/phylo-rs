use num::{Float, Signed, Zero};
use std::{hash::Hash, marker::Sync, fmt::{Debug, Display}, str::FromStr, iter::Sum};

/// Trait bound alias for Edge Weight.
pub trait EdgeWeight: Display + Debug + Sum + FromStr + Float + Signed + Sync + Send{}

/// Trait bound alias for Node Weight.
pub trait NodeWeight: Display + Debug + Sum + FromStr + Float + Zero + Signed + Sync + Send{}

/// Trait bound alias for Node Taxa.
pub trait NodeTaxa: Display + Debug + Clone + FromStr + Ord + Hash + Sync + Send{}
// /// Trait bound alias for Node ID.
// pub trait NodeID: Display + Debug + Hash + Ord + Eq + Copy + Sync{}

impl NodeTaxa for String{}
impl NodeTaxa for u32{}
impl EdgeWeight for f32{}
impl EdgeWeight for f64{}
impl NodeWeight for f32{}
impl NodeWeight for f64{}

/// A trait describing the behaviour of a Node in a n-ary tree
pub trait RootedTreeNode
where
    Self: Clone,
{
    /// Associate type for node identifier. Should be unique within a tree
    type NodeID: Display + Debug + Hash + Ord + PartialEq + Eq + Copy;

    /// Creates a new node with provided id
    fn new(id: Self::NodeID) -> Self;

    /// Returns id of node
    fn get_id(&self) -> Self::NodeID;

    /// Changes id of node
    fn set_id(&mut self, id: Self::NodeID);

    /// Returns id of node parent
    fn get_parent(&self) -> Option<Self::NodeID>;

    /// Sets parent of node
    fn set_parent(&mut self, parent: Option<Self::NodeID>);

    /// Returns slice containing children node ids
    fn get_children(&self) -> &[Self::NodeID];

    /// Add NodeID to node children
    fn add_child(&mut self, child: Self::NodeID);

    /// Remove NodeID from node children
    fn remove_child(&mut self, child: &Self::NodeID);

    /// Checks if node is a leaf node
    fn is_leaf(&self) -> bool {
        self.get_children().is_empty()
    }

    /// Returns Node type as String
    fn node_type(&self) -> String {
        match self.is_leaf() {
            false => "Internal".to_string(),
            true => "Leaf".to_string(),
        }
    }

    /// Adds NodeIDs from Iterator as children
    fn add_children(&mut self, children: impl Iterator<Item = Self::NodeID>) {
        for child in children {
            self.add_child(child);
        }
    }

    /// Removes NodeIDs from Iterator from node children
    fn remove_children(&mut self, children: impl Iterator<Item = Self::NodeID>) {
        for child in children {
            self.remove_child(&child);
        }
    }

    /// Removes all children from node
    fn remove_all_children(&mut self) {
        let children = self.get_children().to_vec();
        for child in children {
            self.remove_child(&child);
        }
    }

    /// Returns number of children of the node.
    fn num_children(&self) -> usize {
        self.get_children().len()
    }

    /// Returns true if node as children.
    fn has_children(&self) -> bool {
        self.num_children() > 0
    }

    /// Returns degree of node.
    fn degree(&self) -> usize {
        match self.get_parent() {
            Some(_) => self.num_children() + 1,
            None => self.num_children(),
        }
    }

    /// Returns ids of all nodes connected to self, including parent if exists.
    fn neighbours(&self) -> impl ExactSizeIterator<Item = Self::NodeID> {
        let mut children = self.get_children().to_vec();
        if let Some(p) = self.get_parent() {
            children.push(p);
        }
        children.into_iter()
    }
}

/// A trait describing the behaviour of a Node in a n-ary tree that carries node annotations
pub trait RootedMetaNode: RootedTreeNode {
    /// Meta annotation of node
    type Meta: NodeTaxa;

    /// Returns node annotation
    fn get_taxa<'a>(&'a self) -> Option<&'a Self::Meta>;

    /// Sets node annotation
    fn set_taxa(&mut self, taxa: Option<Self::Meta>);
}

/// A trait describing the behaviour of a Node in a n-ary tree that has numeric edge annotations
pub trait RootedWeightedNode: RootedTreeNode {
    /// Weight of edge leading into node
    type Weight: EdgeWeight;

    /// Returns weight of edge leading into node
    fn get_weight(&self) -> Option<Self::Weight>;

    /// Sets weight of edge leading into node
    fn set_weight(&mut self, w: Option<Self::Weight>);

    /// Sets weight of edge leading into node as None
    fn unweight(&mut self) {
        self.set_weight(None);
    }

    /// Checks if edge leading into node is weighted
    fn is_weighted(&self) -> bool {
        self.get_weight().is_some()
    }
}

/// A trait describing the behaviour of a Node in a n-ary tree with numeric node annotations
pub trait RootedZetaNode: RootedTreeNode {
    /// Zeta annotation of a node
    type Zeta: NodeWeight;

    /// Returns node annotation
    fn get_zeta(&self) -> Option<Self::Zeta>;

    /// Sets node annotation.
    fn set_zeta(&mut self, w: Option<Self::Zeta>);

    /// Returns true if zeta of node is set
    fn is_zeta_set(&self) -> bool {
        self.get_zeta().is_some()
    }

    /// Sets node zeta to false
    fn remove_zeta(&mut self) {
        self.set_zeta(None);
    }
}
