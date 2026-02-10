#![allow(clippy::needless_lifetimes)]

/// Module with traits of rooted tree nodes
pub mod simple_rnode;

use crate::node::simple_rnode::{
    RootedMetaNode, RootedTreeNode, RootedWeightedNode, RootedZetaNode, EdgeWeight, NodeWeight, NodeTaxa
};
use std::fmt::{Debug, Display};
use std::sync::Arc;

/// Default NodeID type 
pub type NodeID = usize;

/// Default NodeID type 
pub type PhyloNode = Node<String,f32,f32>;

/// Default NodeID type 
pub type DemoNode = Node<u32,f32,f32>;

/// A node structure in an arena-memory managed tree, linking to connected neighbours via NodeID
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    /// A unique identifier for a node
    id: NodeID,
    /// A link to the node parent (set to None for root)
    parent: Option<NodeID>,
    /// Children of node
    children: Vec<NodeID>,
    /// Taxa annotation of node
    taxa: Option<Arc<T>>,
    /// Weight of edge ending in node
    weight: Option<W>,
    /// Real number annotation of node (used by some algorithms)
    zeta: Option<Z>,
}

impl<T,W,Z> RootedTreeNode for Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type NodeID = NodeID;

    fn new(id: Self::NodeID) -> Self {
        Node {
            id,
            parent: None,
            children: vec![],
            taxa: None,
            weight: None,
            zeta: None,
        }
    }

    fn get_id(&self) -> Self::NodeID {
        self.id
    }

    fn set_id(&mut self, id: Self::NodeID) {
        self.id = id
    }

    fn set_parent(&mut self, parent: Option<Self::NodeID>) {
        self.parent = parent;
    }

    fn get_parent(&self) -> Option<Self::NodeID> {
        self.parent
    }

    fn get_children(&self) -> &[Self::NodeID] {
        &self.children
    }

    fn add_child(&mut self, child: Self::NodeID) {
        self.children.push(child);
    }

    fn remove_child(&mut self, child: &Self::NodeID) {
        self.children.retain(|x| x != child);
    }
}

impl<T,W,Z> Node<T,W,Z>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    /// Returns a reference to the inner Arc for shared ownership with the taxa map.
    pub(crate) fn get_taxa_arc(&self) -> Option<&Arc<T>> {
        self.taxa.as_ref()
    }

    /// Sets the taxa field from a pre-built Arc, sharing ownership with the taxa map.
    pub(crate) fn set_taxa_arc(&mut self, taxa: Option<Arc<T>>) {
        self.taxa = taxa;
    }
}

impl<T,W,Z> RootedMetaNode for Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type Meta = T;

    fn get_taxa<'a>(&'a self) -> Option<&'a Self::Meta> {
        self.taxa.as_deref()
    }

    fn set_taxa(&mut self, taxa: Option<Self::Meta>) {
        self.taxa = taxa.map(Arc::new);
    }
}

impl<T,W,Z> RootedWeightedNode for Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type Weight = W;

    fn get_weight(&self) -> Option<Self::Weight> {
        self.weight
    }

    fn set_weight(&mut self, w: Option<Self::Weight>) {
        self.weight = w;
    }
}

impl<T,W,Z> RootedZetaNode for Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type Zeta = Z;

    fn get_zeta(&self) -> Option<Self::Zeta> {
        self.zeta
    }

    fn set_zeta(&mut self, w: Option<Self::Zeta>) {
        self.zeta = w;
    }
}

impl<T,W,Z> Debug for Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}:{}",
            self.get_id(),
            self.node_type(),
            match self.get_taxa() {
                None => "No Taxa".to_string(),
                Some(t) => t.to_string(),
            },
            match self.get_weight() {
                None => "Unweighted".to_string(),
                Some(t) => t.to_string(),
            },
            match self.get_zeta() {
                None => "No Zeta".to_string(),
                Some(z) => z.to_string(),
            }
        )
    }
}

impl<T,W,Z> Display for Node<T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}:{}",
            self.get_id(),
            self.node_type(),
            match self.get_taxa() {
                None => "None".to_string(),
                Some(t) => t.to_string(),
            },
            match self.get_weight() {
                None => "".to_string(),
                Some(t) => t.to_string(),
            },
            match self.get_zeta() {
                None => "No Zeta".to_string(),
                Some(z) => z.to_string(),
            }
        )
    }
}
