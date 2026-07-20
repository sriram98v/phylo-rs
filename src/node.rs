/// Module with traits of rooted tree nodes
pub mod simple_rnode;

use crate::node::simple_rnode::{
    EdgeWeight, NodeTaxa, NodeWeight, RootedAnnotatedNode, RootedMetaNode, RootedTreeNode,
    RootedWeightedNode, RootedZetaNode,
};
use std::fmt::{Debug, Display};
use std::sync::Arc;

/// Default NodeID type
pub type NodeID = usize;

/// Marks a node with no parent, i.e. the root.
///
/// `Option<NodeID>` costs 16 bytes rather than 8: `usize` has no invalid bit
/// pattern for the discriminant to occupy, so the `Option` cannot pack. An
/// arena can never hold `usize::MAX` nodes, so the top value is free to stand
/// for absence.
///
/// This is a storage detail. [`RootedTreeNode::get_parent`] and
/// [`RootedTreeNode::set_parent`] still speak in `Option`, and the serialized
/// form is still `null`.
const NO_PARENT: NodeID = NodeID::MAX;

/// Serializes [`Node::parent`] as the `Option` it used to be.
///
/// The sentinel is how the field is stored, not how it is written down:
/// without this, a root's parent would serialize as `18446744073709551615`
/// instead of `null`, and every tree already on disk would fail to load.
#[cfg(feature = "serde")]
mod parent_as_option {
    use super::{NodeID, NO_PARENT};
    use serde::{Deserialize, Deserializer, Serialize, Serializer};

    pub(super) fn serialize<S: Serializer>(parent: &NodeID, ser: S) -> Result<S::Ok, S::Error> {
        (*parent != NO_PARENT).then_some(*parent).serialize(ser)
    }

    pub(super) fn deserialize<'de, D: Deserializer<'de>>(de: D) -> Result<NodeID, D::Error> {
        Ok(Option::<NodeID>::deserialize(de)?.unwrap_or(NO_PARENT))
    }
}

/// Default NodeID type
pub type PhyloNode = Node<String, f32, f32>;

/// Default NodeID type
pub type DemoNode = Node<u32, f32, f32>;

/// A node structure in an arena-memory managed tree, linking to connected neighbours via NodeID
#[derive(Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Node<T, W, Z>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    /// A unique identifier for a node
    id: NodeID,
    /// A link to the node parent, or [`NO_PARENT`] for the root
    #[cfg_attr(feature = "serde", serde(with = "parent_as_option"))]
    parent: NodeID,
    /// Children of node
    children: Vec<NodeID>,
    /// Taxa annotation of node
    taxa: Option<Arc<T>>,
    /// Weight of edge ending in node
    weight: Option<W>,
    /// Real number annotation of node (used by some algorithms)
    zeta: Option<Z>,
    /// Raw, uninterpreted per-node annotation (e.g. the contents of `[...]`
    /// comments captured while parsing Newick). Stored verbatim so callers can
    /// parse it however they like; `Arc<str>` so cloning a node — or the whole
    /// tree — shares the string rather than reallocating it.
    annotation: Option<Arc<str>>,
}

impl<T, W, Z> RootedTreeNode for Node<T, W, Z>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type NodeID = NodeID;

    fn new(id: Self::NodeID) -> Self {
        Node {
            id,
            parent: NO_PARENT,
            children: vec![],
            taxa: None,
            weight: None,
            zeta: None,
            annotation: None,
        }
    }

    fn get_id(&self) -> Self::NodeID {
        self.id
    }

    fn set_id(&mut self, id: Self::NodeID) {
        self.id = id
    }

    fn set_parent(&mut self, parent: Option<Self::NodeID>) {
        self.parent = parent.unwrap_or(NO_PARENT);
    }

    fn get_parent(&self) -> Option<Self::NodeID> {
        (self.parent != NO_PARENT).then_some(self.parent)
    }

    fn get_children(&self) -> &[Self::NodeID] {
        &self.children
    }

    fn add_child(&mut self, child: Self::NodeID) {
        // A Vec's first push jumps straight to capacity 4. Nodes in a
        // phylogeny are overwhelmingly bifurcating, so that holds twice the
        // memory the children ever use: measured at 15999 of 31999 nodes
        // sitting at len 2, capacity 4, for 250 KiB of a 16000-taxon tree.
        //
        // Ask for exactly what a bifurcation needs. A third child still grows
        // by the usual doubling, so multifurcating trees are unaffected.
        if self.children.capacity() == 0 {
            self.children.reserve_exact(2);
        }
        self.children.push(child);
    }

    fn remove_child(&mut self, child: &Self::NodeID) {
        self.children.retain(|x| x != child);
    }
}

impl<T, W, Z> Node<T, W, Z>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    /// Returns the number of bytes this node has allocated on the heap.
    ///
    /// Counts the children vector's allocation, using capacity rather than
    /// length since that is what is actually reserved.
    ///
    /// Taxa payloads are excluded: they are shared with the tree's taxa map via
    /// [`Arc`], so attributing them to a node would double-count them.
    pub fn heap_size(&self) -> usize {
        self.children.capacity() * std::mem::size_of::<NodeID>()
    }

    /// Returns a reference to the inner Arc for shared ownership with the taxa map.
    ///
    /// Only `SimpleRootedTree` shares taxa this way, so this is dead code
    /// without the feature that defines it.
    #[cfg(feature = "simple_rooted_tree")]
    pub(crate) fn get_taxa_arc(&self) -> Option<&Arc<T>> {
        self.taxa.as_ref()
    }

    /// Sets the taxa field from a pre-built Arc, sharing ownership with the taxa map.
    #[cfg(feature = "simple_rooted_tree")]
    pub(crate) fn set_taxa_arc(&mut self, taxa: Option<Arc<T>>) {
        self.taxa = taxa;
    }
}

impl<T, W, Z> RootedMetaNode for Node<T, W, Z>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type Meta = T;

    fn get_taxa(&self) -> Option<&Self::Meta> {
        self.taxa.as_deref()
    }

    fn set_taxa(&mut self, taxa: Option<Self::Meta>) {
        self.taxa = taxa.map(Arc::new);
    }
}

impl<T, W, Z> RootedWeightedNode for Node<T, W, Z>
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

impl<T, W, Z> RootedZetaNode for Node<T, W, Z>
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

impl<T, W, Z> RootedAnnotatedNode for Node<T, W, Z>
where
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    fn get_annotation(&self) -> Option<&str> {
        self.annotation.as_deref()
    }

    fn set_annotation(&mut self, annotation: Option<Arc<str>>) {
        self.annotation = annotation;
    }
}

impl<T, W, Z> Debug for Node<T, W, Z>
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

impl<T, W, Z> Display for Node<T, W, Z>
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
