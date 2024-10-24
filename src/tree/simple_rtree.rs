use crate::node::simple_rnode::*;
use itertools::Itertools;
use std::fmt::Debug;

/// A type alias for Tree Node ID
pub type TreeNodeID<T> = <<T as RootedTree>::Node as RootedTreeNode>::NodeID;
/// A type alias for Tree Node meta annotation
pub type TreeNodeMeta<T> = <<T as RootedTree>::Node as RootedMetaNode>::Meta;
/// A type alias for Tree edge weight
pub type TreeNodeWeight<T> = <<T as RootedTree>::Node as RootedWeightedNode>::Weight;

/// A trait describing the behaviour of a rooted tree#[allow(clippy::needless_lifetimes)]
#[allow(clippy::needless_lifetimes)]
pub trait RootedTree: Clone + Sync
where
    Self::Node: RootedTreeNode + Debug,
{
    /// An associated node type for a rooted tree
    type Node;

    /// Returns reference to node by ID
    fn get_node<'a>(&'a self, node_id: TreeNodeID<Self>) -> Option<&'a Self::Node>;

    /// Returns a mutable reference to a node
    fn get_node_mut<'a>(&'a mut self, node_id: TreeNodeID<Self>) -> Option<&'a mut Self::Node>;

    /// Reurns an iterator over all NodeID's
    fn get_node_ids(&self) -> impl Iterator<Item = TreeNodeID<Self>>;

    /// Returns an iterator with immutable references to nodes
    fn get_nodes<'a>(&'a self) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        self.get_node_ids()
            .map(|id| self.get_node(id).unwrap())
            .collect_vec()
            .into_iter()
    }

    /// Returns iterator with mutable references to nodes
    fn get_nodes_mut<'a>(&'a mut self) -> impl Iterator<Item = &'a mut Self::Node>;

    /// Returns NodeID of root node
    fn get_root_id(&self) -> TreeNodeID<Self>;

    /// Sets node with NodeID and root node
    fn set_root(&mut self, node_id: TreeNodeID<Self>);

    /// Inserts a floating node into tree.
    fn set_node(&mut self, node: Self::Node);

    /// Adds node as child to an existing node in tree.
    fn add_child(&mut self, parent_id: TreeNodeID<Self>, child: Self::Node) {
        let new_child_id = child.get_id();
        self.set_node(child);
        self.get_node_mut(parent_id)
            .unwrap()
            .add_child(new_child_id);
        self.get_node_mut(new_child_id)
            .unwrap()
            .set_parent(Some(parent_id));
    }

    /// Removes node from tree while deleting any edges if they exist
    fn remove_node(&mut self, node_id: TreeNodeID<Self>) -> Option<Self::Node>;

    /// Removes nodes from tree without deleting any edges that may exist
    fn delete_node(&mut self, node_id: TreeNodeID<Self>);

    /// Returns true if node with node_id exists in tree
    fn contains_node(&self, node_id: TreeNodeID<Self>) -> bool {
        self.get_node(node_id).is_some()
    }

    /// Removes internal nodes of degree 2 and any floating nodes
    fn clean(&mut self) {
        let node_iter = self.get_nodes().cloned().collect_vec();
        for node in &node_iter {
            // remove root with only one child
            let node_id = node.get_id();
            if node.get_id() == self.get_root_id() && node.degree() < 2 {
                let new_root = self.get_root().get_children().next().unwrap();
                self.set_root(new_root);
                self.get_node_mut(self.get_root_id())
                    .unwrap()
                    .set_parent(None);
                self.remove_node(node_id);
            }
            // remove nodes with only one child
            else if !node.is_leaf() && node.get_parent().is_some() && node.degree() < 3 {
                let parent_id = self.get_node_parent_id(node_id);
                let child_id = node.get_children().next().unwrap();
                self.get_node_mut(child_id).unwrap().set_parent(parent_id);
                self.get_node_mut(parent_id.unwrap())
                    .unwrap()
                    .add_child(child_id);
                self.remove_node(node.get_id());
            }
            // Removing dangling references to pruned children
            for chid in node.get_children() {
                if !self.get_nodes().map(|x| x.get_id()).contains(&chid) {
                    self.get_node_mut(node_id).unwrap().remove_child(&chid);
                }
            }
        }
    }

    /// Removes all nodes from tree except root node
    fn clear(&mut self);

    /// Deletes an edge from the tree without deleting an nodes
    fn delete_edge(&mut self, parent_id: TreeNodeID<Self>, child_id: TreeNodeID<Self>) {
        self.get_node_mut(parent_id)
            .unwrap()
            .remove_child(&child_id);
        self.get_node_mut(child_id).unwrap().set_parent(None);
    }

    /// Inserts nodes into tree from iterator. Note: this will overwrite any existing node with a NodeID that already exists in tree.
    fn set_nodes(
        &mut self,
        node_list: impl IntoIterator<
            Item = Self::Node,
            IntoIter = impl ExactSizeIterator<Item = Self::Node>,
        >,
    ) {
        for node in node_list {
            self.set_node(node);
        }
    }

    /// Splits an edge in the tree with provided node.
    fn split_edge(&mut self, edge: (TreeNodeID<Self>, TreeNodeID<Self>), node: Self::Node) {
        let p_id = edge.0;
        let c_id = edge.1;
        let n_id = node.get_id();
        self.set_node(node);
        self.get_node_mut(p_id).unwrap().remove_child(&c_id);
        self.set_child(p_id, n_id);
        self.set_child(n_id, c_id);
    }

    /// Add node as a sibling to the provided NodeID.
    fn add_sibling(
        &mut self,
        node_id: TreeNodeID<Self>,
        split_node: Self::Node,
        sibling_node: Self::Node,
    ) {
        let node_parent_id = self.get_node_parent_id(node_id).unwrap();
        let split_node_id = split_node.get_id();
        self.split_edge((node_parent_id, node_id), split_node);
        self.add_child(split_node_id, sibling_node);
    }

    /// Returns iterator of immutable references to leaf nodes in tree.
    fn get_leaves<'a>(&'a self) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        self.get_nodes()
            .filter(|x| x.is_leaf())
            .collect_vec()
            .into_iter()
    }

    /// Returns an iterator of leaf NodeID's
    fn get_leaf_ids(&self) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        self.get_node_ids()
            .filter(|x| self.is_leaf(*x))
            .collect_vec()
            .into_iter()
    }

    /// Returns an immutable reference to root node
    fn get_root<'a>(&'a self) -> &'a Self::Node {
        self.get_node(self.get_root_id()).unwrap()
    }

    /// Returns a mutable reference to the root node
    fn get_root_mut<'a>(&'a mut self) -> &'a mut Self::Node {
        self.get_node_mut(self.get_root_id()).unwrap()
    }

    /// creates an edge from node with parent ID to child ID. The child node must already exist in tree.
    fn set_child(&mut self, parent_id: TreeNodeID<Self>, child_id: TreeNodeID<Self>) {
        self.get_node_mut(parent_id).unwrap().add_child(child_id);
        self.get_node_mut(child_id)
            .unwrap()
            .set_parent(Some(parent_id));
    }

    /// Removes edge from prant to child without deleting either node.
    fn remove_child<'a>(&'a mut self, parent_id: TreeNodeID<Self>, child_id: TreeNodeID<Self>) {
        self.get_node_mut(parent_id)
            .unwrap()
            .remove_child(&child_id);
    }

    /// Removes set of children from parent node.
    fn remove_children(
        &mut self,
        parent_id: TreeNodeID<Self>,
        child_ids: impl Iterator<Item = TreeNodeID<Self>>,
    ) {
        for child_id in child_ids {
            self.get_node_mut(parent_id)
                .unwrap()
                .remove_child(&child_id);
        }
    }

    /// Removes all children from parent node.
    fn remove_all_children(&mut self, node_id: TreeNodeID<Self>) {
        let node_children_ids = self
            .get_node_children_ids(node_id)
            .collect_vec()
            .into_iter();
        self.remove_children(node_id, node_children_ids);
    }

    /// Returns parent ID of a node in tree
    fn get_node_parent_id(&self, node_id: TreeNodeID<Self>) -> Option<TreeNodeID<Self>> {
        self.get_node(node_id).unwrap().get_parent()
    }

    /// Returns immutable reference to parent for a node
    fn get_node_parent<'a>(&'a self, node_id: TreeNodeID<Self>) -> Option<&'a Self::Node> {
        self.get_node(self.get_node_parent_id(node_id)?)
    }

    /// Returns immutable reference to parent for a node
    fn get_node_parent_mut<'a>(
        &'a mut self,
        node_id: TreeNodeID<Self>,
    ) -> Option<&'a mut Self::Node> {
        self.get_node_mut(self.get_node_parent_id(node_id)?)
    }

    /// Returns an iterator of immutable references to children of a node
    fn get_node_children<'a>(
        &'a self,
        node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = &'a Self::Node> {
        let node = self.get_node(node_id).unwrap();
        node.get_children().map(|x| self.get_node(x).unwrap())
    }

    /// Returns an iterator of node children ids
    fn get_node_children_ids(
        &self,
        node_id: TreeNodeID<Self>,
    ) -> impl ExactSizeIterator<Item = TreeNodeID<Self>> {
        self.get_node(node_id)
            .unwrap()
            .get_children()
            .collect_vec()
            .into_iter()
    }

    /// Returns degree of a node
    fn node_degree<'a>(&'a self, node_id: TreeNodeID<Self>) -> usize {
        self.get_node(node_id).unwrap().degree()
    }

    /// Returns depth of node as number of edges in the path between node and root.
    fn get_node_depth(&self, node_id: TreeNodeID<Self>) -> usize {
        let mut start_id = node_id;
        let mut depth = 0;
        while let Some(parent_id) = self.get_node_parent_id(start_id) {
            depth += 1;
            start_id = parent_id;
        }
        depth
    }

    /// Returns true if tree is binary
    fn is_binary(&self) -> bool {
        for node_id in self.get_node_ids() {
            if node_id == self.get_root_id() {
                if self.node_degree(node_id) != 2 {
                    return false;
                }
            } else if self.node_degree(node_id) % 2 != 1 {
                return false;
            }
        }
        true
    }

    /// Returns true if node with node_id is a leaf node
    fn is_leaf(&self, node_id: TreeNodeID<Self>) -> bool {
        self.get_node(node_id).unwrap().is_leaf()
    }

    /// Returns total number of nodes in tree
    fn num_nodes(&self) -> usize {
        self.get_node_ids().collect_vec().len()
    }

    /// Returns iterator of immutable references to siblings of a node.
    fn get_siblings<'a>(
        &'a self,
        node_id: TreeNodeID<Self>,
    ) -> impl Iterator<Item = &'a Self::Node> {
        let parent_id = self
            .get_node_parent_id(node_id)
            .expect("Node has no siblings!");
        return self
            .get_node_children(parent_id)
            .filter(move |x| x.get_id() != node_id);
    }

    /// Returns iterator of NodeID's of node siblings
    fn get_sibling_ids(&self, node_id: TreeNodeID<Self>) -> impl Iterator<Item = TreeNodeID<Self>> {
        let parent_id = self
            .get_node_parent_id(node_id)
            .expect("Root does not have siblings!");
        let sibling_ids = self
            .get_node_children_ids(parent_id)
            .filter(move |x| x != &node_id);
        sibling_ids
    }

    /// Connects a nodes children to it's parent, then deletes all edges to the node, without deleting the node from the tree
    fn supress_node<'a>(&'a mut self, node_id: TreeNodeID<Self>) -> Option<()> {
        let node_parent_id = self.get_node_parent_id(node_id)?;
        let node_children_ids = self.get_node_children_ids(node_id).collect_vec();
        for child_id in node_children_ids.as_slice() {
            let child = self.get_node_mut(*child_id)?;
            child.set_parent(Some(node_parent_id));
        }
        let node_parent = self.get_node_parent_mut(node_id)?;
        for child_id in node_children_ids {
            node_parent.add_child(child_id);
        }
        self.remove_node(node_id);
        Some(())
    }

    /// Supresses all nodes of degree 2
    fn supress_unifurcations<'a>(&'a mut self);
}

/// A trait describing the behaviour of a rooted tree where some of the nodes have a meta annotation. The terms meta and taxa are used interchangably here.
#[allow(clippy::needless_lifetimes)]
pub trait RootedMetaTree: RootedTree
where
    Self::Node: RootedMetaNode,
{
    ///  Returns an immutable reference to a node with a give meta annotation
    fn get_taxa_node<'a>(&'a self, taxa: &TreeNodeMeta<Self>) -> Option<&'a Self::Node>;

    /// Returns the node id of a node with a meta annotation
    fn get_taxa_node_id(&self, taxa: &TreeNodeMeta<Self>) -> Option<TreeNodeID<Self>> {
        Some(self.get_taxa_node(taxa)?.get_id())
    }

    /// Returns totla number of nodes with a meta annotation
    fn num_taxa(&self) -> usize;

    /// Sets the emta annotation of a node
    fn set_node_taxa<'a>(
        &'a mut self,
        node_id: TreeNodeID<Self>,
        taxa: Option<TreeNodeMeta<Self>>,
    ) {
        self.get_node_mut(node_id).unwrap().set_taxa(taxa)
    }

    /// Returns an immutable reference to the meta annotation of a node, and None is there is no meta annotation
    fn get_node_taxa<'a>(&'a self, node_id: TreeNodeID<Self>) -> Option<&'a TreeNodeMeta<Self>> {
        self.get_node(node_id).unwrap().get_taxa()
    }

    /// Returns a deep copy of the meta annotation of a node
    fn get_node_taxa_cloned(&self, node_id: TreeNodeID<Self>) -> Option<TreeNodeMeta<Self>>;

    /// Returns an iterator with immutable references to all meta annotations in a tree.
    fn get_taxa_space<'a>(&'a self) -> impl ExactSizeIterator<Item = &'a TreeNodeMeta<Self>> {
        self.get_nodes()
            .map(|node| node.get_taxa())
            .filter(|x| x.is_none())
            .map(|x| x.unwrap())
            .collect_vec()
            .into_iter()
    }
}

/// A trait describing the behaviour of a rooted tree where some of the edges are weighted
#[allow(clippy::needless_lifetimes)]
pub trait RootedWeightedTree: RootedTree
where
    Self::Node: RootedWeightedNode,
{
    /// Sets all edge weights to None
    fn unweight<'a>(&'a mut self) {
        let ids = self.get_node_ids().collect_vec();
        for id in ids {
            self.get_node_mut(id).unwrap().set_weight(None);
        }
    }

    /// Sets edge weight to None
    fn set_edge_weight<'a>(
        &'a mut self,
        edge: (TreeNodeID<Self>, TreeNodeID<Self>),
        edge_weight: Option<TreeNodeWeight<Self>>,
    ) {
        self.get_node_mut(edge.1).unwrap().set_weight(edge_weight);
    }

    /// Returns true if edge weight not None
    fn is_weighted<'a>(&'a self) -> bool {
        for node_id in self.get_node_ids() {
            if self.get_node(node_id).unwrap().get_weight().is_none() {
                return false;
            }
        }
        true
    }

    /// Returns weight of edge
    fn get_edge_weight<'a>(
        &'a self,
        _parent_id: TreeNodeID<Self>,
        child_id: TreeNodeID<Self>,
    ) -> Option<TreeNodeWeight<Self>> {
        self.get_node(child_id).unwrap().get_weight()
    }
}
