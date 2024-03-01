pub mod simple_rnode;

use std::{
    fmt::{Debug, Display},
    hash::{Hash, Hasher},
};

use crate::node::simple_rnode::{RootedPhyloNode, RootedTreeNode, WeightedNode};

// use std::rc::Rc;
use std::sync::Arc;

pub type NodeID = Arc<usize>;

#[derive(Clone)]
pub struct Node {
    id: NodeID,
    parent: Option<NodeID>,
    children: Vec<NodeID>,
    taxa: Option<Arc<String>>,
    weight: Option<f64>,
}

impl RootedTreeNode for Node {
    type NodeID = NodeID;

    fn new(id: NodeID) -> Self {
        Node {
            id: Arc::clone(&id),
            parent: None,
            children: vec![],
            taxa: None,
            weight: None,
        }
    }

    fn get_id(&self) -> Self::NodeID {
        Arc::clone(&self.id)
    }

    fn set_id(&mut self, id: Self::NodeID) {
        self.id = Arc::clone(&id)
    }

    fn set_parent(&mut self, parent: Option<Self::NodeID>) {
        self.parent = parent;
    }
    fn get_parent(&self) -> Option<Self::NodeID> {
        self.parent.clone()
    }
    fn get_children(
        &self,
    ) -> impl IntoIterator<Item = Self::NodeID, IntoIter = impl ExactSizeIterator<Item = Self::NodeID>>
    {
        self.children.clone()
    }
    fn add_child(&mut self, child: Self::NodeID) {
        self.children.push(child);
    }
    fn remove_child(&mut self, child: &Self::NodeID) {
        self.children.retain(|x| x != child);
    }
}

impl RootedPhyloNode for Node {
    type Taxa = String;

    fn get_taxa(&self) -> Option<Self::Taxa> {
        match &self.taxa {
            None => None,
            Some(t) => Some(t.to_string()),
        }
    }

    fn set_taxa(&mut self, taxa: Option<String>) {
        self.taxa = match taxa {
            None => None,
            Some(t) => Some(Arc::new(t)),
        };
    }
}

impl WeightedNode for Node {
    type Weight = f64;

    fn get_weight(&self) -> Option<Self::Weight> {
        self.weight.clone()
    }

    fn set_weight(&mut self, w: Option<Self::Weight>) {
        self.weight = w;
    }
}

impl Debug for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}",
            self.get_id(),
            self.node_type(),
            match self.get_taxa() {
                None => "None".to_string(),
                Some(t) => t.to_string(),
            },
            match self.get_weight() {
                None => "".to_string(),
                Some(t) => t.to_string(),
            }
        )
    }
}

impl Display for Node {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}:{}:{}",
            self.get_id(),
            self.node_type(),
            match self.get_taxa() {
                None => "None".to_string(),
                Some(t) => t.to_string(),
            },
            match self.get_weight() {
                None => "".to_string(),
                Some(t) => t.to_string(),
            }
        )
    }
}

impl Hash for Node {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.id.hash(state);
    }
}

impl PartialEq for Node {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Eq for Node {}
