pub mod edge_iter;
pub mod node_iter;

use crate::{
    prelude::*,
    node::Node,
};
use itertools::Itertools;
use std::collections::VecDeque;

#[derive(Clone)]
pub struct BFSIterator {
    stack: VecDeque<usize>,
    nodes: Vec<Option<Node>>,
}

#[derive(Clone)]
pub struct DFSPostOrderIterator {
    stack: VecDeque<Node>,
    nodes: Vec<Option<Node>>,
}

impl BFSIterator {
    pub fn new(
        tree: &impl RootedTree<NodeID = usize, Node = Node>,
        start_id: usize,
    ) -> BFSIterator {
        let max_id = tree.get_node_ids().max().unwrap();
        let mut nodes = vec![None;max_id+1];
        tree.get_nodes()
            .for_each(|node| {
                nodes[node.get_id()] = Some(node.clone())
            });
        BFSIterator {
            stack: vec![start_id].into(),
            nodes,
        }
    }
}

impl DFSPostOrderIterator {
    pub fn new(
        tree: &impl RootedTree<NodeID = usize, Node = Node>,
        start_id: <Node as RootedTreeNode>::NodeID,
    ) -> DFSPostOrderIterator {
        let max_id = tree.get_node_ids().max().unwrap();
        let mut nodes = vec![None;max_id+1];
        tree.get_nodes()
        .for_each(|node| {
            nodes[node.get_id()] = Some(node.clone())
        });
        let start_node = std::mem::replace(&mut nodes[start_id], None).unwrap();
        DFSPostOrderIterator {
            stack: vec![start_node].into(),
            nodes,
        }
    }
}

impl Iterator for BFSIterator {
    type Item = Node;

    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop_front() {
            None => None,
            Some(curr_node_id) => {
                let curr_node = std::mem::replace(&mut self.nodes[curr_node_id], None).unwrap();
                curr_node
                    .get_children()
                    .for_each(|x| self.stack.push_back(x));
                Some(curr_node)
            }
        }
    }
}

impl Iterator for DFSPostOrderIterator {
    type Item = Node;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(node) = self.stack.pop_front() {
            let node_children = node
                .get_children()
                .filter_map(|chid| std::mem::replace(&mut self.nodes[chid], None))
                .collect_vec();
            if !node_children.is_empty() {
                self.stack.push_front(node);
                for child in node_children.into_iter().rev() {
                    self.stack.push_front(child)
                }
            } else {
                return Some(node);
            }
        }
        None
    }
}
