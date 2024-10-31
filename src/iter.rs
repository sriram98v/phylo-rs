/// Module with traits for tree traversals
pub mod node_iter;

use crate::{node::Node, prelude::*};
use itertools::Itertools;
use std::collections::VecDeque;

/// Struct for BFS iteration of nodes in a tree
#[derive(Clone)]
pub struct BFSIterator<'a> {
    stack: VecDeque<usize>,
    nodes: Vec<Option<&'a Node>>,
}

/// Struct for DFS postfix iteration of nodes in a tree
#[derive(Clone)]
pub struct DFSPostOrderIterator<'a> {
    stack: VecDeque<&'a Node>,
    nodes: Vec<Option<&'a Node>>,
}

impl<'a> BFSIterator<'a> {
    /// Creates a new BFS iterator of a tree
    pub fn new(tree: &'a impl RootedTree<Node = Node>, start_id: usize) -> BFSIterator<'a> {
        let max_id = tree.get_node_ids().max().unwrap();
        let mut nodes = vec![None; max_id + 1];
        tree.get_nodes()
            .for_each(|node| nodes[node.get_id()] = Some(node));
        BFSIterator {
            stack: vec![start_id].into(),
            nodes,
        }
    }
}

impl<'a> DFSPostOrderIterator<'a> {
    /// Creates a new DFS postfix iterator of a tree
    pub fn new(
        tree: &'a impl RootedTree<Node = Node>,
        start_id: <Node as RootedTreeNode>::NodeID,
    ) -> DFSPostOrderIterator<'a> {
        let max_id = tree.get_node_ids().max().unwrap();
        let mut nodes = vec![None; max_id + 1];
        tree.get_nodes()
            .for_each(|node| nodes[node.get_id()] = Some(node));
        let start_node = nodes[start_id].take().unwrap();
        DFSPostOrderIterator {
            stack: vec![start_node].into(),
            nodes,
        }
    }
}

impl<'a> Iterator for BFSIterator<'a> {
    type Item = &'a Node;

    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop_front() {
            None => None,
            Some(curr_node_id) => {
                let curr_node = self.nodes[curr_node_id].take().unwrap();
                curr_node
                    .get_children()
                    .for_each(|x| self.stack.push_back(x));
                Some(curr_node)
            }
        }
    }
}

impl<'a> Iterator for DFSPostOrderIterator<'a> {
    type Item = &'a Node;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(node) = self.stack.pop_front() {
            let node_children = node
                .get_children()
                .filter_map(|chid| self.nodes[chid].take())
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
