/// Module with traits for tree traversals
pub mod node_iter;

use crate::{node::Node, prelude::*};
use itertools::Itertools;
use std::collections::VecDeque;

/// Struct for BFS iteration of nodes in a tree
#[derive(Clone)]
pub struct BFSIterator<'a,T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    stack: VecDeque<usize>,
    nodes: Vec<Option<&'a Node<T,W,Z>>>,
}

/// Struct for DFS postfix iteration of nodes in a tree
#[derive(Clone)]
pub struct DFSPostOrderIterator<'a,T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    stack: VecDeque<&'a Node<T,W,Z>>,
    nodes: Vec<Option<&'a Node<T,W,Z>>>,
}

impl<'a,T,W,Z> BFSIterator<'a,T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    /// Creates a new BFS iterator of a tree
    pub fn new(tree: &'a impl RootedTree<Node = Node<T,W,Z>>, start_id: usize) -> BFSIterator<'a,T,W,Z> {
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

impl<'a,T,W,Z> DFSPostOrderIterator<'a,T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    /// Creates a new DFS postfix iterator of a tree
    pub fn new(
        tree: &'a impl RootedTree<Node = Node<T,W,Z>>,
        start_id: usize,
    ) -> DFSPostOrderIterator<'a,T,W,Z> {
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

impl<'a,T,W,Z> Iterator for BFSIterator<'a,T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type Item = &'a Node<T,W,Z>;

    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop_front() {
            None => None,
            Some(curr_node_id) => {
                let curr_node = self.nodes[curr_node_id].take().unwrap();
                for &x in curr_node.get_children() {
                    self.stack.push_back(x);
                }
                Some(curr_node)
            }
        }
    }
}

impl<'a,T,W,Z> Iterator for DFSPostOrderIterator<'a,T,W,Z> 
where 
    T: NodeTaxa,
    W: EdgeWeight,
    Z: NodeWeight,
{
    type Item = &'a Node<T,W,Z>;

    fn next(&mut self) -> Option<Self::Item> {
        while let Some(node) = self.stack.pop_front() {
            let node_children: Vec<&'a Node<T,W,Z>> = node
                .get_children()
                .iter()
                .filter_map(|&chid| self.nodes[chid].take())
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
