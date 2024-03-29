use phylo::node::Node;
use std::sync::Arc;

use itertools::Itertools;
use phylo::node::simple_rnode::{RootedTreeNode, RootedPhyloNode};

#[test]
fn test_set_id(){
let mut n = Node::new(Arc::new(0));
assert_eq!(n.get_id(),Arc::new(0));
n.set_id(Arc::new(10));
assert_eq!(n.get_id(),Arc::new(10));
}


#[test]
fn test_set_taxa(){
let mut n = Node::new(Arc::new(0));
assert_eq!(n.get_taxa(),None);
n.set_taxa(Some(String::from("A")));
assert_eq!(n.get_taxa(),Some(String::from("A")));
}

#[test]
fn test_parent_childs(){
let mut n = Node::new(0.into());
n.add_child(10.into());
n.add_child(20.into());
assert_eq!(n.get_children().into_iter().collect_vec(), vec![10.into(),20.into()]);
n.remove_child(&10.into());
assert_eq!(n.get_children().into_iter().collect_vec(), vec![20.into()]);
n.set_parent(Some(10.into()));
assert_eq!(n.get_parent(),Some(10.into()));

}