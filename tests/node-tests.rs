use phylo::node::PhyloNode;

use phylo::node::simple_rnode::{RootedMetaNode, RootedTreeNode};

#[test]
fn test_set_id() {
    let mut n = PhyloNode::new(0);
    assert_eq!(n.get_id(), 0);
    n.set_id(10);
    assert_eq!(n.get_id(), 10);
}

#[test]
fn test_set_taxa() {
    let mut n = PhyloNode::new(0);
    assert_eq!(n.get_taxa(), None);
    n.set_taxa(Some(String::from("A")));
    assert_eq!(n.get_taxa(), Some(String::from("A")).as_ref());
}

#[test]
fn test_parent_childs() {
    let mut n = PhyloNode::new(0);
    n.add_child(10);
    n.add_child(20);
    assert_eq!(n.get_children(), &[10, 20]);
    n.remove_child(&10);
    assert_eq!(n.get_children(), &[20]);
    n.set_parent(Some(10));
    assert_eq!(n.get_parent(), Some(10));
}
