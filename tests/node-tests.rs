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

#[test]
fn bound_alias_traits_are_open() {
    // Form A: the weight/taxa bound aliases are blanket-implemented, so any type
    // meeting the bounds is a member -- not only the types this crate happens to
    // name. Under the previous sealed impls (only String and u32 for taxa) this
    // test would not compile: i64 was not a NodeTaxa, and a caller could not add
    // it either, since they own neither the trait nor the type.
    fn assert_taxa<T: phylo::prelude::NodeTaxa>() {}
    fn assert_edge_weight<W: phylo::prelude::EdgeWeight>() {}
    fn assert_node_weight<Z: phylo::prelude::NodeWeight>() {}

    // A taxa type this crate never names, now admitted by the bounds alone.
    assert_taxa::<i64>();
    // The types it does name still qualify -- via the bounds, not an explicit impl.
    assert_taxa::<String>();
    assert_edge_weight::<f32>();
    assert_edge_weight::<f64>();
    assert_node_weight::<f64>();
}
