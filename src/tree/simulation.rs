use crate::{
    node::simple_rnode::RootedMetaNode,
    tree::{RootedMetaTree, RootedTree},
};

/// A trait describing generation of a random binary tree under the Yule model.
pub trait Yule: RootedMetaTree
where
    <Self as RootedTree>::Node: RootedMetaNode,
{
    /// Generate a random binary tree under the Yule model with num_taxa
    fn yule(num_taxa: usize) -> Self;
}

/// A trait describing generation of a random binary tree under the Uniform model.
pub trait Uniform: RootedMetaTree
where
    <Self as RootedTree>::Node: RootedMetaNode,
{
    /// Generate a random binary tree under the Uniform model with num_taxa
    fn unif(num_taxa: usize) -> Self;
}
