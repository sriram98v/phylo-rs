use crate::node::simple_rnode::RootedWeightedNode;
use crate::tree::RootedTree;

pub trait PhylogeneticDiversity: RootedTree
where 
    <Self as RootedTree>::Node: RootedWeightedNode
{
    fn precompute_PDs(&mut self);
    fn get_minPD(&self, num_taxa: usize)-><<Self as RootedTree>::Node as RootedWeightedNode>::Weight;
    fn get_norm_minPD(&self, num_taxa: usize)-><<Self as RootedTree>::Node as RootedWeightedNode>::Weight;
}