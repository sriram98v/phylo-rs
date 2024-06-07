use crate::node::simple_rnode::{RootedWeightedNode, RootedTreeNode};
use crate::tree::RootedTree;

pub trait PhylogeneticDiversity: RootedTree
where 
    <Self as RootedTree>::Node: RootedWeightedNode
{
    fn precompute_minPDs(&mut self);

    fn compute_norm_min(&self)->(Vec<Vec<(f32, u32)>>,Vec<Vec<(f32, u32)>>);

    fn get_minPD(&self, num_taxa: usize)-><<Self as RootedTree>::Node as RootedWeightedNode>::Weight;
    fn get_norm_minPD(&self, num_taxa: usize)-><<Self as RootedTree>::Node as RootedWeightedNode>::Weight;

    fn backtrack_min(&self, node_id: <<Self as RootedTree>::Node as RootedTreeNode>::NodeID, num_taxa: usize, taxaset: &mut Vec<<<Self as RootedTree>::Node as RootedTreeNode>::NodeID>);

    fn get_minPD_taxa_set_node(&self, node_id: <<Self as RootedTree>::Node as RootedTreeNode>::NodeID, num_taxa: usize)->impl Iterator<Item = <<Self as RootedTree>::Node as RootedTreeNode>::NodeID>
    {
        let mut taxa_set: Vec<<<Self as RootedTree>::Node as RootedTreeNode>::NodeID> = vec![];
        self.backtrack_min(node_id, num_taxa, &mut taxa_set);
        return taxa_set.into_iter()
    }

    fn get_minPD_taxa_set(&self, num_taxa: usize)->impl Iterator<Item = <<Self as RootedTree>::Node as RootedTreeNode>::NodeID>
    {
        let mut taxa_set: Vec<<<Self as RootedTree>::Node as RootedTreeNode>::NodeID> = vec![];
        self.backtrack_min(self.get_root_id(), num_taxa, &mut taxa_set);
        return taxa_set.into_iter()
    }

    fn get_minPD_node(&self, node_id: <<Self as RootedTree>::Node as RootedTreeNode>::NodeID, num_taxa: usize)-><<Self as RootedTree>::Node as RootedWeightedNode>::Weight;
}