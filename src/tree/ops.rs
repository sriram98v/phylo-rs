use fxhash::FxHashSet as HashSet;
use fxhash::FxHashMap as HashMap;

use itertools::Itertools;
use num::{Float, NumCast, Signed, Zero};
use std::{fmt::{Debug, Display}, hash::Hash};
use std::time::Instant;
use rayon::prelude::*;


use crate::{iter::node_iter::{Ancestors, BFS}, node::simple_rnode::{RootedMetaNode, RootedTreeNode}, tree::simple_rtree::RootedMetaTree};
use super::{distances::PathFunction, Clusters, EulerWalk, RootedTree, RootedZetaNode, DFS};


pub trait SPR: RootedTree + DFS + Sized
{
    /// Attaches input tree to self by spliting an edge
    fn graft(&mut self, tree: Self, edge: (Self::NodeID, Self::NodeID));

    /// Returns subtree starting at given node, while corresponding nodes from self.
    fn prune(&mut self, node_id: Self::NodeID)-> Self;

    /// SPR function
    fn spr(&mut self, edge1: (Self::NodeID, Self::NodeID), edge2: (Self::NodeID, Self::NodeID))
    {
        let pruned_tree = SPR::prune(self, edge1.1);
        SPR::graft(self, pruned_tree, edge2);
    }
}

pub trait NNI
where
    Self: RootedTree + Sized,
{
    fn nni(&mut self, parent_id: Self::NodeID);
}

pub trait Reroot
where
    Self: RootedTree + Sized,
{
    fn reroot_at_node(&mut self, node_id: Self::NodeID);
    fn reroot_at_edge(&mut self, edge: (Self::NodeID, Self::NodeID));
}

pub trait Balance: Clusters + SPR + Sized
where
    Self::NodeID: Display + Debug + Hash + Clone + Ord,
{
    fn balance_subtree(&mut self);
}

pub trait Subtree: Ancestors + DFS + Sized
where
    Self::NodeID: Display + Debug + Hash + Clone + Ord,
{
    fn induce_tree(&self, node_id_list: impl IntoIterator<Item=Self::NodeID, IntoIter = impl ExactSizeIterator<Item = Self::NodeID>>)->Self
    {
        let mut subtree = self.clone();
        subtree.clear();
        for node_id in node_id_list.into_iter()
        {
            let ancestors = self.root_to_node(node_id);
            subtree.set_nodes(ancestors);
        }
        subtree.clean();
        subtree
    }
    fn subtree(&self, node_id: Self::NodeID)->Self
    {
        let mut subtree = self.clone();
        let dfs = self.dfs(node_id);
        subtree.set_nodes(dfs);
        subtree
    }
}

pub trait RobinsonFoulds
where
    Self: RootedTree + Sized,
    Self::NodeID: Display + Debug + Hash + Clone + Ord,
{
    fn rfs(&self, tree: Self)->usize;
}

pub trait ClusterAffinity
where
    Self: RootedTree + Sized,
    Self::NodeID: Display + Debug + Hash + Clone + Ord,
{
    fn ca(&self, tree: Self)->usize;
}

pub trait WeightedRobinsonFoulds
{
    fn wrfs(&self, tree: Self)->usize;
}

pub trait ContractTree: EulerWalk + DFS
{
    fn contracted_tree_nodes(&self, leaf_ids: &Vec<Self::NodeID>)-> impl ExactSizeIterator<Item=Self::Node>
    {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let node_postord_iter = self.postord(new_tree_root_id);
        let mut node_map: HashMap<Self::NodeID, Self::Node> = HashMap::from_iter(vec![(new_tree_root_id, self.get_lca(leaf_ids))].into_iter());
        let mut remove_list = vec![];
        node_postord_iter.for_each(|mut node| {
            match node.is_leaf(){
                true => {
                    match leaf_ids.contains(&node.get_id())
                    {
                        true => {node_map.insert(node.get_id(), node.clone());},
                        false => {},
                    }
                },
                false => {
                    let node_children_ids = node.get_children().collect_vec();
                    for child_id in &node_children_ids{
                        match node_map.get(&child_id).is_some(){
                            true => {},
                            false => node.remove_child(&child_id),
                        }
                    }
                    let node_children_ids = node.get_children().collect_vec();
                    match node_children_ids.len(){
                        0 => {},
                        1 => {
                            // the node is a unifurcation
                            // node should be added to both node_map and remove_list
                            // if child of node is already in remove list, attach node children to node first
                            let child_node_id = node_children_ids[0];

                            match remove_list.contains(&child_node_id){
                                true => {
                                    node.remove_child(&child_node_id);
                                    let grandchildren_ids = node_map.get(&child_node_id).unwrap().get_children().collect_vec();
                                    for grandchild_id in grandchildren_ids{
                                        node_map.get_mut(&grandchild_id).unwrap().set_parent(Some(node.get_id()));
                                        node.add_child(grandchild_id);
                                    }
                                },
                                false => {},
                            }
                            remove_list.push(node.get_id());
                            node_map.insert(node.get_id(), node);
                        },
                        _ => {
                            // node has multiple children
                            // for each child, suppress child if child is in remove list
                            node_children_ids.into_iter()
                                .for_each(|chid| {
                                    match remove_list.contains(&chid){
                                        true => {
                                            // suppress chid 
                                            // remove chid from node children
                                            // children of chid are node grandchildren
                                            // add grandchildren to node children
                                            // set grandchildren parent to node
                                            node.remove_child(&chid);
                                            let node_grandchildren = node_map.get(&chid).unwrap().get_children().collect_vec();
                                            for grandchild in node_grandchildren{
                                                node.add_child(grandchild);
                                                node_map.get_mut(&grandchild).unwrap().set_parent(Some(node.get_id()))
                                            }
                                        },
                                        false => {},
                                    }
                                });
                            if node.get_id()==new_tree_root_id{
                                node.set_parent(None);
                            }
                            node_map.insert(node.get_id(), node);
                        },
                    };
                },
            }
        });
        remove_list.into_iter()
            .for_each(|x| {node_map.remove(&x);});
        return node_map.into_values();
    }

    fn contract_tree(&self, leaf_ids: &Vec<Self::NodeID>)->Self
    {
        let new_tree_root_id = self.get_lca_id(leaf_ids);
        let new_nodes = self.contracted_tree_nodes(leaf_ids).collect_vec();
        let mut new_tree = self.clone();
        new_tree.set_root(new_tree_root_id);
        new_tree.clear();
        new_tree.set_nodes(new_nodes.into_iter());
        return new_tree;
    }
}

pub trait CopheneticDistance: RootedMetaTree<Meta=<Self as CopheneticDistance>::Meta> + EulerWalk + Clusters + Ancestors + ContractTree + PathFunction + Debug + Sized
where
    <Self as RootedTree>::Node: RootedMetaNode<Meta=<Self as CopheneticDistance>::Meta> + RootedZetaNode,
    <<Self as RootedTree>::Node as RootedZetaNode>::Zeta: Signed + Clone + NumCast + std::iter::Sum + Debug + Display + Float + PartialOrd + Copy + Send,
{
    type Meta: Display + Debug + Eq + PartialEq + Clone + Ord + Hash + Send + Sync;

    // helper functions

    /// Reurns the nth norm of an iterator composed of floating point values
    fn compute_norm(vector: impl Iterator<Item=<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, norm: u32)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        if norm==1{
            return vector.sum();
        }
        vector.map(|x| {
                x.powi(norm as i32)
            } )
            .sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>()
            .powf(<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(norm).unwrap().powi(-1))
    }

    /// Reurns the cophenetic distance between two trees using the naive algorithm (\Theta(n^2))
    fn cophen_dist_naive(&self, tree: &Self, norm: u32)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        if !self.is_all_zeta_set() || !tree.is_all_zeta_set(){
            panic!("Zeta values not set");
        }
        let binding1 = self.get_taxa_space().into_iter().collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let binding2 = tree.get_taxa_space().into_iter().collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let taxa_set = binding1.intersection(&binding2).map(|x| x.clone());

        return self.cophen_dist_naive_by_taxa(tree, norm, taxa_set);
    }

    /// Returns the Cophenetic distance between two trees restricted to a taxa set using the \theta(n^2) naive algorithm. 
    fn cophen_dist_naive_by_taxa(&self, tree: &Self, norm: u32, taxa_set: impl Iterator<Item=<Self as RootedMetaTree>::Meta>+Clone)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let dist = Self::compute_norm(taxa_set
            .combinations(2)
            .map(|x| {
                let t_lca_id = self.get_lca_id(&x.iter().map(|a| self.get_taxa_node_id(&a).unwrap()).collect_vec());
                let t_hat_lca_id = tree.get_lca_id(&x.iter().map(|a| tree.get_taxa_node_id(&a).unwrap()).collect_vec());
                let zeta_1 = self.get_zeta(t_lca_id).unwrap();
                let zeta_2 = tree.get_zeta(t_hat_lca_id).unwrap();
                return (zeta_1-zeta_2).abs()
                }), norm);
        return dist;
    }
}
