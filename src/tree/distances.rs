use fxhash::FxHashSet as HashSet;
use itertools::Itertools;
use num::{Float, NumCast, Signed};
use std::fmt::{Debug, Display};

use crate::prelude::*;

/// A type alias for zeta annotation of a node in a tree.
pub type TreeNodeZeta<T> = <<T as RootedTree>::Node as RootedZetaNode>::Zeta;

/// A trait describing the path functions in a tree.
pub trait PathFunction: RootedTree
where
    <Self as RootedTree>::Node: RootedZetaNode,
{
    /// Sets zeta value of all nodes in a tree by function
    fn set_zeta(
        &mut self,
        zeta_func: fn(&Self, TreeNodeID<Self>) -> TreeNodeZeta<Self>,
    ) -> Option<()> {
        let node_ids = self.get_node_ids().collect_vec();
        for node_id in node_ids {
            let zeta = zeta_func(self, node_id);
            self.set_node_zeta(node_id, Some(zeta))?;
        }
        Some(())
    }

    /// Returns zeta value of a node in a tree. None is no zeta value is set
    fn get_zeta(&self, node_id: TreeNodeID<Self>) -> Option<TreeNodeZeta<Self>> {
        self.get_node(node_id)?.get_zeta()
    }

    /// Returns true if node zeta value is not None
    fn is_zeta_set(&self, node_id: TreeNodeID<Self>) -> bool {
        self.get_node(node_id).unwrap().is_zeta_set()
    }

    /// Returns true if all node zeta value is not None    
    fn is_all_zeta_set(&self) -> bool {
        !self.get_nodes().any(|x| !x.is_zeta_set())
    }

    /// Sets zeta value of a node in a tree by value.
    fn set_node_zeta(
        &mut self,
        node_id: TreeNodeID<Self>,
        zeta: Option<TreeNodeZeta<Self>>,
    ) -> Option<()> {
        self.get_node_mut(node_id)?.set_zeta(zeta);
        Some(())
    }
}

/// A trait describing efficient computation of vertex to vertex distances
pub trait DistanceMatrix: RootedWeightedTree
where
    <Self as RootedTree>::Node: RootedWeightedNode,
{
    /// Return the symmetrical pairwise distance matrix.
    fn matrix(&self) -> Vec<Vec<TreeNodeWeight<Self>>>;

    /// Populates a pairwise distance matrix for a subtree
    fn pairwise_distance(
        &self,
        node_id_1: TreeNodeID<Self>,
        node_id_2: TreeNodeID<Self>,
    ) -> TreeNodeWeight<Self>;
}

/// A trait describing naive computation of Robinson Foulds distance
pub trait RobinsonFoulds
where
    Self: RootedTree + RootedMetaTree + Clusters,
    <Self as RootedTree>::Node: RootedMetaNode,
{
    /// Returns Robinson Foulds distance between tree and self.
    fn rfs(&self, tree: &Self) -> usize {
        let self_bps = self
            .get_node_ids()
            .filter(|node_id| node_id == &self.get_root_id())
            .map(|node_id| {
                let node_parent_id = self.get_node_parent_id(node_id).unwrap();
                let bp_ids = self.get_bipartition_ids((node_parent_id, node_id));
                let p0 = bp_ids
                    .0
                    .map(|id| self.get_node_taxa(id).cloned().unwrap())
                    .collect::<HashSet<TreeNodeMeta<Self>>>();
                let p1 = bp_ids
                    .1
                    .map(|id| self.get_node_taxa(id).cloned().unwrap())
                    .collect::<HashSet<TreeNodeMeta<Self>>>();
                (p0, p1)
            })
            .collect_vec();
        let tree_bps = tree
            .get_node_ids()
            .filter(|node_id| node_id == &tree.get_root_id())
            .map(|node_id| {
                let node_parent_id = tree.get_node_parent_id(node_id).unwrap();
                let bp_ids = tree.get_bipartition_ids((node_parent_id, node_id));
                let p0 = bp_ids
                    .0
                    .map(|id| tree.get_node_taxa(id).cloned().unwrap())
                    .collect::<HashSet<TreeNodeMeta<Self>>>();
                let p1 = bp_ids
                    .1
                    .map(|id| tree.get_node_taxa(id).cloned().unwrap())
                    .collect::<HashSet<TreeNodeMeta<Self>>>();
                (p0, p1)
            })
            .collect_vec();

        let mut dist = 0;
        for i in self_bps.iter() {
            if tree_bps.contains(&(i.0.clone(), i.1.clone()))
                || tree_bps.contains(&(i.1.clone(), i.0.clone()))
            {
                dist += 1;
            }
        }
        for i in tree_bps {
            if self_bps.contains(&(i.0.clone(), i.1.clone())) || self_bps.contains(&(i.1, i.0)) {
                dist += 1;
            }
        }

        dist
    }
}

/// A trait describing naive computation of Cluster Affinity distance
pub trait ClusterAffinity
where
    Self: RootedTree + RootedMetaTree + Clusters,
    <Self as RootedTree>::Node: RootedMetaNode,
{
    /// Returns Cluster Affinity distance between tree and self.
    fn ca(&self, tree: &Self) -> usize {
        let self_clusters = self
            .get_node_ids()
            .map(|node_id| {
                self.get_cluster_ids(node_id)
                    .map(|id| self.get_node_taxa(id).unwrap())
                    .collect_vec()
            })
            .collect::<HashSet<_>>();
        let tree_clusters = tree
            .get_node_ids()
            .map(|node_id| {
                tree.get_cluster_ids(node_id)
                    .map(|id| tree.get_node_taxa(id).unwrap())
                    .collect_vec()
            })
            .collect::<HashSet<_>>();

        self_clusters.difference(&tree_clusters).collect_vec().len()
            + tree_clusters.difference(&self_clusters).collect_vec().len()
    }
}

/// A trait describing naive computation of Weighted Robinson Foulds distance
pub trait WeightedRobinsonFoulds
where
    Self: RootedWeightedTree + RootedMetaTree + Clusters,
    <Self as RootedTree>::Node: RootedWeightedNode + RootedMetaNode,
{
    /// Returns weighted Robinson Foulds distance between tree and self.
    fn wrfs(&self, tree: &Self) -> TreeNodeWeight<Self>;
}

/// A trait describing naive computation of cophenetic distance
pub trait CopheneticDistance:
    PathFunction + RootedMetaTree + Clusters + Ancestors + ContractTree + Debug
where
    <Self as RootedTree>::Node: RootedMetaNode + RootedZetaNode,
    <<Self as RootedTree>::Node as RootedZetaNode>::Zeta: Signed
        + Clone
        + NumCast
        + std::iter::Sum
        + Debug
        + Display
        + Float
        + PartialOrd
        + Copy
        + Send,
{
    /// Returns zeta of leaf by taxa
    fn get_zeta_taxa(
        &self,
        taxa: &TreeNodeMeta<Self>,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        self.get_zeta(self.get_taxa_node_id(taxa).unwrap()).unwrap()
    }

    /// Reurns the nth norm of an iterator composed of floating point values
    fn compute_norm(
        vector: impl Iterator<Item = <<Self as RootedTree>::Node as RootedZetaNode>::Zeta>,
        norm: u32,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        if norm == 1 {
            return vector.sum();
        }
        vector
            .map(|x| x.powi(norm as i32))
            .sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>()
            .powf(
                <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(norm)
                    .unwrap()
                    .powi(-1),
            )
    }

    /// Reurns the cophenetic distance between two trees using the naive algorithm (\Theta(n^2))
    fn cophen_dist_naive<'a>(
        &'a self,
        tree: &'a Self,
        norm: u32,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        if !self.is_all_zeta_set() || !tree.is_all_zeta_set() {
            panic!("Zeta values not set");
        }
        let binding1 = self
            .get_taxa_space()
            .collect::<HashSet<&<<Self as RootedTree>::Node as RootedMetaNode>::Meta>>();
        let binding2 = tree
            .get_taxa_space()
            .collect::<HashSet<&<<Self as RootedTree>::Node as RootedMetaNode>::Meta>>();
        let taxa_set = binding1.intersection(&binding2).cloned();

        self.cophen_dist_naive_by_taxa(tree, norm, taxa_set)
    }

    /// Returns the Cophenetic distance between two trees restricted to a taxa set using the \theta(n^2) naive algorithm.
    fn cophen_dist_naive_by_taxa<'a>(
        &'a self,
        tree: &'a Self,
        norm: u32,
        taxa_set: impl Iterator<Item = &'a TreeNodeMeta<Self>> + Clone,
    ) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta {
        let taxa_set = taxa_set.collect_vec();
        let cophen_vec = taxa_set
            .iter()
            .combinations_with_replacement(2)
            .map(|x| match x[0] == x[1] {
                true => vec![x[0]],
                false => x,
            })
            .map(|x| match x.len() {
                1 => {
                    let zeta_1 = self.get_zeta_taxa(x[0]);
                    let zeta_2 = tree.get_zeta_taxa(x[0]);
                    (zeta_1 - zeta_2).abs()
                }
                _ => {
                    let self_ids = x
                        .iter()
                        .map(|a| self.get_taxa_node_id(a).unwrap())
                        .collect_vec();
                    let tree_ids = x
                        .iter()
                        .map(|a| tree.get_taxa_node_id(a).unwrap())
                        .collect_vec();
                    let t_lca_id = self.get_lca_id(self_ids.as_slice());
                    let t_hat_lca_id = tree.get_lca_id(tree_ids.as_slice());
                    let zeta_1 = self.get_zeta(t_lca_id).unwrap();
                    let zeta_2 = tree.get_zeta(t_hat_lca_id).unwrap();
                    (zeta_1 - zeta_2).abs()
                }
            })
            .collect_vec();

        Self::compute_norm(cophen_vec.into_iter(), norm)
    }
}
