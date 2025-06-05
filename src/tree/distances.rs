use itertools::Itertools;
use num::{Zero, Float, NumCast, One};
use std::fmt::Debug;
use vers_vecs::BitVec;

#[cfg(feature = "non_crypto_hash")]
use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
#[cfg(not(feature = "non_crypto_hash"))]
use std::collections::{HashMap, HashSet};

#[cfg(feature = "parallel")]
use rayon::prelude::*;

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
        let mut dist = 0;
        let mut all_taxa: HashSet<&TreeNodeMeta<Self>> = self.get_taxa_space().collect();
        all_taxa.extend(tree.get_taxa_space());
        let num_taxa = all_taxa.len();
        let all_taxa_map: HashMap<&TreeNodeMeta<Self>, usize> = all_taxa
            .into_iter()
            .enumerate()
            .map(|x| (x.1, x.0))
            .collect();


        let mut self_bps: HashMap<TreeNodeID<Self>, BitVec> = vec![].into_iter().collect();
        let mut self_out_bps: HashSet<BitVec> = vec![].into_iter().collect();
        for n_id in self.postord_ids(self.get_root_id()) {
            let mut bp = BitVec::from_zeros(num_taxa);
            match self.is_leaf(n_id) {
                true => {
                    let leaf_meta = self.get_node_taxa(n_id).unwrap();
                    bp.flip_bit(*all_taxa_map.get(leaf_meta).unwrap());
                    self_bps.insert(n_id, bp.clone());  
                }
                false => {
                    if n_id==self.get_root_id(){
                        let chids = self.get_node_children_ids(n_id).collect_vec();
                        let mut skip_root = false;
                        for chid in chids.iter(){
                            if self.is_leaf(*chid){
                                skip_root=true;
                            }
                        }
                        if skip_root{
                            break;
                        }
                        else{
                            let chid = chids[0];
                            let root_bp = self_bps.get(&chid).cloned().unwrap();
                            self_out_bps.insert(root_bp);
                        }
                        break;
                    }
                    self.get_node_children_ids(n_id)
                        .map(|x| self_bps.get(&x).unwrap())
                        .for_each(|x| {let _ = bp.apply_mask_or(x);});
                    if self.get_node_parent_id(n_id)!=Some(self.get_root_id()){
                        self_out_bps.insert(bp.clone());
                    }
                    self_bps.insert(n_id, bp);        
                }
            };
        }

        let mut tree_bps: HashMap<TreeNodeID<Self>, BitVec> = vec![].into_iter().collect();
        let mut tree_out_bps: HashSet<BitVec> = vec![].into_iter().collect();
        for n_id in tree.postord_ids(tree.get_root_id()) {
            let mut bp = BitVec::from_zeros(num_taxa);
            match tree.is_leaf(n_id) {
                true => {
                    let leaf_meta = tree.get_node_taxa(n_id).unwrap();
                    bp.flip_bit(*all_taxa_map.get(leaf_meta).unwrap());
                    tree_bps.insert(n_id, bp.clone());  
                }
                false => {
                    if n_id==tree.get_root_id(){
                        let chids = tree.get_node_children_ids(n_id).collect_vec();

                        let mut skip_root = false;
                        for chid in chids.iter(){
                            if tree.is_leaf(*chid){
                                skip_root=true;
                            }
                        }
                        if skip_root{
                            break;
                        }
                        else{
                            let chid = chids[0];
                            let root_bp = tree_bps.get(&chid).cloned().unwrap();
                            tree_out_bps.insert(root_bp);
                        }
                        break;
                    }
                    tree.get_node_children_ids(n_id)
                        .map(|x| tree_bps.get(&x).unwrap())
                        .for_each(|x| {let _ = bp.apply_mask_or(x);});
                    if tree.get_node_parent_id(n_id)!=Some(tree.get_root_id()){
                        tree_out_bps.insert(bp.clone());
                    }
                    tree_bps.insert(n_id, bp);
                }
            };
        }

        for i in self_out_bps.iter() {
            let mut i_rev = BitVec::from_ones(num_taxa);
            let _ = i_rev.apply_mask_xor(i);
            if !(tree_out_bps.contains(i) || tree_out_bps.contains(&i_rev)) {
                dist += 1;
            }
        }
        for i in tree_out_bps.iter() {
            let mut i_rev = BitVec::from_ones(num_taxa);
            let _ = i_rev.apply_mask_xor(i);
            if !(self_out_bps.contains(i) || self_out_bps.contains(&i_rev)) {
                dist += 1;
            }
        }
        
        dist
    }
}

/// A trait describing naive computation of Cluster Matching distance
pub trait ClusterMatching
where
    Self: RootedTree + RootedMetaTree + Clusters,
    <Self as RootedTree>::Node: RootedMetaNode,
{
    /// Returns Cluster Matching distance between tree and self.
    fn cm(&self, tree: &Self) -> usize {
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

/// A trait describing computation of Cluster Affinity
pub trait ClusterAffinity
where
    Self: RootedTree + RootedMetaTree + Clusters,
    <Self as RootedTree>::Node: RootedMetaNode,
{
    /// Returns Cluster Affinity cost from self to tree..
    fn ca(&self, tree: &Self) -> usize {
        let tree_clusters = tree
            .get_clusters_ids()
            .map(|(id, cluster)| {
                (
                    id,
                    cluster
                        .map(|x| tree.get_node_taxa(x).unwrap())
                        .collect::<HashSet<_>>(),
                )
            })
            .collect::<HashMap<TreeNodeID<Self>, HashSet<_>>>();

        let mut ca_cost = 0;

        for (_, cluster) in self.get_clusters_ids().map(|(id, cluster)| {
            (
                id,
                cluster
                    .map(|x| self.get_node_taxa(x).unwrap())
                    .collect::<HashSet<_>>(),
            )
        }) {
            ca_cost += tree_clusters
                .values()
                .map(|t_cluster| {
                    cluster.difference(t_cluster).collect_vec().len()
                        + t_cluster.difference(&cluster).collect_vec().len()
                })
                .min()
                .unwrap();
        }

        ca_cost
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
    TreeNodeZeta<Self>: NodeWeight,
{
    /// Returns zeta of leaf by taxa
    fn get_zeta_taxa(
        &self,
        taxa: &TreeNodeMeta<Self>,
    ) -> TreeNodeZeta<Self> {
        self.get_zeta(self.get_taxa_node_id(taxa).unwrap()).unwrap()
    }

    /// Reurns the nth norm of an iterator composed of floating point values
    fn compute_norm(
        vector: impl Iterator<Item = TreeNodeZeta<Self>>,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        if norm==0{
            return vector.fold(<TreeNodeZeta<Self>>::zero(), |acc, x| acc.max(x));
        }
        if norm == 1 {
            return vector.sum();
        }
        vector
            .map(|x| {
            let mut out = <TreeNodeZeta<Self>>::one();
                for _ in 0..norm{
                    out = out* x;
                }
                out
            })
            .sum::<TreeNodeZeta<Self>>()
            .powf(
                <TreeNodeZeta<Self> as NumCast>::from(norm)
                    .unwrap()
                    .powi(-1),
            )
    }

    #[cfg(feature = "parallel")]
    /// Returns the vector norm for an iterator
    fn compute_norm_par(
        vector: impl Iterator<Item = TreeNodeZeta<Self>>,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        if norm == 1 {
            return vector.map(|x| x.clone()).sum();
        }
        vector
            .map(|x| {
                x.clone().powi(norm as i32)
            })
            .sum::<TreeNodeZeta<Self>>()
            .powf(
                <TreeNodeZeta<Self> as NumCast>::from(norm)
                    .unwrap()
                    .powi(-1),
            )
    }

    /// Reurns the cophenetic distance between two trees using the naive algorithm (\Theta(n^2))
    fn cophen_dist<'a>(
        &'a self,
        tree: &'a Self,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        if !self.is_all_zeta_set() || !tree.is_all_zeta_set() {
            panic!("Zeta values not set");
        }
        let binding1 = self
            .get_taxa_space()
            .collect::<HashSet<&TreeNodeMeta<Self>>>();
        let binding2 = tree
            .get_taxa_space()
            .collect::<HashSet<&TreeNodeMeta<Self>>>();
        let taxa_set = binding1.intersection(&binding2).cloned();

        self.cophen_dist_by_taxa(tree, norm, taxa_set)
    }

    #[cfg(feature = "parallel")]
    /// Returns the cophenetic distance between two trees using the naive algorithm (\Theta(n^2))
    fn cophen_dist_par<'a>(
        &'a self,
        tree: &'a Self,
        norm: u32,
    ) -> TreeNodeZeta<Self> {
        if !self.is_all_zeta_set() || !tree.is_all_zeta_set() {
            panic!("Zeta values not set");
        }
        let binding1 = self
            .get_taxa_space()
            .collect::<HashSet<&TreeNodeMeta<Self>>>();
        let binding2 = tree
            .get_taxa_space()
            .collect::<HashSet<&TreeNodeMeta<Self>>>();
        let taxa_set = binding1.intersection(&binding2).cloned().map(|x| x.clone()).collect_vec();

        self.cophen_dist_by_taxa_par(tree, norm, taxa_set.iter())
    }

    #[cfg(feature = "parallel")]
    /// Returns the Cophenetic distance between two trees restricted to a taxa set using the \theta(n^2) naive algorithm.
    fn cophen_dist_by_taxa_par<'a>(
        &'a self,
        tree: &'a Self,
        norm: u32,
        taxa_set: impl Iterator<Item = &'a TreeNodeMeta<Self>> + Send,
    ) -> TreeNodeZeta<Self> {
        // let taxa_set = taxa_set.collect_vec();
        let cophen_vec = taxa_set
            .combinations_with_replacement(2)
            .into_iter()
            .map(|x| match x[0] == x[1] {
                true => {
                            let zeta_1 = self.get_zeta_taxa(x[0]);
                            let zeta_2 = tree.get_zeta_taxa(x[0]);
                            (zeta_1 - zeta_2).abs()
                        },
                false => {
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
                        },
            });
            // .collect::<Vec<TreeNodeZeta<Self>>>();

        Self::compute_norm_par(cophen_vec, norm)
    }

    /// Returns the Cophenetic distance between two trees restricted to a taxa set using the \theta(n^2) naive algorithm.
    fn cophen_dist_by_taxa<'a>(
        &'a self,
        tree: &'a Self,
        norm: u32,
        taxa_set: impl Iterator<Item = &'a TreeNodeMeta<Self>>,
    ) -> TreeNodeZeta<Self> {
        // let taxa_set = taxa_set.collect_vec();
        let cophen_vec = taxa_set
            // .iter()
            .combinations_with_replacement(2)
            .map(|x| match x[0] == x[1] {
                true => {
                            let zeta_1 = self.get_zeta_taxa(x[0]);
                            let zeta_2 = tree.get_zeta_taxa(x[0]);
                            (zeta_1 - zeta_2).abs()
                        },
                false => {
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
                        },
            });

        Self::compute_norm(cophen_vec, norm)
    }
}
