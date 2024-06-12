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

    /// Returns taxa present in upper tree.
    fn upper_tree_taxa(&self)->impl Iterator<Item=<Self as CopheneticDistance>::Meta>
    {
        let lower_tree_taxa = self.lower_tree_taxa().collect::<HashSet<_>>();
        self.get_leaves().map(|x| x.get_taxa().unwrap()).filter(move |x| !lower_tree_taxa.contains(x))
    }

    /// Returns taxa present in lower tree.
    fn lower_tree_taxa(&self)->impl Iterator<Item=<Self as CopheneticDistance>::Meta>
    {
        let median_node = self.get_median_node_id();
        self.get_cluster(median_node.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap())
    }

    /// Returns lower tree.
    fn lower_tree(&self)->Self
    {
        let lower_tree_taxa = HashSet::from_iter(self.lower_tree_taxa());
        self.contract_tree(&lower_tree_taxa.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec())

    }

    /// Returns upper tree.
    fn upper_tree(&self)->Self
    {
        let upper_tree_taxa = HashSet::from_iter(self.upper_tree_taxa());
        self.contract_tree(&upper_tree_taxa.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec())
    }

    // Returns zeta of leaf by taxa
    fn get_zeta_taxa(&self, taxa: &<Self as CopheneticDistance>::Meta)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        self.get_zeta(self.get_taxa_node_id(taxa).unwrap()).unwrap()
    }

    /// Returns value of n choose k.
    fn n_choose_k(n: usize, k: usize)->i32
    {
        fn multiplicative_form(n: usize, k:usize)->i32
        {
            (1..k).map(|x| (n+1-x) as f32/(x as f32)).product::<f32>() as i32
        }
    
        match k <= (n as f32/2 as f32) as usize{
            true =>{multiplicative_form(n, k)},
            false =>{multiplicative_form(n-k, n-k)},
        }
    }

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

    /// Returns the Cophenetic distance between two trees in O(pnlog^2n) time. 
    fn cophen_dist(&self, tree: &Self, norm: usize)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let mut ops: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> = vec![];
        let binding1 = self.get_taxa_space().into_iter().collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let binding2 = tree.get_taxa_space().into_iter().collect::<HashSet<<Self as CopheneticDistance>::Meta>>();
        let taxa_set = binding1.intersection(&binding2).map(|x| x.clone()).collect_vec();

        // dbg!(&ops);
        self.populate_op_vec(tree, norm, taxa_set.clone(), &mut ops);
        // dbg!(&ops);

        return ops.into_iter().sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>()
            .powf(<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(norm).unwrap().powi(-1));

    }

    /// Populates vector with divide and conquer distances.
    fn populate_op_vec(&self, tree: &Self, norm: usize, taxa_set: Vec<<Self as RootedMetaTree>::Meta>, op_vec: &mut Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>)
    {
        let double_mix_distance = self.distance_double_mix_type(tree, norm);
        let single_mix_distance = self.distance_single_mix_type(tree, norm);

        op_vec.push(double_mix_distance);
        op_vec.push(single_mix_distance);

        if taxa_set.len()>2{
            let t = self.get_median_node_id();
            let t_hat = tree.get_median_node_id();
            
            let b: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(self.get_cluster(t).filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()).filter(|x| taxa_set.contains(x)));
            let b_hat: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(tree.get_cluster(t_hat).filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()).filter(|x| taxa_set.contains(x)));
    
            let a: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(self.get_taxa_space()).difference(&b).filter(|x| taxa_set.contains(x)).map(|x| x.clone()).collect();
            let a_hat: HashSet<<Self as CopheneticDistance>::Meta> = HashSet::from_iter(tree.get_taxa_space()).difference(&b_hat).filter(|x| taxa_set.contains(x)).map(|x| x.clone()).collect();
    
            let a_int_a_hat = a.intersection(&a_hat).map(|x| x.clone()).collect_vec();
            let a_int_b_hat = a.intersection(&b_hat).map(|x| x.clone()).collect_vec();
            let b_int_a_hat = b.intersection(&a_hat).map(|x| x.clone()).collect_vec();
            let b_int_b_hat = b.intersection(&b_hat).map(|x| x.clone()).collect_vec();
        
            if a_int_a_hat.len()>1{
                let self_tree = self.contract_tree(&a_int_a_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&a_int_a_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, a_int_a_hat, op_vec);
            }
    
            if a_int_b_hat.len()>1{
                let self_tree = self.contract_tree(&a_int_b_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&a_int_b_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, a_int_b_hat, op_vec);
            }
    
            if b_int_b_hat.len()>1{
                let self_tree = self.contract_tree(&b_int_b_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&b_int_b_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, b_int_b_hat, op_vec);
            }
    
            if b_int_a_hat.len()>1{
                let self_tree = self.contract_tree(&b_int_a_hat.iter().map(|x| self.get_taxa_node_id(x).unwrap()).collect_vec());
                let new_tree = tree.contract_tree(&b_int_a_hat.iter().map(|x| tree.get_taxa_node_id(x).unwrap()).collect_vec());
                self_tree.populate_op_vec(&new_tree, norm, b_int_a_hat, op_vec);
            }
        }
    }

    /// Returns ordered iterator used in double mix type cases
    fn get_cntr(&self, leaf_set: HashSet<Self::NodeID>)->Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>
    {
        // line 5 in algo 1
        let mut gamma: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta> = Vec::new();
        let median_node_id = self.get_median_node_id();
        // line 3 in algo 1
        let mut median_path = self.root_to_node(median_node_id.clone()).into_iter().map(|x| (x.get_id(), 0)).collect::<HashMap<_,_>>();
        for node_id in leaf_set{
            // line 4 in algo 1
            median_path.entry(self.get_lca_id(&vec![node_id, median_node_id.clone()])).and_modify(|x| *x+=1);
        }
        for node in self.root_to_node(median_node_id.clone()).into_iter(){
            let c = median_path.get(&node.get_id()).cloned().unwrap();
            for _ in 0..c{
                gamma.push(self.get_zeta(node.get_id()).unwrap())
            }
        }
        gamma
    }

    /// Returns seqPrd^p(\alpha,\beta) used for double mix type cases.
    fn seq_product(mut alpha: Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, mut beta:Vec<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>, norm: usize)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        if alpha.is_empty() || beta.is_empty(){
            return <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as Zero>::zero();
        }
        if alpha.last().unwrap()<&alpha[0]{
            alpha.reverse();
        }
        if beta.last().unwrap()<&beta[0]{
            beta.reverse();
        }
        if alpha.last().unwrap()>beta.last().unwrap(){
            std::mem::swap(&mut alpha, &mut beta);
        }
        let mut out = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as Zero>::zero();
        for a in &alpha{
            for b in &beta{
                out = out + (a.clone()-b.clone()).abs().powi(norm as i32);
            }
        }
        out

    }

    /// This method generates the distance contributed by all taxa pairs 
    /// that are present in different subtrees in both trees(raised to the p^{th} power).
    /// 
    /// This includes the following assignments: AB|A'B', AB|B'A'
    fn distance_double_mix_type(&self, tree: &Self, norm: usize)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let t = self.get_median_node_id();
        let t_hat = tree.get_median_node_id();
        
        let b = self.get_cluster(t.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()).collect::<HashSet<_>>();
        let b_hat = tree.get_cluster(t_hat.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()).collect::<HashSet<_>>();

        let a: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(self.get_cluster(self.get_root_id()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap())).difference(&b).map(|x| x.clone()).collect();
        let a_hat: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(tree.get_cluster(tree.get_root_id()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap())).difference(&b_hat).map(|x| x.clone()).collect();

        // AB|B'A'
        let a_int_b_hat: HashSet<<Self as RootedMetaTree>::Meta> = a.intersection(&b_hat).map(|x| x.clone()).collect();
        let b_int_a_hat: HashSet<<Self as RootedMetaTree>::Meta> = b.intersection(&a_hat).map(|x| x.clone()).collect();

        let alpha = self.get_cntr(a_int_b_hat.iter().map(|x| self.get_taxa_node_id(&x).unwrap()).collect::<HashSet<Self::NodeID>>());
        let beta = tree.get_cntr(b_int_a_hat.iter().map(|x| tree.get_taxa_node_id(&x).unwrap()).collect::<HashSet<Self::NodeID>>());

        // AB|A'B'
        let b_int_b_hat_len = b.intersection(&b_hat).map(|x| x.clone()).collect_vec().len();
        let dd2 = a.intersection(&a_hat).map(|x| x.clone())
            .map(|x| {
                let t_lca_id = self.get_lca_id(&vec![self.get_taxa_node_id(&x).unwrap(), t.clone()]);
                let t_hat_lca_id = tree.get_lca_id(&vec![tree.get_taxa_node_id(&x).unwrap(), t_hat.clone()]);
                let zeta_1 = self.get_zeta(t_lca_id).unwrap();
                let zeta_2 = tree.get_zeta(t_hat_lca_id).unwrap();
                return (zeta_1-zeta_2).abs().powi(norm as i32) * <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(b_int_b_hat_len).unwrap()
            })
            .sum::<<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>();

        return Self::seq_product(alpha, beta, norm) + dd2;
    }

    /// This method generates the distance contributed by all taxa pairs 
    /// that are present in the same subtree in exactly one of the two trees(raised to the p^{th} power).
    /// 
    /// This includes the following assignments: AA|A'B', AA|B'A', BB|A'B', BB|B'A', BA|B'B', BA|A'A', AB|B'B', AB|A'A'.
    fn distance_single_mix_type(&self, tree: &Self, norm: usize)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        if self.num_taxa()<=2{
            return <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta>::zero();
        }
        let d1 = Self::single_mix_xxxy(self, tree, norm);
        let d2 = Self::single_mix_xxxy(tree, self, norm);
        return d1+d2;
    }

    // this method solves AA|A'B'
    fn single_mix_xxxy(t1: &Self, t2: &Self, norm: usize)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let t = t1.get_median_node_id();
        let t_hat = t2.get_median_node_id();

        let self_leaves = t1.get_leaves().map(|x| x.get_taxa().unwrap()).collect::<HashSet<_>>();
        let tree_leaves = t2.get_leaves().map(|x| x.get_taxa().unwrap()).collect::<HashSet<_>>();
        
        let b: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(t1.get_cluster(t.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()));
        let b_hat: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(t2.get_cluster(t_hat.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()));

        let a: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(self_leaves).difference(&b).map(|x| x.clone()).collect();
        let a_hat: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(tree_leaves).difference(&b_hat).map(|x| x.clone()).collect();

        let lower_tree_nodes = t1.postord(t).map(|x| x.get_id()).collect_vec();
        let upper_tree_nodes = t1.postord(t1.get_root_id()).map(|x| x.get_id()).filter(|x| !lower_tree_nodes.contains(x));

        let a_int_b_hat: HashSet<<Self as RootedMetaTree>::Meta> = a.intersection(&b_hat).map(|x| x.clone()).collect();
        let a_int_a_hat: HashSet<<Self as RootedMetaTree>::Meta> = a.intersection(&a_hat).map(|x| x.clone()).collect();

        let mut kappa: HashMap<_,_> = t1.get_node_ids().map(|x| (x, <<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero())).collect();
        match norm%2{
            0 => {
                let self_upper = t1.contract_tree(&a.iter().map(|x| t1.get_taxa_node_id(x).unwrap()).collect_vec());
                // setting sigma to zeros; kappa already set.
                let mut sigma: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1])).collect();        
                for x in a_int_a_hat{
                    let x_node_id = t2.get_taxa_node_id(&x).unwrap();
                    let lca_x_t_hat = t2.get_lca_id(&vec![x_node_id.clone(), t_hat]);
                    let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                    sigma.insert(x_node_id, (0..norm+1).map(|l| beta.powi(l as i32)).collect_vec());
                }
                for v_id in upper_tree_nodes{
                    if v_id!=self_upper.get_root_id(){
                        // calculate v_sigma
                        let v_value = t1.get_node_children_ids(v_id.clone())
                            .map(|x| sigma.get(&x).cloned().unwrap())
                            .fold(vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1], |acc, x| {
                                acc.iter().zip(x).map(|(a,b)| *a+b).collect_vec()
                            });
                        sigma.insert(v_id, v_value);
                        // calculate v_kappa
                        // kappa.insert()
                    }
                }
            }
            _ => {
                let mut sigma_pos: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1])).collect();
                let mut sigma_neg: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1])).collect();        
                let mut delta: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero(); norm+1])).collect();
                for x in &a_int_a_hat{
                    // x node_id in t
                    let x_node_id_t = t1.get_taxa_node_id(&x).unwrap();
                    // x node_id in t_hat
                    let x_node_id_t_hat = t2.get_taxa_node_id(&x).unwrap();
                    let x_parent_id_t = t1.get_node_parent_id(x_node_id_t).unwrap();
                    let lca_x_t_hat = t2.get_lca_id(&vec![x_node_id_t_hat.clone(), t_hat]);
                    let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                    if beta <= t1.get_zeta_taxa(x){
                        // find omega_x
                        let omega_x = t1.node_to_root(x_node_id_t).into_iter()
                            .filter(|w| w.get_zeta().unwrap()<=beta)
                            .min_by(|w, y| {
                                w.get_zeta().unwrap().partial_cmp(&y.get_zeta().unwrap()).unwrap()
                            }).unwrap();
                        // set omega_x.delta
                        delta.entry(omega_x.get_id())
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l] = e[l]+beta.powi(l as i32);
                                }
                            });
                    }
                    if beta <= t1.get_zeta(x_parent_id_t).unwrap(){
                        // set x.sigma_plus
                        sigma_pos.entry(x_node_id_t)
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l] = e[l]+beta.powi(l as i32);
                                }
                            });
                    }
                    else{
                        // set x.sigma_minus
                        sigma_neg.entry(x_node_id_t)
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l] = e[l]+beta.powi(l as i32);
                                }
                            });
                    }
                }
                let self_upper = t1.contract_tree(&a.iter().map(|x| t1.get_taxa_node_id(x).unwrap()).collect_vec());
                for v in self_upper.postord(self_upper.get_root_id()){
                    if v.get_id()!=self_upper.get_root_id(){
                        let v_parent = v.get_parent().unwrap();
                        if !v.is_leaf(){
                            let v_children = v.get_children().collect_vec();
                            let v_left_child_id = v_children[0];
                            let v_right_child_id = v_children[1];
                            let v_delta = sigma_pos.get(&v.get_id()).cloned().unwrap();
                            // calculate v_sigma_pos
                            let v_left_sigma_pos = sigma_pos.get(&v_left_child_id).cloned().unwrap();
                            let v_right_sigma_pos = sigma_pos.get(&v_right_child_id).cloned().unwrap();
                            sigma_pos.entry(v.get_id())
                                .and_modify(|e| {
                                    for l in 0..norm+1{
                                        e[l] = v_left_sigma_pos[l]+v_right_sigma_pos[l]-v_delta[l];
                                    }
                                });
                            // calculate v_sigma_neg
                            let v_left_sigma_neg = sigma_neg.get(&v_left_child_id).cloned().unwrap();
                            let v_right_sigma_neg = sigma_neg.get(&v_right_child_id).cloned().unwrap();
                            sigma_neg.entry(v.get_id())
                                .and_modify(|e| {
                                    for l in 0..norm+1{
                                        e[l] = v_left_sigma_neg[l]+v_right_sigma_neg[l]+v_delta[l];
                                    }
                                });    
                        }
                        // calculate v_kappa
                        let v_sibling = t1.get_sibling_ids(v.get_id()).collect_vec()[0];
                        let v_sibling_leaves = t1.get_cluster_ids(v_sibling).map(|x| t1.get_node_taxa(x).unwrap()).collect::<HashSet<_>>();
                        let Y = a_int_b_hat.intersection(&v_sibling_leaves).collect::<HashSet<_>>();
                        let v_kappa_1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Y.len()).unwrap();
                        let v_kappa_2: <<Self as RootedTree>::Node as RootedZetaNode>::Zeta = (0..norm+1)
                            .map(|l| {
                                let term1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                                let term2 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(-1).unwrap().powi((norm-l) as i32);
                                let term3_1 = t1.get_zeta(v_parent).unwrap().powi(l as i32)*sigma_pos.get(&v.get_id()).unwrap()[norm-l];
                                let term3_2 = t1.get_zeta(v_parent).unwrap().powi((norm-l) as i32)*sigma_neg.get(&v.get_id()).unwrap()[l];
                                let term3 = term3_1+term3_2;
                                term1*term2*term3
                            })
                            .sum();
                        let v_kappa = v_kappa_1*v_kappa_2;
                        kappa.insert(v.get_id(), v_kappa);
                    }
                }
            }
        }
        return kappa.into_values().sum();
    }
    fn single_mix_yyyx(t1: &Self, t2: &Self, norm: usize)-><<Self as RootedTree>::Node as RootedZetaNode>::Zeta
    {
        let t = t1.get_median_node_id();
        let t_hat = t2.get_median_node_id();

        let self_leaves = t1.get_leaves().map(|x| x.get_taxa().unwrap()).collect::<HashSet<_>>();
        let tree_leaves = t2.get_leaves().map(|x| x.get_taxa().unwrap()).collect::<HashSet<_>>();
        
        let b: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(t1.get_cluster(t.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()));
        let b_hat: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(t2.get_cluster(t_hat.clone()).into_iter().filter(|x| x.get_taxa().is_some()).map(|x| x.get_taxa().unwrap()));

        let a: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(self_leaves).difference(&b).map(|x| x.clone()).collect();
        let a_hat: HashSet<<Self as RootedMetaTree>::Meta> = HashSet::from_iter(tree_leaves).difference(&b_hat).map(|x| x.clone()).collect();

        let lower_tree_nodes = t1.postord(t).map(|x| x.get_id()).collect_vec();
        let upper_tree_nodes = t1.postord(t1.get_root_id()).map(|x| x.get_id()).filter(|x| !lower_tree_nodes.contains(x));

        let a_int_b_hat: HashSet<<Self as RootedMetaTree>::Meta> = a.intersection(&b_hat).map(|x| x.clone()).collect();
        let b_int_b_hat: HashSet<<Self as RootedMetaTree>::Meta> = b.intersection(&b_hat).map(|x| x.clone()).collect();

        let mut kappa: HashMap<_,_> = t1.get_node_ids().map(|x| (x, <<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero())).collect();
        match norm%2{
            0 => {
                let self_upper = t1.contract_tree(&a.iter().map(|x| t1.get_taxa_node_id(x).unwrap()).collect_vec());
                // setting sigma to zeros; kappa already set.
                let mut sigma: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1])).collect();        
                for x in b_int_b_hat{
                    let x_node_id = t2.get_taxa_node_id(&x).unwrap();
                    let lca_x_t_hat = t2.get_lca_id(&vec![x_node_id.clone(), t_hat]);
                    let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                    sigma.insert(x_node_id, (0..norm+1).map(|l| beta.powi(l as i32)).collect_vec());
                }
                for v_id in upper_tree_nodes{
                    if v_id!=self_upper.get_root_id(){
                        // calculate v_sigma
                        let v_value = t1.get_node_children_ids(v_id.clone())
                            .map(|x| sigma.get(&x).cloned().unwrap())
                            .fold(vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1], |acc, x| {
                                acc.iter().zip(x).map(|(a,b)| *a+b).collect_vec()
                            });
                        sigma.insert(v_id, v_value);
                        // calculate v_kappa
                        // kappa.insert()
                    }
                }
            }
            _ => {
                let mut sigma_pos: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1])).collect();
                let mut sigma_neg: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero();norm+1])).collect();        
                let mut delta: HashMap<_,_> = t1.get_node_ids().map(|x| (x, vec![<<Self as RootedTree>::Node as RootedZetaNode>::Zeta::zero(); norm+1])).collect();
                for x in &b_int_b_hat{
                    // x node_id in t
                    let x_node_id_t = t1.get_taxa_node_id(&x).unwrap();
                    // x node_id in t_hat
                    let x_node_id_t_hat = t2.get_taxa_node_id(&x).unwrap();
                    let x_parent_id_t = t1.get_node_parent_id(x_node_id_t).unwrap();
                    let lca_x_t_hat = t2.get_lca_id(&vec![x_node_id_t_hat.clone(), t_hat]);
                    let beta = t2.get_zeta(lca_x_t_hat).unwrap();
                    if beta <= t1.get_zeta_taxa(x){
                        // find omega_x
                        let omega_x = t1.node_to_root(x_node_id_t).into_iter()
                            .filter(|w| w.get_zeta().unwrap()<=beta)
                            .min_by(|w, y| {
                                w.get_zeta().unwrap().partial_cmp(&y.get_zeta().unwrap()).unwrap()
                            }).unwrap();
                        // set omega_x.delta
                        delta.entry(omega_x.get_id())
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l] = e[l]+beta.powi(l as i32);
                                }
                            });
                    }
                    if beta <= t1.get_zeta(x_parent_id_t).unwrap(){
                        // set x.sigma_plus
                        sigma_pos.entry(x_node_id_t)
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l] = e[l]+beta.powi(l as i32);
                                }
                            });
                    }
                    else{
                        // set x.sigma_minus
                        sigma_neg.entry(x_node_id_t)
                            .and_modify(|e| {
                                for l in 0..norm+1{
                                    e[l] = e[l]+beta.powi(l as i32);
                                }
                            });
                    }
                }
                let self_upper = t1.contract_tree(&a.iter().map(|x| t1.get_taxa_node_id(x).unwrap()).collect_vec());
                for v in self_upper.postord(self_upper.get_root_id()){
                    if v.get_id()!=self_upper.get_root_id(){
                        let v_parent = v.get_parent().unwrap();
                        if !v.is_leaf(){
                            let v_children = v.get_children().collect_vec();
                            let v_left_child_id = v_children[0];
                            let v_right_child_id = v_children[1];
                            let v_delta = sigma_pos.get(&v.get_id()).cloned().unwrap();
                            // calculate v_sigma_pos
                            let v_left_sigma_pos = sigma_pos.get(&v_left_child_id).cloned().unwrap();
                            let v_right_sigma_pos = sigma_pos.get(&v_right_child_id).cloned().unwrap();
                            sigma_pos.entry(v.get_id())
                                .and_modify(|e| {
                                    for l in 0..norm+1{
                                        e[l] = v_left_sigma_pos[l]+v_right_sigma_pos[l]-v_delta[l];
                                    }
                                });
                            // calculate v_sigma_neg
                            let v_left_sigma_neg = sigma_neg.get(&v_left_child_id).cloned().unwrap();
                            let v_right_sigma_neg = sigma_neg.get(&v_right_child_id).cloned().unwrap();
                            sigma_neg.entry(v.get_id())
                                .and_modify(|e| {
                                    for l in 0..norm+1{
                                        e[l] = v_left_sigma_neg[l]+v_right_sigma_neg[l]+v_delta[l];
                                    }
                                });    
                        }
                        // calculate v_kappa
                        let v_sibling = t1.get_sibling_ids(v.get_id()).collect_vec()[0];
                        let v_sibling_leaves = t1.get_cluster_ids(v_sibling).map(|x| t1.get_node_taxa(x).unwrap()).collect::<HashSet<_>>();
                        let Y = a_int_b_hat.intersection(&v_sibling_leaves).collect::<HashSet<_>>();
                        let v_kappa_1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Y.len()).unwrap();
                        let v_kappa_2: <<Self as RootedTree>::Node as RootedZetaNode>::Zeta = (0..norm+1)
                            .map(|l| {
                                let term1 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Self::n_choose_k(norm, l)).unwrap();
                                let term2 = <<<Self as RootedTree>::Node as RootedZetaNode>::Zeta as NumCast>::from(Y.len()).unwrap().powi((norm-l) as i32);
                                let term3_1 = t1.get_zeta(v_parent).unwrap().powi(l as i32)*sigma_pos.get(&v.get_id()).unwrap()[norm-l];
                                let term3_2 = t1.get_zeta(v_parent).unwrap().powi((norm-l) as i32)*sigma_neg.get(&v.get_id()).unwrap()[l];
                                let term3 = term3_1+term3_2;
                                term1*term2*term3
                            })
                            .sum();
                        let v_kappa = v_kappa_1*v_kappa_2;
                        kappa.insert(v.get_id(), v_kappa);                        
                    }
                }
            }
        }
        return kappa.into_values().sum();
    }
}
