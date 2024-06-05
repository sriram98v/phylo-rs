use std::{collections::VecDeque, ops::Index};

use fxhash::FxHashMap as HashMap;
// use argminmax::ArgMinMax;
use itertools::Itertools;

use crate::{node::{simple_rnode::RootedTreeNode, Node}, tree::simple_rtree::RootedTree};

#[derive(Clone)]
pub struct BFSIterator{
    stack: VecDeque<<Node as RootedTreeNode>::NodeID>,
    nodes: HashMap<<Node as RootedTreeNode>::NodeID,Node>,
}

#[derive(Clone)]
pub struct DFSPostOrderIterator{
    stack: VecDeque<Node>,
    nodes: HashMap<<Node as RootedTreeNode>::NodeID,Node>,
}

impl BFSIterator
{
    pub fn new(tree: &impl RootedTree<NodeID = usize, Node = Node>, start_id: <Node as RootedTreeNode>::NodeID) -> BFSIterator
    {
        BFSIterator{
            stack: vec![start_id].into(),
            nodes: tree.get_nodes().into_iter().map(|x| (x.get_id(), x)).collect::<HashMap<_,_>>(),
        }
    }
}

impl DFSPostOrderIterator
{
    pub fn new(tree: &impl RootedTree<NodeID = usize, Node = Node>, start_id: <Node as RootedTreeNode>::NodeID) -> DFSPostOrderIterator
    {
        let mut nodes = tree.get_nodes().map(|x| (x.get_id(), x)).collect::<HashMap<_,_>>();
        let start_node = nodes.remove(&start_id).unwrap();
        DFSPostOrderIterator{
            stack: vec![start_node].into(),
            nodes: nodes,
        }
    }
}

impl Iterator for BFSIterator
{
    type Item = Node;
    
    fn next(&mut self) -> Option<Self::Item> {
        match self.stack.pop_front(){
            None => return None,
            Some(curr_node_id) => {
                let curr_node = self.nodes.remove(&curr_node_id).unwrap();
                curr_node.get_children()
                    .for_each(|x| {
                        self.stack.push_back(x)
                    });
                return Some(curr_node);

            },
        };
    }
}

impl Iterator for DFSPostOrderIterator
{
    type Item = Node;
    
    fn next(&mut self) -> Option<Self::Item> {

        while let Some(node) = self.stack.pop_front(){
            let node_children = node.get_children()
                .map(|chid| self.nodes.remove(&chid))
                .filter(|chn| chn.is_some())
                .map(|x| x.unwrap())
                .collect_vec();
            if node_children.len()>0{
                self.stack.push_front(node);
                for child in node_children.into_iter().rev(){
                    self.stack.push_front(child)
                }
            }
            else{
                return Some(node);
            }
        }
        return None;
    }
}


pub trait DFS
where
    Self: RootedTree + Sized
{
    fn postord(&self, start_node: Self::NodeID)-> impl Iterator<Item = Self::Node>;

    fn dfs(&self, start_node_id: Self::NodeID)-> impl IntoIterator<Item = Self::Node, IntoIter = impl ExactSizeIterator<Item = Self::Node>>
    {
        let mut stack = VecDeque::from([self.get_node(start_node_id).cloned().unwrap()]);
        let mut out_vec = vec![];
        let mut visited = vec![];
        while let Some(x) = stack.pop_front()
        {
            match visited.contains(&x.get_id())
            {
                false => {
                    visited.push(x.get_id());
                    out_vec.push(x.clone());
                    for child_id in x.get_children().into_iter().collect_vec().into_iter().rev(){
                        stack.push_front(self.get_node(child_id).cloned().unwrap());
                    }
                },
                true => {}
            };
        }
        out_vec
    }
}

pub trait BFS
where
    Self: RootedTree + Sized
{
    fn bfs(&self, start_node_id: Self::NodeID)-> impl Iterator<Item = Self::Node>;
}

pub trait PreOrder
where
    Self: RootedTree + Sized
{
    fn preord(&self, start_node_id: Self::NodeID)-> impl IntoIterator<Item = Self::Node, IntoIter = impl ExactSizeIterator<Item = Self::Node>>
    {
        let mut stack = VecDeque::from([self.get_node(start_node_id).cloned().unwrap()]);
        let mut out_vec = vec![];
        let mut visited = vec![];
        while let Some(x) = stack.pop_front()
        {
            match visited.contains(&x.get_id())
            {
                false => {
                    visited.push(x.get_id());
                    out_vec.push(x.clone());
                    for child_id in x.get_children().into_iter().collect_vec().into_iter().rev(){
                        stack.push_front(self.get_node(child_id).cloned().unwrap());
                    }
                },
                true => {}
            };
        }
        out_vec
    }
}


pub trait Ancestors
where
    Self: RootedTree + Sized
{
    fn root_to_node(&self, start_node_id: Self::NodeID)-> impl IntoIterator<Item = Self::Node, IntoIter = impl ExactSizeIterator<Item = Self::Node>>
    {
        let mut stack = VecDeque::from([self.get_node(start_node_id).cloned().unwrap()]);
        while let Some(x) = stack.pop_front(){
            stack.push_front(x.clone());
            match x.get_parent()
            {
                Some(pid) => {stack.push_front(self.get_node(pid).cloned().unwrap());},
                None => {break;},
            }
        }
        stack
    }

    fn node_to_root(&self, start_node_id: Self::NodeID)-> impl IntoIterator<Item = Self::Node, IntoIter = impl ExactSizeIterator<Item = Self::Node>>
    {
        let mut stack = VecDeque::from([self.get_node(start_node_id).cloned().unwrap()]);
        while let Some(x) = stack.pop_front(){
            stack.push_back(x.clone());
            match x.get_parent()
            {
                Some(pid) => {stack.push_front(self.get_node(pid).cloned().unwrap());},
                None => {break;},
            }
        }
        stack
    }

    fn depth(&self, node_id: Self::NodeID)->usize
    {
        self.node_to_root(node_id).into_iter().len()-1
    }
}

pub trait EulerWalk
where
    Self: RootedTree + Sized
{
    fn precompute_walk(&mut self);

    fn get_precomputed_walk(&self)->Option<&Vec<<Self as RootedTree>::NodeID>>;

    fn precompute_fai(&mut self);

    fn get_precomputed_fai(&self)->Option<impl Index<&Self::NodeID, Output = usize>>;

    fn precompute_da(&mut self);

    fn get_precomputed_da(&self)->Option<&Vec<usize>>;

    fn precompute_constant_time_lca(&mut self)
    {
        self.precompute_walk();
        self.precompute_da();
        self.precompute_fai();
    }

    fn euler_walk(&self, start_node_id: Self::NodeID)->impl ExactSizeIterator<Item = Self::Node>
    {
        let mut stack = VecDeque::from([self.get_node(start_node_id).unwrap()]);
        let mut visited = vec![];
        let mut out_vec = vec![];
        while let Some(node) = stack.pop_front()
        {
            match visited.contains(&node.get_id())
            {
                true => {
                    if let Some(parent_id) = node.get_parent() {out_vec.push(self.get_node(parent_id).cloned().unwrap())}
                },
                false => {
                    visited.push(node.get_id());
                    out_vec.push(node.clone());
                    stack.push_front(node);
                    for child_id in node.get_children().rev(){
                        stack.push_front(self.get_node(child_id).unwrap())
                    }
                }
            }
        }
        out_vec.into_iter()
    }

    fn is_euler_precomputed(&self)->bool
    {
        self.get_precomputed_walk().is_some()
    }

    /// returns the first appearance index of each node in the euler tour
    fn first_appearance(&self)->impl Index<&Self::NodeID, Output = usize>
    {
        let mut fa: HashMap<Self::NodeID, usize> = [].into_iter().collect::<HashMap<_,_>>();
        match self.get_precomputed_walk(){
            Some(walk) => {
                for (idx, node_id) in walk.iter().enumerate(){
                    if fa.contains_key(node_id){}
                    else {
                        fa.insert(node_id.clone(), idx);
                    }
                }
            },
            None => {
                let walk = self.euler_walk(self.get_root_id()).into_iter().map(|x| x.get_id());
                for (idx, node_id) in walk.enumerate(){
                    if fa.contains_key(&node_id){}
                    else {
                        fa.insert(node_id.clone(), idx);
                    }
                }
            }
        }

        fa
    }

    fn is_fai_precomputed(&self)->bool
    {
        self.get_precomputed_fai().is_some()
    }

    fn get_fa_index(&self, node_id: &Self::NodeID)->usize
    {
        match self.get_precomputed_fai(){
            Some(fai) => fai[node_id],
            None => self.first_appearance()[node_id],
        }

    }

    /// depth array for nodes
    fn depth_array(&self)-> Vec<usize>
    {
        let da = match self.get_precomputed_walk() {
            Some(walk) => walk.iter().map(|x| RootedTree::get_node_depth(self, x.clone())).collect_vec(),
            None => {
                self.euler_walk(self.get_root_id()).into_iter().map(|x| RootedTree::get_node_depth(self,x.get_id())).collect_vec()
            }
        };
        da
    }

    fn is_da_precomputed(&self)->bool
    {
        self.get_precomputed_da().is_some()
    }

    fn get_node_depth(&self, node_id: Self::NodeID)->usize
    {
        match self.get_precomputed_da(){
            Some(da) => return da[self.get_fa_index(&node_id)].clone(),
            None => {
                return RootedTree::get_node_depth(self, node_id);
            },
        };
    }

    fn get_euler_slice(&self, start: usize, end: usize)->Vec<<Self as RootedTree>::NodeID>
    {
        match self.get_precomputed_walk(){
            Some(walk) => return walk[start..end].to_vec(),
            None => self.euler_walk(self.get_root_id()).into_iter().map(|x| x.get_id()).collect_vec()[start..end].to_vec()
        }
    }

    fn get_da_slice(&self, start: usize, end: usize)->Vec<usize>
    {
        match self.get_precomputed_da(){
            Some(da) => return da[start..end].to_vec(),
            None => return self.depth_array()[start..end].to_vec(),
        };
    }

    fn get_euler_pos(&self, pos: usize)->Self::NodeID
    {
        match self.get_precomputed_walk(){
            Some(walk) => return walk[pos].clone(),
            None => self.euler_walk(self.get_root_id()).into_iter().map(|x| x.get_id()).nth(pos).unwrap()
        }
    }

    /// constant time lca query
    fn get_lca_id(&self, node_id_vec: &Vec<Self::NodeID>)-> Self::NodeID
    {
        let min_pos = node_id_vec.iter().map(|x| self.get_fa_index(&x)).min().unwrap();
        let max_pos = node_id_vec.iter().map(|x| self.get_fa_index(&x)).max().unwrap();

        let depth_subarray_min_value = self.get_da_slice(min_pos.clone(), max_pos.clone()).into_iter().min().unwrap();
        let depth_subarray_min_pos = self.get_da_slice(min_pos, max_pos).into_iter().position(|x| x==depth_subarray_min_value).unwrap();
        self.get_euler_pos(min_pos+depth_subarray_min_pos).clone()
    }

    fn get_lca(&self, node_id_vec: &Vec<Self::NodeID>)-> Self::Node
    {
        self.get_node(self.get_lca_id(node_id_vec)).cloned().unwrap()
    }
}

pub trait Clusters: DFS + BFS + Sized
{
    fn get_cluster(&self, node_id: Self::NodeID)-> impl ExactSizeIterator<Item = Self::Node>
    {
        self.dfs(node_id).into_iter().filter(|x| x.is_leaf()).collect_vec().into_iter()
    }
    fn get_cluster_ids(&self, node_id: Self::NodeID)-> impl ExactSizeIterator<Item = Self::NodeID>
    {
        self.get_cluster(node_id).map(|x| x.get_id())
    }
    fn get_cluster_size(&self, node_id: Self::NodeID)->usize
    {
        self.get_cluster(node_id).len()
    }

    fn get_bipartition(&self, edge: (Self::NodeID, Self::NodeID))->(impl ExactSizeIterator<Item = Self::Node>, impl ExactSizeIterator<Item = Self::Node>)
    {
        let c2 = self.get_cluster(edge.1.clone());
        let c2_ids = self.get_cluster(edge.1).into_iter().map(|x| x.get_id()).collect_vec();
        let c1 = self.get_cluster(edge.0).into_iter().filter(|x| !c2_ids.contains(&x.get_id())).collect_vec();
        (c1.into_iter(),c2.into_iter())
    }

    fn get_median_node_id_for_leaves(&self, taxa_set : impl ExactSizeIterator<Item=Self::NodeID>)->Self::NodeID
    {
        let mut median_node_id: Self::NodeID = self.get_root_id();
        let leaf_ids = taxa_set.collect_vec();
        loop {
            median_node_id = self.get_node_children_ids(median_node_id)
                .filter(|x| !self.is_leaf(x))
                .max_by(|x, y| {
                    let x_cluster = self.get_cluster_ids(x.clone()).filter(|nids| leaf_ids.contains(nids)).collect_vec();
                    let y_cluster = self.get_cluster_ids(y.clone()).filter(|nids| leaf_ids.contains(nids)).collect_vec();
                    return x_cluster.len().cmp(&y_cluster.len());
                })
                .unwrap();
            let num_leaves = self.get_cluster_ids(median_node_id.clone()).filter(|nids| leaf_ids.contains(nids)).collect_vec().len();
            if (num_leaves)<=(leaf_ids.len()/2){
                break;
            }
        }
        median_node_id
    }

    fn get_median_node_for_leaves(&self, taxa_set : impl ExactSizeIterator<Item=Self::NodeID>)->Self::Node
    {
        self.get_node(self.get_median_node_id_for_leaves(taxa_set)).cloned().unwrap()
    }

    fn get_median_node(&self)->Self::Node
    {
        let leaves = self.get_leaves().map(|x| x.get_id());
        self.get_median_node_for_leaves(leaves)
    }

    fn get_median_node_id(&self)->Self::NodeID
    {
        self.get_median_node().get_id()
    }
}