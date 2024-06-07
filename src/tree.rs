pub mod simple_rtree;
pub mod ops;
pub mod distances;
pub mod io;

use std::ops::Index;

// use std::collections::HashMap;
use fxhash::FxHashMap as HashMap;

use itertools::Itertools;
use anyhow::Result;
use rand::prelude::IteratorRandom;

use crate::node;
use crate::node::simple_rnode::*;
use crate::node::Node;
use crate::tree::{simple_rtree::*, io::*, ops::*};
use crate::iter::node_iter::*;
use vers_vecs::BinaryRmq;

use self::distances::PathFunction;

#[derive(Debug, Clone)]
pub struct SimpleRootedTree{
    root: usize,
    nodes: Vec<Option<Node>>,
    taxa_node_id_map: HashMap<String, usize>,
    precomputed_euler: Option<Vec<usize>>,
    precomputed_fai: Option<Vec<Option<usize>>>,    
    precomputed_da: Option<Vec<usize>>,
    precomputed_rmq: Option<BinaryRmq>,
}

impl SimpleRootedTree{

    pub fn new(root_id: usize)->Self{
        let root_node = Node::new(root_id);
        let mut nodes = vec![None;root_id+1];
        nodes[root_id] = Some(root_node);
        SimpleRootedTree { 
            root: root_id,
            nodes: nodes,
            precomputed_euler: None,
            taxa_node_id_map: [].into_iter().collect::<HashMap<_,_>>(),
            precomputed_fai: None,
            precomputed_da: None,
            precomputed_rmq: None,
        }
    }


    pub fn next_id(&self)->usize
    {
        match &self.nodes.iter().position(|r| r.is_none())
        {
            Some(x) => x.clone(),
            None => self.nodes.len()
        }
    }

    pub fn next_node(&self)->Node
    {
        Node::new(self.next_id())
    }

    pub fn yule(num_taxa: usize)->Result<SimpleRootedTree>
    {   
        let mut tree = SimpleRootedTree::new(0);
        if num_taxa<3{
            return Ok(tree);
        }
        let new_node = Node::new(1);
        tree.add_child(0, new_node);
        tree.set_node_taxa(1, Some("0".to_string()));
        let new_node = Node::new(2);
        tree.add_child(0, new_node);
        tree.set_node_taxa(2, Some("1".to_string()));
        if num_taxa<4{
            return Ok(tree);
        }
        let mut current_leaf_ids = vec![1,2];
        for i in 2..num_taxa{
            let rand_leaf_id = current_leaf_ids.iter().choose(&mut rand::thread_rng()).unwrap();
            let rand_leaf_parent_id = tree.get_node_parent_id(rand_leaf_id.clone()).unwrap();
            let split_node = Node::new(tree.next_id());
            let split_node_id = split_node.get_id();
            tree.split_edge((rand_leaf_parent_id, rand_leaf_id.clone()), split_node);
            let new_leaf = Node::new(tree.next_id());
            let new_leaf_id = new_leaf.get_id();
            tree.add_child(split_node_id, new_leaf);
            tree.set_node_taxa(new_leaf_id.clone(), Some(i.to_string()));
            current_leaf_ids.push(new_leaf_id);
        }
        Ok(tree)
    }
    pub fn unif(num_taxa: usize)->Result<SimpleRootedTree>
    {   
        let mut tree = SimpleRootedTree::new(0);
        if num_taxa<3{
            return Ok(tree);
        }
        let new_node = Node::new(1);
        tree.add_child(0, new_node);
        tree.set_node_taxa(1, Some("0".to_string()));
        let new_node = Node::new(2);
        tree.add_child(0, new_node);
        tree.set_node_taxa(2, Some("1".to_string()));
        if num_taxa<3{
            return Ok(tree);
        }
        let mut current_node_ids = vec![1,2];
        for i in 1..num_taxa{
            let rand_leaf_id = current_node_ids.iter().choose(&mut rand::thread_rng()).unwrap().clone();
            let rand_leaf_parent_id = tree.get_node_parent_id(rand_leaf_id.clone()).unwrap();
            let split_node = Node::new(tree.next_id());
            let split_node_id = split_node.get_id();
            current_node_ids.push(split_node_id.clone());
            tree.split_edge((rand_leaf_parent_id, rand_leaf_id.clone()), split_node);
            let mut new_leaf = Node::new(tree.next_id());
            new_leaf.set_taxa(Some(i.to_string()));
            current_node_ids.push(new_leaf.get_id());
            tree.add_child(split_node_id, new_leaf)
        }
        Ok(tree)
    }
}


impl RootedTree for SimpleRootedTree{

    type NodeID = usize;
    type Node = Node;
    
    /// Returns reference to node by ID
    fn get_node(&self, node_id: Self::NodeID)->Option<&Node>
    {
        self.nodes[node_id].as_ref()
    }

    fn get_node_mut(&mut self, node_id: Self::NodeID)->Option<&mut Node>
    {
        self.nodes[node_id].as_mut()
    }

    fn get_node_ids(&self)->impl Iterator<Item = Self::NodeID>
    {
        (0..self.nodes.len()).filter(|x| self.nodes[x.clone()].is_some())
    }

    fn get_nodes(&self)->impl ExactSizeIterator<Item = Node>
    {
        self.nodes.clone().into_iter().filter(|x| x.is_some()).map(|x| x.unwrap()).collect_vec().into_iter()
    }

    /// Returns reference to node by ID
    fn set_node(&mut self, node: Node)
    {
        let node_id = node.get_id();
        match node.get_taxa(){
            None => {},
            Some(taxa) => {self.taxa_node_id_map.insert(taxa, node.get_id());},
        };
        match self.nodes.len()>node_id{
            true => {self.nodes[node_id] = Some(node)},
            false => {
                let new_len = node.get_id()-self.nodes.len();
                self.nodes.extend((0..new_len+1).map(|_| None));
                self.nodes[node_id] = Some(node);
            }
        };
    }

    fn add_child(&mut self, parent_id: Self::NodeID, child: Node)
    {
        let new_child_id = child.get_id();
        self.set_node(child);
        self.get_node_mut(parent_id).unwrap().add_child(new_child_id);
        self.get_node_mut(new_child_id).unwrap().set_parent(Some(parent_id));
    }

    /// Get root node ID
    fn get_root_id(&self)->Self::NodeID
    {
        self.root
    }

    fn set_root(&mut self, node_id: Self::NodeID) 
    {
        self.root = node_id;
    }

    fn remove_node(&mut self, node_id: Self::NodeID)->Option<Node>
    {
        if let Some(pid) = self.get_node_parent_id(node_id) {self.get_node_mut(pid).unwrap().remove_child(&node_id)}
        std::mem::replace(&mut self.nodes[node_id], None)
    }

    fn delete_node(&mut self, node_id: Self::NodeID) {
        let _ = std::mem::replace(&mut self.nodes[node_id], None);
    }

    fn contains_node(&self, node_id: Self::NodeID)->bool
    {
        self.nodes[node_id].is_some()
    }

    fn delete_edge(&mut self, parent_id: Self::NodeID, child_id: Self::NodeID)
    {
        self.get_node_mut(parent_id).unwrap().remove_child(&child_id);
        self.get_node_mut(child_id).unwrap().set_parent(None);
    }

    fn clean(&mut self)
    {
        let node_iter = self.get_nodes().into_iter().collect_vec();
        for node in &node_iter{
            // remove root with only one child
            let node_id = node.get_id();
            if node.get_id()==self.get_root_id() && node.degree()<2{
                let new_root = self.get_root().get_children().into_iter().next().unwrap();
                self.set_root(new_root);
                self.get_node_mut(self.root).unwrap().set_parent(None);
                self.remove_node(node_id);
            }
            // remove nodes with only one child
            else if !node.is_leaf() && node.get_parent().is_some() && node.degree()<3 {
                let parent_id = self.get_node_parent_id(node_id);
                let child_id = node.get_children().into_iter().next().unwrap();
                self.get_node_mut(child_id).unwrap().set_parent(parent_id);
                self.get_node_mut(parent_id.unwrap()).unwrap().add_child(child_id);
                self.remove_node(node.get_id());
            }
            // Removing dangling references to pruned children
            for chid in node.get_children().into_iter()
            {
                if !self.get_nodes().into_iter().map(|x| x.get_id()).contains(&chid)
                {
                    self.get_node_mut(node_id).unwrap().remove_child(&chid);
                }
            }
        }
    }

    fn clear(&mut self) {
        let root_node = self.get_root().clone();
        let root_node_id = root_node.get_id();
        self.nodes = vec![None;root_node_id+1];
        self.nodes[root_node_id] = Some(root_node);
        self.taxa_node_id_map.clear();
        // self.precomputed_da = None;
        // self.precomputed_euler = None;
        // self.precomputed_fai = None;
        // self.precomputed_rmq = None;
    }
}


impl RootedMetaTree for SimpleRootedTree{

    type Meta = String;

    fn get_taxa_node_id(&self, taxa: &Self::Meta)->Option<Self::NodeID> {
        self.taxa_node_id_map.get(taxa).cloned()
    }

    fn get_taxa_node(&self, taxa: &Self::Meta)->Option<&Self::Node> {
        match self.taxa_node_id_map.get(taxa){
            None => None,
            Some(node_id) => {
                self.get_node(*node_id)
            }
        }
    }

    fn set_node_taxa(&mut self, node_id: Self::NodeID, taxa: Option<Self::Meta>) {
        self.get_node_mut(node_id).unwrap().set_taxa(taxa.clone());
        if let Some(t) = taxa {self.taxa_node_id_map.insert(t, node_id);}        
    }

    fn num_taxa(&self)->usize {
        self.taxa_node_id_map.len()
    }

    fn get_taxa_space(&self)->impl ExactSizeIterator<Item=Self::Meta> + Clone {
        self.taxa_node_id_map.keys().cloned()
    }

}

impl RootedWeightedTree for SimpleRootedTree {

    type Weight = f32;

}

impl PathFunction for SimpleRootedTree{
    fn get_zeta(&self, node_id: Self::NodeID)->Option<f32> {
        self.get_node(node_id).expect(format!("Node with ID {} not found in tree", node_id).as_str()).get_zeta()
    }
    fn set_zeta(&mut self, zeta_func: fn(&Self, Self::NodeID) -> <<Self as RootedTree>::Node as RootedZetaNode>::Zeta) {
        let zetas = self.get_node_ids().map(|node_id| (node_id, Some(zeta_func(self, node_id)))).collect_vec();
        for (node_id, zeta) in zetas{
            self.get_node_mut(node_id).unwrap().set_zeta(zeta)
        }
    }
}

impl Ancestors for SimpleRootedTree{}

impl Subtree for SimpleRootedTree{}

impl PreOrder for SimpleRootedTree{}

impl DFS for SimpleRootedTree{
    fn postord(&self, start_node: Self::NodeID)-> impl Iterator<Item = Self::Node> 
        {
            DFSPostOrderIterator::new(self, start_node)
        }
}

impl BFS for SimpleRootedTree{
    fn bfs(&self, start_node_id: Self::NodeID)-> impl Iterator<Item = Self::Node>
    {
        BFSIterator::new(self, start_node_id)
    }
}

impl ContractTree for SimpleRootedTree{}

impl EulerWalk for SimpleRootedTree{

    fn precompute_fai(&mut self) {
        let max_id = self.nodes.len();
        let index = self.first_appearance();
        let mut fai_vec = vec![None; max_id];
        for id in self.get_node_ids(){
            fai_vec[id] = Some(index[&id]);
        }
        self.precomputed_fai = Some(fai_vec);
    }

    fn precompute_da(&mut self) {
        let da = self.depth_array();
        let rmq = BinaryRmq::from_vec(da.iter().map(|x| x.clone() as u64).collect_vec());
        self.precomputed_rmq = Some(rmq);

        self.precomputed_da = Some(self.depth_array());
    }

    fn precompute_walk(&mut self) {
        self.precomputed_euler = Some(self.euler_walk(self.get_root_id()).into_iter().map(|x| x.get_id()).collect_vec());
    }

    fn get_precomputed_walk(&self)->Option<&Vec<<Self as RootedTree>::NodeID>> {
        self.precomputed_euler.as_ref()
    }

    fn precompute_constant_time_lca(&mut self)
    {
        self.precompute_walk();
        self.precompute_da();
        self.precompute_fai();
    }


    fn get_precomputed_fai(&self)->Option<impl Index<&Self::NodeID, Output = usize>> {
        match &self.precomputed_fai{
            Some(fai_vec) => {
                let hashmap = fai_vec.iter()
                    .enumerate()
                    .filter(|(_, idx)| idx.is_some())
                    .map(|(id, idx)| (id, idx.clone().unwrap()))
                    .collect::<HashMap<Self::NodeID, usize>>();
                return Some(hashmap);
            },
            None => None,
        }
    }

    fn get_precomputed_da(&self)->Option<&Vec<usize>> {
        self.precomputed_da.as_ref()
    }

    fn is_euler_precomputed(&self)->bool {
        self.precomputed_euler.is_some()
    }

    fn first_appearance(&self)->HashMap<Self::NodeID, usize>
    {
        let mut fa: HashMap<Self::NodeID, usize> = [].into_iter().collect::<HashMap<_,_>>();
        match self.get_precomputed_walk(){
            Some(walk) => {
                for node_id in self.get_node_ids(){
                    fa.insert(node_id, walk.iter().position(|x| x==&node_id).unwrap());
                }
                // for (idx, node_id) in walk.iter().enumerate(){
                //     if fa.contains_key(node_id){}
                //     else {
                //         fa.insert(node_id.clone(), idx);
                //     }
                // }
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

    fn get_fa_index(&self, node_id: &Self::NodeID)->usize
    {
        match &self.precomputed_fai{
            Some(fai) => fai[node_id.clone()].unwrap(),
            None => self.first_appearance()[node_id],
        }

    }

    fn get_lca_id(&self, node_id_vec: &Vec<Self::NodeID>)-> Self::NodeID
    {
        let min_pos = node_id_vec.iter().map(|x| self.get_fa_index(&x)).min().unwrap();
        let max_pos = node_id_vec.iter().map(|x| self.get_fa_index(&x)).max().unwrap();

        if min_pos==max_pos{
            return min_pos;
        }
        let dp = self.precomputed_rmq.as_ref().unwrap();
        
        return self.get_euler_pos(dp.range_min(min_pos, max_pos));

    }


}

impl Clusters for SimpleRootedTree{}

impl Newick for SimpleRootedTree
{
    type Weight = f32;
    type Meta = String;
    fn from_newick(newick_str: &[u8])->Self {
                let mut tree = SimpleRootedTree::new(0);
                let mut stack : Vec<Self::NodeID> = Vec::new();
                let mut context : Self::NodeID = tree.get_root_id();
                let mut taxa_str = String::new();
                let mut decimal_str: String = String::new();
                let mut str_ptr: usize = 0;
                let newick_string = String::from_utf8(newick_str.to_vec()).unwrap().chars().filter(|c| !c.is_whitespace()).collect::<Vec<char>>();
                while str_ptr<newick_string.len(){
                    match newick_string[str_ptr]{
                        '(' => {
                            stack.push(context);
                            let new_node = Node::new(tree.next_id());
                            context = new_node.get_id();
                            tree.set_node(new_node);
                            str_ptr +=1;
                        },
                        ')'|',' => {
                            // last context id
                            let last_context = stack.last().expect("Newick string ended abruptly!");
                            // add current context as a child to last context
                            tree.set_child(
                                *last_context,
                                context,
                            );
                            tree.set_edge_weight((*last_context,context), decimal_str.parse::<Self::Weight>().ok());
                            if !taxa_str.is_empty(){
                                tree.set_node_taxa(context, Some(taxa_str.to_string()));
                            }
                            // we clear the strings
                            taxa_str.clear();
                            decimal_str.clear();
        
                            match newick_string[str_ptr] {
                                ',' => {
                                    let new_node = Node::new(tree.next_id());
                                    context = new_node.get_id();
                                    tree.set_node(new_node);
                                    str_ptr += 1;
                                }
                                _ => {
                                    context = stack.pop().expect("Newick string ended abruptly!");
                                    str_ptr += 1;
                                }
                            }
                        },
                        ';'=>{
                            if !taxa_str.is_empty(){
                                tree.set_node_taxa(context, Some(taxa_str));
                            }
                            break;
                        }
                        ':' => {
                                // if the current context had a weight
                                if newick_string[str_ptr]==':'{
                                    str_ptr+=1;
                                    while newick_string[str_ptr].is_ascii_digit() || newick_string[str_ptr]=='.'{
                                        decimal_str.push(newick_string[str_ptr]); 
                                        str_ptr+=1;
                                    }
                                }
                            }
                        _ => {
                            // push taxa characters into taxa string
                            while newick_string[str_ptr]!=':'&&newick_string[str_ptr]!=')'&&newick_string[str_ptr]!=','&&newick_string[str_ptr]!='('&&newick_string[str_ptr]!=';'{
                                taxa_str.push(newick_string[str_ptr]); 
                                str_ptr+=1;
                            }
                        },
                    }
                }
                tree
            }

    fn subtree_to_newick(&self, node_id: Self::NodeID)-> impl std::fmt::Display {
        let node  = self.get_node(node_id).unwrap();
        let mut tmp = String::new();
        if node.get_children().into_iter().len()!=0{
            if node.get_children().into_iter().len()>1{
                tmp.push('(');
            }
            for child_id in node.get_children().into_iter(){
                let child_str = format!("{},", self.subtree_to_newick(child_id));
                tmp.push_str(&child_str);
            }
            tmp.pop();
            if node.get_children().into_iter().collect_vec().len()>1{
                tmp.push(')');
            }
        }
        tmp.push_str(&node.get_taxa().unwrap_or("".to_string()));
        tmp
    }
}

impl SPR for SimpleRootedTree
{
    fn graft(&mut self, tree: Self, edge: (Self::NodeID, Self::NodeID)) {
        let new_node = self.next_node();
        let new_node_id = dbg!(new_node.get_id());
        for node in tree.dfs(tree.get_root_id())
        {
            self.set_node(dbg!(node));
        }
        self.split_edge(edge, new_node);
        self.set_child(dbg!(new_node_id), dbg!(tree.get_root_id()));

    }
    fn prune(&mut self, node_id: Self::NodeID)-> Self {
        let mut pruned_tree = SimpleRootedTree::new(node_id);
        let p_id = self.get_node_parent_id(node_id).unwrap();
        self.get_node_mut(p_id).unwrap().remove_child(&node_id);
        pruned_tree.get_node_mut(pruned_tree.get_root_id()).unwrap().add_children(self.get_node(node_id).unwrap().get_children());
        let dfs = self.dfs(node_id).into_iter().collect_vec();
        for node in dfs
        {    
            self.nodes.remove(node.get_id());
            pruned_tree.set_node(node);
        }
        pruned_tree
    }
}

impl Balance for SimpleRootedTree{
    fn balance_subtree(&mut self) 
    {
        assert!(self.get_cluster(self.get_root_id()).into_iter().collect_vec().len()==4, "Quartets have 4 leaves!");
        assert!(self.is_binary(), "Cannot balance non-binary tree!");
        let root_children = self.get_node_children(self.get_root_id()).into_iter().collect_vec();
        let (child1, child2) = (root_children[0].get_id(), root_children[1].get_id());
        let next_id = self.next_id();
        let split_id = self.next_id()+1;
        match dbg!(((self.get_node(child1).unwrap().is_leaf()), (self.get_node(child2).unwrap().is_leaf()))){
            (false, false) => {},
            (true, false) => {
                let mut leaf_node = self.remove_node(child1).unwrap();
                leaf_node.set_id(next_id);
                let other_leaf_id = &self.get_node_children(child2).into_iter().filter(|node| node.is_leaf()).collect_vec()[0].get_id();
                self.split_edge((child2, *other_leaf_id), Node::new(split_id));
                self.add_child(dbg!(split_id), leaf_node);
            },
            (false, true) => {
                let mut leaf_node = self.remove_node(child2).unwrap();
                leaf_node.set_id(next_id);
                let other_leaf_id = &self.get_node_children(child1).into_iter().filter(|node| node.is_leaf()).collect_vec()[0].get_id();
                self.split_edge((child1, *other_leaf_id), Node::new(split_id));
                self.add_child(split_id, leaf_node);
            },
            _ =>{}
        }
        self.clean();
    }
}

impl CopheneticDistance for SimpleRootedTree{
    type Meta = String;

}