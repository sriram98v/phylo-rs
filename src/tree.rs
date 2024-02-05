pub mod simple_rtree;
pub mod ops;
pub mod distances;
pub mod weighted;
pub mod io;


use std::collections::HashMap;
use std::rc::Rc;
use itertools::Itertools;

use crate::node::simple_rnode::RootedTreeNode;
use crate::node::{Node, NodeID};
use crate::tree::simple_rtree::RootedTree;
use crate::iter::node_iter::*;

use self::io::Newick;
use self::ops::{Balance, Subtree, SPR};
// use crate::tree::ops::SPR;
// use crate::tree::distances::*;
// use crate::iter::{node_iter::*, edge_iter::*};

// pub struct UnrootedPhyloTree{
//     _nodes: HashMap<NodeID, NodeType>,
//     _neighbours: HashMap<NodeID, HashSet<(Option<EdgeWeight>, NodeID)>>,
//     _leaves: HashMap<NodeID, String>,
// }

#[derive(Debug)]
pub struct SimpleRootedTree{
    root: NodeID,
    nodes: HashMap<NodeID, Node>,
}

impl SimpleRootedTree{
    pub fn next_id(&self)->NodeID
    {
        match self.nodes.keys().map(|x| (*Rc::clone(x)).clone()).max()
        {
            Some(x) => Rc::new(x+1),
            None => Rc::new(0)
        }
    }

    pub fn next_node(&self)->Node
    {
        Node::new(self.next_id(), false)
    }
    // pub fn new()->Self{
    //     let root_node = Node::new(Rc::new(0), false);
    //     SimpleRootedTree { 
    //         root: root_node.get_id(),
    //         nodes: HashMap::from([(root_node.get_id(), root_node)]),
    //     }
    // }
}

// impl Default for SimpleRootedTree {
//     fn default() -> Self {
//         Self::new()
//     }
// }

impl RootedTree for SimpleRootedTree{
    
    type NodeID = NodeID;
    type Taxa = String;
    type Node = Node;

    fn new(root_id: Self::NodeID)->Self{
        let root_node = Node::new(root_id, false);
        SimpleRootedTree { 
            root: root_node.get_id(),
            nodes: HashMap::from([(root_node.get_id(), root_node)]),
        }
    }

    /// Returns reference to node by ID
    fn get_node(&self, node_id: Self::NodeID)->Option<&Self::Node>
    {
        self.nodes.get(&node_id)
    }

    fn get_node_mut(&mut self, node_id: Self::NodeID)->Option<&mut Self::Node>
    {
        self.nodes.get_mut(&node_id)
    }

    fn get_node_ids(&self)->impl IntoIterator<Item = Self::NodeID, IntoIter = impl ExactSizeIterator<Item = Self::NodeID>> 
    {
        self.nodes.clone().into_keys()    
    }

    fn get_nodes(&self)->impl IntoIterator<Item = Self::Node, IntoIter = impl ExactSizeIterator<Item = Self::Node>> 
    {
        self.nodes.clone().into_values()    
    }

    /// Returns reference to node by ID
    fn set_node(&mut self, node: Self::Node)
    {
        self.nodes.insert(node.get_id(), node);
    }

    fn add_child(&mut self, parent_id: Self::NodeID, child: Self::Node)
    {
        let new_child_id = child.get_id();
        self.set_node(child);
        self.get_node_mut(Rc::clone(&parent_id)).unwrap().add_child(Rc::clone(&new_child_id));
        self.get_node_mut(new_child_id).unwrap().set_parent(Some(parent_id));
    }

    /// Get root node ID
    fn get_root_id(&self)->Self::NodeID
    {
        Rc::clone(&self.root)
    }

    fn set_root(&mut self, node_id: Self::NodeID) 
    {
        self.root = Rc::clone(&node_id);
    }

    fn remove_node(&mut self, node_id: Self::NodeID)->Option<Self::Node>
    {
        match self.get_node_parent(node_id.clone())
        {
            Some(pid) => self.get_node_mut(pid).unwrap().remove_child(node_id.clone()),
            None => {},
        }
        self.nodes.remove(&node_id)
    }

    fn contains_node(&self, node_id: Self::NodeID)->bool
    {
        self.nodes.contains_key(&node_id)
    }

    fn delete_edge(&mut self, parent_id: Self::NodeID, child_id: Self::NodeID)
    {
        self.get_node_mut(parent_id).unwrap().remove_child(Rc::clone(&child_id));
        self.get_node_mut(child_id).unwrap().set_parent(None);
    }

    fn clean(&mut self)
    {
        // todo!()
        let node_iter = self.get_nodes().into_iter().collect::<Vec<Self::Node>>();
        for node in node_iter.clone(){
            // remove root with only one child
            let node_id = node.get_id();
            if node.get_id()==self.get_root_id() && node.degree()<2{
                let new_root = self.get_root().get_children().into_iter().next().unwrap();
                self.set_root(new_root);
                self.get_node_mut(Rc::clone(&self.root)).unwrap().set_parent(None);
                self.remove_node(node_id.clone());
            }
            // remove nodes with only one child
            else if !node.is_leaf() && node.get_parent()!=None && node.degree()<3 {
                let parent_id = self.get_node_parent(node_id.clone());
                let child_id = node.get_children().into_iter().next().unwrap();
                self.get_node_mut(Rc::clone(&child_id)).unwrap().set_parent(parent_id.clone());
                self.get_node_mut(parent_id.unwrap()).unwrap().add_child(child_id);
                self.remove_node(node.get_id());
            }
            for chid in node.get_children().into_iter()
            {
                if !node_iter.clone().into_iter().map(|x| x.get_id()).contains(&chid)
                {
                    self.get_node_mut(node_id.clone()).unwrap().remove_child(chid);
                }
            }
        }
    }

    fn get_mrca(&self, node_id_list: &Vec<Self::NodeID>)->Self::NodeID
    {
        let euler_tour = self.euler_tour(self.get_root_id()).into_iter().map(|x| x.get_id()).collect_vec();
        let depth_array: Vec<usize> = euler_tour.iter().map(|x| self.get_node_depth(x.clone())).collect_vec();   // todo
        let mut min_pos = euler_tour.len();
        let mut max_pos = 0;
        for node_id in node_id_list
        {
            let pos = euler_tour.iter().position(|r| r == node_id).unwrap();
            match pos<min_pos {
                true => min_pos=pos,
                false => {},
            }
            let pos = euler_tour.iter().rposition(|r| r == node_id).unwrap_or(0);
            match pos>max_pos {
                true => max_pos=pos,
                false => {},
            }
        }
        let depth_subarray_min_value = depth_array[min_pos..max_pos].iter().min().unwrap();
        let depth_subarray_min_pos = depth_array[min_pos..max_pos].iter().position(|x| x==depth_subarray_min_value).unwrap();
        Rc::clone(&euler_tour[min_pos..max_pos][depth_subarray_min_pos])
        // todo!()
    }
}

impl PreOrder for SimpleRootedTree{}

impl PostOrder for SimpleRootedTree{}

impl DFS for SimpleRootedTree{}

impl EulerTour for SimpleRootedTree{}

impl Newick for SimpleRootedTree{
    fn from_newick(newick_str: &[u8])->Self {
                let mut tree = SimpleRootedTree::new(Rc::new(0));
                let mut stack : Vec<NodeID> = Vec::new();
                let mut context : NodeID = tree.get_root_id();
                let mut taxa_str = String::new();
                let mut str_ptr: usize = 0;
                let newick_string = String::from_utf8(newick_str.to_vec()).unwrap().chars().filter(|c| !c.is_whitespace()).collect::<Vec<char>>();
                while str_ptr<newick_string.len(){
                    match newick_string[str_ptr]{
                        '(' => {
                            stack.push(context);
                            let new_node = Node::new(tree.next_id(), false);
                            context = new_node.get_id();
                            tree.set_node(new_node);
                            str_ptr +=1;
                        },
                        ')'|',' => {
                            // last context id
                            let last_context = stack.last().expect("Newick string ended abruptly!");
                            // add current context as a child to last context
                            tree.set_child(
                                Rc::clone(last_context),
                                Rc::clone(&context),
                            );
                            if !taxa_str.is_empty(){
                                tree.set_node_taxa(context, Some(taxa_str.to_string()));
                            }
                            // we clear the strings
                            taxa_str.clear();
        
                            match newick_string[str_ptr] {
                                ',' => {
                                    let new_node = Node::new(tree.next_id(), false);
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
                        _ => {
                            // push taxa characters into taxa string
                            while newick_string[str_ptr]!=')'&&newick_string[str_ptr]!=','&&newick_string[str_ptr]!='('&&newick_string[str_ptr]!=';'{
                                taxa_str.push(newick_string[str_ptr]); 
                                str_ptr+=1;
                            }
                        },
                    }
                }
                let leaf_ids = tree.dfs(tree.get_root_id())
                    .into_iter()
                    .filter(|x| x.get_children().into_iter().collect_vec().is_empty())
                    .map(|x| Rc::clone(&x.get_id()))
                    .collect_vec();
                tree.flip_nodes(leaf_ids.clone().into_iter());
                tree
            }

    fn subtree_to_newick(&self, node_id: Self::NodeID)-> impl std::fmt::Display {
        let node  = self.get_node(dbg!(node_id)).unwrap();
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

        // String::new()
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
        let mut pruned_tree = SimpleRootedTree::new(Rc::clone(&node_id));
        let p_id = self.get_node_parent(Rc::clone(&node_id)).unwrap();
        self.get_node_mut(p_id).unwrap().remove_child(Rc::clone(&node_id));
        pruned_tree.get_node_mut(pruned_tree.get_root_id()).unwrap().add_children(self.get_node(Rc::clone(&node_id)).unwrap().get_children());
        let dfs = self.dfs(node_id.clone()).into_iter().collect_vec().clone();
        for node in dfs
        {    
            self.nodes.remove(&node.get_id());
            pruned_tree.set_node(node);
        }
        pruned_tree
    }
}

impl Clusters for SimpleRootedTree{}

impl Balance for SimpleRootedTree{
    fn balance_subtree(&mut self) 
    {
        assert!(self.get_cluster(self.get_root_id()).into_iter().collect_vec().len()==4, "Quartets have 4 leaves!");
        assert!(self.is_binary(), "Cannot balance non-binary tree!");
        let root_children = self.get_node_children(self.get_root_id()).into_iter().collect_vec();
        let (child1, child2) = (root_children[0].get_id(), root_children[1].get_id());
        let next_id = self.next_id();
        let split_id = Rc::new(*self.next_id()+1);
        match dbg!(((self.get_node(child1.clone()).unwrap().is_leaf()), (self.get_node(child2.clone()).unwrap().is_leaf()))){
            (false, false) => {},
            (true, false) => {
                let mut leaf_node = self.remove_node(child1).unwrap();
                leaf_node.set_id(next_id.clone());
                let other_leaf_id = &self.get_node_children(child2.clone()).into_iter().filter(|node| node.is_leaf()).collect_vec()[0].get_id();
                self.split_edge((child2, Rc::clone(other_leaf_id)), Node::new(split_id.clone(), false));
                self.add_child(dbg!(split_id), leaf_node);
            },
            (false, true) => {
                let mut leaf_node = self.remove_node(child2).unwrap();
                leaf_node.set_id(next_id.clone());
                let other_leaf_id = &self.get_node_children(child1.clone()).into_iter().filter(|node| node.is_leaf()).collect_vec()[0].get_id();
                self.split_edge((child1, Rc::clone(other_leaf_id)), Node::new(split_id.clone(), false));
                self.add_child(split_id, leaf_node);


                // let mut p_t = self.prune(child2);
                // dbg!(self.get_nodes().into_iter().collect_vec());
                // p_t.get_node_mut(p_t.get_root_id()).unwrap().set_id(next_id.clone());
                // p_t.set_root(next_id);
                // let other_leaf = &self.get_node_children(child1.clone()).into_iter().filter(|node| node.is_leaf()).collect_vec()[0].get_id();
                // self.graft(p_t, (child1, Rc::clone(other_leaf)));
                // dbg!(self.get_nodes().into_iter().collect_vec());
            },
            _ =>{}
        }
        self.clean();
    }
}

impl Ancestors for SimpleRootedTree{}

impl Subtree for SimpleRootedTree {}

//     pub fn from_newick(newick_string: String)->Self{
//         let mut tree = RootedPhyloTree::new();
//         let mut stack : Vec<NodeID> = Vec::new();
//         let mut context : NodeID = *tree.get_root();
//         let mut taxa_str = String::new();
//         let mut decimal_str: String = String::new();
//         let mut str_ptr: usize = 0;
//         let newick_string = newick_string.chars().filter(|c| !c.is_whitespace()).collect::<Vec<char>>();
//         while str_ptr<newick_string.len(){
//             match newick_string[str_ptr]{
//                 '(' => {
//                     stack.push(context);
//                     context = tree.add_node();
//                     str_ptr +=1;
//                 },
//                 ')'|',' => {
//                     // last context id
//                     let last_context = stack.last().expect("Newick string ended abruptly!");
//                     // add current context as a child to last context
//                     tree.set_child(
//                         &context,
//                         last_context,
//                         decimal_str.parse::<EdgeWeight>().ok(),
//                         match taxa_str.is_empty(){
//                             true => None,
//                             false => Some(taxa_str.to_string())
//                         }
//                     );
//                     // we clear the strings
//                     taxa_str.clear();
//                     decimal_str.clear();

//                     match newick_string[str_ptr] {
//                         ',' => {
//                             context = tree.add_node();
//                             str_ptr += 1;
//                         }
//                         _ => {
//                             context = stack.pop().expect("Newick string ended abruptly!");
//                             str_ptr += 1;
//                         }
//                     }
//                 },
//                 ';'=>{
//                     if !taxa_str.is_empty(){
//                         tree.assign_taxa(&context, &taxa_str);
//                     }
//                     break;
//                 }
//                 ':' => {
//                     // if the current context had a weight
//                     if newick_string[str_ptr]==':'{
//                         str_ptr+=1;
//                         while newick_string[str_ptr].is_ascii_digit() || newick_string[str_ptr]=='.'{
//                             decimal_str.push(newick_string[str_ptr]); 
//                             str_ptr+=1;
//                         }
//                     }
//                 }
//                 _ => {
//                     // push taxa characters into taxa string
//                     while newick_string[str_ptr]!=':'&&newick_string[str_ptr]!=')'&&newick_string[str_ptr]!=','&&newick_string[str_ptr]!='('&&newick_string[str_ptr]!=';'{
//                         taxa_str.push(newick_string[str_ptr]); 
//                         str_ptr+=1;
//                     }
//                 },
//             }
//         }
//         let mut leaf_ids = Vec::new();
//         tree.leaves_of_node(tree.get_root(), &mut leaf_ids);
//         for leaf_id in leaf_ids{
//             tree.set_leaf(&leaf_id);
//         }
//         tree
//     }

//     fn leaves_of_node(&self, node_id:&NodeID, leaves:&mut Vec<NodeID>){
//         if self.get_node_children(node_id).is_empty(){
//             leaves.push(*node_id);
//         }

//         for (child_node_id, _edge_weight) in self.get_node_children(node_id).iter(){
//             self.leaves_of_node(child_node_id, leaves);
//         }
//     }
// }

// impl SimpleRTree for RootedPhyloTree{
//     fn add_node(&mut self)->NodeID{
//         // New node id
//         let node_id = self.nodes.keys().max().unwrap_or(self.get_root())+&1;
//         // add entry of node in parents and children fields
//         self.nodes.insert(node_id, NodeType::Internal(None));
//         self.parents.insert(node_id, None);
//         self.children.insert(node_id, Vec::new());
//         node_id
//     }

//     fn set_child(&mut self, node_id:&NodeID, parent_id:&NodeID, distance:Option<EdgeWeight>, taxa:Option<String>){
//         self.parents.insert(*node_id, Some(*parent_id));
//         self.children.entry(*parent_id).or_default().push((*node_id, distance));
//         self.nodes.insert(*node_id, NodeType::Internal(taxa));
//     }

//     fn set_leaf(&mut self, node_id: &NodeID) {
//         self.nodes.entry(*node_id).and_modify(|node| node.flip());
//     }

//     fn assign_taxa(&mut self,node:&NodeID, taxa:&str) {
//         self.nodes.insert(*node, NodeType::Internal(Some(taxa.to_string())));
//     }

//     fn set_edge_weight(&mut self, parent:&NodeID, child:&NodeID, edge_weight:Option<EdgeWeight>){
//         self.children.entry(*parent)
//             .and_modify(|children| *children = children.clone().iter()
//                     .map(|(id, w)| {
//                         match id==child{
//                             true => {(*id, edge_weight)},
//                             false => {(*id, *w)},
//                         }
//                     })
//                     .collect()
//         );
//     }

//     fn get_root(&self)->&NodeID{
//         &self.root
//     }

//     fn get_node(&self, node_id: &NodeID)->&NodeType{
//         self.nodes.get(node_id).expect("Invalid NodeID")
//     }

//     fn get_nodes(&self)->&HashMap<NodeID, NodeType>{
//         &self.nodes
//     }

//     fn get_children(&self)->&HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>>{
//         &self.children
//     }

//     fn get_parents(&self)->&HashMap<NodeID, Option<NodeID>>{
//         &self.parents
//     }


//     fn get_node_children(&self, node_id: &NodeID)->&Vec<(NodeID, Option<EdgeWeight>)>{
//         self.children.get(node_id).expect("Invalid NodeID!")
//     }

//     fn get_node_parent(&self, node_id:&NodeID)->Option<&NodeID>{
//         self.parents.get(node_id).expect("Invalid NodeID!").as_ref()
//     }

//     fn get_leaves(&self, node_id: &NodeID)->Vec<(NodeID, NodeType)>{
//         let mut leaf_vec: Vec<NodeID> = Vec::new();
//         self.leaves_of_node(node_id, &mut leaf_vec);
//         leaf_vec.into_iter().map(|leaf_id| (leaf_id, self.nodes.get(&leaf_id).cloned().expect("Invalid NodeID!"))).collect::<Vec<(NodeID, NodeType)>>()
//     }

//     fn get_subtree(&self, node_id: &NodeID)->Box<dyn SimpleRTree>{
//         if self.is_leaf(node_id){
//             panic!("NodeID is a leaf");
//         }
//         let root= *node_id;
//         let mut nodes: HashMap<NodeID, NodeType>= HashMap::new();
//         let mut children: HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>> = HashMap::new();
//         let mut parents: HashMap<NodeID, Option<NodeID>> = HashMap::new();
//         for decsendant_node_id in self.iter_node_pre(node_id){
//             nodes.insert(decsendant_node_id, self.nodes.get(&decsendant_node_id).expect("Invalid NodeID!").clone());
//             children.insert(decsendant_node_id, self.children.get(&decsendant_node_id).expect("Invalid NodeID!").clone());
//             parents.insert(decsendant_node_id, *self.parents.get(&decsendant_node_id).expect("Invalid NodeID!"));
//         }
//         Box::new(
//             RootedPhyloTree{
//                 root,
//                 nodes,
//                 children,
//                 parents,
//             }
//         )
//     }

//     fn get_mrca(&self, node_id_list: Vec<&NodeID>)->NodeID{
//         let ancestor_iter_vec: Vec<std::vec::IntoIter<NodeID>> = node_id_list.iter().map(|x| self.get_ancestors_pre(x).into_iter()).collect();
//         let mut mrca: NodeID = 0;
//         for mut iterator in ancestor_iter_vec{
//             let temp: HashSet<NodeID> = HashSet::new();
//             if let Some(x) = iterator.next() {
//                 match temp.contains(&x){
//                     true => {mrca = x},
//                     false => {
//                         match temp.is_empty(){
//                             true => {},
//                             false => {return mrca}
//                         }
//                     }
//                 }
//             }
//         }
//         mrca
//     }

//     fn is_leaf(&self, node_id: &NodeID)->bool{
//         self.nodes.get(node_id).expect("Invalid NodeID").is_leaf()
//     }

//     fn graft(&mut self, tree: Box<dyn SimpleRTree>, edge: (&NodeID, &NodeID), edge_weights:(Option<EdgeWeight>, Option<EdgeWeight>), graft_edge_weight: Option<EdgeWeight>){
//         let graft_node = self.split_edge(edge, edge_weights);
//         let input_root_id = tree.get_root();
//         for input_node in tree.get_nodes().keys(){
//             if self.get_nodes().contains_key(input_node){
//                 panic!("The NodeIDs in the input tree are already present in the current tree!");
//             }
//         }

//         self.children.extend(tree.get_children().clone().into_iter());
//         self.parents.extend(tree.get_parents().clone().iter());
//         self.nodes.extend(tree.get_nodes().clone().into_iter());
//         self.set_child(input_root_id, &graft_node, graft_edge_weight, Some(tree.get_taxa(input_root_id)))
//     }

//     fn prune(&mut self, node_id: &NodeID)-> Box<dyn SimpleRTree>{
//         let root= *node_id;
//         let root_parent = self.get_node_parent(node_id).expect("Node has no parent! Clean tree first...");
//         self.children.entry(*root_parent).or_default().retain(|(child_id, _w)| *child_id!=root);
//         self.parents.insert(root, None);

//         let mut nodes: HashMap<NodeID, NodeType>= HashMap::new();
//         let mut children: HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>> = HashMap::new();
//         let mut parents: HashMap<NodeID, Option<NodeID>> = HashMap::new();
        
//         for decsendant_node_id in self.iter_node_pre(node_id){
//             nodes.insert(decsendant_node_id, self.nodes.remove(&decsendant_node_id).expect("Invalid NodeID!").clone());
//             children.insert(decsendant_node_id, self.children.remove(&decsendant_node_id).expect("Invalid NodeID!").clone());
//             parents.insert(decsendant_node_id, self.parents.remove(&decsendant_node_id).expect("Invalid NodeID!"));
//             }
//         Box::new(
//             RootedPhyloTree{
//                 root,
//                 nodes,
//                 children,
//                 parents,
//             }
//         )
//     }

//     fn iter_node_pre(&self, start_node_id: &NodeID)->PreOrdNodes{
//         PreOrdNodes::new(start_node_id, &self.children)
//     }

//     fn iter_node_post(&self, start_node_id: &NodeID)->PostOrdNodes{
//         PostOrdNodes::new(start_node_id, &self.children)
//     }

//     fn iter_edges_pre(&self, start_node_id: &NodeID)->PreOrdEdges{
//         PreOrdEdges::new(self, start_node_id)
//     }

//     fn iter_edges_post(&self, start_node_id: &NodeID)->PostOrdEdges{
//         PostOrdEdges::new(self, start_node_id)
//     }

//     fn get_ancestors_pre(&self, node_id: &NodeID)->Vec<NodeID>{
//         let mut node_iter: Vec<NodeID> = Vec::new();
//         let mut curr_node = node_id;
//         while self.parents.get(curr_node).is_some() {
//             match self.parents.get(curr_node).expect("Invalid NodeID!") {
//                 Some(node) => {
//                     node_iter.push(*node);
//                     curr_node = node;
//                 },
//                 None => {
//                     node_iter.push(*self.get_root());
//                     break;
//                 },
//             }
//         }
//         node_iter
//     }

//     fn reroot_at_node(&mut self, node_id: &NodeID){
//         let mut stack: Vec<NodeID> = vec![node_id.clone()];
//         let mut neighbours: HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>> = self.children.clone();
//         let parent_as_edge = self.parents.clone().into_iter()
//         .filter(|(_child_id, parent_id)| parent_id!=&None)
//         .map(|(child_id, parent_id)| (child_id, vec![(parent_id.unwrap(), self.get_edge_weight(parent_id.as_ref().unwrap(), &child_id).cloned())]));
//         for (id, edges) in parent_as_edge{
//             neighbours.entry(id).or_default().extend(edges);
//         }
//         let mut new_children: HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>> = HashMap::new();
//         let mut new_parents: HashMap<NodeID, Option<NodeID>> = HashMap::from([(node_id.clone(), None)]);

//         while !stack.is_empty(){
//             let curr_node = stack.pop().unwrap();
//             if let Some(child) = neighbours.remove(&curr_node){
//                 let curr_node_children = &child.iter().filter(|(id, _w)| !new_parents.keys().contains(id));
//                 new_children.entry(curr_node).or_default().extend(curr_node_children.clone());
//                 for (id, _w) in &child{
//                     new_parents.insert(id.clone(), Some(curr_node.clone()));
//                 }
//                 stack.extend(child.iter().map(|(id, _w)| id.clone()))
//             }
//         }

//         self.children = dbg!(new_children);
//         self.parents = dbg!(new_parents);
//         self.root = *dbg!(node_id);
//     }

//     fn split_edge(&mut self, edge: (&NodeID, &NodeID), edge_weights:(Option<EdgeWeight>, Option<EdgeWeight>))->NodeID{
//         let new_node_id = self.add_node();
//         self.parents.insert(new_node_id, Some(edge.0.clone()));
//         self.children.entry(new_node_id).or_default().push((edge.1.clone(), edge_weights.1));
//         self.parents.insert(edge.1.clone(), Some(new_node_id));
//         self.children.entry(edge.0.clone()).or_default().retain(|(id, _w)| id!=edge.1);
//         self.children.entry(edge.0.clone()).or_default().push((new_node_id, edge_weights.0));
//         new_node_id
//     }

//     fn distance_from_ancestor(&self, node: &NodeID, ancestor: &NodeID, weighted: bool)->f64{
//         let binding = self.get_ancestors_pre(node);
//         let start_idx = binding.iter().position(|&x| x==*ancestor).expect("Provided ancestor is not an ancestor of node!");
//         let mut node_ancestor_pre = binding[start_idx..].iter();
//         let mut curr_parent = node_ancestor_pre.next().unwrap();
//         let mut distance  = 0 as f64;
//         while let Some(node_id) = node_ancestor_pre.next() {
//             let curr_parent_children = self.get_node_children(curr_parent);
//             for (child_id, w) in curr_parent_children{
//                 if child_id==node_id{
//                     match weighted {
//                         true => {distance += w.unwrap_or(0 as f64);}
//                         false => {distance += 1_f64;}
//                     }
//                     curr_parent = node_id;
//                     continue;
//                 }
//                 panic!("Ancestor chain is broken! Clean tree before moving forward...")
//             } 
//         };
//         distance
//     }

//     fn get_bipartition(&self, edge: (&NodeID, &NodeID))->(Vec<(NodeID, NodeType)>, Vec<(NodeID, NodeType)>){
//         let c2 = self.get_cluster(edge.1);
//         (self.nodes.clone().into_iter().filter(|x| !c2.contains(x)).collect_vec(), c2)
//     }

//     fn get_cluster(&self, node_id: &NodeID)-> Vec<(NodeID, NodeType)>{
//         let mut leaves: Vec<NodeID> = Vec::new();
//         self.leaves_of_node(node_id, &mut leaves);
//         leaves.into_iter().map(|leaf_id| (leaf_id, self.get_node(&leaf_id).clone())).collect_vec()
//     }

//     fn clean(&mut self) {
//         let mut remove_list: Vec<&NodeID> = Vec::new();
//         for (node_id, node) in self.nodes.clone().iter(){
//             // remove root with only one child
//             if node_id==self.get_root() && self.get_node_degree(node_id)<2{
//                 let new_root = self.get_node_children(self.get_root())[0].0;
//                 self.root = new_root;
//                 self.parents.entry(new_root).and_modify(|x| *x = None);
//                 remove_list.push(node_id);
//             }
//             // remove nodes with only one child
//             else if !node.is_leaf() &&  self.get_node_degree(node_id)<3{
//                 let parent = self.get_node_parent(node_id).cloned();
//                 let children = self.get_node_children(node_id).clone();
//                 for (child_id, _edge_weight) in children.clone().into_iter(){
//                     self.parents.entry(child_id).and_modify(|x| *x = parent);
//                 }
//                 self.set_children(parent.as_ref().unwrap(), &children);
//             }
//         }
//     }

//     fn get_taxa(&self, node_id:&NodeID)->String {
//         self.get_node(node_id).taxa()
//     }

//     fn incerement_ids(&mut self, value: &usize){
//         self.nodes = self.nodes.clone().into_iter().map(|(node_id, node_type)| (node_id+value, node_type)).collect();
//         self.parents = self.parents.clone().into_iter().map(|(node_id, parent_id)| {
//             (
//                 node_id+value, 
//                 parent_id.map(|id| id + value)
//             )
//         }).collect();
//         self.children = self.children.clone().into_iter().map(|(node_id, children_vec)| {
//             (
//                 node_id+value,
//                 children_vec.into_iter().map(|(child_id, w)| {
//                     (
//                         child_id+value,
//                         w
//                     )
//                 })
//                 .collect()
//             )
//         }).collect();
//     }

// }

// impl RPhyTree for RootedPhyloTree{
//     fn induce_tree(&self, taxa: Vec<String>)->Box<dyn RPhyTree>{
//         let mut nodes: HashMap<NodeID, NodeType> = HashMap::new();
//         let leaf_ids = self.get_nodes()
//             .iter()
//             .filter(|(_id, n_type)| taxa.contains(&n_type.taxa()))
//             .map(|(id, _)| (id));
//         for id in leaf_ids{
//             nodes.insert(id.clone(), self.get_node(id).clone());
//             nodes.extend(self.get_ancestors_pre(id).iter().map(|node_id| (node_id.clone(), self.get_node(node_id).clone())).collect::<HashMap<NodeID, NodeType>>());
//         }
//         let root = self.get_mrca(nodes.keys().collect_vec());
//         let children: HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>> = nodes.keys()
//             .map(|id| (id.clone(), self.get_node_children(id).into_iter().filter(|(child_id, _)| nodes.contains_key(child_id)).map(|i| i.clone()).collect_vec()))
//             .collect();
//         let mut parents: HashMap<NodeID, Option<NodeID>> = nodes.keys()
//             .map(|id| (id.clone(), self.get_node_parent(id).cloned()))
//             .collect();
//         parents.insert(root.clone(), None);
//         Box::new(
//             RootedPhyloTree{
//                 root,
//                 nodes,
//                 children,
//                 parents,
//             }
//         )
//     }

// }