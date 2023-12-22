use std::collections::HashMap;
use itertools::Itertools;

use crate::node::*;
use crate::iter::node_iter::*;
use crate::iter::edge_iter::*;

pub type EdgeWeight = f64;

pub trait SimpleRTree {
    /// Add node to tree
    fn add_node(&mut self)->NodeID;

    /// Sets node_id as child to parent.
    fn set_child(&mut self, node_id:&NodeID, parent_id:&NodeID, distance:Option<EdgeWeight>, taxa: Option<String>);

    /// Sets iterable of node_ids as children to parent
    fn set_children(&mut self, parent: &NodeID, children: &[(NodeID, Option<EdgeWeight>)]){
        for (child_id, edge_weight) in children.iter(){
            self.set_child(child_id, parent, *edge_weight, None);
        }
    }

    /// Converts internal node to leaf_node
    fn set_leaf(&mut self, node_id: &NodeID);

    /// Sets the edge weight between two nodes (None to unweight the edge)
    fn set_edge_weight(&mut self, parent:&NodeID, child:&NodeID, edge_weights:Option<EdgeWeight>);

    /// Returns true of node is child of parent.
    fn node_is_child_of(&self, parent:&NodeID, node_id:&NodeID)->bool{
        self.get_node_children(parent).iter().map(|(id, _weight)| id).contains(node_id)
    }

    /// Assign taxa to leaf node
    fn assign_taxa(&mut self,node_id:&NodeID, taxa:&str);
    
    /// Returns root node id
    fn get_root(&self)->&NodeID;
    
    /// Returns all node ids
    fn get_nodes(&self)->&HashMap<NodeID, NodeType>;

    /// Returns node by ID
    fn get_node(&self, node_id: &NodeID)->&NodeType;

    /// Returns node degree
    fn get_node_degree(&self, node_id:&NodeID)->usize{
        self.get_node_children(node_id).len() + match self.get_node_parent(node_id) {
            Some(_) => 1,
            None => 0
        }
    }

    /// Check if tree is weighted
    fn is_weighted(&self)->bool{
        for (_, _, edge_weight) in self.iter_edges_post(self.get_root()){
            if edge_weight.is_some(){
                return true;
            }
        }
        false
    }

    /// Remove all weights
    fn unweight(&mut self){
        let edge_iter = self.iter_edges_pre(self.get_root());
        for edge in edge_iter{
            self.set_edge_weight(&edge.0, &edge.1, None);
        }
    }

    /// Get all node-child relationships
    fn get_children(&self)->&HashMap<NodeID, Vec<(NodeID, Option<EdgeWeight>)>>;

    /// Get all node-parent relationships
    fn get_parents(&self)->&HashMap<NodeID, Option<NodeID>>;
    
    /// Returns children node ids for given node id 
    fn get_node_children(&self, node_id: &NodeID)->&Vec<(NodeID, Option<EdgeWeight>)>;

    /// Returns node parent
    fn get_node_parent(&self, node_id:&NodeID)->Option<&NodeID>;
    
    /// Returns all leaf node ids
    fn get_leaves(&self, node_id: &NodeID)->Vec<(NodeID, NodeType)>;
    
    /// Returns full subtree rooted at given node
    fn get_subtree(&self, node_id: &NodeID)->Box<dyn SimpleRTree>;
    
    /// Returns most recent common ancestor of give node set
    fn get_mrca(&self, node_id_list: Vec<&NodeID>)->NodeID;
    
    /// Checks if the given node is a leaf node
    fn is_leaf(&self, node_id: &NodeID)->bool;
    
    /// Attaches input tree to self by spliting an edge
    fn graft(&mut self, tree: Box<dyn SimpleRTree>, edge: (&NodeID, &NodeID), edge_weights:(Option<EdgeWeight>, Option<EdgeWeight>), graft_edge_weight: Option<EdgeWeight>);
    
    /// Returns subtree starting at given node, while corresponding nodes from self.
    fn prune(&mut self, node_id: &NodeID)-> Box<dyn SimpleRTree>;

    ///Returns an iterator that iterates over the nodes in Pre-order
    fn iter_node_pre(&self, start_node_id: &NodeID)->PreOrdNodes;
    
    ///Returns an iterator that iterates over the nodes in Post-order
    fn iter_node_post(&self, start_node_id: &NodeID)->PostOrdNodes;
    
    ///Returns an iterator that iterates over the edges in Pre-order
    fn iter_edges_pre(&self, start_node_id: &NodeID)->PreOrdEdges;
    
    ///Returns an iterator that iterates over the edges in Post-order
    fn iter_edges_post(&self, start_node_id: &NodeID)->PostOrdEdges;

    /// Returns all node ids in path from root to given node
    fn get_ancestors_pre(&self, node_id: &NodeID)->Vec<NodeID>;

    /// Returns pairwise distance matrix of the taxa. If weighted is true, then returns sum of edge weights along paths connecting leaves of tree
    fn leaf_distance_matrix(&self, weighted: bool)->HashMap<(NodeID, NodeID), EdgeWeight>{
        let binding = self.get_leaves(self.get_root());
        let leaves = binding.iter().map(|(leaf_id, _taxa)| leaf_id).combinations(2);
        let mut dist_mat: HashMap<(NodeID, NodeID), EdgeWeight> = HashMap::new();
        for node_pair in leaves{
            let w  = self.distance_from_node(node_pair[0], node_pair[1], weighted);
            dist_mat.insert((*node_pair[0], *node_pair[1]), w);
        }
        dist_mat
    }

    /// Returns pairwise distance matrix of all nodes. If weighted is true, then returns sum of edge weights along paths connecting leaves of tree
    fn node_distance_matrix(&self, weighted: bool)->HashMap<(NodeID, NodeID), EdgeWeight>{
        let binding = self.get_nodes();
        let leaves = binding.keys().combinations(2);
        let mut dist_mat: HashMap<(NodeID, NodeID), EdgeWeight> = HashMap::new();
        for node_pair in leaves{
            let w  = self.distance_from_node(node_pair[0], node_pair[1], weighted);
            dist_mat.insert((*node_pair[0], *node_pair[1]), w);
        }
        dist_mat
    }

    /// Rerootes tree at given node.
    fn reroot_at_node(&mut self, node_id: &NodeID);
    
    /// Rerootes tree at edge.
    fn reroot_at_edge(&mut self, edge: (&NodeID, &NodeID), edge_weights: (Option<EdgeWeight>, Option<EdgeWeight>)){
        let split_node_id = self.split_edge(edge, edge_weights);
        self.reroot_at_node(&split_node_id);
    }

    /// Inserts node in the middle of edge given by pair of node ids, and returns the new node id
    fn split_edge(&mut self, edge: (&NodeID, &NodeID), edge_weights:(Option<EdgeWeight>, Option<EdgeWeight>))->NodeID;

    /// Returns distance of node from some ancestor of node. If weighted is true, it returns sum of edges from root to self.
    fn distance_from_ancestor(&self, node_id: &NodeID, ancestor: &NodeID, weighted: bool)->f64;

    /// Returns distance of node from root. If weighted is true, it returns sum of edges from root to self.
    fn distance_from_root(&self, node_id: &NodeID, weighted: bool)->EdgeWeight{
        self.distance_from_ancestor(node_id, self.get_root(), weighted)
    }

    /// Returns distance of node from root. If weighted is true, it returns sum of edges from root to self.
    fn distance_from_node(&self, node1: &NodeID, node2: &NodeID, weighted: bool)->f64{
        let mrca = self.get_mrca(vec![node1, node2]);
        self.distance_from_ancestor(node1, &mrca, weighted) + self.distance_from_ancestor(node2, &mrca, weighted)
    }

    /// Returns bipartition induced by edge
    fn get_bipartition(&self, edge: (&NodeID, &NodeID))->(Vec<(NodeID, NodeType)>, Vec<(NodeID, NodeType)>);

    /// Returns cluster of node
    fn get_cluster(&self, node_id: &NodeID)-> Vec<(NodeID, NodeType)>;

    // removes half-nodes in a subtree starting at input node and returns new root of subtree
    fn clean_subtree(&mut self, node_id: NodeID)->(NodeID, Option<f64>);

    /// Cleans self by removing 1) internal nodes (other than root) with degree 2, 2) Floating root nodes, 3) self loops
    fn clean(&mut self);

    /// Get node taxa
    fn get_taxa(&self, node_id:&NodeID)->String;

    /// Get edge weight
    fn get_edge_weight(&self, parent_id: &NodeID, child_id:&NodeID)->Option<&EdgeWeight>{
        for node_id in self.get_node_children(parent_id).iter(){
            if node_id.0==*child_id{
                return node_id.1.as_ref();
            }
        }
        None
    }

    /// return subtree as newick string
    fn subtree_to_newick(&self, node_id:&NodeID, edge_weight:Option<EdgeWeight>)->String{
        fn print_node(node: &NodeType, weight: Option<EdgeWeight>)->String{
            match weight {
                Some(w) => format!("{}:{}", node.taxa(), w),
                None => node.taxa()
            }
        }

        let node  = self.get_node(node_id);
        let mut tmp = String::new();
        if !self.get_node_children(node_id).is_empty(){
            if self.get_node_children(node_id).len()>1{
                tmp.push('(');
            }
            for (child_id, w) in self.get_node_children(node_id){
                let child_str = format!("{},", self.subtree_to_newick(child_id, *w));
                tmp.push_str(&child_str);
            }
            tmp.pop();
            if self.get_node_children(node_id).len()>1{
                tmp.push(')');
            }
        }
        tmp.push_str(&print_node(node, edge_weight));
        tmp
    }

    /// writes full tree in newick format
    fn to_newick(&self)->String{
        format!("{};", self.subtree_to_newick(self.get_root(), None))
    }

    /// Increment all node_ids
    fn incerement_ids(&mut self, value: &usize);
}

pub trait RPhyTree:SimpleRTree {
    
    /// Checks if tree is binary
    fn is_binary(&self)->bool{
        for (_id, children) in self.get_children().iter(){
            match children.len()<=2{
                true => {},
                false => {return false;}
            }
        };
        return true;
    }

    /// SPR function
    fn spr(&mut self, edge1: (&NodeID, &NodeID), edge2: (&NodeID, &NodeID), edge2_weights: (Option<EdgeWeight>, Option<EdgeWeight>)){
        let graft_edge_weight = self.get_edge_weight(edge1.0, edge1.1).cloned();
        let pruned_tree = self.prune(edge1.1);
        self.graft(pruned_tree, edge2, edge2_weights, graft_edge_weight);
    }

    /// Induce tree
    fn induce_tree(&self, taxa: Vec<String>)->Box<dyn RPhyTree>;

    /// Balance a binary tree of 4 taxa
    fn balance_subtree(&mut self){
        assert!(self.is_binary(), "Tree is not binary!");
        assert!(self.get_leaves(self.get_root()).len()==4, "Quartets have 4 leaves!");
        assert!(!self.is_weighted(), "Cannot balance weighted tree!");
        let (root_children, root) = (self.get_node_children(self.get_root()), self.get_root().clone());
        let (child1, child2) = (root_children[0].0, root_children[1].0);
        match (self.is_leaf(&child1), self.is_leaf(&child2)){
            (false, false) => {},
            (true, false) => {
                let other_leaf = self.get_node_children(&child2).iter().filter(|(id, _)| self.is_leaf(id)).collect_vec()[0].0;
                self.spr((&root, &child1), (&child2, &other_leaf), (None, None));
            },
            (false, true) => {
                let other_leaf = self.get_node_children(&child1).iter().filter(|(id, _)| self.is_leaf(id)).collect_vec()[0].0;
                self.spr((&root, &child2), (&child1, &other_leaf), (None, None));
            },
            _ =>{}
        }
        self.clean()
    }

}
