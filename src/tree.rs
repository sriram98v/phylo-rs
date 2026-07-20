/// Module with traits and structs for ancestral sequence reconstruction
pub mod asr;
/// Module with traits and structs for distance computation
pub mod distances;
/// Module with traits and structs for tree encoding
pub mod io;
/// Module with phylogenetic likelihood under a substitution model
pub mod likelihood;
/// Iterative Newick-format parser
#[cfg(feature = "simple_rooted_tree")]
pub(crate) mod newick;
/// Module with traits and structs for tree operations
pub mod ops;
/// Module with traits and structs for general tree traits
pub mod simple_rtree;
/// Module with traits and structs for tree simulation
pub mod simulation;

#[cfg(feature = "simple_rooted_tree")]
pub use simple_rooted_tree::*;

#[cfg(feature = "simple_rooted_tree")]
mod simple_rooted_tree {
    use super::simulation::{Uniform, Yule};

    use itertools::Itertools;
    use rand::prelude::IteratorRandom;

    use crate::iter::{BFSIterator, DFSPostOrderIterator};
    use crate::node::{Node, NodeID};
    use crate::prelude::*;
    use std::fmt::Debug;
    use std::hash::{Hash, Hasher};
    use std::sync::Arc;

    #[cfg(feature = "non_crypto_hash")]
    use fxhash::{FxHashMap as HashMap, FxHashSet as HashSet};
    #[cfg(not(feature = "non_crypto_hash"))]
    use std::collections::{HashMap, HashSet};

    use crate::tree::asr::{JointAsr, MarginalAsr};

    /// Type alias for Phylogenetic tree.
    pub type PhyloTree = SimpleRootedTree<String, f32, f32>;

    impl MarginalAsr for PhyloTree {
        fn marginal_asr<A: Alphabet>(
            &self,
            model: &GtrModel<A>,
            aln: &Alignment,
            want_posteriors: bool,
        ) -> Result<Reconstruction<A>, AsrError> {
            crate::tree::likelihood::compute_marginal_asr(self, model, aln, want_posteriors)
        }
    }

    impl JointAsr for PhyloTree {
        fn joint_asr<A: Alphabet>(
            &self,
            model: &GtrModel<A>,
            aln: &Alignment,
        ) -> Result<Reconstruction<A>, AsrError> {
            crate::tree::likelihood::compute_joint_asr(self, model, aln)
        }
    }

    /// Pointer-based wrapper around `Arc<T>` for use as HashMap key.
    /// Hashes and compares by Arc pointer identity, avoiding content hashing.
    #[derive(Clone, Debug)]
    pub struct TaxaPtr<T>(pub(crate) Arc<T>);

    impl<T> Hash for TaxaPtr<T> {
        fn hash<H: Hasher>(&self, state: &mut H) {
            Arc::as_ptr(&self.0).hash(state);
        }
    }

    impl<T> PartialEq for TaxaPtr<T> {
        fn eq(&self, other: &Self) -> bool {
            Arc::ptr_eq(&self.0, &other.0)
        }
    }

    impl<T> Eq for TaxaPtr<T> {}

    /// Arena memory-managed tree struct
    #[derive(Debug, Clone)]
    pub struct SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        /// Root NodeID.
        ///
        /// Reachable via [`RootedTree::get_root_id`] and [`RootedTree::set_root`].
        root: NodeID,
        /// Nodes of the tree.
        ///
        /// Reachable via [`RootedTree::get_node`], [`RootedTree::get_node_mut`]
        /// and [`RootedTree::get_node_ids`].
        ///
        /// Private because two separate pieces of derived state -- the
        /// `first_free` cursor and `taxa_node_id_map` -- are only correct if
        /// every write to the arena goes through a method that maintains them.
        /// A direct write cannot, and the failure is silent.
        nodes: Vec<Option<Node<T, W, Z>>>,
        /// Index of nodes by taxa.
        ///
        /// Reachable via [`RootedMetaTree::get_taxa_node_id`] and
        /// [`RootedMetaTree::num_taxa`].
        taxa_node_id_map: HashMap<TaxaPtr<T>, NodeID>,
        /// Lower bound on the first vacant arena slot.
        ///
        /// Invariant: every slot below this index is occupied. That makes it a
        /// hint rather than an answer -- the slot at this index may itself be
        /// occupied, so [`Self::next_id`] still scans, but only from here
        /// rather than from zero.
        ///
        /// Kept private: the invariant is what makes `next_id` cheap, and an
        /// external write to `nodes` could silently break it.
        first_free: NodeID,
    }

    impl<T, W, Z> SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        /// Recomputes [`Self::first_free`] from scratch.
        ///
        /// Only for constructors, which take an arena they did not build.
        fn recompute_first_free(&mut self) {
            self.first_free = self
                .nodes
                .iter()
                .position(Option::is_none)
                .unwrap_or(self.nodes.len());
        }

        /// Restores the invariant after slot `node_id` has been filled.
        fn note_slot_filled(&mut self, node_id: NodeID) {
            if node_id == self.first_free {
                // Walk past everything now occupied. Each slot is stepped over
                // at most once per fill, so this is amortised O(1).
                while self.first_free < self.nodes.len() && self.nodes[self.first_free].is_some() {
                    self.first_free += 1;
                }
            }
        }

        /// Restores the invariant after slot `node_id` has been vacated.
        fn note_slot_vacated(&mut self, node_id: NodeID) {
            if node_id < self.first_free {
                self.first_free = node_id;
            }
        }

        /// Drops `node_id`'s taxon from the lookup map, if it has one.
        ///
        /// The map holds an `Arc` to each taxon, so an entry left behind for a
        /// node that is gone keeps the taxon alive, makes `num_taxa` over-count
        /// and lets a stale name resolve to a vacant slot.
        fn forget_taxa(&mut self, node_id: NodeID) {
            let taxa = self
                .nodes
                .get(node_id)
                .and_then(|slot| slot.as_ref())
                .and_then(|node| node.get_taxa_arc())
                .cloned();
            if let Some(arc) = taxa {
                self.taxa_node_id_map.remove(&TaxaPtr(arc));
            }
        }
    }

    impl<T, W, Z> SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        /// Creates new empty tree
        pub fn new(root_id: NodeID) -> Self {
            let root_node = Node::new(root_id);
            let mut nodes = vec![None; root_id + 1];
            nodes[root_id] = Some(root_node);
            let mut tree = SimpleRootedTree {
                root: root_id,
                nodes,
                taxa_node_id_map: [].into_iter().collect::<HashMap<_, _>>(),
                first_free: 0,
            };
            tree.recompute_first_free();
            tree
        }

        /// Creates tree with specified capacity
        pub fn with_capacity(capacity: usize) -> Self {
            let root_node = Node::new(0);
            let mut nodes = vec![None; capacity];
            nodes[0] = Some(root_node);
            let mut tree = SimpleRootedTree {
                root: 0,
                nodes,
                taxa_node_id_map: [].into_iter().collect::<HashMap<_, _>>(),
                first_free: 0,
            };
            tree.recompute_first_free();
            tree
        }

        /// Returns the lowest vacant arena slot, or the arena length if full.
        ///
        /// Takes `&self`, so it does not consume the slot: calling it twice
        /// without an intervening insert returns the same id both times.
        /// Callers rely on that.
        ///
        /// Scans from [`Self::first_free`] rather than from zero. Since every
        /// slot below that is occupied, the result is identical to a full scan,
        /// but building a tree no longer re-walks the whole arena per node.
        pub fn next_id(&self) -> usize {
            self.nodes[self.first_free..]
                .iter()
                .position(Option::is_none)
                .map_or(self.nodes.len(), |offset| self.first_free + offset)
        }

        /// Creates new node with next NodeID
        pub fn next_node(&self) -> Node<T, W, Z> {
            Node::new(self.next_id())
        }

        /// returns max number of nodes in tree without reallocating node vec
        pub fn get_capacity(&self) -> usize {
            self.nodes.len()
        }

        /// Returns the bytes allocated on the heap by the node arena.
        ///
        /// This counts the arena at capacity, so vacant slots are included:
        /// the arena never shrinks, and a hole costs a full slot.
        pub fn arena_heap_size(&self) -> usize {
            self.nodes.capacity() * std::mem::size_of::<Option<Node<T, W, Z>>>()
                + self
                    .nodes
                    .iter()
                    .flatten()
                    .map(|node| node.heap_size())
                    .sum::<usize>()
        }

        /// Returns the bytes allocated on the heap by the taxa map.
        ///
        /// Excludes the taxa values themselves, which live behind [`Arc`] and
        /// whose size depends on `T`.
        pub fn taxa_map_heap_size(&self) -> usize {
            self.taxa_node_id_map.capacity()
                * (std::mem::size_of::<TaxaPtr<T>>() + std::mem::size_of::<NodeID>())
        }

        /// Returns the total bytes this tree has allocated on the heap.
        ///
        /// The sum of [`Self::arena_heap_size`] and [`Self::taxa_map_heap_size`].
        /// Taxa values are excluded, as described on those methods. The LCA
        /// index is no longer stored on the tree; measure it via
        /// [`LcaOracle::heap_size`] on an oracle built with [`EulerWalk::lca`].
        pub fn heap_size(&self) -> usize {
            self.arena_heap_size() + self.taxa_map_heap_size()
        }
    }

    impl<T, W, Z> RootedTree for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        type Node = Node<T, W, Z>;

        /// Creates new empty tree
        fn new() -> Self {
            let root_node = Node::new(0);
            let mut nodes = vec![None; 1];
            nodes[0] = Some(root_node);
            let mut tree = SimpleRootedTree {
                root: 0,
                nodes,
                taxa_node_id_map: [].into_iter().collect::<HashMap<_, _>>(),
                first_free: 0,
            };
            tree.recompute_first_free();
            tree
        }

        /// Creates tree with specified capacity
        fn with_capacity(capacity: usize) -> Self {
            let root_node = Node::new(0);
            let mut nodes = vec![None; capacity];
            nodes[0] = Some(root_node);
            let mut tree = SimpleRootedTree {
                root: 0,
                nodes,
                taxa_node_id_map: [].into_iter().collect::<HashMap<_, _>>(),
                first_free: 0,
            };
            tree.recompute_first_free();
            tree
        }

        fn from_nodes(nodes: Vec<Option<Self::Node>>, root_id: TreeNodeID<Self>) -> Self {
            // Rebuild taxa_node_id_map from the nodes, exactly as `Deserialize`
            // does. `from_nodes` is public, so a caller passing taxa-bearing
            // nodes must get working taxa lookups (num_taxa, get_taxa_node_id,
            // ...) without a separate re-register pass; leaving the map empty
            // would make those queries silently return nothing.
            let mut taxa_node_id_map: HashMap<TaxaPtr<T>, NodeID> = [].into_iter().collect();
            for node in nodes.iter().flatten() {
                if let Some(arc) = node.get_taxa_arc() {
                    taxa_node_id_map.insert(TaxaPtr(arc.clone()), node.get_id());
                }
            }
            let mut tree = SimpleRootedTree {
                root: root_id,
                nodes,
                taxa_node_id_map,
                first_free: 0,
            };
            tree.recompute_first_free();
            tree
        }

        /// Returns reference to node by ID
        fn get_node(&self, node_id: TreeNodeID<Self>) -> Option<&Node<T, W, Z>> {
            self.nodes[node_id].as_ref()
        }

        fn get_node_mut(&mut self, node_id: TreeNodeID<Self>) -> Option<&mut Node<T, W, Z>> {
            self.nodes[node_id].as_mut()
        }

        fn get_node_ids(&self) -> impl Iterator<Item = TreeNodeID<Self>> {
            (0..self.nodes.len()).filter(|x| self.nodes[*x].is_some())
        }

        fn get_nodes_mut(&mut self) -> impl Iterator<Item = &mut Self::Node> {
            self.nodes.iter_mut().filter_map(|x| x.as_mut())
        }

        fn set_node(&mut self, node: Node<T, W, Z>) {
            let node_id = node.get_id();
            // Whatever was in this slot is about to be replaced; its taxon must
            // not outlive it in the map.
            self.forget_taxa(node_id);
            if let Some(arc) = node.get_taxa_arc() {
                self.taxa_node_id_map
                    .insert(TaxaPtr(arc.clone()), node.get_id());
            }
            match self.nodes.len() > node_id {
                true => self.nodes[node_id] = Some(node),
                false => {
                    let new_len = node.get_id() - self.nodes.len();
                    self.nodes.extend((0..new_len + 1).map(|_| None));
                    self.nodes[node_id] = Some(node);
                }
            };
            self.note_slot_filled(node_id);
        }

        fn get_root_id(&self) -> TreeNodeID<Self> {
            self.root
        }

        fn set_root(&mut self, node_id: TreeNodeID<Self>) {
            self.root = node_id;
        }

        fn remove_node(&mut self, node_id: TreeNodeID<Self>) -> Option<Node<T, W, Z>> {
            if let Some(pid) = self.get_node_parent_id(node_id) {
                // The parent may already be gone: a node keeps recording its
                // parent id after that parent is removed, so unwrapping here
                // panics when a node outlives its parent. If there is no
                // parent, there is no child link to unlink.
                if let Some(parent) = self.get_node_mut(pid) {
                    parent.remove_child(&node_id)
                }
            }
            self.forget_taxa(node_id);
            let removed = self.nodes[node_id].take();
            self.note_slot_vacated(node_id);
            removed
        }

        fn delete_node(&mut self, node_id: TreeNodeID<Self>) {
            self.forget_taxa(node_id);
            let _ = self.nodes[node_id].take();
            self.note_slot_vacated(node_id);
        }

        fn clear(&mut self) {
            let root_node = self.get_root().clone();
            let root_node_id = root_node.get_id();
            self.nodes = vec![None; root_node_id + 1];
            self.nodes[root_node_id] = Some(root_node);
            self.recompute_first_free();
            self.taxa_node_id_map.clear();
        }

        /// Supresses all nodes of degree 2
        fn supress_unifurcations(&mut self) {
            let post_ord_node_ids = self.postord_ids(self.get_root_id()).collect_vec();
            for node_id in post_ord_node_ids {
                if !self.is_leaf(node_id) && node_id != self.root {
                    let node_degree = self.node_degree(node_id);
                    if node_degree == 2 {
                        let node_parent_id = self.get_node_parent_id(node_id).unwrap();
                        let node_child_id = self.get_node_children_ids(node_id).next().unwrap();
                        self.remove_node(node_id);
                        self.set_child(node_parent_id, node_child_id);
                    }
                }
            }
        }
    }

    impl<T, W, Z> RootedMetaTree for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn get_taxa_node(&self, taxa: &TreeNodeMeta<Self>) -> Option<&Self::Node> {
            let node_id = self
                .taxa_node_id_map
                .iter()
                .find(|(tp, _)| tp.0.as_ref() == taxa)
                .map(|(_, id)| *id)?;
            self.get_node(node_id)
        }

        fn set_node_taxa(&mut self, node_id: TreeNodeID<Self>, taxa: Option<TreeNodeMeta<Self>>) {
            // The map is keyed by Arc identity, so a new label inserts a new
            // entry rather than overwriting one. Without dropping the old key
            // first, relabelling grows the map without bound and the previous
            // name goes on resolving.
            self.forget_taxa(node_id);
            if let Some(t) = taxa {
                let arc = Arc::new(t);
                self.get_node_mut(node_id)
                    .unwrap()
                    .set_taxa_arc(Some(arc.clone()));
                self.taxa_node_id_map.insert(TaxaPtr(arc), node_id);
            } else {
                self.get_node_mut(node_id).unwrap().set_taxa(None);
            }
        }

        fn num_taxa(&self) -> usize {
            self.taxa_node_id_map.len()
        }

        fn get_taxa_space(&self) -> impl Iterator<Item = &TreeNodeMeta<Self>> {
            self.taxa_node_id_map.keys().map(|tp| tp.0.as_ref())
        }

        fn get_node_taxa_cloned(&self, node_id: TreeNodeID<Self>) -> Option<TreeNodeMeta<Self>> {
            self.get_node(node_id).unwrap().get_taxa().cloned()
        }
    }

    impl<T, W, Z> Yule for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn yule(num_taxa: usize) -> SimpleRootedTree<T, W, Z> {
            let mut tree = SimpleRootedTree::new(0);
            if num_taxa < 3 {
                return tree;
            }
            let new_node = Node::new(1);
            tree.add_child(0, new_node);
            tree.set_node_taxa(1, T::from_str("0").ok());
            let new_node = Node::new(2);
            tree.add_child(0, new_node);
            tree.set_node_taxa(2, T::from_str("1").ok());
            if num_taxa < 4 {
                return tree;
            }
            let mut current_leaf_ids = vec![1, 2];
            for i in 2..num_taxa {
                let rand_leaf_id = current_leaf_ids
                    .iter()
                    .choose(&mut rand::thread_rng())
                    .unwrap();
                let rand_leaf_parent_id = tree.get_node_parent_id(*rand_leaf_id).unwrap();
                let split_node = Node::new(tree.next_id());
                let split_node_id = split_node.get_id();
                tree.split_edge((rand_leaf_parent_id, *rand_leaf_id), split_node);
                let new_leaf = Node::new(tree.next_id());
                let new_leaf_id = new_leaf.get_id();
                tree.add_child(split_node_id, new_leaf);
                tree.set_node_taxa(new_leaf_id, T::from_str(&i.to_string()).ok());
                current_leaf_ids.push(new_leaf_id);
            }
            tree
        }
    }

    impl<T, W, Z> Uniform for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn unif(num_taxa: usize) -> SimpleRootedTree<T, W, Z> {
            let mut tree = SimpleRootedTree::new(0);
            if num_taxa < 3 {
                return tree;
            }
            let new_node = Node::new(1);
            tree.add_child(0, new_node);
            tree.set_node_taxa(1, T::from_str("0").ok());
            let new_node = Node::new(2);
            tree.add_child(0, new_node);
            tree.set_node_taxa(2, T::from_str("1").ok());
            if num_taxa < 3 {
                return tree;
            }
            let mut current_node_ids = vec![1, 2];
            for i in 2..num_taxa {
                let rand_leaf_id = *current_node_ids
                    .iter()
                    .choose(&mut rand::thread_rng())
                    .unwrap();
                let rand_leaf_parent_id = tree.get_node_parent_id(rand_leaf_id).unwrap();
                let split_node = Node::new(tree.next_id());
                let split_node_id = split_node.get_id();
                current_node_ids.push(split_node_id);
                tree.split_edge((rand_leaf_parent_id, rand_leaf_id), split_node);
                let new_leaf = Node::new(tree.next_id());
                let new_leaf_id = new_leaf.get_id();
                tree.add_child(split_node_id, new_leaf);
                tree.set_node_taxa(new_leaf_id, T::from_str(&i.to_string()).ok());
                current_node_ids.push(new_leaf_id);
            }
            tree
        }
    }

    impl<T, W, Z> RootedWeightedTree for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn unweight(&mut self) {
            // `flatten` yields the occupied slots; the previous body filtered
            // for empty slots and unwrapped them, so it panicked on the first
            // hole and never touched a real node.
            self.nodes
                .iter_mut()
                .flatten()
                .for_each(|node| node.unweight());
        }
    }

    impl<T, W, Z> PathFunction for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> Ancestors for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> Subtree for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> PreOrder for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> ClusterMatching for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> ClusterAffinity for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> RobinsonFoulds for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> DistanceMatrix for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn matrix(&self) -> Vec<Vec<TreeNodeWeight<Self>>> {
            // One euler-tour index, shared across every pair, rather than the
            // naive per-call rebuild the old signature forced.
            let oracle = self.lca();
            let mut out_mat = vec![vec![W::infinity(); self.nodes.len()]; self.nodes.len()];
            for node_ids in self.postord_ids(self.get_root_id()).combinations(2) {
                let n1 = node_ids[0];
                let n2 = node_ids[1];
                out_mat[n1][n1] = W::zero();
                out_mat[n2][n2] = W::zero();
                out_mat[n1][n2] = self.pairwise_distance(&oracle, n1, n2);
                out_mat[n2][n1] = out_mat[n1][n2];
            }
            out_mat
        }

        fn pairwise_distance(
            &self,
            oracle: &LcaOracle<'_, Self>,
            node_id_1: TreeNodeID<Self>,
            node_id_2: TreeNodeID<Self>,
        ) -> TreeNodeWeight<Self> {
            let lca = oracle.get_lca_id(vec![node_id_1, node_id_2].as_slice());
            let d1: TreeNodeWeight<Self> = self
                .node_to_root_ids(node_id_1)
                .map(|x| match x == self.get_root_id() {
                    true => W::zero(),
                    false => self.get_edge_weight(0, x).unwrap_or(W::one()),
                })
                .sum();

            let d2: TreeNodeWeight<Self> = self
                .node_to_root_ids(node_id_2)
                .map(|x| match x == self.get_root_id() {
                    true => W::zero(),
                    false => self.get_edge_weight(0, x).unwrap_or(W::one()),
                })
                .sum();

            let dlca: TreeNodeWeight<Self> = self
                .node_to_root_ids(lca)
                .map(|x| match x == self.get_root_id() {
                    true => W::zero(),
                    false => self.get_edge_weight(0, x).unwrap_or(W::one()),
                })
                .sum();

            d1 + d2 - (W::one() + W::one()) * dlca
        }
    }

    impl<T, W, Z> DFS for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn postord_ids(
            &self,
            start_node: TreeNodeID<Self>,
        ) -> impl Iterator<Item = TreeNodeID<Self>> {
            DFSPostOrderIterator::new(self, start_node).map(|x| x.get_id())
        }

        fn postord_nodes(&self, start_node: TreeNodeID<Self>) -> impl Iterator<Item = &Self::Node> {
            DFSPostOrderIterator::new(self, start_node)
        }
    }

    impl<T, W, Z> BFS for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn bfs_nodes(&self, start_node_id: TreeNodeID<Self>) -> impl Iterator<Item = &Self::Node> {
            BFSIterator::new(self, start_node_id)
        }

        fn bfs_ids(
            &self,
            start_node_id: TreeNodeID<Self>,
        ) -> impl Iterator<Item = TreeNodeID<Self>> {
            BFSIterator::new(self, start_node_id).map(|x| x.get_id())
        }
    }

    impl<T, W, Z> ContractTree for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn contracted_tree_nodes(
            &self,
            leaf_ids: &[TreeNodeID<Self>],
        ) -> impl Iterator<Item = Self::Node> {
            let new_tree_root_id = self.get_lca_id(leaf_ids);
            let node_postord_iter = self.postord_nodes(new_tree_root_id);
            let mut node_map: Vec<Option<Self::Node>> = vec![None; self.nodes.len()];
            node_map[new_tree_root_id] = Some(self.get_lca(leaf_ids).clone());
            let mut leaf_id_set = vec![false; self.nodes.len()];
            for id in leaf_ids {
                leaf_id_set[*id] = true;
            }
            let mut remove_list = vec![false; self.nodes.len()];
            node_postord_iter.for_each(|orig_node| {
                let mut node = orig_node.clone();
                match node.is_leaf() {
                    true => {
                        if leaf_id_set[node.get_id()] {
                            node_map[node.get_id()] = Some(node.clone());
                        }
                    }
                    false => {
                        let node_children_ids = node.get_children().to_vec();
                        for child_id in node_children_ids.iter() {
                            match node_map[*child_id].is_some() {
                                true => {}
                                false => node.remove_child(child_id),
                            }
                        }
                        let node_children_ids = node.get_children().to_vec();
                        match node_children_ids.len() {
                            0 => {}
                            1 => {
                                // the node is a unifurcation
                                // node should be added to both node_map and remove_list
                                // if child of node is already in remove list, attach node children to node first
                                let child_node_id = node_children_ids[0];
                                let child_node_edge_weight = self
                                    .get_node(child_node_id)
                                    .unwrap()
                                    .get_weight()
                                    .unwrap_or(W::zero());

                                if remove_list[child_node_id] {
                                    node.remove_child(&child_node_id);
                                    let grandchildren_ids = node_map[child_node_id]
                                        .as_mut()
                                        .unwrap()
                                        .get_children()
                                        .to_vec();
                                    for grandchild_id in grandchildren_ids {
                                        node_map[grandchild_id]
                                            .as_mut()
                                            .unwrap()
                                            .set_parent(Some(node.get_id()));
                                        let new_edge_weight = node_map[grandchild_id]
                                            .as_ref()
                                            .unwrap()
                                            .get_weight()
                                            .unwrap_or(W::zero())
                                            + child_node_edge_weight;
                                        node_map[grandchild_id]
                                            .as_mut()
                                            .unwrap()
                                            .set_weight(Some(new_edge_weight));
                                        node.add_child(grandchild_id);
                                    }
                                }
                                let n_id = node.get_id();
                                remove_list[n_id] = true;
                                node_map[n_id] = Some(node.clone());
                            }
                            _ => {
                                // node has multiple children
                                // for each child, suppress child if child is in remove list
                                node_children_ids.into_iter().for_each(|chid| {
                                    if remove_list[chid] {
                                        // suppress chid
                                        // remove chid from node children
                                        // children of chid are node grandchildren
                                        // add grandchildren to node children
                                        // set grandchildren parent to node
                                        let chid_weight = self
                                            .get_node(chid)
                                            .unwrap()
                                            .get_weight()
                                            .unwrap_or(W::zero());
                                        node.remove_child(&chid);
                                        let node_grandchildren = node_map[chid]
                                            .as_mut()
                                            .unwrap()
                                            .get_children()
                                            .to_vec();
                                        for grandchild_id in node_grandchildren {
                                            let new_edge_weight = node_map[grandchild_id]
                                                .as_ref()
                                                .unwrap()
                                                .get_weight()
                                                .unwrap_or(W::zero())
                                                + chid_weight;
                                            node.add_child(grandchild_id);
                                            node_map[grandchild_id]
                                                .as_mut()
                                                .unwrap()
                                                .set_parent(Some(node.get_id()));
                                            node_map[grandchild_id]
                                                .as_mut()
                                                .unwrap()
                                                .set_weight(Some(new_edge_weight));
                                        }
                                    }
                                });
                                if node.get_id() == new_tree_root_id {
                                    node.set_parent(None);
                                }
                                node_map[node.get_id()] = Some(node.clone());
                            }
                        };
                    }
                }
            });
            remove_list.into_iter().enumerate().for_each(|(n_id, x)| {
                if x {
                    node_map[n_id] = None;
                }
            });
            node_map.into_iter().flatten()
        }

        fn contract_tree(&self, leaf_ids: &[TreeNodeID<Self>]) -> Result<Self, ()> {
            let new_tree_root_id = self.get_lca_id(leaf_ids);
            let new_nodes = self.contracted_tree_nodes(leaf_ids);
            let mut new_tree = SimpleRootedTree {
                root: new_tree_root_id,
                nodes: vec![None; self.nodes.len()],
                taxa_node_id_map: vec![].into_iter().collect(),
                // The arena starts wholly vacant; `set_nodes` fills it and
                // maintains the invariant from there.
                first_free: 0,
            };
            new_tree.set_nodes(new_nodes);
            Ok(new_tree)
        }

        fn contract_tree_from_iter(
            &self,
            leaf_ids: &[TreeNodeID<Self>],
            node_iter: impl Iterator<Item = TreeNodeID<Self>>,
        ) -> Result<Self, ()> {
            let new_tree_root_id = self.get_lca_id(leaf_ids);
            let new_nodes =
                self.contracted_tree_nodes_from_iter(new_tree_root_id, leaf_ids, node_iter);
            let mut new_tree = SimpleRootedTree {
                root: new_tree_root_id,
                nodes: vec![None; self.nodes.len()],
                taxa_node_id_map: vec![].into_iter().collect(),
                // The arena starts wholly vacant; `set_nodes` fills it and
                // maintains the invariant from there.
                first_free: 0,
            };
            new_tree.set_nodes(new_nodes);
            Ok(new_tree)
        }
    }

    // All `EulerWalk` methods are trait defaults; the euler walks read only
    // topology and the LCA index now lives in a borrowing [`LcaOracle`] built
    // by [`EulerWalk::lca`], so this tree needs no overrides.
    impl<T, W, Z> EulerWalk for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> Clusters for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn get_median_node_id_for_leaves(
            &self,
            taxa_set: impl Iterator<Item = TreeNodeID<Self>>,
        ) -> TreeNodeID<Self> {
            let mut cluster_sizes = vec![0; self.nodes.len()];
            let mut median_node_id: TreeNodeID<Self> = self.get_root_id();
            let leaf_ids: HashSet<TreeNodeID<Self>> = taxa_set.collect();
            for n_id in self.postord_ids(self.get_root_id()) {
                if self.is_leaf(n_id) && leaf_ids.contains(&n_id) {
                    cluster_sizes[n_id] = 1;
                } else {
                    for c_id in self.get_node_children_ids(n_id) {
                        cluster_sizes[n_id] += cluster_sizes[c_id];
                    }
                }
            }
            loop {
                median_node_id = self
                    .get_node_children_ids(median_node_id)
                    .max_by(|x, y| {
                        let x_cluster_size = cluster_sizes[*x];
                        let y_cluster_size = cluster_sizes[*y];
                        x_cluster_size.cmp(&y_cluster_size)
                    })
                    .unwrap();
                if cluster_sizes[median_node_id] <= (leaf_ids.len() / 2) {
                    break;
                }
            }
            median_node_id
        }

        fn get_median_node_for_leaves(
            &self,
            taxa_set: impl Iterator<Item = TreeNodeID<Self>>,
        ) -> &Self::Node {
            self.get_node(self.get_median_node_id_for_leaves(taxa_set))
                .unwrap()
        }

        /// Returns an immutable reference to median node of all leaves in a tree.
        fn get_median_node(&self) -> &Self::Node {
            let leaves = self.get_leaves().map(|x| x.get_id());
            self.get_median_node_for_leaves(leaves)
        }

        /// Returns median Node<T,W,Z>ID of all leaves in a tree.
        fn get_median_node_id(&self) -> TreeNodeID<Self> {
            let leaves = self.get_leaf_ids();
            self.get_median_node_id_for_leaves(leaves)
        }
    }

    impl<T, W, Z> Newick for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn from_newick_with<H: AnnotationHandler>(
            newick_str: &[u8],
            annotations: H,
        ) -> std::io::Result<Self> {
            let input = std::str::from_utf8(newick_str)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))?;
            crate::tree::newick::parse_newick(input, &annotations)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
        }

        fn subtree_to_newick_with<H: AnnotationWriter>(
            &self,
            node_id: TreeNodeID<Self>,
            annotations: H,
        ) -> impl std::fmt::Display {
            // Iterative Euler-tour walk writing into a single buffer: no
            // recursion (so tree height cannot overflow the stack) and no
            // per-subtree string allocation (so the cost is O(n), not the
            // O(n * height) of rebuilding a string at every level). Each stack
            // frame is (node id, index of the next child to emit).
            let mut out = String::new();
            let mut stack: Vec<(TreeNodeID<Self>, usize)> = vec![(node_id, 0)];
            while let Some(&(nid, child_idx)) = stack.last() {
                let node = self.get_node(nid).unwrap();
                let children = node.get_children();
                // Parenthesise only nodes with more than one child, matching the
                // original serialiser (a unary node is flattened).
                if child_idx == 0 && children.len() > 1 {
                    out.push('(');
                }
                if child_idx < children.len() {
                    if child_idx > 0 {
                        out.push(',');
                    }
                    let child = children[child_idx];
                    stack.last_mut().unwrap().1 += 1;
                    stack.push((child, 0));
                } else {
                    // All children emitted: close the subtree, then this node's
                    // own label and branch length.
                    if children.len() > 1 {
                        out.push(')');
                    }
                    if let Some(taxa_str) = &node.get_taxa() {
                        out.push_str(&taxa_str.to_string());
                    }
                    // The writer decides what (if anything) a stored annotation
                    // contributes; the default emits it verbatim after the label.
                    if let Some(annotation) = node.get_annotation() {
                        if let Some(rendered) = annotations.render(annotation) {
                            out.push_str(&rendered);
                        }
                    }
                    if let Some(w) = node.get_weight() {
                        out.push(':');
                        out.push_str(&w.to_string());
                    }
                    stack.pop();
                }
            }
            out
        }
    }

    impl<T, W, Z> Nexus for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    impl<T, W, Z> SPR for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn graft(
            &mut self,
            mut tree: Self,
            edge: (TreeNodeID<Self>, TreeNodeID<Self>),
        ) -> Result<(), ()> {
            let new_node = self.next_node();
            let new_node_id = new_node.get_id();
            for node in tree.get_nodes_mut() {
                node.set_id(self.next_node().get_id());
                self.set_node(node.clone());
            }
            self.split_edge(edge, new_node);
            self.set_child(new_node_id, tree.get_root_id());
            Ok(())
        }
        fn prune(&mut self, node_id: TreeNodeID<Self>) -> Result<Self, ()> {
            let mut pruned_tree = SimpleRootedTree::new(node_id);
            let p_id = self.get_node_parent_id(node_id).unwrap();
            self.get_node_mut(p_id).unwrap().remove_child(&node_id);
            pruned_tree
                .get_node_mut(pruned_tree.get_root_id())
                .unwrap()
                .add_children(
                    self.get_node(node_id)
                        .unwrap()
                        .get_children()
                        .iter()
                        .copied(),
                );
            let dfs = self.dfs(node_id).collect_vec();
            for node in dfs {
                // self.nodes.remove(node.get_id());
                pruned_tree.set_node(node.clone());
            }
            Ok(pruned_tree)
        }
    }

    impl<T, W, Z> NNI for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn nni(&mut self, node_id: TreeNodeID<Self>, left_ch: bool) -> Result<(), ()> {
            if self.is_leaf(node_id) || node_id == self.get_root_id() {
                panic!("NNI cannot be performed at a leaf or root!")
            } else {
                let node_parent_id = self.get_node_parent_id(node_id).unwrap();

                let node_ch_ids = self.get_node_children_ids(node_id).collect_vec();
                let node_ch1 = node_ch_ids[left_ch as usize];
                let node_sibling = self
                    .get_node_children_ids(node_parent_id)
                    .filter(|x| x != &node_id)
                    .collect_vec()[0];

                // set node_ch2 as sibling to parent node
                self.delete_edge(node_id, node_ch1);
                self.delete_edge(node_parent_id, node_sibling);

                self.set_child(node_parent_id, node_ch1);
                self.set_child(node_id, node_sibling);

                Ok(())
            }
        }
    }

    impl<T, W, Z> Balance for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        fn balance_subtree(&mut self) -> Result<(), ()> {
            assert!(
                self.get_cluster(self.get_root_id()).collect_vec().len() == 4,
                "Quartets have 4 leaves!"
            );
            assert!(self.is_binary(), "Cannot balance non-binary tree!");
            let root_children = self.get_node_children(self.get_root_id()).collect_vec();
            let (child1, child2) = (root_children[0].get_id(), root_children[1].get_id());
            let next_id = self.next_id();
            let split_id = self.next_id() + 1;
            match dbg!((
                (self.get_node(child1).unwrap().is_leaf()),
                (self.get_node(child2).unwrap().is_leaf())
            )) {
                (false, false) => {}
                (true, false) => {
                    let mut leaf_node = self.remove_node(child1).unwrap();
                    leaf_node.set_id(next_id);
                    let other_leaf_id = &self
                        .get_node_children(child2)
                        .filter(|node| node.is_leaf())
                        .collect_vec()[0]
                        .get_id();
                    self.split_edge((child2, *other_leaf_id), Node::new(split_id));
                    self.add_child(dbg!(split_id), leaf_node);
                }
                (false, true) => {
                    let mut leaf_node = self.remove_node(child2).unwrap();
                    leaf_node.set_id(next_id);
                    let other_leaf_id = &self
                        .get_node_children(child1)
                        .filter(|node| node.is_leaf())
                        .collect_vec()[0]
                        .get_id();
                    self.split_edge((child1, *other_leaf_id), Node::new(split_id));
                    self.add_child(split_id, leaf_node);
                }
                _ => {}
            }
            self.clean();
            Ok(())
        }
    }

    impl<T, W, Z> CopheneticDistance for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
    }

    #[cfg(feature = "serde")]
    impl<T, W, Z> serde::Serialize for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa + serde::Serialize,
        W: EdgeWeight + serde::Serialize,
        Z: NodeWeight + serde::Serialize,
    {
        fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
        where
            S: serde::Serializer,
        {
            use serde::ser::SerializeStruct;
            let mut state = serializer.serialize_struct("SimpleRootedTree", 2)?;
            state.serialize_field("root", &self.root)?;
            state.serialize_field("nodes", &self.nodes)?;
            state.end()
        }
    }

    #[cfg(feature = "serde")]
    impl<'de, T, W, Z> serde::Deserialize<'de> for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa + serde::Deserialize<'de>,
        W: EdgeWeight + serde::Deserialize<'de>,
        Z: NodeWeight + serde::Deserialize<'de>,
    {
        fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
        where
            D: serde::Deserializer<'de>,
        {
            #[derive(serde::Deserialize)]
            struct Helper<T, W, Z>
            where
                T: NodeTaxa,
                W: EdgeWeight,
                Z: NodeWeight,
            {
                root: NodeID,
                nodes: Vec<Option<Node<T, W, Z>>>,
            }

            let helper: Helper<T, W, Z> = Helper::deserialize(deserializer)?;

            // Rebuild taxa_node_id_map from node data
            let mut taxa_node_id_map: HashMap<TaxaPtr<T>, NodeID> = [].into_iter().collect();
            for node in helper.nodes.iter().flatten() {
                if let Some(arc) = node.get_taxa_arc() {
                    taxa_node_id_map.insert(TaxaPtr(arc.clone()), node.get_id());
                }
            }

            let mut tree = SimpleRootedTree {
                root: helper.root,
                nodes: helper.nodes,
                taxa_node_id_map,
                first_free: 0,
            };
            // Derived from the arena, exactly like `taxa_node_id_map` above, so
            // it stays out of the serialized form and is rebuilt on the way in.
            tree.recompute_first_free();
            Ok(tree)
        }
    }

    impl<T, W, Z> OLA for SimpleRootedTree<T, W, Z>
    where
        T: NodeTaxa,
        W: EdgeWeight,
        Z: NodeWeight,
    {
        /// Decodes an OLATree into a rooted binary tree.
        ///
        /// Leaf ordering σ is taken from `ola.taxa`: leaf l_j has index j.
        /// Each `ola.indices[i-1]` identifies the sibling of l_i in the
        /// restricted tree T^i — a non-negative value is a leaf index, a
        /// negative value is an internal node index.
        fn from_vec(ola: OLATree<T>) -> Self {
            let n = ola.taxa.len();

            if n == 0 {
                return SimpleRootedTree::new(0);
            }

            // For a binary tree on n leaves: n leaves + n-1 internal nodes.
            // Node ID assignment:
            //   Leaf l_j          → NodeID j          (j = 0..n-1)
            //   Internal node I_i → NodeID n + i - 1  (i = 1..n-1)
            let capacity = if n > 1 { 2 * n - 1 } else { 1 };
            let mut nodes: Vec<Option<Node<T, W, Z>>> = vec![None; capacity];

            // Pre-create all leaf nodes (parents set during the loop below)
            #[allow(clippy::needless_range_loop)]
            for j in 0..n {
                nodes[j] = Some(Node::new(j));
            }

            // Build the tree structure by replaying the OLA construction
            let mut root_id: NodeID = 0; // starts as just l₀

            for i in 1..n {
                let e = ola.indices[i - 1];

                // Map OLA entry to the sibling's NodeID
                let sibling_id: NodeID = if e >= 0 {
                    e as NodeID
                } else {
                    // OLA index -k → internal node I_k → NodeID n + k - 1
                    let k = (-e) as usize;
                    n + k - 1
                };

                let internal_id: NodeID = n + i - 1;
                let mut internal_node = Node::new(internal_id);

                // Wire the new internal node in place of the sibling
                match nodes[sibling_id].as_ref().unwrap().get_parent() {
                    None => {
                        // Sibling is the current root; internal becomes the new root
                        root_id = internal_id;
                    }
                    Some(p_id) => {
                        nodes[p_id].as_mut().unwrap().remove_child(&sibling_id);
                        nodes[p_id].as_mut().unwrap().add_child(internal_id);
                        internal_node.set_parent(Some(p_id));
                    }
                }

                // Internal node's children: existing sibling and the new leaf l_i
                internal_node.add_child(sibling_id);
                internal_node.add_child(i);
                nodes[sibling_id]
                    .as_mut()
                    .unwrap()
                    .set_parent(Some(internal_id));
                nodes[i].as_mut().unwrap().set_parent(Some(internal_id));

                nodes[internal_id] = Some(internal_node);
            }

            // Assemble the tree; taxa_node_id_map is populated below
            let mut tree = SimpleRootedTree::from_nodes(nodes, root_id);

            // Register taxa using set_node_taxa so that taxa_node_id_map is populated
            #[allow(clippy::needless_range_loop)]
            for j in 0..n {
                tree.set_node_taxa(j, Some(ola.taxa[j].clone()));
            }

            tree
        }
    }
}
