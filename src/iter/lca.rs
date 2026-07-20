//! Constant-time lowest-common-ancestor queries over an immutably borrowed tree.
//!
//! [`LcaOracle`] owns a precomputed euler-tour/RMQ index and holds a shared
//! reference to the tree it was built from. Because the borrow is shared, the
//! tree cannot be mutated while an oracle is alive -- the staleness that the
//! old in-tree index guarded against at runtime is now a borrow-check error.
//! Build one with [`EulerWalk::lca`](crate::iter::node_iter::EulerWalk::lca),
//! query it, and drop it before mutating the tree again.

use itertools::Itertools;
use vers_vecs::FastRmq;

use crate::iter::node_iter::EulerWalk;
use crate::tree::simple_rtree::TreeNodeID;

/// Marks a node id that never appears in the euler walk, i.e. an arena slot
/// holding no node. The walk is at most `2n - 1` long, so a real position can
/// never reach `u32::MAX`.
const NO_FIRST_APPEARANCE: u32 = u32::MAX;

/// A constant-time LCA oracle borrowing a tree immutably.
///
/// Owns its euler tour, first-appearance index, depth array and range-minimum
/// query; see the module docs for the borrowing contract.
pub struct LcaOracle<'t, Tree: EulerWalk> {
    /// The tree the index was built from. Shared borrow: the tree is frozen for
    /// the oracle's lifetime.
    tree: &'t Tree,
    /// Precomputed euler tour.
    euler: Vec<TreeNodeID<Tree>>,
    /// First-appearance position of each node in `euler`, indexed by node id.
    ///
    /// [`NO_FIRST_APPEARANCE`] marks an id with no node. A `u32` with a sentinel
    /// rather than `Option<usize>`, which costs 16 bytes per entry: `usize` has
    /// no spare bit pattern, so the `Option` cannot pack, and the arena is
    /// indexed by `u32` positions anyway.
    fai: Vec<u32>,
    /// Depth of every entry in `euler`.
    da: Vec<usize>,
    /// Range-minimum-query over `da`.
    ///
    /// `FastRmq` rather than `BinaryRmq`: the latter is a sparse table over
    /// every element, which for a tree of any size is both the largest thing the
    /// library allocates and too big to stay in cache. `FastRmq` keeps a sparse
    /// table over block minima instead.
    rmq: FastRmq,
}

impl<'t, Tree: EulerWalk> LcaOracle<'t, Tree> {
    /// Builds the index from a shared borrow of `tree`.
    ///
    /// Reads only topology, so it takes `&Tree`: the index is materialised into
    /// fresh local vectors and never written back into the tree.
    pub(crate) fn build(tree: &'t Tree) -> Self {
        let euler: Vec<TreeNodeID<Tree>> = tree.euler_walk_ids(tree.get_root_id()).collect_vec();

        let max_id: usize = tree
            .get_node_ids()
            .map(Into::<usize>::into)
            .max()
            .expect("a tree always has at least a root node");

        // First appearance, in one pass over the walk: the earliest position
        // wins. Searching the walk for each node instead would give the same
        // answer, but re-reads the walk once per node.
        let mut fai = vec![NO_FIRST_APPEARANCE; max_id + 1];
        for (pos, node_id) in euler.iter().enumerate() {
            let idx: usize = (*node_id).into();
            if fai[idx] == NO_FIRST_APPEARANCE {
                fai[idx] = pos as u32;
            }
        }

        // Depth of every node in one preorder pass over the parent/child
        // structure. A parent is popped before its children, so each depth is
        // already known when it is needed -- cheaper than asking the tree for
        // the depth of each euler entry, which walks to the root every time
        // over a walk twice the size of the tree.
        let root_id = tree.get_root_id();
        let mut node_depth = vec![0usize; max_id + 1];
        let mut stack = vec![(root_id, 0usize)];
        while let Some((node_id, depth)) = stack.pop() {
            node_depth[Into::<usize>::into(node_id)] = depth;
            for child_id in tree.get_node_children_ids(node_id) {
                stack.push((child_id, depth + 1));
            }
        }
        let da: Vec<usize> = euler
            .iter()
            .map(|x| node_depth[Into::<usize>::into(*x)])
            .collect();

        let rmq = FastRmq::from_vec(da.iter().map(|x| *x as u64).collect_vec());

        LcaOracle {
            tree,
            euler,
            fai,
            da,
            rmq,
        }
    }

    /// Constant-time LCA of a slice of nodes, by NodeID.
    ///
    /// # Panics
    ///
    /// Panics if `node_id_vec` is empty, or if any id in it never appeared in
    /// the euler walk (see [`get_fa_index`](Self::get_fa_index)).
    pub fn get_lca_id(&self, node_id_vec: &[TreeNodeID<Tree>]) -> TreeNodeID<Tree> {
        if node_id_vec.len() == 1 {
            return node_id_vec[0];
        }
        // First-appearance lookups are O(1), so the euler-tour range minimum is
        // the cheap route. One pass for both bounds: `get_fa_index` is not free
        // enough to call twice per node.
        let (min_pos, max_pos) = node_id_vec
            .iter()
            .map(|x| self.get_fa_index(*x))
            .fold((usize::MAX, usize::MIN), |(lo, hi), pos| {
                (lo.min(pos), hi.max(pos))
            });
        self.euler[self.rmq.range_min(min_pos, max_pos)]
    }

    /// Constant-time LCA of a slice of nodes, by immutable reference.
    pub fn get_lca(&self, node_id_vec: &[TreeNodeID<Tree>]) -> &'t Tree::Node {
        self.tree.get_node(self.get_lca_id(node_id_vec)).unwrap()
    }

    /// First-appearance position of `node_id` in the euler tour.
    ///
    /// # Panics
    ///
    /// Panics if `node_id` never appeared in the euler walk, i.e. it is not a
    /// node of the tree this oracle was built from.
    pub fn get_fa_index(&self, node_id: TreeNodeID<Tree>) -> usize {
        let pos = self.fai[Into::<usize>::into(node_id)];
        assert_ne!(
            pos, NO_FIRST_APPEARANCE,
            "node does not appear in the euler walk"
        );
        pos as usize
    }

    /// Depth of `node_id`, read from the depth array.
    pub fn get_node_depth(&self, node_id: TreeNodeID<Tree>) -> usize {
        self.da[self.get_fa_index(node_id)]
    }

    /// NodeID at position `pos` of the euler tour.
    pub fn get_euler_pos(&self, pos: usize) -> TreeNodeID<Tree> {
        self.euler[pos]
    }

    /// Builds an oracle for a subtree extracted from this oracle's tree, reusing
    /// this tour instead of walking the subtree afresh.
    ///
    /// `child` must be the pure subtree rooted at `child.get_root_id()` -- every
    /// descendant kept, node ids and topology preserved, exactly as produced by
    /// [`Subtree::subtree`](crate::tree::ops::Subtree::subtree) or
    /// [`SPR::prune`](crate::tree::ops::SPR::prune). The returned oracle borrows
    /// `child`, not the parent, so the parent tree and oracle may be dropped.
    ///
    /// A rooted subtree occupies one contiguous stretch of the euler tour, so
    /// the euler tour and depth array are sliced out directly; only the RMQ is
    /// rebuilt. Depths are re-based to the subtree root, matching a freshly
    /// built oracle. Passing a `child` that is a subtree of this tree but not a
    /// *pure* one (some descendants dropped, ids or topology altered) yields an
    /// oracle with unspecified (but memory-safe) answers.
    ///
    /// # Panics
    ///
    /// Panics if `child`'s root id does not appear in this oracle's euler tour,
    /// i.e. `child` was not extracted from the tree this oracle was built from.
    pub fn restrict_to_subtree<'c>(&self, child: &'c Tree) -> LcaOracle<'c, Tree> {
        let subtree_root = child.get_root_id();
        let start = self.get_fa_index(subtree_root);
        let root_depth = self.da[start];

        // The subtree is the maximal contiguous euler segment from the root's
        // first appearance whose depths never drop below the root's: stepping
        // out of the subtree means returning to the root's parent, one shallower.
        let mut end = start;
        while end + 1 < self.euler.len() && self.da[end + 1] >= root_depth {
            end += 1;
        }

        let euler: Vec<TreeNodeID<Tree>> = self.euler[start..=end].to_vec();
        // Re-based so the subtree root sits at depth 0.
        let da: Vec<usize> = self.da[start..=end]
            .iter()
            .map(|d| d - root_depth)
            .collect();

        let max_id: usize = euler
            .iter()
            .copied()
            .map(Into::<usize>::into)
            .max()
            .expect("a subtree always contains at least its root");
        let mut fai = vec![NO_FIRST_APPEARANCE; max_id + 1];
        for (pos, node_id) in euler.iter().enumerate() {
            let idx: usize = (*node_id).into();
            if fai[idx] == NO_FIRST_APPEARANCE {
                fai[idx] = pos as u32;
            }
        }

        let rmq = FastRmq::from_vec(da.iter().map(|x| *x as u64).collect_vec());

        LcaOracle {
            tree: child,
            euler,
            fai,
            da,
            rmq,
        }
    }

    /// The full euler tour.
    pub fn euler_slice(&self) -> &[TreeNodeID<Tree>] {
        &self.euler
    }

    /// The full depth array.
    pub fn depth_array(&self) -> &[usize] {
        &self.da
    }

    /// Returns the number of bytes this index has allocated on the heap.
    ///
    /// The RMQ figure comes from the backing implementation's own accounting;
    /// the rest is capacity times element size.
    pub fn heap_size(&self) -> usize {
        self.euler.capacity() * std::mem::size_of::<TreeNodeID<Tree>>()
            + self.fai.capacity() * std::mem::size_of::<u32>()
            + self.da.capacity() * std::mem::size_of::<usize>()
            + self.rmq.heap_size()
    }
}
