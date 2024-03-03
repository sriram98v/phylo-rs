mod dfs_tree;

pub use dfs_tree::{DfsTree, DfsTreeIndex, DfsTreeSiblingsIter};
use std::{collections::VecDeque, ops::ControlFlow};

/// A rooted tree data structure.
///
/// `Data` should be the data that is accessed often. `Meta` should be the data
/// that is stored for every node and is accessed less often.
pub trait RootedTree: Sized {
    type Index: Copy;
    type Data;
    type Meta;

    /// Gets the length of the tree.
    fn len(&self) -> usize;

    /// Gets the root node of the tree.
    fn root(&self) -> Self::Index;

    /// Gets a reference to the data of a node.
    fn get(&self, index: Self::Index) -> &Self::Data;

    /// Gets a mutable reference to the data of a node.
    fn get_mut(&mut self, index: Self::Index) -> &mut Self::Data;

    /// Gets a reference to the meta of a node.
    fn get_meta(&self, index: Self::Index) -> &Self::Meta;

    /// Gets a mutable reference to the meta of a node.
    fn get_meta_mut(&mut self, index: Self::Index) -> &mut Self::Meta;

    /// Returns an iterator over the siblings of a node in any order.
    fn siblings(&self, index: Self::Index) -> impl Iterator<Item = Self::Index>;

    /// Returns the parent of a node.
    fn parent(&self, index: Self::Index) -> Option<Self::Index>;

    /// Returns an iterator over the direct children of a node in any order.
    fn children(&self, index: Self::Index) -> impl Iterator<Item = Self::Index>;

    /// Returns an iterator over the descendants of a node in any order.
    fn descendants(&self, index: Self::Index) -> impl Iterator<Item = Self::Index>;

    /// Returns an iterator over the nodes of the tree in any order.
    fn iter(&self) -> impl Iterator<Item = Self::Index>;

    /// Depth-first search iterator.
    fn dfs(&self, start_node: Self::Index) -> impl Iterator<Item = Self::Index> {
        RootedTreeDfs {
            tree: self,
            stack: vec![start_node],
        }
    }

    /// Depth-first search iteration with mutable access to the tree.
    ///
    /// ## Panics
    ///
    /// - the closure modifies the tree in a way that violates the tree's
    ///   invariants
    fn dfs_mut<T>(
        &mut self,
        start_node: Self::Index,
        mut f: impl FnMut(&mut Self, Self::Index) -> ControlFlow<T, ()>,
    ) -> ControlFlow<T, ()> {
        let mut stack = vec![start_node];
        while let Some(index) = stack.pop() {
            f(self, index)?;
            stack.extend(self.children(index));
        }
        ControlFlow::Continue(())
    }

    /// Breadth-first search iterator.
    fn bfs(&self, start_node: Self::Index) -> impl Iterator<Item = Self::Index> {
        RootedTreeBfs {
            tree: self,
            queue: VecDeque::from([start_node]),
        }
    }

    /// Breadth-first search iteration with mutable access to the tree.
    ///
    /// ## Panics
    ///
    /// - the closure modifies the tree in a way that violates the tree's
    ///   invariants
    fn bfs_mut<T>(
        &mut self,
        start_node: Self::Index,
        mut f: impl FnMut(&mut Self, Self::Index) -> ControlFlow<T, ()>,
    ) -> ControlFlow<T, ()> {
        let mut queue = VecDeque::from([start_node]);
        while let Some(index) = queue.pop_front() {
            f(self, index)?;
            queue.extend(self.children(index));
        }
        ControlFlow::Continue(())
    }
}

struct RootedTreeDfs<'a, T: RootedTree> {
    tree: &'a T,
    stack: Vec<T::Index>,
}

impl<T: RootedTree> Iterator for RootedTreeDfs<'_, T> {
    type Item = T::Index;

    fn next(&mut self) -> Option<Self::Item> {
        self.stack.pop().map(|index| {
            self.stack.extend(self.tree.children(index));
            index
        })
    }
}

struct RootedTreeBfs<'a, T: RootedTree> {
    tree: &'a T,
    queue: VecDeque<T::Index>,
}

impl<T: RootedTree> Iterator for RootedTreeBfs<'_, T> {
    type Item = T::Index;

    fn next(&mut self) -> Option<Self::Item> {
        self.queue.pop_front().map(|index| {
            self.queue.extend(self.tree.children(index));
            index
        })
    }
}

pub trait RootedTreeMut: RootedTree {
    /// Frees all unatached nodes.
    fn clean(&mut self);

    /// Swaps the data and meta of two nodes.
    fn swap(&mut self, a: Self::Index, b: Self::Index);

    /// Adds a new node to the tree.
    fn add(&mut self, parent: Self::Index, data: Self::Data, meta: Self::Meta) -> Self::Index;

    /// Removes a node from the tree.
    fn remove(&mut self, index: Self::Index) -> Option<(Self::Data, Self::Meta)>;
}
