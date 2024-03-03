mod dfs_tree;
mod simple_tree;

pub use dfs_tree::{DfsTree, DfsTreeIndex, DfsTreeSiblingsIter};
pub use simple_tree::SimpleTree;
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
    fn descendants(&self, index: Self::Index) -> impl Iterator<Item = Self::Index> {
        self.dfs(index).skip(1)
    }

    /// Returns an iterator over the nodes of the tree in any order.
    fn iter(&self) -> impl Iterator<Item = Self::Index> {
        self.dfs(self.root())
    }

    /// Depth-first search iterator.
    fn dfs(&self, start_node: Self::Index) -> impl Iterator<Item = Self::Index> {
        let mut stack = vec![start_node];
        std::iter::from_fn(move || {
            stack.pop().inspect(|&index| {
                stack.extend(self.children(index));
            })
        })
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
        let mut queue = VecDeque::from([start_node]);
        std::iter::from_fn(move || {
            queue.pop_front().inspect(|&index| {
                queue.extend(self.children(index));
            })
        })
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

pub trait RootedTreeMut: RootedTree {
    /// Frees all unattached nodes.
    fn clean(&mut self);

    /// Adds a new node to the parent in the tree.
    fn add_child(&mut self, parent: Self::Index, data: Self::Data, meta: Self::Meta)
        -> Self::Index;

    /// Removes a node from the tree.
    fn remove(&mut self, index: Self::Index) -> (Self::Data, Self::Meta);
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn simple() {
        let mut tree = SimpleTree::default();

        let a = tree.add_child(tree.root(), 1, ());
        let b = tree.add_child(a, 2, ());
        let c = tree.add_child(tree.root(), 3, ());
        let d = tree.add_child(c, 4, ());
        let e = tree.add_child(c, 5, ());
        let f = tree.add_child(e, 6, ());
        let g = tree.add_child(tree.root(), 7, ());
        let h = tree.add_child(g, 8, ());
        let i = tree.add_child(h, 9, ());
        let j = tree.add_child(i, 10, ());

        assert_eq!(tree.len(), 11);
        tree.clean();
        assert_eq!(tree.len(), 11);
        tree.remove(j);
        assert_eq!(tree.len(), 10);
        tree.remove(c);
        assert_eq!(tree.len(), 6);
    }
}
