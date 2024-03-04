use std::{iter::Skip, num::NonZeroU32, ops::ControlFlow};

use super::RootedTree;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DfsTreeIndex(NonZeroU32);

impl DfsTreeIndex {
    pub fn from_usize(index: usize) -> Option<DfsTreeIndex> {
        Some(DfsTreeIndex(NonZeroU32::new(index as u32)?))
    }
}

const ROOT: DfsTreeIndex = DfsTreeIndex(unsafe { NonZeroU32::new_unchecked(1) });

type DfsIndexRange = std::ops::Range<DfsTreeIndex>;

/// An immutable tree optimized for depth-first search iteration.
///
/// Internally, the tree is represented in a flat array in DFS order. This
/// allows for fast lookups and certain kinds of iteration, but makes adding and
/// removing nodes slow and requires reordering the array after modifications.
/// So, this tree does not support adding or removing nodes after it has been
/// created.
///
/// Forcing the tree to be immutable also allows for more features such as more
/// specific iterators that allow safe mutable access to the tree data during
/// iteration.
#[derive(Debug, Clone)]
pub struct DfsTree<Data: Default, Meta: Default> {
    siblings: Vec<DfsTreeSiblings>,
    datas: Vec<Data>,
    metas: Vec<Meta>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DfsIterEntryRelation {
    SiblingOfLastEntry,
    ChildOfLastEntry,
    SiblingOfNthParent { depth_offset: usize },
}

#[derive(Debug, Clone)]
pub struct DfsIterEntry<Data, Meta> {
    pub relation: DfsIterEntryRelation,
    pub data: Data,
    pub meta: Meta,
}

impl<Data: Default, Meta: Default, T: Iterator<Item = DfsIterEntry<Data, Meta>>> From<T>
    for DfsTree<Data, Meta>
{
    fn from(iter: T) -> Self {
        let mut tree = Self::default();

        let (lower, upper) = iter.size_hint();
        let reserved = upper.unwrap_or(lower);
        tree.datas.reserve(reserved);
        tree.siblings.reserve(reserved);
        tree.metas.reserve(reserved);

        let mut last_index = tree.root().0.get() as usize;
        let mut parent_stack = vec![];

        for entry in iter {
            let new_index = tree.datas.len();

            tree.datas.push(entry.data);
            tree.metas.push(entry.meta);

            let last_siblings = tree.siblings.get_mut(last_index).unwrap();
            let siblings = match entry.relation {
                DfsIterEntryRelation::SiblingOfLastEntry => {
                    last_siblings.right = Some(DfsTreeIndex::from_usize(new_index).unwrap());
                    DfsTreeSiblings {
                        left: Some(DfsTreeIndex::from_usize(last_index).unwrap()),
                        right: None,
                    }
                }
                DfsIterEntryRelation::ChildOfLastEntry => {
                    parent_stack.push(last_index);
                    DfsTreeSiblings::default()
                }
                DfsIterEntryRelation::SiblingOfNthParent { depth_offset } => {
                    parent_stack.truncate(parent_stack.len() - depth_offset);
                    let parent_index = parent_stack.pop().unwrap();
                    let parent_siblings = tree.siblings.get_mut(parent_index).unwrap();
                    parent_siblings.right = Some(DfsTreeIndex::from_usize(new_index).unwrap());
                    DfsTreeSiblings {
                        left: Some(DfsTreeIndex::from_usize(parent_index).unwrap()),
                        right: None,
                    }
                }
            };
            tree.siblings.push(siblings);

            last_index = new_index;
        }

        tree
    }
}

impl<Data: Default, Meta: Default> Default for DfsTree<Data, Meta> {
    fn default() -> Self {
        // TODO: Consider using MaybeUninit<T> to avoid the Default constraint at the
        // cost of unsafe code.
        Self {
            // [dummy, root]
            datas: vec![Data::default(), Data::default()],
            siblings: vec![DfsTreeSiblings::default(), DfsTreeSiblings::default()],
            metas: vec![Meta::default(), Meta::default()],
        }
    }
}

impl<Data: Default, Meta: Default> RootedTree for DfsTree<Data, Meta> {
    type Data = Data;
    type Meta = Meta;
    type Index = DfsTreeIndex;

    fn len(&self) -> usize {
        self.datas.len() - 1
    }

    fn root(&self) -> DfsTreeIndex {
        ROOT
    }

    fn get(&self, index: DfsTreeIndex) -> &Data {
        self.datas.get(index.0.get() as usize).unwrap()
    }

    fn get_mut(&mut self, index: DfsTreeIndex) -> &mut Data {
        self.datas.get_mut(index.0.get() as usize).unwrap()
    }

    fn get_meta(&self, index: DfsTreeIndex) -> &Meta {
        self.metas.get(index.0.get() as usize).unwrap()
    }

    fn get_meta_mut(&mut self, index: DfsTreeIndex) -> &mut Meta {
        self.metas.get_mut(index.0.get() as usize).unwrap()
    }

    fn parent(&self, index: DfsTreeIndex) -> Option<DfsTreeIndex> {
        let mut left = index;
        while let Some(new_left) = self.siblings.get(left.0.get() as usize).unwrap().left {
            left = new_left;
        }
        (left != ROOT).then(|| DfsTreeIndex((left.0.get() as u32 - 1).try_into().unwrap()))
    }

    fn siblings(&self, index: DfsTreeIndex) -> impl Iterator<Item = Self::Index> {
        DfsTreeSiblingsIter::<'_, true, Data, Meta>::new(
            self,
            self.siblings.get(index.0.get() as usize).unwrap().clone(),
        )
    }

    fn children(&self, index: Self::Index) -> impl Iterator<Item = Self::Index> {
        let next_index = DfsTreeIndex(index.0.get().checked_add(1).unwrap().try_into().unwrap());

        if let Some(next_siblings) = self.siblings.get(next_index.0.get() as usize) {
            if next_siblings.left.is_none() {
                return DfsTreeSiblingsIter::<'_, false, Data, Meta>::new(
                    self,
                    DfsTreeSiblings {
                        left: Some(next_index),
                        right: next_siblings.right,
                    },
                );
            }
        }

        DfsTreeSiblingsIter::new(self, DfsTreeSiblings::default())
    }

    #[allow(refining_impl_trait)]
    fn descendants(&self, index: Self::Index) -> Skip<DfsTreeIter> {
        self.dfs(index).skip(1)
    }

    #[allow(refining_impl_trait)]
    fn iter(&self) -> DfsTreeIter {
        self.dfs(self.root())
    }

    #[allow(refining_impl_trait)]
    fn dfs(&self, start_node: DfsTreeIndex) -> DfsTreeIter {
        let end = self
            .siblings
            .get(start_node.0.get() as usize)
            .unwrap()
            .right
            .unwrap_or_else(|| DfsTreeIndex((self.datas.len() as u32).try_into().unwrap()));
        DfsTreeIter {
            range: start_node..end,
        }
    }

    fn dfs_mut<T>(
        &mut self,
        start_node: DfsTreeIndex,
        mut f: impl FnMut(&mut Self, DfsTreeIndex) -> ControlFlow<T, ()>,
    ) -> ControlFlow<T, ()> {
        for index in self.dfs(start_node) {
            f(self, index)?;
        }
        ControlFlow::Continue(())
    }
}

pub struct DfsTreeIter {
    range: DfsIndexRange,
}

impl Iterator for DfsTreeIter {
    type Item = DfsTreeIndex;

    fn next(&mut self) -> Option<Self::Item> {
        let (DfsTreeIndex(start), DfsTreeIndex(end)) = (self.range.start, self.range.end);
        (start < end).then(|| {
            self.range.start = DfsTreeIndex(start.checked_add(1).unwrap());
            DfsTreeIndex(start)
        })
    }
}

impl ExactSizeIterator for DfsTreeIter {
    fn len(&self) -> usize {
        (self.range.end.0.get() - self.range.start.0.get()) as usize
    }
}

#[derive(Debug, Default, Clone, Copy)]
struct DfsTreeSiblings {
    left: Option<DfsTreeIndex>,
    right: Option<DfsTreeIndex>,
}

pub struct DfsTreeSiblingsIter<'a, const LEFT_FIRST: bool, Data: Default, Meta: Default> {
    tree: &'a DfsTree<Data, Meta>,
    siblings: DfsTreeSiblings,
}

impl<'a, const LEFT_FIRST: bool, Data: Default, Meta: Default>
    DfsTreeSiblingsIter<'a, LEFT_FIRST, Data, Meta>
{
    fn new(tree: &'a DfsTree<Data, Meta>, siblings: DfsTreeSiblings) -> Self {
        Self { tree, siblings }
    }
}

impl<Data: Default, Meta: Default> Iterator for DfsTreeSiblingsIter<'_, true, Data, Meta> {
    type Item = DfsTreeIndex;
    fn next(&mut self) -> Option<Self::Item> {
        self.siblings
            .left
            .inspect(|DfsTreeIndex(left)| {
                self.siblings.left = self.tree.siblings.get(left.get() as usize).unwrap().left;
            })
            .or_else(|| {
                let right = self.siblings.right?;
                self.siblings.right = self
                    .tree
                    .siblings
                    .get(right.0.get() as usize)
                    .unwrap()
                    .right;
                Some(right)
            })
    }
}

impl<Data: Default, Meta: Default> Iterator for DfsTreeSiblingsIter<'_, false, Data, Meta> {
    type Item = DfsTreeIndex;
    fn next(&mut self) -> Option<Self::Item> {
        self.siblings
            .right
            .inspect(|DfsTreeIndex(right)| {
                self.siblings.right = self.tree.siblings.get(right.get() as usize).unwrap().right;
            })
            .or_else(|| {
                let left = self.siblings.left?;
                self.siblings.left = self.tree.siblings.get(left.0.get() as usize).unwrap().left;
                Some(left)
            })
    }
}
