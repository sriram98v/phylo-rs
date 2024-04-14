use std::collections::HashMap;

use rayon::iter::{IntoParallelIterator, ParallelBridge, ParallelIterator};

use super::{ParallelRootedTree, RootedTree, RootedTreeMut};

type SimpleTreeIndex = usize;

#[derive(Debug, Clone)]
struct SimpleTreeNode<Data, Meta> {
    data: Data,
    meta: Meta,
    parent: Option<SimpleTreeIndex>,
    children: Vec<SimpleTreeIndex>,
}

#[derive(Debug, Clone)]
pub struct SimpleTree<Data, Meta> {
    nodes: HashMap<SimpleTreeIndex, SimpleTreeNode<Data, Meta>>,
}

impl<Data: Default, Meta: Default> Default for SimpleTree<Data, Meta> {
    fn default() -> Self {
        let mut nodes = HashMap::new();
        nodes.insert(
            0,
            SimpleTreeNode {
                data: Data::default(),
                meta: Meta::default(),
                parent: None,
                children: vec![],
            },
        );
        Self { nodes }
    }
}

impl<Data, Meta> RootedTree for SimpleTree<Data, Meta> {
    type Index = SimpleTreeIndex;
    type Data = Data;
    type Meta = Meta;

    fn len(&self) -> usize {
        self.iter().count()
    }

    fn root(&self) -> Self::Index {
        0
    }

    fn get(&self, index: Self::Index) -> &Self::Data {
        &self.nodes.get(&index).unwrap().data
    }

    fn get_mut(&mut self, index: Self::Index) -> &mut Self::Data {
        &mut self.nodes.get_mut(&index).unwrap().data
    }

    fn get_meta(&self, index: Self::Index) -> &Self::Meta {
        &self.nodes.get(&index).unwrap().meta
    }

    fn get_meta_mut(&mut self, index: Self::Index) -> &mut Self::Meta {
        &mut self.nodes.get_mut(&index).unwrap().meta
    }

    fn siblings(&self, index: Self::Index) -> impl Iterator<Item = Self::Index> {
        self.nodes
            .get(&index)
            .unwrap()
            .parent
            .map_or_else(
                || [].iter(),
                |parent| self.nodes.get(&parent).unwrap().children.iter(),
            )
            .filter(move |&&node_index| index != node_index)
            .copied()
    }

    fn parent(&self, index: Self::Index) -> Option<Self::Index> {
        self.nodes.get(&index).unwrap().parent
    }

    fn children(&self, index: Self::Index) -> impl Iterator<Item = Self::Index> {
        self.nodes.get(&index).unwrap().children.iter().copied()
    }
}

impl<Data, Meta> RootedTreeMut for SimpleTree<Data, Meta> {
    fn clean(&mut self) {
        let mut marked = vec![false; self.nodes.len()];
        for index in self.iter() {
            marked[index] = true;
        }
        self.nodes.retain(|index, _| marked[*index]);
    }

    fn add_child(
        &mut self,
        parent: Self::Index,
        data: Self::Data,
        meta: Self::Meta,
    ) -> Self::Index {
        let index = self.nodes.len();
        self.nodes.insert(
            index,
            SimpleTreeNode {
                data,
                meta,
                parent: Some(parent),
                children: vec![],
            },
        );
        self.nodes.get_mut(&parent).unwrap().children.push(index);
        index
    }

    fn remove(&mut self, index: Self::Index) -> (Self::Data, Self::Meta) {
        let node = self.nodes.remove(&index).unwrap();
        if let Some(parent) = node.parent {
            let children = &mut self.nodes.get_mut(&parent).unwrap().children;
            children.swap_remove(children.iter().position(|&child| child == index).unwrap());
        }
        (node.data, node.meta)
    }

    fn reserve(&mut self, additional: usize) {
        self.nodes.reserve(additional)
    }
}

impl<Data: Sync, Meta: Sync> ParallelRootedTree for SimpleTree<Data, Meta> {
    fn par_iter(&self) -> impl ParallelIterator<Item = Self::Index> {
        self.nodes.keys().copied().par_bridge().into_par_iter()
    }
}
