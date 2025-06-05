# Phylo


Phylo is a fast, extensible, general-purpose, and WebAssembly-capable library for phylogenetic analysis and inference written in Rust. Phylo-rs leverages a combination of memory safety, speed, and native WebAssembly support offered by Rust to provide a robust set of memory-efficient data structures with basic algorithms ranging from tree manipulation with SPR to computing tree statistics such as phylogenetic diversity.


# A note on implementation


Implementation of tree-like structures in Rust can be complex and time-intensive. Additionally implementing tree traversals and operations on tree structures (recursive or otherwise) can be a substantial task. This crate aims to implement most such methods as easily derivable traits, so you don't have to implement them from scratch where they are not needed.

**We also provide a struct, so you don't have to implement one...**  


# Using `phylo`

Most of the functionality is implemented in [`crate::tree::simple_rtree`]. The

[`crate::tree::ops`] module deals with phylogenetic analysis that requires tree mutations such as SPR, NNI, etc.

[`crate::tree::simulation`] module is used to simulate random trees

[`crate::tree::io`] module is used to read trees from various encodings

[`crate::tree::distances`] module is used to compute various types of distance between nodes in a tree and between trees

[`crate::iter`] is a helper module to provide tree traversals and iterations.


## Building trees

The simplest way to build a tree is to create an empty tree, add a root node, and then add children to the various added nodes:


```rust

use phylo::prelude::*;

let mut tree = PhyloTree::new(1);

let new_node = PhyloNode::new(2);
tree.add_child(tree.get_root_id(), new_node);
let new_node = PhyloNode::new(3);
tree.add_child(tree.get_root_id(), new_node);
let new_node: PhyloNode = PhyloNode::new(4);
tree.add_child(2, new_node);
let new_node: PhyloNode = PhyloNode::new(5);
tree.add_child(2, new_node);
```


## Reading and writing trees

This library can build trees strings (or files) encoded in the [newick](https://en.wikipedia.org/wiki/Newick_format) format:

```rust

use phylo::prelude::*;

let input_str = String::from("((A:0.1,B:0.2),C:0.6);");
let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();
```


## Traversing trees

Several traversals are implemented to visit nodes in a particular order: pre- and post-order. A traversal returns an [`Iterator`] of nodes or NodeIDs in the order in which they are to be visited.

```rust

use phylo::prelude::*;

let input_str = String::from("((A:0.1,B:0.2),C:0.6);");
let tree = PhyloTree::from_newick(input_str.as_bytes()).unwrap();

let dfs_traversal = tree.dfs(tree.get_root_id()).into_iter();
let bfs_traversal = tree.bfs_ids(tree.get_root_id());
let postfix_traversal = tree.postord_ids(tree.get_root_id());

```



## Comparing trees

Several metrics taking into account topology and branch lengths are implemented to compare trees with each other:

```rust

use phylo::prelude::*;

fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
    tree.depth(node_id) as f32
}

let newick_1 = "((A:0.1,B:0.2):0.6,(C:0.3,D:0.4):0.5);";
let newick_2 = "((D:0.3,C:0.4):0.5,(B:0.2,A:0.1):0.6);";

let mut tree_1 = PhyloTree::from_newick(newick_1.as_bytes()).unwrap();
let mut tree_2 = PhyloTree::from_newick(newick_2.as_bytes()).unwrap();

tree_1.precompute_constant_time_lca();
tree_2.precompute_constant_time_lca();

tree_1.set_zeta(depth);
tree_2.set_zeta(depth);


let ca = tree_1.ca(&tree_2);
let cophen = tree_1.cophen_dist(&tree_2, 2);

```


# Examples

The following snippets are code examples of some phylogenetic analyses. You can find these in the `examples` directory of the repository. To visualize the outputs of the examples, please first install the requirements using the following command:


```bash

pip install -r examples/visualization/requirements.txt

```



## Quantifying Phylogenetic Diversity

Quantifying the Phylogenetic Diversity of a set of trees using the Faith Index:

```Rust

#[cfg(feature = "non_crypto_hash")]

use fxhash::FxHashMap as HashMap;

#[cfg(not(feature = "non_crypto_hash"))]

use std::collections::HashMap;


use itertools::Itertools;

use std::fs::{File, read_to_string};

use phylo::prelude::*;

use std::io::Write;


fn main() {

    let paths: HashMap<_, _> = std::fs::read_dir("examples/phylogenetic-diversity/trees")

    .unwrap()

    .map(|x| (x.as_ref().unwrap().file_name().into_string().unwrap(), std::fs::read_dir(x.unwrap().path()).unwrap()

        .map(|f| (f.as_ref().unwrap().file_name().into_string().unwrap().split("-").map(|x| x.to_string()).collect_vec()[0].clone(), PhyloTree::from_newick(read_to_string(f.unwrap().path()).unwrap().as_bytes()).unwrap()))

        .collect::<HashMap<_,_>>()))

    .collect();

    

    for (clade, trees) in paths.iter(){

    println!("Clade: {}", clade);

    let mut pds = vec![];

    for year in 2015..2023{

        let tree = trees.get(&year.to_string());

        match tree{

            Some(t) => {

                println!("{}: {}", year, t.get_nodes().map(|n| n.get_weight().unwrap_or(0.0)).sum::<f32>()); 

                pds.push(t.get_nodes().map(|n| n.get_weight().unwrap_or(0.0)).sum::<f32>());

            },

            _ => {println!("{}: {}", year, 0.0); pds.push(0.0);},

        };

    }

    }

}

```


The variations in Phylogenetic Diversity can be visualized using the python script ```examples/visualization/pd.py```.


## Visualizing Phylogenetic Tree Space

Here, we compare all pairwise RF distances of a set of trees:

```Rust

#[cfg(feature = "non_crypto_hash")]

use fxhash::FxHashMap as HashMap;

#[cfg(not(feature = "non_crypto_hash"))]

use std::collections::HashMap;


use itertools::Itertools;

use std::fs::{File, read_to_string};

use phylo::prelude::*;

use std::io::Write;

use indicatif::{ProgressIterator, ProgressBar, ProgressStyle};


fn main() {

    let trees = (1..11).progress().map(|x| read_to_string(format!("examples/pairwise-distances/r{x}-preprocessed.trees"))

            .unwrap()

            .lines()

            .enumerate()

            .map(|(y,z)| (x,y,PhyloTree::from_newick(z.as_bytes()).unwrap()))

            .collect_vec()

        )

        .flatten()

        .collect_vec();

    

    

    let bar = ProgressBar::new((trees.len()*(trees.len()-1)/2) as u64);

    bar.set_style(ProgressStyle::with_template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg} [eta: {eta}]")

        .unwrap()

        .progress_chars("##-"));

    

    trees.iter().combinations(2).map(|v| (v[0], v[1])).for_each(|(x,y)| {

        let out = format!("{}-{}-{}-{}-{}\n", x.0, y.0, x.1, y.1, x.2.ca(&y.2));

        println!("{}", out);

        bar.inc(1);

    });

    bar.finish();

}

```


The tree space can be visualized using the python script ```examples/visualization/tree-space.py```.




To run the code examples, run the following:

```bash

cargo run --example phylogenetic-diversity


cargo run --example pairwise-distances

```
