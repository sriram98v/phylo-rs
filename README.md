<div align="center">

# 🌳 phylo

**A fast, extensible, WebAssembly-ready phylogenetics library for Rust.**

[![Crates.io](https://img.shields.io/crates/v/phylo.svg)](https://crates.io/crates/phylo)
[![Documentation](https://img.shields.io/docsrs/phylo)](https://docs.rs/phylo)
[![CI](https://github.com/sriram98v/phylo-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sriram98v/phylo-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/crates/l/phylo.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.87-blue.svg)](https://blog.rust-lang.org/)
[![Downloads](https://img.shields.io/crates/d/phylo.svg)](https://crates.io/crates/phylo)

</div>

---

`phylo` provides memory-efficient data structures and algorithms for phylogenetic
analysis and inference — from tree manipulation (SPR, NNI, rerooting) to tree
statistics (phylogenetic diversity, RF distance, cophenetic distance). It leans
on Rust's memory safety, speed, and native WebAssembly support to stay both fast
and portable.

Tree traversals and operations are exposed as **derivable traits**, so you get
DFS/BFS/pre-/post-order, Euler tours, LCA queries, and distance metrics for free
on your own types — and a ready-made `SimpleRootedTree` when you don't want to
implement one.

## Highlights

- **Trait-first design** — compose narrow traits (`RootedTree`, `RootedMetaTree`,
  `EulerWalk`, `DFS`, `Clusters`, …) onto any type, or use the batteries-included
  `PhyloTree`.
- **Arena-allocated trees** — cache-friendly `Vec`-backed storage with `usize`
  node IDs.
- **Constant-time LCA** — an [`LcaOracle`](https://docs.rs/phylo) borrows the tree
  immutably and answers LCA queries in O(1) via an Euler tour + RMQ.
- **Tree comparison** — Robinson-Foulds, weighted RF, cluster affinity, and
  cophenetic distance, with distance-matrix builders.
- **I/O** — Newick and Nexus parsing and serialization.
- **Simulation** — random trees (Yule, uniform).
- **Optional parallelism** — opt into `rayon`-backed computation with the
  `parallel` feature.

## Installation

```sh
cargo add phylo
```

Or add it to `Cargo.toml`:

```toml
[dependencies]
phylo = "5"
```

### Feature flags

| Feature | Default | Description |
| --- | :---: | --- |
| `simple_rooted_tree` | ✅ | The concrete `SimpleRootedTree` / `PhyloTree` implementation. |
| `non_crypto_hash` | ✅ | Use `fxhash` maps/sets instead of `std` for speed. |
| `parallel` | | `rayon`-based parallel computation for the heavy metrics. |
| `serde` | | `Serialize`/`Deserialize` for trees. |

## Quick start

Everything you need is in the prelude:

```rust
use phylo::prelude::*;
```

### Build a tree

Create an empty tree, then attach children to node IDs:

```rust
use phylo::prelude::*;

let mut tree = PhyloTree::new(1);

tree.add_child(tree.get_root_id(), PhyloNode::new(2));
tree.add_child(tree.get_root_id(), PhyloNode::new(3));
tree.add_child(2, PhyloNode::new(4));
tree.add_child(2, PhyloNode::new(5));
```

### Read and write Newick

```rust
use phylo::prelude::*;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2),C:0.6);".as_bytes()).unwrap();
let newick = tree.to_newick();
```

### Traverse

Traversals return an `Iterator` of nodes or node IDs in visiting order:

```rust
use phylo::prelude::*;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2),C:0.6);".as_bytes()).unwrap();

let dfs = tree.dfs(tree.get_root_id());
let bfs = tree.bfs_ids(tree.get_root_id());
let postorder = tree.postord_ids(tree.get_root_id());
```

### Constant-time LCA

Build an `LcaOracle` with `tree.lca()`; it borrows the tree immutably (so
staleness is a compile error, not a runtime bug) and answers queries in O(1):

```rust
use phylo::prelude::*;

let tree = PhyloTree::from_newick("((A,B),(C,D));".as_bytes()).unwrap();

let a = tree.get_taxa_node_id(&"A".to_string()).unwrap();
let b = tree.get_taxa_node_id(&"B".to_string()).unwrap();

let lca = tree.lca();
let ancestor = lca.get_lca_id(&[a, b]);
```

### Compare trees

Metrics account for both topology and branch lengths:

```rust
use phylo::prelude::*;

fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
    tree.depth(node_id) as f32
}

let mut tree_1 = PhyloTree::from_newick("((A:0.1,B:0.2):0.6,(C:0.3,D:0.4):0.5);".as_bytes()).unwrap();
let mut tree_2 = PhyloTree::from_newick("((D:0.3,C:0.4):0.5,(B:0.2,A:0.1):0.6);".as_bytes()).unwrap();

let _ = tree_1.set_zeta(depth);
let _ = tree_2.set_zeta(depth);

let cluster_affinity = tree_1.ca(&tree_2);
let cophenetic = tree_1.cophen_dist(&tree_2, 2);
```

## Module map

| Module | What it does |
| --- | --- |
| [`tree::simple_rtree`](https://docs.rs/phylo/latest/phylo/tree/simple_rtree/) | Core tree traits and `SimpleRootedTree`. |
| [`tree::ops`](https://docs.rs/phylo/latest/phylo/tree/ops/) | Mutating operations: SPR, NNI, reroot, contraction, subtree extraction. |
| [`tree::distances`](https://docs.rs/phylo/latest/phylo/tree/distances/) | RF, weighted RF, cluster affinity, cophenetic distance, distance matrices. |
| [`tree::io`](https://docs.rs/phylo/latest/phylo/tree/io/) | Newick and Nexus reading/writing. |
| [`tree::simulation`](https://docs.rs/phylo/latest/phylo/tree/simulation/) | Random tree generation. |
| [`iter`](https://docs.rs/phylo/latest/phylo/iter/) | Traversals, Euler walks, and the LCA oracle. |

## Examples

Runnable analyses live in the [`examples/`](examples) directory. To visualize
their output, install the Python requirements first:

```sh
pip install -r examples/visualization/requirements.txt
```

**Quantifying phylogenetic diversity** — the Faith index across a set of trees.
Run it, then plot with `examples/visualization/pd.py`:

```sh
cargo run --example phylogenetic-diversity
```

**Visualizing tree space** — all pairwise distances across a set of trees. Run
it, then plot with `examples/visualization/tree-space.py`:

```sh
cargo run --example pairwise-distances
```

## WebAssembly

`phylo` builds for `wasm32` targets out of the box, making it suitable for
in-browser phylogenetics. Use your usual wasm toolchain — e.g. `wasm-pack`, or
`cargo build --target wasm32-unknown-unknown`.

## Citation

If you use `phylo` in your work, please cite
[this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12309125/):

```bibtex
@article{vijendran2025phylo,
  title={Phylo-rs: an extensible phylogenetic analysis library in rust},
  author={Vijendran, Sriram and Anderson, Tavis and Markin, Alexey and Eulenstein, Oliver},
  journal={BMC bioinformatics},
  volume={26},
  pages={197},
  year={2025}
}
```

## License

Licensed under the [MIT License](LICENSE).
