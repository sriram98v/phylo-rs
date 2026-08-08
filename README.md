<div align="center">

# 🌳 phylo

[![Crates.io](https://img.shields.io/crates/v/phylo.svg)](https://crates.io/crates/phylo)
[![Documentation](https://img.shields.io/docsrs/phylo)](https://docs.rs/phylo)
[![CI](https://github.com/sriram98v/phylo-rs/actions/workflows/ci.yml/badge.svg)](https://github.com/sriram98v/phylo-rs/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/crates/l/phylo.svg)](LICENSE)
[![MSRV](https://img.shields.io/badge/MSRV-1.87-blue.svg)](https://blog.rust-lang.org/)
[![Downloads](https://img.shields.io/crates/d/phylo.svg)](https://crates.io/crates/phylo)

</div>

---

<!-- cargo-rdme start -->

A fast, extensible, WebAssembly-ready phylogenetics library for Rust.

`phylo` provides memory-efficient data structures and algorithms for phylogenetic analysis and inference — from tree manipulation (SPR, NNI, rerooting) to tree statistics (phylogenetic diversity, RF distance, cophenetic distance) to maximum-likelihood modelling (GTR+I+G substitution models, Felsenstein pruning, ancestral reconstruction). It leans on Rust's memory safety, speed, and native WebAssembly support to stay both fast and portable.

Tree traversals and operations are exposed as **derivable traits**, so you get DFS/BFS/pre-/post-order, Euler tours, LCA queries, and distance metrics for free on your own types — and a ready-made [`PhyloTree`](https://docs.rs/phylo/latest/phylo/tree/simple_rooted_tree/type.PhyloTree.html) when you don't want to implement one.

## Highlights

- **Trait-first design** — compose narrow traits (`RootedTree`, `RootedMetaTree`, `EulerWalk`, `DFS`, `Clusters`, …) onto any type, or use the batteries-included [`PhyloTree`](https://docs.rs/phylo/latest/phylo/tree/simple_rooted_tree/type.PhyloTree.html).
- **Arena-allocated trees** — cache-friendly `Vec`-backed storage with `usize` node IDs.
- **Constant-time LCA** — an [`LcaOracle`](https://docs.rs/phylo/latest/phylo/iter/lca/struct.LcaOracle.html) borrows the tree immutably and answers LCA queries in O(1) via an Euler tour + RMQ.
- **Tree comparison** — Robinson-Foulds, weighted RF, cluster affinity, and cophenetic distance, with distance-matrix builders.
- **Maximum-likelihood modelling** — GTR+I+G substitution models (JC69 through GTR), Felsenstein-pruning log-likelihood, and marginal/joint ancestral sequence reconstruction.
- **I/O** — Newick and Nexus parsing and serialization.
- **Simulation** — random trees (Yule, uniform).
- **Optional parallelism** — opt into `rayon`-backed computation with the `parallel` feature.
- **Fallible by default** — operations that a caller can misuse return [`Result`](https://doc.rust-lang.org/stable/core/result/enum.Result.html) with a typed [`error::TreeError`](https://docs.rs/phylo/latest/phylo/error/enum.TreeError.html); the library does not panic on bad input.

## Installation

```sh
cargo add phylo
```

## Feature flags

| Feature | Default | Description |
| --- | :---: | --- |
| `simple_rooted_tree` | Yes | The concrete `SimpleRootedTree` / `PhyloTree` implementation. |
| `non_crypto_hash` | Yes | Use `fxhash` maps/sets instead of `std` for speed. |
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

Traversals return an [`Iterator`](https://doc.rust-lang.org/stable/core/iter/traits/iterator/trait.Iterator.html) of nodes or node IDs in visiting order:

```rust
use phylo::prelude::*;

let tree = PhyloTree::from_newick("((A:0.1,B:0.2),C:0.6);".as_bytes()).unwrap();

let dfs = tree.dfs(tree.get_root_id()).unwrap();
let bfs = tree.bfs_ids(tree.get_root_id()).unwrap();
let postorder = tree.postord_ids(tree.get_root_id()).unwrap();
```

### Constant-time LCA

Build an [`LcaOracle`](https://docs.rs/phylo/latest/phylo/iter/lca/struct.LcaOracle.html) with `tree.lca()`; it borrows the tree immutably (so staleness is a compile error, not a runtime bug) and answers queries in O(1):

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

tree_1.set_zeta(depth).unwrap();
tree_2.set_zeta(depth).unwrap();

let cluster_affinity = tree_1.ca(&tree_2);
let cophenetic = tree_1.cophen_dist(&tree_2, 2).unwrap();
```

### Likelihood and ancestral reconstruction

Score an alignment against a tree under a substitution model, or reconstruct ancestral sequences at the internal nodes. The log-likelihood path ([`TreeLikelihood`](https://docs.rs/phylo/latest/phylo/tree/likelihood/trait.TreeLikelihood.html)) runs Felsenstein's pruning algorithm alone — no reconstruction — while marginal/joint ASR ([`MarginalAsr`](https://docs.rs/phylo/latest/phylo/tree/asr/trait.MarginalAsr.html) / [`JointAsr`](https://docs.rs/phylo/latest/phylo/tree/asr/trait.JointAsr.html)) build on the same pruning core:

```rust
use phylo::prelude::*;

let tree =
    PhyloTree::from_newick("((A:0.1,B:0.2):0.15,(C:0.3,D:0.1):0.05);".as_bytes()).unwrap();

// A nucleotide alignment in FASTA — one sequence per leaf taxon.
let fasta = b">A\nACGTACGT\n>B\nACGTATGT\n>C\nACGAACGT\n>D\nTCGTACGA\n";
let aln = Alignment::from_fasta_bytes(fasta).unwrap();

// HKY85 with gamma-distributed rate heterogeneity (+G, 4 categories).
let model = GtrModel::<Nucleotide>::hky85([0.25, 0.25, 0.25, 0.25], 2.0)
    .unwrap()
    .with_gamma(0.5, 4)
    .unwrap();

// Log-likelihood of the alignment given the tree and model (pruning only).
let log_lik = tree.log_likelihood::<Nucleotide>(&model, &aln).unwrap();
assert!(log_lik.is_finite());

// Marginal ancestral sequence reconstruction fills the internal nodes.
let recon = tree.marginal_asr::<Nucleotide>(&model, &aln, false).unwrap();
let root_sequence = recon.sequence_string(tree.get_root_id());
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
| [`models`](https://docs.rs/phylo/latest/phylo/models/) | GTR+I+G substitution models and their named special cases. |
| [`tree::likelihood`](https://docs.rs/phylo/latest/phylo/tree/likelihood/) | Felsenstein-pruning log-likelihood. |
| [`tree::asr`](https://docs.rs/phylo/latest/phylo/tree/asr/) | Marginal and joint ancestral sequence reconstruction. |
| [`error`](https://docs.rs/phylo/latest/phylo/error/) | [`error::TreeError`](https://docs.rs/phylo/latest/phylo/error/enum.TreeError.html) and the parsing/model error types. |

## Examples

Runnable analyses live in the [`examples/`](https://github.com/sriram98v/phylo-rs/tree/main/examples) directory. To visualize their output, install the Python requirements first:

```sh
pip install -r examples/visualization/requirements.txt
```

**Quantifying phylogenetic diversity** — the Faith index across a set of trees. Run it, then plot with `examples/visualization/pd.py`:

```sh
cargo run --example phylogenetic-diversity
```

**Visualizing tree space** — all pairwise distances across a set of trees. Run it, then plot with `examples/visualization/tree-space.py`:

```sh
cargo run --example pairwise-distances
```

## Benchmarks

Benchmarks use [criterion](https://github.com/bheisler/criterion.rs) and are split into four targets by what they measure:

| Target | Covers |
| --- | --- |
| `lca` | `LcaOracle` construction, O(1) queries, and the uncached walk. |
| `traversal` | Post-order traversal, bipartitions, median node. |
| `construction` | Yule simulation, SPR, tree contraction. |
| `distances` | RF, cluster matching, cluster affinity, cophenetic distance. |

```sh
cargo bench                      # everything
cargo bench --bench lca          # one target
cargo bench -- lca_oracle_query  # one group, by regex
cargo bench --features parallel  # adds the parallel cophenetic group
```

Each sweep reports **throughput** in elements (taxa) per second alongside wall time. That is the number to read: a linear routine holds its throughput flat across the sweep, while a quadratic one loses roughly a factor of four per fourfold jump in taxa. Reading the raw times alone hides which is which.

Some groups come in pairs that are only meaningful read together. `contract_tree` builds a throwaway [`LcaOracle`](https://docs.rs/phylo/latest/phylo/iter/lca/struct.LcaOracle.html) per call, while `contract_tree_with_oracle` reuses one the caller hoisted; the gap between them at a given taxa count is what the O(n) index build costs, and it is the difference between the two that tells you whether a change touched the build or the contraction.

Criterion compares every run against the previous one automatically. To pin an explicit reference point:

```sh
git checkout main && cargo bench -- --save-baseline main
git checkout my-branch && cargo bench -- --baseline main
```

An HTML report with plots lands in `target/criterion/report/index.html`.

The quadratic groups (`distances`, and `bipartitions` in `traversal`) run at criterion's minimum sample count so they finish in minutes; their confidence intervals are wide by design. Read them for scaling behaviour, not for single-digit-percent regressions.

## WebAssembly

`phylo` builds for `wasm32` targets out of the box, making it suitable for in-browser phylogenetics — use your usual wasm toolchain (e.g. `wasm-pack`, or `cargo build --target wasm32-unknown-unknown`).

## Citation

If you use `phylo` in your work, please cite [this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12309125/):

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

Licensed under the [MIT License](https://github.com/sriram98v/phylo-rs/blob/main/LICENSE).

<!-- cargo-rdme end -->
