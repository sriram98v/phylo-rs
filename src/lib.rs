#![warn(missing_docs)]

//! A fast, extensible, WebAssembly-ready phylogenetics library for Rust.
//!
//! `phylo` provides memory-efficient data structures and algorithms for
//! phylogenetic analysis and inference — from tree manipulation (SPR, NNI,
//! rerooting) to tree statistics (phylogenetic diversity, RF distance, cophenetic
//! distance) to maximum-likelihood modelling (GTR+I+G substitution models,
//! Felsenstein pruning, ancestral reconstruction). It leans on Rust's memory
//! safety, speed, and native WebAssembly support to stay both fast and portable.
//!
//! Tree traversals and operations are exposed as **derivable traits**, so you get
//! DFS/BFS/pre-/post-order, Euler tours, LCA queries, and distance metrics for
//! free on your own types — and a ready-made [`PhyloTree`](crate::tree::PhyloTree)
//! when you don't want to implement one.
//!
//! # Highlights
//!
//! - **Trait-first design** — compose narrow traits (`RootedTree`,
//!   `RootedMetaTree`, `EulerWalk`, `DFS`, `Clusters`, …) onto any type, or use
//!   the batteries-included [`PhyloTree`](crate::tree::PhyloTree).
//! - **Arena-allocated trees** — cache-friendly `Vec`-backed storage with `usize`
//!   node IDs.
//! - **Constant-time LCA** — an [`LcaOracle`](crate::iter::lca::LcaOracle) borrows
//!   the tree immutably and answers LCA queries in O(1) via an Euler tour + RMQ.
//! - **Tree comparison** — Robinson-Foulds, weighted RF, cluster affinity, and
//!   cophenetic distance, with distance-matrix builders.
//! - **Maximum-likelihood modelling** — GTR+I+G substitution models (JC69 through
//!   GTR), Felsenstein-pruning log-likelihood, and marginal/joint ancestral
//!   sequence reconstruction.
//! - **I/O** — Newick and Nexus parsing and serialization.
//! - **Simulation** — random trees (Yule, uniform).
//! - **Optional parallelism** — opt into `rayon`-backed computation with the
//!   `parallel` feature.
//!
//! # Feature flags
//!
//! | Feature | Default | Description |
//! | --- | :---: | --- |
//! | `simple_rooted_tree` | Yes | The concrete `SimpleRootedTree` / `PhyloTree` implementation. |
//! | `non_crypto_hash` | Yes | Use `fxhash` maps/sets instead of `std` for speed. |
//! | `parallel` | | `rayon`-based parallel computation for the heavy metrics. |
//! | `serde` | | `Serialize`/`Deserialize` for trees. |
//!
//! # Quick start
//!
//! Everything you need is in the prelude:
//!
//! ```
//! use phylo::prelude::*;
//! ```
//!
//! ## Build a tree
//!
//! Create an empty tree, then attach children to node IDs:
//!
//! ```
//! use phylo::prelude::*;
//!
//! let mut tree = PhyloTree::new(1);
//!
//! tree.add_child(tree.get_root_id(), PhyloNode::new(2));
//! tree.add_child(tree.get_root_id(), PhyloNode::new(3));
//! tree.add_child(2, PhyloNode::new(4));
//! tree.add_child(2, PhyloNode::new(5));
//! ```
//!
//! ## Read and write Newick
//!
//! ```
//! use phylo::prelude::*;
//!
//! let tree = PhyloTree::from_newick("((A:0.1,B:0.2),C:0.6);".as_bytes()).unwrap();
//! let newick = tree.to_newick();
//! ```
//!
//! ## Traverse
//!
//! Traversals return an [`Iterator`] of nodes or node IDs in visiting order:
//!
//! ```
//! use phylo::prelude::*;
//!
//! let tree = PhyloTree::from_newick("((A:0.1,B:0.2),C:0.6);".as_bytes()).unwrap();
//!
//! let dfs = tree.dfs(tree.get_root_id());
//! let bfs = tree.bfs_ids(tree.get_root_id());
//! let postorder = tree.postord_ids(tree.get_root_id());
//! ```
//!
//! ## Constant-time LCA
//!
//! Build an [`LcaOracle`](crate::iter::lca::LcaOracle) with `tree.lca()`; it
//! borrows the tree immutably (so staleness is a compile error, not a runtime bug)
//! and answers queries in O(1):
//!
//! ```
//! use phylo::prelude::*;
//!
//! let tree = PhyloTree::from_newick("((A,B),(C,D));".as_bytes()).unwrap();
//!
//! let a = tree.get_taxa_node_id(&"A".to_string()).unwrap();
//! let b = tree.get_taxa_node_id(&"B".to_string()).unwrap();
//!
//! let lca = tree.lca();
//! let ancestor = lca.get_lca_id(&[a, b]);
//! ```
//!
//! ## Compare trees
//!
//! Metrics account for both topology and branch lengths:
//!
//! ```
//! use phylo::prelude::*;
//!
//! fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
//!     tree.depth(node_id) as f32
//! }
//!
//! let mut tree_1 = PhyloTree::from_newick("((A:0.1,B:0.2):0.6,(C:0.3,D:0.4):0.5);".as_bytes()).unwrap();
//! let mut tree_2 = PhyloTree::from_newick("((D:0.3,C:0.4):0.5,(B:0.2,A:0.1):0.6);".as_bytes()).unwrap();
//!
//! let _ = tree_1.set_zeta(depth);
//! let _ = tree_2.set_zeta(depth);
//!
//! let cluster_affinity = tree_1.ca(&tree_2);
//! let cophenetic = tree_1.cophen_dist(&tree_2, 2);
//! ```
//!
//! ## Likelihood and ancestral reconstruction
//!
//! Score an alignment against a tree under a substitution model, or reconstruct
//! ancestral sequences at the internal nodes. The log-likelihood path
//! ([`TreeLikelihood`](crate::tree::likelihood::TreeLikelihood)) runs Felsenstein's
//! pruning algorithm alone — no reconstruction — while marginal/joint ASR
//! ([`MarginalAsr`](crate::tree::asr::MarginalAsr) /
//! [`JointAsr`](crate::tree::asr::JointAsr)) build on the same pruning core:
//!
//! ```
//! use phylo::prelude::*;
//!
//! let tree =
//!     PhyloTree::from_newick("((A:0.1,B:0.2):0.15,(C:0.3,D:0.1):0.05);".as_bytes()).unwrap();
//!
//! // A nucleotide alignment in FASTA — one sequence per leaf taxon.
//! let fasta = b">A\nACGTACGT\n>B\nACGTATGT\n>C\nACGAACGT\n>D\nTCGTACGA\n";
//! let aln = Alignment::from_fasta_bytes(fasta).unwrap();
//!
//! // HKY85 with gamma-distributed rate heterogeneity (+G, 4 categories).
//! let model = GtrModel::<Nucleotide>::hky85([0.25, 0.25, 0.25, 0.25], 2.0)
//!     .unwrap()
//!     .with_gamma(0.5, 4)
//!     .unwrap();
//!
//! // Log-likelihood of the alignment given the tree and model (pruning only).
//! let log_lik = tree.log_likelihood::<Nucleotide>(&model, &aln).unwrap();
//! assert!(log_lik.is_finite());
//!
//! // Marginal ancestral sequence reconstruction fills the internal nodes.
//! let recon = tree.marginal_asr::<Nucleotide>(&model, &aln, false).unwrap();
//! let root_sequence = recon.sequence_string(tree.get_root_id());
//! ```
//!
//! # Module map
//!
//! | Module | What it does |
//! | --- | --- |
//! | [`tree::simple_rtree`] | Core tree traits and `SimpleRootedTree`. |
//! | [`tree::ops`] | Mutating operations: SPR, NNI, reroot, contraction, subtree extraction. |
//! | [`tree::distances`] | RF, weighted RF, cluster affinity, cophenetic distance, distance matrices. |
//! | [`tree::io`] | Newick and Nexus reading/writing. |
//! | [`tree::simulation`] | Random tree generation. |
//! | [`iter`] | Traversals, Euler walks, and the LCA oracle. |
//! | [`models`] | GTR+I+G substitution models and their named special cases. |
//! | [`tree::likelihood`] | Felsenstein-pruning log-likelihood. |
//! | [`tree::asr`] | Marginal and joint ancestral sequence reconstruction. |
//!
//! # Examples
//!
//! Runnable analyses live in the `examples/` directory of the repository (e.g.
//! `cargo run --example phylogenetic-diversity`, `cargo run --example
//! pairwise-distances`). See the repository README for how to visualize their
//! output.
//!
//! # WebAssembly
//!
//! `phylo` builds for `wasm32` targets out of the box, making it suitable for
//! in-browser phylogenetics — use your usual wasm toolchain (e.g. `wasm-pack`, or
//! `cargo build --target wasm32-unknown-unknown`).
//!
//! # Citation
//!
//! If you use `phylo` in your work, please cite
//! [this paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12309125/):
//!
//! ```bibtex
//! @article{vijendran2025phylo,
//!   title={Phylo-rs: an extensible phylogenetic analysis library in rust},
//!   author={Vijendran, Sriram and Anderson, Tavis and Markin, Alexey and Eulenstein, Oliver},
//!   journal={BMC bioinformatics},
//!   volume={26},
//!   pages={197},
//!   year={2025}
//! }
//! ```
//!
//! # License
//!
//! Licensed under the MIT License.

/// Module with multiple sequence alignments and column compression.
pub mod alignment;
/// Module with sequence alphabets (nucleotides, amino acids).
pub mod alphabet;
/// Module with errors.
pub mod error;
/// Module with tree traversal iterator traits and structs
pub mod iter;
/// Module with tree node traits and structs
pub mod node;
/// Module with tree traits and structs
pub mod tree;

/// Module with substitution models for molecular evolution (GTR+I+G and special cases).
pub mod models;

/// Prelude module that imports all active and tested traits along with any required struct and type alias.
pub mod prelude {
    #[doc(no_inline)]
    pub use crate::alignment::*;
    #[doc(no_inline)]
    pub use crate::alphabet::*;
    #[doc(no_inline)]
    pub use crate::error::*;
    #[doc(no_inline)]
    pub use crate::iter::lca::*;
    #[doc(no_inline)]
    pub use crate::iter::node_iter::*;
    pub use crate::models::*;
    #[doc(no_inline)]
    pub use crate::node::{simple_rnode::*, Node, PhyloNode};
    #[doc(no_inline)]
    pub use crate::tree::asr::*;
    #[doc(no_inline)]
    pub use crate::tree::distances::*;
    #[doc(no_inline)]
    pub use crate::tree::io::*;
    #[doc(no_inline)]
    pub use crate::tree::likelihood::*;
    #[doc(no_inline)]
    pub use crate::tree::ops::*;
    #[doc(no_inline)]
    pub use crate::tree::simple_rtree::*;
    #[doc(no_inline)]
    pub use crate::tree::simulation::*;

    #[cfg(feature = "simple_rooted_tree")]
    pub use crate::tree::{PhyloTree, SimpleRootedTree};
}
