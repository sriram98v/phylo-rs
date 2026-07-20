//! Heap usage report for [`PhyloTree`].
//!
//! This is a reporting tool rather than a pass/fail test, so it is `#[ignore]`d
//! and does not run in CI. Run it explicitly:
//!
//! ```text
//! cargo test --release --test memory-report -- --ignored --nocapture
//! ```
//!
//! It exists to make the effect of a layout change legible: run it before and
//! after, and diff. `divan` measures time, so it is the wrong tool for bytes.

use phylo::node::PhyloNode;
use phylo::prelude::*;
use phylo::tree::PhyloTree;

const SIZES: &[usize] = &[1000, 4000, 16000];

fn kib(bytes: usize) -> String {
    format!("{:.1} KiB", bytes as f64 / 1024.0)
}

#[test]
#[ignore = "reporting tool: run with --ignored --nocapture"]
fn memory_report() {
    println!();
    println!(
        "size_of::<PhyloNode>()         = {} B",
        std::mem::size_of::<PhyloNode>()
    );
    println!(
        "size_of::<Option<PhyloNode>>() = {} B  <- the arena slot size",
        std::mem::size_of::<Option<PhyloNode>>()
    );
    println!();
    println!(
        "{:>7} {:>7} {:>12} {:>12} {:>10} {:>12} {:>8} {:>7}",
        "taxa", "nodes", "arena", "index", "taxa_map", "total", "B/node", "idx/tree"
    );

    for &n in SIZES {
        let tree = PhyloTree::yule(n);
        let nodes = tree.num_nodes();

        // The tree proper: arena + taxa map. The LCA index no longer lives on
        // the tree, so `heap_size` never includes it.
        let arena = tree.arena_heap_size();
        let taxa_map = tree.taxa_map_heap_size();
        let bare = tree.heap_size();

        // The borrowing oracle any LCA-using algorithm builds on demand.
        let index = tree.lca().heap_size();
        let total = bare + index;

        println!(
            "{:>7} {:>7} {:>12} {:>12} {:>10} {:>12} {:>8.1} {:>6.1}x",
            n,
            nodes,
            kib(arena),
            kib(index),
            kib(taxa_map),
            kib(total),
            arena as f64 / nodes as f64,
            index as f64 / bare as f64,
        );
    }
    println!();
}
