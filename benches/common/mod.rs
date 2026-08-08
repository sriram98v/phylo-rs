//! Shared fixtures and sweep sizes for the criterion benchmark targets.
//!
//! Each target declares `mod common;` and pulls in the whole module, so any
//! given target leaves some of these unused.
#![allow(dead_code)]

use itertools::Itertools;
use phylo::prelude::*;
use phylo::tree::PhyloTree;
use rand::rngs::StdRng;
use rand::seq::IteratorRandom;
use rand::SeedableRng;

/// Taxa counts for benchmarks that scale near-linearly.
///
/// Geometric rather than linear: the point of sweeping at all is to expose the
/// growth rate, and a linear sweep spends most of its time re-measuring the
/// same complexity class.
pub const TAXA: &[usize] = &[1_000, 4_000, 16_000];

/// Taxa counts for benchmarks that are quadratic in the number of taxa.
///
/// These compare every node against every node, so they need their own, much
/// smaller sweep to stay runnable.
pub const QUADRATIC_TAXA: &[usize] = &[250, 500, 1_000];

/// Normalisation exponent used by the cophenetic distance benchmarks.
pub const NORM: u32 = 1;

/// Seed for every random choice a benchmark makes.
///
/// Benchmarks are compared against saved baselines, so an input that changes
/// between runs shows up as a performance change. Everything random here is
/// drawn from this fixed seed instead of the thread RNG.
const SEED: u64 = 0x0070_796C_6F5F_7273;

/// A random Yule tree with `taxa` leaves.
pub fn yule(taxa: usize) -> PhyloTree {
    PhyloTree::yule(taxa)
}

/// Two independent Yule trees of the same size, for the comparison metrics.
pub fn yule_pair(taxa: usize) -> (PhyloTree, PhyloTree) {
    (PhyloTree::yule(taxa), PhyloTree::yule(taxa))
}

/// Two Yule trees with node depths written into zeta.
///
/// Cophenetic distance reads zeta, so it must be populated before timing.
pub fn yule_pair_with_zeta(taxa: usize) -> (PhyloTree, PhyloTree) {
    fn depth(tree: &PhyloTree, node_id: usize) -> f32 {
        tree.depth(node_id) as f32
    }

    let (mut t1, mut t2) = yule_pair(taxa);
    let _ = t1.set_zeta(depth);
    let _ = t2.set_zeta(depth);
    (t1, t2)
}

/// The `(parent, leaf)` edges of the first two leaves — a pair of SPR endpoints.
pub fn first_two_leaf_edges(tree: &PhyloTree) -> ((usize, usize), (usize, usize)) {
    let edges = tree
        .get_leaf_ids()
        .map(|leaf| (tree.get_node_parent_id(leaf).unwrap(), leaf))
        .take(2)
        .collect_vec();
    (edges[0], edges[1])
}

/// A deterministic random sample of `fraction` of the tree's leaves.
///
/// Samples actual leaf IDs. The divan benchmarks sampled from `0..taxa_size`,
/// which is a range of arena slots rather than leaves, so the subset handed to
/// `contract_tree` was mostly internal nodes.
pub fn leaf_subset(tree: &PhyloTree, fraction: f64) -> Vec<usize> {
    let mut rng = StdRng::seed_from_u64(SEED);
    let n = ((tree.num_taxa() as f64) * fraction) as usize;
    tree.get_leaf_ids().choose_multiple(&mut rng, n.max(2))
}
