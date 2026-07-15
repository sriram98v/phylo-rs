//! Ancestral Sequence Reconstruction (ASR) module.
//!
//! This module provides tools for reconstructing ancestral states at internal nodes
//! of a phylogenetic tree given observed sequences at the leaves.

/// Alphabet definitions for ASR (e.g., Nucleotides, Amino Acids).
pub mod alphabet;

/// Likelihood profiles and scaling for numerical stability.
pub mod profile;

/// GTR (General Time Reversible) substitution models.
pub mod gtr;

/// Multiple sequence alignment handling and column compression.
pub mod alignment;

/// Core algorithms for marginal and joint ancestral reconstruction.
pub mod reconstruction;

#[cfg(test)]
mod integration_test;

pub use self::alignment::Alignment;
pub use self::alphabet::Alphabet;
pub use self::gtr::GtrModel;
pub use self::reconstruction::Reconstruction;

use crate::error::AsrError;

// Used only by the concrete PhyloTree reconstructions below. The traits
// themselves are usable without `simple_rooted_tree`.
#[cfg(feature = "simple_rooted_tree")]
use {
    self::profile::Profile, crate::node::NodeID, crate::prelude::*, crate::tree::PhyloTree,
    nalgebra::DVector, num_traits::NumCast, std::collections::HashMap,
};

/// Trait for performing marginal ancestral sequence reconstruction.
pub trait MarginalAsr {
    /// Performs marginal ML reconstruction of ancestral sequences.
    fn marginal_asr<A: Alphabet>(
        &self,
        model: &GtrModel<A>,
        aln: &Alignment,
        want_posteriors: bool,
    ) -> Result<Reconstruction<A>, AsrError>;
}

/// Trait for performing joint ancestral sequence reconstruction.
pub trait JointAsr {
    /// Performs joint ML reconstruction of ancestral sequences (Viterbi).
    fn joint_asr<A: Alphabet>(
        &self,
        model: &GtrModel<A>,
        aln: &Alignment,
    ) -> Result<Reconstruction<A>, AsrError>;
}

/// Internal implementation of Marginal ASR logic.
///
/// Concrete in `PhyloTree`, so it is gated on the feature that defines it. The
/// `MarginalAsr` trait itself stays available without that feature, for callers
/// bringing their own tree type.
#[cfg(feature = "simple_rooted_tree")]
pub fn compute_marginal_asr<A>(
    tree: &PhyloTree,
    model: &GtrModel<A>,
    aln: &Alignment,
    want_posteriors: bool,
) -> Result<Reconstruction<A>, AsrError>
where
    A: Alphabet,
{
    let comp = aln.compress_columns();
    let root = tree.get_root_id();
    let n_states = A::N_STATES;
    let pi = model.equilibrium();

    // Map alignment leaf order to NodeIDs
    let mut leaf_id_map = Vec::with_capacity(comp.leaf_order.len());
    for name in &comp.leaf_order {
        let node_id = tree.get_taxa_node_id(name).ok_or_else(|| {
            AsrError::AlphabetMismatch(format!("Taxon {} in alignment not found in tree", name))
        })?;
        leaf_id_map.push(node_id);
    }

    let mut total_log_likelihood = 0.0;
    let mut final_sequences = HashMap::new();
    let mut final_posteriors = if want_posteriors {
        let mut map = HashMap::new();
        for node_id in tree.get_node_ids() {
            map.insert(node_id, vec![vec![0.0; n_states]; aln.width]);
        }
        Some(map)
    } else {
        None
    };

    // Initialize sequences
    for node_id in tree.get_node_ids() {
        final_sequences.insert(node_id, vec![0; aln.width]);
    }

    for (p_idx, pattern) in comp.patterns.iter().enumerate() {
        let multiplicity = comp.multiplicity[p_idx] as f64;

        // UP pass: Rooted post-order traversal
        let postord = tree.postord_ids(root).collect::<Vec<_>>();
        let mut profiles: HashMap<NodeID, Profile> = HashMap::new();

        for v in &postord {
            if tree.is_leaf(*v) {
                let pos = leaf_id_map.iter().position(|&id| id == *v).ok_or_else(|| {
                    AsrError::InvalidAlignment(
                        "Leaf in tree not found in alignment leaf order".to_string(),
                    )
                })?;
                let char_val = pattern[pos];
                let prof_vals = A::profile(char_val).ok_or_else(|| {
                    AsrError::AlphabetMismatch("Invalid char in alignment".to_string())
                })?;
                profiles.insert(*v, Profile::new(prof_vals, 0.0).scale());
            } else {
                let mut v_vals = DVector::from_element(n_states, 1.0);
                let mut sum_log_scale = 0.0;

                for c in tree.get_node_children_ids(*v) {
                    let prof_c = profiles.get(&c).ok_or(AsrError::NumericalInstability)?;
                    let weight = tree
                        .get_edge_weight(*v, c)
                        .and_then(NumCast::from)
                        .unwrap_or(0.0);
                    let p_t = model.transition(weight);

                    let child_contrib = p_t * DVector::from_vec(prof_c.values.clone());

                    for i in 0..n_states {
                        v_vals[i] *= child_contrib[i];
                    }
                    sum_log_scale += prof_c.log_scale;
                }
                profiles.insert(
                    *v,
                    Profile::new(v_vals.as_slice().to_vec(), sum_log_scale).scale(),
                );
            }
        }

        let root_prof = profiles.get(&root).unwrap();
        let mut root_mass = 0.0;
        for i in 0..n_states {
            root_mass += pi[i] * root_prof.values[i];
        }
        total_log_likelihood += multiplicity * (root_mass.ln() + root_prof.log_scale);

        // DOWN pass: Pre-order traversal
        let preord = tree.preord_ids(root).collect::<Vec<_>>();
        let mut node_posteriors: HashMap<NodeID, Vec<f64>> = HashMap::new();

        for v in &preord {
            if *v == root {
                let mut post = vec![0.0; n_states];
                let mut sum = 0.0;
                for i in 0..n_states {
                    post[i] = pi[i] * root_prof.values[i];
                    sum += post[i];
                }
                if sum > 0.0 {
                    for val in &mut post {
                        *val /= sum;
                    }
                }
                node_posteriors.insert(*v, post);
            } else {
                let p = tree.get_node_parent_id(*v).unwrap();
                let post_p = node_posteriors.get(&p).unwrap();

                let weight = tree
                    .get_edge_weight(p, *v)
                    .and_then(NumCast::from)
                    .unwrap_or(0.0);
                let p_t = model.transition(weight);

                let post_p_vec = DVector::from_vec(post_p.clone());
                let msg_vec = p_t.transpose() * post_p_vec;

                let prof_v = profiles.get(v).unwrap();
                let mut post_v = vec![0.0; n_states];
                let mut sum = 0.0;
                for i in 0..n_states {
                    post_v[i] = prof_v.values[i] * msg_vec[i];
                    sum += post_v[i];
                }
                if sum > 0.0 {
                    for val in &mut post_v {
                        *val /= sum;
                    }
                }
                node_posteriors.insert(*v, post_v);
            }
        }

        for (v, post) in node_posteriors {
            let best_state = post
                .iter()
                .enumerate()
                .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                .unwrap()
                .0;
            for site in 0..aln.width {
                if comp.site_to_pattern[site] == p_idx {
                    final_sequences.get_mut(&v).unwrap()[site] = best_state;
                    if let Some(ref mut map) = final_posteriors {
                        map.get_mut(&v).unwrap()[site] = post.clone();
                    }
                }
            }
        }
    }

    Ok(Reconstruction {
        sequences: final_sequences,
        posteriors: final_posteriors,
        log_likelihood: total_log_likelihood,
        alphabet: std::marker::PhantomData,
    })
}

/// Implementation of Joint ASR logic (Viterbi).
///
/// Concrete in `PhyloTree`, so it is gated on the feature that defines it. The
/// `JointAsr` trait itself stays available without that feature.
#[cfg(feature = "simple_rooted_tree")]
pub fn compute_joint_asr<A>(
    tree: &PhyloTree,
    model: &GtrModel<A>,
    aln: &Alignment,
) -> Result<Reconstruction<A>, AsrError>
where
    A: Alphabet,
{
    let comp = aln.compress_columns();
    let root = tree.get_root_id();
    let n_states = A::N_STATES;
    let pi = model.equilibrium();

    let mut leaf_id_map = Vec::with_capacity(comp.leaf_order.len());
    for name in &comp.leaf_order {
        let node_id = tree.get_taxa_node_id(name).ok_or_else(|| {
            AsrError::AlphabetMismatch(format!("Taxon {} in alignment not found in tree", name))
        })?;
        leaf_id_map.push(node_id);
    }

    let mut total_log_likelihood = 0.0;
    let mut final_sequences = HashMap::new();
    for node_id in tree.get_node_ids() {
        final_sequences.insert(node_id, vec![0; aln.width]);
    }

    for (p_idx, pattern) in comp.patterns.iter().enumerate() {
        let multiplicity = comp.multiplicity[p_idx] as f64;

        let postord = tree.postord_ids(root).collect::<Vec<_>>();
        let mut c_values: HashMap<NodeID, Vec<f64>> = HashMap::new();
        let mut ptrs: HashMap<(NodeID, NodeID), Vec<usize>> = HashMap::new();

        for v in &postord {
            if tree.is_leaf(*v) {
                let pos = leaf_id_map.iter().position(|&id| id == *v).ok_or_else(|| {
                    AsrError::InvalidAlignment(
                        "Leaf in tree not found in alignment leaf order".to_string(),
                    )
                })?;
                let char_val = pattern[pos];
                let prof = A::profile(char_val)
                    .ok_or_else(|| AsrError::AlphabetMismatch("Invalid char".to_string()))?;
                let c_v = prof
                    .iter()
                    .map(|&p| if p > 0.0 { p.ln() } else { f64::NEG_INFINITY })
                    .collect();
                c_values.insert(*v, c_v);
            } else {
                let mut c_v = vec![0.0; n_states];
                for c in tree.get_node_children_ids(*v) {
                    let weight = tree
                        .get_edge_weight(*v, c)
                        .and_then(NumCast::from)
                        .unwrap_or(0.0);
                    let p_t = model.transition(weight);
                    let c_c = c_values.get(&c).ok_or(AsrError::NumericalInstability)?;

                    let mut current_ptrs = vec![0; n_states];

                    for i in 0..n_states {
                        let mut local_max = f64::NEG_INFINITY;
                        let mut local_argmax = 0;
                        for j in 0..n_states {
                            let val = p_t[(i, j)].ln() + c_c[j];
                            if val > local_max {
                                local_max = val;
                                local_argmax = j;
                            }
                        }
                        c_v[i] += local_max;
                        current_ptrs[i] = local_argmax;
                    }
                    ptrs.insert((*v, c), current_ptrs);
                }
                c_values.insert(*v, c_v);
            }
        }

        let root_c = c_values.get(&root).unwrap();
        let mut root_best_state = 0;
        let mut root_max_ll = f64::NEG_INFINITY;
        for i in 0..n_states {
            let ll = pi[i].ln() + root_c[i];
            if ll > root_max_ll {
                root_max_ll = ll;
                root_best_state = i;
            }
        }
        total_log_likelihood += multiplicity * root_max_ll;

        let preord = tree.preord_ids(root).collect::<Vec<_>>();
        let mut states = HashMap::new();
        states.insert(root, root_best_state);

        for v in &preord {
            let s_v = states[v];
            for c in tree.get_node_children_ids(*v) {
                let ptr_vc = ptrs.get(&(*v, c)).unwrap();
                states.insert(c, ptr_vc[s_v]);
            }
        }

        for (v, s_v) in states {
            for site in 0..aln.width {
                if comp.site_to_pattern[site] == p_idx {
                    final_sequences.get_mut(&v).unwrap()[site] = s_v;
                }
            }
        }
    }

    Ok(Reconstruction {
        sequences: final_sequences,
        posteriors: None,
        log_likelihood: total_log_likelihood,
        alphabet: std::marker::PhantomData,
    })
}
