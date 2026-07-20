//! Phylogenetic likelihood under a substitution model.
//!
//! Felsenstein's pruning algorithm over a tree and an [`Alignment`], generalized
//! across the rate categories of a [`GtrModel`].
//!
//! [`compute_log_likelihood`] returns just the tree's log-likelihood — the
//! Felsenstein up pass, no ancestral states. [`compute_marginal_asr`] reuses the
//! same per-pattern pruning core (`prune_pattern_category`) and adds a pre-order
//! down pass to reconstruct states and posteriors, so the two never drift. The
//! [`TreeLikelihood`] trait exposes the likelihood-only path.
//!
//! The joint (Viterbi) engine keeps its own recursion: it maximizes rather than
//! sums over states (a different semiring), so it cannot share the marginal core.
//!
//! Making pruning generic over `RootedTree` (rather than concrete in `PhyloTree`)
//! is a deliberate non-goal here.

/// Likelihood profiles and scaling for numerical stability.
pub mod profile;

/// Result type carrying reconstructed states and the tree log-likelihood.
pub mod reconstruction;

#[cfg(test)]
mod integration_test;

pub use self::reconstruction::Reconstruction;

/// Log-likelihood of an alignment given a tree and a substitution model.
///
/// Feature-free (like [`crate::tree::asr::MarginalAsr`]) so a caller bringing its
/// own tree type can implement it without `simple_rooted_tree`.
pub trait TreeLikelihood {
    /// Natural-log likelihood of `aln` given this tree and `model`, computed with
    /// Felsenstein's pruning algorithm. No ancestral states are reconstructed.
    fn log_likelihood<A: crate::alphabet::Alphabet>(
        &self,
        model: &crate::models::GtrModel<A>,
        aln: &crate::alignment::Alignment,
    ) -> Result<f64, crate::error::AsrError>;
}

// Every engine here is concrete in PhyloTree, so the module's imports gate as a
// block. What stays available without the feature is `crate::tree::asr`, which
// is trait-level and depends on none of this.
#[cfg(feature = "simple_rooted_tree")]
use {
    self::profile::Profile, crate::alignment::Alignment, crate::alphabet::Alphabet,
    crate::error::AsrError, crate::models::GtrModel, crate::node::NodeID, crate::prelude::*,
    crate::tree::PhyloTree, nalgebra::DVector, num_traits::NumCast, std::collections::HashMap,
};

/// Numerically stable log-sum-exp: `ln(sum_i exp(xs[i]))`.
///
/// Only the gated engines below mix rate categories, so this is dead
/// code without `simple_rooted_tree`.
#[cfg(feature = "simple_rooted_tree")]
fn log_sum_exp(xs: &[f64]) -> f64 {
    let max = xs.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() {
        return max;
    }
    let sum: f64 = xs.iter().map(|&x| (x - max).exp()).sum();
    max + sum.ln()
}

/// Felsenstein's pruning up pass for a single compressed pattern under a single
/// rate category.
///
/// Runs the rooted post-order traversal, building a scaled likelihood [`Profile`]
/// at every node (branch lengths scaled by category `cat_idx`'s rate), and returns
/// the profile map together with this category's log-likelihood contribution
/// `ln(weight) + ln(root_mass) + root_log_scale`.
///
/// This is the shared core: [`compute_log_likelihood`] keeps only the returned
/// `cat_ll` (dropping the profiles), while [`compute_marginal_asr`] retains the
/// profiles to seed its down pass. There is exactly one pruning recursion.
#[cfg(feature = "simple_rooted_tree")]
#[allow(clippy::too_many_arguments)]
fn prune_pattern_category<A: Alphabet>(
    tree: &PhyloTree,
    model: &GtrModel<A>,
    cat_idx: usize,
    category_weight: f64,
    pattern: &[u8],
    leaf_id_map: &[NodeID],
    postord: &[NodeID],
    pi: &DVector<f64>,
    n_states: usize,
) -> Result<(HashMap<NodeID, Profile>, f64), AsrError> {
    let root = tree.get_root_id();
    let mut profiles: HashMap<NodeID, Profile> = HashMap::new();

    for v in postord {
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
                let p_t = model.category_transition(cat_idx, weight);

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

    let root_prof = profiles.get(&root).ok_or(AsrError::NumericalInstability)?;
    let mut root_mass = 0.0;
    for i in 0..n_states {
        root_mass += pi[i] * root_prof.values[i];
    }
    let cat_ll = category_weight.ln() + root_mass.ln() + root_prof.log_scale;

    Ok((profiles, cat_ll))
}

/// Log-likelihood of an alignment given a tree and a substitution model.
///
/// Felsenstein's pruning up pass only: for each compressed alignment pattern the
/// per-category log-likelihoods (`prune_pattern_category`) are mixed with
/// [`log_sum_exp`] and accumulated by pattern multiplicity. No down pass, and no
/// per-node sequence/posterior allocation — unlike [`compute_marginal_asr`], whose
/// `log_likelihood` field this reproduces exactly.
///
/// Concrete in `PhyloTree`, so it is gated on the feature that defines it. The
/// [`TreeLikelihood`] trait itself stays available without that feature.
#[cfg(feature = "simple_rooted_tree")]
pub fn compute_log_likelihood<A>(
    tree: &PhyloTree,
    model: &GtrModel<A>,
    aln: &Alignment,
) -> Result<f64, AsrError>
where
    A: Alphabet,
{
    let comp = aln.compress_columns();
    let root = tree.get_root_id();
    let n_states = A::N_STATES;
    let pi = model.equilibrium();
    let categories = model.categories();
    let n_categories = categories.len();

    // Map alignment leaf order to NodeIDs.
    let mut leaf_id_map = Vec::with_capacity(comp.leaf_order.len());
    for name in &comp.leaf_order {
        let node_id = tree.get_taxa_node_id(name).ok_or_else(|| {
            AsrError::AlphabetMismatch(format!("Taxon {} in alignment not found in tree", name))
        })?;
        leaf_id_map.push(node_id);
    }

    let postord = tree.postord_ids(root).collect::<Vec<_>>();

    let mut total_log_likelihood = 0.0;
    for (p_idx, pattern) in comp.patterns.iter().enumerate() {
        let multiplicity = comp.multiplicity[p_idx] as f64;

        let mut cat_log_likelihoods = Vec::with_capacity(n_categories);
        for (cat_idx, category) in categories.iter().enumerate() {
            let (_profiles, cat_ll) = prune_pattern_category(
                tree,
                model,
                cat_idx,
                category.weight,
                pattern,
                &leaf_id_map,
                &postord,
                pi,
                n_states,
            )?;
            cat_log_likelihoods.push(cat_ll);
        }

        total_log_likelihood += multiplicity * log_sum_exp(&cat_log_likelihoods);
    }

    Ok(total_log_likelihood)
}

/// Internal implementation of Marginal ASR logic.
///
/// Generalized over the model's rate categories (see [`GtrModel::categories`]): each
/// compressed alignment pattern's likelihood is computed independently per category
/// (branch lengths scaled by that category's rate) and combined via a weighted mixture
/// (log-sum-exp for the site log-likelihood; a category- and likelihood-weighted average
/// of posteriors for the marginal reconstruction). A model with a single unit-rate
/// category (the default, no `+I`/`+G`) reduces exactly to plain single-rate GTR.
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
    let categories = model.categories();
    let n_categories = categories.len();

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

    let postord = tree.postord_ids(root).collect::<Vec<_>>();
    let preord = tree.preord_ids(root).collect::<Vec<_>>();

    for (p_idx, pattern) in comp.patterns.iter().enumerate() {
        let multiplicity = comp.multiplicity[p_idx] as f64;

        // Per-category site log-likelihood (including that category's ln(weight)), and
        // per-category, per-node marginal posteriors.
        let mut cat_log_likelihoods = Vec::with_capacity(n_categories);
        let mut cat_posteriors: Vec<HashMap<NodeID, Vec<f64>>> = Vec::with_capacity(n_categories);

        for (cat_idx, category) in categories.iter().enumerate() {
            // UP pass: shared Felsenstein pruning core (see `prune_pattern_category`).
            // The profiles are retained here to seed the marginal down pass below.
            let (profiles, cat_ll) = prune_pattern_category(
                tree,
                model,
                cat_idx,
                category.weight,
                pattern,
                &leaf_id_map,
                &postord,
                pi,
                n_states,
            )?;
            cat_log_likelihoods.push(cat_ll);
            let root_prof = profiles.get(&root).unwrap();

            // DOWN pass: Pre-order traversal, marginal posteriors for this category.
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
                    let p_t = model.category_transition(cat_idx, weight);

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

            cat_posteriors.push(node_posteriors);
        }

        let site_ll = log_sum_exp(&cat_log_likelihoods);
        total_log_likelihood += multiplicity * site_ll;

        // Mixture weight for each category = its share of the total site likelihood.
        // When every category has zero likelihood (e.g. contradictory data forced by
        // zero-length branches), `site_ll` is -inf; fall back to a uniform mixture so the
        // subtraction below doesn't produce NaN (-inf - -inf).
        let cat_mix_weights: Vec<f64> = if site_ll.is_finite() {
            cat_log_likelihoods
                .iter()
                .map(|&ll| (ll - site_ll).exp())
                .collect()
        } else {
            vec![1.0 / n_categories as f64; n_categories]
        };

        let mut mixed_posteriors: HashMap<NodeID, Vec<f64>> = HashMap::new();
        for v in &preord {
            let v = *v;
            let mut mixed = vec![0.0; n_states];
            for (cat_idx, mix_w) in cat_mix_weights.iter().enumerate() {
                let post = cat_posteriors[cat_idx].get(&v).unwrap();
                for i in 0..n_states {
                    mixed[i] += mix_w * post[i];
                }
            }
            let sum: f64 = mixed.iter().sum();
            if sum > 0.0 {
                for val in &mut mixed {
                    *val /= sum;
                }
            }
            mixed_posteriors.insert(v, mixed);
        }

        for (v, post) in mixed_posteriors {
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
/// Generalized over the model's rate categories: for each compressed pattern, runs the
/// Viterbi recursion independently per category (branch lengths scaled by that category's
/// rate), then jointly maximizes over (states, rate category) by picking whichever
/// category yields the highest `ln(weight) + root log-likelihood` and tracing back its
/// pointers. A model with a single unit-rate category reduces exactly to plain joint ASR.
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
    let categories = model.categories();

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

    let postord = tree.postord_ids(root).collect::<Vec<_>>();
    let preord = tree.preord_ids(root).collect::<Vec<_>>();

    for (p_idx, pattern) in comp.patterns.iter().enumerate() {
        let multiplicity = comp.multiplicity[p_idx] as f64;

        // Best (state, root log-likelihood, pointers) per category.
        let mut best_overall_ll = f64::NEG_INFINITY;
        let mut best_overall_states: HashMap<NodeID, usize> = HashMap::new();

        for (cat_idx, category) in categories.iter().enumerate() {
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
                        let p_t = model.category_transition(cat_idx, weight);
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
            let cat_ll = category.weight.ln() + root_max_ll;

            if cat_ll > best_overall_ll {
                best_overall_ll = cat_ll;

                let mut states = HashMap::new();
                states.insert(root, root_best_state);
                for v in &preord {
                    let s_v = states[v];
                    for c in tree.get_node_children_ids(*v) {
                        let ptr_vc = ptrs.get(&(*v, c)).unwrap();
                        states.insert(c, ptr_vc[s_v]);
                    }
                }
                best_overall_states = states;
            }
        }

        total_log_likelihood += multiplicity * best_overall_ll;

        for (v, s_v) in best_overall_states {
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
