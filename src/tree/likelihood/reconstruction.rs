use crate::node::NodeID;
use std::collections::HashMap;

/// The result of an ancestral sequence reconstruction.
pub struct Reconstruction<A> {
    /// The inferred state index for each node at each site.
    /// Keyed by NodeID, value is a vector of state indices from the alphabet.
    pub sequences: HashMap<NodeID, Vec<usize>>,
    /// Optional marginal posterior probabilities for each state at each site for each node.
    pub posteriors: Option<HashMap<NodeID, Vec<Vec<f64>>>>,
    /// The total log-likelihood of the reconstruction.
    pub log_likelihood: f64,
    /// The alphabet used for the reconstruction.
    pub alphabet: std::marker::PhantomData<A>,
}

impl<A: crate::alphabet::Alphabet> Reconstruction<A> {
    /// Returns the sequence for a specific node as a String of canonical characters.
    pub fn sequence_string(&self, node: NodeID) -> Option<String> {
        let seq = self.sequences.get(&node)?;
        let mut s = String::with_capacity(seq.len());
        for &idx in seq {
            s.push(A::char_of(idx) as char);
        }
        Some(s)
    }
}
