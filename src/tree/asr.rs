//! Ancestral Sequence Reconstruction (ASR).
//!
//! Reconstructs ancestral states at internal nodes of a phylogenetic tree given
//! observed sequences at the leaves.
//!
//! This is a thin layer over [`crate::tree::likelihood`], which holds the
//! substance: Felsenstein's pruning and the Viterbi recursion, generalized over
//! a model's rate categories. Reconstruction is the argmax (and, for the
//! marginal case, the posteriors) taken from that computation.

use crate::alignment::Alignment;
use crate::alphabet::Alphabet;
use crate::error::AsrError;
use crate::models::GtrModel;

pub use crate::tree::likelihood::Reconstruction;

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
