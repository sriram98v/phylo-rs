use thiserror::Error;

/// A type for errors when parsing newick strings
#[derive(Error, Debug)]
pub enum NewickError {
    /// Invalid character in source
    #[error("invalid character at {idx}")]
    InvalidCharacter {
        /// Position of the invalid character in the source
        idx: usize,
    },
}

/// A type for errors when parsing Nexus files
#[derive(Error, Debug)]
pub enum NexusError {
    /// Invalid header format
    #[error("expected \"#NEXUS\" at the start of the input")]
    InvalidHeader,
}

/// A type for errors during ancestral sequence reconstruction
#[derive(Error, Debug)]
pub enum AsrError {
    /// A branch length was missing or non-positive
    #[error("missing or non-positive branch length")]
    MissingBranchLength,
    /// The provided alphabet is incompatible with the sequence data
    #[error("alphabet mismatch: {0}")]
    AlphabetMismatch(String),
    /// The alignment is ragged or contains invalid identifiers
    #[error("invalid alignment: {0}")]
    InvalidAlignment(String),
    /// Symmetric eigendecomposition failed to converge
    #[error("eigendecomposition failed to converge")]
    EigendecompFailure,
    /// Numerical instability encountered during scaling
    #[error("numerical instability encountered during ASR scaling")]
    NumericalInstability,
}
