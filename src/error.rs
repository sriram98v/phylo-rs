use thiserror::Error;

/// A type for errors when parsing newick strings
#[derive(Error, Debug)]
pub enum NewickError {
    /// Invalid character in source
    #[error("invalid character at byte {idx}")]
    InvalidCharacter {
        /// Byte offset of the invalid character in the source
        idx: usize,
    },
    /// A `)` or `,` appeared with no matching `(`, or a `(` was never closed
    #[error("unbalanced parentheses at byte {idx}")]
    UnbalancedParens {
        /// Byte offset at which the imbalance was detected
        idx: usize,
    },
    /// A quoted label (`'...'`) was never closed
    #[error("unterminated quoted label starting at byte {idx}")]
    UnterminatedQuote {
        /// Byte offset of the opening quote
        idx: usize,
    },
    /// A `[...]` comment was never closed
    #[error("unterminated comment starting at byte {idx}")]
    UnterminatedComment {
        /// Byte offset of the opening bracket
        idx: usize,
    },
    /// The input contained no tree
    #[error("empty input: no tree found")]
    Empty,
}

/// A type for errors when parsing Nexus files
#[derive(Error, Debug)]
pub enum NexusError {
    /// Invalid header format
    #[error("expected \"#NEXUS\" at the start of the input")]
    InvalidHeader,
    /// No parseable tree definition was found
    #[error("no tree definition (expected a \"... = <newick>;\" entry in a TREES block)")]
    MissingTreeBlock,
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
    /// A substitution model parameter (e.g. kappa, alpha, p_inv) was out of range
    #[error("invalid model parameter: {0}")]
    InvalidModelParameter(String),
}
