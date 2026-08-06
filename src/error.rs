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

/// A type for errors from tree traversal, queries and structural operations.
///
/// Every variant is reachable from caller input — an id that is not in the tree,
/// an empty query set, an annotation that was never set. States that the
/// algorithms themselves guarantee are asserted with `expect` at the point of
/// use rather than surfaced here, so no variant of this enum is unreachable.
///
/// Node ids are reported as `usize`, which every
/// [`NodeID`](crate::node::simple_rnode::RootedTreeNode::NodeID) converts into.
#[derive(Error, Debug, PartialEq, Eq, Clone)]
pub enum TreeError {
    /// The given id does not name a node of this tree
    #[error("no node with id {0} in this tree")]
    UnknownNode(usize),
    /// An operation needing at least one node was given none
    #[error("expected at least one node id, got none")]
    EmptyNodeSet,
    /// The given pair is not an edge of this tree
    #[error("({0}, {1}) is not an edge of this tree")]
    UnknownEdge(usize, usize),
    /// The node has no parent, but the operation requires one
    #[error("node {0} has no parent")]
    NoParent(usize),
    /// The operation reads branch weights, and this node has none set
    #[error("node {0} has no branch weight set")]
    MissingWeight(usize),
    /// The operation reads zeta annotations, and this node has none set
    #[error("node {0} has no zeta annotation set; call set_zeta first")]
    MissingZeta(usize),
    /// The tree does not have enough taxa for the operation
    #[error("operation needs at least {expected} taxa, tree has {actual}")]
    TooFewTaxa {
        /// Minimum number of taxa the operation requires
        expected: usize,
        /// Number of taxa the tree actually has
        actual: usize,
    },
    /// The two trees do not span the same taxa, so they cannot be compared
    #[error("trees do not span the same taxa set")]
    TaxaSetMismatch,
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
