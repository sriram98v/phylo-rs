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
