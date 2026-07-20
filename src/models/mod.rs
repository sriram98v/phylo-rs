//! Substitution models for molecular evolution.
//!
//! The base model is **GTR+I+G** ([`GtrModel`](crate::models::gtr::GtrModel)): the General Time Reversible model
//! with a proportion of invariant sites (+I) and discrete-gamma rate heterogeneity (+G).
//! Every other named nucleotide model (JC69, K80, F81, HKY85, TN93, K81, TIM, TVM, SYM,
//! GTR) is exposed as a special case via a named constructor on [`GtrModel`](crate::models::gtr::GtrModel), per the
//! nested model hierarchy in Posada & Crandall (2001), Sysbio 50(4):580.

/// Eigendecomposition-based rate matrix core shared by all substitution models.
pub mod rate_matrix;

/// Discrete-gamma rate heterogeneity (Yang 1994).
pub mod gamma;

/// GTR+I+G base model and its named special cases (JC69, K80, F81, HKY85, TN93, ...).
pub mod gtr;

pub use self::gamma::discrete_gamma;
pub use self::gtr::{GtrModel, RateCategory};
pub use self::rate_matrix::RateMatrix;
