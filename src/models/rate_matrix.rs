use crate::error::AsrError;
use nalgebra::{DMatrix, DVector};

/// Core eigendecomposition engine for a reversible substitution rate matrix.
///
/// Given equilibrium frequencies `pi` and a symmetric exchangeability matrix `w`,
/// builds the instantaneous rate matrix `Q`, symmetrizes it for a stable
/// eigendecomposition, and exposes `transition(t)` to compute `P(t) = exp(Qt)`.
///
/// This is the numerical core shared by all nucleotide substitution models
/// (JC69, K80, F81, HKY85, TN93, SYM, GTR, ...); see [`crate::models::gtr::GtrModel`]
/// for the higher-level model that wraps this with rate heterogeneity (+G) and a
/// proportion of invariant sites (+I).
#[derive(Clone, Debug)]
pub struct RateMatrix {
    /// Equilibrium state frequencies ($\pi$).
    pi: DVector<f64>,
    /// Eigenvalues of the symmetrized generator matrix $S = V \Lambda V^T$.
    eigenvalues: DVector<f64>,
    /// Eigenvectors of the symmetrized generator matrix.
    eigenvectors: DMatrix<f64>,
    sqrt_pi: DVector<f64>,
    inv_sqrt_pi: DVector<f64>,
    n_states: usize,
}

impl RateMatrix {
    /// Creates a new rate matrix from equilibrium frequencies `pi` and exchangeability
    /// matrix `w`. `w` must be symmetric and positive. `pi` need not be pre-normalized.
    ///
    /// When `normalize` is `true`, the resulting rate matrix is scaled so that the mean
    /// substitution rate (at equilibrium) is 1, i.e. one unit of branch length
    /// corresponds to one expected substitution per site.
    pub fn new(
        n_states: usize,
        pi: Vec<f64>,
        w: DMatrix<f64>,
        normalize: bool,
    ) -> Result<Self, AsrError> {
        if pi.len() != n_states {
            return Err(AsrError::AlphabetMismatch(
                "pi length does not match alphabet states".to_string(),
            ));
        }
        if w.nrows() != n_states || w.ncols() != n_states {
            return Err(AsrError::AlphabetMismatch(
                "W matrix dimensions do not match alphabet states".to_string(),
            ));
        }

        let pi_vec = DVector::from_vec(pi);
        let pi_sum: f64 = pi_vec.sum();
        if pi_sum <= 0.0 {
            return Err(AsrError::InvalidModelParameter(
                "equilibrium frequencies must sum to a positive value".to_string(),
            ));
        }
        let pi_norm = pi_vec / pi_sum;

        // Q_ij = pi_j * W_ij for i != j
        let mut q = DMatrix::zeros(n_states, n_states);
        for i in 0..n_states {
            for j in 0..n_states {
                if i != j {
                    q[(i, j)] = pi_norm[j] * w[(i, j)];
                }
            }
        }

        // Diagonal elements: Q_ii = -sum_{j!=i} Q_ij
        for i in 0..n_states {
            let row_sum: f64 = q.row(i).iter().sum();
            q[(i, i)] = -row_sum;
        }

        // Normalization: mu = -sum(pi_i * Q_ii) = 1
        if normalize {
            let mut mu = 0.0;
            for i in 0..n_states {
                mu -= pi_norm[i] * q[(i, i)];
            }
            if mu <= 0.0 {
                return Err(AsrError::NumericalInstability);
            }
            q /= mu;
        }

        // Symmetrization for eigendecomposition: S = diag(sqrt(pi)) * Q * diag(1/sqrt(pi))
        let sqrt_pi = pi_norm.map(|x| x.sqrt());
        let inv_sqrt_pi = sqrt_pi.map(|x| 1.0 / x);

        let mut s = DMatrix::zeros(n_states, n_states);
        for i in 0..n_states {
            for j in 0..n_states {
                s[(i, j)] = sqrt_pi[i] * q[(i, j)] * inv_sqrt_pi[j];
            }
        }

        // SymmetricEigen is stable and wasm-safe
        let eigen = s.symmetric_eigen();

        Ok(Self {
            pi: pi_norm,
            eigenvalues: eigen.eigenvalues,
            eigenvectors: eigen.eigenvectors,
            sqrt_pi,
            inv_sqrt_pi,
            n_states,
        })
    }

    /// Returns the transition probability matrix P(t) = exp(Qt).
    pub fn transition(&self, t: f64) -> DMatrix<f64> {
        let n = self.n_states;
        if t == 0.0 {
            return DMatrix::identity(n, n);
        }

        // exp(S*t) = V * exp(Lambda*t) * V^T
        let mut exp_lambda_t = DMatrix::zeros(n, n);
        for i in 0..n {
            exp_lambda_t[(i, i)] = (self.eigenvalues[i] * t).exp();
        }

        let s_t = &self.eigenvectors * exp_lambda_t * self.eigenvectors.transpose();

        // P(t) = diag(1/sqrt(pi)) * exp(S*t) * diag(sqrt(pi))
        let mut p_t = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                p_t[(i, j)] = self.inv_sqrt_pi[i] * s_t[(i, j)] * self.sqrt_pi[j];
            }
        }
        p_t
    }

    /// Returns the equilibrium frequencies.
    pub fn equilibrium(&self) -> &DVector<f64> {
        &self.pi
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_jc_transition() {
        let n = 4;
        let pi = vec![0.25; 4];
        let w = DMatrix::from_element(n, n, 1.0);
        let model = RateMatrix::new(n, pi, w, true).unwrap();
        let p_t = model.transition(1.0);

        let diag = p_t[(0, 0)];
        let off_diag = p_t[(0, 1)];
        assert!((diag - off_diag).abs() > 0.0);

        for i in 0..4 {
            let sum: f64 = p_t.row(i).iter().sum();
            assert!((sum - 1.0).abs() < 1e-10);
        }
    }
}
