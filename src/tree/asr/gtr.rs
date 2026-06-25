use nalgebra::{DMatrix, DVector};
use crate::error::AsrError;
use crate::tree::asr::alphabet::Alphabet;

/// General Time Reversible (GTR) substitution model.
pub struct GtrModel<A: Alphabet> {
    /// Equilibrium state frequencies ($\pi$).
    pi: DVector<f64>,
    /// Symmetrized generator matrix $S$ such that $S = \text{diag}(\sqrt{\pi}) Q \text{diag}(1/\sqrt{\pi})$.
    /// $S = V \Lambda V^T$.
    eigenvalues: DVector<f64>,
    eigenvectors: DMatrix<f64>,
    sqrt_pi: DVector<f64>,
    inv_sqrt_pi: DVector<f64>,
    _phantom: std::marker::PhantomData<A>,
}

impl<A: Alphabet> GtrModel<A> {
    /// Creates a new GTR model from equilibrium frequencies `pi` and exchangeability matrix `w`.
    /// `w` must be symmetric and positive.
    pub fn new(pi: Vec<f64>, w: DMatrix<f64>, normalize: bool) -> Result<Self, AsrError> {
        let n = A::N_STATES;
        if pi.len() != n {
            return Err(AsrError::AlphabetMismatch("pi length does not match alphabet states".to_string()));
        }
        if w.nrows() != n || w.ncols() != n {
            return Err(AsrError::AlphabetMismatch("W matrix dimensions do not match alphabet states".to_string()));
        }

        let pi_vec = DVector::from_vec(pi);
        let pi_sum: f64 = pi_vec.sum();
        let pi_norm = pi_vec / pi_sum;

        // Q_ij = pi_j * W_ij for i != j
        let mut q = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
                if i != j {
                    q[(i, j)] = pi_norm[j] * w[(i, j)];
                }
            }
        }

        // Diagonal elements: Q_ii = -sum_{j!=i} Q_ij
        for i in 0..n {
            let row_sum: f64 = q.row(i).iter().sum();
            q[(i, i)] = -row_sum;
        }

        // Normalization: mu = -sum(pi_i * Q_ii) = 1
        if normalize {
            let mut mu = 0.0;
            for i in 0..n {
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

        let mut s = DMatrix::zeros(n, n);
        for i in 0..n {
            for j in 0..n {
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
            _phantom: std::marker::PhantomData,
        })
    }

    /// Returns the transition probability matrix P(t) = exp(Qt).
    pub fn transition(&self, t: f64) -> DMatrix<f64> {
        let n = A::N_STATES;
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

    /// Jukes-Cantor model: uniform pi, all W = 1.
    pub fn jukes_cantor() -> Result<Self, AsrError> {
        let n = A::N_STATES;
        let pi = vec![1.0 / (n as f64); n];
        let w = DMatrix::from_element(n, n, 1.0);
        Self::new(pi, w, true)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::tree::asr::alphabet::Nucleotide;

    #[test]
    fn test_jc_transition() {
        let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
        let p_t = model.transition(1.0);

        // In JC, all off-diagonals are equal, all diagonals are equal
        let diag = p_t[(0, 0)];
        let off_diag = p_t[(0, 1)];

        assert!((diag - off_diag).abs() > 0.0);

        // Row sums must be 1
        for i in 0..4 {
            let sum: f64 = p_t.row(i).iter().sum();
            assert!((sum - 1.0).abs() < 1e-10);
        }
    }
}
