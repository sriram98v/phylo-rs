/// Numerical core for ancestral sequence reconstruction.
///
/// Handles scaling of likelihood vectors to prevent underflow during
/// the Felsenstein pruning algorithm.
pub struct Profile {
    /// Linear-space likelihoods for each state.
    pub values: Vec<f64>,
    /// Accumulated log-scaling factor for this profile.
    pub log_scale: f64,
}

impl Profile {
    /// Creates a new profile from a raw vector of likelihoods.
    pub fn new(values: Vec<f64>, log_scale: f64) -> Self {
        Self { values, log_scale }
    }

    /// Normalizes the profile by dividing by the maximum value and adding the log of that max
    /// to the log_scale. This is the standard "scaling" approach in phylogenetic ML.
    pub fn scale(mut self) -> Self {
        let max_val = self.values.iter().copied().fold(f64::NEG_INFINITY, f64::max);

        if max_val <= 0.0 {
            // All zero or negative; nothing to scale, but we keep the structure
            return self;
        }

        let scale_factor = max_val;
        for val in &mut self.values {
            *val /= scale_factor;
        }

        self.log_scale += scale_factor.ln();
        self
    }

    /// Returns the true log-likelihood of the state distribution.
    pub fn total_log_likelihood(&self) -> f64 {
        let sum: f64 = self.values.iter().sum();
        if sum <= 0.0 {
            f64::NEG_INFINITY
        } else {
            sum.ln() + self.log_scale
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_scaling() {
        // Small values that would underflow if multiplied repeatedly
        let vals = vec![1e-20, 2e-20, 3e-20];
        let profile = Profile::new(vals, 0.0).scale();

        // Use tolerance for f64 comparisons (1.0/3.0 ≠ exact binary)
        assert!((profile.values[0] - 1.0 / 3.0).abs() < 1e-15);
        assert!((profile.values[1] - 2.0 / 3.0).abs() < 1e-15);
        assert!(profile.values[2] - 1.0 < 1e-15);
        assert!((profile.log_scale - 3e-20f64.ln()).abs() < 1e-10);
    }

    #[test]
    fn test_log_likelihood() {
        let vals = vec![0.1, 0.2, 0.3, 0.4];
        let scale = 10.0f64.ln();
        let profile = Profile::new(vals, scale);

        // log(1.0) + ln(10) = 2.3025...
        assert!((profile.total_log_likelihood() - scale).abs() < 1e-10);
    }
}
