use crate::alphabet::{Alphabet, Nucleotide};
use crate::error::AsrError;
use crate::models::gamma::discrete_gamma;
use crate::models::rate_matrix::RateMatrix;
use nalgebra::DMatrix;

/// A single rate category in a rate-heterogeneous substitution model: a relative
/// substitution rate and its mixture weight. Weights across all categories of a model
/// sum to 1.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RateCategory {
    /// Relative substitution rate for this category (branch lengths are scaled by this
    /// before computing the transition matrix). A rate of `0.0` models an invariant site.
    pub rate: f64,
    /// Mixture weight of this category.
    pub weight: f64,
}

/// General Time Reversible model with a proportion of invariant sites and discrete-gamma
/// rate heterogeneity across variable sites: **GTR+I+G**.
///
/// This is the most general nucleotide substitution model in the Posada & Crandall (2001)
/// hierarchy (Sysbio 50(4):580); every other named nucleotide model (JC69, K80, F81,
/// HKY85, TN93, K81, TIM, TVM, SYM, GTR) is a special case obtainable via a named
/// constructor that constrains its equilibrium frequencies and/or exchangeability matrix.
/// `+I` and `+G` are optional decorations applied via [`GtrModel::with_invariant`] and
/// [`GtrModel::with_gamma`]; a freshly-constructed model has neither (a single rate
/// category with rate 1.0 and weight 1.0), which reduces exactly to the plain substitution
/// model.
pub struct GtrModel<A: Alphabet> {
    matrix: RateMatrix,
    categories: Vec<RateCategory>,
    _phantom: std::marker::PhantomData<A>,
}

impl<A: Alphabet> GtrModel<A> {
    /// Creates a new GTR model (no rate heterogeneity) from equilibrium frequencies `pi`
    /// and exchangeability matrix `w`. `w` must be symmetric and positive.
    pub fn new(pi: Vec<f64>, w: DMatrix<f64>, normalize: bool) -> Result<Self, AsrError> {
        let matrix = RateMatrix::new(A::N_STATES, pi, w, normalize)?;
        Ok(Self {
            matrix,
            categories: vec![RateCategory {
                rate: 1.0,
                weight: 1.0,
            }],
            _phantom: std::marker::PhantomData,
        })
    }

    /// Jukes-Cantor (JC69) model: uniform equilibrium frequencies, all exchangeabilities equal.
    pub fn jukes_cantor() -> Result<Self, AsrError> {
        let n = A::N_STATES;
        let pi = vec![1.0 / (n as f64); n];
        let w = DMatrix::from_element(n, n, 1.0);
        Self::new(pi, w, true)
    }

    /// Returns the transition probability matrix `P(t)` for the given rate category,
    /// under the plain (unscaled by category) generator.
    pub fn category_transition(&self, category: usize, t: f64) -> DMatrix<f64> {
        let rate = self.categories[category].rate;
        self.matrix.transition(rate * t)
    }

    /// Returns the equilibrium frequencies.
    pub fn equilibrium(&self) -> &nalgebra::DVector<f64> {
        self.matrix.equilibrium()
    }

    /// Returns the rate categories (rate, weight) making up this model's rate heterogeneity.
    pub fn categories(&self) -> &[RateCategory] {
        &self.categories
    }

    /// Number of rate categories.
    pub fn n_categories(&self) -> usize {
        self.categories.len()
    }

    /// Returns a copy of this model with discrete-gamma rate heterogeneity across `k`
    /// categories with shape parameter `alpha` (the "+G" decoration). Replaces any
    /// previously-set gamma categories; preserves an existing proportion of invariant
    /// sites set via [`GtrModel::with_invariant`].
    pub fn with_gamma(&self, alpha: f64, k: usize) -> Result<Self, AsrError>
    where
        Self: Sized,
    {
        let p_inv = self.proportion_invariant();
        let gamma_rates = discrete_gamma(alpha, k)?;
        let variable_weight = 1.0 - p_inv;

        let mut categories = Vec::with_capacity(k + if p_inv > 0.0 { 1 } else { 0 });
        if p_inv > 0.0 {
            categories.push(RateCategory {
                rate: 0.0,
                weight: p_inv,
            });
        }
        for rate in gamma_rates {
            categories.push(RateCategory {
                rate,
                weight: variable_weight / (k as f64),
            });
        }

        Ok(Self {
            matrix: self.matrix.clone(),
            categories,
            _phantom: std::marker::PhantomData,
        })
    }

    /// Returns a copy of this model with a proportion `p_inv` of invariant sites (the "+I"
    /// decoration). `p_inv` must be in `[0, 1)`. Reweights any existing rate categories
    /// (e.g. from [`GtrModel::with_gamma`]) so all weights still sum to 1.
    pub fn with_invariant(&self, p_inv: f64) -> Result<Self, AsrError>
    where
        Self: Sized,
    {
        if !(0.0..1.0).contains(&p_inv) {
            return Err(AsrError::InvalidModelParameter(
                "proportion of invariant sites must be in [0, 1)".to_string(),
            ));
        }

        let old_p_inv = self.proportion_invariant();
        let variable_categories: Vec<&RateCategory> =
            self.categories.iter().filter(|c| c.rate != 0.0).collect();
        let old_variable_weight = 1.0 - old_p_inv;
        let new_variable_weight = 1.0 - p_inv;

        let mut categories = Vec::with_capacity(variable_categories.len() + 1);
        if p_inv > 0.0 {
            categories.push(RateCategory {
                rate: 0.0,
                weight: p_inv,
            });
        }
        for c in variable_categories {
            let scale = if old_variable_weight > 0.0 {
                new_variable_weight / old_variable_weight
            } else {
                0.0
            };
            categories.push(RateCategory {
                rate: c.rate,
                weight: c.weight * scale,
            });
        }
        if categories.iter().all(|c| c.rate == 0.0) {
            // No variable categories existed (e.g. fresh model): fall back to a single
            // unit-rate variable category carrying the full variable weight.
            categories.push(RateCategory {
                rate: 1.0,
                weight: new_variable_weight,
            });
        }

        Ok(Self {
            matrix: self.matrix.clone(),
            categories,
            _phantom: std::marker::PhantomData,
        })
    }

    /// Returns the total weight currently assigned to the invariant (rate == 0) category,
    /// or 0.0 if none is present.
    fn proportion_invariant(&self) -> f64 {
        self.categories
            .iter()
            .filter(|c| c.rate == 0.0)
            .map(|c| c.weight)
            .sum()
    }
}

/// Named constructors for special cases of GTR+I+G, nested per Posada & Crandall (2001)
/// Fig 5. All constructors accept nucleotide alphabets only.
impl GtrModel<Nucleotide> {
    /// K80 (Kimura 1980): equal equilibrium frequencies, transition/transversion ratio `kappa`.
    pub fn k80(kappa: f64) -> Result<Self, AsrError> {
        validate_positive(kappa, "kappa")?;
        let pi = vec![0.25; 4];
        let w = ti_tv_matrix(kappa, 1.0);
        Self::new(pi, w, true)
    }

    /// F81 (Felsenstein 1981): free equilibrium frequencies, all exchangeabilities equal.
    pub fn f81(pi: [f64; 4]) -> Result<Self, AsrError> {
        let w = DMatrix::from_element(4, 4, 1.0);
        Self::new(pi.to_vec(), w, true)
    }

    /// HKY85 (Hasegawa, Kishino, Yano 1985): free equilibrium frequencies, transition/
    /// transversion ratio `kappa`.
    pub fn hky85(pi: [f64; 4], kappa: f64) -> Result<Self, AsrError> {
        validate_positive(kappa, "kappa")?;
        let w = ti_tv_matrix(kappa, 1.0);
        Self::new(pi.to_vec(), w, true)
    }

    /// TN93 (Tamura & Nei 1993): free equilibrium frequencies, separate rates `kappa1`
    /// (A<->G transitions) and `kappa2` (C<->T transitions), transversions rate 1.
    pub fn tn93(pi: [f64; 4], kappa1: f64, kappa2: f64) -> Result<Self, AsrError> {
        validate_positive(kappa1, "kappa1")?;
        validate_positive(kappa2, "kappa2")?;
        // State order A=0, C=1, G=2, T=3.
        let w = symmetric_matrix(1.0, kappa1, 1.0, 1.0, kappa2, 1.0);
        Self::new(pi.to_vec(), w, true)
    }

    /// K81 (Kimura 1981, a.k.a. TPM/K3ST): equal equilibrium frequencies, three
    /// transversion/transition rate classes: `a` (A<->G, C<->T), `b` (A<->C, G<->T),
    /// `c` (A<->T, C<->G).
    pub fn k81(a: f64, b: f64, c: f64) -> Result<Self, AsrError> {
        validate_positive(a, "a")?;
        validate_positive(b, "b")?;
        validate_positive(c, "c")?;
        let pi = vec![0.25; 4];
        // AC=b, AG=a, AT=c, CG=c, CT=a, GT=b
        let w = symmetric_matrix(b, a, c, c, a, b);
        Self::new(pi, w, true)
    }

    /// TIM (Transition model): free equilibrium frequencies, two transition rates
    /// (`kappa1` for A<->G, `kappa2` for C<->T) and two transversion rates (`a` for
    /// A<->C and G<->T, `b` for A<->T and C<->G).
    pub fn tim(pi: [f64; 4], kappa1: f64, kappa2: f64, a: f64, b: f64) -> Result<Self, AsrError> {
        validate_positive(kappa1, "kappa1")?;
        validate_positive(kappa2, "kappa2")?;
        validate_positive(a, "a")?;
        validate_positive(b, "b")?;
        // AC=a, AG=kappa1, AT=b, CG=b, CT=kappa2, GT=a
        let w = symmetric_matrix(a, kappa1, b, b, kappa2, a);
        Self::new(pi.to_vec(), w, true)
    }

    /// TVM (Transversion model): free equilibrium frequencies, transitions share a single
    /// rate `kappa`, transversions split into three free rates `a`, `b`, `c`.
    pub fn tvm(pi: [f64; 4], kappa: f64, a: f64, b: f64, c: f64) -> Result<Self, AsrError> {
        validate_positive(kappa, "kappa")?;
        validate_positive(a, "a")?;
        validate_positive(b, "b")?;
        validate_positive(c, "c")?;
        // AC=a, AG=kappa, AT=b, CG=c, CT=kappa, GT=... (TVM has AC, AT, CG, GT free; here
        // parameterized as: AC=a, AT=b, CG=c, GT=(a*b*c-consistent free 4th), transitions=kappa)
        let w = symmetric_matrix(a, kappa, b, c, kappa, (a * b * c).cbrt());
        Self::new(pi.to_vec(), w, true)
    }

    /// SYM (Symmetric model, Zharkikh 1994): equal equilibrium frequencies, all six
    /// exchangeabilities free.
    pub fn sym(rates: [f64; 6]) -> Result<Self, AsrError> {
        for (name, r) in ["ac", "ag", "at", "cg", "ct", "gt"]
            .iter()
            .zip(rates.iter())
        {
            validate_positive(*r, name)?;
        }
        let pi = vec![0.25; 4];
        let w = symmetric_matrix(rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]);
        Self::new(pi, w, true)
    }

    /// GTR (General Time Reversible, Tavare 1986): free equilibrium frequencies, all six
    /// exchangeabilities free. This is the fully general base model (before +I/+G).
    pub fn gtr(pi: [f64; 4], rates: [f64; 6]) -> Result<Self, AsrError> {
        for (name, r) in ["ac", "ag", "at", "cg", "ct", "gt"]
            .iter()
            .zip(rates.iter())
        {
            validate_positive(*r, name)?;
        }
        let w = symmetric_matrix(rates[0], rates[1], rates[2], rates[3], rates[4], rates[5]);
        Self::new(pi.to_vec(), w, true)
    }
}

fn validate_positive(value: f64, name: &str) -> Result<(), AsrError> {
    if value <= 0.0 || !value.is_finite() {
        return Err(AsrError::InvalidModelParameter(format!(
            "{name} must be positive and finite"
        )));
    }
    Ok(())
}

/// Builds the symmetric exchangeability matrix for state order A=0, C=1, G=2, T=3 from
/// the six pairwise rates: AC, AG, AT, CG, CT, GT.
fn symmetric_matrix(ac: f64, ag: f64, at: f64, cg: f64, ct: f64, gt: f64) -> DMatrix<f64> {
    let mut w = DMatrix::from_element(4, 4, 0.0);
    let pairs = [
        (0, 1, ac),
        (0, 2, ag),
        (0, 3, at),
        (1, 2, cg),
        (1, 3, ct),
        (2, 3, gt),
    ];
    for (i, j, r) in pairs {
        w[(i, j)] = r;
        w[(j, i)] = r;
    }
    w
}

/// Builds the exchangeability matrix for a transition/transversion model: transitions
/// (A<->G, C<->T) get rate `kappa`, transversions get rate `tv`.
fn ti_tv_matrix(kappa: f64, tv: f64) -> DMatrix<f64> {
    symmetric_matrix(tv, kappa, tv, tv, kappa, tv)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabet::Nucleotide;

    #[test]
    fn test_jc_transition() {
        let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
        let p_t = model.category_transition(0, 1.0);

        let diag = p_t[(0, 0)];
        let off_diag = p_t[(0, 1)];
        assert!((diag - off_diag).abs() > 0.0);

        for i in 0..4 {
            let sum: f64 = p_t.row(i).iter().sum();
            assert!((sum - 1.0).abs() < 1e-10);
        }
    }

    #[test]
    fn test_default_model_has_single_unit_category() {
        let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
        assert_eq!(
            model.categories(),
            &[RateCategory {
                rate: 1.0,
                weight: 1.0
            }]
        );
    }

    #[test]
    fn test_hky85_with_kappa_one_matches_jc69() {
        let jc = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
        let hky = GtrModel::<Nucleotide>::hky85([0.25, 0.25, 0.25, 0.25], 1.0).unwrap();

        let p_jc = jc.category_transition(0, 0.5);
        let p_hky = hky.category_transition(0, 0.5);

        for i in 0..4 {
            for j in 0..4 {
                assert!((p_jc[(i, j)] - p_hky[(i, j)]).abs() < 1e-10);
            }
        }
    }

    #[test]
    fn test_with_gamma_categories_sum_to_one_weight() {
        let model = GtrModel::<Nucleotide>::jukes_cantor()
            .unwrap()
            .with_gamma(0.5, 4)
            .unwrap();
        assert_eq!(model.n_categories(), 4);
        let total_weight: f64 = model.categories().iter().map(|c| c.weight).sum();
        assert!((total_weight - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_with_invariant_adds_rate_zero_category() {
        let model = GtrModel::<Nucleotide>::jukes_cantor()
            .unwrap()
            .with_invariant(0.3)
            .unwrap();
        assert_eq!(model.n_categories(), 2);
        let inv_cat = model.categories().iter().find(|c| c.rate == 0.0).unwrap();
        assert!((inv_cat.weight - 0.3).abs() < 1e-10);
        let total_weight: f64 = model.categories().iter().map(|c| c.weight).sum();
        assert!((total_weight - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_with_gamma_then_invariant_combined() {
        let model = GtrModel::<Nucleotide>::jukes_cantor()
            .unwrap()
            .with_gamma(0.5, 4)
            .unwrap()
            .with_invariant(0.2)
            .unwrap();
        assert_eq!(model.n_categories(), 5);
        let total_weight: f64 = model.categories().iter().map(|c| c.weight).sum();
        assert!((total_weight - 1.0).abs() < 1e-10);
        let inv_weight: f64 = model
            .categories()
            .iter()
            .filter(|c| c.rate == 0.0)
            .map(|c| c.weight)
            .sum();
        assert!((inv_weight - 0.2).abs() < 1e-10);
    }

    #[test]
    fn test_invalid_invariant_proportion() {
        let model = GtrModel::<Nucleotide>::jukes_cantor().unwrap();
        assert!(model.with_invariant(1.0).is_err());
        assert!(model.with_invariant(-0.1).is_err());
    }

    #[test]
    fn test_named_constructors_validate_params() {
        assert!(GtrModel::<Nucleotide>::k80(-1.0).is_err());
        assert!(GtrModel::<Nucleotide>::hky85([0.25; 4], 0.0).is_err());
    }
}
