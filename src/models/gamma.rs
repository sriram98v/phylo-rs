//! Discrete-gamma rate heterogeneity (Yang 1994).
//!
//! Approximates a continuous Gamma(alpha, alpha) rate distribution (mean 1, shape = rate =
//! `alpha`) with `k` equiprobable discrete categories, each represented by the conditional
//! mean of its bin. This is the standard PAML/RAxML/IQ-TREE "+G" discretization.

use crate::error::AsrError;

/// Returns `k` discrete gamma rate categories approximating Gamma(alpha, alpha) (mean 1).
///
/// Each of the `k` returned rates is the mean of the Gamma(alpha, alpha) distribution
/// conditional on falling in the `i`-th of `k` equiprobable quantile bins. The rates are
/// unweighted (each category carries weight `1/k`); callers combine this with a proportion
/// of invariant sites separately.
pub fn discrete_gamma(alpha: f64, k: usize) -> Result<Vec<f64>, AsrError> {
    if alpha <= 0.0 {
        return Err(AsrError::InvalidModelParameter(
            "gamma shape parameter alpha must be positive".to_string(),
        ));
    }
    if k == 0 {
        return Err(AsrError::InvalidModelParameter(
            "number of gamma categories must be at least 1".to_string(),
        ));
    }
    if k == 1 {
        return Ok(vec![1.0]);
    }

    let kf = k as f64;
    // Quantile boundaries of Gamma(alpha, alpha) at i/k, i = 1..k-1
    let mut bounds = Vec::with_capacity(k - 1);
    for i in 1..k {
        let p = i as f64 / kf;
        bounds.push(point_gamma(p, alpha, alpha)?);
    }

    // Conditional mean of each bin using the incomplete gamma function with shape alpha+1,
    // per Yang (1994): mean of bin i = k * [G(alpha+1, bound_i * alpha) - G(alpha+1, bound_{i-1} * alpha)] / alpha
    // where G is the regularized lower incomplete gamma function.
    let mut cdf_bounds = Vec::with_capacity(k + 1);
    cdf_bounds.push(0.0);
    for &b in &bounds {
        cdf_bounds.push(incomplete_gamma(b * alpha, alpha + 1.0));
    }
    cdf_bounds.push(1.0);

    let mut rates = Vec::with_capacity(k);
    for i in 0..k {
        let rate = kf * (cdf_bounds[i + 1] - cdf_bounds[i]);
        rates.push(rate);
    }

    Ok(rates)
}

/// Quantile function (inverse CDF) of the Gamma(shape, rate) distribution, via
/// Newton-Raphson refinement of the Wilson-Hilferty approximation. `prob` in (0, 1).
fn point_gamma(prob: f64, shape: f64, rate: f64) -> Result<f64, AsrError> {
    if !(0.0..1.0).contains(&prob) {
        return Err(AsrError::InvalidModelParameter(
            "quantile probability must be in [0, 1)".to_string(),
        ));
    }

    // Wilson-Hilferty cube-root normal approximation for the initial guess.
    let ln_gamma_shape = ln_gamma(shape);
    let ch = shape
        * (1.0 - 1.0 / (9.0 * shape) + inverse_normal_cdf(prob) * (1.0 / (9.0 * shape)).sqrt())
            .powi(3);
    // `ch` approximates the quantile in the Gamma(shape, 1) ("x*rate") space; divide by
    // `rate` to bring it back into the x-space that Newton-Raphson iterates over below.
    let mut x = if ch > 0.0 { ch / rate } else { 1e-10 };

    // Newton-Raphson on the regularized incomplete gamma function. The derivative of
    // incomplete_gamma(x*rate, shape) w.r.t. x is the Gamma(shape, rate) density at x.
    for _ in 0..100 {
        let f = incomplete_gamma(x * rate, shape) - prob;
        let ln_df = (shape - 1.0) * (x * rate).ln() - x * rate - ln_gamma_shape + rate.ln();
        let df = ln_df.exp();
        if !df.is_finite() || df.abs() < 1e-300 {
            break;
        }
        let step = f / df;
        let x_new = x - step;
        if !x_new.is_finite() || x_new <= 0.0 {
            break;
        }
        if (x_new - x).abs() < 1e-12 * x.max(1e-12) {
            x = x_new;
            break;
        }
        x = x_new;
    }

    Ok(x)
}

/// Natural log of the gamma function (Lanczos approximation).
fn ln_gamma(x: f64) -> f64 {
    const G: f64 = 7.0;
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.5203681218851,
        -1259.1392167224028,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507343278686905,
        -0.13857109526572012,
        9.984_369_578_019_572e-6,
        1.5056327351493116e-7,
    ];

    if x < 0.5 {
        // Reflection formula
        return (std::f64::consts::PI / (std::f64::consts::PI * x).sin()).ln() - ln_gamma(1.0 - x);
    }

    let x = x - 1.0;
    let mut a = COEFFS[0];
    let t = x + G + 0.5;
    for (i, c) in COEFFS.iter().enumerate().skip(1) {
        a += c / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

/// Regularized lower incomplete gamma function P(shape, x) = gamma(shape, x) / Gamma(shape).
fn incomplete_gamma(x: f64, shape: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let ln_gamma_shape = ln_gamma(shape);

    if x < shape + 1.0 {
        // Series expansion
        let mut term = 1.0 / shape;
        let mut sum = term;
        let mut n = shape;
        for _ in 0..500 {
            n += 1.0;
            term *= x / n;
            sum += term;
            if term.abs() < sum.abs() * 1e-15 {
                break;
            }
        }
        (sum * (-x + shape * x.ln() - ln_gamma_shape).exp()).clamp(0.0, 1.0)
    } else {
        // Continued fraction for the upper incomplete gamma function, then Q = 1 - P
        let mut b = x + 1.0 - shape;
        let mut c = 1e300;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..500 {
            let an = -(i as f64) * (i as f64 - shape);
            b += 2.0;
            d = an * d + b;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            c = b + an / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-15 {
                break;
            }
        }
        let q = (-x + shape * x.ln() - ln_gamma_shape).exp() * h;
        (1.0 - q).clamp(0.0, 1.0)
    }
}

/// Approximate inverse standard normal CDF (Acklam's algorithm), used only to seed
/// Newton-Raphson in [`point_gamma`].
fn inverse_normal_cdf(p: f64) -> f64 {
    // Beasley-Springer-Moro approximation.
    let a = [
        -3.969683028665376e+01,
        2.209460984245205e+02,
        -2.759285104469687e+02,
        1.383_577_518_672_69e2,
        -3.066479806614716e+01,
        2.506628277459239e+00,
    ];
    let b = [
        -5.447609879822406e+01,
        1.615858368580409e+02,
        -1.556989798598866e+02,
        6.680131188771972e+01,
        -1.328068155288572e+01,
    ];
    let c = [
        -7.784894002430293e-03,
        -3.223964580411365e-01,
        -2.400758277161838e+00,
        -2.549732539343734e+00,
        4.374664141464968e+00,
        2.938163982698783e+00,
    ];
    let d = [
        7.784695709041462e-03,
        3.224671290700398e-01,
        2.445134137142996e+00,
        3.754408661907416e+00,
    ];

    let p_low = 0.02425;
    let p_high = 1.0 - p_low;

    if p < p_low {
        let q = (-2.0 * p.ln()).sqrt();
        (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    } else if p <= p_high {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p).ln()).sqrt();
        -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_single_category_is_unit_rate() {
        let rates = discrete_gamma(0.5, 1).unwrap();
        assert_eq!(rates, vec![1.0]);
    }

    #[test]
    fn test_weighted_mean_is_one() {
        for &alpha in &[0.1, 0.5, 1.0, 2.0, 5.0] {
            let rates = discrete_gamma(alpha, 4).unwrap();
            let mean: f64 = rates.iter().sum::<f64>() / rates.len() as f64;
            assert!((mean - 1.0).abs() < 1e-3, "alpha={alpha}, mean={mean}");
        }
    }

    #[test]
    fn test_rates_are_increasing() {
        let rates = discrete_gamma(1.0, 4).unwrap();
        for i in 1..rates.len() {
            assert!(rates[i] > rates[i - 1]);
        }
    }

    #[test]
    fn test_large_alpha_reduces_heterogeneity() {
        // As alpha -> infinity, the gamma distribution concentrates around 1: rates converge.
        let rates_small_alpha = discrete_gamma(0.1, 4).unwrap();
        let rates_large_alpha = discrete_gamma(50.0, 4).unwrap();

        let spread = |r: &[f64]| {
            r.iter().cloned().fold(f64::MIN, f64::max) - r.iter().cloned().fold(f64::MAX, f64::min)
        };
        assert!(spread(&rates_large_alpha) < spread(&rates_small_alpha));
    }

    #[test]
    fn test_invalid_params() {
        assert!(discrete_gamma(0.0, 4).is_err());
        assert!(discrete_gamma(-1.0, 4).is_err());
        assert!(discrete_gamma(1.0, 0).is_err());
    }
}
