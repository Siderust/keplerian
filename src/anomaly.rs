// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed anomaly solvers and anomaly conversions.
//!
//! ## Scientific scope
//! This module implements the standard elliptic, parabolic, and hyperbolic
//! anomaly relations used in two-body Keplerian motion. The numerical methods
//! target ordinary astrodynamics workloads and do not model perturbations.
//!
//! ## Technical scope
//! Public APIs use [`qtty::angular::Radians`] for angular quantities and
//! [`crate::eccentricity::Eccentricity`] for the conic parameter. Private
//! numeric kernels remain raw `f64`.
//!
//! ## References
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

use core::f64::consts::PI;

use qtty::angular::Radians;

use crate::eccentricity::Eccentricity;

/// Options controlling iterative Kepler-equation solves.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AnomalyOptions {
    /// Maximum Newton or bisection iterations.
    pub max_iter: u32,
    /// Absolute residual tolerance in raw radians.
    pub tol: f64,
}

impl Default for AnomalyOptions {
    fn default() -> Self {
        Self {
            max_iter: 64,
            tol: 1.0e-12,
        }
    }
}

/// Errors returned by anomaly solvers.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum AnomalyError {
    /// Iteration did not converge within the configured limit.
    #[error(
        "Kepler solver did not converge after {iterations} iterations (residual {residual:e})"
    )]
    NotConverged {
        /// Iterations completed.
        iterations: u32,
        /// Absolute residual at the final iterate.
        residual: f64,
    },
    /// Eccentricity is outside the valid range for the selected conic regime.
    #[error("invalid eccentricity {0}")]
    InvalidEccentricity(f64),
    /// Mean anomaly is not finite.
    #[error("invalid mean anomaly {0}")]
    InvalidMeanAnomaly(f64),
}

/// Solves the elliptic Kepler equation `M = E - e sin(E)`.
///
/// # Errors
///
/// Returns [`AnomalyError::InvalidEccentricity`] if `ecc` is not in `[0, 1)`,
/// or [`AnomalyError::NotConverged`] if the solver fails to reach tolerance.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{kepler_elliptic, AnomalyOptions};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let e = kepler_elliptic(Radians::new(0.5), Eccentricity::new_unchecked(0.1), AnomalyOptions::default()).unwrap();
/// assert!((e.value() - 0.5524799869).abs() < 1e-9);
/// ```
pub fn kepler_elliptic(
    mean_anomaly: Radians,
    ecc: Eccentricity,
    opts: AnomalyOptions,
) -> Result<Radians, AnomalyError> {
    let ecc_value = ecc.value();
    let mean_value = mean_anomaly.value();
    if !(0.0..1.0).contains(&ecc_value) || !ecc_value.is_finite() {
        return Err(AnomalyError::InvalidEccentricity(ecc_value));
    }
    if !mean_value.is_finite() {
        return Err(AnomalyError::InvalidMeanAnomaly(mean_value));
    }
    if ecc_value == 0.0 {
        return Ok(mean_anomaly);
    }
    let mut e_anom = mean_value + ecc_value * mean_value.sin();
    for i in 0..opts.max_iter {
        let mut residual = elliptic_residual(e_anom, ecc_value, mean_value);
        if residual.abs() <= opts.tol {
            return Ok(Radians::new(e_anom));
        }
        let fp = 1.0 - ecc_value * e_anom.cos();
        e_anom -= residual / fp;
        if !e_anom.is_finite() {
            break;
        }
        if i + 1 == opts.max_iter {
            residual = elliptic_residual(e_anom, ecc_value, mean_value);
            if residual.abs() <= opts.tol {
                return Ok(Radians::new(e_anom));
            }
        }
    }
    elliptic_bisection(mean_value, ecc_value, opts).map(Radians::new)
}

/// Solves the hyperbolic Kepler equation `M = e sinh(F) - F`.
///
/// # Errors
///
/// Returns [`AnomalyError::InvalidEccentricity`] if `ecc ≤ 1`,
/// or [`AnomalyError::NotConverged`] if the solver fails to reach tolerance.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{kepler_hyperbolic, mean_from_hyperbolic, AnomalyOptions};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ecc = Eccentricity::new_unchecked(1.5);
/// let f = kepler_hyperbolic(Radians::new(1.0), ecc, AnomalyOptions::default()).unwrap();
/// assert!((mean_from_hyperbolic(f.value(), ecc).value() - 1.0).abs() < 1e-12);
/// ```
pub fn kepler_hyperbolic(
    mean_anomaly: Radians,
    ecc: Eccentricity,
    opts: AnomalyOptions,
) -> Result<Radians, AnomalyError> {
    let ecc_value = ecc.value();
    let mean_value = mean_anomaly.value();
    if ecc_value <= 1.0 || !ecc_value.is_finite() {
        return Err(AnomalyError::InvalidEccentricity(ecc_value));
    }
    if !mean_value.is_finite() {
        return Err(AnomalyError::InvalidMeanAnomaly(mean_value));
    }
    if mean_value == 0.0 {
        return Ok(Radians::new(0.0));
    }
    let abs_mean = mean_value.abs();
    let mut f_anom = if abs_mean > 50.0 * ecc_value {
        mean_value.signum() * (2.0 * abs_mean / ecc_value).ln()
    } else {
        (mean_value / ecc_value).asinh()
    };
    for _ in 0..opts.max_iter {
        let residual = hyperbolic_residual(f_anom, ecc_value, mean_value);
        if residual.abs() <= opts.tol {
            return Ok(Radians::new(f_anom));
        }
        let fp = ecc_value * f_anom.cosh() - 1.0;
        f_anom -= residual / fp;
        if !f_anom.is_finite() {
            break;
        }
    }
    hyperbolic_bisection(mean_value, ecc_value, opts).map(Radians::new)
}

/// Solves Barker's parabolic equation `M = D + D³/3` analytically.
///
/// # Examples
///
/// ```
/// let d = keplerian::anomaly::kepler_parabolic(0.25);
/// assert!((d + d.powi(3) / 3.0 - 0.25).abs() < 1e-12);
/// ```
#[must_use]
pub fn kepler_parabolic(mean_anomaly_d: f64) -> f64 {
    let a = 1.5 * mean_anomaly_d;
    (a + (a * a + 1.0).sqrt()).cbrt() - ((a * a + 1.0).sqrt() - a).cbrt()
}

/// Converts eccentric anomaly to true anomaly for elliptic motion.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::true_from_eccentric;
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let nu = true_from_eccentric(Radians::new(0.5), Eccentricity::new_unchecked(0.1));
/// assert!(nu.value().is_finite());
/// ```
#[must_use]
pub fn true_from_eccentric(ea: Radians, ecc: Eccentricity) -> Radians {
    let ea_value = ea.value();
    let ecc_value = ecc.value();
    let s = ((1.0 + ecc_value).sqrt() * (0.5 * ea_value).sin())
        .atan2((1.0 - ecc_value).sqrt() * (0.5 * ea_value).cos());
    Radians::new(wrap_two_pi_raw(2.0 * s))
}

/// Converts true anomaly to eccentric anomaly for elliptic motion.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::eccentric_from_true;
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ea = eccentric_from_true(Radians::new(0.5), Eccentricity::new_unchecked(0.1));
/// assert!(ea.value().is_finite());
/// ```
#[must_use]
pub fn eccentric_from_true(nu: Radians, ecc: Eccentricity) -> Radians {
    let nu_value = nu.value();
    let ecc_value = ecc.value();
    let e = 2.0
        * (((1.0 - ecc_value).sqrt() * (0.5 * nu_value).sin())
            .atan2((1.0 + ecc_value).sqrt() * (0.5 * nu_value).cos()));
    Radians::new(wrap_two_pi_raw(e))
}

/// Converts eccentric anomaly to mean anomaly for elliptic motion.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::mean_from_eccentric;
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let m = mean_from_eccentric(Radians::new(0.5), Eccentricity::new_unchecked(0.1));
/// assert!((m.value() - (0.5 - 0.1_f64 * 0.5_f64.sin())).abs() < 1e-12);
/// ```
#[must_use]
pub fn mean_from_eccentric(ea: Radians, ecc: Eccentricity) -> Radians {
    let ea_value = ea.value();
    Radians::new(ea_value - ecc.value() * ea_value.sin())
}

/// Converts mean anomaly to eccentric anomaly for elliptic motion.
///
/// # Errors
///
/// Returns [`AnomalyError`] when the solver fails (see [`kepler_elliptic`]).
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{eccentric_from_mean, mean_from_eccentric, AnomalyOptions};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ecc = Eccentricity::new_unchecked(0.2);
/// let m = Radians::new(1.0);
/// let ea = eccentric_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
/// let m2 = mean_from_eccentric(ea, ecc);
/// assert!((m2.value() - m.value()).abs() < 1e-12);
/// ```
pub fn eccentric_from_mean(
    m: Radians,
    ecc: Eccentricity,
    opts: AnomalyOptions,
) -> Result<Radians, AnomalyError> {
    kepler_elliptic(m, ecc, opts)
}

/// Converts mean anomaly to true anomaly for elliptic motion.
///
/// # Errors
///
/// Returns [`AnomalyError`] when the solver fails (see [`kepler_elliptic`]).
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{true_from_mean, mean_from_true, AnomalyOptions};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ecc = Eccentricity::new_unchecked(0.2);
/// let m = Radians::new(1.0);
/// let nu = true_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
/// let m2 = mean_from_true(nu, ecc);
/// assert!((m2.value() - m.value()).abs() < 1e-12);
/// ```
pub fn true_from_mean(
    m: Radians,
    ecc: Eccentricity,
    opts: AnomalyOptions,
) -> Result<Radians, AnomalyError> {
    Ok(true_from_eccentric(eccentric_from_mean(m, ecc, opts)?, ecc))
}

/// Converts true anomaly to mean anomaly for elliptic motion.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{mean_from_true, true_from_mean, AnomalyOptions};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ecc = Eccentricity::new_unchecked(0.15);
/// let nu = Radians::new(0.8);
/// let m = mean_from_true(nu, ecc);
/// let nu2 = true_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
/// assert!((nu2.value() - nu.value()).abs() < 1e-12);
/// ```
#[must_use]
pub fn mean_from_true(nu: Radians, ecc: Eccentricity) -> Radians {
    mean_from_eccentric(eccentric_from_true(nu, ecc), ecc)
}

/// Converts hyperbolic anomaly to true anomaly.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{true_from_hyperbolic, hyperbolic_from_true};
/// use keplerian::Eccentricity;
/// let ecc = keplerian::Eccentricity::new_unchecked(2.0);
/// let f = 0.5_f64;
/// let nu = true_from_hyperbolic(f, ecc);
/// let f2 = hyperbolic_from_true(nu, ecc);
/// assert!((f2 - f).abs() < 1e-12);
/// ```
#[must_use]
pub fn true_from_hyperbolic(fa: f64, ecc: Eccentricity) -> Radians {
    let ecc_value = ecc.value();
    Radians::new(
        2.0 * (((ecc_value + 1.0).sqrt() * (0.5 * fa).sinh())
            .atan2((ecc_value - 1.0).sqrt() * (0.5 * fa).cosh())),
    )
}

/// Converts true anomaly to hyperbolic anomaly.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{hyperbolic_from_true, true_from_hyperbolic};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ecc = Eccentricity::new_unchecked(2.0);
/// let nu = Radians::new(0.4);
/// let f = hyperbolic_from_true(nu, ecc);
/// let nu2 = true_from_hyperbolic(f, ecc);
/// assert!((nu2.value() - nu.value()).abs() < 1e-12);
/// ```
#[must_use]
pub fn hyperbolic_from_true(nu: Radians, ecc: Eccentricity) -> f64 {
    let ecc_value = ecc.value();
    let t = (0.5 * nu.value()).tan() * ((ecc_value - 1.0) / (ecc_value + 1.0)).sqrt();
    2.0 * t.atanh()
}

/// Converts hyperbolic anomaly to mean anomaly.
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::mean_from_hyperbolic;
/// use keplerian::Eccentricity;
/// let ecc = keplerian::Eccentricity::new_unchecked(1.5);
/// let f = 0.3_f64;
/// let m = mean_from_hyperbolic(f, ecc);
/// assert!((m.value() - (1.5 * 0.3_f64.sinh() - 0.3)).abs() < 1e-12);
/// ```
#[must_use]
pub fn mean_from_hyperbolic(fa: f64, ecc: Eccentricity) -> Radians {
    Radians::new(ecc.value() * fa.sinh() - fa)
}

/// Converts mean anomaly to hyperbolic anomaly.
///
/// # Errors
///
/// Returns [`AnomalyError`] when the solver fails (see [`kepler_hyperbolic`]).
///
/// # Examples
///
/// ```
/// use keplerian::anomaly::{hyperbolic_from_mean, mean_from_hyperbolic, AnomalyOptions};
/// use keplerian::Eccentricity;
/// use qtty::angular::Radians;
/// let ecc = Eccentricity::new_unchecked(1.5);
/// let m = Radians::new(1.0);
/// let f = hyperbolic_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
/// assert!((mean_from_hyperbolic(f, ecc).value() - m.value()).abs() < 1e-12);
/// ```
pub fn hyperbolic_from_mean(
    m: Radians,
    ecc: Eccentricity,
    opts: AnomalyOptions,
) -> Result<f64, AnomalyError> {
    kepler_hyperbolic(m, ecc, opts).map(|value| value.value())
}

/// Wraps a raw angle to `[0, 2π)` for internal numeric use.
#[must_use]
pub(crate) fn wrap_two_pi_raw(x: f64) -> f64 {
    x.rem_euclid(2.0 * PI)
}

#[inline]
fn elliptic_residual(ea: f64, ecc: f64, mean_anomaly: f64) -> f64 {
    ea - ecc * ea.sin() - mean_anomaly
}

fn elliptic_bisection(
    mean_anomaly: f64,
    ecc: f64,
    opts: AnomalyOptions,
) -> Result<f64, AnomalyError> {
    let mut lower = mean_anomaly - PI;
    let mut upper = mean_anomaly + PI;
    let mut residual = elliptic_residual(0.5 * (lower + upper), ecc, mean_anomaly);

    for _ in 0..opts.max_iter {
        let mid = 0.5 * (lower + upper);
        residual = elliptic_residual(mid, ecc, mean_anomaly);
        if residual.abs() <= opts.tol || (upper - lower).abs() <= opts.tol {
            return Ok(mid);
        }
        if residual.is_sign_positive() {
            upper = mid;
        } else {
            lower = mid;
        }
    }

    Err(AnomalyError::NotConverged {
        iterations: opts.max_iter,
        residual: residual.abs(),
    })
}

#[inline]
fn hyperbolic_residual(fa: f64, ecc: f64, mean_anomaly: f64) -> f64 {
    ecc * fa.sinh() - fa - mean_anomaly
}

fn hyperbolic_bisection(
    mean_anomaly: f64,
    ecc: f64,
    opts: AnomalyOptions,
) -> Result<f64, AnomalyError> {
    let sign = mean_anomaly.signum();
    let positive_mean_anomaly = mean_anomaly.abs();
    let mut lower = 0.0;
    let mut upper = (positive_mean_anomaly / ecc).asinh().max(1.0);
    let mut upper_residual = hyperbolic_residual(upper, ecc, positive_mean_anomaly);

    for _ in 0..opts.max_iter {
        if upper_residual.is_sign_positive() || upper_residual == 0.0 {
            break;
        }
        upper *= 2.0;
        upper_residual = hyperbolic_residual(upper, ecc, positive_mean_anomaly);
    }

    if upper_residual.is_sign_negative() {
        return Err(AnomalyError::NotConverged {
            iterations: opts.max_iter,
            residual: upper_residual.abs(),
        });
    }

    let mut residual = upper_residual;
    for _ in 0..opts.max_iter {
        let mid = 0.5 * (lower + upper);
        residual = hyperbolic_residual(mid, ecc, positive_mean_anomaly);
        if residual.abs() <= opts.tol || (upper - lower).abs() <= opts.tol {
            return Ok(sign * mid);
        }
        if residual.is_sign_positive() {
            upper = mid;
        } else {
            lower = mid;
        }
    }

    Err(AnomalyError::NotConverged {
        iterations: opts.max_iter,
        residual: residual.abs(),
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::TAU;

    #[test]
    fn solvers_converge() {
        let zero = Eccentricity::new(0.0).unwrap();
        assert_eq!(
            kepler_elliptic(Radians::new(0.4), zero, AnomalyOptions::default()).unwrap(),
            Radians::new(0.4)
        );
        let elliptic_ecc = Eccentricity::new(0.4).unwrap();
        let e =
            kepler_elliptic(Radians::new(1.0), elliptic_ecc, AnomalyOptions::default()).unwrap();
        assert!((mean_from_eccentric(e, elliptic_ecc).value() - 1.0).abs() < 1e-12);
        let hyperbolic_ecc = Eccentricity::new(1.4).unwrap();
        let f = kepler_hyperbolic(Radians::new(1.0), hyperbolic_ecc, AnomalyOptions::default())
            .unwrap();
        assert!((mean_from_hyperbolic(f.value(), hyperbolic_ecc).value() - 1.0).abs() < 1e-12);
    }

    #[test]
    fn round_trip_elliptic_anomalies() {
        let nu = Radians::new(1.2);
        let ecc = Eccentricity::new(0.3).unwrap();
        let ea = eccentric_from_true(nu, ecc);
        let m = mean_from_eccentric(ea, ecc);
        let nu2 = true_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
        let diff = (nu2.value() - nu.value() + core::f64::consts::PI).rem_euclid(TAU)
            - core::f64::consts::PI;
        assert!(diff.abs() < 1e-12);
    }

    #[test]
    fn round_trip_hyperbolic_anomalies() {
        let nu = Radians::new(0.8);
        let ecc = Eccentricity::new(1.7).unwrap();
        let f = hyperbolic_from_true(nu, ecc);
        let m = mean_from_hyperbolic(f, ecc);
        let f2 = hyperbolic_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
        assert!((f2 - f).abs() < 1e-12);
    }

    #[test]
    fn invalid_eccentricity() {
        assert!(matches!(
            kepler_elliptic(
                Radians::new(0.0),
                Eccentricity::new_unchecked(1.0),
                AnomalyOptions::default()
            ),
            Err(AnomalyError::InvalidEccentricity(_))
        ));
        assert!(matches!(
            kepler_hyperbolic(
                Radians::new(0.0),
                Eccentricity::new_unchecked(0.9),
                AnomalyOptions::default()
            ),
            Err(AnomalyError::InvalidEccentricity(_))
        ));
    }

    #[test]
    fn detects_non_convergence() {
        let err = kepler_elliptic(
            Radians::new(1.0),
            Eccentricity::new(0.9).unwrap(),
            AnomalyOptions {
                max_iter: 0,
                tol: 1e-16,
            },
        )
        .unwrap_err();
        assert!(matches!(err, AnomalyError::NotConverged { .. }));
    }

    #[test]
    fn near_parabolic_hyperbolic_solver_converges() {
        let ecc = Eccentricity::new(1.0000001).unwrap();
        let mean = Radians::new(0.001_f64.to_radians() + 0.01720209895 * 0.5);
        let f = kepler_hyperbolic(
            mean,
            ecc,
            AnomalyOptions {
                max_iter: 100,
                tol: 1e-14,
            },
        )
        .unwrap();

        assert!(f.value().is_finite());
        assert!(hyperbolic_residual(f.value(), ecc.value(), mean.value()).abs() < 1e-13);
    }
}
