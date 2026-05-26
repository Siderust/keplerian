// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Analytic transfer and invariant helpers for two-body motion.
//!
//! ## Scientific scope
//! This module provides textbook two-body helpers such as Hohmann transfer
//! estimates, vis-viva speed, escape speed, and the classical energy and
//! angular-momentum invariants.
//!
//! ## Technical scope
//! Public functions consume typed states and typed quantities. Returned orbital
//! invariants use dedicated `qtty::dynamics` aliases for specific energy and
//! specific angular momentum.
//!
//! ## References
//! - Hohmann, W. (1925). *Die Erreichbarkeit der Himmelskörper*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::dynamics::{
    GravitationalParameter, KmPerSeconds, SpecificAngularMomentum, SpecificOrbitalEnergy,
};
use qtty::length::Kilometers;
use qtty::Second;

use crate::problem::KeplerProblem;
use crate::state::CartesianState;
use crate::vec3::{cross, norm};

/// Errors returned by fallible transfer helpers.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum TransferError {
    /// Gravitational parameter is not strictly positive or not finite.
    #[error("invalid gravitational parameter: {0}")]
    InvalidGravitationalParameter(f64),
    /// A radius argument is not strictly positive or not finite.
    #[error("invalid {0} radius: {1}")]
    InvalidRadius(&'static str, f64),
    /// Semi-major axis is not strictly positive or not finite.
    #[error("invalid semi-major axis: {0}")]
    InvalidSemiMajorAxis(f64),
    /// Computation produced a non-finite result.
    #[error("non-finite result")]
    NonFiniteResult,
}

/// Result of an ideal coplanar Hohmann transfer between circular orbits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HohmannResult {
    /// Departure impulse magnitude in km/s.
    pub dv_departure: KmPerSeconds,
    /// Arrival circularization impulse magnitude in km/s.
    pub dv_arrival: KmPerSeconds,
    /// Sum of departure and arrival impulse magnitudes in km/s.
    pub total: KmPerSeconds,
    /// Half-period of the transfer ellipse in seconds.
    pub transfer_time: Second,
}

/// Returns the elliptic orbital period for `a > 0`, or `None` for non-elliptic `a`.
///
/// # Examples
///
/// ```
/// use keplerian::problem::KeplerProblem;
/// use keplerian::transfer::orbital_period;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
/// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
///
/// let mu = GravitationalParameter::new(398600.4418);
/// let p = KeplerProblem::<C, F>::new(mu);
/// let t = orbital_period(&p, Kilometers::new(6678.0)).unwrap();
/// assert!((t.value() - 5431.0).abs() < 10.0);
/// ```
#[must_use]
pub fn orbital_period<C: ReferenceCenter, F: ReferenceFrame>(
    problem: &KeplerProblem<C, F>,
    semi_major_axis: Kilometers,
) -> Option<Second> {
    let a = semi_major_axis.value();
    (a > 0.0).then(|| {
        Second::new(2.0 * core::f64::consts::PI * (a * a * a / problem.mu().value()).sqrt())
    })
}

/// Computes specific orbital energy `ε = v²/2 - μ/r` in km²/s².
///
/// # Examples
///
/// ```
/// use keplerian::transfer::specific_orbital_energy;
/// use keplerian::state::CartesianState;
/// use qtty::dynamics::{GravitationalParameter, KmPerSecond};
/// use qtty::length::Kilometer;
/// use affn::cartesian::{Position, Velocity};
///
/// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
/// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
///
/// let pos = Position::new(6678.0, 0.0, 0.0);
/// let vel = Velocity::new(0.0, 7.72, 0.0);
/// let state = CartesianState::<C, F>::new(pos, vel);
/// let mu = GravitationalParameter::new(398600.4418);
/// let e = specific_orbital_energy(&state, mu);
/// assert!(e.value() < 0.0); // bound orbit
/// ```
#[must_use]
pub fn specific_orbital_energy<C: ReferenceCenter, F: ReferenceFrame>(
    state: &CartesianState<C, F>,
    mu: GravitationalParameter,
) -> SpecificOrbitalEnergy {
    let r = [
        state.position().x().value(),
        state.position().y().value(),
        state.position().z().value(),
    ];
    let v = [
        state.velocity().x().value(),
        state.velocity().y().value(),
        state.velocity().z().value(),
    ];
    SpecificOrbitalEnergy::new(
        0.5 * (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) - mu.value() / norm(r),
    )
}

/// Computes specific angular momentum magnitude `|r × v|` in km²/s.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::specific_angular_momentum;
/// use keplerian::state::CartesianState;
/// use affn::cartesian::{Position, Velocity};
///
/// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
/// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
///
/// let pos = Position::new(6678.0, 0.0, 0.0);
/// let vel = Velocity::new(0.0, 7.72, 0.0);
/// let state = CartesianState::<C, F>::new(pos, vel);
/// let h = specific_angular_momentum(&state);
/// assert!((h.value() - 6678.0 * 7.72).abs() < 0.01);
/// ```
#[must_use]
pub fn specific_angular_momentum<C: ReferenceCenter, F: ReferenceFrame>(
    state: &CartesianState<C, F>,
) -> SpecificAngularMomentum {
    let r = [
        state.position().x().value(),
        state.position().y().value(),
        state.position().z().value(),
    ];
    let v = [
        state.velocity().x().value(),
        state.velocity().y().value(),
        state.velocity().z().value(),
    ];
    SpecificAngularMomentum::new(norm(cross(r, v)))
}

/// Computes ideal Hohmann transfer impulses between coplanar circular orbits.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::hohmann_delta_v;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// let mu = GravitationalParameter::new(398600.4418);
/// let h = hohmann_delta_v(mu, Kilometers::new(6678.0), Kilometers::new(42164.0));
/// assert!((h.total.value() - 3.91).abs() < 0.1);
/// ```
#[must_use]
pub fn hohmann_delta_v(
    mu: GravitationalParameter,
    r1: Kilometers,
    r2: Kilometers,
) -> HohmannResult {
    let mu = mu.value();
    let r1 = r1.value();
    let r2 = r2.value();
    let a_t = 0.5 * (r1 + r2);
    let v1 = (mu / r1).sqrt();
    let v2 = (mu / r2).sqrt();
    let vt1 = (mu * (2.0 / r1 - 1.0 / a_t)).sqrt();
    let vt2 = (mu * (2.0 / r2 - 1.0 / a_t)).sqrt();
    let d1 = (vt1 - v1).abs();
    let d2 = (v2 - vt2).abs();
    HohmannResult {
        dv_departure: KmPerSeconds::new(d1),
        dv_arrival: KmPerSeconds::new(d2),
        total: KmPerSeconds::new(d1 + d2),
        transfer_time: Second::new(core::f64::consts::PI * (a_t * a_t * a_t / mu).sqrt()),
    }
}

/// Computes vis-viva speed `sqrt(μ(2/r - 1/a))` in km/s.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::vis_viva_speed;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// let mu = GravitationalParameter::new(398600.4418);
/// let v = vis_viva_speed(mu, Kilometers::new(6678.0), Kilometers::new(6678.0));
/// assert!((v.value() - 7.726).abs() < 0.01);
/// ```
#[must_use]
pub fn vis_viva_speed(mu: GravitationalParameter, r: Kilometers, a: Kilometers) -> KmPerSeconds {
    KmPerSeconds::new((mu.value() * (2.0 / r.value() - 1.0 / a.value())).sqrt())
}

/// Computes two-body escape speed `sqrt(2μ/r)` in km/s.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::{escape_speed, vis_viva_speed};
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// let mu = GravitationalParameter::new(398600.4418);
/// let r = Kilometers::new(6678.0);
/// let v_esc = escape_speed(mu, r);
/// assert!(v_esc.value() > vis_viva_speed(mu, r, r).value());
/// ```
#[must_use]
pub fn escape_speed(mu: GravitationalParameter, r: Kilometers) -> KmPerSeconds {
    KmPerSeconds::new((2.0 * mu.value() / r.value()).sqrt())
}

/// Fallible version of [`orbital_period`].
///
/// # Errors
///
/// Returns [`TransferError`] if `mu ≤ 0`, `a ≤ 0`, or the result is not finite.
///
/// # Examples
///
/// ```
/// use keplerian::problem::KeplerProblem;
/// use keplerian::transfer::try_orbital_period;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
/// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
///
/// let p = KeplerProblem::<C, F>::new(GravitationalParameter::new(398600.4418));
/// let t = try_orbital_period(&p, Kilometers::new(7000.0)).unwrap();
/// assert!((t.value() - 5840.0).abs() < 200.0);
/// ```
pub fn try_orbital_period<C: ReferenceCenter, F: ReferenceFrame>(
    problem: &KeplerProblem<C, F>,
    semi_major_axis: Kilometers,
) -> Result<Second, TransferError> {
    let mu = problem.mu().value();
    if !mu.is_finite() || mu <= 0.0 {
        return Err(TransferError::InvalidGravitationalParameter(mu));
    }
    let a = semi_major_axis.value();
    if !a.is_finite() || a <= 0.0 {
        return Err(TransferError::InvalidSemiMajorAxis(a));
    }
    let t = 2.0 * core::f64::consts::PI * (a * a * a / mu).sqrt();
    if !t.is_finite() {
        return Err(TransferError::NonFiniteResult);
    }
    Ok(Second::new(t))
}

/// Fallible version of [`hohmann_delta_v`].
///
/// # Errors
///
/// Returns [`TransferError`] if inputs are non-positive or the result is not finite.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::try_hohmann_delta_v;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// let h = try_hohmann_delta_v(
///     GravitationalParameter::new(398600.4418),
///     Kilometers::new(6678.0),
///     Kilometers::new(42164.0),
/// ).unwrap();
/// assert!((h.total.value() - 3.91).abs() < 0.1);
/// ```
pub fn try_hohmann_delta_v(
    mu: GravitationalParameter,
    r1: Kilometers,
    r2: Kilometers,
) -> Result<HohmannResult, TransferError> {
    let mu_val = mu.value();
    if !mu_val.is_finite() || mu_val <= 0.0 {
        return Err(TransferError::InvalidGravitationalParameter(mu_val));
    }
    let r1_val = r1.value();
    if !r1_val.is_finite() || r1_val <= 0.0 {
        return Err(TransferError::InvalidRadius("r1", r1_val));
    }
    let r2_val = r2.value();
    if !r2_val.is_finite() || r2_val <= 0.0 {
        return Err(TransferError::InvalidRadius("r2", r2_val));
    }
    let result = hohmann_delta_v(mu, r1, r2);
    if !result.total.value().is_finite() {
        return Err(TransferError::NonFiniteResult);
    }
    Ok(result)
}

/// Fallible version of [`vis_viva_speed`].
///
/// # Errors
///
/// Returns [`TransferError`] if inputs are non-positive or the result is not finite.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::try_vis_viva_speed;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// let v = try_vis_viva_speed(
///     GravitationalParameter::new(398600.4418),
///     Kilometers::new(7000.0),
///     Kilometers::new(7000.0),
/// ).unwrap();
/// assert!((v.value() - 7.546).abs() < 0.01);
/// ```
pub fn try_vis_viva_speed(
    mu: GravitationalParameter,
    r: Kilometers,
    a: Kilometers,
) -> Result<KmPerSeconds, TransferError> {
    let mu_val = mu.value();
    if !mu_val.is_finite() || mu_val <= 0.0 {
        return Err(TransferError::InvalidGravitationalParameter(mu_val));
    }
    let r_val = r.value();
    if !r_val.is_finite() || r_val <= 0.0 {
        return Err(TransferError::InvalidRadius("r", r_val));
    }
    let a_val = a.value();
    if !a_val.is_finite() || a_val == 0.0 {
        return Err(TransferError::InvalidSemiMajorAxis(a_val));
    }
    let v = vis_viva_speed(mu, r, a);
    if !v.value().is_finite() {
        return Err(TransferError::NonFiniteResult);
    }
    Ok(v)
}

/// Fallible version of [`escape_speed`].
///
/// # Errors
///
/// Returns [`TransferError`] if inputs are non-positive or the result is not finite.
///
/// # Examples
///
/// ```
/// use keplerian::transfer::try_escape_speed;
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::length::Kilometers;
///
/// let v = try_escape_speed(
///     GravitationalParameter::new(398600.4418),
///     Kilometers::new(6378.0),
/// ).unwrap();
/// assert!((v.value() - 11.18).abs() < 0.1);
/// ```
pub fn try_escape_speed(
    mu: GravitationalParameter,
    r: Kilometers,
) -> Result<KmPerSeconds, TransferError> {
    let mu_val = mu.value();
    if !mu_val.is_finite() || mu_val <= 0.0 {
        return Err(TransferError::InvalidGravitationalParameter(mu_val));
    }
    let r_val = r.value();
    if !r_val.is_finite() || r_val <= 0.0 {
        return Err(TransferError::InvalidRadius("r", r_val));
    }
    let v = escape_speed(mu, r);
    if !v.value().is_finite() {
        return Err(TransferError::NonFiniteResult);
    }
    Ok(v)
}

#[cfg(test)]
mod tests {
    use super::*;
    use affn::cartesian::{Position, Velocity};
    use qtty::dynamics::KmPerSecond;
    use qtty::length::Kilometer;

    #[derive(Debug, Copy, Clone)]
    struct C;
    impl ReferenceCenter for C {
        type Params = ();
        fn center_name() -> &'static str {
            "C"
        }
    }
    #[derive(Debug, Copy, Clone)]
    struct F;
    impl ReferenceFrame for F {
        fn frame_name() -> &'static str {
            "F"
        }
    }

    #[test]
    fn helper_values_are_reasonable() {
        let mu = GravitationalParameter::new(398600.4418);
        let h = hohmann_delta_v(mu, Kilometers::new(6678.0), Kilometers::new(42164.0));
        assert!((h.total.value() - 3.91).abs() < 0.2);
        assert!((escape_speed(mu, Kilometers::new(6378.0)).value() - 11.18).abs() < 0.1);
        assert!(
            (vis_viva_speed(mu, Kilometers::new(7000.0), Kilometers::new(7000.0)).value() - 7.546)
                .abs()
                < 0.01
        );
    }

    #[test]
    fn invariants_compute() {
        let s = CartesianState::<C, F>::new(
            Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0),
            Velocity::<F, KmPerSecond>::new(0.0, 7.5, 0.0),
        );
        assert!(
            specific_orbital_energy(&s, GravitationalParameter::new(398600.4418)).value() < 0.0
        );
        assert!((specific_angular_momentum(&s).value() - 52500.0).abs() < 1e-12);
    }
}
