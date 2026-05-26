// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed central-force two-body propagation context.
//!
//! ## Scientific scope
//! This module propagates Cartesian states under the unperturbed two-body
//! Kepler problem. It supports elliptic and hyperbolic motion and rejects the
//! parabolic limit.
//!
//! ## Technical scope
//! [`KeplerProblem`] stores a typed gravitational parameter and propagates
//! [`crate::state::CartesianState`] by a typed elapsed [`qtty::Second`].
//!
//! ## References
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

use core::marker::PhantomData;

use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::dynamics::GravitationalParameter;
use qtty::Second;

use crate::anomaly::{
    hyperbolic_from_mean, hyperbolic_from_true, mean_from_hyperbolic, mean_from_true,
    true_from_hyperbolic, true_from_mean, AnomalyError, AnomalyOptions, MeanAnomaly, TrueAnomaly,
};
use crate::elements::{ConicRegime, ConversionError, KeplerianElements};
use crate::state::CartesianState;

/// Errors returned by two-body propagation.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum PropagationError {
    /// Eccentricity is invalid for propagation.
    #[error("invalid eccentricity {0}")]
    InvalidEccentricity(f64),
    /// Parabolic propagation is not currently supported.
    #[error("parabolic propagation is unsupported")]
    ParabolicUnsupported,
    /// Anomaly solver failed.
    #[error(transparent)]
    AnomalyError(#[from] AnomalyError),
    /// Cartesian/element conversion failed.
    #[error(transparent)]
    ConversionError(#[from] ConversionError),
}

/// A central-force Kepler problem for a center/frame pair.
#[derive(Debug, Clone, Copy)]
pub struct KeplerProblem<C: ReferenceCenter, F: ReferenceFrame> {
    mu: GravitationalParameter,
    _marker: PhantomData<(C, F)>,
}

impl<C: ReferenceCenter, F: ReferenceFrame> KeplerProblem<C, F> {
    /// Creates a new central-force problem.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::problem::KeplerProblem;
    /// use qtty::dynamics::GravitationalParameter;
    ///
    /// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
    /// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
    ///
    /// let mu = GravitationalParameter::new(398600.4418);
    /// let p = KeplerProblem::<C, F>::new(mu);
    /// assert_eq!(p.mu().value(), 398600.4418);
    /// ```
    #[must_use]
    pub fn new(mu: GravitationalParameter) -> Self {
        Self {
            mu,
            _marker: PhantomData,
        }
    }

    /// Returns the problem gravitational parameter.
    #[must_use]
    pub fn mu(&self) -> GravitationalParameter {
        self.mu
    }

    /// Propagates a state by a typed elapsed duration using Kepler's equation.
    ///
    /// # Errors
    ///
    /// Returns [`PropagationError`] for parabolic orbits, invalid eccentricity,
    /// anomaly solver failures, or degenerate Cartesian/element conversions.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::problem::KeplerProblem;
    /// use keplerian::state::CartesianState;
    /// use qtty::dynamics::GravitationalParameter;
    /// use qtty::Second;
    /// use affn::cartesian::{Position, Velocity};
    ///
    /// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
    /// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
    ///
    /// let mu = GravitationalParameter::new(398600.4418);
    /// // Exact circular speed at 7000 km ensures e ≈ 0 and a clean round-trip.
    /// let v_circ = (398600.4418_f64 / 7000.0_f64).sqrt();
    /// let pos = Position::new(7000.0, 0.0, 0.0);
    /// let vel = Velocity::new(0.0, v_circ, 0.0);
    /// let state = CartesianState::<C, F>::new(pos, vel);
    /// let problem = KeplerProblem::<C, F>::new(mu);
    /// let s2 = problem.propagate(&state, Second::new(0.0)).unwrap();
    /// let p = s2.position();
    /// assert!((p.x().value() - 7000.0).abs() < 1.0);
    /// ```
    pub fn propagate(
        &self,
        state: &CartesianState<C, F>,
        dt: Second,
    ) -> Result<CartesianState<C, F>, PropagationError>
    where
        C: ReferenceCenter<Params = ()>,
    {
        let el = KeplerianElements::<F>::from_cartesian(state, self.mu)?;
        let a = el.semi_major_axis.value();
        let ecc_value = el.eccentricity.value();
        if !ecc_value.is_finite() || ecc_value < 0.0 {
            return Err(PropagationError::InvalidEccentricity(ecc_value));
        }
        let true_anomaly = match el.conic_kind() {
            ConicRegime::Elliptic => {
                let n = (self.mu.value() / (a * a * a)).sqrt();
                let m0 = mean_from_true(TrueAnomaly::new(el.true_anomaly), el.eccentricity);
                let m = MeanAnomaly::from_value(m0.value() + n * dt.value());
                true_from_mean(m, el.eccentricity, AnomalyOptions::default())?
            }
            ConicRegime::Hyperbolic => {
                let n = (-self.mu.value() / (a * a * a)).sqrt();
                let f0 = hyperbolic_from_true(TrueAnomaly::new(el.true_anomaly), el.eccentricity);
                let m0 = mean_from_hyperbolic(f0, el.eccentricity);
                let m = MeanAnomaly::from_value(m0.value() + n * dt.value());
                true_from_hyperbolic(
                    hyperbolic_from_mean(m, el.eccentricity, AnomalyOptions::default())?,
                    el.eccentricity,
                )
            }
            ConicRegime::Parabolic => return Err(PropagationError::ParabolicUnsupported),
        };
        let next = KeplerianElements::<F>::new(
            el.semi_major_axis,
            el.eccentricity,
            el.inclination,
            el.raan,
            el.arg_periapsis,
            true_anomaly.radians(),
        )?;
        Ok(next.try_to_cartesian::<C>(self.mu)?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use affn::cartesian::{Position, Velocity};
    use qtty::dynamics::KmPerSecond;
    use qtty::length::Kilometer;

    use crate::state::CartesianState;
    use crate::transfer::{specific_angular_momentum, specific_orbital_energy};

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
    fn circular_returns_after_period() {
        let mu = GravitationalParameter::new(398600.4418);
        let r = 7000.0;
        let state = CartesianState::<C, F>::new(
            Position::<C, F, Kilometer>::new(r, 0.0, 0.0),
            Velocity::<F, KmPerSecond>::new(0.0, (mu.value() / r).sqrt(), 0.0),
        );
        let period = 2.0 * core::f64::consts::PI * (r * r * r / mu.value()).sqrt();
        let out = KeplerProblem::<C, F>::new(mu)
            .propagate(&state, Second::new(period))
            .unwrap();
        assert!((out.position().x().value() - r).abs() < 1e-6);
        assert!(out.position().y().value().abs() < 1e-5);
    }

    #[test]
    fn invariants_are_conserved() {
        let mu = GravitationalParameter::new(398600.4418);
        let state = CartesianState::<C, F>::new(
            Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0),
            Velocity::<F, KmPerSecond>::new(0.0, 7.2, 1.0),
        );
        let out = KeplerProblem::<C, F>::new(mu)
            .propagate(&state, Second::new(1000.0))
            .unwrap();
        assert!(
            (specific_orbital_energy(&state, mu) - specific_orbital_energy(&out, mu))
                .value()
                .abs()
                < 1e-8
        );
        assert!(
            (specific_angular_momentum(&state) - specific_angular_momentum(&out))
                .value()
                .abs()
                < 1e-8
        );
    }
}
