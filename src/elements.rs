// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed Keplerian orbital elements and Cartesian conversions.
//!
//! ## Scientific scope
//! This module models the classical six-element representation of a two-body
//! conic and its conversion to and from Cartesian state vectors. The formulas
//! assume a point-mass central field.
//!
//! ## Technical scope
//! [`KeplerianElements`] stores typed distances and angles, plus the semantic
//! [`crate::eccentricity::Eccentricity`] newtype. Conversions operate on
//! [`crate::state::CartesianState`] and typed gravitational parameters.
//!
//! ## References
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

use core::f64::consts::PI;
use core::marker::PhantomData;

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::angular::Radians;
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
use qtty::length::{Kilometer, Kilometers};

use crate::anomaly::wrap_two_pi_raw;
use crate::eccentricity::Eccentricity;
use crate::state::CartesianState;
use crate::vec3::{cross, dot, norm, scale, sub};

pub use crate::eccentricity::ConicRegime;

const EPS: f64 = 1.0e-10;

/// Errors returned by element validation and Cartesian conversion.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
pub enum ConversionError {
    /// Eccentricity is negative, non-finite, or otherwise invalid.
    #[error("invalid eccentricity {0}")]
    InvalidEccentricity(f64),
    /// Inclination is outside `[0, π]` or non-finite.
    #[error("invalid inclination {0}")]
    InvalidInclination(f64),
    /// A named scalar field is not finite.
    #[error("non-finite {field}: {value}")]
    NonFiniteValue {
        /// Field name.
        field: &'static str,
        /// Rejected value.
        value: f64,
    },
    /// The conversion encountered a degenerate geometry.
    #[error("degenerate orbital geometry: {0}")]
    Degenerate(&'static str),
    /// Semi-major axis sign is inconsistent with the eccentricity regime.
    #[error("semi-major axis {sma} km is incoherent with eccentricity {ecc}")]
    IncoherentRegime {
        /// Semi-major axis value.
        sma: f64,
        /// Eccentricity value.
        ecc: f64,
    },
}

/// Classical Keplerian elements with typed distance and angular quantities.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KeplerianElements<F: ReferenceFrame> {
    /// Semi-major axis in kilometres; negative for hyperbolic trajectories.
    pub semi_major_axis: Kilometers,
    /// Dimensionless eccentricity.
    pub eccentricity: Eccentricity,
    /// Inclination in radians.
    pub inclination: Radians,
    /// Right ascension of ascending node in radians.
    pub raan: Radians,
    /// Argument of periapsis in radians.
    pub arg_periapsis: Radians,
    /// True anomaly in radians.
    pub true_anomaly: Radians,
    _frame: PhantomData<F>,
}

impl<F: ReferenceFrame> KeplerianElements<F> {
    /// Creates validated Keplerian elements.
    ///
    /// # Errors
    ///
    /// Returns [`ConversionError`] if any field is non-finite, `eccentricity < 0`,
    /// or `inclination ∉ [0, π]`.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::elements::KeplerianElements;
    /// use keplerian::Eccentricity;
    /// use qtty::angular::Radians;
    /// use qtty::length::Kilometers;
    ///
    /// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
    ///
    /// let el = KeplerianElements::<F>::new(
    ///     Kilometers::new(7000.0),
    ///     Eccentricity::new_unchecked(0.01),
    ///     Radians::new(0.5),
    ///     Radians::new(0.0),
    ///     Radians::new(0.0),
    ///     Radians::new(0.0),
    /// ).unwrap();
    /// assert_eq!(el.semi_major_axis.value(), 7000.0);
    /// ```
    pub fn new(
        semi_major_axis: Kilometers,
        eccentricity: Eccentricity,
        inclination: Radians,
        raan: Radians,
        arg_periapsis: Radians,
        true_anomaly: Radians,
    ) -> Result<Self, ConversionError> {
        validate_finite("semi_major_axis", semi_major_axis.value())?;
        validate_finite("eccentricity", eccentricity.value())?;
        validate_finite("inclination", inclination.value())?;
        validate_finite("raan", raan.value())?;
        validate_finite("arg_periapsis", arg_periapsis.value())?;
        validate_finite("true_anomaly", true_anomaly.value())?;
        if eccentricity.value() < 0.0 {
            return Err(ConversionError::InvalidEccentricity(eccentricity.value()));
        }
        if !(0.0..=PI).contains(&inclination.value()) {
            return Err(ConversionError::InvalidInclination(inclination.value()));
        }
        let e = eccentricity.value();
        let a = semi_major_axis.value();
        if e < 1.0 - EPS && a <= 0.0 {
            return Err(ConversionError::IncoherentRegime { sma: a, ecc: e });
        }
        if e > 1.0 + EPS && a >= 0.0 {
            return Err(ConversionError::IncoherentRegime { sma: a, ecc: e });
        }
        Ok(Self {
            semi_major_axis,
            eccentricity,
            inclination,
            raan: Radians::new(wrap_two_pi_raw(raan.value())),
            arg_periapsis: Radians::new(wrap_two_pi_raw(arg_periapsis.value())),
            true_anomaly: Radians::new(wrap_two_pi_raw(true_anomaly.value())),
            _frame: PhantomData,
        })
    }

    /// Returns the eccentricity-based conic regime.
    #[must_use]
    pub fn conic_kind(&self) -> ConicRegime {
        if self.eccentricity.is_parabolic(EPS) {
            ConicRegime::Parabolic
        } else if self.eccentricity.is_elliptic() {
            ConicRegime::Elliptic
        } else {
            ConicRegime::Hyperbolic
        }
    }

    /// Converts elements to a typed Cartesian state.
    ///
    /// # Errors
    ///
    /// Returns [`ConversionError`] if `mu ≤ 0`, semi-latus rectum `p ≤ 0`, or
    /// the orbit denominator `1 + e cos(ν) ≤ ε`.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::elements::KeplerianElements;
    /// use keplerian::Eccentricity;
    /// use qtty::angular::Radians;
    /// use qtty::length::Kilometers;
    /// use qtty::dynamics::GravitationalParameter;
    ///
    /// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
    /// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
    ///
    /// // Circular orbit: radius equals semi-major axis at every true anomaly.
    /// let el = KeplerianElements::<F>::new(
    ///     Kilometers::new(7000.0),
    ///     Eccentricity::new_unchecked(0.0),
    ///     Radians::new(0.5),
    ///     Radians::new(0.0),
    ///     Radians::new(0.0),
    ///     Radians::new(0.0),
    /// ).unwrap();
    /// let state = el.try_to_cartesian::<C>(GravitationalParameter::new(398600.4418)).unwrap();
    /// let r = state.position();
    /// assert!((r.x().value().hypot(r.y().value()).hypot(r.z().value()) - 7000.0).abs() < 1.0);
    /// ```
    pub fn try_to_cartesian<C: ReferenceCenter<Params = ()>>(
        &self,
        mu: GravitationalParameter,
    ) -> Result<CartesianState<C, F>, ConversionError> {
        let mu_val = mu.value();
        if !mu_val.is_finite() || mu_val <= 0.0 {
            return Err(ConversionError::Degenerate(
                "non-positive or non-finite gravitational parameter",
            ));
        }
        let a = self.semi_major_axis.value();
        let e = self.eccentricity.value();
        let nu = self.true_anomaly.value();
        let p = a * (1.0 - e * e);
        if p <= 0.0 {
            return Err(ConversionError::Degenerate(
                "non-positive semi-latus rectum",
            ));
        }
        let denom = 1.0 + e * nu.cos();
        if denom <= EPS {
            return Err(ConversionError::Degenerate("degenerate orbit denominator"));
        }
        let r = p / denom;
        let root = (mu_val / p).sqrt();
        let r_pqw = [r * nu.cos(), r * nu.sin(), 0.0];
        let v_pqw = [-root * nu.sin(), root * (e + nu.cos()), 0.0];
        let r_ijk = rotate_pqw(
            r_pqw,
            self.raan.value(),
            self.inclination.value(),
            self.arg_periapsis.value(),
        );
        let v_ijk = rotate_pqw(
            v_pqw,
            self.raan.value(),
            self.inclination.value(),
            self.arg_periapsis.value(),
        );
        Ok(CartesianState::new(
            Position::<C, F, Kilometer>::new(r_ijk[0], r_ijk[1], r_ijk[2]),
            Velocity::<F, KmPerSecond>::new(v_ijk[0], v_ijk[1], v_ijk[2]),
        ))
    }

    /// Converts a typed Cartesian state into classical Keplerian elements.
    ///
    /// # Errors
    ///
    /// Returns [`ConversionError`] for degenerate geometry (zero r, zero h,
    /// parabolic orbit) or non-finite inputs.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::elements::KeplerianElements;
    /// use keplerian::state::CartesianState;
    /// use keplerian::Eccentricity;
    /// use qtty::angular::Radians;
    /// use qtty::length::Kilometers;
    /// use qtty::dynamics::GravitationalParameter;
    /// use affn::cartesian::{Position, Velocity};
    ///
    /// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
    /// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
    ///
    /// let pos = Position::new(7000.0, 0.0, 0.0);
    /// let vel = Velocity::new(0.0, 7.546, 0.0);
    /// let state = CartesianState::<C, F>::new(pos, vel);
    /// let mu = GravitationalParameter::new(398600.4418);
    /// let el = KeplerianElements::<F>::from_cartesian(&state, mu).unwrap();
    /// assert!((el.semi_major_axis.value() - 7000.0).abs() < 10.0);
    /// ```
    pub fn from_cartesian<C: ReferenceCenter>(
        state: &CartesianState<C, F>,
        mu: GravitationalParameter,
    ) -> Result<Self, ConversionError> {
        let r = vec3_from_pos(state.position());
        let v = vec3_from_vel(state.velocity());
        let mu = mu.value();
        validate_finite("mu", mu)?;
        if mu <= 0.0 {
            return Err(ConversionError::Degenerate(
                "non-positive gravitational parameter",
            ));
        }
        let rmag = norm(r);
        let vmag = norm(v);
        if rmag <= EPS {
            return Err(ConversionError::Degenerate("zero position"));
        }
        let h = cross(r, v);
        let hmag = norm(h);
        if hmag <= EPS {
            return Err(ConversionError::Degenerate("zero angular momentum"));
        }
        let n = [-h[1], h[0], 0.0];
        let nmag = norm(n);
        let e_vec = sub(scale(cross(v, h), 1.0 / mu), scale(r, 1.0 / rmag));
        let ecc_value = norm(e_vec);
        let energy = 0.5 * vmag * vmag - mu / rmag;
        if energy.abs() <= EPS {
            return Err(ConversionError::Degenerate("parabolic orbit"));
        }
        let a = -mu / (2.0 * energy);
        let inc = (h[2] / hmag).clamp(-1.0, 1.0).acos();
        let raan = if nmag > EPS {
            wrap_two_pi_raw(n[1].atan2(n[0]))
        } else {
            0.0
        };
        let hhat = scale(h, 1.0 / hmag);
        let argp = if nmag > EPS && ecc_value > EPS {
            wrap_two_pi_raw(dot(cross(n, e_vec), hhat).atan2(dot(n, e_vec)))
        } else if ecc_value > EPS {
            // Equatorial orbit: periapsis direction is e_vec projected to the equatorial plane.
            wrap_two_pi_raw(e_vec[1].atan2(e_vec[0]))
        } else {
            0.0
        };
        let nu = if ecc_value > EPS {
            wrap_two_pi_raw(
                (dot(cross(e_vec, r), hhat) / (ecc_value * rmag))
                    .atan2(dot(e_vec, r) / (ecc_value * rmag)),
            )
        } else if nmag > EPS {
            wrap_two_pi_raw(
                (dot(cross(n, r), hhat) / (nmag * rmag)).atan2(dot(n, r) / (nmag * rmag)),
            )
        } else {
            wrap_two_pi_raw(r[1].atan2(r[0]))
        };
        let eccentricity =
            Eccentricity::new(ecc_value).ok_or(ConversionError::InvalidEccentricity(ecc_value))?;
        Self::new(
            Kilometers::new(a),
            eccentricity,
            Radians::new(inc),
            Radians::new(raan),
            Radians::new(argp),
            Radians::new(nu),
        )
    }
}

fn validate_finite(field: &'static str, value: f64) -> Result<(), ConversionError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(ConversionError::NonFiniteValue { field, value })
    }
}

fn vec3_from_pos<C: ReferenceCenter, F: ReferenceFrame>(p: &Position<C, F, Kilometer>) -> [f64; 3] {
    [p.x().value(), p.y().value(), p.z().value()]
}

fn vec3_from_vel<F: ReferenceFrame>(v: &Velocity<F, KmPerSecond>) -> [f64; 3] {
    [v.x().value(), v.y().value(), v.z().value()]
}

fn rotate_pqw(v: [f64; 3], raan: f64, inc: f64, argp: f64) -> [f64; 3] {
    let (co, so) = (raan.cos(), raan.sin());
    let (ci, si) = (inc.cos(), inc.sin());
    let (cw, sw) = (argp.cos(), argp.sin());
    let r11 = co * cw - so * sw * ci;
    let r12 = -co * sw - so * cw * ci;
    let r21 = so * cw + co * sw * ci;
    let r22 = -so * sw + co * cw * ci;
    let r31 = sw * si;
    let r32 = cw * si;
    [
        r11 * v[0] + r12 * v[1],
        r21 * v[0] + r22 * v[1],
        r31 * v[0] + r32 * v[1],
    ]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Copy, Clone)]
    struct Center;
    impl ReferenceCenter for Center {
        type Params = ();
        fn center_name() -> &'static str {
            "C"
        }
    }
    #[derive(Debug, Copy, Clone)]
    struct Frame;
    impl ReferenceFrame for Frame {
        fn frame_name() -> &'static str {
            "F"
        }
    }

    #[test]
    fn circular_equatorial_case() {
        let state = CartesianState::<Center, Frame>::new(
            Position::new(7000.0, 0.0, 0.0),
            Velocity::new(0.0, (398600.4418_f64 / 7000.0).sqrt(), 0.0),
        );
        let el =
            KeplerianElements::from_cartesian(&state, GravitationalParameter::new(398600.4418))
                .unwrap();
        assert!(el.eccentricity.value() < 1e-12);
        assert!(el.inclination.value().abs() < 1e-12);
        assert!((el.semi_major_axis.value() - 7000.0).abs() < 1e-8);
    }

    #[test]
    fn round_trip_moderate_eccentricity() {
        let el = KeplerianElements::<Frame>::new(
            Kilometers::new(12000.0),
            Eccentricity::new(0.2).unwrap(),
            Radians::new(0.4),
            Radians::new(0.3),
            Radians::new(0.7),
            Radians::new(1.1),
        )
        .unwrap();
        let st = el
            .try_to_cartesian::<Center>(GravitationalParameter::new(398600.4418))
            .unwrap();
        let back = KeplerianElements::from_cartesian(&st, GravitationalParameter::new(398600.4418))
            .unwrap();
        assert!((back.semi_major_axis.value() - el.semi_major_axis.value()).abs() < 1e-8);
        assert!((back.eccentricity.value() - el.eccentricity.value()).abs() < 1e-12);
    }
}
