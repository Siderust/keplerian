// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Domain-semantic eccentricity scalar and conic regime classifier.
//!
//! ## Scientific scope
//! Eccentricity classifies conic sections in central-force motion. This module
//! models the scalar `e` itself and the resulting [`ConicRegime`]; it does not
//! attach any frame, epoch, or body semantics.
//!
//! ## Technical scope
//! [`Eccentricity`] wraps a raw `f64` so public APIs can distinguish orbital
//! eccentricity from unrelated dimensionless diagnostics such as tolerances or
//! residuals. [`ConicRegime`] is derived via [`Eccentricity::classify`].
//!
//! ## References
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

/// Conic regime classified by eccentricity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ConicRegime {
    /// Bounded ellipse, `0 ≤ e < 1 − ε`.
    Elliptic,
    /// Parabola, `|e − 1| ≤ ε`.
    Parabolic,
    /// Unbounded hyperbola, `e > 1 + ε`.
    Hyperbolic,
}

/// Errors returned by fallible eccentricity constructors.
#[derive(Debug, Clone, Copy, PartialEq, thiserror::Error)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EccentricityError {
    /// The value is negative.
    #[error("eccentricity must be non-negative, got {0}")]
    Negative(f64),
    /// The value is not finite.
    #[error("eccentricity must be finite, got {0}")]
    NotFinite(f64),
    /// An elliptic constructor received `e ≥ 1`.
    #[error("elliptic eccentricity requires e < 1, got {0}")]
    NotElliptic(f64),
    /// A hyperbolic constructor received `e ≤ 1`.
    #[error("hyperbolic eccentricity requires e > 1, got {0}")]
    NotHyperbolic(f64),
}

/// Dimensionless eccentricity scalar (domain-semantic newtype).
///
/// Prevents accidentally mixing eccentricity with tolerances, residuals,
/// or other dimensionless scalars. Valid range is `[0, ∞)`.
///
/// # Examples
///
/// ```
/// use keplerian::Eccentricity;
/// let e = Eccentricity::new(0.3).unwrap();
/// assert!(e.is_elliptic());
/// assert!(!e.is_hyperbolic());
/// ```
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct Eccentricity(f64);

impl Eccentricity {
    /// Creates a new eccentricity, returning `None` if `e` is negative or non-finite.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::Eccentricity;
    /// assert!(Eccentricity::new(0.5).is_some());
    /// assert!(Eccentricity::new(-0.1).is_none());
    /// ```
    #[must_use]
    pub fn new(e: f64) -> Option<Self> {
        (e.is_finite() && e >= 0.0).then_some(Self(e))
    }

    /// Creates a new eccentricity, returning [`EccentricityError`] if invalid.
    ///
    /// # Errors
    ///
    /// Returns [`EccentricityError::NotFinite`] for NaN/±∞, or
    /// [`EccentricityError::Negative`] for `e < 0`.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::eccentricity::{Eccentricity, EccentricityError};
    /// assert!(Eccentricity::try_new(0.5).is_ok());
    /// assert!(matches!(Eccentricity::try_new(-0.1), Err(EccentricityError::Negative(_))));
    /// ```
    pub fn try_new(e: f64) -> Result<Self, EccentricityError> {
        if !e.is_finite() {
            return Err(EccentricityError::NotFinite(e));
        }
        if e < 0.0 {
            return Err(EccentricityError::Negative(e));
        }
        Ok(Self(e))
    }

    /// Creates an elliptic eccentricity (`0 ≤ e < 1`).
    ///
    /// # Errors
    ///
    /// Returns [`EccentricityError`] if `e < 0`, non-finite, or `e ≥ 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::eccentricity::{Eccentricity, EccentricityError};
    /// assert!(Eccentricity::new_elliptic(0.5).is_ok());
    /// assert!(matches!(Eccentricity::new_elliptic(1.5), Err(EccentricityError::NotElliptic(_))));
    /// ```
    pub fn new_elliptic(e: f64) -> Result<Self, EccentricityError> {
        let ec = Self::try_new(e)?;
        if e >= 1.0 {
            return Err(EccentricityError::NotElliptic(e));
        }
        Ok(ec)
    }

    /// Creates a hyperbolic eccentricity (`e > 1`).
    ///
    /// # Errors
    ///
    /// Returns [`EccentricityError`] if `e < 0`, non-finite, or `e ≤ 1`.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::eccentricity::{Eccentricity, EccentricityError};
    /// assert!(Eccentricity::new_hyperbolic(1.5).is_ok());
    /// assert!(matches!(Eccentricity::new_hyperbolic(0.5), Err(EccentricityError::NotHyperbolic(_))));
    /// ```
    pub fn new_hyperbolic(e: f64) -> Result<Self, EccentricityError> {
        let ec = Self::try_new(e)?;
        if e <= 1.0 {
            return Err(EccentricityError::NotHyperbolic(e));
        }
        Ok(ec)
    }

    /// Returns the parabolic eccentricity `e = 1.0` exactly.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::eccentricity::Eccentricity;
    /// assert_eq!(Eccentricity::parabolic().value(), 1.0);
    /// ```
    #[must_use]
    pub const fn parabolic() -> Self {
        Self(1.0)
    }

    /// Creates a new eccentricity without checking.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::Eccentricity;
    /// let e = Eccentricity::new_unchecked(0.7);
    /// assert_eq!(e.value(), 0.7);
    /// ```
    #[must_use]
    pub const fn new_unchecked(e: f64) -> Self {
        Self(e)
    }

    /// Returns the raw dimensionless value.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::Eccentricity;
    /// assert_eq!(Eccentricity::new_unchecked(0.1).value(), 0.1);
    /// ```
    #[must_use]
    pub const fn value(self) -> f64 {
        self.0
    }

    /// Classifies this eccentricity into a [`ConicRegime`] using tolerance `eps`.
    ///
    /// `|e − 1| ≤ eps` → [`ConicRegime::Parabolic`]; `e < 1 − eps` → Elliptic;
    /// `e > 1 + eps` → Hyperbolic.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::eccentricity::{ConicRegime, Eccentricity};
    /// assert_eq!(Eccentricity::new_unchecked(0.5).classify(1e-10), ConicRegime::Elliptic);
    /// assert_eq!(Eccentricity::new_unchecked(1.0).classify(1e-10), ConicRegime::Parabolic);
    /// assert_eq!(Eccentricity::new_unchecked(1.5).classify(1e-10), ConicRegime::Hyperbolic);
    /// ```
    #[must_use]
    pub fn classify(self, eps: f64) -> ConicRegime {
        if (self.0 - 1.0).abs() <= eps {
            ConicRegime::Parabolic
        } else if self.0 < 1.0 {
            ConicRegime::Elliptic
        } else {
            ConicRegime::Hyperbolic
        }
    }

    /// Returns `true` for an elliptic orbit (`0 ≤ e < 1`).
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::Eccentricity;
    /// assert!(Eccentricity::new_unchecked(0.5).is_elliptic());
    /// assert!(!Eccentricity::new_unchecked(1.5).is_elliptic());
    /// ```
    #[must_use]
    pub fn is_elliptic(self) -> bool {
        self.0 < 1.0
    }

    /// Returns `true` for a parabolic orbit within tolerance `eps`.
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::Eccentricity;
    /// assert!(Eccentricity::new_unchecked(1.0).is_parabolic(1e-10));
    /// ```
    #[must_use]
    pub fn is_parabolic(self, eps: f64) -> bool {
        (self.0 - 1.0).abs() <= eps
    }

    /// Returns `true` for a hyperbolic orbit (`e > 1`).
    ///
    /// # Examples
    ///
    /// ```
    /// use keplerian::Eccentricity;
    /// assert!(Eccentricity::new_unchecked(1.5).is_hyperbolic());
    /// ```
    #[must_use]
    pub fn is_hyperbolic(self) -> bool {
        self.0 > 1.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn is_hyperbolic_classifies_strictly_above_one() {
        assert!(Eccentricity::new_unchecked(1.5).is_hyperbolic());
        assert!(!Eccentricity::new_unchecked(0.5).is_hyperbolic());
    }

    #[test]
    fn try_new_accepts_valid_values() {
        let e = Eccentricity::try_new(0.5).unwrap();
        assert!((e.value() - 0.5).abs() < 1e-15);
        assert_eq!(Eccentricity::try_new(0.0).unwrap().value(), 0.0);
    }

    #[test]
    fn try_new_rejects_invalid_values() {
        assert!(matches!(
            Eccentricity::try_new(-0.1),
            Err(EccentricityError::Negative(_))
        ));
        assert!(matches!(
            Eccentricity::try_new(f64::NAN),
            Err(EccentricityError::NotFinite(_))
        ));
        assert!(matches!(
            Eccentricity::try_new(f64::INFINITY),
            Err(EccentricityError::NotFinite(_))
        ));
    }

    #[test]
    fn new_elliptic_validates_regime() {
        let e = Eccentricity::new_elliptic(0.3).unwrap();
        assert!((e.value() - 0.3).abs() < 1e-15);
        assert!(matches!(
            Eccentricity::new_elliptic(1.0),
            Err(EccentricityError::NotElliptic(_))
        ));
        assert!(matches!(
            Eccentricity::new_elliptic(1.5),
            Err(EccentricityError::NotElliptic(_))
        ));
        assert!(matches!(
            Eccentricity::new_elliptic(-0.1),
            Err(EccentricityError::Negative(_))
        ));
    }

    #[test]
    fn new_hyperbolic_validates_regime() {
        let e = Eccentricity::new_hyperbolic(1.5).unwrap();
        assert!((e.value() - 1.5).abs() < 1e-15);
        assert!(matches!(
            Eccentricity::new_hyperbolic(1.0),
            Err(EccentricityError::NotHyperbolic(_))
        ));
        assert!(matches!(
            Eccentricity::new_hyperbolic(0.5),
            Err(EccentricityError::NotHyperbolic(_))
        ));
        assert!(matches!(
            Eccentricity::new_hyperbolic(f64::NAN),
            Err(EccentricityError::NotFinite(_))
        ));
    }

    #[test]
    fn parabolic_constructor_is_exactly_one() {
        assert_eq!(Eccentricity::parabolic().value(), 1.0);
    }

    #[test]
    fn classify_uses_tolerance_band() {
        assert_eq!(
            Eccentricity::new_unchecked(0.5).classify(1e-10),
            ConicRegime::Elliptic
        );
        assert_eq!(
            Eccentricity::new_unchecked(1.0).classify(1e-10),
            ConicRegime::Parabolic
        );
        assert_eq!(
            Eccentricity::new_unchecked(1.0 + 5e-11).classify(1e-10),
            ConicRegime::Parabolic
        );
        assert_eq!(
            Eccentricity::new_unchecked(1.5).classify(1e-10),
            ConicRegime::Hyperbolic
        );
    }
}
