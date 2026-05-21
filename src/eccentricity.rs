// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Domain-semantic eccentricity scalar.
//!
//! ## Scientific scope
//! Eccentricity classifies conic sections in central-force motion. This module
//! models only the scalar `e` itself; it does not attach any frame, epoch, or
//! body semantics.
//!
//! ## Technical scope
//! [`Eccentricity`] wraps a raw `f64` so public APIs can distinguish orbital
//! eccentricity from unrelated dimensionless diagnostics such as tolerances or
//! residuals.
//!
//! ## References
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

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
    /// Creates a new eccentricity.
    ///
    /// Returns `None` if `e` is negative or non-finite.
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
