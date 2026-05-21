//! Unified error type for Lambert solving.
//!
//! ## Scientific scope
//! This module classifies numerical and geometric failure modes of Lambert's
//! two-point boundary-value problem.
//!
//! ## Technical scope
//! [`LambertError`] is shared by both typed public entry points and the
//! crate-private raw numeric backend.
//!
//! ## References
//! - Izzo, D. (2014). *Revisiting Lambert's Problem*.

/// Errors returned by the Lambert solver.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum LambertError {
    /// Gravitational parameter must be strictly positive (km³/s²).
    #[error("non-positive gravitational parameter ({0})")]
    NonPositiveMu(f64),
    /// Time of flight must be a finite, strictly positive duration.
    #[error("non-positive time of flight ({0} s)")]
    NonPositiveTof(f64),
    /// Initial or final position vector has effectively zero magnitude.
    #[error("zero-magnitude position vector")]
    ZeroPosition,
    /// `r1`, `r2` and the origin are collinear, so the chord/transfer plane is undefined.
    #[error("collinear positions: chord cannot be resolved unambiguously")]
    Collinear,
    /// The requested revolution count exceeds the maximum permitted by the supplied time of flight.
    #[error("requested {requested} revolutions exceeds N_max = {max} for the given TOF")]
    RevolutionsExceedNMax {
        /// Number of revolutions requested by the caller.
        requested: u32,
        /// Largest physically admissible revolution count.
        max: u32,
    },
    /// Householder iteration failed to converge within the iteration cap.
    #[error("Householder iteration failed to converge (residual = {residual:e})")]
    DidNotConverge {
        /// Final residual `|T(x) − T*|` in non-dimensional units.
        residual: f64,
    },
}
