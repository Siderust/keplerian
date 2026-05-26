// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed Lambert boundary-value solving.
//!
//! ## Scientific scope
//! Solves Lambert's two-point boundary-value problem for zero- and
//! multi-revolution Keplerian transfers using Izzo's reformulation.
//!
//! ## Technical scope
//! Public callers use the typed entry points [`lambert`] and [`lambert_n_rev`].
//! The raw array solver remains crate-private and backs the typed wrappers.
//!
//! ## References
//! - Izzo, D. (2014). *Revisiting Lambert's Problem*.
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

mod error;
mod izzo;
mod typed;

pub use error::LambertError;
pub use izzo::{LambertBranch, LambertDiagnostics, NRevBranch};
pub use typed::{lambert, lambert_n_rev, TypedLambertSolution};
