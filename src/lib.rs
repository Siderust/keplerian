// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Domain-agnostic Keplerian dynamics on typed quantities.
//!
//! ## Scientific scope
//! This crate provides reusable two-body astrodynamics primitives: anomaly
//! solving, Keplerian elements, Lambert transfers, and transfer invariants.
//!
//! ## Technical scope
//! Public APIs depend only on `qtty` and `affn`; no time-scale or astronomy
//! orchestration types are introduced here.
//!
//! ## References
//! - Battin, R. H. (1999). *An Introduction to the Mathematics and Methods of
//!   Astrodynamics*.
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

#![forbid(unsafe_code)]
#![deny(missing_docs)]
#![cfg_attr(not(feature = "std"), no_std)]
#![doc = include_str!("../README.md")]

#[cfg(feature = "alloc")]
extern crate alloc;

pub mod anomaly;
pub mod eccentricity;
pub mod elements;
pub mod error;
pub mod lambert;
pub mod prelude;
pub mod problem;
pub mod state;
pub mod transfer;

#[cfg(feature = "alloc")]
pub mod search;

mod vec3;

pub use eccentricity::Eccentricity;
pub use error::KeplerError;
