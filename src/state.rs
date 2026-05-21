// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Typed Cartesian two-body state vectors.
//!
//! ## Scientific scope
//! This module stores the instantaneous Cartesian state needed for two-body
//! propagation and boundary-value problems.
//!
//! ## Technical scope
//! [`CartesianState`] combines an `affn` typed position and velocity while
//! preserving center and frame tags.
//!
//! ## References
//! - Vallado, D. A. (2013). *Fundamentals of Astrodynamics and Applications*.

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::dynamics::KmPerSecond;
use qtty::length::Kilometer;

/// A Cartesian state with typed center and frame markers.
#[derive(Debug)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
#[cfg_attr(
    feature = "serde",
    serde(bound(
        serialize = "Position<C, F, Kilometer>: serde::Serialize, Velocity<F, KmPerSecond>: serde::Serialize",
        deserialize = "Position<C, F, Kilometer>: serde::Deserialize<'de>, Velocity<F, KmPerSecond>: serde::Deserialize<'de>"
    ))
)]
pub struct CartesianState<C: ReferenceCenter, F: ReferenceFrame> {
    /// Position measured from center `C` in frame `F`, kilometres.
    pub position: Position<C, F, Kilometer>,
    /// Free velocity vector in frame `F`, kilometres per second.
    pub velocity: Velocity<F, KmPerSecond>,
}

impl<C: ReferenceCenter, F: ReferenceFrame> Clone for CartesianState<C, F>
where
    Position<C, F, Kilometer>: Clone,
    Velocity<F, KmPerSecond>: Clone,
{
    fn clone(&self) -> Self {
        Self {
            position: self.position.clone(),
            velocity: self.velocity,
        }
    }
}

impl<C: ReferenceCenter, F: ReferenceFrame> Copy for CartesianState<C, F>
where
    Position<C, F, Kilometer>: Copy,
    Velocity<F, KmPerSecond>: Copy,
{
}

impl<C: ReferenceCenter, F: ReferenceFrame> CartesianState<C, F> {
    /// Creates a new Cartesian state.
    #[must_use]
    pub fn new(position: Position<C, F, Kilometer>, velocity: Velocity<F, KmPerSecond>) -> Self {
        Self { position, velocity }
    }

    /// Returns the typed position.
    #[must_use]
    pub fn position(&self) -> &Position<C, F, Kilometer> {
        &self.position
    }

    /// Returns the typed velocity.
    #[must_use]
    pub fn velocity(&self) -> &Velocity<F, KmPerSecond> {
        &self.velocity
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[derive(Debug, Copy, Clone)]
    struct TestCenter;
    impl ReferenceCenter for TestCenter {
        type Params = ();
        fn center_name() -> &'static str {
            "TestCenter"
        }
    }

    #[derive(Debug, Copy, Clone)]
    struct TestFrame;
    impl ReferenceFrame for TestFrame {
        fn frame_name() -> &'static str {
            "TestFrame"
        }
    }

    #[test]
    fn preserves_center_frame_tags() {
        let s = CartesianState::<TestCenter, TestFrame>::new(
            Position::<TestCenter, TestFrame, Kilometer>::new(1.0, 2.0, 3.0),
            Velocity::<TestFrame, KmPerSecond>::new(4.0, 5.0, 6.0),
        );
        assert_eq!(s.position().z().value(), 3.0);
        assert_eq!(s.velocity().x().value(), 4.0);
    }
}
