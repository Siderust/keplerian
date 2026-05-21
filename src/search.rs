// SPDX-License-Identifier: AGPL-3.0-or-later
// Copyright (C) 2026 Vallés Puig, Ramon

//! Caller-driven Lambert transfer-grid search.
//!
//! ## Scientific scope
//! This module performs brute-force grids of Lambert boundary-value solves over
//! departure epochs and flight times supplied by the caller.
//!
//! ## Technical scope
//! Search inputs and outputs remain typed through `affn` and `qtty`. The module
//! is enabled only with the `alloc` feature because it stores `Vec`-backed grids.
//!
//! ## References
//! - Izzo, D. (2014). *Revisiting Lambert's Problem*.

use alloc::vec::Vec;

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::dynamics::{GravitationalParameter, KmPerSecond, KmPerSeconds};
use qtty::length::Kilometer;
use qtty::Second;

use crate::lambert::{lambert, LambertBranch, LambertError};

/// Provider of typed positions as a function of elapsed seconds.
pub trait TrajectoryProvider<C: ReferenceCenter, F: ReferenceFrame> {
    /// Provider-specific error type.
    type Error: core::fmt::Debug;

    /// Returns the typed position at elapsed time `t`.
    fn position_at(&self, t: Second) -> Result<Position<C, F, Kilometer>, Self::Error>;
}

/// Cartesian product of departure epochs and flight times.
#[derive(Debug, Clone)]
pub struct SearchGrid {
    /// Departure elapsed seconds.
    pub departures: Vec<Second>,
    /// Flight times in seconds.
    pub flight_times: Vec<Second>,
}

/// Successful Lambert candidate for one grid cell.
#[derive(Debug, Clone, Copy)]
pub struct TransferCandidate<F: ReferenceFrame> {
    /// Departure elapsed seconds.
    pub departure: Second,
    /// Lambert time of flight.
    pub flight_time: Second,
    /// Departure velocity from Lambert.
    pub v1: Velocity<F, KmPerSecond>,
    /// Arrival velocity from Lambert.
    pub v2: Velocity<F, KmPerSecond>,
    /// Sum of Lambert endpoint speed magnitudes.
    pub total_dv: KmPerSeconds,
}

/// Outcome stored for each search grid cell.
#[derive(Debug)]
pub enum CellOutcome<F: ReferenceFrame, E> {
    /// Lambert solve succeeded.
    Success(TransferCandidate<F>),
    /// Lambert solve failed after both positions were available.
    LambertFailed(LambertError),
    /// A trajectory provider failed.
    ProviderFailed(E),
}

/// Full row-major candidate grid.
#[derive(Debug)]
pub struct CandidateGrid<F: ReferenceFrame, E> {
    /// The searched grid.
    pub grid: SearchGrid,
    /// Cells indexed as `[departure_index][flight_time_index]`.
    pub cells: Vec<Vec<CellOutcome<F, E>>>,
}

/// Runs Lambert on every source/target position pair in a grid.
#[must_use]
pub fn lambert_search<C, F, P1, P2>(
    source: &P1,
    target: &P2,
    grid: SearchGrid,
    mu: GravitationalParameter,
    branch: LambertBranch,
) -> CandidateGrid<F, P1::Error>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
    P1: TrajectoryProvider<C, F>,
    P2: TrajectoryProvider<C, F, Error = P1::Error>,
{
    let mut cells = Vec::with_capacity(grid.departures.len());
    for dep in &grid.departures {
        let mut row = Vec::with_capacity(grid.flight_times.len());
        for tof in &grid.flight_times {
            let outcome = match source.position_at(*dep) {
                Err(e) => CellOutcome::ProviderFailed(e),
                Ok(r1) => match target.position_at(Second::new(dep.value() + tof.value())) {
                    Err(e) => CellOutcome::ProviderFailed(e),
                    Ok(r2) => match lambert(r1, r2, *tof, mu, branch) {
                        Err(e) => CellOutcome::LambertFailed(e),
                        Ok(sol) => {
                            let total = speed(&sol.v1) + speed(&sol.v2);
                            CellOutcome::Success(TransferCandidate {
                                departure: *dep,
                                flight_time: *tof,
                                v1: sol.v1,
                                v2: sol.v2,
                                total_dv: KmPerSeconds::new(total),
                            })
                        }
                    },
                },
            };
            row.push(outcome);
        }
        cells.push(row);
    }
    CandidateGrid { grid, cells }
}

fn speed<F: ReferenceFrame>(v: &Velocity<F, KmPerSecond>) -> f64 {
    let x = v.x().value();
    let y = v.y().value();
    let z = v.z().value();
    (x * x + y * y + z * z).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

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
    #[derive(Debug)]
    struct Provider {
        x: f64,
        fail: bool,
    }
    impl TrajectoryProvider<C, F> for Provider {
        type Error = &'static str;
        fn position_at(&self, _: Second) -> Result<Position<C, F, Kilometer>, Self::Error> {
            if self.fail {
                Err("failed")
            } else {
                Ok(Position::<C, F, Kilometer>::new(self.x, 0.0, 0.0))
            }
        }
    }

    #[test]
    fn grid_dimensions_and_failures() {
        let grid = SearchGrid {
            departures: alloc::vec![Second::new(0.0), Second::new(1.0)],
            flight_times: alloc::vec![Second::new(1000.0)],
        };
        let out = lambert_search(
            &Provider {
                x: 7000.0,
                fail: true,
            },
            &Provider {
                x: 8000.0,
                fail: false,
            },
            grid,
            GravitationalParameter::new(398600.4418),
            LambertBranch::Prograde,
        );
        assert_eq!(out.cells.len(), 2);
        assert_eq!(out.cells[0].len(), 1);
        assert!(matches!(out.cells[0][0], CellOutcome::ProviderFailed(_)));
    }
}
