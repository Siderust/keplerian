//! Example binary for caller-driven Lambert search.
#![cfg(feature = "alloc")]
#![allow(clippy::print_stdout)]

use affn::cartesian::Position;
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use keplerian::lambert::LambertBranch;
use keplerian::search::{lambert_search, CellOutcome, SearchGrid, TrajectoryProvider};
use qtty::dynamics::GravitationalParameter;
use qtty::length::Kilometer;
use qtty::Second;

#[derive(Debug, Copy, Clone)]
struct Center;
impl ReferenceCenter for Center {
    type Params = ();
    fn center_name() -> &'static str {
        "Center"
    }
}
#[derive(Debug, Copy, Clone)]
struct Frame;
impl ReferenceFrame for Frame {
    fn frame_name() -> &'static str {
        "Frame"
    }
}
struct Circle {
    radius: f64,
    phase: f64,
}
impl TrajectoryProvider<Center, Frame> for Circle {
    type Error = &'static str;
    fn position_at(&self, t: Second) -> Result<Position<Center, Frame, Kilometer>, Self::Error> {
        let theta = self.phase + t.value() / 86_400.0 * 0.01;
        Ok(Position::<Center, Frame, Kilometer>::new(
            self.radius * theta.cos(),
            self.radius * theta.sin(),
            0.0,
        ))
    }
}
fn main() {
    let grid = SearchGrid {
        departures: (0..5).map(|i| Second::new(i as f64 * 86_400.0)).collect(),
        flight_times: (1..=5).map(|i| Second::new(i as f64 * 3600.0)).collect(),
    };
    let out = lambert_search(
        &Circle {
            radius: 7000.0,
            phase: 0.0,
        },
        &Circle {
            radius: 9000.0,
            phase: 0.8,
        },
        grid,
        GravitationalParameter::new(398600.4418),
        LambertBranch::Prograde,
    );
    let successes = out
        .cells
        .iter()
        .flatten()
        .filter(|c| matches!(c, CellOutcome::Success(_)))
        .count();
    println!(
        "{}x{} grid, {successes} successful cells",
        out.grid.departures.len(),
        out.grid.flight_times.len()
    );
}
