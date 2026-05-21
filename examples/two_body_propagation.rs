//! Example binary for two-body propagation.
#![allow(clippy::print_stdout)]

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use keplerian::{problem::KeplerProblem, state::CartesianState};
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
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

fn main() {
    let mu = GravitationalParameter::new(398_600.441_8);
    let r = 7000.0;
    let state = CartesianState::<Center, Frame>::new(
        Position::<Center, Frame, Kilometer>::new(r, 0.0, 0.0),
        Velocity::<Frame, KmPerSecond>::new(0.0, (mu.value() / r).sqrt(), 0.0),
    );
    let period = 2.0 * core::f64::consts::PI * (r * r * r / mu.value()).sqrt();
    let out = KeplerProblem::<Center, Frame>::new(mu)
        .propagate(&state, Second::new(period))
        .unwrap();
    assert!((out.position().x().value() - r).abs() < 1e-6);
    println!("returned to x = {:.6} km", out.position().x().value());
}
