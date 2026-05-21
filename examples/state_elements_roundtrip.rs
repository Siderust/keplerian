//! Example binary for Cartesian/element round-trips.
#![allow(clippy::print_stdout)]

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use keplerian::elements::KeplerianElements;
use keplerian::state::CartesianState;
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
use qtty::length::Kilometer;

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
    let state = CartesianState::<Center, Frame>::new(
        Position::<Center, Frame, Kilometer>::new(7000.0, 0.0, 0.0),
        Velocity::<Frame, KmPerSecond>::new(0.0, 7.2, 1.0),
    );
    let elements = KeplerianElements::from_cartesian(&state, mu).unwrap();
    let back = elements.to_cartesian::<Center>(mu);
    assert!((back.position().x().value() - state.position().x().value()).abs() < 1e-8);
    println!(
        "a = {:.3} km, e = {:.6}",
        elements.semi_major_axis.value(),
        elements.eccentricity.value()
    );
}
