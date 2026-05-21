//! Integration tests for propagation invariants.

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use keplerian::problem::KeplerProblem;
use keplerian::state::CartesianState;
use keplerian::transfer::{specific_angular_momentum, specific_orbital_energy};
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
use qtty::length::Kilometer;
use qtty::Second;

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

#[test]
fn propagation_preserves_invariants_and_reverses() {
    let mu = GravitationalParameter::new(398600.4418);
    let state = CartesianState::<C, F>::new(
        Position::<C, F, Kilometer>::new(8000.0, 1000.0, 200.0),
        Velocity::<F, KmPerSecond>::new(-0.5, 7.0, 1.0),
    );
    let problem = KeplerProblem::<C, F>::new(mu);
    let fwd = problem.propagate(&state, Second::new(1200.0)).unwrap();
    let back = problem.propagate(&fwd, Second::new(-1200.0)).unwrap();
    assert!(
        (specific_orbital_energy(&state, mu) - specific_orbital_energy(&fwd, mu))
            .value()
            .abs()
            < 1e-8
    );
    assert!(
        (specific_angular_momentum(&state) - specific_angular_momentum(&fwd))
            .value()
            .abs()
            < 1e-8
    );
    assert!((back.position().x().value() - state.position().x().value()).abs() < 1e-6);
    assert!((back.position().y().value() - state.position().y().value()).abs() < 1e-6);
}
