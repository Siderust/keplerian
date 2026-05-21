//! Integration tests for typed Lambert validation cases.

use affn::cartesian::Position;
use affn::frames::ICRS;
use keplerian::lambert::{lambert, LambertBranch};
use qtty::dynamics::GravitationalParameter;
use qtty::length::Kilometer;
use qtty::Second;

#[test]
fn vallado_example_7_5() {
    let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
    let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
    let sol = lambert(
        r1,
        r2,
        Second::new(4560.0),
        GravitationalParameter::new(398600.4418),
        LambertBranch::Prograde,
    )
    .unwrap();
    assert!((sol.v1.x().value() - 2.058913).abs() < 1e-3);
    assert!((sol.v1.y().value() - 2.915965).abs() < 1e-3);
    assert!((sol.v2.x().value() - -3.451565).abs() < 1e-3);
    assert!((sol.v2.y().value() - 0.910315).abs() < 1e-3);
}
