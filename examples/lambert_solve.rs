//! Example binary for a typed Lambert solve.
#![allow(clippy::print_stdout)]

use affn::cartesian::Position;
use affn::frames::ICRS;
use keplerian::lambert::{lambert, LambertBranch};
use qtty::dynamics::GravitationalParameter;
use qtty::length::Kilometer;
use qtty::Second;

const AU_KM: f64 = 1.495_978_707e8;
const MU_SUN: f64 = 1.327_124_400_18e11;

fn main() {
    let r_earth_km = AU_KM;
    let r_mars_km = 1.524 * AU_KM;
    let phase = 60.0_f64.to_radians();
    let r1 = Position::<(), ICRS, Kilometer>::new(r_earth_km, 0.0, 0.0);
    let r2 =
        Position::<(), ICRS, Kilometer>::new(r_mars_km * phase.cos(), r_mars_km * phase.sin(), 0.0);
    let solution = lambert(
        r1,
        r2,
        Second::new(258.0 * 86_400.0),
        GravitationalParameter::new(MU_SUN),
        LambertBranch::Prograde,
    )
    .unwrap();
    println!(
        "v1 = ({:.4}, {:.4}, {:.4}) km/s",
        solution.v1.x().value(),
        solution.v1.y().value(),
        solution.v1.z().value()
    );
}
