//! Example binary for elliptic Kepler solves.
#![allow(clippy::print_stdout)]

use keplerian::anomaly::{kepler_elliptic, AnomalyOptions, MeanAnomaly};
use keplerian::Eccentricity;

fn main() {
    let ecc = Eccentricity::new(0.2).unwrap();
    for m in [0.0, 0.5, 1.0, 1.5] {
        let e =
            kepler_elliptic(MeanAnomaly::from_value(m), ecc, AnomalyOptions::default()).unwrap();
        println!("M = {m:.3} rad -> E = {:.6} rad", e.value());
    }
}
