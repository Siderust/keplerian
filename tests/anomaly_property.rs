//! Property tests for typed anomaly conversions.

use keplerian::anomaly::{eccentric_from_mean, mean_from_eccentric, AnomalyOptions, MeanAnomaly};
use keplerian::Eccentricity;
use proptest::prelude::*;

proptest! {
    #[test]
    fn elliptic_mean_round_trips(m in 0.0_f64..core::f64::consts::TAU, e in 0.0_f64..0.9) {
        let ecc = Eccentricity::new(e).unwrap();
        let ea = eccentric_from_mean(MeanAnomaly::from_value(m), ecc, AnomalyOptions::default()).unwrap();
        let m2 = mean_from_eccentric(ea, ecc);
        prop_assert!((m2.value() - m).abs() < 1e-10);
    }
}
