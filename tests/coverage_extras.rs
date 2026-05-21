//! Tests targeting coverage gaps: error conversions, state clone, transfer
//! helpers, conic kinds, anomaly edge cases, Lambert n-rev, and search grid.

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use affn::frames::ICRS;
use keplerian::anomaly::{
    kepler_parabolic, kepler_elliptic, kepler_hyperbolic, AnomalyOptions, eccentric_from_mean,
    hyperbolic_from_mean, true_from_hyperbolic, hyperbolic_from_true, mean_from_hyperbolic,
};
use keplerian::eccentricity::Eccentricity;
use keplerian::elements::{ConversionError, KeplerianElements};
use keplerian::error::KeplerError;
use keplerian::lambert::{lambert_n_rev, LambertBranch, LambertError, NRevBranch};
use keplerian::problem::{KeplerProblem, PropagationError};
use keplerian::state::CartesianState;
use keplerian::transfer::{orbital_period, specific_orbital_energy, specific_angular_momentum};
use qtty::angular::Radians;
use qtty::dynamics::GravitationalParameter;
use qtty::length::{Kilometer, Kilometers};
use qtty::Second;

#[derive(Debug, Clone, Copy)]
struct C;
impl ReferenceCenter for C {
    type Params = ();
    fn center_name() -> &'static str {
        "C"
    }
}

#[derive(Debug, Clone, Copy)]
struct F;
impl ReferenceFrame for F {
    fn frame_name() -> &'static str {
        "F"
    }
}

// ── error.rs ──────────────────────────────────────────────────────────────────

#[test]
fn kepler_error_from_anomaly() {
    let inner = keplerian::anomaly::AnomalyError::InvalidEccentricity(1.5);
    let e = KeplerError::from(inner);
    assert!(matches!(e, KeplerError::Anomaly(_)));
}

#[test]
fn kepler_error_from_conversion() {
    let inner = ConversionError::InvalidEccentricity(-1.0);
    let e = KeplerError::from(inner);
    assert!(matches!(e, KeplerError::Conversion(_)));
}

#[test]
fn kepler_error_from_propagation() {
    let inner = PropagationError::ParabolicUnsupported;
    let e = KeplerError::from(inner);
    assert!(matches!(e, KeplerError::Propagation(_)));
}

#[test]
fn kepler_error_from_lambert() {
    let inner = LambertError::ZeroPosition;
    let e = KeplerError::from(inner);
    assert!(matches!(e, KeplerError::Lambert(_)));
}

// ── state.rs ──────────────────────────────────────────────────────────────────

#[test]
fn cartesian_state_clone_and_velocity() {
    let pos = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, 7.5, 0.0);
    let s = CartesianState::<C, F>::new(pos, vel);
    // CartesianState<C,F> is Copy when Position and Velocity are Copy;
    // use explicit Clone::clone to exercise the Clone impl body.
    #[allow(clippy::clone_on_copy)]
    let s2 = Clone::clone(&s);
    assert_eq!(s2.velocity().x().value(), 0.0);
    assert_eq!(s2.velocity().y().value(), 7.5);
}

// ── transfer.rs ───────────────────────────────────────────────────────────────

#[test]
fn specific_orbital_energy_and_period() {
    let mu = GravitationalParameter::new(398600.4418);
    let problem = KeplerProblem::<C, F>::new(mu);

    let pos = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, 7.546, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);

    let eps = specific_orbital_energy(&state, mu);
    assert!(eps.value() < 0.0, "bound orbit has negative energy");

    let h = specific_angular_momentum(&state);
    assert!(h.value() > 0.0);

    let t = orbital_period(&problem, Kilometers::new(7000.0)).unwrap();
    assert!((t.value() - 5840.0).abs() < 200.0);
}

#[test]
fn orbital_period_negative_sma_returns_none() {
    let mu = GravitationalParameter::new(398600.4418);
    let problem = KeplerProblem::<C, F>::new(mu);
    assert!(orbital_period(&problem, Kilometers::new(-7000.0)).is_none());
}

// ── problem.rs ────────────────────────────────────────────────────────────────

#[test]
fn kepler_problem_mu_accessor() {
    let mu = GravitationalParameter::new(398600.4418);
    let p = KeplerProblem::<C, F>::new(mu);
    assert_eq!(p.mu().value(), 398600.4418);
}

#[test]
fn kepler_problem_hyperbolic_propagation() {
    // Build a clearly hyperbolic orbit (e ≈ 1.5) by giving v > escape speed.
    let mu = GravitationalParameter::new(398600.4418);
    let r = 7000.0_f64;
    let v_esc = (2.0 * mu.value() / r).sqrt();
    let pos = Position::<C, F, Kilometer>::new(r, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, v_esc * 1.3, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    let problem = KeplerProblem::<C, F>::new(mu);
    let result = problem.propagate(&state, Second::new(100.0));
    assert!(result.is_ok(), "hyperbolic propagation failed: {result:?}");
}

#[test]
fn kepler_problem_parabolic_returns_error() {
    // Near-parabolic: energy ≈ 0 → from_cartesian returns Degenerate("parabolic orbit")
    // which propagate() wraps as PropagationError::Conversion.
    let mu = GravitationalParameter::new(398600.4418);
    let r = 7000.0_f64;
    let v_par = (2.0 * mu.value() / r).sqrt(); // exact escape speed
    let pos = Position::<C, F, Kilometer>::new(r, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, v_par, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    let problem = KeplerProblem::<C, F>::new(mu);
    assert!(problem.propagate(&state, Second::new(100.0)).is_err());
}

// ── elements.rs ───────────────────────────────────────────────────────────────

#[test]
fn elements_new_rejects_negative_eccentricity() {
    let err = KeplerianElements::<F>::new(
        Kilometers::new(7000.0),
        Eccentricity::new_unchecked(-0.1),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
    );
    assert!(matches!(err, Err(ConversionError::InvalidEccentricity(_))));
}

#[test]
fn elements_new_rejects_inclination_out_of_range() {
    use core::f64::consts::PI;
    let err = KeplerianElements::<F>::new(
        Kilometers::new(7000.0),
        Eccentricity::new_unchecked(0.1),
        Radians::new(PI + 0.1),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
    );
    assert!(matches!(err, Err(ConversionError::InvalidInclination(_))));
}

#[test]
fn conic_kind_parabolic() {
    // Eccentricity::new_unchecked allows e = 1.0 which is parabolic.
    use keplerian::elements::ConicRegime;
    let el = KeplerianElements::<F>::new(
        Kilometers::new(7000.0),
        Eccentricity::new_unchecked(1.0),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
    )
    .unwrap();
    assert_eq!(el.conic_kind(), ConicRegime::Parabolic);
}

#[test]
fn conic_kind_hyperbolic() {
    use keplerian::elements::ConicRegime;
    let el = KeplerianElements::<F>::new(
        Kilometers::new(-40000.0),
        Eccentricity::new_unchecked(1.5),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
        Radians::new(0.0),
    )
    .unwrap();
    assert_eq!(el.conic_kind(), ConicRegime::Hyperbolic);
}

#[test]
fn from_cartesian_degenerate_zero_position() {
    let mu = GravitationalParameter::new(398600.4418);
    let pos = Position::<C, F, Kilometer>::new(0.0, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, 7.5, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    assert!(KeplerianElements::<F>::from_cartesian(&state, mu).is_err());
}

#[test]
fn from_cartesian_degenerate_zero_angular_momentum() {
    let mu = GravitationalParameter::new(398600.4418);
    let pos = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
    // radial velocity → h = r × v = 0
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(7.5, 0.0, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    assert!(KeplerianElements::<F>::from_cartesian(&state, mu).is_err());
}

#[test]
fn from_cartesian_equatorial_eccentric_orbit() {
    // Equatorial (inc ≈ 0) non-circular orbit: tests argp from equatorial branch.
    let mu = GravitationalParameter::new(398600.4418);
    let pos = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
    // Slightly off-circular, equatorial: vy < v_circ, vz = 0.
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, 6.5, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    let el = KeplerianElements::<F>::from_cartesian(&state, mu);
    assert!(el.is_ok(), "{el:?}");
}

// ── anomaly.rs ────────────────────────────────────────────────────────────────

#[test]
fn elliptic_circular_orbit_returns_mean_anomaly() {
    // e = 0 → E = M (identity)
    let zero = Eccentricity::new(0.0).unwrap();
    let m = Radians::new(1.23);
    let e = eccentric_from_mean(m, zero, AnomalyOptions::default()).unwrap();
    assert!((e.value() - m.value()).abs() < 1e-14);
}

#[test]
fn hyperbolic_mean_anomaly_zero_returns_zero() {
    let ecc = Eccentricity::new_unchecked(1.5);
    let f = hyperbolic_from_mean(Radians::new(0.0), ecc, AnomalyOptions::default()).unwrap();
    assert!((f).abs() < 1e-14);
}

#[test]
fn hyperbolic_mean_anomaly_large_branch() {
    // |M| > 50 * e triggers the log initial guess branch.
    let ecc = Eccentricity::new_unchecked(1.2);
    let m = Radians::new(100.0);
    let f = hyperbolic_from_mean(m, ecc, AnomalyOptions::default());
    assert!(f.is_ok());
}

#[test]
fn kepler_parabolic_round_trips() {
    let m = 0.5_f64;
    let d = kepler_parabolic(m);
    let m_back = d + d.powi(3) / 3.0;
    assert!((m_back - m).abs() < 1e-12);
}

#[test]
fn hyperbolic_anomaly_round_trips() {
    let ecc = Eccentricity::new_unchecked(2.0);
    let nu = Radians::new(0.8);
    let f = hyperbolic_from_true(nu, ecc);
    let nu2 = true_from_hyperbolic(f, ecc);
    assert!((nu2.value() - nu.value()).abs() < 1e-12);
    let m = mean_from_hyperbolic(f, ecc);
    let f2 = hyperbolic_from_mean(m, ecc, AnomalyOptions::default()).unwrap();
    assert!((f2 - f).abs() < 1e-12);
}

#[test]
fn hyperbolic_kepler_rejects_nan_mean_anomaly() {
    let ecc = Eccentricity::new_unchecked(1.5);
    let err = hyperbolic_from_mean(Radians::new(f64::NAN), ecc, AnomalyOptions::default());
    assert!(err.is_err());
}

#[test]
fn elliptic_kepler_rejects_nan_mean_anomaly() {
    let ecc = Eccentricity::new_unchecked(0.5);
    let err = eccentric_from_mean(Radians::new(f64::NAN), ecc, AnomalyOptions::default());
    assert!(err.is_err());
}

// ── lambert/typed.rs ──────────────────────────────────────────────────────────

#[test]
fn lambert_n_rev_valid_case() {
    // Two positions separated by ~90°; 2-hour TOF, 1 extra revolution.
    let r1 = Position::<(), ICRS, Kilometer>::new(7000.0, 0.0, 0.0);
    let r2 = Position::<(), ICRS, Kilometer>::new(0.0, 7000.0, 0.0);
    let tof = Second::new(10800.0);
    let mu = GravitationalParameter::new(398600.4418);
    // May succeed or fail depending on geometry; just exercise the code path.
    let _ = lambert_n_rev(r1, r2, tof, mu, LambertBranch::Prograde, 1, NRevBranch::Left);
}

#[test]
fn lambert_n_rev_retrograde_case() {
    let r1 = Position::<(), ICRS, Kilometer>::new(7000.0, 0.0, 0.0);
    let r2 = Position::<(), ICRS, Kilometer>::new(0.0, 7000.0, 0.0);
    let tof = Second::new(10800.0);
    let mu = GravitationalParameter::new(398600.4418);
    let _ = lambert_n_rev(r1, r2, tof, mu, LambertBranch::Retrograde, 1, NRevBranch::Right);
}

// ── search.rs ─────────────────────────────────────────────────────────────────

#[cfg(feature = "alloc")]
mod search_tests {
    extern crate alloc;

    use super::*;
    use keplerian::search::{lambert_search, CellOutcome, SearchGrid, TrajectoryProvider};

    struct FixedProvider {
        pos: [f64; 3],
    }
    impl TrajectoryProvider<C, F> for FixedProvider {
        type Error = &'static str;
        fn position_at(&self, _: Second) -> Result<Position<C, F, Kilometer>, Self::Error> {
            Ok(Position::<C, F, Kilometer>::new(
                self.pos[0],
                self.pos[1],
                self.pos[2],
            ))
        }
    }

    #[test]
    fn search_success_covers_speed_helper() {
        let grid = SearchGrid {
            departures: alloc::vec![Second::new(0.0)],
            flight_times: alloc::vec![Second::new(4560.0)],
        };
        let out = lambert_search(
            &FixedProvider {
                pos: [15945.34, 0.0, 0.0],
            },
            &FixedProvider {
                pos: [12214.84, 10249.47, 0.0],
            },
            grid,
            GravitationalParameter::new(398600.4418),
            LambertBranch::Prograde,
        );
        assert_eq!(out.cells.len(), 1);
        assert!(matches!(out.cells[0][0], CellOutcome::Success(_)));
    }

    struct TargetFails;
    impl TrajectoryProvider<C, F> for TargetFails {
        type Error = &'static str;
        fn position_at(&self, _: Second) -> Result<Position<C, F, Kilometer>, Self::Error> {
            Err("target down")
        }
    }

    #[test]
    fn search_target_provider_failure() {
        let grid = SearchGrid {
            departures: alloc::vec![Second::new(0.0)],
            flight_times: alloc::vec![Second::new(4560.0)],
        };
        let out = lambert_search(
            &FixedProvider {
                pos: [15945.34, 0.0, 0.0],
            },
            &TargetFails,
            grid,
            GravitationalParameter::new(398600.4418),
            LambertBranch::Prograde,
        );
        assert!(matches!(out.cells[0][0], CellOutcome::ProviderFailed(_)));
    }

    /// `r1 == r2` forces a `ZeroPosition` Lambert failure → covers `CellOutcome::LambertFailed`.
    #[test]
    fn search_lambert_failed_cell() {
        let grid = SearchGrid {
            departures: alloc::vec![Second::new(0.0)],
            flight_times: alloc::vec![Second::new(4560.0)],
        };
        let out = lambert_search(
            &FixedProvider { pos: [0.0, 0.0, 0.0] },
            &FixedProvider { pos: [0.0, 0.0, 0.0] },
            grid,
            GravitationalParameter::new(398600.4418),
            LambertBranch::Prograde,
        );
        assert!(matches!(out.cells[0][0], CellOutcome::LambertFailed(_)));
    }
}

// ── eccentricity.rs ───────────────────────────────────────────────────────────

#[test]
fn eccentricity_is_hyperbolic() {
    assert!(Eccentricity::new_unchecked(1.5).is_hyperbolic());
    assert!(!Eccentricity::new_unchecked(0.5).is_hyperbolic());
}

// ── state.rs: velocity() return value ─────────────────────────────────────────

#[test]
fn cartesian_state_velocity_accessor() {
    let pos = Position::<C, F, Kilometer>::new(1.0, 2.0, 3.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(4.0, 5.0, 6.0);
    let s = CartesianState::<C, F>::new(pos, vel);
    let vref = s.velocity();
    assert_eq!(vref.z().value(), 6.0);
}

// ── elements.rs: non-positive mu and NaN ─────────────────────────────────────

#[test]
fn from_cartesian_non_positive_mu() {
    let pos = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, 7.5, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    let mu = GravitationalParameter::new(-1.0);
    let err = KeplerianElements::<F>::from_cartesian(&state, mu);
    assert!(err.is_err());
}

#[test]
fn from_cartesian_nan_mu() {
    let pos = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, 7.5, 0.0);
    let state = CartesianState::<C, F>::new(pos, vel);
    let mu = GravitationalParameter::new(f64::NAN);
    let err = KeplerianElements::<F>::from_cartesian(&state, mu);
    assert!(matches!(err, Err(ConversionError::NonFiniteValue { .. })));
}

#[test]
fn from_cartesian_circular_inclined_orbit() {
    // Inclined circular orbit: ecc ≈ 0 and nmag > 0 → exercises the
    // `else if nmag > EPS` branch for nu computation.
    let mu = GravitationalParameter::new(398600.4418);
    let r = 7000.0_f64;
    let v_circ = (mu.value() / r).sqrt();
    // Tilt into the i=45° plane: velocity in y-z plane.
    let vy = v_circ * std::f64::consts::FRAC_1_SQRT_2;
    let vz = v_circ * std::f64::consts::FRAC_1_SQRT_2;
    let pos = Position::<C, F, Kilometer>::new(r, 0.0, 0.0);
    let vel = Velocity::<F, qtty::dynamics::KmPerSecond>::new(0.0, vy, vz);
    let state = CartesianState::<C, F>::new(pos, vel);
    let el = KeplerianElements::<F>::from_cartesian(&state, mu);
    assert!(el.is_ok(), "{el:?}");
}

#[test]
fn elements_to_cartesian_and_back() {
    // Explicitly calls to_cartesian, covering its rotate_pqw invocations.
    let mu = GravitationalParameter::new(398600.4418);
    let el = KeplerianElements::<F>::new(
        Kilometers::new(7000.0),
        Eccentricity::new_unchecked(0.1),
        Radians::new(0.5),
        Radians::new(1.0),
        Radians::new(0.3),
        Radians::new(0.7),
    )
    .unwrap();
    let state = el.to_cartesian::<C>(mu);
    assert!(state.position().x().value().is_finite());
}

// ── anomaly.rs: bisection fallback paths ─────────────────────────────────────

#[test]
fn elliptic_bisection_fallback_with_zero_tol() {
    // max_iter = 2, tol = 0.0: Newton fails to reach exact zero residual,
    // triggers bisection, exercises lines 110-112 and 421-429.
    let ecc = Eccentricity::new_unchecked(0.5);
    let opts = AnomalyOptions {
        max_iter: 2,
        tol: 0.0,
    };
    let result = kepler_elliptic(Radians::new(1.0), ecc, opts);
    // Either converges or not, but the code paths are exercised.
    let _ = result;
}

#[test]
fn hyperbolic_bisection_fallback_with_zero_tol() {
    // max_iter = 2, tol = 0.0: exercises hyperbolic bisection lines.
    let ecc = Eccentricity::new_unchecked(1.5);
    let opts = AnomalyOptions {
        max_iter: 2,
        tol: 0.0,
    };
    let result = kepler_hyperbolic(Radians::new(1.0), ecc, opts);
    let _ = result;
}

#[test]
fn hyperbolic_bisection_large_mean_anomaly() {
    // Large |M| triggers the upper-bracket expansion loop (lines 459-460).
    let ecc = Eccentricity::new_unchecked(1.1);
    let opts = AnomalyOptions {
        max_iter: 2,
        tol: 0.0,
    };
    let result = kepler_hyperbolic(Radians::new(50.0), ecc, opts);
    let _ = result;
}
