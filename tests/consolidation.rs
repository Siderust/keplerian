//! Regression tests for conic geometry projection and canonical anomalies.
use affn::{conic::ConicKind, frames::ReferenceFrame};
use keplerian::{
    anomaly::{eccentric_from_true, true_from_eccentric, EccentricAnomaly},
    elements::KeplerianElements,
    Eccentricity,
};
use qtty::{angular::Radians, length::Kilometers};
#[derive(Clone, Copy, Debug)]
struct Frame;
impl ReferenceFrame for Frame {
    fn frame_name() -> &'static str {
        "Frame"
    }
}

#[test]
fn anomaly_conversion_matches_independent_cartesian_geometry() {
    for e in [0.0, 1e-12, 0.8, 0.999999, 1.0 - f64::EPSILON] {
        for angle in [
            -1e12,
            -20.0,
            -core::f64::consts::PI,
            -1e-8,
            0.0,
            1e-8,
            core::f64::consts::PI,
            20.0,
            1e12,
        ] {
            let ecc = Eccentricity::new(e).unwrap();
            let nu = true_from_eccentric(EccentricAnomaly::from_value(angle), ecc);
            assert!((0.0..core::f64::consts::TAU).contains(&nu.value()));
            // Focus-centered ellipse coordinates provide an independent oracle.
            let x = (1.0 - e) - 2.0 * (angle * 0.5).sin().powi(2);
            let y = ((1.0 - e) * (1.0 + e)).sqrt() * angle.sin();
            let expected = y.atan2(x);
            assert!(
                (nu.value().sin() - expected.sin()).abs() < 3e-8,
                "e={e}, E={angle}"
            );
            assert!((nu.value().cos() - expected.cos()).abs() < 3e-8);
            // Inversion is ill-conditioned as e approaches one; use regular e.
            if e < 0.99 {
                let back = eccentric_from_true(nu, ecc).value();
                assert!((back.sin() - angle.sin()).abs() < 1e-12);
                assert!((back.cos() - angle.cos()).abs() < 1e-12);
            }
        }
    }
}

#[test]
fn geometry_projection_preserves_conventions_and_revalidates_public_fields() {
    for (a, e, kind) in [
        (7000.0, 0.0, ConicKind::Elliptic),
        (-9000.0, 1.5, ConicKind::Hyperbolic),
    ] {
        let mut el = KeplerianElements::<Frame>::new(
            Kilometers::new(a),
            Eccentricity::new(e).unwrap(),
            Radians::new(0.7),
            Radians::new(1.2),
            Radians::new(0.4),
            Radians::new(0.3),
        )
        .unwrap();
        let geometry = el.to_geometry().unwrap();
        assert_eq!(geometry.kind(), kind);
        assert_eq!(geometry.shape().semi_major_axis().value(), a.abs());
        assert_eq!(geometry.shape().eccentricity(), e);
        assert!(
            (geometry
                .orientation()
                .inclination()
                .to::<qtty::units::Radian>()
                .value()
                - 0.7)
                .abs()
                < 1e-15
        );
        el.semi_major_axis = Kilometers::new(-a);
        assert!(el.to_geometry().is_err());
    }
    for e in [1.0 - 1e-12, 1.0, 1.0 + 1e-12] {
        let a = if e > 1.0 { -7000.0 } else { 7000.0 };
        let el = KeplerianElements::<Frame>::new(
            Kilometers::new(a),
            Eccentricity::new(e).unwrap(),
            Radians::new(0.0),
            Radians::new(0.0),
            Radians::new(0.0),
            Radians::new(0.0),
        )
        .unwrap();
        assert_eq!(el.to_geometry().is_err(), e == 1.0);
    }
}

#[test]
fn cartesian_rotation_matches_independent_classical_reference() {
    const MU: f64 = 398_600.441_8;

    // The first case is retrograde and elliptic; the second is hyperbolic.
    // The reference below expands Rz(raan) Rx(inclination) Rz(arg_periapsis)
    // directly instead of using affn's Euler-angle implementation.
    for (a, e, inclination, raan, arg_periapsis, true_anomaly) in [
        (12_000.0, 0.35, 2.2, 1.1, 0.65, 1.2),
        (-18_000.0, 1.6, 0.8, 2.4, 0.35, 0.7),
    ] {
        let elements = KeplerianElements::<Frame>::new(
            Kilometers::new(a),
            Eccentricity::new(e).unwrap(),
            Radians::new(inclination),
            Radians::new(raan),
            Radians::new(arg_periapsis),
            Radians::new(true_anomaly),
        )
        .unwrap();
        let state = elements
            .try_to_cartesian::<TestCenter>(qtty::dynamics::GravitationalParameter::new(MU))
            .unwrap();

        let p = a * (1.0 - e * e);
        let denominator = 1.0 + e * true_anomaly.cos();
        let radius = p / denominator;
        let root = (MU / p).sqrt();
        let r_pqw = [
            radius * true_anomaly.cos(),
            radius * true_anomaly.sin(),
            0.0,
        ];
        let v_pqw = [
            -root * true_anomaly.sin(),
            root * (e + true_anomaly.cos()),
            0.0,
        ];
        let (sw, cw) = arg_periapsis.sin_cos();
        let (si, ci) = inclination.sin_cos();
        let (so, co) = raan.sin_cos();
        let rotate = |vector: [f64; 3]| {
            let x = cw * vector[0] - sw * vector[1];
            let y = sw * vector[0] + cw * vector[1];
            let z = 0.0;
            [
                co * x - so * (ci * y - si * z),
                so * x + co * (ci * y - si * z),
                si * y + ci * z,
            ]
        };
        let expected_position = rotate(r_pqw);
        let expected_velocity = rotate(v_pqw);

        let position = state.position();
        let velocity = state.velocity();
        for (actual, expected) in [
            (position.x().value(), expected_position[0]),
            (position.y().value(), expected_position[1]),
            (position.z().value(), expected_position[2]),
        ] {
            assert!(
                (actual - expected).abs() < 1e-8,
                "position mismatch: actual={actual}, expected={expected}"
            );
        }
        for (actual, expected) in [
            (velocity.x().value(), expected_velocity[0]),
            (velocity.y().value(), expected_velocity[1]),
            (velocity.z().value(), expected_velocity[2]),
        ] {
            assert!(
                (actual - expected).abs() < 1e-11,
                "velocity mismatch: actual={actual}, expected={expected}"
            );
        }

        let position_norm = position
            .x()
            .value()
            .hypot(position.y().value())
            .hypot(position.z().value());
        let velocity_norm = velocity
            .x()
            .value()
            .hypot(velocity.y().value())
            .hypot(velocity.z().value());
        assert!((position_norm - radius).abs() < 1e-8);
        assert!((velocity_norm * velocity_norm - (2.0 * MU / radius - MU / a)).abs() < 1e-12);

        let angular_momentum = [
            position.y().value() * velocity.z().value()
                - position.z().value() * velocity.y().value(),
            position.z().value() * velocity.x().value()
                - position.x().value() * velocity.z().value(),
            position.x().value() * velocity.y().value()
                - position.y().value() * velocity.x().value(),
        ];
        let expected_angular_momentum = [si * so, -si * co, ci];
        let angular_momentum_norm = angular_momentum
            .iter()
            .map(|component| component * component)
            .sum::<f64>()
            .sqrt();
        for (actual, expected) in angular_momentum
            .iter()
            .map(|component| component / angular_momentum_norm)
            .zip(expected_angular_momentum)
        {
            assert!((actual - expected).abs() < 1e-12);
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct TestCenter;
impl affn::centers::ReferenceCenter for TestCenter {
    type Params = ();

    fn center_name() -> &'static str {
        "TestCenter"
    }
}
