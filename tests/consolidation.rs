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
