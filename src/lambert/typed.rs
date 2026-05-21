//! Typed entry-points for Lambert's problem.
//!
//! ## Scientific scope
//! This module solves Lambert's two-point boundary-value problem between two
//! positions in a fixed frame and central field.
//!
//! ## Technical scope
//! Public functions consume typed `affn` positions, typed `qtty::Second` time
//! of flight, and a typed gravitational parameter. The raw numeric backend is
//! kept crate-private.
//!
//! ## References
//! - Izzo, D. (2014). *Revisiting Lambert's Problem*.

use affn::cartesian::{Position, Velocity};
use affn::centers::ReferenceCenter;
use affn::frames::ReferenceFrame;
use qtty::dynamics::{GravitationalParameter, KmPerSecond};
use qtty::length::Kilometer;
use qtty::Second;

use super::error::LambertError;
use super::izzo::{
    solve_lambert as solve_lambert_arr, solve_lambert_n_rev as solve_lambert_n_rev_arr,
    LambertBranch, LambertDiagnostics, LambertSolution, NRevBranch,
};

/// Typed Lambert solution — departure / arrival velocities plus diagnostics.
#[derive(Debug, Clone, Copy)]
pub struct TypedLambertSolution<F: ReferenceFrame> {
    /// Departure velocity at `r1`, in the same frame as the input positions.
    pub v1: Velocity<F, KmPerSecond>,
    /// Arrival velocity at `r2`, in the same frame as the input positions.
    pub v2: Velocity<F, KmPerSecond>,
    /// Householder iteration diagnostics from the underlying numeric solver.
    pub diagnostics: LambertDiagnostics,
}

/// Solve Lambert's problem (single revolution) on typed inputs.
///
/// # Errors
///
/// Returns [`LambertError`] if the solver fails to converge or the geometry is
/// degenerate (e.g., anti-parallel positions).
///
/// # Examples
///
/// ```
/// use keplerian::lambert::{lambert, LambertBranch};
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::Second;
/// use qtty::length::Kilometer;
/// use affn::cartesian::Position;
///
/// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
/// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
///
/// let r1 = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
/// let r2 = Position::<C, F, Kilometer>::new(0.0, 7000.0, 0.0);
/// let tof = Second::new(3600.0);
/// let mu = GravitationalParameter::new(398600.4418);
/// let sol = lambert(r1, r2, tof, mu, LambertBranch::Prograde).unwrap();
/// assert!(sol.v1.x().value().is_finite());
/// ```
pub fn lambert<C, F>(
    r1: Position<C, F, Kilometer>,
    r2: Position<C, F, Kilometer>,
    tof: Second,
    mu: GravitationalParameter,
    branch: LambertBranch,
) -> Result<TypedLambertSolution<F>, LambertError>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let LambertSolution {
        v1,
        v2,
        diagnostics,
    } = solve_lambert_arr(
        position_to_array(&r1),
        position_to_array(&r2),
        tof.value(),
        mu.value(),
        branch,
    )?;

    Ok(TypedLambertSolution {
        v1: Velocity::<F, KmPerSecond>::new(v1[0], v1[1], v1[2]),
        v2: Velocity::<F, KmPerSecond>::new(v2[0], v2[1], v2[2]),
        diagnostics,
    })
}

/// Solve Lambert's problem with `N ≥ 1` complete revolutions, on typed inputs.
///
/// # Errors
///
/// Returns [`LambertError`] if the solver fails to converge, `revolutions == 0`,
/// or the geometry is degenerate.
///
/// # Examples
///
/// ```
/// use keplerian::lambert::{lambert_n_rev, LambertBranch, NRevBranch};
/// use qtty::dynamics::GravitationalParameter;
/// use qtty::Second;
/// use qtty::length::Kilometer;
/// use affn::cartesian::Position;
///
/// #[derive(Debug, Clone, Copy)] struct C; impl affn::centers::ReferenceCenter for C { type Params = (); fn center_name() -> &'static str { "C" } }
/// #[derive(Debug, Clone, Copy)] struct F; impl affn::frames::ReferenceFrame for F { fn frame_name() -> &'static str { "F" } }
///
/// let r1 = Position::<C, F, Kilometer>::new(7000.0, 0.0, 0.0);
/// let r2 = Position::<C, F, Kilometer>::new(0.0, 7000.0, 0.0);
/// let tof = Second::new(7200.0);
/// let mu = GravitationalParameter::new(398600.4418);
/// // Multi-rev may fail if TOF is too short; just check it compiles/runs.
/// let _ = lambert_n_rev(r1, r2, tof, mu, LambertBranch::Prograde, 1, NRevBranch::Left);
/// ```
pub fn lambert_n_rev<C, F>(
    r1: Position<C, F, Kilometer>,
    r2: Position<C, F, Kilometer>,
    tof: Second,
    mu: GravitationalParameter,
    branch: LambertBranch,
    revolutions: u32,
    side: NRevBranch,
) -> Result<TypedLambertSolution<F>, LambertError>
where
    C: ReferenceCenter<Params = ()>,
    F: ReferenceFrame,
{
    let LambertSolution {
        v1,
        v2,
        diagnostics,
    } = solve_lambert_n_rev_arr(
        position_to_array(&r1),
        position_to_array(&r2),
        tof.value(),
        mu.value(),
        branch,
        revolutions,
        side,
    )?;

    Ok(TypedLambertSolution {
        v1: Velocity::<F, KmPerSecond>::new(v1[0], v1[1], v1[2]),
        v2: Velocity::<F, KmPerSecond>::new(v2[0], v2[1], v2[2]),
        diagnostics,
    })
}

fn position_to_array<C, F>(p: &Position<C, F, Kilometer>) -> [f64; 3]
where
    C: ReferenceCenter,
    F: ReferenceFrame,
{
    [p.x().value(), p.y().value(), p.z().value()]
}

#[cfg(test)]
mod tests {
    use super::*;
    use affn::frames::ICRS;

    #[test]
    fn typed_zero_rev_matches_array_kernel() {
        let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
        let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
        let tof = Second::new(4_560.0);
        let mu = GravitationalParameter::new(398_600.441_8);

        let typed = lambert(r1, r2, tof, mu, LambertBranch::Prograde).unwrap();
        let raw = solve_lambert_arr(
            [15945.34, 0.0, 0.0],
            [12214.83899, 10249.46731, 0.0],
            4_560.0,
            398_600.441_8,
            LambertBranch::Prograde,
        )
        .unwrap();
        for i in 0..3 {
            assert!((typed.v1.as_array()[i].value() - raw.v1[i]).abs() < 1e-12);
            assert!((typed.v2.as_array()[i].value() - raw.v2[i]).abs() < 1e-12);
        }
    }

    #[test]
    fn typed_n_rev_propagates_error() {
        let r1 = Position::<(), ICRS, Kilometer>::new(15945.34, 0.0, 0.0);
        let r2 = Position::<(), ICRS, Kilometer>::new(12214.83899, 10249.46731, 0.0);
        let tof = Second::new(4_560.0);
        let mu = GravitationalParameter::new(398_600.441_8);
        let err = lambert_n_rev(
            r1,
            r2,
            tof,
            mu,
            LambertBranch::Prograde,
            1,
            NRevBranch::Left,
        )
        .unwrap_err();
        match err {
            LambertError::RevolutionsExceedNMax { requested, .. } => {
                assert_eq!(requested, 1);
            }
            other => panic!("unexpected error: {other:?}"),
        }
    }
}
