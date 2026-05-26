# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.2.0] - 2026-05-26

### Added

- `MeanAnomaly`, `TrueAnomaly`, `EccentricAnomaly`, `HyperbolicAnomaly`,
  `ParabolicAnomaly` typed anomaly newtypes; all anomaly solver functions now
  accept and return these typed values instead of raw `Radians` or `f64`.
- `AnomalyOptions::try_new(max_iter, tol)` fallible constructor; fields are now
  private; `max_iter()` and `tol()` accessors added.
- `AnomalyError::InvalidOptions` variant for invalid solver configuration.
- `ConicRegime` and `EccentricityError` moved to `eccentricity` module;
  `elements::ConicRegime` is a re-export for backwards compatibility.
- `Eccentricity::try_new`, `new_elliptic`, `new_hyperbolic`, `parabolic`,
  `classify` constructors.
- `parabolic_from_mean`, `true_from_parabolic`, `parabolic_from_true` parabolic
  anomaly helpers (replaces the old `kepler_parabolic` raw-f64 function).
- `KeplerianElements::try_to_cartesian` fallible Cartesian conversion (replaces
  infallible `to_cartesian`).
- `ConversionError::IncoherentRegime` — returned when the semi-major axis sign
  is inconsistent with the eccentricity regime.
- `a`/`e` coherence check in `KeplerianElements::new`: elliptic requires `a > 0`;
  hyperbolic requires `a < 0`.
- `TransferError` enum for fallible transfer helpers.
- `try_orbital_period`, `try_hohmann_delta_v`, `try_vis_viva_speed`,
  `try_escape_speed` fallible variants of transfer functions.
- `prelude` module re-exporting the most commonly used public items.
- `#![cfg_attr(not(feature = "std"), no_std)]` — crate is now no-std capable
  with `alloc` feature for heap-using APIs.
- `Clone` derive on `LambertError`; `PartialEq` and serde derives on
  `LambertDiagnostics`, `LambertBranch`, `NRevBranch`, `TypedLambertSolution`.
- CI: additional checks for `--features serde` and `--features alloc,serde`
  no-default-features combinations.

### Changed

- `KeplerProblem.mu` field is now private; use `mu()` accessor.
- `TransferCandidate.total_dv` renamed to `endpoint_speed_sum`.
- `kepler_parabolic` renamed to `parabolic_from_mean`; accepts and returns
  `ParabolicAnomaly` instead of raw `f64`.

### Removed

- Infallible `KeplerianElements::to_cartesian`; replaced by
  `try_to_cartesian`.
- `AnomalyOptions` struct literal construction (fields privatised); use
  `AnomalyOptions::try_new`.

## [0.1.0] - 2026-05-22

### Added

- `Eccentricity` newtype with `new_unchecked` constructor; compile-time
  validation enforced.
- `AnomalyOptions` control struct and `AnomalyError` for all iterative solvers.
- Elliptic anomaly solvers: `eccentric_from_mean`, `true_from_eccentric`,
  `eccentric_from_true`, `mean_from_eccentric`.
- Hyperbolic anomaly solvers: `hyperbolic_from_mean`, `true_from_hyperbolic`.
- `KeplerianElements<F>` — six classical elements typed over an `affn` frame;
  conversions `from_cartesian`.
- `KeplerProblem<C, F>` — two-body IVP solver; `new` + `propagate`.
- Lambert boundary-value solver (`lambert`, `lambert_n_rev`) with multi-rev
  support.
- Transfer invariants: `specific_orbital_energy`, `specific_angular_momentum`,
  `vis_viva_speed`, `orbital_period`, `escape_speed`, `hohmann_delta_v`.
- `CartesianState` Cartesian state type.
- `KeplerError` crate-level error family.
- `alloc`-gated `search` module: Lambert transfer-search grids.
- Optional `serde` feature for public data types.
- CI workflow with fmt, Clippy, no-std/alloc checks, tests, and doc-test jobs.
- Audit, deny, coverage (llvm-cov ≥ 90 % line rate), and Miri jobs.
- Publish workflow triggered on `v*.*.*` tags.

[0.2.0]: https://github.com/siderust/keplerian/releases/tag/v0.2.0
[0.1.0]: https://github.com/siderust/keplerian/releases/tag/v0.1.0
