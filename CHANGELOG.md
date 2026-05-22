# Changelog

All notable changes to this project will be documented in this file.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2026-05-22

### Added

- `Eccentricity` newtype with `new_unchecked`, `new_elliptic`, `new_hyperbolic`,
  and `classify` constructors; compile-time validation enforced.
- `AnomalyOptions` control struct (max iterations, tolerance) and
  `AnomalyError` (non-convergence) for all iterative solvers.
- Elliptic anomaly solvers: `eccentric_from_mean`, `true_from_eccentric`,
  `eccentric_from_true`, `mean_from_eccentric`.
- Parabolic anomaly solver: `kepler_parabolic` (Barker's equation).
- Hyperbolic anomaly solvers: `hyperbolic_from_mean`, `true_from_hyperbolic`.
- `KeplerianElements<F>` — six classical elements typed over an `affn` frame;
  conversions `to_cartesian` / `from_cartesian`.
- `KeplerProblem<C, F>` — two-body IVP solver; `new` + `propagate`.
- Lambert boundary-value solver (`lambert`, `lambert_n_rev`) with multi-rev
  support.
- Transfer invariants: `specific_orbital_energy`, `specific_angular_momentum`,
  `vis_viva_speed`, `orbital_period`, `escape_speed`, `hohmann_delta_v`.
- `OrbitalState` and `StateDerivative` Cartesian state types.
- `KeplerError` crate-level error family.
- `alloc`-gated `search` module: Lambert transfer-search grids.
- Optional `serde` feature for all public data types.
- CI workflow with fmt, Clippy (default + all features), no-std/alloc checks,
  tests (default + all features), and doc-test jobs.
- Audit, deny, coverage (llvm-cov ≥ 90 % line rate), and Miri jobs.
- Publish workflow triggered on `v*.*.*` tags.

[0.1.0]: https://github.com/siderust/keplerian/releases/tag/v0.1.0
