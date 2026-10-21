# keplerian

[![Crates.io](https://img.shields.io/crates/v/keplerian)](https://crates.io/crates/keplerian)
[![docs.rs](https://img.shields.io/docsrs/keplerian)](https://docs.rs/keplerian)
[![CI](https://github.com/Siderust/keplerian/actions/workflows/ci.yml/badge.svg)](https://github.com/Siderust/keplerian/actions/workflows/ci.yml)
[![License: AGPL-3.0-only](https://img.shields.io/badge/license-AGPL--3.0--only-blue)](LICENSE)

Domain-agnostic Keplerian dynamics in Rust: typed anomaly solvers,
classical elements, two-body propagation, Lambert transfers, and
transfer/search helpers built only on [`qtty`] and [`affn`].

## Scope

This crate owns reusable **central-force / two-body** dynamics that do
not require astronomy-specific semantics.

In scope:

- elliptic, parabolic, and hyperbolic anomaly conversions and Kepler solves,
- typed Keplerian elements and typed Cartesian states,
- analytic two-body propagation under a fixed gravitational parameter,
- Lambert boundary-value solvers,
- transfer invariants and caller-driven Lambert search grids.

Out of scope:

- epochs, time scales, and calendars,
- ephemerides, bodies, observatories, or frame pipelines,
- non-Keplerian perturbation models.

## Crate boundary

```text
qtty   ─── typed quantities and units
affn   ─── typed geometry, frames, centers
   │
   └──> keplerian
```

## Features

| Feature | Default | Effect |
|---|:---:|---|
| `std` | yes | Standard library (implies `alloc`); forwards `qtty/std` and `affn/std`. |
| `alloc` | no | Heap-backed APIs such as `Vec` search grids; forwards `qtty/alloc` and `affn/alloc`. |
| `serde` | no | Serde derives for public data types (implies `alloc`). |

Pure `core`-only builds use `--no-default-features`. Combine with `--features alloc` and/or `--features serde` as needed.

## License

AGPL-3.0-only.

[`qtty`]: https://crates.io/crates/qtty
[`affn`]: https://crates.io/crates/affn
