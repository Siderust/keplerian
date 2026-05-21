# keplerian

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
| `std` | yes | Standard-library support via `qtty`. |
| `alloc` | no | Enables `Vec`-backed search grids. |
| `serde` | no | Serde derives for public data types. |

## License

AGPL-3.0-only.

[`qtty`]: https://crates.io/crates/qtty
[`affn`]: https://crates.io/crates/affn
