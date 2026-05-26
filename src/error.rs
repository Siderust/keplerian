// SPDX-License-Identifier: AGPL-3.0-only
// Copyright (C) 2026 Vallés Puig, Ramon

//! Top-level error aggregation for crate workflows.
//!
//! ## Scientific scope
//! This module does not add new orbital mechanics. It only bundles the error
//! taxonomies of anomaly conversion, element conversion, propagation, and
//! Lambert solving.
//!
//! ## Technical scope
//! [`KeplerError`] is a convenience enum for callers that want one top-level
//! error type across multiple subsystems.
//!
//! ## References
//! - This module aggregates the references documented by the submodules it wraps.

use crate::{anomaly, elements, lambert, problem};

/// Unified error type for high-level Keplerian workflows.
#[derive(Debug, thiserror::Error)]
pub enum KeplerError {
    /// An anomaly conversion or solver failed.
    #[error(transparent)]
    Anomaly(anomaly::AnomalyError),
    /// Cartesian/element conversion failed.
    #[error(transparent)]
    Conversion(elements::ConversionError),
    /// Two-body propagation failed.
    #[error(transparent)]
    Propagation(problem::PropagationError),
    /// Lambert solving failed.
    #[error(transparent)]
    Lambert(lambert::LambertError),
}

impl From<anomaly::AnomalyError> for KeplerError {
    fn from(value: anomaly::AnomalyError) -> Self {
        Self::Anomaly(value)
    }
}

impl From<elements::ConversionError> for KeplerError {
    fn from(value: elements::ConversionError) -> Self {
        Self::Conversion(value)
    }
}

impl From<problem::PropagationError> for KeplerError {
    fn from(value: problem::PropagationError) -> Self {
        Self::Propagation(value)
    }
}

impl From<lambert::LambertError> for KeplerError {
    fn from(value: lambert::LambertError) -> Self {
        Self::Lambert(value)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kepler_error_from_anomaly() {
        let inner = anomaly::AnomalyError::InvalidEccentricity(1.5);
        let e = KeplerError::from(inner);
        assert!(matches!(e, KeplerError::Anomaly(_)));
    }

    #[test]
    fn kepler_error_from_conversion() {
        let inner = elements::ConversionError::InvalidEccentricity(-1.0);
        let e = KeplerError::from(inner);
        assert!(matches!(e, KeplerError::Conversion(_)));
    }

    #[test]
    fn kepler_error_from_propagation() {
        let inner = problem::PropagationError::ParabolicUnsupported;
        let e = KeplerError::from(inner);
        assert!(matches!(e, KeplerError::Propagation(_)));
    }

    #[test]
    fn kepler_error_from_lambert() {
        let inner = lambert::LambertError::ZeroPosition;
        let e = KeplerError::from(inner);
        assert!(matches!(e, KeplerError::Lambert(_)));
    }
}
