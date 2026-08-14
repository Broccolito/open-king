//! KING relatedness estimators.
//!
//! Clean-room implementation of the estimators described in Manichaikul et al.
//! 2010 (Bioinformatics 26:2867-2873) and the KING 2.x documentation.
//!
//! Filled in by the implementation phase; see `docs/SPEC.md` section 4.

#![forbid(unsafe_code)]

pub mod counts;
pub mod kinship;
pub mod infer;
