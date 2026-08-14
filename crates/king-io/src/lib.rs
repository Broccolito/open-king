//! PLINK 1 binary fileset I/O.
//!
//! Reads `.bed` / `.bim` / `.fam` triples into the packed, SNP-major bit-plane
//! representation that the relatedness kernels in `king-core` consume.
//!
//! Filled in by the implementation phase; see `docs/SPEC.md` section 3.

#![forbid(unsafe_code)]

pub mod bed;
pub mod bim;
pub mod fam;
