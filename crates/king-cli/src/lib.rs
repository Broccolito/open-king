//! Command line and console layer for `king`.
//!
//! Two modules, both reproducing the reference binary's behaviour byte for byte:
//!
//! * [`cli`] — the argument parser and the option table.
//! * [`console`] — every line the program prints, from the banner to the run-time
//!   progress ticks the analysis engines emit.
//!
//! The binary in `main.rs` is a thin wiring layer over these.

#![forbid(unsafe_code)]

pub mod cli;
pub mod console;
