//! Safety contract for deliberately excluded product-scope requests.
//!
//! These cases intentionally depart from KING's exit-0 behavior. A recognized analysis
//! must never appear to succeed after producing none of the files its name promises.

use std::process::{Command, Output};

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_open-king"))
        .args(args)
        .output()
        .expect("open-king binary runs")
}

#[track_caller]
fn assert_rejected(args: &[&str], requests: &str) {
    let out = run(args);
    let stdout = String::from_utf8(out.stdout).expect("stdout is UTF-8");
    assert_eq!(out.status.code(), Some(1));
    assert!(out.stderr.is_empty(), "all diagnostics belong on stdout");
    assert!(
        stdout.contains(&format!(
            "open-king's minimal relatedness product does not implement: {requests}."
        )),
        "missing product-scope diagnostic in:\n{stdout}"
    );
    assert!(stdout.contains("Supported analyses: --related"));
    assert!(stdout.contains("See docs/SCOPE.md"));
    assert!(
        !stdout.contains("KING starts at"),
        "scope validation must happen before analysis startup"
    );
    assert!(
        !stdout.contains("cannot be opened"),
        "scope validation must happen before probing input files"
    );
}

#[test]
fn excluded_analysis_is_fatal_before_input_is_opened() {
    assert_rejected(&["-b", "does-not-exist.bed", "--pca"], "--pca");
}

#[test]
fn excluded_analysis_cannot_hide_behind_a_supported_analysis() {
    assert_rejected(&["-b", "does-not-exist.bed", "--kinship", "--roh"], "--roh");
}

#[test]
fn excluded_parameters_and_output_modes_are_fatal_together() {
    assert_rejected(
        &[
            "-b",
            "does-not-exist.bed",
            "--kinship",
            "--plink",
            "--rpath",
            "/tmp/R",
            "--pcs",
            "0",
        ],
        "--plink, --rpath, --pcs",
    );
}

#[test]
fn comma_separated_multifileset_input_is_fatal_before_file_lookup() {
    assert_rejected(
        &["-b", "a.bed,b.bed", "--kinship"],
        "comma-separated multi-fileset input",
    );
}
