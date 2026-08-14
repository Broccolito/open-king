//! Byte-for-byte parity against console output captured from KING 2.3.2.
//!
//! Each fixture in `tests/fixtures/` is the verbatim stdout of the reference binary for
//! one command line. The tests run our binary with the same arguments and compare the
//! bytes. The **only** normalisation is on `KING starts at` / `KING ends at` lines, whose
//! content is a wall-clock timestamp; everything else — including the BEL that opens the
//! WARNING block and the tab/space mix in the closing citation — must match exactly.
//!
//! Captured with:
//!
//! ```text
//! king                              > noargs.stdout
//! king --related                    > related.stdout
//! king -b nonexistent.bed --related > bednone.stdout
//! king --bogusflag                  > bogus.stdout
//! ```
//!
//! All four exit 1 and write nothing to stderr.

use std::process::{Command, Output};

const NOARGS: &str = include_str!("fixtures/noargs.stdout.txt");
const RELATED: &str = include_str!("fixtures/related.stdout.txt");
const BEDNONE: &str = include_str!("fixtures/bednone.stdout.txt");
const BOGUS: &str = include_str!("fixtures/bogus.stdout.txt");

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_king"))
        .args(args)
        .output()
        .expect("king binary runs")
}

/// Replace the wall-clock part of the two timestamp lines, and nothing else.
///
/// Splitting on `\n` rather than using `lines()` keeps trailing newlines significant —
/// the captures end with a blank line and that has to be compared too.
fn normalise(s: &str) -> String {
    s.split('\n')
        .map(|line| {
            for prefix in ["KING starts at ", "KING ends at "] {
                if let Some(rest) = line.strip_prefix(prefix) {
                    assert_eq!(rest.len(), 24, "ctime is 24 characters: {rest:?}");
                    return format!("{prefix}<TIME>");
                }
            }
            line.to_string()
        })
        .collect::<Vec<_>>()
        .join("\n")
}

#[track_caller]
fn assert_parity(args: &[&str], expected: &str) {
    let out = run(args);
    let stdout = String::from_utf8(out.stdout).expect("stdout is UTF-8");
    assert!(
        out.stderr.is_empty(),
        "the reference writes nothing to stderr for `king {}`",
        args.join(" ")
    );
    assert_eq!(
        out.status.code(),
        Some(1),
        "the reference exits 1 for `king {}`",
        args.join(" ")
    );
    assert_eq!(normalise(&stdout), normalise(expected));
}

#[test]
fn no_arguments() {
    assert_parity(&[], NOARGS);
}

#[test]
fn related_without_a_fileset() {
    assert_parity(&["--related"], RELATED);
}

#[test]
fn related_with_a_missing_bed() {
    assert_parity(&["-b", "nonexistent.bed", "--related"], BEDNONE);
}

#[test]
fn undefined_option() {
    assert_parity(&["--bogusflag"], BOGUS);
}

/// [`normalise`] must touch the timestamp lines and nothing else, or the four tests
/// above would be quietly comparing less than they claim.
#[test]
fn normalisation_touches_only_timestamps() {
    assert_eq!(
        normalise("KING ends at Thu Aug 13 17:29:25 2026"),
        "KING ends at <TIME>"
    );
    // These three fixtures contain no timestamp at all, so normalising is the identity —
    // trailing blank lines and the BEL included.
    for fixture in [NOARGS, RELATED, BOGUS] {
        assert_eq!(normalise(fixture), fixture);
    }
    assert!(
        BOGUS.contains('\u{7}'),
        "the WARNING block opens with a BEL"
    );
    assert_eq!(BEDNONE.matches("KING starts at").count(), 1);
}
