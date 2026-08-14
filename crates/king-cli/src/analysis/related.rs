//! `--related`: the close-relative pass, and the small-sample downgrade that replaces it.
//!
//! `--related` is **not** a synonym for `--kinship`. On a fileset of ten or more samples
//! it writes a sixteen-column `.kin` and a fourteen-column `.kin0` whose extra columns
//! (`HetConc`, `HomIBS0`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`) come from the
//! IBD-segment engine, plus `<prefix>allsegs.txt`. That path is unimplemented; see
//! `docs/PARITY.md` §11.
//!
//! Below ten samples the reference downgrades the whole pass to `--kinship`, and that
//! path *is* implemented here. The downgrade is total: the ten-column `.kin` and
//! eight-column `.kin0` come out byte-identical to a plain `--kinship` run on the same
//! fileset, no `allsegs.txt` is written, and `--degree` is discarded — echoed nowhere and
//! applied nowhere.

use crate::cli::{Opt, Options};

/// Fewest samples the full `--related` pass will run on.
///
/// Ten. Established by a ladder of filesets: nine samples print the replacement notice
/// and emit the ten-column `.kin`, ten samples run the real pass and emit the sixteen-
/// column one. The corpus agrees — `dups` and `sexchr` (ten samples each) take the full
/// path while `missing` and `nuclear` (six) do not.
const MIN_SAMPLES: usize = 10;

/// Whether this sample count sends `--related` down the `--kinship` path.
pub fn downgrades_to_kinship(n_samples: usize) -> bool {
    n_samples < MIN_SAMPLES
}

/// `--related is replaced with --kinship for a small sample size.`, with the blank line
/// the reference puts above it.
///
/// Note the wording differs from `--ibdseg`'s own downgrade notice
/// (`--kinship analysis carried out instead for such a small sample size.`) — a
/// `--related --ibdseg` run on a three-sample fileset prints both, in that order, and
/// runs the kinship pass twice.
pub fn small_sample_notice() -> String {
    "\n--related is replaced with --kinship for a small sample size.\n".to_string()
}

/// The options the full pass echoes under `Options in effect:`.
///
/// `--related` then `--degree <d>`, matching the `core/*__related_degree*` captures. Not
/// reachable yet — the downgrade path echoes a bare `--kinship` instead — but it is the
/// shape the unimplemented pass owes.
pub fn options_in_effect(opts: &Options) -> Vec<String> {
    let mut v = vec!["--related".to_string()];
    let degree = opts.int(Opt::Degree);
    if degree != 0 {
        v.push(format!("--degree {degree}"));
    }
    let cpus = opts.int(Opt::Cpus);
    if cpus != 0 {
        v.push(format!("--cpus {cpus}"));
    }
    let prefix = opts.string(Opt::Prefix);
    if prefix != "king" {
        v.push(format!("--prefix {prefix}"));
    }
    v
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli;

    fn parse(args: &[&str]) -> Options {
        let owned: Vec<String> = args.iter().map(|s| (*s).to_string()).collect();
        cli::parse(&owned).options
    }

    #[test]
    fn the_downgrade_boundary_is_ten_samples() {
        assert!(downgrades_to_kinship(1));
        assert!(downgrades_to_kinship(9));
        assert!(!downgrades_to_kinship(10));
        assert!(!downgrades_to_kinship(200));
    }

    #[test]
    fn the_notice_carries_its_own_leading_blank_line() {
        assert_eq!(
            small_sample_notice(),
            "\n--related is replaced with --kinship for a small sample size.\n"
        );
    }

    #[test]
    fn the_downgrade_discards_degree() {
        // `--ibdseg --degree 2` on `singleton` prints the *unfiltered* between-family
        // line; the filtered one is what `--degree` would have produced.
        let opts = parse(&["--related", "--degree", "2"]);
        assert_eq!(opts.int(Opt::Degree), 2);
        assert_eq!(opts.without_degree().int(Opt::Degree), 0);
    }

    #[test]
    fn the_full_pass_echoes_related_then_degree() {
        assert_eq!(
            options_in_effect(&parse(&["--related", "--degree", "3"])),
            ["--related", "--degree 3"]
        );
        assert_eq!(options_in_effect(&parse(&["--related"])), ["--related"]);
    }
}
