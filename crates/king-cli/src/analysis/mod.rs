//! The analysis engines, one module per `--flag`, plus the pieces they share.
//!
//! Every analysis follows the same shape, which is the reference's own:
//!
//! ```text
//! [Autosome genotypes stored in <w> words for each of <n> individuals.]   (load::prints_preamble)
//! Options in effect:
//! \t--<analysis>
//! \t[--degree <d>] [--minConc <c>] [--cpus <n>] [--prefix <p>]
//!
//! <the analysis body>
//! ```
//!
//! and `KING ends at <time>` closes the run once, after the last pass.
//!
//! This module owns only what more than one analysis needs: the sample-ID comparator,
//! the three row orders, the `Options in effect` list and the C-compatible number
//! formatting. Anything used by exactly one analysis lives in that analysis's module.

pub mod autoqc;
pub mod build;
pub mod cluster;
pub mod duplicate;
pub mod ibdseg;
pub mod ibs;
pub mod kinship;
pub mod qc;
pub mod related;
pub mod segments;
pub mod splitped;
pub mod unrelated;
pub mod xkinship;

use std::cmp::Ordering;
use std::collections::HashSet;

use king_io::Sample;

use crate::cli::{Kind, Opt, Options};

// ---------------------------------------------------------------------------
// Sample identifiers
// ---------------------------------------------------------------------------

/// The comparator the reference orders `.kin`/`.ibs` FID blocks and family members by.
///
/// Fitted against 44 probe families in `docs/BEHAVIOR.md` §Q6, which is also where the
/// raw input→output orderings live. Four rules, applied while walking the two IDs
/// together:
///
/// * a **run of digits** compares against another digit run by **length first, then
///   bytes**, so `7 < 70 < 007` — no integer is ever parsed, which is why 26-digit IDs
///   order correctly and there is no overflow behaviour to emulate;
/// * a **non-digit** compares one character at a time, ASCII-uppercase folded, so
///   `z1 < [1` (`Z` is 0x5A, `[` is 0x5B) and `ab < aB` compare equal;
/// * a non-digit sorts **before** a digit: `b1 < 1a`, `_1 < 1`;
/// * a shorter ID that is a prefix of a longer one sorts first: `a < a1`, `x < x0`.
///
/// Plain lexicographic order and plain "natural sort" each get part of this wrong; both
/// were tried against the reference and rejected.
pub fn king_id_cmp(a: &[u8], b: &[u8]) -> Ordering {
    let (mut i, mut j) = (0usize, 0usize);
    loop {
        match (i == a.len(), j == b.len()) {
            (true, true) => return Ordering::Equal,
            (true, false) => return Ordering::Less,
            (false, true) => return Ordering::Greater,
            _ => {}
        }
        let (ca, cb) = (a[i], b[j]);
        let (da, db) = (ca.is_ascii_digit(), cb.is_ascii_digit());
        if da != db {
            return if da {
                Ordering::Greater
            } else {
                Ordering::Less
            };
        }
        if da {
            let ra = digit_run(&a[i..]);
            let rb = digit_run(&b[j..]);
            match ra.len().cmp(&rb.len()).then_with(|| ra.cmp(rb)) {
                Ordering::Equal => {
                    i += ra.len();
                    j += rb.len();
                }
                ord => return ord,
            }
        } else {
            match ca.to_ascii_uppercase().cmp(&cb.to_ascii_uppercase()) {
                Ordering::Equal => {
                    i += 1;
                    j += 1;
                }
                ord => return ord,
            }
        }
    }
}

/// The leading run of ASCII digits.
fn digit_run(s: &[u8]) -> &[u8] {
    let n = s.iter().take_while(|c| c.is_ascii_digit()).count();
    &s[..n]
}

// ---------------------------------------------------------------------------
// Row orders
// ---------------------------------------------------------------------------

/// The families of a sample set, ordered the way the *sorted* writers emit them.
///
/// Families are ordered by [`king_id_cmp`] of the FID and their members by
/// [`king_id_cmp`] of the IID — **both independent of `.fam` order**, which is the trap:
/// the same `.bed` analysed with its rows shuffled emits identical `.kin`/`.ibs`.
pub fn family_blocks(samples: &[Sample]) -> Vec<Vec<usize>> {
    let mut fids: Vec<&str> = samples.iter().map(|s| s.fid.as_str()).collect();
    fids.sort_by(|a, b| king_id_cmp(a.as_bytes(), b.as_bytes()));
    fids.dedup_by(|a, b| king_id_cmp(a.as_bytes(), b.as_bytes()) == Ordering::Equal);

    fids.into_iter()
        .map(|fid| {
            let mut members: Vec<usize> = samples
                .iter()
                .enumerate()
                .filter(|(_, s)| king_id_cmp(s.fid.as_bytes(), fid.as_bytes()) == Ordering::Equal)
                .map(|(i, _)| i)
                .collect();
            members.sort_by(|&x, &y| {
                king_id_cmp(samples[x].iid.as_bytes(), samples[y].iid.as_bytes())
            });
            members
        })
        .collect()
}

/// Add a founder row for every `PAT`/`MAT` that names someone the `.fam` does not list.
///
/// The reference resolves parent references **by string, within the family**, and
/// materialises the ones it cannot find: two rows of family `F` whose `PAT`/`MAT` are
/// `PA`/`MA` are full siblings even though neither parent is genotyped, and the reference
/// duly prints `Z0 0.250  Phi 0.2500` for them. It does this even when an id of that name
/// exists in a *different* family.
///
/// `king_core::infer::Pedigree` treats an unresolvable parent as unknown, so without this
/// the `Phi`/`Z0` columns collapse to `0.0000`/`1.000` on any `.fam` whose founders are
/// not themselves genotyped — and, because the PO/FS cutoff is calibrated on declared
/// sibships, the summary's `PO` column collapses with them.
///
/// A phantom's **sex is inferred from the slot that named it**: a `PAT` reference makes a
/// male, a `MAT` reference a female. Only the X pass reads it — the X kinship recurrence
/// branches on sex — but it costs nothing to get right here. Phantoms are appended after
/// the real rows, so every real sample keeps its index.
pub fn with_phantom_parents(samples: &[Sample]) -> Vec<Sample> {
    let mut known: HashSet<(&str, &str)> = samples
        .iter()
        .map(|s| (s.fid.as_str(), s.iid.as_str()))
        .collect();
    let mut phantoms: Vec<Sample> = Vec::new();
    for s in samples {
        for (parent, sex) in [(&s.pat, 1), (&s.mat, 2)] {
            if parent == "0" || !known.insert((&s.fid, parent)) {
                continue;
            }
            phantoms.push(Sample {
                fid: s.fid.clone(),
                iid: parent.clone(),
                pat: "0".to_string(),
                mat: "0".to_string(),
                sex,
                pheno: "-9".to_string(),
            });
        }
    }
    let mut out = samples.to_vec();
    out.append(&mut phantoms);
    out
}

/// Within-family pairs in `.kin` / `.ibs` order: the `i < j` upper triangle of each
/// family's sorted member list, family blocks in sorted order.
pub fn within_family_pairs(samples: &[Sample]) -> Vec<(usize, usize)> {
    let mut out = Vec::new();
    for members in family_blocks(samples) {
        for (k, &a) in members.iter().enumerate() {
            for &b in &members[k + 1..] {
                out.push((a, b));
            }
        }
    }
    out
}

/// Between-family pairs in the **square-tiled** order the reference emits.
///
/// Pairs are `i < j` over `.fam` index order, but sorted by `(i / block, j / block, i, j)`
/// — the writer walks a block-tiled loop. `block` is **8** for `.ibs0` and 32 for
/// `.kin0`. Plain ascending `(i, j)` coincides with the tiled order whenever
/// `n <= block`, so it looks right on every small fixture and silently diverges on real
/// data; verified against `bigish` (200 samples, 19 327 rows) and `unrelated` (30
/// samples, 390 rows), both of which match `block = 8` and nothing else.
pub fn between_family_pairs(samples: &[Sample], block: usize) -> Vec<(usize, usize)> {
    let n = samples.len();
    let mut out = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            if samples[i].fid != samples[j].fid {
                out.push((i, j));
            }
        }
    }
    out.sort_by_key(|&(i, j)| (i / block, j / block, i, j));
    out
}

/// Every pair in `.fam` serial order — plain row-major `i < j`.
///
/// This is `.con`'s order, and it is **not** tiled: a 296-row `.con` from a 200-sample
/// run is in exact row-major order, which rules out the `.ibs0`/`.kin0` tiling.
pub fn serial_pairs(n: usize) -> Vec<(usize, usize)> {
    let mut out = Vec::with_capacity(n * n / 2);
    for i in 0..n {
        for j in i + 1..n {
            out.push((i, j));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Console / naming helpers
// ---------------------------------------------------------------------------

/// The options the reference echoes under `Options in effect:` for an analysis pass.
///
/// The analysis flag itself, then every **explicitly given** option from a short
/// whitelist, in option-table order. The whitelist is what the corpus shows: `--degree`,
/// `--minConc`, `--cpus`, `--prefix`. It is emphatically not "every option that was
/// given" — `--seglength 10` is accepted, changes the output, and is **not** echoed
/// (`--ibdseg --seglength 10` prints a bare `--ibdseg`), so the block is a fixed list
/// rather than a dump of the parse state.
pub fn options_in_effect(opts: &Options, analysis: Opt) -> Vec<String> {
    let mut out = vec![format!("--{}", crate::cli::echo_name(analysis))];
    for o in [Opt::Degree, Opt::MinConc, Opt::Cpus, Opt::Prefix] {
        let value = match o.kind() {
            // An integer option carries its own "unset": `--degree 0` and `--cpus 0` are
            // not echoed, exactly as an unmentioned one is not.
            Kind::Int => match opts.int(o) {
                0 => continue,
                v => v.to_string(),
            },
            // A double is echoed whenever it was explicitly given, `--minConc 0`
            // included, so the value alone cannot decide it. `%G`, so `1` prints as `1`
            // and `0.9` as `0.9` — never `1.000000`.
            Kind::Double if opts.was_given(o) => format!("{}", opts.double(o)),
            // `--prefix` is echoed only when it differs from the default it always holds.
            Kind::Str if opts.string(o) != Options::new().string(o) => opts.string(o).to_string(),
            _ => continue,
        };
        out.push(format!("--{} {}", o.name(), value));
    }
    out
}

/// An output file name: `--prefix` is plain **concatenation**, not a stem plus a dot.
///
/// `--prefix ZZ_` gives `ZZ_.ibs` and `ZZ_allsegs.txt`; the default prefix `king` gives
/// `king.ibs` and `kingallsegs.txt`. Note there is no dot before `allsegs.txt`, which is
/// why the suffix has to carry its own.
pub fn out_path(opts: &Options, suffix: &str) -> String {
    format!("{}{suffix}", opts.string(Opt::Prefix))
}

/// How many CPU cores to claim in the console lines.
///
/// The count is normalised away before any parity diff (it is a property of the host,
/// not of the data), so this only has to be a plausible positive integer; `--cpus` wins
/// when it was given, exactly as the reference's own line follows the flag.
pub fn cpu_count(opts: &Options) -> usize {
    let n = opts.int(Opt::Cpus);
    if n > 0 {
        return n as usize;
    }
    std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1)
}

/// `printf("%.*lf")` for a `double`, including the spellings C uses for non-finite
/// values.
///
/// Rust prints `NaN` and Rust's `{:.4}` of an infinity is `inf`; C prints `nan` and
/// `inf`/`-inf`. The reference reaches all three on degenerate pairs (an all-homozygous
/// sample gives `HomConc = nan`, a within-family pair with no heterozygotes and a
/// non-zero IBS0 gives `Kinship = -inf`), so the spelling matters.
pub fn f(x: f64, precision: usize) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-inf" } else { "inf" }.to_string();
    }
    format!("{x:.precision$}")
}

/// `printf("%G")` — the conversion the `Error` column and `MI_Removal` use.
///
/// Only ever reached with `0`, `0.5` and `1`, which it renders without trailing zeros.
pub fn g(x: f64) -> String {
    if x.is_nan() {
        return "nan".to_string();
    }
    if x.is_infinite() {
        return if x < 0.0 { "-INF" } else { "INF" }.to_string();
    }
    format!("{x}")
}

// ---------------------------------------------------------------------------
// Relationship bands
// ---------------------------------------------------------------------------

/// The kinship boundaries between inference classes.
///
/// The **exact** `2^-(k + 3/2)` grid, not the rounded constants the manual prints: a
/// designed pair at `0.1767775` — above `2^-2.5` but below the printed `0.17678` — is
/// treated as 1st degree by the reference.
pub mod band {
    /// At or above this: duplicate / MZ twin. `2^-1.5`.
    pub const MZ: f64 = 0.353_553_390_593_273_8;
    /// At or above this: 1st degree. `2^-2.5`.
    pub const FIRST: f64 = 0.176_776_695_296_636_9;
    /// At or above this: 2nd degree. `2^-3.5`.
    pub const SECOND: f64 = 0.088_388_347_648_318_45;
    /// At or above this: 3rd degree. `2^-4.5`.
    pub const THIRD: f64 = 0.044_194_173_824_159_22;
    /// At or above this: 4th degree; below it, unrelated. `2^-5.5`.
    pub const FOURTH: f64 = 0.022_097_086_912_079_61;
}

/// One column of the console's relationship-summary table.
///
/// `Other` is the table's `OTHER`: everything below 3rd degree, i.e. 4th-degree and
/// unrelated pairs lumped together, plus every pedigree relationship that is not one of
/// the five named classes.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    Mz,
    Po,
    Fs,
    Second,
    Third,
    Other,
}

impl Class {
    /// Whether the class counts towards the summary's `total relatives`, which is every
    /// column except `OTHER`.
    pub fn is_relative(self) -> bool {
        self != Class::Other
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(fid: &str, iid: &str) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: "0".to_string(),
            mat: "0".to_string(),
            sex: 0,
            pheno: "0".to_string(),
        }
    }

    fn order(ids: &[&str]) -> Vec<String> {
        let mut v: Vec<String> = ids.iter().map(|s| s.to_string()).collect();
        v.sort_by(|a, b| king_id_cmp(a.as_bytes(), b.as_bytes()));
        v
    }

    /// Every row is a capture from the reference: input order → emitted `.kin` order.
    #[test]
    fn id_comparator_matches_the_probe_families() {
        assert_eq!(order(&["007", "7", "70"]), ["7", "70", "007"]);
        assert_eq!(
            order(&["00", "000", "1", "01", "001"]),
            ["1", "00", "01", "000", "001"]
        );
        assert_eq!(
            order(&["9", "10", "11", "100", "99"]),
            ["9", "10", "11", "99", "100"]
        );
        assert_eq!(
            order(&["1", "a", "1a", "a1", "10", "z"]),
            ["a", "a1", "z", "1", "1a", "10"]
        );
        assert_eq!(order(&["a", "A9", "a10", "A2"]), ["a", "A2", "A9", "a10"]);
        assert_eq!(
            order(&["A", "B", "Z", "aa", "bb"]),
            ["A", "aa", "B", "bb", "Z"]
        );
        assert_eq!(
            order(&["[1", "A1", "`1", "^1", "{1", "~1", "z1"]),
            ["A1", "z1", "[1", "^1", "`1", "{1", "~1"]
        );
        assert_eq!(order(&["ab", "a1", "a2", "ac"]), ["ab", "ac", "a1", "a2"]);
        assert_eq!(
            order(&["aa1", "a1a", "aab", "a1"]),
            ["aab", "aa1", "a1", "a1a"]
        );
        assert_eq!(order(&["x", "x0", "x00", "x1"]), ["x", "x0", "x1", "x00"]);
        assert_eq!(order(&["9z", "z9", "a9", "9a"]), ["a9", "z9", "9a", "9z"]);
        // No integer parse: 26-digit strings still order by length then bytes.
        assert_eq!(
            order(&[
                "99999999999999999999999999",
                "1",
                "99999999999999999999999998"
            ]),
            [
                "1",
                "99999999999999999999999998",
                "99999999999999999999999999"
            ]
        );
    }

    #[test]
    fn within_family_order_ignores_fam_order() {
        // Captured: the same family emits `2 9 10` however the `.fam` lists it.
        let scrambled = [sample("F", "2"), sample("F", "10"), sample("F", "9")];
        assert_eq!(within_family_pairs(&scrambled), [(0, 2), (0, 1), (2, 1)]);
        // ...and FID blocks are sorted by the same comparator.
        let two = [
            sample("10", "a"),
            sample("2", "a"),
            sample("10", "b"),
            sample("2", "b"),
        ];
        assert_eq!(within_family_pairs(&two), [(1, 3), (0, 2)]);
    }

    #[test]
    fn between_family_order_is_block_tiled() {
        let samples: Vec<Sample> = (0..10).map(|i| sample(&format!("F{i}"), "x")).collect();
        let pairs = between_family_pairs(&samples, 8);
        // The first tile is the whole 0..8 triangle; only then does column tile 1 open.
        assert_eq!(&pairs[..3], [(0, 1), (0, 2), (0, 3)]);
        let first_cross = pairs.iter().position(|&(_, j)| j >= 8).unwrap();
        assert_eq!(pairs[first_cross], (0, 8));
        assert_eq!(first_cross, 28); // C(8,2) pairs precede it
                                     // Plain row-major would have put (0,8) at index 7.
        assert_ne!(pairs, serial_pairs(10));
    }

    #[test]
    fn fixed_formatting_uses_c_spellings() {
        assert_eq!(f(0.5, 4), "0.5000");
        assert_eq!(f(-0.000_01, 4), "-0.0000");
        assert_eq!(f(f64::NAN, 4), "nan");
        assert_eq!(f(f64::NEG_INFINITY, 4), "-inf");
        assert_eq!(f(41_542_807.0, 3), "41542807.000");
        assert_eq!(g(0.0), "0");
        assert_eq!(g(0.5), "0.5");
        assert_eq!(g(1.0), "1");
    }

    #[test]
    fn prefix_is_concatenated_not_joined() {
        let mut opts = Options::new();
        assert_eq!(out_path(&opts, ".ibs"), "king.ibs");
        assert_eq!(out_path(&opts, "allsegs.txt"), "kingallsegs.txt");
        opts = crate::cli::parse(&["--prefix".to_string(), "ZZ_".to_string()]).options;
        assert_eq!(out_path(&opts, ".ibs"), "ZZ_.ibs");
        assert_eq!(out_path(&opts, "allsegs.txt"), "ZZ_allsegs.txt");
    }

    #[test]
    fn options_in_effect_is_a_whitelist_in_table_order() {
        let parsed = |args: &[&str]| {
            let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
            crate::cli::parse(&owned).options
        };
        let o = parsed(&["--duplicate"]);
        assert_eq!(options_in_effect(&o, Opt::Duplicate), ["--duplicate"]);
        let o = parsed(&["--duplicate", "--minConc", "0.9", "--cpus", "1"]);
        assert_eq!(
            options_in_effect(&o, Opt::Duplicate),
            ["--duplicate", "--minConc 0.9", "--cpus 1"]
        );
        // `%G`, so an integral value loses its decimal point.
        let o = parsed(&["--duplicate", "--minConc", "1"]);
        assert_eq!(options_in_effect(&o, Opt::Duplicate)[1], "--minConc 1");
        let o = parsed(&["--duplicate", "--minConc", "0"]);
        assert_eq!(options_in_effect(&o, Opt::Duplicate)[1], "--minConc 0");
        // --seglength is accepted and changes behaviour but is never echoed.
        let o = parsed(&["--ibs", "--seglength", "10"]);
        assert_eq!(options_in_effect(&o, Opt::Ibs), ["--ibs"]);
    }
}
