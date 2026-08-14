//! `<prefix>X.seg` — the X-chromosome half of `--ibdseg`.
//!
//! One extra file and one extra console line, both appended to the end of the `--ibdseg`
//! body:
//!
//! ```text
//! Summary statistics of IBD segments for individual pairs saved in file <p>.seg
//! Additional summary statistics of X-Chr IBD segments saved in file <p>X.seg
//! ```
//!
//! # When it is written
//!
//! Two conditions, and neither is the one the file's name suggests.
//!
//! * **`--degree` must be given**, with a non-zero value. This is the reference's own
//!   oddity: plain `--ibdseg` on a fileset with 1 500 X markers writes no `X.seg` at all,
//!   `--degree 1`, `2`, `3`, `5` and `10` all write it, and `--degree 0` — which is how an
//!   absent `--degree` is spelled internally — does not. Even `--degree -1`, which filters
//!   every row out of `<prefix>.seg`, still writes a header-only `X.seg`. The *value*
//!   changes only which rows survive; its presence is what creates the file.
//! * **the X map must yield a usable segment**, by exactly the construction
//!   [`king_core::ibdseg::usable_segments`] applies to the autosomes: at least five
//!   complete 64-marker words, spanning more than 10 Mb after word alignment, with no
//!   marker gap over 1 Mb. This is the same test that prints `In addition to autosomes,
//!   <n> segments of length <x> Mb on X-chr can be further used.`, and the two always
//!   agree — that line and this file appear together or not at all.
//!
//! **There is no 512-marker threshold here.** `--kinship`'s X pass has one
//! ([`crate::load::X_PASS_MIN_SNPS`]) and this pass does not: 320 X markers over 30 Mb
//! write `X.seg`, 319 do not, and 640 markers packed into 5 Mb do not either. Swept
//! against the reference at 320/384/448/511/512/513/576/640/1000 markers — the count
//! matters only through the usable-segment construction. There is no family-count
//! condition either: a single-family fileset writes `X.seg` just as `--kinship`'s X pass
//! refuses to run on one.
//!
//! # Rows
//!
//! **Exactly the rows of `<prefix>.seg`, in exactly its order** — the same pairs, the same
//! 16-sample tiling ([`super::ibdseg`]), including the `--degree` filter, which is applied
//! to the *autosomal* `PropIBD` and never recomputed on X. The X table is a second set of
//! columns for the pairs the autosomal pass already chose to report, so:
//!
//! * a pair with no X sharing at all still gets a row, printed `0.0000 0.0000 0.0000` —
//!   `sexchr`'s father–son pairs are exactly that;
//! * the `>10 Mb` long-segment pair filter and the `--seglength` floor are the autosomal
//!   pass's; nothing re-screens a pair on its X evidence;
//! * an empty `<prefix>.seg` gives a header-only `X.seg`, still announced;
//! * samples whose `.fam` sex is neither 1 nor 2 are **not** excluded, unlike in
//!   `--kinship`'s `X.kin`. Their rows appear with the raw code in the `Sex` column.
//!
//! # The estimates
//!
//! The autosomal segment caller, run unchanged over the X marker array and the X planes.
//! There is no male/female branch anywhere in it, and none is needed: a hemizygous male is
//! stored **homozygous**, so he is never heterozygous and never IBS0 against a
//! heterozygote. That single representational fact is what makes the reference's X output
//! read the way it does — a father–son pair scores 0 (no shared X), a father–daughter pair
//! 1.0000 IBD1, and two brothers who drew the same maternal X score IBD2.
//!
//! `--seglength` applies on X exactly as on the autosomes (checked at 3/5/10 Mb: the X
//! columns move), and the denominator is the X usable-segment total, computed over the X
//! array alone — the autosomal denominator never enters.
//!
//! # Columns
//!
//! ```text
//! FID1\tID1\tFID2\tID2\tSex1\tSex2\tMaxIBD1\tMaxIBD2\tIBD1Seg\tIBD2Seg\tPropIBD
//! ```
//!
//! **The header is wrong and must be copied wrong.** It names eleven columns; every data
//! row carries nine values and a trailing tab. `MaxIBD1` and `MaxIBD2` are never written,
//! so the three numbers a row does carry — `IBD1Seg`, `IBD2Seg`, `PropIBD` — sit under the
//! headings `MaxIBD1`, `MaxIBD2`, `IBD1Seg`. Reproducing the misalignment is the point:
//! it is what the reference emits.
//!
//! `PropIBD` here is the **full-precision** `IBD2Seg + IBD1Seg/2`, not the
//! printed-columns recombination [`king_core::ibdseg::seg_prop_ibd`] that
//! `<prefix>.seg` uses one file away. The `sexchr` capture decides it outright: the pair
//! at `IBD1Seg 0.4257 / IBD2Seg 0.0000` prints `0.2128`, and the printed-column rule gives
//! `0.2129` (the double nearest `4257 * 5e-5` lands just above the decimal half); the pair
//! at `0.9067` prints `0.4533` against that rule's `0.4534`. Two refutations in fourteen
//! rows. `X.seg` and `X.kin` agree with each other value for value, which is the same
//! thing said the other way round — `X.kin` was already known to take full precision.

use std::fmt::Write as _;

use king_core::ibdseg::{self, Usable};
use king_io::{Genotypes, Sample};

use crate::analysis::f;

/// The eleven-name header, two of whose columns are never filled in. See the module docs.
const HEADER: &str =
    "FID1\tID1\tFID2\tID2\tSex1\tSex2\tMaxIBD1\tMaxIBD2\tIBD1Seg\tIBD2Seg\tPropIBD\n";

/// `Additional summary statistics of X-Chr IBD segments saved in file <path>` — the one
/// line this pass adds to the console, after `--ibdseg`'s own `.seg` line.
pub fn saved_line(path: &str) -> String {
    format!("Additional summary statistics of X-Chr IBD segments saved in file {path}\n")
}

/// Whether `--ibdseg` writes `<prefix>X.seg` on this run.
///
/// `degree` is the effective `--degree` (0 when unset), `x_segs` the usable segments of
/// the X marker array. The planes themselves are the caller's to supply; they are decoded
/// whenever the map holds an X marker, so a non-empty `x_segs` always has them.
pub fn runs(degree: i32, x_segs: &[Usable]) -> bool {
    degree != 0 && !x_segs.is_empty()
}

/// The X marker array prepared for scanning, alongside the planes it indexes.
///
/// Built from the pieces `--ibdseg` has already computed for `<prefix>allsegs.txt`, so the
/// usable-segment construction runs once per invocation rather than once per file.
pub struct XSegments<'a> {
    pos: &'a [i64],
    segs: &'a [Usable],
    genotypes: &'a Genotypes,
    denom: i64,
    seglength_bp: i64,
}

impl<'a> XSegments<'a> {
    pub fn new(
        pos: &'a [i64],
        segs: &'a [Usable],
        genotypes: &'a Genotypes,
        seglength_bp: i64,
    ) -> Self {
        XSegments {
            pos,
            segs,
            genotypes,
            denom: ibdseg::denominator(segs, pos),
            seglength_bp,
        }
    }

    /// `(IBD1Seg, IBD2Seg, PropIBD)` on X for one pair, unrounded.
    fn pair(&self, i: usize, j: usize) -> (f64, f64, f64) {
        let s = ibdseg::pair_segments(self.genotypes, self.pos, self.segs, i, j, self.seglength_bp);
        (
            s.ibd1_seg(self.denom),
            s.ibd2_seg(self.denom),
            s.prop_ibd(self.denom),
        )
    }
}

/// Render `<prefix>X.seg` for `pairs`, which must be `<prefix>.seg`'s reported pairs in
/// its own order.
pub fn text(samples: &[Sample], pairs: &[(usize, usize)], x: &XSegments) -> String {
    let mut s = String::from(HEADER);
    for &(i, j) in pairs {
        let (pi1, pi2, prop) = x.pair(i, j);
        s.push_str(&row(samples, i, j, pi1, pi2, prop));
    }
    s
}

/// One data row: nine values and a trailing tab, against an eleven-name header — the
/// reference's own misalignment, reproduced deliberately.
///
/// The `Sex` columns are the `.fam` codes printed raw, with no exclusion and no
/// normalisation; `%.4lf` on the three estimates.
fn row(samples: &[Sample], i: usize, j: usize, pi1: f64, pi2: f64, prop: f64) -> String {
    let mut s = String::new();
    let _ = writeln!(
        s,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
        samples[i].fid,
        samples[i].iid,
        samples[j].fid,
        samples[j].iid,
        samples[i].sex,
        samples[j].sex,
        f(pi1, 4),
        f(pi2, 4),
        f(prop, 4),
    );
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seg(lo: usize, hi: usize) -> Usable {
        Usable { chr: 23, lo, hi }
    }

    /// Only the *presence* of `--degree` creates the file, and 0 is how an absent
    /// `--degree` is spelled. `-1` is present, and writes a header-only file.
    #[test]
    fn degree_zero_is_an_absent_degree() {
        let segs = [seg(0, 639)];
        assert!(!runs(0, &segs));
        for d in [-1, 1, 2, 3, 5, 10] {
            assert!(runs(d, &segs), "--degree {d}");
        }
    }

    /// No usable X segment, no file — whatever the marker count and whatever the degree.
    #[test]
    fn an_empty_x_map_never_runs() {
        for d in [-1, 0, 1, 2, 10] {
            assert!(!runs(d, &[]), "--degree {d}");
        }
    }

    #[test]
    fn the_header_names_two_columns_no_row_fills_in() {
        let names: Vec<&str> = HEADER.trim_end().split('\t').collect();
        assert_eq!(names.len(), 11);
        assert_eq!(&names[6..8], ["MaxIBD1", "MaxIBD2"]);
    }

    fn sample(fid: &str, iid: &str, sex: u8) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: "0".to_string(),
            mat: "0".to_string(),
            sex,
            pheno: "-9".to_string(),
        }
    }

    /// Two rows copied verbatim out of the `sexchr --ibdseg --degree 2` capture — the
    /// misaligned trailing tab, the raw sex codes, and the full-precision `PropIBD` that
    /// `.seg`'s printed-column rule would render `0.2129` instead.
    #[test]
    fn a_row_matches_the_reference_byte_for_byte() {
        let s = vec![
            sample("SEX", "S_SON1", 1),
            sample("SEX", "S_DAU1", 2),
            sample("SEX", "S_F", 1),
        ];
        assert_eq!(
            // `IBD1Seg` prints 0.4257 and half of it prints 0.2128 — the printed-column
            // rule would give 0.2129 off those same digits, which is what rules it out.
            row(&s, 0, 1, 0.4256996, 0.0, 0.4256996 / 2.0),
            "SEX\tS_SON1\tSEX\tS_DAU1\t1\t2\t0.4257\t0.0000\t0.2128\t\n"
        );
        assert_eq!(
            row(&s, 2, 0, 0.0, 0.0, 0.0),
            "SEX\tS_F\tSEX\tS_SON1\t1\t1\t0.0000\t0.0000\t0.0000\t\n"
        );
    }

    /// A `.fam` sex outside 1/2 is printed, not filtered — `--kinship`'s `X.kin` drops
    /// such a sample, this file does not.
    #[test]
    fn an_unknown_sex_still_gets_a_row() {
        let s = vec![sample("F1", "A", 0), sample("F1", "B", 2)];
        assert_eq!(
            row(&s, 0, 1, 1.0, 0.0, 0.5),
            "F1\tA\tF1\tB\t0\t2\t1.0000\t0.0000\t0.5000\t\n"
        );
    }

    /// The console line the pass appends, prefix included.
    #[test]
    fn the_saved_line_carries_the_prefixed_name() {
        assert_eq!(
            saved_line("ZZ_X.seg"),
            "Additional summary statistics of X-Chr IBD segments saved in file ZZ_X.seg\n"
        );
    }
}
