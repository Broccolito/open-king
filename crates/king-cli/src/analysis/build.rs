//! `--build`: family clustering followed by pedigree reconstruction.
//!
//! The pass opens with the same clustering prologue as `--cluster` and `--unrelated`
//! ([`unrelated::clustering_prologue`]) and then reconstructs pedigrees inside the
//! families that clustering built:
//!
//! ```text
//! Pedigree reconstruction starts at <t>
//! Reconstructing pedigree...
//! Age information not provided.
//! Total length of <k> chromosomal segments usable for IBD segment analysis is <d> Mb.
//!   Information of these chromosomal segments can be found in file <p>allsegs.txt
//!
//! <the reconstruction log, echoed to the console as well as to the file>
//!
//! Details of pedigree reconstruction are available in log file <p>build.log
//! Update-ID information is saved in file <p>updateids.txt
//! No pedigrees can be reconstructed.
//! Pedigree reconstruction ends at <t>
//! ```
//!
//! Note the segment total is printed **twice** — once by the clustering prologue and
//! again here — and that the two blank lines above `Details of …` are the segment
//! block's own blank plus the one that closes an empty log.
//!
//! # What gets reconstructed
//!
//! Only the **newly clustered** families do. Every capture in which clustering joined
//! nothing — nine of the thirteen, `threegen`'s single twelve-member three-generation
//! family included — produces a **zero-byte** `<p>build.log`, a **zero-byte**
//! `<p>updateparents.txt`, no `<p>updateids.txt` at all despite the console line saying
//! otherwise, and the closing `No pedigrees can be reconstructed.`. `bigish`, the one
//! fileset large enough for the hundred-sample clustering gate to fire, is also the one
//! that logs rules and inferences and writes all three files.
//!
//! **The reconstruction itself is unimplemented**; this module reproduces the
//! no-reconstruction outcome and stops. See `docs/PARITY.md` §6.2.
//!
//! # `INFERENCE AV.FS` and its `Join3/Join2` — measured, not implemented
//!
//! The one statistic the log needs was open when §6.2 was written. It is now identified,
//! and it is **segment-derived**, which is why implementing the surrounding rules would
//! still not make `apps/bigish__build` pass.
//!
//! For an ordered triple `(R; N1, N2)` write `IBD(x, y)` for the union of that pair's
//! called IBD1 and IBD2 segments, as a set of base pairs on the usable-segment map. Then
//!
//! ```text
//! Join2 = | IBD(R,N1) ∩ IBD(R,N2) |
//! Join3 = | IBD(R,N1) ∩ IBD(R,N2) ∩ IBD(N1,N2) |
//! ```
//!
//! and the log prints `Join3/Join2` at `%.3lf`. The genetics behind it: where `R` is IBD
//! to both sibs, a *grandparent* forces the sibs to have inherited the same parental
//! haplotype, so the ratio is 1; an *avuncular* does not, so it sits near 2/3.
//!
//! Scored against **53 `AV.FS` values the reference emitted over 19 filesets** — the
//! corpus `bigish` plus eighteen purpose-built two- and three-family fixtures with
//! sibships of 2…5 — the formula reproduces every one to a mean of **+0.0035**, range
//! **−0.0001 … +0.0118**, one-sided high exactly like the IBD1 residual. Only 5 of the
//! 53 round to the printed three decimals, and **none of `bigish`'s five do**.
//!
//! Two further rules, measured the same way:
//!
//! * **Which two sibs are named is a property of the sibship, not of `R`.** Every
//!   `AV.FS` line raised against one sibship names the same pair whatever `R` is (three
//!   distinct `R` in one fixture, two in another), and that pair is the sibship's first
//!   two members in the order the `RULE FS1` line prints them — a fixture whose log
//!   reads `joins in sibship (A_C2 A_C3 A_C1)` names `A_C2 and A_C3` in all three of its
//!   `AV.FS` lines. What orders a sibship is *not* identified; it is data-dependent, not
//!   a fixed permutation (two-member sibships in one fixture print both `(C1 C2)` and
//!   `(C2 C1)`).
//! * **The verdict is a cut on the ratio.** Below it the line reads `<R> is uncle|aunt of
//!   N1 and N2`; above it, `<R> is grandfather|grandmother, HS, or nephew|niece of N1 and
//!   N2`, the word pairs following `R`'s sex. Over the same 53 values the cut is bracketed
//!   to **(0.848, 0.901)** — largest `uncle` 0.848, smallest ambiguous 0.901 — which does
//!   not separate 0.85, 0.875 and 0.9.
//!
//! `<p>updateparents.txt` needs none of this: it is `RULE FS0`/`FS1` only — every
//! clustered member on one tab-separated `FID IID FATHER MOTHER` row, in `updateids.txt`
//! order, keeping its `.fam` parents, except that each mutually-full-sib group declaring no
//! parents is given the next unused pair of synthetic parent ids (`1 2`, then `3 4`, …,
//! counted across the whole run, one pair per group however many members it has —
//! checked on a three-father sibship). It is left unwritten all the same: the only shape
//! the corpus exercises is that one, the neighbouring `RULE PO.*` family and the phantom
//! materialisation a non-sibship merge triggers are unrecovered, and a writer that
//! handled only the `bigish` shape would be a rule fitted to a single file.

use std::io::Write;
use std::path::Path;

use king_core::infer::{pedigree_kinship, KinshipCache, Pedigree};

use crate::analysis::{band, ibdseg, out_path, unrelated, with_phantom_parents};
use crate::cli::Options;
use crate::console;
use crate::load::Loaded;

/// `PropIBD` below which a pedigree 1st-degree pair is reported as not looking like one:
/// the 1st-degree band edge on `PropIBD`, `2^-1.5`.
///
/// `PropIBD` is twice a kinship, so this is [`band::FIRST`] doubled — the same cut-point
/// the segment `InfType` uses to open its `PO`/`FS` band.
const FIRST_DEGREE_PROP_IBD: f64 = 2.0 * band::FIRST;

/// Run the pass. The caller has already printed `Options in effect:`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let clustering = unrelated::clustering_prologue(opts, loaded, out);
    if clustering.tiny {
        // Under ten samples the reference stops at the disabled notice: no
        // reconstruction block, and not one file.
        return;
    }

    if clustering.any_merged {
        // `--cluster` writes and announces the same file at the same point; `--unrelated`
        // shares the prologue above and neither writes nor announces it.
        let ids_path = out_path(opts, "updateids.txt");
        let _ = std::fs::write(
            &ids_path,
            clustering.updateids_text(&loaded.fileset.samples),
        );
        let _ = out
            .write_all(format!("Update-ID information is saved in file {ids_path}\n\n").as_bytes());
    }

    let _ = out.write_all(
        format!(
            "Pedigree reconstruction starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(b"Reconstructing pedigree...\n");
    // No `.fam` carries ages, and the corpus never provides one, so this notice is
    // unconditional here.
    let _ = out.write_all(b"Age information not provided.\n");
    let _ = out.write_all(ibdseg::segment_prepass(opts, loaded).as_bytes());

    // Console only: `monomorphic`'s warned run leaves a zero-byte `build.log`.
    let _ = out.write_all(first_degree_warnings(opts, loaded).as_bytes());

    // The log, echoed to stdout and written to the file. Empty unless clustering merged
    // families, which is the branch this module does not implement.
    let log = String::new();
    let _ = out.write_all(log.as_bytes());
    let _ = out.write_all(b"\n");

    let log_path = out_path(opts, "build.log");
    let parents_path = out_path(opts, "updateparents.txt");
    let _ = std::fs::write(Path::new(&log_path), log.as_bytes());
    let _ = std::fs::write(Path::new(&parents_path), b"");

    let _ = out.write_all(
        format!("Details of pedigree reconstruction are available in log file {log_path}\n")
            .as_bytes(),
    );
    // Announced whether or not the file is written: on an unmerged run the reference
    // prints this line and leaves no `updateids.txt` behind.
    let _ = out.write_all(
        format!(
            "Update-ID information is saved in file {}\n",
            out_path(opts, "updateids.txt")
        )
        .as_bytes(),
    );
    if clustering.any_merged {
        let _ = out.write_all(
            format!("Update-parent information is saved in file {parents_path}\n").as_bytes(),
        );
    } else {
        let _ = out.write_all(b"No pedigrees can be reconstructed.\n");
    }
    let _ = out.write_all(
        format!(
            "Pedigree reconstruction ends at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
}

/// The `does not look like 1st-degree relatives` block, one per offending pair.
///
/// Reconstruction refuses to trust a family whose declared 1st-degree pairs the genotypes
/// contradict, and says so before it starts:
///
/// ```text
/// Warning: (P_C3 P_C4) does not look like 1st-degree relatives.
/// please fix within-family errors first before pedigree recontruction.
/// ```
///
/// (the misspelling is the reference's, and both lines repeat per pair).
///
/// # The predicate
///
/// A pair the **pedigree** puts in the 1st-degree band whose **segment** `PropIBD` falls
/// below [`FIRST_DEGREE_PROP_IBD`]. It is neither the kinship estimate nor `Error`:
///
/// * `monomorphic`'s `P_C3`/`P_C4` warns at `Kinship 0.1477`, and bisecting it up to
///   0.2169 never silences the warning, while `multifam`'s declared sib pair `B_C1`/`B_C2`
///   at 0.1708 never raises it — so it is not the estimate;
/// * both of those pairs print `Error 0.5` in `.kin` — so it is not `Error` either;
/// * `monomorphic`'s `P_C1`/`P_C4` is `InfType 2nd` and does **not** warn — so it is not
///   `InfType`.
///
/// What separates them is `PropIBD`: 0.2406 for the warned pair against 0.3564 and 0.4487
/// for the two silent ones.
///
/// The cut-point was then bisected against the reference on a fixture built for it:
/// `missing` (6 samples, below the reconstruction gate) padded to 14 with unrelated
/// singletons so that reconstruction runs, and its `M_C2`/`M_C3` pair walked down through
/// the boundary by forcing an opposite homozygote at the first *k* markers. The warning
/// appears between `k = 1124` and `k = 1125`, where the reported `PropIBD` is 0.3541 on
/// both sides, and `InfType` reads `FS` throughout — so the test is on `PropIBD` and not
/// on the label. Two decimals of context: 0.3543 is silent, 0.3538 warns.
///
/// That measured cut sits ~0.0006 **above** the `2^-1.5` used here, which is the one
/// `PropIBD` band edge the binary is otherwise built on (`inf_type`'s `D1`, verified to
/// four decimals). The likeliest explanation is that the internal test and the reported
/// column divide by denominators differing by about 0.15 %, which would put the internal
/// comparison exactly on `2^-1.5`; the alternative is an unexplained 0.3541. The constant
/// is left principled, and the disputed window (0.35355, 0.35415) is 0.0006 wide and
/// empty in the whole corpus — the nearest pair either way is 0.3538 and 0.3564.
///
/// Across all thirteen corpus datasets the rule fires on exactly the one pair the
/// reference warns about. It reaches the right answer only when the segment engine does:
/// `multifam`'s `B_C2`/`B_C3` is 0.3583 for the reference and 0.3526 here, so this warning
/// currently fires there too. That is `docs/PARITY.md` §4.1, not a second rule.
fn first_degree_warnings(opts: &Options, loaded: &Loaded) -> String {
    let Some(segments) = ibdseg::Segments::new(opts, loaded) else {
        return String::new();
    };
    let samples = &loaded.fileset.samples;
    let pedigree = Pedigree::from_samples(&with_phantom_parents(samples));
    let mut cache = KinshipCache::default();
    let mut s = String::new();
    for block in crate::analysis::family_blocks(samples) {
        for (n, &a) in block.iter().enumerate() {
            for &b in &block[n + 1..] {
                let phi = pedigree_kinship(&pedigree, &mut cache, a, b);
                if !(band::FIRST..band::MZ).contains(&phi) {
                    continue;
                }
                if segments.of(loaded, a, b).2 < FIRST_DEGREE_PROP_IBD {
                    s.push_str(&format!(
                        "Warning: ({} {}) does not look like 1st-degree relatives.\n\
                         please fix within-family errors first before pedigree recontruction.\n",
                        samples[a].iid, samples[b].iid
                    ));
                }
            }
        }
    }
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The cut-point is the 1st-degree band edge on `PropIBD` — a doubled kinship — and it
    /// separates every pair the reference was observed to judge either way.
    #[test]
    fn the_first_degree_cut_is_the_doubled_band_edge() {
        assert_eq!(FIRST_DEGREE_PROP_IBD, band::MZ);
        // Warned by the reference: `monomorphic` P_C3/P_C4, and `missing`'s three sib
        // pairs once the fileset is padded past the reconstruction gate.
        for warned in [0.2406, 0.3450, 0.2624, 0.2284] {
            assert!(warned < FIRST_DEGREE_PROP_IBD, "{warned}");
        }
        // Silent: `multifam`'s two sib pairs and `monomorphic`'s P_C1/P_C4, which is
        // `InfType 2nd` and still not warned about.
        for silent in [0.3564, 0.3583, 0.4487] {
            assert!(silent >= FIRST_DEGREE_PROP_IBD, "{silent}");
        }
    }
}
