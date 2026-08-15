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
//! Of the three, **`<p>updateparents.txt` is implemented and byte-identical**, and
//! `<p>build.log`'s **rule half is too** — see the two sections on them below. What the
//! log still omits is its `INFERENCE` half, which needs both an ordering rule this module
//! has *falsified every candidate for* and a segment statistic this engine reproduces only
//! to about 0.005. That is what remains of `docs/PARITY.md` §6.2.
//!
//! # `<p>build.log` — the format, and which half is written
//!
//! The file is a per-cluster narration of the reconstruction. It is written to disk and
//! echoed to **stdout byte for byte**: the console block between the segment prepass and
//! `Details of pedigree reconstruction …` is the file's contents verbatim, plus one blank
//! line after it that the file does not carry (`bigish`'s capture ends `0.803\n`, 806
//! bytes, and its stdout repeats those 18 lines exactly). A cluster contributes nothing at
//! all unless it raises a line — a duplicate-joined or `PO`-joined merge leaves a
//! zero-byte log — and when it does, its first line is the header `Family KING<k>:`.
//!
//! Inside a cluster the order is fixed: the `RULE` lines that build sibships, then one or
//! more blank lines, then the `INFERENCE` lines, grouped by the sibship they are about.
//! The templates, all of them mined from the reference binary's `__cstring` section and
//! then confirmed against runs, are
//!
//! ```text
//! Family KING1:
//!   Family KING1 RULE FS0: Sibship (B01_F B02_F)'s parents are (1 2)
//!   Family KING1 RULE FS1: C_F joins in sibship (A_F B_F)
//!   Family KING1 RULE FS2: Sibship (…) and sibship (…) are combined
//!   Reconstruct parent-offspring pair (C_F, D_F)...
//!
//!   Family KING1 INFERENCE AV.FS: B02_F is uncle of B01_C2 and B01_C3, Join3/Join2=0.778
//!   Family KING1 INFERENCE AV.FS: A_F is grandfather, HS, or nephew of B_C2 and B_C3, …
//!     HS B02_C4 unrelated to B01_M
//!   Family KING1 INFERENCE HS.UN2: B01_C3 and B02_C4 are HS
//! ```
//!
//! `FS0` **creates** a sibship out of a full-sib pair and names the parents it takes;
//! `FS1` adds one more member to a sibship that already exists; `FS2` merges two declared
//! sibships and has never been observed to fire. [`sibship_parents`] documents which of
//! the three a component raises and why. All three word pairs in the `INFERENCE` lines —
//! `uncle`/`aunt`, `grandfather`/`grandmother`, `nephew`/`niece` — follow `R`'s `.fam` sex
//! code, and the choice between the two templates is a cut on the ratio: over 53 measured
//! values the largest `uncle` is 0.846 and the smallest ambiguous one 0.902, which still
//! does not separate 0.85, 0.875 and 0.9.
//!
//! **This module writes the header and the `RULE FS0`/`FS1` lines and stops there.**
//! Scored against the reference on the 30 held-out shapes of
//! `docs/research/fixtures/buildlog.py` (`rules`, over `build_shapes.py`'s twenty merge
//! shapes and ten `avfs.py` ones), the rule lines are
//! byte-identical on **23**; of the other seven, three differ only in which cluster is
//! called `KING1` (the numbering bug at the bottom of this doc), two are the `<FID>-><IID>`
//! renaming shapes that are out of scope for this binary, and two differ only in the
//! *order* the `FS1` line lists an existing sibship's members — the same unidentified
//! ordering that picks the `AV.FS` pair, below. On `bigish` all six lines this module
//! writes are byte-identical to the capture, 243 of its 806 bytes, and every byte written
//! is a byte the reference wrote.
//!
//! The blank lines are **not** written, because their count is a function of the
//! `INFERENCE` half: one is printed before each sibship's block *until the family prints
//! its first inference*, and if it never prints one, every block still prints its blank.
//! That reproduces 1, 1, 2 for `bigish`'s three clusters and 3, 2, 1 for `three_clusters`'s
//! (whose first cluster infers nothing at all), and 12 of the 13 other measured clusters;
//! the one it misses is `three_fs`, which prints three where the rule says two.
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
//! Re-scored on **34 `AV.FS` values over 16 filesets** — the corpus `bigish` plus fifteen
//! purpose-built two- and three-family fixtures with sibships of 2…6, via
//! `docs/research/fixtures/avfs_score.py` — the formula lands one-sided high on every one
//! of them: mean **+0.0039**, range **+0.0003 … +0.0102**, and only **1 of 34** rounds to
//! the printed three decimals. None of `bigish`'s five does; re-measured, they are
//!
//! | triple (as the log names it) | reference | ours | residual |
//! |---|---:|---:|---:|
//! | `B02_F` uncle of `B01_C2`, `B01_C3` | 0.778 | 0.7828 | +0.0048 |
//! | `B01_F` uncle of `B02_C3`, `B02_C4` | 0.801 | 0.8062 | +0.0052 |
//! | `B14_F` uncle of `B13_C2`, `B13_C1` | 0.779 | 0.7828 | +0.0038 |
//! | `B13_F` uncle of `B14_C1`, `B14_C2` | 0.827 | 0.8278 | +0.0008 |
//! | `B25_F` uncle of `B26_C3`, `B26_C1` | 0.803 | 0.8065 | +0.0035 |
//!
//! — so four of the five need the segment caller to move by 0.0035…0.0052 and the fifth by
//! 0.0008, and writing the log today would turn `apps/bigish__build`'s 18 missing lines
//! into 5 wrong numbers rather than a pass. Note also that `avfs_score.py` still prints
//! the **retracted** `[0, dU_sib/(J2/D)]` bound as `OK` on all five: `bigish` is where that
//! bound was measured in sample, and it is 11 of 34 out of it, so those `OK`s are not
//! evidence.
//!
//! ## The residual is not the formula, and the old bound argument is now dead outright
//!
//! An earlier reading of this module claimed the residual was **entirely** accounted for
//! by our sib-pair union over-call `ΔS`, since that can inflate `Join3` by at most `ΔS`
//! and so the ratio by at most `ΔS / Join2` — quoted as *39 of 39 triples inside*
//! `[0, ΔS / Join2]`. **That does not survive fresh shapes and seeds: it is 11 of 34.**
//! It is now refuted on `bigish` itself. Re-measured against the *current* segment engine,
//! all fifteen pair totals behind those five triples — every `IBD1Seg + IBD2Seg` for
//! `(R,N1)`, `(R,N2)` and `(N1,N2)` — match the reference's to four decimals, so `ΔS` is
//! **0.0000** and the bound is `0`, while the residual is unchanged at `+0.0008 … +0.0052`.
//! Five of five now score `OUTSIDE`. The bound was always too crude to be evidence — `dU`
//! is a *rounded total*, and `Join3/Join2` is an intersection of three sets, so it reads
//! segment **placement**, which a matching total says nothing about — and the placement is
//! where the residual lives. Its sign is the giveaway: dilating all three sets slightly
//! grows a triple intersection proportionally more than a double one, so ours is high.
//!
//! What does survive is the conclusion, reached independently: **no variant of the
//! formula removes the residual**, so it is an input problem, not an arithmetic one.
//! Swept over the same 34 triples, every candidate is worse or no better —
//!
//! | variant | exact/34 | mean |
//! |---|---|---|
//! | as above (base pairs, refined endpoints, `IBD1 ∪ IBD2`) | 1 | +0.0039 |
//! | marker counts instead of base pairs | 1 | +0.0039 |
//! | minimum piece length on `Join3` (0.1…3 Mb) | ≤4 | +0.0015…+0.0039 |
//! | minimum piece length on both (0.25…5 Mb) | ≤4 | +0.0014…+0.0038 |
//! | eroding every set by 1…63 markers | ≤10 | +0.0034…−0.0287 |
//!
//! — and the residual is *heteroscedastic*, spanning thirtyfold across triples, which is
//! the signature of data-dependent caller error rather than of a constant the formula is
//! missing. (Eroding by 6 markers zeroes the *mean* at 10 of 34 exact; it is an
//! unprincipled knob, it does not approach exactness, and it is not landed. Two variants
//! measured earlier are also still worse: word-aligning the intervals gives mean −0.025,
//! and re-calling at a different minimum segment length changes nothing below 5 Mb and
//! hurts at 10.)
//!
//! So `apps/bigish__build` is blocked by the segment caller, but note the weaker claim:
//! `docs/PARITY.md` §4.1 closing is *necessary*, and no longer demonstrably sufficient.
//!
//! Three further rules, measured the same way. The second of them is the one that keeps
//! the `INFERENCE` lines unwritten even where the ratio would be close enough to look
//! right, and it is stated here as an open problem with its candidate space *closed*, not
//! as a guess.
//!
//! * **Which `R` an `AV.FS` line can be raised for.** `R` must be an inferred
//!   **2nd-degree** relative of *both* named members. That, and nothing weaker, reproduces
//!   the candidate set every time it was checked: the three-father shape names exactly the
//!   two children of the third family against the father sibship `(A_F B_F)`, the
//!   four-father one exactly the six children of the third and fourth, and the shape whose
//!   families have a single child each names exactly the one candidate — each family's own
//!   children are 1st-degree to their own father and so are excluded. One `R` may print
//!   the *same* line two to four times; the repeat count is per `(R, sibship)`, is not the
//!   number of sib pairs, and is not identified.
//! * **Which two sibs are named — unidentified, and here is the closed candidate space.**
//!   The pair is a property of the sibship: every line raised against one sibship names the
//!   same two whatever `R` is, verified on four sibships against three distinct `R` each
//!   and two more against two each, including cases where the verdicts differ (`uncle` for
//!   one `R`, `grandfather, HS, or nephew` for another). Where the sibship is one a
//!   `RULE FS0`/`FS1` built, the named pair is its **first two members in the order the
//!   rule line prints** — `RULE FS1: B_X joins in sibship (A_C2 A_C3 A_C1)` and the
//!   `AV.FS` line naming `A_C2 and A_C3` come from the same fileset. So the open question
//!   is one ordering, shared by both lines. What it is **not**:
//!   - **not genotype-derived.** Four fresh seeds — complete genotype reseeding, the
//!     sibship's own kinships moving over a 0.10 range — give byte-identical `FS1` orders
//!     at each of three sibship sizes: `(C2 C3 C1)`, `(C3 C4 C2 C1)`, `(C4 C1 C5 C3 C2)`.
//!   - **not the `.fam` row order.** Permuting a sibship's three rows inside the `.fam`
//!     (genotypes moved with them) leaves the named pair on the same two *individuals*,
//!     now at different positions.
//!   - **not the absolute sample index.** Moving all 80 padding singletons of a
//!     four-family fixture to the front of the `.fam` leaves the whole log byte-identical.
//!   - **not the sibship's size or position.** Four three-child sibships in one cluster
//!     print four different orders, and `bigish`'s `B01` and `B13` — structurally
//!     identical three-child sibships, same cluster shape, same sexes, same phenotypes —
//!     name `(C2, C3)` and `(C2, C1)`.
//!   - **not any pairwise statistic.** Over the 27 measured sibships of three or more
//!     children, no `argmin` or `argmax` of `HetHet`, `IBS0`, `HetConc`, `HomIBS0`,
//!     `Kinship`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `N_SNP`, `Z0` or `Error` picks the
//!     named pair more than **11 of 27** times, against a chance baseline of one in three
//!     or worse; and neither does any of ten segment-level statistics computed here
//!     (`|IBD1|`, `|IBD2|`, their union and complement, segment counts, longest segment):
//!     best **6 of 20** on the subset those were run over. The named pair's rank on the
//!     ratio itself runs from first to last — one sibship names the pair with the *lowest*
//!     `Join3/Join2` of its three and another the *highest* of its fifteen.
//!
//!   The earlier reading in this doc — "a function of the pedigree shape alone",
//!   with a positional map — is **withdrawn**: it was measured only on the first family of
//!   two-family fixtures, where the answer happens to be constant, and the four-sibship
//!   fixture refutes it directly.
//! * **The verdict is a cut on the ratio.** Below it the line reads `<R> is uncle|aunt of
//!   N1 and N2`; above it, `<R> is grandfather|grandmother, HS, or nephew|niece of N1 and
//!   N2`, the word pairs following `R`'s sex. Over 53 values the cut is bracketed to
//!   **(0.846, 0.902)** — largest `uncle` 0.846, smallest ambiguous 0.902 — which still
//!   does not separate 0.85, 0.875 and 0.9.
//!
//! # `<p>updateparents.txt` — implemented, and what pinned each clause
//!
//! This file needs none of the segment statistic: it carries only what the `RULE FS*`
//! lines decided. [`sibship_parents`] implements it and [`reconstruct`] writes it, in the
//! same walk that writes those rule lines, so the two cannot disagree.
//! Every clause below was measured against the reference on **twenty held-out merge
//! shapes** built for it, not on `bigish`; the rig is
//! `docs/research/fixtures/build_shapes.py`, which re-runs the scorecard end to end.
//! Eighteen of the twenty are in scope, **fifteen of them byte-identical on the file and
//! on the console tail**, and the other three differ only in which cluster is called
//! `KING1` — the separate numbering bug at the bottom of this doc.
//!
//! * **Rows.** Every member of every merged cluster, `FID IID FATHER MOTHER`, tab
//!   separated, in `updateids.txt` order (cluster, then the ID comparator).
//! * **Parents.** Each member keeps its `.fam` parents, except that a *sibship the
//!   inference touched* takes one parent pair for all of its members. A sibship is a
//!   connected component of `inferred FS` ∪ `declares the same non-missing couple`, and
//!   it takes:
//!   - the couple one of its members already declares, scanning in member order — the
//!     `RULE FS0: Sibship (A_F B_F)'s parents are (A_G_F A_G_M)` form; failing that,
//!   - the next unused pair of synthetic ids, `1 2`, then `3 4`, … — one pair per
//!     sibship however many members it has, counted across the whole run and in cluster
//!     order. Two sibships inside one cluster take consecutive pairs, ordered by their
//!     first member under the ID comparator and **not** by `.fam` row order (the shape
//!     that separates the two writes each family's mother first and names her so she
//!     sorts last; the fathers still take `1 2`).
//! * **A cluster holding an inferred duplicate contributes no rows at all** — not even
//!   identity ones. A cluster merged by `PO` alone does contribute them, with nobody's
//!   parents changed, because without ages the reference will not orient the pair.
//! * **Nothing is written unless some sibship got parents.** A run whose only merges are
//!   `PO` or duplicate leaves a zero-byte file, no `Update-parent information is saved…`
//!   line, and the closing `No pedigrees can be reconstructed.` — even though those
//!   clusters *are* in `updateids.txt`. The old code keyed that tail on `any_merged`,
//!   which is wrong on any fileset whose merges are all `PO`.
//!
//! Two shapes are deliberately **out of scope** and excluded from the scorecard: when a
//! `.fam` names a parent that lives in another family, the reference materialises a
//! phantom for it and, because the id is then no longer unique, renames every individual
//! to `<FID>-><IID>`. Nothing in this binary implements that renaming, so those two
//! captures would fail on `updateids.txt` first.
//!
//! # A held-out clustering bug this rig also found — *not* fixed here
//!
//! Merged clusters are **not** numbered in family order. They are numbered by the
//! *relationship type* of the pair that joins them: every `Dup/MZ`-joined cluster first,
//! then the `PO`-joined ones, then the `FS`-joined ones, ties broken by family order. A
//! fixture whose three clusters are joined, in family order, by `FS`, `PO` and a
//! duplicate is numbered `KING3`, `KING2`, `KING1`. `unrelated::clusters` uses family
//! order alone, which agrees on all thirteen corpus filesets — every merge in `bigish`
//! is `FS` — so no capture moves either way. The fix belongs in `unrelated.rs`.

use std::io::Write;
use std::path::Path;

use king_core::infer::{pedigree_kinship, KinshipCache, Pedigree};
use king_core::{counts, kinship, Scope};
use king_io::Sample;

use crate::analysis::{band, ibdseg, out_path, related, unrelated, with_phantom_parents};
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

    // The log, echoed to stdout and written to the file byte for byte. Only its `RULE`
    // half is built: the `INFERENCE AV.FS` lines need two rules this engine has not
    // identified as well as a statistic it reproduces only to ~0.005, so they are left
    // out rather than invented. See the module doc.
    let Reconstruction { parents, log } = reconstruct(opts, loaded, &clustering);
    let _ = out.write_all(log.as_bytes());
    let _ = out.write_all(b"\n");

    let log_path = out_path(opts, "build.log");
    let parents_path = out_path(opts, "updateparents.txt");
    let _ = std::fs::write(Path::new(&log_path), log.as_bytes());
    let _ = std::fs::write(
        Path::new(&parents_path),
        parents.as_deref().unwrap_or("").as_bytes(),
    );

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
    // Keyed on whether reconstruction actually assigned a parent, not on whether
    // clustering merged anything: a run whose only merges are `PO` or duplicate lands in
    // the `No pedigrees` branch with its clusters still listed in `updateids.txt`.
    if parents.is_some() {
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

/// What one `--build` pass reconstructed: the update file, and the log that narrates it.
struct Reconstruction {
    /// `<prefix>updateparents.txt`, or `None` when reconstruction assigned nobody a parent.
    ///
    /// `None` is the whole-run verdict the console tail prints as `No pedigrees can be
    /// reconstructed.`; the caller still writes the (empty) file, exactly as the reference
    /// leaves a zero-byte one behind.
    parents: Option<String>,
    /// `<prefix>build.log`, which is also echoed to stdout unchanged.
    log: String,
}

/// Reconstruct every merged cluster, building both outputs in one walk.
///
/// The two are two views of the same decisions — `updateparents.txt` records the parent a
/// sibship ended up with, the log records the rule that gave it one — so they are built
/// together and cannot drift apart. The clause-by-clause evidence for both is in the
/// module doc.
fn reconstruct(
    opts: &Options,
    loaded: &Loaded,
    clustering: &unrelated::Clustering,
) -> Reconstruction {
    let samples = &loaded.fileset.samples;
    let types = InfTypes::new(opts, loaded);

    let mut text = String::new();
    let mut log = String::new();
    let mut next_synthetic = 1u32;
    let mut reconstructed = false;

    for (key, members) in clustering.merged() {
        // An unresolved duplicate inside the cluster suppresses it whole: the reference
        // emits neither a rule nor a single identity row for such a family.
        if pairs(members).any(|(a, b)| types.of(loaded, a, b) == Some("Dup/MZ")) {
            continue;
        }
        let is_fs = |i: usize, j: usize| types.of(loaded, i, j) == Some("FS");
        let assigned = sibship_parents(samples, members, &is_fs, &mut next_synthetic);
        reconstructed |= !assigned.parents.is_empty();
        // The header is printed once, immediately before the cluster's first line — a
        // cluster that raises no rule at all does not appear in the log.
        if !assigned.rules.is_empty() {
            log.push_str(&console::build_family_header(key));
            let name = |n: &usize| samples[members[*n]].iid.as_str();
            for rule in &assigned.rules {
                let names: Vec<&str> = rule.members.iter().map(name).collect();
                log.push_str(&match rule.joiner {
                    None => console::build_rule_fs0(key, &names, &rule.pat, &rule.mat),
                    Some(j) => console::build_rule_fs1(key, name(&j), &names),
                });
            }
        }
        for (n, &i) in members.iter().enumerate() {
            let (pat, mat) = assigned
                .parents
                .iter()
                .find(|(m, _)| *m == n)
                .map(|(_, p)| (p.0.as_str(), p.1.as_str()))
                .unwrap_or((samples[i].pat.as_str(), samples[i].mat.as_str()));
            text.push_str(&format!("{key}\t{}\t{pat}\t{mat}\n", samples[i].iid));
        }
    }
    Reconstruction {
        parents: reconstructed.then_some(text),
        log,
    }
}

/// One `RULE FS0` or `RULE FS1` line, with its people held as positions in the cluster.
///
/// `joiner` is `None` for `FS0` — the rule that *creates* a sibship and names its parents —
/// and `Some(who)` for `FS1`, where `members` is the sibship as it stood before the join
/// and `pat`/`mat` are unused.
struct Rule {
    joiner: Option<usize>,
    members: Vec<usize>,
    pat: String,
    mat: String,
}

/// The parent pairs and the rule lines one cluster's reconstruction produces.
struct Assignment {
    /// `(position in the cluster, (father, mother))`.
    parents: Vec<(usize, (String, String))>,
    rules: Vec<Rule>,
}

/// Every unordered pair of positions in `members`, as index pairs into it.
fn pairs(members: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    (0..members.len())
        .flat_map(move |a| (a + 1..members.len()).map(move |b| (members[a], members[b])))
}

/// The parent pair each sibship in one cluster takes, and the rules that gave it one.
///
/// A sibship is a connected component of `inferred FS` ∪ `declares the same non-missing
/// couple`, and only a component the *inference* touched is given parents — a declared
/// sibship nobody was joined to keeps what the `.fam` already says. `next_synthetic`
/// carries the `1 2`, `3 4`, … counter across clusters, so it advances only where a
/// sibship actually needed inventing.
///
/// # Which rule a component raises
///
/// A component that already contains a **declared sibship** — two or more members naming
/// the same couple — did not have to be created, so every other member of the component
/// `FS1`-joins it, one line each, and no `FS0` is raised: `fs_kids`'s `B_X joins in
/// sibship (A_C2 A_C3 A_C1)` is exactly that. Every other component *is* created, by the
/// full-sib pair that opens it: `FS0` names those two and the couple they take, and the
/// third and later members `FS1`-join in turn — `m43`'s `FS0 (A_F B_F)`, then `C_F joins
/// in sibship (A_F B_F)`, then `D_F joins in sibship (A_F B_F C_F)`. One member declaring
/// a couple is not a declared sibship: `fs_one_declared` raises `FS0` and the pair takes
/// that member's parents rather than a synthetic one.
///
/// A component holding **two** declared sibships is the reference's `RULE FS2: Sibship
/// (…) and sibship (…) are combined`. Nothing here has ever raised it and no capture
/// contains one, so no rule is emitted for such a component — its parents are still
/// assigned, which is what `updateparents.txt` needs.
fn sibship_parents(
    samples: &[Sample],
    members: &[usize],
    is_fs: &dyn Fn(usize, usize) -> bool,
    next_synthetic: &mut u32,
) -> Assignment {
    let n = members.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    let union = |parent: &mut Vec<usize>, a: usize, b: usize| {
        let (ra, rb) = (find(parent, a), find(parent, b));
        if ra != rb {
            parent[rb] = ra;
        }
    };

    let mut joined = vec![false; n];
    for a in 0..n {
        for b in a + 1..n {
            let (i, j) = (members[a], members[b]);
            if is_fs(i, j) {
                union(&mut parent, a, b);
                joined[a] = true;
                joined[b] = true;
            }
            // Declared full sibs: the same couple, both named.
            if samples[i].pat != "0"
                && samples[i].pat == samples[j].pat
                && samples[i].mat == samples[j].mat
            {
                union(&mut parent, a, b);
            }
        }
    }

    // Components in the order of their first member, which is the ID comparator order the
    // cluster is already held in.
    let mut out = Assignment {
        parents: Vec::new(),
        rules: Vec::new(),
    };
    let mut seen: Vec<usize> = Vec::new();
    for a in 0..n {
        let root = find(&mut parent, a);
        if seen.contains(&root) {
            continue;
        }
        seen.push(root);
        let group: Vec<usize> = (0..n).filter(|&b| find(&mut parent, b) == root).collect();
        if !group.iter().any(|&b| joined[b]) {
            continue;
        }
        let couple = group
            .iter()
            .map(|&b| &samples[members[b]])
            .find(|s| s.pat != "0")
            .map(|s| (s.pat.clone(), s.mat.clone()))
            .unwrap_or_else(|| {
                let pair = (
                    next_synthetic.to_string(),
                    (*next_synthetic + 1).to_string(),
                );
                *next_synthetic += 2;
                pair
            });
        out.rules
            .extend(component_rules(samples, members, &group, &couple));
        for b in group {
            out.parents.push((b, couple.clone()));
        }
    }
    out
}

/// The `RULE FS0`/`FS1` lines one component raises, in the order the log prints them.
///
/// Empty for the `FS2` shape — two declared sibships in one component — which is
/// unimplemented and documented as such on [`sibship_parents`].
fn component_rules(
    samples: &[Sample],
    members: &[usize],
    group: &[usize],
    couple: &(String, String),
) -> Vec<Rule> {
    // The declared sibships inside the component: each named couple, with the members
    // that name it, in cluster order.
    let mut declared: Vec<(&str, &str, Vec<usize>)> = Vec::new();
    for &b in group {
        let s = &samples[members[b]];
        if s.pat == "0" {
            continue;
        }
        match declared
            .iter_mut()
            .find(|(pat, mat, _)| *pat == s.pat && *mat == s.mat)
        {
            Some((_, _, who)) => who.push(b),
            None => declared.push((&s.pat, &s.mat, vec![b])),
        }
    }
    declared.retain(|(_, _, who)| who.len() > 1);
    if declared.len() > 1 {
        return Vec::new();
    }

    let mut rules = Vec::new();
    let mut sibship: Vec<usize> = match declared.first() {
        Some((_, _, who)) => who.clone(),
        None => {
            let opening: Vec<usize> = group.iter().take(2).copied().collect();
            rules.push(Rule {
                joiner: None,
                members: opening.clone(),
                pat: couple.0.clone(),
                mat: couple.1.clone(),
            });
            opening
        }
    };
    let joiners: Vec<usize> = group
        .iter()
        .filter(|b| !sibship.contains(b))
        .copied()
        .collect();
    for b in joiners {
        rules.push(Rule {
            joiner: Some(b),
            members: sibship.clone(),
            pat: String::new(),
            mat: String::new(),
        });
        sibship.push(b);
    }
    rules
}

/// The `InfType` of a cross-family pair that clustering treated as 1st-degree.
///
/// Same two steps the screening summary uses — the kinship gate that merged the families,
/// then the segment `InfType` that labels the pair — kept in one place so reconstruction
/// and the summary cannot disagree about what joined a cluster. Within-family pairs are
/// `None`: their relationship is declared, not inferred.
struct InfTypes {
    engine: related::Engine,
    edge: f64,
}

impl InfTypes {
    fn new(opts: &Options, loaded: &Loaded) -> Self {
        InfTypes {
            engine: related::Engine::autosomes(
                &loaded.fileset.variants,
                &loaded.fileset.kept,
                i64::from(opts.int(crate::cli::Opt::Sexchr)),
                ibdseg::seglength_bp(opts),
            ),
            edge: band::FIRST,
        }
    }

    fn of(&self, loaded: &Loaded, i: usize, j: usize) -> Option<&'static str> {
        let samples = &loaded.fileset.samples;
        if samples[i].fid == samples[j].fid {
            return None;
        }
        let genotypes = &loaded.fileset.genotypes;
        let counts = counts::pair_counts(genotypes, i, j);
        if counts.n_snp == 0 || kinship::kinship(&counts, Scope::BetweenFamily) <= self.edge {
            return None;
        }
        Some(
            self.engine
                .pair(genotypes, i, j)
                .inf_type(kinship::het_concordance(&counts)),
        )
    }
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

    fn person(fid: &str, iid: &str, pat: &str, mat: &str) -> Sample {
        Sample {
            fid: fid.into(),
            iid: iid.into(),
            pat: pat.into(),
            mat: mat.into(),
            sex: 0,
            pheno: "-9".into(),
        }
    }

    /// `bigish`'s shape: two nuclear families whose fathers turn out to be full sibs.
    /// Members are held in the ID comparator's order, as a cluster always is.
    fn two_families() -> Vec<Sample> {
        vec![
            person("FA", "A_C1", "A_F", "A_M"),
            person("FA", "A_C2", "A_F", "A_M"),
            person("FA", "A_F", "0", "0"),
            person("FA", "A_M", "0", "0"),
            person("FB", "B_C1", "B_F", "B_M"),
            person("FB", "B_C2", "B_F", "B_M"),
            person("FB", "B_F", "0", "0"),
            person("FB", "B_M", "0", "0"),
        ]
    }

    fn assign(
        samples: &[Sample],
        fs: &[(usize, usize)],
        next: &mut u32,
    ) -> Vec<(String, String, String)> {
        let members: Vec<usize> = (0..samples.len()).collect();
        let is_fs = |i: usize, j: usize| fs.contains(&(i, j)) || fs.contains(&(j, i));
        let got = sibship_parents(samples, &members, &is_fs, next);
        members
            .iter()
            .enumerate()
            .map(|(n, &i)| {
                let (pat, mat) = got
                    .parents
                    .iter()
                    .find(|(m, _)| *m == n)
                    .map(|(_, p)| (p.0.clone(), p.1.clone()))
                    .unwrap_or((samples[i].pat.clone(), samples[i].mat.clone()));
                (samples[i].iid.clone(), pat, mat)
            })
            .collect()
    }

    /// A sibship the inference joined and nobody declared parents for takes the next
    /// synthetic pair, one pair for the whole sibship; every other member keeps its
    /// `.fam` parents, the declared sibship of children included.
    #[test]
    fn an_undeclared_sibship_takes_the_next_synthetic_pair() {
        let samples = two_families();
        let mut next = 1;
        // A_F (2) and B_F (6) are the inferred full sibs.
        let rows = assign(&samples, &[(2, 6)], &mut next);
        assert_eq!(rows[2], ("A_F".into(), "1".into(), "2".into()));
        assert_eq!(rows[6], ("B_F".into(), "1".into(), "2".into()));
        // The mothers are in no sibship, the children in a declared one nobody joined.
        assert_eq!(rows[3], ("A_M".into(), "0".into(), "0".into()));
        assert_eq!(rows[0], ("A_C1".into(), "A_F".into(), "A_M".into()));
        // One pair consumed, however many members the sibship had.
        assert_eq!(next, 3);
    }

    /// Two sibships in one cluster take consecutive pairs, ordered by their first member
    /// under the ID comparator the cluster is already sorted by.
    #[test]
    fn two_sibships_take_consecutive_pairs_in_member_order() {
        let samples = two_families();
        let mut next = 1;
        // Fathers (2, 6) and mothers (3, 7) are each an inferred sibship.
        let rows = assign(&samples, &[(2, 6), (3, 7)], &mut next);
        assert_eq!(rows[2].1, "1");
        assert_eq!(rows[6].1, "1");
        assert_eq!(rows[3], ("A_M".into(), "3".into(), "4".into()));
        assert_eq!(rows[7], ("B_M".into(), "3".into(), "4".into()));
        assert_eq!(next, 5);
    }

    /// When a member of the sibship already declares a couple, the sibship inherits it
    /// and no synthetic id is spent — the `RULE FS0: … parents are (A_G_F A_G_M)` form.
    #[test]
    fn a_declared_couple_wins_over_a_synthetic_pair() {
        let mut samples = two_families();
        samples[2] = person("FA", "A_F", "A_G_F", "A_G_M");
        let mut next = 1;
        let rows = assign(&samples, &[(2, 6)], &mut next);
        assert_eq!(rows[2], ("A_F".into(), "A_G_F".into(), "A_G_M".into()));
        assert_eq!(rows[6], ("B_F".into(), "A_G_F".into(), "A_G_M".into()));
        assert_eq!(next, 1);
    }

    /// An outsider who is a full sib of a *declared* sibship joins it and takes its
    /// parents — `RULE FS1: B_X joins in sibship (…)` — and the sibship keeps them.
    #[test]
    fn joining_a_declared_sibship_inherits_its_parents() {
        let mut samples = two_families();
        samples.push(person("FB", "B_X", "0", "0"));
        let mut next = 1;
        // B_X (8) is a full sib of A_C1 (0) and A_C2 (1), who declare A_F/A_M.
        let rows = assign(&samples, &[(0, 8), (1, 8)], &mut next);
        assert_eq!(rows[8], ("B_X".into(), "A_F".into(), "A_M".into()));
        assert_eq!(rows[0], ("A_C1".into(), "A_F".into(), "A_M".into()));
        assert_eq!(next, 1);
    }

    /// The rule lines one cluster raises, rendered the way the log renders them.
    fn rule_lines(samples: &[Sample], fs: &[(usize, usize)], next: &mut u32) -> Vec<String> {
        let members: Vec<usize> = (0..samples.len()).collect();
        let is_fs = |i: usize, j: usize| fs.contains(&(i, j)) || fs.contains(&(j, i));
        let got = sibship_parents(samples, &members, &is_fs, next);
        let name = |n: &usize| samples[members[*n]].iid.as_str();
        got.rules
            .iter()
            .map(|r| {
                let names: Vec<&str> = r.members.iter().map(name).collect();
                match r.joiner {
                    None => console::build_rule_fs0("KING1", &names, &r.pat, &r.mat),
                    Some(j) => console::build_rule_fs1("KING1", name(&j), &names),
                }
            })
            .collect()
    }

    /// A sibship the inference creates opens with `FS0` naming the pair and the couple it
    /// takes; the third and later members `FS1`-join it one line at a time, each naming the
    /// sibship as it stood. `m43`'s capture is exactly these three lines.
    #[test]
    fn a_created_sibship_opens_with_fs0_and_grows_by_fs1() {
        let mut samples = two_families();
        samples.push(person("FC", "C_F", "0", "0"));
        let mut next = 1;
        // The three fathers — A_F (2), B_F (6), C_F (8) — are mutual full sibs.
        let got = rule_lines(&samples, &[(2, 6), (2, 8), (6, 8)], &mut next);
        assert_eq!(
            got,
            vec![
                "  Family KING1 RULE FS0: Sibship (A_F B_F)'s parents are (1 2)\n",
                "  Family KING1 RULE FS1: C_F joins in sibship (A_F B_F)\n",
            ]
        );
    }

    /// A component that already holds a declared sibship was never created, so it raises no
    /// `FS0` — the outsider just `FS1`-joins what the `.fam` already declared.
    #[test]
    fn joining_a_declared_sibship_raises_only_fs1() {
        let mut samples = two_families();
        samples.push(person("FB", "B_X", "0", "0"));
        let mut next = 1;
        let got = rule_lines(&samples, &[(0, 8), (1, 8)], &mut next);
        assert_eq!(
            got,
            vec!["  Family KING1 RULE FS1: B_X joins in sibship (A_C1 A_C2)\n"]
        );
        assert_eq!(next, 1, "a declared couple costs no synthetic pair");
    }

    /// `FS2` — two declared sibships combined — is unimplemented, and the component still
    /// gets its parents so `updateparents.txt` is unaffected.
    #[test]
    fn two_declared_sibships_in_one_component_raise_no_rule() {
        let samples = two_families();
        let mut next = 1;
        // A_C1 (0) and B_C1 (4) are inferred full sibs, so both declared sibships merge.
        let got = rule_lines(&samples, &[(0, 4)], &mut next);
        assert!(got.is_empty(), "{got:?}");
        let rows = assign(&samples, &[(0, 4)], &mut 1);
        assert_eq!(rows[4], ("B_C1".into(), "A_F".into(), "A_M".into()));
    }

    /// Nobody joined, nobody reassigned: a cluster merged by something other than a full
    /// sib pair leaves every parent column exactly as the `.fam` had it.
    #[test]
    fn a_cluster_with_no_inferred_sibship_assigns_nothing() {
        let samples = two_families();
        let members: Vec<usize> = (0..samples.len()).collect();
        let mut next = 1;
        let got = sibship_parents(&samples, &members, &|_, _| false, &mut next);
        assert!(got.parents.is_empty());
        assert!(got.rules.is_empty());
        assert_eq!(next, 1);
    }

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
