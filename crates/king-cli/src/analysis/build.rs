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
//! Of the three, **`<p>updateparents.txt` is implemented and byte-identical** — see the
//! section on it below. `<p>build.log` is not: it is left empty, because the `INFERENCE`
//! lines that dominate it carry a segment-derived statistic this engine cannot yet
//! reproduce, and a half-written log would be no closer to the capture than an empty one.
//! That is what remains of `docs/PARITY.md` §6.2.
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
//! ## The residual is not the formula, and the old bound argument for it is wrong
//!
//! An earlier reading of this module claimed the residual was **entirely** accounted for
//! by our sib-pair union over-call `ΔS`, since that can inflate `Join3` by at most `ΔS`
//! and so the ratio by at most `ΔS / Join2` — quoted as *39 of 39 triples inside*
//! `[0, ΔS / Join2]`. **That does not survive fresh shapes and seeds: it is 11 of 34.**
//! Worse, on the **8 triples where all three pair unions match the reference's reported
//! `IBD1Seg + IBD2Seg` to four decimals** the residual is still `+0.0003 … +0.0049`, mean
//! `+0.0031`. The bound is simply too crude to be evidence — `dU` is a *rounded total*,
//! and `Join3/Join2` is an intersection of three sets, so it reads segment **placement**,
//! which a matching total says nothing about.
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
//! Two further rules, measured the same way:
//!
//! * **Which two sibs are named is a property of the sibship, not of `R`.** Every
//!   `AV.FS` line raised against one sibship names the same pair whatever `R` is (three
//!   distinct `R` in one fixture, two in another), and where the sibship is the
//!   `RULE FS0`/`FS1` one it is that sibship's first two members in the order the rule
//!   line prints them.
//!
//!   For a *declared* sibship — the children of one `.fam` couple, which is what `bigish`
//!   names — the order is **not** the `.fam` order, and the earlier reading that it is
//!   "data-dependent" is wrong. It survives complete genotype reseeding: the 4:4 shape
//!   gives `(A_C3 A_C4)` and `(B_C1 B_C2)` on **all nine** seeds tried, and seven other
//!   shapes agree across three seeds each, though every printed `Join3/Join2` differs.
//!   It is also unchanged by each child's sex (five sex patterns, including all-male and
//!   all-female) and by sliding the whole pedigree down the `.fam` behind 0…8 extra
//!   singletons — so it is a *position* inside the sibship, not a sample index.
//!   What it is **not** a function of: the sibship's size alone, nor any pairwise
//!   statistic — over 19 measured triples the named pair's rank on `Join2`, `Join3`, the
//!   ratio itself, the sibs' mutual `PropIBD` and `Kinship`, and `PropIBD` to `R` each
//!   range from first to last. A four-child second family names positions `(1,2)`,
//!   `(1,3)`, `(3,4)` or `(3,2)` according to how many children the *first* family has
//!   and how many unrelated singletons pad the cohort — the one input that moves it
//!   while genotypes do not. `avfs_score.py --pairs` re-measures the map.
//! * **The verdict is a cut on the ratio.** Below it the line reads `<R> is uncle|aunt of
//!   N1 and N2`; above it, `<R> is grandfather|grandmother, HS, or nephew|niece of N1 and
//!   N2`, the word pairs following `R`'s sex. Over the same 53 values the cut is bracketed
//!   to **(0.848, 0.901)** — largest `uncle` 0.848, smallest ambiguous 0.901 — which does
//!   not separate 0.85, 0.875 and 0.9.
//!
//! # `<p>updateparents.txt` — implemented, and what pinned each clause
//!
//! This file needs none of the segment statistic: it carries only what the `RULE FS*`
//! lines decided. [`sibship_parents`] implements it and [`updateparents_text`] writes it.
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

    // The log, echoed to stdout and written to the file. Its `RULE` lines are known but
    // its `INFERENCE AV.FS` ones carry a segment-derived statistic this engine cannot yet
    // reproduce, so it is left empty rather than half written; see the module doc.
    let log = String::new();
    let _ = out.write_all(log.as_bytes());
    let _ = out.write_all(b"\n");

    let parents = updateparents_text(opts, loaded, &clustering);

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

/// `<prefix>updateparents.txt`, or `None` when reconstruction assigned nobody a parent.
///
/// `None` is the whole-run verdict the console tail prints as `No pedigrees can be
/// reconstructed.`; the caller still writes the (empty) file, exactly as the reference
/// leaves a zero-byte one behind. The clause-by-clause evidence is in the module doc.
fn updateparents_text(
    opts: &Options,
    loaded: &Loaded,
    clustering: &unrelated::Clustering,
) -> Option<String> {
    let samples = &loaded.fileset.samples;
    let types = InfTypes::new(opts, loaded);

    let mut text = String::new();
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
        reconstructed |= !assigned.is_empty();
        for (n, &i) in members.iter().enumerate() {
            let (pat, mat) = assigned
                .iter()
                .find(|(m, _)| *m == n)
                .map(|(_, p)| (p.0.as_str(), p.1.as_str()))
                .unwrap_or((samples[i].pat.as_str(), samples[i].mat.as_str()));
            text.push_str(&format!("{key}\t{}\t{pat}\t{mat}\n", samples[i].iid));
        }
    }
    reconstructed.then_some(text)
}

/// Every unordered pair of positions in `members`, as index pairs into it.
fn pairs(members: &[usize]) -> impl Iterator<Item = (usize, usize)> + '_ {
    (0..members.len())
        .flat_map(move |a| (a + 1..members.len()).map(move |b| (members[a], members[b])))
}

/// The parent pair each sibship in one cluster takes, as `(position in members, pair)`.
///
/// A sibship is a connected component of `inferred FS` ∪ `declares the same non-missing
/// couple`, and only a component the *inference* touched is given parents — a declared
/// sibship nobody was joined to keeps what the `.fam` already says. `next_synthetic`
/// carries the `1 2`, `3 4`, … counter across clusters, so it advances only where a
/// sibship actually needed inventing.
fn sibship_parents(
    samples: &[Sample],
    members: &[usize],
    is_fs: &dyn Fn(usize, usize) -> bool,
    next_synthetic: &mut u32,
) -> Vec<(usize, (String, String))> {
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
    let mut out: Vec<(usize, (String, String))> = Vec::new();
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
        for b in group {
            out.push((b, couple.clone()));
        }
    }
    out
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

    /// Nobody joined, nobody reassigned: a cluster merged by something other than a full
    /// sib pair leaves every parent column exactly as the `.fam` had it.
    #[test]
    fn a_cluster_with_no_inferred_sibship_assigns_nothing() {
        let samples = two_families();
        let members: Vec<usize> = (0..samples.len()).collect();
        let mut next = 1;
        assert!(sibship_parents(&samples, &members, &|_, _| false, &mut next).is_empty());
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
