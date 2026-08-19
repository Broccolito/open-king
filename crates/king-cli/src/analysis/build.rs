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
//! All three outputs are implemented. The primary `bigish` parity case is byte-identical
//! across stdout, `build.log`, `updateids.txt`, and `updateparents.txt`. A wider research
//! replay is intentionally retained below to identify harmless reference repetition and
//! rare pedigree shapes beyond that shipped corpus case.
//!
//! # `<p>build.log` — the format, and which half is written
//!
//! The file is a per-cluster narration of the reconstruction. It is written to disk and
//! echoed to **stdout byte for byte**: the console block between the segment prepass and
//! `Details of pedigree reconstruction …` is the file's contents verbatim, plus one blank
//! line after it that the file does not carry (`bigish`'s capture ends `0.803\n`, 806
//! bytes, and its stdout repeats those 18 lines exactly). A cluster contributes nothing at
//! all unless it raises a line, and when it does, its first line is the header
//! `Family KING<k>:`.
//!
//! Inside a cluster the order is fixed: the `RULE` lines that build sibships, then one or
//! more blank lines, then the `INFERENCE` block. The templates, all of them mined from the
//! reference binary's `__cstring` section and then confirmed against runs, are
//!
//! ```text
//! Family KING1:
//!   Duplicate QBC_D (of QBB_C1) is removed.
//!   Family KING1 RULE FS0: Sibship (B01_F B02_F)'s parents are (1 2)
//!   Family KING1 RULE FS1: C_F joins in sibship (A_F B_F)
//!   Family KING1 RULE FS2: Sibship (…) and sibship (…) are combined
//!   Reconstruct parent-offspring pair (C_F, D_F)...
//!   RBA_F's sibship is used to determine the parent/offspring
//!   Family KING1 RULE PO.S: RBA_F is now father of RBC_F
//!     3 is created as RBC_F's mother.
//!
//!   Family KING1 INFERENCE AV.FS: B02_F is uncle of B01_C2 and B01_C3, Join3/Join2=0.778
//!   Family KING1 INFERENCE AV.FS: A_F is grandfather, HS, or nephew of B_C2 and B_C3, …
//!   Family KING1 INFERENCE AV.HS: A_F is uncle of (B_C1 D_C3), Join3/Join2=0.699
//!     HS B02_C4 unrelated to B01_M
//!   Family KING1 INFERENCE HS.UN2: B01_C3 and B02_C4 are HS
//! ```
//!
//! ## Which template is in which half, and what triggers it
//!
//! The split still matters to emission order and to the condition that suppresses a whole
//! PO-only block, even though this module now writes both halves.
//!
//! | template | half | trigger |
//! |---|---|---|
//! | `Family KING<k>:` | — | printed once, before the cluster's first line, whatever raises it |
//! | `Duplicate <a> (of <b>) is removed.` | **rule** | an inferred `Dup/MZ` pair inside the cluster, *provided the cluster raises something else* |
//! | `RULE FS0` | **rule** | a component of *inferred FS* ∪ *declares the same couple* that the inference **created** |
//! | `RULE FS1` | **rule** | one more member joining a component that already contained a declared sibship, or one `FS0` just made |
//! | `RULE FS2` | rule | two *declared* sibships in one component; needs a mis-declaring `.fam`, and does fire — see below |
//! | `Reconstruct parent-offspring pair (X, Y)...` | **inference** | an inferred `PO` pair, but only inside a cluster whose inference block also speaks |
//! | `<X>'s sibship is used to determine the parent/offspring`, `RULE PO.S`, `<n> is created as <Y>'s mother.` | inference | the `PO` pair above, when one member's sibship orients it |
//! | `INFERENCE AV.FS`, `uncle\|aunt` form | inference | an `R` inferred 2nd-degree to **both** named members of a sibship, with `Join3/Join2 < 0.85` |
//! | `INFERENCE AV.FS`, `grandfather\|…, HS, or nephew\|…` form | inference | the same candidate with `Join3/Join2 > 0.90` |
//! | *(no line at all)* | — | the same candidate with `0.85 <= Join3/Join2 <= 0.90` — the dead band, see below |
//! | `INFERENCE AV.HS: <R> is <verdict> of (<a> <b>), …` | inference | the same statistic and the same band, run against a candidate **HS pair** instead of a sibship |
//! | `HS <a> unrelated to <b>` | inference | a candidate HS pair (`PropIBD > 0.1875`), one line per declared parent the other member is unrelated to |
//! | `INFERENCE HS.UN2` | inference | that pair, when both sides produced such a line |
//!
//! **`Reconstruct parent-offspring pair` is not a rule line**, which is the correction
//! that matters most here. It looks like one — no `Family KING<k>` prefix, sitting above
//! the blank lines, raised by a `PO` merge — and this module briefly wrote it on that
//! reading. Two measurements say otherwise: across every reference log these rigs have
//! produced, **42 of 42** clusters that print it also print an `INFERENCE` line; and a `PO`
//! merge between two families with **no sibship anywhere** — two one-person families, and
//! two childless couples, three seeds each — prints nothing at all, not even a header
//! (`work/poprobe`). A `PO`-joined cluster still contributes its identity rows to
//! `<p>updateparents.txt`; it just does not narrate them here.
//!
//! `Duplicate … is removed.` is the opposite case and **is** written. `dupkeep.py` scores
//! it over ten shapes × three seeds: **23 of 27** runs print it in a file with no
//! `INFERENCE` line anywhere, so it is rule-half beyond doubt. See [`duplicate_verdict`]
//! for which of the two copies goes.
//!
//! ## What this module writes, and how it scores
//!
//! The complete cached research replay currently has **277 of 347** logs byte-identical.
//! Another **52** have the same distinct lines and differ only in repetition/count residue
//! from the reference's unstable inference loop. Of the remaining 18, two cached reference
//! logs were truncated during debugger work and two exercise the deliberately excluded
//! cross-FID parent renaming (`<FID>-><IID>`). The remaining rare shapes are documented as
//! held-out behavioral residuals rather than hidden behind the primary corpus score.
//!
//! ## `bigish`'s eighteen lines, one at a time
//!
//! All eighteen lines now land byte for byte, including the five AV.FS statistics, the
//! HS.UN2 block, and the four blank lines. This is pinned by the parity harness rather
//! than inferred from a reduced line set.
//!
//! # `INFERENCE AV.FS` and its `Join3/Join2` — solved
//!
//! The one statistic the log needs is **segment-derived**, and the set it is computed
//! over is the one the pair **reports**, not the one the caller finds. For an ordered
//! triple `(R; N1, N2)` write `S(x, y)` for the base pairs covered by that pair's
//! *reported* segments — every IBD2 call, plus the IBD1 pieces left when the IBD2 calls
//! are cut out at marker granularity, each piece having cleared `--seglength` on its own.
//! That is exactly the interval set `<prefix>.seg` would list, and exactly what
//! [`reported_intervals`] builds. Then
//!
//! ```text
//! Join2 = | S(R,N1) INTERSECT S(R,N2) |
//! Join3 = | S(R,N1) INTERSECT S(R,N2) INTERSECT S(N1,N2) |
//! ```
//!
//! and the log prints `Join3/Join2` at `%.3lf` ([`join_ratio`]). The genetics behind it:
//! where `R` is IBD to both sibs, a *grandparent* forces the sibs to have inherited the
//! same parental haplotype, so the ratio is 1; an *avuncular* does not, so it sits near
//! 2/3.
//!
//! ## What the correction was, and how it scores
//!
//! The previous reading of this doc used the **raw calls** — every IBD1 call and every
//! IBD2 call, unified — and that set is one-sided large: it keeps the sub-`--seglength`
//! IBD1 fringe the pair never reports. It landed **1 of 34** triples on the printed three
//! decimals, mean `+0.0039`, and the residual was written up here as a caller error that
//! no variant of the formula would remove. It was a formula error. Cutting the IBD2 calls
//! out of the IBD1 ones and floor-testing each piece — the same two clauses `IBD1Seg`
//! itself is summed under — removes it completely:
//!
//! | interval set | exact at `%.3lf` |
//! |---|---|
//! | raw IBD1 and IBD2 calls unified (the old reading) | 13 of 297 |
//! | **reported segments** (this one) | **296 of 297** |
//! | reported segments, endpoints pushed to the marker-gap midpoint | 0 of 5 |
//! | IBD2 calls alone | undefined, `Join2` is 0 |
//!
//! Those 297 are every `AV.FS` line every rig under `docs/research/fixtures/` has
//! captured — `bigish` plus roughly 120 purpose-built filesets across `avfs.py`,
//! `clusternum.py`, `build_shapes.py`, `dupkeep.py`, the genotype-surgery sweeps and the
//! fresh-seed battery. The set was chosen on `bigish`'s five and then held out on the
//! other 292. `bigish`'s five now read 0.7784, 0.8007, 0.7793, 0.8266, 0.8027 against the
//! capture's 0.778, 0.801, 0.779, 0.827, 0.803 — **five of five**, and no triple anywhere
//! is off by more than 0.0005.
//!
//! The single miss is not ours. `clusternum.py`'s `po_fs` prints `Join3/Join2=2.555` — a
//! ratio above 1, which no intersection of `Join3` inside `Join2` can produce. It is the
//! only capture whose `AV.FS` line sits under a `Reconstruct parent-offspring pair (…)`
//! rather than a `RULE FS0`, so whatever that branch divides by, it is not this `Join2`.
//! One case in 297; recorded, not chased.
//!
//! # The verdict is a **band**, not a cut — and that is why lines go missing
//!
//! The rule this doc used to carry, "below the cut it reads `uncle`, above it reads
//! `grandfather, HS, or nephew`", is wrong in a way that hid an entire third branch:
//!
//! ```text
//! ratio < 0.85   ->  "<R> is uncle|aunt of N1 and N2"
//! ratio > 0.90   ->  "<R> is grandfather|grandmother, HS, or nephew|niece of N1 and N2"
//! otherwise      ->  nothing is printed at all
//! ```
//!
//! ([`av_verdict`].) The word pairs follow `R`'s sex. The dead band is what explains every
//! `AV.FS` line a candidate ought to have raised and did not — `bigish`'s `KING3`, where
//! `B26_F` is 2nd-degree to all three `B25` children and stays silent, is a ratio of 0.889
//! falling inside it.
//!
//! ## How both edges were bisected
//!
//! Not by collecting values and looking at the extremes — by walking one triple through
//! each edge with a continuous knob. `docs/research/fixtures/work/bandcut.py` overwrites
//! `N2`'s calls with `N1`'s over the first `f` of the markers: inside that stretch the sib
//! pair is IBD everywhere, so it joins `Join3` wherever `R` already shared with `N1`, and
//! the ratio walks up while the `.fam`, the ids and the pedigree stay put. Because
//! [`join_ratio`] is now exact, the bracket carries **no rounding slop at all**:
//!
//! | edge | last value one side | first value the other | bracket |
//! |---|---|---|---|
//! | `uncle` / silent | 0.848718 prints `uncle 0.849` | 0.851164 is silent | **(0.848718, 0.851164]** |
//! | silent / ambiguous | 0.896895 is silent | 0.900106 prints `… 0.900` | **(0.896895, 0.900106]** |
//!
//! `0.85` and `0.90` are the only round constants inside, and the previous single-cut
//! bracket `(0.846, 0.902)` is superseded twice over: it was one interval where there are
//! two, and each of the two is now about 0.003 wide instead of 0.056.
//!
//! Held out afterwards on 14 fresh-seed filesets built for it (`battery.py band`), each
//! carrying two two-child sibships so the named pair is forced: **27 of 27** verdicts
//! agree with the band — 14 `uncle`, 5 ambiguous, 8 silent, every silent one inside
//! `[0.8532, 0.8902]`. (The 28th is a cluster whose only `FS` pair reconstruction
//! rejected, so it raises no candidate at all and carries no evidence either way.)
//!
//! # The `HS` half — trigger found
//!
//! The three remaining templates are one pass over **candidate half-sib pairs**, and a
//! candidate is a cross-family pair whose segments put it *well inside* 2nd degree:
//!
//! ```text
//! PropIBD > 0.1875          (equivalently IBD1Seg > 0.375; every candidate has IBD2Seg 0)
//! ```
//!
//! Bracketed to **(0.1868, 0.1878)** on 13 held-out candidate pairs from fresh five-child
//! shapes plus the 13 the older rigs had left behind — 26 pairs, 0 refutations. It is not
//! a kinship test: over the same 26 the `Kinship` column **overlaps** (silent as high as
//! 0.0932, firing as low as 0.0897), so `2^-3.5` and every other kinship constant is
//! refuted outright. `bigish` is decided by this gate: `B01_C3`/`B02_C4` at 0.1901 fires
//! and `B01_C2`/`B02_C2` at 0.1799 does not, which is why `KING1` carries exactly one
//! `HS.UN2` and not two.
//!
//! For each candidate `(p1, p2)` in id order the pass walks `p1`'s declared parents and
//! prints one line per parent the *other* member is inferred unrelated to, then does the
//! same for `p2`:
//!
//! ```text
//!     HS <p2> unrelated to <a parent of p1>
//!     HS <p1> unrelated to <a parent of p2>
//!   Family KING1 INFERENCE HS.UN2: <p1> and <p2> are HS
//! ```
//!
//! The fathers never appear because they are the cluster's `FS` pair and so *are* related
//! to the other's children; the mothers do. `HS.UN2` closes the block when both sides
//! produced a line — "UN2" is the two `unrelated` findings. The clean confirmation is a
//! constructed one: `battery.py hsrel` makes the mothers full sibs too (double first
//! cousins), so **no** parent check can pass, and the log carries nine candidate pairs
//! above the gate, two `RULE FS0` lines, three `AV.FS` lines — and not one `HS` line.
//!
//! ## `INFERENCE AV.HS` is `AV.FS` with the pair in place of the sibship
//!
//! Ahead of that block, every `R` inferred 2nd-degree to **both** members of the candidate
//! pair raises
//!
//! ```text
//!   Family KING1 INFERENCE AV.HS: A_F is uncle of (B_C1 D_C3), Join3/Join2=0.699
//! ```
//!
//! — the pair in parentheses, in the same order the `HS.UN2` line will name it, and the
//! same [`join_ratio`] under the same [`av_verdict`] band. `buildlog/m43`'s two read
//! 0.699 and 0.587 and [`join_ratio`] returns 0.699224 and 0.586998, so this template needs
//! nothing new: the pair it names is *given* by the HS candidate, not by a sibship order,
//! which makes `AV.HS` the one inference template that is fully implementable today.
//!
//! That fileset also carries the tightest **printed** pair of verdicts anywhere —
//! `uncle … 0.846` and `grandfather, HS, or nephew … 0.903` — which brackets the two cuts
//! to `(0.8455, …)` and `(…, 0.9035)` from printed values alone, independently of the
//! genotype-surgery bisection and agreeing with it.
//!
//! # `RULE FS2` — made to fire
//!
//! Never seen until now. It needs **two declared sibships in one component**, which no
//! honest pedigree can produce: an inferred `FS` pair with one end in each sibship means
//! the two ends share both parents, so at least one family's `.fam` must mis-declare them.
//! Simulate four true full sibs, declare two under one couple and two under another, and
//! the reference prints
//!
//! ```text
//! Family KING1:
//!   Family KING1 RULE FS2: Sibship (K_C1 K_C2) and sibship (K_C4 K_C3) are combined
//! ```
//!
//! — and nothing else: no parents line, no `FS0`, no inference. Note `(K_C4 K_C3)`, which
//! is the same internal order the next section is about.
//!
//! # The internal sibling order — solved from the actual library code
//!
//! KING first sorts the reconstructed pedigree by `(FID, IID)`, then
//! `Pedigree::MakeSibships` applies libStatGen's deliberately unstable `QuickSort` by
//! `(FID, father, mother)`. Equal-parent siblings retain the deterministic swap residue of
//! that second sort. [`reconstructed_pedigree_order`] and [`unstable_king_sort`] reproduce
//! it exactly; the five-child sentinel is `(A_C4 A_C1 A_C5 A_C3 A_C2)`.
//!
//! # The blank lines — the `reject` rule is refuted
//!
//! Their count was previously fitted to `1 + the number of candidates examined and turned
//! down before the first line prints` (107 of 113). That rule is now **refuted by a
//! matched pair**: `battery.py band`'s `bnd07` and `bnd09` are the same shape, the same
//! two two-child sibships, one candidate `R` each, and *both candidates silent in both
//! filesets* — and they print **3** blanks and **2**. No function of the verdicts can
//! separate them, and no threshold on the ratio can either, since `bnd07`'s pair straddles
//! `bnd09`'s: 0.8532, 0.8902 against 0.8810, 0.8711.
//!
//! Whatever counts the blanks is the same unidentified iteration that decides how often a
//! line **repeats**: the reference prints an identical `AV.FS` line 1 to 6 times, always
//! the last candidate of a sibship block, and the count is not the number of sib pairs,
//! not the number of candidates and not the number of families. `bigish` happens to sit at
//! 1 everywhere, so this unstable repetition does not affect the primary parity case.
//!
//! # `<p>updateparents.txt` — implemented, and what pinned each clause
//!
//! This file needs none of the segment statistic: it carries only what the `RULE FS*`
//! lines decided. [`sibship_parents`] implements it and [`reconstruct`] writes it, in the
//! same walk that writes those rule lines, so the two cannot disagree.
//! Every clause below was measured against the reference on **twenty held-out merge
//! shapes** built for it, not on `bigish`; the rig is
//! `docs/research/fixtures/build_shapes.py`, which re-runs the scorecard end to end.
//! Eighteen of the twenty are in scope and **all eighteen are byte-identical on the file
//! and on the console tail** — the three that used to differ were the cluster-numbering
//! shapes, and that bug is fixed in [`unrelated::clusters`].
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
//! * **A duplicate is removed, not fatal to the cluster.** One copy drops out of every
//!   sibship — [`duplicate_verdict`] decides which — and the rest of the cluster
//!   reconstructs around it, the dropped copy keeping its `.fam` parents in its own row.
//!   The rule this doc used to carry, "a cluster holding an inferred duplicate contributes
//!   no rows at all", is the *special case* where removing the duplicate leaves nothing
//!   behind: a pure `Dup/MZ` merge has no `FS` or `PO` pair left afterwards, so it raises
//!   no line and contributes no row, the removal line included. A cluster that carries a
//!   duplicate *and* an `FS` pair contributes rows for all of its members, the removed
//!   copy included (`clusternum.py`'s `mixed_cluster`, 13 rows). A cluster merged by `PO`
//!   alone also contributes them, with nobody's parents changed, because without ages the
//!   reference will not orient the pair.
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
//! # The clustering bug this rig found — now fixed in `unrelated.rs`
//!
//! Merged clusters are **not** numbered in family order but in the order a **staged merge
//! queue** creates them, duplicates first, then parent–offspring, then full sibs. The rule
//! and the nineteen held-out shapes that pinned it are documented on
//! [`unrelated::clusters`]; the corpus is indifferent to it — every merge in `bigish` is
//! `FS` and no other corpus fileset merges at all — so no capture moved either way, but
//! six of the shapes here did.

use std::collections::HashMap;
use std::io::Write;
use std::path::Path;

use king_core::ibdseg::Called;
use king_core::infer::{pedigree_kinship, KinshipCache, Pedigree};
use king_io::Sample;

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

/// `Join3/Join2` at or above which an `INFERENCE AV.FS` candidate stops reading `uncle`.
///
/// Bisected to `(0.848718, 0.851164]` with a continuous genotype knob on one triple, with
/// no rounding slop because [`join_ratio`] reproduces the reference's own value — see the
/// module doc. `0.85` is the only round constant inside.
pub const AV_UNCLE_MAX: f64 = 0.85;

/// `Join3/Join2` above which an `INFERENCE AV.FS` candidate reads
/// `grandfather|grandmother, HS, or nephew|niece`.
///
/// Bisected to `(0.896895, 0.900106]` the same way. Between this and [`AV_UNCLE_MAX`] the
/// reference prints **nothing**, which is the branch that had been missed entirely.
pub const AV_AMBIGUOUS_MIN: f64 = 0.90;

/// `PropIBD` above which a cross-family 2nd-degree pair becomes a **half-sib candidate**
/// — the pair that can raise `HS <a> unrelated to <b>` and `INFERENCE HS.UN2`.
///
/// `3/16`, bracketed to `(0.1868, 0.1878)` on 26 candidate pairs, 13 of them from fresh
/// five-child shapes built for it. Equivalently `IBD1Seg > 0.375`: every candidate seen
/// carries `IBD2Seg 0`, so nothing here separates the two readings. The `Kinship` column
/// is refuted as the test — over the same 26 its silent and firing ranges overlap.
pub const HS_CANDIDATE_PROP_IBD: f64 = 0.1875;

/// The base pairs a pair **reports** as IBD, which is the set `Join3/Join2` is built over.
///
/// Every IBD2 call, plus the pieces of each IBD1 call left once the IBD2 calls are cut out
/// — at marker granularity, excluding the IBD2 call's own end markers — with each piece
/// facing `seglength_bp` on its own. Those are the two clauses `IBD1Seg` is summed under,
/// so this is exactly the interval set `<prefix>.seg` lists for the pair.
///
/// Using the **raw** calls instead is the error that made the old `Join3/Join2` read
/// one-sided high by ~0.004: the raw set keeps the sub-`--seglength` IBD1 fringe that the
/// pair never reports. See the module doc for the 297-triple scorecard.
///
/// `ibd1` and `ibd2` are one pair's merged calls over one usable segment, as
/// [`king_core::ibdseg::Scan`] returns them; `pos` is the autosomal position array the
/// scan was run against. The result is sorted and non-overlapping.
pub fn reported_intervals(
    pos: &[i64],
    ibd1: &[Called],
    ibd2: &[Called],
    seglength_bp: i64,
) -> Vec<(i64, i64)> {
    let mut out: Vec<(i64, i64)> = ibd2.iter().map(|c| (pos[c.lo], pos[c.hi])).collect();
    for call in ibd1 {
        let mut cur = call.lo;
        let piece = |lo: usize, hi: usize, out: &mut Vec<(i64, i64)>| {
            if lo <= hi && pos[hi] - pos[lo] >= seglength_bp {
                out.push((pos[lo], pos[hi]));
            }
        };
        for o in ibd2 {
            if o.hi < call.lo || o.lo > call.hi {
                continue;
            }
            if o.lo > cur {
                piece(cur, o.lo - 1, &mut out);
            }
            cur = cur.max(o.hi + 1);
        }
        if cur <= call.hi {
            piece(cur, call.hi, &mut out);
        }
    }
    out.sort_unstable();
    let mut merged: Vec<(i64, i64)> = Vec::new();
    for iv in out {
        match merged.last_mut() {
            Some(last) if iv.0 <= last.1 => last.1 = last.1.max(iv.1),
            _ => merged.push(iv),
        }
    }
    merged
}

/// Intersect two sorted, non-overlapping interval sets.
fn intersect(a: &[(i64, i64)], b: &[(i64, i64)]) -> Vec<(i64, i64)> {
    let (mut x, mut y, mut out) = (0, 0, Vec::new());
    while x < a.len() && y < b.len() {
        let (lo, hi) = (a[x].0.max(b[y].0), a[x].1.min(b[y].1));
        if lo < hi {
            out.push((lo, hi));
        }
        if a[x].1 < b[y].1 {
            x += 1;
        } else {
            y += 1;
        }
    }
    out
}

/// `Join3/Join2` for the triple `(R; N1, N2)`, or `None` when `Join2` is empty.
///
/// The three arguments are [`reported_intervals`] for `(R,N1)`, `(R,N2)` and `(N1,N2)`,
/// concatenated across every usable segment. Positions must be comparable across
/// chromosomes — lay them end to end first, since a `.bim` restarts at 1 on every one.
pub fn join_ratio(r_n1: &[(i64, i64)], r_n2: &[(i64, i64)], n1_n2: &[(i64, i64)]) -> Option<f64> {
    let total = |v: &[(i64, i64)]| v.iter().map(|(a, b)| b - a).sum::<i64>();
    let join2 = intersect(r_n1, r_n2);
    let denom = total(&join2);
    (denom != 0).then(|| total(&intersect(&join2, n1_n2)) as f64 / denom as f64)
}

/// What an `INFERENCE AV.FS` candidate reads, given its ratio and `R`'s sex — or `None`
/// when the ratio lands in the dead band and the reference prints no line at all.
///
/// `sex` is the `.fam` code: 2 is female, anything else takes the masculine wording (the
/// reference has no third form, and every capture with `sex 0` reads masculine).
pub fn av_verdict(ratio: f64, sex: u8) -> Option<&'static str> {
    if ratio < AV_UNCLE_MAX {
        Some(if sex == 2 { "aunt" } else { "uncle" })
    } else if ratio > AV_AMBIGUOUS_MIN {
        Some(if sex == 2 {
            "grandmother, HS, or niece"
        } else {
            "grandfather, HS, or nephew"
        })
    } else {
        None
    }
}

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
    let types = unrelated::InfTypes::new(opts, loaded);
    let segments = ibdseg::Segments::new(opts, loaded);
    let sibship_order = reconstructed_pedigree_order(samples, clustering);

    let mut text = String::new();
    let mut log = String::new();
    let mut next_synthetic = 1u32;
    let mut reconstructed = false;

    for (key, members) in clustering.merged() {
        // The cluster's own 1st-degree pairs, as positions in `members`, split by what
        // reconstruction does with each kind.
        let (mut dups, mut po, mut fs) = (Vec::new(), Vec::new(), Vec::new());
        for a in 0..members.len() {
            for b in a + 1..members.len() {
                match types.first_degree(loaded, members[a], members[b]) {
                    Some("Dup/MZ") => dups.push((a, b)),
                    Some("PO") if clustering.is_joining_pair(members[a], members[b]) => {
                        po.push((a, b));
                    }
                    Some("FS") => fs.push((a, b)),
                    _ => {}
                }
            }
        }
        // A duplicate is *removed*, not a reason to abandon the cluster: one copy drops
        // out of every sibship and the rest reconstructs around it. What a pure duplicate
        // merge leaves behind is nothing to reconstruct, which is why such a cluster
        // raises no line and contributes no row — the removal line goes with it.
        let dups: Vec<(usize, usize)> = dups
            .into_iter()
            .map(|(a, b)| duplicate_verdict(samples, members[a], members[b]))
            .collect();
        let gone: Vec<usize> = dups.iter().map(|&(_, r)| r).collect();
        let live =
            |&(a, b): &(usize, usize)| !gone.contains(&members[a]) && !gone.contains(&members[b]);
        fs.retain(live);
        po.retain(live);
        // `PO` alone raises no rule line — the narration of a parent-offspring pair is
        // part of the inference half, see the module doc — but it does still make the
        // cluster contribute its identity rows.
        if fs.is_empty() && po.is_empty() {
            continue;
        }
        let is_fs = |i: usize, j: usize| {
            !gone.contains(&i)
                && !gone.contains(&j)
                && types.first_degree(loaded, i, j) == Some("FS")
        };
        let mut assigned = sibship_parents(
            samples,
            members,
            &is_fs,
            &sibship_order,
            &mut next_synthetic,
        );
        // A joining PO edge can already be explained by a declared parent of the FS
        // component it helped connect (for example a grandfather against the sibling
        // who just joined that parent's child). KING leaves that consistent pedigree
        // alone instead of orienting the edge a second time as PO.S.
        po.retain(|&(a, b)| {
            let pa = assigned.parents.iter().find(|(m, _)| *m == a);
            let pb = assigned.parents.iter().find(|(m, _)| *m == b);
            let (assigned_couple, outsider) = match (pa, pb) {
                (Some((_, p)), None) => (p, b),
                (None, Some((_, p))) => (p, a),
                _ => return true,
            };
            let outsider_id = &samples[members[outsider]].iid;
            !assigned.parents.iter().any(|(position, parents)| {
                parents == assigned_couple
                    && (samples[members[*position]].pat == *outsider_id
                        || samples[members[*position]].mat == *outsider_id)
            })
        });
        let mut po_orientations = Vec::new();
        for &(a, b) in &po {
            let assigned_a = assigned.parents.iter().any(|(m, _)| *m == a);
            let assigned_b = assigned.parents.iter().any(|(m, _)| *m == b);
            let (parent, child) = match (assigned_a, assigned_b) {
                (true, false) => (a, b),
                (false, true) => (b, a),
                _ => continue,
            };
            let created = next_synthetic.to_string();
            next_synthetic += 1;
            let (pat, mat, parent_role, created_role) = match samples[members[parent]].sex {
                1 => (
                    samples[members[parent]].iid.clone(),
                    created.clone(),
                    "father",
                    "mother",
                ),
                2 => (
                    created.clone(),
                    samples[members[parent]].iid.clone(),
                    "mother",
                    "father",
                ),
                _ => continue,
            };
            assigned.parents.push((child, (pat, mat)));
            po_orientations.push((a, b, parent, child, parent_role, created, created_role));
        }
        reconstructed |= !assigned.parents.is_empty();
        let name = |n: &usize| samples[members[*n]].iid.as_str();
        let iid = |i: usize| samples[i].iid.as_str();
        let mut lines = String::new();
        for &(kept, removed) in &dups {
            lines.push_str(&console::build_duplicate_removed(iid(removed), iid(kept)));
        }
        for rule in &assigned.rules {
            let names: Vec<&str> = rule.members.iter().map(name).collect();
            lines.push_str(&match (&rule.second, rule.joiner) {
                (Some(second), _) => {
                    let second_names: Vec<&str> = second.iter().map(name).collect();
                    console::build_rule_fs2(key, &names, &second_names)
                }
                (None, None) => console::build_rule_fs0(key, &names, &rule.pat, &rule.mat),
                (None, Some(j)) => console::build_rule_fs1(key, name(&j), &names),
            });
        }
        let effective_parents: HashMap<usize, (String, String)> = assigned
            .parents
            .iter()
            .map(|(position, parents)| (members[*position], parents.clone()))
            .collect();
        let inference = segments.as_ref().map_or_else(String::new, |segments| {
            let mut reconstructed_sibships: Vec<Vec<usize>> = Vec::new();
            for rule in &assigned.rules {
                if rule.second.is_some() {
                    continue;
                }
                match rule.joiner {
                    None => reconstructed_sibships
                        .push(rule.members.iter().map(|p| members[*p]).collect()),
                    Some(joiner) => {
                        if let Some(group) = reconstructed_sibships
                            .iter_mut()
                            .find(|group| rule.members.iter().all(|p| group.contains(&members[*p])))
                        {
                            group.push(members[joiner]);
                        }
                    }
                }
            }
            inference_lines(
                key,
                members,
                loaded,
                &types,
                segments,
                &InferencePedigree {
                    sibship_order: &sibship_order,
                    reconstructed_sibships: &reconstructed_sibships,
                    allow_hs_without_av: !assigned.rules.is_empty(),
                    effective_parents: &effective_parents,
                },
            )
        });
        reconstructed |= !po.is_empty() && !inference.trim().is_empty();
        if !inference.trim().is_empty() {
            for &(a, b) in &po {
                lines.push_str(&console::build_reconstruct_po(name(&a), name(&b)));
                if let Some((_, _, parent, child, parent_role, created, created_role)) =
                    po_orientations.iter().find(|(x, y, ..)| *x == a && *y == b)
                {
                    lines.push_str(&console::build_rule_po_s(
                        key,
                        name(parent),
                        name(child),
                        parent_role,
                        created,
                        created_role,
                    ));
                }
            }
        }
        if !lines.is_empty() || !inference.trim().is_empty() {
            log.push_str(&console::build_family_header(key));
            log.push_str(&lines);
            log.push_str(&inference);
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

/// KING's internal person order after clustering has rewritten connected FIDs.
///
/// The reference loads a temporary pedigree whose merged families are already named
/// `KING<n>`, sorts its persons by `(FID, IID)`, then `Pedigree::MakeSibships` applies
/// libStatGen's deliberately unstable quicksort by `(FID, father, mother)`.  Equal-parent
/// members therefore acquire a deterministic but non-obvious order.  It is the order
/// shared by `RULE FS1`, `RULE FS2`, and the named pair in `INFERENCE AV.FS`.
fn reconstructed_pedigree_order(
    samples: &[Sample],
    clustering: &unrelated::Clustering,
) -> Vec<usize> {
    let mut fid: Vec<&str> = samples.iter().map(|s| s.fid.as_str()).collect();
    for (key, members) in clustering.merged() {
        for &i in members {
            fid[i] = key;
        }
    }
    let cmp = |a: &[u8], b: &[u8]| crate::analysis::king_id_cmp(a, b);
    let mut order: Vec<usize> = (0..samples.len()).collect();
    unstable_king_sort(&mut order, |a, b| {
        cmp(fid[a].as_bytes(), fid[b].as_bytes())
            .then_with(|| cmp(samples[a].iid.as_bytes(), samples[b].iid.as_bytes()))
    });
    unstable_king_sort(&mut order, |a, b| {
        cmp(fid[a].as_bytes(), fid[b].as_bytes())
            .then_with(|| cmp(samples[a].pat.as_bytes(), samples[b].pat.as_bytes()))
            .then_with(|| cmp(samples[a].mat.as_bytes(), samples[b].mat.as_bytes()))
    });
    order
}

/// libStatGen's `QuickSort`, including its instability for equal keys.
///
/// A stable Rust sort is observably wrong here: it changes the sibling named in pedigree
/// reconstruction even though every relationship and statistic is unchanged.
fn unstable_king_sort(values: &mut [usize], cmp: impl Fn(usize, usize) -> std::cmp::Ordering) {
    if values.len() <= 1 {
        return;
    }
    let mut stack: Vec<(usize, usize)> = Vec::new();
    let (mut l, mut r) = (0usize, values.len() - 1);
    loop {
        while r > l {
            let (pivot, mut scanl, mut scanr);
            if r - l > 7 {
                let mid = (r + l) / 2;
                if cmp(values[mid], values[l]).is_lt() {
                    values.swap(mid, l);
                }
                if cmp(values[r], values[l]).is_lt() {
                    values.swap(r, l);
                }
                if cmp(values[r], values[mid]).is_lt() {
                    values.swap(r, mid);
                }
                pivot = r - 1;
                values.swap(mid, pivot);
                scanl = l + 1;
                scanr = r - 2;
            } else {
                pivot = r;
                scanl = l;
                scanr = r - 1;
            }
            loop {
                while scanl < r && cmp(values[scanl], values[pivot]).is_lt() {
                    scanl += 1;
                }
                while scanr > l && cmp(values[pivot], values[scanr]).is_lt() {
                    scanr -= 1;
                }
                if scanl >= scanr {
                    break;
                }
                values.swap(scanl, scanr);
                if scanl < r {
                    scanl += 1;
                }
                if scanr > l {
                    scanr -= 1;
                }
            }
            values.swap(pivot, scanl);
            let (lsize, rsize) = (scanl - l, r - scanl);
            if lsize > rsize {
                stack.push((l, scanl - 1));
                if rsize != 0 {
                    l = scanl + 1;
                } else {
                    break;
                }
            } else {
                stack.push((scanl + 1, r));
                if lsize != 0 {
                    r = scanl - 1;
                } else {
                    break;
                }
            }
        }
        if let Some((left, right)) = stack.pop() {
            l = left;
            r = right;
        } else {
            break;
        }
    }
}

/// The inference half of one merged family's reconstruction log.
///
/// This is deliberately driven from the same unrounded segment calls and inferred labels
/// as clustering.  Recomputing from formatted `.kin0`/`.seg` rows would move candidates
/// at the measured boundaries and would make `--build` depend on files another pass may
/// not have written.
struct InferencePedigree<'a> {
    sibship_order: &'a [usize],
    reconstructed_sibships: &'a [Vec<usize>],
    allow_hs_without_av: bool,
    effective_parents: &'a HashMap<usize, (String, String)>,
}

fn inference_lines(
    key: &str,
    members: &[usize],
    loaded: &Loaded,
    types: &unrelated::InfTypes,
    segments: &ibdseg::Segments,
    reconstructed: &InferencePedigree<'_>,
) -> String {
    let samples = &loaded.fileset.samples;
    let parents = |i: usize| {
        reconstructed
            .effective_parents
            .get(&i)
            .map(|p| (p.0.as_str(), p.1.as_str()))
            .unwrap_or((samples[i].pat.as_str(), samples[i].mat.as_str()))
    };
    let pedigree = Pedigree::from_samples(&with_phantom_parents(samples));
    let mut pedigree_cache = KinshipCache::default();
    let mut second_degree: HashMap<(usize, usize), bool> = HashMap::new();
    for (p, &a) in members.iter().enumerate() {
        for &b in &members[p + 1..] {
            let pair = if a < b { (a, b) } else { (b, a) };
            let is_second = if samples[a].fid == samples[b].fid {
                let phi = pedigree_kinship(&pedigree, &mut pedigree_cache, a, b);
                (band::SECOND..band::FIRST).contains(&phi)
            } else {
                types.inferred(loaded, a, b) == Some("2nd")
            };
            second_degree.insert(pair, is_second);
        }
    }
    let is_second = |a: usize, b: usize| {
        let pair = if a < b { (a, b) } else { (b, a) };
        second_degree.get(&pair).copied().unwrap_or(false)
    };
    let mut cache: HashMap<(usize, usize), Vec<(i64, i64)>> = HashMap::new();
    let intervals = |a: usize, b: usize, cache: &mut HashMap<(usize, usize), Vec<(i64, i64)>>| {
        let pair = if a < b { (a, b) } else { (b, a) };
        cache
            .entry(pair)
            .or_insert_with(|| segments.reported_intervals(loaded, pair.0, pair.1))
            .clone()
    };
    let ratio =
        |r: usize, a: usize, b: usize, cache: &mut HashMap<(usize, usize), Vec<(i64, i64)>>| {
            let ra = intervals(r, a, cache);
            let rb = intervals(r, b, cache);
            let ab = intervals(a, b, cache);
            join_ratio(&ra, &rb, &ab)
        };

    // Declared sibships in cluster order.  The first two entries are the pair the
    // inference pass tests; two-child fixtures pin every trigger independent of the
    // otherwise unstable ordering of larger sibships.
    let mut sibships: Vec<(&str, &str, Vec<usize>)> = Vec::new();
    for &i in reconstructed.sibship_order {
        if !members.contains(&i) {
            continue;
        }
        let s = &samples[i];
        if s.pat == "0" {
            continue;
        }
        match sibships
            .iter_mut()
            .find(|(pat, mat, _)| *pat == s.pat && *mat == s.mat)
        {
            Some((_, _, who)) => who.push(i),
            None => sibships.push((&s.pat, &s.mat, vec![i])),
        }
    }
    sibships.retain(|(_, _, who)| who.len() > 1);
    for group in reconstructed.reconstructed_sibships {
        let who = group.clone();
        if who.len() > 1
            && !sibships
                .iter()
                .any(|(_, _, existing)| existing.iter().all(|i| who.contains(i)))
        {
            sibships.push(("", "", who));
        }
    }

    let mut body = String::new();
    let mut rejected_before_first = 0usize;
    let mut printed = false;
    for (_, _, sibs) in &sibships {
        let (a, b) = (sibs[0], sibs[1]);
        for &r in members {
            if sibs.contains(&r)
                || !is_second(r, a)
                || !is_second(r, b)
                || segments.of(loaded, r, a).2 <= HS_CANDIDATE_PROP_IBD
                || segments.of(loaded, r, b).2 <= HS_CANDIDATE_PROP_IBD
            {
                continue;
            }
            let Some(value) = ratio(r, a, b, &mut cache) else {
                continue;
            };
            let Some(verdict) = av_verdict(value, samples[r].sex) else {
                if !printed {
                    rejected_before_first += 1;
                }
                continue;
            };
            if !printed {
                body.push_str(&"\n".repeat(1 + rejected_before_first));
                printed = true;
            }
            body.push_str(&console::build_inference_av_fs(
                key,
                &samples[r].iid,
                verdict,
                &samples[a].iid,
                &samples[b].iid,
                value,
            ));
        }
    }

    // Candidate half-sibs are cross-family pairs above the segment gate.  Their avuncular
    // checks precede the two parent-unrelated findings and the closing HS.UN2 line. A
    // PO-only family does not open this half unless an AV line already made it speak.
    if printed || reconstructed.allow_hs_without_av {
        for (p, &a) in members.iter().enumerate() {
            for &b in &members[p + 1..] {
                if samples[a].fid == samples[b].fid
                    || reconstructed
                        .reconstructed_sibships
                        .iter()
                        .any(|sibship| sibship.contains(&a) || sibship.contains(&b))
                    || parents(a).0 == "0"
                    || parents(b).0 == "0"
                    || !is_second(a, b)
                    || segments.of(loaded, a, b).2 <= HS_CANDIDATE_PROP_IBD
                {
                    continue;
                }
                let mut block = String::new();
                for &r in members {
                    if r == a || r == b || !is_second(r, a) || !is_second(r, b) {
                        continue;
                    }
                    let Some(value) = ratio(r, a, b, &mut cache) else {
                        continue;
                    };
                    if let Some(verdict) = av_verdict(value, samples[r].sex) {
                        block.push_str(&console::build_inference_av_hs(
                            key,
                            &samples[r].iid,
                            verdict,
                            &samples[a].iid,
                            &samples[b].iid,
                            value,
                        ));
                    }
                }
                let unrelated_parents = |child: usize, other: usize| {
                    let p = parents(child);
                    [p.0, p.1]
                        .into_iter()
                        .filter(|id| *id != "0")
                        .filter_map(|id| {
                            members
                                .iter()
                                .copied()
                                .find(|&i| samples[i].iid == id)
                                .filter(|&parent| {
                                    types.inferred(loaded, other, parent) == Some("UN")
                                })
                        })
                        .collect::<Vec<_>>()
                };
                let ua = unrelated_parents(a, b);
                let ub = unrelated_parents(b, a);
                if !ua.is_empty() && !ub.is_empty() {
                    for parent in &ua {
                        block.push_str(&console::build_hs_unrelated(
                            &samples[b].iid,
                            &samples[*parent].iid,
                        ));
                    }
                    for parent in &ub {
                        block.push_str(&console::build_hs_unrelated(
                            &samples[a].iid,
                            &samples[*parent].iid,
                        ));
                    }
                    block.push_str(&console::build_inference_hs_un2(
                        key,
                        &samples[a].iid,
                        &samples[b].iid,
                    ));
                }
                if !block.is_empty() {
                    if !printed {
                        body.push_str(&"\n".repeat(1 + rejected_before_first));
                        printed = true;
                    }
                    body.push_str(&block);
                }
            }
        }
    }
    if !printed && !sibships.is_empty() {
        body.push_str(&"\n".repeat(1 + rejected_before_first));
    }
    body
}

/// Of an inferred duplicate pair, `(the copy that stays, the copy that is removed)`.
///
/// The reference keeps the copy the `.fam` connects to more people — its **declared
/// 1st-degree relatives that the fileset actually carries**, counting named parents, full
/// sibs naming the same couple, and children naming it — and breaks a tie on the ID
/// comparator, keeping the later id.
///
/// Both halves are measured, out of sample, on ten shapes × three seeds
/// (`docs/research/fixtures/dupkeep.py`): 27 of 27 for this rule, against 21 of 27 for
/// "keep the later id" and 6 of 27 for "keep the earlier". The shapes that separate them
/// are the ones where the better-connected copy sorts *first* — an unparented singleton
/// against a child of a declared couple, two declared children one of which also has a
/// declared sib, and a lone copy against one that is itself a declared father.
fn duplicate_verdict(samples: &[Sample], i: usize, j: usize) -> (usize, usize) {
    let rank = |x: usize| declared_first_degree(samples, x);
    let (a, b) = match rank(i).cmp(&rank(j)) {
        std::cmp::Ordering::Greater => (i, j),
        std::cmp::Ordering::Less => (j, i),
        // Tie: the later id stays.
        std::cmp::Ordering::Equal => {
            if crate::analysis::king_id_cmp(samples[i].iid.as_bytes(), samples[j].iid.as_bytes())
                == std::cmp::Ordering::Greater
            {
                (i, j)
            } else {
                (j, i)
            }
        }
    };
    (a, b)
}

/// How many declared 1st-degree relatives of `i` the fileset actually carries.
fn declared_first_degree(samples: &[Sample], i: usize) -> usize {
    let me = &samples[i];
    let present = |id: &str| id != "0" && samples.iter().any(|s| s.iid == id);
    let mut n = usize::from(present(&me.pat)) + usize::from(present(&me.mat));
    if me.pat != "0" {
        n += samples
            .iter()
            .filter(|s| s.iid != me.iid && s.pat == me.pat && s.mat == me.mat)
            .count();
    }
    n + samples
        .iter()
        .filter(|s| s.pat == me.iid || s.mat == me.iid)
        .count()
}

/// One `RULE FS0`, `FS1`, or `FS2` line, with its people held as cluster positions.
///
/// `joiner` is `None` for `FS0` — the rule that *creates* a sibship and names its parents —
/// and `Some(who)` for `FS1`, where `members` is the sibship as it stood before the join.
/// `second` distinguishes `FS2` and holds its second declared sibship. `pat`/`mat` are
/// used only by `FS0`.
struct Rule {
    joiner: Option<usize>,
    members: Vec<usize>,
    second: Option<Vec<usize>>,
    pat: String,
    mat: String,
}

/// The parent pairs and the rule lines one cluster's reconstruction produces.
struct Assignment {
    /// `(position in the cluster, (father, mother))`.
    parents: Vec<(usize, (String, String))>,
    rules: Vec<Rule>,
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
/// A component holding **two** declared sibships raises `RULE FS2: Sibship (…) and
/// sibship (…) are combined`; the same libStatGen sibling order supplies both lists.
fn sibship_parents(
    samples: &[Sample],
    members: &[usize],
    is_fs: &dyn Fn(usize, usize) -> bool,
    sibship_order: &[usize],
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
        out.rules.extend(component_rules(
            samples,
            members,
            &group,
            &couple,
            sibship_order,
        ));
        for b in group {
            out.parents.push((b, couple.clone()));
        }
    }
    out
}

/// The `RULE FS0`/`FS1` lines one component raises, in the order the log prints them.
///
fn component_rules(
    samples: &[Sample],
    members: &[usize],
    group: &[usize],
    couple: &(String, String),
    sibship_order: &[usize],
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
    let rank = |position: &usize| {
        sibship_order
            .iter()
            .position(|&i| i == members[*position])
            .unwrap_or(usize::MAX)
    };
    for (_, _, who) in &mut declared {
        who.sort_by_key(rank);
    }
    if declared.len() > 1 {
        return vec![Rule {
            joiner: None,
            members: declared[0].2.clone(),
            second: Some(declared[1].2.clone()),
            pat: String::new(),
            mat: String::new(),
        }];
    }

    let mut rules = Vec::new();
    let mut sibship: Vec<usize> = match declared.first() {
        Some((_, _, who)) => who.clone(),
        None => {
            let opening: Vec<usize> = group.iter().take(2).copied().collect();
            rules.push(Rule {
                joiner: None,
                members: opening.clone(),
                second: None,
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
            second: None,
            pat: String::new(),
            mat: String::new(),
        });
        sibship.push(b);
    }
    rules
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

    /// The IBD1 fringe an IBD2 call cuts off is dropped when it is under the floor and
    /// kept when it is over — the clause that took `Join3/Join2` from 13 of 297 to 296.
    #[test]
    fn reported_intervals_cuts_ibd2_out_and_floors_each_piece() {
        // Markers every 1 Mb, so index and megabase agree.
        let pos: Vec<i64> = (0..40).map(|k| k * 1_000_000).collect();
        let ibd1 = [Called { lo: 0, hi: 30 }];
        let ibd2 = [Called { lo: 4, hi: 20 }];
        // Floor 3 Mb: the left piece [0, 3] is 3 Mb and stays, the right [21, 30] is 9.
        assert_eq!(
            reported_intervals(&pos, &ibd1, &ibd2, 3_000_000),
            vec![
                (0, 3_000_000),
                (4_000_000, 20_000_000),
                (21_000_000, 30_000_000)
            ]
        );
        // Floor 5 Mb: the left piece is now under it and is not reported at all, while
        // the IBD2 call it was cut by is untouched.
        assert_eq!(
            reported_intervals(&pos, &ibd1, &ibd2, 5_000_000),
            vec![(4_000_000, 20_000_000), (21_000_000, 30_000_000)]
        );
        // With no IBD2 call the whole IBD1 call is one piece.
        assert_eq!(
            reported_intervals(&pos, &ibd1, &[], 3_000_000),
            vec![(0, 30_000_000)]
        );
    }

    #[test]
    fn join_ratio_is_the_triple_over_the_double_intersection() {
        let a = [(0, 100), (200, 300)];
        let b = [(50, 250)];
        // Join2 = [50,100] + [200,250] = 100.
        let s = [(60, 90)];
        assert_eq!(join_ratio(&a, &b, &s), Some(0.3));
        // Disjoint from R: no denominator, so no line.
        assert_eq!(join_ratio(&a, &[(400, 500)], &s), None);
    }

    /// The three branches, including the one the old single-cut reading did not have.
    #[test]
    fn av_verdict_has_a_dead_band() {
        assert_eq!(av_verdict(0.8487, 1), Some("uncle"));
        assert_eq!(av_verdict(0.8487, 2), Some("aunt"));
        assert_eq!(av_verdict(0.8512, 1), None);
        assert_eq!(av_verdict(0.8969, 1), None);
        assert_eq!(av_verdict(0.9001, 1), Some("grandfather, HS, or nephew"));
        assert_eq!(av_verdict(0.9001, 2), Some("grandmother, HS, or niece"));
        // Both edges are measured as closed on the silent side.
        assert_eq!(av_verdict(AV_UNCLE_MAX, 1), None);
        assert_eq!(av_verdict(AV_AMBIGUOUS_MIN, 1), None);
    }

    /// `bigish`'s two `KING1` cousin pairs, which is why it prints one `HS.UN2` not two.
    #[test]
    fn hs_candidate_gate_splits_bigish_two_cousin_pairs() {
        let candidate = |prop_ibd: f64| prop_ibd > HS_CANDIDATE_PROP_IBD;
        assert!(candidate(0.1901), "B01_C3/B02_C4 raises HS.UN2");
        assert!(!candidate(0.1799), "B01_C2/B02_C2 does not");
        // The bracket the 26 measured pairs leave the constant in.
        assert!(candidate(0.1878) && !candidate(0.1868));
    }

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
        let got = sibship_parents(samples, &members, &is_fs, &members, next);
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

    /// Of a duplicate pair, the copy the `.fam` connects to more people stays — and that
    /// beats "the later id", which is the rule the first shape to show a duplicate
    /// suggested. `dupkeep.py`'s `lone_first` in miniature: `A1` sorts first and declares
    /// nothing, `Z_C1` sorts last and has two parents and a sib.
    #[test]
    fn the_better_connected_copy_of_a_duplicate_stays() {
        let samples = vec![
            person("FA1", "A1", "0", "0"),
            person("FZ", "Z_C1", "Z_F", "Z_M"),
            person("FZ", "Z_C2", "Z_F", "Z_M"),
            person("FZ", "Z_F", "0", "0"),
            person("FZ", "Z_M", "0", "0"),
        ];
        assert_eq!(declared_first_degree(&samples, 0), 0);
        assert_eq!(declared_first_degree(&samples, 1), 3);
        assert_eq!(duplicate_verdict(&samples, 0, 1), (1, 0));
        assert_eq!(duplicate_verdict(&samples, 1, 0), (1, 0));
    }

    /// With nothing to separate them the later id stays: three singleton full sibs and a
    /// duplicate of the first removes `Q1` and keeps `Q4`.
    #[test]
    fn a_tied_duplicate_keeps_the_later_id() {
        let samples = vec![person("FQ1", "Q1", "0", "0"), person("FQ4", "Q4", "0", "0")];
        assert_eq!(duplicate_verdict(&samples, 0, 1), (1, 0));
    }

    /// Children count towards the connectivity, not just parents and sibs.
    #[test]
    fn a_declared_parent_outranks_a_lone_copy() {
        let samples = vec![
            person("FD", "D_C1", "D_F", "D_M"),
            person("FD", "D_F", "0", "0"),
            person("FD", "D_M", "0", "0"),
            person("FZL", "ZL1", "0", "0"),
        ];
        assert_eq!(declared_first_degree(&samples, 1), 1);
        assert_eq!(declared_first_degree(&samples, 3), 0);
        assert_eq!(duplicate_verdict(&samples, 1, 3), (1, 3));
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
        let got = sibship_parents(samples, &members, &is_fs, &members, next);
        let name = |n: &usize| samples[members[*n]].iid.as_str();
        got.rules
            .iter()
            .map(|r| {
                let names: Vec<&str> = r.members.iter().map(name).collect();
                match (&r.second, r.joiner) {
                    (Some(second), _) => {
                        let second_names: Vec<&str> = second.iter().map(name).collect();
                        console::build_rule_fs2("KING1", &names, &second_names)
                    }
                    (None, None) => console::build_rule_fs0("KING1", &names, &r.pat, &r.mat),
                    (None, Some(j)) => console::build_rule_fs1("KING1", name(&j), &names),
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

    /// Two declared sibships tied by an inferred full-sib pair raise the one `FS2` line.
    #[test]
    fn two_declared_sibships_in_one_component_raise_fs2() {
        let samples = two_families();
        let mut next = 1;
        // A_C1 (0) and B_C1 (4) are inferred full sibs, so both declared sibships merge.
        let got = rule_lines(&samples, &[(0, 4)], &mut next);
        assert_eq!(
            got,
            vec![
                "  Family KING1 RULE FS2: Sibship (A_C1 A_C2) and sibship (B_C1 B_C2) are combined\n"
            ]
        );
        let rows = assign(&samples, &[(0, 4)], &mut 1);
        assert_eq!(rows[4], ("B_C1".into(), "A_F".into(), "A_M".into()));
    }

    /// libStatGen's equal-key partition residue is intentionally unstable and differs
    /// from either Rust sort. The end-to-end `bigish` parity case pins the global pedigree
    /// context; this unit sentinel pins the low-level five-element permutation.
    #[test]
    fn libstatgen_equal_parent_order_matches_the_reference() {
        let ids = ["A_C1", "A_C2", "A_C3", "A_C4", "A_C5"];
        let mut order: Vec<usize> = (0..ids.len()).collect();
        unstable_king_sort(&mut order, |a, b| {
            crate::analysis::king_id_cmp(ids[a].as_bytes(), ids[b].as_bytes())
        });
        unstable_king_sort(&mut order, |_, _| std::cmp::Ordering::Equal);
        assert_eq!(
            order.into_iter().map(|i| ids[i]).collect::<Vec<_>>(),
            ["A_C3", "A_C4", "A_C5", "A_C2", "A_C1"]
        );
    }

    /// Nobody joined, nobody reassigned: a cluster merged by something other than a full
    /// sib pair leaves every parent column exactly as the `.fam` had it.
    #[test]
    fn a_cluster_with_no_inferred_sibship_assigns_nothing() {
        let samples = two_families();
        let members: Vec<usize> = (0..samples.len()).collect();
        let mut next = 1;
        let got = sibship_parents(&samples, &members, &|_, _| false, &members, &mut next);
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
