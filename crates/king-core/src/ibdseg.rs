//! IBD-segment calling — the engine behind `--ibdseg` and the segment columns of
//! `--related`.
//!
//! # Provenance
//!
//! KING's IBD-segment algorithm is **unpublished**: the manual says the manuscript is
//! "yet to be published", and the citation the binary prints ("Chen et al. 2024") does
//! not exist in any index. Nothing here comes from KING's source. Every rule below was
//! established by running the reference binary on filesets whose genotypes were
//! constructed so that the answer is forced, and the experiment that fixes each rule is
//! named where the rule is stated. `docs/BEHAVIOR.md` carries the raw sweeps.
//!
//! # The one structural fact everything rests on
//!
//! KING scans genotypes in **64-marker words of the global marker array** — the same
//! words the loader reports as "Autosome genotypes stored in N words". Word `w` covers
//! marker indices `64w ..= 64w+63` counting from the first retained autosomal marker of
//! the whole `.bim`, **not** from the start of a chromosome. That is why a chromosome's
//! fate can depend on how many markers precede it: see [`usable_segments`].
//!
//! # What is verified and what is not
//!
//! * [`usable_segments`] — **verified**. Reproduces `<prefix>allsegs.txt` byte for byte
//!   on all ten corpus datasets that emit one, including the two cases (`dups`
//!   chr14/chr15, `threegen` chr21/chr22) where a *longer* chromosome is dropped and a
//!   shorter one kept, which no per-chromosome rule can explain.
//! * The **boundary convention** — verified. Segments live on the word grid; an IBD1 call
//!   is then refined at both ends by the flanking word's *last* IBS0, an IBD2 call only at
//!   its right end. Inverting the `MaxIBD2` column of `--ibs` locates 154 of 158 corpus
//!   segments to exactly one word interval `[64u, 64v+63]`, every one of them aligned, and
//!   the convention is what makes a parent–offspring pair read exactly
//!   `1.0000 / 0.0000 / 0.5000`. Confirmed by forced-IBS0 sweeps on constructed filesets.
//! * [`inf_type`] — verified. Reproduces all 8 722 `InfType` values in the captured
//!   corpus from each row's own printed columns.
//! * [`Scan::ibd1`]'s word rule (no IBS0 tolerance at all) — verified by forced-IBS0
//!   sweeps.
//! * [`MIN_INFORMATIVE`] — **measured** on hand-written-genotype fixtures where the count
//!   is exact, then validated on data that had no part in choosing it: the corpus
//!   separates on it with no overlap at all, and a 512-invocation word-grid sweep agrees
//!   511 times. It is what [`MIN_RUN1`] used to stand in for.
//! * [`Scan::ibd2`], [`IBD2_HET_DIRTY`] and [`IBD2_REACH`] — **measured on a `.seg`-native
//!   canvas** (`docs/research/17-seg-caller.md`), which is the instrument the rest of this
//!   header used to ask for: every constant below is a bisection read off the reference
//!   through a marker ruler that recovers the number of calls and the number of words
//!   exactly, not an inversion of some other pass's column. One clause — the bridging
//!   condition — is fitted rather than bisected and is flagged in place. See the section
//!   below for what is left.
//! * [`reported_at_degree`] — **measured**, over 38 298 differential cases plus a
//!   constructed fixture for the clause the corpus cannot reach.
//! * [`Scan::ibd2_words`] — **solved for the corpus.** The `--ibs` IBD2 caller is a
//!   quantised confirmation scan (`docs/research/16-segment-extension.md`), and both of
//!   `--ibs`'s IBD2 columns — `MaxIBD2` and `Pr_IBD2` — are exact on **all 21 560**
//!   `.ibs`/`.ibs0` corpus rows and on all **158** pairs the reference grades. Out of
//!   sample on constructed random word canvases it is right about 93 % of the time, which
//!   is stated at [`Scan::ibd2_words`] rather than rounded up to "done".
//!
//! * [`Scan::ibd1`] and [`ibd1_pieces`] — **measured on an IBD1-native canvas**
//!   (`docs/research/18-ibd1-caller.md`), the mirror of the one above: the painted region
//!   is held IBD2-free so the printed `IBD1Seg` reads back marker intervals directly. Its
//!   words, its gate, the absence of a push and the absence of bridging are all bisections,
//!   and the `IBD1Seg` column is now exact on **all 982** corpus rows at the default floor.
//!
//! * [`Scan::join_runs`], [`Scan::join_runs2`] and [`Scan::merge_ok`] — the
//!   **`--seglength` run merge** (`docs/research/20-seglength-floor.md`, corrected by
//!   `docs/research/21-push-merge.md`), the clause `docs/research/18-ibd1-caller.md` §9
//!   measured on two of its five conditions and deliberately left out because that version
//!   made the corpus much worse. Both passes join two runs across a short interruption
//!   when a budget `cost * (bad - 2) <= informative` passes over the interrupting words,
//!   and the gate is asked *first*, so a refused run lies inside an interruption instead
//!   of ending one. **The two passes are otherwise not the same rule.** IBD1 joins across
//!   at most two unusable words and measures the gap run-to-run; IBD2 has no cap at all,
//!   measures everything between the two runs' **gate windows** rather than between the
//!   runs, and its informative count is HetHet with a switch to A1A1/A1A1 below
//!   [`MIN_INFORMATIVE`]. The conditioned merged calls also feed the ">10 Mb" pair
//!   filter; a held-out IBD1 pair and an independent IBD2 canvas require that ordering.
//!   Neither can fire at the default floor on real spacings; at 5 and 10 Mb they take
//!   `IBD1Seg` from 910/844 to **982/970** and `IBD2Seg` from 946/937 to **982/972**,
//!   worth **11 parity cases** between them.
//!
//! * The **one-word push** — every call after the first in a usable segment starting one
//!   word later (`docs/research/17-seg-caller.md` §6) — is **conditional**
//!   (`docs/research/21-push-merge.md` §2): a call arms it only when it reaches half the
//!   floor, measured from its own gate-start word. At the default floor that is almost
//!   always true, which is why §6 read it as unconditional.
//!
//! * The **fringe** — the partial word beyond a usable segment's word grid, read by both
//!   passes — is **measured** (`docs/research/19-ibd2seg-residual.md`) on a third canvas,
//!   `docs/research/fixtures/fringecanvas.py`, which builds a segment that does *not* start
//!   on a word boundary by shortening the carrier chromosome by `f` markers. The other two
//!   rigs are word-aligned and cannot reach any of it, so this clause was a corpus fit
//!   until that rig existed. Both stops are bisected at 16 positions a side, and the
//!   asymmetry is real and would not have been guessed: **an opposite homozygote in a
//!   fringe does not stop an IBD2 call**, though one anywhere in a complete word inside the
//!   grid disqualifies that word outright. Each pass stops at its own breaking marker — an
//!   opposite homozygote for IBD1, a het-vs-hom mismatch for IBD2. See `fringe_tests`.
//!
//! # What is still not right — and none of it is at the default floor
//!
//! Against the captured reference `.seg` files at the default 3 Mb floor the caller
//! reproduces **all 982 rows** with all four printed columns identical: `IBD1Seg` 982,
//! `IBD2Seg` 982, `PropIBD` 982, `InfType` 982, mean and worst absolute error 0.0000, and
//! the reported pair set exactly right (0 extra, 0 missing) on all ten datasets. There is
//! no residual left to point at here.
//!
//! Two of those four columns were closed by a *writer* rule rather than by this function,
//! and the distinction matters for anyone continuing. `<prefix>.seg` computes `PropIBD`
//! from the four decimals it is about to print rather than from the totals
//! ([`seg_prop_ibd`]), and lists its rows in 16-sample blocks rather than by sample index
//! (`analysis::ibdseg::seg_pair_order`). Neither touches a segment. Before they landed,
//! this caller scored 806 of 982 exact rows with both estimate columns already at 982, and
//! `docs/research/19-ibd2seg-residual.md` §9 concluded from those 176 rows that the IBD1
//! pass was "systematically about a marker short". It was not — the diagnosis was an
//! artefact of comparing two different formulas, and `20-seg-writer.md` is the correction.
//! **Do not re-derive the sub-ulp IBD1 argument; it is dead.**
//!
//! **All three captured floors are now exact.** `--seglength` 3, 5 and 10 each reproduce
//! **982 of 982** rows on all four printed columns, MAE 0.000000, worst row 0.0000, 0
//! extra and 0 missing pairs (`tests/parity/fit/scorecard.py`). The last two `.seg` parity
//! cases closed with `docs/research/23-gap-bound.md`, which found that the 10 Mb residual
//! was neither the merge's gap nor an invented merge — the diagnosis `…/21-…` §8.1 handed
//! on, and which `chrprobe.py` refuted by reading the reference one chromosome at a time:
//!
//! * **The floor is asked twice.** A run is reported only if its **gate window** spans
//!   [`WINDOW_FRACTION`]-th of `--seglength`, whatever the reported call measures. On the
//!   corpus an 11.2066 Mb IBD2 call is reported at 6 290 751 bp and gone at 6 290 752,
//!   which is twice its 3.1 Mb one-word window. Both passes have it, with different
//!   strictness, and it is asked at emit — a run the bound refuses still merges.
//! * **The IBD1 merge's budget is summed over every word between the two runs**, a
//!   gate-refused run included, while [`MERGE_MAX_WORDS`] still counts only the unusable
//!   ones. That is `20-…` §11 item 4, left undecided there and bisected in `23-…` §5.
//!
//! **What is still wrong is out of sample, and measured.** `docs/PARITY.md` §4.6:
//! `docs/research/fixtures/oosseg.py` runs 24 fresh filesets on unused seeds through both
//! binaries at all three floors and gets **68 of 72** byte-identical — 4 value-differing
//! rows in 6 713, with 0 extra and 0 missing. All four occur at exactly 40 000 markers;
//! 39 999 and 40 001-marker controls are exact. This is KING's measured uninitialised
//! exact-multiple-of-64 tail read, intentionally not reproduced by safe Rust. The formerly
//! missing distant pair is fixed by letting the measured merged calls feed the filter.
//!
//! **How to grade work on it.** Not with the `.seg` row counts, which are now saturated at
//! all three floors, and not with `--ibs`, whose IBD2 columns have been exact under every
//! candidate since [`Scan::ibd2_words`]. Grade out of sample — `oosseg.py` for the whole
//! program, the canvases (`docs/research/fixtures/window1.py` §7, `mergelab.py`,
//! `push1.py`) for one clause at a time — and, when a real row does go wrong, localise it
//! with `chrprobe.py` before theorising: every campaign before `23-…` had to guess which
//! segment of which pair a wrong row came from, and the two guesses `21-…` §8.1 recorded
//! were both wrong. Do **not** re-sweep this function's
//! constants: forty single-knob perturbations and all 32 combinations of the two IBD1
//! endpoint rules crossed with the two IBD1 fringe rules were scored in the final pass, and
//! the committed values are the unique maximum on both exact rows and mean error
//! (`20-seg-writer.md` §6).
//!
//! One further open item: the **100 Mb usable-total floor** (`…/17-…` §2) — the reference
//! refuses a fileset whose usable total is under 100 000 000 bp, bisected to the base pair,
//! and this crate does not model it. No corpus dataset is anywhere near it.
//!
//! **What is *not* the answer, measured** (`…/16-segment-extension.md` §9,
//! `…/17-seg-caller.md` §8): porting [`Scan::ibd2_words`]'s chunk geometry to this pass.
//! `.seg` is **not a quantised confirmation scan** — no chunk quantum, no confirmation
//! count, HetHet and A1A1/A1A1 interchangeable where `--ibs` ignores the latter, and
//! nothing ever cut. The port scored 709 exact rows against 705 but nearly tripled mean
//! `PropIBD` error (0.00138 → 0.00356); `tests/parity/fit/segtry.py` keeps that measurement
//! reproducible so nobody re-runs the experiment blind.
//!
//! Read the per-dataset split with the caveat that four of the ten filesets report only
//! the 14 within-family pairs of one six-person nuclear family, over 5 000 to 10 000
//! markers — and in `monomorphic`'s case half of those markers are monomorphic or
//! ultra-rare. The reference's own numbers there are nowhere near the pedigree truth, so
//! they grade nothing. `bigish` 582/763, `multifam` 77/104, `threegen` 30/39, `admixed`
//! 13/16, `monomorphic` 13/14, `dups` 2/3, `unrelated` 1/1; `nuclear` 9/14, `missing`
//! 10/14, `sexchr` 10/14. `docs/PARITY.md` §5 carries the evidence.

use king_io::Genotypes;

/// Markers per scan word. The whole engine is quantised to this.
pub const WORD: usize = 64;

/// Complete words a run must span before it can become a reported IBD1 segment.
///
/// **One, not two.** A two-word floor is what the corpus wants while
/// [`MIN_INFORMATIVE`] is missing — raising it from one to two took the `.seg` rows
/// agreeing at all four printed decimals from 315/982 to 626/982 — but it was the
/// informativeness gate in disguise: a lone clean word usually falls short of ten
/// informative markers, and when it does not, the reference calls it. A deterministic
/// one-word fixture settles it directly (`docs/research/13-informativeness-gate.md` §6):
/// at 9 informative markers the reference reports nothing, at 10 it reports the word's
/// full 127 marker intervals.
const MIN_RUN1: usize = 1;

/// The same floor for IBD2 runs.
///
/// Measured, not fitted: on a constructed fixture a single IBD2-clean word is reported as
/// a segment of exactly 63 marker intervals (`docs/research/10-segment-rule-fixtures.md`
/// §3), and `MaxIBD2` agrees. The IBD1 floor reached the same value later and by a
/// different route, so the two constants stay separate: they were established by
/// different experiments and nothing says they must move together.
const MIN_RUN2: usize = 1;

/// Unusable words an **IBD1** `--seglength` merge may bridge — `20-seglength-floor.md`.
///
/// Bisected on `mergelab.py` (§3): a one- and a two-word interruption are joined at every
/// floor above their own gap, a three-word one at none, however little the interrupting
/// words carry and however short the gap. Runs the gate refused do **not** count here —
/// they lie inside the interruption and are stepped over (§6).
///
/// The IBD2 pass has **no** such cap — `21-push-merge.md` §4 re-measures it on a fixture
/// built for the purpose and joins fifteen unusable words. The two passes really do
/// differ: the same fixture on the IBD1 pass still refuses three.
const MERGE_MAX_WORDS: usize = 2;

/// The fraction of `--seglength` a call must reach to arm the one-word push.
///
/// `21-push-merge.md` §2, bisected to the base pair on three spacings: a call arms the
/// push iff `pos[hi] - pos[64 * gs] >= seglength / PUSH_FRACTION`, where `gs` is its
/// gate-start word. The integer division is the reference's own: at a floor of 5 080 001
/// bp a 2 540 000 bp call still arms it, at 5 080 100 it does not.
const PUSH_FRACTION: i64 = 2;

/// The fraction of `--seglength` a run's **gate window** must span to be reported.
///
/// `23-gap-bound.md`: the floor is asked twice, and the second question is not about
/// the reported call at all. A call whose window is one word is dropped once the floor
/// passes twice that window's span, however long the call itself measures — on the corpus
/// an 11.2066 Mb IBD2 call is reported at `--seglength 6.290751` and gone at 6.290752,
/// which is `2 * 3 145 375 + 1`, its window being one 3.1 Mb word. Bisected to the base
/// pair on two independent corpus calls (§1) and on purpose-built canvases at four
/// spacings (§2), with the same integer division [`PUSH_FRACTION`] uses.
///
/// The two passes ask it with different strictness, which is bisected and not a guess:
/// IBD2 keeps a window of exactly `min_bp / 2` and IBD1 does not (§4). Equivalently the
/// IBD1 span is measured one base pair shorter; no fixture can separate the two readings.
const WINDOW_FRACTION: i64 = 2;

/// Bad markers a merge gets for nothing, and what each further one costs.
///
/// The merge test is `cost * (bad - MERGE_FREE) <= informative`, summed over the
/// interrupting words. `MERGE_FREE` is bisected from the all-A2A2 interruption, which
/// carries no informative marker at all and joins at two opposite homozygotes but not at
/// three, on **both** passes; the costs are bisected off the boundary line of the
/// `(bad, informative)` grid — ten values of `bad` on the IBD1 pass (§4) and the
/// HetHet/A1A1 sweep on the IBD2 one (§7).
const MERGE_FREE: u32 = 2;
/// IBD1: one further opposite homozygote costs four informative markers.
const MERGE_COST1: u32 = 4;
/// IBD2: one further bad marker costs three. The IBD2 pass counts a het-vs-hom mismatch
/// as bad alongside an opposite homozygote; the IBD1 pass does not count it at all.
const MERGE_COST2: u32 = 3;

/// Informative markers a run must carry over its **own complete words** to be called.
///
/// The absence of a contradiction is only evidence where a contradiction had the chance
/// to appear, and this is that test: a run `[u..v]` is reported only if at least ten of
/// the markers `64u ..= 64(v+1)-1` are informative for the pair, in the sense of
/// `WordDiff::inf1` (IBD1) or `WordDiff::inf2` (IBD2). Failing runs are dropped
/// outright — not shortened, not merged, not re-scored.
///
/// **Measured, then validated out of sample** (`docs/research/13-informativeness-gate.md`).
/// The constant comes from hand-written-genotype fixtures where the count is exact: ten
/// passes and nine fails at every run width from 1 to 14 words and for three different
/// placements of the informative markers inside the run. The corpus then separates on it
/// without having chosen it — over 1 170 pairs every one the reference refuses has at
/// most 9 and every one it reports has at least 10, with 62 refusals sitting at exactly 9
/// and 60 acceptances at exactly 10, so 9 costs 62 extra pairs and 11 costs 60 missing
/// ones. A 512-invocation word-grid sweep (shifting the grid under fixed genotypes)
/// agrees with the reference 511 times with no false accepts.
const MIN_INFORMATIVE: u32 = 10;

/// Het-vs-hom disagreements that make a word too dirty to sit inside a `.seg` IBD2 run.
///
/// **Two: one mismatch leaves a word usable, two do not.** Measured directly on the
/// `.seg`-native canvas of `docs/research/17-seg-caller.md` §3 — a block of eight
/// pure-HetHet words with `j` consecutive words replaced by one carrying `m` mismatches,
/// read back through a marker ruler that recovers the number of words *and* the number of
/// calls. `m = 1` never splits the block at any `j`; `m = 2` always does. The corpus agrees
/// and its gradient is sharp: `IBD2Seg` is exact on 896 of 982 rows at 2, against
/// 867 / 868 / 837 at 1 / 3 / 5.
///
/// The **5** this used to be was inverted out of `--ibs`'s `MaxIBD2`, which is a different
/// caller's column (`docs/research/15-ibs-ibd2-rules.md` §3) — see [`IBS_IBD2_DIRTY`],
/// which is that caller's threshold and is still 5.
const IBD2_HET_DIRTY: u32 = 2;

/// Markers a `.seg` IBD2 call reaches past the nearest het-vs-hom mismatch bounding it.
///
/// The endpoints are **not** word-aligned. Swept over all 64 bit positions of the flanking
/// word (`docs/research/17-seg-caller.md` §5): the extension is `127 − lastbit` to the left
/// and `64 + firstbit` to the right, which is the same statement twice — the call runs 63
/// markers past the nearest mismatch on that side. Bit 31 gives exactly +96, so 62 (→95)
/// and 64 (→97) are both excluded. Corpus: `IBD2Seg` exact on 896 rows at 63, against
/// 824 / 825 / 830 / 826 at 0 / 32 / 62 / 64.
const IBD2_REACH: usize = 63;

/// Whether scan word `k` is too dirty to sit inside an IBD2 run.
///
/// An opposite homozygote is disqualifying on its own — the two samples cannot be IBD2
/// where they share no allele — and het-vs-hom disagreements disqualify from
/// [`IBD2_HET_DIRTY`] up.
fn ibd2_dirty(scan: &Scan, k: usize) -> bool {
    scan.ibs0_at(k) != 0 || scan.ibs1[k].count_ones() >= IBD2_HET_DIRTY
}

/// A usable segment is cut wherever two consecutive markers are further apart than this.
///
/// Verified to the byte: a gap of exactly 1 000 000 does **not** cut, 1 000 001 does.
pub const MAX_MARKER_GAP: i64 = 1_000_000;

/// ...and also wherever one whole scan word spans more than this.
///
/// Measured by sweeping uniform marker spacing `s`: `s = 156_250` (so `64s` is exactly
/// 10 000 000) leaves the chromosome whole; `s = 156_251` shatters it into single words.
/// The span compared is `pos[64(w+1)] - pos[64w]`, i.e. **64** gaps, not 63 — spacing
/// 157 000 cuts even though `63s` is still under 10 Mb.
pub const MAX_WORD_SPAN: i64 = 10_000_000;

/// A usable segment must contain at least this many complete words...
pub const MIN_WORDS: usize = 5;

/// ...and its word-aligned span must exceed this.
///
/// "Word-aligned" is the trap: the length tested is
/// `pos[last marker of the last complete word] - pos[first marker of the first complete
/// word]`, not the piece's full span. A 2000-marker chromosome spanning 10.075 Mb is
/// dropped (word-aligned 9.994 Mb) while one spanning 10.081 Mb is kept (10.000 Mb).
pub const MIN_SEGMENT_BP: i64 = 10_000_000;

/// A pair is reported only if it has at least one segment longer than this.
///
/// The reference announces it verbatim on every run: "Sample pairs without any long IBD
/// segments (>10Mb) are excluded." Not tunable — `--seglength` does not move it.
pub const LONG_SEGMENT_BP: i64 = 10_000_000;

/// Default `--seglength`, in base pairs.
pub const DEFAULT_SEGLENGTH_BP: i64 = 3_000_000;

/// One usable chromosomal segment: the denominator's unit, one row of `allsegs.txt`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Usable {
    /// Chromosome code, as printed.
    pub chr: i64,
    /// First marker index, in the array this segment was cut from.
    pub lo: usize,
    /// Last marker index, inclusive.
    pub hi: usize,
}

impl Usable {
    /// Markers in the segment.
    pub fn n_snp(self) -> usize {
        self.hi - self.lo + 1
    }

    /// First complete word of the global grid lying entirely inside the segment.
    pub fn first_word(self) -> usize {
        self.lo.div_ceil(WORD)
    }

    /// Last complete word of the global grid lying entirely inside the segment.
    ///
    /// Only meaningful when [`Usable::words`] is non-zero.
    pub fn last_word(self) -> usize {
        (self.hi + 1) / WORD - 1
    }

    /// How many complete words the segment contains.
    pub fn words(self) -> usize {
        ((self.hi + 1) / WORD).saturating_sub(self.lo.div_ceil(WORD))
    }
}

/// Cut a marker array into the segments KING considers usable.
///
/// `chr` and `pos` are parallel arrays over **one** analysis array — the retained
/// autosomal markers, or the X markers — in `.bim` order. Indices are the global word
/// grid for that array.
///
/// Three cuts and two filters, in this order:
///
/// 1. cut at every chromosome change and every marker gap over [`MAX_MARKER_GAP`];
/// 2. inside a piece, cut between complete words `w` and `w+1` whenever
///    `pos[64(w+1)] - pos[64w]` exceeds [`MAX_WORD_SPAN`] — the sub-piece boundary lands
///    on the word boundary, so the left part ends at marker `64(w+1)-1`;
/// 3. keep a sub-piece only if it holds at least [`MIN_WORDS`] complete words **and** its
///    word-aligned span exceeds [`MIN_SEGMENT_BP`].
///
/// The word-count filter is what makes this dataset-global rather than per-chromosome: on
/// the `dups` corpus fileset chromosome 14 (372 markers, 18.6 Mb) is dropped while
/// chromosome 15 (355 markers, 17.7 Mb) is kept, purely because chromosome 14 starts at
/// global index 7622 and so straddles the grid badly enough to contain only four complete
/// words, while chromosome 15 contains five.
pub fn usable_segments(chr: &[i64], pos: &[i64]) -> Vec<Usable> {
    assert_eq!(chr.len(), pos.len());
    let mut out = Vec::new();
    let n = chr.len();
    if n == 0 {
        return out;
    }
    let mut start = 0usize;
    for i in 1..=n {
        let cut = i == n || chr[i] != chr[i - 1] || pos[i] - pos[i - 1] > MAX_MARKER_GAP;
        if cut {
            split_by_word_span(chr, pos, start, i - 1, &mut out);
            start = i;
        }
    }
    out
}

/// Apply cut (2) and filter (3) of [`usable_segments`] to one gap-free piece.
fn split_by_word_span(chr: &[i64], pos: &[i64], lo: usize, hi: usize, out: &mut Vec<Usable>) {
    let piece = Usable {
        chr: chr[lo],
        lo,
        hi,
    };
    let mut sub_lo = lo;
    if piece.words() >= 2 {
        for w in piece.first_word()..piece.last_word() {
            if pos[WORD * (w + 1)] - pos[WORD * w] > MAX_WORD_SPAN {
                push_if_usable(chr, pos, sub_lo, WORD * (w + 1) - 1, out);
                sub_lo = WORD * (w + 1);
            }
        }
    }
    push_if_usable(chr, pos, sub_lo, hi, out);
}

/// Filter (3): keep a sub-piece only if it is both wide enough and word-rich enough.
fn push_if_usable(chr: &[i64], pos: &[i64], lo: usize, hi: usize, out: &mut Vec<Usable>) {
    let seg = Usable {
        chr: chr[lo],
        lo,
        hi,
    };
    if seg.words() < MIN_WORDS {
        return;
    }
    let aligned = pos[WORD * seg.last_word() + WORD - 1] - pos[WORD * seg.first_word()];
    if aligned > MIN_SEGMENT_BP {
        out.push(seg);
    }
}

/// Total length in base pairs of a set of usable segments — the IBD denominator.
pub fn denominator(segs: &[Usable], pos: &[i64]) -> i64 {
    segs.iter().map(|s| pos[s.hi] - pos[s.lo]).sum()
}

// ---------------------------------------------------------------------------
// Per-pair segment calling
// ---------------------------------------------------------------------------

/// A called segment, as a closed marker-index range.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Called {
    pub lo: usize,
    pub hi: usize,
}

/// The two genotype disagreements the scan counts, evaluated one word at a time.
///
/// Bit `k` of the returned masks is marker `64w + k`.
///
/// * `ibs0` — opposite homozygotes. With the loader's plane encoding (`plane0` = "is
///   homozygous", `plane1` = "carries A1" for a homozygote) that is
///   `hom_i & hom_j & (a1_i ^ a1_j)`: a heterozygote or a missing call has `plane0 = 0`
///   and so can never contribute.
/// * `ibs1` — one heterozygote against one homozygote, both called.
/// * `inf1` — markers that count towards [`MIN_INFORMATIVE`] for an IBD1 run: both
///   samples carry the A1 allele and at least one is homozygous for it. Equivalently,
///   inside a word with no IBS0, "at least one of the pair is A1A1 and both are called" —
///   exactly the markers at which an IBS0 *could* have been seen. The two readings differ
///   only where an IBS0 is present, which disqualifies the word anyway, so no experiment
///   inside a run can tell them apart.
/// * `inf2` — the same count for an IBD2 run: both samples carry A1 **and neither is a
///   heterozygote against the other's A1A1**, i.e. HetHet plus A1A1/A1A1. HetHet is worth
///   1 to an IBD2 run and 0 to an IBD1 one, verified in both directions on fixtures
///   (`docs/research/13-informativeness-gate.md` §5); a pair A2A2/A2A2 is worth 0 to
///   either. The `& !ibs1` clause is the correction of `…/17-seg-caller.md` §14.3: twenty
///   words carrying one het-vs-A1A1 marker each are refused by the `.seg` gate where ten
///   HetHet or ten A1A1/A1A1 markers pass it, so a het-vs-A1A1 marker — which *is*
///   `p1 & p1` — is **not** informative. It is the only mask that changed, and only
///   [`Scan::ibd2`] reads it; the `--ibs` pass counts `hethet`.
///
/// `A1` is the `.bim`'s **first allele column**, taken literally — not the minor allele
/// and not any cohort frequency. A pair homozygous for the minor allele counts zero when
/// that allele sits in the A2 column, and sliding the rest of the cohort's genotypes
/// (hence the marker's MAF) from 2/12 to 6/12 with the pair held fixed does not move a
/// single call. KING's insistence that A1 be the minor allele is what makes reading the
/// column behave like a frequency filter.
#[derive(Clone, Copy, Debug, Default)]
struct WordDiff {
    ibs0: u64,
    ibs1: u64,
    inf1: u64,
    inf2: u64,
    hethet: u64,
}

fn word_diff(g: &Genotypes, i: usize, j: usize, w: usize) -> WordDiff {
    let (p0i, p1i) = (g.plane0[i][w], g.plane1[i][w]);
    let (p0j, p1j) = (g.plane0[j][w], g.plane1[j][w]);
    let het_i = !p0i & p1i;
    let het_j = !p0j & p1j;
    // `plane1` is "carries A1" — set for A1A1 and for a heterozygote, clear for a
    // missing call — so `share` already excludes missing genotypes on both sides.
    let share = p1i & p1j;
    let ibs1 = (het_i & p0j) | (p0i & het_j);
    WordDiff {
        ibs0: p0i & p0j & (p1i ^ p1j),
        ibs1,
        inf1: share & (p0i | p0j),
        inf2: share & !ibs1,
        hethet: het_i & het_j,
    }
}

/// One pair's scan over one usable segment.
///
/// Held as a struct so the IBD1 and IBD2 passes can share the per-word masks: computing
/// them is the whole cost of the analysis.
pub struct Scan {
    /// Per-word IBS0 masks, indexed from [`Usable::first_word`].
    ibs0: Vec<u64>,
    /// Per-word IBS1 masks, same indexing.
    ibs1: Vec<u64>,
    /// Per-word IBD1-informative masks, same indexing. See [`MIN_INFORMATIVE`].
    inf1: Vec<u64>,
    /// Per-word IBD2-informative masks, same indexing.
    inf2: Vec<u64>,
    /// Per-word "both heterozygous" masks, same indexing. The `--ibs` IBD2 pass counts
    /// these and nothing else — see [`Scan::ibd2_words`].
    hethet: Vec<u64>,
    seg: Usable,
    /// Head fringe: markers `seg.lo ..< 64*first_word`, as a mask of IBS0 positions
    /// relative to `64*first_word - 64`.
    head_ibs0: u64,
    /// Tail fringe: markers `64*(last_word+1) ..= seg.hi`, as a mask relative to
    /// `64*(last_word+1)`.
    tail_ibs0: u64,
    /// The same two fringes, as masks of **het-vs-hom mismatch** positions.
    ///
    /// The IBD2 pass reads these where the IBD1 pass reads the IBS0 ones: past the
    /// segment's word grid the scan is marker by marker, and each pass stops at its own
    /// breaking marker — an opposite homozygote for IBD1, a mismatch for IBD2. See
    /// [`Scan::ibd2`] and `docs/research/19-ibd2seg-residual.md` §3.
    head_ibs1: u64,
    tail_ibs1: u64,
}

impl Scan {
    /// Compute the per-word disagreement masks for one pair over one usable segment.
    pub fn new(g: &Genotypes, i: usize, j: usize, seg: Usable) -> Scan {
        let (w0, w1) = (seg.first_word(), seg.last_word());
        let nwords = seg.words();
        let mut ibs0 = Vec::with_capacity(nwords);
        let mut ibs1 = Vec::with_capacity(nwords);
        let mut inf1 = Vec::with_capacity(nwords);
        let mut inf2 = Vec::with_capacity(nwords);
        let mut hethet = Vec::with_capacity(nwords);
        for w in w0..w0 + nwords {
            let d = word_diff(g, i, j, w);
            ibs0.push(d.ibs0);
            ibs1.push(d.ibs1);
            inf1.push(d.inf1);
            inf2.push(d.inf2);
            hethet.push(d.hethet);
        }
        // The fringes are the markers of the segment that fall in a word the segment does
        // not wholly own. They take no part in the word scan but they do bound the
        // boundary refinement, so their IBS0 pattern is kept.
        let (head, head_mis) = if nwords == 0 || seg.lo == WORD * w0 {
            (0, 0)
        } else {
            let d = word_diff(g, i, j, w0 - 1);
            let drop = !mask_low(seg.lo - WORD * (w0 - 1));
            (d.ibs0 & drop, d.ibs1 & drop)
        };
        let (tail, tail_mis) = if nwords == 0 || seg.hi == WORD * (w1 + 1) - 1 {
            (0, 0)
        } else {
            let d = word_diff(g, i, j, w1 + 1);
            let keep = mask_low(seg.hi - WORD * (w1 + 1) + 1);
            (d.ibs0 & keep, d.ibs1 & keep)
        };
        Scan {
            ibs0,
            ibs1,
            inf1,
            inf2,
            hethet,
            seg,
            head_ibs0: head,
            tail_ibs0: tail,
            head_ibs1: head_mis,
            tail_ibs1: tail_mis,
        }
    }

    /// How far left an IBD2 call may creep into the partial word before the segment's
    /// word grid: one marker past the **last** mismatch among the segment's own markers
    /// there, or the segment's first marker when that fringe carries none.
    fn head_stop(&self) -> usize {
        if self.head_ibs1 == 0 {
            return self.seg.lo;
        }
        let w0 = self.seg.first_word();
        let last = WORD * (w0 - 1) + (63 - self.head_ibs1.leading_zeros()) as usize;
        (last + 1).max(self.seg.lo)
    }

    /// The mirror: one marker before the **first** mismatch in the partial word after the
    /// grid, or the segment's last marker when that fringe carries none.
    fn tail_stop(&self) -> usize {
        if self.tail_ibs1 == 0 {
            return self.seg.hi;
        }
        let w1 = self.seg.last_word();
        let first = WORD * (w1 + 1) + self.tail_ibs1.trailing_zeros() as usize;
        (first - 1).min(self.seg.hi)
    }

    /// Whether the run of scan words `k0..=k1` carries [`MIN_INFORMATIVE`] markers of
    /// `inf`.
    ///
    /// The window is the run's **own complete words** and nothing else. Markers in the
    /// flanking words the reported segment reaches into lengthen the call but contribute
    /// nothing to the count, and a sub-threshold run is never rescued by them: with 8 or
    /// 9 informative markers in the core the fixture reports nothing however many are
    /// added to the flanks, and with 10 it reports, the length growing independently.
    /// There is no per-word component either — ten markers packed into the first word of
    /// a fourteen-word run pass.
    fn informative(inf: &[u64], k0: usize, k1: usize) -> bool {
        let mut n = 0u32;
        for &m in &inf[k0..=k1] {
            n += m.count_ones();
            if n >= MIN_INFORMATIVE {
                return true;
            }
        }
        false
    }

    /// Whether two runs may be joined across the words `mid` — `20-seglength-floor.md` §4.
    ///
    /// `mid` lists the words the budget is summed over. On the IBD1 pass that is **every**
    /// word between the two runs, gate-refused runs included: `20-…` §11 item 4 left it
    /// undecided and `23-gap-bound.md` §5 settles it, sweeping a refused run's own
    /// het-vs-A1A1 load from 0 to 9 with the budget on its boundary and finding the merge
    /// turn on at 2, which is where that run's markers take `V` from 8 to 10. The
    /// [`MERGE_MAX_WORDS`] cap still counts only the *unusable* words, so a refused run is
    /// stepped over by the cap and paid for by the budget. On the IBD2 pass `mid` is what
    /// `21-…` §3 measured: the words no gate window covers.
    ///
    /// The budget is `cost * (bad - MERGE_FREE) <= X`:
    ///
    /// * IBD1 — `bad` is the opposite homozygotes; `X` is the A1A1/A1A1 markers, unless
    ///   the het-vs-A1A1 ones alone reach [`MIN_INFORMATIVE`], in which case it is those.
    ///   The switch is a bisection, not a guess: with 16 to 40 A1A1/A1A1 markers in the
    ///   interruption, 9 het-vs-A1A1 markers join and 10 do not (§5).
    /// * IBD2 — `bad` is the opposite homozygotes **plus** the het-vs-hom mismatches, and
    ///   `X` is `inf2` (HetHet + A1A1/A1A1), the very count the gate uses (§7).
    fn merge_ok(&self, mid: &[usize], pass2: bool) -> bool {
        let (mut bad, mut x, mut v) = (0u32, 0u32, 0u32);
        for &k in mid {
            bad += self.ibs0[k].count_ones();
            if pass2 {
                // `21-push-merge.md` §5: `X` is the **HetHet** count, with the same
                // switch at `MIN_INFORMATIVE` the IBD1 pass has — `20-…` §7 read it as
                // `inf2` because a HetHet filler was the only one it varied.
                // `inf1 & !ibs1` is the A1A1/A1A1 half of `inf2`; the rest is HetHet.
                bad += self.ibs1[k].count_ones();
                v += (self.inf2[k] & !self.inf1[k]).count_ones();
                x += (self.inf1[k] & !self.ibs1[k]).count_ones();
            } else {
                x += (self.inf1[k] & !self.ibs1[k]).count_ones();
                v += (self.inf1[k] & self.ibs1[k]).count_ones();
            }
        }
        if v >= MIN_INFORMATIVE {
            x = v;
        }
        let cost = if pass2 { MERGE_COST2 } else { MERGE_COST1 };
        cost * bad.saturating_sub(MERGE_FREE) <= x
    }

    /// Join adjacent gate-passing IBD2 runs — `21-push-merge.md` §3 and §4.
    ///
    /// The IBD2 pass merges on the same budget as the IBD1 one but over a different
    /// interruption and with **no cap on its width** (§4: fifteen unusable words join
    /// when the gap and the budget allow, where the IBD1 pass refuses three at any
    /// floor). What separates the two runs is the space between their **gate windows**,
    /// not between the runs themselves (§3): the earlier window ends at `ge_of(b)` — the
    /// one word its right end reaches into — and the later one opens at its gate-start
    /// word `gs`. A word covered by a window is not part of the interruption, and that
    /// holds after *any* usable word, so a gate-refused run's own reach word is skipped
    /// too.
    fn join_runs2(
        &self,
        runs: Vec<(usize, usize)>,
        ok: &[bool],
        mis: &[u32],
        pos: &[i64],
        min_bp: i64,
    ) -> Vec<(usize, usize)> {
        let n = self.nwords();
        let ge_of = |b: usize| {
            if b + 1 < n && self.ibs0_at(b + 1) == 0 && mis[b + 1] != 0 {
                b + 1
            } else {
                b
            }
        };
        let mut out: Vec<(usize, usize)> = Vec::with_capacity(runs.len());
        for (a, b) in runs {
            if let Some(&(pa, pb)) = out.last() {
                let q = ge_of(pb);
                let g2 = (a..=b).find(|&t| mis[t] == 0).unwrap_or(a);
                let covered = |k: usize| ok[k] || (k > 0 && ok[k - 1] && self.ibs0_at(k) == 0);
                let mid: Vec<usize> = (q + 1..a).filter(|&k| !covered(k)).collect();
                if (pb + 1..a).any(|k| !ok[k])
                    && pos[self.marker(g2, 0)] - pos[self.marker(q + 1, 0) - 1] < min_bp
                    && self.merge_ok(&mid, true)
                {
                    *out.last_mut().unwrap() = (pa, b);
                    continue;
                }
            }
            out.push((a, b));
        }
        out
    }

    /// Join adjacent gate-passing runs across a short interruption — `20-…` §2.
    ///
    /// `runs` are the runs that already cleared the gate, in order; `usable` marks every
    /// word of the scan. Two runs join iff at most [`MERGE_MAX_WORDS`] unusable words lie
    /// between them, the gap from the earlier run's last marker to the later run's first
    /// is **strictly** under `--seglength`, and [`Scan::merge_ok`] passes.
    ///
    /// The cap and the budget do **not** read the same words: the cap counts only the
    /// unusable ones, so a gate-refused run between them is stepped over (`20-…` §6),
    /// while the budget is summed over every word in the interruption, that run included
    /// (`23-gap-bound.md` §5). The merged run then takes the gate, the endpoints and the
    /// floor exactly as an unmerged one does.
    fn join_runs(
        &self,
        runs: Vec<(usize, usize)>,
        usable: &[bool],
        pos: &[i64],
        min_bp: i64,
        pass2: bool,
    ) -> Vec<(usize, usize)> {
        let mut out: Vec<(usize, usize)> = Vec::with_capacity(runs.len());
        for (a, b) in runs {
            if let Some(&(pa, pb)) = out.last() {
                let bad_words = (pb + 1..a).filter(|&k| !usable[k]).count();
                let mid: Vec<usize> = (pb + 1..a).collect();
                if bad_words > 0
                    && bad_words <= MERGE_MAX_WORDS
                    && pos[self.marker(a, 0)] - pos[self.marker(pb + 1, 0) - 1] < min_bp
                    && self.merge_ok(&mid, pass2)
                {
                    *out.last_mut().unwrap() = (pa, b);
                    continue;
                }
            }
            out.push((a, b));
        }
        out
    }

    fn nwords(&self) -> usize {
        self.ibs0.len()
    }

    /// IBS0 mask of scan word `k` (0-based within the segment).
    fn ibs0_at(&self, k: usize) -> u64 {
        self.ibs0[k]
    }

    /// Marker index of scan word `k`'s bit `b`.
    fn marker(&self, k: usize, b: u32) -> usize {
        WORD * (self.seg.first_word() + k) + b as usize
    }

    /// IBD1 segments: maximal runs of words with **no** IBS0 at all, refined at the ends.
    ///
    /// Every clause below was **re-measured on an IBD1-native canvas**
    /// (`docs/research/18-ibd1-caller.md`) — the mirror of the `.seg` canvas of
    /// `…/17-seg-caller.md`, with the painted region held IBD2-free by giving every word
    /// thirty-four het-vs-hom mismatches, so the printed `IBD1Seg` reads back the number
    /// of marker intervals called and `IBD2Seg` reads 0.0000 on every fixture. The
    /// geometry that campaign found is exactly the one already here; what it changed is
    /// [`ibd1_pieces`], which is a different function.
    ///
    /// * **The word rule has no error tolerance whatsoever** — a single opposite
    ///   homozygote anywhere in a word breaks it, and nothing else does. On a block of
    ///   eight callable words with `j` consecutive words replaced by one carrying `z`
    ///   opposite homozygotes, `z = 1` splits the block at every `j`, and 64 het-vs-hom
    ///   mismatches, 64 missing calls, 64 HetHet or 64 A1A1/A1A1 never do (§1).
    /// * **A lone bad word is never absorbed** — `j = 1` reads back as two calls, not one,
    ///   whatever the bad word carries (§1, §5). (There is one exception, and it is a
    ///   `--seglength` effect rather than a property of the pass: see §9 and the note
    ///   below.)
    /// * The refinement is **asymmetric**, which is the part no amount of intuition
    ///   produces: the right end runs to the **last** IBS0 marker inside the *next* word
    ///   — not to the first, and not to the run's own last marker — while the left end
    ///   starts one marker **after** the last IBS0 before the run. Swept over thirteen bit
    ///   patterns in the flanking word (§2), the extension is `1 + lastbit` on the right
    ///   and `63 − lastbit` on the left, and `{0, 63}` reads the same as `{63}`, which is
    ///   what makes it the *last*. The scan reads the immediately flanking word and no
    ///   further: an IBS0 two words out moves nothing.
    /// * **There is no push.** [`Scan::ibd2`]'s "every call after the first starts one
    ///   word late" has no analogue here, at any number of calls (§3).
    /// * The gate is [`MIN_INFORMATIVE`] markers of `inf1` over the run's **own** complete
    ///   words: 9 refused / 10 accepted on A1A1/A1A1 and on het-vs-A1A1 alike, split
    ///   across two words, and mixed; HetHet is worth nothing, and a run carrying 9 is
    ///   still refused when the flanking word its call reaches into carries 40 (§4).
    ///
    /// Consecutive segments are then clipped so they cannot overlap, earlier one wins.
    /// Where the run reaches the last complete word, the segment instead creeps into the
    /// trailing fringe marker by marker and stops just before the first IBS0 there.
    ///
    /// **The `--seglength` run merge** (`docs/research/20-seglength-floor.md`) is applied
    /// by [`Scan::runs`] before any of the above: the gate is asked first, and two
    /// surviving runs are joined when at most two unusable words lie between them, the
    /// run-to-run gap is strictly under `--seglength`, and
    /// `4 * (opposite homozygotes - 2) <= X` over those words — where `X` is the
    /// A1A1/A1A1 count unless the het-vs-A1A1 markers alone reach [`MIN_INFORMATIVE`].
    /// See [`Scan::merge_ok`]. It cannot fire at the default floor on real spacings, which
    /// is why `IBD1Seg` was already exact on all 982 corpus rows there; at 5 and 10 Mb this
    /// pass alone takes that column from 910 and 844 to 959 and 960, and the IBD2 side of
    /// the merge (`Scan::join_runs2`, `21-push-merge.md`) carries it the rest of the way to
    /// **982 and 970** by removing IBD2 territory this column would otherwise keep.
    pub fn ibd1(&self, pos: &[i64], min_bp: i64, merge: bool) -> Vec<Called> {
        self.runs(|k| self.ibs0_at(k) == 0, MIN_RUN1, pos, min_bp, merge)
    }

    /// Maximal runs of `good` words, at least `min_run` long, turned into segments.
    ///
    /// Ordering matters and is not obvious: a segment that falls under the `--seglength`
    /// floor is dropped **before** it can clip its successor's start, so a short call
    /// never eats the beginning of the long one behind it.
    fn runs(
        &self,
        good: impl Fn(usize) -> bool,
        min_run: usize,
        pos: &[i64],
        min_bp: i64,
        merge: bool,
    ) -> Vec<Called> {
        let n = self.nwords();
        let usable: Vec<bool> = (0..n).map(&good).collect();
        // The gate is asked **first** (`20-…` §6): a run under [`MIN_INFORMATIVE`] is
        // refused outright, and then lies inside a later interruption rather than ending
        // one — it can never be the endpoint of a merged segment, but does not stop one.
        let mut kept: Vec<(usize, usize)> = Vec::new();
        let mut k = 0usize;
        while k < n {
            if !usable[k] {
                k += 1;
                continue;
            }
            let k0 = k;
            while k < n && usable[k] {
                k += 1;
            }
            let k1 = k - 1;
            if k1 + 1 - k0 >= min_run && Scan::informative(&self.inf1, k0, k1) {
                kept.push((k0, k1));
            }
        }
        if merge {
            kept = self.join_runs(kept, &usable, pos, min_bp, false);
        }
        let mut out: Vec<Called> = Vec::new();
        for (k0, k1) in kept {
            // The window bound (`23-gap-bound.md` §4). The IBD1 window is the run's own
            // complete words — the very span [`Scan::informative`] counts over — and the
            // comparison is **strict**, one unit of `min_bp / 2` tighter than the IBD2
            // pass's. Asked here, after the merge, so a merged run is measured whole.
            let window = pos[self.marker(k1, 63)] - pos[self.marker(k0, 0)];
            if window <= min_bp / WINDOW_FRACTION {
                continue;
            }
            let hi = self.right_end(k1);
            let mut lo = self.left_end(k0);
            if let Some(prev) = out.last() {
                lo = lo.max(prev.hi + 1);
            }
            if lo <= hi && pos[hi] - pos[lo] >= min_bp {
                out.push(Called { lo, hi });
            }
        }
        out
    }

    /// Right end of a run finishing at scan word `k1`.
    ///
    /// The run always reaches **into the word that ended it**: out to that word's *last*
    /// IBS0, or — when the word that ended it carries no IBS0 at all, which only happens
    /// on the IBD2 pass, where runs end on a heterozygote mismatch — all the way through
    /// it. Inverting the `MaxIBD2` column of `--ibs` over the corpus resolves 154 of 158
    /// segments to an interval of the form `[64u, 64v+63]`, which is what forces the
    /// second case; taking the run's own last word instead costs 0.0015 of mean `PropIBD`
    /// error.
    fn right_end(&self, k1: usize) -> usize {
        if k1 + 1 < self.nwords() {
            match self.ibs0_at(k1 + 1) {
                0 => self.marker(k1 + 1, 63).min(self.seg.hi),
                m => self.marker(k1 + 1, 63 - m.leading_zeros()),
            }
        } else if self.tail_ibs0 != 0 {
            WORD * (self.seg.last_word() + 1) + self.tail_ibs0.trailing_zeros() as usize - 1
        } else {
            self.seg.hi
        }
    }

    /// Left end of a run starting at scan word `k0`, one marker past the last IBS0 before
    /// it — or the word boundary when the word before holds no IBS0 at all.
    fn left_end(&self, k0: usize) -> usize {
        if k0 > 0 {
            match self.ibs0_at(k0 - 1) {
                0 => self.marker(k0, 0),
                m => self.marker(k0 - 1, 63 - m.leading_zeros()) + 1,
            }
        } else if self.head_ibs0 != 0 {
            WORD * (self.seg.first_word() - 1) + (63 - self.head_ibs0.leading_zeros()) as usize + 1
        } else {
            self.seg.lo
        }
    }

    /// IBD2 segments: stretches of words that are not [`ibd2_dirty`], with their own
    /// geometry — **not** [`Scan::ibd1`]'s and **not** [`Scan::ibd2_words`]'s.
    ///
    /// Every constant below is a reading off the `.seg`-native canvas of
    /// `docs/research/17-seg-caller.md`: chromosome 1 carries one IBD1 segment so the pair
    /// earns a `.seg` row at all, chromosome 2 is painted one complete word at a time from
    /// an explicit composition and walled with all-IBS0 words, and its uniform spacing is
    /// chosen so that one ulp of the printed `IBD2Seg` is a fifth of a marker gap. The
    /// printed column then reads back **the number of marker intervals called** — and,
    /// because a word-aligned call over `n` words measures `64n − 1`, a total `M` over `w`
    /// words from `c` calls is `64w − c`, so `c = (−M) mod 64` recovers the number of calls
    /// and the number of words exactly. That is what makes the rules below bisections
    /// rather than fits.
    ///
    /// # The rule
    ///
    /// * **A word is usable iff it carries no opposite homozygote and at most
    ///   [`IBD2_HET_DIRTY`] − 1 het-vs-hom mismatches.** One IBS0 anywhere in a word
    ///   disqualifies it at any HetHet density; missing calls, A2A2/A2A2 and A1A1/A1A1
    ///   markers are all irrelevant, and a het-vs-A1A1 mismatch counts exactly like a
    ///   het-vs-A2A2 one (§3).
    /// * **A lone unusable word carrying no IBS0 is absorbed iff the gate passes on both
    ///   sides of it** — the run so far, counted from its gate-start through *this* word,
    ///   and the continuation, counted from the **very next** word (which must therefore be
    ///   mismatch-free) through the words its own right end reaches. The bridge carries no
    ///   constant of its own: it is the ordinary gate below, asked twice (§14). `CyC` is one
    ///   call over three words where `Cyz` is one over two and `CyzC` is one over four
    ///   (§3, §7); `Q(9)dCC` splits where `Q(10)dCC` bridges, and so does `Cy Q(9)` against
    ///   `Cy Q(10)`, which is the same 9/10 bisection on each half independently. Two
    ///   unusable words in a row are never absorbed, under any of 64 compositions.
    /// * **The gate is [`MIN_INFORMATIVE`] markers of `inf2`, counted from the run's first
    ///   *mismatch-free* word** through the words the right end reaches into. Where the
    ///   count starts is not a detail: `zx` and `xz` are the same two words with the same
    ///   total and only `zx` is called (§4). A run with no mismatch-free word at all is
    ///   refused outright, which is why a uniform block at one mismatch per word reports
    ///   nothing however wide it is. `inf2` here is HetHet + A1A1/A1A1 — a het-vs-A1A1
    ///   marker is not informative, bisected at 10 against 20 in §14.3.
    /// * **The endpoints are marker-level, not word-aligned**: a call reaches
    ///   [`IBD2_REACH`] markers past the nearest het-vs-hom mismatch in the word that
    ///   bounds it — right past the *first*, left before the *last*. An opposite
    ///   homozygote blocks that reach **whole-word**: a flanking word carrying one stops
    ///   the call on the run's own boundary whatever bit it sits at, and one in the second
    ///   word out caps the reach at that word's boundary (§5). This asymmetry —
    ///   marker-level for mismatches, word-level for IBS0 — is what made
    ///   `docs/research/14-ibd2-geometry.md` §6.2's "two words and one marker" look
    ///   unexplainable.
    /// * **Every call after the first in a usable segment starts one word late**, counted
    ///   from its gate-start word. It is the *emitted* call that pushes, not the break: a
    ///   run refused by the gate pushes nothing, and the push survives its cause being
    ///   dropped by `--seglength`, so the clip is applied before the length filter (§6).
    /// * **The fringe.** Once an end lands on the grid's own edge the word scan is over and
    ///   a **marker** scan takes over across the partial word beyond it — outwards to the
    ///   segment's own first or last marker, or only as far as the nearest **het-vs-hom
    ///   mismatch** there, stopping one marker short of it (the *last* such marker on the
    ///   left, the *first* on the right). An opposite homozygote in that partial word does
    ///   **not** stop the call, though one anywhere in a complete word inside the grid
    ///   disqualifies that word outright. It is the computed end that snaps, not the run:
    ///   a call whose run opens at the *second* complete word still reaches the fringe when
    ///   the reach above carries its left end back to the edge. All of it is bisected on
    ///   `docs/research/fixtures/fringecanvas.py` — the rig that finally builds a segment
    ///   which does not start on a word boundary — and written up in
    ///   `docs/research/19-ibd2seg-residual.md` §1-§4; `fringe_tests` holds the cases. This
    ///   clause is why `dups`' duplicate pair reads `IBD2Seg 1.0000` and not the
    ///   word-aligned 0.8984, and it is worth 86 exact `IBD2Seg` rows.
    ///
    /// Requiring a whole *word* to be usable rather than only cutting at boundaries is what
    /// makes a parent–offspring pair print `IBD2Seg 0.0000`: PO genotypes disagree far too
    /// often in every word for any of them to qualify.
    ///
    /// # The two length measures — and, separately, the two callers
    ///
    /// `--ibs` measures a call **word-aligned**, from `64u` to `64e+63`, while `.seg`
    /// measures it to the usable segment's own ends. That much is only a ruler, and the
    /// ruler is verified: a duplicate pair in `dups` prints `IBD2Seg 1.0000` and
    /// `Pr_IBD2 0.8984`, and 0.8984 is, to the last digit, the word-aligned total over that
    /// fileset's usable segments divided by the same `D` (357 701 908 / 398 163 465). This
    /// function returns the `.seg` measure.
    ///
    /// The rulers were once thought to be the *whole* difference. They are not, and
    /// `docs/research/17-seg-caller.md` §8 states the negative in full: **`.seg` is not a
    /// quantised confirmation scan.** There is no chunk quantum (sweeping a uniform block's
    /// width gives exactly `64W − 1` at every `W`, with no staircase at any composition), no
    /// confirmation count (a uniform block at one mismatch per word is refused at every
    /// width, yet *one* clean word in front of it buys an unbounded tail), HetHet and
    /// A1A1/A1A1 are interchangeable here where [`Scan::ibd2_words`] ignores the latter
    /// entirely, and nothing is ever cut — a run is called whole or refused. Porting the
    /// chunk geometry over wholesale (`tests/parity/fit/segtry.py`) was tried and rejected:
    /// 709 exact rows against 705 but mean `PropIBD` error 0.00356 against 0.00138. It was
    /// not mis-tuned, it was the wrong kind of rule.
    ///
    /// # Accuracy
    ///
    /// Out of sample on the captured corpus — nothing in it chose a constant — `IBD2Seg` is
    /// exact on **all 982** primary rows, and with the writer rules of [`seg_prop_ibd`] and
    /// `analysis::ibdseg::seg_pair_order` alongside it the whole file is byte-identical on
    /// all thirteen datasets at the default floor. This function's own contribution, scored
    /// on the unrounded scale so the generations are comparable: `IBD2Seg` 822 → 896 → 982,
    /// exact rows 705 → 747 → 806, mean `PropIBD` error 0.00138 → 0.000067 → 0.000023, worst
    /// row 0.2109 → 0.0042 → 0.0001. `tests/parity/fit/engine.py` pins all three as named
    /// parameter bundles (`RETIRED`, `FRINGE18`, `PROP19`), so every one of those numbers
    /// re-runs.
    ///
    /// Neither the **bridge** nor the **fringe** is a fit any longer (`…/17-…` §14,
    /// `…/19-…` §1-§4). The binary is graded against the reference on the canvas batteries
    /// directly: over 6 000 word-aligned canvases — the exhaustive sequences of length ≤ 4
    /// and 5 over `{clean, quiet, 1-mismatch, 2-mismatch}`, of length 4 over an
    /// eight-letter alphabet, and six families of random and "rich" random canvases — this
    /// rule reproduces the reference on **6 000**, where the fitted lookahead it replaces
    /// reproduced 5 723; and on **504** further canvases whose segment does not start on a
    /// word boundary (two unused random seeds plus an exhaustive composition sweep crossed
    /// with six fringe shapes) it reproduces the reference on **504**.
    pub fn ibd2(&self, pos: &[i64], min_bp: i64, merge: bool) -> Vec<Called> {
        let n = self.nwords();
        let (w0, w1) = (self.seg.first_word(), self.seg.last_word());
        if n == 0 {
            return Vec::new();
        }
        let mis: Vec<u32> = self.ibs1.iter().map(|m| m.count_ones()).collect();
        let inf2: Vec<u32> = self.inf2.iter().map(|m| m.count_ones()).collect();
        let usable: Vec<bool> = (0..n).map(|k| !ibd2_dirty(self, k)).collect();
        // The two marker-level stops beyond the word grid — see `Scan::head_stop`.
        let (head_stop, tail_stop) = (self.head_stop(), self.tail_stop());

        // The last word a run ending at `b` reaches into — the whole-word form of the
        // right-hand reach below, which is the window the gate is counted over.
        let ge_of = |b: usize| {
            if b + 1 < n && self.ibs0_at(b + 1) == 0 && mis[b + 1] != 0 {
                b + 1
            } else {
                b
            }
        };
        let gate_ok =
            |g: usize, b: usize| inf2[g..=ge_of(b)].iter().sum::<u32>() >= MIN_INFORMATIVE;

        // **The bridge, and it is the gate asked twice.** A lone unusable word carrying no
        // opposite homozygote is absorbed iff both halves would pass the gate on their own:
        // the run so far, from its gate-start through *this* word (the unusable word's own
        // `inf2` counts — `zyCC` bridges where `zdCC` does not), and the continuation, from
        // the very next word — which must therefore be mismatch-free — through the words
        // its own right end reaches. Read from `usable`, never from the running copy, so
        // two unusable words in a row cannot chain their way in.
        let mut ok = usable.clone();
        let mut gate_start: Option<usize> = None;
        for k in 0..n {
            if usable[k] {
                if gate_start.is_none() && mis[k] == 0 {
                    gate_start = Some(k);
                }
                continue;
            }
            let bridged = match gate_start {
                Some(g)
                    if k > 0
                        && self.ibs0_at(k) == 0
                        && k + 1 < n
                        && usable[k + 1]
                        && mis[k + 1] == 0 =>
                {
                    let mut b = k + 1;
                    while b + 1 < n && usable[b + 1] {
                        b += 1;
                    }
                    gate_ok(g, k - 1) && gate_ok(k + 1, b)
                }
                _ => false,
            };
            if bridged {
                ok[k] = true;
            } else {
                gate_start = None;
            }
        }

        // The gate is asked **first**, so a run it refuses can sit inside a later merge's
        // interruption without ending it (`20-seglength-floor.md` §6). This is the same
        // `gs`/`gate_ok` pair the emit loop below re-asks on the merged run.
        let mut kept: Vec<(usize, usize)> = Vec::new();
        let mut k = 0usize;
        while k < n {
            if !ok[k] {
                k += 1;
                continue;
            }
            let a = k;
            while k < n && ok[k] {
                k += 1;
            }
            let b = k - 1;
            if b + 1 - a < MIN_RUN2 {
                continue;
            }
            match (a..=b).find(|&t| mis[t] == 0) {
                Some(gs) if gate_ok(gs, b) => kept.push((a, b)),
                _ => continue,
            }
        }
        if merge {
            kept = self.join_runs2(kept, &ok, &mis, pos, min_bp);
        }

        let mut out: Vec<Called> = Vec::new();
        // Whether the one-word push is armed. `17-…` §6 read it as "every call after the
        // first", which is what it looks like at the default floor; `21-push-merge.md` §2
        // bisects the condition: a call arms the push only when it is at least **half**
        // the floor long, measured from its own gate-start word rather than from its left
        // end. Whether it survived `--seglength` itself does not enter. Once armed it
        // stays armed for the rest of the usable segment.
        let mut armed = false;
        for (a, b) in kept {
            let (u, v) = (w0 + a, w0 + b);

            // Left: `IBD2_REACH` markers before the *last* mismatch of the word that
            // opened the run — but only if that word is inside the segment and carries no
            // opposite homozygote, and only as far as the next word out, which blocks the
            // reach whole-word when it is outside the segment or carries an IBS0.
            let mut left = WORD * u;
            if a > 0 && self.ibs0_at(a - 1) == 0 && self.ibs1[a - 1] != 0 {
                let last = self.marker(a - 1, 63 - self.ibs1[a - 1].leading_zeros());
                left = last.saturating_sub(IBD2_REACH);
                if a < 2 || self.ibs0_at(a - 2) != 0 {
                    left = left.max(WORD * (u - 1));
                }
            }
            // Once the end reaches the grid's own first marker the word scan is over and
            // the marker scan takes over, whether that moves the end out (the segment's
            // fringe carries no mismatch) or pulls it back in (it does).
            if left <= WORD * w0 {
                left = head_stop;
            }
            // Right: the mirror image, past the *first* mismatch of the word that ended it.
            let mut right = WORD * v + WORD - 1;
            if b + 1 < n && self.ibs0_at(b + 1) == 0 && self.ibs1[b + 1] != 0 {
                let first = self.marker(b + 1, self.ibs1[b + 1].trailing_zeros());
                right = first + IBD2_REACH;
                if b + 2 >= n || self.ibs0_at(b + 2) != 0 {
                    right = right.min(WORD * (v + 2) - 1);
                }
            }
            if right >= WORD * (w1 + 1) - 1 {
                right = tail_stop;
            }

            // The gate, from the run's first mismatch-free word through the one word the
            // right end reaches into. A run without a mismatch-free word is refused
            // outright. This is the same `gate_ok` the bridge asks, on the same window:
            // the reach can spill *markers* into the word after next, but the gate counts
            // whole words and stops at the first of them (§14.2).
            let Some(gs) = (a..=b).find(|&t| mis[t] == 0) else {
                continue;
            };
            if !gate_ok(gs, b) {
                continue;
            }
            // The same window, asked a second time — for its **length** rather than its
            // informative content (`23-gap-bound.md`). It is asked here, at emit and
            // after the merge, not with the gate: a run the bound refuses still merges
            // with its neighbour, and the merged window is then measured whole (§6).
            let window = pos[self.marker(ge_of(b), 63)] - pos[self.marker(gs, 0)];
            if window < min_bp / WINDOW_FRACTION {
                continue;
            }

            if armed {
                left = left.max(WORD * (w0 + gs + 1));
            }
            let mut lo = left.max(self.seg.lo);
            let hi = right.min(self.seg.hi);
            // Consecutive calls may touch, but not overlap.
            if let Some(prev) = out.last() {
                lo = lo.max(prev.hi);
            }
            if lo > hi {
                continue;
            }
            let from_gs = (WORD * (w0 + gs)).clamp(self.seg.lo, hi);
            armed = armed || pos[hi] - pos[from_gs] >= min_bp / PUSH_FRACTION;
            if pos[hi] - pos[lo] >= min_bp {
                out.push(Called { lo, hi });
            }
        }
        out
    }
}

// ---------------------------------------------------------------------------
// The `--ibs` IBD2 pass
// ---------------------------------------------------------------------------

/// Het-vs-hom disagreements that make a word too dirty for the `--ibs` IBD2 scan.
///
/// The same 5 as [`IBD2_HET_DIRTY`], now measured directly rather than inverted out of
/// `MaxIBD2`: on a fixture whose canvas is IBD1 (so no word carries an opposite
/// homozygote) two clean blocks separated by two words holding *n* het-vs-hom mismatches
/// each are reported as one segment for `n <= 4` and as two for `n >= 5`
/// (`docs/research/15-ibs-ibd2-rules.md` §2).
const IBS_IBD2_DIRTY: u32 = 5;

/// HetHet markers one confirmation chunk needs before the `--ibs` pass will keep it.
///
/// **Measured, not inverted, and bisected four separate ways.** The first measurement was
/// a block of complete words in which exactly `k` markers are heterozygous in both samples
/// and the rest are homozygous-reference in everybody: reported iff `k >= 95`, invariant to
/// the block's width (2, 3, 4 and 5 words), its position, the marker spacing, the sample
/// count and the carrier chromosome's length — eight independent bisections all landing on
/// 95 (`docs/research/15-ibs-ibd2-rules.md` §5). `docs/research/16-segment-extension.md` §5
/// then found the same integer from the other side, against an already-established run: a
/// trailing tail of words at `(mismatches, HetHet) = (m, h)` per word is absorbed at
/// `(1, 19)`, `(2, 24)`, `(3, 32)` and `(4, 48)` and refused at `(1, 18)`, `(2, 23)`,
/// `(3, 31)` and `(4, 47)` — that is `5×19 = 95`, `4×24 = 96`, `3×32 = 96`, `2×48 = 96`
/// against `90, 92, 93, 94`, whose intersection is the single integer `(94, 95]`.
///
/// Markers where both samples are homozygous for A1 do not count at all: the same fixture
/// with 200 of them and no HetHet reports nothing, which is what separates this from
/// [`MIN_INFORMATIVE`]'s `inf2`.
const IBS_IBD2_HETHET: u32 = 95;

/// Het-vs-hom mismatches that close one confirmation chunk.
///
/// The quantum the whole scan is built out of. Measured off the staircase of
/// `docs/research/16-segment-extension.md` §4: sweeping the width `W` of a uniform block of
/// words at `(m, h)` per word, the reported end lands on `min(W - 1, 5⌊W/5⌋)` at `m = 1`,
/// and the staircase's period is `⌈5/m⌉` words for every `m` in `1..=5` — five mismatches,
/// however they are spread. §5 sees the same quantum directly: a 20-word all-HetHet prefix
/// buys a trailing tail *exactly five* mismatches and never ten, so the counters are reset
/// at each chunk rather than carried.
const IBS_IBD2_CHUNK_MIS: u32 = 5;

/// Complete words one confirmation chunk must span before it can be confirmed.
///
/// A HetHet count alone does not decide a chunk: 2-word chunks are refused at 122 HetHet
/// while 3-word chunks are confirmed at 96 (`…/16-segment-extension.md` §5, the `m = 3` and
/// `m = 4` rows, refused at *every* `h`). A grid search over `{1, 2, 3, 4, 5}` scores
/// 614/658 constructed filesets at 3 against 599 at 1–2 and 516 at 4.
const IBS_IBD2_CHUNK_WORDS: usize = 3;

/// Mismatches the measured interval may pick up past the last confirmed chunk.
///
/// The reported interval does not stop dead on the confirmed word: it reaches on through
/// the run while it picks up at most this many further mismatches. Visible as the `+1`
/// in §5's covered-word figure, which steps `…, 5, 6, then 10, 11, then 15, 16, …` — every
/// fifth mismatch (the chunk) plus one. A grid search over `{0, 1, 2}` scores 614/658 at 1
/// against 577 at 0 and 563 at 2.
const IBS_IBD2_EXT_MIS: u32 = 1;

/// Complete words the reported interval must span.
///
/// A two-word interval is refused however informative it is — the only way to produce one
/// is a run against the usable segment's last word, and a 128-marker block of pure HetHet
/// there reports nothing while the same block one word earlier (which measures three
/// words) reports.
const IBS_IBD2_MIN_WORDS: usize = 3;

/// One `--ibs` IBD2 call, as a closed **word** interval of the global grid.
///
/// `--ibs` measures a call from `64 * lo` through `64 * hi + 63` — the ruler that makes
/// `dups`' MZ pair read `Pr_IBD2 0.8984` against `IBD2Seg 1.0000`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct WordCall {
    pub lo: usize,
    pub hi: usize,
}

impl Scan {
    /// The IBD2 calls `--ibs` reports, as word intervals of the global grid.
    ///
    /// # Why this is not [`Scan::ibd2`]
    ///
    /// It was supposed to be: one caller, two rulers. The rulers *are* different and that
    /// part holds — `.seg` measures a call to the usable segment's own ends, `--ibs`
    /// measures it word-aligned — but the calls themselves are not the same set, and a
    /// constructed fixture says so outright. Take a `W`-word all-HetHet block bounded by
    /// words in which **every** marker is an opposite homozygote, on a chromosome that is
    /// otherwise all such words. `--ibdseg` reports an `IBD2Seg` worth exactly the block;
    /// `--ibs` reports a `MaxIBD2` worth the **whole usable segment**, IBS0 words and all
    /// (`docs/research/15-ibs-ibd2-rules.md` §3). Whatever the `.seg` pass does with an
    /// opposite homozygote, this pass does not stop at one.
    ///
    /// # The shape of the rule: a quantised confirmation scan
    ///
    /// A run is *not* accepted or refused whole, and it is not extended word by word
    /// either. It is confirmed in **chunks of [`IBS_IBD2_CHUNK_MIS`] het-vs-hom
    /// mismatches**, each of which must independently carry [`IBS_IBD2_HETHET`] HetHet
    /// markers over at least [`IBS_IBD2_CHUNK_WORDS`] words, and the reported interval
    /// stops at the last chunk that was confirmed. `docs/research/16-segment-extension.md`
    /// derives this and rules out the three shapes that were guessed before it:
    ///
    /// * **Not a running score / X-drop.** Credit does not accumulate. A 20-word
    ///   all-HetHet prefix — 1 280 HetHet markers — buys the tail *exactly five* trailing
    ///   mismatches and never ten, at every prefix length. The counters reset at each
    ///   chunk boundary (§5).
    /// * **Not an HMM or a Viterbi path.** Both are aggregate-optimal over the words they
    ///   see, so a uniform block bounded by walls admits no interior optimum; the
    ///   reference reports **partial** intervals of uniform blocks, and where it cuts them
    ///   depends on the block's width (§3, §4).
    /// * **Not greedy word-at-a-time seed extension.** The scan *is* one forward pass, but
    ///   the unit of extension is the chunk, accepted or refused whole.
    ///
    /// The two facts that made the residual look non-causal both fall out of the quantum:
    /// the staircase `e = min(W - 1, 5⌊W/5⌋)` is the chunk period, and the
    /// order-dependence (`8 clean + 8 dirty-ish` reported whole, the reverse not) is that
    /// HetHet is counted *within* a chunk, from wherever that chunk happens to start.
    ///
    /// # The rules, each with the fixture that fixes it
    ///
    /// * **A word breaks the run iff it carries [`IBS_IBD2_DIRTY`] het-vs-hom
    ///   mismatches.** Opposite homozygotes are irrelevant — a gap word of 64 IBS0
    ///   markers and nothing else is scanned straight through, and so is one of 64 missing
    ///   calls, while five het-vs-hom mismatches split the run.
    /// * **One dirty word between two clean ones is absorbed**, two in a row are not.
    /// * **The scan runs one word past the run's last clean word** — that word is what
    ///   makes the mismatch counter fire, and its HetHet counts towards the chunk it
    ///   closes — but never past the usable segment's own last word.
    /// * **A refused chunk ends the segment at the last confirmation and starts a new
    ///   one.** Where exactly it restarts is the least-supported clause here: see
    ///   `restart` below.
    /// * **The interval reaches past the last confirmed chunk while picking up at most
    ///   [`IBS_IBD2_EXT_MIS`] further mismatch**, and snaps to the segment's last word
    ///   when the run ends within two words of it.
    /// * **...and the whole HetHet test is waived where the run reaches the segment's last
    ///   two words**: a block of 384 markers with no HetHet at all is reported when it
    ///   ends on `w1` or `w1 - 1` and refused one word earlier. The exemption follows the
    ///   usable segment, not the genotype array — it is the same whether the segment is
    ///   the first chromosome or the last.
    /// * **The interval must span [`IBS_IBD2_MIN_WORDS`] words.**
    ///
    /// Consecutive calls are clipped so they cannot share a word, earlier one wins.
    ///
    /// # Accuracy, and the one clause that is fitted rather than measured
    ///
    /// Out of sample this reproduces **`MaxIBD2` 158/158 and `Pr_IBD2` 158/158** — both
    /// `--ibs` IBD2 columns exactly, on every corpus pair the reference grades, over ten
    /// datasets none of which chose a constant. The rule it replaced scored 148 and 100.
    /// On 658 constructed filesets, 200 of them *random* word sequences the fit never saw,
    /// it reproduces the reference's exact interval 614 times (93.3 %); the residual is
    /// concentrated on sequences that alternate 20- and 64-mismatch words with near-empty
    /// ones, which nothing in the corpus resembles.
    ///
    /// The **restart** clause is the exception and is flagged in place below: the corpus
    /// cannot see it (restarting after the refusing word scores identically on all 316
    /// corpus targets), and it is fitted from fourteen irregular constructed patterns.
    pub fn ibd2_words(&self) -> Vec<WordCall> {
        let n = self.nwords();
        if n == 0 {
            return Vec::new();
        }
        let w0 = self.seg.first_word();
        // Everything below is in *local* word indices, `0 ..= last`; `last` is the usable
        // segment's own final complete word, which several rules key off.
        let last = n - 1;
        let mis: Vec<u32> = self.ibs1.iter().map(|m| m.count_ones()).collect();
        let hh: Vec<u32> = self.hethet.iter().map(|m| m.count_ones()).collect();

        let clean: Vec<bool> = mis.iter().map(|&m| m < IBS_IBD2_DIRTY).collect();
        // A lone dirty word between two clean ones is absorbed. Read from `clean`, never
        // from the running copy, so two dirty words in a row cannot chain their way in.
        let mut ok = clean.clone();
        for k in 1..n.saturating_sub(1) {
            if !clean[k] && clean[k - 1] && clean[k + 1] {
                ok[k] = true;
            }
        }

        /// How far past the last confirmed word `conf` the measured interval reaches:
        /// on through the run's words until the mismatches picked up would exceed
        /// [`IBS_IBD2_EXT_MIS`]. `b` is the run's own last clean word, never crossed.
        fn extend(mis: &[u32], conf: usize, b: usize) -> usize {
            let mut cum = 0u32;
            let mut e = conf;
            for (k, &m) in mis.iter().enumerate().take(b + 1).skip(conf + 1) {
                cum += m;
                if cum > IBS_IBD2_EXT_MIS {
                    break;
                }
                e = k;
            }
            e
        }

        let mut raw: Vec<(usize, usize)> = Vec::new();
        let mut k = 0usize;
        while k < n {
            if !ok[k] {
                k += 1;
                continue;
            }
            let a = k;
            while k < n && ok[k] {
                k += 1;
            }
            let b = k - 1;
            // One word past the run's last clean word, clamped to the segment.
            let scan_last = (b + 1).min(last);
            // The run reaches the usable segment's own last two words.
            let exempt = b + 1 >= last;

            let mut u = a; // where the segment currently being built started
            let mut acc_mis = 0u32; // mismatches in the open chunk
            let mut acc_het = 0u32; // HetHet in the open chunk
            let mut cstart = a; // first word of the open chunk
            let mut conf: Option<usize> = None; // last confirmed word
            let mut last_mis: Option<usize> = None; // last word holding a mismatch
            let mut s = a;
            while s <= scan_last {
                let (m, h) = (mis[s], hh[s]);
                let before = acc_mis;
                acc_mis += m;
                acc_het += h;
                if acc_mis >= IBS_IBD2_CHUNK_MIS {
                    // The chunk closes on this word. A chunk closed by the usable
                    // segment's own last word is exempt from the HetHet test, the same
                    // tail exemption that waives it for a run reaching `last`.
                    let confirmed = (acc_het >= IBS_IBD2_HETHET
                        && s + 1 - cstart >= IBS_IBD2_CHUNK_WORDS)
                        || (exempt && s >= scan_last);
                    if confirmed {
                        conf = Some(s);
                        acc_mis = 0;
                        acc_het = 0;
                        last_mis = None;
                        cstart = s + 1;
                    } else {
                        if let Some(c) = conf {
                            raw.push((u, extend(&mis, c, b)));
                        }
                        // **The least-supported clause in this function.** Where the word
                        // that closed the chunk holds only the chunk's fifth mismatch, the
                        // next segment opens after the *fourth* mismatch's word instead of
                        // after the refusing word. Fourteen irregular constructed patterns
                        // separate the two (`[2,2,1] -> 2` against `[2,2,2] -> 3`,
                        // `[4,1] -> 1` against `[4,4] -> 2`); the corpus cannot, and scores
                        // identically either way. Treat it as fitted, not measured.
                        u = match last_mis {
                            Some(lm) if before + 1 == IBS_IBD2_CHUNK_MIS && m == 1 => lm + 1,
                            _ => s + 1,
                        };
                        acc_mis = 0;
                        acc_het = 0;
                        conf = None;
                        last_mis = None;
                        cstart = u;
                        s = u;
                        continue;
                    }
                }
                if m > 0 {
                    last_mis = Some(s);
                }
                s += 1;
            }
            if exempt {
                raw.push((u, last));
            } else if let Some(c) = conf {
                // A run ending exactly two words short of the segment's last word still
                // takes it; otherwise the interval stops at the confirmed end plus its
                // one-mismatch overhang.
                raw.push((
                    u,
                    if b + 2 >= last {
                        last
                    } else {
                        extend(&mis, c, b)
                    },
                ));
            }
        }

        let mut out: Vec<WordCall> = Vec::new();
        let mut prev_hi: Option<usize> = None;
        for (mut lo, hi) in raw {
            if lo > hi {
                continue;
            }
            if let Some(p) = prev_hi {
                lo = lo.max(p + 1);
            }
            if lo > hi || hi + 1 - lo < IBS_IBD2_MIN_WORDS {
                continue;
            }
            prev_hi = Some(hi);
            out.push(WordCall {
                lo: w0 + lo,
                hi: w0 + hi,
            });
        }
        out
    }
}

/// `(1 << n) - 1`, saturating at all-ones so `n == 64` is not UB-adjacent.
fn mask_low(n: usize) -> u64 {
    if n >= 64 {
        u64::MAX
    } else {
        (1u64 << n) - 1
    }
}

// ---------------------------------------------------------------------------
// Aggregation
// ---------------------------------------------------------------------------

/// One pair's segment summary — the four columns of `<prefix>.seg`.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PairSegments {
    /// Total base pairs called IBD1 but not IBD2.
    pub ibd1_bp: i64,
    /// Total base pairs called IBD2.
    pub ibd2_bp: i64,
    /// Longest single reported segment, for the ">10Mb" pair filter.
    pub longest_bp: i64,
}

impl PairSegments {
    /// `IBD1Seg` — π1.
    pub fn ibd1_seg(&self, denom: i64) -> f64 {
        if denom == 0 {
            0.0
        } else {
            self.ibd1_bp as f64 / denom as f64
        }
    }

    /// `IBD2Seg` — π2.
    pub fn ibd2_seg(&self, denom: i64) -> f64 {
        if denom == 0 {
            0.0
        } else {
            self.ibd2_bp as f64 / denom as f64
        }
    }

    /// `PropIBD = IBD2Seg + IBD1Seg/2`, in full precision.
    ///
    /// This is the value the **`.kin` family** prints and the value every *decision* is
    /// made on — [`inf_type`], [`reported_at_degree`], `--unrelated`'s greedy and
    /// `--related`'s `Error` grader all take it. `<prefix>.seg` prints a different one;
    /// see [`seg_prop_ibd`], which is a formatting rule and nothing more.
    pub fn prop_ibd(&self, denom: i64) -> f64 {
        self.ibd2_seg(denom) + self.ibd1_seg(denom) / 2.0
    }

    /// Whether the pair survives the fixed ">10Mb" filter and so gets a row at all.
    ///
    /// The console text says `>10Mb`; the binary means `>=`. Bisected on a fixture whose
    /// only segment could be sized to the base pair: 9 990 000 bp is absent,
    /// 10 000 000 bp is present.
    pub fn reported(&self) -> bool {
        self.longest_bp >= LONG_SEGMENT_BP
    }
}

/// Scan one pair across every usable segment and aggregate.
///
/// `seglength_bp` is `--seglength` in base pairs; segments shorter than it are neither
/// reported nor counted, but they still bound the pair filter's "longest segment" the
/// same way the reference's console text implies ("not reported/utilized"). It is
/// [`ibd1_pieces`], not the IBD1 call, that the floor is applied to on the IBD1 side.
pub fn pair_segments(
    g: &Genotypes,
    pos: &[i64],
    segs: &[Usable],
    i: usize,
    j: usize,
    seglength_bp: i64,
) -> PairSegments {
    let mut acc = PairSegments::default();
    for &seg in segs {
        if seg.words() == 0 {
            continue;
        }
        let scan = Scan::new(g, i, j, seg);
        let ibd2 = scan.ibd2(pos, seglength_bp, true);
        let ibd1 = scan.ibd1(pos, seglength_bp, true);
        // The pair-reporting filter reads the conditioned, floor-dependent merged calls.
        // A held-out distant pair has two sub-10 Mb IBD1 calls at 3 Mb, then one 14.6 Mb
        // merged call at 5/10 Mb; KING reports it only at the raised floors. A separate
        // canvas pins the same rule on IBD2. The older unconditioned merge candidate did
        // invent hundreds of pairs, but the measured merge above does not.
        for c in &ibd2 {
            acc.longest_bp = acc.longest_bp.max(pos[c.hi] - pos[c.lo]);
        }
        for c in &ibd1 {
            acc.longest_bp = acc.longest_bp.max(pos[c.hi] - pos[c.lo]);
        }
        for c in &ibd2 {
            acc.ibd2_bp += pos[c.hi] - pos[c.lo];
        }
        for c in &ibd1 {
            // IBD1 is reported as the part of an IBD1 call that is not already IBD2:
            // a duplicate pair is IBD1 everywhere by the IBS0 rule yet reports
            // `IBD1Seg 0.0000`. The pieces are separate segments and each faces the
            // `--seglength` floor on its own — see [`ibd1_pieces`].
            for f in ibd1_pieces(*c, &ibd2) {
                let len = pos[f.hi] - pos[f.lo];
                if len >= seglength_bp {
                    acc.ibd1_bp += len;
                }
            }
        }
    }
    acc
}

/// An IBD1 call with the IBD2 calls cut out of it — the pieces `IBD1Seg` actually sums.
///
/// **Measured on the IBD1-native canvas** (`docs/research/18-ibd1-caller.md` §6), and the
/// two clauses below are what took `IBD1Seg` from 826 exact corpus rows to all 982:
///
/// * the cut is at **marker** granularity and **excludes the IBD2 call's own end
///   markers** — a call `[lo, hi]` cut by an IBD2 call `[a, b]` leaves `[lo, a-1]` and
///   `[b+1, hi]`, so each piece is one marker gap shorter than the naive
///   "length minus overlap". On the canvas `K·B³·K` the reference reports 63 marker
///   intervals where the naive subtraction gives 64, and `K⁴·B⁴·K` reports 161 + 63 where
///   it gives 162 + 64; a graded (non-uniform) ruler picks the same endpoints out of five
///   candidate conventions.
/// * every piece then faces the **`--seglength` floor on its own**. On `K·B³·K` (whose
///   piece is exactly 4 410 000 bp) the floor is bisected to the base pair: at
///   `--seglength 4.410000` the piece is counted, at `4.410001` it is dropped and
///   `IBD1Seg` reads zero even though the IBD1 call it came from spans 26.8 Mb. So this
///   floor is applied to the *pieces*, not to the call.
///
/// `others` are the IBD2 calls that **survived** `--seglength`: a dropped IBD2 call is not
/// subtracted at all (canvas `K·B·K`, where raising the floor past the IBD2 call's length
/// takes `IBD2Seg` to zero and hands the whole call back to `IBD1Seg`).
///
/// They are ordered and may touch but never overlap, which is what lets this walk them
/// once.
fn ibd1_pieces(c: Called, others: &[Called]) -> Vec<Called> {
    let mut out = Vec::new();
    let mut cur = c.lo;
    for o in others {
        if o.hi < c.lo || o.lo > c.hi {
            continue;
        }
        if o.lo > cur {
            out.push(Called {
                lo: cur,
                hi: o.lo - 1,
            });
        }
        cur = cur.max(o.hi + 1);
    }
    if cur <= c.hi {
        out.push(Called { lo: cur, hi: c.hi });
    }
    out
}

// ---------------------------------------------------------------------------
// The `.seg` writer's own PropIBD
// ---------------------------------------------------------------------------

/// The `PropIBD` **`<prefix>.seg` prints** — computed from the two columns beside it,
/// after they have been rounded to the four decimals the file shows.
///
/// ```text
/// i1 = the integer <prefix>.seg prints as IBD1Seg, scaled by 10 000
/// i2 = the same for IBD2Seg
/// PropIBD = i2 * 1e-4 + i1 * 5e-5          , printed "%.4lf"
/// ```
///
/// # Why this is not the same number as [`Segments::prop_ibd`]
///
/// The reference contradicts itself, and that is the whole finding. Run it once:
///
/// ```text
/// king -b bigish.bed --related --degree 2 --ibdseg --cpus 1 --prefix r
/// ```
///
/// **147** pairs land in both `r.kin` and `r.seg`; **all 147** carry identical `IBD1Seg`
/// and `IBD2Seg` in the two files, and **43** carry a different `PropIBD` — e.g.
/// `IBD1Seg 0.4885 / IBD2Seg 0.2974` prints `0.5417` in `.kin` and `0.5416` in `.seg`,
/// and `0.3852 / 0.3123` prints `0.5048` against `0.5049`. Same invocation, same pair,
/// same printed inputs, two answers. No single expression can match both files, so the
/// two writers get one each: `.kin` and its relatives take `prop_ibd` (full precision,
/// byte-exact on all 4 805 corpus rows), and `.seg` takes this.
///
/// # How it was pinned, and how strong the evidence is
///
/// Entirely from the reference's own captured output — the rule is a function of two
/// printed columns, so it needs no genotypes to test:
///
/// * Over **all 4 172** `.seg` rows the corpus captures, `PropIBD` is consistent with
///   *some* rounding of `IBD2Seg + IBD1Seg/2` read off the printed columns: **0
///   refutations**. 2 859 of those rows determine the answer outright; the other 1 313
///   land on an exact decimal half, where the reference rounds up 1 099 times and down
///   214, so the tie is not a convention — it is arithmetic.
/// * That the inputs are the *printed* columns and not the underlying totals is what
///   1.7 % of rows decide: on those the full-precision value sits at least half a printed
///   ulp away from the printed-column combination, so the full-precision hypothesis would
///   have been refuted on ≈ 71 of the 4 172. It was refuted on **0**.
/// * The expression is not free. `(i1 + 2·i2)/20000`, `i2/10000 + i1/20000`,
///   `(i1/2 + i2)/10000`, `i2 + i1·0.5` on the printed doubles, `(i1+2·i2)·5e-5`, and
///   integer round-half-up all reproduce between 3 804 and 4 086 of the 4 172. **This one
///   reproduces 4 172.** The 1 313 ties are decided by whether `i2 * 1e-4 + i1 * 5e-5`
///   lands infinitesimally above or below the exact decimal half, and it agrees with the
///   reference on every one of them.
/// * The rule is `.seg`'s alone. Applied to the reference's own `.kin` it reproduces
///   3 714 of 4 248 rows, `.kin0` 236 of 302, `X.kin` 78 of 90, `cluster.kin` 139 of 165
///   — all four of which open-king already reproduces byte for byte with `prop_ibd`.
///
/// # What it is worth
///
/// On the primary 3 Mb capture, byte-exact `.seg` rows go **806 → 982 of 982**, mean
/// `PropIBD` error 0.000018 → **0.000000**, worst 0.0001 → **0.0000**. Held out at the
/// two floors that had no part in finding it, row-exactness rises to exactly the number
/// of rows whose two estimate columns are already right — 755 → **900** of 900 at 5 Mb
/// and 713 → **832** of 832 at 10 Mb. `PropIBD` now contributes **no error at all** at
/// any floor: whatever `.seg` still gets wrong is `IBD1Seg` or `IBD2Seg`.
///
/// It changes no decision anywhere. `InfType`, the `--degree` filter, `--unrelated` and
/// `--related`'s `Error` all keep reading [`Segments::prop_ibd`]; this value reaches one
/// column of one file.
pub fn seg_prop_ibd(ibd1_seg: f64, ibd2_seg: f64) -> f64 {
    printed_units(ibd2_seg) as f64 * 1e-4 + printed_units(ibd1_seg) as f64 * 5e-5
}

/// The integer a `{:.4}` column shows, scaled by 10 000 — read back off the formatter
/// itself so it cannot drift from the digits actually written beside it.
fn printed_units(x: f64) -> i64 {
    let s = format!("{x:.4}");
    let digits: i64 = s
        .bytes()
        .filter(u8::is_ascii_digit)
        .fold(0i64, |a, b| a * 10 + i64::from(b - b'0'));
    if s.starts_with('-') {
        -digits
    } else {
        digits
    }
}

// ---------------------------------------------------------------------------
// InfType
// ---------------------------------------------------------------------------

/// Relationship label from the segment estimates — the `InfType` column.
///
/// First match wins, on the **unrounded** f64 estimates. The degree cut-points are the
/// ones KING's own emitted R script states (`<prefix>_ibd1vsibd2.R`); the six literal
/// decimals were bracketed against the reference binary on synthetic pairs with
/// prescribed (π1, π2) — `0.32` to within (0.3199, 0.3201), `0.15` to (0.1500, 0.1508],
/// `0.96` to (0.9599, 0.9601].
///
/// There are **two** full-sib clauses, and missing the second is the easy mistake: a pair
/// at π = 0.33 with π2 = 0.20 is `FS`, not `2nd`, even though it is below the 2^-1.5 line
/// that clause A tests. Together the clauses reproduce every one of the 8 722 `InfType`
/// values in the captured corpus from that row's own printed columns.
///
/// Two further traps: the binary writes `Dup/MZ` where the manual says `Dup/MZTwin`, and
/// the `2nd` bucket has **no upper bound**, so a pair at π = 0.45 with π2 = 0 is `2nd`.
pub fn inf_type(pi1: f64, pi2: f64, prop: f64) -> &'static str {
    const D1: f64 = 0.353_553_390_593_273_8; // 2^-1.5
    const D2: f64 = 0.176_776_695_296_636_9; // 2^-2.5
    const D3: f64 = 0.088_388_347_648_318_45; // 2^-3.5
    const D4: f64 = 0.044_194_173_824_159_22; // 2^-4.5
    if pi2 > 0.7 {
        "Dup/MZ"
    } else if pi1 + pi2 > 0.96 || (pi1 + pi2 > 0.9 && pi2 < 0.08) {
        "PO"
    } else if (prop > D1 && pi2 >= 0.08) || (prop > 0.32 && pi2 > 0.15) {
        "FS"
    } else if prop > D2 {
        "2nd"
    } else if prop > D3 {
        "3rd"
    } else if prop > D4 {
        "4th"
    } else {
        "UN"
    }
}

/// `2^-(|d| + 0.5)` — the `PropIBD` band edge `--degree d` compares against.
///
/// The magnitude only; [`reported_at_degree`] owns the direction of the comparison, which
/// is not the same for a negative `d`.
pub fn degree_cutoff(degree: i32) -> f64 {
    2f64.powf(-(f64::from(degree.unsigned_abs()) + 0.5))
}

/// `IBD2Seg` at or above which `--degree 1` reports a pair whatever its `PropIBD`.
///
/// The same 0.08 the binary's own emitted R script draws as a horizontal line and uses in
/// its `d1.FS` predicate.
pub const FIRST_DEGREE_IBD2: f64 = 0.08;

/// Whether `--degree d` reports a pair with these segment estimates.
///
/// Three branches, all measured against the reference over 38 298
/// (dataset, `--seglength`, `--degree`, pair) cases spanning `d` from −6 to 6, with no
/// disagreement. Reproduce with
/// `python3 tests/parity/probes/degree_filter.py --ref <reference king>`:
///
/// * **`d == 0`** — no filter. An integer option carries its own "unset", so an absent
///   `--degree` and an explicit `--degree 0` are the same thing here as they are in the
///   banner.
/// * **`d > 0`** — `PropIBD > 2^-(d+0.5)`, *or*, **at `d == 1` only**,
///   `IBD2Seg >= `[`FIRST_DEGREE_IBD2`]. The second clause is not a rounding of the
///   first and it does not generalise: a constructed pair sharing one IBD2 block and
///   nothing else, at `IBD2Seg = PropIBD = 0.0981`, is reported at `--degree 1`,
///   **not** reported at `--degree 2`, and reported again at `--degree 3`, where
///   `0.0981 > 2^-3.5` carries it on the `PropIBD` clause alone
///   (`docs/research/fixtures/gate8.py`). The corpus cannot see this clause at all —
///   a real first-degree pair has `IBD1Seg ≈ 0.5`, so its `PropIBD` clears `2^-1.5`
///   anyway, and over 52 974 corpus cases no pair has `IBD2Seg` strictly between 0 and
///   0.1089. Sweeping the fixture's block length brackets the constant to
///   (0.0789, 0.0812], which contains the R script's own literal `0.08`.
/// * **`d < 0`** — `PropIBD <= 2^-(|d|+0.5)`: the comparison inverts, so `--degree -2`
///   reports the complement of `--degree 2`. On `bigish` the two report 321 and 442 of
///   763 pairs, and 321 + 442 = 763 exactly; `multifam` is 62 and 43 of 104, one over,
///   the overlap being the single pair the `IBD2Seg` clause admits at `d == 1`.
pub fn reported_at_degree(degree: i32, ibd2_seg: f64, prop_ibd: f64) -> bool {
    match degree {
        0 => true,
        d if d < 0 => prop_ibd <= degree_cutoff(d),
        1 => prop_ibd > degree_cutoff(1) || ibd2_seg >= FIRST_DEGREE_IBD2,
        d => prop_ibd > degree_cutoff(d),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn uniform(n: usize, chr: i64, step: i64) -> (Vec<i64>, Vec<i64>) {
        (
            vec![chr; n],
            (0..n).map(|i| 1_000_000 + i as i64 * step).collect(),
        )
    }

    #[test]
    fn ibd1_pieces_exclude_the_ibd2_call_s_own_end_markers() {
        // `docs/research/18-ibd1-caller.md` §6.1: the canvas `K B3 K` reports 63 marker
        // intervals where "length minus overlap" gives 64, and `K4 B4 K` 161 + 63 where
        // it gives 162 + 64.
        let c = Called { lo: 0, hi: 383 };
        let ibd2 = [Called { lo: 0, hi: 319 }];
        assert_eq!(ibd1_pieces(c, &ibd2), vec![Called { lo: 320, hi: 383 }]);
        let c = Called { lo: 0, hi: 639 };
        let ibd2 = [Called { lo: 162, hi: 575 }];
        assert_eq!(
            ibd1_pieces(c, &ibd2),
            vec![Called { lo: 0, hi: 161 }, Called { lo: 576, hi: 639 }]
        );
    }

    #[test]
    fn ibd1_pieces_survive_touching_and_enclosing_ibd2_calls() {
        // `Scan::ibd2` clips consecutive calls to `lo = max(lo, prev.hi)`, so two IBD2
        // calls may share a marker; and a call may swallow the IBD1 one whole.
        let c = Called { lo: 0, hi: 255 };
        let touching = [Called { lo: 64, hi: 128 }, Called { lo: 128, hi: 191 }];
        assert_eq!(
            ibd1_pieces(c, &touching),
            vec![Called { lo: 0, hi: 63 }, Called { lo: 192, hi: 255 }]
        );
        assert!(ibd1_pieces(c, &[Called { lo: 0, hi: 255 }]).is_empty());
        assert_eq!(ibd1_pieces(c, &[]), vec![c]);
    }

    #[test]
    fn a_gap_of_exactly_one_megabase_does_not_cut() {
        let (chr, mut pos) = uniform(1280, 1, 50_000);
        for p in pos.iter_mut().skip(641) {
            *p += 1_000_000 - 50_000;
        }
        assert_eq!(usable_segments(&chr, &pos).len(), 1);
    }

    #[test]
    fn a_gap_over_one_megabase_cuts() {
        let (chr, mut pos) = uniform(1280, 1, 50_000);
        for p in pos.iter_mut().skip(641) {
            *p += 1_000_001 - 50_000;
        }
        let segs = usable_segments(&chr, &pos);
        assert_eq!(segs.len(), 2);
        assert_eq!(segs[0].n_snp(), 641);
        assert_eq!(segs[1].n_snp(), 639);
    }

    #[test]
    fn a_word_spanning_over_ten_megabases_cuts() {
        // 64 * 156_250 == 10_000_000 exactly, which is not "over".
        let (chr, pos) = uniform(1280, 1, 156_250);
        assert_eq!(usable_segments(&chr, &pos).len(), 1);
        let (chr, pos) = uniform(1280, 1, 156_251);
        assert!(usable_segments(&chr, &pos).is_empty());
    }

    #[test]
    fn four_complete_words_is_not_enough() {
        // 256 markers starting at index 0 are exactly four words.
        let (chr, pos) = uniform(256, 1, 100_000);
        assert!(usable_segments(&chr, &pos).is_empty());
        let (chr, pos) = uniform(320, 1, 100_000);
        assert_eq!(usable_segments(&chr, &pos).len(), 1);
    }

    #[test]
    fn the_length_test_is_word_aligned() {
        // 2000 markers: complete words 0..=30, so the tested span stops at marker 1983.
        let (chr, pos) = uniform(2000, 1, 5_040); // aligned 9.994 Mb, full 10.075 Mb
        assert!(usable_segments(&chr, &pos).is_empty());
        let (chr, pos) = uniform(2000, 1, 5_043); // aligned 10.000 Mb
        assert_eq!(usable_segments(&chr, &pos).len(), 1);
    }

    #[test]
    fn word_alignment_decides_between_two_chromosomes() {
        // The `dups` shape: chr14 holds four complete words, chr15 five, even though
        // chr14 has more markers. Indices 7622..7993 and 7994..8348.
        let mut chr = vec![13i64; 7622];
        let mut pos: Vec<i64> = (0..7622).map(|i| i as i64 * 50_000).collect();
        chr.extend(std::iter::repeat_n(14, 372));
        pos.extend((0..372).map(|i| 1_000_000 + i as i64 * 50_000));
        chr.extend(std::iter::repeat_n(15, 355));
        pos.extend((0..355).map(|i| 1_000_000 + i as i64 * 50_000));
        let segs = usable_segments(&chr, &pos);
        let chrs: Vec<i64> = segs.iter().map(|s| s.chr).collect();
        assert!(chrs.contains(&15), "chr15 has five complete words");
        assert!(!chrs.contains(&14), "chr14 has only four");
    }

    /// Every case here is a real row of the reference's own `<prefix>.seg`, and every one
    /// of them is an exact decimal tie — `IBD2Seg + IBD1Seg/2` lands on `x.xxxx5`. They
    /// do **not** all go the same way, which is the point: the direction is arithmetic,
    /// not a rounding convention, and any "round half up"/"half even" rule fails half of
    /// this list. See `seg_prop_ibd`.
    #[test]
    fn the_seg_writer_resolves_exact_ties_in_both_directions() {
        let p = |a: f64, b: f64| format!("{:.4}", seg_prop_ibd(a, b));
        // rounds DOWN
        assert_eq!(p(0.4885, 0.2974), "0.5416"); // .kin prints 0.5417 for this pair
                                                 // rounds UP
        assert_eq!(p(0.7151, 0.0000), "0.3576");
        assert_eq!(p(0.5207, 0.1808), "0.4412"); // .kin prints 0.4411
        assert_eq!(p(0.4987, 0.0000), "0.2494");
        assert_eq!(p(0.4461, 0.2703), "0.4934");
        assert_eq!(p(0.4135, 0.0000), "0.2068");
        // and a row that is not a tie at all, so every candidate rule agrees
        assert_eq!(p(0.3852, 0.3123), "0.5049"); // .kin prints 0.5048
    }

    /// The full-precision value is a *different number*, and the difference is what the
    /// two writers disagree about. Keep them distinct.
    #[test]
    fn the_two_prop_ibd_rules_are_not_the_same_function() {
        // 0.4885 / 0.2974: full precision says 0.54165 -> "0.5417" only if the underlying
        // totals push it there; the printed-column rule pins it at 0.5416.
        assert_ne!(
            format!("{:.4}", seg_prop_ibd(0.4885, 0.2974)),
            format!("{:.4}", 0.2974_f64 + 0.4885_f64 / 2.0 + 1e-9)
        );
        // ...and it is a pure function of the two printed columns: anything inside the
        // same 4-dp cell gives the same answer.
        assert_eq!(
            seg_prop_ibd(0.48851, 0.29739),
            seg_prop_ibd(0.48849, 0.29741)
        );
    }

    #[test]
    fn inf_type_bands() {
        assert_eq!(inf_type(0.0, 1.0, 1.0), "Dup/MZ");
        assert_eq!(inf_type(1.0, 0.0, 0.5), "PO");
        assert_eq!(inf_type(0.4002, 0.3238, 0.5239), "FS");
        // The 2nd bucket is not bounded above.
        assert_eq!(inf_type(0.8962, 0.0, 0.4481), "2nd");
        assert_eq!(inf_type(0.2097, 0.0, 0.1048), "3rd");
        assert_eq!(inf_type(0.1, 0.0, 0.05), "4th");
        assert_eq!(inf_type(0.004, 0.0, 0.002), "UN");
    }

    #[test]
    fn the_second_full_sib_clause_catches_what_the_first_misses() {
        // Below 2^-1.5 = 0.35355, so clause A says "2nd"; the reference says FS.
        assert_eq!(inf_type(0.36, 0.16, 0.34), "FS");
        // ...but only above 0.32 and only with IBD2Seg over 0.15. Both bounds bracketed
        // against the reference binary.
        assert_eq!(inf_type(0.36, 0.16, 0.31), "2nd");
        assert_eq!(inf_type(0.36, 0.15, 0.34), "2nd");
    }

    #[test]
    fn the_ten_megabase_pair_filter_is_inclusive() {
        // The console says ">10Mb"; the binary reports a pair whose longest segment is
        // exactly 10 000 000 bp and drops one at 9 990 000.
        let at = |bp| PairSegments {
            ibd1_bp: bp,
            ibd2_bp: 0,
            longest_bp: bp,
        };
        assert!(at(LONG_SEGMENT_BP).reported());
        assert!(!at(LONG_SEGMENT_BP - 10_000).reported());
    }

    #[test]
    fn degree_cutoffs_are_powers_of_two() {
        assert!((degree_cutoff(1) - 0.353_553_390_6).abs() < 1e-9);
        assert!((degree_cutoff(3) - 0.088_388_347_6).abs() < 1e-9);
        // The magnitude only: a negative degree names the same band edge.
        assert_eq!(degree_cutoff(-3), degree_cutoff(3));
    }

    #[test]
    fn no_degree_reports_every_pair() {
        assert!(reported_at_degree(0, 0.0, 0.000_1));
    }

    #[test]
    fn a_positive_degree_keeps_the_band_and_everything_above_it() {
        assert!(reported_at_degree(2, 0.0, 0.1778));
        assert!(!reported_at_degree(2, 0.0, 0.1756));
        assert!(reported_at_degree(3, 0.0, 0.0893));
        assert!(!reported_at_degree(3, 0.0, 0.0866));
    }

    /// The constructed pair of `docs/research/fixtures/gate8.py`: one IBD2 block and no
    /// other sharing, so `PropIBD` equals `IBD2Seg` and the two clauses can be told apart.
    /// Reported at degree 1 by the IBD2 clause, **not** at degree 2, and again at degree 3
    /// where `PropIBD` alone carries it.
    #[test]
    fn the_first_degree_ibd2_clause_applies_at_degree_one_only() {
        assert!(reported_at_degree(1, 0.0981, 0.0981));
        assert!(!reported_at_degree(2, 0.0981, 0.0981));
        assert!(reported_at_degree(3, 0.0981, 0.0981));
        // Bracketed to (0.0789, 0.0812] by sweeping the block length.
        assert!(reported_at_degree(1, 0.0812, 0.0812));
        assert!(!reported_at_degree(1, 0.0789, 0.0789));
    }

    #[test]
    fn a_negative_degree_reports_the_complement() {
        for prop in [0.0020, 0.1000, 0.3465, 0.3600, 0.9000] {
            assert_ne!(
                reported_at_degree(2, 0.0, prop),
                reported_at_degree(-2, 0.0, prop),
                "prop {prop}"
            );
        }
    }

    // ----------------------------------------------------------------- IBD2 rule
    //
    // A two-sample fileset in which both samples are homozygous for A1 at every marker,
    // so every word is IBD2-clean, and the test then dirties named words. `WORDS` is
    // large enough that the tail rule and the ordinary one-word extension can be told
    // apart: the guard only reaches back two words from `w1`.

    const WORDS: usize = 10;

    /// Genotype codes: 0 hom A1, 1 het, 2 hom A2, 3 missing.
    fn genotypes(a: &[u8], b: &[u8]) -> Genotypes {
        let n = a.len();
        let w = n.div_ceil(WORD);
        let mut plane0 = vec![vec![0u64; w]; 2];
        let mut plane1 = vec![vec![0u64; w]; 2];
        for (s, codes) in [a, b].iter().enumerate() {
            for (m, &g) in codes.iter().enumerate() {
                let (b0, b1) = match g {
                    0 => (1, 1),
                    1 => (0, 1),
                    2 => (1, 0),
                    _ => (0, 0),
                };
                plane0[s][m / WORD] |= b0 << (m % WORD);
                plane1[s][m / WORD] |= b1 << (m % WORD);
            }
        }
        Genotypes {
            plane0,
            plane1,
            n_samples: 2,
            n_variants: n,
        }
    }

    /// The whole canvas as one usable segment, word-aligned at both ends.
    fn whole() -> Usable {
        Usable {
            chr: 1,
            lo: 0,
            hi: WORDS * WORD - 1,
        }
    }

    /// One scan over an all-hom-A1 canvas with `edits` — `(sample, marker, code)`.
    fn scan_over(edits: &[(usize, usize, u8)], seg: Usable) -> (Scan, Vec<i64>) {
        let n = WORDS * WORD;
        let mut g = [vec![0u8; n], vec![0u8; n]];
        for &(s, m, code) in edits {
            g[s][m] = code;
        }
        let gt = genotypes(&g[0], &g[1]);
        let pos: Vec<i64> = (0..n).map(|i| i as i64 * 1_000).collect();
        (Scan::new(&gt, 0, 1, seg), pos)
    }

    /// IBD2 calls over an all-hom-A1 canvas with `edits` — `(sample, marker, code)` —
    /// measured over `seg`.
    pub(super) fn ibd2_calls_over(edits: &[(usize, usize, u8)], seg: Usable) -> Vec<Called> {
        let (scan, pos) = scan_over(edits, seg);
        scan.ibd2(&pos, 0, true)
    }

    /// The `--seglength` run merge — `docs/research/20-seglength-floor.md` §2, §4.
    ///
    /// The canvas is A1A1/A1A1 everywhere, so an untouched word carries 64 A1A1/A1A1
    /// markers and no opposite homozygote; word 4 is given `z` opposite homozygotes and
    /// therefore holds `64 - z` informative ones. The budget `4 * (z - 2) <= 64 - z` puts
    /// the threshold at `z = 14`, and markers are 1 kb apart so the run-to-run gap is
    /// exactly 65 000 bp.
    #[test]
    fn seglength_merge_joins_two_runs_across_one_word() {
        let calls = |z: usize, min_bp: i64| {
            let edits: Vec<(usize, usize, u8)> =
                (0..z).map(|b| (1usize, 4 * WORD + b, 2u8)).collect();
            let (scan, pos) = scan_over(&edits, whole());
            (
                scan.ibd1(&pos, min_bp, true).len(),
                scan.ibd1(&pos, min_bp, false).len(),
            )
        };
        // Above the gap the budget decides, and it is bisected at 14 / 15.
        for z in [1usize, 2, 8, 13, 14] {
            assert_eq!(calls(z, 66_000), (1, 2), "z = {z} should merge");
        }
        for z in [15usize, 20, 64] {
            assert_eq!(calls(z, 66_000), (2, 2), "z = {z} should not merge");
        }
        // ...and the gap comparison is strict: 65 000 bp is not under 65 000.
        assert_eq!(calls(1, 65_000), (2, 2));
        assert_eq!(calls(1, 65_001), (1, 2));
        // At the default floor with real spacing the clause is live; with none, dead.
        assert_eq!(calls(1, 0), (2, 2));
    }

    /// The gate window's length bound on the IBD2 pass — `23-gap-bound.md` §2.
    ///
    /// Words 0..=3 carry two het-vs-hom mismatches each, so they are unusable and none of
    /// them can bridge; word 4 is untouched and is a one-word run; word 5 carries an
    /// opposite homozygote, which stops the right reach and keeps `ge_of(4) = 4`. The
    /// call therefore reaches back over the whole of word 3 and measures 189 000 bp while
    /// its window is one word, 63 000. Markers are 1 kb apart, so the bound cuts at
    /// `2 * 63 000 + 1`, far under the call's own length.
    #[test]
    fn an_ibd2_call_is_dropped_when_its_gate_window_is_under_half_the_floor() {
        let mut edits: Vec<(usize, usize, u8)> = Vec::new();
        for w in 0..4 {
            edits.push((1, w * WORD, 1));
            edits.push((1, w * WORD + 1, 1));
        }
        edits.push((1, 5 * WORD, 2));
        let (scan, pos) = scan_over(&edits, whole());
        let call = |min_bp| scan.ibd2(&pos, min_bp, false);

        // The call itself clears both floors: 130 .. 319 is 189 000 bp.
        let kept = call(126_001);
        assert_eq!(kept.first(), Some(&Called { lo: 130, hi: 319 }));
        assert_eq!(pos[319] - pos[130], 189_000);
        assert_eq!(pos[319] - pos[256], 63_000); // ...and its window is 63 000
        assert_eq!(kept.len(), 2);

        // One base pair of `--seglength` more and that window is under half the floor.
        assert_eq!(call(126_002), vec![Called { lo: 384, hi: 639 }]);

        // §6: the bound is asked at emit, after the merge. With the merge live the two
        // runs join across word 5 and the merged window — 383 000 bp — passes, so the
        // very call the bound refuses on its own comes back as one segment.
        assert_eq!(
            scan.ibd2(&pos, 126_002, true),
            vec![Called { lo: 130, hi: 639 }]
        );
    }

    /// The same bound on the IBD1 pass, over the run's own words and with the **strict**
    /// comparison — one unit of `min_bp / 2` tighter than the IBD2 pass's
    /// (`23-gap-bound.md` §4). Word 3 and word 5 carry one opposite homozygote each, so
    /// word 4 is a one-word run whose call reaches over both flanks; the merge is off so
    /// the three runs stay separate.
    #[test]
    fn the_ibd1_window_bound_is_strict_where_the_ibd2_one_is_not() {
        let edits = [(1usize, 3 * WORD + 1, 2u8), (1, 5 * WORD + 62, 2)];
        let (scan, pos) = scan_over(&edits, whole());
        let call = |min_bp| scan.ibd1(&pos, min_bp, false);

        let kept = call(125_999);
        assert!(kept.contains(&Called { lo: 194, hi: 382 }));
        assert_eq!(pos[382] - pos[194], 188_000);
        assert_eq!(pos[319] - pos[256], 63_000); // the run's own window

        // `63 000 > 126 000 / 2` is false, where the IBD2 pass keeps the equal case.
        assert!(!call(126_000).contains(&Called { lo: 194, hi: 382 }));
    }

    /// The IBD1 merge budget is summed over **every** word between the two runs, a
    /// gate-refused run included — `23-gap-bound.md` §5, which is `20-…` §11 item 4.
    ///
    /// Words 3 and 5 are the interruption: two opposite homozygotes and four het-vs-A1A1
    /// markers each, the rest A2A2/A2A2 and so worth nothing. Word 4 has no opposite
    /// homozygote, so it is a run of its own, and the gate refuses it while it carries
    /// under [`MIN_INFORMATIVE`] informative markers. `bad = 4` and `4 * (4 - 2) = 8`,
    /// against `X = V if V >= 10 else U` with `U = 0`: over the unusable words alone
    /// `V = 8` and nothing merges, and word 4's own markers are what take it to 10.
    #[test]
    fn a_gate_refused_run_pays_into_the_ibd1_merge_budget() {
        let calls = |v: usize| {
            let mut edits: Vec<(usize, usize, u8)> = Vec::new();
            for w in [3usize, 5] {
                edits.push((1, w * WORD, 2));
                edits.push((1, w * WORD + 1, 2));
                for b in 2..6 {
                    edits.push((1, w * WORD + b, 1));
                }
                for b in 6..WORD {
                    edits.push((0, w * WORD + b, 2));
                    edits.push((1, w * WORD + b, 2));
                }
            }
            for b in 0..v {
                edits.push((1, 4 * WORD + b, 1));
            }
            for b in v..WORD {
                edits.push((0, 4 * WORD + b, 2));
                edits.push((1, 4 * WORD + b, 2));
            }
            let (scan, pos) = scan_over(&edits, whole());
            scan.ibd1(&pos, 200_000, true)
        };
        // V = 8 + v. Under ten the budget has nothing to spend and the runs stay apart,
        // so only the right-hand one clears the floor.
        for v in [0usize, 1] {
            assert_eq!(calls(v), vec![Called { lo: 322, hi: 639 }], "v = {v}");
        }
        // At ten `X` switches to those markers and the whole canvas is one call.
        for v in [2usize, 4, 9] {
            assert_eq!(calls(v), vec![Called { lo: 0, hi: 639 }], "v = {v}");
        }
    }

    /// The same canvas read by the IBD1 pass, so the two fringe predicates can be
    /// compared side by side (`docs/research/19-ibd2seg-residual.md` §5).
    pub(super) fn ibd1_calls_over(edits: &[(usize, usize, u8)], seg: Usable) -> Vec<Called> {
        let (scan, pos) = scan_over(edits, seg);
        scan.ibd1(&pos, 0, true)
    }

    /// IBD2 calls for a pair that is identical everywhere except at `edits`, which set
    /// the second sample's genotype at a marker (1 = het, so an IBS1; 2 = hom A2, an IBS0).
    fn ibd2_calls(edits: &[(usize, u8)]) -> Vec<Called> {
        let e: Vec<(usize, usize, u8)> = edits.iter().map(|&(m, g)| (1, m, g)).collect();
        ibd2_calls_over(&e, whole())
    }

    /// `k` het-vs-hom disagreements inside word `w`.
    fn het(w: usize, k: usize) -> Vec<(usize, u8)> {
        (0..k).map(|i| (WORD * w + i, 1u8)).collect()
    }

    /// Het-vs-hom mismatches at *named bits* of word `w`, so the reach can be swept.
    fn mismatch_bits(w: usize, bits: &[usize]) -> Vec<(usize, usize, u8)> {
        bits.iter().map(|&b| (1, WORD * w + b, 1u8)).collect()
    }

    /// Every marker of word `w` an opposite homozygote — a wall no run crosses and no
    /// reach passes.
    fn ibs0_wall(w: usize) -> Vec<(usize, usize, u8)> {
        (0..WORD).map(|i| (1, WORD * w + i, 2u8)).collect()
    }

    /// A `z` word (`docs/research/17-seg-caller.md` §4): both samples homozygous for A2,
    /// so the word is usable and mismatch-free but carries no `inf2` at all.
    fn quiet_word(w: usize) -> Vec<(usize, usize, u8)> {
        (0..WORD)
            .flat_map(|i| [(0, WORD * w + i, 2u8), (1, WORD * w + i, 2u8)])
            .collect()
    }

    #[test]
    fn a_lone_dirty_word_is_absorbed_when_the_run_picks_up_cleanly() {
        // Word 4 is unusable and word 5 carries no mismatch at all, so the bridge holds
        // whatever word 4 carries — 40 mismatches or all 64.
        for k in [40, 64] {
            assert_eq!(
                ibd2_calls(&het(4, k)),
                vec![Called {
                    lo: 0,
                    hi: WORDS * WORD - 1
                }],
                "{k} mismatches in a lone word"
            );
        }
    }

    /// The conditional half of the bridge (`…/17-seg-caller.md` §3, §7): word 5 now
    /// carries a single mismatch, which leaves it *usable* but not mismatch-free, and the
    /// lone unusable word 4 is no longer absorbed. `Cyz` against `Cyx` in the doc's
    /// notation.
    #[test]
    fn a_lone_dirty_word_is_not_absorbed_when_the_next_word_carries_a_mismatch() {
        let mut e = het(4, 40);
        e.extend(het(5, 1));
        assert_eq!(
            ibd2_calls(&e),
            vec![
                Called { lo: 0, hi: 319 },
                Called {
                    lo: WORD * 7,
                    hi: WORDS * WORD - 1
                },
            ]
        );
    }

    /// `j` markers of word `w` A1A1 in both samples and the rest A2A2 in both — a
    /// mismatch-free word worth exactly `j` to the gate. `Q(j)` in `…/17-seg-caller.md`
    /// §14's notation.
    fn gate_word(w: usize, j: usize) -> Vec<(usize, usize, u8)> {
        (j..WORD)
            .flat_map(|i| [(0, WORD * w + i, 2u8), (1, WORD * w + i, 2u8)])
            .collect()
    }

    /// Words `ws` walled off with opposite homozygotes everywhere else.
    fn seg_walled(ws: &[usize], body: Vec<(usize, usize, u8)>) -> Vec<Called> {
        let mut e: Vec<(usize, usize, u8)> = (0..WORDS)
            .filter(|w| !ws.contains(w))
            .flat_map(ibs0_wall)
            .collect();
        e.extend(body);
        ibd2_calls_over(&e, whole())
    }

    /// **The bridge's left half faces the gate** (`…/17-seg-caller.md` §14.1). `[Q(j), d,
    /// C, C]`: nine informative markers before the lone unusable word and the run does not
    /// reach across it, ten and it does. Bisected on the reference at exactly 9/10 — the
    /// ordinary gate's own constant, which is what makes the bridge carry none of its own.
    #[test]
    fn the_bridge_needs_the_gate_on_the_run_it_continues() {
        // Word 2 is `Q(j)`, word 3 is unusable and worth nothing, words 4 and 5 are clean.
        let body = |j: usize| {
            let mut e = gate_word(2, j);
            e.extend(gate_word(3, 0));
            e.extend(mismatch_bits(3, &[0, 1]));
            seg_walled(&[2, 3, 4, 5], e)
        };
        // Ten: one call over words 2..=5, opening on word 2 (the wall before it blocks the
        // reach whole-word) and ending on word 5's last marker.
        assert_eq!(body(10), vec![Called { lo: 128, hi: 383 }]);
        // Nine: the run stops, is refused by the gate, and only the continuation is
        // called — which opens 63 markers before word 3's *last* mismatch, at 130.
        assert_eq!(body(9), vec![Called { lo: 130, hi: 383 }]);
    }

    /// **...and so does the continuation** (`…/17-seg-caller.md` §14.2). `[C, y, Q(j)]`:
    /// the same 9/10 bisection on the other side, with the left half held passing.
    #[test]
    fn the_bridge_needs_the_gate_on_the_continuation_too() {
        let body = |j: usize| {
            let mut e = mismatch_bits(3, &[0, 1]);
            e.extend(gate_word(3, 0));
            e.extend(mismatch_bits(3, &[0, 1]));
            e.extend(gate_word(4, j));
            seg_walled(&[2, 3, 4], e)
        };
        // Ten: the bridge holds and the call runs to word 4's last marker.
        assert_eq!(body(10), vec![Called { lo: 128, hi: 319 }]);
        // Nine: word 2 is called on its own, reaching 63 markers past word 3's *first*
        // mismatch, and the continuation is refused by the gate.
        assert_eq!(body(9), vec![Called { lo: 128, hi: 255 }]);
    }

    /// **A het-vs-A1A1 marker is not informative for the `.seg` gate**
    /// (`…/17-seg-caller.md` §14.3). It is `p1 & p1` — both samples carry A1 — so the
    /// retired statistic counted it; the reference does not. Word 2 is the quiet
    /// gate-start and word 3 carries nine HetHet plus one more informative-looking marker:
    /// ten HetHet clears the gate, nine does not, and nine plus a het-vs-A1A1 does not
    /// either, though it is ten under `p1 & p1`.
    #[test]
    fn a_het_against_a_homozygote_for_a1_is_not_informative() {
        let hethet = |w: usize, k: usize| -> Vec<(usize, usize, u8)> {
            (0..k)
                .flat_map(|i| [(0, WORD * w + i, 1u8), (1, WORD * w + i, 1u8)])
                .collect()
        };
        let body = |k: usize, ibs1b: bool| {
            let mut e = quiet_word(2);
            e.extend(gate_word(3, 0));
            e.extend(hethet(3, k));
            if ibs1b {
                e.extend(mismatch_bits(3, &[k]));
            }
            seg_walled(&[2, 3], e)
        };
        assert_eq!(body(10, false), vec![Called { lo: 128, hi: 255 }]);
        assert_eq!(body(9, false), vec![]);
        assert_eq!(body(9, true), vec![]);
    }

    #[test]
    fn two_consecutive_dirty_words_break_an_ibd2_run() {
        let mut e = het(4, 40);
        e.extend(het(5, 40));
        // The first call reaches 63 markers past word 4's *first* mismatch (marker 256),
        // so it ends at 319 — deep inside the word that ended it, not on a word boundary.
        // The second run opens on word 6 and is pushed one word: its gate-start word is 6,
        // so it starts at 64 * 7 rather than the 296 the left reach alone would give.
        assert_eq!(
            ibd2_calls(&e),
            vec![
                Called { lo: 0, hi: 319 },
                Called {
                    lo: WORD * 7,
                    hi: WORDS * WORD - 1
                },
            ]
        );
    }

    /// [`IBD2_HET_DIRTY`] bisected: one mismatch per word never splits a run, two always
    /// do (`…/17-seg-caller.md` §3). This is the constant that used to be 5, inverted out
    /// of a different caller's column.
    #[test]
    fn the_het_mismatch_threshold_is_two_per_word() {
        let mut one = het(4, 1);
        one.extend(het(5, 1));
        assert_eq!(
            ibd2_calls(&one).len(),
            1,
            "one mismatch leaves a word usable"
        );
        let mut two = het(4, 2);
        two.extend(het(5, 2));
        assert_eq!(ibd2_calls(&two).len(), 2, "two make it unusable");
    }

    /// A call reaches [`IBD2_REACH`] markers past the nearest het-vs-hom mismatch, so
    /// moving that mismatch one bit moves the call's end one marker
    /// (`…/17-seg-caller.md` §5). Words 4 and 5 are both unusable, so the run is [0, 3]
    /// and its right end is `64 * 4 + firstbit + 63`.
    #[test]
    fn a_call_reaches_sixtythree_markers_past_the_bounding_mismatch() {
        for (bits, want) in [
            (vec![0, 1], 256 + IBD2_REACH),
            (vec![30, 31], 256 + 30 + IBD2_REACH),
            (vec![62, 63], 256 + 62 + IBD2_REACH),
        ] {
            let mut e = mismatch_bits(4, &bits);
            e.extend(mismatch_bits(5, &bits));
            assert_eq!(
                ibd2_calls_over(&e, whole())[0],
                Called { lo: 0, hi: want },
                "mismatches at bits {bits:?}"
            );
        }
    }

    #[test]
    fn an_opposite_homozygote_is_never_bridged() {
        // One IBS0 in word 4, with usable words either side. It is not absorbed, and it
        // also blocks the reach **whole-word**: the first call stops on the run's own last
        // marker (255) rather than reaching into word 4, whatever bit the IBS0 sits at.
        // The second run opens on word 5 and is pushed one word, to 64 * 6.
        let calls = ibd2_calls(&[(WORD * 4 + 10, 2)]);
        assert_eq!(
            calls,
            vec![
                Called {
                    lo: 0,
                    hi: WORD * 4 - 1
                },
                Called {
                    lo: WORD * 6,
                    hi: WORDS * WORD - 1
                },
            ]
        );
    }

    /// The gate counts `inf2` from the run's first **mismatch-free** word, not from the
    /// run's first word — which is what makes `zx` a call and `xz` not, at the same two
    /// words, the same run and the same total (`…/17-seg-caller.md` §4).
    #[test]
    fn the_gate_counts_from_the_runs_first_mismatch_free_word() {
        let walls: Vec<(usize, usize, u8)> = [0usize, 1, 2, 5, 6, 7, 8, 9]
            .iter()
            .flat_map(|&w| ibs0_wall(w))
            .collect();
        let body = |e: Vec<(usize, usize, u8)>| {
            let mut all = walls.clone();
            all.extend(e);
            ibd2_calls_over(&all, whole())
        };
        let called = vec![Called {
            lo: WORD * 3,
            hi: WORD * 5 - 1,
        }];
        // `zx`: the count opens on the quiet word and picks up the mismatch word's 64.
        let mut zx = quiet_word(3);
        zx.extend(mismatch_bits(4, &[0]));
        assert_eq!(body(zx), called);
        // `xz`: the same two words the other way round. The count opens on word 4, which
        // carries nothing, and the run is refused.
        let mut xz = mismatch_bits(3, &[0]);
        xz.extend(quiet_word(4));
        assert_eq!(body(xz), vec![]);
        // ...and two quiet words are refused on their own, so it is the `x` that pays.
        let mut zz = quiet_word(3);
        zz.extend(quiet_word(4));
        assert_eq!(body(zz), vec![]);
    }

    /// A call touching the usable segment's own first or last complete word runs on to the
    /// segment's first or last **marker** — the fringe the word grid does not cover, and
    /// what makes `dups`' duplicate pair read `IBD2Seg 1.0000` rather than 0.8984.
    #[test]
    fn a_call_touching_the_segments_ends_takes_its_own_first_and_last_marker() {
        let seg = Usable {
            chr: 1,
            lo: 5,
            hi: 630,
        };
        assert_eq!(ibd2_calls_over(&[], seg), vec![Called { lo: 5, hi: 630 }]);
    }

    /// **The tail snap is gone.** The rule this replaced took a run stopping within two
    /// words of `w1` all the way to the segment's end; the `.seg`-native canvas says a run
    /// ends where its reach ends and nowhere else (`…/17-seg-caller.md` §7, §8 — nothing is
    /// ever extended or cut to a boundary).
    #[test]
    fn a_run_stopping_short_of_the_segments_last_word_does_not_snap_to_it() {
        let mut e = het(8, 40);
        e.extend(het(9, 40));
        assert_eq!(ibd2_calls(&e), vec![Called { lo: 0, hi: 575 }]);
    }

    /// Three one- and two-word runs, and the push in three places at once: without it the
    /// second call would open at 104 (word 2's last mismatch, less the reach) and the third
    /// at 360, where the reference opens them a whole word later — at `64 * 4` and `64 * 8`,
    /// one past each run's gate-start word (`…/17-seg-caller.md` §6).
    #[test]
    fn an_ibd2_run_of_one_word_is_a_segment_and_later_calls_start_one_word_late() {
        let mut e = het(1, 40);
        e.extend(het(2, 40));
        e.extend(het(5, 40));
        e.extend(het(6, 40));
        assert_eq!(
            ibd2_calls(&e),
            vec![
                Called { lo: 0, hi: 127 },
                Called {
                    lo: WORD * 4,
                    hi: 383
                },
                Called {
                    lo: WORD * 8,
                    hi: WORDS * WORD - 1
                },
            ]
        );
    }

    // -----------------------------------------------------------------------
    // The `--ibs` pass: `Scan::ibd2_words`
    // -----------------------------------------------------------------------

    /// Word calls for a pair that is heterozygous at every marker — 64 HetHet per word —
    /// with `edits` overriding one marker of one sample (`(sample, marker, code)`).
    fn word_calls(edits: &[(usize, usize, u8)]) -> Vec<WordCall> {
        let n = WORDS * WORD;
        let mut g = [vec![1u8; n], vec![1u8; n]];
        for &(s, m, code) in edits {
            g[s][m] = code;
        }
        let gt = genotypes(&g[0], &g[1]);
        let seg = Usable {
            chr: 1,
            lo: 0,
            hi: n - 1,
        };
        Scan::new(&gt, 0, 1, seg).ibd2_words()
    }

    /// Every marker of words `ws` an opposite homozygote: no heterozygote in sight.
    fn ibs0_words(ws: &[usize]) -> Vec<(usize, usize, u8)> {
        let mut e = Vec::new();
        for &w in ws {
            for i in 0..WORD {
                e.push((0, WORD * w + i, 0u8));
                e.push((1, WORD * w + i, 2u8));
            }
        }
        e
    }

    /// `k` het-vs-hom mismatches in word `w`, on the second sample.
    fn mismatches(w: usize, k: usize) -> Vec<(usize, usize, u8)> {
        (0..k).map(|i| (1, WORD * w + i, 0u8)).collect()
    }

    /// Two whole words of opposite homozygotes do not break an `--ibs` IBD2 run — the
    /// rule that separates this pass from [`Scan::ibd2`], measured on a fixture whose
    /// `--ibs` `MaxIBD2` spans the entire usable segment while `--ibdseg`'s `IBD2Seg`
    /// stops at the block.
    #[test]
    fn opposite_homozygotes_do_not_break_an_ibs_ibd2_run() {
        let calls = word_calls(&ibs0_words(&[4, 5]));
        assert_eq!(calls, vec![WordCall { lo: 0, hi: 9 }]);
    }

    /// Four het-vs-hom mismatches leave a word usable; five split the run, and the call
    /// reaches one word past it.
    #[test]
    fn five_het_mismatches_break_an_ibs_ibd2_run_and_four_do_not() {
        let mut four = mismatches(5, 4);
        four.extend(mismatches(6, 4));
        assert_eq!(word_calls(&four), vec![WordCall { lo: 0, hi: 9 }]);

        let mut five = mismatches(5, 5);
        five.extend(mismatches(6, 5));
        assert_eq!(
            word_calls(&five),
            vec![WordCall { lo: 0, hi: 5 }, WordCall { lo: 7, hi: 9 }]
        );
    }

    /// An interior call needs 95 HetHet markers over the words it is measured across;
    /// a call against the segment's own tail needs none.
    #[test]
    fn the_hethet_count_is_waived_only_at_the_segments_tail() {
        // Words 0..2 carry `h` HetHet markers each and are otherwise homozygous-concordant;
        // words 3 and 4 are dirty, so the first call is [0, 3] and its window is those
        // four words. Words 5..9 have no HetHet at all and reach the segment's last word.
        let build = |h: usize| {
            let mut e = Vec::new();
            for w in 0..WORDS {
                let keep = if w < 3 { h } else { 0 };
                for i in keep..WORD {
                    e.push((0, WORD * w + i, 0u8));
                    e.push((1, WORD * w + i, 0u8));
                }
            }
            for w in [3, 4] {
                for i in 0..5 {
                    e.push((0, WORD * w + i, 1u8));
                    e.push((1, WORD * w + i, 0u8));
                }
            }
            e
        };
        let tail = WordCall { lo: 5, hi: 9 };
        assert_eq!(
            word_calls(&build(31)),
            vec![tail],
            "93 HetHet is not enough"
        );
        assert_eq!(
            word_calls(&build(32)),
            vec![WordCall { lo: 0, hi: 3 }, tail],
            "96 HetHet is"
        );
    }

    /// A two-word interval is refused however informative it is; three words are kept.
    #[test]
    fn an_ibs_ibd2_call_must_span_three_words() {
        let mut e = mismatches(6, 5);
        e.extend(mismatches(7, 5));
        // Words 8..9 are clean and reach `w1`, so the interval is [8, 9] — two words.
        assert_eq!(word_calls(&e), vec![WordCall { lo: 0, hi: 6 }]);
    }

    // -----------------------------------------------------------------------
    // The chunk scan (`docs/research/16-segment-extension.md`)
    //
    // A second canvas, painted one *whole word* at a time from an explicit composition,
    // because everything the chunk rule does is a function of the per-word pair
    // (mismatches, HetHet) and nothing else. This is the rig `segfit.py` drives against
    // the reference binary, reproduced in-process.
    // -----------------------------------------------------------------------

    /// Word calls for a canvas in which word `w` carries `comp[w].0` het-vs-hom mismatches
    /// and `comp[w].1` HetHet markers and is homozygous-concordant everywhere else.
    fn composed(comp: &[(usize, usize)]) -> Vec<WordCall> {
        let n = comp.len() * WORD;
        let (mut a, mut b) = (vec![0u8; n], vec![0u8; n]);
        for (w, &(m, h)) in comp.iter().enumerate() {
            assert!(m + h <= WORD, "a word holds 64 markers");
            for i in 0..m {
                a[WORD * w + i] = 1; // het against hom A1: one het-vs-hom mismatch
            }
            for i in m..m + h {
                a[WORD * w + i] = 1;
                b[WORD * w + i] = 1; // both heterozygous: one HetHet
            }
        }
        let gt = genotypes(&a, &b);
        let seg = Usable {
            chr: 1,
            lo: 0,
            hi: n - 1,
        };
        Scan::new(&gt, 0, 1, seg).ibd2_words()
    }

    /// `body` walled off at both ends by two words of 64 mismatches each — two, so the
    /// bridging rule cannot absorb them — with three clean words after it so that `body`
    /// is never inside the usable segment's own exempt tail. `body` starts at word 2, and
    /// the trailing clean run is always reported as the segment's last three words.
    fn walled(body: &[(usize, usize)]) -> Vec<WordCall> {
        const WALL: (usize, usize) = (WORD, 0);
        let mut comp = vec![WALL, WALL];
        comp.extend_from_slice(body);
        comp.extend([WALL, WALL, (0, WORD), (0, WORD), (0, WORD)]);
        composed(&comp)
    }

    /// The reported interval stops at the last **confirmed chunk**, not at the run's end —
    /// and where that is depends on words the interval does not reach, which is what no
    /// word-local rule can produce (`…/16-segment-extension.md` §3, §4).
    #[test]
    fn the_ibs_ibd2_interval_is_cut_at_the_last_confirmed_chunk() {
        // Nine words at one mismatch and 19 HetHet each. Words 2..6 accumulate five
        // mismatches carrying 5 x 19 = 95 HetHet, which confirms them; words 7..10 only
        // reach four, and the wall that finally closes their chunk brings no HetHet at
        // all, so that chunk is refused and the interval stops at the confirmation plus
        // its one-mismatch overhang.
        let tail = WordCall { lo: 13, hi: 15 };
        assert_eq!(walled(&[(1, 19); 9]), vec![WordCall { lo: 2, hi: 7 }, tail]);
        // One further identical word, and the same nine words are no longer cut: the
        // tenth mismatch now closes a second chunk inside the block, and 5 x 19 confirms
        // it. The staircase of §4 in two lines.
        assert_eq!(
            walled(&[(1, 19); 10]),
            vec![WordCall { lo: 2, hi: 11 }, WordCall { lo: 14, hi: 16 }]
        );
    }

    /// [`IBS_IBD2_HETHET`] bisected against an established run: 5 x 19 confirms, 5 x 18
    /// does not, and an unconfirmed block is not shortened but dropped entirely.
    #[test]
    fn a_chunk_of_five_mismatches_needs_ninetyfive_hethet() {
        let tail = WordCall { lo: 13, hi: 15 };
        assert_eq!(
            walled(&[(1, 18); 9]),
            vec![tail],
            "5 x 18 = 90 confirms nothing"
        );
        assert_eq!(
            walled(&[(1, 19); 9]),
            vec![WordCall { lo: 2, hi: 7 }, tail],
            "5 x 19 = 95 confirms the first chunk"
        );
    }

    /// A chunk must also span [`IBS_IBD2_CHUNK_WORDS`] words: HetHet alone does not buy
    /// one. Two words carrying 122 HetHet between them are refused where three carrying
    /// 144 are confirmed, so the floor is on words and not on the count.
    #[test]
    fn a_chunk_must_span_three_words_however_much_hethet_it_carries() {
        let tail = WordCall { lo: 13, hi: 15 };
        // (3, 61): five mismatches accumulate in two words, at 122 HetHet — refused.
        assert_eq!(walled(&[(3, 61); 9]), vec![tail]);
        // (2, 48): five accumulate in three words, at 144 HetHet — confirmed, three times
        // over, and the whole block is reported.
        assert_eq!(
            walled(&[(2, 48); 9]),
            vec![WordCall { lo: 2, hi: 10 }, tail]
        );
    }

    /// The counters are **reset** at each chunk, not carried: a 20-word all-HetHet prefix
    /// buys the trailing words exactly five mismatches and never ten, at any tail length.
    /// This is the measurement that rules out a running score or an X-drop
    /// (`…/16-segment-extension.md` §5).
    #[test]
    fn a_hethet_prefix_buys_one_chunk_of_five_mismatches_and_no_more() {
        let cut = |j: usize| {
            let mut body = vec![(0usize, WORD); 20];
            body.extend(vec![(1usize, 0usize); j]);
            walled(&body)[0]
        };
        // 1 280 HetHet markers in front of the tail. Ten trailing mismatches or twenty,
        // the interval covers the same six trailing words (22..=27).
        assert_eq!(cut(10), WordCall { lo: 2, hi: 27 });
        assert_eq!(cut(20), WordCall { lo: 2, hi: 27 });
    }
}

/// The **fringe**: what happens in the partial word beyond a usable segment's word grid.
///
/// Every case here is a bisection off the reference on
/// `docs/research/fixtures/fringecanvas.py`, the rig that builds a usable segment which
/// does not start on a word boundary — see `docs/research/19-ibd2seg-residual.md`. The two
/// canvases the other campaigns used are word-aligned and cannot reach any of it, so these
/// were fitted to the corpus until §1-§5 measured them.
#[cfg(test)]
mod fringe_tests {
    use super::tests::*;
    use super::*;

    /// A segment opening 32 markers before word 1 and closing on word 8's last marker:
    /// head fringe = word 0's bits 32..=63, no tail fringe.
    fn head_seg() -> Usable {
        Usable {
            chr: 1,
            lo: WORD - 32,
            hi: 9 * WORD - 1,
        }
    }

    /// The mirror: word-aligned at the front, closing 32 markers into word 9.
    fn tail_seg() -> Usable {
        Usable {
            chr: 1,
            lo: WORD,
            hi: 9 * WORD + 31,
        }
    }

    /// One edit on the second sample: 1 = het (a het-vs-hom mismatch), 2 = hom A2 (IBS0).
    fn at(m: usize, code: u8) -> (usize, usize, u8) {
        (1, m, code)
    }

    #[test]
    fn a_call_runs_into_a_clean_fringe_to_the_segments_own_end() {
        // `19-…` §0.2: eight widths on each side, all `extend`.
        assert_eq!(
            ibd2_calls_over(&[], head_seg()),
            vec![Called { lo: 32, hi: 575 }]
        );
        assert_eq!(
            ibd2_calls_over(&[], tail_seg()),
            vec![Called { lo: 64, hi: 607 }]
        );
    }

    #[test]
    fn a_mismatch_in_the_head_fringe_stops_the_call_one_marker_past_it() {
        // `19-…` §1, bisected at 16 positions: `off = 1`, never 0.
        for q in [32usize, 40, 50, 63] {
            assert_eq!(
                ibd2_calls_over(&[at(q, 1)], head_seg()),
                vec![Called { lo: q + 1, hi: 575 }],
                "mismatch at head marker {q}"
            );
        }
    }

    #[test]
    fn the_head_fringe_stop_reads_the_last_mismatch() {
        // `19-…` §2: two breakers, the one nearest the grid wins.
        assert_eq!(
            ibd2_calls_over(&[at(40, 1), at(50, 1)], head_seg()),
            vec![Called { lo: 51, hi: 575 }]
        );
    }

    #[test]
    fn a_mismatch_in_the_tail_fringe_stops_the_call_one_marker_before_it() {
        for q in [576usize, 586, 600, 607] {
            assert_eq!(
                ibd2_calls_over(&[at(q, 1)], tail_seg()),
                vec![Called { lo: 64, hi: q - 1 }],
                "mismatch at tail marker {q}"
            );
        }
        // ...and it is the *first* one that stops it.
        assert_eq!(
            ibd2_calls_over(&[at(586, 1), at(600, 1)], tail_seg()),
            vec![Called { lo: 64, hi: 585 }]
        );
    }

    /// The result of `19-…` §2, and the one no symmetry argument would have produced: an
    /// opposite homozygote **anywhere in a complete word** disqualifies that word at any
    /// HetHet density, and the same marker in a *fringe* is invisible to the IBD2 pass.
    #[test]
    fn an_opposite_homozygote_in_a_fringe_does_not_stop_an_ibd2_call() {
        assert_eq!(
            ibd2_calls_over(&[at(50, 2)], head_seg()),
            vec![Called { lo: 32, hi: 575 }]
        );
        assert_eq!(
            ibd2_calls_over(&[at(586, 2)], tail_seg()),
            vec![Called { lo: 64, hi: 607 }]
        );
        // An IBS0 further out than a mismatch changes nothing: the mismatch still stops it.
        assert_eq!(
            ibd2_calls_over(&[at(40, 1), at(50, 2)], head_seg()),
            vec![Called { lo: 41, hi: 575 }]
        );
    }

    /// The fringe is the *segment's own* markers. Word 0's bits 0..=31 belong to whatever
    /// precedes this segment — on the corpus, another chromosome — and must not be read.
    #[test]
    fn markers_before_the_segment_are_not_part_of_its_fringe() {
        for m in [0usize, 20, 31] {
            assert_eq!(
                ibd2_calls_over(&[at(m, 1)], head_seg()),
                vec![Called { lo: 32, hi: 575 }],
                "mismatch at marker {m}, outside the segment"
            );
        }
    }

    /// `19-…` §4. Word 1 is unusable but IBS0-free, so the run is words 2..=8 and the
    /// reach carries its left end back over word 1 to exactly the grid's edge — where the
    /// marker scan takes over. Stating the clause as "a call whose *run* starts at `w0`"
    /// instead would stop this call 32 markers short.
    #[test]
    fn the_reach_snaps_out_to_the_fringe_when_it_lands_on_the_grid_edge() {
        let dirty = [at(WORD, 1), at(WORD + 1, 1)]; // two mismatches in word 1
        assert_eq!(
            ibd2_calls_over(&dirty, head_seg()),
            vec![Called { lo: 32, hi: 575 }]
        );
        // ...and the fringe's own breaker still applies once it has snapped there.
        let mut e = dirty.to_vec();
        e.push(at(50, 1));
        assert_eq!(
            ibd2_calls_over(&e, head_seg()),
            vec![Called { lo: 51, hi: 575 }]
        );
    }

    /// `19-…` §5: the IBD1 pass reads the same two partial words with the *other*
    /// predicate. Each pass stops at its own breaking marker and cannot see the other's.
    #[test]
    fn the_ibd1_fringe_breaks_on_ibs0_and_ignores_mismatches() {
        assert_eq!(
            ibd1_calls_over(&[at(50, 2)], head_seg()),
            vec![Called { lo: 51, hi: 575 }]
        );
        assert_eq!(
            ibd1_calls_over(&[at(50, 1)], head_seg()),
            vec![Called { lo: 32, hi: 575 }]
        );
        assert_eq!(
            ibd1_calls_over(&[at(586, 2)], tail_seg()),
            vec![Called { lo: 64, hi: 585 }]
        );
        assert_eq!(
            ibd1_calls_over(&[at(586, 1)], tail_seg()),
            vec![Called { lo: 64, hi: 607 }]
        );
        // Head: the last IBS0 wins; tail: the first.
        assert_eq!(
            ibd1_calls_over(&[at(40, 2), at(50, 2)], head_seg()),
            vec![Called { lo: 51, hi: 575 }]
        );
        assert_eq!(
            ibd1_calls_over(&[at(586, 2), at(600, 2)], tail_seg()),
            vec![Called { lo: 64, hi: 585 }]
        );
    }
}
