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
//! # What is still not right — and it is all [`Scan::ibd2`]
//!
//! Against the captured reference `.seg` files at the default 3 Mb floor the caller
//! reproduces **709 of 982 rows** with all four printed columns identical, **all 982**
//! `InfType` labels, and a mean absolute `PropIBD` error of **0.00037** (worst 0.0089).
//! The set of pairs reported is exactly right — 0 extra, 0 missing, on all ten datasets —
//! so what is left is the *length* of calls that are already found.
//!
//! Split those rows by whether the reference reports any IBD2 at all and the residual
//! localises completely (`docs/research/14-ibd2-geometry.md` §2 drew this table for the
//! rule before this one; the second column is what moved):
//!
//! | reference row | rows | both estimate columns exact |
//! | --- | ---: | ---: |
//! | `IBD2Seg == 0.0000` | 823 | **823** (was 819) |
//! | `IBD2Seg > 0` | 159 | **2** (was 1) |
//!
//! The first line is now clean: the four exceptions were all `monomorphic`, where the
//! reference prints `IBD2Seg 0.0000` and the retired rule invented 0.42, 0.27, 0.03 and
//! 0.08 of the genome out of words carrying two to four mismatches. So [`Scan::ibd1`] and
//! its refinement are **finished**, the "is there any IBD2 here at all" question is
//! finished, and every failing `--ibdseg` parity case is failing on the *length* of the
//! 159 rows that do carry IBD2. Three notes for anyone continuing this:
//!
//! * **`.seg` exact-row count is a poor gradient for IBD2 work** — it is dominated by the
//!   823 rows IBD2 barely touches, which is why a rule that cuts mean error by a factor of
//!   3.7 and the worst row by a factor of 24 moves it only from 705 to 709. Grade with the
//!   `IBD2Seg` column (822 → 896) and the mean, as `tests/parity/fit/seg17.py` does.
//! * **`--ibs` is not a gradient either.** It used to be the sharp one, but
//!   [`Scan::ibd2_words`] made both its IBD2 columns exact, so every candidate `.seg` rule
//!   now scores an identical 158/158 on them. The grader that works is the `.seg`-native
//!   canvas of `docs/research/17-seg-caller.md` §1 — `docs/research/fixtures/segcanvas.py`,
//!   with 872 reference answers cached, so it re-runs in seconds without the binary.
//! * **What is left is one clause.** All 11 misses of the exhaustive length-≤4 battery and
//!   the whole 8 % random-canvas residual are the bridging condition of [`Scan::ibd2`],
//!   whose lookahead is guessed: it stops at the next unusable word, and the canvas `Cyzy`
//!   says it should not (`…/17-seg-caller.md` §11.1, §11.2). Two further open items are
//!   recorded there: a two-marker left-end effect on runs that open on a mismatch word, and
//!   the **100 Mb usable-total floor** (§2) — the reference refuses a fileset whose usable
//!   total is under 100 000 000 bp, bisected to the base pair, and this crate does not
//!   model it. No corpus dataset is anywhere near it.
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
//! they grade nothing. `bigish` 557/763, `multifam` 73/104, `threegen` 28/39, `admixed`
//! 12/16, `monomorphic` 12/14, `dups` 2/3, `unrelated` 1/1; `nuclear` 8/14, `missing`
//! 8/14, `sexchr` 8/14. `docs/PARITY.md` §5 carries the evidence.

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
/// * `inf2` — the same count for an IBD2 run, which drops the homozygosity clause: both
///   samples carry A1. HetHet is worth 1 to an IBD2 run and 0 to an IBD1 one, verified in
///   both directions on fixtures (`docs/research/13-informativeness-gate.md` §5); a pair
///   A2A2/A2A2 is worth 0 to either.
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
    WordDiff {
        ibs0: p0i & p0j & (p1i ^ p1j),
        ibs1: (het_i & p0j) | (p0i & het_j),
        inf1: share & (p0i | p0j),
        inf2: share,
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
        let head = if nwords == 0 || seg.lo == WORD * w0 {
            0
        } else {
            let d = word_diff(g, i, j, w0 - 1);
            let keep = seg.lo - WORD * (w0 - 1);
            d.ibs0 & !((1u64 << keep) - 1)
        };
        let tail = if nwords == 0 || seg.hi == WORD * (w1 + 1) - 1 {
            0
        } else {
            let d = word_diff(g, i, j, w1 + 1);
            let keep = seg.hi - WORD * (w1 + 1) + 1;
            d.ibs0 & mask_low(keep)
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
        }
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
    /// The word rule has no error tolerance whatsoever — a single opposite homozygote
    /// anywhere in a word breaks it. Pinned by putting one forced IBS0 at marker `j` in
    /// an otherwise perfectly IBD1 pair and sweeping `j`: every `j` costs exactly one
    /// marker interval of reported length, which is the signature of a split, and the
    /// split point moves with `j`'s word.
    ///
    /// The refinement is **asymmetric**, which is the part no amount of intuition
    /// produces:
    ///
    /// * the right end runs to the **last** IBS0 marker inside the *next* word — not to
    ///   the first, and not to the run's own last marker. A pair that is IBD1 over
    ///   markers 0..63 and opposite-homozygous from 64 onwards reports its segment out to
    ///   marker 127, deep inside all-IBS0 territory;
    /// * the left end starts one marker **after** the last IBS0 before the run.
    ///
    /// Consecutive segments are then clipped so they cannot overlap, earlier one wins.
    /// Where the run reaches the last complete word, the segment instead creeps into the
    /// trailing fringe marker by marker and stops just before the first IBS0 there.
    pub fn ibd1(&self, pos: &[i64], min_bp: i64) -> Vec<Called> {
        self.runs(|k| self.ibs0_at(k) == 0, MIN_RUN1, pos, min_bp)
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
    ) -> Vec<Called> {
        let mut out: Vec<Called> = Vec::new();
        let n = self.nwords();
        let mut k = 0usize;
        while k < n {
            if !good(k) {
                k += 1;
                continue;
            }
            let k0 = k;
            while k < n && good(k) {
                k += 1;
            }
            let k1 = k - 1;
            if k1 + 1 - k0 < min_run {
                continue;
            }
            if !Scan::informative(&self.inf1, k0, k1) {
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
    /// * **A lone unusable word carrying no IBS0 is absorbed** — but conditionally, which
    ///   is the part `docs/research/10-segment-rule-fixtures.md` §3 could not see: the word
    ///   after it must carry *no* mismatch at all, and the usable words from there on must
    ///   carry [`MIN_INFORMATIVE`] `inf2` between them. `CyC` is one call over three words
    ///   where `Cyz` is one over two and `CyzC` is one over four (§3, §7).
    /// * **The gate is [`MIN_INFORMATIVE`] markers of `inf2`, counted from the run's first
    ///   *mismatch-free* word** through the words the right end reaches into. Where the
    ///   count starts is not a detail: `zx` and `xz` are the same two words with the same
    ///   total and only `zx` is called (§4). A run with no mismatch-free word at all is
    ///   refused outright, which is why a uniform block at one mismatch per word reports
    ///   nothing however wide it is.
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
    /// * **A call touching the usable segment's own first or last complete word runs on to
    ///   the segment's first or last marker.** The canvas cannot see this (its chromosome 2
    ///   is word-aligned); the corpus can — it is worth 37 exact `IBD2Seg` values, and it is
    ///   why `dups`' duplicate pair reads `IBD2Seg 1.0000` and not the word-aligned 0.8984.
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
    /// # Accuracy, and the one clause that is fitted rather than bisected
    ///
    /// Out of sample on the captured corpus — nothing in it chose a constant — this scores
    /// **709 of 982 exact rows** against the previous rule's 705, `IBD2Seg` exact on
    /// **896 rows against 822**, mean `PropIBD` error **0.00037 against 0.00138**, worst row
    /// **0.0089 against 0.2109**, with the reported pair set still exactly right (0 extra,
    /// 0 missing). On 160 random canvases from two seeds that had no part in the fit it
    /// reproduces the reference **147 times (92 %)**, where the geometry it replaces
    /// reproduces 18 %; on the exhaustive battery of every word sequence of length ≤ 4 over
    /// `{clean, quiet, 1-mismatch, 2-mismatch}` it scores **329 of 340** against 4.
    ///
    /// The **bridging condition** is the exception and is flagged in place below: its
    /// lookahead is guessed rather than bisected, it accounts for all 11 exhaustive misses
    /// and the whole 8 % random residual, and on the corpus it buys `IBD2Seg` 896 against
    /// 874 while costing a little mean (0.00037 against 0.00033) and tail (0.0089 against
    /// 0.0071) — the signature of a clause that is *nearly* right
    /// (`docs/research/17-seg-caller.md` §11.1). Treat it as fitted, not measured.
    pub fn ibd2(&self, pos: &[i64], min_bp: i64) -> Vec<Called> {
        let n = self.nwords();
        let (w0, w1) = (self.seg.first_word(), self.seg.last_word());
        if n == 0 {
            return Vec::new();
        }
        let mis: Vec<u32> = self.ibs1.iter().map(|m| m.count_ones()).collect();
        let inf2: Vec<u32> = self.inf2.iter().map(|m| m.count_ones()).collect();
        let usable: Vec<bool> = (0..n).map(|k| !ibd2_dirty(self, k)).collect();

        // **The one fitted clause.** A lone unusable word carrying no opposite homozygote
        // is absorbed only when the run picks up cleanly on the other side: the next word
        // must carry no mismatch at all, and the usable words from there on must carry the
        // gate's worth of `inf2` between them. Read from `usable`, never from the running
        // copy, so two unusable words in a row cannot chain their way in.
        let mut ok = usable.clone();
        for k in 1..n.saturating_sub(1) {
            if usable[k] || self.ibs0_at(k) != 0 || !usable[k - 1] || !usable[k + 1] {
                continue;
            }
            if mis[k + 1] != 0 {
                continue;
            }
            let mut acc = 0u32;
            let mut t = k + 1;
            while t < n && usable[t] && acc < MIN_INFORMATIVE {
                acc += inf2[t];
                t += 1;
            }
            if acc >= MIN_INFORMATIVE {
                ok[k] = true;
            }
        }

        let mut out: Vec<Called> = Vec::new();
        // Calls that cleared the gate in *this* usable segment, whether or not they
        // survived `--seglength`. The push counts these, not the breaks between runs.
        let mut emitted = 0usize;
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
            // Right: the mirror image, past the *first* mismatch of the word that ended it.
            let mut right = WORD * v + WORD - 1;
            if b + 1 < n && self.ibs0_at(b + 1) == 0 && self.ibs1[b + 1] != 0 {
                let first = self.marker(b + 1, self.ibs1[b + 1].trailing_zeros());
                right = first + IBD2_REACH;
                if b + 2 >= n || self.ibs0_at(b + 2) != 0 {
                    right = right.min(WORD * (v + 2) - 1);
                }
            }

            // The gate, from the run's first mismatch-free word through the words the
            // right end reaches into. A run without one is refused outright.
            let Some(gs) = (a..=b).find(|&t| mis[t] == 0) else {
                continue;
            };
            let ge = if right > WORD * v + WORD - 1 {
                (n - 1).min(right / WORD - w0)
            } else {
                b
            };
            if inf2[gs..=ge].iter().sum::<u32>() < MIN_INFORMATIVE {
                continue;
            }

            if emitted > 0 {
                left = left.max(WORD * (w0 + gs + 1));
            }
            emitted += 1;
            // The segment's own fringes: a call touching `w0`/`w1` runs to the segment's
            // own first/last marker, which the word grid does not cover.
            if u == w0 {
                left = left.min(self.seg.lo);
            }
            if v == w1 {
                right = right.max(self.seg.hi);
            }
            let mut lo = left.max(self.seg.lo);
            let hi = right.min(self.seg.hi);
            // Consecutive calls may touch, but not overlap.
            if let Some(prev) = out.last() {
                lo = lo.max(prev.hi);
            }
            if lo <= hi && pos[hi] - pos[lo] >= min_bp {
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

    /// `PropIBD = IBD2Seg + IBD1Seg/2`, in full precision — never from the printed 4 dp.
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
/// same way the reference's console text implies ("not reported/utilized").
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
        let ibd2 = scan.ibd2(pos, seglength_bp);
        let ibd1 = scan.ibd1(pos, seglength_bp);
        for c in &ibd2 {
            let len = pos[c.hi] - pos[c.lo];
            acc.ibd2_bp += len;
            acc.longest_bp = acc.longest_bp.max(len);
        }
        for c in &ibd1 {
            let len = pos[c.hi] - pos[c.lo];
            acc.longest_bp = acc.longest_bp.max(len);
            // IBD1 is reported as the part of an IBD1 call that is not already IBD2:
            // a duplicate pair is IBD1 everywhere by the IBS0 rule yet reports
            // `IBD1Seg 0.0000`.
            acc.ibd1_bp += len - overlap(*c, &ibd2, pos);
        }
    }
    acc
}

/// Base pairs of `c` also covered by `others` (which are disjoint and ordered).
fn overlap(c: Called, others: &[Called], pos: &[i64]) -> i64 {
    let mut total = 0;
    for o in others {
        let lo = c.lo.max(o.lo);
        let hi = c.hi.min(o.hi);
        if lo < hi {
            total += pos[hi] - pos[lo];
        }
    }
    total
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

    /// IBD2 calls over an all-hom-A1 canvas with `edits` — `(sample, marker, code)` —
    /// measured over `seg`.
    fn ibd2_calls_over(edits: &[(usize, usize, u8)], seg: Usable) -> Vec<Called> {
        let n = WORDS * WORD;
        let mut g = [vec![0u8; n], vec![0u8; n]];
        for &(s, m, code) in edits {
            g[s][m] = code;
        }
        let gt = genotypes(&g[0], &g[1]);
        let pos: Vec<i64> = (0..n).map(|i| i as i64 * 1_000).collect();
        Scan::new(&gt, 0, 1, seg).ibd2(&pos, 0)
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
