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
//! * [`Scan::ibd2`] and [`IBD2_HET_DIRTY`] — **measured against `--ibs`'s `MaxIBD2`**,
//!   which grades one exact segment length per pair rather than an aggregate. The rule
//!   reproduces **145 of the 158** corpus values, where the previous rule reached 95 and
//!   the best rule the earlier recon found reached 51.
//! * [`reported_at_degree`] — **measured**, over 38 298 differential cases plus a
//!   constructed fixture for the clause the corpus cannot reach.
//!
//! # What is still not right
//!
//! The **±1 word at a segment's ends**. Against the captured reference `.seg` files at
//! the default 3 Mb floor the caller reproduces **705 of 982 rows** with all four printed
//! columns identical, **981 of 982** `InfType` labels, and a mean absolute `PropIBD`
//! error of **0.00137** (worst 0.2109). The set of pairs reported is exactly right — 0
//! extra, 0 missing, on all ten datasets — so what is left is the length of calls that
//! are already found, and it goes both ways: of the 277 inexact rows, `IBD1Seg` is too
//! high on 139 and too low on 21, `IBD2Seg` too low on 121 and too high on 39.
//!
//! Read the per-dataset split with the caveat that four of the ten filesets report only
//! the 14 within-family pairs of one six-person nuclear family, over 5 000 to 10 000
//! markers — and in `monomorphic`'s case half of those markers are monomorphic or
//! ultra-rare. The reference's own numbers there are nowhere near the pedigree truth, so
//! they grade nothing. `bigish` 557/763, `multifam` 73/104, `threegen` 28/39, `admixed`
//! 12/16, `dups` 2/3, `unrelated` 1/1; `nuclear` 8/14, `missing` 8/14, `monomorphic`
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

/// Het-vs-hom disagreements that make a word too dirty to sit inside an IBD2 run.
///
/// **Measured against `--ibs`'s `MaxIBD2`**, which prints the length in base pairs of one
/// single IBD2 segment for 158 corpus pairs and so grades a candidate rule pass/fail per
/// segment instead of through an aggregate. Inverting those 158 numbers locates 154 of
/// them to exactly one word interval; over the 92 whose start is interior to a usable
/// segment, the first word of the segment has an IBS1 count of **at most 4** in 92 cases
/// out of 92, and the word before it has **at least 5** in 92 out of 92. The threshold is
/// not a fit, it is a gap in the data with nothing in it.
const IBD2_HET_DIRTY: u32 = 5;

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
        for w in w0..w0 + nwords {
            let d = word_diff(g, i, j, w);
            ibs0.push(d.ibs0);
            ibs1.push(d.ibs1);
            inf1.push(d.inf1);
            inf2.push(d.inf2);
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
    /// geometry — **not** [`Scan::ibd1`]'s.
    ///
    /// Everything here is measured against `--ibs`'s `MaxIBD2` column, which gives one
    /// exact segment length per pair for 158 corpus pairs; the fitted rule reproduces
    /// **145** of them, where the previous rule reproduced 51. Three things separate an
    /// IBD2 run from an IBD1 one:
    ///
    /// * **A lone dirty word does not break a run.** Two consecutive ones always do. Over
    ///   the 154 located segments, 103 contain a dirty word in their interior and only
    ///   three contain two in a row, and a constructed fixture makes the same point
    ///   without statistics: an interior word in which *all 64* markers are het-vs-hom
    ///   mismatches leaves the call untouched
    ///   (`docs/research/10-segment-rule-fixtures.md` §3). A word carrying an opposite
    ///   homozygote is never bridged.
    /// * **The run reaches one word past its last clean word, and not one marker past its
    ///   first.** Of the 92 located segments whose end is interior to a usable segment, 88
    ///   end on a word the rule calls dirty; every one of the 92 has a dirty word
    ///   immediately after it. The start is the opposite: all 92 begin on a clean word.
    /// * **The last word boundary of a usable segment never breaks a run.** Where the two
    ///   dirty words that would end a run are the segment's own last two, the call runs to
    ///   the segment's end instead — `nuclear`'s `N_C1`/`N_C4`, `multifam`'s 17/19 and
    ///   `bigish`'s 157/158 are all exactly this, and no other rule reaches them.
    ///
    /// Requiring a whole *word* to be clean rather than only cutting at boundaries is what
    /// makes a parent–offspring pair print `IBD2Seg 0.0000`: PO genotypes disagree far too
    /// often in every word for any of them to qualify.
    ///
    /// # The two length measures
    ///
    /// `--ibs` reports these same segments **word-aligned**, from `64u` to `64e+63`, while
    /// `.seg` measures them to the usable segment's own ends. That is not a contradiction
    /// and it is not two callers: a duplicate pair in `dups` prints `IBD2Seg 1.0000` and
    /// `Pr_IBD2 0.8984`, and 0.8984 is, to the last digit, the word-aligned total over
    /// that fileset's usable segments divided by the same `D` (357 701 908 / 398 163 465).
    /// This function returns the `.seg` measure.
    pub fn ibd2(&self, pos: &[i64], min_bp: i64) -> Vec<Called> {
        let n = self.nwords();
        let (w0, w1) = (self.seg.first_word(), self.seg.last_word());
        if n == 0 {
            return Vec::new();
        }
        let clean: Vec<bool> = (0..n).map(|k| !ibd2_dirty(self, k)).collect();
        // A single dirty word with a clean word on either side is absorbed; an opposite
        // homozygote is never absorbed. Read from `clean`, never from the running copy,
        // so two dirty words in a row can not chain their way in.
        let mut ok = clean.clone();
        for k in 1..n.saturating_sub(1) {
            if !clean[k] && clean[k - 1] && clean[k + 1] && self.ibs0_at(k) == 0 {
                ok[k] = true;
            }
        }

        let mut out: Vec<Called> = Vec::new();
        let mut k = 0usize;
        while k < n {
            if !ok[k] {
                k += 1;
                continue;
            }
            let k0 = k;
            while k < n && ok[k] {
                k += 1;
            }
            let k1 = k - 1;
            if k1 + 1 - k0 < MIN_RUN2 {
                continue;
            }
            if !Scan::informative(&self.inf2, k0, k1) {
                continue;
            }
            let (u, v) = (w0 + k0, w0 + k1);
            // One word past the run — except that the segment's own last boundary is
            // never a break, so a run stopping two words short of `w1` still takes `w1`.
            let e = if v + 2 >= w1 { w1 } else { v + 1 };
            // The start is the run's own first marker — an IBD2 call never reaches back
            // into the word that opened it, where an IBD1 call does. All 92 located
            // segments whose start is interior to a usable segment begin exactly on
            // `64u`, and refining the start by the flanking word's last IBS0 the way
            // [`Scan::left_end`] does moves no row of the 982-row corpus either way.
            let mut lo = if u == w0 { self.seg.lo } else { WORD * u };
            // The end does reach into the word that ended the run, and stops on that
            // word's last IBS0 if it has one — the same asymmetry [`Scan::right_end`]
            // applies to IBD1. Worth two IBD2 columns on the corpus; `MaxIBD2` cannot see
            // it, since every located segment's flanking words happen to be IBS0-free.
            let hi = if e == w1 {
                self.seg.hi
            } else {
                match self.ibs0_at(e - w0) {
                    0 => WORD * e + WORD - 1,
                    m => self.marker(e - w0, 63 - m.leading_zeros()),
                }
            };
            if let Some(prev) = out.last() {
                lo = lo.max(prev.hi + 1);
            }
            if lo <= hi && pos[hi] - pos[lo] >= min_bp {
                out.push(Called { lo, hi });
            }
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

    /// IBD2 calls for a pair that is identical everywhere except at `edits`, which set
    /// the second sample's genotype at a marker (1 = het, so an IBS1; 2 = hom A2, an IBS0).
    fn ibd2_calls(edits: &[(usize, u8)]) -> Vec<Called> {
        let n = WORDS * WORD;
        let a = vec![0u8; n];
        let mut b = vec![0u8; n];
        for &(m, g) in edits {
            b[m] = g;
        }
        let g = genotypes(&a, &b);
        let pos: Vec<i64> = (0..n).map(|i| i as i64 * 1_000).collect();
        let seg = Usable {
            chr: 1,
            lo: 0,
            hi: n - 1,
        };
        Scan::new(&g, 0, 1, seg).ibd2(&pos, 0)
    }

    /// `k` het-vs-hom disagreements inside word `w`.
    fn het(w: usize, k: usize) -> Vec<(usize, u8)> {
        (0..k).map(|i| (WORD * w + i, 1u8)).collect()
    }

    #[test]
    fn a_lone_dirty_word_does_not_break_an_ibd2_run() {
        let calls = ibd2_calls(&het(4, 40));
        assert_eq!(
            calls,
            vec![Called {
                lo: 0,
                hi: WORDS * WORD - 1
            }]
        );
    }

    #[test]
    fn two_consecutive_dirty_words_break_an_ibd2_run() {
        let mut e = het(4, 40);
        e.extend(het(5, 40));
        // The first call reaches one word past its last clean word — into word 4, the
        // word that ended it — and the second opens on the next clean word.
        assert_eq!(
            ibd2_calls(&e),
            vec![
                Called {
                    lo: 0,
                    hi: WORD * 5 - 1
                },
                Called {
                    lo: WORD * 6,
                    hi: WORDS * WORD - 1
                },
            ]
        );
    }

    #[test]
    fn the_het_mismatch_threshold_is_five_per_word() {
        // Four disagreements leave both words clean, so the run is never broken...
        let mut four = het(4, 4);
        four.extend(het(5, 4));
        assert_eq!(ibd2_calls(&four).len(), 1);
        // ...and five break it, because two dirty words now sit side by side.
        let mut five = het(4, 5);
        five.extend(het(5, 5));
        assert_eq!(ibd2_calls(&five).len(), 2);
    }

    #[test]
    fn an_opposite_homozygote_is_never_bridged() {
        // One IBS0 in word 4, with clean words either side: an IBD1 run would bridge
        // nothing here either, but the point is that the IBD2 pass must not absorb it.
        let calls = ibd2_calls(&[(WORD * 4 + 10, 2)]);
        assert_eq!(
            calls,
            vec![
                // the run ends on that word's last (and only) IBS0, not at the word end
                Called {
                    lo: 0,
                    hi: WORD * 4 + 10
                },
                Called {
                    lo: WORD * 5,
                    hi: WORDS * WORD - 1
                },
            ]
        );
    }

    #[test]
    fn a_usable_segments_last_word_boundary_never_breaks_a_run() {
        // Words 8 and 9 are the segment's last two and both dirty. The boundary between
        // them is never tested, so the call runs to the segment's end rather than
        // stopping one word past word 7.
        let mut e = het(8, 40);
        e.extend(het(9, 40));
        assert_eq!(
            ibd2_calls(&e),
            vec![Called {
                lo: 0,
                hi: WORDS * WORD - 1
            }]
        );
    }

    #[test]
    fn an_ibd2_run_of_one_word_is_a_segment() {
        // Words 1..2 and 5..6 dirty leaves a single clean word at 3 and 4 — one run.
        let mut e = het(1, 40);
        e.extend(het(2, 40));
        e.extend(het(5, 40));
        e.extend(het(6, 40));
        let calls = ibd2_calls(&e);
        assert_eq!(
            calls[0],
            Called {
                lo: 0,
                hi: WORD * 2 - 1
            }
        );
        assert_eq!(
            calls[1],
            Called {
                lo: WORD * 3,
                hi: WORD * 6 - 1
            }
        );
        assert_eq!(
            calls[2],
            Called {
                lo: WORD * 7,
                hi: WORDS * WORD - 1
            }
        );
    }
}
