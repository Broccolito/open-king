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
//! * [`Scan::ibd1`] — the word rule (no error tolerance at all) and the asymmetric
//!   boundary refinement are pinned by forced-IBS0 sweeps; the experiments are named in
//!   its doc comment.
//! * [`Scan::ibd2`] — the break condition reproduces a measured two-word contingency
//!   table exactly, but its boundary placement is **fitted**.
//! * [`MIN_RUN1`] — **fitted**, and the weakest link. It is the one constant chosen to
//!   fit the corpus rather than derived from a forced experiment, and one controlled
//!   probe contradicts it.
//!
//! Measured against the captured reference `.seg` files (`tests/ibdseg_parity.rs`), the
//! caller reproduces **625 of 982 rows** with both `IBD1Seg` and `IBD2Seg` identical at
//! the printed four decimals, **972 of 982** `InfType` labels, and a mean absolute
//! `PropIBD` error of **0.0049** (worst 0.168). It reports every pair the reference
//! reports and 188 extra weak pairs. It is **not** byte-identical.

use king_io::Genotypes;

/// Markers per scan word. The whole engine is quantised to this.
pub const WORD: usize = 64;

/// Complete words a run must span before it can become a reported IBD1 segment.
///
/// **Fitted, not derived.** Two is what the corpus wants: raising the IBD1 floor from one
/// word to two takes the number of `.seg` rows that match the reference at all four
/// printed decimals from 317/982 to 624/982, and raising it to three starts losing rows
/// the reference reports (257 pairs disappear). One controlled probe disagrees with it —
/// an isolated single good word flanked by solid IBS0 *is* reported by the reference —
/// so the true rule is probably a length or density test that this approximates rather
/// than a word count. See the module header.
const MIN_RUN1: usize = 2;

/// The same floor for IBD2 runs. The corpus barely distinguishes 2 from 3 here
/// (625 against 625 exact rows), so it is set to match [`MIN_RUN1`].
///
/// Dropping the boundary extension instead of raising this floor is *not* the fix: with
/// segments clipped to their run's own words the caller loses 269 of the reference's 982
/// rows outright, which is the strongest single confirmation that the asymmetric
/// extension in [`Scan::right_end`] is real.
const MIN_RUN2: usize = 2;

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
#[derive(Clone, Copy, Debug, Default)]
struct WordDiff {
    ibs0: u64,
    ibs1: u64,
}

fn word_diff(g: &Genotypes, i: usize, j: usize, w: usize) -> WordDiff {
    let (p0i, p1i) = (g.plane0[i][w], g.plane1[i][w]);
    let (p0j, p1j) = (g.plane0[j][w], g.plane1[j][w]);
    let het_i = !p0i & p1i;
    let het_j = !p0j & p1j;
    WordDiff {
        ibs0: p0i & p0j & (p1i ^ p1j),
        ibs1: (het_i & p0j) | (p0i & het_j),
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
        for w in w0..w0 + nwords {
            let d = word_diff(g, i, j, w);
            ibs0.push(d.ibs0);
            ibs1.push(d.ibs1);
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
            seg,
            head_ibs0: head,
            tail_ibs0: tail,
        }
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
    pub fn ibd1(&self) -> Vec<Called> {
        let mut out: Vec<Called> = Vec::new();
        let n = self.nwords();
        let mut k = 0usize;
        while k < n {
            if self.ibs0_at(k) != 0 {
                k += 1;
                continue;
            }
            let k0 = k;
            while k < n && self.ibs0_at(k) == 0 {
                k += 1;
            }
            let k1 = k - 1;
            if k1 + 1 - k0 < MIN_RUN1 {
                continue;
            }
            let hi = self.right_end(k1);
            let mut lo = self.left_end(k0);
            if let Some(prev) = out.last() {
                lo = lo.max(prev.hi + 1);
            }
            if lo <= hi {
                out.push(Called { lo, hi });
            }
        }
        out
    }

    /// Right end of a run finishing at scan word `k1`.
    ///
    /// A run ended by an IBS0 reaches out to that word's **last** IBS0. A run ended for
    /// any other reason — the IBD2 pass's heterozygote break, or simply the end of the
    /// segment's complete words — stops on the word boundary instead.
    fn right_end(&self, k1: usize) -> usize {
        if k1 + 1 < self.nwords() {
            match self.ibs0_at(k1 + 1) {
                0 => self.marker(k1, 63),
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

    /// IBD2 segments: runs of IBS0-free words that are also not split by a heterozygote
    /// disagreement.
    ///
    /// **This rule is fitted, not pinned.** What is solid is the contingency table: with
    /// `a` forced het-vs-hom disagreements in word `w` and `b` in word `w+1`, an
    /// otherwise perfect duplicate pair keeps its single genome-wide IBD2 segment for
    /// every `(a, b)` with `a < 2` or `b = 0` — including `(64, 0)`, sixty-four
    /// consecutive disagreements — and loses it for every `(a, b)` with `a >= 2` and
    /// `b >= 1`. Sixty-four disagreements inside one word are tolerated; three straddling
    /// a word boundary are not. The break test below reproduces that table exactly.
    ///
    /// What is *not* pinned is where the resulting boundaries land: the same forced-IBS0
    /// probe that fixed [`Scan::ibd1`] shows IBD2 pulling back about a word further than
    /// this code does. IBD2 lengths should be read as close, not exact.
    pub fn ibd2(&self) -> Vec<Called> {
        let n = self.nwords();
        let mut out: Vec<Called> = Vec::new();
        let mut k = 0usize;
        while k < n {
            if !self.ibd2_word(k) {
                k += 1;
                continue;
            }
            let k0 = k;
            while k < n && self.ibd2_word(k) {
                k += 1;
            }
            let k1 = k - 1;
            if k1 + 1 - k0 < MIN_RUN2 {
                continue;
            }
            let hi = self.right_end(k1);
            let mut lo = self.left_end(k0);
            if let Some(prev) = out.last() {
                lo = lo.max(prev.hi + 1);
            }
            if lo <= hi {
                out.push(Called { lo, hi });
            }
        }
        out
    }

    /// Whether scan word `k` can be part of an IBD2 segment.
    ///
    /// It must be free of opposite homozygotes, and neither of the two boundaries it
    /// touches may carry a heterozygote break. Marking the *word* rather than only the
    /// boundary is what makes a parent-offspring pair report `IBD2Seg 0.0000`: PO
    /// genotypes disagree in every word, so every word is disqualified, whereas a rule
    /// that only cut at boundaries would leave a chain of single-word IBD2 segments each
    /// comfortably over the 3 Mb reporting floor.
    fn ibd2_word(&self, k: usize) -> bool {
        self.ibs0_at(k) == 0 && !self.het_break(k) && !(k > 0 && self.het_break(k - 1))
    }

    /// Whether a heterozygote break sits between scan words `k` and `k+1`.
    ///
    /// Two disagreements on the left of the boundary and one on the right. The asymmetry
    /// is real and was measured, not assumed: `(a, b) = (2, 1)` breaks and `(1, 2)` does
    /// not.
    fn het_break(&self, k: usize) -> bool {
        k + 1 < self.nwords() && self.ibs1[k].count_ones() >= 2 && self.ibs1[k + 1] != 0
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
    pub fn reported(&self) -> bool {
        self.longest_bp > LONG_SEGMENT_BP
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
        let ibd2 = keep_long(scan.ibd2(), pos, seglength_bp);
        let ibd1 = keep_long(scan.ibd1(), pos, seglength_bp);
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

/// Drop calls shorter than `min_bp`.
fn keep_long(mut v: Vec<Called>, pos: &[i64], min_bp: i64) -> Vec<Called> {
    v.retain(|c| pos[c.hi] - pos[c.lo] >= min_bp);
    v
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
/// First match wins. The cut-points are the ones KING's own emitted R script states
/// (`<prefix>_ibd1vsibd2.R`), cross-checked against 225 synthetic pairs with prescribed
/// (π1, π2). Two traps: the binary writes `Dup/MZ` where the manual says `Dup/MZTwin`,
/// and the `2nd` bucket has **no upper bound**, so a pair at π = 0.45 with π2 = 0 is
/// `2nd`, not `FS`.
pub fn inf_type(pi1: f64, pi2: f64, prop: f64) -> &'static str {
    const D1: f64 = 0.353_553_390_593_273_8; // 2^-1.5
    const D2: f64 = 0.176_776_695_296_636_9; // 2^-2.5
    const D3: f64 = 0.088_388_347_648_318_45; // 2^-3.5
    const D4: f64 = 0.044_194_173_824_159_22; // 2^-4.5
    if pi2 > 0.7 {
        "Dup/MZ"
    } else if pi1 + pi2 > 0.96 || (pi1 + pi2 > 0.9 && pi2 < 0.08) {
        "PO"
    } else if prop > D1 && pi2 >= 0.08 {
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

/// `--degree d` keeps rows with `PropIBD > 2^-(d+0.5)`; no `--degree` keeps everything.
pub fn degree_cutoff(degree: i32) -> f64 {
    2f64.powf(-(f64::from(degree) + 0.5))
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
    fn degree_cutoffs_are_powers_of_two() {
        assert!((degree_cutoff(1) - 0.353_553_390_6).abs() < 1e-9);
        assert!((degree_cutoff(3) - 0.088_388_347_6).abs() < 1e-9);
    }
}
