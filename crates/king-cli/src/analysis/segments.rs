//! The IBD-segment pre-pass: which stretches of the map are usable, and the IBD2
//! sharing the `--ibs` writers report on top of them.
//!
//! Two very different confidence levels live in this file, and the split matters:
//!
//! * [`usable`] and [`Segments::total_length`] are **established**. They reproduce
//!   `<prefix>allsegs.txt` byte for byte on all 13 corpus datasets — every column, every
//!   row, including the two datasets whose kept-chromosome set is not monotone in SNP
//!   count — and the `Total length of …` console figure on all 13.
//! * [`Segments::ibd2`] is **fitted, not established**. The IBD2 calling rule is
//!   unpublished (`docs/SPEC.md` §8 item 16); what is proven here is its *shape* (see
//!   below), not its constants. Read the doc comment before trusting a number.

use std::fmt::Write as _;

use king_io::{Genotypes, Variant};

use crate::load;

/// SNPs per 64-bit word of a bit plane. Every boundary in this file is a multiple of it.
const SNPS_PER_WORD: usize = 64;

/// Largest base-pair gap that does **not** cut a chromosome into two segments.
///
/// `10 Mb / 64 = 156 250`, and the boundary is inclusive: a uniform 156 250 bp spacing
/// stays one usable segment while 156 251 yields `No informative IBD segments.` The
/// `.bim` cM column is ignored entirely — re-running the same fileset with `cM = 0` and
/// with `cM = 10 × Mb` gives identical output.
const MAX_GAP_BP: i64 = 156_250;

/// Complete 64-marker words a run needs before it counts as usable.
///
/// The words are aligned to the **global retained-marker index**, not to the run, which
/// is the whole subtlety: a chromosome's usable-word count depends on where its markers
/// happen to fall on that grid, so the rule is not monotone in SNP count. In `dups`,
/// chromosome 14 (372 markers, ~18.6 Mb) is dropped while chromosome 15 (355 markers,
/// ~17.7 Mb) is kept — 4 aligned words against 5. A "span > 10 Mb and ≥ 320 markers"
/// rule gets that backwards, and `docs/SPEC.md` §8 item 17 records it as an unexplained
/// counter-example; the alignment is the explanation.
const MIN_WORDS: usize = 5;

/// Total usable length at or above which the segment machinery is considered informative.
///
/// Below it the reference prints `Segments too short.` and drops `MaxIBD2`/`Pr_IBD2`
/// from both `.ibs` and `.ibs0`. Bracketed by a span sweep to `(99 990 000, 100 000 000]`.
const INFORMATIVE_BP: i64 = 100_000_000;

/// Smallest kinship for which a pair's IBD2 columns are computed at all.
///
/// `2^-3.5`; a designed dataset puts the gate between kinship `0.0880` (`-9`) and
/// `0.0900` (computed).
pub const IBD2_KINSHIP_GATE: f64 = 0.088_388_347_648_318_45;

/// One maximal run of markers dense enough for segment analysis.
///
/// `first` and `last` index the **retained** marker array — the same indexing the bit
/// planes use — and are inclusive.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Segment {
    /// Chromosome code, as the reference prints it in `allsegs.txt`.
    pub chrom: i64,
    pub first: usize,
    pub last: usize,
}

impl Segment {
    /// Markers in the run.
    pub fn snps(&self) -> usize {
        self.last - self.first + 1
    }

    /// The first and last **complete, grid-aligned** words the run contains.
    ///
    /// `None` when it contains none, which cannot happen for a segment that survived
    /// [`usable`] but can while classifying.
    fn word_span(&self) -> Option<(usize, usize)> {
        let first = self.first.div_ceil(SNPS_PER_WORD);
        let last_end = (self.last + 1) / SNPS_PER_WORD;
        (last_end > first).then(|| (first, last_end - 1))
    }
}

/// The usable segments of a map, plus the markers they are measured against.
pub struct Segments {
    pub list: Vec<Segment>,
    /// Base-pair position of every marker of the analysed set, in order.
    positions: Vec<i64>,
    /// Marker names, parallel to `positions`; `allsegs.txt` prints a segment's endpoints.
    names: Vec<String>,
}

impl Segments {
    /// Length of a segment in base pairs: last marker minus first marker.
    pub fn length(&self, s: &Segment) -> i64 {
        self.positions[s.last] - self.positions[s.first]
    }

    /// `D` — the sum of every usable segment's length, and the denominator of `Pr_IBD2`.
    ///
    /// This is exactly the figure the console reports as
    /// `Total length of <n> chromosomal segments … is <D/1e6> Mb.`; verified against all
    /// 13 datasets.
    pub fn total_length(&self) -> i64 {
        self.list.iter().map(|s| self.length(s)).sum()
    }

    /// Whether the map carries enough usable length for the IBD2 columns to appear.
    ///
    /// The same predicate drives the `Segments too short.` console line, so the two can
    /// never disagree.
    pub fn informative(&self) -> bool {
        self.total_length() >= INFORMATIVE_BP
    }

    /// The rows an X-chromosome segment set contributes to an `allsegs.txt` that already
    /// holds `preceding` autosomal rows.
    ///
    /// X segments are appended to the same file with the running `Segment` index
    /// continuing across them and `Chr` printing the X code, but they are **not** part of
    /// the `Total length of <n> chromosomal segments …` figure — the console reports them
    /// separately as `In addition to autosomes, …`. `sexchr`'s capture shows all three:
    /// two autosomal rows, a total of 199.9 Mb, and a third row with `Chr` 23.
    pub fn allsegs_continued(&self, preceding: usize) -> String {
        let mut s = String::new();
        for (n, seg) in self.list.iter().enumerate() {
            let start = self.positions[seg.first] as f64 / 1e6;
            let stop = self.positions[seg.last] as f64 / 1e6;
            let _ = writeln!(
                s,
                "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}",
                preceding + n + 1,
                seg.chrom,
                start,
                stop,
                stop - start,
                seg.snps(),
                self.names[seg.first],
                self.names[seg.last],
            );
        }
        s
    }

    /// `<prefix>allsegs.txt`, header included.
    ///
    /// ```text
    /// Segment→Chr→StartMB→StopMB→Length→N_SNP→StartSNP→StopSNP
    /// 1→1→1.009→44.257→43.248→866→rs1_1008530→rs1_44256507
    /// ```
    /// `Segment` is a 1-based running index across the whole genome, `Length` is exactly
    /// `StopMB − StartMB`, and the marker names are the run's own endpoints.
    pub fn allsegs(&self) -> String {
        let mut s =
            String::from("Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP\n");
        s.push_str(&self.allsegs_continued(0));
        s
    }

    /// IBD2 sharing for one pair: the longest IBD2 segment in base pairs and the
    /// proportion of `D` called IBD2.
    ///
    /// # What is established
    ///
    /// * **The reported endpoints are word-aligned.** Searching every marker pair in a
    ///   dataset for the exact `MaxIBD2` value the reference printed locates a unique
    ///   pair every time, and it is always `(64·a, 64·(b+1) − 1)` for some words `a ≤ b`
    ///   — first marker of a word to last marker of a word. Six independent hits in
    ///   `nuclear` alone (e.g. `108752454 = bp[2303] − bp[128]`, words 2..35).
    /// * **The denominator is `D`**, [`Segments::total_length`], not the genome span and
    ///   not a marker count: the marker-count form misprints every value.
    /// * **The gate is [`IBD2_KINSHIP_GATE`]**, and un-gated pairs print `-9`/`-9` in
    ///   `.ibs0` but `0.000`/`0.0000` in `.ibs`.
    ///
    /// # `--ibs` does not share `--ibdseg`'s caller, and must not be made to
    ///
    /// The two disagree **in the reference's own output**. On `nuclear`, the pair
    /// `N_C1`/`N_C2` is `Pr_IBD2 0.2173` in `king.ibs` and `IBD2Seg 0.2626` in
    /// `king.seg`; all six of that dataset's IBD2-sharing pairs differ the same way, and
    /// `.ibs` is the smaller every time (0.2749/0.3144, 0.4669/0.5095, 0.4604/0.5194,
    /// 0.2942/0.3531, 0.2812/0.3108). Same binary, same fileset, same denominator `D` —
    /// so `--ibs` runs a *different*, tighter IBD2 rule than `--ibdseg` does. Wiring this
    /// function to `king_core::ibdseg` is the obvious tidy-up, and it is wrong: it would
    /// make `Pr_IBD2` systematically too large. The duplication is the finding.
    ///
    /// # Alternatives already measured, so they need not be tried again
    ///
    /// Scored over the 152 gated, non-zero `MaxIBD2` rows of the golden `.ibs`/`.ibs0`
    /// corpus (`monomorphic` and `sexchr` are the two datasets the rule below already
    /// reproduces byte for byte, and any replacement has to keep them):
    ///
    /// | word test | run rule | right end | exact `MaxIBD2` |
    /// | --- | --- | --- | ---: |
    /// | no IBS0, ≤4 het mismatches | maximal runs | into the next word | **83** |
    /// | no IBS0, ≤3 / ≤5 / ≤8 | maximal runs | into the next word | 80 / 78 / 62 |
    /// | no IBS0, ≤4 | maximal runs | word-aligned | 46 |
    /// | `king_core::ibdseg`'s het-break test | maximal runs | word-aligned | 45 |
    /// | no IBS0, any het count | maximal runs | word-aligned | 35 |
    ///
    /// The top row scores nearly twice the rule below on `MaxIBD2` and still **fails the
    /// two datasets that currently pass**, so it is not an improvement in the only
    /// metric that counts. It also reproduces `Pr_IBD2` on zero rows, which localises
    /// the residual: the reference sometimes reports one segment where every variant
    /// here reports two (`nuclear`'s `N_C1`/`N_C4` is 108 752 454 bp against a longest
    /// run of 57 547 501), and a sum over runs feels every such split while a maximum
    /// does not. Whatever bridges those runs is the missing rule.
    ///
    /// # What is fitted and unverified
    ///
    /// The word-compatibility test and the run state machine below. A grid search over
    /// eight datasets and 171 gated pairs puts the best setting at "a word is IBD2 iff it
    /// carries no het/hom mismatch and at most one opposite-homozygote call; a run opens
    /// on three consecutive compatible words, survives up to two consecutive
    /// incompatible ones, and ends one word past its last compatible word". That
    /// reproduces `MaxIBD2` on 159/171 and `Pr_IBD2` on 143/171 — good enough to show
    /// the shape is right and **not** good enough for byte parity.
    ///
    /// Known-wrong classes, each a lead for whoever finishes this: `monomorphic` and
    /// `missing` pairs where uninformative words (all-homozygous or all-missing) look
    /// compatible and the reference still calls no segment, and long runs where the
    /// reference's boundary sits one word further out than the fitted rule puts it. Until
    /// those are explained, `--ibs` parity holds only for datasets in which no pair
    /// reaches the gate.
    pub fn ibd2(&self, g: &Genotypes, i: usize, j: usize) -> Ibd2 {
        let mut lengths: Vec<i64> = Vec::new();
        for seg in &self.list {
            let Some((first_word, last_word)) = seg.word_span() else {
                continue;
            };
            let compatible: Vec<bool> = (first_word..=last_word)
                .map(|w| word_is_ibd2(g, i, j, w))
                .collect();
            for (a, b) in runs(&compatible) {
                let (a, b) = (first_word + a, first_word + b);
                let start = self.positions[a * SNPS_PER_WORD];
                let stop = self.positions[(b + 1) * SNPS_PER_WORD - 1];
                lengths.push(stop - start);
            }
        }
        let total = self.total_length();
        Ibd2 {
            max_bp: lengths.iter().copied().max().unwrap_or(0) as f64,
            proportion: if total == 0 {
                0.0
            } else {
                lengths.iter().sum::<i64>() as f64 / total as f64
            },
        }
    }
}

/// A pair's IBD2 summary, in the units the `.ibs` columns print.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Ibd2 {
    /// `MaxIBD2` — longest IBD2 segment in base pairs, printed `%.3lf`.
    pub max_bp: f64,
    /// `Pr_IBD2` — IBD2 length over `D`, printed `%.4lf`.
    pub proportion: f64,
}

/// Words a run must open on, and consecutive incompatible words it survives.
///
/// Fitted; see [`Segments::ibd2`].
const RUN_OPEN_WORDS: usize = 3;
const RUN_BREAK_WORDS: usize = 3;

/// Whether one 64-marker word looks IBD2 for a pair.
///
/// IBD2 means the two share both alleles, so *any* het/hom mismatch or opposite
/// homozygote contradicts it; the single tolerated opposite homozygote is the fitted
/// allowance for genotyping error. Missing calls contribute to neither count, which is
/// the known weakness recorded in [`Segments::ibd2`].
fn word_is_ibd2(g: &Genotypes, i: usize, j: usize, w: usize) -> bool {
    let (x0, x1) = (g.plane0[i][w], g.plane1[i][w]);
    let (y0, y1) = (g.plane0[j][w], g.plane1[j][w]);
    let nm_i = x0 | x1;
    let nm_j = y0 | y1;
    let het_i = !x0 & x1;
    let het_j = !y0 & y1;
    let both_hom = x0 & y0;
    let mismatch = ((het_i & nm_j) | (het_j & nm_i)) & !(het_i & het_j);
    let ibs0 = both_hom & (x1 ^ y1);
    mismatch == 0 && ibs0.count_ones() <= 1
}

/// Maximal IBD2 runs over a segment's compatibility flags, as inclusive word indices.
///
/// A run opens on [`RUN_OPEN_WORDS`] consecutive compatible words, survives up to
/// `RUN_BREAK_WORDS - 1` consecutive incompatible ones, and stops one word past its last
/// compatible word — the asymmetry is what the reference's boundaries show, not a
/// simplification.
fn runs(compatible: &[bool]) -> Vec<(usize, usize)> {
    let n = compatible.len();
    let mut out = Vec::new();
    let mut start: Option<usize> = None;
    let mut last_ok = 0usize;
    let mut broken = 0usize;
    for k in 0..n {
        match start {
            None => {
                let opens = compatible[k]
                    && k + RUN_OPEN_WORDS <= n
                    && compatible[k..k + RUN_OPEN_WORDS].iter().all(|&c| c);
                if opens {
                    start = Some(k);
                    last_ok = k;
                    broken = 0;
                }
            }
            Some(a) => {
                if compatible[k] {
                    last_ok = k;
                    broken = 0;
                } else {
                    broken += 1;
                    if broken >= RUN_BREAK_WORDS {
                        out.push((a, (last_ok + 1).min(n - 1)));
                        start = None;
                    }
                }
            }
        }
    }
    if let Some(a) = start {
        out.push((a, (last_ok + 1).min(n - 1)));
    }
    out
}

/// Cut the retained map into usable segments.
///
/// Two steps, both read off the reference:
///
/// 1. **Cut** wherever the chromosome changes or consecutive markers are more than
///    [`MAX_GAP_BP`] apart.
/// 2. **Drop** every run that does not contain [`MIN_WORDS`] complete words of the
///    global 64-marker grid.
///
/// `sexchr` is needed only to resolve chromosome labels to codes, so that `01` and `1`
/// are one chromosome and the printed `Chr` column matches the reference's.
pub fn usable(variants: &[Variant], kept: &[usize], sexchr: i64) -> Segments {
    let selected: Vec<&Variant> = kept.iter().map(|&k| &variants[k]).collect();
    cut(&selected, sexchr)
}

/// The same construction over the X chromosome alone.
///
/// The reference packs X genotypes into their own word array, so the 64-marker grid runs
/// over the X markers by themselves rather than over their `.bim` positions. The result
/// never enters `Total length …` or `Pr_IBD2`; it is reported by its own console line and
/// appended to `allsegs.txt`.
pub fn usable_x(variants: &[Variant], sexchr: i64) -> Segments {
    let selected: Vec<&Variant> = variants
        .iter()
        .filter(|v| load::chromosome_code(&v.chrom, sexchr) == sexchr)
        .collect();
    cut(&selected, sexchr)
}

/// The shared cut-and-filter pass over an already-selected marker list.
fn cut(selected: &[&Variant], sexchr: i64) -> Segments {
    let positions: Vec<i64> = selected.iter().map(|v| v.bp).collect();
    let names: Vec<String> = selected.iter().map(|v| v.id.clone()).collect();
    let codes: Vec<i64> = selected
        .iter()
        .map(|v| load::chromosome_code(&v.chrom, sexchr))
        .collect();

    let mut list = Vec::new();
    let mut start = 0usize;
    for k in 1..=positions.len() {
        let cut = k == positions.len()
            || codes[k] != codes[k - 1]
            || positions[k] - positions[k - 1] > MAX_GAP_BP
            // A decreasing position is not a gap; treat it as a cut rather than let a
            // negative length into the totals. The reference refuses such a map outright
            // (`Chromosomes unsorted: …`), which is a separate, unimplemented path.
            || positions[k] < positions[k - 1];
        if !cut {
            continue;
        }
        let seg = Segment {
            chrom: codes[start],
            first: start,
            last: k - 1,
        };
        if seg.word_span().is_some_and(|(a, b)| b + 1 - a >= MIN_WORDS) {
            list.push(seg);
        }
        start = k;
    }
    Segments {
        list,
        positions,
        names,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn map(spec: &[(&str, i64, usize)]) -> (Vec<Variant>, Vec<usize>) {
        let mut variants = Vec::new();
        for (chrom, start, n) in spec {
            for k in 0..*n {
                variants.push(Variant {
                    chrom: chrom.to_string(),
                    id: format!("rs{}_{}", chrom, k),
                    cm: 0.0,
                    bp: start + 50_000 * k as i64,
                    a1: "A".into(),
                    a2: "G".into(),
                });
            }
        }
        let kept = (0..variants.len()).collect();
        (variants, kept)
    }

    /// The `dups` puzzle: chromosome 14 has *more* markers than 15 and is dropped, because
    /// the global word grid gives it one complete word fewer.
    #[test]
    fn usable_segments_count_grid_aligned_words_not_markers() {
        // 7622 markers before chr14 reproduces the captured offsets exactly.
        let (variants, kept) = map(&[
            ("1", 1_000_000, 7622),
            ("14", 1_000_000, 372),
            ("15", 1_000_000, 355),
        ]);
        let segs = usable(&variants, &kept, 23);
        let chroms: Vec<i64> = segs.list.iter().map(|s| s.chrom).collect();
        assert_eq!(chroms, [1, 15], "chr14 has 4 aligned words, chr15 has 5");
    }

    #[test]
    fn a_gap_of_exactly_the_limit_does_not_cut() {
        let mut variants = Vec::new();
        for k in 0..400 {
            variants.push(Variant {
                chrom: "1".into(),
                id: format!("rs{k}"),
                cm: 0.0,
                // One 156 250 bp gap in the middle: usable.
                bp: 1_000_000
                    + if k < 200 {
                        50_000 * k
                    } else {
                        50_000 * k + 106_250
                    },
                a1: "A".into(),
                a2: "G".into(),
            });
        }
        let kept: Vec<usize> = (0..variants.len()).collect();
        assert_eq!(usable(&variants, &kept, 23).list.len(), 1);

        // One base pair more and the chromosome splits into two 200-marker pieces, each
        // with fewer than five aligned words, so nothing survives.
        variants[200].bp += 1;
        for v in &mut variants[201..] {
            v.bp += 1;
        }
        assert!(usable(&variants, &kept, 23).list.is_empty());
    }

    #[test]
    fn allsegs_renders_the_captured_columns() {
        let (variants, kept) = map(&[("1", 1_008_530, 866)]);
        let segs = usable(&variants, &kept, 23);
        let text = segs.allsegs();
        let mut lines = text.lines();
        assert_eq!(
            lines.next().unwrap(),
            "Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP"
        );
        assert_eq!(
            lines.next().unwrap(),
            "1\t1\t1.009\t44.259\t43.250\t866\trs1_0\trs1_865"
        );
    }

    #[test]
    fn informative_at_exactly_one_hundred_megabases() {
        // 2001 markers at 50 kb spacing span exactly 100 Mb.
        let (variants, kept) = map(&[("1", 1_000_000, 2001)]);
        let segs = usable(&variants, &kept, 23);
        assert_eq!(segs.total_length(), 100_000_000);
        assert!(segs.informative());

        let (variants, kept) = map(&[("1", 1_000_000, 2000)]);
        assert!(!usable(&variants, &kept, 23).informative());
    }

    #[test]
    fn runs_open_late_and_close_one_word_past_the_last_match() {
        let f = |s: &str| -> Vec<bool> { s.bytes().map(|c| c == b'.').collect() };
        // Neither of the two lone compatible words can open a run; the third position
        // with three in a row does, and the trailing mismatch is inside the run.
        assert_eq!(runs(&f("x.x...........x")), [(3, 14)]);
        // Up to two consecutive mismatching words are bridged.
        assert_eq!(runs(&f("....xx.....")), [(0, 10)]);
        // Three end the run one word past its last match, and a new one may open after.
        assert_eq!(runs(&f("....xxx....")), [(0, 4), (7, 10)]);
        // Nothing at all when the segment never gets three consecutive matches.
        assert!(runs(&f("..x..x..x")).is_empty());
    }
}
