//! The IBD-segment pre-pass: which stretches of the map are usable, and the IBD2
//! sharing the `--ibs` writers report on top of them.
//!
//! [`usable`] and [`Segments::total_length`] are **established**: they reproduce
//! `<prefix>allsegs.txt` byte for byte on all 13 corpus datasets — every column, every
//! row, including the two datasets whose kept-chromosome set is not monotone in SNP
//! count — and the `Total length of …` console figure on all 13.
//!
//! [`Segments::ibd2`] carries no caller of its own: it delegates to
//! [`open_king_core::ibdseg::Scan::ibd2_words`] and only measures and totals what comes back —
//! see its doc comment and `docs/research/15-ibs-ibd2-rules.md`.

use std::fmt::Write as _;

use open_king_core::ibdseg::{self, Usable};
use open_king_io::{Genotypes, Variant};

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
pub(crate) const INFORMATIVE_BP: i64 = 100_000_000;

/// Shortest IBD2 call that lets a pair's `Pr_IBD2` be non-zero.
///
/// Bracketed to `(9 999 423, 10 000 187]` by sweeping the marker spacing under a fixed
/// three-word call: below it the pair prints `Pr_IBD2 0.0000` while still printing the
/// call's length in `MaxIBD2`. It gates the pair, not the call — see
/// [`Segments::ibd2`].
const LONG_SEGMENT_BP: i64 = 10_000_000;

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
    pub fn empty() -> Self {
        Segments {
            list: Vec::new(),
            positions: Vec::new(),
            names: Vec::new(),
        }
    }

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
    /// # The caller lives in `open_king_core`
    ///
    /// The calls come from [`open_king_core::ibdseg::Scan::ibd2_words`], which owns every rule
    /// and states the fixture behind each one. This function only measures them and
    /// applies the two rules that belong to the `--ibs` writer rather than to the scan:
    ///
    /// * **The ruler.** A call is measured word-aligned, `64 * lo` through
    ///   `64 * hi + 63`, over the denominator [`Segments::total_length`]. That is what
    ///   makes `dups`' MZ pair print `Pr_IBD2 0.8984` where `.seg` prints `IBD2Seg
    ///   1.0000`: 357 701 908 word-aligned base pairs over the same `D` of 398 163 465.
    /// * **The 10 Mb long-segment rule, which gates the *pair*, not the call.** If no
    ///   call reaches [`LONG_SEGMENT_BP`], `Pr_IBD2` is `0.0000` however much shorter
    ///   sharing there is; if one does, *every* call counts towards the total, the short
    ///   ones included. Both halves are fixture-measured: a lone 9 999 423 bp call prints
    ///   `MaxIBD2 9999423.000` with `Pr_IBD2 0.0000` and one of 10 000 187 bp prints
    ///   `0.0498`, while adding a 19.15 Mb call beside a 9.55 Mb one makes the total
    ///   28.7 Mb — the pair of them (`docs/research/15-ibs-ibd2-rules.md` §4).
    ///
    /// `MaxIBD2` is never gated: it reports the longest call whatever its length.
    ///
    /// # The rest, unchanged
    ///
    /// * **The kinship gate is [`IBD2_KINSHIP_GATE`]**, and un-gated pairs print `-9`/`-9`
    ///   in `.ibs0` but `0.000`/`0.0000` in `.ibs`.
    /// * `--seglength` never reaches here: `--ibs` has no such option.
    pub fn ibd2(&self, g: &Genotypes, i: usize, j: usize) -> Ibd2 {
        let mut max_bp = 0i64;
        let mut called_bp = 0i64;
        for seg in &self.list {
            let usable = Usable {
                chr: seg.chrom,
                lo: seg.first,
                hi: seg.last,
            };
            if usable.words() == 0 {
                continue;
            }
            for c in ibdseg::Scan::new(g, i, j, usable).ibd2_words() {
                let len = self.positions[SNPS_PER_WORD * c.hi + SNPS_PER_WORD - 1]
                    - self.positions[SNPS_PER_WORD * c.lo];
                max_bp = max_bp.max(len);
                called_bp += len;
            }
        }
        let total = self.total_length();
        Ibd2 {
            max_bp: max_bp as f64,
            proportion: if total == 0 || max_bp < LONG_SEGMENT_BP {
                0.0
            } else {
                called_bp as f64 / total as f64
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

    /// Two samples homozygous for the same allele everywhere are IBD2 across the board,
    /// and `--ibs` measures that on the word grid: `64·w0` through `64·w1 + 63`, not the
    /// usable segment's own first and last marker.
    #[test]
    fn a_full_ibd2_pair_is_measured_on_the_word_grid() {
        let (variants, kept) = map(&[("1", 1_000_000, 2001)]);
        let segs = usable(&variants, &kept, 23);
        let words = variants.len().div_ceil(SNPS_PER_WORD);
        let g = Genotypes {
            plane0: vec![vec![u64::MAX; words]; 2],
            plane1: vec![vec![u64::MAX; words]; 2],
            n_samples: 2,
            n_variants: variants.len(),
        };
        let seg = &segs.list[0];
        let (w0, w1) = (
            seg.first.div_ceil(SNPS_PER_WORD),
            (seg.last + 1) / SNPS_PER_WORD - 1,
        );
        let want = segs.positions[SNPS_PER_WORD * w1 + SNPS_PER_WORD - 1]
            - segs.positions[SNPS_PER_WORD * w0];
        let got = segs.ibd2(&g, 0, 1);
        assert_eq!(got.max_bp, want as f64);
        // Strictly under 1.0: `D` is measured to the segment's own ends, the call to the
        // word grid inside them.
        assert!(got.proportion < 1.0);
        assert_eq!(got.proportion, want as f64 / segs.total_length() as f64);
    }

    /// `Pr_IBD2` is gated on the pair having one call of at least [`LONG_SEGMENT_BP`],
    /// and `MaxIBD2` is not gated at all. A 200-marker map at 50 kb spacing carries a
    /// single usable segment barely 10 Mb wide, so a fully-IBD2 pair over it reports its
    /// length and a proportion of zero.
    #[test]
    fn a_call_under_ten_megabases_prints_a_length_but_no_proportion() {
        let mut variants = Vec::new();
        for k in 0..320 {
            variants.push(Variant {
                chrom: "1".into(),
                id: format!("rs{k}"),
                cm: 0.0,
                bp: 1_000_000 + 30_000 * k,
                a1: "A".into(),
                a2: "G".into(),
            });
        }
        let kept: Vec<usize> = (0..variants.len()).collect();
        let segs = usable(&variants, &kept, 23);
        let words = variants.len().div_ceil(SNPS_PER_WORD);
        let g = Genotypes {
            plane0: vec![vec![u64::MAX; words]; 2],
            plane1: vec![vec![u64::MAX; words]; 2],
            n_samples: 2,
            n_variants: variants.len(),
        };
        let got = segs.ibd2(&g, 0, 1);
        assert!(got.max_bp > 0.0 && got.max_bp < LONG_SEGMENT_BP as f64);
        assert_eq!(got.proportion, 0.0);
    }
}
