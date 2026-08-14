//! The pairwise counting kernel.
//!
//! This is the hot loop of the whole program: it is `O(n_pairs * n_words)` and every
//! other stage is negligible beside it. Each count is obtained as a popcount over a
//! word-wise boolean combination of the two samples' bit planes, so a pair costs one
//! pass over `2 * ceil(n_variants / 64)` words per sample with no branching.
//!
//! # Deriving the expressions from the plane encoding
//!
//! [`Genotypes`] stores, for sample `s` and variant `m`, bit `m` of `plane0[s]` and
//! `plane1[s]`:
//!
//! | `plane0` | `plane1` | genotype |
//! |----------|----------|----------|
//! | 1        | 1        | hom A1/A1 |
//! | 0        | 1        | het |
//! | 1        | 0        | hom A2/A2 |
//! | 0        | 0        | missing |
//!
//! Three masks fall straight out of that table:
//!
//! * `nonmissing = plane0 | plane1` — missing is the *only* code with both bits clear,
//!   so the OR is set exactly on called genotypes.
//! * `hom = plane0` — `plane0` is set for both homozygous codes and clear for het and
//!   for missing, so it is already an "is a homozygous call" mask. Note `hom` implies
//!   `nonmissing`.
//! * `het = !plane0 & plane1` — het is the unique code with `plane0` clear and `plane1`
//!   set. Equivalently `nonmissing & !hom`. Note `het` also implies `nonmissing`.
//!
//! Every count in [`PairCounts`] is taken over the **pairwise** non-missing set
//! `M_ij = nonmissing_i & nonmissing_j` (see `docs/VERIFIED_FORMULAS.md`), which gives:
//!
//! * `n_snp = popcount(nonmissing_i & nonmissing_j)` — the size of `M_ij` itself.
//! * `het_i = popcount(het_i_bits & nonmissing_j)`. The mask by *`j`'s* callability is
//!   the load-bearing part: `het_i_bits` already implies `nonmissing_i`, so ANDing with
//!   `nonmissing_j` restricts the count to `M_ij`. A per-sample heterozygote count
//!   computed once and reused across pairs is a different, wrong number whenever the two
//!   samples differ in missingness.
//! * `het_j = popcount(het_j_bits & nonmissing_i)` — the same by symmetry.
//! * `het_het = popcount(het_i_bits & het_j_bits)` — each operand implies its own
//!   sample's callability, so their AND already lies in `M_ij`; no extra mask is needed.
//! * `hom_hom = popcount(plane0[i] & plane0[j])` — likewise self-masking.
//! * `ibs0 = popcount(plane0[i] & plane0[j] & (plane1[i] ^ plane1[j]))`. IBS0 means
//!   *opposite* homozygotes. `plane0[i] & plane0[j]` restricts to sites where both are
//!   homozygous; among homozygotes `plane1` is 1 for A1/A1 and 0 for A2/A2, so it is a
//!   1-bit label of *which* homozygote, and the XOR selects exactly the sites where the
//!   two labels disagree. The `plane0` factors are what keep het and missing sites out
//!   of the XOR.
//!
//! All six expressions were checked black-box against KING 2.3.2's `.ibs`/`.ibs0` output
//! on generated filesets — including one whose per-sample missing rates ranged from 2%
//! to 70%, where the pairwise mask changes `N_Het1`/`N_Het2` on every pair — with zero
//! discrepancies over all six columns.
//!
//! The exact identity `hom_hom = n_snp - het_i - het_j + het_het` holds (partition `M_ij`
//! by which member is heterozygous), and the tests assert it, but `hom_hom` is counted
//! directly rather than derived so that the assertion stays an independent check.
//!
//! Bits past `n_variants` in the last word of each plane are required by the
//! [`Genotypes`] contract to be zero, so no tail masking happens here.

use king_io::Genotypes;
use rayon::prelude::*;

use crate::PairCounts;

/// Count the six pairwise quantities for samples `i` and `j`.
///
/// Both orders are supported: `pair_counts(g, j, i)` returns the same counts with
/// `het_i` and `het_j` exchanged.
///
/// # Panics
///
/// Panics if `i` or `j` is not a valid sample index, if either sample's planes are
/// shorter than [`Genotypes::words_per_sample`], or if `n_variants` exceeds `u32::MAX`
/// (in which case the counts could not be represented in [`PairCounts`]).
#[inline]
pub fn pair_counts(g: &Genotypes, i: usize, j: usize) -> PairCounts {
    assert!(
        g.n_variants <= u32::MAX as usize,
        "n_variants {} exceeds the u32 counters in PairCounts",
        g.n_variants
    );

    let w = g.words_per_sample();
    // Slicing to a common length up front is what lets the zipped loop below run
    // without per-word bounds checks, and turns a short plane into a clear panic here
    // rather than a silently truncated count.
    let a0 = &g.plane0[i][..w];
    let a1 = &g.plane1[i][..w];
    let b0 = &g.plane0[j][..w];
    let b1 = &g.plane1[j][..w];

    debug_assert!(
        tail_is_clean(a0, a1, g.n_variants),
        "sample {i} has dirty tail bits"
    );
    debug_assert!(
        tail_is_clean(b0, b1, g.n_variants),
        "sample {j} has dirty tail bits"
    );

    // Each accumulator is bounded by the number of variants, which the assertion above
    // pins below u32::MAX, so none of them can overflow.
    let mut n_snp: u32 = 0;
    let mut het_i: u32 = 0;
    let mut het_j: u32 = 0;
    let mut het_het: u32 = 0;
    let mut hom_hom: u32 = 0;
    let mut ibs0: u32 = 0;

    for (((&x0, &x1), &y0), &y1) in a0.iter().zip(a1).zip(b0).zip(b1) {
        let nm_i = x0 | x1; // called in i
        let nm_j = y0 | y1; // called in j
        let het_bits_i = !x0 & x1; // i is heterozygous
        let het_bits_j = !y0 & y1; // j is heterozygous
        let both_hom = x0 & y0; // both homozygous (implies both called)

        n_snp += (nm_i & nm_j).count_ones();
        het_i += (het_bits_i & nm_j).count_ones();
        het_j += (het_bits_j & nm_i).count_ones();
        het_het += (het_bits_i & het_bits_j).count_ones();
        hom_hom += both_hom.count_ones();
        ibs0 += (both_hom & (x1 ^ y1)).count_ones();
    }

    PairCounts {
        n_snp,
        het_i,
        het_j,
        het_het,
        ibs0,
        hom_hom,
    }
}

/// Count every pair in `pairs`, in parallel, returning results in the same order.
///
/// The output is **independent of the thread count and of how rayon happens to split
/// the work**: each element is a pure function of `g` and its own `(i, j)`, nothing is
/// accumulated across pairs, and collecting an indexed parallel iterator writes each
/// result to its own slot. There is no floating point anywhere in the kernel, so there
/// is not even an associativity question to get wrong.
pub fn all_pairs(g: &Genotypes, pairs: &[(usize, usize)]) -> Vec<PairCounts> {
    pairs
        .par_iter()
        .map(|&(i, j)| pair_counts(g, i, j))
        .collect()
}

/// Counts for a **single** sample: `(called, heterozygous)` over its own non-missing set.
///
/// Distinct from the `het_i` of [`pair_counts`], which is restricted to the *pairwise*
/// non-missing set. Both are needed: the pairwise form drives every printed estimator,
/// while this one is what the X-chromosome pass medians over to impute a heterozygosity
/// for hemizygous males (`analysis::xkinship`).
///
/// # Panics
///
/// Panics if `i` is not a valid sample index or if its planes are shorter than
/// [`Genotypes::words_per_sample`].
pub fn sample_counts(g: &Genotypes, i: usize) -> (u32, u32) {
    let w = g.words_per_sample();
    let p0 = &g.plane0[i][..w];
    let p1 = &g.plane1[i][..w];
    let mut called: u32 = 0;
    let mut het: u32 = 0;
    for (&x0, &x1) in p0.iter().zip(p1) {
        called += (x0 | x1).count_ones();
        het += (!x0 & x1).count_ones();
    }
    (called, het)
}

/// Whether the unused high bits of the last plane word are zero, as the [`Genotypes`]
/// contract requires. Used only by `debug_assert!`.
fn tail_is_clean(p0: &[u64], p1: &[u64], n_variants: usize) -> bool {
    let used = n_variants % 64;
    if used == 0 {
        return true;
    }
    let tail = !0u64 << used;
    match (p0.last(), p1.last()) {
        (Some(&x), Some(&y)) => (x | y) & tail == 0,
        _ => true,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- test-side genotype model -------------------------------------------------

    #[derive(Debug, Clone, Copy, PartialEq, Eq)]
    enum Geno {
        HomA1,
        Het,
        HomA2,
        Missing,
    }
    use Geno::*;

    /// Every genotype, in a fixed order, for exhaustive enumeration.
    const ALL: [Geno; 4] = [HomA1, Het, HomA2, Missing];

    /// Build bit planes from a `[sample][variant]` genotype matrix, honouring the
    /// documented encoding and leaving the tail bits of the last word clear.
    fn planes(rows: &[Vec<Geno>]) -> Genotypes {
        let n_variants = rows.first().map_or(0, |r| r.len());
        assert!(rows.iter().all(|r| r.len() == n_variants));
        let words = n_variants.div_ceil(64);
        let mut plane0 = vec![vec![0u64; words]; rows.len()];
        let mut plane1 = vec![vec![0u64; words]; rows.len()];
        for (s, row) in rows.iter().enumerate() {
            for (m, g) in row.iter().enumerate() {
                let (bit0, bit1) = match g {
                    HomA1 => (true, true),
                    Het => (false, true),
                    HomA2 => (true, false),
                    Missing => (false, false),
                };
                if bit0 {
                    plane0[s][m / 64] |= 1 << (m % 64);
                }
                if bit1 {
                    plane1[s][m / 64] |= 1 << (m % 64);
                }
            }
        }
        Genotypes {
            plane0,
            plane1,
            n_samples: rows.len(),
            n_variants,
        }
    }

    /// Decode one genotype back out of the planes — the inverse of [`planes`], used to
    /// keep the naive reference honest about reading the same bits the kernel reads.
    fn decode(g: &Genotypes, s: usize, m: usize) -> Geno {
        let bit = |plane: &Vec<Vec<u64>>| plane[s][m / 64] >> (m % 64) & 1 == 1;
        match (bit(&g.plane0), bit(&g.plane1)) {
            (true, true) => HomA1,
            (false, true) => Het,
            (true, false) => HomA2,
            (false, false) => Missing,
        }
    }

    /// Naive per-variant reference: decode both genotypes, count in a match. Deliberately
    /// shares no logic with the popcount kernel.
    fn naive(g: &Genotypes, i: usize, j: usize) -> PairCounts {
        let mut c = PairCounts::default();
        for m in 0..g.n_variants {
            match (decode(g, i, m), decode(g, j, m)) {
                (Missing, _) | (_, Missing) => continue,
                (Het, Het) => {
                    c.het_i += 1;
                    c.het_j += 1;
                    c.het_het += 1;
                }
                (Het, HomA1) | (Het, HomA2) => c.het_i += 1,
                (HomA1, Het) | (HomA2, Het) => c.het_j += 1,
                (HomA1, HomA1) | (HomA2, HomA2) => c.hom_hom += 1,
                (HomA1, HomA2) | (HomA2, HomA1) => {
                    c.hom_hom += 1;
                    c.ibs0 += 1;
                }
            }
            c.n_snp += 1;
        }
        c
    }

    fn counts(
        n_snp: u32,
        het_i: u32,
        het_j: u32,
        het_het: u32,
        ibs0: u32,
        hom_hom: u32,
    ) -> PairCounts {
        PairCounts {
            n_snp,
            het_i,
            het_j,
            het_het,
            ibs0,
            hom_hom,
        }
    }

    // ---- exhaustive 4x4 genotype-pair table ---------------------------------------

    /// The hand-written truth table for a single variant, one row per element of
    /// `{HomA1, Het, HomA2, Missing}^2`. Written out from the *definitions* of the
    /// counts, not from the kernel or the naive implementation.
    fn expected_single(gi: Geno, gj: Geno) -> PairCounts {
        //                                  n_snp het_i het_j hh ibs0 homhom
        match (gi, gj) {
            (HomA1, HomA1) => counts(1, 0, 0, 0, 0, 1),
            (HomA1, Het) => counts(1, 0, 1, 0, 0, 0),
            (HomA1, HomA2) => counts(1, 0, 0, 0, 1, 1),
            (Het, HomA1) => counts(1, 1, 0, 0, 0, 0),
            (Het, Het) => counts(1, 1, 1, 1, 0, 0),
            (Het, HomA2) => counts(1, 1, 0, 0, 0, 0),
            (HomA2, HomA1) => counts(1, 0, 0, 0, 1, 1),
            (HomA2, Het) => counts(1, 0, 1, 0, 0, 0),
            (HomA2, HomA2) => counts(1, 0, 0, 0, 0, 1),
            // A missing call on either side removes the site from M_ij entirely.
            (Missing, _) | (_, Missing) => counts(0, 0, 0, 0, 0, 0),
        }
    }

    #[test]
    fn exhaustive_genotype_pair_table() {
        let mut seen = 0;
        for gi in ALL {
            for gj in ALL {
                let g = planes(&[vec![gi], vec![gj]]);
                let got = pair_counts(&g, 0, 1);
                assert_eq!(
                    got,
                    expected_single(gi, gj),
                    "kernel wrong for ({gi:?}, {gj:?})"
                );
                assert_eq!(
                    naive(&g, 0, 1),
                    expected_single(gi, gj),
                    "naive wrong for ({gi:?}, {gj:?})"
                );
                seen += 1;
            }
        }
        assert_eq!(
            seen, 16,
            "the enumeration must cover all 4x4 genotype pairs"
        );
    }

    #[test]
    fn exhaustive_table_summed_in_one_pass() {
        // All 16 combinations as 16 variants of a single pair. The expected totals are
        // the 3x3 non-missing block: 9 called sites, 3 het rows, 3 het columns, 1
        // het/het cell, 2 opposite-homozygote cells, 4 hom/hom cells.
        let (mut ri, mut rj) = (Vec::new(), Vec::new());
        for gi in ALL {
            for gj in ALL {
                ri.push(gi);
                rj.push(gj);
            }
        }
        let g = planes(&[ri, rj]);
        let got = pair_counts(&g, 0, 1);
        assert_eq!(got, counts(9, 3, 3, 1, 2, 4));
        assert_eq!(got, naive(&g, 0, 1));
        assert_eq!(got.ibs1(), 3 + 3 - 2, "IBS1 = het_i + het_j - 2*het_het");
        assert_eq!(got.ibs2(), 9 - 2 - 4);
        assert_eq!(got.hom_hom + got.het_i + got.het_j, got.n_snp + got.het_het);
    }

    #[test]
    fn exhaustive_table_is_position_independent() {
        // Slide the same 16-variant table across word boundaries: the counts must not
        // depend on which 64-bit word the informative variants land in.
        for offset in [0usize, 1, 48, 62, 63, 64, 65, 100, 127, 128, 190] {
            let (mut ri, mut rj) = (vec![Missing; offset], vec![Missing; offset]);
            for gi in ALL {
                for gj in ALL {
                    ri.push(gi);
                    rj.push(gj);
                }
            }
            let g = planes(&[ri, rj]);
            assert_eq!(
                pair_counts(&g, 0, 1),
                counts(9, 3, 3, 1, 2, 4),
                "counts changed at offset {offset}"
            );
        }
    }

    #[test]
    fn het_counts_are_masked_by_the_other_sample() {
        // i is het everywhere; j is called at only the first two of four variants. If
        // het_i were taken over i's own non-missing set it would be 4, not 2. This is
        // the single most likely way to get kinship subtly wrong.
        let g = planes(&[vec![Het, Het, Het, Het], vec![Het, HomA1, Missing, Missing]]);
        let got = pair_counts(&g, 0, 1);
        assert_eq!(got, counts(2, 2, 1, 1, 0, 0));
        assert_eq!(got, naive(&g, 0, 1));
    }

    #[test]
    fn ibs0_requires_opposite_homozygotes() {
        // Same homozygote, opposite homozygotes, het against each homozygote, and a
        // homozygote against a missing call.
        let g = planes(&[
            vec![HomA1, HomA1, HomA2, Het, Het, HomA1, HomA2],
            vec![HomA1, HomA2, HomA1, HomA1, HomA2, Missing, Missing],
        ]);
        let got = pair_counts(&g, 0, 1);
        assert_eq!(
            got.ibs0, 2,
            "only the two opposite-homozygote sites are IBS0"
        );
        assert_eq!(got.hom_hom, 3);
        assert_eq!(got.n_snp, 5);
        assert_eq!(got, naive(&g, 0, 1));
    }

    #[test]
    fn self_pair_is_the_samples_own_profile() {
        let row = vec![HomA1, Het, HomA2, Missing, Het, HomA1, Missing, Het];
        let g = planes(&[row]);
        let got = pair_counts(&g, 0, 0);
        assert_eq!(got.n_snp, 6, "own non-missing count");
        assert_eq!(got.het_i, 3);
        assert_eq!(got.het_j, 3);
        assert_eq!(
            got.het_het, 3,
            "a sample is het with itself exactly where it is het"
        );
        assert_eq!(got.hom_hom, 3);
        assert_eq!(
            got.ibs0, 0,
            "a sample can never be an opposite homozygote to itself"
        );
        assert_eq!(got.ibs1(), 0);
        assert_eq!(got.ibs2(), 6);
    }

    #[test]
    fn swapping_the_pair_transposes_het_i_and_het_j() {
        let g = random_genotypes(&mut XorShift::new(0xA5A5_1234), 6, 300, 0.2);
        for (i, j) in [(0usize, 1usize), (2, 5), (3, 4)] {
            let a = pair_counts(&g, i, j);
            let b = pair_counts(&g, j, i);
            assert_eq!(a.het_i, b.het_j);
            assert_eq!(a.het_j, b.het_i);
            assert_eq!(a.n_snp, b.n_snp);
            assert_eq!(a.het_het, b.het_het);
            assert_eq!(a.hom_hom, b.hom_hom);
            assert_eq!(a.ibs0, b.ibs0);
            assert_ne!(
                a.het_i, a.het_j,
                "the fixture must not make the check vacuous"
            );
        }
    }

    #[test]
    fn empty_and_all_missing_inputs() {
        let empty = planes(&[vec![], vec![]]);
        assert_eq!(empty.words_per_sample(), 0);
        assert_eq!(pair_counts(&empty, 0, 1), PairCounts::default());

        let gone = planes(&[vec![Missing; 130], vec![HomA1; 130]]);
        assert_eq!(pair_counts(&gone, 0, 1), PairCounts::default());
        assert_eq!(pair_counts(&gone, 1, 0), PairCounts::default());
    }

    // ---- randomized differential test ---------------------------------------------

    /// Marsaglia xorshift64, so the fixtures are reproducible without a dependency.
    struct XorShift(u64);

    impl XorShift {
        fn new(seed: u64) -> Self {
            assert!(seed != 0, "xorshift cannot be seeded with zero");
            XorShift(seed)
        }

        fn next_u64(&mut self) -> u64 {
            let mut x = self.0;
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            self.0 = x;
            x
        }

        /// Uniform in `[0, 1)`.
        fn next_f64(&mut self) -> f64 {
            (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn random_genotypes(
        rng: &mut XorShift,
        n_samples: usize,
        n_variants: usize,
        missing: f64,
    ) -> Genotypes {
        let rows: Vec<Vec<Geno>> = (0..n_samples)
            .map(|_| {
                (0..n_variants)
                    .map(|_| {
                        if rng.next_f64() < missing {
                            Missing
                        } else {
                            match rng.next_u64() % 4 {
                                0 => HomA1,
                                1 | 2 => Het,
                                _ => HomA2,
                            }
                        }
                    })
                    .collect()
            })
            .collect();
        planes(&rows)
    }

    #[test]
    fn xorshift_is_deterministic_and_not_degenerate() {
        let a: Vec<u64> = (0..8).map(|_| XorShift::new(1).next_u64()).collect();
        assert!(
            a.windows(2).all(|w| w[0] == w[1]),
            "same seed must give the same value"
        );
        let mut rng = XorShift::new(1);
        let seq: Vec<u64> = (0..64).map(|_| rng.next_u64()).collect();
        assert_eq!(
            seq.len(),
            seq.iter().collect::<std::collections::HashSet<_>>().len()
        );
        let mut rng = XorShift::new(0xDEAD_BEEF);
        let mean: f64 = (0..4096).map(|_| rng.next_f64()).sum::<f64>() / 4096.0;
        assert!((mean - 0.5).abs() < 0.03, "next_f64 mean was {mean}");
    }

    #[test]
    fn kernel_matches_naive_on_random_matrices() {
        let mut rng = XorShift::new(0x0FED_C0DE_1234_5678);
        let mut checked = 0;
        // Variant counts deliberately straddle word boundaries; missing rates run from
        // none to almost everything.
        for &n_variants in &[1usize, 7, 63, 64, 65, 127, 128, 129, 200, 511, 512, 1000] {
            for &missing in &[0.0f64, 0.02, 0.25, 0.6, 0.95] {
                let g = random_genotypes(&mut rng, 5, n_variants, missing);
                for i in 0..g.n_samples {
                    for j in 0..g.n_samples {
                        let got = pair_counts(&g, i, j);
                        assert_eq!(
                            got,
                            naive(&g, i, j),
                            "mismatch at n_variants={n_variants} missing={missing} pair=({i},{j})"
                        );
                        // Structural invariants that must hold for any input.
                        assert!(got.n_snp <= n_variants as u32);
                        assert!(got.het_het <= got.het_i.min(got.het_j));
                        assert!(got.ibs0 <= got.hom_hom);
                        // hom_hom = n_snp - het_i - het_j + het_het, rearranged so the
                        // u32 arithmetic cannot underflow on an all-heterozygous pair.
                        assert_eq!(got.hom_hom + got.het_i + got.het_j, got.n_snp + got.het_het);
                        assert_eq!(got.n_snp, got.ibs0 + got.ibs1() + got.ibs2());
                        checked += 1;
                    }
                }
            }
        }
        assert_eq!(checked, 12 * 5 * 25);
    }

    #[test]
    fn kernel_matches_naive_with_extreme_per_sample_missingness() {
        // Per-sample missing rates that differ wildly, which is exactly the regime where
        // a non-pairwise heterozygote count would diverge.
        let mut rng = XorShift::new(0x5EED_0001);
        let rates = [0.0f64, 0.02, 0.3, 0.7, 0.98];
        let n_variants = 777;
        let rows: Vec<Vec<Geno>> = rates
            .iter()
            .map(|&r| {
                (0..n_variants)
                    .map(|_| {
                        if rng.next_f64() < r {
                            Missing
                        } else {
                            match rng.next_u64() % 3 {
                                0 => HomA1,
                                1 => Het,
                                _ => HomA2,
                            }
                        }
                    })
                    .collect()
            })
            .collect();
        let g = planes(&rows);
        let mut differed = 0;
        for i in 0..g.n_samples {
            for j in 0..g.n_samples {
                let got = pair_counts(&g, i, j);
                assert_eq!(got, naive(&g, i, j), "pair ({i},{j})");
                if i != j {
                    let own_het =
                        (0..n_variants).filter(|&m| decode(&g, i, m) == Het).count() as u32;
                    if got.het_i != own_het {
                        differed += 1;
                    }
                }
            }
        }
        assert!(
            differed >= 15,
            "fixture should exercise the pairwise mask, differed={differed}"
        );
    }

    // ---- parallel determinism -------------------------------------------------------

    fn upper_triangle(n: usize) -> Vec<(usize, usize)> {
        (0..n)
            .flat_map(|i| (i + 1..n).map(move |j| (i, j)))
            .collect()
    }

    #[test]
    fn all_pairs_is_independent_of_thread_count() {
        let g = random_genotypes(&mut XorShift::new(0x7E57_1010), 40, 1500, 0.15);
        let pairs = upper_triangle(g.n_samples);
        assert_eq!(pairs.len(), 780);

        let sequential: Vec<PairCounts> =
            pairs.iter().map(|&(i, j)| pair_counts(&g, i, j)).collect();

        for threads in [1usize, 2, 3, 8, 17] {
            let pool = rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("thread pool");
            let parallel = pool.install(|| all_pairs(&g, &pairs));
            assert_eq!(
                parallel, sequential,
                "all_pairs differed with {threads} threads"
            );
        }

        // Not vacuous: the fixture must actually contain varied counts.
        assert!(sequential.iter().any(|c| c.ibs0 > 0));
        assert!(sequential.windows(2).any(|w| w[0] != w[1]));
    }

    #[test]
    fn all_pairs_preserves_the_requested_order() {
        let g = random_genotypes(&mut XorShift::new(0x0DD_0DD), 8, 500, 0.1);
        // Out of order, repeated, reversed, and self-pairs.
        let pairs = [
            (5usize, 2usize),
            (0, 0),
            (2, 5),
            (7, 1),
            (5, 2),
            (3, 3),
            (1, 7),
        ];
        let got = all_pairs(&g, &pairs);
        assert_eq!(got.len(), pairs.len());
        for (k, &(i, j)) in pairs.iter().enumerate() {
            assert_eq!(got[k], pair_counts(&g, i, j), "slot {k}");
        }
        assert_eq!(got[0], got[4], "the same pair must give the same counts");
        assert_eq!(got[0].het_i, got[2].het_j, "(5,2) and (2,5) are transposes");
        assert!(all_pairs(&g, &[]).is_empty());
    }

    // ---- throughput -----------------------------------------------------------------

    /// Timing probe. Run with:
    /// `cargo test -p king-core --release -- --ignored --nocapture bench_throughput`
    #[test]
    #[ignore = "timing probe, not a correctness test"]
    fn bench_throughput() {
        use std::time::Instant;

        let n_samples = 200;
        let n_variants = 50_000;
        let g = random_genotypes(&mut XorShift::new(0xBEEF_0001), n_samples, n_variants, 0.05);
        let pairs = upper_triangle(n_samples);

        // Warm up the pool and the caches, and confirm the timed path is still correct.
        let warm = all_pairs(&g, &pairs[..64]);
        assert_eq!(warm[7], naive(&g, pairs[7].0, pairs[7].1));

        let t0 = Instant::now();
        let out = all_pairs(&g, &pairs);
        let elapsed = t0.elapsed();

        let cells = pairs.len() as f64 * n_variants as f64;
        let secs = elapsed.as_secs_f64();
        println!(
            "{} samples x {} variants = {} pairs in {:.3} ms => {:.2} G pair-genotypes/s ({:.1} M pairs/s)",
            n_samples,
            n_variants,
            pairs.len(),
            secs * 1e3,
            cells / secs / 1e9,
            pairs.len() as f64 / secs / 1e6,
        );

        assert_eq!(out.len(), pairs.len());
        assert!(out.iter().all(|c| c.n_snp > 0));
        // Generous enough to survive a loaded machine while still failing on an
        // accidental order-of-magnitude regression.
        assert!(secs < 5.0, "20k pairs x 50k SNPs took {secs:.3} s");
    }
}
