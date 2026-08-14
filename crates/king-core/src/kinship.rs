//! Pairwise relatedness estimators computed from [`PairCounts`].
//!
//! Every formula here is transcribed from `docs/VERIFIED_FORMULAS.md`, which states the
//! numeric check against the reference KING 2.3.2 binary. The counts arrive as exact
//! integers; each estimator converts them to `f64` once and performs a **single**
//! division, so the result is the correctly-rounded value of an exact rational. Counts
//! below `2^53` convert and combine losslessly, so there is no accumulation error to
//! worry about.
//!
//! # Degenerate denominators
//!
//! Real data reaches every one of these. Each case below was produced deliberately by
//! constructing a PLINK fileset and running the reference binary on it; the behaviour
//! recorded is what the reference actually printed, not what IEEE happens to give.
//!
//! | Situation | Reference output | Handling here |
//! | --- | --- | --- |
//! | `n_snp == 0` (pair with no overlapping calls) | `nan` for `IBS`, `Dist`, `HetConc`, `Het2\|1`, `Het1\|2`, `HomConc`, and for the `HetHet`/`IBS0` proportion columns | plain IEEE `0.0 / 0.0` → `NaN` |
//! | `hom_hom == 0` (both samples heterozygous everywhere) | `HomConc` = `nan` | plain IEEE → `NaN` |
//! | `het_i == 0` | `Het2\|1` = `nan` | plain IEEE → `NaN` |
//! | within-family, `het_i + het_j == 0`, `ibs0 == 0` | `Kinship` = `nan` | plain IEEE → `NaN` |
//! | within-family, `het_i + het_j == 0`, `ibs0 > 0` | `Kinship` = `-inf` | plain IEEE → `NEG_INFINITY` |
//! | **between-family, `min(het_i, het_j) == 0`** | `Kinship` = **`0.0000`** | **explicit guard returning `0.0`** |
//!
//! The between-family guard is the one place where the reference is *not* IEEE: it
//! returns exactly zero rather than `nan`/`-inf`, for every combination of `ibs0` and of
//! which sample has no heterozygotes, and regardless of `n_snp`. Third-party
//! implementations disagree here (plink2 yields `-inf`, SNPRelate `NaN`); the reference
//! is the authority and it prints `0.0000`.
//!
//! Non-finite values are *not* clamped and are *not* filtered — a within-family pair with
//! one heterozygote each and a thousand opposite-homozygote sites prints `-1000.0000`
//! verbatim. Row filtering is the writer's job, and the two writers differ: `.kin` omits
//! a pair whose kinship is not finite, while `.ibs` omits only `nan` and prints `-inf`.

use crate::{PairCounts, Scope};

/// Lossless widening of a count to `f64`.
#[inline]
fn c(x: u32) -> f64 {
    f64::from(x)
}

/// Kinship coefficient for one pair, using the estimator selected by `scope`.
///
/// * [`Scope::WithinFamily`] — `(HetHet - 2*IBS0) / (Het_i + Het_j)`
/// * [`Scope::BetweenFamily`] — `0.5 + (2*HetHet - 4*IBS0 - Het_i - Het_j) / (4*min(Het_i, Het_j))`
///
/// The two coincide exactly when `het_i == het_j` and diverge otherwise; which one applies
/// is decided solely by whether the pair shares an `FID`, never guessed.
///
/// See the module docs for the degenerate-denominator behaviour, which is matched to the
/// reference binary rather than left to IEEE in the between-family case.
pub fn kinship(counts: &PairCounts, scope: Scope) -> f64 {
    let het_i = c(counts.het_i);
    let het_j = c(counts.het_j);
    let het_het = c(counts.het_het);
    let ibs0 = c(counts.ibs0);

    match scope {
        Scope::WithinFamily => (het_het - 2.0 * ibs0) / (het_i + het_j),
        Scope::BetweenFamily => {
            // Observed: the reference prints 0.0000 whenever the smaller heterozygote
            // count is zero, instead of the nan / -inf that the raw expression yields.
            if counts.het_i == 0 || counts.het_j == 0 {
                return 0.0;
            }
            let smaller = if het_i < het_j { het_i } else { het_j };
            0.5 + (2.0 * het_het - 4.0 * ibs0 - het_i - het_j) / (4.0 * smaller)
        }
    }
}

/// `N_IBS1 = Het_i + Het_j - 2*HetHet`, as `f64`.
///
/// Computed in floating point rather than through [`PairCounts::ibs1`] so that a
/// malformed count set cannot panic on unsigned underflow; for valid counts the two are
/// identical and both exact.
#[inline]
fn ibs1(counts: &PairCounts) -> f64 {
    c(counts.het_i) + c(counts.het_j) - 2.0 * c(counts.het_het)
}

/// `N_IBS2 = N_SNP - IBS0 - N_IBS1`, as `f64`.
#[inline]
fn ibs2(counts: &PairCounts) -> f64 {
    c(counts.n_snp) - c(counts.ibs0) - ibs1(counts)
}

/// Mean IBS allele sharing: `(N_IBS1 + 2*N_IBS2) / N_SNP`. Ranges 0..=2.
pub fn ibs_mean(counts: &PairCounts) -> f64 {
    (ibs1(counts) + 2.0 * ibs2(counts)) / c(counts.n_snp)
}

/// Mean squared genotype distance: `(N_IBS1 + 4*IBS0) / N_SNP`.
///
/// This is **not** `2 - ibs_mean`; that identity holds only when `ibs0 == 0`.
pub fn dist(counts: &PairCounts) -> f64 {
    (ibs1(counts) + 4.0 * c(counts.ibs0)) / c(counts.n_snp)
}

/// Heterozygote Jaccard concordance: `HetHet / (Het_i + Het_j - HetHet)`.
pub fn het_concordance(counts: &PairCounts) -> f64 {
    c(counts.het_het) / (c(counts.het_i) + c(counts.het_j) - c(counts.het_het))
}

/// `Het2|1` — `HetHet / Het_i`: fraction of sample `i`'s heterozygotes that `j` shares.
pub fn het_2given1(counts: &PairCounts) -> f64 {
    c(counts.het_het) / c(counts.het_i)
}

/// `Het1|2` — `HetHet / Het_j`.
pub fn het_1given2(counts: &PairCounts) -> f64 {
    c(counts.het_het) / c(counts.het_j)
}

/// Concordance among sites where both samples are homozygous:
/// `(HomHom - IBS0) / HomHom`.
pub fn hom_concordance(counts: &PairCounts) -> f64 {
    (c(counts.hom_hom) - c(counts.ibs0)) / c(counts.hom_hom)
}

/// The `HetHet` column of `.kin` / `.kin0`: `HetHet / N_SNP`, a proportion, not a count.
pub fn het_het_prop(counts: &PairCounts) -> f64 {
    c(counts.het_het) / c(counts.n_snp)
}

/// The `IBS0` column of `.kin` / `.kin0`: `IBS0 / N_SNP`, a proportion, not a count.
pub fn ibs0_prop(counts: &PairCounts) -> f64 {
    c(counts.ibs0) / c(counts.n_snp)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Round to 4 decimals the way `%.4f` does, so tests compare what would be printed.
    fn r4(x: f64) -> f64 {
        (x * 10_000.0).round() / 10_000.0
    }

    /// The four pairs tabulated in `docs/VERIFIED_FORMULAS.md`, with the raw counts read
    /// out of `docs/reference-captures/tiny/out_ibs/king.ibs{,0}` and the printed values
    /// read out of the same capture (Kinship cross-checked against `out_kinship`).
    struct RefPair {
        name: &'static str,
        scope: Scope,
        counts: PairCounts,
        ibs: f64,
        dist: f64,
        het_conc: f64,
        het_2given1: f64,
        het_1given2: f64,
        hom_conc: f64,
        kinship: f64,
        n_ibs1: u32,
        n_ibs2: u32,
    }

    const REFERENCE_PAIRS: [RefPair; 4] = [
        RefPair {
            name: "f1dad/f1kid1",
            scope: Scope::WithinFamily,
            counts: PairCounts {
                n_snp: 1970,
                het_i: 602,
                het_j: 595,
                het_het: 302,
                ibs0: 0,
                hom_hom: 1075,
            },
            ibs: 1.6990,
            dist: 0.3010,
            het_conc: 0.3374,
            het_2given1: 0.5017,
            het_1given2: 0.5076,
            hom_conc: 1.0000,
            kinship: 0.2523,
            n_ibs1: 593,
            n_ibs2: 1377,
        },
        RefPair {
            name: "f1dad/f1kid2",
            scope: Scope::WithinFamily,
            counts: PairCounts {
                n_snp: 1942,
                het_i: 595,
                het_j: 537,
                het_het: 269,
                ibs0: 0,
                hom_hom: 1079,
            },
            ibs: 1.6941,
            dist: 0.3059,
            het_conc: 0.3117,
            het_2given1: 0.4521,
            het_1given2: 0.5009,
            hom_conc: 1.0000,
            kinship: 0.2376,
            n_ibs1: 594,
            n_ibs2: 1348,
        },
        RefPair {
            name: "f1dad/f2dad",
            scope: Scope::BetweenFamily,
            counts: PairCounts {
                n_snp: 1960,
                het_i: 595,
                het_j: 558,
                het_het: 337,
                ibs0: 28,
                hom_hom: 1144,
            },
            ibs: 1.7270,
            dist: 0.3015,
            het_conc: 0.4130,
            het_2given1: 0.5664,
            het_1given2: 0.6039,
            hom_conc: 0.9755,
            kinship: 0.2352,
            n_ibs1: 479,
            n_ibs2: 1453,
        },
        RefPair {
            name: "f1dad/f2mom",
            scope: Scope::BetweenFamily,
            counts: PairCounts {
                n_snp: 1803,
                het_i: 549,
                het_j: 566,
                het_het: 221,
                ibs0: 110,
                hom_hom: 909,
            },
            ibs: 1.5047,
            dist: 0.6173,
            het_conc: 0.2472,
            het_2given1: 0.4026,
            het_1given2: 0.3905,
            hom_conc: 0.8790,
            kinship: -0.0068,
            n_ibs1: 673,
            n_ibs2: 1020,
        },
    ];

    #[test]
    fn reproduces_the_verified_reference_pairs_to_four_decimals() {
        for p in &REFERENCE_PAIRS {
            let k = &p.counts;
            assert_eq!(k.ibs1(), p.n_ibs1, "{}: N_IBS1", p.name);
            assert_eq!(k.ibs2(), p.n_ibs2, "{}: N_IBS2", p.name);
            assert_eq!(r4(ibs_mean(k)), p.ibs, "{}: IBS", p.name);
            assert_eq!(r4(dist(k)), p.dist, "{}: Dist", p.name);
            assert_eq!(r4(het_concordance(k)), p.het_conc, "{}: HetConc", p.name);
            assert_eq!(r4(het_2given1(k)), p.het_2given1, "{}: Het2|1", p.name);
            assert_eq!(r4(het_1given2(k)), p.het_1given2, "{}: Het1|2", p.name);
            assert_eq!(r4(hom_concordance(k)), p.hom_conc, "{}: HomConc", p.name);
            assert_eq!(r4(kinship(k, p.scope)), p.kinship, "{}: Kinship", p.name);
        }
    }

    #[test]
    fn proportion_columns_match_the_capture() {
        // docs/VERIFIED_FORMULAS.md: HetHet 302 / N_SNP 1970 -> 0.1533,
        // IBS0 28 / N_SNP 1960 -> 0.0143.
        assert_eq!(r4(het_het_prop(&REFERENCE_PAIRS[0].counts)), 0.1533);
        assert_eq!(r4(ibs0_prop(&REFERENCE_PAIRS[2].counts)), 0.0143);
        assert_eq!(r4(het_het_prop(&REFERENCE_PAIRS[2].counts)), 0.1719);
        assert_eq!(r4(ibs0_prop(&REFERENCE_PAIRS[0].counts)), 0.0000);
    }

    #[test]
    fn choosing_the_wrong_scope_changes_the_answer() {
        // The whole point of the two estimators: they must not be interchangeable when
        // the heterozygote counts differ. f1dad/f2mom is 549 vs 566.
        let k = &REFERENCE_PAIRS[3].counts;
        assert_ne!(
            r4(kinship(k, Scope::WithinFamily)),
            r4(kinship(k, Scope::BetweenFamily))
        );
    }

    #[test]
    fn the_two_estimators_coincide_when_het_counts_are_equal() {
        let k = PairCounts {
            n_snp: 5000,
            het_i: 1200,
            het_j: 1200,
            het_het: 500,
            ibs0: 37,
            hom_hom: 2600,
        };
        let within = kinship(&k, Scope::WithinFamily);
        let between = kinship(&k, Scope::BetweenFamily);
        assert!(
            (within - between).abs() < 1e-15,
            "within {within} vs between {between}"
        );
    }

    #[test]
    fn duplicate_pair_gives_one_half_under_both_estimators() {
        let k = PairCounts {
            n_snp: 2000,
            het_i: 640,
            het_j: 640,
            het_het: 640,
            ibs0: 0,
            hom_hom: 1360,
        };
        assert_eq!(kinship(&k, Scope::WithinFamily), 0.5);
        assert_eq!(kinship(&k, Scope::BetweenFamily), 0.5);
    }

    // ---------------------------------------------------------------------------------
    // Degenerate denominators. Every expectation below was read off the reference binary
    // run on a purpose-built fileset; see the module docs.
    // ---------------------------------------------------------------------------------

    /// Reference `.ibs` row `famD D1 D2`: two samples heterozygous at all 2000 SNPs.
    /// Printed `... NHomHom 0 ... HomConc nan  Kinship 0.5000`.
    #[test]
    fn zero_homozygous_sites_gives_nan_hom_concordance() {
        let k = PairCounts {
            n_snp: 2000,
            het_i: 2000,
            het_j: 2000,
            het_het: 2000,
            ibs0: 0,
            hom_hom: 0,
        };
        assert!(hom_concordance(&k).is_nan());
        assert_eq!(kinship(&k, Scope::WithinFamily), 0.5);
        assert_eq!(het_concordance(&k), 1.0);
        assert_eq!(ibs_mean(&k), 2.0);
        assert_eq!(dist(&k), 0.0);
    }

    /// Reference `.ibs` row `famB B1 B2`: both samples homozygous everywhere, 1000 of the
    /// 2000 sites opposite homozygotes. Printed
    /// `HetConc nan  Het2|1 nan  Het1|2 nan  HomConc 0.5000  Kinship -inf`,
    /// and the pair is absent from `.kin` because `.kin` drops non-finite kinship.
    #[test]
    fn within_family_zero_heterozygotes_with_ibs0_gives_negative_infinity() {
        let k = PairCounts {
            n_snp: 2000,
            het_i: 0,
            het_j: 0,
            het_het: 0,
            ibs0: 1000,
            hom_hom: 2000,
        };
        assert_eq!(kinship(&k, Scope::WithinFamily), f64::NEG_INFINITY);
        assert!(het_concordance(&k).is_nan());
        assert!(het_2given1(&k).is_nan());
        assert!(het_1given2(&k).is_nan());
        assert_eq!(hom_concordance(&k), 0.5);
        assert_eq!(ibs_mean(&k), 1.0);
        assert_eq!(dist(&k), 2.0);
    }

    /// Two identical all-homozygous samples in one family: numerator and denominator are
    /// both zero. The reference emitted no `.kin` and no `.ibs` row for such a pair, which
    /// is consistent with a `nan` kinship being filtered out by both writers.
    #[test]
    fn within_family_zero_heterozygotes_without_ibs0_gives_nan() {
        let k = PairCounts {
            n_snp: 2000,
            het_i: 0,
            het_j: 0,
            het_het: 0,
            ibs0: 0,
            hom_hom: 2000,
        };
        assert!(kinship(&k, Scope::WithinFamily).is_nan());
    }

    /// The between-family guard. All four of these were printed as `0.0000` by the
    /// reference in `.kin0` and `.ibs0` — never `nan`, never `-inf`.
    #[test]
    fn between_family_zero_min_heterozygotes_gives_exactly_zero() {
        // famA A1 / famC C1: identical all-homozygous samples, no IBS0.
        let identical_homs = PairCounts {
            n_snp: 2000,
            het_i: 0,
            het_j: 0,
            het_het: 0,
            ibs0: 0,
            hom_hom: 2000,
        };
        // famA A1 / famB B1: opposite homozygotes at every site.
        let opposite_homs = PairCounts {
            n_snp: 2000,
            het_i: 0,
            het_j: 0,
            het_het: 0,
            ibs0: 2000,
            hom_hom: 2000,
        };
        // famA A1 / famD D1: only one of the two has heterozygotes.
        let one_sided = PairCounts {
            n_snp: 2000,
            het_i: 0,
            het_j: 2000,
            het_het: 0,
            ibs0: 0,
            hom_hom: 0,
        };
        // famX X1 / famY Y2: no overlapping calls at all.
        let no_overlap = PairCounts::default();

        for k in [identical_homs, opposite_homs, one_sided, no_overlap] {
            let v = kinship(&k, Scope::BetweenFamily);
            assert_eq!(v, 0.0, "{k:?}");
            assert!(v.is_sign_positive(), "must be +0.0, not -0.0: {k:?}");
        }
        // The guard is exactly `min == 0`: with a single heterozygote it must vanish and
        // the raw expression must come through, however extreme.
        let one_het = PairCounts {
            n_snp: 2000,
            het_i: 1,
            het_j: 1,
            het_het: 0,
            ibs0: 1000,
            hom_hom: 1998,
        };
        assert_eq!(kinship(&one_het, Scope::BetweenFamily), -1000.0);
    }

    /// Reference row `famX X1 famY Y2` in `.kin0`/`.ibs0`: a pair with no site called in
    /// both. `N_SNP 0`, every derived column `nan`, `Kinship 0.0000`.
    #[test]
    fn zero_overlapping_snps_gives_nan_everywhere_but_kinship() {
        let k = PairCounts::default();
        assert!(ibs_mean(&k).is_nan());
        assert!(dist(&k).is_nan());
        assert!(het_concordance(&k).is_nan());
        assert!(het_2given1(&k).is_nan());
        assert!(het_1given2(&k).is_nan());
        assert!(hom_concordance(&k).is_nan());
        assert!(het_het_prop(&k).is_nan());
        assert!(ibs0_prop(&k).is_nan());
        assert!(kinship(&k, Scope::WithinFamily).is_nan());
        assert_eq!(kinship(&k, Scope::BetweenFamily), 0.0);
    }

    /// One sample homozygous everywhere, the other heterozygous everywhere: `Het2|1` is
    /// `0/0` and prints `nan`, while `Het1|2` is `0/2000` and prints `0.0000`. The
    /// reference row `famX X1 famY Y1` shows exactly this asymmetry, which is the check
    /// that the two are not accidentally the same expression.
    #[test]
    fn conditional_heterozygosity_is_asymmetric() {
        let k = PairCounts {
            n_snp: 1000,
            het_i: 0,
            het_j: 1000,
            het_het: 0,
            ibs0: 0,
            hom_hom: 0,
        };
        assert!(het_2given1(&k).is_nan());
        assert_eq!(het_1given2(&k), 0.0);
        assert_eq!(het_concordance(&k), 0.0);
    }

    /// Reference rows `famG G1 G2` (`.ibs`, within) and `famG G1 famH H1` (`.ibs0`,
    /// between) both printed `-1000.0000`: extreme negative kinship passes through
    /// unmodified. Nothing is clamped to 0 or to -1.
    #[test]
    fn extreme_negative_kinship_is_not_clamped() {
        let k = PairCounts {
            n_snp: 2000,
            het_i: 1,
            het_j: 1,
            het_het: 0,
            ibs0: 1000,
            hom_hom: 1998,
        };
        assert_eq!(kinship(&k, Scope::WithinFamily), -1000.0);
        assert_eq!(kinship(&k, Scope::BetweenFamily), -1000.0);
        // and the accompanying columns from the same captured row
        assert_eq!(r4(ibs_mean(&k)), 0.9990);
        assert_eq!(r4(dist(&k)), 2.0010);
        assert_eq!(r4(hom_concordance(&k)), 0.4995);
    }

    /// Counts near the top of `u32` must still convert and combine losslessly, so the
    /// division is of exact integers.
    #[test]
    fn large_counts_stay_exact() {
        let k = PairCounts {
            n_snp: 4_000_000_000,
            het_i: 2_000_000_000,
            het_j: 2_000_000_000,
            het_het: 1_000_000_000,
            ibs0: 0,
            hom_hom: 1_000_000_000,
        };
        assert_eq!(kinship(&k, Scope::WithinFamily), 0.25);
        assert_eq!(kinship(&k, Scope::BetweenFamily), 0.25);
        assert_eq!(het_het_prop(&k), 0.25);
    }
}
