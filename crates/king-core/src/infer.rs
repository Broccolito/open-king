//! Relationship inference: turning estimates into classes, and reading the expected
//! values back out of a declared pedigree.
//!
//! Two independent jobs live here.
//!
//! * [`classify`] maps an estimated kinship (plus the IBS0 proportion, which is the only
//!   thing that separates parent–offspring from full siblings) onto a [`Relationship`].
//! * [`Pedigree`], [`pedigree_kinship`] and [`pedigree_z0`] compute the *expected* values
//!   from the declared `PAT`/`MAT` graph. These populate the `Phi` and `Z0` columns of
//!   `.kin`, which are pedigree-derived and owe nothing to the genotypes.
//!
//! # Boundary semantics
//!
//! Cutoffs are **inclusive at the lower edge of the higher class**: a kinship of exactly
//! [`cutoff::FIRST`] is 1st degree, not 2nd. This was chosen by construction, not by
//! taste — filesets were built whose within-family kinship lands exactly on each printed
//! constant, and the reference put every one of them in the higher class:
//!
//! | Constructed kinship | Reference class |
//! | --- | --- |
//! | `0.3540` | Dup/MZ |
//! | `0.1770` | 1st degree |
//! | `0.0884` | 2nd degree |
//!
//! Non-finite input falls through every comparison and lands on [`Relationship::Unrelated`],
//! which keeps the function total; a `NaN` kinship means a degenerate pair whose row the
//! reference omits from `.kin` entirely, so the caller should be filtering it out anyway.
//!
//! # The reference's real cutoffs are the exact powers, not the printed decimals
//!
//! [`cutoff`] holds the rounded values published in the paper and the manual. Bisecting
//! the reference binary with constructed kinships shows its actual class boundaries are
//! the *unrounded* `2^(-k-3/2)` grid:
//!
//! | Boundary | Constructed kinship → class | Bracket | Exact power |
//! | --- | --- | --- | --- |
//! | Dup/MZ vs 1st | `0.353553` → 1st, `0.353554` → Dup/MZ | `(0.353553, 0.353554]` | `2^(-3/2) = 0.35355339…` |
//! | 1st vs 2nd | `0.176776` → 2nd, `0.176777` → 1st | `(0.176776, 0.176777]` | `2^(-5/2) = 0.17677670…` |
//! | 2nd vs 3rd | `0.088388` → 3rd, `0.088389` → 2nd | `(0.088388, 0.088389]` | `2^(-7/2) = 0.08838835…` |
//!
//! So [`Cutoffs::PUBLISHED`] and [`Cutoffs::REFERENCE`] disagree on kinships inside three
//! narrow bands — `[0.3535534, 0.354)`, `[0.1767767, 0.177)`, `[0.0883883, 0.0884)`.
//! [`classify`] uses `PUBLISHED` because those are the constants the crate publishes in
//! [`cutoff`]; [`classify_with`] lets a caller chasing byte parity select `REFERENCE`.

use std::collections::HashMap;

use king_io::Sample;

use crate::{cutoff, Relationship};

/// The kinship boundaries between relationship classes.
///
/// Compared with `>=`, so each field is the smallest kinship belonging to the class it
/// names.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Cutoffs {
    /// At or above this: duplicate / MZ twin.
    pub dup_mz: f64,
    /// At or above this: 1st degree.
    pub first: f64,
    /// At or above this: 2nd degree.
    pub second: f64,
    /// At or above this: 3rd degree.
    pub third: f64,
    /// At or above this: 4th degree. Below it: unrelated.
    pub fourth: f64,
}

impl Cutoffs {
    /// The rounded constants published in the paper and the KING manual, as re-exported
    /// by [`cutoff`]. This is what [`classify`] uses.
    pub const PUBLISHED: Self = Self {
        dup_mz: cutoff::DUP_MZ,
        first: cutoff::FIRST,
        second: cutoff::SECOND,
        third: cutoff::THIRD,
        fourth: cutoff::FOURTH,
    };

    /// The boundaries the reference binary actually uses: the unrounded `2^(-k-3/2)`
    /// grid, bracketed to six decimal places by construction (see the module docs). The
    /// top three were measured directly; `third` and `fourth` continue the grid and are
    /// **not** separately confirmed, because the reference's within-family summary lumps
    /// 4th-degree and unrelated pairs together.
    pub const REFERENCE: Self = Self {
        dup_mz: 0.353_553_390_593_273_8,  // 2^(-3/2)
        first: 0.176_776_695_296_636_9,   // 2^(-5/2)
        second: 0.088_388_347_648_318_45, // 2^(-7/2)
        third: 0.044_194_173_824_159_22,  // 2^(-9/2)
        fourth: 0.022_097_086_912_079_61, // 2^(-11/2)
    };
}

/// Classify a pair from its estimated kinship and IBS0 proportion.
///
/// `po_threshold` is the IBS0-proportion cutoff separating parent–offspring from full
/// siblings among 1st-degree pairs; the reference derives it from the data and announces
/// it as `1st-degree relatives are treated as parent-offspring if IBS0 < <x>`, so the
/// comparison here is **strictly less than** — a pair sitting exactly on the threshold is
/// a full sibling pair.
///
/// Uses [`Cutoffs::PUBLISHED`]; see [`classify_with`] and the module docs.
pub fn classify(kinship: f64, ibs0_prop: f64, po_threshold: f64) -> Relationship {
    classify_with(kinship, ibs0_prop, po_threshold, &Cutoffs::PUBLISHED)
}

/// [`classify`] against a caller-chosen set of boundaries.
pub fn classify_with(
    kinship: f64,
    ibs0_prop: f64,
    po_threshold: f64,
    cutoffs: &Cutoffs,
) -> Relationship {
    if kinship >= cutoffs.dup_mz {
        Relationship::DupMz
    } else if kinship >= cutoffs.first {
        // phi cannot separate PO from FS — both sit at 0.25. IBS0 can: a true PO pair
        // shares an allele IBD at every autosomal site, so opposite homozygotes only
        // arise from genotyping error.
        if ibs0_prop < po_threshold {
            Relationship::ParentOffspring
        } else {
            Relationship::FullSib
        }
    } else if kinship >= cutoffs.second {
        Relationship::Second
    } else if kinship >= cutoffs.third {
        Relationship::Third
    } else if kinship >= cutoffs.fourth {
        Relationship::Fourth
    } else {
        Relationship::Unrelated
    }
}

// -------------------------------------------------------------------------------------
// Pedigree side
// -------------------------------------------------------------------------------------

/// The declared parent graph of a sample set, indexed exactly like the `.fam` rows it was
/// built from.
///
/// Parent references are resolved by IID **within the sample's own FID**, which is what
/// PLINK's format means. Any reference that does not resolve — `"0"`, an unknown IID, an
/// IID in another family — is treated as an unknown parent.
///
/// Construction breaks pedigree cycles: if a `PAT`/`MAT` edge would make an individual its
/// own ancestor, that single edge is dropped, so the graph handed to the recurrence is
/// always acyclic and the recursion always terminates. Real `.fam` files should never need
/// this, but a corrupt one must not hang or overflow the stack.
#[derive(Debug, Clone, Default)]
pub struct Pedigree {
    pat: Vec<Option<usize>>,
    mat: Vec<Option<usize>>,
    /// Longest path to a founder. An ancestor always has strictly smaller depth, which is
    /// what lets the recurrence pick the individual it is legal to expand.
    depth: Vec<u32>,
}

impl Pedigree {
    /// Build from `.fam` rows, in file order.
    pub fn from_samples(samples: &[Sample]) -> Self {
        let mut index: HashMap<(&str, &str), usize> = HashMap::with_capacity(samples.len());
        for (i, s) in samples.iter().enumerate() {
            index.entry((s.fid.as_str(), s.iid.as_str())).or_insert(i);
        }
        let resolve = |fid: &str, id: &str, self_idx: usize| -> Option<usize> {
            if id == "0" {
                return None;
            }
            index.get(&(fid, id)).copied().filter(|&p| p != self_idx)
        };

        let mut pat = Vec::with_capacity(samples.len());
        let mut mat = Vec::with_capacity(samples.len());
        for (i, s) in samples.iter().enumerate() {
            pat.push(resolve(&s.fid, &s.pat, i));
            mat.push(resolve(&s.fid, &s.mat, i));
        }

        let depth = break_cycles_and_measure_depth(&mut pat, &mut mat);
        Self { pat, mat, depth }
    }

    /// Number of individuals.
    pub fn len(&self) -> usize {
        self.pat.len()
    }

    /// Whether the pedigree holds no individuals.
    pub fn is_empty(&self) -> bool {
        self.pat.is_empty()
    }

    /// Resolved `(PAT, MAT)` indices for an individual.
    pub fn parents(&self, i: usize) -> (Option<usize>, Option<usize>) {
        (self.pat[i], self.mat[i])
    }

    /// Whether two distinct individuals declare the *same two known* parents.
    ///
    /// Both parents must be known: the reference treats a pair that shares only a father,
    /// with both mothers unknown, as half siblings (`Z0` `0.500`), not as full siblings.
    pub fn is_full_sib_pair(&self, a: usize, b: usize) -> bool {
        a != b
            && self.pat[a].is_some()
            && self.mat[a].is_some()
            && self.pat[a] == self.pat[b]
            && self.mat[a] == self.mat[b]
    }

    /// [`pedigree_kinship`] with a throwaway cache. Convenient for one-off queries; use
    /// the free function with a shared [`KinshipCache`] when computing many pairs.
    pub fn kinship(&self, a: usize, b: usize) -> f64 {
        pedigree_kinship(self, &mut KinshipCache::default(), a, b)
    }

    /// [`pedigree_z0`] with a throwaway cache.
    pub fn z0(&self, a: usize, b: usize) -> f64 {
        pedigree_z0(self, &mut KinshipCache::default(), a, b)
    }
}

/// Depth-first pass that drops back edges and records each individual's generation depth.
fn break_cycles_and_measure_depth(
    pat: &mut [Option<usize>],
    mat: &mut [Option<usize>],
) -> Vec<u32> {
    const WHITE: u8 = 0;
    const GRAY: u8 = 1;
    const BLACK: u8 = 2;

    let n = pat.len();
    let mut mark = vec![WHITE; n];
    let mut depth = vec![0u32; n];
    let mut stack: Vec<(usize, u8)> = Vec::new();

    for start in 0..n {
        if mark[start] != WHITE {
            continue;
        }
        mark[start] = GRAY;
        stack.push((start, 0));
        while let Some(&(u, step)) = stack.last() {
            if step == 2 {
                let from_pat = pat[u].map_or(0, |p| depth[p] + 1);
                let from_mat = mat[u].map_or(0, |m| depth[m] + 1);
                depth[u] = from_pat.max(from_mat);
                mark[u] = BLACK;
                stack.pop();
                continue;
            }
            if let Some(top) = stack.last_mut() {
                top.1 = step + 1;
            }
            let parent = if step == 0 { pat[u] } else { mat[u] };
            let Some(v) = parent else { continue };
            match mark[v] {
                WHITE => {
                    mark[v] = GRAY;
                    stack.push((v, 0));
                }
                // v is still on the stack, so u is one of v's descendants: this edge
                // closes a cycle. Drop it rather than looping forever.
                GRAY => {
                    if step == 0 {
                        pat[u] = None;
                    } else {
                        mat[u] = None;
                    }
                }
                _ => {}
            }
        }
    }
    depth
}

/// Memo table for the kinship recurrence. Reuse one across all pairs of a family.
#[derive(Debug, Clone, Default)]
pub struct KinshipCache {
    memo: HashMap<(usize, usize), f64>,
}

impl KinshipCache {
    /// Discard everything memoised so far.
    pub fn clear(&mut self) {
        self.memo.clear();
    }
}

/// Expected kinship coefficient `phi` from the declared pedigree — the `Phi` column of
/// `.kin`.
///
/// The standard recurrence, memoised:
///
/// ```text
///     phi(i, i) = (1 + F_i) / 2,   F_i = phi(pat_i, mat_i)
///     phi(i, j) = [ phi(pat_i, j) + phi(mat_i, j) ] / 2
/// ```
///
/// where the individual expanded in the second line is always the one that cannot be an
/// ancestor of the other, and an unknown parent contributes zero. `F_i` is zero unless
/// *both* parents are known.
///
/// This is a real traversal, so it gives the right answer on inbred pedigrees where no
/// lookup table would: the child of first cousins has `phi = 0.28125` with each parent
/// and `F = 0.0625`, both of which the reference prints.
pub fn pedigree_kinship(ped: &Pedigree, cache: &mut KinshipCache, a: usize, b: usize) -> f64 {
    // Expand the deeper individual; ties broken by index so the memo key is canonical.
    let (x, y) = if (ped.depth[a], a) >= (ped.depth[b], b) {
        (a, b)
    } else {
        (b, a)
    };
    if let Some(&v) = cache.memo.get(&(x, y)) {
        return v;
    }
    let v = if x == y {
        let inbreeding = match (ped.pat[x], ped.mat[x]) {
            (Some(p), Some(m)) => pedigree_kinship(ped, cache, p, m),
            _ => 0.0,
        };
        0.5 * (1.0 + inbreeding)
    } else {
        // x is at least as deep as y and differs from it, so x is not an ancestor of y
        // and expanding x is legal. A founder x therefore contributes nothing.
        let via_pat = match ped.pat[x] {
            Some(p) => pedigree_kinship(ped, cache, p, y),
            None => 0.0,
        };
        let via_mat = match ped.mat[x] {
            Some(m) => pedigree_kinship(ped, cache, m, y),
            None => 0.0,
        };
        0.5 * (via_pat + via_mat)
    };
    cache.memo.insert((x, y), v);
    v
}

/// The largest pedigree kinship at which the reference still reports the full-sibling
/// `Z0` of `0.25`.
///
/// Bracketed by construction to `(0.296875, 0.3046875]`: sibling pairs built with
/// `phi` = 0.25, 0.28125, 0.2890625 and 0.296875 all printed `0.250`, while 0.3046875,
/// 0.3125, 0.3515625 and 0.375 all printed `0.000`. The conservative end of the bracket is
/// used here, so every measured case is reproduced and only the unobserved gap is at risk.
/// Reaching that gap needs siblings whose parents are themselves closer than first
/// cousins, so no ordinary consanguineous pedigree goes near it.
pub const FS_Z0_KINSHIP_LIMIT: f64 = 0.296875;

/// Expected `Pr[IBD = 0]` from the declared pedigree — the `Z0` column of `.kin`.
///
/// The reference does **not** compute the true identity coefficients; gene-dropping
/// simulation disagrees with it by 0.1 or more on inbred pedigrees. What it computes,
/// reproduced here, is:
///
/// ```text
///     full siblings, phi <= FS_Z0_KINSHIP_LIMIT   ->  0.25
///     otherwise, with z = 1 - 4*phi:
///         z >= 0.25                               ->  z
///         z <  0.25                               ->  0
/// ```
///
/// `1 - 4*phi` is exact for the canonical relationships — 0 for parent–offspring, 0.5 for
/// 2nd degree, 0.75 for 3rd, 1 for unrelated — and the reference carries it through to
/// non-canonical inbred values too (`phi` 0.09375 prints `0.625`). The floor is the part
/// that has to be observed: anything that would land strictly between 0 and the
/// full-sibling value of 0.25 is reported as `0.000` instead, verified at `phi` 0.203125
/// and 0.21875 (both printed `0.000`, not `0.187`/`0.125`) against `phi` 0.1875, which
/// printed `0.250`.
///
/// Reproduces all 437 `Z0`/`Phi` rows the reference emitted across six purpose-built
/// pedigree filesets, including three levels of consanguinity.
pub fn pedigree_z0(ped: &Pedigree, cache: &mut KinshipCache, a: usize, b: usize) -> f64 {
    let phi = pedigree_kinship(ped, cache, a, b);
    if ped.is_full_sib_pair(a, b) && phi <= FS_Z0_KINSHIP_LIMIT {
        return 0.25;
    }
    let z = 1.0 - 4.0 * phi;
    if z >= 0.25 {
        z
    } else {
        0.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The next representable `f64` below a positive, finite `x`.
    fn just_below(x: f64) -> f64 {
        assert!(x > 0.0 && x.is_finite());
        f64::from_bits(x.to_bits() - 1)
    }

    const NO_IBS0: f64 = 0.0;
    const PO_THR: f64 = 0.1;

    // ---------------------------------------------------------------------------------
    // classify
    // ---------------------------------------------------------------------------------

    #[test]
    fn every_published_cutoff_is_inclusive_at_its_lower_edge() {
        let cases = [
            (cutoff::DUP_MZ, Relationship::DupMz),
            (cutoff::FIRST, Relationship::ParentOffspring),
            (cutoff::SECOND, Relationship::Second),
            (cutoff::THIRD, Relationship::Third),
            (cutoff::FOURTH, Relationship::Fourth),
        ];
        for (k, expected) in cases {
            assert_eq!(classify(k, NO_IBS0, PO_THR), expected, "exactly on {k}");
        }
    }

    #[test]
    fn one_ulp_below_each_cutoff_drops_to_the_next_class() {
        let cases = [
            (cutoff::DUP_MZ, Relationship::ParentOffspring),
            (cutoff::FIRST, Relationship::Second),
            (cutoff::SECOND, Relationship::Third),
            (cutoff::THIRD, Relationship::Fourth),
            (cutoff::FOURTH, Relationship::Unrelated),
        ];
        for (k, expected) in cases {
            let below = just_below(k);
            assert!(below < k);
            assert_eq!(classify(below, NO_IBS0, PO_THR), expected, "just below {k}");
        }
    }

    #[test]
    fn parent_offspring_split_is_strictly_below_the_threshold() {
        let first_degree = 0.25;
        assert_eq!(
            classify(first_degree, PO_THR, PO_THR),
            Relationship::FullSib,
            "exactly on the IBS0 threshold is a full sib, not PO"
        );
        assert_eq!(
            classify(first_degree, just_below(PO_THR), PO_THR),
            Relationship::ParentOffspring
        );
        assert_eq!(
            classify(first_degree, 0.0, PO_THR),
            Relationship::ParentOffspring
        );
        assert_eq!(classify(first_degree, 0.9, PO_THR), Relationship::FullSib);
    }

    #[test]
    fn the_ibs0_split_only_applies_to_first_degree_pairs() {
        // A high IBS0 proportion must not perturb any other class.
        assert_eq!(classify(0.5, 0.9, PO_THR), Relationship::DupMz);
        assert_eq!(classify(0.12, 0.9, PO_THR), Relationship::Second);
        assert_eq!(classify(0.06, 0.9, PO_THR), Relationship::Third);
        assert_eq!(classify(0.03, 0.9, PO_THR), Relationship::Fourth);
        assert_eq!(classify(0.01, 0.9, PO_THR), Relationship::Unrelated);
    }

    #[test]
    fn non_finite_input_is_total_and_lands_somewhere_sensible() {
        // A degenerate pair (nan kinship) is one the reference drops from .kin outright;
        // classify must not panic and must not claim relatedness.
        assert_eq!(classify(f64::NAN, NO_IBS0, PO_THR), Relationship::Unrelated);
        assert_eq!(
            classify(f64::NEG_INFINITY, NO_IBS0, PO_THR),
            Relationship::Unrelated
        );
        assert_eq!(
            classify(f64::INFINITY, NO_IBS0, PO_THR),
            Relationship::DupMz
        );
        // A nan IBS0 proportion fails `<`, so the pair is a full sib rather than PO.
        assert_eq!(classify(0.25, f64::NAN, PO_THR), Relationship::FullSib);
    }

    #[test]
    fn labels_and_degrees_line_up_with_the_classes() {
        assert_eq!(classify(0.4, NO_IBS0, PO_THR).label(), "Dup/MZ");
        assert_eq!(classify(0.25, 0.0, PO_THR).label(), "PO");
        assert_eq!(classify(0.25, 0.5, PO_THR).label(), "FS");
        assert_eq!(classify(0.25, 0.5, PO_THR).degree(), Some(1));
        assert_eq!(classify(0.25, 0.0, PO_THR).degree(), Some(1));
        assert_eq!(classify(0.1, NO_IBS0, PO_THR).degree(), Some(2));
        assert_eq!(classify(0.0, NO_IBS0, PO_THR).degree(), None);
    }

    /// The bisection results from the reference binary. Each pair of kinships brackets one
    /// boundary to six decimals; the values are exactly the ones the binary was run on.
    #[test]
    fn measured_reference_boundaries_reproduce_under_the_reference_cutoffs() {
        let c = &Cutoffs::REFERENCE;
        let cases = [
            (0.353_553, Relationship::ParentOffspring),
            (0.353_554, Relationship::DupMz),
            (0.176_776, Relationship::Second),
            (0.176_777, Relationship::ParentOffspring),
            (0.088_388, Relationship::Third),
            (0.088_389, Relationship::Second),
        ];
        for (k, expected) in cases {
            assert_eq!(classify_with(k, NO_IBS0, PO_THR, c), expected, "{k}");
        }
    }

    #[test]
    fn published_and_reference_cutoffs_disagree_only_inside_the_known_bands() {
        let c = &Cutoffs::REFERENCE;
        // Inside a band: the reference promotes, the published constants do not.
        assert_eq!(classify(0.088_389, NO_IBS0, PO_THR), Relationship::Third);
        assert_eq!(
            classify_with(0.088_389, NO_IBS0, PO_THR, c),
            Relationship::Second
        );
        // Outside every band the two agree.
        for k in [0.5, 0.36, 0.2, 0.1, 0.07, 0.05, 0.03, 0.01, 0.0] {
            assert_eq!(
                classify(k, NO_IBS0, PO_THR),
                classify_with(k, NO_IBS0, PO_THR, c),
                "{k}"
            );
        }
        // Every REFERENCE cutoff sits strictly below its PUBLISHED counterpart, so the
        // disagreement is always "reference promotes", never the other way round.
        let (r, p) = (Cutoffs::REFERENCE, Cutoffs::PUBLISHED);
        for (reference, published) in [
            (r.dup_mz, p.dup_mz),
            (r.first, p.first),
            (r.second, p.second),
            (r.third, p.third),
            (r.fourth, p.fourth),
        ] {
            assert!(reference < published, "{reference} vs {published}");
        }
    }

    #[test]
    fn reference_cutoffs_are_the_powers_they_claim_to_be() {
        let c = Cutoffs::REFERENCE;
        for (got, k) in [
            (c.dup_mz, 3.0f64),
            (c.first, 5.0),
            (c.second, 7.0),
            (c.third, 9.0),
            (c.fourth, 11.0),
        ] {
            assert_eq!(got, 2f64.powf(-k / 2.0), "2^(-{k}/2)");
        }
    }

    // ---------------------------------------------------------------------------------
    // Pedigree
    // ---------------------------------------------------------------------------------

    fn ped_of(fid: &str, rows: &[(&str, &str, &str)]) -> Pedigree {
        let samples: Vec<Sample> = rows
            .iter()
            .map(|(iid, pat, mat)| Sample {
                fid: fid.to_string(),
                iid: (*iid).to_string(),
                pat: (*pat).to_string(),
                mat: (*mat).to_string(),
                sex: 0,
                pheno: "1".to_string(),
            })
            .collect();
        Pedigree::from_samples(&samples)
    }

    /// A three-generation outbred pedigree covering every relationship the `.kin` columns
    /// are specified against.
    ///
    /// ```text
    ///     gpa x gma -> dad, aunt          (dad and aunt are full sibs)
    ///     dad x mom -> kid1, kid2         (full sibs)
    ///     dad x mom2 -> half              (half sib of kid1/kid2)
    ///     aunt x unc -> cous              (first cousin of kid1/kid2)
    /// ```
    const OUTBRED: [(&str, &str, &str); 10] = [
        ("gpa", "0", "0"),
        ("gma", "0", "0"),
        ("mom", "0", "0"),
        ("mom2", "0", "0"),
        ("unc", "0", "0"),
        ("dad", "gpa", "gma"),
        ("aunt", "gpa", "gma"),
        ("kid1", "dad", "mom"),
        ("kid2", "dad", "mom"),
        ("half", "dad", "mom2"),
    ];

    fn idx(rows: &[(&str, &str, &str)], iid: &str) -> usize {
        rows.iter().position(|r| r.0 == iid).expect("unknown iid")
    }

    #[test]
    fn canonical_relationships_give_the_documented_phi_and_z0() {
        let mut rows = OUTBRED.to_vec();
        rows.push(("cous", "unc", "aunt"));
        let p = ped_of("f", &rows);
        let i = |name: &str| idx(&rows, name);

        // (a, b, phi, z0)
        let cases: [(&str, &str, f64, f64); 9] = [
            ("dad", "kid1", 0.25, 0.000),    // parent-offspring
            ("mom", "kid1", 0.25, 0.000),    // parent-offspring
            ("kid1", "kid2", 0.25, 0.250),   // full sibs
            ("kid1", "half", 0.125, 0.500),  // half sibs
            ("gpa", "kid1", 0.125, 0.500),   // grandparent
            ("aunt", "kid1", 0.125, 0.500),  // avuncular
            ("kid1", "cous", 0.0625, 0.750), // first cousins
            ("dad", "mom", 0.0, 1.000),      // unrelated within a family
            ("gpa", "unc", 0.0, 1.000),      // unrelated founders
        ];
        for (a, b, phi, z0) in cases {
            let (x, y) = (i(a), i(b));
            assert_eq!(p.kinship(x, y), phi, "phi({a},{b})");
            assert_eq!(p.z0(x, y), z0, "z0({a},{b})");
            // symmetric
            assert_eq!(p.kinship(y, x), phi, "phi({b},{a})");
            assert_eq!(p.z0(y, x), z0, "z0({b},{a})");
        }
        // self
        assert_eq!(p.kinship(i("kid1"), i("kid1")), 0.5);
        assert_eq!(p.z0(i("kid1"), i("kid1")), 0.0);
    }

    /// First-cousin mating. `inb` is the child of two first cousins, so `F = 1/16` and the
    /// recurrence produces values no lookup table would hold. Every number below was
    /// printed by the reference binary for this exact pedigree.
    ///
    /// ```text
    ///     gpa x gma -> dad, aunt
    ///     dad x mom -> kid1                (kid1 and cous are first cousins)
    ///     unc x aunt -> cous
    ///     cous x kid1 -> inb, inb2         (inb, inb2 are full sibs, F = 0.0625)
    /// ```
    const INBRED: [(&str, &str, &str); 12] = [
        ("gpa", "0", "0"),
        ("gma", "0", "0"),
        ("mom", "0", "0"),
        ("unc", "0", "0"),
        ("mom2", "0", "0"),
        ("dad", "gpa", "gma"),
        ("aunt", "gpa", "gma"),
        ("kid1", "dad", "mom"),
        ("cous", "unc", "aunt"),
        ("half", "dad", "mom2"),
        ("inb", "cous", "kid1"),
        ("inb2", "cous", "kid1"),
    ];

    #[test]
    fn first_cousin_mating_matches_the_reference() {
        let rows = INBRED.to_vec();
        let p = ped_of("ped", &rows);
        let i = |name: &str| idx(&rows, name);
        let mut cache = KinshipCache::default();

        // Inbreeding coefficient of the child of first cousins is 1/16, so phi(inb, inb)
        // is (1 + 1/16)/2 rather than 1/2.
        assert_eq!(p.kinship(i("cous"), i("kid1")), 0.0625);
        assert_eq!(p.kinship(i("inb"), i("inb")), 0.53125);

        // (a, b, phi, z0) exactly as printed by the reference for this pedigree.
        let cases: [(&str, &str, f64, f64); 8] = [
            ("kid1", "inb", 0.28125, 0.000), // parent-offspring, inbred child
            ("cous", "inb", 0.28125, 0.000), // parent-offspring, inbred child
            ("inb", "inb2", 0.28125, 0.250), // full sibs, inbred
            ("dad", "inb", 0.1875, 0.250),   // grandparent + great-uncle path
            ("half", "inb", 0.09375, 0.625), // no canonical class at all
            ("gpa", "inb", 0.125, 0.500),
            ("mom", "inb", 0.125, 0.500),
            ("kid1", "cous", 0.0625, 0.750),
        ];
        for (a, b, phi, z0) in cases {
            let (x, y) = (i(a), i(b));
            assert_eq!(pedigree_kinship(&p, &mut cache, x, y), phi, "phi({a},{b})");
            assert_eq!(pedigree_z0(&p, &mut cache, x, y), z0, "z0({a},{b})");
        }
    }

    /// Deeper consanguinity, again all measured. These are the cases that pin the two
    /// pieces of `pedigree_z0` that are not simply `1 - 4*phi`.
    #[test]
    fn deeper_inbreeding_matches_the_reference() {
        // Full-sib mating twice over: sibs -> sibs -> child.
        let rows = vec![
            ("a1", "0", "0"),
            ("a2", "0", "0"),
            ("b1", "a1", "a2"),
            ("b2", "a1", "a2"),
            ("k1", "b1", "b2"),
            ("k2", "b1", "b2"),
            ("l1", "k1", "k2"),
            ("h1", "a1", "0"),
            ("h2", "a1", "0"),
        ];
        let p = ped_of("px", &rows);
        let i = |name: &str| idx(&rows, name);

        // Full sibs whose parents are themselves full sibs: phi 0.375, and the reference
        // reports Z0 0.000 rather than the usual full-sib 0.250.
        assert_eq!(p.kinship(i("k1"), i("k2")), 0.375);
        assert_eq!(p.z0(i("k1"), i("k2")), 0.0);
        // Parent and inbred child.
        assert_eq!(p.kinship(i("b1"), i("k1")), 0.375);
        assert_eq!(p.z0(i("b1"), i("k1")), 0.0);
        assert_eq!(p.kinship(i("k1"), i("l1")), 0.5);
        assert_eq!(p.z0(i("k1"), i("l1")), 0.0);
        // Shared father, both mothers unknown: half sibs, NOT full sibs.
        assert!(!p.is_full_sib_pair(i("h1"), i("h2")));
        assert_eq!(p.kinship(i("h1"), i("h2")), 0.125);
        assert_eq!(p.z0(i("h1"), i("h2")), 0.5);
    }

    /// The full-sibling `Z0` limit, bracketed against the reference. Sibling pairs were
    /// built at each `phi` below and the reference printed exactly these `Z0` values.
    #[test]
    fn full_sib_z0_switches_off_above_the_measured_limit() {
        // Reference: 0.25, 0.28125, 0.2890625, 0.296875 -> 0.250;
        //            0.3046875, 0.3125, 0.3515625, 0.375 -> 0.000.
        for phi in [0.25, 0.28125, 0.2890625, 0.296875] {
            assert!(phi <= FS_Z0_KINSHIP_LIMIT, "{phi} should keep Z0 = 0.25");
        }
        for phi in [0.3046875, 0.3125, 0.3515625, 0.375] {
            assert!(phi > FS_Z0_KINSHIP_LIMIT, "{phi} should drop Z0 to 0");
        }
    }

    /// `1 - 4*phi` is reported down to 0.25 and then falls off a cliff to zero. Measured:
    /// `phi` 0.1875 printed `0.250`, while 0.203125 and 0.21875 printed `0.000` rather
    /// than `0.187` and `0.125`.
    #[test]
    fn z0_floors_values_that_fall_between_zero_and_the_full_sib_value() {
        let rows = vec![
            ("gpa", "0", "0"),
            ("gma", "0", "0"),
            ("mom", "0", "0"),
            ("unc", "0", "0"),
            ("dad", "gpa", "gma"),
            ("aunt", "gpa", "gma"),
            ("kid1", "dad", "mom"),
            ("cous", "unc", "aunt"),
            ("inb", "cous", "kid1"),
            ("gk", "dad", "inb"),
        ];
        let p = ped_of("pz", &rows);
        let i = |name: &str| idx(&rows, name);

        // phi 0.1875 -> 1 - 4*phi = 0.25, which is reported as-is.
        assert_eq!(p.kinship(i("dad"), i("inb")), 0.1875);
        assert_eq!(p.z0(i("dad"), i("inb")), 0.25);
        assert_eq!(p.kinship(i("gma"), i("gk")), 0.1875);
        assert_eq!(p.z0(i("gma"), i("gk")), 0.25);
        // phi 0.203125 -> 1 - 4*phi = 0.1875, which is floored to 0.
        assert_eq!(p.kinship(i("cous"), i("gk")), 0.203125);
        assert_eq!(p.z0(i("cous"), i("gk")), 0.0);
        // phi 0.21875 -> 1 - 4*phi = 0.125, also floored to 0.
        assert_eq!(p.kinship(i("aunt"), i("gk")), 0.21875);
        assert_eq!(p.z0(i("aunt"), i("gk")), 0.0);
        // Well past the floor, where 1 - 4*phi is negative anyway.
        assert_eq!(p.kinship(i("kid1"), i("gk")), 0.265625);
        assert_eq!(p.z0(i("kid1"), i("gk")), 0.0);
        assert_eq!(p.kinship(i("inb"), i("gk")), 0.359375);
        assert_eq!(p.z0(i("inb"), i("gk")), 0.0);
        assert_eq!(p.kinship(i("mom"), i("gk")), 0.0625);
        assert_eq!(p.z0(i("mom"), i("gk")), 0.75);
    }

    #[test]
    fn declaration_order_does_not_change_the_answer() {
        // Same pedigree as OUTBRED but with children listed before their parents, which a
        // hand-edited .fam can easily do.
        let forward = OUTBRED.to_vec();
        let mut reversed = OUTBRED.to_vec();
        reversed.reverse();
        let pf = ped_of("f", &forward);
        let pr = ped_of("f", &reversed);
        for (a, b) in [
            ("kid1", "kid2"),
            ("gpa", "kid1"),
            ("dad", "kid1"),
            ("aunt", "half"),
        ] {
            assert_eq!(
                pf.kinship(idx(&forward, a), idx(&forward, b)),
                pr.kinship(idx(&reversed, a), idx(&reversed, b)),
                "phi({a},{b})"
            );
            assert_eq!(
                pf.z0(idx(&forward, a), idx(&forward, b)),
                pr.z0(idx(&reversed, a), idx(&reversed, b)),
                "z0({a},{b})"
            );
        }
    }

    #[test]
    fn parents_are_resolved_within_the_family_only() {
        let samples = vec![
            Sample {
                fid: "f1".into(),
                iid: "x".into(),
                pat: "0".into(),
                mat: "0".into(),
                sex: 1,
                pheno: "1".into(),
            },
            Sample {
                fid: "f2".into(),
                iid: "y".into(),
                // "x" exists, but in another family: not a parent.
                pat: "x".into(),
                mat: "0".into(),
                sex: 1,
                pheno: "1".into(),
            },
        ];
        let p = Pedigree::from_samples(&samples);
        assert_eq!(p.parents(1), (None, None));
        assert_eq!(p.kinship(0, 1), 0.0);
        assert_eq!(p.z0(0, 1), 1.0);
    }

    #[test]
    fn cycles_are_broken_instead_of_hanging() {
        // c is its own grandparent, and d is its own father.
        let rows = vec![
            ("a", "0", "0"),
            ("b", "c", "a"),
            ("c", "b", "a"),
            ("d", "d", "0"),
        ];
        let p = ped_of("cyc", &rows);
        assert_eq!(
            p.parents(3),
            (None, None),
            "self-parent edge must be dropped"
        );
        // Whatever the exact numbers, the traversal must terminate and stay in range.
        for a in 0..p.len() {
            for b in 0..p.len() {
                let phi = p.kinship(a, b);
                assert!(
                    phi.is_finite() && (0.0..=1.0).contains(&phi),
                    "phi({a},{b}) = {phi}"
                );
                let z0 = p.z0(a, b);
                assert!(
                    z0.is_finite() && (0.0..=1.0).contains(&z0),
                    "z0({a},{b}) = {z0}"
                );
            }
        }
    }

    #[test]
    fn empty_and_single_sample_pedigrees_are_well_formed() {
        let p = Pedigree::from_samples(&[]);
        assert!(p.is_empty());
        assert_eq!(p.len(), 0);

        let one = ped_of("f", &[("solo", "0", "0")]);
        assert_eq!(one.len(), 1);
        assert!(!one.is_empty());
        assert_eq!(one.kinship(0, 0), 0.5);
        assert!(!one.is_full_sib_pair(0, 0));
    }

    #[test]
    fn a_shared_cache_gives_the_same_values_as_fresh_ones() {
        let rows = INBRED.to_vec();
        let p = ped_of("ped", &rows);
        let mut cache = KinshipCache::default();
        for a in 0..p.len() {
            for b in 0..p.len() {
                assert_eq!(pedigree_kinship(&p, &mut cache, a, b), p.kinship(a, b));
                assert_eq!(pedigree_z0(&p, &mut cache, a, b), p.z0(a, b));
            }
        }
        cache.clear();
        assert_eq!(
            pedigree_kinship(&p, &mut cache, idx(&rows, "inb"), idx(&rows, "inb2")),
            0.28125
        );
    }

    /// The pedigree classes feed straight into the same cutoffs `classify` uses, so the
    /// expected kinships must land in the classes their names promise.
    #[test]
    fn pedigree_kinships_classify_as_their_relationship_names() {
        assert_eq!(classify(0.5, 0.0, PO_THR), Relationship::DupMz);
        assert_eq!(classify(0.25, 0.0, PO_THR), Relationship::ParentOffspring);
        assert_eq!(classify(0.25, 0.02, 0.01), Relationship::FullSib);
        assert_eq!(classify(0.125, NO_IBS0, PO_THR), Relationship::Second);
        assert_eq!(classify(0.0625, NO_IBS0, PO_THR), Relationship::Third);
        assert_eq!(classify(0.03125, NO_IBS0, PO_THR), Relationship::Fourth);
        assert_eq!(classify(0.0, NO_IBS0, PO_THR), Relationship::Unrelated);
    }
}
