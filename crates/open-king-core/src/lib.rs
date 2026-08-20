//! KING relatedness estimators.
//!
//! Clean-room implementation of the estimators described in Manichaikul et al. 2010
//! (*Bioinformatics* 26:2867–2873) and the KING 2.x documentation. Every formula
//! implemented here is stated, with its verification against the reference binary, in
//! `docs/VERIFIED_FORMULAS.md`.
//!
//! These type definitions are the contract between the submodules and `open-king-cli`.

#![forbid(unsafe_code)]

pub mod counts;
pub mod ibdseg;
pub mod infer;
pub mod kinship;
pub mod report;

pub use report::{
    Bundle, BundleError, PairStatistics, PoThreshold, RelatedPair, RelatednessOptions,
    RelatednessReport, SegmentStatistics,
};

/// Raw pairwise genotype counts for one pair of samples.
///
/// Every count is taken over the **pairwise non-missing** variant set — the set of
/// variants at which *both* samples have a call. `n_snp` is the size of that set.
///
/// Counting `het_i` / `het_j` over each sample's own non-missing set instead is the
/// single most likely way to get subtly wrong kinship; see `docs/VERIFIED_FORMULAS.md`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct PairCounts {
    /// Variants called in both samples (printed as `N_SNP`).
    pub n_snp: u32,
    /// Variants where sample `i` is heterozygous (printed as `N_Het1`).
    pub het_i: u32,
    /// Variants where sample `j` is heterozygous (printed as `N_Het2`).
    pub het_j: u32,
    /// Variants where both are heterozygous (printed as `NHetHet`).
    pub het_het: u32,
    /// Variants where the two are opposite homozygotes (printed as `N_IBS0`).
    pub ibs0: u32,
    /// Variants where both are homozygous (printed as `NHomHom`).
    pub hom_hom: u32,
    /// Variants at which **either** sample is homozygous for A1 — the union, not the
    /// intersection. It is the denominator of `--related`'s `HomIBS0` column and of
    /// nothing else, which is why it has no `.ibs` counterpart to check it against;
    /// `docs/VERIFIED_FORMULAS.md` records the re-derivation from raw `.bed` that fixed
    /// it.
    pub hom_a1_union: u32,
}

impl PairCounts {
    /// `N_IBS1 = het_i + het_j - 2*het_het` — exactly one of the pair is heterozygous.
    pub fn ibs1(&self) -> u32 {
        self.het_i + self.het_j - 2 * self.het_het
    }

    /// `N_IBS2 = n_snp - ibs0 - ibs1`.
    pub fn ibs2(&self) -> u32 {
        self.n_snp - self.ibs0 - self.ibs1()
    }
}

/// Which kinship estimator applies to a pair.
///
/// Selected purely by whether the two samples share an `FID`. The two forms coincide
/// when `het_i == het_j` and diverge otherwise, so this must never be guessed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Scope {
    /// Same family — `.kin`.
    WithinFamily,
    /// Different families — `.kin0`. Robust to population structure.
    BetweenFamily,
}

/// Inferred relationship class, ordered from closest to most distant.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Relationship {
    /// Duplicate sample or monozygotic twin.
    DupMz,
    /// Parent–offspring (1st degree, separated from `FullSib` by IBS0).
    ParentOffspring,
    /// Full siblings (1st degree).
    FullSib,
    Second,
    Third,
    Fourth,
    Unrelated,
}

impl Relationship {
    /// The label the reference binary prints for this class.
    pub fn label(&self) -> &'static str {
        match self {
            Relationship::DupMz => "Dup/MZ",
            Relationship::ParentOffspring => "PO",
            Relationship::FullSib => "FS",
            Relationship::Second => "2nd",
            Relationship::Third => "3rd",
            Relationship::Fourth => "4th",
            Relationship::Unrelated => "UN",
        }
    }

    /// Degree of relationship, or `None` for duplicates and unrelated pairs.
    pub fn degree(&self) -> Option<u8> {
        match self {
            Relationship::ParentOffspring | Relationship::FullSib => Some(1),
            Relationship::Second => Some(2),
            Relationship::Third => Some(3),
            Relationship::Fourth => Some(4),
            _ => None,
        }
    }
}

/// Kinship coefficient boundaries between relationship classes.
///
/// Successive halvings on the `2^(-k/2)` grid, per the paper and the KING manual.
pub mod cutoff {
    /// Above this is a duplicate / MZ twin: `2^(-3/2)`.
    pub const DUP_MZ: f64 = 0.354;
    /// Above this is 1st degree.
    pub const FIRST: f64 = 0.177;
    /// Above this is 2nd degree.
    pub const SECOND: f64 = 0.0884;
    /// Above this is 3rd degree.
    pub const THIRD: f64 = 0.0442;
    /// Above this is 4th degree; below it, unrelated.
    pub const FOURTH: f64 = 0.0221;
}
