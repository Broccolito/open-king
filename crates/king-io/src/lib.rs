//! PLINK 1 binary fileset I/O.
//!
//! Reads `.bed` / `.bim` / `.fam` triples into the packed bit-plane representation the
//! relatedness kernels in `king-core` consume.
//!
//! These type definitions are the contract between `king-io`, `king-core` and
//! `king-cli`. Implementations live in the submodules.

#![forbid(unsafe_code)]

pub mod bed;
pub mod bim;
pub mod fam;

use std::fmt;
use std::path::PathBuf;

/// One row of a `.fam` file.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sample {
    pub fid: String,
    pub iid: String,
    /// Paternal IID, or `"0"` when absent.
    pub pat: String,
    /// Maternal IID, or `"0"` when absent.
    pub mat: String,
    /// 1 = male, 2 = female, 0 = unknown.
    pub sex: u8,
    /// Phenotype, kept verbatim as text so it can be round-tripped.
    pub pheno: String,
}

/// One row of a `.bim` file.
#[derive(Debug, Clone, PartialEq)]
pub struct Variant {
    /// Chromosome exactly as written in the `.bim`.
    pub chrom: String,
    pub id: String,
    /// Genetic distance in centimorgans.
    pub cm: f64,
    /// Base-pair coordinate.
    pub bp: i64,
    /// Allele 1 (PLINK's minor/A1 column). May be `"0"` for a missing allele.
    pub a1: String,
    /// Allele 2 (PLINK's major/A2 column).
    pub a2: String,
}

impl Variant {
    /// Numeric chromosome code, using PLINK's convention
    /// (X = 23, Y = 24, XY = 25, MT = 26). `None` for unrecognised codes.
    pub fn chrom_code(&self) -> Option<u8> {
        chrom_code(&self.chrom)
    }

    /// Whether this variant is on an autosome (codes 1..=22).
    pub fn is_autosome(&self) -> bool {
        matches!(self.chrom_code(), Some(1..=22))
    }
}

/// Map a `.bim` chromosome string to PLINK's numeric code.
pub fn chrom_code(s: &str) -> Option<u8> {
    let t = s.trim();
    let t = t.strip_prefix("chr").unwrap_or(t);
    if let Ok(n) = t.parse::<u8>() {
        return if (1..=26).contains(&n) { Some(n) } else { None };
    }
    match t.to_ascii_uppercase().as_str() {
        "X" => Some(23),
        "Y" => Some(24),
        "XY" | "PAR" => Some(25),
        "MT" | "M" => Some(26),
        _ => None,
    }
}

/// Genotypes stored as two bit planes per sample, SNP-major within a sample.
///
/// For sample `s` and variant index `m`, bit `m` of `plane0[s]` and `plane1[s]` encode:
///
/// | `plane0` | `plane1` | genotype           |
/// |----------|----------|--------------------|
/// | 1        | 1        | homozygous A1/A1   |
/// | 0        | 1        | heterozygous       |
/// | 1        | 0        | homozygous A2/A2   |
/// | 0        | 0        | **missing**        |
///
/// This is chosen so that `plane0` alone is a "is homozygous" mask and
/// `plane0 | plane1` is a "is non-missing" mask, which makes every count the
/// relatedness kernels need a single `popcount` over ANDs — see
/// [`king_core::counts`](../king_core/counts/index.html).
///
/// Bits beyond `n_variants` in the final word of each plane MUST be zero, so that
/// popcounts never need a tail mask.
#[derive(Debug, Clone, Default)]
pub struct Genotypes {
    /// `plane0[s]` is the word array for sample `s`; length `words_per_sample`.
    pub plane0: Vec<Vec<u64>>,
    pub plane1: Vec<Vec<u64>>,
    pub n_samples: usize,
    pub n_variants: usize,
}

impl Genotypes {
    /// Number of 64-bit words per sample: `ceil(n_variants / 64)`.
    ///
    /// This is the figure the reference binary reports as
    /// "Autosome genotypes stored in N words for each of M individuals."
    pub fn words_per_sample(&self) -> usize {
        self.n_variants.div_ceil(64)
    }
}

/// A loaded `.bed`/`.bim`/`.fam` triple.
#[derive(Debug, Clone)]
pub struct Fileset {
    pub samples: Vec<Sample>,
    pub variants: Vec<Variant>,
    /// Genotypes restricted to the variants retained for analysis.
    pub genotypes: Genotypes,
    /// Indices into `variants` that `genotypes` actually carries, in order.
    pub kept: Vec<usize>,
}

/// Which variants to load into the bit planes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariantFilter {
    /// Autosomes only (chromosome codes 1..=22) — the default relatedness path.
    Autosomes,
    /// A single chromosome code, e.g. 23 for X.
    Chromosome(u8),
    /// Everything in the `.bim`.
    All,
}

/// Failure modes shared across the readers, worded to match the reference binary's
/// console messages where those are known.
#[derive(Debug)]
pub enum IoError {
    Open { path: PathBuf, source: std::io::Error },
    Read { path: PathBuf, source: std::io::Error },
    /// `.bed` did not start with the 0x6c 0x1b magic.
    BadMagic { path: PathBuf },
    /// `.bed` mode byte was neither 0x00 nor 0x01.
    BadMode { path: PathBuf, mode: u8 },
    /// Individual-major `.bed` (mode 0x00), which KING does not accept.
    IndividualMajor { path: PathBuf },
    /// `.bed` length disagreed with `n_samples` x `n_variants`.
    SizeMismatch { path: PathBuf, expected: u64, found: u64 },
    /// A line had the wrong number of fields.
    Fields { path: PathBuf, line: usize, expected: usize, found: usize },
    /// A numeric field failed to parse.
    Parse { path: PathBuf, line: usize, field: &'static str, value: String },
    /// Two samples shared an (FID, IID) key.
    DuplicateSample { fid: String, iid: String },
}

impl fmt::Display for IoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            IoError::Open { path, source } => {
                write!(f, "cannot open {}: {source}", path.display())
            }
            IoError::Read { path, source } => {
                write!(f, "error reading {}: {source}", path.display())
            }
            IoError::BadMagic { path } => {
                write!(f, "{} is not a PLINK binary genotype file", path.display())
            }
            IoError::BadMode { path, mode } => {
                write!(f, "{} has unrecognized mode byte {mode:#04x}", path.display())
            }
            IoError::IndividualMajor { path } => write!(
                f,
                "{} is in individual-major order, which is not supported",
                path.display()
            ),
            IoError::SizeMismatch { path, expected, found } => write!(
                f,
                "{} is {found} bytes but {expected} were expected",
                path.display()
            ),
            IoError::Fields { path, line, expected, found } => write!(
                f,
                "{}: line {line} has {found} fields, expected {expected}",
                path.display()
            ),
            IoError::Parse { path, line, field, value } => write!(
                f,
                "{}: line {line}: cannot parse {field} from {value:?}",
                path.display()
            ),
            IoError::DuplicateSample { fid, iid } => {
                write!(f, "duplicate sample {fid} {iid}")
            }
        }
    }
}

impl std::error::Error for IoError {}

pub type Result<T> = std::result::Result<T, IoError>;
