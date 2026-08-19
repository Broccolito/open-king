//! Owned, typed entry point for relatedness analysis.
//!
//! [`Bundle`] validates a loaded PLINK fileset once. [`Bundle::relatedness`] then returns
//! every pair as structured Rust data: exact integer counts, the KING estimators,
//! pedigree expectations, inferred relationship, and optional IBD-segment metrics. No
//! console output or text-file round trip is involved.

use std::fmt;
use std::path::Path;

use king_io::{bed, Fileset, IoError, Sample, VariantFilter};

use crate::counts;
use crate::ibdseg::{self, PairSegments, Usable, DEFAULT_SEGLENGTH_BP};
use crate::infer::{self, Cutoffs, Pedigree};
use crate::kinship;
use crate::{PairCounts, Relationship, Scope};

/// A validated genotype/map/sample bundle ready for relatedness analysis.
#[derive(Debug, Clone)]
pub struct Bundle {
    fileset: Fileset,
}

impl Bundle {
    /// Validate and adopt an already loaded fileset.
    pub fn new(fileset: Fileset) -> Result<Self, BundleError> {
        validate(&fileset)?;
        Ok(Self { fileset })
    }

    /// Load the autosomal PLINK fileset at `prefix_or_bed` and validate it.
    pub fn from_plink(prefix_or_bed: impl AsRef<Path>) -> Result<Self, BundleError> {
        let fileset =
            bed::read_fileset(prefix_or_bed.as_ref(), VariantFilter::Autosomes, None, None)?;
        Self::new(fileset)
    }

    /// Borrow the validated low-level fileset.
    pub fn fileset(&self) -> &Fileset {
        &self.fileset
    }

    /// Recover the owned low-level fileset.
    pub fn into_fileset(self) -> Fileset {
        self.fileset
    }

    /// Compute an owned relatedness report without going through the CLI or text files.
    pub fn relatedness(&self, options: &RelatednessOptions) -> RelatednessReport {
        analyze(self, options)
    }
}

impl TryFrom<Fileset> for Bundle {
    type Error = BundleError;

    fn try_from(fileset: Fileset) -> Result<Self, Self::Error> {
        Self::new(fileset)
    }
}

/// How the IBS0 boundary between parent–offspring and full siblings is selected.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PoThreshold {
    /// Match the relatedness workflow: use half the mean IBS0 proportion among declared
    /// full-sibling pairs whose estimate is first-degree. If there are no such anchors,
    /// only an exactly zero IBS0 proportion is parent–offspring.
    FromPedigree,
    /// Use this caller-supplied, strict upper bound (`IBS0 < threshold`).
    Fixed(f64),
}

/// Controls for the typed relatedness analysis.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RelatednessOptions {
    /// Relationship-class boundaries. The default uses the measured KING 2.3.2 powers.
    pub cutoffs: Cutoffs,
    /// Parent–offspring/full-sibling calibration policy.
    pub po_threshold: PoThreshold,
    /// Minimum called IBD segment length in base pairs. `None` skips segment scanning.
    pub segment_length_bp: Option<i64>,
}

impl Default for RelatednessOptions {
    fn default() -> Self {
        Self {
            cutoffs: Cutoffs::REFERENCE,
            po_threshold: PoThreshold::FromPedigree,
            segment_length_bp: Some(DEFAULT_SEGLENGTH_BP),
        }
    }
}

/// All scalar genotype statistics exposed by the relatedness and IBS reports.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairStatistics {
    pub kinship: f64,
    pub ibs_mean: f64,
    pub distance: f64,
    pub het_concordance: f64,
    pub het_2_given_1: f64,
    pub het_1_given_2: f64,
    pub hom_concordance: f64,
    pub het_het_proportion: f64,
    pub ibs0_proportion: f64,
    /// The `HomIBS0` statistic used by the close-relative report.
    pub hom_ibs0_proportion: f64,
}

/// Structured IBD-segment results for one pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SegmentStatistics {
    /// Raw called lengths, including the longest-call reporting gate.
    pub counts: PairSegments,
    pub denominator_bp: i64,
    pub ibd1_proportion: f64,
    pub ibd2_proportion: f64,
    pub proportion_ibd: f64,
    /// KING's `.seg` display formula, which rounds π1/π2 before combining them.
    pub displayed_proportion_ibd: f64,
    /// Whether the pair survives KING's fixed 10 Mb long-segment gate.
    pub reported: bool,
    /// Segment-based relationship label used in `.seg`/`.kin` output.
    pub inference: &'static str,
}

/// One unordered sample pair. `first < second` always holds.
#[derive(Debug, Clone, PartialEq)]
pub struct RelatedPair {
    pub first: usize,
    pub second: usize,
    pub scope: Scope,
    pub counts: PairCounts,
    pub statistics: PairStatistics,
    pub pedigree_kinship: f64,
    pub pedigree_z0: f64,
    pub relationship: Relationship,
    pub segments: Option<SegmentStatistics>,
}

/// Owned result of analyzing every unordered sample pair.
#[derive(Debug, Clone, PartialEq)]
pub struct RelatednessReport {
    pub samples: Vec<Sample>,
    pub pairs: Vec<RelatedPair>,
    /// Calibrated cutoff, or `None` when no declared full-sibling anchors existed.
    pub po_threshold: Option<f64>,
    /// IBD denominator for the retained map, or `None` when scanning was disabled.
    pub segment_denominator_bp: Option<i64>,
    /// The map intervals used in the segment denominator.
    pub usable_segments: Vec<Usable>,
}

impl RelatednessReport {
    /// Find an unordered pair by sample indices.
    pub fn pair(&self, a: usize, b: usize) -> Option<&RelatedPair> {
        let (first, second) = if a < b { (a, b) } else { (b, a) };
        self.pairs
            .iter()
            .find(|pair| pair.first == first && pair.second == second)
    }
}

fn analyze(bundle: &Bundle, options: &RelatednessOptions) -> RelatednessReport {
    let fileset = bundle.fileset();
    let n = fileset.samples.len();
    let pair_indices: Vec<(usize, usize)> = (0..n)
        .flat_map(|first| ((first + 1)..n).map(move |second| (first, second)))
        .collect();
    let pair_counts = counts::all_pairs(&fileset.genotypes, &pair_indices);
    let pedigree = Pedigree::from_samples(&fileset.samples);

    let preliminary: Vec<(Scope, PairStatistics)> = pair_indices
        .iter()
        .zip(&pair_counts)
        .map(|(&(first, second), counts)| {
            let scope = if fileset.samples[first].fid == fileset.samples[second].fid {
                Scope::WithinFamily
            } else {
                Scope::BetweenFamily
            };
            (scope, statistics(counts, scope))
        })
        .collect();

    let calibrated_po = match options.po_threshold {
        PoThreshold::Fixed(value) => Some(value),
        PoThreshold::FromPedigree => {
            let mut sum = 0.0;
            let mut anchors = 0usize;
            for ((&(first, second), (_, stats)), _) in
                pair_indices.iter().zip(&preliminary).zip(&pair_counts)
            {
                if pedigree.is_full_sib_pair(first, second)
                    && stats.kinship >= options.cutoffs.first
                    && stats.kinship < options.cutoffs.dup_mz
                {
                    sum += stats.ibs0_proportion;
                    anchors += 1;
                }
            }
            (anchors > 0).then(|| 0.5 * sum / anchors as f64)
        }
    };
    // No anchors is not equivalent to a strict cutoff of zero: KING calls exact IBS0=0
    // parent–offspring in that case. The smallest positive f64 encodes that boundary.
    let effective_po = calibrated_po.unwrap_or(f64::MIN_POSITIVE);

    let segment_context = options.segment_length_bp.map(|minimum_bp| {
        let (chromosomes, positions) = retained_map(fileset);
        let usable = ibdseg::usable_segments(&chromosomes, &positions);
        let denominator = ibdseg::denominator(&usable, &positions);
        (minimum_bp, positions, usable, denominator)
    });

    let mut pedigree_cache = infer::KinshipCache::default();
    let pairs = pair_indices
        .into_iter()
        .zip(pair_counts)
        .zip(preliminary)
        .map(|(((first, second), counts), (scope, statistics))| {
            let relationship = infer::classify_with(
                statistics.kinship,
                statistics.ibs0_proportion,
                effective_po,
                &options.cutoffs,
            );
            let segments =
                segment_context
                    .as_ref()
                    .map(|(minimum_bp, positions, usable, denominator)| {
                        let raw = ibdseg::pair_segments(
                            &fileset.genotypes,
                            positions,
                            usable,
                            first,
                            second,
                            *minimum_bp,
                        );
                        let ibd1 = raw.ibd1_seg(*denominator);
                        let ibd2 = raw.ibd2_seg(*denominator);
                        let prop = raw.prop_ibd(*denominator);
                        SegmentStatistics {
                            counts: raw,
                            denominator_bp: *denominator,
                            ibd1_proportion: ibd1,
                            ibd2_proportion: ibd2,
                            proportion_ibd: prop,
                            displayed_proportion_ibd: ibdseg::seg_prop_ibd(ibd1, ibd2),
                            reported: raw.reported(),
                            inference: ibdseg::inf_type(ibd1, ibd2, prop),
                        }
                    });
            let pedigree_kinship =
                infer::pedigree_kinship(&pedigree, &mut pedigree_cache, first, second);
            let pedigree_z0 = infer::pedigree_z0(&pedigree, &mut pedigree_cache, first, second);
            RelatedPair {
                first,
                second,
                scope,
                counts,
                statistics,
                pedigree_kinship,
                pedigree_z0,
                relationship,
                segments,
            }
        })
        .collect();

    let (segment_denominator_bp, usable_segments) = match segment_context {
        Some((_, _, usable, denominator)) => (Some(denominator), usable),
        None => (None, Vec::new()),
    };

    RelatednessReport {
        samples: fileset.samples.clone(),
        pairs,
        po_threshold: calibrated_po,
        segment_denominator_bp,
        usable_segments,
    }
}

fn statistics(counts: &PairCounts, scope: Scope) -> PairStatistics {
    PairStatistics {
        kinship: kinship::kinship(counts, scope),
        ibs_mean: kinship::ibs_mean(counts),
        distance: kinship::dist(counts),
        het_concordance: kinship::het_concordance(counts),
        het_2_given_1: kinship::het_2given1(counts),
        het_1_given_2: kinship::het_1given2(counts),
        hom_concordance: kinship::hom_concordance(counts),
        het_het_proportion: kinship::het_het_prop(counts),
        ibs0_proportion: kinship::ibs0_prop(counts),
        hom_ibs0_proportion: f64::from(counts.ibs0) / f64::from(counts.hom_a1_union),
    }
}

fn retained_map(fileset: &Fileset) -> (Vec<i64>, Vec<i64>) {
    fileset
        .kept
        .iter()
        .map(|&index| {
            let variant = &fileset.variants[index];
            (
                i64::from(
                    bed::king_chrom_code(&variant.chrom)
                        .expect("Bundle validation guarantees a recognized chromosome"),
                ),
                variant.bp,
            )
        })
        .unzip()
}

fn validate(fileset: &Fileset) -> Result<(), BundleError> {
    king_io::fam::check_duplicates(&fileset.samples)?;
    let g = &fileset.genotypes;
    if g.n_samples != fileset.samples.len() {
        return Err(BundleError::SampleCount {
            samples: fileset.samples.len(),
            genotypes: g.n_samples,
        });
    }
    if g.n_variants != fileset.kept.len() {
        return Err(BundleError::VariantCount {
            kept: fileset.kept.len(),
            genotypes: g.n_variants,
        });
    }
    if g.n_variants > u32::MAX as usize {
        return Err(BundleError::TooManyVariants {
            variants: g.n_variants,
        });
    }
    if g.plane0.len() != g.n_samples || g.plane1.len() != g.n_samples {
        return Err(BundleError::PlaneCount {
            samples: g.n_samples,
            plane0: g.plane0.len(),
            plane1: g.plane1.len(),
        });
    }
    let words = g.words_per_sample();
    for sample in 0..g.n_samples {
        for (plane, values) in [(0u8, &g.plane0[sample]), (1u8, &g.plane1[sample])] {
            if values.len() != words {
                return Err(BundleError::PlaneWords {
                    sample,
                    plane,
                    expected: words,
                    found: values.len(),
                });
            }
            if g.n_variants % 64 != 0 {
                let tail = !0u64 << (g.n_variants % 64);
                if values.last().is_some_and(|word| word & tail != 0) {
                    return Err(BundleError::DirtyTail { sample, plane });
                }
            }
        }
    }

    let mut previous: Option<(u8, i64, usize)> = None;
    for (position, &index) in fileset.kept.iter().enumerate() {
        let Some(variant) = fileset.variants.get(index) else {
            return Err(BundleError::KeptIndex {
                position,
                index,
                variants: fileset.variants.len(),
            });
        };
        let Some(chromosome) = bed::king_chrom_code(&variant.chrom) else {
            return Err(BundleError::UnknownChromosome {
                variant: index,
                chromosome: variant.chrom.clone(),
            });
        };
        if let Some((previous_chromosome, previous_bp, previous_index)) = previous {
            if chromosome < previous_chromosome
                || (chromosome == previous_chromosome && variant.bp < previous_bp)
            {
                return Err(BundleError::MapOrder {
                    previous: previous_index,
                    current: index,
                });
            }
        }
        previous = Some((chromosome, variant.bp, index));
    }
    Ok(())
}

/// A malformed or unreadable [`Bundle`].
#[derive(Debug)]
pub enum BundleError {
    Io(IoError),
    SampleCount {
        samples: usize,
        genotypes: usize,
    },
    VariantCount {
        kept: usize,
        genotypes: usize,
    },
    TooManyVariants {
        variants: usize,
    },
    PlaneCount {
        samples: usize,
        plane0: usize,
        plane1: usize,
    },
    PlaneWords {
        sample: usize,
        plane: u8,
        expected: usize,
        found: usize,
    },
    DirtyTail {
        sample: usize,
        plane: u8,
    },
    KeptIndex {
        position: usize,
        index: usize,
        variants: usize,
    },
    UnknownChromosome {
        variant: usize,
        chromosome: String,
    },
    MapOrder {
        previous: usize,
        current: usize,
    },
}

impl From<IoError> for BundleError {
    fn from(error: IoError) -> Self {
        Self::Io(error)
    }
}

impl fmt::Display for BundleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Io(error) => error.fmt(f),
            Self::SampleCount { samples, genotypes } => write!(
                f,
                "sample table has {samples} rows but genotypes declare {genotypes} samples"
            ),
            Self::VariantCount { kept, genotypes } => write!(
                f,
                "kept map has {kept} variants but genotypes declare {genotypes} variants"
            ),
            Self::TooManyVariants { variants } => {
                write!(f, "{variants} variants exceed the u32 pair-count capacity")
            }
            Self::PlaneCount {
                samples,
                plane0,
                plane1,
            } => write!(
                f,
                "genotypes declare {samples} samples but planes contain {plane0} and {plane1} rows"
            ),
            Self::PlaneWords {
                sample,
                plane,
                expected,
                found,
            } => write!(
                f,
                "sample {sample} plane {plane} has {found} words; expected {expected}"
            ),
            Self::DirtyTail { sample, plane } => {
                write!(f, "sample {sample} plane {plane} has set bits past the map")
            }
            Self::KeptIndex {
                position,
                index,
                variants,
            } => write!(
                f,
                "kept map position {position} points to variant {index}, but only {variants} exist"
            ),
            Self::UnknownChromosome {
                variant,
                chromosome,
            } => write!(
                f,
                "variant {variant} has unrecognized chromosome {chromosome:?}"
            ),
            Self::MapOrder { previous, current } => write!(
                f,
                "retained map regresses between variants {previous} and {current}"
            ),
        }
    }
}

impl std::error::Error for BundleError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Io(error) => Some(error),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use king_io::{Fileset, Genotypes, Sample, Variant};

    use super::*;

    fn sample(fid: &str, iid: &str) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: "0".to_string(),
            mat: "0".to_string(),
            sex: 0,
            pheno: "-9".to_string(),
        }
    }

    fn variants(n: usize) -> Vec<Variant> {
        (0..n)
            .map(|index| Variant {
                chrom: "1".to_string(),
                id: format!("v{index}"),
                cm: 0.0,
                bp: 1 + index as i64 * 50_000,
                a1: "A".to_string(),
                a2: "G".to_string(),
            })
            .collect()
    }

    #[test]
    fn typed_report_exposes_counts_estimators_relationships_and_lookup() {
        let fileset = Fileset {
            samples: vec![sample("F", "A"), sample("F", "B"), sample("G", "C")],
            variants: variants(4),
            genotypes: Genotypes {
                // A/B: A1, het, A2, missing. C reverses the two homozygotes.
                plane0: vec![vec![0b0101], vec![0b0101], vec![0b0101]],
                plane1: vec![vec![0b0011], vec![0b0011], vec![0b0110]],
                n_samples: 3,
                n_variants: 4,
            },
            kept: (0..4).collect(),
        };
        let report = Bundle::new(fileset)
            .expect("valid bundle")
            .relatedness(&RelatednessOptions {
                segment_length_bp: None,
                ..RelatednessOptions::default()
            });

        assert_eq!(report.samples.len(), 3);
        assert_eq!(report.pairs.len(), 3);
        assert_eq!(report.segment_denominator_bp, None);
        let duplicate = report.pair(1, 0).expect("symmetric lookup");
        assert_eq!((duplicate.first, duplicate.second), (0, 1));
        assert_eq!(duplicate.scope, Scope::WithinFamily);
        assert_eq!(duplicate.counts.n_snp, 3);
        assert_eq!(duplicate.counts.het_het, 1);
        assert_eq!(duplicate.statistics.kinship, 0.5);
        assert_eq!(duplicate.relationship, Relationship::DupMz);
        assert!(duplicate.segments.is_none());

        let cross_family = report.pair(0, 2).expect("cross-family pair");
        assert_eq!(cross_family.scope, Scope::BetweenFamily);
        assert_eq!(cross_family.counts.ibs0, 2);
        assert_eq!(cross_family.statistics.kinship, -1.5);
        assert_eq!(cross_family.relationship, Relationship::Unrelated);
    }

    #[test]
    fn plink_entrypoint_matches_a_committed_reference_row_without_text_roundtrip() {
        let bed =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/parity/golden/multifam.bed");
        let report = Bundle::from_plink(bed)
            .expect("committed PLINK fixture loads")
            .relatedness(&RelatednessOptions {
                segment_length_bp: None,
                ..RelatednessOptions::default()
            });

        assert_eq!(report.samples.len(), 20);
        assert_eq!(report.pairs.len(), 190);
        assert!(report.po_threshold.is_some());
        // FAM1 A_C1/A_C2 is the first row of the committed KING 2.3.2 `.kin` and
        // `.ibs` captures. The API reaches the same integers and printed digits directly.
        let pair = report.pair(2, 3).expect("A_C1/A_C2");
        assert_eq!(pair.counts.n_snp, 15_000);
        assert_eq!(pair.counts.ibs0, 230);
        assert_eq!(pair.counts.ibs1(), 3_885);
        assert_eq!(pair.counts.ibs2(), 10_885);
        assert_eq!(pair.counts.het_het, 3_329);
        assert_eq!(pair.counts.hom_hom, 7_786);
        assert_eq!(pair.counts.het_i, 5_257);
        assert_eq!(pair.counts.het_j, 5_286);
        assert_eq!(format!("{:.4}", pair.statistics.ibs_mean), "1.7103");
        assert_eq!(format!("{:.4}", pair.statistics.distance), "0.3203");
        assert_eq!(format!("{:.4}", pair.statistics.kinship), "0.2721");
        assert_eq!(pair.pedigree_z0, 0.25);
        assert_eq!(pair.pedigree_kinship, 0.25);
        assert_eq!(pair.relationship, Relationship::FullSib);
    }

    #[test]
    fn report_wires_the_usable_map_into_segment_metrics() {
        let n = 384usize;
        let words = n.div_ceil(64);
        let all_heterozygous = vec![u64::MAX; words];
        let fileset = Fileset {
            samples: vec![sample("F", "A"), sample("F", "B")],
            variants: variants(n),
            genotypes: Genotypes {
                plane0: vec![vec![0; words], vec![0; words]],
                plane1: vec![all_heterozygous.clone(), all_heterozygous],
                n_samples: 2,
                n_variants: n,
            },
            kept: (0..n).collect(),
        };
        let report = Bundle::new(fileset)
            .expect("valid bundle")
            .relatedness(&RelatednessOptions::default());
        assert_eq!(report.usable_segments.len(), 1);
        assert!(report
            .segment_denominator_bp
            .is_some_and(|bp| bp > 10_000_000));
        let segments = report.pairs[0].segments.expect("segment scan enabled");
        assert_eq!(
            segments.denominator_bp,
            report.segment_denominator_bp.unwrap()
        );
        assert!(segments.ibd2_proportion > 0.0);
        assert!(segments.reported);
    }

    #[test]
    fn malformed_planes_and_dirty_tail_bits_are_typed_errors() {
        let base = Fileset {
            samples: vec![sample("F", "A")],
            variants: variants(1),
            genotypes: Genotypes {
                plane0: vec![vec![1]],
                plane1: vec![vec![1]],
                n_samples: 1,
                n_variants: 1,
            },
            kept: vec![0],
        };

        let mut short = base.clone();
        short.genotypes.plane1[0].clear();
        assert!(matches!(
            Bundle::new(short),
            Err(BundleError::PlaneWords {
                sample: 0,
                plane: 1,
                expected: 1,
                found: 0
            })
        ));

        let mut dirty = base;
        dirty.genotypes.plane0[0][0] |= 1 << 63;
        assert!(matches!(
            Bundle::new(dirty),
            Err(BundleError::DirtyTail {
                sample: 0,
                plane: 0
            })
        ));
    }
}
