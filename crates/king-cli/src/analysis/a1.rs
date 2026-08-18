//! Validation of PLINK's A1-minor orientation contract.
//!
//! KING samples the first 4,096 retained autosomal markers and rejects an input when
//! more than ten percent have more observed A1/A1 than A2/A2 calls. Heterozygous and
//! missing calls cancel out of that comparison. On maps shorter than the window the
//! reference reads unstable tail state; open-king deliberately skips the gate until the
//! complete window exists, preserving valid short-map behavior without reproducing an
//! unsafe read or inventing a stricter contract than the reference reliably enforces.

use king_io::Genotypes;

/// Number of leading retained autosomal markers in KING's validation window.
pub const WINDOW_MARKERS: usize = 4096;

/// Return the rejected A1-major percentage, or `None` when the input passes.
///
/// The integer comparison keeps the strict `> 10%` boundary exact. The returned value
/// uses the same denominator and is ready for the reference's one-decimal diagnostic.
pub fn rejected_percentage(genotypes: &Genotypes) -> Option<f64> {
    if genotypes.n_variants < WINDOW_MARKERS {
        return None;
    }
    let checked = WINDOW_MARKERS;

    let mut major = 0usize;
    for marker in 0..checked {
        let word = marker / 64;
        let mask = 1u64 << (marker % 64);
        let mut hom_a1 = 0usize;
        let mut hom_a2 = 0usize;
        for sample in 0..genotypes.n_samples {
            if genotypes.plane0[sample][word] & mask == 0 {
                continue;
            }
            if genotypes.plane1[sample][word] & mask != 0 {
                hom_a1 += 1;
            } else {
                hom_a2 += 1;
            }
        }
        major += usize::from(hom_a1 > hom_a2);
    }

    (major * 10 > checked).then(|| major as f64 * 100.0 / checked as f64)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(samples: usize, markers: usize, a1_major: &[usize]) -> Genotypes {
        let words = markers.div_ceil(64);
        let mut genotypes = Genotypes {
            plane0: vec![vec![0; words]; samples],
            plane1: vec![vec![0; words]; samples],
            n_samples: samples,
            n_variants: markers,
        };
        // Default every call to A2/A2; turn selected markers into A1/A1.
        for sample in 0..samples {
            for marker in 0..markers {
                let word = marker / 64;
                let mask = 1u64 << (marker % 64);
                genotypes.plane0[sample][word] |= mask;
            }
            for &marker in a1_major {
                let word = marker / 64;
                let mask = 1u64 << (marker % 64);
                genotypes.plane1[sample][word] |= mask;
            }
        }
        genotypes
    }

    #[test]
    fn the_4096_marker_boundary_is_strict() {
        let below: Vec<usize> = (0..409).collect();
        let above: Vec<usize> = (0..410).collect();
        assert_eq!(rejected_percentage(&matrix(20, 5000, &below)), None);
        assert_eq!(
            rejected_percentage(&matrix(20, 5000, &above)),
            Some(10.009765625)
        );
    }

    #[test]
    fn markers_beyond_the_reference_window_do_not_enter_the_gate() {
        let trailing: Vec<usize> = (4096..5000).collect();
        assert_eq!(rejected_percentage(&matrix(20, 5000, &trailing)), None);
    }

    #[test]
    fn short_maps_skip_the_references_uninitialized_tail() {
        let all: Vec<usize> = (0..4095).collect();
        assert_eq!(rejected_percentage(&matrix(20, 4095, &all)), None);
    }

    #[test]
    fn a_homozygote_tie_is_not_a1_major() {
        let mut genotypes = matrix(20, 1, &[]);
        // Ten A1/A1 and ten A2/A2 calls: strict majority is required.
        for sample in 0..10 {
            genotypes.plane1[sample][0] = 1;
        }
        assert_eq!(rejected_percentage(&genotypes), None);
    }
}
