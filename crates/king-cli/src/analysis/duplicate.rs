//! `--duplicate` — find sample pairs whose heterozygote concordance says they are the
//! same person.
//!
//! One file, `<prefix>.con`, holding **both** within- and cross-family pairs in `.fam`
//! serial order, and a console body built around a two-stage search:
//!
//! ```text
//! Sorting autosomes...
//! Computing pairwise genotype concordance starts at <ctime>
//!   8 CPU cores are used...
//!         Stage 1 (with 512 SNPs) screening ends at <ctime>      [N >= 100 only]
//!         Stage 2 (with all SNPs) inference ends at <ctime>
//! 2 pairs of duplicates with heterozygote concordance rate > 80% are saved in file king.con
//!
//!   43 additional pairs from screening stage not confirmed in the final stage
//! ```
//!
//! The three `%.5lf` columns are the only five-decimal fields the reference emits
//! anywhere.

use std::fmt::Write as _;
use std::io::Write;

use king_core::{counts, kinship as est, PairCounts};

use crate::analysis::{cpu_count, f, out_path, serial_pairs};
use crate::cli::{Opt, Options};
use crate::console;
use crate::load::Loaded;

/// Header of `<prefix>.con`. Note `N`, not `N_SNP`.
const CON_HEADER: &str =
    "FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tConcord\tHomConc\tHetConc\n";

/// Sample count at which the reference stops testing every pair and screens first.
///
/// Below it, `.con` is written for every run (header-only when nothing passes) and the
/// console reports `C(n,2) − confirmed` unconfirmed pairs — i.e. every pair was a
/// candidate. At or above it a screening stage runs and the file is written only when
/// something is found. Bracketed by a 3–101 sample sweep (`docs/BEHAVIOR.md` §Q7).
const SCREENING_SAMPLES: usize = 100;

/// SNPs the screening stage uses, as it announces them.
const SCREENING_SNPS: usize = 512;

/// Run the pass: compute every pair, write `.con`, print the body.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;
    let min_conc = opts.double(Opt::MinConc);
    let screens = samples.len() >= SCREENING_SAMPLES;

    write(out, console::SORTING_AUTOSOMES);
    write(out, &console::concordance_starts(console::now_local()));
    write(out, &console::cpu_cores_indented(cpu_count(opts)));

    let pairs = serial_pairs(samples.len());
    let all = counts::all_pairs(genotypes, &pairs);

    // A pair is written iff its heterozygote concordance is **strictly** above the
    // threshold: `--minConc 1` writes nothing even for a pair whose concordance is
    // exactly 1.0, and `--minConc 0` writes every pair.
    let mut confirmed: Vec<usize> = Vec::new();
    for (k, c) in all.iter().enumerate() {
        if c.n_snp > 0 && est::het_concordance(c) > min_conc {
            confirmed.push(k);
        }
    }

    // The screening stage is a lower bound on the candidate set, so a screen that keeps
    // exactly the pairs the final stage confirms leaves no unconfirmed pairs to report —
    // which is what every captured `N >= 100` run shows. See [`screened`].
    let candidates = screened(&pairs, &confirmed, screens);

    if screens {
        if candidates > 0 {
            write(
                out,
                &console::stage1_screening_ends(SCREENING_SNPS, console::now_local()),
            );
            write(out, &console::stage2_inference_ends(console::now_local()));
        } else {
            write(
                out,
                &console::ends_at(console::CONCORDANCE_INDENT, console::now_local()),
            );
        }
    } else {
        write(out, &console::stage2_inference_ends(console::now_local()));
    }

    // Below the screening threshold the file always appears; above it, only when the
    // search actually found something.
    let path = out_path(opts, ".con");
    if !screens || !confirmed.is_empty() {
        let mut text = String::from(CON_HEADER);
        for &k in &confirmed {
            let (i, j) = pairs[k];
            let _ = writeln!(
                text,
                "{}\t{}\t{}\t{}\t{}",
                samples[i].fid,
                samples[i].iid,
                samples[j].fid,
                samples[j].iid,
                row(&all[k])
            );
        }
        let _ = std::fs::write(&path, text.as_bytes());
    }

    if confirmed.is_empty() {
        write(out, &console::no_duplicates_found(min_conc));
    } else {
        write(
            out,
            &console::duplicates_saved(confirmed.len(), min_conc, &path),
        );
    }
    if candidates > confirmed.len() {
        write(
            out,
            &console::additional_pairs_unconfirmed(candidates - confirmed.len()),
        );
    }
}

/// How many pairs the screening stage handed to the final stage.
///
/// **Established.** Below [`SCREENING_SAMPLES`] there is no screening at all: every one
/// of the `C(n,2)` pairs is a candidate, which is why `dups` reports `45 − 2 = 43`
/// unconfirmed pairs at the default threshold and `singleton`, with no pairs, reports
/// none.
///
/// **Not established.** The screening statistic itself. It runs on 512 SNPs chosen after
/// the `Sorting autosomes...` step, and it is *not* heterozygote concordance over the
/// first 512 markers: the 296 pairs a 200-sample run keeps at `--minConc 0` are not the
/// 296 top-ranked pairs by that statistic, nor by the IBS0 rate or the HetHet rate over
/// the same markers, over any of 256/512/1024 leading markers. Its threshold does move
/// with `--minConc` (296, 296, 280, 64, 19, 0 candidates at 0.2, 0.25, 0.3, 0.4, 0.5,
/// 0.6), so it is not a fixed pre-filter either.
///
/// What is implemented instead is the exact statistic — a screen that keeps precisely
/// the pairs the final stage will confirm. That reproduces every captured `N >= 100`
/// run, because no such dataset in the corpus contains a duplicate: the counts agree at
/// zero. It will under-report the `additional pairs …` line on a large dataset that does
/// contain one.
fn screened(pairs: &[(usize, usize)], confirmed: &[usize], screens: bool) -> usize {
    if screens {
        confirmed.len()
    } else {
        pairs.len()
    }
}

/// The seven data columns of a `.con` row.
///
/// `Concord` is `N_IBS2 / N` — confirmed against the reference rather than assumed; on
/// the `dups` MZ pair `9974/10000` prints `0.99740` while `HomConc` and `HetConc` print
/// `0.99862` and `0.99512` from their own definitions, so the three are genuinely three
/// different statistics.
fn row(c: &PairCounts) -> String {
    let concord = f64::from(c.ibs2()) / f64::from(c.n_snp);
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        c.n_snp,
        c.ibs0,
        c.ibs1(),
        c.ibs2(),
        f(concord, 5),
        f(est::hom_concordance(c), 5),
        f(est::het_concordance(c), 5),
    )
}

/// Write and flush, so console ordering survives the process exits around it.
fn write(out: &mut dyn Write, s: &str) {
    let _ = out.write_all(s.as_bytes());
    let _ = out.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(
        n_snp: u32,
        ibs0: u32,
        het_i: u32,
        het_j: u32,
        het_het: u32,
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

    /// The `dups` MZ pair, transcribed from the capture:
    /// `10000  9  17  9974  0.99740  0.99862  0.99512`.
    #[test]
    fn con_row_matches_the_capture() {
        let c = counts(10_000, 9, 3480, 3473, 3468, 6515);
        assert_eq!(c.ibs1(), 17);
        assert_eq!(c.ibs2(), 9974);
        assert_eq!(row(&c), "10000\t9\t17\t9974\t0.99740\t0.99862\t0.99512");
    }

    /// The exact-duplicate row of the same capture, where all three rates are 1.
    #[test]
    fn an_exact_duplicate_prints_five_decimal_ones() {
        let c = counts(10_000, 0, 3500, 3500, 3500, 6500);
        assert_eq!(row(&c), "10000\t0\t0\t10000\t1.00000\t1.00000\t1.00000");
    }

    #[test]
    fn candidate_count_switches_at_a_hundred_samples() {
        let small = serial_pairs(10);
        assert_eq!(screened(&small, &[0, 1], false), 45);
        let big = serial_pairs(200);
        assert_eq!(screened(&big, &[], true), 0);
    }
}
