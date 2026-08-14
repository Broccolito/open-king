//! The X-chromosome half of `--kinship`: `<prefix>X.kin` and `<prefix>X.kin0`.
//!
//! The reference runs this as a second pass, after the autosomal `.kin0` has been
//! written, and announces it with its own short block:
//!
//! ```text
//!
//! X-chromosome analysis...
//! X-chromosome genotypes stored in <w> 64-bit words for each of <n> individuals.
//! Within-family kinship data saved in file <p>X.kin
//! Relationship inference across families starts at <t>
//!                                          ends at <t>
//! Between-family kinship data saved in file <p>X.kin0
//! ```
//!
//! No CPU-core line, no progress ticks, no relationship summary, no `Note --kinship
//! --degree` hint, and no blank lines inside the block. Both files are written
//! unconditionally once the pass starts — the `sexchr__kinship_cpus1_altfam` capture has
//! a header-only `X.kin` (its `.fam` override leaves no family with two members) and
//! still prints the "Within-family…" line, while the *autosomal* `.kin` is absent from
//! that same run.
//!
//! # When the pass runs
//!
//! Three conditions, all established against the reference:
//!
//! * **at least [`load::X_PASS_MIN_SNPS`] X markers.** Bisected: 511 → nothing, 512 →
//!   both files. On the `sexchr` fixture that is why `--sexchr 24` (300 X SNPs) and
//!   `--sexchr 25` (150) are silent while `--sexchr 2` (2 000) is not.
//! * **no `--degree`.** `--kinship --degree 2` on a fileset with 1 500 X markers writes
//!   `king.kin`/`king.kin0` and nothing else, and prints no X block.
//! * **the autosomal between-family stage ran.** A single-family fileset stops at
//!   `There is only one family.` and never reaches the X pass.
//!
//! # Who takes part
//!
//! Samples whose `.fam` sex is neither 1 nor 2 are **excluded outright** — `sexchr`'s
//! `S_U0A`/`S_U0B` appear in no X row and in no derived statistic. Everyone else is in,
//! and the `<n> individuals` of the console line still counts the whole `.fam`.
//!
//! # The estimators
//!
//! Each pair is classified by the two sexes, and the `Sex` column spells them in the
//! row's own order (`FM` is a female `ID1` and a male `ID2`). Counts are over the
//! pairwise non-missing X set, exactly as on the autosomes.
//!
//! ```text
//! FF   Het = (Het_i + Het_j) / 2N        phi = (HetHet - 2*IBS0) / (Het_i + Het_j)
//! FM   Het = Het_female / N              phi = 0.5  - IBS0 / Het_female
//! MM   Het = H                           phi = 0.75 - IBS0 / (N * H)
//! ```
//!
//! Two things about that table are worth stating outright, because neither is guessable:
//!
//! * **`.kin0` uses the same three forms as `.kin`.** There is no `min(Het_i, Het_j)`
//!   population-structure-robust variant on X — the between-family FF rows of
//!   `sexchr`'s `kingX.kin0` match the within-family formula to four decimals and the
//!   robust one to none of them.
//! * **`H` is imputed for hemizygous males.** A male has no heterozygotes, so `MM` has
//!   no denominator of its own. The reference gives every male the **lower median of the
//!   female heterozygosity rates in his own family** — `sorted[(n-1)/2]`, verified for
//!   family sizes 1 to 7 — computed over each female's *own* non-missing X set, not the
//!   pair's. An `MM` pair then uses **the smaller of the two males' values**, which is
//!   how a cross-family pair between a family whose median is 0.50 and one whose median
//!   is 0.10 comes out at 0.10. A male whose family has no genotyped female has no value
//!   at all; he simply does not take part in the minimum, so `sexchr`'s lone male
//!   singleton `S_UM` pairs off at the big family's 0.331. Only an `MM` pair with no
//!   value on *either* side is dropped, silently and with no gap in anything else.
//!
//! `FM` deliberately does *not* take the minimum: `sexchr`'s `S_M`/`S_SON1` prints the
//! female's 0.339 even though the male's imputed value is 0.331.
//!
//! # `PhiX`
//!
//! The pedigree X-kinship coefficient, from the standard recurrence — males carry one X
//! and inherit it from their mother alone:
//!
//! ```text
//! phi_X(i, i) = 1                    i male
//! phi_X(i, i) = (1 + phi_X(pat_i, mat_i)) / 2      i female
//! phi_X(i, j) = phi_X(i, mat_j)                    j male
//! phi_X(i, j) = (phi_X(i, pat_j) + phi_X(i, mat_j)) / 2   j female
//! ```
//!
//! which gives 0.5 for brothers and for a father–daughter pair, 0.25 for mother–
//! daughter, 0.375 for sisters and **0** for father–son — all four confirmed against
//! `sexchr`'s `PhiX` column. Phantom parents count, so a sibship whose parents are not
//! genotyped still gets the right value.

use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

use king_core::counts;
use king_core::infer::Pedigree;
use king_io::{Genotypes, Sample};

use crate::analysis::{f, family_blocks, out_path, with_phantom_parents};
use crate::cli::{Opt, Options};
use crate::console;
use crate::load::Loaded;

/// Header of `<prefix>X.kin`. No `Error` column, and `Sex`/`PhiX`/`KinshipX` in place of
/// the autosomal `Z0`/`Phi`/`Kinship`.
const XKIN_HEADER: &str = "FID\tID1\tID2\tSex\tN_SNP\tPhiX\tHet\tIBS0\tKinshipX\n";
/// Header of `<prefix>X.kin0`. Same columns minus `PhiX`, plus the second FID.
const XKIN0_HEADER: &str = "FID1\tID1\tFID2\tID2\tSex\tN_SNP\tHet\tIBS0\tKinshipX\n";

/// Sex codes as the `.fam` spells them.
const MALE: u8 = 1;
const FEMALE: u8 = 2;

/// Whether the X pass runs at all.
///
/// `between_family_ran` is the caller's answer to the third condition in the module
/// docs; the other two are read from here.
pub fn runs(opts: &Options, loaded: &Loaded, between_family_ran: bool) -> bool {
    between_family_ran && opts.int(Opt::Degree) == 0 && loaded.x_genotypes.is_some()
}

/// Write both X files and print the block.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let Some(genotypes) = loaded.x_genotypes.as_ref() else {
        return;
    };
    let samples = &loaded.fileset.samples;

    write(out, console::X_CHROMOSOME_ANALYSIS);
    write(
        out,
        &console::x_chromosome_words(genotypes.words_per_sample(), samples.len()),
    );

    let blocks = family_blocks(samples);
    let imputed = imputed_male_het(samples, genotypes, &blocks);

    // --- within family ---------------------------------------------------
    let mut text = String::from(XKIN_HEADER);
    for members in &blocks {
        for (k, &i) in members.iter().enumerate() {
            for &j in &members[k + 1..] {
                if let Some(row) = row(samples, genotypes, &imputed, i, j) {
                    let phi = phi_x(&pedigree_of(samples), i, j);
                    let _ = writeln!(
                        text,
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        samples[i].fid,
                        samples[i].iid,
                        samples[j].iid,
                        row.sex,
                        row.n_snp,
                        f(phi, 4),
                        f(row.het, 3),
                        f(row.ibs0, 4),
                        f(row.kinship, 4),
                    );
                }
            }
        }
    }
    let path = out_path(opts, "X.kin");
    let _ = std::fs::write(Path::new(&path), text.as_bytes());
    write(out, &console::within_family_kinship_saved(&path));

    // --- between families ------------------------------------------------
    write(
        out,
        &console::relationship_inference_starts(console::now_local()),
    );

    // Family-major order — families in sorted order, members sorted inside them — and
    // the plain `i < j` upper triangle over it. This is *not* the autosomal `.kin0`
    // order, which is `.fam` order tiled at 32: a `.fam` listing `MA1, FA, MB1, FB`
    // emits X rows in the order `FA-FB, FA-MB1, MA1-FB, MA1-MB1`.
    let order: Vec<usize> = blocks.iter().flatten().copied().collect();
    let mut text = String::from(XKIN0_HEADER);
    for (k, &i) in order.iter().enumerate() {
        for &j in &order[k + 1..] {
            if samples[i].fid == samples[j].fid {
                continue;
            }
            if let Some(row) = row(samples, genotypes, &imputed, i, j) {
                let _ = writeln!(
                    text,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    samples[i].fid,
                    samples[i].iid,
                    samples[j].fid,
                    samples[j].iid,
                    row.sex,
                    row.n_snp,
                    f(row.het, 3),
                    f(row.ibs0, 4),
                    f(row.kinship, 4),
                );
            }
        }
    }

    write(
        out,
        &console::ends_at(console::RELATIONSHIP_INFERENCE_INDENT, console::now_local()),
    );
    let path = out_path(opts, "X.kin0");
    let _ = std::fs::write(Path::new(&path), text.as_bytes());
    write(out, &console::between_family_kinship_saved(&path));
}

// ---------------------------------------------------------------------------
// Rows
// ---------------------------------------------------------------------------

/// The four numeric columns both X files share, plus the `Sex` spelling.
struct Row {
    sex: &'static str,
    n_snp: u32,
    het: f64,
    ibs0: f64,
    kinship: f64,
}

/// Build one row, or `None` when the reference emits none for this pair.
fn row(
    samples: &[Sample],
    genotypes: &Genotypes,
    imputed: &[Option<f64>],
    i: usize,
    j: usize,
) -> Option<Row> {
    let (si, sj) = (samples[i].sex, samples[j].sex);
    if !matches!(si, MALE | FEMALE) || !matches!(sj, MALE | FEMALE) {
        return None;
    }
    let c = counts::pair_counts(genotypes, i, j);
    if c.n_snp == 0 {
        return None;
    }
    let n = f64::from(c.n_snp);
    let ibs0 = f64::from(c.ibs0);

    let (sex, het, kinship) = match (si, sj) {
        (FEMALE, FEMALE) => {
            let denom = f64::from(c.het_i + c.het_j);
            if denom == 0.0 {
                return None;
            }
            (
                "FF",
                denom / (2.0 * n),
                (f64::from(c.het_het) - 2.0 * ibs0) / denom,
            )
        }
        (FEMALE, MALE) | (MALE, FEMALE) => {
            let female = f64::from(if si == FEMALE { c.het_i } else { c.het_j });
            if female == 0.0 {
                return None;
            }
            let sex = if si == FEMALE { "FM" } else { "MF" };
            (sex, female / n, 0.5 - ibs0 / female)
        }
        _ => {
            // Both hemizygous: the denominator is imputed. The **smaller** of the two
            // imputed values wins, but only defined values take part — a male whose
            // family has no genotyped female contributes nothing rather than a zero, so
            // `sexchr`'s singleton male `S_UM` still pairs off against the big family at
            // *its* median. Only a pair with no defined value on either side is dropped.
            let h = match (imputed[i], imputed[j]) {
                (Some(a), Some(b)) => a.min(b),
                (Some(a), None) | (None, Some(a)) => a,
                (None, None) => return None,
            };
            if h == 0.0 {
                return None;
            }
            ("MM", h, 0.75 - ibs0 / (n * h))
        }
    };
    Some(Row {
        sex,
        n_snp: c.n_snp,
        het,
        ibs0: ibs0 / n,
        kinship,
    })
}

/// Every male's imputed heterozygosity: the lower median of the female heterozygosity
/// rates in his family, or `None` when his family has no genotyped female.
///
/// Females get `None` too — nothing reads it for them.
fn imputed_male_het(
    samples: &[Sample],
    genotypes: &Genotypes,
    blocks: &[Vec<usize>],
) -> Vec<Option<f64>> {
    let mut out = vec![None; samples.len()];
    for members in blocks {
        let mut rates: Vec<f64> = members
            .iter()
            .filter(|&&m| samples[m].sex == FEMALE)
            .filter_map(|&m| {
                let (called, het) = counts::sample_counts(genotypes, m);
                (called > 0).then(|| f64::from(het) / f64::from(called))
            })
            .collect();
        if rates.is_empty() {
            continue;
        }
        rates.sort_by(|a, b| a.total_cmp(b));
        let median = rates[(rates.len() - 1) / 2];
        for &m in members {
            if samples[m].sex == MALE {
                out[m] = Some(median);
            }
        }
    }
    out
}

// ---------------------------------------------------------------------------
// PhiX
// ---------------------------------------------------------------------------

/// Pedigree plus the sexes it needs, phantom founders included.
struct XPedigree {
    ped: Pedigree,
    sex: Vec<u8>,
}

fn pedigree_of(samples: &[Sample]) -> XPedigree {
    let augmented = with_phantom_parents(samples);
    XPedigree {
        sex: augmented.iter().map(|s| s.sex).collect(),
        ped: Pedigree::from_samples(&augmented),
    }
}

/// `phi_X(a, b)` by the recurrence in the module docs.
///
/// Expansion always picks the individual whose ancestry is deeper, which terminates
/// because `Pedigree` has already broken any cycles and a parent is strictly shallower.
fn phi_x(p: &XPedigree, a: usize, b: usize) -> f64 {
    fn go(p: &XPedigree, a: usize, b: usize, fuel: u32) -> f64 {
        if fuel == 0 {
            return 0.0;
        }
        if a == b {
            return if p.sex[a] == MALE {
                1.0
            } else {
                let (pat, mat) = p.ped.parents(a);
                match (pat, mat) {
                    (Some(f), Some(m)) => (1.0 + go(p, f, m, fuel - 1)) / 2.0,
                    _ => 0.5,
                }
            };
        }
        // Expand `b`, chosen as the deeper of the two so founders are reached.
        let (a, b) = if p.ped.depth_of(a) > p.ped.depth_of(b) {
            (b, a)
        } else {
            (a, b)
        };
        let (pat, mat) = p.ped.parents(b);
        if p.sex[b] == MALE {
            // One X, from the mother; an unknown mother is an unrelated founder.
            mat.map_or(0.0, |m| go(p, a, m, fuel - 1))
        } else {
            let from_pat = pat.map_or(0.0, |f| go(p, a, f, fuel - 1));
            let from_mat = mat.map_or(0.0, |m| go(p, a, m, fuel - 1));
            (from_pat + from_mat) / 2.0
        }
    }
    go(p, a, b, 64)
}

fn write(out: &mut dyn Write, s: &str) {
    let _ = out.write_all(s.as_bytes());
    let _ = out.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(fid: &str, iid: &str, pat: &str, mat: &str, sex: u8) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: pat.to_string(),
            mat: mat.to_string(),
            sex,
            pheno: "-9".to_string(),
        }
    }

    /// The `sexchr` nuclear family, whose whole `PhiX` column the reference printed.
    fn nuclear() -> Vec<Sample> {
        vec![
            sample("SEX", "S_F", "0", "0", 1),
            sample("SEX", "S_M", "0", "0", 2),
            sample("SEX", "S_SON1", "S_F", "S_M", 1),
            sample("SEX", "S_SON2", "S_F", "S_M", 1),
            sample("SEX", "S_DAU1", "S_F", "S_M", 2),
            sample("SEX", "S_DAU2", "S_F", "S_M", 2),
        ]
    }

    #[test]
    fn phi_x_reproduces_the_reference_column() {
        let s = nuclear();
        let p = pedigree_of(&s);
        // Indices: 0 father, 1 mother, 2/3 sons, 4/5 daughters.
        assert_eq!(phi_x(&p, 0, 1), 0.0, "father-mother");
        assert_eq!(phi_x(&p, 0, 2), 0.0, "father-son");
        assert_eq!(phi_x(&p, 1, 2), 0.5, "mother-son");
        assert_eq!(phi_x(&p, 2, 3), 0.5, "brothers");
        assert_eq!(phi_x(&p, 0, 4), 0.5, "father-daughter");
        assert_eq!(phi_x(&p, 1, 4), 0.25, "mother-daughter");
        assert_eq!(phi_x(&p, 4, 5), 0.375, "sisters");
        assert_eq!(phi_x(&p, 2, 4), 0.25, "brother-sister");
    }

    #[test]
    fn phantom_parents_make_a_sibship() {
        // Neither parent is genotyped, and the reference still calls these brothers.
        let s = vec![
            sample("F", "B1", "PA", "MA", 1),
            sample("F", "B2", "PA", "MA", 1),
        ];
        let p = pedigree_of(&s);
        assert_eq!(phi_x(&p, 0, 1), 0.5);
    }

    #[test]
    fn unrelated_founders_are_zero() {
        let s = vec![sample("A", "X", "0", "0", 1), sample("B", "Y", "0", "0", 2)];
        let p = pedigree_of(&s);
        assert_eq!(phi_x(&p, 0, 1), 0.0);
    }

    #[test]
    fn the_imputed_median_is_the_lower_one() {
        // Rates 0.1/0.2/0.3/0.4 -> index (4-1)/2 = 1 -> 0.2, as the reference prints.
        let mut rates: Vec<f64> = vec![0.4, 0.1, 0.3, 0.2];
        rates.sort_by(|a, b| a.total_cmp(b));
        assert_eq!(rates[(rates.len() - 1) / 2], 0.2);
        let mut two: Vec<f64> = vec![0.4, 0.1];
        two.sort_by(|a, b| a.total_cmp(b));
        assert_eq!(two[(two.len() - 1) / 2], 0.1);
    }
}
