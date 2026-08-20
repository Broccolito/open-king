//! `--ibs` — the full IBS / concordance table for every pair.
//!
//! Two files and one console body:
//!
//! * `<prefix>.ibs` — within-family pairs, sorted order, with the pedigree's `Z0`/`Phi`
//!   ahead of the counts. **Always** created, header-only when no family has two
//!   members, and never truncated (unlike `.kin`).
//! * `<prefix>.ibs0` — between-family pairs, block-tiled order, no `Z0`/`Phi`, and the
//!   between-family kinship estimator. Created iff the `.fam` names two or more
//!   families.
//! * `<prefix>allsegs.txt` — the usable-segment table, from the segment pre-pass that
//!   `--ibs` runs before anything else.
//!
//! ```text
//! Total length of 14 chromosomal segments usable for IBD segment analysis is 398.2 Mb.
//!   Information of these chromosomal segments can be found in file kingallsegs.txt
//!
//! Within-family IBS data saved in file king.ibs
//!
//! Relationship summary (total relatives: 1 by pedigree, 2 by inference)
//!   ...
//!
//! IBS and relationship inference across families starts at <ctime>
//! 8 CPU cores are used.
//!                                          ends at <ctime>
//! Between-family IBS data saved in file king.ibs0
//! ```
//!
//! The two trailing columns `MaxIBD2`/`Pr_IBD2` are present iff the usable segments
//! total at least 100 Mb — the same condition that suppresses the `Segments too short.`
//! line, so the header and the console can never disagree. Their *values* come from
//! [`segments::Segments::ibd2`], which is fitted rather than established; see its doc
//! comment for exactly how far it is trusted.

use std::fmt::Write as _;
use std::io::Write;

use open_king_core::infer::{pedigree_kinship, pedigree_z0, KinshipCache, Pedigree};
use open_king_core::{counts, kinship as est, PairCounts, Scope};
use open_king_io::Sample;

use crate::analysis::segments::{self, Segments};
use crate::analysis::{
    band, between_family_pairs, f, family_blocks, out_path, within_family_pairs, Class,
};
use crate::cli::{Opt, Options};
use crate::console::{self, RelationshipCounts};
use crate::load::Loaded;

/// Tile width of the `.ibs0` row order — half `.kin0`'s.
const IBS0_BLOCK: usize = 8;

/// The columns shared by both files, from `N_SNP` onwards.
const SHARED_COLUMNS: &str = concat!(
    "N_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\t",
    "IBS\tDist\tHetConc\tHet2|1\tHet1|2\tHomConc\tKinship"
);
/// The two columns the segment pre-pass adds when it has at least 100 Mb to work with.
const SEGMENT_COLUMNS: &str = "\tMaxIBD2\tPr_IBD2";

/// Run the pass: the segment pre-pass, then the two writers.
///
/// The caller has already emitted the preamble and the `Options in effect:` block, and
/// closes the run with `KING ends at`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let sexchr = i64::from(opts.int(Opt::Sexchr));
    let map_warning = crate::analysis::ibdseg::map_order_warning(&loaded.fileset.variants, sexchr);
    let (segs, x_segs) = if map_warning.is_some() {
        (segments::Segments::empty(), segments::Segments::empty())
    } else {
        (
            segments::usable(&loaded.fileset.variants, &loaded.fileset.kept, sexchr),
            segments::usable_x(&loaded.fileset.variants, sexchr),
        )
    };

    let allsegs = out_path(opts, "allsegs.txt");
    if let Some(warning) = map_warning.as_ref() {
        write(out, warning);
    } else if segs.list.is_empty() {
        write(out, console::NO_INFORMATIVE_SEGMENTS);
    } else {
        write(
            out,
            &console::segments_usable(segs.list.len(), segs.total_length()),
        );
        if !x_segs.list.is_empty() {
            write(
                out,
                &console::x_segments_usable(x_segs.list.len(), x_segs.total_length()),
            );
        }
        write(out, &console::segments_file(&allsegs));
        let mut text = segs.allsegs();
        text.push_str(&x_segs.allsegs_continued(segs.list.len()));
        let _ = std::fs::write(&allsegs, text.as_bytes());
    }
    if map_warning.is_none() && !segs.informative() {
        write(out, console::SEGMENTS_TOO_SHORT);
    }

    within_family(opts, loaded, &segs, out);

    if family_blocks(samples).len() < 2 {
        write(out, console::ONLY_ONE_FAMILY);
        return;
    }
    between_family(opts, loaded, &segs, out);
}

// ---------------------------------------------------------------------------
// Within family
// ---------------------------------------------------------------------------

/// One emitted `.ibs` row, kept whole so the relationship summary is computed from the
/// same numbers the file carries.
struct Row {
    /// `.fam` indices, kept so the IBD2 pass can re-read the pair's bit planes.
    i: usize,
    j: usize,
    fid: String,
    id1: String,
    id2: String,
    counts: PairCounts,
    z0: f64,
    phi: f64,
    kinship: f64,
    pedigree: Class,
}

impl Row {
    fn ibs0_prop(&self) -> f64 {
        est::ibs0_prop(&self.counts)
    }
}

fn within_family(opts: &Options, loaded: &Loaded, segs: &Segments, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;
    let pairs = within_family_pairs(samples);

    let pedigree = Pedigree::from_samples(&with_phantom_parents(samples));
    let mut cache = KinshipCache::default();
    let all = counts::all_pairs(genotypes, &pairs);

    let mut rows: Vec<Row> = Vec::with_capacity(pairs.len());
    for (&(i, j), &c) in pairs.iter().zip(all.iter()) {
        // A pair with no jointly-called variant is omitted from the file entirely, and
        // so is one whose estimate is not a number: the reference prints `-inf` but
        // never `nan`.
        let kinship = est::kinship(&c, Scope::WithinFamily);
        if c.n_snp == 0 || kinship.is_nan() {
            continue;
        }
        let phi = pedigree_kinship(&pedigree, &mut cache, i, j);
        let z0 = pedigree_z0(&pedigree, &mut cache, i, j);
        rows.push(Row {
            i,
            j,
            fid: samples[i].fid.clone(),
            id1: samples[i].iid.clone(),
            id2: samples[j].iid.clone(),
            counts: c,
            z0,
            phi,
            kinship,
            pedigree: pedigree_class(phi, z0),
        });
    }

    let path = out_path(opts, ".ibs");
    let mut text = format!("FID\tID1\tID2\tZ0\tPhi\t{SHARED_COLUMNS}");
    if segs.informative() {
        text.push_str(SEGMENT_COLUMNS);
    }
    text.push('\n');
    for r in &rows {
        let _ = write!(
            text,
            "{}\t{}\t{}\t{}\t{}\t",
            r.fid,
            r.id1,
            r.id2,
            f(r.z0, 3),
            f(r.phi, 4)
        );
        text.push_str(&shared_columns(&r.counts, r.kinship));
        if segs.informative() {
            // Un-analysed pairs are `0.000`/`0.0000` here, where `.ibs0` writes `-9`.
            let v = if r.kinship >= segments::IBD2_KINSHIP_GATE {
                segs.ibd2(genotypes, r.i, r.j)
            } else {
                segments::Ibd2::default()
            };
            let _ = write!(text, "\t{}\t{}", f(v.max_bp, 3), f(v.proportion, 4));
        }
        text.push('\n');
    }
    let _ = std::fs::write(&path, text.as_bytes());

    if rows.is_empty() {
        write(out, console::EACH_FAMILY_ONE_INDIVIDUAL);
        return;
    }
    write(out, &console::within_family_ibs_saved(&path));
    if let Some(table) = summary(&rows) {
        write(out, &table);
    }
}

// ---------------------------------------------------------------------------
// Between family
// ---------------------------------------------------------------------------

fn between_family(opts: &Options, loaded: &Loaded, segs: &Segments, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;
    let pairs = between_family_pairs(samples, IBS0_BLOCK);

    write(
        out,
        &console::ibs_across_families_starts(console::now_local()),
    );
    write(out, &console::cpu_cores(crate::analysis::cpu_count(opts)));
    let all = counts::all_pairs(genotypes, &pairs);

    let path = out_path(opts, ".ibs0");
    let mut text = format!("FID1\tID1\tFID2\tID2\t{SHARED_COLUMNS}");
    if segs.informative() {
        text.push_str(SEGMENT_COLUMNS);
    }
    text.push('\n');
    for (&(i, j), &c) in pairs.iter().zip(all.iter()) {
        let kinship = est::kinship(&c, Scope::BetweenFamily);
        if c.n_snp == 0 || kinship.is_nan() {
            continue;
        }
        let _ = write!(
            text,
            "{}\t{}\t{}\t{}\t",
            samples[i].fid, samples[i].iid, samples[j].fid, samples[j].iid
        );
        text.push_str(&shared_columns(&c, kinship));
        if segs.informative() {
            if kinship >= segments::IBD2_KINSHIP_GATE {
                let v = segs.ibd2(genotypes, i, j);
                let _ = write!(text, "\t{}\t{}", f(v.max_bp, 3), f(v.proportion, 4));
            } else {
                // The bare token `-9`, not `-9.000`/`-9.0000`.
                text.push_str("\t-9\t-9");
            }
        }
        text.push('\n');
    }
    let _ = std::fs::write(&path, text.as_bytes());

    write(
        out,
        &console::ends_at(console::IBS_ACROSS_FAMILIES_INDENT, console::now_local()),
    );
    write(out, &console::between_family_ibs_saved(&path));
}

// ---------------------------------------------------------------------------
// Shared row body
// ---------------------------------------------------------------------------

/// The 15 columns from `N_SNP` to `Kinship`, identical in both files.
///
/// Eight exact integer counts, then six derived rates and the estimate. `N_IBS1` and
/// `N_IBS2` are the exact identities `Het_i + Het_j - 2·HetHet` and
/// `N_SNP - IBS0 - N_IBS1`, verified against the reference's own integer columns with
/// zero mismatches.
fn shared_columns(c: &PairCounts, kinship: f64) -> String {
    format!(
        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
        c.n_snp,
        c.ibs0,
        c.ibs1(),
        c.ibs2(),
        c.het_het,
        c.hom_hom,
        c.het_i,
        c.het_j,
        f(est::ibs_mean(c), 4),
        f(est::dist(c), 4),
        f(est::het_concordance(c), 4),
        f(est::het_2given1(c), 4),
        f(est::het_1given2(c), 4),
        f(est::hom_concordance(c), 4),
        f(kinship, 4),
    )
}

// ---------------------------------------------------------------------------
// Pedigree and the relationship summary
// ---------------------------------------------------------------------------

/// Add a founder row for every `PAT`/`MAT` the `.fam` names but does not list.
///
/// The reference resolves parent references by string within the family and
/// materialises the ones it cannot find, so two rows whose `PAT`/`MAT` are `PA`/`MA` are
/// full siblings even when neither parent is genotyped — it prints `Z0 0.250 Phi 0.2500`
/// for them. `open_king_core::infer::Pedigree` treats an unresolvable parent as unknown, so
/// without this every such pair collapses to `Z0 1.000 Phi 0.0000`. Phantoms are
/// appended after the real rows, so every real sample keeps its index.
fn with_phantom_parents(samples: &[Sample]) -> Vec<Sample> {
    let mut known: Vec<(String, String)> = samples
        .iter()
        .map(|s| (s.fid.clone(), s.iid.clone()))
        .collect();
    let mut out = samples.to_vec();
    for s in samples {
        for parent in [&s.pat, &s.mat] {
            if parent == "0" || known.iter().any(|(f, i)| f == &s.fid && i == parent) {
                continue;
            }
            known.push((s.fid.clone(), parent.clone()));
            out.push(Sample {
                fid: s.fid.clone(),
                iid: parent.clone(),
                pat: "0".to_string(),
                mat: "0".to_string(),
                sex: 0,
                pheno: "-9".to_string(),
            });
        }
    }
    out
}

/// Bucket a pedigree-expected `(Phi, Z0)` into a summary column.
fn pedigree_class(phi: f64, z0: f64) -> Class {
    if phi >= band::MZ {
        Class::Mz
    } else if phi >= band::FIRST {
        if z0 == 0.0 {
            Class::Po
        } else {
            Class::Fs
        }
    } else if phi >= band::SECOND {
        Class::Second
    } else if phi >= band::THIRD {
        Class::Third
    } else {
        Class::Other
    }
}

/// The IBS0 proportion below which `--ibs` calls a 1st-degree pair parent–offspring.
///
/// **`--ibs` does not use `--kinship`'s cutoff.** The two summaries are identical on 8 of
/// the 10 corpus datasets and disagree on the other two — `nuclear` is `8 PO / 6 FS`
/// under `--kinship` and `9 PO / 5 FS` under `--ibs`, `monomorphic` is `9/4` against
/// `12/1` — so whatever `--ibs` does here, it is a different rule and must not be shared.
///
/// It is a **constant**, not the data-derived value `--kinship` computes: every
/// data-derived candidate tried (half the mean, half the max, the median of the non-zero
/// IBS0 proportions, each over several pair sets) is refuted by `dups` and `trio`, whose
/// only 1st-degree pairs sit at exactly `IBS0 = 0` and are called `PO` — a rule scaled to
/// the observed spread yields a cutoff of 0 there and calls them `FS`.
///
/// Solving each dataset's summary for the cutoff interval that reproduces it and
/// intersecting gives `(0.0068, 0.0085]`, from `nuclear`'s `(0.0068, 0.0099]` on the low
/// side and `threegen`'s `(0, 0.0085]` on the high side. Any value in that bracket
/// reproduces all 10; `2^-7` is picked because KING's other thresholds sit on the
/// `2^-k` grid, and the bracket — not this particular point in it — is what is
/// established.
const PO_IBS0_CUTOFF: f64 = 0.0078125;

/// Bucket an estimate into a summary column.
fn inferred_class(kinship: f64, ibs0_prop: f64) -> Class {
    if kinship >= band::MZ {
        Class::Mz
    } else if kinship >= band::FIRST {
        if ibs0_prop < PO_IBS0_CUTOFF {
            Class::Po
        } else {
            Class::Fs
        }
    } else if kinship >= band::SECOND {
        Class::Second
    } else if kinship >= band::THIRD {
        Class::Third
    } else {
        Class::Other
    }
}

/// The relationship summary block, or `None` when the reference prints nothing.
///
/// Suppressed when neither row of the table has a single relative in it: `unrelated`'s
/// 45 within-family pairs are all unrelated by both measures and print no table at all,
/// going straight from the `Within-family IBS data saved …` line into the
/// between-family stage with no blank line between them.
fn summary(rows: &[Row]) -> Option<String> {
    let mut pedigree = RelationshipCounts::default();
    let mut inference = RelationshipCounts::default();
    let mut any = false;
    for r in rows {
        let inferred = inferred_class(r.kinship, r.ibs0_prop());
        any |= r.pedigree.is_relative() || inferred.is_relative();
        bump(&mut pedigree, r.pedigree);
        bump(&mut inference, inferred);
    }
    any.then(|| console::relationship_summary(pedigree, inference))
}

fn bump(c: &mut RelationshipCounts, class: Class) {
    let slot = match class {
        Class::Mz => &mut c.mz,
        Class::Po => &mut c.po,
        Class::Fs => &mut c.fs,
        Class::Second => &mut c.second,
        Class::Third => &mut c.third,
        Class::Other => &mut c.other,
    };
    *slot += 1;
}

/// Write and flush, so console ordering survives the process exits around it.
fn write(out: &mut dyn Write, s: &str) {
    let _ = out.write_all(s.as_bytes());
    let _ = out.flush();
}
