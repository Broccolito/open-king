//! `--kinship`: the within-family `.kin` and between-family `.kin0` writers, and the
//! console body that describes them.
//!
//! The estimators come from `docs/VERIFIED_FORMULAS.md` and the file shapes from
//! `docs/SPEC.md` §6.2/6.3; every sequencing rule below was established by running the
//! reference binary, and the probe that pins each one is named where it is used.
//!
//! # The console body
//!
//! ```text
//! Autosome genotypes stored in <w> words for each of <n> individuals.
//!
//! Options in effect:
//! <tab>--kinship
//! [<tab>--degree <d>][<tab>--cpus <n>][<tab>--prefix <p>]
//!
//! [Within-family kinship data saved in file <p>.kin]
//! [<relationship summary>]
//! [There is only one family.]                     -- or the between-family stage:
//! Relationship inference across families starts at <t>
//! <n> CPU cores are used.
//!                                          ends at <t>
//! Between-family kinship data saved in file <p>.kin0
//! [Note --kinship --degree <n> can filter & speed up the kinship computing.]
//! ```
//!
//! Three conditionals, each pinned against the reference:
//!
//! * **The within-family block** appears iff some family has at least two members. The
//!   `singleton` and `pair` captures — one sample per family — have no `.kin`, no
//!   "Within-family…" line and no summary.
//! * **`There is only one family.`** replaces the entire between-family stage, `.kin0`
//!   included, when the `.fam` names one FID *and* the within-family block ran. A
//!   one-sample, one-family `.fam` does **not** print it: `singleton` goes straight on
//!   to write a header-only `.kin0`. Probed directly with a two-member single family,
//!   which does print it.
//! * **The `Note …` line** appears only when `--degree` is absent; the filtered form of
//!   the "Between-family…" line replaces both.
//!
//! # `--degree`
//!
//! Filters `.kin0` only — never `.kin`, never its columns — on the kinship **estimate**,
//! keeping `kinship >= 2^-(d+1.5)` against the exact double rather than the `%.5lf` the
//! console prints (`docs/BEHAVIOR.md` §Q5). `--degree 0` means unset.
//!
//! # The PO/FS split in the summary
//!
//! `docs/SPEC.md` §7.4 concluded that "for Tier 1 there is no PO/FS discrimination at
//! all". That is wrong for the summary table, which splits the 1st-degree band into `PO`
//! and `FS` and does so from the data. The rule, and the sweep that fixes it, is
//! [`po_cutoff`].

use std::collections::HashSet;
use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

use king_core::infer::{pedigree_kinship, pedigree_z0, KinshipCache, Pedigree};
use king_core::{counts, kinship as est, PairCounts, Scope};
use king_io::{Genotypes, Sample};

use crate::analysis::{
    band, between_family_pairs, cpu_count, f, family_blocks, g, out_path, with_phantom_parents,
    xkinship, Class,
};
use crate::cli::{Opt, Options};
use crate::console::{self, RelationshipCounts};
use crate::load::Loaded;

/// Header of `<prefix>.kin`, tab separated like every row under it.
const KIN_HEADER: &str = "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n";
/// Header of `<prefix>.kin0`. No `Error` column — that one is `.kin` only.
const KIN0_HEADER: &str = "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n";

/// Tile width of the `.kin0` row order. `.ibs0` tiles at 8; this writer is the one that
/// tiles at 32, which is why `bigish` is the only fixture that can tell them apart.
const KIN0_BLOCK: usize = 32;

/// Buffer size at which the reference's `.kin` writer flushes, and so the granularity of
/// the truncation bug that follows from it. See [`flushed_prefix`].
const FLUSH_BYTES: usize = 65_536;

/// A sample needs at least this many called autosomal variants to enter the
/// between-family stage.
///
/// The console announces `M<512` and `docs/SPEC.md` §5.2 records the boundary as
/// `M < 513`; neither is the number. A ladder of samples truncated to 505, 508, 510, 511,
/// 512, 513, 514, 515, 516, 520, 530, 531, 535, 540, 543, **544**, **545**, 546, 548,
/// 550, 552, 556 and 558 calls puts the flip at **544 excluded / 545 included**, which is
/// also where the independent sweep in §5.2's "unresolved wrinkle" landed from the other
/// direction (varying the total variant count rather than the per-sample one). Two
/// constructions agreeing on 545 is what makes it a constant rather than an artefact.
const MIN_CALLS: u32 = 545;

// ---------------------------------------------------------------------------
// Entry point
// ---------------------------------------------------------------------------

/// The options this pass echoes under `Options in effect:`.
///
/// The order is the reference's own and does **not** follow the command line: probing
/// `--kinship --cpus 1 --degree 1` emits `--kinship`, `--degree 1`, `--cpus 1`. A
/// `--degree 0` or `--cpus 0` counts as unset and is not echoed at all — the
/// `trio__kinship_cpus0` capture prints a bare `--kinship` — and `--prefix` appears only
/// when it differs from the default.
pub fn options_in_effect(opts: &Options) -> Vec<String> {
    let mut v = vec!["--kinship".to_string()];
    let degree = opts.int(Opt::Degree);
    if degree != 0 {
        v.push(format!("--degree {degree}"));
    }
    let cpus = opts.int(Opt::Cpus);
    if cpus != 0 {
        v.push(format!("--cpus {cpus}"));
    }
    let prefix = opts.string(Opt::Prefix);
    if prefix != "king" {
        v.push(format!("--prefix {prefix}"));
    }
    v
}

/// Run the pass: write the files, print the body.
///
/// The caller has already emitted the preamble and the `Options in effect:` block, and
/// closes the run with `KING ends at`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;

    let blocks = family_blocks(samples);
    let within_ran = blocks.iter().any(|m| m.len() >= 2);

    if within_ran {
        let rows = within_family_rows(samples, genotypes, &blocks);
        let path = out_path(opts, ".kin");
        write_kin(&path, &rows, blocks.len() == 1);
        write(out, &console::within_family_kinship_saved(&path));
        if let Some(table) = summary(&rows) {
            write(out, &table);
        }
    }

    // One family, and a within-family stage that ran: the reference stops here and
    // writes no `.kin0` at all — and never reaches the X pass either.
    if within_ran && blocks.len() == 1 {
        write(out, console::ONLY_ONE_FAMILY);
        return;
    }

    between_family(opts, loaded, out);

    // The X pass is a second, self-contained stage that follows the autosomal one; it
    // decides for itself whether the map and the flags allow it to run.
    if xkinship::runs(opts, loaded, true) {
        xkinship::run(opts, loaded, out);
    }
}

// ---------------------------------------------------------------------------
// Within family
// ---------------------------------------------------------------------------

/// One emitted `.kin` row, kept whole so the summary is computed from exactly the values
/// the file carries.
struct KinRow {
    fid: String,
    id1: String,
    id2: String,
    counts: PairCounts,
    z0: f64,
    phi: f64,
    kinship: f64,
    ibs0_prop: f64,
    pedigree: Class,
}

fn within_family_rows(
    samples: &[Sample],
    genotypes: &Genotypes,
    blocks: &[Vec<usize>],
) -> Vec<KinRow> {
    let pedigree = Pedigree::from_samples(&with_phantom_parents(samples));
    let mut cache = KinshipCache::default();

    let mut rows = Vec::new();
    for members in blocks {
        for (k, &i) in members.iter().enumerate() {
            for &j in &members[k + 1..] {
                let c = counts::pair_counts(genotypes, i, j);
                // Two skips, both verified against the reference: a pair with no shared
                // calls is omitted from every file, and a within-family pair with no
                // heterozygote on either side is omitted from `.kin`, its estimator
                // having a zero denominator.
                if c.n_snp == 0 || c.het_i + c.het_j == 0 {
                    continue;
                }
                let phi = pedigree_kinship(&pedigree, &mut cache, i, j);
                let z0 = pedigree_z0(&pedigree, &mut cache, i, j);
                rows.push(KinRow {
                    fid: samples[i].fid.clone(),
                    id1: samples[i].iid.clone(),
                    id2: samples[j].iid.clone(),
                    z0,
                    phi,
                    kinship: est::kinship(&c, Scope::WithinFamily),
                    ibs0_prop: est::ibs0_prop(&c),
                    pedigree: pedigree_class(phi, z0),
                    counts: c,
                });
            }
        }
    }
    rows
}

/// Render `.kin` and write it, honouring the reference's truncation bug.
fn write_kin(path: &str, rows: &[KinRow], single_family: bool) {
    let mut text = String::from(KIN_HEADER);
    for r in rows {
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.fid,
            r.id1,
            r.id2,
            r.counts.n_snp,
            f(r.z0, 3),
            f(r.phi, 4),
            f(est::het_het_prop(&r.counts), 4),
            f(r.ibs0_prop, 4),
            f(r.kinship, 4),
            g(error_flag(r.kinship, r.phi)),
        );
    }
    let body = if single_family {
        flushed_prefix(&text)
    } else {
        &text
    };
    let _ = std::fs::write(Path::new(path), body.as_bytes());
}

/// The prefix of `text` that actually reaches disk when the dataset is a single family.
///
/// `.kin` rows go into a buffer that is written out and cleared every time it reaches
/// 64 KiB, and on the one-family path the file is never closed, so the final partial
/// buffer is lost. Sizes measured across a 50–200 sample sweep are all `k × ~65 545`
/// bytes and always end on a line boundary; padding the FID moves the row count while
/// pinning the byte count, which is what identifies a buffer rather than a row limit as
/// the mechanism (`docs/BEHAVIOR.md` §Q7).
///
/// Every dataset small enough to buffer under 64 KiB therefore yields a **zero-byte**
/// file — which is all the parity corpus ever shows, and why the rule was first recorded
/// as "one family ⇒ empty `.kin`".
fn flushed_prefix(text: &str) -> &str {
    let mut pending = 0usize;
    let mut flushed = 0usize;
    let mut at = 0usize;
    for line in text.split_inclusive('\n') {
        at += line.len();
        pending += line.len();
        if pending >= FLUSH_BYTES {
            flushed = at;
            pending = 0;
        }
    }
    &text[..flushed]
}

/// `Error`, the graded pedigree-vs-estimate disagreement flag.
///
/// Not a boolean and not an integer: `0` agrees, `0.5` warns, `1` is an error
/// (`docs/SPEC.md` §4.8, fitted over 247 rows). The bands are multiplicative in
/// `Kinship / Phi` — within a factor of `sqrt(2)` is fine, within a factor of 2 is a
/// warning — and they are **not** the adjacent degree classes.
///
/// Boundary inclusivity was left open by the spec and is settled here by construction: a
/// pair built at `Kinship 0.1250` against `Phi 0.2500` (ratio exactly `1/2`) prints
/// `0.5`, and one at `Kinship 0.5000` (ratio exactly `2`) also prints `0.5`, so both
/// ends of the warning band are closed; `0.1249` prints `1`.
///
/// **Shared with `--related`**, which feeds it the segment kinship `PropIBD / 2` for the
/// `2nd`/`3rd`/`4th` rows of its own `Error` column and decides the rest by label — see
/// [`crate::analysis::related::error_flag`]. The two columns still disagree on nine
/// corpus rows, but only because the two estimates do.
pub(crate) fn error_flag(kinship: f64, phi: f64) -> f64 {
    const SQRT2: f64 = std::f64::consts::SQRT_2;
    if phi > 0.0 {
        let r = kinship / phi;
        if (1.0 / SQRT2..=SQRT2).contains(&r) {
            0.0
        } else if (0.5..=2.0).contains(&r) {
            0.5
        } else {
            1.0
        }
    } else if kinship <= band::FOURTH {
        0.0
    } else if kinship <= band::THIRD {
        0.5
    } else {
        1.0
    }
}

// ---------------------------------------------------------------------------
// Relationship summary
// ---------------------------------------------------------------------------

/// Bucket a pedigree-expected `(Phi, Z0)` into a summary column.
///
/// `Phi` lands on the same bands as an estimate; the 1st-degree band is then split by
/// `Z0`, which is exactly `0.000` for parent–offspring and `0.250` for full siblings.
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

/// Bucket an estimate into a summary column. `po_cut` comes from [`po_cutoff`].
pub(crate) fn inferred_class(kinship: f64, ibs0_prop: f64, po_cut: Option<f64>) -> Class {
    if kinship >= band::MZ {
        Class::Mz
    } else if kinship >= band::FIRST {
        let po = match po_cut {
            Some(cut) => ibs0_prop < cut,
            None => ibs0_prop <= 0.0,
        };
        if po {
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

/// The IBS0 cutoff that splits the 1st-degree band into `PO` and `FS`.
///
/// **Half the mean IBS0 proportion of the pedigree-declared full-sibling pairs that the
/// estimate also places in the 1st degree** — the midpoint between the two classes'
/// expectations, a true parent–offspring pair having no opposite-homozygote site at all.
/// When that set is empty there is nothing to calibrate against and only an exactly zero
/// IBS0 counts as `PO`; that case returns `None` rather than a cutoff of zero, because
/// the two behave differently at the boundary (below).
///
/// This is **not** the `Cutoff value for IBS0 between FS and PO is set at %.4f` of the
/// `--build`/`--cluster` path, whose value `docs/BEHAVIOR.md` §Q2 leaves open: that one
/// sits near `0.0055` on ordinary data, while this one is whatever the declared sibships
/// say. On the corpus they differ by more than an order of magnitude, and a constant
/// `0.0055` mis-classifies `monomorphic`.
///
/// # How it was established
///
/// A designed sweep, every pair at kinship `0.25` over 20 000 SNPs, nine "anchor"
/// sibling pairs pinned at `IBS0 = 0.0100` and one test pair swept:
///
/// | test pair `IBS0` | 0.0040 | 0.004 70 | 0.004 75 | 0.0050 |
/// | --- | --- | --- | --- | --- |
/// | reference calls it | `PO` | `PO` | `FS` | `FS` |
///
/// The flip lands between `94/20000` and `95/20000`, and `0.5·(9·0.01 + t)/10` crosses
/// `t` at `t = 0.09/19 = 0.004 7368`: predicted flip 94→95, observed flip 94→95.
///
/// Three further probes fix the shape:
///
/// * **The comparison is strict.** Ten pairs whose only pedigree-FS 1st-degree member
///   sits at `IBS0 = 0` give a cutoff of `0`, and the reference calls that pair `FS` —
///   so `ibs0 < cutoff`, not `<=`.
/// * **The empty case is different.** The same genotypes with the sibships *not*
///   declared — every sample a founder — call the `IBS0 = 0` pair `PO`, while `5e-6` and
///   up stay `FS`. An empty set is therefore not "cutoff 0 with `<`".
/// * **Only 1st-degree declared sibs count.** Nine declared sib pairs whose estimates
///   land in the 2nd-degree band contribute nothing: with them the calibration set is
///   just the test pair and the sweep never produces a `PO`.
///
/// The rule reproduces every summary row in the parity corpus, including two that no
/// fixed threshold can: in `monomorphic` a declared sib pair at `IBS0 = 0.0008` is `PO`
/// while `0.0044` is not, and *that same pair* is `FS` once its siblings are taken out
/// of the family.
fn po_cutoff(rows: &[KinRow]) -> Option<f64> {
    let mut sum = 0.0;
    let mut n = 0usize;
    for r in rows {
        if r.pedigree == Class::Fs && r.kinship >= band::FIRST && r.kinship < band::MZ {
            sum += r.ibs0_prop;
            n += 1;
        }
    }
    (n > 0).then(|| 0.5 * sum / n as f64)
}

/// The relationship summary block, or `None` when the reference prints nothing.
///
/// It is suppressed when neither row has a single relative in it: the `unrelated`
/// capture has 45 within-family pairs, all pedigree-unrelated and all inferred
/// unrelated, and prints no table at all — going straight from the "Within-family…" line
/// to the between-family stage with no blank line between them.
fn summary(rows: &[KinRow]) -> Option<String> {
    let po_cut = po_cutoff(rows);
    let mut pedigree = RelationshipCounts::default();
    let mut inference = RelationshipCounts::default();
    let mut any = false;
    for r in rows {
        let inferred = inferred_class(r.kinship, r.ibs0_prop, po_cut);
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

// ---------------------------------------------------------------------------
// Between families
// ---------------------------------------------------------------------------

fn between_family(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;

    let dropped = screened_out(genotypes);
    // The screen always applies, but it is only *announced* on the unfiltered path: with
    // any `--degree` the reference drops the same samples silently. Probed on a fileset
    // with one 100-call sample — `--degree 9` keeps every pair above its threshold and
    // still omits that sample's, without a word of explanation.
    if !dropped.is_empty() && opts.int(Opt::Degree) == 0 {
        // The reference's own display bug: the count is right, the names are the first
        // `count` rows of the `.fam` in file order, which are generally not the samples
        // that were dropped. Reproduced deliberately.
        let names: Vec<(&str, &str)> = samples[..dropped.len()]
            .iter()
            .map(|s| (s.fid.as_str(), s.iid.as_str()))
            .collect();
        write(
            out,
            &console::samples_excluded_from_kinship(dropped.len(), &names),
        );
    }

    write(
        out,
        &console::relationship_inference_starts(console::now_local()),
    );
    write(out, &console::cpu_cores(cpu_count(opts)));

    let pairs: Vec<(usize, usize)> = between_family_pairs(samples, KIN0_BLOCK)
        .into_iter()
        .filter(|(i, j)| !dropped.contains(i) && !dropped.contains(j))
        .collect();
    let all = counts::all_pairs(genotypes, &pairs);
    let threshold = degree_threshold(opts.int(Opt::Degree));

    let mut text = String::from(KIN0_HEADER);
    let mut written = 0u64;
    for (&(i, j), c) in pairs.iter().zip(&all) {
        if c.n_snp == 0 {
            continue;
        }
        let kinship = est::kinship(c, Scope::BetweenFamily);
        // `>=` against the exact double, and a `NaN` estimate fails it — which is the
        // right way round: a filtered `.kin0` should not carry a pair whose estimate
        // does not exist.
        if threshold.is_some_and(|t| kinship < t || kinship.is_nan()) {
            continue;
        }
        written += 1;
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            samples[i].fid,
            samples[i].iid,
            samples[j].fid,
            samples[j].iid,
            c.n_snp,
            f(est::het_het_prop(c), 4),
            f(est::ibs0_prop(c), 4),
            f(kinship, 4),
        );
    }

    write(
        out,
        &console::ends_at(console::RELATIONSHIP_INFERENCE_INDENT, console::now_local()),
    );

    let path = out_path(opts, ".kin0");
    let _ = std::fs::write(Path::new(&path), text.as_bytes());
    match opts.int(Opt::Degree) {
        0 => {
            write(out, &console::between_family_kinship_saved(&path));
            write(out, console::KINSHIP_DEGREE_NOTE);
        }
        d => write(
            out,
            &console::between_family_kinship_saved_degree(d, written, &path),
        ),
    }
}

/// `--degree d` keeps `kinship >= 2^-(d+1.5)`, compared against the exact double.
///
/// `--degree 0` is "unset" and keeps everything — the reference does not even echo the
/// flag. Negative degrees diverge between `--kinship` and `--related` and are documented
/// as unspecified; the formula is applied to them unchanged.
fn degree_threshold(degree: i32) -> Option<f64> {
    (degree != 0).then(|| 2f64.powf(-(f64::from(degree) + 1.5)))
}

/// Sample indices the between-family stage drops.
///
/// The reference screens out anyone with fewer than [`MIN_CALLS`] called autosomal
/// variants before the cross-family pass. Screened-out samples still appear in `.kin`,
/// and the pairs they would have contributed are simply absent from `.kin0` — verified
/// on a 12-sample fileset where the three truncated samples vanish from the file while
/// the other nine keep every pair.
///
/// No capture in the parity corpus reaches this: the smallest map there carries 4 150
/// autosomal variants and no sample is missing more than 40 % of them.
fn screened_out(g: &Genotypes) -> HashSet<usize> {
    (0..g.n_samples)
        .filter(|&i| called(g, i) < MIN_CALLS)
        .collect()
}

/// How many non-missing calls a sample has. Missing is the only code with both plane
/// bits clear, so the OR of the planes is the callability mask.
fn called(g: &Genotypes, i: usize) -> u32 {
    let w = g.words_per_sample();
    g.plane0[i][..w]
        .iter()
        .zip(&g.plane1[i][..w])
        .map(|(a, b)| (a | b).count_ones())
        .sum()
}

fn write(out: &mut dyn Write, s: &str) {
    let _ = out.write_all(s.as_bytes());
    let _ = out.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(fid: &str, iid: &str, pat: &str, mat: &str) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: pat.to_string(),
            mat: mat.to_string(),
            sex: 1,
            pheno: "-9".to_string(),
        }
    }

    fn row(pedigree: Class, kinship: f64, ibs0_prop: f64) -> KinRow {
        KinRow {
            fid: "F".into(),
            id1: "A".into(),
            id2: "B".into(),
            counts: PairCounts::default(),
            z0: 0.0,
            phi: 0.0,
            kinship,
            ibs0_prop,
            pedigree,
        }
    }

    #[test]
    fn phantom_parents_make_declared_sibs_full_sibs() {
        // Neither PA nor MA is genotyped; the reference still calls these two full sibs
        // and prints Z0 0.250 / Phi 0.2500 for them.
        let s = vec![sample("F", "A", "PA", "MA"), sample("F", "B", "PA", "MA")];
        let ped = Pedigree::from_samples(&with_phantom_parents(&s));
        assert_eq!(ped.kinship(0, 1), 0.25);
        assert_eq!(ped.z0(0, 1), 0.25);
        // Without the phantoms the same `.fam` looks like two unrelated founders.
        assert_eq!(Pedigree::from_samples(&s).kinship(0, 1), 0.0);
    }

    #[test]
    fn phantom_parents_are_scoped_to_their_family() {
        // `P_F` exists as a sample of another family; MONO's children must not resolve
        // to it, but they must still come out siblings of each other.
        let s = vec![
            sample("OTHER", "P_F", "0", "0"),
            sample("MONO", "C1", "P_F", "P_M"),
            sample("MONO", "C2", "P_F", "P_M"),
        ];
        let ped = Pedigree::from_samples(&with_phantom_parents(&s));
        assert_eq!(ped.kinship(1, 2), 0.25);
        assert_eq!(ped.kinship(0, 1), 0.0);
    }

    #[test]
    fn pedigree_classes_follow_phi_then_z0() {
        assert_eq!(pedigree_class(0.25, 0.0), Class::Po);
        assert_eq!(pedigree_class(0.25, 0.25), Class::Fs);
        assert_eq!(pedigree_class(0.125, 0.5), Class::Second);
        assert_eq!(pedigree_class(0.0625, 0.75), Class::Third);
        // 4th degree has no column of its own in the within-family table.
        assert_eq!(pedigree_class(0.03125, 0.875), Class::Other);
        assert_eq!(pedigree_class(0.0, 1.0), Class::Other);
    }

    #[test]
    fn error_flag_bands_are_closed_at_the_factor_of_two() {
        // Reference probe, Phi = 0.25 throughout.
        assert_eq!(error_flag(0.5000, 0.25), 0.5); // ratio exactly 2
        assert_eq!(error_flag(0.1250, 0.25), 0.5); // ratio exactly 1/2
        assert_eq!(error_flag(0.1249, 0.25), 1.0);
        assert_eq!(error_flag(0.2500, 0.25), 0.0);
        assert_eq!(error_flag(0.1768, 0.25), 0.0); // just above 2^-0.5
        assert_eq!(error_flag(0.3535, 0.25), 0.0); // just below 2^0.5
                                                   // Golden rows from the corpus.
        assert_eq!(error_flag(0.1708, 0.25), 0.5); // multifam B_C1/B_C2
        assert_eq!(error_flag(0.4962, 0.0), 1.0); // dups MZ_1/MZ_2
        assert_eq!(error_flag(-0.0085, 0.0), 0.0); // unrelated P01/P02
        assert_eq!(g(0.0), "0");
        assert_eq!(g(0.5), "0.5");
        assert_eq!(g(1.0), "1");
    }

    #[test]
    fn degree_thresholds_are_the_exact_powers() {
        assert_eq!(degree_threshold(0), None);
        assert_eq!(degree_threshold(1), Some(band::FIRST));
        assert_eq!(degree_threshold(2), Some(band::SECOND));
        // A pair at 0.1767775 sits above 2^-2.5 but below the printed 0.17678 and the
        // reference keeps it, which is what pins the comparison to the exact double.
        assert!(0.176_777_5 >= degree_threshold(1).unwrap());
    }

    #[test]
    fn single_family_kin_is_truncated_to_flushed_chunks() {
        // 200 000 bytes of two-byte lines: three whole buffers reach disk and the
        // 3 392-byte remainder dies with the unclosed file.
        let text: String = "X\n".repeat(100_000);
        assert_eq!(flushed_prefix(&text).len(), 3 * FLUSH_BYTES);
        // A line length that does not divide 64 KiB makes each flush overshoot to the
        // next line boundary, so the cut lands on a line rather than on the mark — which
        // is what the reference's byte counts show (65 540, 65 566, 65 611 …).
        let line = "XY\n";
        let text: String = line.repeat(70_000);
        let kept = flushed_prefix(&text);
        assert_eq!(kept.len() % line.len(), 0);
        assert_eq!(
            kept.len(),
            3 * FLUSH_BYTES.div_ceil(line.len()) * line.len()
        );
        assert!(kept.len() > 3 * FLUSH_BYTES);
        // Anything under one buffer never reaches disk at all — the whole parity corpus
        // is in this regime, which is why the rule first looked like "always empty".
        assert_eq!(flushed_prefix("FID\tID1\n a\tb\n"), "");
    }

    #[test]
    fn po_cutoff_is_half_the_declared_sibs_mean() {
        // Nine anchors at 0.0100 plus one test pair: the reference flips at 94/20000.
        let anchors = || (0..9).map(|_| row(Class::Fs, 0.25, 0.01));
        let rows: Vec<KinRow> = anchors()
            .chain([row(Class::Fs, 0.25, 94.0 / 20000.0)])
            .collect();
        assert!(94.0 / 20000.0 < po_cutoff(&rows).unwrap());
        let rows: Vec<KinRow> = anchors()
            .chain([row(Class::Fs, 0.25, 95.0 / 20000.0)])
            .collect();
        assert!(95.0 / 20000.0 >= po_cutoff(&rows).unwrap());
    }

    #[test]
    fn po_cutoff_ignores_declared_sibs_outside_the_first_degree() {
        let rows: Vec<KinRow> = (0..9)
            .map(|_| row(Class::Fs, 0.125, 0.01))
            .chain([row(Class::Fs, 0.25, 0.0)])
            .collect();
        // The 2nd-degree sibs contribute nothing, so the cutoff is zero — and a zero
        // cutoff makes even an IBS0 of zero a full sibling, because the test is strict.
        assert_eq!(po_cutoff(&rows), Some(0.0));
        assert_eq!(inferred_class(0.25, 0.0, Some(0.0)), Class::Fs);
    }

    #[test]
    fn without_declared_sibs_only_a_zero_ibs0_is_parent_offspring() {
        let rows = vec![row(Class::Other, 0.25, 0.0), row(Class::Other, 0.25, 5e-6)];
        assert_eq!(po_cutoff(&rows), None);
        assert_eq!(inferred_class(0.25, 0.0, None), Class::Po);
        assert_eq!(inferred_class(0.25, 5e-6, None), Class::Fs);
    }

    #[test]
    fn summary_is_suppressed_when_nothing_is_related() {
        // The `unrelated` capture: 45 within-family pairs, none related either way.
        let rows: Vec<KinRow> = (0..45).map(|_| row(Class::Other, 0.001, 0.068)).collect();
        assert!(summary(&rows).is_none());
        // One relative anywhere brings the whole table back.
        let rows: Vec<KinRow> = std::iter::once(row(Class::Other, 0.25, 0.0))
            .chain((0..44).map(|_| row(Class::Other, 0.001, 0.068)))
            .collect();
        let table = summary(&rows).unwrap();
        assert!(table.contains("total relatives: 0 by pedigree, 1 by inference"));
    }

    #[test]
    fn options_are_echoed_in_the_reference_order() {
        let echo = |args: &[&str]| {
            let argv: Vec<String> = args.iter().map(|s| s.to_string()).collect();
            options_in_effect(&crate::cli::parse(&argv).options)
        };
        assert_eq!(echo(&["--kinship"]), ["--kinship"]);
        // The block order is the reference's own, not the command line's: this exact
        // invocation emits --kinship, --degree 1, --cpus 1.
        assert_eq!(
            echo(&["--kinship", "--cpus", "1", "--degree", "1"]),
            ["--kinship", "--degree 1", "--cpus 1"]
        );
        assert_eq!(
            echo(&["--kinship", "--prefix", "custom"]),
            ["--kinship", "--prefix custom"]
        );
        // Zero means unset for both integers, and the default prefix is silent — the
        // `trio__kinship_cpus0` capture prints a bare `--kinship`.
        assert_eq!(echo(&["--kinship", "--cpus", "0"]), ["--kinship"]);
        assert_eq!(echo(&["--kinship", "--degree", "0"]), ["--kinship"]);
        assert_eq!(echo(&["--kinship", "--prefix", "king"]), ["--kinship"]);
    }
}
