//! `--related`: the close-relative pass, and the small-sample downgrade that replaces it.
//!
//! `--related` is **not** a synonym for `--kinship`. On a fileset of ten or more samples
//! it runs the segment pre-pass, writes a sixteen-column `.kin` and, when it finds any
//! candidate, a fourteen-column `.kin0` — both carrying the IBD-segment columns
//! `IBD1Seg`/`IBD2Seg`/`PropIBD`/`InfType` on top of the kinship ones. Below ten samples
//! the reference downgrades the whole pass to `--kinship`.
//!
//! # The console body
//!
//! ```text
//! Autosome genotypes stored in <w> words for each of <n> individuals.
//!
//! Options in effect:
//! <tab>--related
//! [<tab>--degree <d>]
//!
//! Total length of <n> chromosomal segments usable for IBD segment analysis is <x> Mb.
//! [  In addition to autosomes, ...]
//!   Information of these chromosomal segments can be found in file <p>allsegs.txt
//!
//! [Within-family kinship data saved in file <p>.kin]
//! [<relationship summary>]
//! [Within-family X-chr IBD-sharing inference saved in file <p>X.kin]
//! [There is only one family.]                          -- or the between-family stage
//! ```
//!
//! # The between-family stage
//!
//! Three flows, and which one runs is the subtlest thing in this module. Everything
//! below was bracketed against the reference binary on constructed filesets; the probes
//! are named where each constant is defined.
//!
//! ```text
//! [A subset of informative SNPs will be used to screen close relatives.]  -- degree <= 2
//! [Sorting autosomes...]
//! Relationship inference across families starts at <ctime>
//! <c> CPU cores are used...
//! ```
//! then one of
//! ```text
//! <41sp>ends at <ctime>                                   -- nothing to look at
//! No close relatives are inferred.
//!
//! ```
//! ```text
//!   Stages 1&2 (with <s> SNPs): <d> pairs of relatives are detected (with kinship > <t>)
//! <31sp>Screening ends at <ctime>                         -- degree <= 2, n >= 100
//!   Final Stage (with <m> SNPs): <c> pairs of relatives (up to <k>-degree) are confirmed
//! <31sp>Inference ends at <ctime>
//! ```
//! ```text
//! <31sp>Inference ends at <ctime>                         -- degree >= 3
//!   <c> pairs of relatives (up to <k>-degree) are identified
//! ```
//! and finally, when `<c>` is zero, `No cryptic relatedness (up to the <d>-degree) is
//! found.` and nothing else; otherwise the between-family summary table, the
//! `Between-family relatives (kinship >= <t>) saved in file <p>.kin0` line and — at
//! degree 1 only — the two-line `Note only duplicates …` advertisement.

use std::collections::HashSet;
use std::fmt::Write as _;
use std::io::Write;
use std::path::Path;

use king_core::ibdseg::{self, Usable};
use king_core::infer::{pedigree_kinship, pedigree_z0, KinshipCache, Pedigree};
use king_core::{counts, kinship as est, PairCounts, Scope};
use king_io::{Genotypes, Sample, Variant};

use crate::analysis::{
    band, between_family_pairs, cpu_count, f, family_blocks, g, out_path, with_phantom_parents,
    xkinship, Class,
};
use crate::cli::{Opt, Options};
use crate::console::{self, RelationshipCounts};
use crate::load::{self, Loaded};

/// Fewest samples the full `--related` pass will run on.
///
/// Ten. Established by a ladder of filesets: nine samples print the replacement notice
/// and emit the ten-column `.kin`, ten samples run the real pass and emit the sixteen-
/// column one. The corpus agrees — `dups` and `sexchr` (ten samples each) take the full
/// path while `missing` and `nuclear` (six) do not.
const MIN_SAMPLES: usize = 10;

/// Fewest samples at which the two-stage screening path runs at all.
///
/// **100**, bisected on sample-count ladders cut from `bigish`: 99 samples print
/// `No close relatives are inferred.` and write no `.kin0`, 100 print the
/// `Stages 1&2 …` block. The gate is unconditional, and it is a reference *bug* worth
/// knowing about: a ten-sample fileset holding a duplicate pair, an MZ pair and a
/// parent–offspring pair across families — every one of which `--kinship --degree 1`
/// reports — still comes out of `--related --degree 1` as "No close relatives are
/// inferred." with no `.kin0` at all. Every corpus dataset but `bigish` is below the
/// gate, which is why they all take the quiet path at degree 1 and 2.
const SCREEN_MIN_SAMPLES: usize = 100;

/// Largest SNP subset the screening stage will use.
///
/// `2^15`. `Stages 1&2 (with <s> SNPs)` prints `min(m, 32768)` — 5 000, 10 000, 20 000
/// and 30 000 on `bigish` truncated to those maps, then 32 768 at both 40 000 and
/// 50 000 — while `Final Stage (with <m> SNPs)` always prints the whole map.
const SCREEN_SNPS: usize = 32_768;

/// Buffer size at which the reference's `.kin` writer flushes, and so the granularity of
/// the truncation bug that follows from it (`kinship::flushed_prefix` documents it).
const FLUSH_BYTES: usize = 65_536;

/// A sample needs at least this many called autosomal variants to enter the
/// between-family stage — the same screen `--kinship` applies.
const MIN_CALLS: u32 = 545;

/// Header of `<prefix>.kin`: the ten `--kinship` columns plus `HetConc`, `HomIBS0` and
/// the four segment columns, with `Error` still last.
const KIN_HEADER: &str = concat!(
    "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\t",
    "IBD1Seg\tIBD2Seg\tPropIBD\tInfType\tError\n"
);
/// Header of `<prefix>.kin0`. No `Z0`/`Phi` and no `Error` — those are `.kin` only.
const KIN0_HEADER: &str = concat!(
    "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\t",
    "IBD1Seg\tIBD2Seg\tPropIBD\tInfType\n"
);
/// Header of `<prefix>X.kin` — the X pass reports IBD sharing, not kinship.
const XKIN_HEADER: &str = "FID\tID1\tID2\tSex1\tSex2\tPhiX\tIBD1Seg\tIBD2Seg\tPropIBD\n";
/// Header of `<prefix>X.kin0`.
const XKIN0_HEADER: &str = "FID1\tID1\tFID2\tID2\tSex1\tSex2\tIBD1Seg\tIBD2Seg\tPropIBD\n";

/// Tile width of the `.kin0` row order, shared with `--kinship`.
const KIN0_BLOCK: usize = 32;

/// Indent of the between-family stage's `Screening ends at` / `Inference ends at` lines.
const STAGE_INDENT: usize = 31;

// ---------------------------------------------------------------------------
// Entry points the dispatcher uses
// ---------------------------------------------------------------------------

/// Whether this sample count sends `--related` down the `--kinship` path.
pub fn downgrades_to_kinship(n_samples: usize) -> bool {
    n_samples < MIN_SAMPLES
}

/// `--related is replaced with --kinship for a small sample size.`, with the blank line
/// the reference puts above it.
///
/// Note the wording differs from `--ibdseg`'s own downgrade notice
/// (`--kinship analysis carried out instead for such a small sample size.`) — a
/// `--related --ibdseg` run on a three-sample fileset prints both, in that order, and
/// runs the kinship pass twice.
pub fn small_sample_notice() -> String {
    "\n--related is replaced with --kinship for a small sample size.\n".to_string()
}

/// The options the full pass echoes under `Options in effect:`.
pub fn options_in_effect(opts: &Options) -> Vec<String> {
    let mut v = vec!["--related".to_string()];
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

// ---------------------------------------------------------------------------
// The narrow interface to the IBD-segment engine
// ---------------------------------------------------------------------------

/// One pair's IBD-segment summary — everything the relatedness writers need from the
/// segment engine, and the only thing they take from it.
///
/// `king_core::ibdseg` owns the calling rule; this module owns the columns. Keeping the
/// contract down to these four numbers is what lets the two be developed against each
/// other: improve the caller and every `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType` and
/// `MaxIBD2` column in the program moves with it, with no other coordination.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PairIbd {
    /// `IBD1Seg` — the share of the usable genome called IBD1 and not IBD2.
    pub ibd1_seg: f64,
    /// `IBD2Seg` — the share called IBD2.
    pub ibd2_seg: f64,
    /// `PropIBD = IBD2Seg + IBD1Seg / 2`, in full precision.
    pub prop_ibd: f64,
    /// `MaxIBD2` — the longest single IBD2 segment, in base pairs.
    pub max_ibd2: f64,
}

impl PairIbd {
    /// The `InfType` label, from the unrounded proportions.
    pub fn inf_type(&self) -> &'static str {
        ibdseg::inf_type(self.ibd1_seg, self.ibd2_seg, self.prop_ibd)
    }
}

/// One marker array — the autosomes or the X chromosome — prepared for the engine.
///
/// KING keeps the two in separate matrices, so their 64-marker word grids are
/// independent and each carries its own usable-segment set and its own denominator.
pub struct Engine {
    /// Positions of the markers this array covers, in bit-plane order.
    pos: Vec<i64>,
    /// The usable segments over them.
    segs: Vec<Usable>,
    /// `D` — the sum of the usable segments' lengths, and the denominator of every
    /// proportion this engine reports.
    denom: i64,
    /// `--seglength` in base pairs.
    seglength_bp: i64,
}

impl Engine {
    /// Build the autosomal engine from the retained marker set — the same markers, in
    /// the same order, that the autosomal bit planes carry.
    pub fn autosomes(
        variants: &[Variant],
        kept: &[usize],
        sexchr: i64,
        seglength_bp: i64,
    ) -> Engine {
        let chr: Vec<i64> = kept
            .iter()
            .map(|&k| load::chromosome_code(&variants[k].chrom, sexchr))
            .collect();
        let pos: Vec<i64> = kept.iter().map(|&k| variants[k].bp).collect();
        Engine::new(chr, pos, seglength_bp)
    }

    /// Build the X engine over the X markers alone, in `.bim` order — which is the order
    /// the X bit planes are packed in.
    pub fn x_chromosome(variants: &[Variant], sexchr: i64, seglength_bp: i64) -> Engine {
        let mut chr = Vec::new();
        let mut pos = Vec::new();
        for v in variants {
            if load::chromosome_code(&v.chrom, sexchr) == sexchr {
                chr.push(sexchr);
                pos.push(v.bp);
            }
        }
        Engine::new(chr, pos, seglength_bp)
    }

    fn new(chr: Vec<i64>, pos: Vec<i64>, seglength_bp: i64) -> Engine {
        let segs = ibdseg::usable_segments(&chr, &pos);
        let denom = ibdseg::denominator(&segs, &pos);
        Engine {
            pos,
            segs,
            denom,
            seglength_bp,
        }
    }

    /// Whether the array yielded any usable segment at all.
    pub fn is_empty(&self) -> bool {
        self.segs.is_empty()
    }

    /// `D`, in base pairs.
    pub fn total_bp(&self) -> i64 {
        self.denom
    }

    /// Scan one pair and reduce it to the four reported numbers.
    ///
    /// `ibd1_seg`/`ibd2_seg`/`prop_ibd` come straight from
    /// [`king_core::ibdseg::pair_segments`]; `max_ibd2` is the longest IBD2 call of the
    /// same scan, which that function aggregates away.
    pub fn pair(&self, genotypes: &Genotypes, i: usize, j: usize) -> PairIbd {
        let seg = ibdseg::pair_segments(genotypes, &self.pos, &self.segs, i, j, self.seglength_bp);
        PairIbd {
            ibd1_seg: seg.ibd1_seg(self.denom),
            ibd2_seg: seg.ibd2_seg(self.denom),
            prop_ibd: seg.prop_ibd(self.denom),
            max_ibd2: self.max_ibd2(genotypes, i, j) as f64,
        }
    }

    /// Longest single IBD2 call, in base pairs.
    ///
    /// Re-runs the engine's own scan and takes the maximum over the same calls
    /// `pair_segments` sums, so `MaxIBD2` can never describe a different segment set
    /// from `IBD2Seg`.
    fn max_ibd2(&self, genotypes: &Genotypes, i: usize, j: usize) -> i64 {
        let mut best = 0;
        for &seg in &self.segs {
            if seg.words() == 0 {
                continue;
            }
            for c in ibdseg::Scan::new(genotypes, i, j, seg).ibd2(&self.pos, self.seglength_bp) {
                best = best.max(self.pos[c.hi] - self.pos[c.lo]);
            }
        }
        best
    }
}

// ---------------------------------------------------------------------------
// The pass
// ---------------------------------------------------------------------------

/// Run the pass: write the files, print the body.
///
/// The caller has already emitted the preamble and the `Options in effect:` block, and
/// closes the run with `KING ends at`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let seglength_bp = crate::analysis::ibdseg::seglength_bp(opts);
    let sexchr = i64::from(opts.int(Opt::Sexchr));
    let engine = Engine::autosomes(
        &loaded.fileset.variants,
        &loaded.fileset.kept,
        sexchr,
        seglength_bp,
    );

    // The segment pre-pass owns the `Total length …` block and `allsegs.txt`; it is
    // byte-identical to the one `--ibdseg` and the QC reports emit.
    write(out, &crate::analysis::ibdseg::segment_prepass(opts, loaded));

    let blocks = family_blocks(samples);
    let within_ran = blocks.iter().any(|m| m.len() >= 2);
    if within_ran {
        let rows = within_family_rows(loaded, &engine, &blocks);
        let path = out_path(opts, ".kin");
        write_kin(&path, &rows, blocks.len() == 1);
        write(out, &console::within_family_kinship_saved(&path));
        if let Some(table) = summary(&rows) {
            write(out, &table);
        }
    }

    // The X pass is a self-contained stage between the two autosomal ones. Unlike
    // `--kinship`'s, it is not suppressed by `--degree`: `sexchr` emits `X.kin` at every
    // degree from 0 to 4.
    let x = x_engine(opts, loaded, seglength_bp);
    if let Some((xengine, xgenotypes)) = x.as_ref() {
        let path = out_path(opts, "X.kin");
        write_x_kin(&path, loaded, xengine, xgenotypes, &blocks);
        write(out, &x_within_saved(&path));
    }

    // One family, and a within-family stage that ran: the reference stops here.
    if within_ran && blocks.len() == 1 {
        write(out, console::ONLY_ONE_FAMILY);
        return;
    }
    between_family(opts, loaded, &engine, x.as_ref(), out);
}

/// The X engine and planes, when the map carries enough X markers to run the X pass.
fn x_engine<'a>(
    opts: &Options,
    loaded: &'a Loaded,
    seglength_bp: i64,
) -> Option<(Engine, &'a Genotypes)> {
    let genotypes = loaded.x_genotypes.as_ref()?;
    let sexchr = i64::from(opts.int(Opt::Sexchr));
    let engine = Engine::x_chromosome(&loaded.fileset.variants, sexchr, seglength_bp);
    Some((engine, genotypes))
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
    ibd: PairIbd,
    pedigree: Class,
}

fn within_family_rows(loaded: &Loaded, engine: &Engine, blocks: &[Vec<usize>]) -> Vec<KinRow> {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;
    let pedigree = Pedigree::from_samples(&with_phantom_parents(samples));
    let mut cache = KinshipCache::default();

    let mut rows = Vec::new();
    for members in blocks {
        for (k, &i) in members.iter().enumerate() {
            for &j in &members[k + 1..] {
                let c = counts::pair_counts(genotypes, i, j);
                // The same two skips `--kinship` applies, and the row sets match on all
                // ten corpus datasets that reach this pass.
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
                    ibd: engine.pair(genotypes, i, j),
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
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            r.fid,
            r.id1,
            r.id2,
            r.counts.n_snp,
            f(r.z0, 3),
            f(r.phi, 4),
            f(est::het_het_prop(&r.counts), 4),
            f(est::ibs0_prop(&r.counts), 4),
            f(est::het_concordance(&r.counts), 4),
            f(hom_ibs0(&r.counts), 4),
            f(r.kinship, 4),
            f(r.ibd.ibd1_seg, 4),
            f(r.ibd.ibd2_seg, 4),
            f(r.ibd.prop_ibd, 4),
            r.ibd.inf_type(),
            g(error_flag(r.pedigree_label(), r.ibd)),
        );
    }
    let body = if single_family {
        flushed_prefix(&text)
    } else {
        &text
    };
    let _ = std::fs::write(Path::new(path), body.as_bytes());
}

impl KinRow {
    /// The pedigree's own relationship label, for the `Error` comparison.
    fn pedigree_label(&self) -> &'static str {
        pedigree_label(self.phi, self.z0)
    }
}

/// The prefix of `text` that actually reaches disk when the dataset is a single family.
///
/// `.kin` rows go into a 64 KiB buffer that is never closed on the one-family path, so
/// the final partial buffer is lost; `threegen`'s 66-row `.kin` comes out zero bytes
/// under both `--kinship` and `--related`. See `kinship::flushed_prefix` for the sweep
/// that identified a buffer rather than a row limit as the mechanism.
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

/// `HomIBS0` — `N_IBS0` over the number of variants at which **either** sample is
/// homozygous for A1, counted over the pairwise non-missing set.
///
/// Undocumented and not guessable from the name: it is neither `IBS0 / HomHom` nor
/// `1 - HomConc`. Re-derived from the raw `.bed` on the 727 within-family rows of
/// `dups`, `multifam`, `monomorphic` and `admixed` with zero mismatches.
///
/// A pair with no A1 homozygote on either side has a zero denominator, and the reference
/// prints the `nan` that follows — the same spelling `HomConc` already reaches.
fn hom_ibs0(c: &PairCounts) -> f64 {
    f64::from(c.ibs0) / f64::from(c.hom_a1_union)
}

/// `Error`, the pedigree-versus-segment disagreement flag.
///
/// **Not `--kinship`'s rule.** That one compares the kinship *estimate* against `Phi` on
/// a multiplicative scale; this one compares the pedigree's relationship *label* against
/// the segment-inferred one, and the two disagree on nine corpus rows.
///
/// The comparison is by label first — an exact match is `0` — then by degree: `0.5` when
/// the two degrees differ by exactly one **and both are second degree or more distant**,
/// `1` otherwise. So a pedigree `PO` inferred as `FS` is a full error even though both
/// are first degree, while `2nd` against `3rd` is only a warning.
///
/// The `InfType` this compares against is *not* the one printed: clause A of the `FS`
/// test drops its `IBD2Seg >= 0.08` requirement here. Validated over the 4 550
/// `InfType`-carrying rows of the golden `.kin`/`.kin0` corpus with zero mismatches.
fn error_flag(pedigree: &str, ibd: PairIbd) -> f64 {
    let inferred = inf_type_for_error(ibd);
    if pedigree == inferred {
        return 0.0;
    }
    let (a, b) = (label_degree(pedigree), label_degree(inferred));
    if a.abs_diff(b) == 1 && a >= 2 && b >= 2 {
        0.5
    } else {
        1.0
    }
}

/// [`king_core::ibdseg::inf_type`] with the `IBD2Seg >= 0.08` guard on the first `FS`
/// clause removed — the variant the `Error` column compares against.
fn inf_type_for_error(ibd: PairIbd) -> &'static str {
    let (p1, p2, pi) = (ibd.ibd1_seg, ibd.ibd2_seg, ibd.prop_ibd);
    if p2 > 0.7 {
        "Dup/MZ"
    } else if p1 + p2 > 0.96 || (p1 + p2 > 0.9 && p2 < 0.08) {
        "PO"
    } else if pi > band::MZ || (pi > 0.32 && p2 > 0.15) {
        "FS"
    } else if pi > band::FIRST {
        "2nd"
    } else if pi > band::SECOND {
        "3rd"
    } else if pi > band::THIRD {
        "4th"
    } else {
        "UN"
    }
}

/// The pedigree's relationship label, on the same alphabet `InfType` uses.
fn pedigree_label(phi: f64, z0: f64) -> &'static str {
    if phi >= band::MZ {
        "Dup/MZ"
    } else if phi >= band::FIRST {
        if z0 == 0.0 {
            "PO"
        } else {
            "FS"
        }
    } else if phi >= band::SECOND {
        "2nd"
    } else if phi >= band::THIRD {
        "3rd"
    } else if phi >= band::FOURTH {
        "4th"
    } else {
        "UN"
    }
}

/// Degree of a relationship label; `Dup/MZ` is zero and `UN` is beyond every band.
fn label_degree(label: &str) -> u32 {
    match label {
        "Dup/MZ" => 0,
        "PO" | "FS" => 1,
        "2nd" => 2,
        "3rd" => 3,
        "4th" => 4,
        _ => 9,
    }
}

/// The summary column a `Dup/MZ`…`UN` label falls in. `4th` and `UN` share `OTHER`.
fn label_class(label: &str) -> Class {
    match label {
        "Dup/MZ" => Class::Mz,
        "PO" => Class::Po,
        "FS" => Class::Fs,
        "2nd" => Class::Second,
        "3rd" => Class::Third,
        _ => Class::Other,
    }
}

/// Bucket a pedigree-expected `(Phi, Z0)` into a summary column.
fn pedigree_class(phi: f64, z0: f64) -> Class {
    label_class(pedigree_label(phi, z0))
}

/// The within-family relationship summary, or `None` when the reference prints nothing.
///
/// The `Pedigree` row is `--kinship`'s; the `Inference` row is **not**. It is a tally of
/// the `InfType` column, not of the kinship estimate — which is why `--related`'s
/// summary differs from `--kinship`'s on five of the seven corpus datasets that reach
/// both (`multifam` is `24/11/1` by kinship and `24/12/0` by segment). Verified column
/// for column against the `InfType` tallies of all seven golden `.kin` files.
fn summary(rows: &[KinRow]) -> Option<String> {
    let mut pedigree = RelationshipCounts::default();
    let mut inference = RelationshipCounts::default();
    let mut any = false;
    for r in rows {
        let inferred = label_class(r.ibd.inf_type());
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

fn between_family(
    opts: &Options,
    loaded: &Loaded,
    engine: &Engine,
    x: Option<&(Engine, &Genotypes)>,
    out: &mut dyn Write,
) {
    let samples = &loaded.fileset.samples;
    let genotypes = &loaded.fileset.genotypes;
    let degree = effective_degree(opts);
    let screening = degree <= 2;

    if screening {
        write(out, SCREENING_HEADER);
    }
    write(
        out,
        &console::relationship_inference_starts(console::now_local()),
    );
    write(out, &cpu_cores(cpu_count(opts)));

    let dropped = screened_out(genotypes);
    let pairs: Vec<(usize, usize)> = between_family_pairs(samples, KIN0_BLOCK)
        .into_iter()
        .filter(|(i, j)| !dropped.contains(i) && !dropped.contains(j))
        .collect();
    let all = counts::all_pairs(genotypes, &pairs);
    let kinships: Vec<f64> = all
        .iter()
        .map(|c| est::kinship(c, Scope::BetweenFamily))
        .collect();

    // The screening stage estimates on a marker subset; every other stage uses the whole
    // map. Below `SCREEN_SNPS` markers the two coincide and the prefix is skipped.
    let screen_snps = genotypes.n_variants.min(SCREEN_SNPS);
    let screen_kinships = match screening_planes(genotypes, screen_snps) {
        None => kinships.clone(),
        Some(prefix) => counts::all_pairs(&prefix, &pairs)
            .iter()
            .map(|c| est::kinship(c, Scope::BetweenFamily))
            .collect(),
    };
    let detected = detected_pairs(
        if screening {
            &screen_kinships
        } else {
            &kinships
        },
        samples.len(),
        degree,
        screening,
    );
    if detected == 0 {
        write(
            out,
            &console::ends_at(console::RELATIONSHIP_INFERENCE_INDENT, console::now_local()),
        );
        write(out, NO_CLOSE_RELATIVES);
        return;
    }

    // The rows themselves. Inclusion is a **disjunction**: a pair is reported if its
    // kinship reaches the degree's band *or* its segment sharing does. The kinship half
    // alone loses six corpus rows — `bigish --degree 2` reports a pair at 0.0870, just
    // under `2^-3.5`, on the strength of its `PropIBD`.
    let kin_cut = 2f64.powf(-(f64::from(degree) + 1.5));
    let prop_cut = 2f64.powf(-(f64::from(degree) + 0.5));
    let mut rows = Vec::new();
    for ((&(i, j), c), &kinship) in pairs.iter().zip(&all).zip(&kinships) {
        if c.n_snp == 0 {
            continue;
        }
        let ibd = engine.pair(genotypes, i, j);
        if !(kinship >= kin_cut || ibd.prop_ibd > prop_cut) {
            continue;
        }
        rows.push((i, j, *c, kinship, ibd));
    }

    // `N pairs … are identified` counts the summary table, and the table never
    // increments its own `4th` column — so a run whose only rows are fourth-degree
    // reports zero pairs while still writing them.
    let mut tally = RelationshipCounts::default();
    for (_, _, _, _, ibd) in &rows {
        let class = label_class(ibd.inf_type());
        if class != Class::Other {
            bump(&mut tally, class);
        }
    }
    let confirmed = tally.mz + tally.po + tally.fs + tally.second + tally.third;

    if screening {
        let snps = loaded.fileset.genotypes.n_variants;
        write(
            out,
            &stages_line(screen_snps, detected, screen_cutoff(degree)),
        );
        write(out, &stage_ends("Screening"));
        write(out, &final_stage_line(snps, confirmed, degree));
        write(out, &stage_ends("Inference"));
    } else {
        write(out, &stage_ends("Inference"));
        write(out, &identified_line(confirmed, degree));
    }

    write_kin0(&out_path(opts, ".kin0"), samples, &rows);
    if let Some((xengine, xgenotypes)) = x {
        write_x_kin0(
            &out_path(opts, "X.kin0"),
            samples,
            xengine,
            xgenotypes,
            &rows,
        );
    }

    if confirmed == 0 {
        write(out, &no_cryptic_relatedness(degree));
        return;
    }
    write(out, &between_family_summary(tally, confirmed));
    write(out, &saved_line(kin_cut, &out_path(opts, ".kin0")));
    if degree == 1 {
        write(out, DEGREE_NOTE);
    }
}

/// `--degree`, with the unset value resolved to the 1 the reference actually applies.
fn effective_degree(opts: &Options) -> i32 {
    match opts.int(Opt::Degree) {
        0 => 1,
        d => d,
    }
}

/// How many pairs the stage counts as candidates — and, when that is zero, the whole
/// stage collapses into `No close relatives are inferred.`
///
/// Two rules, one per flow:
///
/// * **Screening (degree ≤ 2).** Below [`SCREEN_MIN_SAMPLES`] the answer is
///   unconditionally zero, however related the samples are. At or above it, the count is
///   the number of pairs whose kinship exceeds [`screen_cutoff`], and it is exactly the
///   number the `Stages 1&2` line prints: on `bigish` truncated to 32 768 markers — the
///   size at which no subsetting can happen — the reference printed 18 at degree 1 and
///   50 at degree 2, against 18 and 50 pairs over `2^-3` and `2^-4` in that fileset's
///   own `.kin0`.
/// * **Exhaustive (degree ≥ 3).** A candidate is a pair whose kinship exceeds
///   `2^-(degree + 2.5)`, i.e. the reporting threshold of one degree further out.
///   Bracketed to `(0.0209, 0.0228]` at degree 3 by sweeping `unrelated` over 17 marker
///   subsets — `2^-5.5 = 0.02210` — and consistent with every corpus case: `admixed`
///   goes on to write a `.kin0` at degree 3 on the strength of a single 0.0254 pair that
///   the 0.04419 reporting threshold then rejects, leaving the file header-only.
///
/// # Known gap
///
/// The screening count is exact only while `m <= 32768`. Above it the reference screens
/// on a subset of that size and this counts on the whole map, which reproduces `bigish`
/// at degree 1 (18) but not at degree 2 (50 against the reference's 36). Which 32 768
/// markers it picks is unresolved: the first, last, evenly spaced and
/// highest-minor-allele-frequency subsets were each built and scored against the
/// reference, and none reproduces both degrees.
fn detected_pairs(kinships: &[f64], n_samples: usize, degree: i32, screening: bool) -> usize {
    if screening && n_samples < SCREEN_MIN_SAMPLES {
        return 0;
    }
    let cut = if screening {
        screen_cutoff(degree)
    } else {
        2f64.powf(-(f64::from(degree) + 2.5))
    };
    kinships.iter().filter(|&&k| k > cut).count()
}

/// The kinship the `Stages 1&2` line screens on: `2^-(degree + 2)`, printed `%.4lf`.
fn screen_cutoff(degree: i32) -> f64 {
    2f64.powf(-(f64::from(degree) + 2.0))
}

/// The bit planes the screening stage estimates on: the first `snps` markers.
///
/// `None` when the whole map is already that short, which is every corpus dataset but
/// `bigish`. `snps` is always a multiple of 64 there, so the truncation lands on a word
/// boundary and no tail masking is needed; the general case masks anyway rather than
/// leave a `Genotypes` that breaks its own clean-tail contract.
///
/// **Which** markers is the open question. The reference calls them "a subset of
/// informative SNPs" and this takes the map's own prefix, which reproduces `bigish` at
/// degree 1 (18 detected) and not at degree 2 (50 against 36) — see [`detected_pairs`].
fn screening_planes(g: &Genotypes, snps: usize) -> Option<Genotypes> {
    if snps >= g.n_variants {
        return None;
    }
    let words = snps.div_ceil(64);
    let tail = snps % 64;
    let cut = |plane: &[Vec<u64>]| -> Vec<Vec<u64>> {
        plane
            .iter()
            .map(|s| {
                let mut w = s[..words].to_vec();
                if tail != 0 {
                    w[words - 1] &= (1u64 << tail) - 1;
                }
                w
            })
            .collect()
    };
    Some(Genotypes {
        plane0: cut(&g.plane0),
        plane1: cut(&g.plane1),
        n_samples: g.n_samples,
        n_variants: snps,
    })
}

/// Sample indices the between-family stage drops, the same screen `--kinship` applies.
fn screened_out(g: &Genotypes) -> HashSet<usize> {
    (0..g.n_samples)
        .filter(|&i| called(g, i) < MIN_CALLS)
        .collect()
}

/// How many non-missing calls a sample has.
fn called(g: &Genotypes, i: usize) -> u32 {
    let w = g.words_per_sample();
    g.plane0[i][..w]
        .iter()
        .zip(&g.plane1[i][..w])
        .map(|(a, b)| (a | b).count_ones())
        .sum()
}

type Kin0Row = (usize, usize, PairCounts, f64, PairIbd);

fn write_kin0(path: &str, samples: &[Sample], rows: &[Kin0Row]) {
    let mut text = String::from(KIN0_HEADER);
    for (i, j, c, kinship, ibd) in rows {
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            samples[*i].fid,
            samples[*i].iid,
            samples[*j].fid,
            samples[*j].iid,
            c.n_snp,
            f(est::het_het_prop(c), 4),
            f(est::ibs0_prop(c), 4),
            f(est::het_concordance(c), 4),
            f(hom_ibs0(c), 4),
            f(*kinship, 4),
            f(ibd.ibd1_seg, 4),
            f(ibd.ibd2_seg, 4),
            f(ibd.prop_ibd, 4),
            ibd.inf_type(),
        );
    }
    let _ = std::fs::write(Path::new(path), text.as_bytes());
}

// ---------------------------------------------------------------------------
// The X-chromosome stage
// ---------------------------------------------------------------------------

/// `<prefix>X.kin` — X IBD sharing for every within-family pair.
///
/// The columns are the pedigree's X kinship and the segment engine's proportions over
/// the X marker array; there is no kinship estimate and no `InfType`. Sexes are the
/// `.fam` codes, printed raw. The engine needs no sex logic of its own: a hemizygous
/// male is stored homozygous, which is what makes a mother–son pair read `IBD1Seg
/// 1.0000` and two brothers sharing a maternal X read a large `IBD2Seg`, both of which
/// the `sexchr` capture shows.
fn write_x_kin(
    path: &str,
    loaded: &Loaded,
    engine: &Engine,
    genotypes: &Genotypes,
    blocks: &[Vec<usize>],
) {
    let samples = &loaded.fileset.samples;
    let pedigree = xkinship::pedigree_of(samples);
    let mut text = String::from(XKIN_HEADER);
    for members in blocks {
        for (k, &i) in members.iter().enumerate() {
            for &j in &members[k + 1..] {
                let ibd = engine.pair(genotypes, i, j);
                let _ = writeln!(
                    text,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    samples[i].fid,
                    samples[i].iid,
                    samples[j].iid,
                    samples[i].sex,
                    samples[j].sex,
                    f(xkinship::phi_x(&pedigree, i, j), 4),
                    f(ibd.ibd1_seg, 4),
                    f(ibd.ibd2_seg, 4),
                    f(ibd.prop_ibd, 4),
                );
            }
        }
    }
    let _ = std::fs::write(Path::new(path), text.as_bytes());
}

/// `<prefix>X.kin0` — the same X columns for the between-family rows that were reported.
///
/// Written alongside `.kin0` and never announced: the `sexchr --degree 3` capture emits
/// a header-only `X.kin0` without a word about it on the console.
fn write_x_kin0(
    path: &str,
    samples: &[Sample],
    engine: &Engine,
    genotypes: &Genotypes,
    rows: &[Kin0Row],
) {
    let mut text = String::from(XKIN0_HEADER);
    for (i, j, _, _, _) in rows {
        let ibd = engine.pair(genotypes, *i, *j);
        let _ = writeln!(
            text,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            samples[*i].fid,
            samples[*i].iid,
            samples[*j].fid,
            samples[*j].iid,
            samples[*i].sex,
            samples[*j].sex,
            f(ibd.ibd1_seg, 4),
            f(ibd.ibd2_seg, 4),
            f(ibd.prop_ibd, 4),
        );
    }
    let _ = std::fs::write(Path::new(path), text.as_bytes());
}

// ---------------------------------------------------------------------------
// Console lines this pass owns
// ---------------------------------------------------------------------------

/// The two lines that open the screening flow, printed whenever the effective degree is
/// 1 or 2 — before the sample-count gate, and whatever that gate then decides.
const SCREENING_HEADER: &str = concat!(
    "A subset of informative SNPs will be used to screen close relatives.\n",
    "Sorting autosomes...\n",
);

/// What the stage prints when it has no candidate at all. The trailing blank line is the
/// reference's own.
const NO_CLOSE_RELATIVES: &str = "No close relatives are inferred.\n\n";

/// The advertisement that closes a degree-1 run.
const DEGREE_NOTE: &str = concat!(
    "\nNote only duplicates and 1st-degree relatives are included in the inference.\n",
    "  Specifying '--degree 2' if a higher degree relationship inference is needed.\n\n",
);

/// `<c> CPU cores are used...` — three dots, where `--kinship`'s line has one.
fn cpu_cores(n: usize) -> String {
    format!("{n} CPU cores are used...\n")
}

fn stages_line(snps: usize, detected: usize, cutoff: f64) -> String {
    format!(
        "  Stages 1&2 (with {snps} SNPs): {detected} pairs of relatives are detected (with kinship > {})\n",
        f(cutoff, 4)
    )
}

fn final_stage_line(snps: usize, confirmed: u64, degree: i32) -> String {
    format!(
        "  Final Stage (with {snps} SNPs): {confirmed} pairs of relatives (up to {}-degree) are confirmed\n",
        ordinal(degree)
    )
}

fn identified_line(confirmed: u64, degree: i32) -> String {
    format!(
        "  {confirmed} pairs of relatives (up to {}-degree) are identified\n",
        ordinal(degree)
    )
}

fn stage_ends(what: &str) -> String {
    format!(
        "{:STAGE_INDENT$}{what} ends at {}\n",
        "",
        console::ctime(console::now_local())
    )
}

fn no_cryptic_relatedness(degree: i32) -> String {
    format!("No cryptic relatedness (up to the {degree}-degree) is found.\n")
}

fn saved_line(cutoff: f64, path: &str) -> String {
    format!(
        "\nBetween-family relatives (kinship >= {}) saved in file {path}\n",
        f(cutoff, 5)
    )
}

fn x_within_saved(path: &str) -> String {
    format!("Within-family X-chr IBD-sharing inference saved in file {path}\n")
}

/// The reference's ordinals, typos included: 1st, 2nd, then `3nd` and `4nd`.
fn ordinal(degree: i32) -> String {
    match degree {
        1 => "1st".to_string(),
        2 => "2nd".to_string(),
        d => format!("{d}nd"),
    }
}

/// The between-family summary table — one `Inference` row, six columns ending in `4th`
/// rather than `OTHER`, and a rule two characters shorter than the within-family one.
///
/// The `4th` column is never incremented, which is what lets `N pairs … identified` come
/// out one short of the `.kin0` row count on `bigish --degree 4` (59 against 60).
fn between_family_summary(tally: RelationshipCounts, confirmed: u64) -> String {
    let mut s = format!(
        "\nRelationship summary (total relatives: 0 by pedigree, {confirmed} by inference)\n"
    );
    s.push_str("        \tMZ\tPO\tFS\t2nd\t3rd\t4th\n");
    s.push_str("  =========================================================\n");
    let _ = writeln!(
        s,
        "  Inference\t{}\t{}\t{}\t{}\t{}\t0",
        tally.mz, tally.po, tally.fs, tally.second, tally.third
    );
    s.push('\n');
    s
}

fn write(out: &mut dyn Write, s: &str) {
    let _ = out.write_all(s.as_bytes());
    let _ = out.flush();
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli;

    fn parse(args: &[&str]) -> Options {
        let owned: Vec<String> = args.iter().map(|s| (*s).to_string()).collect();
        cli::parse(&owned).options
    }

    fn ibd(p1: f64, p2: f64) -> PairIbd {
        PairIbd {
            ibd1_seg: p1,
            ibd2_seg: p2,
            prop_ibd: p2 + p1 / 2.0,
            max_ibd2: 0.0,
        }
    }

    #[test]
    fn the_downgrade_boundary_is_ten_samples() {
        assert!(downgrades_to_kinship(1));
        assert!(downgrades_to_kinship(9));
        assert!(!downgrades_to_kinship(10));
        assert!(!downgrades_to_kinship(200));
    }

    #[test]
    fn the_notice_carries_its_own_leading_blank_line() {
        assert_eq!(
            small_sample_notice(),
            "\n--related is replaced with --kinship for a small sample size.\n"
        );
    }

    #[test]
    fn the_downgrade_discards_degree() {
        let opts = parse(&["--related", "--degree", "2"]);
        assert_eq!(opts.int(Opt::Degree), 2);
        assert_eq!(opts.without_degree().int(Opt::Degree), 0);
    }

    #[test]
    fn the_full_pass_echoes_related_then_degree() {
        assert_eq!(
            options_in_effect(&parse(&["--related", "--degree", "3"])),
            ["--related", "--degree 3"]
        );
        assert_eq!(options_in_effect(&parse(&["--related"])), ["--related"]);
    }

    #[test]
    fn an_absent_degree_behaves_as_degree_one() {
        assert_eq!(effective_degree(&parse(&["--related"])), 1);
        assert_eq!(effective_degree(&parse(&["--related", "--degree", "4"])), 4);
    }

    /// Golden rows, one per label the corpus reaches.
    #[test]
    fn inf_type_matches_the_golden_columns() {
        assert_eq!(ibd(0.0436, 0.9223).inf_type(), "Dup/MZ");
        assert_eq!(ibd(1.0, 0.0).inf_type(), "PO");
        assert_eq!(ibd(0.5328, 0.2569).inf_type(), "FS");
        assert_eq!(ibd(0.5268, 0.0).inf_type(), "2nd");
        assert_eq!(ibd(0.0, 0.0).inf_type(), "UN");
    }

    #[test]
    fn error_compares_labels_then_degrees() {
        // dups: an undeclared MZ pair, pedigree UN against inferred Dup/MZ.
        assert_eq!(error_flag("UN", ibd(0.0436, 0.9223)), 1.0);
        // dups: a declared PO pair the segments agree with.
        assert_eq!(error_flag("PO", ibd(1.0, 0.0)), 0.0);
        // multifam: declared founders, no sharing.
        assert_eq!(error_flag("UN", ibd(0.0, 0.0)), 0.0);
        // Adjacent distant degrees are a warning — a declared 2nd inferred 3rd...
        assert_eq!(inf_type_for_error(ibd(0.30, 0.0)), "3rd");
        assert_eq!(error_flag("2nd", ibd(0.30, 0.0)), 0.5);
        // ...but the PO/FS split inside the first degree is a full error.
        assert_eq!(error_flag("PO", ibd(0.53, 0.26)), 1.0);
    }

    #[test]
    fn the_error_column_drops_the_ibd2_guard_on_the_first_fs_clause() {
        // PropIBD over 2^-1.5 with almost no IBD2: printed as `2nd`, compared as `FS`.
        let row = ibd(0.8962, 0.0);
        assert_eq!(row.inf_type(), "2nd");
        assert_eq!(inf_type_for_error(row), "FS");
    }

    #[test]
    fn the_screening_flow_is_gated_on_a_hundred_samples() {
        // A hundred pairs far above any cutoff, but only 99 samples: no candidates.
        let kinships = vec![0.5; 100];
        assert_eq!(detected_pairs(&kinships, 99, 1, true), 0);
        assert_eq!(detected_pairs(&kinships, 100, 1, true), 100);
        // The exhaustive flow has no such gate.
        assert_eq!(detected_pairs(&kinships, 10, 3, false), 100);
    }

    #[test]
    fn the_exhaustive_candidate_cutoff_is_one_degree_looser_than_the_row_cutoff() {
        // `unrelated --degree 3`: the sweep put the flip between these two.
        assert_eq!(detected_pairs(&[0.0209], 30, 3, false), 0);
        assert_eq!(detected_pairs(&[0.0228], 30, 3, false), 1);
        // `admixed --degree 3` keeps its 0.0254 pair as a candidate and then rejects it
        // from the file, leaving a header-only `.kin0`.
        assert_eq!(detected_pairs(&[0.0254], 40, 3, false), 1);
        assert!(0.0254 < 2f64.powf(-4.5));
    }

    #[test]
    fn the_screening_cutoffs_are_the_printed_ones() {
        assert_eq!(f(screen_cutoff(1), 4), "0.1250");
        assert_eq!(f(screen_cutoff(2), 4), "0.0625");
    }

    #[test]
    fn the_saved_line_prints_five_decimals_of_the_band() {
        for (degree, want) in [
            (1, "0.17678"),
            (2, "0.08839"),
            (3, "0.04419"),
            (4, "0.02210"),
        ] {
            let cut = 2f64.powf(-(f64::from(degree) + 1.5));
            assert!(
                saved_line(cut, "king.kin0").contains(want),
                "degree {degree}"
            );
        }
    }

    #[test]
    fn ordinals_reproduce_the_reference_typos() {
        assert_eq!(ordinal(1), "1st");
        assert_eq!(ordinal(2), "2nd");
        assert_eq!(ordinal(3), "3nd");
        assert_eq!(ordinal(4), "4nd");
    }

    #[test]
    fn the_between_family_table_never_fills_its_fourth_degree_column() {
        let tally = RelationshipCounts {
            fs: 3,
            second: 23,
            ..RelationshipCounts::default()
        };
        let table = between_family_summary(tally, 26);
        assert!(table.contains("total relatives: 0 by pedigree, 26 by inference"));
        assert!(table.contains("  Inference\t0\t0\t3\t23\t0\t0\n"));
        assert!(table.contains("        \tMZ\tPO\tFS\t2nd\t3rd\t4th\n"));
        assert!(table.contains("\n  =========================================================\n"));
    }

    #[test]
    fn single_family_kin_is_truncated_to_flushed_chunks() {
        let text: String = "X\n".repeat(100_000);
        assert_eq!(flushed_prefix(&text).len(), 3 * FLUSH_BYTES);
        assert_eq!(flushed_prefix("FID\tID1\n a\tb\n"), "");
    }
}
