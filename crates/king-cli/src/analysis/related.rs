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
    band, between_family_pairs, cpu_count, f, family_blocks, g, kinship, out_path,
    with_phantom_parents, xkinship, Class,
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

/// Heterozygote concordance above which `--related` is willing to call a pair `Dup/MZ`.
///
/// **`0.8`, and it is not `--minConc`.** See [`PairIbd::inf_type`].
const DUP_HET_CONCORDANCE: f64 = 0.8;

impl PairIbd {
    /// The `InfType` label as `.kin`/`.kin0` print it: the segment table of
    /// [`king_core::ibdseg::inf_type`] with its `Dup/MZ` clause additionally gated on the
    /// pair's **heterozygote concordance**.
    ///
    /// `.seg` and `.kin` disagree, and this is the whole of the disagreement. Run
    /// `--ibdseg` and `--related` over the same fileset and a pair at
    /// `IBD1Seg 0.0182 / IBD2Seg 0.8128 / PropIBD 0.8219` comes out `Dup/MZ` in `.seg`
    /// and `FS` in `.kin`: `IBD2Seg > 0.7` fires in both, but `.kin` also requires
    /// `HetConc > 0.8` and falls through to the `FS` clause without it.
    ///
    /// Established on a 72-pair ladder sweeping the IBD2 fraction through the boundary
    /// (`tests/parity/probes/mkpairs.py`, targets `(0.0, 0.760…0.828)` × four offsets):
    /// the highest non-`Dup/MZ` row is `HetConc 0.7986` and the lowest `Dup/MZ` row is
    /// `0.8004`, while `Kinship` is *non-monotone* across the same boundary (`0.4082`
    /// prints `FS`, `0.4086` prints `Dup/MZ`), which rules the kinship estimate out as
    /// the gate. `--minConc 0.7` and `--minConc 0.9` leave every label unchanged, so the
    /// `0.8` is hard-coded rather than the duplicate pass's option. Whether the
    /// comparison is `>` or `>=` is unobservable — `HetConc` is a ratio of counts that
    /// lands on `0.8` only by coincidence — and `>` is assumed, as everywhere else in
    /// the table.
    ///
    /// Checked over 6 371 rows — every `InfType`-carrying `.kin` and `.kin0` in the
    /// golden corpus plus 2 123 purpose-built probe rows — with **0** mismatches, against
    /// 181 for the ungated `.seg` rule. The corpus alone cannot see it: its only pairs
    /// over `IBD2Seg 0.7` are true duplicates, which clear `0.8` comfortably.
    pub fn inf_type(&self, het_conc: f64) -> &'static str {
        if self.ibd2_seg > 0.7 && het_conc <= DUP_HET_CONCORDANCE {
            // The `Dup/MZ` clause is suppressed; the rest of the table decides. With
            // `IBD2Seg > 0.7` it collapses: the `PO` clause's second disjunct needs
            // `IBD2Seg < 0.08`, and `PropIBD >= IBD2Seg > 0.7 > 2^-1.5` makes the first
            // `FS` clause fire on everything the `PO` one leaves. So the fall-through is
            // `PO` above the `0.96` sum and `FS` below it — every probe row that reaches
            // here is an `FS`; the `PO` branch follows from the clause order rather than
            // from an observation.
            return if self.ibd1_seg + self.ibd2_seg > 0.96 {
                "PO"
            } else {
                "FS"
            };
        }
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
            for c in
                ibdseg::Scan::new(genotypes, i, j, seg).ibd2(&self.pos, self.seglength_bp, true)
            {
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
    // degree from 0 to 4. It *is* suppressed by the same thing that suppresses the
    // autosomal `.kin` — a fileset with no family of two writes neither, and prints
    // neither line. (`--kinship`'s X pass differs here too: it emits a header-only
    // `X.kin` in exactly that case. Checked on a fileset of twelve singleton families,
    // against the same fileset with one two-member family added, which writes both.)
    let x = x_engine(opts, loaded, seglength_bp);
    if let (true, Some((xengine, xgenotypes))) = (within_ran, x.as_ref()) {
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

/// The X engine and planes, when the X map yields a usable segment to measure over.
///
/// **That construction is the whole gate** — not the 512-marker count `--kinship`'s X
/// pass uses. Both directions were checked against the reference on built filesets: 320 X
/// markers over 30 Mb (one usable segment, well under 512) write `X.kin`, and 640 markers
/// packed into 5 Mb (no usable segment, well over 512) do not. The corpus cannot tell the
/// two rules apart — `sexchr` clears both — which is why the count stood in for the
/// segment test here until it was measured.
fn x_engine<'a>(
    opts: &Options,
    loaded: &'a Loaded,
    seglength_bp: i64,
) -> Option<(Engine, &'a Genotypes)> {
    let genotypes = loaded.x_genotypes.as_ref()?;
    let sexchr = i64::from(opts.int(Opt::Sexchr));
    let engine = Engine::x_chromosome(&loaded.fileset.variants, sexchr, seglength_bp);
    (!engine.is_empty()).then_some((engine, genotypes))
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
        let het_conc = est::het_concordance(&r.counts);
        let inferred = r.ibd.inf_type(het_conc);
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
            f(het_conc, 4),
            f(hom_ibs0(&r.counts), 4),
            f(r.kinship, 4),
            f(r.ibd.ibd1_seg, 4),
            f(r.ibd.ibd2_seg, 4),
            f(r.ibd.prop_ibd, 4),
            inferred,
            g(error_flag(r.pedigree_label(), r.phi, inferred, r.ibd)),
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
///
/// # Known gap: exact four-decimal ties
///
/// The *counts* are right everywhere; the *rounding* of a ratio that lands exactly on a
/// four-decimal tie is not, and the reference's tie-breaking is not reproduced by any
/// evaluation order of this ratio. A fileset with hand-placed genotypes
/// (`probes/pederr.py homibs0`) puts nine pairs on nine exact ties:
///
/// | `N_IBS0` / union | exact value | reference | `a/b` |
/// | --- | --- | --- | --- |
/// | 3 / 32 | 0.09375 | 0.0938 | 0.0938 |
/// | 31 / 32 | 0.96875 | 0.9688 | 0.9688 |
/// | 15 / 96 | 0.15625 | **0.1563** | 0.1562 |
/// | 51 / 96 | 0.53125 | 0.5312 | 0.5312 |
/// | 9 / 160 | 0.05625 | **0.0562** | 0.0563 |
/// | 17 / 160 | 0.10625 | 0.1062 | 0.1062 |
/// | 13 / 160 | 0.08125 | **0.0812** | 0.0813 |
/// | 7 / 160 | 0.04375 | 0.0437 | 0.0437 |
/// | 414 / 960 | 0.43125 | **0.4312** | 0.4313 |
///
/// `15/96` and `51/96` are both exactly representable doubles and both exact ties, and
/// the reference rounds the first up and the second down — so whatever it prints is not
/// the correctly-rounded value of `a/b` under any single tie rule, and its inputs must be
/// perturbed by however it accumulates them. Nineteen algebraically equal forms were
/// scored against the nine rows (`a/b`, `a/(hom_i + hom_j - both)`, `(a/N)/(b/N)`,
/// `a*(1/b)`, `1 - (b-a)/b`, single-precision variants, …); the best scores 8 of 9 and
/// none scores 9, so **none is adopted** — a form fitted to eight hand-made ties would be
/// worse than an honest ratio. `a/b` is kept because it is the algebraically correct one.
///
/// Incidence: **0 of the 4 248** golden `.kin` rows and **2 of 1 189** rows of a random
/// pedigree probe. Ties are rare because the union denominator is rarely a small multiple
/// of a power of two.
fn hom_ibs0(c: &PairCounts) -> f64 {
    f64::from(c.ibs0) / f64::from(c.hom_a1_union)
}

/// The three `InfType` labels `Error` grades numerically rather than by name.
///
/// Everything else — `Dup/MZ`, `PO`, `FS`, `UN` — is a class with a signature of its own,
/// and `Error` scores it by exact agreement. See [`error_flag`].
const GRADED_LABELS: [&str; 3] = ["2nd", "3rd", "4th"];

/// `Error`, the pedigree-versus-inference disagreement flag.
///
/// It is `--kinship`'s [`kinship::error_flag`] fed the **segment** kinship
/// `PropIBD / 2` instead of the kinship estimate — but only for the middle degrees. The
/// full rule, which is the one thing on this page that is neither documented nor
/// guessable:
///
/// ```text
/// InfType in {2nd, 3rd, 4th}  ->  kinship::error_flag(PropIBD / 2, Phi)
/// otherwise                   ->  0 if InfType == the pedigree's own label, else 1
/// ```
///
/// So `Dup/MZ`, `PO`, `FS` and `UN` are all-or-nothing — a declared `FS` inferred `PO`
/// scores `1` however close `PropIBD / 2` lands to `0.25`, and a declared `FS` inferred
/// `FS` scores `0` even at `PropIBD 0.8352`, where the ratio alone would say `0.5` — while
/// `2nd`/`3rd`/`4th` are graded on the same multiplicative `sqrt(2)` / `2` bands
/// `--kinship` uses. That is why two rows with the *same* pair of labels can score
/// differently: `bigish` `B15_C2`/`B15_G_F` is a declared 2nd inferred `3rd` at
/// `PropIBD 0.1756` and scores `0.5`, and a probe pair with the same labels at
/// `PropIBD 0.1215` scores `1` — ratio `0.7024` against `0.4860`, either side of `0.5`.
///
/// The `Phi == 0` fall-back in `kinship::error_flag` is what makes a declared-unrelated
/// pair inferred `4th` score `0.5`: `PropIBD / 2` then lands in `(2^-5.5, 2^-4.5]` by
/// construction.
///
/// Fitted and checked over 5 813 rows — all 4 248 `InfType`-carrying rows of the golden
/// `.kin` corpus plus 1 565 rows from six purpose-built pedigree probes crossing declared
/// `PO`/`FS`/2nd/3rd/4th/unrelated with prescribed `(IBD1Seg, IBD2Seg)` right across the
/// simplex — with **0** mismatches. The label-degree rule this replaced misses 11 of them.
fn error_flag(pedigree: &str, phi: f64, inferred: &str, ibd: PairIbd) -> f64 {
    if !GRADED_LABELS.contains(&inferred) {
        return if inferred == pedigree { 0.0 } else { 1.0 };
    }
    kinship::error_flag(ibd.prop_ibd / 2.0, phi)
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
        let inferred = label_class(r.ibd.inf_type(est::het_concordance(&r.counts)));
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
    for (_, _, c, _, ibd) in &rows {
        let class = label_class(ibd.inf_type(est::het_concordance(c)));
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
/// The screening count is exact only while `m <= 32768`; above it the reference screens
/// on something this does not model, and this counts on the map's first 32 768 markers.
/// It reproduces `bigish` at degree 1 (18) but not at degree 2 (50 against the
/// reference's 36) — the one stdout line in the corpus that `--related` still gets wrong.
///
/// Full record: **`docs/research/22-screen.md`**, instrument
/// `docs/research/fixtures/screendeflate.py` (`facts` re-measures everything below;
/// `screencanvas.py` and `screenweight.py` are the two earlier rigs it builds on). The
/// headline is that the reference's screen is **not the kinship over any subset of
/// markers**, so the placeholder below cannot be repaired by choosing a better subset —
/// and neither can anything else of that shape.
///
/// ## The law the screen obeys
///
/// A **dilution bisection** — replace one member of a real `bigish` pair, at a growing
/// random marker set, with genotypes drawn from the fileset's own allele frequencies, so
/// the pair's kinship falls continuously while no new relative appears — locates the
/// acceptance boundary to one marker, i.e. to ~1e-5 in kinship. 36 bisections at each
/// cutoff on `bigish` (m = 50 000, n = 169):
///
/// ```text
/// cutoff 0.0625:  R = 1.02257 ± 0.00065   boundary kinship 0.07216
/// cutoff 0.1250:  R = 1.02079 ± 0.00062   boundary kinship 0.13264
/// ```
///
/// with `k_screen = 0.5 + R*(k - 0.5)`. The two agree to 0.2 %; a multiplicative rule
/// would need 0.866 and 0.943, an additive offset 0.0097 and 0.0076. On a synthetic
/// flat-MAF fileset, where the deflation is four times larger and so the lever arm 25
/// times longer, the same holds (`R` = 1.0798 / 1.0838 while `cut/k*` moves 0.659 →
/// 0.812). **The law is affine about 0.5, and `R` is exactly 1 whenever `m <= 32768`**
/// (0.99999 ± 0.00001 on `bigish`; 1.00000 on three synthetic MAF spectra).
///
/// The deflation is **systematic, not sampling noise**: realisation spread of the
/// boundary is 0.0018 against a deflation of 0.0089, and the per-pair labels are a sharp
/// threshold — every `bigish` pair above 0.0731 accepted, every one below 0.0718
/// rejected, one inversion inside a 0.0009 window. The per-pair labels also re-confirm
/// that the stage is per-pair: they sum to 36 at degree 2 (and to 17 against 18 at degree
/// 1, so "per-pair" is exact only to ±1).
///
/// ## Why no marker subset can be the answer
///
/// For a pair with IBD probabilities `(k0, k1, k2)`, at a marker of frequency `p`,
/// `E[N_l] = 4pq(1 - 2φ)` and `E[het_l] = 2pq` — the `p²q²` terms cancel exactly. Both the
/// numerator `N = het_i + het_j + 4·IBS0 - 2·HetHet` and the denominator
/// `min(het_i, het_j)` are therefore proportional to `Σ pq` over whatever index set they
/// are summed on, so **the KING robust kinship over any marker subset, under any
/// non-negative per-marker weighting, is unbiased for the same φ**. That is a proof, not
/// a failed search, and three measurements agree with it:
///
/// * top-K-by-MAF subsets of `bigish` count 47/45/44/48 pairs over `2^-4` at
///   K = 50 000/32 768/25 000/16 384, and 41/41/41/40 on its first 16 384 markers — flat,
///   where the reference gives 36;
/// * **replicating a map `r` times** leaves every kinship bit-identical (KING's own
///   `.kin0` confirms it) and still moves the count: 41 → 36 → 33 → 29 → 27 at
///   r = 2…6 from `bigish`'s first 16 384 markers, i.e. `R` = 1.000, 1.021, 1.037, 1.055,
///   1.065. Every sub-multiset of a replicated map is a weighting of the base map;
/// * the one loophole — a subset chosen from data that includes the pair — is simulated
///   and closed: top-32 768 by in-sample MAF gives `R` = 0.995 ± 0.002 (no bias), by
///   in-sample heterozygote count 0.916 ± 0.003 (bias of the wrong sign).
///
/// Nine permutations of `bigish`'s marker order print 36/18 every time, which retires
/// prefixes, strides and word decimations; the boundary bisection, forty times finer,
/// moves by 0.0004 (5 % of the deflation), which is the size of a tie-break inside an
/// informativeness ranking and no more.
///
/// ## What the deflation does track
///
/// It needs the markers that overflow the 32 768 budget to be *informative*, and it
/// grows with how much equally-informative material overflows. Appending 17 232 markers
/// at MAF 0.02 to `bigish`'s first 32 768 leaves the printed count at exactly its
/// m = 32 768 value (50/18, where a `bigish`-sized deflation would read 42); appending the
/// real tail gives 36/18. Two-point MAF maps at m = 65 536 put `R`'s minimum (1.008)
/// exactly where the budget need not split a tied group, climbing to 1.060 at K = 40 000
/// and 1.078 at K = 50 000. Across spectra at m = 50 000: flat 1.080, uniform 1.033,
/// `bigish` 1.022, a low-MAF-heavy beta 1.007 — and one beta point sits *below* one
/// (0.9980 ± 0.0003), so "never below 1" is retired too. `R` is not a function of `(m, n)`
/// alone: `bigish` at m = 50 000 reads 1.0216 while its first 25 000 markers replicated
/// twice — same `m`, same `n` — read 1.0280.
///
/// What that leaves is a shape, recorded as a shape and **not** as a rule: when a map
/// holds more equally-informative markers than 32 768, the reference reaches its budget by
/// something lossy applied *uniformly across markers* rather than by keeping some and
/// dropping others. The uniformity is measured directly — on a flat-MAF map a contiguous
/// clone block at markers `[0,·)`, `[20000,·)`, `[32768,·)` or at the tail gives boundary
/// kinships 0.0957/0.0957/0.0930/0.0940, with no preference for the head of the file.
/// `22-screen.md` §5 lists the two candidate mechanisms that survive.
///
/// Landing the affine law with a fitted `R` would reproduce `bigish` and nothing else, so
/// it is deliberately not landed. The consequence is contained: the count reaches stdout
/// and nothing else. `.kin0`'s row set comes from the exhaustive re-estimate below and is
/// byte-correct at every degree, including on the two `bigish` cases whose stdout this
/// line spoils — the 14 pairs the reference's screen drops all sit below the 0.08839
/// reporting threshold, so no reported row depends on it.
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
/// The reference calls them "a subset of informative SNPs"; this takes the map's own
/// prefix, and **that is a placeholder, not a finding**. Clone-canvas runs against the
/// reference ([`detected_pairs`], `docs/research/fixtures/screencanvas.py`) accept a
/// clone window lying entirely past marker 32 768 — where this prefix estimate reads
/// 0.0020 — at the same true kinship as one at the head of the map, so the reference's
/// subset is not a prefix, a stride or a word decimation, and its screening statistic
/// moves with the *other samples'* genotypes, which no fixed marker choice can do. The
/// prefix is kept only because it is the cheapest subset that reproduces the printed
/// count at degree 1 on every map tried; swapping it for the whole map would trade the
/// degree-1 case (18, correct today) for no gain at degree 2 (47 against 36).
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
            ibd.inf_type(est::het_concordance(c)),
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

    /// Golden rows, one per label the corpus reaches. The corpus's only pairs over
    /// `IBD2Seg 0.7` are true duplicates, so its concordances are all far above the gate.
    #[test]
    fn inf_type_matches_the_golden_columns() {
        assert_eq!(ibd(0.0436, 0.9223).inf_type(0.94), "Dup/MZ");
        assert_eq!(ibd(1.0, 0.0).inf_type(0.34), "PO");
        assert_eq!(ibd(0.5328, 0.2569).inf_type(0.45), "FS");
        assert_eq!(ibd(0.5268, 0.0).inf_type(0.30), "2nd");
        assert_eq!(ibd(0.0, 0.0).inf_type(0.26), "UN");
    }

    /// The `HetConc > 0.8` gate on the `Dup/MZ` clause, at the bracket the ladder gives:
    /// probe rows at `0.7986` and `0.8004` land on opposite sides with `IBD2Seg` in the
    /// same place, which is why the gate cannot be the kinship estimate.
    #[test]
    fn the_dup_clause_needs_het_concordance_over_four_fifths() {
        let row = ibd(0.0, 0.7822);
        assert_eq!(row.inf_type(0.8004), "Dup/MZ");
        assert_eq!(row.inf_type(0.7986), "FS");
        // `.seg` has no such gate: the same numbers there are `Dup/MZ` either way.
        assert_eq!(ibdseg::inf_type(0.0, 0.7822, 0.7822), "Dup/MZ");
        // Below the segment clause the concordance is inert.
        assert_eq!(ibd(0.5328, 0.2569).inf_type(0.99), "FS");
    }

    #[test]
    fn error_grades_the_middle_degrees_and_matches_the_rest_by_name() {
        // dups: an undeclared MZ pair, pedigree UN against inferred Dup/MZ.
        assert_eq!(error_flag("UN", 0.0, "Dup/MZ", ibd(0.0436, 0.9223)), 1.0);
        // dups: a declared PO pair the segments agree with.
        assert_eq!(error_flag("PO", 0.25, "PO", ibd(1.0, 0.0)), 0.0);
        // multifam: declared founders, no sharing.
        assert_eq!(error_flag("UN", 0.0, "UN", ibd(0.0, 0.0)), 0.0);
        // The PO/FS split inside the first degree is a full error however close the
        // segment kinship lands to `Phi` — here it is 0.2625 against 0.25.
        assert_eq!(error_flag("PO", 0.25, "FS", ibd(0.53, 0.26)), 1.0);
        // ...and an agreeing name is a clean 0 even when the ratio alone says 0.5:
        // probe `FS02`, PropIBD 0.8352 against Phi 0.25, ratio 1.67.
        assert_eq!(error_flag("FS", 0.25, "FS", ibd(0.1448, 0.7628)), 0.0);
    }

    /// Two rows with the *same* pair of labels and different scores — the observation
    /// that rules a label-degree rule out. Both are a declared 2nd inferred `3rd`.
    #[test]
    fn the_middle_degrees_are_graded_on_the_kinship_ratio() {
        // bigish B15_C2/B15_G_F: PropIBD 0.1756, ratio 0.7024, just under 1/sqrt(2).
        assert_eq!(error_flag("2nd", 0.125, "3rd", ibd(0.3512, 0.0)), 0.5);
        // probe HA07: PropIBD 0.1215, ratio 0.4860, just under 1/2.
        assert_eq!(error_flag("2nd", 0.125, "3rd", ibd(0.2430, 0.0)), 1.0);
        // A declared-unrelated pair inferred `4th` takes the `Phi == 0` fall-back, which
        // puts PropIBD/2 in the fourth-degree band by construction: probe `un4` P000.
        assert_eq!(error_flag("UN", 0.0, "4th", ibd(0.1374, 0.0)), 0.5);
        assert_eq!(error_flag("UN", 0.0, "3rd", ibd(0.2557, 0.0)), 1.0);
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
