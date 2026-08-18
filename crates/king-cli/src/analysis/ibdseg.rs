//! `--ibdseg` — pairwise IBD-segment inference.
//!
//! The pass owns three files unconditionally, a fourth ([`xseg`]) when the run has an X
//! chromosome and a `--degree`, and one console body:
//!
//! ```text
//! kingsplitped.txt is generated for certain pedigree plot applications.
//!
//! Options in effect:
//! \t--ibdseg
//! [\t--degree <d>]
//!
//! Total length of <n> chromosomal segments usable for IBD segment analysis is <x> Mb.
//! [  In addition to autosomes, <n> segments of length <x> Mb on X-chr can be further used.]
//!   Information of these chromosomal segments can be found in file <p>allsegs.txt
//!
//! IBD segment analysis starts at <ctime>
//! <c> CPU cores are used for autosome inference...
//!                        ends at <ctime>
//!
//! Note with relationship inference as the primary goal, the following filters are applied:
//!   Sample pairs without any long IBD segments (>10Mb) are excluded.
//!   Short IBD segments (<3Mb) are not reported/utilized.
//! Summary statistics of IBD segments for individual pairs saved in file <p>.seg
//! [Additional summary statistics of X-Chr IBD segments saved in file <p>X.seg]
//! ```
//!
//! Two lines in there are hard-coded in the reference and must be copied verbatim even
//! when they are wrong: the `(<3Mb)` still says `3` under `--seglength 10`, and the
//! `(>10Mb)` pair filter is not tunable at all.
//!
//! `<prefix>splitped.txt` is written from [`crate::analysis::splitped`] before any
//! segment work; it depends only on the `.fam`, which is why it is byte-identical under
//! every `--degree` / `--seglength` / `--related` variant.
//!
//! # `<prefix>.seg` is not the `.kin` table with different columns
//!
//! It differs from every other pairwise file this crate writes in two ways, and both were
//! found the same way — by diffing bytes after the numbers stopped being wrong:
//!
//! * **`PropIBD` is computed from the two columns printed beside it**, not from the
//!   underlying totals. The reference prints two different `PropIBD` values for the same
//!   pair in the same run — one in `.kin`, one in `.seg`. [`ibdseg::seg_prop_ibd`].
//! * **The rows are ordered by 16-sample block**, not by sample index. [`seg_pair_order`].
//!
//! Neither changes a value or a row: the same pairs with the same estimates, printed
//! differently. Together they took this pass from 18 of 52 captured `--ibdseg`
//! invocations byte-identical to 41 of 52, and the default 3 Mb floor to all of them.

use std::io::Write;
use std::path::PathBuf;

use king_core::ibdseg::{self, Usable};
use king_io::Variant;

use crate::analysis::{splitped, xseg};
use crate::cli::{Opt, Options};
use crate::console;
use crate::load::{self, Class, Loaded};

const SORT_NOTE: &str =
    "  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n";
const NO_SEGMENTS: &str = concat!(
    "No informative IBD segments.\n",
    "  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n",
);

/// The first map-order problem KING reports before any segment analysis.
pub fn map_order_warning(variants: &[Variant], sexchr: i64) -> Option<String> {
    let mut previous: Option<(&Variant, i64)> = None;
    for variant in variants {
        let class = load::classify(&variant.chrom, sexchr);
        if !class.is_autosomal() && class != Class::X {
            continue;
        }
        let code = load::chromosome_code(&variant.chrom, sexchr);
        if let Some((before, before_code)) = previous {
            if code < before_code {
                return Some(format!(
                    "Chromosomes unsorted: {} on chr {}, {} on chr {}.\n",
                    before.id, before_code, variant.id, code
                ));
            }
            if code == before_code && variant.bp < before.bp {
                return Some(format!(
                    "Positions unsorted: {} at {}, {} at {}.\n",
                    before.id, before.bp, variant.id, variant.bp
                ));
            }
        }
        previous = Some((variant, code));
    }
    None
}

/// Below this many samples the reference silently runs `--kinship` instead.
///
/// Measured by sweeping `n` over identical filesets: 2, 3 and 4 print
/// `--kinship analysis carried out instead for such a small sample size.` and take the
/// kinship path; 5 and up run the segment engine. This is why the corpus's `pair`,
/// `trio` and `singleton` captures contain no `.seg` at all.
pub const MIN_SAMPLES: usize = 5;

/// Whether `--ibdseg` on this fileset is really going to be a `--kinship` run.
pub fn downgrades_to_kinship(n_samples: usize) -> bool {
    n_samples < MIN_SAMPLES
}

/// `--kinship analysis carried out instead for such a small sample size.`, with the blank
/// line the reference puts above it.
pub fn small_sample_notice() -> String {
    "\n--kinship analysis carried out instead for such a small sample size.\n".to_string()
}

/// The `Options in effect:` entries for this pass.
///
/// `--seglength` is deliberately absent: the reference echoes it earlier, in its own
/// `Minimum segment length is set as <n> bp` line, and never in this block.
pub fn options_in_effect(opts: &Options) -> Vec<String> {
    let mut v = vec!["--ibdseg".to_string()];
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

/// `--seglength` in base pairs.
///
/// The flag's units are Mb and its accepted range is 1..=10; anything outside silently
/// reverts to 3 Mb (`--seglength 0` and `--seglength 11` both produce output that is
/// md5-identical to the default).
pub fn seglength_bp(opts: &Options) -> i64 {
    let mb = opts.double(Opt::Seglength);
    if mb > console::SEGLENGTH_MIN && mb < console::SEGLENGTH_MAX {
        (mb * 1e6).round() as i64
    } else {
        ibdseg::DEFAULT_SEGLENGTH_BP
    }
}

/// The two marker arrays the engine works over, each with its own global word grid.
///
/// KING keeps autosomes and X in separate matrices, so their word grids are independent:
/// an X segment's word alignment does not depend on how many autosomal markers precede
/// it. Both are built here in `.bim` order.
struct Arrays {
    /// Chromosome codes of the retained autosomal markers.
    auto_chr: Vec<i64>,
    /// Positions of the retained autosomal markers.
    auto_pos: Vec<i64>,
    /// Indices into the `.bim` for those markers, for `StartSNP`/`StopSNP`.
    auto_idx: Vec<usize>,
    x_chr: Vec<i64>,
    x_pos: Vec<i64>,
    x_idx: Vec<usize>,
}

fn arrays(variants: &[Variant], sexchr: i64) -> Arrays {
    let mut a = Arrays {
        auto_chr: Vec::new(),
        auto_pos: Vec::new(),
        auto_idx: Vec::new(),
        x_chr: Vec::new(),
        x_pos: Vec::new(),
        x_idx: Vec::new(),
    };
    for (i, v) in variants.iter().enumerate() {
        let class = load::classify(&v.chrom, sexchr);
        let code = load::chromosome_code(&v.chrom, sexchr);
        if class.is_autosomal() {
            a.auto_chr.push(code);
            a.auto_pos.push(v.bp);
            a.auto_idx.push(i);
        } else if class == Class::X {
            a.x_chr.push(code);
            a.x_pos.push(v.bp);
            a.x_idx.push(i);
        }
    }
    a
}

/// Render `<prefix>allsegs.txt`.
///
/// Tab separated; `%.3lf` on the three Mb columns; `Segment` is a 1-based running index
/// across the **whole file**, autosomes then X, not per chromosome.
fn allsegs_text(variants: &[Variant], a: &Arrays, auto: &[Usable], xseg: &[Usable]) -> String {
    let mut s = String::from("Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP\n");
    let mut n = 0usize;
    for (segs, pos, idx) in [(auto, &a.auto_pos, &a.auto_idx), (xseg, &a.x_pos, &a.x_idx)] {
        for seg in segs {
            n += 1;
            let start = pos[seg.lo] as f64 / 1e6;
            let stop = pos[seg.hi] as f64 / 1e6;
            s.push_str(&format!(
                "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}\n",
                n,
                seg.chr,
                start,
                stop,
                stop - start,
                seg.n_snp(),
                variants[idx[seg.lo]].id,
                variants[idx[seg.hi]].id,
            ));
        }
    }
    s
}

/// One reported pair.
struct Row {
    /// `.fam` indices of the pair, kept so [`xseg`] can re-scan exactly these pairs, in
    /// exactly this order, over the X array.
    i: usize,
    j: usize,
    fid1: String,
    id1: String,
    fid2: String,
    id2: String,
    pi1: f64,
    pi2: f64,
    prop: f64,
}

/// Run the pass: write the files, print the body.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let prefix = opts.string(Opt::Prefix).to_string();
    let sexchr = i64::from(opts.int(Opt::Sexchr));
    // Announced before `Options in effect:` and written before any segment work, so it
    // survives even the `No informative IBD segments.` early exit below.
    write_file(
        &format!("{prefix}splitped.txt"),
        &splitped::text(&loaded.fileset.samples),
    );
    if let Some(warning) = map_order_warning(&loaded.fileset.variants, sexchr) {
        let _ = out.write_all(warning.as_bytes());
        let _ = out.write_all(SORT_NOTE.as_bytes());
        return;
    }
    let a = arrays(&loaded.fileset.variants, sexchr);
    let auto = ibdseg::usable_segments(&a.auto_chr, &a.auto_pos);
    let x_usable = ibdseg::usable_segments(&a.x_chr, &a.x_pos);

    if auto.is_empty() {
        // No denominator, so nothing downstream can be computed and no file is written.
        let _ = out.write_all(NO_SEGMENTS.as_bytes());
        return;
    }

    let denom = ibdseg::denominator(&auto, &a.auto_pos);
    let _ = out.write_all(
        format!(
            "Total length of {} chromosomal segments usable for IBD segment analysis is {:.1} Mb.\n",
            auto.len(),
            denom as f64 / 1e6
        )
        .as_bytes(),
    );
    if !x_usable.is_empty() {
        let xlen = ibdseg::denominator(&x_usable, &a.x_pos);
        let _ = out.write_all(
            format!(
                "  In addition to autosomes, {} segments of length {:.1} Mb on X-chr can be further used.\n",
                x_usable.len(),
                xlen as f64 / 1e6
            )
            .as_bytes(),
        );
    }
    let allsegs_path = format!("{prefix}allsegs.txt");
    let _ = out.write_all(
        format!(
            "  Information of these chromosomal segments can be found in file {allsegs_path}\n\n"
        )
        .as_bytes(),
    );
    write_file(
        &allsegs_path,
        &allsegs_text(&loaded.fileset.variants, &a, &auto, &x_usable),
    );

    let _ = out.write_all(
        format!(
            "IBD segment analysis starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(
        format!(
            "{} CPU cores are used for autosome inference...\n",
            cpus(opts)
        )
        .as_bytes(),
    );

    // ---- the analysis itself -------------------------------------------
    let seglen = seglength_bp(opts);
    // An integer option carries its own "unset": `--degree 0` is not echoed in the banner
    // and does not filter, exactly as an unmentioned `--degree` does not. Asking
    // [`Options::was_given`] instead is silently always-false — it is only ever set for a
    // `Kind::Double`, which is why this filter did nothing at all until it was measured:
    // `--ibdseg --degree 2` on `bigish` reported all 763 pairs against the reference's
    // 442. The rule itself is [`ibdseg::reported_at_degree`].
    let degree = opts.int(Opt::Degree);
    let samples = &loaded.fileset.samples;
    let g = &loaded.fileset.genotypes;
    let mut rows = Vec::new();
    for (i, j) in seg_pair_order(samples.len()) {
        let seg = ibdseg::pair_segments(g, &a.auto_pos, &auto, i, j, seglen);
        if !seg.reported() {
            continue;
        }
        let (pi1, pi2, prop) = (
            seg.ibd1_seg(denom),
            seg.ibd2_seg(denom),
            seg.prop_ibd(denom),
        );
        if !ibdseg::reported_at_degree(degree, pi2, prop) {
            continue;
        }
        rows.push(Row {
            i,
            j,
            fid1: samples[i].fid.clone(),
            id1: samples[i].iid.clone(),
            fid2: samples[j].fid.clone(),
            id2: samples[j].iid.clone(),
            pi1,
            pi2,
            prop,
        });
    }

    let _ = out.write_all(console::ends_at(ENDS_AT_INDENT, console::now_local()).as_bytes());
    let _ = out.write_all(b"\n");

    let seg_path = format!("{prefix}.seg");
    let mut text = String::from("FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n");
    for r in &rows {
        text.push_str(&format!(
            "{}\t{}\t{}\t{}\t{:.4}\t{:.4}\t{:.4}\t{}\n",
            r.fid1,
            r.id1,
            r.fid2,
            r.id2,
            r.pi1,
            r.pi2,
            // `.seg` alone prints PropIBD from the two columns above it rather than from
            // the underlying totals — the reference's own two writers disagree on this
            // and only one rule can match each. See `ibdseg::seg_prop_ibd`.
            ibdseg::seg_prop_ibd(r.pi1, r.pi2),
            ibdseg::inf_type(r.pi1, r.pi2, r.prop),
        ));
    }
    write_file(&seg_path, &text);

    let _ = out.write_all(FILTER_NOTE.as_bytes());
    let _ = out.write_all(
        format!(
            "Summary statistics of IBD segments for individual pairs saved in file {seg_path}\n"
        )
        .as_bytes(),
    );

    // The X table, when this run has one: the same reported pairs measured a second time
    // over the X array. Gate and rules in [`xseg`].
    if let (true, Some(xg)) = (xseg::runs(degree, &x_usable), loaded.x_genotypes.as_ref()) {
        let x = xseg::XSegments::new(&a.x_pos, &x_usable, xg, seglen);
        let x_path = format!("{prefix}X.seg");
        write_file(&x_path, &xseg::text(samples, &pairs(&rows), &x));
        let _ = out.write_all(xseg::saved_line(&x_path).as_bytes());
    }
}

/// The `.fam` index pairs of the reported rows, in the order `<prefix>.seg` printed them.
fn pairs(rows: &[Row]) -> Vec<(usize, usize)> {
    rows.iter().map(|r| (r.i, r.j)).collect()
}

/// Indent that lines this pass's `ends at` up under `IBD segment analysis starts at`.
const ENDS_AT_INDENT: usize = 23;

/// The two filter lines. Both are hard-coded in the reference: the `3` does **not** track
/// `--seglength`, and the `10` is not tunable at all.
const FILTER_NOTE: &str = concat!(
    "Note with relationship inference as the primary goal, the following filters are applied:\n",
    "  Sample pairs without any long IBD segments (>10Mb) are excluded.\n",
    "  Short IBD segments (<3Mb) are not reported/utilized.\n",
);

/// `kingsplitped.txt is generated for certain pedigree plot applications.` — printed
/// *before* the `Options in effect:` block, unlike everything else this pass says.
pub fn splitped_notice(prefix: &str) -> String {
    format!("{prefix}splitped.txt is generated for certain pedigree plot applications.\n\n")
}

/// The autosomal segment estimates for one pair at a time.
///
/// `--ibdseg` sweeps every pair and writes `<prefix>.seg`; `--build` and `--cluster` want
/// the same three numbers for a handful of named pairs and write no `.seg` at all, so the
/// usable-segment construction is done once here and each pair scanned on demand.
pub struct Segments {
    pos: Vec<i64>,
    segs: Vec<Usable>,
    denom: i64,
    seglen: i64,
}

impl Segments {
    /// `None` when the map yields no usable autosomal segment — the case in which the
    /// reference computes nothing either and `--ibdseg` prints `No informative IBD
    /// segments.`.
    pub fn new(opts: &Options, loaded: &Loaded) -> Option<Self> {
        let a = arrays(&loaded.fileset.variants, i64::from(opts.int(Opt::Sexchr)));
        let segs = ibdseg::usable_segments(&a.auto_chr, &a.auto_pos);
        if segs.is_empty() {
            return None;
        }
        let denom = ibdseg::denominator(&segs, &a.auto_pos);
        Some(Segments {
            pos: a.auto_pos,
            segs,
            denom,
            seglen: seglength_bp(opts),
        })
    }

    /// `(IBD1Seg, IBD2Seg, PropIBD)` for one pair, unrounded.
    pub fn of(&self, loaded: &Loaded, i: usize, j: usize) -> (f64, f64, f64) {
        let s = ibdseg::pair_segments(
            &loaded.fileset.genotypes,
            &self.pos,
            &self.segs,
            i,
            j,
            self.seglen,
        );
        (
            s.ibd1_seg(self.denom),
            s.ibd2_seg(self.denom),
            s.prop_ibd(self.denom),
        )
    }
}

/// Samples per block in the `<prefix>.seg` row order — see [`seg_pair_order`].
const SEG_BLOCK: usize = 16;

/// The order `<prefix>.seg` lists its pairs in: **by 16-sample block, then by index**.
///
/// Not `i` ascending then `j` ascending, which is what every other pairwise file uses and
/// what this pass used until it was measured. The reference walks blocks of
/// [`SEG_BLOCK`] samples — for each block `b1`, for each block `b2 >= b1`, every reported
/// pair with `i` in `b1` and `j` in `b2`, `i` then `j` ascending. On a fileset of 16
/// samples or fewer there is one block and the two orders coincide, which is why nine of
/// the corpus's thirteen datasets never showed it.
///
/// `multifam` (20 samples) does: after the within-block pairs run out at `(13, 14)` the
/// reference jumps to `(11, 16)`, deferring every pair that reaches into the second block
/// until the first block is exhausted. `bigish` (200 samples, 13 blocks) shows the same
/// shape thirteen times over.
///
/// **The block size is 16 and nothing else.** Sweeping it over 2..80 against the row
/// order of **all 50 captured `.seg` files**, exactly one value reproduces every one:
/// `threegen` (12 samples) rules out everything below 12, `multifam` rules out everything
/// from 20 up, and inside that window only 16 survives `bigish`.
///
/// It is a pure ordering: the same rows with the same values, and `measure_gaps.py`
/// (which matches rows on their identifier columns before comparing) reported **0 extra
/// and 0 missing** both before and after. Only a byte-for-byte diff can see it, which is
/// why it stayed hidden behind the `PropIBD` residual until `ibdseg::seg_prop_ibd` closed
/// that: with the numbers finally exact on `multifam` and `bigish`, the order was all
/// that was left.
fn seg_pair_order(n: usize) -> impl Iterator<Item = (usize, usize)> {
    let blocks = n.div_ceil(SEG_BLOCK);
    (0..blocks).flat_map(move |b1| {
        (b1..blocks).flat_map(move |b2| {
            let i_hi = ((b1 + 1) * SEG_BLOCK).min(n);
            (b1 * SEG_BLOCK..i_hi).flat_map(move |i| {
                let j_lo = if b1 == b2 { i + 1 } else { b2 * SEG_BLOCK };
                let j_hi = ((b2 + 1) * SEG_BLOCK).min(n);
                (j_lo..j_hi).map(move |j| (i, j))
            })
        })
    })
}

/// `--cpus` when given, otherwise however many cores are visible.
fn cpus(opts: &Options) -> usize {
    match opts.int(Opt::Cpus) {
        n if n > 0 => n as usize,
        _ => std::thread::available_parallelism().map_or(1, |n| n.get()),
    }
}

fn write_file(path: &str, text: &str) {
    let _ = std::fs::write(PathBuf::from(path), text);
}

/// The segment pre-pass `--unrelated`, `--bysample` and `--bySNP` share.
///
/// All three run the same usable-segment construction `--ibdseg` does and write
/// `<prefix>allsegs.txt` from it — the two QC reports silently, `--unrelated` with the
/// console block this returns. A run whose map yields no usable segment writes no file
/// and reports nothing.
pub fn segment_prepass(opts: &Options, loaded: &Loaded) -> String {
    let prefix = opts.string(Opt::Prefix).to_string();
    let sexchr = i64::from(opts.int(Opt::Sexchr));
    if map_order_warning(&loaded.fileset.variants, sexchr).is_some() {
        return String::new();
    }
    let a = arrays(&loaded.fileset.variants, sexchr);
    let auto = ibdseg::usable_segments(&a.auto_chr, &a.auto_pos);
    let xseg = ibdseg::usable_segments(&a.x_chr, &a.x_pos);
    if auto.is_empty() {
        return String::new();
    }
    let path = format!("{prefix}allsegs.txt");
    write_file(
        &path,
        &allsegs_text(&loaded.fileset.variants, &a, &auto, &xseg),
    );
    let denom = ibdseg::denominator(&auto, &a.auto_pos);
    let mut s = format!(
        "Total length of {} chromosomal segments usable for IBD segment analysis is {:.1} Mb.\n",
        auto.len(),
        denom as f64 / 1e6
    );
    if !xseg.is_empty() {
        let xlen = ibdseg::denominator(&xseg, &a.x_pos);
        s.push_str(&format!(
            "  In addition to autosomes, {} segments of length {:.1} Mb on X-chr can be further used.\n",
            xseg.len(),
            xlen as f64 / 1e6
        ));
    }
    s.push_str(&format!(
        "  Information of these chromosomal segments can be found in file {path}\n\n"
    ));
    s
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli;

    fn parse(args: &[&str]) -> Options {
        let v: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        cli::parse(&v).options
    }

    fn marker(chrom: &str, id: &str, bp: i64) -> Variant {
        Variant {
            chrom: chrom.to_string(),
            id: id.to_string(),
            cm: 0.0,
            bp,
            a1: "A".to_string(),
            a2: "G".to_string(),
        }
    }

    #[test]
    fn map_order_reports_the_first_position_or_chromosome_regression() {
        let positions = [marker("1", "a", 20), marker("1", "b", 10)];
        assert_eq!(
            map_order_warning(&positions, 23).as_deref(),
            Some("Positions unsorted: a at 20, b at 10.\n")
        );

        let chromosomes = [marker("22", "z", 20), marker("21", "y", 10)];
        assert_eq!(
            map_order_warning(&chromosomes, 23).as_deref(),
            Some("Chromosomes unsorted: z on chr 22, y on chr 21.\n")
        );

        let sorted = [
            marker("1", "a", 10),
            marker("1", "b", 20),
            marker("2", "c", 1),
        ];
        assert_eq!(map_order_warning(&sorted, 23), None);
    }

    #[test]
    fn the_small_sample_gate_is_five() {
        assert!(downgrades_to_kinship(4));
        assert!(!downgrades_to_kinship(5));
    }

    #[test]
    fn no_segments_note_matches_the_reference() {
        assert_eq!(
            NO_SEGMENTS,
            concat!(
                "No informative IBD segments.\n",
                "  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.\n",
            )
        );
    }

    #[test]
    fn seglength_clamps_to_three_megabases_outside_its_range() {
        assert_eq!(seglength_bp(&parse(&["--ibdseg"])), 3_000_000);
        assert_eq!(
            seglength_bp(&parse(&["--ibdseg", "--seglength", "5"])),
            5_000_000
        );
        assert_eq!(
            seglength_bp(&parse(&["--ibdseg", "--seglength", "0"])),
            3_000_000
        );
        assert_eq!(
            seglength_bp(&parse(&["--ibdseg", "--seglength", "11"])),
            3_000_000
        );
        assert_eq!(
            seglength_bp(&parse(&["--ibdseg", "--seglength", "10"])),
            10_000_000
        );
    }

    #[test]
    fn the_seg_pair_order_is_plain_index_order_within_one_block() {
        let got: Vec<_> = seg_pair_order(4).collect();
        assert_eq!(
            got,
            vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)],
            "16 samples or fewer is one block, so nothing is deferred"
        );
        // Every pair, exactly once, at the block boundary too.
        for n in [1, 15, 16, 17, 33] {
            let v: Vec<_> = seg_pair_order(n).collect();
            assert_eq!(v.len(), n * n.saturating_sub(1) / 2, "n={n}");
            let mut s = v.clone();
            s.sort_unstable();
            s.dedup();
            assert_eq!(s.len(), v.len(), "n={n}: a pair was emitted twice");
            assert!(v.iter().all(|&(i, j)| i < j), "n={n}");
        }
    }

    /// The shape `multifam` exposed: with 20 samples the first block is finished before
    /// any pair reaching into the second one is written.
    #[test]
    fn the_seg_pair_order_defers_cross_block_pairs() {
        let v: Vec<_> = seg_pair_order(20).collect();
        let at = |p: (usize, usize)| v.iter().position(|&q| q == p).unwrap();
        assert!(at((13, 14)) < at((11, 16)), "within-block pairs come first");
        assert!(at((11, 19)) < at((12, 16)), "cross-block runs i then j");
        assert!(
            at((15, 19)) < at((16, 17)),
            "the second block's own pairs are last"
        );
        assert_eq!(v[0], (0, 1));
        assert_eq!(*v.last().unwrap(), (18, 19));
    }

    #[test]
    fn seglength_is_not_echoed_in_options_in_effect() {
        assert_eq!(
            options_in_effect(&parse(&["--ibdseg", "--seglength", "5"])),
            vec!["--ibdseg".to_string()]
        );
        assert_eq!(
            options_in_effect(&parse(&["--ibdseg", "--degree", "2"])),
            vec!["--ibdseg".to_string(), "--degree 2".to_string()]
        );
    }
}
