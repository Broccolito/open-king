//! `--ibdseg` — pairwise IBD-segment inference.
//!
//! The pass owns three files and one console body:
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
//! ```
//!
//! Two lines in there are hard-coded in the reference and must be copied verbatim even
//! when they are wrong: the `(<3Mb)` still says `3` under `--seglength 10`, and the
//! `(>10Mb)` pair filter is not tunable at all.
//!
//! # Known gap
//!
//! `<prefix>splitped.txt` is **announced but not written**. It is a pedigree-splitting
//! artefact, not a segment artefact: it renames families that turn out to be several
//! disconnected pedigrees (`POOL` → `POOL_S1`, `POOL_S2`, …), imports a founder's
//! genotyped parents into the family that references them, drops families that contribute
//! no informative pair, and reorders members. Reproducing it is a pedigree-reconstruction
//! task shared with `--build`, and it is byte-identical under every `--degree` /
//! `--seglength` / `--related` variant, so it can be added later without touching this
//! module's numbers.

use std::io::Write;
use std::path::PathBuf;

use king_core::ibdseg::{self, Usable};
use king_io::Variant;

use crate::cli::{Opt, Options};
use crate::console;
use crate::load::{self, Class, Loaded};

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
    let a = arrays(&loaded.fileset.variants, sexchr);
    let auto = ibdseg::usable_segments(&a.auto_chr, &a.auto_pos);
    let xseg = ibdseg::usable_segments(&a.x_chr, &a.x_pos);

    if auto.is_empty() {
        // No denominator, so nothing downstream can be computed and no file is written.
        let _ = out.write_all(b"No informative IBD segments.\n");
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
    if !xseg.is_empty() {
        let xlen = ibdseg::denominator(&xseg, &a.x_pos);
        let _ = out.write_all(
            format!(
                "  In addition to autosomes, {} segments of length {:.1} Mb on X-chr can be further used.\n",
                xseg.len(),
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
        &allsegs_text(&loaded.fileset.variants, &a, &auto, &xseg),
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
    let cutoff = if opts.was_given(Opt::Degree) {
        Some(ibdseg::degree_cutoff(opts.int(Opt::Degree)))
    } else {
        None
    };
    let samples = &loaded.fileset.samples;
    let g = &loaded.fileset.genotypes;
    let mut rows = Vec::new();
    for i in 0..samples.len() {
        for j in i + 1..samples.len() {
            let seg = ibdseg::pair_segments(g, &a.auto_pos, &auto, i, j, seglen);
            if !seg.reported() {
                continue;
            }
            let (pi1, pi2, prop) = (
                seg.ibd1_seg(denom),
                seg.ibd2_seg(denom),
                seg.prop_ibd(denom),
            );
            if let Some(c) = cutoff {
                if prop <= c {
                    continue;
                }
            }
            rows.push(Row {
                fid1: samples[i].fid.clone(),
                id1: samples[i].iid.clone(),
                fid2: samples[j].fid.clone(),
                id2: samples[j].iid.clone(),
                pi1,
                pi2,
                prop,
            });
        }
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
            r.prop,
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

    #[test]
    fn the_small_sample_gate_is_five() {
        assert!(downgrades_to_kinship(4));
        assert!(!downgrades_to_kinship(5));
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
