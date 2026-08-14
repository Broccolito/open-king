//! Loading a PLINK binary fileset, and the console sequence the reference emits while it
//! happens.
//!
//! Everything here was established by running the reference binary, never by reading its
//! source. The module owns three things the readers in `king-io` deliberately do not:
//!
//! 1. **The chromosome partition**, which depends on `--sexchr` and therefore cannot be
//!    expressed as a [`VariantFilter`] — see [`classify`].
//! 2. **The order of the checks**, which is observable: some failures print before
//!    `Loading genotype data …` and some after `Read in PLINK bed file …`.
//! 3. **The progress ticks**, whose count and placement are a pure function of the
//!    retained-SNP count — see [`Ticks`].
//!
//! # The sequence
//!
//! ```text
//! [probe the .bed — skipped entirely on the --risk --model path]
//! Loading genotype data in PLINK binary format...
//! Read in PLINK fam file <fam>...
//!   PLINK pedigrees loaded: <n> samples
//! Read in PLINK bim file <bim>...
//!   Genotype data consist of <a> autosome SNPs[ (including <xy> XY SNPs)][, …]
//! [  <r> other SNPs are removed.]
//!   PLINK maps loaded: <a+x+y+mt> SNPs
//! Read in PLINK bed file <bed>...
//! 0%\r6%\r…88%\r  PLINK binary genotypes loaded.
//! 94%\r  KING format genotype data successfully converted.
//! ```
//!
//! The last two lines are the trap: the ticks are **not** all emitted before
//! `loaded.`. See [`Ticks`] for the rule and the sweep that pinned it.
//!
//! `Autosome genotypes stored in <w> words for each of <n> individuals.` is *not* part of
//! this sequence — the reference reprints it once per analysis pass, so it belongs to the
//! dispatcher. [`Loaded::preamble`] renders it.

use std::fs::File;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};

use king_io::{bed, bim, fam, Fileset, IoError, Variant, VariantFilter};

use crate::cli::{Opt, Options};
use crate::console;

/// The body of a FATAL ERROR, ready for `console::fatal_block`.
pub type Fatal = String;

// ---------------------------------------------------------------------------
// Chromosome partition
// ---------------------------------------------------------------------------

/// Which class of the reference's chromosome partition a `.bim` row falls into.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Class {
    /// Chromosomes `1 ..= sexchr-1`. Analysed, and counted in the autosome total.
    Autosome,
    /// Chromosome `sexchr+2`, the pseudo-autosomal `XY`. Also analysed as an autosome,
    /// but additionally reported in the `(including <n> XY SNPs)` parenthetical.
    Xy,
    /// Chromosome `sexchr`.
    X,
    /// Chromosome `sexchr+1`.
    Y,
    /// Chromosome `sexchr+3`.
    Mt,
    /// Everything else, including chromosome `0` and unparseable contig names. Dropped,
    /// and reported by `<n> other SNPs are removed.`
    Removed,
}

impl Class {
    /// Whether the class enters the autosomal genotype matrix.
    pub fn is_autosomal(self) -> bool {
        matches!(self, Class::Autosome | Class::Xy)
    }
}

/// Numeric chromosome code for a `.bim` label, under a given `--sexchr`.
///
/// Two rules, in this order, both read off the reference:
///
/// * `X`, `Y`, `XY` and `MT` (case insensitive, whole label) bind to the **class**, not to
///   a fixed number: under `--sexchr 5` a row labelled `X` is chromosome 5, while a row
///   labelled `23` is removed. A label merely *starting* with one of these names does not
///   match — `X1` is not `X`.
/// * Otherwise the leading run of ASCII digits is read, with **no sign accepted**. `1abc`,
///   `01`, `1.5` and `2x` are all chromosome 1, 1, 1 and 2; `+1`, `-1`, `chr1` and
///   `unplaced` all yield 0 and are removed.
///
/// Note in particular that a `chr` prefix is **not** stripped — `chr1` is dropped, which a
/// tolerant reader would silently get wrong.
pub fn chromosome_code(label: &str, sexchr: i64) -> i64 {
    if label.eq_ignore_ascii_case("X") {
        return sexchr;
    }
    if label.eq_ignore_ascii_case("Y") {
        return sexchr + 1;
    }
    if label.eq_ignore_ascii_case("XY") {
        return sexchr + 2;
    }
    if label.eq_ignore_ascii_case("MT") {
        return sexchr + 3;
    }
    leading_digits(label)
}

/// The leading run of ASCII digits as a number, or 0 when there is none.
///
/// Saturates rather than wrapping, so an absurdly long digit string cannot alias onto a
/// real chromosome code.
fn leading_digits(label: &str) -> i64 {
    let mut n: i64 = 0;
    for c in label.bytes() {
        match c.checked_sub(b'0').filter(|d| *d < 10) {
            Some(d) => n = n.saturating_mul(10).saturating_add(i64::from(d)),
            None => break,
        }
    }
    n
}

/// Classify a `.bim` label under a given `--sexchr`.
///
/// The partition, fully determined by sweeping `--sexchr` over a fileset carrying every
/// chromosome at a distinct count:
///
/// ```text
/// autosomes = 1 ..= sexchr-1,  plus sexchr+2 (XY)
/// X = sexchr      Y = sexchr+1      XY = sexchr+2      MT = sexchr+3
/// everything else, chromosome 0 included, is removed
/// ```
///
/// At the default `--sexchr 23` that is `1..=22` **plus 25** — the easy half of the rule
/// to miss, and it changes the autosomal SNP count on any real fileset with PAR markers.
pub fn classify(label: &str, sexchr: i64) -> Class {
    let code = chromosome_code(label, sexchr);
    if code == sexchr + 2 {
        Class::Xy
    } else if code == sexchr {
        Class::X
    } else if code == sexchr + 1 {
        Class::Y
    } else if code == sexchr + 3 {
        Class::Mt
    } else if code >= 1 && code < sexchr {
        Class::Autosome
    } else {
        Class::Removed
    }
}

/// How a map broke down across the partition.
///
/// `autosome` is the figure the reference prints and **includes** `xy`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct MapCounts {
    pub autosome: usize,
    pub xy: usize,
    pub x: usize,
    pub y: usize,
    pub mt: usize,
    pub removed: usize,
}

impl MapCounts {
    /// Count every classified variant — what `PLINK maps loaded: <n> SNPs` reports.
    ///
    /// This is **not** the number of lines in the `.bim`: removed variants are excluded,
    /// so a map with unplaced contigs reports fewer SNPs loaded than it has rows.
    pub fn loaded(&self) -> usize {
        self.autosome + self.x + self.y + self.mt
    }

    fn tally(&mut self, class: Class) {
        match class {
            Class::Autosome => self.autosome += 1,
            Class::Xy => {
                self.autosome += 1;
                self.xy += 1;
            }
            Class::X => self.x += 1,
            Class::Y => self.y += 1,
            Class::Mt => self.mt += 1,
            Class::Removed => self.removed += 1,
        }
    }
}

/// Classify a whole map, returning the per-variant keep flags and the tallies.
pub fn partition(variants: &[Variant], sexchr: i64) -> (Vec<bool>, MapCounts) {
    let mut counts = MapCounts::default();
    let mut keep = Vec::with_capacity(variants.len());
    for v in variants {
        let class = classify(&v.chrom, sexchr);
        counts.tally(class);
        keep.push(class.is_autosomal());
    }
    (keep, counts)
}

// ---------------------------------------------------------------------------
// Progress ticks
// ---------------------------------------------------------------------------

/// The tick schedule for a load of `m` retained SNPs.
///
/// The reference reads the `.bed` in blocks of whole 64-SNP words and prints one
/// `<pct>%\r` tick **after** each block it completes. Three constants fall out of a sweep
/// over `m` from 1 to 50 000, cross-checked against the `sexchr` fixture:
///
/// * the block is [`Ticks::words_per_block`] `= ceil(w / 16)` words, where `w` is the word
///   count — so the run aims for 16 ticks and rounds the block size up;
/// * the percentage is always `round(100 * i / 16)`, with the denominator **fixed at 16**
///   however many ticks actually occur — a 25-word load prints `0 6 13 … 75`, not
///   `0 8 15 … 92`;
/// * the last tick lands on the far side of `PLINK binary genotypes loaded.` whenever the
///   final block is short, i.e. whenever `m` is not a multiple of the block.
///
/// That last point is the one worth restating: for `m = 5000` the reference emits fifteen
/// ticks, then `loaded.`, then `94%`, then `converted.`. For `m = 960` — a whole number of
/// blocks — all fifteen ticks precede `loaded.` and no tick follows it. For `m = 1` there
/// are no ticks before `loaded.` at all and a lone `0%` after it. Each case was captured
/// from the reference.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Ticks {
    /// Retained SNPs per block.
    pub snps_per_block: usize,
    /// Blocks that complete during the read, i.e. ticks printed before `loaded.`.
    pub before: usize,
    /// Whether a short final block adds one tick after `loaded.`.
    pub after: bool,
}

/// The tick denominator. Constant, whatever the real block count turns out to be.
const TICK_DENOMINATOR: usize = 16;
/// SNPs in one 64-bit word of a bit plane.
const SNPS_PER_WORD: usize = 64;

impl Ticks {
    /// Work out the schedule for `m` retained SNPs.
    pub fn plan(m: usize) -> Ticks {
        let words = m.div_ceil(SNPS_PER_WORD);
        let snps_per_block = SNPS_PER_WORD * Self::words_per_block(words);
        Ticks {
            snps_per_block,
            before: m / snps_per_block,
            after: m % snps_per_block != 0,
        }
    }

    /// `ceil(w / 16)`, floored at one word so a zero-SNP map cannot divide by zero.
    fn words_per_block(words: usize) -> usize {
        words.div_ceil(TICK_DENOMINATOR).max(1)
    }

    /// The `i`-th tick, `<pct>%\r`.
    fn tick(i: usize) -> String {
        console::progress(i, TICK_DENOMINATOR)
    }

    /// Emit the whole schedule, `loaded.` and `converted.` included.
    fn emit(self, out: &mut dyn Write) {
        for i in 0..self.before {
            write(out, &Self::tick(i));
        }
        write(out, &console::binary_genotypes_loaded());
        if self.after {
            write(out, &Self::tick(self.before));
        }
        write(out, &console::genotypes_converted());
    }

    /// Emit only the ticks that would have completed before `retained` SNPs had been
    /// packed — what a truncated `.bed` shows before its fatal.
    fn emit_partial(self, retained: usize, out: &mut dyn Write) {
        for i in 0..(retained / self.snps_per_block) {
            write(out, &Self::tick(i));
        }
    }
}

// ---------------------------------------------------------------------------
// Loading
// ---------------------------------------------------------------------------

/// A fileset plus the map tallies the console reported for it.
pub struct Loaded {
    pub fileset: Fileset,
    pub counts: MapCounts,
    /// The X-chromosome bit planes, decoded only when the map carries enough X markers
    /// for the X pass to run at all ([`X_PASS_MIN_SNPS`]).
    pub x_genotypes: Option<king_io::Genotypes>,
}

/// Fewest X-chromosome SNPs that will start the X-chromosome pass.
///
/// **512**, bisected against the reference on a fileset whose X marker count was swept:
/// 511 produces no `X.kin`/`X.kin0` and no `X-chromosome analysis...` block, 512 produces
/// both. It is a count of *markers*, not of calls — 512 markers of which 200 are missing
/// in every sample still runs. This is what keeps `--sexchr 24` (300 X SNPs on the
/// `sexchr` fixture) and `--sexchr 25` (150) silent while `--sexchr 2` (2 000) is not.
pub const X_PASS_MIN_SNPS: usize = 512;

impl Loaded {
    /// Words per sample in the autosomal bit planes.
    pub fn words(&self) -> usize {
        self.fileset.genotypes.words_per_sample()
    }

    /// `Autosome genotypes stored in <w> words for each of <n> individuals.` plus the
    /// blank line under it.
    ///
    /// This is **per analysis pass, not per load**, and only for the analyses in
    /// [`prints_preamble`]. A `--kinship --ibs` run emits it twice — once before each
    /// pass's own `Options in effect:` block — while a `--build` run never emits it at
    /// all. Rendering it unconditionally after the load looks right on `--kinship` and is
    /// wrong on five of the analyses.
    pub fn preamble(&self) -> String {
        format!(
            "{}\n",
            console::autosome_words(self.words(), self.fileset.samples.len())
        )
    }
}

/// Whether an analysis opens its pass with [`Loaded::preamble`].
///
/// Probed one flag at a time against the reference on the same fileset:
///
/// | Prints it | Does not |
/// | --- | --- |
/// | `--related` `--kinship` `--ibs` `--duplicate` | `--ibdseg` `--mds` `--pca` `--cluster` `--build` `--bysample` `--bysnp` `--tdt` `--unrelated` `--roh` `--gdt` `--makeGRM` `--lmm` `--risk` `--rplot` `--pngplot` |
///
/// `--autoQC` prints a line of the same shape but **not this one**: its word count is
/// `ceil(m / 16)`, not `ceil(m / 64)` — 313 against 79 on a 5 000-SNP map, 625 against
/// 157 on 10 000 — so it is reporting its own denser packing and must render the line
/// itself rather than borrow the loader's.
pub fn prints_preamble(opt: Opt) -> bool {
    matches!(opt, Opt::Related | Opt::Kinship | Opt::Ibs | Opt::Duplicate)
}

/// Whether the reference probes the `.bed` before it starts loading.
///
/// It does, for every analysis — except `--risk` with a `--model`, which goes straight to
/// the `.fam` and fails there instead. That is the reference's own inconsistency, not
/// ours: `king -b nonexistent.bed --risk --model m.txt` reports
/// `Pedigree file nonexistent.fam cannot be opened`, and it keeps doing so when
/// `--kinship` is added alongside.
pub fn probes_bed(opts: &Options) -> bool {
    !opts.flag(Opt::Risk) || opts.string(Opt::Model).is_empty()
}

/// The three fileset paths, honouring `--fam` and `--bim`.
///
/// Both overrides are used **verbatim**: `--fam alt` reads a file literally named `alt`
/// and reports `Pedigree file alt cannot be opened` when it is missing. Only the derived
/// names get a suffix appended, and the `.bed` name itself is always `opts.bed` exactly as
/// typed.
pub fn paths(opts: &Options) -> (PathBuf, PathBuf, PathBuf) {
    let (bed_path, default_bim, default_fam) = bed::paths_for(Path::new(&opts.bed));
    let fam_path = match opts.string(Opt::Fam) {
        "" => default_fam,
        s => PathBuf::from(s),
    };
    let bim_path = match opts.string(Opt::Bim) {
        "" => default_bim,
        s => PathBuf::from(s),
    };
    (bed_path, bim_path, fam_path)
}

/// Load the fileset named on the command line, printing what the reference prints.
///
/// Returns the fatal message body on failure; the caller wraps it in
/// `console::fatal_block` and exits 1.
pub fn load(opts: &Options, out: &mut dyn Write) -> Result<Loaded, Fatal> {
    let (bed_path, bim_path, fam_path) = paths(opts);
    let sexchr = i64::from(opts.int(Opt::Sexchr));

    if probes_bed(opts) {
        probe_bed(&opts.bed)?;
    }

    write(out, &console::loading_plink_binary());

    // --- .fam ------------------------------------------------------------
    // The reference opens first and announces second: a missing `.fam` produces the fatal
    // with no `Read in PLINK fam file …` line above it.
    openable(&fam_path).map_err(|_| console::pedigree_file_unopenable(&display(&fam_path)))?;
    write(out, &console::read_fam(&display(&fam_path)));
    // The reference rewrites the pedigree through `<prefix>$TMP$.ped` and opens that file
    // here, before it has parsed a single row: with `--prefix sub/pre` and a duplicated
    // sample in the `.fam`, the unwritable-prefix fatal is what comes out, so the open
    // precedes even the duplicate check. It is a real writability test rather than a
    // "does the directory exist" one — a read-only directory fails it too — and the
    // reference removes the file again, leaving nothing behind on a successful run.
    probe_tmp_ped(opts.string(Opt::Prefix))?;
    let samples = fam::read_fam(&fam_path).map_err(io_fatal)?;
    if let Err(IoError::DuplicateSample { fid, iid }) = fam::check_duplicates(&samples) {
        // Printed before the fatal, and `PLINK pedigrees loaded` never appears.
        write(out, &console::person_duplicated(&fid, &iid));
        return Err(console::PEDIGREE_STRUCTURE_PROBLEMS.to_string());
    }
    write(out, &console::pedigrees_loaded(samples.len()));

    // --- .bim ------------------------------------------------------------
    openable(&bim_path).map_err(|_| console::map_file_unopenable(&display(&bim_path)))?;
    write(out, &console::read_bim(&display(&bim_path)));
    let variants = bim::read_bim(&bim_path).map_err(io_fatal)?;
    let (keep, counts) = partition(&variants, sexchr);
    write(
        out,
        &console::genotype_data_counts(counts.autosome, counts.xy, counts.x, counts.y, counts.mt),
    );
    if counts.removed > 0 {
        write(out, &console::other_snps_removed(counts.removed));
    }
    write(out, &console::maps_loaded(counts.loaded()));
    if counts.autosome == 0 {
        // Reported only after the whole map has been announced.
        return Err(console::NO_AUTOSOME_SNPS.to_string());
    }

    // --- .bed ------------------------------------------------------------
    write(out, &console::read_bed(&display(&bed_path)));
    let ticks = Ticks::plan(counts.autosome);
    let decoded = decode(&bed_path, samples.len(), &variants, &keep);
    let (genotypes, kept) = match decoded {
        Ok(v) => v,
        Err(DecodeError::Truncated { marker }) => {
            // The reference has already printed a tick for every block it finished, so
            // replay exactly those before failing.
            let retained = keep[..marker].iter().filter(|k| **k).count();
            ticks.emit_partial(retained, out);
            return Err(console::not_enough_genotypes(marker));
        }
        Err(DecodeError::NotSnpMajor) => return Err(console::SNP_MAJOR_ONLY.to_string()),
        Err(DecodeError::NotBinary) => {
            return Err(console::PLINK_OR_KING_BINARY_REQUIRED.to_string())
        }
        Err(DecodeError::Io(e)) => return Err(e),
    };
    ticks.emit(out);

    // The X planes cost a second pass over the `.bed`, so take it only when the X pass
    // could run. The autosomal decode has already validated the file's length and header,
    // so this one cannot introduce a new failure mode.
    let x_genotypes = (counts.x >= X_PASS_MIN_SNPS)
        .then(|| {
            let x_keep: Vec<bool> = variants
                .iter()
                .map(|v| classify(&v.chrom, sexchr) == Class::X)
                .collect();
            decode(&bed_path, samples.len(), &variants, &x_keep)
                .ok()
                .map(|(g, _)| g)
        })
        .flatten();

    Ok(Loaded {
        fileset: Fileset {
            samples,
            variants,
            genotypes,
            kept,
        },
        counts,
        x_genotypes,
    })
}

/// The `.bed` checks that run *before* `Loading genotype data …`.
///
/// Three of them, in this order — the order is observable, because each produces a
/// different message on the same input:
///
/// 1. the file must open (`-b missing.bed` → `Genotype file missing.bed cannot be opened`);
/// 2. the argument must end in a literal lowercase `.bed` (a readable file named `plain`
///    → `Please use PLINK binary format as input.`);
/// 3. the first two bytes must be the PLINK magic (a zero-length or garbage-headed file
///    → `Please use either PLINK or KING binary format as input.`).
///
/// The **mode** byte is deliberately not checked here: an individual-major `.bed` gets all
/// the way to `Read in PLINK bed file …` before failing.
fn probe_bed(arg: &str) -> Result<(), Fatal> {
    let mut file = File::open(arg).map_err(|_| console::genotype_file_unopenable(arg))?;
    if !arg.ends_with(".bed") {
        return Err(console::PLINK_BINARY_REQUIRED.to_string());
    }
    let mut magic = [0u8; 2];
    match file.read_exact(&mut magic) {
        Ok(()) if magic == bed::MAGIC => Ok(()),
        _ => Err(console::PLINK_OR_KING_BINARY_REQUIRED.to_string()),
    }
}

/// Why a `.bed` could not be decoded, in the reference's own terms.
enum DecodeError {
    /// Fewer bytes than the map and sample table call for; `marker` is the 0-based index
    /// of the first row that is not fully present.
    Truncated {
        marker: usize,
    },
    NotSnpMajor,
    NotBinary,
    Io(Fatal),
}

/// Decode the `.bed` into bit planes carrying only the variants `keep` marks.
///
/// Two length cases the readers in `king-io` treat as errors are handled here first,
/// because the reference does not treat them alike: a **short** `.bed` is fatal with
/// `Not enough genotypes at the <k>th marker`, while a **long** one is read happily and
/// its trailing bytes ignored. Only the over-long case pays for buffering.
fn decode(
    path: &Path,
    n_samples: usize,
    variants: &[Variant],
    keep: &[bool],
) -> Result<(king_io::Genotypes, Vec<usize>), DecodeError> {
    let probe = keep_list(keep);
    let expected = bed::expected_bed_len(n_samples, variants.len());
    let found = std::fs::metadata(path).map(|m| m.len()).unwrap_or(0);

    let result = if found < expected {
        // Distinguish "not a .bed at all" from "a .bed that stops early": the header has
        // to be intact before a marker index means anything.
        if found < bed::HEADER_LEN as u64 {
            return Err(DecodeError::NotBinary);
        }
        let row = bed::bytes_per_variant(n_samples) as u64;
        let complete = (found - bed::HEADER_LEN as u64) / row.max(1);
        return Err(DecodeError::Truncated {
            marker: complete as usize,
        });
    } else if found > expected {
        let mut bytes = vec![0u8; expected as usize];
        File::open(path)
            .and_then(|mut f| f.read_exact(&mut bytes))
            .map_err(|e| DecodeError::Io(e.to_string()))?;
        bed::read_bed_bytes(&bytes, path, n_samples, &probe, VariantFilter::Autosomes)
    } else {
        bed::read_bed(path, n_samples, &probe, VariantFilter::Autosomes)
    };

    result.map_err(|e| match e {
        IoError::IndividualMajor { .. } | IoError::BadMode { .. } => DecodeError::NotSnpMajor,
        IoError::BadMagic { .. } => DecodeError::NotBinary,
        other => DecodeError::Io(other.to_string()),
    })
}

/// Encode a keep list as the stand-in `.bim` column `king-io` filters on.
///
/// `king-io` selects variants through [`VariantFilter`], whose chromosome model stops at
/// code 26 and strips a `chr` prefix. Neither matches the reference — `--sexchr 30` puts
/// MT at 33, and `chr1` is a *removed* contig, not chromosome 1 — so the partition is
/// decided here (see [`classify`]) and `king-io` is handed the resolved answer: chromosome
/// `1` for the rows to keep and `0`, which is in no class, for the rows to drop. Only the
/// chromosome field and the vector's length are read back out, so the rest is left empty.
fn keep_list(keep: &[bool]) -> Vec<Variant> {
    keep.iter()
        .map(|&k| Variant {
            chrom: if k { "1" } else { "0" }.to_string(),
            id: String::new(),
            cm: 0.0,
            bp: 0,
            a1: String::new(),
            a2: String::new(),
        })
        .collect()
}

/// Whether a path can be opened for reading, which is the only thing the reference tests
/// before it announces the file.
fn openable(path: &Path) -> std::io::Result<()> {
    File::open(path).map(|_| ())
}

/// Reproduce the reference's `<prefix>$TMP$.ped` writability probe.
///
/// Create the file the reference creates, then remove it. Creating rather than merely
/// testing the directory is what makes `--prefix ../read-only-dir/x` fail the way the
/// reference fails it; removing it again is what keeps a successful run from leaving an
/// output file the reference never leaves.
fn probe_tmp_ped(prefix: &str) -> Result<(), Fatal> {
    let path = format!("{prefix}$TMP$.ped");
    match File::create(&path) {
        Ok(_) => {
            let _ = std::fs::remove_file(&path);
            Ok(())
        }
        Err(_) => Err(console::cannot_open_tmp_ped(prefix)),
    }
}

/// A path as the reference prints it: the string it was given, unresolved.
fn display(path: &Path) -> String {
    path.to_string_lossy().into_owned()
}

/// Malformed-input errors `king-io` reports and the reference does not, surfaced as a
/// fatal rather than silently mis-parsed.
///
/// The reference reads a short `.fam` line or a ragged `.bim` without complaint and
/// analyses the wreckage; `docs/SPEC.md` §3.7 records that as a deliberate divergence, so
/// these paths trade a silent wrong answer for a diagnosable stop.
fn io_fatal(e: IoError) -> Fatal {
    e.to_string()
}

/// Write and flush, so ordering survives the process exits around it.
fn write(out: &mut dyn Write, s: &str) {
    let _ = out.write_all(s.as_bytes());
    let _ = out.flush();
}

#[cfg(test)]
mod tests {
    use super::*;

    fn classes(labels: &[&str], sexchr: i64) -> Vec<Class> {
        labels.iter().map(|l| classify(l, sexchr)).collect()
    }

    #[test]
    fn default_partition_folds_xy_into_the_autosomes() {
        // Captured: a map of 2000+2000 on chr1/chr2, 1500 X, 300 Y, 150 XY, 50 MT reports
        // "4150 autosome SNPs (including 150 XY SNPs), 1500 X…, 300 Y…, 50 mitochondrial".
        assert_eq!(
            classes(&["1", "22", "23", "24", "25", "26", "0", "27"], 23),
            [
                Class::Autosome,
                Class::Autosome,
                Class::X,
                Class::Y,
                Class::Xy,
                Class::Mt,
                Class::Removed,
                Class::Removed,
            ]
        );
    }

    #[test]
    fn sexchr_shifts_the_whole_partition() {
        // --sexchr 5 over {0:7,1:11,2:13,3:17,20:19,21:23,22:29,23:31,24:37,25:41,26:43,
        // 27:47,30:53} reported 41 autosome SNPs and 330 removed: only chr1-3 survive.
        assert_eq!(
            classes(&["1", "4", "5", "6", "7", "8", "9", "23", "25"], 5),
            [
                Class::Autosome,
                Class::Autosome,
                Class::X,
                Class::Y,
                Class::Xy,
                Class::Mt,
                Class::Removed,
                Class::Removed,
                Class::Removed,
            ]
        );
        // --sexchr 22 put 37 XY SNPs (chr24) inside a 120-SNP autosome total.
        assert_eq!(
            classes(&["21", "24", "22", "23", "25"], 22),
            [Class::Autosome, Class::Xy, Class::X, Class::Y, Class::Mt,]
        );
    }

    #[test]
    fn names_bind_to_the_class_not_to_a_number() {
        // Under --sexchr 5 the reference put "X"/"x" in the X count (8 = 3+5) and dropped
        // the rows labelled "23" entirely.
        assert_eq!(classify("X", 5), Class::X);
        assert_eq!(classify("x", 5), Class::X);
        assert_eq!(classify("23", 5), Class::Removed);
        assert_eq!(classify("XY", 5), Class::Xy);
        assert_eq!(classify("mt", 5), Class::Mt);
        // ...and at the default they coincide with their usual numbers.
        assert_eq!(classify("X", 23), classify("23", 23));
        assert_eq!(classify("mt", 23), classify("26", 23));
    }

    #[test]
    fn labels_parse_like_the_reference() {
        // Captured totals: "1"(2) "1abc"(3) "01"(5) "1.5"(13) "2x"(17) "0002"(19) summed
        // to the 59 autosome SNPs reported, while "+1"(7) "-1"(11) "X1"(29) made up the
        // 47 removed.
        assert_eq!(chromosome_code("1abc", 23), 1);
        assert_eq!(chromosome_code("01", 23), 1);
        assert_eq!(chromosome_code("1.5", 23), 1);
        assert_eq!(chromosome_code("2x", 23), 2);
        assert_eq!(chromosome_code("0002", 23), 2);
        assert_eq!(chromosome_code("+1", 23), 0);
        assert_eq!(chromosome_code("-1", 23), 0);
        assert_eq!(chromosome_code("X1", 23), 0);
        // A `chr` prefix is a contig name, not chromosome 1.
        assert_eq!(classify("chr1", 23), Class::Removed);
        assert_eq!(classify("chr2", 23), Class::Removed);
        assert_eq!(classify("unplaced", 23), Class::Removed);
        // "M" alone is not mitochondrial; only "MT" is.
        assert_eq!(classify("M", 23), Class::Removed);
        assert_eq!(classify("m", 23), Class::Removed);
    }

    #[test]
    fn maps_loaded_excludes_removed_variants() {
        let variants: Vec<Variant> = ["1", "25", "23", "24", "26", "0", "chr1"]
            .iter()
            .map(|c| Variant {
                chrom: c.to_string(),
                id: String::new(),
                cm: 0.0,
                bp: 0,
                a1: String::new(),
                a2: String::new(),
            })
            .collect();
        let (keep, counts) = partition(&variants, 23);
        assert_eq!(keep, [true, true, false, false, false, false, false]);
        assert_eq!(counts.autosome, 2);
        assert_eq!(counts.xy, 1);
        assert_eq!(counts.removed, 2);
        assert_eq!(counts.loaded(), 5);
    }

    /// Every row is a capture from the reference: `(m, ticks before `loaded.`, tick after)`.
    #[test]
    fn tick_schedule_matches_the_reference_sweep() {
        let cases = [
            (1, 0, true),
            (2, 0, true),
            (63, 0, true),
            (64, 1, false),
            (65, 1, true),
            (100, 1, true),
            (128, 2, false),
            (320, 5, false),
            (321, 5, true),
            (640, 10, false),
            (960, 15, false),
            (1280, 10, false),
            (1600, 12, true),
            (1920, 15, false),
            (2240, 11, true),
            (2560, 13, true),
            (4150, 12, true), // the `sexchr` fixture
            (5000, 15, true),
            (10000, 15, true),
            (20000, 15, true),
            (50000, 15, true), // `bigish`
        ];
        for (m, before, after) in cases {
            let t = Ticks::plan(m);
            assert_eq!((t.before, t.after), (before, after), "m = {m}");
        }
    }

    #[test]
    fn tick_percentages_use_a_fixed_denominator_of_16() {
        // A 25-word load emits thirteen ticks, and they still run 0 6 13 … 75 — not the
        // 0 8 15 … 92 a denominator of thirteen would give.
        let t = Ticks::plan(1600);
        let pct: Vec<String> = (0..=t.before).map(Ticks::tick).collect();
        let want = [
            "0%", "6%", "13%", "19%", "25%", "31%", "38%", "44%", "50%", "56%", "63%", "69%", "75%",
        ];
        assert_eq!(pct.len(), want.len());
        for (got, w) in pct.iter().zip(want) {
            assert_eq!(got, &format!("{w}\r"));
        }
    }

    #[test]
    fn truncation_replays_only_completed_blocks() {
        // Captured: 5000 markers, the first 2000 relabelled to chr23, `.bed` cut after
        // 3000 rows. 1000 retained SNPs had been packed, the block is 192, so five ticks
        // print and the fatal names marker 3000.
        let t = Ticks::plan(3000);
        assert_eq!(t.snps_per_block, 192);
        assert_eq!(1000 / t.snps_per_block, 5);
        assert_eq!(
            console::not_enough_genotypes(3000),
            "Not enough genotypes at the 3000th marker\n"
        );
    }

    #[test]
    fn risk_with_a_model_skips_the_bed_probe() {
        let mut opts = Options::new();
        assert!(probes_bed(&opts));
        opts = parsed(&["--risk", "--model", "m.txt"]);
        assert!(!probes_bed(&opts));
        // Adding another analysis does not bring the probe back.
        opts = parsed(&["--risk", "--model", "m.txt", "--kinship"]);
        assert!(!probes_bed(&opts));
        // `--risk` without a model never reaches the loader, but is not the skip case.
        opts = parsed(&["--risk"]);
        assert!(probes_bed(&opts));
    }

    fn parsed(args: &[&str]) -> Options {
        let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        crate::cli::parse(&owned).options
    }

    #[test]
    fn overrides_are_used_verbatim() {
        let opts = parsed(&["-b", "sub/study.v2.bed"]);
        let (bed, bim, fam) = paths(&opts);
        assert_eq!(bed, PathBuf::from("sub/study.v2.bed"));
        assert_eq!(bim, PathBuf::from("sub/study.v2.bim"));
        assert_eq!(fam, PathBuf::from("sub/study.v2.fam"));

        // `--fam alt` is a file called `alt`, not `alt.fam`.
        let opts = parsed(&["-b", "t.bed", "--fam", "alt", "--bim", "other.bim"]);
        let (_, bim, fam) = paths(&opts);
        assert_eq!(fam, PathBuf::from("alt"));
        assert_eq!(bim, PathBuf::from("other.bim"));
    }

    #[test]
    fn counts_line_degenerates_to_the_simple_form() {
        // The general renderer must agree with the single-clause helper.
        for n in [0, 1, 5000] {
            assert_eq!(
                console::genotype_data_counts(n, 0, 0, 0, 0),
                console::autosome_snps(n)
            );
        }
    }
}
