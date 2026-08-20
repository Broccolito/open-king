//! Byte-for-byte parity for the fileset-loading console sequence.
//!
//! Every expectation below is a transcript of the reference binary's own output for a
//! fileset this test writes from scratch, so the suite pins the sequence without needing
//! KING 2.3.2 on the machine. The filesets are tiny and synthetic; the interesting ones
//! were chosen to exercise a specific rule:
//!
//! * [`simple`] — the tick interleave, whose two halves straddle
//!   `PLINK binary genotypes loaded.`
//! * [`mixed_chromosomes`] / [`sexchr_shifts_the_partition`] — the chromosome partition
//!   and the shape of the map-composition line
//! * the `fails_*` tests — the check order, which is observable because each failure
//!   stops at a different point in the sequence
//!
//! Output is compared from `Loading genotype data …` onwards: the banner and parameter
//! block above it are `console.rs`'s business and carry the fileset's absolute path.

use std::path::{Path, PathBuf};
use std::process::Command;

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Genotype codes as PLINK packs them, two bits per sample.
const HOM1: u8 = 0b00;
const MISSING: u8 = 0b01;
const HET: u8 = 0b10;
const HOM2: u8 = 0b11;

/// A scratch directory that cleans itself up.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let dir = std::env::temp_dir().join(format!("open-king-load-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch dir");
        Scratch(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Write a `.bed`/`.bim`/`.fam` triple.
///
/// `chroms` gives one chromosome label per variant, so the map's composition is spelled
/// out at the call site. Genotypes cycle through a fixed pattern — the loader never looks
/// at them, only at how many there are.
fn write_fileset(scratch: &Scratch, stem: &str, families: &[(&str, &str)], chroms: &[&str]) {
    let fam: String = families
        .iter()
        .map(|(fid, iid)| format!("{fid} {iid} 0 0 1 1\n"))
        .collect();
    std::fs::write(scratch.path(&format!("{stem}.fam")), fam).expect("write fam");

    let bim: String = chroms
        .iter()
        .enumerate()
        .map(|(i, c)| format!("{c} rs{i} 0 {} A G\n", i + 1))
        .collect();
    std::fs::write(scratch.path(&format!("{stem}.bim")), bim).expect("write bim");

    let pattern = [HOM1, HET, HOM2, HOM1, MISSING, HET];
    let mut bed = vec![0x6c, 0x1b, 0x01];
    for v in 0..chroms.len() {
        let mut row = vec![0u8; families.len().div_ceil(4)];
        for s in 0..families.len() {
            let code = pattern[(v + s) % pattern.len()];
            row[s / 4] |= code << (2 * (s % 4));
        }
        bed.extend_from_slice(&row);
    }
    std::fs::write(scratch.path(&format!("{stem}.bed")), bed).expect("write bed");
}

/// Six samples in two families — enough for `.kin` and `.kin0` both to be reachable.
fn two_families() -> Vec<(&'static str, &'static str)> {
    vec![
        ("F1", "S1"),
        ("F1", "S2"),
        ("F1", "S3"),
        ("F2", "S4"),
        ("F2", "S5"),
        ("F2", "S6"),
    ]
}

/// `n` copies of one chromosome label.
fn chrom(label: &'static str, n: usize) -> Vec<&'static str> {
    vec![label; n]
}

// ---------------------------------------------------------------------------
// Harness
// ---------------------------------------------------------------------------

/// Run the binary and return stdout from `Loading genotype data …` onwards, with the
/// scratch directory replaced by `<DIR>` and the timestamp lines blanked.
fn load_output(scratch: &Scratch, args: &[&Path]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_open-king"))
        .args(args)
        // Every analysis writes its output files relative to the working directory. The
        // fileset paths passed in are absolute, so running inside the scratch directory
        // costs nothing and keeps `king.kin` and friends out of the source tree — the
        // scratch directory is removed on drop.
        .current_dir(&scratch.0)
        .output()
        .expect("open-king binary runs");
    assert!(
        out.stderr.is_empty(),
        "the reference writes nothing to stderr"
    );
    let stdout = String::from_utf8(out.stdout).expect("stdout is UTF-8");

    // Everything from the load onwards. A run that dies before the loader speaks — a
    // missing `.bed`, say — has only the FATAL block left after the `KING starts at` line.
    let start = stdout
        .find("Loading genotype data")
        .or_else(|| stdout.find("\nFATAL ERROR"))
        .unwrap_or_else(|| panic!("no load section in:\n{stdout}"));
    let body = stdout[start..].replace(&scratch.0.to_string_lossy().to_string(), "<DIR>");
    // The expectations below spell the scratch paths `<DIR>/t.fam`. On Windows the path
    // the binary echoes is joined with a backslash, so normalise the separator in the
    // part we substituted. This is a property of the fixture, not of the program: the
    // reference prints whatever separator its platform uses, and so do we.
    if std::path::MAIN_SEPARATOR == '\\' {
        body.replace("<DIR>\\", "<DIR>/")
    } else {
        body
    }
}

/// The loading sequence must *begin* with `expected`; what an analysis prints after it is
/// the dispatcher's business, so this stays green as the engines land.
#[track_caller]
fn assert_loads(scratch: &Scratch, args: &[&Path], expected: &str) {
    let got = load_output(scratch, args);
    assert!(
        got.starts_with(expected),
        "load sequence diverged.\n--- want prefix ---\n{expected}\n--- got ---\n{got}"
    );
}

/// A run that ends in a FATAL prints nothing after it, so compare the whole tail.
#[track_caller]
fn assert_fatal(scratch: &Scratch, args: &[&Path], expected: &str) {
    let got = load_output(scratch, args);
    assert_eq!(got, expected, "fatal sequence diverged");
}

fn args<'a>(items: &'a [&'a Path]) -> Vec<&'a Path> {
    items.to_vec()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// 100 autosomal SNPs: two words, one 64-SNP block, and a 36-SNP remainder.
///
/// The remainder is what pushes the second tick past `PLINK binary genotypes loaded.` —
/// the reference prints `0%`, then `loaded.`, then `6%`, then `converted.`, and that
/// ordering is the single easiest thing to get wrong in this whole sequence.
#[test]
fn simple() {
    let s = Scratch::new("simple");
    write_fileset(&s, "t", &two_families(), &chrom("1", 100));
    assert_loads(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 100 autosome SNPs\n",
            "  PLINK maps loaded: 100 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            "0%\r  PLINK binary genotypes loaded.\n",
            "6%\r  KING format genotype data successfully converted.\n",
            "Autosome genotypes stored in 2 words for each of 6 individuals.\n\n",
        ),
    );
}

/// A whole number of blocks emits every tick *before* `loaded.` and none after it.
#[test]
fn exact_block_emits_no_trailing_tick() {
    let s = Scratch::new("exact");
    write_fileset(&s, "t", &two_families(), &chrom("1", 128));
    assert_loads(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 128 autosome SNPs\n",
            "  PLINK maps loaded: 128 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            "0%\r6%\r  PLINK binary genotypes loaded.\n",
            "  KING format genotype data successfully converted.\n",
            "Autosome genotypes stored in 2 words for each of 6 individuals.\n\n",
        ),
    );
}

/// Fewer SNPs than one word: no tick at all before `loaded.`, and a lone `0%` after it.
#[test]
fn sub_word_map_emits_its_only_tick_after_loaded() {
    let s = Scratch::new("subword");
    write_fileset(&s, "t", &two_families(), &chrom("1", 40));
    assert_loads(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 40 autosome SNPs\n",
            "  PLINK maps loaded: 40 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            "  PLINK binary genotypes loaded.\n",
            "0%\r  KING format genotype data successfully converted.\n",
            "Autosome genotypes stored in 1 words for each of 6 individuals.\n\n",
        ),
    );
}

/// The map-composition line in full, and `maps loaded` counting only classified variants.
///
/// XY is folded into the autosome total *and* repeated in the parenthetical, `chr1` is an
/// unrecognised contig rather than chromosome 1, and the two unplaced rows plus `chr1` make
/// up the removed count — so 20 `.bim` rows report as 17 SNPs loaded.
#[test]
fn mixed_chromosomes() {
    let s = Scratch::new("mixed");
    let mut chroms = chrom("1", 5);
    chroms.extend(chrom("22", 3));
    chroms.extend(chrom("25", 2)); // XY
    chroms.extend(chrom("23", 4)); // X
    chroms.extend(chrom("24", 2)); // Y
    chroms.extend(chrom("26", 1)); // MT
    chroms.extend(chrom("0", 2)); // unplaced
    chroms.extend(chrom("chr1", 1)); // a contig name, not chromosome 1
    write_fileset(&s, "t", &two_families(), &chroms);
    assert_loads(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 10 autosome SNPs (including 2 XY SNPs), ",
            "4 X-chromosome SNPs, 2 Y-chromosome SNPs, 1 mitochondrial SNPs\n",
            "  3 other SNPs are removed.\n",
            "  PLINK maps loaded: 17 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
        ),
    );
}

/// `--sexchr 5` moves the whole partition: autosomes become 1–4, X is 5, XY is 7.
///
/// The rows labelled `23`/`25` — X and XY at the default — are now outside every class and
/// are removed, while the rows labelled `X` follow the class rather than the number.
#[test]
fn sexchr_shifts_the_partition() {
    let s = Scratch::new("sexchr");
    let mut chroms = chrom("1", 4);
    chroms.extend(chrom("4", 2));
    chroms.extend(chrom("7", 3)); // XY under --sexchr 5
    chroms.extend(chrom("X", 5)); // binds to the class, i.e. chromosome 5
    chroms.extend(chrom("23", 6)); // no longer X: removed
    chroms.extend(chrom("25", 1)); // no longer XY: removed
    write_fileset(&s, "t", &two_families(), &chroms);
    assert_loads(
        &s,
        &args(&[
            Path::new("-b"),
            &s.path("t.bed"),
            Path::new("--kinship"),
            Path::new("--sexchr"),
            Path::new("5"),
        ]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 9 autosome SNPs (including 3 XY SNPs), ",
            "5 X-chromosome SNPs\n",
            "  7 other SNPs are removed.\n",
            "  PLINK maps loaded: 14 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
        ),
    );
}

/// `--fam` and `--bim` are used verbatim, and the `.bed` still comes from `-b`.
#[test]
fn fam_and_bim_overrides() {
    let s = Scratch::new("override");
    write_fileset(&s, "t", &two_families(), &chrom("1", 64));
    // An alternate map that moves half the variants onto X, so the override is visible in
    // the counts rather than only in the echoed path.
    let alt: String = (0..64)
        .map(|i| {
            format!(
                "{} rs{i} 0 {} A G\n",
                if i < 32 { "23" } else { "1" },
                i + 1
            )
        })
        .collect();
    std::fs::write(s.path("alt.bim"), alt).expect("write alt bim");
    assert_loads(
        &s,
        &args(&[
            Path::new("-b"),
            &s.path("t.bed"),
            Path::new("--kinship"),
            Path::new("--bim"),
            &s.path("alt.bim"),
        ]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/alt.bim...\n",
            "  Genotype data consist of 32 autosome SNPs, 32 X-chromosome SNPs\n",
            "  PLINK maps loaded: 64 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            // Only 32 SNPs survive the override, so this is a sub-word map: no tick
            // before `loaded.` and a lone `0%` after it, exactly as in
            // `sub_word_map_emits_its_only_tick_after_loaded`. The tick schedule follows
            // the *retained* count, not the `.bim`'s 64 rows.
            "  PLINK binary genotypes loaded.\n",
            "0%\r  KING format genotype data successfully converted.\n",
            "Autosome genotypes stored in 1 words for each of 6 individuals.\n\n",
        ),
    );
}

/// A missing `.bed` stops before the loader announces itself.
#[test]
fn fails_before_loading_when_the_bed_is_missing() {
    let s = Scratch::new("nobed");
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("gone.bed"), Path::new("--kinship")]),
        "\nFATAL ERROR - \nGenotype file <DIR>/gone.bed cannot be opened\n\n",
    );
}

/// The parser/loader preserves KING's `--risk --model` probe quirk (unit-tested in
/// `load.rs`), but the executable's product-scope gate now rejects this excluded analysis
/// before that reference path can open any file.
#[test]
fn excluded_risk_is_rejected_before_the_loader_probe_quirk() {
    let s = Scratch::new("risk");
    assert_fatal(
        &s,
        &args(&[
            Path::new("-b"),
            &s.path("gone.bed"),
            Path::new("--risk"),
            Path::new("--model"),
            &s.path("m.txt"),
        ]),
        concat!(
            "\nFATAL ERROR - \n",
            "open-king's minimal relatedness product does not implement: --risk, --model.\n",
            "Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, ",
            "--unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.\n",
            "See docs/SCOPE.md for the product-scope contract.\n\n",
        ),
    );
}

/// A `.bim` with no autosomes still announces the whole map before it gives up.
#[test]
fn fails_after_announcing_a_map_with_no_autosomes() {
    let s = Scratch::new("noauto");
    let mut chroms = chrom("23", 10);
    chroms.extend(chrom("24", 4));
    write_fileset(&s, "t", &two_families(), &chroms);
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 0 autosome SNPs, 10 X-chromosome SNPs, ",
            "4 Y-chromosome SNPs\n",
            "  PLINK maps loaded: 14 SNPs\n",
            "\nFATAL ERROR - \nNo autosome SNPs are available. Please check your map file.\n\n",
        ),
    );
}

/// A truncated `.bed` replays the ticks for the blocks it finished, then names the marker.
///
/// 200 SNPs is four words, so one 64-SNP block; cutting the file after 100 rows leaves one
/// completed block and hence one tick. The marker index is 0-based and the suffix is a
/// literal `th`.
#[test]
fn truncated_bed_replays_completed_blocks() {
    let s = Scratch::new("trunc");
    write_fileset(&s, "t", &two_families(), &chrom("1", 200));
    let bed = std::fs::read(s.path("t.bed")).expect("read bed");
    // 6 samples is 2 bytes per variant.
    std::fs::write(s.path("t.bed"), &bed[..3 + 100 * 2]).expect("truncate bed");
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 200 autosome SNPs\n",
            "  PLINK maps loaded: 200 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            "0%\r",
            "\nFATAL ERROR - \nNot enough genotypes at the 100th marker\n\n\n",
        ),
    );
}

/// An over-long `.bed` is read happily and its trailing bytes ignored.
#[test]
fn over_long_bed_loads() {
    let s = Scratch::new("long");
    write_fileset(&s, "t", &two_families(), &chrom("1", 64));
    let mut bed = std::fs::read(s.path("t.bed")).expect("read bed");
    bed.extend_from_slice(b"trailing garbage");
    std::fs::write(s.path("t.bed"), bed).expect("extend bed");
    assert_loads(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 64 autosome SNPs\n",
            "  PLINK maps loaded: 64 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            "0%\r  PLINK binary genotypes loaded.\n",
            "  KING format genotype data successfully converted.\n",
            "Autosome genotypes stored in 1 words for each of 6 individuals.\n\n",
        ),
    );
}

/// An individual-major `.bed` gets all the way to `Read in PLINK bed file …` first.
///
/// Contrast [`fails_on_bad_magic`], which stops before the loader speaks at all: the magic
/// and the mode byte are checked at different points, and only the magic is checked early.
#[test]
fn individual_major_bed_fails_late() {
    let s = Scratch::new("indiv");
    write_fileset(&s, "t", &two_families(), &chrom("1", 64));
    let mut bed = std::fs::read(s.path("t.bed")).expect("read bed");
    bed[2] = 0x00;
    std::fs::write(s.path("t.bed"), bed).expect("rewrite bed");
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "  PLINK pedigrees loaded: 6 samples\n",
            "Read in PLINK bim file <DIR>/t.bim...\n",
            "  Genotype data consist of 64 autosome SNPs\n",
            "  PLINK maps loaded: 64 SNPs\n",
            "Read in PLINK bed file <DIR>/t.bed...\n",
            "\nFATAL ERROR - \nCurrently only SNP-major mode can be analyzed.\n\n",
        ),
    );
}

#[test]
fn fails_on_bad_magic() {
    let s = Scratch::new("magic");
    write_fileset(&s, "t", &two_families(), &chrom("1", 64));
    let mut bed = std::fs::read(s.path("t.bed")).expect("read bed");
    bed[0] = b'X';
    std::fs::write(s.path("t.bed"), bed).expect("rewrite bed");
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        "\nFATAL ERROR - \nPlease use either PLINK or KING binary format as input.\n\n",
    );
}

/// A readable file whose name does not end in `.bed` is rejected on the name alone.
#[test]
fn fails_on_a_name_that_is_not_a_bed() {
    let s = Scratch::new("suffix");
    write_fileset(&s, "t", &two_families(), &chrom("1", 64));
    std::fs::copy(s.path("t.bed"), s.path("plain")).expect("copy bed");
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("plain"), Path::new("--kinship")]),
        "\nFATAL ERROR - \nPlease use PLINK binary format as input.\n\n",
    );
}

/// A duplicated `(FID, IID)` is named, and `PLINK pedigrees loaded` never appears.
#[test]
fn duplicate_sample_is_named_before_the_fatal() {
    let s = Scratch::new("dup");
    let families = vec![("F1", "S1"), ("F1", "S1"), ("F2", "S2")];
    write_fileset(&s, "t", &families, &chrom("1", 64));
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "Family F1: Person S1 is duplicated\n",
            "\nFATAL ERROR - \nPlease correct problems with pedigree structure\n\n\n",
        ),
    );
}

#[test]
fn case_colliding_sample_is_named_before_the_fatal() {
    let s = Scratch::new("case-duplicate");
    write_fileset(&s, "t", &[("F1", "A_F"), ("f1", "a_f")], &chrom("1", 64));
    assert_fatal(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--kinship")]),
        concat!(
            "Loading genotype data in PLINK binary format...\n",
            "Read in PLINK fam file <DIR>/t.fam...\n",
            "Family f1: Person a_f is duplicated\n",
            "\nFATAL ERROR - \nPlease correct problems with pedigree structure\n\n\n",
        ),
    );
}

/// Analyses that do not open with the preamble must not be given one.
#[test]
fn build_does_not_print_the_preamble() {
    let s = Scratch::new("build");
    write_fileset(&s, "t", &two_families(), &chrom("1", 64));
    let got = load_output(
        &s,
        &args(&[Path::new("-b"), &s.path("t.bed"), Path::new("--build")]),
    );
    let tail = got
        .split_once("  KING format genotype data successfully converted.\n")
        .expect("the load ends with the conversion line")
        .1;
    assert!(
        tail.starts_with("\nOptions in effect:\n\t--build\n"),
        "--build goes straight from the conversion line to its own block, with no \
         `Autosome genotypes stored in …` preamble between them:\n{got}"
    );
}
