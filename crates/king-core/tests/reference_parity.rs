//! End-to-end numeric parity against the reference KING 2.3.2 binary.
//!
//! Each test builds a PLINK 1 fileset from scratch (the `.bed`/`.bim`/`.fam` bytes are
//! written here, so nothing depends on PLINK being installed), runs the reference binary
//! over it with `--ibs --kinship`, then runs `king-io` + `king-core` over the same files
//! and compares every column of all four output tables:
//!
//! * counts (`N_SNP`, `N_IBS0`, `N_IBS1`, `N_IBS2`, `NHetHet`, `NHomHom`, `N_Het1`,
//!   `N_Het2`) must be **exactly** equal;
//! * everything printed with `%.4f` or `%.3f` must agree **as printed**, i.e. the string
//!   our value formats to must equal the string the reference wrote.
//!
//! Gated on the `KING_REF` environment variable so a machine without the reference binary
//! still passes:
//!
//! ```text
//! KING_REF=/path/to/king cargo test -p king-core --test reference_parity -- --nocapture
//! ```
//!
//! # The SNP-inclusion rule these tests encode
//!
//! Established black-box (see `snp_inclusion` below, which re-derives it from the binary
//! on every run): the reference keeps a `.bim` row for the default relatedness analysis
//! **iff its chromosome field maps to PLINK code 1..=22 or 25 (XY)**. Nothing else is
//! filtered — not monomorphic markers, not markers with an allele written `0`, not
//! markers at any missingness rate, and there is no MAF cutoff. Missing calls are handled
//! purely pairwise, at counting time.

#![forbid(unsafe_code)]

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

use king_core::infer::{KinshipCache, Pedigree};
use king_core::{counts, kinship, PairCounts, Scope};
use king_io::{bed, Fileset, VariantFilter};

// ---------------------------------------------------------------------------
// Reference binary discovery
// ---------------------------------------------------------------------------

/// Path to the reference binary, or `None` when `KING_REF` is unset.
///
/// A set-but-wrong `KING_REF` is a hard error: silently skipping there would turn a
/// typo into a green test run.
fn reference_binary() -> Option<PathBuf> {
    let raw = std::env::var_os("KING_REF")?;
    if raw.is_empty() {
        return None;
    }
    let path = PathBuf::from(raw);
    assert!(
        path.is_file(),
        "KING_REF points at {} which is not a file",
        path.display()
    );
    Some(path)
}

/// Bind the reference binary or return from the test with a skip notice.
macro_rules! reference_or_skip {
    ($name:expr) => {
        match reference_binary() {
            Some(p) => p,
            None => {
                eprintln!(
                    "SKIP {}: set KING_REF=/path/to/king (KING 2.3.2) to run reference-parity tests",
                    $name
                );
                return;
            }
        }
    };
}

/// A scratch directory under the crate's integration-test temp dir.
fn workdir(name: &str) -> PathBuf {
    let dir = Path::new(env!("CARGO_TARGET_TMPDIR"))
        .join("parity")
        .join(name);
    if dir.exists() {
        fs::remove_dir_all(&dir).expect("clear workdir");
    }
    fs::create_dir_all(&dir).expect("create workdir");
    dir
}

// ---------------------------------------------------------------------------
// Deterministic RNG
// ---------------------------------------------------------------------------

/// xorshift64. Fixed seeds keep every fileset reproducible run to run.
struct Rng(u64);

impl Rng {
    fn new(seed: u64) -> Self {
        Self(seed | 1)
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.0;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.0 = x;
        x
    }

    /// Uniform in `[0, 1)`.
    fn unit(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    fn chance(&mut self, p: f64) -> bool {
        self.unit() < p
    }

    fn range(&mut self, lo: f64, hi: f64) -> f64 {
        lo + (hi - lo) * self.unit()
    }
}

// ---------------------------------------------------------------------------
// PLINK 1 fileset writer
// ---------------------------------------------------------------------------

/// Dosage of allele 1: 2 = hom A1/A1, 1 = het, 0 = hom A2/A2, `None` = no call.
type Geno = Option<u8>;

/// The two-bit code PLINK writes for a genotype, low-order pair = first sample in a byte.
fn plink_code(g: Geno) -> u8 {
    match g {
        Some(2) => 0b00,
        None => 0b01,
        Some(1) => 0b10,
        Some(0) => 0b11,
        Some(other) => panic!("dosage {other} is not 0, 1 or 2"),
    }
}

#[derive(Clone)]
struct Person {
    fid: String,
    iid: String,
    pat: String,
    mat: String,
    sex: u8,
}

impl Person {
    fn founder(fid: &str, iid: &str, sex: u8) -> Self {
        Self {
            fid: fid.into(),
            iid: iid.into(),
            pat: "0".into(),
            mat: "0".into(),
            sex,
        }
    }

    fn child(fid: &str, iid: &str, sex: u8, pat: &str, mat: &str) -> Self {
        Self {
            fid: fid.into(),
            iid: iid.into(),
            pat: pat.into(),
            mat: mat.into(),
            sex,
        }
    }
}

#[derive(Clone)]
struct Marker {
    chrom: String,
    id: String,
    bp: u32,
    a1: String,
    a2: String,
}

impl Marker {
    fn new(chrom: &str, id: String, bp: u32) -> Self {
        Self {
            chrom: chrom.into(),
            id,
            bp,
            a1: "A".into(),
            a2: "G".into(),
        }
    }
}

/// A complete fileset: people in `.fam` order, markers in `.bim` order, and
/// `geno[marker][person]`.
struct Dataset {
    name: &'static str,
    people: Vec<Person>,
    markers: Vec<Marker>,
    geno: Vec<Vec<Geno>>,
}

impl Dataset {
    /// Write `.bed`/`.bim`/`.fam` into `dir` and return the prefix.
    fn write(&self, dir: &Path) -> PathBuf {
        let prefix = dir.join(self.name);
        let n_samples = self.people.len();

        let mut fam = String::new();
        for p in &self.people {
            fam.push_str(&format!(
                "{} {} {} {} {} 1\n",
                p.fid, p.iid, p.pat, p.mat, p.sex
            ));
        }
        fs::write(prefix.with_extension("fam"), fam).expect("write .fam");

        let mut bim = String::new();
        for m in &self.markers {
            bim.push_str(&format!(
                "{}\t{}\t0\t{}\t{}\t{}\n",
                m.chrom, m.id, m.bp, m.a1, m.a2
            ));
        }
        fs::write(prefix.with_extension("bim"), bim).expect("write .bim");

        let mut bed = vec![0x6c, 0x1b, 0x01];
        for row in &self.geno {
            assert_eq!(row.len(), n_samples, "genotype row has the wrong width");
            for chunk in row.chunks(4) {
                let mut byte = 0u8;
                for (k, &g) in chunk.iter().enumerate() {
                    byte |= plink_code(g) << (2 * k);
                }
                bed.push(byte);
            }
        }
        fs::write(prefix.with_extension("bed"), bed).expect("write .bed");

        prefix
    }

    /// Every unordered pair, split by whether the two share an `FID`.
    fn pair_split(&self) -> (usize, usize) {
        let (mut within, mut between) = (0, 0);
        for i in 0..self.people.len() {
            for j in i + 1..self.people.len() {
                if self.people[i].fid == self.people[j].fid {
                    within += 1;
                } else {
                    between += 1;
                }
            }
        }
        (within, between)
    }
}

// ---------------------------------------------------------------------------
// Genotype simulation
// ---------------------------------------------------------------------------

/// Where a person's alleles come from. Parents must precede their children.
#[derive(Clone, Copy)]
enum Origin {
    Founder,
    Child {
        pat: usize,
        mat: usize,
    },
    /// A byte-exact copy of another sample, missing calls included.
    Duplicate(usize),
}

/// Draw unlinked biallelic markers down a pedigree, then knock out calls per sample.
///
/// `missing[s]` is the per-call dropout probability for sample `s`; a `Duplicate` copies
/// its source's final genotype, so the copy's dropout is identical rather than independent
/// (that is what makes it an *exact* duplicate).
fn simulate(
    rng: &mut Rng,
    origins: &[Origin],
    n_markers: usize,
    missing: &[f64],
) -> Vec<Vec<Geno>> {
    let n = origins.len();
    assert_eq!(missing.len(), n);
    let mut out = Vec::with_capacity(n_markers);

    for _ in 0..n_markers {
        let freq = rng.range(0.10, 0.90);
        let mut hap: Vec<(u8, u8)> = vec![(0, 0); n];
        for s in 0..n {
            hap[s] = match origins[s] {
                Origin::Founder => (u8::from(rng.chance(freq)), u8::from(rng.chance(freq))),
                Origin::Child { pat, mat } => {
                    let from_pat = if rng.chance(0.5) {
                        hap[pat].0
                    } else {
                        hap[pat].1
                    };
                    let from_mat = if rng.chance(0.5) {
                        hap[mat].0
                    } else {
                        hap[mat].1
                    };
                    (from_pat, from_mat)
                }
                Origin::Duplicate(src) => hap[src],
            };
        }

        let mut row: Vec<Geno> = hap.iter().map(|&(a, b)| Some(a + b)).collect();
        for s in 0..n {
            match origins[s] {
                Origin::Duplicate(src) => row[s] = row[src],
                _ => {
                    if rng.chance(missing[s]) {
                        row[s] = None;
                    }
                }
            }
        }
        out.push(row);
    }
    out
}

/// `n` autosomal markers on chromosome 1, evenly spaced.
fn autosomal_markers(n: usize) -> Vec<Marker> {
    (0..n)
        .map(|k| Marker::new("1", format!("rs{k}"), (k as u32 + 1) * 1000))
        .collect()
}

// ---------------------------------------------------------------------------
// Reference binary invocation and output parsing
// ---------------------------------------------------------------------------

/// Run the reference binary and return its combined stdout + stderr.
fn run_reference(king: &Path, prefix: &Path, outdir: &Path, args: &[&str]) -> String {
    fs::create_dir_all(outdir).expect("create output dir");
    let out = Command::new(king)
        .arg("-b")
        .arg(prefix.with_extension("bed"))
        .arg("--prefix")
        .arg(outdir.join("king"))
        .args(args)
        .output()
        .unwrap_or_else(|e| panic!("cannot run {}: {e}", king.display()));
    let text =
        String::from_utf8_lossy(&out.stdout).into_owned() + &String::from_utf8_lossy(&out.stderr);
    assert!(
        !text.contains("FATAL ERROR"),
        "reference binary failed on {}:\n{text}",
        prefix.display()
    );
    text
}

/// One row of a reference output table, keyed by header name.
type Row = BTreeMap<String, String>;

/// Parse a tab-separated reference table. `None` when the file was not created; an empty
/// vector when it is zero-byte or header-only (both of which the reference does emit).
fn read_table(path: &Path) -> Option<Vec<Row>> {
    let text = fs::read_to_string(path).ok()?;
    let mut lines = text.lines().filter(|l| !l.trim().is_empty());
    let Some(header) = lines.next() else {
        return Some(Vec::new());
    };
    let header: Vec<&str> = header.split('\t').collect();
    let rows = lines
        .map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                fields.len(),
                header.len(),
                "{}: row has {} fields, header has {}",
                path.display(),
                fields.len(),
                header.len()
            );
            header
                .iter()
                .zip(fields)
                .map(|(h, v)| ((*h).to_string(), v.to_string()))
                .collect()
        })
        .collect();
    Some(rows)
}

/// Fetch a column, failing loudly rather than silently skipping a comparison.
fn col<'a>(row: &'a Row, name: &str) -> &'a str {
    row.get(name)
        .unwrap_or_else(|| panic!("reference row has no column {name:?}: {row:?}"))
}

/// The console's `Genotype data consist of N autosome SNPs...` figure.
fn autosome_snps_reported(console: &str) -> usize {
    let tail = console
        .split("Genotype data consist of ")
        .nth(1)
        .unwrap_or_else(|| panic!("console has no SNP-count line:\n{console}"));
    tail.split_whitespace()
        .next()
        .and_then(|n| n.parse().ok())
        .unwrap_or_else(|| panic!("cannot parse the SNP count from {tail:?}"))
}

// ---------------------------------------------------------------------------
// Formatting, matched to the reference's printf
// ---------------------------------------------------------------------------

/// `%.*f` as C prints it: non-finite values are `nan` / `inf` / `-inf`, not Rust's `NaN`.
fn fmt_c(x: f64, places: usize) -> String {
    if x.is_nan() {
        "nan".to_string()
    } else if x.is_infinite() {
        if x > 0.0 {
            "inf".to_string()
        } else {
            "-inf".to_string()
        }
    } else {
        format!("{x:.places$}")
    }
}

// ---------------------------------------------------------------------------
// Comparison bookkeeping
// ---------------------------------------------------------------------------

/// Per-column pass/fail tally for one dataset, so a failure names the column and the pair.
#[derive(Default)]
struct Report {
    order: Vec<String>,
    checked: BTreeMap<String, usize>,
    failed: BTreeMap<String, usize>,
    detail: Vec<String>,
}

impl Report {
    fn record(&mut self, column: &str, ok: bool, context: &str, ours: &str, theirs: &str) {
        if !self.checked.contains_key(column) {
            self.order.push(column.to_string());
            self.checked.insert(column.to_string(), 0);
            self.failed.insert(column.to_string(), 0);
        }
        *self.checked.get_mut(column).unwrap() += 1;
        if !ok {
            *self.failed.get_mut(column).unwrap() += 1;
            if self.detail.len() < 12 {
                self.detail.push(format!(
                    "  {context} {column}: ours={ours} reference={theirs}"
                ));
            }
        }
    }

    fn int(&mut self, column: &str, context: &str, ours: u32, theirs: &str) {
        let want: u32 = theirs
            .parse()
            .unwrap_or_else(|_| panic!("{context}: {column} {theirs:?} is not an integer"));
        self.record(column, ours == want, context, &ours.to_string(), theirs);
    }

    fn float(&mut self, column: &str, context: &str, ours: f64, theirs: &str, places: usize) {
        let printed = fmt_c(ours, places);
        self.record(column, printed == theirs, context, &printed, theirs);
    }

    /// Print the per-column tally and panic if anything mismatched.
    fn finish(self, dataset: &str, pairs: usize) {
        println!("\n[{dataset}] {pairs} pairs compared");
        let width = self.order.iter().map(String::len).max().unwrap_or(0).max(8);
        for column in &self.order {
            let n = self.checked[column];
            let bad = self.failed[column];
            let verdict = if bad == 0 {
                "PASS".to_string()
            } else {
                format!("FAIL ({bad} mismatched)")
            };
            println!("  {column:<width$}  {n:>6} checks  {verdict}");
        }
        let total_bad: usize = self.failed.values().sum();
        assert!(
            total_bad == 0,
            "[{dataset}] {total_bad} column values differ from the reference:\n{}",
            self.detail.join("\n")
        );
    }
}

// ---------------------------------------------------------------------------
// Our engine
// ---------------------------------------------------------------------------

struct Engine {
    fileset: Fileset,
    index: BTreeMap<(String, String), usize>,
    pedigree: Pedigree,
}

impl Engine {
    fn load(prefix: &Path) -> Self {
        let fileset = bed::read_fileset(prefix, VariantFilter::Autosomes, None, None)
            .unwrap_or_else(|e| panic!("king-io failed on {}: {e}", prefix.display()));
        let index = fileset
            .samples
            .iter()
            .enumerate()
            .map(|(i, s)| ((s.fid.clone(), s.iid.clone()), i))
            .collect();
        let pedigree = Pedigree::from_samples(&fileset.samples);
        Self {
            fileset,
            index,
            pedigree,
        }
    }

    fn lookup(&self, fid: &str, iid: &str) -> usize {
        *self
            .index
            .get(&(fid.to_string(), iid.to_string()))
            .unwrap_or_else(|| panic!("reference names a sample we did not load: {fid} {iid}"))
    }

    fn counts_for(&self, pairs: &[(usize, usize)]) -> Vec<PairCounts> {
        counts::all_pairs(&self.fileset.genotypes, pairs)
    }
}

// ---------------------------------------------------------------------------
// The comparison itself
// ---------------------------------------------------------------------------

/// Columns of `.ibs` / `.ibs0` that carry raw counts.
fn check_counts(report: &mut Report, context: &str, c: &PairCounts, row: &Row) {
    report.int("N_SNP", context, c.n_snp, col(row, "N_SNP"));
    report.int("N_IBS0", context, c.ibs0, col(row, "N_IBS0"));
    report.int("N_IBS1", context, c.ibs1(), col(row, "N_IBS1"));
    report.int("N_IBS2", context, c.ibs2(), col(row, "N_IBS2"));
    report.int("NHetHet", context, c.het_het, col(row, "NHetHet"));
    report.int("NHomHom", context, c.hom_hom, col(row, "NHomHom"));
    report.int("N_Het1", context, c.het_i, col(row, "N_Het1"));
    report.int("N_Het2", context, c.het_j, col(row, "N_Het2"));
}

/// Columns of `.ibs` / `.ibs0` derived from the counts.
fn check_ibs_floats(report: &mut Report, context: &str, c: &PairCounts, row: &Row, scope: Scope) {
    report.float("IBS", context, kinship::ibs_mean(c), col(row, "IBS"), 4);
    report.float("Dist", context, kinship::dist(c), col(row, "Dist"), 4);
    report.float(
        "HetConc",
        context,
        kinship::het_concordance(c),
        col(row, "HetConc"),
        4,
    );
    report.float(
        "Het2|1",
        context,
        kinship::het_2given1(c),
        col(row, "Het2|1"),
        4,
    );
    report.float(
        "Het1|2",
        context,
        kinship::het_1given2(c),
        col(row, "Het1|2"),
        4,
    );
    report.float(
        "HomConc",
        context,
        kinship::hom_concordance(c),
        col(row, "HomConc"),
        4,
    );
    report.float(
        "Kinship",
        context,
        kinship::kinship(c, scope),
        col(row, "Kinship"),
        4,
    );
}

/// `Z0` / `Phi` are pedigree expectations, not estimates; they appear in `.kin` and `.ibs`.
fn check_pedigree(
    report: &mut Report,
    context: &str,
    engine: &Engine,
    cache: &mut KinshipCache,
    a: usize,
    b: usize,
    row: &Row,
) {
    let phi = king_core::infer::pedigree_kinship(&engine.pedigree, cache, a, b);
    let z0 = king_core::infer::pedigree_z0(&engine.pedigree, cache, a, b);
    report.float("Z0", context, z0, col(row, "Z0"), 3);
    report.float("Phi", context, phi, col(row, "Phi"), 4);
}

/// Build, run, compare. Returns the number of pairs the reference actually emitted.
fn parity(king: &Path, ds: &Dataset) -> usize {
    let dir = workdir(ds.name);
    let prefix = ds.write(&dir);
    let outdir = dir.join("ref");
    let console = run_reference(king, &prefix, &outdir, &["--ibs", "--kinship"]);

    let engine = Engine::load(&prefix);
    let mut report = Report::default();
    let mut cache = KinshipCache::default();

    // The console's own SNP figure is a check on the inclusion rule before any pair math.
    let expected_autosomes = ds
        .markers
        .iter()
        .filter(|m| reference_keeps(&m.chrom))
        .count();
    assert_eq!(
        autosome_snps_reported(&console),
        expected_autosomes,
        "[{}] reference kept a different marker set than the established rule predicts",
        ds.name
    );
    assert_eq!(
        engine.fileset.kept.len(),
        expected_autosomes,
        "[{}] king-io kept {} markers, the reference kept {}",
        ds.name,
        engine.fileset.kept.len(),
        expected_autosomes
    );

    let (want_within, want_between) = ds.pair_split();

    // ---- within-family: .ibs (raw counts) and .kin (proportions) ----
    let ibs = read_table(&outdir.join("king.ibs")).unwrap_or_default();
    let kin = read_table(&outdir.join("king.kin")).unwrap_or_default();
    let within_pairs: Vec<(usize, usize)> = ibs
        .iter()
        .map(|r| {
            let fid = col(r, "FID");
            (
                engine.lookup(fid, col(r, "ID1")),
                engine.lookup(fid, col(r, "ID2")),
            )
        })
        .collect();
    let within_counts = engine.counts_for(&within_pairs);
    for ((row, &(a, b)), c) in ibs.iter().zip(&within_pairs).zip(&within_counts) {
        let context = format!(
            "{}:{}/{}",
            col(row, "FID"),
            col(row, "ID1"),
            col(row, "ID2")
        );
        check_counts(&mut report, &context, c, row);
        check_ibs_floats(&mut report, &context, c, row, Scope::WithinFamily);
        check_pedigree(&mut report, &context, &engine, &mut cache, a, b, row);
    }

    // `.kin` is written zero-byte when the dataset holds a single family; when it is
    // populated it must cover exactly the same pairs as `.ibs`.
    if !kin.is_empty() {
        assert_eq!(kin.len(), ibs.len(), "[{}] .kin and .ibs disagree", ds.name);
        let keyed: BTreeMap<(String, String, String), &Row> = kin
            .iter()
            .map(|r| {
                (
                    (
                        col(r, "FID").to_string(),
                        col(r, "ID1").to_string(),
                        col(r, "ID2").to_string(),
                    ),
                    r,
                )
            })
            .collect();
        for ((ibs_row, &(a, b)), c) in ibs.iter().zip(&within_pairs).zip(&within_counts) {
            let key = (
                col(ibs_row, "FID").to_string(),
                col(ibs_row, "ID1").to_string(),
                col(ibs_row, "ID2").to_string(),
            );
            let row = keyed
                .get(&key)
                .unwrap_or_else(|| panic!("[{}] .kin has no row for {key:?}", ds.name));
            let context = format!("kin {}:{}/{}", key.0, key.1, key.2);
            report.int("N_SNP", &context, c.n_snp, col(row, "N_SNP"));
            report.float(
                "HetHet(prop)",
                &context,
                kinship::het_het_prop(c),
                col(row, "HetHet"),
                4,
            );
            report.float(
                "IBS0(prop)",
                &context,
                kinship::ibs0_prop(c),
                col(row, "IBS0"),
                4,
            );
            report.float(
                "Kinship",
                &context,
                kinship::kinship(c, Scope::WithinFamily),
                col(row, "Kinship"),
                4,
            );
            check_pedigree(&mut report, &context, &engine, &mut cache, a, b, row);
        }
    }

    // ---- between-family: .ibs0 (raw counts) and .kin0 (proportions) ----
    let ibs0 = read_table(&outdir.join("king.ibs0")).unwrap_or_default();
    let kin0 = read_table(&outdir.join("king.kin0")).unwrap_or_default();
    let between_pairs: Vec<(usize, usize)> = ibs0
        .iter()
        .map(|r| {
            (
                engine.lookup(col(r, "FID1"), col(r, "ID1")),
                engine.lookup(col(r, "FID2"), col(r, "ID2")),
            )
        })
        .collect();
    let between_counts = engine.counts_for(&between_pairs);
    for (row, c) in ibs0.iter().zip(&between_counts) {
        let context = format!(
            "{}:{}/{}:{}",
            col(row, "FID1"),
            col(row, "ID1"),
            col(row, "FID2"),
            col(row, "ID2")
        );
        check_counts(&mut report, &context, c, row);
        check_ibs_floats(&mut report, &context, c, row, Scope::BetweenFamily);
    }

    if !kin0.is_empty() {
        let keyed: BTreeMap<(String, String, String, String), &Row> = kin0
            .iter()
            .map(|r| {
                (
                    (
                        col(r, "FID1").to_string(),
                        col(r, "ID1").to_string(),
                        col(r, "FID2").to_string(),
                        col(r, "ID2").to_string(),
                    ),
                    r,
                )
            })
            .collect();
        assert_eq!(
            kin0.len(),
            ibs0.len(),
            "[{}] .kin0 and .ibs0 disagree",
            ds.name
        );
        for (ibs0_row, c) in ibs0.iter().zip(&between_counts) {
            let key = (
                col(ibs0_row, "FID1").to_string(),
                col(ibs0_row, "ID1").to_string(),
                col(ibs0_row, "FID2").to_string(),
                col(ibs0_row, "ID2").to_string(),
            );
            let row = keyed
                .get(&key)
                .unwrap_or_else(|| panic!("[{}] .kin0 has no row for {key:?}", ds.name));
            let context = format!("kin0 {}:{}/{}:{}", key.0, key.1, key.2, key.3);
            report.int("N_SNP", &context, c.n_snp, col(row, "N_SNP"));
            report.float(
                "HetHet(prop)",
                &context,
                kinship::het_het_prop(c),
                col(row, "HetHet"),
                4,
            );
            report.float(
                "IBS0(prop)",
                &context,
                kinship::ibs0_prop(c),
                col(row, "IBS0"),
                4,
            );
            report.float(
                "Kinship",
                &context,
                kinship::kinship(c, Scope::BetweenFamily),
                col(row, "Kinship"),
                4,
            );
        }
    }

    // Coverage: the reference must have emitted every pair we expected it to.
    assert_eq!(
        ibs.len(),
        want_within,
        "[{}] .ibs covered {} of {want_within} within-family pairs",
        ds.name,
        ibs.len()
    );
    assert_eq!(
        ibs0.len(),
        want_between,
        "[{}] .ibs0 covered {} of {want_between} between-family pairs",
        ds.name,
        ibs0.len()
    );

    let pairs = ibs.len() + ibs0.len();
    report.finish(ds.name, pairs);
    pairs
}

/// The established inclusion rule: PLINK chromosome codes 1..=22 and 25 (XY).
///
/// Only the spellings the reference itself recognises map to a code — bare integers
/// (leading zeros and whitespace tolerated) and the bare letters `X`, `Y`, `XY`, `MT`.
/// A `chr` prefix, `M`, `PAR1` and anything else make the reference drop the marker.
fn reference_keeps(chrom: &str) -> bool {
    matches!(reference_chrom_code(chrom), Some(1..=22) | Some(25))
}

/// PLINK code the reference assigns to a `.bim` chromosome field, or `None` when it
/// refuses the marker entirely. Derived black-box; `snp_inclusion_chromosome_set` below
/// re-checks every arm against the binary.
fn reference_chrom_code(chrom: &str) -> Option<u8> {
    let t = chrom.trim();
    // `23.0` parses as 23: the integer prefix is what counts.
    let digits: String = t.chars().take_while(char::is_ascii_digit).collect();
    if !digits.is_empty() {
        let n: u32 = digits.parse().ok()?;
        return if (1..=26).contains(&n) {
            Some(n as u8)
        } else {
            None
        };
    }
    match t.to_ascii_uppercase().as_str() {
        "X" => Some(23),
        "Y" => Some(24),
        "XY" => Some(25),
        "MT" => Some(26),
        _ => None,
    }
}

// ---------------------------------------------------------------------------
// Datasets
// ---------------------------------------------------------------------------

const N_MARKERS: usize = 3000;

/// Two independent trios. Two families are needed for `.kin` to be written at all.
fn dataset_trio() -> Dataset {
    let mut people = Vec::new();
    let mut origins = Vec::new();
    for f in 1..=2 {
        let fid = format!("T{f}");
        let base = people.len();
        people.push(Person::founder(&fid, &format!("T{f}_dad"), 1));
        people.push(Person::founder(&fid, &format!("T{f}_mom"), 2));
        people.push(Person::child(
            &fid,
            &format!("T{f}_kid"),
            1,
            &format!("T{f}_dad"),
            &format!("T{f}_mom"),
        ));
        origins.push(Origin::Founder);
        origins.push(Origin::Founder);
        origins.push(Origin::Child {
            pat: base,
            mat: base + 1,
        });
    }
    let missing = vec![0.0; people.len()];
    let mut rng = Rng::new(0x5EED_0001);
    let geno = simulate(&mut rng, &origins, N_MARKERS, &missing);
    Dataset {
        name: "trio",
        people,
        markers: autosomal_markers(N_MARKERS),
        geno,
    }
}

/// One nuclear family with four full sibs, plus an unrelated two-person family.
fn dataset_nuclear_four_sibs() -> Dataset {
    let mut people = vec![
        Person::founder("N1", "N1_dad", 1),
        Person::founder("N1", "N1_mom", 2),
    ];
    let mut origins = vec![Origin::Founder, Origin::Founder];
    for k in 1..=4 {
        people.push(Person::child(
            "N1",
            &format!("N1_kid{k}"),
            1 + (k % 2) as u8,
            "N1_dad",
            "N1_mom",
        ));
        origins.push(Origin::Child { pat: 0, mat: 1 });
    }
    people.push(Person::founder("N2", "N2_a", 1));
    people.push(Person::founder("N2", "N2_b", 2));
    origins.push(Origin::Founder);
    origins.push(Origin::Founder);

    let missing = vec![0.0; people.len()];
    let mut rng = Rng::new(0x5EED_0002);
    let geno = simulate(&mut rng, &origins, N_MARKERS, &missing);
    Dataset {
        name: "nuclear4",
        people,
        markers: autosomal_markers(N_MARKERS),
        geno,
    }
}

/// Wildly unequal per-sample call rates, so the pairwise-missingness rule is load-bearing
/// on every pair and `N_Het1`/`N_Het2` cannot come from a per-sample cache.
fn dataset_asymmetric_missing() -> Dataset {
    let people = vec![
        Person::founder("A1", "A1_dad", 1),
        Person::founder("A1", "A1_mom", 2),
        Person::child("A1", "A1_kid", 1, "A1_dad", "A1_mom"),
        Person::founder("A2", "A2_a", 1),
        Person::founder("A2", "A2_b", 2),
        Person::founder("A3", "A3_a", 1),
        Person::founder("A3", "A3_b", 2),
    ];
    let origins = vec![
        Origin::Founder,
        Origin::Founder,
        Origin::Child { pat: 0, mat: 1 },
        Origin::Founder,
        Origin::Founder,
        Origin::Founder,
        Origin::Founder,
    ];
    // The trio's child is the 30%-missing sample, so a *related* pair is asymmetric too.
    let missing = vec![0.00, 0.02, 0.30, 0.00, 0.45, 0.15, 0.60];
    let mut rng = Rng::new(0x5EED_0003);
    let geno = simulate(&mut rng, &origins, N_MARKERS, &missing);
    Dataset {
        name: "asym_missing",
        people,
        markers: autosomal_markers(N_MARKERS),
        geno,
    }
}

/// Monomorphic markers of every flavour, `0`-allele markers, and all-missing markers,
/// mixed into a polymorphic background.
fn dataset_monomorphic() -> Dataset {
    let people = vec![
        Person::founder("M1", "M1_dad", 1),
        Person::child("M1", "M1_kid", 2, "M1_dad", "0"),
        Person::founder("M2", "M2_a", 1),
        Person::founder("M2", "M2_b", 2),
        Person::founder("M3", "M3_a", 1),
        Person::founder("M3", "M3_b", 2),
    ];
    let origins = vec![
        Origin::Founder,
        Origin::Child { pat: 0, mat: 0 },
        Origin::Founder,
        Origin::Founder,
        Origin::Founder,
        Origin::Founder,
    ];
    let n = people.len();
    let missing = vec![0.03; n];
    let mut rng = Rng::new(0x5EED_0004);

    // 2000 polymorphic + 1000 degenerate markers.
    let polymorphic = 2000;
    let mut geno = simulate(&mut rng, &origins, polymorphic, &missing);
    let mut markers = autosomal_markers(polymorphic);

    let push = |markers: &mut Vec<Marker>,
                geno: &mut Vec<Vec<Geno>>,
                count: usize,
                tag: &str,
                a1: &str,
                row: Vec<Geno>| {
        for k in 0..count {
            let idx = markers.len();
            let mut m = Marker::new("1", format!("{tag}{k}"), (idx as u32 + 1) * 1000);
            m.a1 = a1.to_string();
            markers.push(m);
            geno.push(row.clone());
        }
    };
    push(&mut markers, &mut geno, 400, "monA2", "A", vec![Some(0); n]);
    push(&mut markers, &mut geno, 200, "monA1", "A", vec![Some(2); n]);
    push(
        &mut markers,
        &mut geno,
        100,
        "allhet",
        "A",
        vec![Some(1); n],
    );
    push(
        &mut markers,
        &mut geno,
        200,
        "zeroa1",
        "0",
        vec![Some(0); n],
    );
    push(&mut markers, &mut geno, 100, "nocall", "A", vec![None; n]);

    Dataset {
        name: "monomorphic",
        people,
        markers,
        geno,
    }
}

/// Byte-exact duplicate samples, one pair inside a family and one across families, so
/// both kinship estimators are exercised at their `phi = 0.5` ceiling.
fn dataset_duplicate() -> Dataset {
    let people = vec![
        Person::founder("D1", "D1_a", 1),
        Person::founder("D1", "D1_copy", 1),
        Person::founder("D2", "D2_a", 2),
        Person::founder("D3", "D3_copy", 2),
        Person::founder("D4", "D4_a", 1),
        Person::founder("D4", "D4_b", 2),
    ];
    let origins = vec![
        Origin::Founder,
        Origin::Duplicate(0),
        Origin::Founder,
        Origin::Duplicate(2),
        Origin::Founder,
        Origin::Founder,
    ];
    let missing = vec![0.03, 0.0, 0.03, 0.0, 0.03, 0.03];
    let mut rng = Rng::new(0x5EED_0005);
    let geno = simulate(&mut rng, &origins, N_MARKERS, &missing);
    Dataset {
        name: "dup_pair",
        people,
        markers: autosomal_markers(N_MARKERS),
        geno,
    }
}

/// Sixty unrelated samples, one per family: 1770 between-family pairs.
fn dataset_unrelated60() -> Dataset {
    let people: Vec<Person> = (1..=60)
        .map(|k| Person::founder(&format!("U{k}"), &format!("U{k}"), 1 + (k % 2) as u8))
        .collect();
    let origins = vec![Origin::Founder; people.len()];
    let missing = vec![0.02; people.len()];
    let mut rng = Rng::new(0x5EED_0006);
    let geno = simulate(&mut rng, &origins, N_MARKERS, &missing);
    Dataset {
        name: "unrelated60",
        people,
        markers: autosomal_markers(N_MARKERS),
        geno,
    }
}

/// Autosomes mixed with every interesting kind of non-autosome label, to prove the
/// inclusion rule end to end rather than only through the console's SNP tally.
///
/// The `25`/`XY` blocks are the ones that matter: the reference counts XY *as an
/// autosome* and analyses it, so getting this wrong changes `N_SNP` on every pair. The
/// `chr7` block is the mirror image — a label a permissive parser would happily map to
/// chromosome 7, which the reference refuses outright.
fn dataset_mixed_chromosomes() -> Dataset {
    let mut ds = dataset_trio();
    ds.name = "mixedchrom";
    let n = ds.people.len();
    let mut rng = Rng::new(0x5EED_0007);
    for (label, count) in [
        ("23", 300),
        ("25", 300),
        ("XY", 100),
        ("26", 100),
        ("chr7", 100),
        ("GL000192.1", 100),
    ] {
        for _ in 0..count {
            let idx = ds.markers.len();
            ds.markers.push(Marker::new(
                label,
                format!("m{idx}"),
                (idx as u32 + 1) * 1000,
            ));
            ds.geno.push(
                (0..n)
                    .map(|_| Some(u8::from(rng.chance(0.5)) + u8::from(rng.chance(0.5))))
                    .collect(),
            );
        }
    }
    ds
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[test]
fn parity_trio() {
    let king = reference_or_skip!("parity_trio");
    assert_eq!(parity(&king, &dataset_trio()), 6 + 9);
}

#[test]
fn parity_nuclear_four_sibs() {
    let king = reference_or_skip!("parity_nuclear_four_sibs");
    assert_eq!(parity(&king, &dataset_nuclear_four_sibs()), 16 + 12);
}

#[test]
fn parity_asymmetric_missingness() {
    let king = reference_or_skip!("parity_asymmetric_missingness");
    assert_eq!(parity(&king, &dataset_asymmetric_missing()), 5 + 16);
}

#[test]
fn parity_monomorphic_markers() {
    let king = reference_or_skip!("parity_monomorphic_markers");
    assert_eq!(parity(&king, &dataset_monomorphic()), 3 + 12);
}

#[test]
fn parity_duplicate_pair() {
    let king = reference_or_skip!("parity_duplicate_pair");
    assert_eq!(parity(&king, &dataset_duplicate()), 2 + 13);
}

#[test]
fn parity_sixty_unrelated() {
    let king = reference_or_skip!("parity_sixty_unrelated");
    assert_eq!(parity(&king, &dataset_unrelated60()), 60 * 59 / 2);
}

#[test]
fn parity_mixed_chromosome_map() {
    let king = reference_or_skip!("parity_mixed_chromosome_map");
    let ds = dataset_mixed_chromosomes();
    // 3000 chr1 + 300 chr25 + 100 XY analysed; 300 X, 100 MT, 100 chr7 and 100 scaffold
    // markers dropped. `parity` asserts this against both the reference's console tally
    // and king-io's kept set before it compares a single pair.
    assert_eq!(
        ds.markers
            .iter()
            .filter(|m| reference_keeps(&m.chrom))
            .count(),
        3400
    );
    assert_eq!(ds.markers.len(), 4000);
    assert_eq!(parity(&king, &ds), 6 + 9);
}

/// A duplicate pair must actually land at the estimator's ceiling, otherwise
/// `parity_duplicate_pair` would be checking a pair of near-unrelated samples and pass
/// vacuously.
#[test]
fn duplicate_pair_is_really_a_duplicate() {
    let ds = dataset_duplicate();
    let dir = workdir("dup_selfcheck");
    let prefix = ds.write(&dir);
    let engine = Engine::load(&prefix);
    let within = engine.counts_for(&[(0, 1)])[0];
    let between = engine.counts_for(&[(2, 3)])[0];

    assert_eq!(
        within.ibs0, 0,
        "an exact copy cannot have opposite homozygotes"
    );
    assert_eq!(between.ibs0, 0);
    assert_eq!(within.het_i, within.het_het);
    assert_eq!(between.het_j, between.het_het);
    let k_within = kinship::kinship(&within, Scope::WithinFamily);
    let k_between = kinship::kinship(&between, Scope::BetweenFamily);
    assert!(
        (k_within - 0.5).abs() < 1e-12,
        "within-family duplicate kinship was {k_within}"
    );
    assert!(
        (k_between - 0.5).abs() < 1e-12,
        "between-family duplicate kinship was {k_between}"
    );
    // The copy is byte-exact, so its call count matches the source exactly.
    assert!(within.n_snp > 2500, "n_snp collapsed to {}", within.n_snp);
}

/// The trio's parent–offspring pairs must have essentially no IBS0 and the sibs must
/// have plenty, otherwise the pedigree datasets are not exercising what they claim to.
#[test]
fn pedigree_datasets_have_the_structure_they_claim() {
    let ds = dataset_nuclear_four_sibs();
    let dir = workdir("nuclear_selfcheck");
    let prefix = ds.write(&dir);
    let engine = Engine::load(&prefix);

    // dad(0) x kid1(2) is parent-offspring; kid1(2) x kid2(3) are full sibs.
    let c = engine.counts_for(&[(0, 2), (2, 3), (0, 1)]);
    let (po, fs, spouses) = (c[0], c[1], c[2]);
    assert_eq!(po.ibs0, 0, "simulated parent-offspring must have IBS0 = 0");
    // Full sibs are IBD-0 at a quarter of markers, and an IBD-0 pair is an opposite
    // homozygote with probability 2p^2q^2, which averages ~0.06 over the U(0.1, 0.9)
    // frequencies used here: ~45 of 3000 markers. The point of the assertion is that it is
    // firmly non-zero while the parent-offspring pair above is exactly zero, which is the
    // fact the PO/FS split relies on.
    assert!(
        fs.ibs0 > 25,
        "full sibs should show many opposite homozygotes, saw {}",
        fs.ibs0
    );
    let k_po = kinship::kinship(&po, Scope::WithinFamily);
    let k_fs = kinship::kinship(&fs, Scope::WithinFamily);
    let k_un = kinship::kinship(&spouses, Scope::WithinFamily);
    assert!((0.20..0.30).contains(&k_po), "PO kinship was {k_po}");
    assert!((0.18..0.32).contains(&k_fs), "FS kinship was {k_fs}");
    assert!(k_un.abs() < 0.05, "unrelated founders had kinship {k_un}");
}

// ---------------------------------------------------------------------------
// SNP inclusion: re-derive the rule from the binary
// ---------------------------------------------------------------------------

/// What one chromosome-label probe measured, on both sides.
#[derive(Debug, PartialEq, Eq)]
struct Probe {
    /// Markers the reference reported as its autosomal analysis set.
    reference_kept: usize,
    /// `N_SNP` the reference printed for the first pair.
    reference_n_snp: usize,
    /// Markers `king-io` loaded under `VariantFilter::Autosomes`.
    our_kept: usize,
    /// `N_SNP` `king-core` computed for the same pair.
    our_n_snp: usize,
}

/// Build a fileset of `anchor` chr-1 markers plus `probe` markers on `chrom`, run both
/// the reference and our engine over it, and report what each one analysed.
fn probe_chromosome(king: &Path, tag: &str, chrom: &str, anchor: usize, probe: usize) -> Probe {
    let people: Vec<Person> = (1..=4)
        .flat_map(|f| {
            [
                Person::founder(&format!("P{f}"), &format!("P{f}_a"), 1),
                Person::founder(&format!("P{f}"), &format!("P{f}_b"), 2),
            ]
        })
        .collect();
    let origins = vec![Origin::Founder; people.len()];
    let mut rng = Rng::new(0x5EED_1000 + tag.len() as u64);
    let geno = simulate(&mut rng, &origins, anchor + probe, &vec![0.0; people.len()]);
    let mut markers = autosomal_markers(anchor);
    for k in 0..probe {
        let idx = markers.len();
        markers.push(Marker::new(chrom, format!("p{k}"), (idx as u32 + 1) * 1000));
    }
    let ds = Dataset {
        name: "probe",
        people,
        markers,
        geno,
    };
    let dir = workdir(&format!("probe_{tag}"));
    let prefix = ds.write(&dir);
    let outdir = dir.join("ref");
    let console = run_reference(king, &prefix, &outdir, &["--ibs"]);
    let ibs = read_table(&outdir.join("king.ibs")).unwrap_or_default();
    let engine = Engine::load(&prefix);
    Probe {
        reference_kept: autosome_snps_reported(&console),
        reference_n_snp: col(&ibs[0], "N_SNP").parse().expect("N_SNP"),
        our_kept: engine.fileset.kept.len(),
        our_n_snp: engine.counts_for(&[(0, 1)])[0].n_snp as usize,
    }
}

/// Which chromosome labels reach the default relatedness analysis.
///
/// Every expectation here was measured against KING 2.3.2. Two of them are the reason
/// this test exists: `25`/`XY` **is** part of the autosome set, and a `chr` prefix makes
/// the reference discard the marker outright.
#[test]
fn snp_inclusion_chromosome_set() {
    let king = reference_or_skip!("snp_inclusion_chromosome_set");
    const ANCHOR: usize = 700;
    const PROBE: usize = 300;

    // (label, whether the probe block joins the autosomal analysis set)
    let cases: &[(&str, bool)] = &[
        ("1", true),
        ("22", true),
        ("25", true),
        ("XY", true),
        ("23", false),
        ("X", false),
        ("x", false),
        ("24", false),
        ("Y", false),
        ("26", false),
        ("MT", false),
        ("27", false),
        ("0", false),
        ("chr1", false),
        ("chrX", false),
        ("M", false),
        ("PAR1", false),
        ("GL000192.1", false),
        ("-1", false),
        ("99", false),
    ];

    for (label, kept) in cases {
        let tag: String = label
            .chars()
            .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
            .collect();
        let expected = if *kept { ANCHOR + PROBE } else { ANCHOR };
        assert_eq!(
            probe_chromosome(&king, &tag, label, ANCHOR, PROBE),
            Probe {
                reference_kept: expected,
                reference_n_snp: expected,
                our_kept: expected,
                our_n_snp: expected,
            },
            "chromosome label {label:?}"
        );
        assert_eq!(
            reference_keeps(label),
            *kept,
            "reference_keeps({label:?}) disagrees with the binary"
        );
    }
}

/// Nothing but the chromosome is filtered: no MAF cutoff, no monomorphic drop, no
/// `0`-allele drop, no per-marker call-rate cutoff.
///
/// Each case adds a 300-marker block of one degenerate kind to a 700-marker polymorphic
/// background; the reference's own SNP count and the pair `N_SNP` are both asserted, so a
/// filter would show up either as a smaller analysis set or as a smaller overlap.
#[test]
fn snp_inclusion_no_content_based_filtering() {
    let king = reference_or_skip!("snp_inclusion_no_content_based_filtering");
    const ANCHOR: usize = 700;
    const BLOCK: usize = 300;

    let people: Vec<Person> = (1..=4)
        .flat_map(|f| {
            [
                Person::founder(&format!("P{f}"), &format!("P{f}_a"), 1),
                Person::founder(&format!("P{f}"), &format!("P{f}_b"), 2),
            ]
        })
        .collect();
    let n = people.len();
    let origins = vec![Origin::Founder; n];

    // (name, A1 allele, the block's genotype row, N_SNP contribution of the block)
    let mut singleton = vec![Some(0); n];
    singleton[0] = Some(1);
    let mut half_called = vec![Some(1); n];
    for g in half_called.iter_mut().skip(n / 2) {
        *g = None;
    }
    let cases: Vec<(&str, &str, Vec<Geno>, usize)> = vec![
        ("monomorphic_hom_a2", "A", vec![Some(0); n], BLOCK),
        ("monomorphic_hom_a1", "A", vec![Some(2); n], BLOCK),
        ("monomorphic_all_het", "A", vec![Some(1); n], BLOCK),
        ("zero_a1_allele", "0", vec![Some(0); n], BLOCK),
        ("zero_a1_but_polymorphic", "0", singleton.clone(), BLOCK),
        ("maf_one_in_sixteen", "A", singleton, BLOCK),
        ("all_calls_missing", "A", vec![None; n], 0),
    ];

    for (name, a1, row, contribution) in cases {
        let mut rng = Rng::new(0x5EED_2000 + name.len() as u64);
        let mut geno = simulate(&mut rng, &origins, ANCHOR, &vec![0.0; n]);
        let mut markers = autosomal_markers(ANCHOR);
        for k in 0..BLOCK {
            let idx = markers.len();
            let mut m = Marker::new("1", format!("d{k}"), (idx as u32 + 1) * 1000);
            m.a1 = a1.to_string();
            markers.push(m);
            geno.push(row.clone());
        }
        let ds = Dataset {
            name: "degenerate",
            people: people.clone(),
            markers,
            geno,
        };
        let dir = workdir(&format!("degen_{name}"));
        let prefix = ds.write(&dir);
        let outdir = dir.join("ref");
        let console = run_reference(&king, &prefix, &outdir, &["--ibs"]);
        let ibs = read_table(&outdir.join("king.ibs")).unwrap_or_default();
        let n_snp: usize = col(&ibs[0], "N_SNP").parse().expect("N_SNP");

        assert_eq!(
            autosome_snps_reported(&console),
            ANCHOR + BLOCK,
            "{name}: the reference dropped markers from the analysis set"
        );
        assert_eq!(
            n_snp,
            ANCHOR + contribution,
            "{name}: pairwise N_SNP disagrees with the no-content-filter rule"
        );

        // And our loader agrees on the analysis set.
        let engine = Engine::load(&prefix);
        assert_eq!(engine.fileset.kept.len(), ANCHOR + BLOCK, "{name}: king-io");
        assert_eq!(
            engine.counts_for(&[(0, 1)])[0].n_snp as usize,
            ANCHOR + contribution,
            "{name}: king-core N_SNP"
        );
    }
}

/// `reference_chrom_code` is a hand-written transcription of measured behaviour; keep its
/// odd corners pinned so a future edit cannot quietly widen it.
#[test]
fn reference_chrom_code_corner_cases() {
    assert_eq!(reference_chrom_code("1"), Some(1));
    assert_eq!(reference_chrom_code("022"), Some(22));
    assert_eq!(reference_chrom_code(" 23"), Some(23));
    assert_eq!(reference_chrom_code("23.0"), Some(23));
    assert_eq!(reference_chrom_code("X"), Some(23));
    assert_eq!(reference_chrom_code("xy"), Some(25));
    assert_eq!(reference_chrom_code("MT"), Some(26));
    assert_eq!(reference_chrom_code("27"), None);
    assert_eq!(reference_chrom_code("0"), None);
    assert_eq!(reference_chrom_code("chr1"), None);
    assert_eq!(reference_chrom_code("M"), None);
    assert_eq!(reference_chrom_code("PAR1"), None);
    assert_eq!(reference_chrom_code("-1"), None);

    assert!(reference_keeps("25") && reference_keeps("XY"));
    assert!(!reference_keeps("23") && !reference_keeps("24") && !reference_keeps("26"));
}

// ---------------------------------------------------------------------------
// Sample exclusion from the between-family analysis
// ---------------------------------------------------------------------------

/// Build `n_markers` fully-called autosomal markers over four two-person families, with
/// sample 0's first `masked` markers set to no-call. Returns whether the reference
/// excluded anybody, and how many between-family rows it emitted.
fn exclusion_probe(king: &Path, tag: &str, n_markers: usize, masked: usize) -> (bool, usize) {
    let people: Vec<Person> = (1..=4)
        .flat_map(|f| {
            [
                Person::founder(&format!("E{f}"), &format!("E{f}_a"), 1),
                Person::founder(&format!("E{f}"), &format!("E{f}_b"), 2),
            ]
        })
        .collect();
    let n = people.len();
    let mut rng = Rng::new(0x5EED_3000);
    let mut geno = simulate(
        &mut rng,
        &vec![Origin::Founder; n],
        n_markers,
        &vec![0.0; n],
    );
    for row in geno.iter_mut().take(masked) {
        row[0] = None;
    }
    let ds = Dataset {
        name: "excl",
        people,
        markers: autosomal_markers(n_markers),
        geno,
    };
    let dir = workdir(&format!("excl_{tag}"));
    let prefix = ds.write(&dir);
    let outdir = dir.join("ref");
    let console = run_reference(king, &prefix, &outdir, &["--kinship"]);
    let rows = read_table(&outdir.join("king.kin0"))
        .unwrap_or_default()
        .len();
    (
        console.contains("are excluded from the kinship analysis"),
        rows,
    )
}

/// The reference drops samples from the **between-family** analysis when they carry too
/// few called autosomal markers, printing `The following N samples are excluded from the
/// kinship analysis (M<512)`. Within-family output is unaffected.
///
/// This is row selection, not arithmetic, so it belongs to whoever writes the output
/// files rather than to the estimator — but it decides which pairs exist at all, so the
/// parity corpus deliberately stays well above the cutoff and the cutoff itself is
/// pinned here.
///
/// `M` is **not** simply the sample's called-marker count, and the exact quantity is
/// still open. Measured (each figure reproducible across independent genotype draws, so
/// the boundary is structural, not statistical):
///
/// | analysed markers | fewest calls still included |
/// |------------------|-----------------------------|
/// | 545 (all called) | 545 — 544 excludes everyone |
/// | 600              | 553                         |
/// | 1000             | 537                         |
/// | 3000             | 521                         |
///
/// So the requirement tightens as the map shrinks and approaches the printed 512 from
/// above as it grows. What this test pins is the part that is unambiguous and that a
/// caller can rely on: the cutoff exists, it is evaluated per sample against that
/// sample's own calls, and excluding a sample removes exactly that sample's
/// between-family pairs. Anything near the boundary is deliberately left to the two
/// full-call points, which are exact.
#[test]
fn between_family_analysis_needs_enough_called_markers() {
    let king = reference_or_skip!("between_family_analysis_needs_enough_called_markers");

    // Marker count alone, no missingness.
    assert_eq!(exclusion_probe(&king, "n400", 400, 0), (true, 0));
    assert_eq!(exclusion_probe(&king, "n544", 544, 0), (true, 0));
    assert_eq!(exclusion_probe(&king, "n545", 545, 0), (false, 24));
    assert_eq!(exclusion_probe(&king, "n3000", 3000, 0), (false, 24));

    // It is per-sample *called* markers, not the size of the map: masking one sample down
    // below the cutoff drops exactly that sample's six between-family pairs, 24 -> 18.
    let (excluded, rows) = exclusion_probe(&king, "masked", 3000, 2600);
    assert!(excluded, "a 400-call sample should have been excluded");
    assert_eq!(rows, 18, "only the low-call sample's pairs should vanish");

    // Just above the cutoff, nobody is dropped.
    let (excluded, rows) = exclusion_probe(&king, "unmasked", 3000, 2000);
    assert!(!excluded, "a 1000-call sample should be kept");
    assert_eq!(rows, 24);
}

// ---------------------------------------------------------------------------
// Non-vacuity: the corpus must be able to tell right rules from plausible wrong ones
// ---------------------------------------------------------------------------

/// Heterozygote count over a sample's **own** non-missing set — the tempting per-sample
/// cache that `docs/VERIFIED_FORMULAS.md` warns against. Used only to prove the parity
/// datasets can detect it.
fn own_het_count(g: &king_io::Genotypes, s: usize) -> u32 {
    g.plane0[s]
        .iter()
        .zip(&g.plane1[s])
        .map(|(&p0, &p1)| (!p0 & p1).count_ones())
        .sum()
}

/// How many markers a sample has a call at.
fn own_call_count(g: &king_io::Genotypes, s: usize) -> u32 {
    g.plane0[s]
        .iter()
        .zip(&g.plane1[s])
        .map(|(&p0, &p1)| (p0 | p1).count_ones())
        .sum()
}

/// If `N_Het1`/`N_Het2` were counted per sample instead of pairwise, the reference's
/// numbers would be unreachable on `asym_missing` — so that dataset's PASS is real
/// evidence for the pairwise rule, not an artifact of everything being fully called.
#[test]
fn asymmetric_dataset_would_expose_a_per_sample_het_count() {
    let ds = dataset_asymmetric_missing();
    let dir = workdir("asym_nonvacuity");
    let prefix = ds.write(&dir);
    let engine = Engine::load(&prefix);
    let n = engine.fileset.samples.len();

    let mut pairs = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            pairs.push((i, j));
        }
    }
    let counts = engine.counts_for(&pairs);
    // The two rules can only agree on a pair where *both* samples are fully called, and
    // this dataset has exactly two such samples, hence exactly one such pair.
    let full: Vec<bool> = (0..n)
        .map(|s| {
            own_call_count(&engine.fileset.genotypes, s) as usize
                == engine.fileset.genotypes.n_variants
        })
        .collect();
    let expected_differ = pairs
        .iter()
        .filter(|&&(i, j)| !(full[i] && full[j]))
        .count();
    assert_eq!(
        full.iter().filter(|&&f| f).count(),
        2,
        "dataset changed: rebuild the expectation below"
    );

    let mut differ = 0;
    for (&(i, j), c) in pairs.iter().zip(&counts) {
        let naive_i = own_het_count(&engine.fileset.genotypes, i);
        let naive_j = own_het_count(&engine.fileset.genotypes, j);
        assert!(c.het_i <= naive_i && c.het_j <= naive_j, "pairwise <= own");
        if c.het_i != naive_i || c.het_j != naive_j {
            differ += 1;
        }
    }
    assert_eq!(
        differ, expected_differ,
        "every pair with a missing call anywhere should distinguish the pairwise count \
         from the per-sample count"
    );
    assert_eq!(expected_differ, 20, "20 of the 21 pairs are discriminating");

    // And the wrong count moves kinship far past the 4th-decimal comparison threshold.
    let (i, j) = (0, 2); // fully-called father vs 30%-missing child
    let c = counts[pairs.iter().position(|&p| p == (i, j)).unwrap()];
    let mut wrong = c;
    wrong.het_i = own_het_count(&engine.fileset.genotypes, i);
    wrong.het_j = own_het_count(&engine.fileset.genotypes, j);
    let right = kinship::kinship(&c, Scope::WithinFamily);
    let bad = kinship::kinship(&wrong, Scope::WithinFamily);
    assert!(
        (right - bad).abs() > 1e-3,
        "per-sample het counts changed kinship by only {}",
        (right - bad).abs()
    );
}

/// The within/between estimator split must matter on this corpus too: applying the
/// within-family formula to between-family pairs has to change the printed value, or
/// `Scope` would be untested by the parity runs.
#[test]
fn corpus_distinguishes_the_two_kinship_estimators() {
    let ds = dataset_asymmetric_missing();
    let dir = workdir("scope_nonvacuity");
    let prefix = ds.write(&dir);
    let engine = Engine::load(&prefix);

    // A1_dad (0) and A3_b (6) are in different families and differ sharply in call rate.
    let c = engine.counts_for(&[(0, 6)])[0];
    assert_ne!(
        c.het_i, c.het_j,
        "the two forms coincide when het_i == het_j"
    );
    let within = fmt_c(kinship::kinship(&c, Scope::WithinFamily), 4);
    let between = fmt_c(kinship::kinship(&c, Scope::BetweenFamily), 4);
    assert_ne!(
        within, between,
        "the two estimators print the same value, so the corpus cannot test the split"
    );
}
