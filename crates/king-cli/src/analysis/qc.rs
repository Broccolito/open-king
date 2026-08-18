//! `--bysample` and `--bySNP` — the two pedigree-driven QC reports.
//!
//! Both walk the same three pedigree relations, both are **space** separated (every
//! relatedness file is tab separated; these two are not), and both grow and shrink their
//! column list with the data. Everything here was established by running the reference
//! binary; the experiments are named at the rule they fix.
//!
//! # What counts as a Mendelian inconsistency
//!
//! The reference splits the two checks, and the split is *not* "the trio check subsumes
//! the pair check":
//!
//! * a **parent–offspring pair** is inconsistent when the two are opposite homozygotes;
//! * a **trio** is inconsistent when the offspring is heterozygous and **both parents are
//!   homozygous for the same allele**, and in no other case.
//!
//! Enumerating all 27 (father, mother, offspring) genotype combinations against the
//! reference — one combination per SNP, read back out of `bySNP.txt` — shows exactly two
//! trio errors, `(aa, aa, Aa)` and `(AA, AA, Aa)`. A combination like `(aa, aa, AA)` is
//! flatly impossible yet `N_errTrio` stays 0 for it: the reference books it as two *pair*
//! errors instead. Implementing the textbook trio rule therefore over-counts, and would
//! be invisible on any correctly simulated pedigree — where neither kind occurs.
//!
//! `N_HetOff`, the denominator of `Err_InHetTrio`, is the tell: the trio statistic is
//! defined over heterozygous offspring, because that is the only offspring genotype whose
//! inconsistency the pair check cannot see.
//!
//! # Which chromosomes take part
//!
//! `bySample.txt` is autosomal throughout (`XY` included, as everywhere else): `sexchr`'s
//! 4 000 autosomal + 150 XY markers give `N_SNP = 4150` and `N_pair = 4 × 4150`. X, Y and
//! MT contribute only their own count/heterozygosity block.
//!
//! `bySNP.txt` covers every retained marker, and the three non-autosomal classes each get
//! their own treatment — all three verified by feeding the same 27 genotype combinations
//! to chromosomes 1, X, Y and MT at once:
//!
//! | Class | PO pairs | Trios |
//! | --- | --- | --- |
//! | autosome, `XY` | counted | counted |
//! | X | counted, same rule | counted, same rule |
//! | Y | counted, same rule | **never**; the five trio columns print as bare `0` |
//! | MT | **never**; the five PO columns print as bare `0` | **never**, likewise |
//!
//! The bare `0` is not a rounding artefact: an autosomal row with no data at all prints
//! `0 0 0 0.0000 0.0000`, so the zero-denominator case is `%.4lf`. Y and MT print `%d`
//! zeros because a different writer emits them.

use std::fmt::Write as _;
use std::io::Write;

use king_io::{bed, Genotypes, Sample, VariantFilter};

use crate::analysis::{f, out_path};
use crate::cli::{Opt, Options};
use crate::console;
use crate::load::{self, Class, Loaded};

/// A parent–offspring pair whose Mendelian error rate must exceed this before the sample
/// is proposed for removal in `MI_Removal`.
///
/// Bracketed by probe filesets at known error rates: a 0.48 % pair rate leaves the column
/// 0, a 1.1 % rate sets it. The reference's exact predicate is unresolved — see the note
/// on [`mi_removal`].
const MI_REMOVAL_RATE: f64 = 0.01;

// ---------------------------------------------------------------------------
// Pedigree
// ---------------------------------------------------------------------------

/// The pedigree relations the two reports are built on, plus the printed parent columns.
pub struct Pedigree {
    /// `(parent, offspring)`, in `.fam` order, father before mother for each offspring.
    pub po: Vec<(usize, usize)>,
    /// `(father, mother, offspring)` for every offspring with both parents genotyped.
    pub trios: Vec<(usize, usize, usize)>,
    /// Full-sibling pairs — reported on the console, used nowhere else.
    pub full_sibs: usize,
    /// The `FA` and `MO` columns as the reference prints them, phantoms included.
    fa_label: Vec<String>,
    mo_label: Vec<String>,
}

impl Pedigree {
    /// Read the relations out of a `.fam`.
    ///
    /// The one surprise is the parent columns. When a sample names exactly one parent the
    /// reference invents the other and prints it as `KING<k>`, `k` counting up from 1 in
    /// `.fam` row order — `dups`' `PO_C`, whose mother is `0`, is printed
    /// `PO_C PO_P KING1`. Reordering the families renumbers them, so the counter follows
    /// the rows and not the IDs. The invented parent gets no row of its own and forms no
    /// trio: `dups` has one PO pair and no trio block.
    pub fn build(samples: &[Sample]) -> Self {
        let n = samples.len();
        let index = |fid: &str, iid: &str| {
            samples
                .iter()
                .position(|s| s.fid == fid && s.iid == iid)
                .filter(|_| iid != "0")
        };
        let mut father = vec![None; n];
        let mut mother = vec![None; n];
        let mut fa_label: Vec<String> = samples.iter().map(|s| s.pat.clone()).collect();
        let mut mo_label: Vec<String> = samples.iter().map(|s| s.mat.clone()).collect();
        let mut phantom = 0usize;
        for (i, s) in samples.iter().enumerate() {
            father[i] = index(&s.fid, &s.pat);
            mother[i] = index(&s.fid, &s.mat);
            let (has_pat, has_mat) = (s.pat != "0", s.mat != "0");
            if has_pat != has_mat {
                phantom += 1;
                let name = format!("KING{phantom}");
                if has_pat {
                    mo_label[i] = name;
                } else {
                    fa_label[i] = name;
                }
            }
        }
        let mut po = Vec::new();
        let mut trios = Vec::new();
        for i in 0..n {
            if let Some(p) = father[i] {
                po.push((p, i));
            }
            if let Some(p) = mother[i] {
                po.push((p, i));
            }
            if let (Some(fa), Some(mo)) = (father[i], mother[i]) {
                trios.push((fa, mo, i));
            }
        }
        let mut full_sibs = 0;
        for a in 0..n {
            for b in a + 1..n {
                if father[a].is_some() && father[a] == father[b] && mother[a] == mother[b] {
                    full_sibs += 1;
                }
            }
        }
        Pedigree {
            po,
            trios,
            full_sibs,
            fa_label,
            mo_label,
        }
    }

    /// `There are N parent-offspring pairs and M trios, and K full-sibling pairs …`.
    ///
    /// No pluralisation anywhere — the `trio` capture prints `1 trios`.
    pub fn console_line(&self) -> String {
        format!(
            "There are {} parent-offspring pairs and {} trios, and {} full-sibling pairs according to the pedigree.\n",
            self.po.len(),
            self.trios.len(),
            self.full_sibs
        )
    }
}

// ---------------------------------------------------------------------------
// Genotype access
// ---------------------------------------------------------------------------

/// Every retained marker's genotypes, which the loader's autosome-only planes do not
/// carry.
///
/// `Loaded::fileset.genotypes` holds the autosomal bit planes alone, because that is what
/// the relatedness kernels want. `bySNP.txt` reports X, Y and MT as well — and so does
/// `--autoQC`, which shares this reader — so this reads the `.bed` a second time under
/// [`VariantFilter::All`] and keeps the classification beside it. The re-read cannot fail:
/// the same file loaded moments ago.
pub(crate) struct Calls {
    planes: Genotypes,
    /// Class of every `.bim` row, parallel to `variants`.
    pub(crate) class: Vec<Class>,
}

impl Calls {
    pub(crate) fn read(opts: &Options, loaded: &Loaded) -> Option<Self> {
        let (bed_path, _, _) = load::paths(opts);
        let sexchr = i64::from(opts.int(Opt::Sexchr));
        let variants = &loaded.fileset.variants;
        let (planes, kept) = bed::read_bed(
            &bed_path,
            loaded.fileset.samples.len(),
            variants,
            VariantFilter::All,
        )
        .ok()?;
        // `All` keeps every row in order, so plane bit `i` is `.bim` row `i`.
        debug_assert!(kept.iter().enumerate().all(|(i, k)| i == *k));
        let class = variants
            .iter()
            .map(|v| load::classify(&v.chrom, sexchr))
            .collect();
        Some(Calls { planes, class })
    }

    /// Dosage of the `.bim` A1 allele at `(sample, variant)`, or `None` when uncalled.
    ///
    /// The planes encode `(1,1)` hom-A1, `(0,1)` het, `(1,0)` hom-A2, `(0,0)` missing.
    pub(crate) fn get(&self, sample: usize, variant: usize) -> Option<u8> {
        let (word, bit) = (variant / 64, variant % 64);
        let p0 = (self.planes.plane0[sample][word] >> bit) & 1;
        let p1 = (self.planes.plane1[sample][word] >> bit) & 1;
        match (p0, p1) {
            (1, 1) => Some(2),
            (0, 1) => Some(1),
            (1, 0) => Some(0),
            _ => None,
        }
    }

    /// Indices of every `.bim` row in the class, in file order.
    pub(crate) fn rows(&self, want: impl Fn(Class) -> bool) -> Vec<usize> {
        (0..self.class.len())
            .filter(|&i| want(self.class[i]))
            .collect()
    }
}

/// Whether a parent–offspring genotype pair is Mendelian-inconsistent: opposite
/// homozygotes, and nothing else.
fn pair_error(parent: u8, offspring: u8) -> bool {
    (parent == 0 && offspring == 2) || (parent == 2 && offspring == 0)
}

/// Whether a trio is Mendelian-inconsistent *in the sense the reference counts*: a
/// heterozygous offspring of two parents homozygous for the same allele.
fn trio_error(father: u8, mother: u8, offspring: u8) -> bool {
    offspring == 1 && father == mother && father != 1
}

// ---------------------------------------------------------------------------
// --bysample
// ---------------------------------------------------------------------------

/// Per-sample tallies over one marker class.
#[derive(Default, Clone)]
struct ClassTally {
    called: Vec<usize>,
    het: Vec<usize>,
}

fn tally(calls: &Calls, rows: &[usize], n: usize) -> ClassTally {
    let mut t = ClassTally {
        called: vec![0; n],
        het: vec![0; n],
    };
    for &v in rows {
        for s in 0..n {
            if let Some(d) = calls.get(s, v) {
                t.called[s] += 1;
                if d == 1 {
                    t.het[s] += 1;
                }
            }
        }
    }
    t
}

/// Run the `--bysample` pass: console body plus `<prefix>bySample.txt`.
pub fn run_bysample(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let ped = Pedigree::build(samples);
    let path = out_path(opts, "bySample.txt");

    let _ = out.write_all(
        format!(
            "QC-by-sample starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(ped.console_line().as_bytes());
    let _ = out.write_all(b"QC starts...\n");

    if let Some(calls) = Calls::read(opts, loaded) {
        write_file(&path, &by_sample_text(samples, &ped, &calls));
    }

    let _ = out.write_all(
        format!(
            "  QC-by-sample ends at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(format!("QC statistics by samples saved in file {path}\n\n").as_bytes());
}

/// The whole of `bySample.txt`, header included.
fn by_sample_text(samples: &[Sample], ped: &Pedigree, calls: &Calls) -> String {
    let n = samples.len();
    let auto = calls.rows(Class::is_autosomal);
    let xs = calls.rows(|c| c == Class::X);
    let ys = calls.rows(|c| c == Class::Y);
    let mts = calls.rows(|c| c == Class::Mt);

    let a = tally(calls, &auto, n);
    let x = tally(calls, &xs, n);
    let y = tally(calls, &ys, n);
    let mt = tally(calls, &mts, n);

    // Mendelian accounting is autosomal only.
    let mut n_pair = vec![0usize; n];
    let mut n_mip = vec![0usize; n];
    let mut n_trio = vec![0usize; n];
    let mut n_mit = vec![0usize; n];
    let mut pair_snps = vec![0usize; ped.po.len()];
    let mut pair_errs = vec![0usize; ped.po.len()];
    for &v in &auto {
        for (k, &(p, c)) in ped.po.iter().enumerate() {
            let (Some(gp), Some(gc)) = (calls.get(p, v), calls.get(c, v)) else {
                continue;
            };
            n_pair[p] += 1;
            n_pair[c] += 1;
            pair_snps[k] += 1;
            if pair_error(gp, gc) {
                n_mip[p] += 1;
                n_mip[c] += 1;
                pair_errs[k] += 1;
            }
        }
        for &(fa, mo, c) in &ped.trios {
            let (Some(gf), Some(gm), Some(gc)) =
                (calls.get(fa, v), calls.get(mo, v), calls.get(c, v))
            else {
                continue;
            };
            for i in [fa, mo, c] {
                n_trio[i] += 1;
            }
            if trio_error(gf, gm, gc) {
                for i in [fa, mo, c] {
                    n_mit[i] += 1;
                }
            }
        }
    }
    let removal = mi_removal(ped, &pair_snps, &pair_errs, n);

    let has_po = !ped.po.is_empty();
    let has_trio = !ped.trios.is_empty();
    let mut text = String::new();
    let mut header = vec![
        "FID",
        "IID",
        "FA",
        "MO",
        "SEX",
        "N_SNP",
        "Missing",
        "Heterozygosity",
    ];
    if !xs.is_empty() {
        header.extend(["N_xSNP", "xHeterozygosity"]);
    }
    if !ys.is_empty() {
        header.extend(["N_ySNP", "N_yHetero"]);
    }
    if !mts.is_empty() {
        header.extend(["N_mtSNP", "N_mtHetero"]);
    }
    if has_po {
        header.extend(["N_pair", "N_MIp", "Err_MIp"]);
    }
    if has_trio {
        header.extend(["N_trio", "N_MIt", "Err_MIt"]);
    }
    if has_po {
        header.push("MI_Removal");
    }
    let _ = writeln!(text, "{}", header.join(" "));

    for i in 0..n {
        let s = &samples[i];
        let mut row = vec![
            s.fid.clone(),
            s.iid.clone(),
            ped.fa_label[i].clone(),
            ped.mo_label[i].clone(),
            s.sex.to_string(),
            a.called[i].to_string(),
            f(rate(auto.len() - a.called[i], auto.len()), 4),
            f(rate(a.het[i], a.called[i]), 4),
        ];
        if !xs.is_empty() {
            row.push(x.called[i].to_string());
            row.push(f(rate(x.het[i], x.called[i]), 4));
        }
        if !ys.is_empty() {
            row.push(y.called[i].to_string());
            // A count, not a rate — the reference's own asymmetry with the X block.
            row.push(y.het[i].to_string());
        }
        if !mts.is_empty() {
            row.push(mt.called[i].to_string());
            row.push(mt.het[i].to_string());
        }
        if has_po {
            row.push(n_pair[i].to_string());
            row.push(n_mip[i].to_string());
            row.push(f(rate(n_mip[i], n_pair[i]), 4));
        }
        if has_trio {
            row.push(n_trio[i].to_string());
            row.push(n_mit[i].to_string());
            row.push(f(rate(n_mit[i], n_trio[i]), 4));
        }
        if has_po {
            row.push(crate::analysis::g(f64::from(u8::from(removal[i]))));
        }
        let _ = writeln!(text, "{}", row.join(" "));
    }
    text
}

/// `x / y`, with the reference's zero-denominator answer of `0`.
fn rate(x: usize, y: usize) -> f64 {
    if y == 0 {
        0.0
    } else {
        x as f64 / y as f64
    }
}

/// The `MI_Removal` flag: the samples whose removal would clear the Mendelian errors.
///
/// **Approximate — the reference's exact predicate is unresolved.** What the probe
/// filesets do establish:
///
/// * it is driven by the *pair* error rate, never the trio one: a fileset whose only
///   errors are trio-type leaves `MI_Removal` 0 on the offspring that carries them and
///   sets it on its two error-free-by-trio siblings, whose pair rates are the high ones;
/// * it is a cover, not a threshold on each sample: with the father carrying every error,
///   only the father is flagged even though each of his three children is above the rate
///   that flags a sample elsewhere;
/// * it is silent below roughly a 1 % pair rate.
///
/// So: repeatedly flag the sample with the highest error rate over its still-uncovered PO
/// pairs while that rate exceeds [`MI_REMOVAL_RATE`], then drop its pairs. That reproduces
/// six of seven probes; the seventh flagged a sample at a 0.18 % rate, which no monotone
/// rule explains and which the known non-determinism in this code path may account for.
/// Every capture in the golden corpus is `0`.
fn mi_removal(ped: &Pedigree, pair_snps: &[usize], pair_errs: &[usize], n: usize) -> Vec<bool> {
    let mut flagged = vec![false; n];
    let mut live: Vec<bool> = vec![true; ped.po.len()];
    loop {
        let mut snps = vec![0usize; n];
        let mut errs = vec![0usize; n];
        for (k, &(p, c)) in ped.po.iter().enumerate() {
            if !live[k] {
                continue;
            }
            for i in [p, c] {
                snps[i] += pair_snps[k];
                errs[i] += pair_errs[k];
            }
        }
        let worst = (0..n)
            .filter(|&i| !flagged[i])
            .max_by(|&i, &j| rate(errs[i], snps[i]).total_cmp(&rate(errs[j], snps[j])));
        let Some(i) = worst.filter(|&i| rate(errs[i], snps[i]) > MI_REMOVAL_RATE) else {
            return flagged;
        };
        flagged[i] = true;
        for (k, &(p, c)) in ped.po.iter().enumerate() {
            if p == i || c == i {
                live[k] = false;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// --bySNP
// ---------------------------------------------------------------------------

/// Run the `--bySNP` pass: console body plus `<prefix>bySNP.txt`.
pub fn run_bysnp(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let ped = Pedigree::build(samples);
    let path = out_path(opts, "bySNP.txt");

    let _ = out.write_all(
        format!(
            "QC-by-SNP starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(ped.console_line().as_bytes());
    let _ = out.write_all(
        format!(
            "Scanning autosomes for QC-by-SNP with {} CPU cores...\n",
            crate::analysis::cpu_count(opts)
        )
        .as_bytes(),
    );

    if let Some(calls) = Calls::read(opts, loaded) {
        // One line per non-autosomal class present, in partition order — the same order
        // their rows take in the file.
        for (class, name) in [(Class::X, "X"), (Class::Y, "Y"), (Class::Mt, "MT")] {
            if !calls.rows(|c| c == class).is_empty() {
                let _ = out
                    .write_all(format!("Scanning chromosome {name} for QC-by-SNP...\n").as_bytes());
            }
        }
        let sexchr = i64::from(opts.int(Opt::Sexchr));
        write_file(&path, &by_snp_text(loaded, &ped, &calls, sexchr));
    }

    let _ = out.write_all(
        format!(
            "QC-by-SNP ends at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(format!("QC statistics by SNPs saved in file {path}\n\n").as_bytes());
}

/// The two QC lines emitted before their shared A1-major orientation gate.
pub fn a1_gate_prelude(opt: Opt, loaded: &Loaded, out: &mut dyn Write) {
    let name = if opt == Opt::Bysample {
        "QC-by-sample"
    } else {
        "QC-by-SNP"
    };
    let _ = out.write_all(
        format!(
            "{name} starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
    let _ = out.write_all(
        Pedigree::build(&loaded.fileset.samples)
            .console_line()
            .as_bytes(),
    );
}

/// The `Chr` column: the class symbol for X, Y and MT, the numeric code for everything
/// else — `XY` prints as its number (`25` under the default `--sexchr`), not as `XY`.
fn chr_symbol(label: &str, sexchr: i64) -> String {
    match load::classify(label, sexchr) {
        Class::X => "X".to_string(),
        Class::Y => "Y".to_string(),
        Class::Mt => "MT".to_string(),
        _ => load::chromosome_code(label, sexchr).to_string(),
    }
}

/// The whole of `bySNP.txt`, header included.
///
/// Rows are the autosomal block in `.bim` order, then X, then Y, then MT — the loader's
/// own partition order, not `.bim` order overall, which `sexchr` shows: its chromosome-25
/// markers print before X even though the `.bim` lists 23 and 24 first.
fn by_snp_text(loaded: &Loaded, ped: &Pedigree, calls: &Calls, sexchr: i64) -> String {
    let samples = &loaded.fileset.samples;
    let variants = &loaded.fileset.variants;
    let n = samples.len();
    let has_po = !ped.po.is_empty();
    let has_trio = !ped.trios.is_empty();

    let mut header = vec![
        "SNP", "Chr", "Pos", "Label_A", "Label_a", "Freq_A", "N", "N_AA", "N_Aa", "N_aa",
        "CallRate",
    ];
    if has_po {
        header.extend(["N_PO", "N_HomPO", "N_errPO", "Err_InPO", "Err_InHomPO"]);
    }
    if has_trio {
        header.extend([
            "N_trio",
            "N_HetOff",
            "N_errTrio",
            "Err_InTrio",
            "Err_InHetTrio",
        ]);
    }
    let mut text = String::new();
    let _ = writeln!(text, "{}", header.join(" "));

    let order = calls
        .rows(Class::is_autosomal)
        .into_iter()
        .chain(calls.rows(|c| c == Class::X))
        .chain(calls.rows(|c| c == Class::Y))
        .chain(calls.rows(|c| c == Class::Mt));

    for v in order {
        let class = calls.class[v];
        let var = &variants[v];
        let g: Vec<Option<u8>> = (0..n).map(|s| calls.get(s, v)).collect();
        let count = |d: u8| g.iter().filter(|x| **x == Some(d)).count();
        let (n_aa, n_ab, n_bb) = (count(2), count(1), count(0));
        let called = n_aa + n_ab + n_bb;
        let mut row = vec![
            var.id.clone(),
            chr_symbol(&var.chrom, sexchr),
            var.bp.to_string(),
            var.a1.clone(),
            var.a2.clone(),
            f(rate(2 * n_aa + n_ab, 2 * called), 4),
            called.to_string(),
            n_aa.to_string(),
            n_ab.to_string(),
            n_bb.to_string(),
            f(rate(called, n), 4),
        ];
        if has_po {
            // MT carries no parent–offspring statistic at all, and prints integer zeros
            // rather than the `%.4lf` an empty autosomal row would print.
            if class == Class::Mt {
                row.extend(["0", "0", "0", "0", "0"].map(str::to_string));
            } else {
                let (mut n_po, mut n_hom, mut n_err) = (0usize, 0usize, 0usize);
                for &(p, c) in &ped.po {
                    let (Some(gp), Some(gc)) = (g[p], g[c]) else {
                        continue;
                    };
                    n_po += 1;
                    if gp != 1 && gc != 1 {
                        n_hom += 1;
                    }
                    if pair_error(gp, gc) {
                        n_err += 1;
                    }
                }
                row.push(n_po.to_string());
                row.push(n_hom.to_string());
                row.push(n_err.to_string());
                row.push(f(rate(n_err, n_po), 4));
                row.push(f(rate(n_err, n_hom), 4));
            }
        }
        if has_trio {
            if matches!(class, Class::Y | Class::Mt) {
                row.extend(["0", "0", "0", "0", "0"].map(str::to_string));
            } else {
                let (mut n_t, mut n_het, mut n_err) = (0usize, 0usize, 0usize);
                for &(fa, mo, c) in &ped.trios {
                    let (Some(gf), Some(gm), Some(gc)) = (g[fa], g[mo], g[c]) else {
                        continue;
                    };
                    n_t += 1;
                    if gc == 1 {
                        n_het += 1;
                    }
                    if trio_error(gf, gm, gc) {
                        n_err += 1;
                    }
                }
                row.push(n_t.to_string());
                row.push(n_het.to_string());
                row.push(n_err.to_string());
                row.push(f(rate(n_err, n_t), 4));
                row.push(f(rate(n_err, n_het), 4));
            }
        }
        let _ = writeln!(text, "{}", row.join(" "));
    }
    text
}

/// Write an output file, ignoring failure exactly as the reference does.
fn write_file(path: &str, text: &str) {
    let _ = std::fs::write(path, text);
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

    #[test]
    fn pair_rule_is_opposite_homozygotes_only() {
        assert!(pair_error(0, 2));
        assert!(pair_error(2, 0));
        for (p, c) in [(0, 0), (0, 1), (1, 0), (1, 1), (1, 2), (2, 1), (2, 2)] {
            assert!(!pair_error(p, c), "({p},{c}) is not a pair error");
        }
    }

    /// The reference's own answer for all 27 combinations, read out of `bySNP.txt`.
    #[test]
    fn trio_rule_fires_on_exactly_two_combinations() {
        let mut hits = Vec::new();
        for fa in 0..3u8 {
            for mo in 0..3u8 {
                for c in 0..3u8 {
                    if trio_error(fa, mo, c) {
                        hits.push((fa, mo, c));
                    }
                }
            }
        }
        assert_eq!(hits, [(0, 0, 1), (2, 2, 1)]);
    }

    #[test]
    fn phantom_parents_number_in_fam_order() {
        let ped = Pedigree::build(&[
            sample("FA", "P1", "0", "0"),
            sample("FA", "C1", "P1", "0"),
            sample("FB", "M2", "0", "0"),
            sample("FB", "C2", "0", "M2"),
            sample("FC", "U", "0", "0"),
        ]);
        assert_eq!(ped.mo_label[1], "KING1");
        assert_eq!(ped.fa_label[3], "KING2");
        // A founder keeps its literal zeros, and neither phantom makes a trio.
        assert_eq!(
            (ped.fa_label[4].as_str(), ped.mo_label[4].as_str()),
            ("0", "0")
        );
        assert!(ped.trios.is_empty());
        assert_eq!(ped.po, [(0, 1), (2, 3)]);
    }

    #[test]
    fn console_line_does_not_pluralise() {
        let ped = Pedigree::build(&[
            sample("T", "F", "0", "0"),
            sample("T", "M", "0", "0"),
            sample("T", "C", "F", "M"),
        ]);
        assert_eq!(
            ped.console_line(),
            "There are 2 parent-offspring pairs and 1 trios, and 0 full-sibling pairs according to the pedigree.\n"
        );
    }

    #[test]
    fn full_sibs_need_both_parents_shared() {
        let ped = Pedigree::build(&[
            sample("F", "FA", "0", "0"),
            sample("F", "MO", "0", "0"),
            sample("F", "M2", "0", "0"),
            sample("F", "C1", "FA", "MO"),
            sample("F", "C2", "FA", "MO"),
            sample("F", "H1", "FA", "M2"),
        ]);
        assert_eq!(ped.full_sibs, 1);
    }

    #[test]
    fn zero_denominators_render_as_zero() {
        assert_eq!(f(rate(0, 0), 4), "0.0000");
        assert_eq!(f(rate(3, 4), 4), "0.7500");
    }
}
