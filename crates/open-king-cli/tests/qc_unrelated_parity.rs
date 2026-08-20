//! Byte-for-byte parity for `--bysample`, `--bySNP` and `--unrelated`.
//!
//! Every expectation is a transcript of the reference binary's answer for a fileset this
//! test writes from scratch, so the suite pins the rules without needing KING 2.3.2 on
//! the machine. Each fixture exists for a rule that the golden corpus cannot show,
//! because a correctly simulated pedigree never produces a Mendelian error:
//!
//! * [`mendelian_rules`] — the 27 (father, mother, offspring) genotype combinations, and
//!   with them the split between the pair check and the trio check;
//! * [`non_autosomal_blocks`] — Y with no trio statistic and MT with neither, both
//!   printing bare `%d` zeros rather than the `%.4lf` an empty autosomal row prints;
//! * [`phantom_parent_is_numbered`] — the invented `KING<k>` parent;
//! * [`tiny_dataset_takes_the_pedigree_path`] and [`unrelated_visit_order`] — the two
//!   `--unrelated` regimes.

use std::path::{Path, PathBuf};
use std::process::Command;

/// PLINK's two-bit codes, by A1 dosage.
const HOM_A1: u8 = 0b00;
const HET: u8 = 0b10;
const HOM_A2: u8 = 0b11;

struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let dir = std::env::temp_dir().join(format!("open-king-qc-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch dir");
        Scratch(dir)
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A `.fam` row: `FID IID FA MO SEX`.
struct Row<'a>(&'a str, &'a str, &'a str, &'a str, u8);

/// Write a fileset whose genotypes are given per variant, one code per sample.
///
/// `chroms` names one chromosome per variant; markers are 1 Mb apart within a chromosome
/// so the map is never degenerate.
fn write_fileset(dir: &Path, stem: &str, rows: &[Row], chroms: &[&str], codes: &[Vec<u8>]) {
    let fam: String = rows
        .iter()
        .map(|r| format!("{} {} {} {} {} -9\n", r.0, r.1, r.2, r.3, r.4))
        .collect();
    std::fs::write(dir.join(format!("{stem}.fam")), fam).unwrap();

    let mut bim = String::new();
    let mut seen: Vec<(&str, i64)> = Vec::new();
    for (i, c) in chroms.iter().enumerate() {
        let slot = match seen.iter_mut().find(|(name, _)| name == c) {
            Some(s) => s,
            None => {
                seen.push((c, 0));
                seen.last_mut().unwrap()
            }
        };
        slot.1 += 1;
        let bp = slot.1 * 1_000_000;
        bim.push_str(&format!("{c}\trs{i}\t{}\t{bp}\tA\tG\n", bp as f64 / 1e6));
    }
    std::fs::write(dir.join(format!("{stem}.bim")), bim).unwrap();

    let mut bed = vec![0x6c, 0x1b, 0x01];
    for row in codes {
        assert_eq!(row.len(), rows.len());
        let mut packed = vec![0u8; rows.len().div_ceil(4)];
        for (i, code) in row.iter().enumerate() {
            packed[i / 4] |= code << (2 * (i % 4));
        }
        bed.extend_from_slice(&packed);
    }
    std::fs::write(dir.join(format!("{stem}.bed")), bed).unwrap();
}

/// Run the binary in `dir` and return its stdout.
fn run(dir: &Path, stem: &str, args: &[&str]) -> String {
    let bed = dir.join(format!("{stem}.bed"));
    let out = Command::new(env!("CARGO_BIN_EXE_open-king"))
        .current_dir(dir)
        .arg("-b")
        .arg(&bed)
        .args(args)
        .output()
        .expect("run open-king");
    String::from_utf8_lossy(&out.stdout).into_owned()
}

fn read(dir: &Path, name: &str) -> String {
    std::fs::read_to_string(dir.join(name)).unwrap_or_default()
}

/// The columns of one `bySNP.txt` row, whitespace split.
fn snp_row<'a>(text: &'a str, id: &str) -> Vec<&'a str> {
    text.lines()
        .find(|l| l.split(' ').next() == Some(id))
        .unwrap_or_else(|| panic!("no row for {id}"))
        .split(' ')
        .collect()
}

/// One marker's genotype codes plus the (father, mother, offspring) dosages behind them.
struct Combination {
    codes: Vec<u8>,
    dosages: (u8, u8, u8),
}

/// All 27 (father, mother, offspring) combinations, one per marker.
fn combinations() -> Vec<Combination> {
    let dose = [HOM_A2, HET, HOM_A1]; // 0, 1, 2 copies of A1
    let mut out = Vec::new();
    for f in 0..3u8 {
        for m in 0..3u8 {
            for c in 0..3u8 {
                out.push(Combination {
                    codes: vec![
                        dose[f as usize],
                        dose[m as usize],
                        dose[c as usize],
                        HOM_A2,
                        HOM_A2,
                    ],
                    dosages: (f, m, c),
                });
            }
        }
    }
    out
}

/// A parent–offspring pair is inconsistent when the two are opposite homozygotes; a trio
/// only when the offspring is heterozygous and both parents are homozygous for the same
/// allele. `(aa, aa, AA)` is impossible yet books as two pair errors and no trio error.
#[test]
fn mendelian_rules() {
    let s = Scratch::new("mi");
    let rows = [
        Row("F", "FA", "0", "0", 1),
        Row("F", "MO", "0", "0", 2),
        Row("F", "KID", "FA", "MO", 1),
        Row("G", "U1", "0", "0", 1),
        Row("G", "U2", "0", "0", 2),
    ];
    let combos = combinations();
    let codes: Vec<Vec<u8>> = combos.iter().map(|c| c.codes.clone()).collect();
    let chroms = vec!["1"; codes.len()];
    write_fileset(&s.0, "mi", &rows, &chroms, &codes);
    run(&s.0, "mi", &["--bySNP"]);
    let text = read(&s.0, "kingbySNP.txt");

    for (i, combo) in combos.iter().enumerate() {
        let (f, m, c) = (&combo.dosages.0, &combo.dosages.1, &combo.dosages.2);
        let row = snp_row(&text, &format!("rs{i}"));
        let n_err_po: u32 = row[13].parse().unwrap();
        let n_err_trio: u32 = row[18].parse().unwrap();
        let n_het_off: u32 = row[17].parse().unwrap();

        let expect_po = u32::from(*f == 0 && *c == 2 || *f == 2 && *c == 0)
            + u32::from(*m == 0 && *c == 2 || *m == 2 && *c == 0);
        let expect_trio = u32::from(*c == 1 && f == m && *f != 1);
        assert_eq!(n_err_po, expect_po, "N_errPO for ({f},{m},{c})");
        assert_eq!(n_err_trio, expect_trio, "N_errTrio for ({f},{m},{c})");
        assert_eq!(n_het_off, u32::from(*c == 1), "N_HetOff for ({f},{m},{c})");
    }
}

/// Y carries the pair statistic but never the trio one; MT carries neither. Both print
/// the columns they do not compute as bare `0`, which an empty autosomal row does not.
#[test]
fn non_autosomal_blocks() {
    let s = Scratch::new("chr");
    let rows = [
        Row("F", "FA", "0", "0", 1),
        Row("F", "MO", "0", "0", 2),
        Row("F", "KID", "FA", "MO", 1),
        Row("G", "U1", "0", "0", 1),
        Row("G", "U2", "0", "0", 2),
    ];
    // One marker per class, all five samples hom A2, plus autosomal padding so the map
    // is loadable.
    let chroms = ["1", "1", "1", "1", "23", "24", "26"];
    let codes: Vec<Vec<u8>> = chroms.iter().map(|_| vec![HOM_A2; 5]).collect();
    write_fileset(&s.0, "chr", &rows, &chroms, &codes);
    run(&s.0, "chr", &["--bySNP"]);
    let text = read(&s.0, "kingbySNP.txt");

    let auto = snp_row(&text, "rs0");
    assert_eq!(&auto[1..3], ["1", "1000000"]);
    assert_eq!(
        &auto[11..21],
        ["2", "2", "0", "0.0000", "0.0000", "1", "0", "0", "0.0000", "0.0000"]
    );

    let x = snp_row(&text, "rs4");
    assert_eq!(x[1], "X");
    assert_eq!(
        &x[11..21],
        ["2", "2", "0", "0.0000", "0.0000", "1", "0", "0", "0.0000", "0.0000"]
    );

    let y = snp_row(&text, "rs5");
    assert_eq!(y[1], "Y");
    assert_eq!(
        &y[11..16],
        ["2", "2", "0", "0.0000", "0.0000"],
        "Y keeps the pair block"
    );
    assert_eq!(
        &y[16..21],
        ["0", "0", "0", "0", "0"],
        "Y prints integer zeros for trios"
    );

    let mt = snp_row(&text, "rs6");
    assert_eq!(mt[1], "MT");
    assert_eq!(
        &mt[11..21],
        ["0", "0", "0", "0", "0", "0", "0", "0", "0", "0"]
    );

    // Rows are grouped by class, not left in `.bim` order.
    let order: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|l| l.split(' ').nth(1).unwrap())
        .collect();
    assert_eq!(order, ["1", "1", "1", "1", "X", "Y", "MT"]);
}

/// A sample naming exactly one parent gets the other invented as `KING<k>`, numbered in
/// `.fam` row order. It forms no trio, so the trio block never appears.
#[test]
fn phantom_parent_is_numbered() {
    let s = Scratch::new("phantom");
    let rows = [
        Row("A", "P1", "0", "0", 1),
        Row("A", "C1", "P1", "0", 2),
        Row("B", "M2", "0", "0", 2),
        Row("B", "C2", "0", "M2", 1),
        Row("C", "U", "0", "0", 1),
    ];
    let chroms = ["1"; 8];
    let codes: Vec<Vec<u8>> = chroms.iter().map(|_| vec![HOM_A2; 5]).collect();
    write_fileset(&s.0, "phantom", &rows, &chroms, &codes);
    let stdout = run(&s.0, "phantom", &["--bysample"]);
    let text = read(&s.0, "kingbySample.txt");

    assert!(stdout.contains(
        "There are 2 parent-offspring pairs and 0 trios, and 0 full-sibling pairs according to the pedigree.\n"
    ));
    let header = text.lines().next().unwrap();
    assert!(
        header.ends_with("N_pair N_MIp Err_MIp MI_Removal"),
        "{header}"
    );
    assert!(!header.contains("N_trio"));
    assert!(text.contains("A C1 P1 KING1 2 "), "{text}");
    assert!(text.contains("B C2 KING2 M2 1 "), "{text}");
    assert!(text.contains("C U 0 0 1 "), "{text}");
}

/// Under ten samples the reference disables clustering, prints so, and still writes both
/// files from the pedigree alone: a declared chain loses everything but its ends.
#[test]
fn tiny_dataset_takes_the_pedigree_path() {
    let s = Scratch::new("tiny");
    let rows = [
        Row("F", "A", "0", "0", 1),
        Row("F", "B", "A", "0", 1),
        Row("F", "C", "B", "0", 1),
        Row("S0", "S0", "0", "0", 1),
        Row("S1", "S1", "0", "0", 1),
    ];
    let chroms = ["1"; 8];
    let codes: Vec<Vec<u8>> = chroms.iter().map(|_| vec![HOM_A2; 5]).collect();
    write_fileset(&s.0, "tiny", &rows, &chroms, &codes);
    let stdout = run(&s.0, "tiny", &["--unrelated"]);

    assert!(stdout
        .contains("This function is currently disabled for tiny dataset with sample size < 10.\n"));
    assert!(!stdout.contains("Sorting autosomes"));
    // The grandchild is 2nd degree to the kept founder, so both descendants go.
    assert_eq!(read(&s.0, "kingunrelated.txt"), "F\tA\nS0\tS0\nS1\tS1\n");
    assert_eq!(read(&s.0, "kingunrelated_toberemoved.txt"), "F\tB\nF\tC\n");
    assert!(stdout.contains("A list of 3 unrelated individuals saved in file kingunrelated.txt\n"));
    assert!(stdout.contains(
        "An alternative list of 2 to-be-removed individuals saved in file kingunrelated_toberemoved.txt\n"
    ));
}

/// Ten mutually unrelated members of one FID are all kept, in the measured visit order
/// rather than in ID order.
#[test]
fn unrelated_visit_order() {
    let s = Scratch::new("order");
    let rows: Vec<Row> = (1..=10)
        .map(|i| {
            let iid: &'static str = Box::leak(format!("I{i:02}").into_boxed_str());
            Row("POOL", iid, "0", "0", 1)
        })
        .collect();
    // Every sample its own private heterozygous marker: no two share anything, so no pair
    // is anywhere near the 4th-degree edge.
    let mut codes = Vec::new();
    for i in 0..10 {
        let mut row = vec![HOM_A2; 10];
        row[i] = HET;
        codes.push(row.clone());
        codes.push(row);
    }
    let chroms = vec!["1"; codes.len()];
    write_fileset(&s.0, "order", &rows, &chroms, &codes);
    let stdout = run(&s.0, "order", &["--unrelated"]);

    assert!(stdout.contains("Clustering up to 1st-degree relatives in families...\n"));
    assert!(stdout.contains("Individual IDs are unique across all families.\n"));
    assert!(stdout.contains("No families were found to be connected.\n"));
    let kept: Vec<String> = read(&s.0, "kingunrelated.txt")
        .lines()
        .map(|l| l.split('\t').nth(1).unwrap().to_string())
        .collect();
    let want: Vec<String> = [4, 3, 10, 9, 2, 5, 1, 8, 6, 7]
        .iter()
        .map(|i| format!("I{i:02}"))
        .collect();
    assert_eq!(kept, want);
    assert_eq!(read(&s.0, "kingunrelated_toberemoved.txt"), "");
    // `--degree` moves the console line and nothing else.
    let two = run(&s.0, "order", &["--unrelated", "--degree", "2"]);
    assert!(two.contains("Clustering up to 2nd-degree relatives in families...\n"));
    let after: Vec<String> = read(&s.0, "kingunrelated.txt")
        .lines()
        .map(|l| l.split('\t').nth(1).unwrap().to_string())
        .collect();
    assert_eq!(after, want);
}
