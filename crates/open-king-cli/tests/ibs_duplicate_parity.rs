//! Byte-for-byte parity for the `--ibs` and `--duplicate` passes.
//!
//! The expectations below are transcripts of rules established against the reference
//! binary (`docs/BEHAVIOR.md`, plus the captures under `tests/parity/golden`), pinned
//! here on filesets this test writes from scratch so the suite runs without KING 2.3.2
//! on the machine. The parity harness covers the reference's own datasets; these tests
//! cover the *decisions* — file existence, column sets, row order and the selection
//! threshold — on inputs small enough to reason about by hand.

use std::path::PathBuf;
use std::process::Command;

// ---------------------------------------------------------------------------
// Fixtures
// ---------------------------------------------------------------------------

/// Genotype codes as PLINK packs them, two bits per sample.
const HOM1: u8 = 0b00;
const HET: u8 = 0b10;
const HOM2: u8 = 0b11;

/// A scratch directory that cleans itself up.
struct Scratch(PathBuf);

impl Scratch {
    fn new(tag: &str) -> Scratch {
        let dir = std::env::temp_dir().join(format!("king-ibs-{}-{tag}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(&dir).expect("scratch dir");
        Scratch(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    fn read(&self, name: &str) -> String {
        std::fs::read_to_string(self.path(name))
            .unwrap_or_else(|_| panic!("{name} was not written"))
    }

    fn exists(&self, name: &str) -> bool {
        self.path(name).exists()
    }
}

impl Drop for Scratch {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// Write a `.bed`/`.bim`/`.fam` triple from an explicit `[sample][variant]` matrix.
///
/// Markers are placed 50 kb apart on chromosome 1 starting at 1 Mb, which is dense enough
/// for the segment cut (`<= 156 250 bp`) never to fire, so the only thing that decides
/// whether a segment survives is how many complete 64-marker words the run holds.
fn write_fileset(scratch: &Scratch, stem: &str, families: &[(&str, &str)], genotypes: &[Vec<u8>]) {
    let n_variants = genotypes[0].len();
    assert!(genotypes.iter().all(|g| g.len() == n_variants));

    let fam: String = families
        .iter()
        .map(|(fid, iid)| format!("{fid} {iid} 0 0 1 1\n"))
        .collect();
    std::fs::write(scratch.path(&format!("{stem}.fam")), fam).expect("write fam");

    let bim: String = (0..n_variants)
        .map(|i| format!("1 rs{i} 0 {} A G\n", 1_000_000 + 50_000 * i))
        .collect();
    std::fs::write(scratch.path(&format!("{stem}.bim")), bim).expect("write bim");

    let mut bed = vec![0x6c, 0x1b, 0x01];
    for v in 0..n_variants {
        let mut row = vec![0u8; families.len().div_ceil(4)];
        for (s, g) in genotypes.iter().enumerate() {
            row[s / 4] |= g[v] << (2 * (s % 4));
        }
        bed.extend_from_slice(&row);
    }
    std::fs::write(scratch.path(&format!("{stem}.bed")), bed).expect("write bed");
}

/// Pseudo-random genotypes for one sample, at Hardy-Weinberg proportions for `p = 0.5`.
///
/// Two samples drawn with different seeds are unrelated to within sampling noise, which
/// keeps every synthetic pair below the IBD2 gate; drawing twice with the same seed makes
/// an exact duplicate. A short repeating pattern will not do: it correlates the samples
/// enough to push ordinary pairs over the gate.
fn pattern(seed: usize, n: usize) -> Vec<u8> {
    let mut state = 0x2545_F491_4F6C_DD1Du64.wrapping_mul(seed as u64 + 1);
    (0..n)
        .map(|_| {
            state = state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1_442_695_040_888_963_407);
            match (state >> 33) % 4 {
                0 => HOM1,
                3 => HOM2,
                _ => HET,
            }
        })
        .collect()
}

/// Run the binary in the scratch directory and return its stdout from `Options in
/// effect:` onwards, with timestamps blanked.
fn run(scratch: &Scratch, stem: &str, extra: &[&str]) -> String {
    let bed = scratch.path(&format!("{stem}.bed"));
    let mut args: Vec<String> = vec!["-b".into(), bed.to_string_lossy().into_owned()];
    args.extend(extra.iter().map(|s| (*s).to_string()));
    let out = Command::new(env!("CARGO_BIN_EXE_open-king"))
        .args(&args)
        .current_dir(&scratch.0)
        .output()
        .expect("king binary runs");
    assert!(out.stderr.is_empty(), "nothing is written to stderr");
    let stdout = String::from_utf8(out.stdout).expect("stdout is UTF-8");
    let start = stdout
        .find("Options in effect:")
        .unwrap_or_else(|| panic!("no analysis body in:\n{stdout}"));
    blank_timestamps(&stdout[start..])
}

/// Replace every ` at <ctime>` tail with a fixed token.
fn blank_timestamps(s: &str) -> String {
    s.lines()
        .map(|line| match line.rfind(" at ") {
            Some(i) if line.len() > i + 4 => format!("{} at <TS>", &line[..i]),
            _ => line.to_string(),
        })
        .collect::<Vec<_>>()
        .join("\n")
}

fn header(text: &str) -> Vec<&str> {
    text.lines()
        .next()
        .expect("a header line")
        .split('\t')
        .collect()
}

fn rows(text: &str) -> Vec<Vec<&str>> {
    text.lines()
        .skip(1)
        .map(|l| l.split('\t').collect())
        .collect()
}

// ---------------------------------------------------------------------------
// --duplicate
// ---------------------------------------------------------------------------

/// Four samples, two families, the first pair an exact copy of each other.
fn duplicate_fixture(scratch: &Scratch) {
    let n = 400;
    let g = vec![pattern(0, n), pattern(0, n), pattern(2, n), pattern(4, n)];
    write_fileset(
        scratch,
        "d",
        &[("F1", "A"), ("F1", "ACOPY"), ("F2", "B"), ("F2", "C")],
        &g,
    );
}

/// The whole `--duplicate` body, including the two blank lines that frame the
/// unconfirmed-pairs count and the `Stage 2` line a below-threshold sample count always
/// prints.
#[test]
fn duplicate_console_body() {
    let s = Scratch::new("dup-console");
    duplicate_fixture(&s);
    assert_eq!(
        run(&s, "d", &["--duplicate"]),
        concat!(
            "Options in effect:\n",
            "\t--duplicate\n",
            "\n",
            "Sorting autosomes...\n",
            "Computing pairwise genotype concordance starts at <TS>\n",
            "  <N> CPU cores are used...\n",
            "        Stage 2 (with all SNPs) inference ends at <TS>\n",
            "1 pairs of duplicates with heterozygote concordance rate > 80% are saved in file king.con\n",
            "\n",
            "  5 additional pairs from screening stage not confirmed in the final stage\n",
            "\n",
            "KING ends at <TS>\n",
        )
        .replace("<N>", &cpu_line(&run(&s, "d", &["--duplicate"])))
    );
}

/// The CPU count in the body is a property of the host, so read it back rather than
/// pinning it.
fn cpu_line(body: &str) -> String {
    body.lines()
        .find(|l| l.contains("CPU cores are used"))
        .and_then(|l| l.split_whitespace().next())
        .expect("a CPU line")
        .to_string()
}

/// `.con` carries both within- and cross-family pairs, and its three concordance columns
/// are the only `%.5lf` fields the program emits.
#[test]
fn con_layout_and_five_decimals() {
    let s = Scratch::new("dup-con");
    duplicate_fixture(&s);
    run(&s, "d", &["--duplicate"]);
    let con = s.read("king.con");
    assert_eq!(
        header(&con),
        [
            "FID1", "ID1", "FID2", "ID2", "N", "N_IBS0", "N_IBS1", "N_IBS2", "Concord", "HomConc",
            "HetConc"
        ]
    );
    // An exact copy: every site IBS2, so all three rates are 1 at five decimals.
    assert_eq!(
        rows(&con),
        [["F1", "A", "F1", "ACOPY", "400", "0", "0", "400", "1.00000", "1.00000", "1.00000"]]
    );
}

/// Selection is **strictly** above `--minConc`, so an exact duplicate fails `1`.
///
/// The file is still written, with its header alone, and every pair is reported as an
/// unconfirmed candidate.
#[test]
fn minconc_is_a_strict_threshold() {
    let s = Scratch::new("dup-minconc");
    duplicate_fixture(&s);
    let body = run(&s, "d", &["--duplicate", "--minConc", "1"]);
    assert!(
        body.contains("\t--minConc 1\n"),
        "the value is echoed as %G:\n{body}"
    );
    assert!(
        body.contains("No duplicates are found with heterozygote concordance rate > 100%.\n"),
        "{body}"
    );
    assert!(
        body.contains(
            "  6 additional pairs from screening stage not confirmed in the final stage\n"
        ),
        "{body}"
    );
    assert_eq!(s.read("king.con").lines().count(), 1, "header only");

    // ...and a threshold of 0 keeps every pair, with no unconfirmed ones left to report.
    let s = Scratch::new("dup-minconc0");
    duplicate_fixture(&s);
    let body = run(&s, "d", &["--duplicate", "--minConc", "0"]);
    assert!(body.contains("6 pairs of duplicates"), "{body}");
    assert!(!body.contains("additional pairs"), "{body}");
}

/// `--prefix` is concatenated, not joined with a dot.
#[test]
fn duplicate_honours_the_prefix() {
    let s = Scratch::new("dup-prefix");
    duplicate_fixture(&s);
    let body = run(&s, "d", &["--duplicate", "--prefix", "ZZ_"]);
    assert!(body.contains("saved in file ZZ_.con"), "{body}");
    assert!(s.exists("ZZ_.con"));
}

// ---------------------------------------------------------------------------
// --ibs
// ---------------------------------------------------------------------------

/// Ten samples, each its own family, over a map too short for segment analysis.
fn singleton_families(scratch: &Scratch, n_variants: usize) {
    let families: Vec<(String, String)> = (0..10)
        .map(|i| (format!("F{i:02}"), format!("S{i:02}")))
        .collect();
    let borrowed: Vec<(&str, &str)> = families
        .iter()
        .map(|(f, i)| (f.as_str(), i.as_str()))
        .collect();
    let g: Vec<Vec<u8>> = (0..10).map(|i| pattern(i, n_variants)).collect();
    write_fileset(scratch, "s", &borrowed, &g);
}

/// A map whose usable segments total under 100 Mb loses `MaxIBD2`/`Pr_IBD2` from **both**
/// files, and says so on the console. 400 markers 50 kb apart is one 6-word segment of
/// ~20 Mb: usable, but far short.
#[test]
fn short_segments_drop_the_trailing_columns() {
    let s = Scratch::new("ibs-short");
    singleton_families(&s, 400);
    let body = run(&s, "s", &["--ibs"]);
    assert!(body.contains("Segments too short.\n"), "{body}");
    assert!(
        body.contains(
            "Total length of 1 chromosomal segments usable for IBD segment analysis is 19.9 Mb.\n"
        ),
        "{body}"
    );
    assert_eq!(*header(&s.read("king.ibs")).last().unwrap(), "Kinship");
    assert_eq!(*header(&s.read("king.ibs0")).last().unwrap(), "Kinship");
    assert_eq!(header(&s.read("king.ibs")).len(), 20);
    assert_eq!(header(&s.read("king.ibs0")).len(), 19);
}

/// At 100 Mb the columns appear on both files. 2001 markers 50 kb apart span exactly
/// 100 000 000 bp, which is the inclusive end of the bracketed boundary.
#[test]
fn a_hundred_megabases_adds_the_trailing_columns() {
    let s = Scratch::new("ibs-long");
    singleton_families(&s, 2001);
    let body = run(&s, "s", &["--ibs"]);
    assert!(!body.contains("Segments too short."), "{body}");
    assert_eq!(
        header(&s.read("king.ibs"))
            .iter()
            .rev()
            .take(2)
            .collect::<Vec<_>>(),
        [&"Pr_IBD2", &"MaxIBD2"]
    );
    assert_eq!(
        header(&s.read("king.ibs0"))
            .iter()
            .rev()
            .take(2)
            .collect::<Vec<_>>(),
        [&"Pr_IBD2", &"MaxIBD2"]
    );
    // Un-analysed cross-family pairs carry the bare token, not a formatted -9.
    let ibs0 = s.read("king.ibs0");
    for row in rows(&ibs0) {
        assert_eq!(&row[row.len() - 2..], ["-9", "-9"]);
    }
}

/// `.ibs` is always created — header-only when no family has two members — while `.ibs0`
/// needs two families. `allsegs.txt` lists the segments the console counted.
#[test]
fn file_existence_and_the_one_family_lines() {
    let s = Scratch::new("ibs-exist");
    singleton_families(&s, 400);
    let body = run(&s, "s", &["--ibs"]);
    assert!(
        body.contains("Each family consists of one individual.\n"),
        "{body}"
    );
    assert!(!body.contains("Within-family IBS data saved"), "{body}");
    assert_eq!(s.read("king.ibs").lines().count(), 1, "header only");
    assert!(s.exists("king.ibs0"));
    assert_eq!(
        s.read("kingallsegs.txt"),
        concat!(
            "Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP\n",
            "1\t1\t1.000\t20.950\t19.950\t400\trs0\trs399\n",
        )
    );

    // One family, two members: `.ibs` gets rows and there is no between-family stage.
    let s = Scratch::new("ibs-one-family");
    let g = vec![pattern(0, 400), pattern(3, 400)];
    write_fileset(&s, "s", &[("F", "A"), ("F", "B")], &g);
    let body = run(&s, "s", &["--ibs"]);
    assert!(
        body.contains("Within-family IBS data saved in file king.ibs\n"),
        "{body}"
    );
    assert!(body.contains("There is only one family.\n"), "{body}");
    assert_eq!(s.read("king.ibs").lines().count(), 2);
    assert!(!s.exists("king.ibs0"), ".ibs0 needs two families");
}

/// `.ibs0` rows are square-tiled with a block of 8, not row-major.
///
/// With ten singleton families the two orders differ from the eighth row on: row-major
/// would put `(S00, S08)` immediately after `(S00, S07)`, while the tiled order finishes
/// the whole `0..8` triangle first.
#[test]
fn ibs0_rows_are_tiled_by_eight() {
    let s = Scratch::new("ibs-tiled");
    singleton_families(&s, 400);
    run(&s, "s", &["--ibs"]);
    let ibs0 = s.read("king.ibs0");
    let pairs: Vec<(String, String)> = rows(&ibs0)
        .iter()
        .map(|r| (r[1].to_string(), r[3].to_string()))
        .collect();
    assert_eq!(pairs.len(), 45);
    assert_eq!(pairs[7], ("S01".into(), "S02".into()));
    assert_eq!(pairs[27], ("S06".into(), "S07".into()));
    assert_eq!(pairs[28], ("S00".into(), "S08".into()));
}

/// Within a family the rows follow the ID comparator, not `.fam` order: `2 < 9 < 10`.
#[test]
fn ibs_rows_ignore_fam_order() {
    let s = Scratch::new("ibs-order");
    let g = vec![pattern(0, 400), pattern(2, 400), pattern(4, 400)];
    write_fileset(&s, "s", &[("F", "10"), ("F", "2"), ("F", "9")], &g);
    run(&s, "s", &["--ibs"]);
    let ibs = s.read("king.ibs");
    let pairs: Vec<(String, String)> = rows(&ibs)
        .iter()
        .map(|r| (r[1].to_string(), r[2].to_string()))
        .collect();
    assert_eq!(
        pairs,
        [
            ("2".to_string(), "9".to_string()),
            ("2".to_string(), "10".to_string()),
            ("9".to_string(), "10".to_string()),
        ]
    );
}

/// The `.ibs` column set, in the reference's own order, pipes included.
#[test]
fn ibs_header_layout() {
    let s = Scratch::new("ibs-header");
    let g = vec![pattern(0, 400), pattern(3, 400)];
    write_fileset(&s, "s", &[("F", "A"), ("F", "B")], &g);
    run(&s, "s", &["--ibs"]);
    assert_eq!(
        header(&s.read("king.ibs")),
        [
            "FID", "ID1", "ID2", "Z0", "Phi", "N_SNP", "N_IBS0", "N_IBS1", "N_IBS2", "NHetHet",
            "NHomHom", "N_Het1", "N_Het2", "IBS", "Dist", "HetConc", "Het2|1", "Het1|2", "HomConc",
            "Kinship"
        ]
    );
}
