//! Numeric parity of the IBD-segment engine against captured reference output.
//!
//! Gated on `KING_GOLDEN`, which must point at `tests/parity/golden` in this repo:
//!
//! ```text
//! KING_GOLDEN=tests/parity/golden cargo test -p king-core --test ibdseg_parity -- --nocapture
//! ```
//!
//! Two checks per dataset:
//!
//! * `<ds>__ibdseg/kingallsegs.txt` must be reproduced **byte for byte** — this is the
//!   IBD denominator, so any error here rescales every number downstream;
//! * `<ds>__ibdseg/king.seg` is compared row by row, and the report prints how many rows
//!   agree at the printed four decimals and how far the worst one is off. That file is
//!   the honest measure of where the segment caller stands.
//!
//! **Two row counts are printed and both matter.** `king.seg` grades four fields per row:
//!
//! ```text
//! IBD1Seg  IBD2Seg  PropIBD  InfType
//! ```
//!
//! `est=` counts rows whose two **estimate** columns both round to the reference's, and
//! `row=` counts rows where all four fields do — the latter is what the parity harness
//! needs, since it diffs the file byte for byte. Reporting only `est=` would overstate the
//! engine, and for two generations of the caller the two counts were far apart: `PropIBD`
//! was computed as `IBD2Seg + IBD1Seg/2` at full precision, and 176 of 982 rows had both
//! estimates exact and still printed a different `PropIBD`.
//!
//! They now agree, because that was never a caller problem: `<prefix>.seg` computes
//! `PropIBD` from the two columns *after* rounding them to the four decimals it prints
//! (`king_core::ibdseg::seg_prop_ibd`), while `<prefix>.kin` uses full precision — the
//! reference prints two different `PropIBD` values for the same pair in the same run.
//! With each writer given its own rule, `row=` and `est=` are both 982 at the default
//! floor and `PropIBD` contributes no error at any floor. `docs/PARITY.md` §4.3 and §4.4.

#![forbid(unsafe_code)]

use std::fs;
use std::path::{Path, PathBuf};

use king_core::ibdseg;
use king_io::{bed, Variant, VariantFilter};

const DATASETS: &[&str] = &[
    "nuclear",
    "threegen",
    "multifam",
    "dups",
    "missing",
    "monomorphic",
    "sexchr",
    "unrelated",
    "admixed",
    "bigish",
];

fn golden() -> Option<PathBuf> {
    let raw = std::env::var("KING_GOLDEN").ok()?;
    let p = PathBuf::from(&raw);
    let p = if p.is_absolute() {
        p
    } else {
        Path::new(env!("CARGO_MANIFEST_DIR")).join("../..").join(p)
    };
    p.is_dir().then_some(p)
}

/// One analysis array: chromosome codes, positions, and `.bim` row indices.
#[derive(Default)]
struct Arr {
    chr: Vec<i64>,
    pos: Vec<i64>,
    idx: Vec<usize>,
}

/// Split a `.bim` into the autosomal and X analysis arrays, as the engine sees them.
fn arrays(variants: &[Variant]) -> (Arr, Arr) {
    let (mut ac, mut ap, mut ai) = (vec![], vec![], vec![]);
    let (mut xc, mut xp, mut xi) = (vec![], vec![], vec![]);
    for (i, v) in variants.iter().enumerate() {
        let code: i64 = v.chrom.parse().unwrap_or(0);
        if (1..=22).contains(&code) || code == 25 {
            ac.push(code);
            ap.push(v.bp);
            ai.push(i);
        } else if code == 23 {
            xc.push(code);
            xp.push(v.bp);
            xi.push(i);
        }
    }
    (
        Arr {
            chr: ac,
            pos: ap,
            idx: ai,
        },
        Arr {
            chr: xc,
            pos: xp,
            idx: xi,
        },
    )
}

#[test]
fn allsegs_is_byte_identical() {
    let Some(root) = golden() else {
        eprintln!("KING_GOLDEN not set; skipping");
        return;
    };
    let mut checked = 0;
    for ds in DATASETS {
        let want_path = root.join(format!("ibdseg/{ds}__ibdseg/kingallsegs.txt"));
        let Ok(want) = fs::read_to_string(&want_path) else {
            continue;
        };
        // The genotype files are deliberately never committed (see .gitignore); a fresh
        // checkout has the golden text but not the corpus. Skip rather than panic, and
        // let the `checked` assertion below report the shortfall.
        let Ok(variants) = king_io::bim::read_bim(&root.join(format!("{ds}.bim"))) else {
            continue;
        };
        let (au, x) = arrays(&variants);
        let auto = ibdseg::usable_segments(&au.chr, &au.pos);
        let xseg = ibdseg::usable_segments(&x.chr, &x.pos);

        let mut got =
            String::from("Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP\n");
        let mut n = 0;
        for (segs, arr) in [(&auto, &au), (&xseg, &x)] {
            for s in segs {
                n += 1;
                let (start, stop) = (arr.pos[s.lo] as f64 / 1e6, arr.pos[s.hi] as f64 / 1e6);
                got.push_str(&format!(
                    "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{}\t{}\n",
                    n,
                    s.chr,
                    start,
                    stop,
                    stop - start,
                    s.n_snp(),
                    variants[arr.idx[s.lo]].id,
                    variants[arr.idx[s.hi]].id
                ));
            }
        }
        assert_eq!(got, want, "{ds}: allsegs.txt differs");
        checked += 1;
    }
    assert!(
        checked >= 10,
        "only {checked} of {} datasets checked -- generate the corpus first:\n  \
         python3 tests/parity/generate_corpus.py --outdir tests/parity/golden",
        DATASETS.len()
    );
}

#[test]
fn seg_rows_report() {
    let Some(root) = golden() else {
        eprintln!("KING_GOLDEN not set; skipping");
        return;
    };
    let mut datasets_measured = 0usize;
    let (mut tot_rows, mut tot_exact, mut tot_extra, mut tot_missing) = (0, 0, 0, 0);
    let (mut tot_types, mut tot_dsum, mut tot_full) = (0usize, 0.0f64, 0usize);
    let mut worst: f64 = 0.0;
    for ds in DATASETS {
        let want_path = root.join(format!("ibdseg/{ds}__ibdseg/king.seg"));
        let Ok(want) = fs::read_to_string(&want_path) else {
            continue;
        };
        // Same as above: no corpus, no check. Never a panic.
        let Ok(fs_) = bed::read_fileset(
            &root.join(format!("{ds}.bed")),
            VariantFilter::Autosomes,
            None,
            None,
        ) else {
            continue;
        };
        datasets_measured += 1;
        let (au, _) = arrays(&fs_.variants);
        let ap = au.pos;
        let auto = ibdseg::usable_segments(&au.chr, &ap);
        let denom = ibdseg::denominator(&auto, &ap);

        let mut got = std::collections::BTreeMap::new();
        for i in 0..fs_.samples.len() {
            for j in i + 1..fs_.samples.len() {
                let s = ibdseg::pair_segments(
                    &fs_.genotypes,
                    &ap,
                    &auto,
                    i,
                    j,
                    ibdseg::DEFAULT_SEGLENGTH_BP,
                );
                if !s.reported() {
                    continue;
                }
                let (p1, p2) = (s.ibd1_seg(denom), s.ibd2_seg(denom));
                got.insert(
                    (fs_.samples[i].iid.clone(), fs_.samples[j].iid.clone()),
                    // `PropIBD` as `<prefix>.seg` prints it, not the `.kin` value —
                    // this grades the file the parity harness diffs.
                    (p1, p2, ibdseg::seg_prop_ibd(p1, p2)),
                );
            }
        }

        let mut want_rows = std::collections::BTreeMap::new();
        for line in want.lines().skip(1) {
            let f: Vec<&str> = line.split('\t').collect();
            want_rows.insert(
                (f[1].to_string(), f[3].to_string()),
                (
                    f[4].parse::<f64>().unwrap(),
                    f[5].parse::<f64>().unwrap(),
                    f[6].parse::<f64>().unwrap(),
                ),
            );
        }

        let (mut exact, mut dmax, mut dsum, mut types) = (0usize, 0.0f64, 0.0f64, 0usize);
        let mut full = 0usize;
        for (k, w) in &want_rows {
            if let Some(g) = got.get(k) {
                let same = format!("{:.4}", g.0) == format!("{:.4}", w.0)
                    && format!("{:.4}", g.1) == format!("{:.4}", w.1);
                if same {
                    exact += 1;
                }
                // `InfType` reads the full-precision value on both sides.
                let same_type =
                    ibdseg::inf_type(g.0, g.1, g.1 + g.0 / 2.0) == ibdseg::inf_type(w.0, w.1, w.2);
                if same_type {
                    types += 1;
                }
                // The harness diffs the file, so a row only counts when every printed
                // field matches — `PropIBD` included. See the module doc.
                if same && same_type && format!("{:.4}", g.2) == format!("{:.4}", w.2) {
                    full += 1;
                }
                dmax = dmax.max((g.2 - w.2).abs());
                dsum += (g.2 - w.2).abs();
            }
        }
        tot_types += types;
        tot_dsum += dsum;
        tot_full += full;
        let extra = got.keys().filter(|k| !want_rows.contains_key(*k)).count();
        let missing = want_rows.keys().filter(|k| !got.contains_key(*k)).count();
        eprintln!(
            "{ds:<12} gold={:<5} ours={:<5} row={:<5} est={:<5} infType={:<5} missing={:<4} extra={:<4} maxdPropIBD={:.6}",
            want_rows.len(),
            got.len(),
            full,
            exact,
            types,
            missing,
            extra,
            dmax
        );
        tot_rows += want_rows.len();
        tot_exact += exact;
        tot_extra += extra;
        tot_missing += missing;
        worst = worst.max(dmax);
    }
    eprintln!(
        "TOTAL gold={tot_rows} row={tot_full} est={tot_exact} infType={tot_types} missing={tot_missing} extra={tot_extra} meandPropIBD={:.6} worst={worst:.4}",
        tot_dsum / tot_rows as f64
    );

    // A report that measured nothing must not read as a pass: the corpus is
    // gitignored, so an absent one would otherwise make this test silently green.
    assert!(
        datasets_measured >= 10,
        "only {datasets_measured} of {} datasets measured -- generate the corpus first:\n  \\
         python3 tests/parity/generate_corpus.py --outdir tests/parity/golden",
        DATASETS.len()
    );
}
