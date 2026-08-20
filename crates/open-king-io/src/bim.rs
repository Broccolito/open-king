//! `.bim` — the PLINK variant table.
//!
//! Six whitespace-separated columns:
//!
//! ```text
//! CHR  ID  CM  BP  A1  A2
//! ```
//!
//! Fewer than six fields is an error; extra trailing fields are ignored, matching the
//! reference (see [`parse_bim`]'s tests).
//!
//! * `CHR` is kept verbatim (`1`, `chr1`, `X`, `MT`, or a contig name we do not
//!   recognise); [`crate::chrom_code`] maps it to PLINK's numeric code when asked, and
//!   [`crate::bed::king_chrom_code`] — the stricter mapping the reference itself uses —
//!   decides what the analysis filter keeps.
//! * `CM` is a genetic position in centimorgans. It is usually the literal `0`, and is
//!   frequently written as an integer, so it is parsed as a float either way.
//! * `BP` is signed: PLINK itself emits negative base-pair positions to mark variants it
//!   has been told to exclude from position-based operations, and such files must still
//!   load.
//! * `A1`/`A2` are kept as strings. They may be multi-character (indels) and either may be
//!   `0`, PLINK's missing-allele convention — a monomorphic variant gets `0` for the
//!   allele that was never observed. Such variants are *not* dropped; the reference binary
//!   counts them in its "autosome SNPs" total.

use std::path::Path;

use crate::fam::read_text;
use crate::{IoError, Result, Variant};

/// Number of columns every `.bim` line must have.
pub const N_FIELDS: usize = 6;

/// Read and parse a `.bim` file.
pub fn read_bim(path: &Path) -> Result<Vec<Variant>> {
    let text = read_text(path)?;
    parse_bim(&text, path)
}

/// Parse `.bim` text that has already been read into memory.
///
/// `path` is used only to label errors.
pub fn parse_bim(text: &str, path: &Path) -> Result<Vec<Variant>> {
    let mut variants = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let lineno = i + 1;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_ascii_whitespace().collect();
        // Too few fields is fatal (the reference reports `0 autosome SNPs` then a FATAL
        // ERROR); extra trailing fields are ignored, as the reference ignores them.
        if f.len() < N_FIELDS {
            return Err(IoError::Fields {
                path: path.to_path_buf(),
                line: lineno,
                expected: N_FIELDS,
                found: f.len(),
            });
        }
        let cm = parse_cm(f[2]).ok_or_else(|| IoError::Parse {
            path: path.to_path_buf(),
            line: lineno,
            field: "CM",
            value: f[2].to_string(),
        })?;
        let bp = f[3].parse::<i64>().map_err(|_| IoError::Parse {
            path: path.to_path_buf(),
            line: lineno,
            field: "BP",
            value: f[3].to_string(),
        })?;
        variants.push(Variant {
            chrom: f[0].to_string(),
            id: f[1].to_string(),
            cm,
            bp,
            a1: f[4].to_string(),
            a2: f[5].to_string(),
        });
    }
    Ok(variants)
}

/// Parse a centimorgan field, rejecting the non-finite values `f64::from_str` accepts.
///
/// `"0"`, `"0.0"`, `"1.25"` and `"1e-3"` are all fine; `"nan"` and `"inf"` are not, since
/// they would silently poison any later position arithmetic.
fn parse_cm(field: &str) -> Option<f64> {
    match field.parse::<f64>() {
        Ok(v) if v.is_finite() => Some(v),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::VariantFilter;
    use std::path::PathBuf;

    fn p() -> PathBuf {
        PathBuf::from("test.bim")
    }

    #[test]
    fn parses_six_columns_with_integer_and_fractional_cm() {
        let text = "\
1\trs1\t0\t1000\tA\tG
1 rs2 1.25 2000 C T
2\trs3\t2\t3000\tA\tG
";
        let v = parse_bim(text, &p()).unwrap();
        assert_eq!(v.len(), 3);
        assert_eq!(v[0].chrom, "1");
        assert_eq!(v[0].id, "rs1");
        assert_eq!(v[0].cm, 0.0);
        assert_eq!(v[0].bp, 1000);
        assert_eq!(v[0].a1, "A");
        assert_eq!(v[0].a2, "G");
        assert_eq!(v[1].cm, 1.25);
        // an integer cm is the same value as its float spelling
        assert_eq!(v[2].cm, 2.0);
    }

    #[test]
    fn negative_base_pair_positions_are_accepted() {
        let v = parse_bim("1 rs1 0 -12345 A G\n", &p()).unwrap();
        assert_eq!(v[0].bp, -12_345);
    }

    #[test]
    fn missing_and_multi_character_alleles_survive() {
        // a monomorphic site (A1 == "0") and an indel
        let text = "1 mono 0 100 0 G\n1 indel 0 200 GATTACA AT\n";
        let v = parse_bim(text, &p()).unwrap();
        assert_eq!((v[0].a1.as_str(), v[0].a2.as_str()), ("0", "G"));
        assert_eq!((v[1].a1.as_str(), v[1].a2.as_str()), ("GATTACA", "AT"));
        // and they are still ordinary autosomal variants
        assert!(v[0].is_autosome());
        assert!(v[1].is_autosome());
    }

    #[test]
    fn chromosome_text_is_kept_verbatim_but_still_codes() {
        let text = "chr1 a 0 1 A G\nX b 0 1 A G\nMT c 0 1 A G\nGL000226 d 0 1 A G\n";
        let v = parse_bim(text, &p()).unwrap();
        assert_eq!(v[0].chrom, "chr1");
        assert_eq!(v[0].chrom_code(), Some(1));
        assert!(v[0].is_autosome());
        assert_eq!(v[1].chrom, "X");
        assert_eq!(v[1].chrom_code(), Some(23));
        assert!(!v[1].is_autosome());
        assert_eq!(v[2].chrom_code(), Some(26));
        // an unplaced contig has no code and is not an autosome
        assert_eq!(v[3].chrom, "GL000226");
        assert_eq!(v[3].chrom_code(), None);
        assert!(!v[3].is_autosome());
        // ... which also means the default filter never keeps it
        assert!(!crate::bed::is_kept(&v[3], VariantFilter::Autosomes));
        assert!(crate::bed::is_kept(&v[3], VariantFilter::All));
    }

    #[test]
    fn too_few_fields_is_rejected() {
        match parse_bim("1 rs1 0 1000 A\n", &p()) {
            Err(IoError::Fields {
                line,
                expected,
                found,
                ..
            }) => {
                assert_eq!((line, expected, found), (1, 6, 5));
            }
            other => panic!("expected Fields error, got {other:?}"),
        }
    }

    /// The reference accepts extra trailing columns and ignores them.
    ///
    /// Probed directly: `king -b t.bed --bim seven.bim --ibs`, with a seventh field
    /// appended to every `.bim` line, printed the same `70 autosome SNPs` and wrote a
    /// `king.ibs` byte-identical to the six-column run. (A *five*-column `.bim` gives
    /// `0 autosome SNPs` and a FATAL ERROR there, which is why too-few stays an error
    /// here.)
    #[test]
    fn extra_trailing_columns_are_ignored_as_the_reference_ignores_them() {
        let six = parse_bim("1\trs1\t0\t1000\tA\tG\n2\trs2\t1.5\t20\tC\tT\n", &p()).unwrap();
        let seven =
            parse_bim("1\trs1\t0\t1000\tA\tG\tX\n2\trs2\t1.5\t20\tC\tT\tY\n", &p()).unwrap();
        assert_eq!(seven, six);
        assert_eq!(seven[1].a2, "T", "A2 must still come from column 6");
    }

    #[test]
    fn unparsable_numeric_fields_name_the_column_and_the_text() {
        match parse_bim("1 rs1 zzz 1000 A G\n", &p()) {
            Err(IoError::Parse {
                line, field, value, ..
            }) => {
                assert_eq!((line, field, value.as_str()), (1, "CM", "zzz"));
            }
            other => panic!("expected Parse error, got {other:?}"),
        }
        match parse_bim("1 rs1 0 1e9 A G\n", &p()) {
            // BP must be an integer; float syntax is not silently truncated
            Err(IoError::Parse { field, value, .. }) => {
                assert_eq!((field, value.as_str()), ("BP", "1e9"));
            }
            other => panic!("expected Parse error, got {other:?}"),
        }
    }

    #[test]
    fn non_finite_cm_is_rejected() {
        for bad in ["nan", "NaN", "inf", "-inf", "infinity"] {
            let line = format!("1 rs1 {bad} 1000 A G\n");
            match parse_bim(&line, &p()) {
                Err(IoError::Parse { field, .. }) => assert_eq!(field, "CM"),
                other => panic!("expected Parse error for {bad:?}, got {other:?}"),
            }
        }
        assert_eq!(parse_cm("1e-3"), Some(0.001));
        assert_eq!(parse_cm("-0.5"), Some(-0.5));
    }

    #[test]
    fn blank_lines_are_skipped_but_still_counted() {
        let v = parse_bim("\n1 rs1 0 1 A G\n\n1 rs2 0 2 A G\n", &p()).unwrap();
        assert_eq!(v.len(), 2);
        match parse_bim("1 rs1 0 1 A G\n\n1 rs2 0 zz A G\n", &p()) {
            Err(IoError::Parse { line, field, .. }) => assert_eq!((line, field), (3, "BP")),
            other => panic!("expected Parse error, got {other:?}"),
        }
    }

    #[test]
    fn missing_file_is_an_open_error() {
        let path = Path::new("/nonexistent-dir-open-king/none.bim");
        match read_bim(path) {
            Err(IoError::Open { path: p, .. }) => assert_eq!(p, path),
            other => panic!("expected Open error, got {other:?}"),
        }
    }
}
