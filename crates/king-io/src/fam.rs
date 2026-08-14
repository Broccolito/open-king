//! `.fam` — the PLINK sample table.
//!
//! Six whitespace-separated columns:
//!
//! ```text
//! FID  IID  PAT  MAT  SEX  PHENO
//! ```
//!
//! Every field except `SEX` is kept verbatim as text so that a fileset can be written
//! back out unchanged; `SEX` is normalised to `1`/`2`/`0` because the relatedness code
//! needs to compare it. Blank lines are ignored; a line with **fewer** than six fields is
//! an error. Extra trailing fields are read past and dropped rather than rejected: the
//! reference binary does the same, and a `.fam` given a seventh column analysed
//! byte-identically to the same file without it, so refusing it here would turn a file the
//! reference accepts into a FATAL ERROR.

use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::path::Path;

use crate::{IoError, Result, Sample};

/// Number of columns every `.fam` line must have.
pub const N_FIELDS: usize = 6;

/// Parse a `SEX` field.
///
/// PLINK writes `1` for male and `2` for female. Everything else — `0`, `-9`, `M`, `F`,
/// an empty field — is "unknown" and becomes `0`, matching PLINK's own leniency here.
pub fn parse_sex(field: &str) -> u8 {
    match field {
        "1" => 1,
        "2" => 2,
        _ => 0,
    }
}

/// Read and parse a `.fam` file.
///
/// Duplicate `(FID, IID)` pairs are *not* rejected here — use [`check_duplicates`] or
/// [`SampleIndex::build`] for that, so that callers which only want to inspect a malformed
/// file still can.
pub fn read_fam(path: &Path) -> Result<Vec<Sample>> {
    let text = read_text(path)?;
    parse_fam(&text, path)
}

/// Parse `.fam` text that has already been read into memory.
///
/// `path` is used only to label errors.
pub fn parse_fam(text: &str, path: &Path) -> Result<Vec<Sample>> {
    let mut samples = Vec::new();
    for (i, line) in text.lines().enumerate() {
        let lineno = i + 1;
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split_ascii_whitespace().collect();
        // Too few fields is fatal; extra trailing fields are ignored, because that is
        // what the reference does — a `.fam` with a seventh column analysed byte-identically
        // to the same file without it. See the tests for the probe.
        if f.len() < N_FIELDS {
            return Err(IoError::Fields {
                path: path.to_path_buf(),
                line: lineno,
                expected: N_FIELDS,
                found: f.len(),
            });
        }
        samples.push(Sample {
            fid: f[0].to_string(),
            iid: f[1].to_string(),
            pat: f[2].to_string(),
            mat: f[3].to_string(),
            sex: parse_sex(f[4]),
            pheno: f[5].to_string(),
        });
    }
    Ok(samples)
}

/// Return the first duplicated `(FID, IID)` pair, if any.
pub fn find_duplicate(samples: &[Sample]) -> Option<(&str, &str)> {
    let mut seen: HashMap<(&str, &str), usize> = HashMap::with_capacity(samples.len());
    for s in samples {
        if seen.insert((s.fid.as_str(), s.iid.as_str()), 0).is_some() {
            return Some((s.fid.as_str(), s.iid.as_str()));
        }
    }
    None
}

/// Error if any `(FID, IID)` pair occurs twice.
///
/// The reference binary treats this as fatal (`Family <FID>: Person <IID> is duplicated`),
/// so the loader does too.
pub fn check_duplicates(samples: &[Sample]) -> Result<()> {
    match find_duplicate(samples) {
        Some((fid, iid)) => Err(IoError::DuplicateSample {
            fid: fid.to_string(),
            iid: iid.to_string(),
        }),
        None => Ok(()),
    }
}

/// Lookup from `(FID, IID)` to position in the sample vector.
///
/// Borrows the samples so that lookups allocate nothing.
#[derive(Debug, Clone, Default)]
pub struct SampleIndex<'a> {
    by_key: HashMap<(&'a str, &'a str), usize>,
}

impl<'a> SampleIndex<'a> {
    /// Build the index, rejecting duplicate `(FID, IID)` keys.
    pub fn build(samples: &'a [Sample]) -> Result<Self> {
        let mut by_key = HashMap::with_capacity(samples.len());
        for (i, s) in samples.iter().enumerate() {
            if by_key.insert((s.fid.as_str(), s.iid.as_str()), i).is_some() {
                return Err(IoError::DuplicateSample {
                    fid: s.fid.clone(),
                    iid: s.iid.clone(),
                });
            }
        }
        Ok(Self { by_key })
    }

    /// Index of the sample with this `(FID, IID)`, or `None`.
    pub fn get(&self, fid: &str, iid: &str) -> Option<usize> {
        self.by_key.get(&(fid, iid)).copied()
    }

    /// Number of indexed samples.
    pub fn len(&self) -> usize {
        self.by_key.len()
    }

    /// Whether the index is empty.
    pub fn is_empty(&self) -> bool {
        self.by_key.is_empty()
    }
}

/// Open a file and slurp it, distinguishing "could not open" from "could not read".
pub(crate) fn read_text(path: &Path) -> Result<String> {
    let mut file = File::open(path).map_err(|source| IoError::Open {
        path: path.to_path_buf(),
        source,
    })?;
    let mut text = String::new();
    file.read_to_string(&mut text)
        .map_err(|source| IoError::Read {
            path: path.to_path_buf(),
            source,
        })?;
    Ok(text)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    fn p() -> PathBuf {
        PathBuf::from("test.fam")
    }

    const SIX: &str = "\
fam1 s1 0 0 1 -9
fam1 s2 s1 0 2 1.5
fam2\ts3\t0\t0\t0\t2
";

    #[test]
    fn parses_six_columns_and_keeps_text_verbatim() {
        let s = parse_fam(SIX, &p()).unwrap();
        assert_eq!(s.len(), 3);
        assert_eq!(
            s[1],
            Sample {
                fid: "fam1".into(),
                iid: "s2".into(),
                pat: "s1".into(),
                mat: "0".into(),
                sex: 2,
                // kept as text, not rounded through a float
                pheno: "1.5".into(),
            }
        );
        // tab-separated lines parse identically to space-separated ones
        assert_eq!(s[2].fid, "fam2");
        assert_eq!(s[2].iid, "s3");
        assert_eq!(s[2].pheno, "2");
        assert_eq!(s[0].sex, 1);
    }

    #[test]
    fn unrecognised_sex_becomes_zero() {
        for field in ["0", "-9", "M", "F", "3", "", "01", " 1"] {
            assert_eq!(parse_sex(field), 0, "sex field {field:?}");
        }
        assert_eq!(parse_sex("1"), 1);
        assert_eq!(parse_sex("2"), 2);
    }

    #[test]
    fn too_few_fields_is_rejected_with_the_offending_line() {
        let short = "fam1 s1 0 0 1\n";
        match parse_fam(short, &p()) {
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
        // ... and it is the *first* bad line that is named
        match parse_fam("f a 0 0 1 -9\nf b 0 0\n", &p()) {
            Err(IoError::Fields { line, found, .. }) => assert_eq!((line, found), (2, 4)),
            other => panic!("expected Fields error, got {other:?}"),
        }
    }

    /// The reference accepts extra trailing columns and ignores them.
    ///
    /// Probed directly: `king -b t.bed --fam seven.fam --ibs`, where `seven.fam` is the
    /// six-column `.fam` with a seventh field appended to every line, produced a
    /// `king.ibs` **byte-identical** to the six-column run — no warning, no error. Failing
    /// the load here would turn a file the reference analyses into a FATAL ERROR.
    #[test]
    fn extra_trailing_columns_are_ignored_as_the_reference_ignores_them() {
        let six = parse_fam(SIX, &p()).unwrap();
        let seven: String = SIX
            .lines()
            .map(|l| format!("{l} EXTRA\n"))
            .collect::<String>();
        let extra = parse_fam(&seven, &p()).unwrap();
        assert_eq!(extra, six, "the seventh column must not change anything");

        // a `.ped`-shaped line is read as its first six fields, exactly as the reference
        // reads it
        let ped = parse_fam("fam1 s1 0 0 1 -9 A A C T G G\n", &p()).unwrap();
        assert_eq!(ped.len(), 1);
        assert_eq!(ped[0].pheno, "-9");
        assert_eq!(ped[0].sex, 1);
    }

    #[test]
    fn blank_lines_are_skipped_but_still_counted() {
        let text = "fam1 s1 0 0 1 -9\n\n   \nfam1 s2 0 0\n";
        match parse_fam(text, &p()) {
            // the bad row is physical line 4, not "line 2 of the samples"
            Err(IoError::Fields { line, found, .. }) => assert_eq!((line, found), (4, 4)),
            other => panic!("expected Fields error, got {other:?}"),
        }
        let ok = parse_fam("\nfam1 s1 0 0 1 -9\n\n", &p()).unwrap();
        assert_eq!(ok.len(), 1);
    }

    #[test]
    fn crlf_line_endings_do_not_leak_into_fields() {
        let s = parse_fam("fam1 s1 0 0 1 -9\r\nfam1 s2 0 0 2 -9\r\n", &p()).unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!(s[0].pheno, "-9");
        assert_eq!(s[1].iid, "s2");
    }

    #[test]
    fn duplicate_fid_iid_pairs_are_detected() {
        let dup = parse_fam("f1 a 0 0 1 -9\nf2 b 0 0 1 -9\nf1 a 0 0 2 -9\n", &p()).unwrap();
        assert_eq!(find_duplicate(&dup), Some(("f1", "a")));
        match check_duplicates(&dup) {
            Err(IoError::DuplicateSample { fid, iid }) => {
                assert_eq!((fid.as_str(), iid.as_str()), ("f1", "a"));
            }
            other => panic!("expected DuplicateSample, got {other:?}"),
        }
        assert!(SampleIndex::build(&dup).is_err());
    }

    #[test]
    fn same_iid_in_different_families_is_not_a_duplicate() {
        let ok = parse_fam("f1 a 0 0 1 -9\nf2 a 0 0 1 -9\n", &p()).unwrap();
        assert_eq!(find_duplicate(&ok), None);
        check_duplicates(&ok).unwrap();
        let idx = SampleIndex::build(&ok).unwrap();
        assert_eq!(idx.get("f1", "a"), Some(0));
        assert_eq!(idx.get("f2", "a"), Some(1));
    }

    #[test]
    fn sample_index_maps_keys_to_fam_order() {
        let s = parse_fam(SIX, &p()).unwrap();
        let idx = SampleIndex::build(&s).unwrap();
        assert_eq!(idx.len(), 3);
        assert!(!idx.is_empty());
        assert_eq!(idx.get("fam1", "s1"), Some(0));
        assert_eq!(idx.get("fam1", "s2"), Some(1));
        assert_eq!(idx.get("fam2", "s3"), Some(2));
        // keys are the pair, not either half on its own
        assert_eq!(idx.get("fam2", "s1"), None);
        assert_eq!(idx.get("nope", "s1"), None);
        assert_eq!(idx.get("fam1", "nope"), None);
    }

    #[test]
    fn empty_index_reports_empty() {
        let none: Vec<Sample> = Vec::new();
        let idx = SampleIndex::build(&none).unwrap();
        assert!(idx.is_empty());
        assert_eq!(idx.len(), 0);
    }

    #[test]
    fn missing_file_is_an_open_error() {
        let path = Path::new("/nonexistent-dir-open-king/none.fam");
        match read_fam(path) {
            Err(IoError::Open { path: p, .. }) => assert_eq!(p, path),
            other => panic!("expected Open error, got {other:?}"),
        }
    }
}
