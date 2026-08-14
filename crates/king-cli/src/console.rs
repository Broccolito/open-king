//! Byte-for-byte reproduction of KING 2.3.2's console output.
//!
//! Every literal and every column here was derived from captured output of the reference
//! binary, not from its source. Each block documents the capture it reproduces.
//!
//! # Layout of the parameters block
//!
//! ```text
//!  <---------- 30 ---------->|   |<---- 15 ---->|
//!                 Binary File : nonexistent.bed (-bname)
//!
//!  <------------- 33 ------------>|   |
//!               Inference Parameter : --degree, --noscreen [-1717986816],
//!                                     --seglength, --minConc [0.80]
//!  <-------------- 36 --------------->|
//! ```
//!
//! * The "Binary File" row is `{label:>30} : {value:>15} (-bname)` — both fields are
//!   right-aligned and both grow if their content is longer.
//! * Every "Additional Options" row is `{section:>33} : ` followed by the section's
//!   options separated by `", "`.
//! * A row wraps when the next option would end past **column 78**: the comma is written
//!   on the current line, then a newline and [`CONTINUATION_INDENT`] spaces. The first
//!   option after `" : "` and the first on a continuation line are never wrap-tested, so
//!   an over-long option simply overflows.
//!
//! The wrap column was pinned by two captured rows: `Genetic Risk Score` ends at column
//! 78 without wrapping, and `Optional Input` wraps a `--covfile []` that would have
//! ended at column 82 — and `--sexchr [23]`, which would have ended at 79.

use std::fmt::Write as _;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::cli::{Kind, Opt, Options, GROUPS};

/// Width of the label column in the "Binary File" row.
const FILE_LABEL_WIDTH: usize = 30;
/// Width of the value column in the "Binary File" row.
const FILE_VALUE_WIDTH: usize = 15;
/// Width of the section-header column under "Additional Options".
const SECTION_WIDTH: usize = 33;
/// Column at which a wrapped option list resumes: [`SECTION_WIDTH`] plus `" : "`.
const CONTINUATION_INDENT: usize = SECTION_WIDTH + 3;
/// An option may end at this column but not past it.
const WRAP_COLUMN: usize = 78;

/// The version banner, without its newline.
pub const BANNER: &str = "KING 2.3.2 - (c) 2010-2023 Wei-Min Chen";

/// Body of the FATAL ERROR raised when no genotype file was given.
///
/// The closing paragraph mixes five tabs on one line with ten spaces on the next; that is
/// what the reference emits and it is reproduced verbatim.
pub const GENOTYPE_FILES_REQUIRED: &str = concat!(
    "Genotype files are required. e.g.,\n",
    "  king -b ex.bed --related\n",
    "\n",
    "Please check the reference paper Manichaikul et al. 2010 Bioinformatics,\n",
    "\t\t\t\t\tChen et al. 2024,\n",
    "          or the KING website at kingrelatedness.com",
);

/// The notice printed when a fileset was given but no analysis was requested.
///
/// The trailing `-- ` entries are in the reference's own string; it prints 24 slots and
/// two of them are empty.
pub const NO_ANALYSIS_NOTICE: &str = concat!(
    "Please specify one of the following 24 options: ",
    "--related --kinship --autoQC --mtscore --risk --ibs --homog --ibdseg --mds --pca ",
    "--cluster --build --bysample --bysnp --tdt --unrelated --duplicate --roh --grm ",
    "--gdt --pc -- --pcgdt --",
);

// ---------------------------------------------------------------------------
// Startup blocks
// ---------------------------------------------------------------------------

/// Banner plus the full "parameters in effect" block, ending with a blank line.
pub fn startup(opts: &Options) -> String {
    let mut s = String::with_capacity(2048);
    let _ = writeln!(s, "{BANNER}\n");
    s.push_str(&parameters_block(opts));
    s
}

/// The "The following parameters are in effect:" block, ending with a blank line.
pub fn parameters_block(opts: &Options) -> String {
    let mut s = String::with_capacity(2048);
    s.push_str("The following parameters are in effect:\n");
    let _ = writeln!(
        s,
        "{label:>lw$} : {value:>vw$} (-bname)",
        label = "Binary File",
        value = opts.bed,
        lw = FILE_LABEL_WIDTH,
        vw = FILE_VALUE_WIDTH,
    );
    s.push('\n');
    s.push_str("Additional Options\n");
    for (section, opts_in_section) in GROUPS {
        let _ = write!(s, "{section:>w$} : ", w = SECTION_WIDTH);
        let mut col = CONTINUATION_INDENT;
        for (i, &opt) in opts_in_section.iter().enumerate() {
            let item = render(opts, opt);
            if i > 0 {
                s.push(',');
                col += 1;
                if col + 1 + item.chars().count() > WRAP_COLUMN {
                    s.push('\n');
                    for _ in 0..CONTINUATION_INDENT {
                        s.push(' ');
                    }
                    col = CONTINUATION_INDENT;
                } else {
                    s.push(' ');
                    col += 1;
                }
            }
            s.push_str(&item);
            col += item.chars().count();
        }
        s.push('\n');
    }
    s.push('\n');
    s
}

/// One option as it appears in the block: `--kinship`, `--kinship [ON]`, `--prefix [king]`.
fn render(opts: &Options, opt: Opt) -> String {
    let name = opt.name();
    match opt.kind() {
        Kind::Flag => {
            if opts.flag(opt) {
                format!("--{name} [ON]")
            } else {
                format!("--{name}")
            }
        }
        // Integers are hidden entirely when zero, which is why the default `--degree`
        // shows nothing and `--degree 0` also shows nothing.
        Kind::Int => {
            let v = opts.int(opt);
            if v == 0 {
                format!("--{name}")
            } else {
                format!("--{name} [{v}]")
            }
        }
        // Doubles show when non-zero *or* when explicitly given: `--maxP 0` prints
        // `[0.00]` while an untouched `--maxP` prints nothing.
        Kind::Double => {
            let v = opts.double(opt);
            if v == 0.0 && !opts.was_given(opt) {
                format!("--{name}")
            } else {
                format!("--{name} [{}]", format_double(v))
            }
        }
        Kind::Str => format!("--{name} [{}]", opts.string(opt)),
    }
}

/// The reference's double format: `%.2f`, except that a non-zero value below 0.01
/// (negatives included) switches to `%.1e`.
///
/// Confirmed on the binary: `0.01` prints `0.01`, `0.0099999` prints `1.0e-02`,
/// `-2` prints `-2.0e+00`, `-0.0` prints `-0.00`, `1e20` prints
/// `100000000000000000000.00`.
fn format_double(v: f64) -> String {
    if v == 0.0 || v >= 0.01 {
        format!("{v:.2}")
    } else {
        format_scientific(v)
    }
}

/// C's `%.1e`: one fractional digit and a signed, at-least-two-digit exponent.
///
/// Rust writes `-5e-1` where C writes `-5.0e-01`, so the exponent is re-padded here.
fn format_scientific(v: f64) -> String {
    let raw = format!("{v:.1e}");
    let Some((mantissa, exp)) = raw.split_once('e') else {
        return raw;
    };
    let (sign, digits) = match exp.strip_prefix('-') {
        Some(d) => ('-', d),
        None => ('+', exp.strip_prefix('+').unwrap_or(exp)),
    };
    format!("{mantissa}e{sign}{digits:0>2}")
}

// ---------------------------------------------------------------------------
// Diagnostics
// ---------------------------------------------------------------------------

/// The WARNING block. Starts with a BEL, exactly as the reference does.
///
/// Returns an empty string when there is nothing to report.
pub fn warning_block(messages: &[String]) -> String {
    if messages.is_empty() {
        return String::new();
    }
    let mut s = String::new();
    s.push_str("\n\x07WARNING - \nProblems encountered parsing command line:\n\n");
    for m in messages {
        let _ = writeln!(s, "{m}");
    }
    s.push('\n');
    s
}

/// The FATAL ERROR block, including the blank line that always precedes and follows it.
pub fn fatal_block(message: &str) -> String {
    format!("\nFATAL ERROR - \n{message}\n\n")
}

/// The notice printed when a fileset was given but no analysis was selected.
pub fn no_analysis_block() -> String {
    format!("\n{NO_ANALYSIS_NOTICE}\n\n")
}

/// The notice printed when more than one analysis was selected.
///
/// Two spaces follow the colon, because the reference writes the colon and then one
/// space per entry: `separately:  --related --duplicate`.
pub fn separate_analyses_block(analyses: &[String]) -> String {
    let mut s = String::from("The following analyses will run separately: ");
    for a in analyses {
        let _ = write!(s, " {a}");
    }
    s.push_str("\n\n");
    s
}

/// `Genotype file <name> cannot be opened`.
pub fn genotype_file_unopenable(path: &str) -> String {
    format!("Genotype file {path} cannot be opened")
}

/// The lowest `--sexchr` the reference accepts; 1 and below are fatal.
pub const MIN_SEX_CHROMOSOME: i32 = 2;

/// The `--sexchr` value that means "human", and so prints no notice.
pub const HUMAN_SEX_CHROMOSOME: i32 = 23;

/// The line printed right after `KING starts at` whenever `--sexchr` is not 23.
pub fn non_human_notice(sexchr: i32) -> String {
    format!("Non-human samples are analyzed, with {sexchr} pairs of chromosomes\n")
}

/// `Sex chromosome <n> out of range.`
///
/// The trailing newline is in the reference's own message, so this fatal block ends with
/// three newlines where every other one ends with two. Verified with `od`.
pub fn sex_chromosome_out_of_range(sexchr: i32) -> String {
    format!("Sex chromosome {sexchr} out of range.\n")
}

/// Body of the FATAL ERROR raised by `--risk` without a `--model`.
pub const RISK_MODEL_REQUIRED: &str = "Please use --model <file> to specify a risk model.";

/// `--seglength` is accepted strictly between these two, in Mb. Both boundaries were
/// bisected against the binary: 0.99 and 10.01 are rejected, and the next representable
/// doubles either side of them are accepted.
pub const SEGLENGTH_MIN: f64 = 0.99;
/// See [`SEGLENGTH_MIN`].
pub const SEGLENGTH_MAX: f64 = 10.01;

/// `Minimum segment length is set as <n> bp`, for an in-range `--seglength`.
///
/// The reference follows it with a lone `.` and **no** newline, so whatever prints next
/// continues on that line — reproduced here.
pub fn segment_length_notice(seglength_mb: f64) -> String {
    let bp = (seglength_mb * 1_000_000.0) as i64;
    format!("Minimum segment length is set as {bp} bp\n.")
}

/// What the reference prints for an out-of-range `--seglength`.
pub const SEGLENGTH_OUT_OF_RANGE: &str = concat!(
    "KING supports minimum segment length from 1 to 10 Mb at the moment.\n",
    "Default seglength of 3Mb is used.\n",
);

/// What the reference prints when `--minConc` is outside 0..=1.
pub const MINCONC_OUT_OF_RANGE: &str = "minConc value is out of range and not specified.\n";

/// Body of the FATAL ERROR raised by a `--maxP` outside `0 < p < 2`.
///
/// The reference reports the tail that failed — `p/2` for the lower one, `1 - p/2` for
/// the upper — formatted with C's `%g`.
pub fn p_value_out_of_range(tail: f64) -> String {
    format!("p-value [{}] outside range in ninv()", format_g(tail))
}

/// C's `%g` with the default precision of 6: `%f` or `%e` depending on the exponent,
/// with trailing zeros and a bare decimal point removed.
fn format_g(v: f64) -> String {
    const PRECISION: i32 = 6;
    if v.is_nan() {
        return "nan".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "inf" } else { "-inf" }.to_string();
    }
    let sci = format!("{v:.*e}", (PRECISION - 1) as usize);
    let (mantissa, exponent) = sci.split_once('e').expect("Rust always writes an exponent");
    let x: i32 = exponent.parse().expect("exponent is an integer");
    if (-4..PRECISION).contains(&x) {
        let decimals = (PRECISION - 1 - x).max(0) as usize;
        strip_trailing_zeros(&format!("{v:.decimals$}"))
    } else {
        let m = strip_trailing_zeros(mantissa);
        let sign = if x < 0 { '-' } else { '+' };
        format!("{m}e{sign}{:02}", x.abs())
    }
}

fn strip_trailing_zeros(s: &str) -> String {
    if !s.contains('.') {
        return s.to_string();
    }
    s.trim_end_matches('0').trim_end_matches('.').to_string()
}

// ---------------------------------------------------------------------------
// Run-time lines (for the analysis engines)
// ---------------------------------------------------------------------------

/// `KING starts at Thu Aug 13 17:29:25 2026`, with its newline.
pub fn king_starts_at(unix_secs: i64) -> String {
    format!("KING starts at {}\n", ctime(unix_secs))
}

/// `KING ends at Thu Aug 13 17:29:25 2026`, with its newline.
pub fn king_ends_at(unix_secs: i64) -> String {
    format!("KING ends at {}\n", ctime(unix_secs))
}

/// A phase's `... ends at <time>` line, right-aligned under the matching `starts at`.
///
/// The reference pads with spaces so the timestamps line up, e.g. the 41 spaces that
/// precede `ends at` under `Relationship inference across families starts at`.
pub fn ends_at(indent: usize, unix_secs: i64) -> String {
    format!("{:indent$}ends at {}\n", "", ctime(unix_secs))
}

/// `Loading genotype data in PLINK binary format...`
pub fn loading_plink_binary() -> String {
    "Loading genotype data in PLINK binary format...\n".to_string()
}

/// `Read in PLINK fam file <path>...`
pub fn read_fam(path: &str) -> String {
    format!("Read in PLINK fam file {path}...\n")
}

/// `  PLINK pedigrees loaded: <n> samples`
pub fn pedigrees_loaded(samples: usize) -> String {
    format!("  PLINK pedigrees loaded: {samples} samples\n")
}

/// `Read in PLINK bim file <path>...`
pub fn read_bim(path: &str) -> String {
    format!("Read in PLINK bim file {path}...\n")
}

/// `  Genotype data consist of <n> autosome SNPs`
pub fn autosome_snps(snps: usize) -> String {
    format!("  Genotype data consist of {snps} autosome SNPs\n")
}

/// `  PLINK maps loaded: <n> SNPs`
pub fn maps_loaded(snps: usize) -> String {
    format!("  PLINK maps loaded: {snps} SNPs\n")
}

/// `Read in PLINK bed file <path>...`
pub fn read_bed(path: &str) -> String {
    format!("Read in PLINK bed file {path}...\n")
}

/// One `\r`-terminated progress tick, e.g. `13%\r`.
///
/// The percentage is `done / total` rounded half-up, matching the captured sequence
/// `0% 6% 13% 19% 25% 31% 38% 44% 50% 56% 63% 69% 75% 81% 88% 94%` for 16 steps.
pub fn progress(done: usize, total: usize) -> String {
    if total == 0 {
        return "0%\r".to_string();
    }
    let pct = (200 * done + total) / (2 * total);
    format!("{pct}%\r")
}

/// `  PLINK binary genotypes loaded.`
pub fn binary_genotypes_loaded() -> String {
    "  PLINK binary genotypes loaded.\n".to_string()
}

/// `  KING format genotype data successfully converted.`
pub fn genotypes_converted() -> String {
    "  KING format genotype data successfully converted.\n".to_string()
}

/// `Autosome genotypes stored in <w> words for each of <n> individuals.`
pub fn autosome_words(words: usize, individuals: usize) -> String {
    format!("Autosome genotypes stored in {words} words for each of {individuals} individuals.\n")
}

/// The tab-indented `Options in effect:` block, ending with a blank line.
pub fn options_in_effect(flags: &[String]) -> String {
    let mut s = String::from("Options in effect:\n");
    for f in flags {
        let _ = writeln!(s, "\t{f}");
    }
    s.push('\n');
    s
}

/// `<n> CPU cores are used.`
pub fn cpu_cores(n: usize) -> String {
    format!("{n} CPU cores are used.\n")
}

/// `  <n> CPU cores are used...` — the indented, ellipsised variant.
pub fn cpu_cores_indented(n: usize) -> String {
    format!("  {n} CPU cores are used...\n")
}

/// One row of the relationship summary table.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RelationshipCounts {
    pub mz: u64,
    pub po: u64,
    pub fs: u64,
    pub second: u64,
    pub third: u64,
    pub other: u64,
}

impl RelationshipCounts {
    /// "Total relatives" counts every class except OTHER.
    fn relatives(&self) -> u64 {
        self.mz + self.po + self.fs + self.second + self.third
    }
}

/// The relationship summary table, with the blank line before and after it that the
/// reference always emits.
///
/// ```text
///
/// Relationship summary (total relatives: 7 by pedigree, 9 by inference)
///   Source  MZ  PO  FS  2nd 3rd OTHER
///   ===========================================================
///   Pedigree    0   6   1   0   0   5
///   Inference   0   6   1   1   1   3
///
/// ```
/// (the columns are separated by single tabs).
pub fn relationship_summary(pedigree: RelationshipCounts, inference: RelationshipCounts) -> String {
    let row = |name: &str, c: RelationshipCounts| {
        format!(
            "  {name}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            c.mz, c.po, c.fs, c.second, c.third, c.other
        )
    };
    let mut s = String::new();
    let _ = write!(
        s,
        "\nRelationship summary (total relatives: {} by pedigree, {} by inference)\n",
        pedigree.relatives(),
        inference.relatives()
    );
    s.push_str("  Source\tMZ\tPO\tFS\t2nd\t3rd\tOTHER\n");
    s.push_str("  ===========================================================\n");
    s.push_str(&row("Pedigree", pedigree));
    s.push_str(&row("Inference", inference));
    s.push('\n');
    s
}

// ---------------------------------------------------------------------------
// Time
// ---------------------------------------------------------------------------

/// Seconds since the Unix epoch, right now.
pub fn now() -> i64 {
    SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map(|d| d.as_secs() as i64)
        .unwrap_or(0)
}

/// Current local time, ready for [`king_starts_at`] and friends.
pub fn now_local() -> i64 {
    let t = now();
    t + tz::local_offset_secs(t)
}

/// C's `ctime()` without its trailing newline: `Thu Aug 13 17:29:25 2026`.
///
/// `asctime` formats the day of month with `%3d`, so single-digit days are space padded
/// and produce two spaces after the month: `Thu Aug  3 09:00:00 2026`.
pub fn ctime(unix_secs: i64) -> String {
    const DAYS: [&str; 7] = ["Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat"];
    const MONTHS: [&str; 12] = [
        "Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec",
    ];
    let days = unix_secs.div_euclid(86_400);
    let secs = unix_secs.rem_euclid(86_400);
    let (year, month, day) = civil_from_days(days);
    let weekday = (days + 4).rem_euclid(7) as usize;
    format!(
        "{} {} {:>2} {:02}:{:02}:{:02} {}",
        DAYS[weekday],
        MONTHS[(month - 1) as usize],
        day,
        secs / 3600,
        (secs / 60) % 60,
        secs % 60,
        year
    )
}

/// Proleptic-Gregorian date for a day count relative to 1970-01-01 (Hinnant's algorithm).
fn civil_from_days(z: i64) -> (i64, u32, u32) {
    let z = z + 719_468;
    let era = if z >= 0 { z } else { z - 146_096 } / 146_097;
    let doe = (z - era * 146_097) as u64; // [0, 146096]
    let yoe = (doe - doe / 1460 + doe / 36_524 - doe / 146_096) / 365; // [0, 399]
    let y = yoe as i64 + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100); // [0, 365]
    let mp = (5 * doy + 2) / 153; // [0, 11]
    let d = (doy - (153 * mp + 2) / 5 + 1) as u32; // [1, 31]
    let m = if mp < 10 { mp + 3 } else { mp - 9 } as u32; // [1, 12]
    (if m <= 2 { y + 1 } else { y }, m, d)
}

/// Local UTC offset, read from the system's zoneinfo database.
///
/// Deliberately dependency-free. On a platform without a TZif database (Windows) the
/// offset falls back to 0, i.e. timestamps print in UTC; a platform seam for that lives
/// here and nowhere else.
mod tz {
    use std::path::PathBuf;

    /// Offset in seconds to add to a UTC timestamp to get local time at that instant.
    pub fn local_offset_secs(unix_secs: i64) -> i64 {
        tzfile_path()
            .and_then(|p| std::fs::read(p).ok())
            .and_then(|d| offset_from_tzif(&d, unix_secs))
            .unwrap_or(0)
    }

    fn tzfile_path() -> Option<PathBuf> {
        if let Ok(tz) = std::env::var("TZ") {
            let tz = tz.strip_prefix(':').unwrap_or(&tz);
            if !tz.is_empty() && !tz.starts_with('/') {
                let p = PathBuf::from("/usr/share/zoneinfo").join(tz);
                if p.is_file() {
                    return Some(p);
                }
            } else if tz.starts_with('/') && PathBuf::from(tz).is_file() {
                return Some(PathBuf::from(tz));
            }
        }
        let p = PathBuf::from("/etc/localtime");
        p.is_file().then_some(p)
    }

    fn be32(b: &[u8]) -> i64 {
        i32::from_be_bytes([b[0], b[1], b[2], b[3]]) as i64
    }
    fn be64(b: &[u8]) -> i64 {
        i64::from_be_bytes([b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7]])
    }

    /// Parse a TZif file far enough to answer "what was the offset at `t`".
    fn offset_from_tzif(data: &[u8], t: i64) -> Option<i64> {
        let mut pos = 0usize;
        let mut time_width = 4usize;
        let mut header = parse_header(data, pos)?;
        if data.get(4).copied().unwrap_or(0) >= b'2' {
            // Skip the 32-bit block and re-read the 64-bit one behind it.
            pos = 44 + block_len(&header, 4);
            header = parse_header(data, pos)?;
            time_width = 8;
        }
        let body = pos + 44;
        let times = body;
        let idx = times + header.timecnt * time_width;
        let types = idx + header.timecnt;
        if data.len() < types + header.typecnt * 6 {
            return None;
        }

        let read_time = |i: usize| -> i64 {
            let at = times + i * time_width;
            if time_width == 8 {
                be64(&data[at..at + 8])
            } else {
                be32(&data[at..at + 4])
            }
        };
        let gmtoff = |type_index: usize| -> i64 { be32(&data[types + type_index * 6..]) };

        // Last transition at or before `t`.
        let mut chosen: Option<usize> = None;
        for i in 0..header.timecnt {
            if read_time(i) <= t {
                chosen = Some(i);
            } else {
                break;
            }
        }
        match chosen {
            Some(i) => Some(gmtoff(data[idx + i] as usize)),
            // Before the first transition: first non-DST type, else type 0.
            None => {
                let first_std = (0..header.typecnt).find(|&k| data[types + k * 6 + 4] == 0);
                (header.typecnt > 0).then(|| gmtoff(first_std.unwrap_or(0)))
            }
        }
    }

    struct Header {
        isutcnt: usize,
        isstdcnt: usize,
        leapcnt: usize,
        timecnt: usize,
        typecnt: usize,
        charcnt: usize,
    }

    fn parse_header(data: &[u8], pos: usize) -> Option<Header> {
        if data.len() < pos + 44 || &data[pos..pos + 4] != b"TZif" {
            return None;
        }
        let n = |k: usize| -> usize { be32(&data[pos + 20 + k * 4..]) as usize };
        Some(Header {
            isutcnt: n(0),
            isstdcnt: n(1),
            leapcnt: n(2),
            timecnt: n(3),
            typecnt: n(4),
            charcnt: n(5),
        })
    }

    fn block_len(h: &Header, time_width: usize) -> usize {
        h.timecnt * (time_width + 1)
            + h.typecnt * 6
            + h.charcnt
            + h.leapcnt * (time_width + 4)
            + h.isstdcnt
            + h.isutcnt
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ctime_matches_c() {
        // Cross-checked against C's asctime() for each timestamp.
        assert_eq!(ctime(1_786_642_165), "Thu Aug 13 17:29:25 2026");
        // The %3d quirk: a single-digit day is space padded.
        assert_eq!(ctime(1_785_747_907), "Mon Aug  3 09:05:07 2026");
        assert_eq!(ctime(0), "Thu Jan  1 00:00:00 1970");
        assert_eq!(ctime(951_868_799), "Tue Feb 29 23:59:59 2000");
        assert_eq!(ctime(946_684_799), "Fri Dec 31 23:59:59 1999");
    }

    #[test]
    fn local_offset_is_a_whole_number_of_minutes() {
        let off = tz::local_offset_secs(now());
        assert_eq!(off % 60, 0);
        assert!(off.abs() <= 14 * 3600);
    }

    #[test]
    fn scientific_matches_c() {
        assert_eq!(format_double(-1.5), "-1.5e+00");
        assert_eq!(format_double(-2.0), "-2.0e+00");
        assert_eq!(format_double(-0.5), "-5.0e-01");
        assert_eq!(format_double(0.001), "1.0e-03");
        assert_eq!(format_double(1e-7), "1.0e-07");
        assert_eq!(format_double(0.005), "5.0e-03");
        assert_eq!(format_double(-1e9), "-1.0e+09");
        assert_eq!(format_double(0.0099999), "1.0e-02");
        assert_eq!(format_double(1e-30), "1.0e-30");
    }

    #[test]
    fn fixed_matches_c() {
        assert_eq!(format_double(0.0), "0.00");
        assert_eq!(format_double(-0.0), "-0.00");
        assert_eq!(format_double(0.01), "0.01");
        assert_eq!(format_double(0.25), "0.25");
        assert_eq!(format_double(0.8), "0.80");
        assert_eq!(format_double(2.5), "2.50");
        assert_eq!(format_double(3.0), "3.00");
        assert_eq!(format_double(9.995), "9.99");
        assert_eq!(format_double(0.09999), "0.10");
        assert_eq!(format_double(1234567.0), "1234567.00");
        assert_eq!(format_double(1e20), "100000000000000000000.00");
    }

    #[test]
    fn progress_ticks_match_capture() {
        let ticks: Vec<String> = (0..16).map(|i| progress(i, 16)).collect();
        let expected = [
            "0%", "6%", "13%", "19%", "25%", "31%", "38%", "44%", "50%", "56%", "63%", "69%",
            "75%", "81%", "88%", "94%",
        ];
        for (got, want) in ticks.iter().zip(expected) {
            assert_eq!(got, &format!("{want}\r"));
        }
    }
}
