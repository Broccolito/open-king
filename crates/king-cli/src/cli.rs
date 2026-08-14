//! Hand-written command-line parser reproducing KING 2.3.2's quirks.
//!
//! Deliberately not `clap`: the reference binary's parser has behaviour no general
//! argument library reproduces — case-insensitive **prefix** matching on long names,
//! value tokens that are only consumed when they *look* like a number, integer options
//! that toggle when given no value, and a shared-storage bug between `--noscreen` and
//! `--minConc`. Every rule below was derived by black-box probing of the reference
//! binary; the probes are quoted in the doc comment of each rule.
//!
//! Nothing here reads or is derived from KING's source code.

use std::collections::BTreeSet;

/// How an option carries its value, which also decides how it renders in the
/// "parameters in effect" block.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    /// `--kinship` — no value; renders `--kinship [ON]` when set.
    Flag,
    /// `--degree 2` — renders `--degree [2]`, hidden entirely when the value is 0.
    Int,
    /// `--minConc 0.9` — renders `--minConc [0.90]`.
    Double,
    /// `--prefix out` — renders `--prefix [out]`, always shown, even when empty.
    Str,
}

/// Every option the reference binary lists in its parameters-in-effect block.
///
/// The discriminants are storage slots, and the declaration order is the order the
/// block prints them in, so do not reorder without also reordering [`GROUPS`].
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Debug)]
#[repr(usize)]
pub enum Opt {
    Related = 0,
    Duplicate,
    Kinship,
    Ibdseg,
    Ibs,
    MakeGrm,
    Degree,
    Noscreen,
    Seglength,
    MinConc,
    Unrelated,
    Cluster,
    Build,
    Bysample,
    BySnp,
    Roh,
    AutoQc,
    CallrateN,
    CallrateM,
    Pca,
    Mds,
    Projection,
    Pcs,
    Lmm,
    Tdt,
    Gdt,
    Trait,
    Covariate,
    MaxP,
    Invnorm,
    Risk,
    Model,
    Prevalence,
    Noflip,
    Cpus,
    Fam,
    Bim,
    Phefile,
    Covfile,
    Prunedsnp,
    Sexchr,
    Rplot,
    Pngplot,
    Plink,
    Prefix,
    Rpath,
}

impl Opt {
    /// Number of options in the table.
    pub const COUNT: usize = 46;

    /// The option's spelling, without the leading `--`. Case matters for display only;
    /// matching is case-insensitive.
    pub fn name(self) -> &'static str {
        SPECS[self as usize].0
    }

    /// How the option carries its value.
    pub fn kind(self) -> Kind {
        SPECS[self as usize].1
    }
}

/// `(spelling, kind)` indexed by `Opt as usize`.
static SPECS: [(&str, Kind); Opt::COUNT] = [
    ("related", Kind::Flag),
    ("duplicate", Kind::Flag),
    ("kinship", Kind::Flag),
    ("ibdseg", Kind::Flag),
    ("ibs", Kind::Flag),
    ("makeGRM", Kind::Flag),
    ("degree", Kind::Int),
    ("noscreen", Kind::Int),
    ("seglength", Kind::Double),
    ("minConc", Kind::Double),
    ("unrelated", Kind::Flag),
    ("cluster", Kind::Flag),
    ("build", Kind::Flag),
    ("bysample", Kind::Flag),
    ("bySNP", Kind::Flag),
    ("roh", Kind::Flag),
    ("autoQC", Kind::Flag),
    ("callrateN", Kind::Double),
    ("callrateM", Kind::Double),
    ("pca", Kind::Flag),
    ("mds", Kind::Flag),
    ("projection", Kind::Int),
    ("pcs", Kind::Int),
    ("lmm", Kind::Flag),
    ("tdt", Kind::Flag),
    ("gdt", Kind::Flag),
    ("trait", Kind::Str),
    ("covariate", Kind::Str),
    ("maxP", Kind::Double),
    ("invnorm", Kind::Flag),
    ("risk", Kind::Flag),
    ("model", Kind::Str),
    ("prevalence", Kind::Double),
    ("noflip", Kind::Flag),
    ("cpus", Kind::Int),
    ("fam", Kind::Str),
    ("bim", Kind::Str),
    ("phefile", Kind::Str),
    ("covfile", Kind::Str),
    ("prunedsnp", Kind::Str),
    ("sexchr", Kind::Int),
    ("rplot", Kind::Flag),
    ("pngplot", Kind::Flag),
    ("plink", Kind::Flag),
    ("prefix", Kind::Str),
    ("rpath", Kind::Str),
];

/// The "Additional Options" sections, in print order, with their right-aligned headers.
pub static GROUPS: &[(&str, &[Opt])] = &[
    ("Close Relative Inference", &[Opt::Related, Opt::Duplicate]),
    (
        "Pairwise Relatedness Inference",
        &[Opt::Kinship, Opt::Ibdseg, Opt::Ibs, Opt::MakeGrm],
    ),
    (
        "Inference Parameter",
        &[Opt::Degree, Opt::Noscreen, Opt::Seglength, Opt::MinConc],
    ),
    (
        "Relationship Application",
        &[Opt::Unrelated, Opt::Cluster, Opt::Build],
    ),
    (
        "QC Report",
        &[Opt::Bysample, Opt::BySnp, Opt::Roh, Opt::AutoQc],
    ),
    ("QC Parameter", &[Opt::CallrateN, Opt::CallrateM]),
    ("Population Structure", &[Opt::Pca, Opt::Mds]),
    ("Structure Parameter", &[Opt::Projection, Opt::Pcs]),
    ("Quantitative Trait GWAS", &[Opt::Lmm]),
    ("Binary Trait GWAS", &[Opt::Tdt, Opt::Gdt]),
    (
        "Association Model",
        &[Opt::Trait, Opt::Covariate, Opt::MaxP],
    ),
    ("Association Method Parameter", &[Opt::Invnorm]),
    (
        "Genetic Risk Score",
        &[Opt::Risk, Opt::Model, Opt::Prevalence, Opt::Noflip],
    ),
    ("Computing Parameter", &[Opt::Cpus]),
    (
        "Optional Input",
        &[
            Opt::Fam,
            Opt::Bim,
            Opt::Phefile,
            Opt::Covfile,
            Opt::Prunedsnp,
            Opt::Sexchr,
        ],
    ),
    ("Output", &[Opt::Rplot, Opt::Pngplot, Opt::Plink]),
    ("Output Parameter", &[Opt::Prefix, Opt::Rpath]),
];

/// Options whose presence means "an analysis was requested".
///
/// Derived by running `king -b <file> --<opt>` for every option and checking which ones
/// suppress the "Please specify one of the following 24 options" notice. Note that
/// `--rplot` and `--pngplot` count and `--plink` does not.
static ANALYSES: &[Opt] = &[
    Opt::Related,
    Opt::Duplicate,
    Opt::Kinship,
    Opt::Ibdseg,
    Opt::Ibs,
    Opt::MakeGrm,
    Opt::Unrelated,
    Opt::Cluster,
    Opt::Build,
    Opt::Bysample,
    Opt::BySnp,
    Opt::Roh,
    Opt::AutoQc,
    Opt::Pca,
    Opt::Mds,
    Opt::Lmm,
    Opt::Tdt,
    Opt::Gdt,
    Opt::Risk,
    Opt::Rplot,
    Opt::Pngplot,
];

/// Analyses that run as separate passes, with the names the reference uses for them and
/// in the order it lists them.
///
/// This is the order of the reference's own "Please specify one of the following 24
/// options" list, minus the four entries that have no option in the displayed table
/// (`--homog`, `--pc`, `--pcgdt` and two empty slots) and minus `--cluster`, which
/// counts as an analysis but never appears in the "will run separately" line. Note the
/// renames: `--makeGRM` is listed as `--grm`, `--bySNP` as `--bysnp` and `--lmm` as
/// `--mtscore`.
static SEPARATE_ANALYSES: &[(Opt, &str)] = &[
    (Opt::Related, "related"),
    (Opt::Kinship, "kinship"),
    (Opt::AutoQc, "autoQC"),
    (Opt::Lmm, "mtscore"),
    (Opt::Risk, "risk"),
    (Opt::Ibs, "ibs"),
    (Opt::Ibdseg, "ibdseg"),
    (Opt::Mds, "mds"),
    (Opt::Pca, "pca"),
    (Opt::Build, "build"),
    (Opt::Bysample, "bysample"),
    (Opt::BySnp, "bysnp"),
    (Opt::Tdt, "tdt"),
    (Opt::Unrelated, "unrelated"),
    (Opt::Duplicate, "duplicate"),
    (Opt::Roh, "roh"),
    (Opt::MakeGrm, "grm"),
    (Opt::Gdt, "gdt"),
];

/// The name the reference echoes for an analysis under `Options in effect:`.
///
/// Three analyses are echoed under a different spelling from the one that invokes them —
/// `--bySNP` prints as `--bysnp`, `--makeGRM` as `--grm` and `--lmm` as `--mtscore` —
/// exactly as they are listed in the "will run separately" line.
pub fn echo_name(opt: Opt) -> &'static str {
    SEPARATE_ANALYSES
        .iter()
        .find(|(o, _)| *o == opt)
        .map(|(_, name)| *name)
        .unwrap_or_else(|| opt.name())
}

/// Every option, in print order.
pub fn all() -> impl Iterator<Item = Opt> {
    GROUPS.iter().flat_map(|(_, opts)| opts.iter().copied())
}

/// `--noscreen` and `--minConc` overlap in the reference binary's memory.
///
/// Probing shows an `i32` and an `f64` that overlap with a **one byte** offset:
///
/// ```text
///   byte:      0    1    2    3    4    5    6    7    8
///   noscreen [-------i32-------]
///   minConc       [---------------f64------------------]
/// ```
///
/// The buffer starts zeroed and `minConc`'s default `0.8` is written into it, so the
/// `i32` reads back the low three mantissa bytes of `0.8` shifted up by one byte —
/// which is exactly the famous `-1717986816` that the reference prints as `--noscreen`'s
/// default. The model was confirmed against the binary on every combination probed:
///
/// ```text
///   king                                 --noscreen [-1717986816]
///   king --minConc 0.9                   --noscreen [-858993408]
///   king --noscreen 7 --minConc 0.9      --noscreen [-858993401]   (0xCCCCCD07)
///   king --noscreen 256 --minConc 0.9    --noscreen [-858993408]
///   king --minConc 0.9 --noscreen 7      --noscreen [7]
/// ```
#[derive(Clone, Copy, Debug)]
struct Overlap([u8; 9]);

impl Overlap {
    fn new(min_conc: f64) -> Self {
        let mut o = Overlap([0u8; 9]);
        o.set_min_conc(min_conc);
        o
    }
    fn noscreen(&self) -> i32 {
        i32::from_le_bytes([self.0[0], self.0[1], self.0[2], self.0[3]])
    }
    fn set_noscreen(&mut self, v: i32) {
        self.0[0..4].copy_from_slice(&v.to_le_bytes());
    }
    fn min_conc(&self) -> f64 {
        let mut b = [0u8; 8];
        b.copy_from_slice(&self.0[1..9]);
        f64::from_le_bytes(b)
    }
    fn set_min_conc(&mut self, v: f64) {
        self.0[1..9].copy_from_slice(&v.to_le_bytes());
    }
}

/// Parsed command line: every option's value plus the binary fileset name.
#[derive(Clone, Debug)]
pub struct Options {
    /// `-b NAME` / `-bNAME`. Used verbatim as the `.bed` path; the block labels it
    /// "Binary File" and annotates it `(-bname)`.
    pub bed: String,
    flags: [bool; Opt::COUNT],
    ints: [i32; Opt::COUNT],
    doubles: [f64; Opt::COUNT],
    strings: Vec<String>,
    /// Set when a `Kind::Double` option actually consumed a value token. The reference
    /// prints a double whenever it is non-zero **or** was explicitly given, which is why
    /// `--maxP 0` prints `[0.00]` while an untouched `--maxP` prints nothing at all.
    touched: [bool; Opt::COUNT],
    overlap: Overlap,
}

impl Default for Options {
    fn default() -> Self {
        let mut o = Options {
            bed: String::new(),
            flags: [false; Opt::COUNT],
            ints: [0; Opt::COUNT],
            doubles: [0.0; Opt::COUNT],
            strings: vec![String::new(); Opt::COUNT],
            touched: [false; Opt::COUNT],
            overlap: Overlap::new(0.80),
        };
        o.ints[Opt::Sexchr as usize] = 23;
        o.strings[Opt::Prefix as usize] = "king".to_string();
        o
    }
}

impl Options {
    /// KING's defaults.
    pub fn new() -> Self {
        Self::default()
    }

    /// Value of a [`Kind::Flag`] option.
    pub fn flag(&self, o: Opt) -> bool {
        debug_assert_eq!(o.kind(), Kind::Flag);
        self.flags[o as usize]
    }

    /// Value of a [`Kind::Int`] option.
    pub fn int(&self, o: Opt) -> i32 {
        debug_assert_eq!(o.kind(), Kind::Int);
        if o == Opt::Noscreen {
            self.overlap.noscreen()
        } else {
            self.ints[o as usize]
        }
    }

    /// Value of a [`Kind::Double`] option.
    pub fn double(&self, o: Opt) -> f64 {
        debug_assert_eq!(o.kind(), Kind::Double);
        if o == Opt::MinConc {
            self.overlap.min_conc()
        } else {
            self.doubles[o as usize]
        }
    }

    /// Value of a [`Kind::Str`] option.
    pub fn string(&self, o: Opt) -> &str {
        debug_assert_eq!(o.kind(), Kind::Str);
        &self.strings[o as usize]
    }

    /// Whether a double option consumed an explicit value; see [`Options::touched`].
    pub fn was_given(&self, o: Opt) -> bool {
        self.touched[o as usize]
    }

    /// Whether any option that counts as an analysis was requested.
    pub fn any_analysis(&self) -> bool {
        ANALYSES.iter().any(|&o| self.flags[o as usize])
    }

    /// The analyses requested, in table order, spelled with their leading `--`.
    ///
    /// This is what the reference lists under "Options in effect:".
    pub fn analyses_in_effect(&self) -> Vec<String> {
        ANALYSES
            .iter()
            .filter(|&&o| self.flags[o as usize])
            .map(|o| format!("--{}", o.name()))
            .collect()
    }

    /// The analyses that will each run as their own pass, named and ordered the way the
    /// reference's "The following analyses will run separately" line does.
    ///
    /// Empty when `--cluster` is on: the reference never prints that line then, whatever
    /// else was asked for.
    pub fn separate_analyses(&self) -> Vec<String> {
        if self.flags[Opt::Cluster as usize] {
            return Vec::new();
        }
        SEPARATE_ANALYSES
            .iter()
            .filter(|(o, _)| self.flags[*o as usize])
            .map(|(_, name)| format!("--{name}"))
            .collect()
    }

    /// A copy of this command line with `--degree` unset.
    ///
    /// The two "too few samples" downgrades — `--related` below ten samples and
    /// `--ibdseg` below five — hand the run to the `--kinship` pass, and the reference
    /// runs that pass as if `--degree` had never been typed: `--ibdseg --degree 2` on
    /// `singleton` prints the unfiltered `Between-family kinship data saved in file
    /// king.kin0` plus the `Note --kinship --degree <n> …` hint, not the filtered form.
    /// The echoed `Options in effect:` block is a bare `--kinship` for the same reason.
    pub fn without_degree(&self) -> Options {
        let mut o = self.clone();
        o.ints[Opt::Degree as usize] = 0;
        o
    }

    fn toggle_flag(&mut self, o: Opt) {
        self.flags[o as usize] = !self.flags[o as usize];
    }
    fn set_int(&mut self, o: Opt, v: i32) {
        if o == Opt::Noscreen {
            self.overlap.set_noscreen(v);
        } else {
            self.ints[o as usize] = v;
        }
    }
    fn set_double(&mut self, o: Opt, v: f64) {
        if o == Opt::MinConc {
            self.overlap.set_min_conc(v);
        } else {
            self.doubles[o as usize] = v;
        }
        self.touched[o as usize] = true;
    }
    fn set_string(&mut self, o: Opt, v: &str) {
        self.strings[o as usize] = v.to_string();
    }
}

/// Result of parsing a command line.
#[derive(Clone, Debug)]
pub struct Parsed {
    pub options: Options,
    /// Messages for the WARNING block, in the order they were produced. Parsing always
    /// continues after a warning.
    pub warnings: Vec<String>,
}

/// Outcome of matching a `--name` against the option table.
enum Match {
    Found(Opt),
    Ambiguous,
    Undefined,
}

/// Case-insensitive prefix match, exactly like the reference:
/// `--rel` sets `--related`, `--RELATED` works, `--r` is ambiguous, `--duplicatex` is
/// undefined.
fn lookup(name: &str) -> Match {
    let needle = name.to_ascii_lowercase();
    let hits: BTreeSet<Opt> = all()
        .filter(|o| o.name().to_ascii_lowercase().starts_with(&needle))
        .collect();
    match hits.len() {
        0 => Match::Undefined,
        1 => Match::Found(*hits.iter().next().expect("one hit")),
        _ => Match::Ambiguous,
    }
}

/// Does this token look like an integer to the reference binary?
///
/// Probed: `-3`, `+3`, `007`, `4294967296` and even a bare `-` are consumed;
/// `3x`, `0x10`, `.5`, `3.0`, `3-4`, `" 5"` and `""` are not (they end up in the
/// WARNING block as "ignored").
fn looks_like_int(s: &str) -> bool {
    let mut it = s.chars();
    match it.next() {
        Some(c) if c == '+' || c == '-' => {}
        Some(c) if c.is_ascii_digit() => {}
        _ => return false,
    }
    it.all(|c| c.is_ascii_digit())
}

/// Does this token look like a double to the reference binary?
///
/// Two rules, both derived by probing and both needed:
///
/// 1. the first character must be a digit, `.`, `+` or `-` — so `inf`, `nan`, `e5`,
///    `" 5"` and `""` are rejected while `1x`, `-x`, `.x` and `1.2.3` are accepted;
/// 2. if the token contains an `e`/`E`, what follows the first one must be non-empty and
///    consist of an optional sign and digits only.
///
/// Rule 2 explains a whole family of results at once: `1e` and `12e` are rejected but
/// `1e+` is accepted, `1e3` is accepted but `1e5x` and `1e2e` are not — and, crucially,
/// `--seglength --sexchr` leaves `--sexchr` alone (the `e` is followed by `xchr`) while
/// `--seglength --minConc` swallows `--minConc` as a value of 0.0, because `minConc` has
/// no `e` in it at all. Both were confirmed against the binary.
fn looks_like_double(s: &str) -> bool {
    match s.chars().next() {
        Some(c) if c.is_ascii_digit() || c == '.' || c == '+' || c == '-' => {}
        _ => return false,
    }
    let Some(exp) = s.find(['e', 'E']) else {
        return true;
    };
    let after = &s[exp + 1..];
    if after.is_empty() {
        return false;
    }
    let digits = after.strip_prefix(['+', '-']).unwrap_or(after);
    digits.bytes().all(|b| b.is_ascii_digit())
}

/// C `atoi` semantics: leading sign, decimal digits, wrap into `int`.
///
/// `strtol` saturates at 64 bits and the cast to `int` truncates, which is why
/// `--degree 4294967296` reads back as 0.
fn c_atoi(s: &str) -> i32 {
    let (neg, digits) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };
    let mut acc: i64 = 0;
    for c in digits.chars() {
        let Some(d) = c.to_digit(10) else { break };
        acc = match acc.checked_mul(10).and_then(|a| a.checked_add(d as i64)) {
            Some(v) => v,
            None => {
                // strtol clamps, then the narrowing cast truncates.
                return if neg {
                    i64::MIN as i32
                } else {
                    i64::MAX as i32
                };
            }
        };
    }
    (if neg { -acc } else { acc }) as i32
}

/// C `atof` semantics: parse the longest numeric prefix, 0.0 when there is none.
///
/// `1.2.3` reads as 1.2 and `-x` reads as 0.0. C99 hex literals count too, which is how
/// the reference turns `--callrateN 0x10` into `[16.00]`, `0X1A` into `[26.00]`,
/// `0x1p4` into `[16.00]` and `0x.8p1` into `[1.00]`.
fn c_atof(s: &str) -> f64 {
    let b = s.as_bytes();
    let mut i = 0;
    if i < b.len() && (b[i] == b'+' || b[i] == b'-') {
        i += 1;
    }
    if let Some(v) = c_atof_hex(&s[i..]) {
        return if b.first() == Some(&b'-') { -v } else { v };
    }
    let mut digits = 0;
    while i < b.len() && b[i].is_ascii_digit() {
        i += 1;
        digits += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        while i < b.len() && b[i].is_ascii_digit() {
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return 0.0;
    }
    // Optional exponent, only when it is complete.
    if i < b.len() && (b[i] == b'e' || b[i] == b'E') {
        let mut j = i + 1;
        if j < b.len() && (b[j] == b'+' || b[j] == b'-') {
            j += 1;
        }
        let digits_start = j;
        while j < b.len() && b[j].is_ascii_digit() {
            j += 1;
        }
        if j > digits_start {
            i = j;
        }
    }
    s[..i].parse::<f64>().unwrap_or(0.0)
}

/// The unsigned `0x…` branch of [`c_atof`]. `None` when `s` is not a hex literal, in
/// which case the decimal scanner takes over (`0x` alone reads as the plain `0`).
fn c_atof_hex(s: &str) -> Option<f64> {
    let rest = s.strip_prefix("0x").or_else(|| s.strip_prefix("0X"))?;
    let b = rest.as_bytes();
    let mut i = 0;
    let mut value = 0.0f64;
    let mut digits = 0;
    while i < b.len() && b[i].is_ascii_hexdigit() {
        value = value * 16.0 + f64::from((b[i] as char).to_digit(16).expect("hex digit"));
        i += 1;
        digits += 1;
    }
    if i < b.len() && b[i] == b'.' {
        i += 1;
        let mut scale = 1.0 / 16.0;
        while i < b.len() && b[i].is_ascii_hexdigit() {
            value += f64::from((b[i] as char).to_digit(16).expect("hex digit")) * scale;
            scale /= 16.0;
            i += 1;
            digits += 1;
        }
    }
    if digits == 0 {
        return None;
    }
    // Binary exponent, only when it actually has digits behind it.
    if i < b.len() && (b[i] == b'p' || b[i] == b'P') {
        let mut j = i + 1;
        let neg = b.get(j) == Some(&b'-');
        if neg || b.get(j) == Some(&b'+') {
            j += 1;
        }
        let start = j;
        let mut exp: i32 = 0;
        while j < b.len() && b[j].is_ascii_digit() {
            exp = exp
                .saturating_mul(10)
                .saturating_add(i32::from(b[j] - b'0'));
            j += 1;
        }
        if j > start {
            value *= 2f64.powi(if neg { -exp } else { exp });
        }
    }
    Some(value)
}

/// Resolve `name` and apply it. `spelling` is the text to quote in a warning, `next` the
/// following token when the option is allowed to take a value.
///
/// Returns whether `next` was consumed.
fn apply_long(
    options: &mut Options,
    warnings: &mut Vec<String>,
    spelling: &str,
    name: &str,
    next: Option<&str>,
) -> bool {
    let opt = match lookup(name) {
        Match::Undefined => {
            warnings.push(format!("Command line parameter {spelling} is undefined"));
            return false;
        }
        Match::Ambiguous => {
            warnings.push(format!("Command line parameter {spelling} is ambiguous"));
            return false;
        }
        Match::Found(opt) => opt,
    };

    match opt.kind() {
        // Flags never consume a value; a value left behind is warned about on its own
        // turn round the loop. Repeating a flag toggles it: `--related --related`
        // leaves it off.
        Kind::Flag => {
            options.toggle_flag(opt);
            false
        }
        Kind::Int => match next {
            Some(v) if looks_like_int(v) => {
                options.set_int(opt, c_atoi(v));
                true
            }
            // No value: the reference toggles the integer, which is how `--degree`
            // alone becomes 1 and `--sexchr` alone becomes 0.
            _ => {
                let cur = options.int(opt);
                options.set_int(opt, i32::from(cur == 0));
                false
            }
        },
        Kind::Double => match next {
            Some(v) if looks_like_double(v) => {
                options.set_double(opt, c_atof(v));
                true
            }
            // No value: left untouched, unlike the integer case.
            _ => false,
        },
        // Strings swallow whatever comes next, even another option:
        // `--trait --related` sets trait to the literal "--related".
        Kind::Str => match next {
            Some(v) => {
                options.set_string(opt, v);
                true
            }
            None => false,
        },
    }
}

/// Parse a command line. `args` must exclude the program name; warning messages number
/// arguments from 1, matching the reference's `(#N)`.
pub fn parse(args: &[String]) -> Parsed {
    let mut options = Options::new();
    let mut warnings = Vec::new();
    let mut i = 0usize;

    while i < args.len() {
        let tok = args[i].as_str();
        let argno = i + 1;

        // A lone `--` borrows the next token as its name — `king -- related` turns
        // --related on — but only when that token is not itself dash-led, and the
        // resulting option never takes a value: `king -- degree 3` leaves the 3 behind.
        let merged = tok == "--" && args.get(i + 1).is_some_and(|next| !next.starts_with('-'));

        if merged {
            let name = args[i + 1].as_str();
            apply_long(
                &mut options,
                &mut warnings,
                &format!("--{name}"),
                name,
                None,
            );
            i += 2;
            continue;
        }

        if let Some(name) = tok.strip_prefix("--") {
            let consumed = apply_long(
                &mut options,
                &mut warnings,
                tok,
                name,
                args.get(i + 1).map(String::as_str),
            );
            if consumed {
                i += 1;
            }
        } else if tok.len() >= 2 && tok.starts_with('-') && matches!(&tok[1..2], "b" | "B") {
            if tok.len() > 2 {
                options.bed = tok[2..].to_string();
            } else {
                // `-b` takes the next token unless it looks like another option, and
                // clears the fileset when it cannot: `-bx.bed -b` ends up with no
                // fileset at all.
                match args.get(i + 1) {
                    Some(v) if !v.starts_with('-') => {
                        options.bed = v.clone();
                        i += 1;
                    }
                    _ => options.bed.clear(),
                }
            }
        } else {
            warnings.push(format!("Command line parameter {tok} (#{argno}) ignored"));
        }

        i += 1;
    }

    Parsed { options, warnings }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse_str(args: &[&str]) -> Parsed {
        let owned: Vec<String> = args.iter().map(|s| s.to_string()).collect();
        parse(&owned)
    }

    #[test]
    fn table_is_complete_and_ordered() {
        let listed: Vec<Opt> = all().collect();
        assert_eq!(listed.len(), Opt::COUNT);
        for (i, o) in listed.iter().enumerate() {
            assert_eq!(*o as usize, i, "GROUPS order must match Opt discriminants");
        }
    }

    #[test]
    fn no_name_is_a_prefix_of_another() {
        // Needed for prefix matching to be unambiguous on full spellings.
        for a in all() {
            for b in all() {
                if a != b {
                    assert!(
                        !b.name()
                            .to_ascii_lowercase()
                            .starts_with(&a.name().to_ascii_lowercase()),
                        "{} is a prefix of {}",
                        a.name(),
                        b.name()
                    );
                }
            }
        }
    }

    #[test]
    fn defaults_match_reference() {
        let o = Options::new();
        assert_eq!(o.int(Opt::Sexchr), 23);
        assert_eq!(o.string(Opt::Prefix), "king");
        assert_eq!(o.double(Opt::MinConc), 0.80);
        assert_eq!(o.int(Opt::Noscreen), -1_717_986_816);
        assert_eq!(o.int(Opt::Degree), 0);
        assert!(!o.flag(Opt::Related));
        assert!(o.bed.is_empty());
    }

    #[test]
    fn noscreen_overlaps_minconc() {
        // Every line here was verified against the reference binary.
        assert_eq!(
            parse_str(&["--minConc", "0.9"]).options.int(Opt::Noscreen),
            -858_993_408
        );
        assert_eq!(
            parse_str(&["--noscreen", "7", "--minConc", "0.9"])
                .options
                .int(Opt::Noscreen),
            -858_993_401
        );
        assert_eq!(
            parse_str(&["--noscreen", "256", "--minConc", "0.9"])
                .options
                .int(Opt::Noscreen),
            -858_993_408
        );
        assert_eq!(
            parse_str(&["--minConc", "0.9", "--noscreen", "7"])
                .options
                .int(Opt::Noscreen),
            7
        );
        // Writing noscreen corrupts minConc's low mantissa bytes but not its value at
        // two decimals, which is why the reference still prints [0.80].
        let p = parse_str(&["--noscreen", "7"]);
        assert!((p.options.double(Opt::MinConc) - 0.8).abs() < 1e-6);
    }

    #[test]
    fn prefix_matching() {
        assert!(parse_str(&["--rel"]).options.flag(Opt::Related));
        assert!(parse_str(&["--RELATED"]).options.flag(Opt::Related));
        assert!(parse_str(&["--duplicat"]).options.flag(Opt::Duplicate));
        assert_eq!(
            parse_str(&["--r"]).warnings,
            vec!["Command line parameter --r is ambiguous"]
        );
        assert_eq!(
            parse_str(&["--duplicatex"]).warnings,
            vec!["Command line parameter --duplicatex is undefined"]
        );
        assert_eq!(
            parse_str(&["--"]).warnings,
            vec!["Command line parameter -- is ambiguous"]
        );
    }

    #[test]
    fn int_value_predicate() {
        assert_eq!(parse_str(&["--degree", "3"]).options.int(Opt::Degree), 3);
        assert_eq!(parse_str(&["--degree", "-3"]).options.int(Opt::Degree), -3);
        assert_eq!(parse_str(&["--degree", "+3"]).options.int(Opt::Degree), 3);
        assert_eq!(parse_str(&["--degree", "007"]).options.int(Opt::Degree), 7);
        assert_eq!(parse_str(&["--degree", "-"]).options.int(Opt::Degree), 0);
        assert_eq!(
            parse_str(&["--degree", "4294967296"])
                .options
                .int(Opt::Degree),
            0
        );
        // Rejected value tokens toggle the option and are reported as ignored.
        for bad in ["3x", "0x10", ".5", "3.0", "3-4", " 5", ""] {
            let p = parse_str(&["--degree", bad]);
            assert_eq!(p.options.int(Opt::Degree), 1, "{bad}");
            assert_eq!(
                p.warnings,
                vec![format!("Command line parameter {bad} (#2) ignored")],
                "{bad}"
            );
        }
        // Bare integer options toggle: 0 becomes 1, non-zero becomes 0.
        assert_eq!(parse_str(&["--degree"]).options.int(Opt::Degree), 1);
        assert_eq!(parse_str(&["--sexchr"]).options.int(Opt::Sexchr), 0);
        assert_eq!(parse_str(&["--noscreen"]).options.int(Opt::Noscreen), 0);
    }

    #[test]
    fn double_value_predicate() {
        assert_eq!(
            parse_str(&["--seglength", "1.2.3"])
                .options
                .double(Opt::Seglength),
            1.2
        );
        assert_eq!(
            parse_str(&["--seglength", "1x"])
                .options
                .double(Opt::Seglength),
            1.0
        );
        assert_eq!(
            parse_str(&["--seglength", "-x"])
                .options
                .double(Opt::Seglength),
            0.0
        );
        assert_eq!(
            parse_str(&["--seglength", "1e3"])
                .options
                .double(Opt::Seglength),
            1000.0
        );
        assert_eq!(
            parse_str(&["--seglength", ".5"])
                .options
                .double(Opt::Seglength),
            0.5
        );
        assert_eq!(
            parse_str(&["--seglength", "5."])
                .options
                .double(Opt::Seglength),
            5.0
        );
        for bad in ["inf", "nan", "e5", " 5", ""] {
            let p = parse_str(&["--seglength", bad]);
            assert!(!p.options.was_given(Opt::Seglength), "{bad}");
            assert_eq!(
                p.warnings,
                vec![format!("Command line parameter {bad} (#2) ignored")],
                "{bad}"
            );
        }
        // A bare double is left alone, and stays untouched.
        let p = parse_str(&["--seglength"]);
        assert!(!p.options.was_given(Opt::Seglength));
    }

    #[test]
    fn strings_swallow_anything() {
        let p = parse_str(&["--trait", "--related"]);
        assert_eq!(p.options.string(Opt::Trait), "--related");
        assert!(!p.options.flag(Opt::Related));
        assert!(p.warnings.is_empty());
        // Nothing left to take: unchanged, no warning.
        assert_eq!(parse_str(&["--trait"]).options.string(Opt::Trait), "");
    }

    #[test]
    fn binary_fileset() {
        assert_eq!(parse_str(&["-b", "x.bed"]).options.bed, "x.bed");
        assert_eq!(parse_str(&["-bx.bed"]).options.bed, "x.bed");
        assert_eq!(parse_str(&["-B", "x.bed"]).options.bed, "x.bed");
        // A dash-leading token is not a file name.
        let p = parse_str(&["-b", "--related"]);
        assert_eq!(p.options.bed, "");
        assert!(p.options.flag(Opt::Related));
        let p = parse_str(&["-b", "-3"]);
        assert_eq!(p.options.bed, "");
        assert_eq!(p.warnings, vec!["Command line parameter -3 (#2) ignored"]);
    }

    #[test]
    fn unknown_tokens() {
        assert_eq!(
            parse_str(&["--unknownopt", "value"]).warnings,
            vec![
                "Command line parameter --unknownopt is undefined",
                "Command line parameter value (#2) ignored",
            ]
        );
        assert_eq!(
            parse_str(&["-x"]).warnings,
            vec!["Command line parameter -x (#1) ignored"]
        );
        assert_eq!(
            parse_str(&["foo"]).warnings,
            vec!["Command line parameter foo (#1) ignored"]
        );
        assert_eq!(
            parse_str(&["--related=1"]).warnings,
            vec!["Command line parameter --related=1 is undefined"]
        );
    }

    #[test]
    fn flags_toggle() {
        // Found by differential fuzzing: a repeated flag is not idempotent.
        assert!(parse_str(&["--related"]).options.flag(Opt::Related));
        assert!(!parse_str(&["--related", "--related"])
            .options
            .flag(Opt::Related));
        assert!(parse_str(&["--rel", "--RELATED", "--related"])
            .options
            .flag(Opt::Related));
    }

    #[test]
    fn lone_dash_dash_borrows_the_next_token_as_a_name() {
        assert!(parse_str(&["--", "related"]).options.flag(Opt::Related));
        assert!(parse_str(&["--", "rel"]).options.flag(Opt::Related));
        // The borrowed option never takes a value of its own.
        let p = parse_str(&["--", "degree", "3"]);
        assert_eq!(p.options.int(Opt::Degree), 1);
        assert_eq!(p.warnings, vec!["Command line parameter 3 (#3) ignored"]);
        // A dash-led follower is not borrowed.
        let p = parse_str(&["--", "--", "related"]);
        assert!(p.options.flag(Opt::Related));
        assert_eq!(p.warnings, vec!["Command line parameter -- is ambiguous"]);
    }

    #[test]
    fn value_predicates_ignore_option_lookalikes_by_spelling_only() {
        // `--sexchr` has an `e` followed by junk so it is not a number; `--minConc`
        // has no `e` at all, so the reference eats it as a value of 0.0.
        let p = parse_str(&["--seglength", "--sexchr"]);
        assert!(!p.options.was_given(Opt::Seglength));
        assert_eq!(p.options.int(Opt::Sexchr), 0, "--sexchr still got its turn");
        let p = parse_str(&["--seglength", "--minConc"]);
        assert!(p.options.was_given(Opt::Seglength));
        assert_eq!(p.options.double(Opt::Seglength), 0.0);
    }

    #[test]
    fn hex_values() {
        assert_eq!(
            parse_str(&["--callrateN", "0x10"])
                .options
                .double(Opt::CallrateN),
            16.0
        );
        assert_eq!(
            parse_str(&["--callrateN", "0X1A"])
                .options
                .double(Opt::CallrateN),
            26.0
        );
        assert_eq!(
            parse_str(&["--callrateN", "0x1p4"])
                .options
                .double(Opt::CallrateN),
            16.0
        );
        assert_eq!(
            parse_str(&["--callrateN", "0x.8p1"])
                .options
                .double(Opt::CallrateN),
            1.0
        );
        // `0x` alone falls back to the plain decimal 0.
        assert_eq!(
            parse_str(&["--callrateN", "0x"])
                .options
                .double(Opt::CallrateN),
            0.0
        );
        // Integers do not take hex at all.
        let p = parse_str(&["--degree", "0x10"]);
        assert_eq!(p.options.int(Opt::Degree), 1);
    }

    #[test]
    fn separate_analyses_use_the_reference_names_and_order() {
        assert_eq!(
            parse_str(&["--duplicate", "--related"])
                .options
                .separate_analyses(),
            vec!["--related", "--duplicate"]
        );
        assert_eq!(
            parse_str(&["--tdt", "--build", "--autoQC"])
                .options
                .separate_analyses(),
            vec!["--autoQC", "--build", "--tdt"]
        );
        assert_eq!(
            parse_str(&["--makeGRM", "--pca", "--lmm", "--bySNP"])
                .options
                .separate_analyses(),
            vec!["--mtscore", "--pca", "--bysnp", "--grm"]
        );
        // --cluster silences the line entirely.
        assert!(parse_str(&["--unrelated", "--cluster", "--build"])
            .options
            .separate_analyses()
            .is_empty());
        // ...and these three never appear in it, though they do count as analyses.
        assert!(parse_str(&["--rplot", "--pngplot"])
            .options
            .separate_analyses()
            .is_empty());
    }

    #[test]
    fn analyses_in_effect() {
        assert!(!Options::new().any_analysis());
        assert!(parse_str(&["--kinship"]).options.any_analysis());
        assert!(parse_str(&["--rplot"]).options.any_analysis());
        // Parameters are not analyses.
        assert!(!parse_str(&["--plink"]).options.any_analysis());
        assert!(!parse_str(&["--pcs", "2"]).options.any_analysis());
        assert_eq!(
            parse_str(&["--ibs", "--related"])
                .options
                .analyses_in_effect(),
            vec!["--related", "--ibs"]
        );
    }
}
