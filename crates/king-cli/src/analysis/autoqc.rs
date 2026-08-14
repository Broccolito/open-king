//! `--autoQC` — the automated QC pipeline, and the four report files it writes.
//!
//! Eight numbered steps, three of them conditional on the map carrying sex chromosomes:
//!
//! ```text
//! 1  SNP call rate (autosomes + X), then monomorphic markers (autosomes + X + Y)
//! 2  sample call rate, over the surviving autosomes
//! 3  SNP call rate again (autosomes + X), over the surviving samples
//! 4  Y call rate in males / X heterozygosity in males / Y call rate in females  [X and Y]
//! 5  gender check against the Y-marker count, and inference for unreported sexes [X and Y]
//! 6  the same three filters at tighter thresholds                               [X and Y]
//! 7  final counts
//! 8  the summary table, which is also written to `<prefix>_autoQC_Summary.txt`
//! ```
//!
//! Everything below was established by running the reference binary; each rule names the
//! experiment that fixes it. The pass is **class-local**: autosomes (`XY` included, as
//! everywhere else), X and Y each carry their own marker array, and MT is loaded, counted
//! and then never looked at again — `sexchr`'s 11 monomorphic MT markers survive QC
//! untouched, and MT appears in no console line at all.
//!
//! # The dropped word
//!
//! `--autoQC` packs genotypes **16 markers to a word** — it says so itself, reporting
//! `ceil(m / 16)` words where every other analysis reports `ceil(m / 64)` — and its
//! scanning loops walk `16 * (words - 1)` markers plus `m % 16`. When `m` is a multiple of
//! 16 that last term is zero and **the final 16 markers of the class are never examined**:
//! not filtered, not listed in the removal file, not counted in the final tally. This is
//! the reference's own arithmetic bug, and reproducing it is load-bearing for parity — it
//! is why `bigish` reports 50 000 raw markers, removes none, and finishes with 49 984.
//!
//! Swept directly: a 1..40 sweep of a clean autosomal fileset finishes with `m` markers at
//! every `m` except 16, 32 and 48, where it finishes with `m - 16`. The same holds
//! independently for X (32 X markers with the last 8 monomorphic → `0 SNPs are
//! monomorphic` and 16 X markers left) and for Y. MT is exempt because nothing scans it.
//! The truncation is decided **once, from the class's marker count at load**, and is not
//! recomputed as markers are removed: 100 markers of which 20 are monomorphic finish at
//! 80, not at 64.
//!
//! The console headers and the gender cutoff use the *untruncated* count — `step 2 …
//! (with 5232 SNPs)` on `missing` includes the 16 markers the pass cannot see, and 32 Y
//! markers give a cutoff of `16.0` although no sample's Y count can exceed 16.
//!
//! # Thresholds
//!
//! `--callrateM` (default 0.95) is the SNP threshold and `--callrateN` (default 0.95) the
//! sample threshold. Step 1 uses a looser threshold derived from `--callrateM`; swept:
//!
//! ```text
//! 0.5 → 40%   0.7 → 60%   0.8 → 70%   0.85 → 70%   0.9 → 80%   0.95 → 80%   1.2 → 80%
//! ```
//!
//! which is [`step1_threshold`]. Every call-rate test is written against the **missing**
//! rate rather than the call rate: `missing / n > 1 - threshold`. The difference shows at
//! exact boundaries, because `1 - 0.9` is `0.09999999999999998` in binary — a marker
//! missing in exactly 2 of 20 samples *is* removed under `--callrateM 0.9`, while a marker
//! at exactly 95% survives the default (`0.05 < 0.05000000000000004`). Writing the test as
//! `call_rate < threshold` keeps the second case and loses the first.
//!
//! # Removal reasons
//!
//! Reason strings carry the threshold as a rounded integer percentage, and both Y
//! call-rate filters (steps 4a and 6a) write the **`--callrateM`** percentage even though
//! step 4a applies the looser step-1 threshold: under `--callrateM 0.9`, a marker removed
//! by 4a's 80% test is still labelled `CallRateLessThan90`.

use std::fmt::Write as _;
use std::io::Write;

use king_io::{Sample, Variant};

use crate::analysis::out_path;
use crate::analysis::qc::Calls;
use crate::cli::{Opt, Options};
use crate::console;
use crate::load::{Class, Loaded};

/// Markers per word in `--autoQC`'s packing, and so the size of the tail it drops.
const WORD: usize = 16;

/// Ceiling of the step-1 call-rate threshold, whatever `--callrateM` says.
const STEP1_CEILING: f64 = 0.8;

/// X heterozygosity above this, in males, removes the marker in step 4b.
const X_HET_LOOSE: f64 = 0.05;
/// The same filter in step 6b.
const X_HET_TIGHT: f64 = 0.01;
/// Y call rate above this, in females, removes the marker in step 4c.
const Y_IN_FEMALES_LOOSE: f64 = 0.10;
/// The same filter in step 6c.
const Y_IN_FEMALES_TIGHT: f64 = 0.02;

/// Default X-heterozygosity cut for the gender check, and the ceiling of the fitted one.
const X_HET_CUT: f64 = 0.10;
/// A female's heterozygosity has to clear this before it can lower the cut.
const X_HET_FLOOR: f64 = 0.05;

/// `--callrateN` and `--callrateM` when the command line did not give them.
///
/// The default lives here rather than in [`Options`] because the banner distinguishes the
/// two: `--callrateN` prints bare by default and `[0.90]` when given, so the parse state
/// has to keep holding zero. An explicit `--callrateM 0` is therefore taken literally —
/// the reference duly applies a 0% filter and a −10% one in step 1.
const CALLRATE_DEFAULT: f64 = 0.95;

/// A male is unsuspicious only while his share of the Y panel is **above two thirds**, and
/// a female only while hers is **below one third**.
///
/// Both boundaries were bisected on panels of 2, 20, 30, 32, 48, 64 and 100 markers, and
/// only the exact thirds fit all of them: the male flip lands on 14 of 20, 21 of 30, 22 of
/// 32, 33 of 48 and 43 of 64, and the female flip on a 100-marker panel is at 34 — not at
/// the 31 a 0.3 threshold would give. Equality is suspicious on both sides: 32 of 48 and
/// 16 of 48 are both flagged.
///
/// The comparison is on the **share**, not on a scaled count, which is what an emptied Y
/// panel shows: `0 / 0` is NaN, every comparison against it is false, and the reference
/// duly flags nobody at all when step 4 has removed every Y marker.
const Y_MALE_SHARE: f64 = 2.0 / 3.0;
const Y_FEMALE_SHARE: f64 = 1.0 / 3.0;

/// `.fam` sex codes.
const MALE: u8 = 1;
const FEMALE: u8 = 2;
const UNREPORTED: u8 = 0;

// ---------------------------------------------------------------------------
// Marker classes
// ---------------------------------------------------------------------------

/// The three classes `--autoQC` filters. MT is not one of them.
#[derive(Clone, Copy, PartialEq, Eq)]
enum MarkerClass {
    Auto,
    X,
    Y,
}

/// One class's marker array: the rows the pass can see, and the count it reports.
struct Markers {
    /// `.bim` rows of the class in file order, **truncated to the scanned prefix**.
    rows: Vec<usize>,
    /// Whether each entry of `rows` is still in the study.
    alive: Vec<bool>,
    /// The class's marker count as loaded, dropped tail included.
    raw: usize,
    /// `raw` minus everything removed so far — the figure the console headers print, and
    /// therefore *not* the number of markers the pass can actually see.
    logical: usize,
    /// Whether the class lost a whole word to the scan loops.
    truncated: bool,
}

impl Markers {
    fn new(rows: Vec<usize>) -> Markers {
        let raw = rows.len();
        // The dropped word: a whole-word class loses its final word to the scan loops.
        let scanned = if raw >= WORD && raw % WORD == 0 {
            raw - WORD
        } else {
            raw
        };
        Markers {
            rows: rows[..scanned].to_vec(),
            alive: vec![true; scanned],
            raw,
            logical: raw,
            truncated: scanned != raw,
        }
    }

    /// Markers the step-5 per-sample loops count but never read.
    ///
    /// The dropped word is invisible to every *marker* loop and to the step-2 sample
    /// call rate, but the two per-sample statistics of step 5 walk it anyway and take all
    /// 16 of its slots as **called and homozygous**. Two probes pin it: with a 32-marker
    /// Y panel, a female with no Y calls at all is treated as having 16 — one real call
    /// puts her over the 16.0 cutoff and she is reported mislabeled — and with a
    /// 32-marker X panel, females heterozygous at 3 of the 16 visible markers print
    /// `X-chr heterozygosity is set as 0.09`, which is 3/32 and not 3/16.
    fn phantom(&self) -> usize {
        if self.truncated {
            WORD
        } else {
            0
        }
    }

    /// `(slot, .bim row)` for every marker still in the study, in file order.
    fn live(&self) -> Vec<(usize, usize)> {
        self.rows
            .iter()
            .enumerate()
            .filter(|(i, _)| self.alive[*i])
            .map(|(i, &row)| (i, row))
            .collect()
    }

    /// The `.bim` rows still in the study.
    fn live_rows(&self) -> Vec<usize> {
        self.live().into_iter().map(|(_, row)| row).collect()
    }

    /// Markers the pass can still see — what step 7 counts.
    fn visible(&self) -> usize {
        self.alive.iter().filter(|a| **a).count()
    }
}

// ---------------------------------------------------------------------------
// The pass
// ---------------------------------------------------------------------------

/// Everything one `--autoQC` run reads and mutates.
struct AutoQc<'a> {
    calls: &'a Calls,
    samples: &'a [Sample],
    variants: &'a [Variant],
    /// Whether each sample is still in the study.
    alive: Vec<bool>,
    /// `.fam` sex, verbatim — step 5 infers the unreported ones into the report without
    /// writing them back here.
    sex: Vec<u8>,
    auto: Markers,
    x: Markers,
    y: Markers,
    /// MT markers are counted and nothing else, so only the count is kept.
    mt: usize,
    /// `<prefix>_autoQC_snptoberemoved.txt`, in the order the steps append to it.
    snp_file: String,
    /// `<prefix>_autoQC_sampletoberemoved.txt`, likewise.
    sample_file: String,
    /// `<prefix>_autoQC_updatesex.txt`; written only when a sex was inferred.
    updatesex: String,
    /// The summary table's counters, in its own numbering.
    counts: Counts,
}

/// The summary table's counters, named as it numbers them.
#[derive(Default)]
struct Counts {
    low_call_rate: usize,
    monomorphic: usize,
    sample_call_rate: usize,
    snp_call_rate: usize,
    y_in_men_loose: usize,
    x_het_loose: usize,
    y_in_women_loose: usize,
    mislabeled_male: usize,
    mislabeled_female: usize,
    suspicious: usize,
    y_in_men_tight: usize,
    x_het_tight: usize,
    y_in_women_tight: usize,
}

/// Run the `--autoQC` pass: console body plus up to four report files.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let Some(calls) = Calls::read(opts, loaded) else {
        return;
    };
    let mut qc = AutoQc::new(&calls, &loaded.fileset.samples, &loaded.fileset.variants);
    let log = qc.pipeline(opts);
    let _ = out.write_all(log.as_bytes());
}

/// `Autosome genotypes stored in <w> words for each of <n> individuals.`, with
/// `--autoQC`'s own denser packing: `ceil(m / 16)`, not the loader's `ceil(m / 64)`.
pub fn preamble(loaded: &Loaded) -> String {
    format!(
        "{}\n",
        console::autosome_words(
            loaded.counts.autosome.div_ceil(WORD),
            loaded.fileset.samples.len()
        )
    )
}

impl<'a> AutoQc<'a> {
    fn new(calls: &'a Calls, samples: &'a [Sample], variants: &'a [Variant]) -> AutoQc<'a> {
        AutoQc {
            calls,
            samples,
            variants,
            alive: vec![true; samples.len()],
            sex: samples.iter().map(|s| s.sex).collect(),
            auto: Markers::new(calls.rows(Class::is_autosomal)),
            x: Markers::new(calls.rows(|c| c == Class::X)),
            y: Markers::new(calls.rows(|c| c == Class::Y)),
            mt: calls.rows(|c| c == Class::Mt).len(),
            snp_file: "SNP\tREASON\n".to_string(),
            sample_file: "FID\tIID\tREASON\n".to_string(),
            updatesex: String::new(),
            counts: Counts::default(),
        }
    }

    // -- primitives ------------------------------------------------------

    /// Samples still in the study, in `.fam` order.
    fn subjects(&self) -> Vec<usize> {
        (0..self.samples.len()).filter(|&s| self.alive[s]).collect()
    }

    /// Live samples of one **reported** sex — the denominators of every step 4/6 filter.
    ///
    /// A sample whose `.fam` sex is 0 is in neither set, in step 6 as much as in step 4:
    /// that is what keeps `sexchr`'s two sex-0 samples from dragging 40% X heterozygosity
    /// into the male filter and wiping out the X panel, and it is why four samples inferred
    /// in step 5 leave 6b and 6c reporting zero.
    fn of_sex(&self, sex: u8) -> Vec<usize> {
        (0..self.samples.len())
            .filter(|&s| self.alive[s] && self.sex[s] == sex)
            .collect()
    }

    /// Fraction of `subjects` with no call at `row`.
    ///
    /// Every call-rate test compares this against `1 - threshold` rather than comparing a
    /// call rate against the threshold; see the module note on exact boundaries.
    fn missing_rate(&self, row: usize, subjects: &[usize]) -> f64 {
        let missing = subjects
            .iter()
            .filter(|&&s| self.calls.get(s, row).is_none())
            .count();
        missing as f64 / subjects.len() as f64
    }

    /// Heterozygosity at `row` over the `subjects` that have a call there.
    ///
    /// The denominator is the **called** subjects: a marker heterozygous in 1 of 18 called
    /// males (2 of the 20 uncalled) is 5.6% and step 4b removes it, where 1/20 would be 5%
    /// and would survive.
    fn het_rate(&self, row: usize, subjects: &[usize]) -> f64 {
        let mut called = 0usize;
        let mut het = 0usize;
        for &s in subjects {
            if let Some(dosage) = self.calls.get(s, row) {
                called += 1;
                het += usize::from(dosage == 1);
            }
        }
        if called == 0 {
            0.0
        } else {
            het as f64 / called as f64
        }
    }

    /// One sample's heterozygosity over the live markers of a class.
    ///
    /// Called markers only — a female missing 50 of 100 X markers and heterozygous at 5
    /// of the 50 she has is at 0.10, not 0.05 — plus the class's phantom word, which
    /// enlarges the denominator without contributing a heterozygote.
    ///
    /// With nothing called at all the quotient is `0 / 0`, and C's NaN makes every
    /// comparison against it false: a fileset whose whole X panel went monomorphic in
    /// step 1 flags nobody, where a 0.0 here would flag every female.
    fn sample_het(&self, sample: usize, rows: &[usize], phantom: usize) -> f64 {
        let mut called = phantom;
        let mut het = 0usize;
        for &row in rows {
            if let Some(dosage) = self.calls.get(sample, row) {
                called += 1;
                het += usize::from(dosage == 1);
            }
        }
        het as f64 / called as f64
    }

    /// Whether every call at `row` is the same homozygote — the reference's "monomorphic".
    ///
    /// A marker whose calls are *all heterozygous* is polymorphic and the reference agrees:
    /// an all-het column survives while an all-hom one is removed. A marker with no calls
    /// at all counts as monomorphic.
    fn monomorphic(&self, row: usize, subjects: &[usize]) -> bool {
        let mut seen_a1 = false;
        let mut seen_a2 = false;
        for &s in subjects {
            match self.calls.get(s, row) {
                Some(2) => seen_a1 = true,
                Some(0) => seen_a2 = true,
                Some(_) => return false,
                None => {}
            }
            if seen_a1 && seen_a2 {
                return false;
            }
        }
        true
    }

    /// How many of `rows` a sample has a call at — the gender check's Y count.
    fn call_count(&self, sample: usize, rows: &[usize]) -> usize {
        rows.iter()
            .filter(|&&row| self.calls.get(sample, row).is_some())
            .count()
    }

    fn markers(&mut self, class: MarkerClass) -> &mut Markers {
        match class {
            MarkerClass::Auto => &mut self.auto,
            MarkerClass::X => &mut self.x,
            MarkerClass::Y => &mut self.y,
        }
    }

    /// Drop a marker and record it in the SNP report.
    fn drop_marker(&mut self, class: MarkerClass, slot: usize, reason: &str) {
        let markers = self.markers(class);
        markers.alive[slot] = false;
        markers.logical -= 1;
        let row = markers.rows[slot];
        let id = self.variants[row].id.clone();
        let _ = writeln!(self.snp_file, "{id}\t{reason}");
    }

    /// Drop a sample and record it in the sample report.
    fn drop_sample(&mut self, sample: usize, reason: &str) {
        self.alive[sample] = false;
        let (fid, iid) = (&self.samples[sample].fid, &self.samples[sample].iid);
        let _ = writeln!(self.sample_file, "{fid}\t{iid}\t{reason}");
    }

    // -- the eight steps -------------------------------------------------

    fn pipeline(&mut self, opts: &Options) -> String {
        let callrate_m = callrate(opts, Opt::CallrateM);
        let callrate_n = callrate(opts, Opt::CallrateN);
        let loose = step1_threshold(callrate_m);
        let raw_snps = self.auto.raw + self.x.raw + self.y.raw + self.mt;
        let raw_samples = self.samples.len();

        let mut log = String::with_capacity(4096);
        let _ = writeln!(
            log,
            "Auto-QC starts at {}\n",
            console::ctime(console::now_local())
        );

        self.step1(&mut log, loose);
        self.step2(&mut log, callrate_n);
        self.step3(&mut log, callrate_m);
        if self.gender_qc() {
            self.step4(&mut log, loose, callrate_m);
            self.step5(&mut log);
            self.step6(&mut log, callrate_m);
        } else {
            // Each half of the gate reports itself: a map with X but no Y prints one line.
            if self.x.raw == 0 {
                log.push_str("X-Chr SNPs are not available. Gender QC is skipped.\n");
            }
            if self.y.raw == 0 {
                log.push_str("Y-Chr SNPs are not available. Gender QC is skipped.\n");
            }
        }
        self.step7(&mut log);
        self.step8(&mut log, opts, raw_samples, raw_snps);
        log
    }

    /// Whether steps 4–6 run at all.
    ///
    /// On the **map's** marker counts, not on what survives step 1: a fileset whose every
    /// X marker is monomorphic still runs the gender block, with an empty X panel.
    fn gender_qc(&self) -> bool {
        self.x.raw > 0 && self.y.raw > 0
    }

    /// Step 1 — the loose SNP call-rate filter, then monomorphic markers.
    ///
    /// The two are sequential, not parallel: the monomorphic scan only sees the markers
    /// the call-rate scan left behind. On `missing`, 629 of the 1569 low-call-rate markers
    /// are also monomorphic and the reference reports 3199 monomorphic, not 3828.
    fn step1(&mut self, log: &mut String, loose: f64) {
        let _ = writeln!(
            log,
            "Auto-QC step 1: Apply SNP call rate filter {} on {} SNPs (in {} samples)",
            percent(loose),
            self.auto.raw,
            self.samples.len()
        );
        let subjects = self.subjects();
        let reason = call_rate_reason(loose);
        let low_auto = self.filter_call_rate(MarkerClass::Auto, &subjects, loose, &reason);
        let _ = writeln!(
            log,
            "  {low_auto} autosome SNPs have call rate < {}",
            percent(loose)
        );
        let low_x = self.filter_call_rate(MarkerClass::X, &subjects, loose, &reason);
        let _ = writeln!(
            log,
            "  {low_x} X-chr SNPs have call rate < {}",
            percent(loose)
        );
        self.counts.low_call_rate = low_auto + low_x;

        // Monomorphic markers: autosomes, then X, then Y, in that order in the report.
        let mut mono = 0;
        for class in [MarkerClass::Auto, MarkerClass::X, MarkerClass::Y] {
            let live = self.markers(class).live();
            for (slot, row) in live {
                if self.monomorphic(row, &subjects) {
                    self.drop_marker(class, slot, "Monomorphic");
                    mono += 1;
                }
            }
        }
        let _ = writeln!(log, "  {mono} SNPs are monomorphic");
        self.counts.monomorphic = mono;
    }

    /// Remove every marker of `class` whose missing rate clears `threshold`.
    fn filter_call_rate(
        &mut self,
        class: MarkerClass,
        subjects: &[usize],
        threshold: f64,
        reason: &str,
    ) -> usize {
        let mut removed = 0;
        let live = self.markers(class).live();
        for (slot, row) in live {
            if self.missing_rate(row, subjects) > 1.0 - threshold {
                self.drop_marker(class, slot, reason);
                removed += 1;
            }
        }
        removed
    }

    /// Step 2 — the sample call-rate filter, over the surviving **autosomes** only.
    ///
    /// A sample missing every X marker but complete on the autosomes is not removed, and
    /// neither is one whose missingness all falls inside the dropped word.
    fn step2(&mut self, log: &mut String, callrate_n: f64) {
        let _ = writeln!(
            log,
            "\nAuto-QC step 2: Apply sample call rate filter {} on {} samples (with {} SNPs)",
            percent(callrate_n),
            self.samples.len(),
            self.auto.logical
        );
        let rows = self.auto.live_rows();
        let reason = format!("MissingMoreThan{}", int_percent(1.0 - callrate_n));
        let mut dropped = 0;
        for s in self.subjects() {
            let missing = rows.len() - self.call_count(s, &rows);
            if missing as f64 / rows.len() as f64 > 1.0 - callrate_n {
                self.drop_sample(s, &reason);
                dropped += 1;
            }
        }
        let _ = writeln!(
            log,
            "  {dropped} samples have call rate < {}",
            percent(callrate_n)
        );
        self.counts.sample_call_rate = dropped;
    }

    /// Step 3 — the SNP call-rate filter again, at `--callrateM`, over the samples left.
    fn step3(&mut self, log: &mut String, callrate_m: f64) {
        let subjects = self.subjects();
        let _ = writeln!(
            log,
            "\nAuto-QC step 3: Apply SNP call rate filter {} on {} SNPs (in {} samples)",
            percent(callrate_m),
            self.auto.logical,
            subjects.len()
        );
        let reason = call_rate_reason(callrate_m);
        let low_auto = self.filter_call_rate(MarkerClass::Auto, &subjects, callrate_m, &reason);
        let _ = writeln!(
            log,
            "  {low_auto} SNPs have call rate < {}",
            percent(callrate_m)
        );
        let low_x = self.filter_call_rate(MarkerClass::X, &subjects, callrate_m, &reason);
        let _ = writeln!(
            log,
            "  {low_x} chr-X SNPs have call rate < {}",
            percent(callrate_m)
        );
        self.counts.snp_call_rate = low_auto + low_x;
    }

    /// Step 4 — the three sex-chromosome filters at their loose thresholds.
    ///
    /// The header prints the class's **raw** Y count (300 on `sexchr`, where 205 Y markers
    /// had already gone as monomorphic); step 6's header prints the current one.
    fn step4(&mut self, log: &mut String, loose: f64, callrate_m: f64) {
        let _ = writeln!(
            log,
            "\nAuto-QC step 4: Apply call rate filters on {} Y-chr SNPs",
            self.y.raw
        );
        let (a, b, c) =
            self.sex_filters(log, '4', loose, X_HET_LOOSE, Y_IN_FEMALES_LOOSE, callrate_m);
        self.counts.y_in_men_loose = a;
        self.counts.x_het_loose = b;
        self.counts.y_in_women_loose = c;
    }

    /// Step 6 — the same three filters, tightened, over the samples step 5 left behind.
    fn step6(&mut self, log: &mut String, callrate_m: f64) {
        let _ = writeln!(
            log,
            "\nAuto-QC step 6: Apply call rate filters on {} Y-chr SNPs",
            self.y.logical
        );
        let (a, b, c) = self.sex_filters(
            log,
            '6',
            callrate_m,
            X_HET_TIGHT,
            Y_IN_FEMALES_TIGHT,
            callrate_m,
        );
        self.counts.y_in_men_tight = a;
        self.counts.x_het_tight = b;
        self.counts.y_in_women_tight = c;
    }

    /// The body shared by steps 4 and 6, which differ only in their thresholds and in the
    /// sample set they inherit.
    fn sex_filters(
        &mut self,
        log: &mut String,
        step: char,
        y_in_males: f64,
        x_het: f64,
        y_in_females: f64,
        callrate_m: f64,
    ) -> (usize, usize, usize) {
        let males = self.of_sex(MALE);
        let females = self.of_sex(FEMALE);
        // Both Y call-rate filters label their removals with `--callrateM`, step 4a
        // included, although it applies the looser step-1 threshold.
        let reason = call_rate_reason(callrate_m);

        let _ = writeln!(
            log,
            "\n  Step {step}a: Apply Y-chr call rate filter {} in males",
            percent(y_in_males)
        );
        let mut a = 0;
        let live = self.y.live();
        for (slot, row) in live {
            if self.missing_rate(row, &males) > 1.0 - y_in_males {
                self.drop_marker(MarkerClass::Y, slot, &reason);
                a += 1;
            }
        }
        let _ = writeln!(
            log,
            "  {a} chr-Y SNPs have call rate < {} in males",
            percent(y_in_males)
        );

        let _ = writeln!(
            log,
            "\n  Step {step}b: Apply X-chr heterozygosity filter {}% in males",
            int_percent(x_het)
        );
        let mut b = 0;
        let live = self.x.live();
        for (slot, row) in live {
            if self.het_rate(row, &males) > x_het {
                self.drop_marker(MarkerClass::X, slot, "xHeterozygosityInMale");
                b += 1;
            }
        }
        let _ = writeln!(
            log,
            "  {b} X-chr SNPs have heterozygosity > {}% in males",
            int_percent(x_het)
        );

        let _ = writeln!(
            log,
            "\n  Step {step}c: Apply Y-chr call rate filter {}% in females",
            int_percent(y_in_females)
        );
        let mut c = 0;
        let live = self.y.live();
        for (slot, row) in live {
            if 1.0 - self.missing_rate(row, &females) > y_in_females {
                self.drop_marker(MarkerClass::Y, slot, "YSNPInFemales");
                c += 1;
            }
        }
        let _ = writeln!(
            log,
            "  {c} chr-Y SNPs have call rate > {}% in females",
            int_percent(y_in_females)
        );
        (a, b, c)
    }

    /// Step 5 — the gender check.
    ///
    /// Every sample gets a Y-marker count, and the panel gets a cutoff of **half** the
    /// class's current marker count — a fixed fraction, not a fitted one: 95 Y markers
    /// give 47.5, 31 give 15.5 and 2 give 1.0, whatever the counts themselves look like.
    ///
    /// A reported sex is then wrong in one of two ways, each its own counter:
    ///
    /// * **mislabeled** — a male below the cutoff, a female above it;
    /// * **suspicious** — otherwise, a male below `0.7 x panel` or above the X-het cut, a
    ///   female above `0.3 x panel` or below it. Both fractions were bisected: with 20 Y
    ///   markers a male is clean at 14 and suspicious at 13; with 30, clean at 21 and
    ///   suspicious at 20. The same sweep on females flips at `0.3 x panel`.
    ///
    /// The X-het cut starts at 0.10, is lowered by any **reported female** whose own
    /// heterozygosity lies in `(0.05, 0.10)`, and is then truncated to whole percent.
    /// Males and unreported samples never move it, and a female at or below 0.05 does not
    /// either — she is simply flagged. Ten females at 0.09 print `set as 0.09` and none is
    /// flagged; ten at 0.05 print `0.10` and all ten are.
    fn step5(&mut self, log: &mut String) {
        let subjects = self.subjects();
        let _ = writeln!(
            log,
            "\nAuto-QC step 5: Gender QC on {} samples",
            subjects.len()
        );
        let panel = self.y.logical;
        let cutoff = panel as f64 / 2.0;
        let _ = writeln!(
            log,
            "\n  Step 5a: Determine thresholds in {panel} Y-chr SNPs for gender checking"
        );
        let _ = writeln!(
            log,
            "  {cutoff:.1} is used as the cutoff value for the chr-Y SNP count between males vs females"
        );

        let y_rows = self.y.live_rows();
        let x_rows = self.x.live_rows();
        let y_count: Vec<usize> = (0..self.samples.len())
            .map(|s| self.call_count(s, &y_rows) + self.y.phantom())
            .collect();
        let x_het: Vec<f64> = (0..self.samples.len())
            .map(|s| self.sample_het(s, &x_rows, self.x.phantom()))
            .collect();

        let mut cut = X_HET_CUT;
        for &s in &subjects {
            if self.sex[s] == FEMALE && x_het[s] > X_HET_FLOOR && x_het[s] < cut {
                cut = x_het[s];
            }
        }
        // ...and then truncated to whole percent, which the printed value gives away: ten
        // females at 3/40 print `filter 0.070`, not `0.075`. The floor test happens
        // before the truncation — a female at 0.055 does lower the cut, to 0.05.
        let cut = whole_percent(cut);

        let _ = writeln!(log, "\n  Step 5b: Gender checking");
        let _ = writeln!(log, "  X-chr heterozygosity is set as {cut:.2}");

        let mut mislabeled_male = Vec::new();
        let mut mislabeled_female = Vec::new();
        for &s in &subjects {
            match self.sex[s] {
                MALE if (y_count[s] as f64) < cutoff => mislabeled_male.push(s),
                FEMALE if (y_count[s] as f64) > cutoff => mislabeled_female.push(s),
                _ => {}
            }
        }
        let mut flagged = Vec::new();
        for &s in &subjects {
            if mislabeled_male.contains(&s) || mislabeled_female.contains(&s) {
                continue;
            }
            if suspicious(self.sex[s], y_count[s], panel, x_het[s], cut) {
                flagged.push(s);
            }
        }
        let unreported: Vec<usize> = subjects
            .iter()
            .copied()
            .filter(|&s| self.sex[s] == UNREPORTED)
            .collect();

        let _ = writeln!(
            log,
            "  {} (reported) males have gender errors (less than {cutoff:.1} Y-chr SNPs)",
            mislabeled_male.len()
        );
        let _ = writeln!(
            log,
            "  {} (reported) females have gender errors (more than {cutoff:.1} Y-chr SNPs)",
            mislabeled_female.len()
        );
        let _ = writeln!(
            log,
            "  {} samples have additional gender errors (according to X-Chr heterozygosity filter {cut:.3})",
            flagged.len()
        );
        let _ = writeln!(
            log,
            "  Genders of {} samples are not reported but will be inferred now",
            unreported.len()
        );

        self.counts.mislabeled_male = mislabeled_male.len();
        self.counts.mislabeled_female = mislabeled_female.len();
        self.counts.suspicious = flagged.len();
        // The three groups are appended in check order, not in `.fam` order: a suspicious
        // sample listed first in the `.fam` still follows the mislabeled ones.
        for s in mislabeled_male {
            self.drop_sample(s, "MislabeledAsMale");
        }
        for s in mislabeled_female {
            self.drop_sample(s, "MislabeledAsFemale");
        }
        for s in flagged {
            self.drop_sample(s, "GenderQC");
        }
        // Inference is by Y-marker count alone — X heterozygosity has no say, which an
        // unreported sample with a female-looking X and a full Y panel shows — and the
        // test is **not strict**: with every Y marker gone the cutoff is 0.0 and the
        // reference infers male.
        //
        // The inferred sex goes into the report and **nowhere else**: step 6 still counts
        // only the reported males and females, so four samples inferred here as two of
        // each leave 6b and 6c reporting zero.
        for s in unreported {
            let sex = if y_count[s] as f64 >= cutoff {
                MALE
            } else {
                FEMALE
            };
            let (fid, iid) = (&self.samples[s].fid, &self.samples[s].iid);
            let _ = writeln!(self.updatesex, "{fid}\t{iid}\t{sex}");
        }
    }

    /// Step 7 — the final counts. X and Y are named only when they have markers left, and
    /// MT never is.
    fn step7(&mut self, log: &mut String) {
        let _ = writeln!(log, "\nAuto-QC step 7: Final check");
        let _ = write!(
            log,
            "  {} samples, {} autosome SNPs",
            self.subjects().len(),
            self.auto.visible()
        );
        if self.x.visible() > 0 {
            let _ = write!(log, ", {} X-chr SNPs", self.x.visible());
        }
        if self.y.visible() > 0 {
            let _ = write!(log, ", {} Y-chr SNPs", self.y.visible());
        }
        log.push('\n');
    }

    /// Step 8 — the summary table, which goes to the console and to the summary file
    /// verbatim, then the "saved in" lines and the closing timestamp.
    fn step8(&mut self, log: &mut String, opts: &Options, raw_samples: usize, raw_snps: usize) {
        let _ = writeln!(log, "\nAuto-QC step 8: QC Summary Report\n");

        let c = &self.counts;
        let mut table = String::with_capacity(1024);
        table.push_str(&count_row("Step", "Description", "Subjects", "SNPs"));
        table.push_str(&count_row(
            "1",
            "Raw data counts",
            &raw_samples.to_string(),
            &raw_snps.to_string(),
        ));
        table.push_str(&snp_row(
            "1.1",
            "SNPs with very low call rate < 80% (removed)",
            c.low_call_rate,
        ));
        table.push_str(&snp_row("1.2", "Monomorphic SNPs (removed)", c.monomorphic));
        table.push_str(&subject_row(
            "1.3",
            "Sample call rate < 95% (removed)",
            c.sample_call_rate,
        ));
        table.push_str(&snp_row(
            "1.4",
            "SNPs with call rate < 95% (removed)",
            c.snp_call_rate,
        ));
        if self.gender_qc() {
            // The step-2 row is the raw counts less what steps 1–3 removed, dropped word
            // included — it is arithmetic on the counters, not a recount.
            let left = raw_snps - c.low_call_rate - c.monomorphic - c.snp_call_rate;
            table.push_str(&count_row(
                "2",
                "data counts for gender error checking",
                &(raw_samples - c.sample_call_rate).to_string(),
                &left.to_string(),
            ));
            table.push_str(&snp_row(
                "2.1",
                "Y-chr SNPs with call rate < 80% in men (removed)",
                c.y_in_men_loose,
            ));
            table.push_str(&snp_row(
                "2.2",
                "X-chr SNPs with heterozygosity > 5% in men (removed)",
                c.x_het_loose,
            ));
            table.push_str(&snp_row(
                "2.3",
                "Y-chr SNPs with genotypes in >10% women (removed)",
                c.y_in_women_loose,
            ));
            table.push_str(&subject_row(
                "2.4",
                "Mislabeled as male (removed)",
                c.mislabeled_male,
            ));
            table.push_str(&subject_row(
                "2.5",
                "Mislabeled as female (removed)",
                c.mislabeled_female,
            ));
            table.push_str(&subject_row(
                "2.6",
                "Suspicious gender error (removed)",
                c.suspicious,
            ));
            table.push_str(&snp_row(
                "2.7",
                "Y-chr SNPs with call rate < 95% in men (removed)",
                c.y_in_men_tight,
            ));
            table.push_str(&snp_row(
                "2.8",
                "X-chr SNPs with heterozygosity > 1% in men (removed)",
                c.x_het_tight,
            ));
            table.push_str(&snp_row(
                "2.9",
                "Y-chr SNPs with genotypes in >2% women (removed)",
                c.y_in_women_tight,
            ));
        }
        table.push_str(&format!("{:<5}{:<55}\n", "3", "Generate Final Study Files"));
        // The final SNP count is a recount rather than raw-minus-counters, because the
        // dropped word is in neither the raw count nor any counter.
        let final_snps = self.auto.visible() + self.x.visible() + self.y.visible() + self.mt;
        table.push_str(&count_row(
            "",
            "Final QC'ed data",
            &self.subjects().len().to_string(),
            &final_snps.to_string(),
        ));
        log.push_str(&table);

        let summary = out_path(opts, "_autoQC_Summary.txt");
        let snps = out_path(opts, "_autoQC_snptoberemoved.txt");
        let samples = out_path(opts, "_autoQC_sampletoberemoved.txt");
        write_file(&summary, &table);
        write_file(&snps, &self.snp_file);
        write_file(&samples, &self.sample_file);
        let _ = writeln!(log, "\nQC summary report saved in {summary}");
        let _ = writeln!(log, "SNP-removal QC file saved in {snps}");
        let _ = writeln!(log, "Sample-removal QC file saved in {samples}");
        if !self.updatesex.is_empty() {
            let path = out_path(opts, "_autoQC_updatesex.txt");
            write_file(&path, &self.updatesex);
            let _ = writeln!(log, "Update-sex QC file saved in {path}");
        }
        let _ = writeln!(
            log,
            "\nAuto-QC ends at {}",
            console::ctime(console::now_local())
        );
    }
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

/// One of the two call-rate options, defaulted when the command line was silent.
fn callrate(opts: &Options, o: Opt) -> f64 {
    if opts.was_given(o) {
        opts.double(o)
    } else {
        CALLRATE_DEFAULT
    }
}

/// The step-1 SNP threshold derived from `--callrateM`.
///
/// `min(0.8, 0.1 * (trunc(callrateM * 10) - 1))` — truncate toward zero, then cap at 80%.
/// Swept against the reference: 0.85 gives 70% (not 75%), 0.9 and everything above give
/// 80%, and −0.5 gives −60%, which is what pins the truncation as toward zero rather than
/// downward.
///
/// The **multiplication** matters: `0.1 * 7` is `0.7000000000000001` while `7.0 / 10.0` is
/// `0.6999999999999999`, and the two disagree on a Y marker called in exactly 70% of the
/// males — the reference removes it, so it is holding the larger of the two.
fn step1_threshold(callrate_m: f64) -> f64 {
    (0.1 * ((callrate_m * 10.0).trunc() - 1.0)).min(STEP1_CEILING)
}

/// Whether a sample with a reported sex fails the second-tier gender check.
///
/// Two tests, either of which flags: the sample's share of the Y panel has to be strictly
/// above `2/3` for a male and strictly below `1/3` for a female, and its X heterozygosity
/// has to be at most the cut for a male and at least the cut for a female. `share` is a
/// division, so an emptied Y panel makes it NaN and both comparisons false — which is
/// exactly what the reference does when step 4 has removed every Y marker. An X
/// heterozygosity computed over nothing is NaN for the same reason.
fn suspicious(sex: u8, y_count: usize, panel: usize, x_het: f64, cut: f64) -> bool {
    let share = y_count as f64 / panel as f64;
    match sex {
        MALE => share <= Y_MALE_SHARE || x_het > cut,
        FEMALE => share >= Y_FEMALE_SHARE || x_het < cut,
        _ => false,
    }
}

/// A heterozygosity truncated to whole percent, the form the gender check stores its cut
/// in: 0.075 becomes 0.07, not 0.08.
fn whole_percent(x: f64) -> f64 {
    (x * 100.0).trunc() / 100.0
}

/// A threshold as the console spells it: `%.1f%%` of the percentage.
fn percent(threshold: f64) -> String {
    format!("{:.1}%", threshold * 100.0)
}

/// A threshold as the reason strings and the 4b/4c headers spell it: a **rounded** integer
/// percentage. `(int)(0.85 * 100)` would be 84; the reference prints 85.
fn int_percent(threshold: f64) -> i64 {
    (threshold * 100.0).round() as i64
}

/// The reason string a call-rate filter writes into the SNP report.
fn call_rate_reason(threshold: f64) -> String {
    format!("CallRateLessThan{}", int_percent(threshold))
}

/// A summary row carrying both counts: `%-5s%-55s%-10s%-10s`.
fn count_row(step: &str, description: &str, subjects: &str, snps: &str) -> String {
    format!("{step:<5}{description:<55}{subjects:<10}{snps:<10}\n")
}

/// A summary row whose count belongs in the SNPs column: `%-5s%-65s(%d)`.
fn snp_row(step: &str, description: &str, n: usize) -> String {
    format!("{step:<5}{description:<65}({n})\n")
}

/// A summary row whose count belongs in the Subjects column: `%-5s%-55s(%d)`.
fn subject_row(step: &str, description: &str, n: usize) -> String {
    format!("{step:<5}{description:<55}({n})\n")
}

fn write_file(path: &str, text: &str) {
    let _ = std::fs::write(path, text);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn step1_threshold_matches_the_sweep() {
        for (given, want) in [
            (0.5, 0.4),
            (0.6, 0.5),
            (0.7, 0.6),
            (0.8, 0.7),
            (0.85, 0.7),
            (0.9, 0.8),
            (0.95, 0.8),
            (0.99, 0.8),
            (1.0, 0.8),
            (1.2, 0.8),
            (-0.5, -0.6),
        ] {
            assert_eq!(percent(step1_threshold(given)), percent(want), "{given}");
        }
        assert_eq!(percent(step1_threshold(0.0)), "-10.0%");
    }

    #[test]
    fn reason_strings_carry_a_rounded_percentage() {
        assert_eq!(call_rate_reason(0.95), "CallRateLessThan95");
        assert_eq!(call_rate_reason(0.85), "CallRateLessThan85");
        assert_eq!(call_rate_reason(0.87), "CallRateLessThan87");
        assert_eq!(call_rate_reason(0.8), "CallRateLessThan80");
        assert_eq!(int_percent(1.0 - 0.95), 5);
        assert_eq!(int_percent(1.0 - 0.96), 4);
    }

    /// Every summary row is 80 columns wide up to its count.
    #[test]
    fn summary_rows_are_column_aligned() {
        assert_eq!(
            count_row("Step", "Description", "Subjects", "SNPs").len(),
            81
        );
        assert_eq!(count_row("1", "Raw data counts", "3", "5000").len(), 81);
        assert_eq!(
            snp_row("1.1", "SNPs with very low call rate < 80% (removed)", 0).len(),
            74
        );
        assert_eq!(
            subject_row("1.3", "Sample call rate < 95% (removed)", 0).len(),
            64
        );
    }

    /// Both gender-check boundaries, as bisected against the reference.
    #[test]
    fn gender_suspicion_boundaries() {
        let male = |count, panel| suspicious(MALE, count, panel, 0.0, 0.10);
        // (panel, first clean count) — everything below it is suspicious.
        for (panel, clean) in [(2, 2), (20, 14), (30, 21), (32, 22), (48, 33), (64, 43)] {
            assert!(male(clean - 1, panel), "{panel}/{clean}");
            assert!(!male(clean, panel), "{panel}/{clean}");
        }
        let female = |count, panel| suspicious(FEMALE, count, panel, 0.5, 0.10);
        for (panel, flagged) in [(20, 7), (48, 16), (100, 34)] {
            assert!(female(flagged, panel), "{panel}/{flagged}");
            assert!(!female(flagged - 1, panel), "{panel}/{flagged}");
        }
        // An emptied Y panel flags nobody: 0/0 is NaN and every comparison against it is
        // false. Neither does an X heterozygosity computed over no markers.
        assert!(!suspicious(MALE, 0, 0, 0.0, 0.10));
        assert!(!suspicious(FEMALE, 0, 0, 0.5, 0.10));
        assert!(!suspicious(FEMALE, 0, 20, f64::NAN, 0.10));
        // X heterozygosity alone is enough, in either direction.
        assert!(suspicious(MALE, 20, 20, 0.11, 0.10));
        assert!(!suspicious(MALE, 20, 20, 0.10, 0.10));
        assert!(suspicious(FEMALE, 0, 20, 0.09, 0.10));
        assert!(!suspicious(FEMALE, 0, 20, 0.10, 0.10));
        // An unreported sex is inferred, never flagged.
        assert!(!suspicious(UNREPORTED, 0, 20, 0.0, 0.10));
    }

    #[test]
    fn the_het_cut_is_truncated_to_whole_percent() {
        assert_eq!(whole_percent(0.075), 0.07);
        assert_eq!(whole_percent(0.09375), 0.09);
        assert_eq!(whole_percent(0.0625), 0.06);
        assert_eq!(whole_percent(0.055), 0.05);
        assert_eq!(whole_percent(0.10), 0.10);
    }

    /// The dropped word: a class whose marker count is a multiple of 16 loses its last 16.
    #[test]
    fn whole_word_classes_lose_their_last_word() {
        assert_eq!(Markers::new((0..32).collect()).rows.len(), 16);
        assert_eq!(Markers::new((0..33).collect()).rows.len(), 33);
        assert_eq!(Markers::new((0..16).collect()).rows.len(), 0);
        assert_eq!(Markers::new((0..15).collect()).rows.len(), 15);
        assert_eq!(Markers::new((0..50_000).collect()).rows.len(), 49_984);
        assert_eq!(Markers::new(Vec::new()).rows.len(), 0);
        // ...while the count the headers print is untouched by it.
        assert_eq!(Markers::new((0..32).collect()).logical, 32);
        // ...and step 5's per-sample loops count its 16 slots as called.
        assert_eq!(Markers::new((0..32).collect()).phantom(), 16);
        assert_eq!(Markers::new((0..33).collect()).phantom(), 0);
        assert_eq!(Markers::new(Vec::new()).phantom(), 0);
    }
}
