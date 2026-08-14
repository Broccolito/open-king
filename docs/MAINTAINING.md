# Maintaining open-king

Everything you need to change this project safely: the one rule that cannot be broken, how
the tree is laid out, how the test corpus and its golden output are produced, how to run and
extend the parity suite, **the fixture technique the whole segment engine rests on (§8)**, and
how to add an analysis.

Read `docs/PARITY.md` first if you want to know what currently works.

If you are picking this up cold and intend to work on the IBD-segment caller — which is where
all the remaining parity gap lives — read in this order: `docs/PARITY.md` §5.0 (what is
solved, what is not, what to run next), then §8 below (the instrument), then
`docs/research/17-seg-caller.md` and `18-ibd1-caller.md`.

---

## 1. The clean-room constraint

**Never read, download, decompile or copy KING's source code.** Not to check a formula,
not to settle an argument, not "just to look".

open-king exists so that a permissively licensed implementation of KING's file formats and
estimators can be embedded in other software. That is only defensible if the implementation
is genuinely independent. The moment a contributor has read the original source, everything
they touch afterwards is derivative, and the fact is not recoverable after the fact — a
reviewer cannot tell by looking at a diff. So the rule is absolute rather than
best-effort, and it is stated at the top of `crates/king-core/src/ibdseg.rs` as well as
here.

**What you may use instead**, and what everything currently in the tree came from:

1. **Peer-reviewed publications** — Manichaikul *et al.* 2010 for the kinship and IBS
   estimators. Cite the paper where you implement a formula.
2. **The public manual and website**, including whatever the binary itself emits: the
   `--help` text, the banner, the error strings, and the R scripts `--rplot` writes
   (`<prefix>_ibd1vsibd2.R` is where the relationship-degree cut-points are stated in the
   binary's own words).
3. **Black-box observation of the reference binary.** Run it on inputs you constructed,
   and read the outputs. This is the workhorse, and it is what §5 is about.

The IBD-segment algorithm is not published at all — KING's manual says the manuscript is
"yet to be published" — so every rule in `king-core::ibdseg` was established the third way.
`docs/research/` is the log of those experiments; each rule's doc comment names the
experiment that fixed it.

**Writing it down.** When you establish a rule, record *how you know*, not just what the
rule is. A constant with no provenance is indistinguishable from a guess, and the next
person will not be able to tell whether changing it is safe. The existing doc comments are
the house style: state the rule, then the measurement that pins it, then the alternative
that was ruled out.

**Do not fit to dataset names.** If a conditional anywhere in `crates/*/src/` mentions
`bigish` or `nuclear` or any other corpus dataset, the rule it guards is wrong — you have
encoded the corpus rather than the reference's behaviour. Branch on properties of the data
instead. (Dataset names appear legitimately in `crates/*/tests/`, as the list a scorecard
iterates over.)

---

## 2. Layout

```
crates/
  king-io/        PLINK 1 fileset I/O: .bed / .bim / .fam -> packed bit planes.
                  lib.rs states the types; the submodules implement them.
  king-core/      The estimators, with no I/O and no console.
                  counts.rs   the per-pair counting kernel
                  kinship.rs  the kinship / IBS estimators
                  ibdseg.rs   the IBD-segment engine (usable segments, per-pair
                              calling, aggregation, InfType, the --degree filter)
                  infer.rs    relationship inference from the estimates
  king-cli/       The command line, the console, and one module per analysis.
                  cli.rs      option table + parser (the banner depends on it)
                  console.rs  every line the program prints, C-compatible formatting
                  load.rs     fileset loading, chromosome partition, progress ticks
                  analysis/   one file per --flag; mod.rs owns only what is shared
                  main.rs     dispatch: which passes run, in which order

docs/
  PARITY.md              what is byte-identical and what is not — authoritative
  MAINTAINING.md         this file
  SPEC.md                the reference's observable behaviour, flag by flag
  BEHAVIOR.md            raw sweeps behind the rules
  VERIFIED_FORMULAS.md   the estimators and the experiment that fixed each
  research/              the investigation log, numbered in the order it happened;
                         18-ibd1-caller.md is the newest; it and 17-seg-caller.md
                         (whose §14 is the bridge, landed) cover the .seg caller
  research/fixtures/     the fixture rigs (see §8 — read it before touching the
                         segment caller): filesets whose answer is forced by
                         construction, used to pin constants the corpus cannot see.
                         fixlab.py builds a fileset and drives the reference;
                         gate8.py brackets the --degree 1 IBD2 clause;
                         segcanvas.py is the .seg canvas (+ its measured cache);
                         ibd1canvas.py is the same canvas built IBD1-side up;
                         gradebinary.py replays those canvases with OUR binary and
                         grades it against the cached reference answers;
                         avfs.py / avfs_score.py are the --build AV.FS rig

tests/parity/
  generate_corpus.py     builds the 13 input datasets from a fixed seed
  make_alt_inputs.py     builds the alternate --fam / --bim inputs
  capture_params.py      captures the `params` golden group from the reference
  run_parity.py          the acceptance test: replay every capture, diff the bytes
  measure_gaps.py        how *big* each remaining difference is
  verify_formulas.py     checks docs/VERIFIED_FORMULAS.md against captured output
  verify_row_order.py    checks the .kin0 / .ibs0 block-tiled row ordering
  golden/<group>/<case>/ one captured reference invocation per directory
  work/                  generated inputs (gitignored) — data/ and alt/
  fit/                   analysis scripts used while fitting rules
  probes/                small one-off reference probes worth keeping
```

**One file in `fit/` is not a scratch script.** `fit/engine.py` is a line-for-line Python
mirror of `crates/king-core/src/ibdseg.rs` with every disputed rule exposed as a `Params`
field, so a candidate rule can be scored over the whole corpus in a second instead of a
rebuild-and-replay. It is a *mirror, not a second source of truth*: `fit/check_mirror.py`
asserts that with default `Params` it reproduces the built binary's own `.seg` columns and
`MaxIBD2` on every corpus row. **If you change a rule in `ibdseg.rs`, either update
`engine.py` to match or expect `check_mirror.py` to fail** — and when the two disagree,
the Rust is right and the mirror has the bug. `fit/seg17.py` and `fit/seg18.py` are the
scorecards built on it: each prints the committed rule and the retired one side by side over
all 982 primary rows — `seg17.py` for the `.seg` IBD2 caller, `seg18.py` for the `IBD1Seg`
overlap rule at 3, 5 and 10 Mb — and `R17(...)` / `R18(...)` expose every knob for a
candidate. `seg17.py grid19` and `seg18.py grid` sweep them all.

**Two committed files are measurement caches, and a non-reference binary will silently
corrupt them.** `docs/research/fixtures/segcanvas_measured.json` (6 416 answers) and
`ibd1canvas_measured.json` (1 013) hold readings *measured from the reference*, which is what
lets both rigs re-run in under a second without it. `segcanvas.py` writes whatever it measures
back into its cache and takes the binary from `$KING`, so **never point `$KING` at our build
while running `segcanvas.py` or `ibd1canvas.py` in the tree.** To grade our own build, use
`gradebinary.py` (§8), which reads those caches, never writes them, and keeps its own answers
in `$TMPDIR`. `git diff --stat docs/research/fixtures/` before committing; those two files
should only ever change when the reference was the thing being run.

A second hazard in that rig, worth knowing before you trust a single probe: the reference
has a **major-allele QC check seeded from the clock**, which aborts a run with
`FATAL ERROR - Too many first alleles as the major allele` at random on small constructed
filesets. `segcanvas.py` retries up to 24 times with a sleep between attempts for exactly
this reason. A one-shot probe of a fixture is not evidence — run it a dozen times and count.

Three rules of thumb the tree already follows:

* **`lib.rs` and `mod.rs` state contracts**; implementations live in submodules. If you
  are adding behaviour, it usually belongs in a submodule, not in the module that
  declares the shape.
* **`king-core` never prints and never touches the filesystem.** Formatting a number the
  way C's `printf` would is `king-cli::console`'s job.
* **Anything used by exactly one analysis lives in that analysis's module.**
  `analysis/mod.rs` is for what two or more of them share.

All four crate roots carry `#![forbid(unsafe_code)]`. Keep it that way.

---

## 3. Build and check

```bash
cargo build --release                                   # -> target/release/king
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check
```

CI (`.github/workflows/ci.yml`) runs exactly these four on Linux, macOS and Windows. All
four must be clean before the parity suite is worth looking at.

CI does **not** run the parity suite, and there are two separate reasons — only the first
of which is usually remembered:

1. **It cannot gate as-is.** `run_parity.py` exits 1 whenever any case fails, and the suite
   sits at 408/480 by design, so a bare invocation would fail every build. Gating it needs
   a regression *floor* (`--json`, then assert the pass count has not dropped), not a
   pass/fail.
2. **It has only ever been run on macOS.** The goldens were captured there. Nothing in them
   is obviously host-dependent — timestamps, thread counts, progress tokens and absolute
   paths are all normalized (`docs/PARITY.md` §7) — but that is an argument, not a
   measurement. Run the suite on Linux and confirm 408/480 *before* wiring it into CI, and
   restrict the step to one OS until Windows is checked too (the binary is `king.exe`
   there).

The goldens are committed and the corpus regenerates from a seed, so neither reason is
about missing inputs; only the 480/480 self-check needs the reference binary.

### Before tagging a release

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace                                             # 290 passed, 1 ignored
cargo clean && cargo build --release --offline   # must work from a clean checkout, offline
python3 tests/parity/run_parity.py --impl target/release/king      # record the exact count
python3 tests/parity/run_parity.py --impl "<reference>"            # must be 480/480
python3 tests/parity/measure_gaps.py --impl target/release/king -q # the numbers PARITY.md §4 quotes
KING_GOLDEN=tests/parity/golden \
  cargo test -p king-core --test ibdseg_parity -- --nocapture      # the row-level .seg scorecard
cd tests/parity/fit && python3 check_mirror.py                     # must print MIRROR OK
python3 docs/research/fixtures/gradebinary.py target/release/king          # 6000/6000
python3 docs/research/fixtures/gradebinary.py target/release/king --ibd1   # 540/540 closed
git ls-files | grep -E '\.(bed|bim|fam|vcf|bcf)$'                  # must print nothing
git diff --stat docs/research/fixtures/*_measured.json             # the caches must be untouched
```

The last release measured **408 PASS / 72 FAIL / 480**, self-check **480/480**, 290 tests
passing, a clean offline build in **8.07 s** from a pristine copy of the tree, and that
clean-tree binary re-measured at 408/480. Do not publish a count you have not just re-run: the
parity number is the project's entire claim, and it is cheap to check (the suite takes about
two seconds).

**Publish both counts, not just the headline.** The 408 is a *whole-file* number — a case
turns `PASS` only when every row of every file it writes is byte-exact — and the residual is
spread thinly enough that a large row-level gain can move it by nothing. `docs/PARITY.md` §4.4
is the row-level scoreboard and §3 the file-level one; a release note that quotes one without
the other misleads in one direction or the other. The `17-seg-caller.md` §14 correction moved
**zero** corpus rows and zero cases and was still right (§8 below is how that was shown).

The last one is not decoration. `.gitignore` excludes `/docs/research/fixtures/work/` and
`/tests/parity/work/`, but **`.gitignore` does not untrack files already committed** — a
fixture fileset added before the rule stays in the index forever until someone runs
`git rm --cached`. Fixture filesets are regenerated by the rigs in
`docs/research/fixtures/`; none of them is evidence that needs preserving, and none should
ever be committed. Keep `docs/PARITY.md` and `README.md` in step with the recorded count:
the README must claim nothing `PARITY.md` does not measure.

Two test suites want the goldens and skip themselves without them:

```bash
KING_GOLDEN=tests/parity/golden cargo test -p king-core --test ibdseg_parity -- --nocapture
```

prints the per-dataset `.seg` scorecard rather than just passing.

---

## 4. Regenerating the input corpus

The datasets are **not committed**; the seeded generator is. Anything that needs them
regenerates them automatically, so normally you do nothing. To do it by hand:

```bash
python3 tests/parity/generate_corpus.py --outdir tests/parity/work/data
python3 tests/parity/make_alt_inputs.py --datadir tests/parity/work/data \
                                        --outdir  tests/parity/work/alt
```

About 20 seconds for all 13 datasets. `--only <name> ...` regenerates a subset;
`--seed` changes the seed, which **invalidates every golden** — do not.

The 13 datasets are deliberately varied and each exists to exercise something:
`trio`, `nuclear`, `threegen`, `multifam` (pedigree shapes), `dups` (a duplicate pair),
`missing` (per-sample and per-SNP missingness), `monomorphic` (monomorphic, ultra-rare and
missing-allele markers), `sexchr` (X, Y, XY and MT markers), `unrelated`, `admixed`
(two-population structure), `singleton`, `pair` (degenerate sizes), `bigish` (200 samples —
the only one large enough for the between-family and screening paths).

`generate_corpus.py` is pure and seeded: rerunning it produces byte-identical files, which
is what makes the goldens meaningful. If you change it, every capture that touches the
dataset you changed must be re-captured (§5).

---

## 5. Capturing golden output from the reference

A golden case is a directory under `tests/parity/golden/<group>/<case>/`:

```
cmd.txt        argv, one token per line or one shell line, with {KING} {DATA} {ALT}
exitcode.txt   the reference's exit status
stdout.txt     verbatim, including the \r progress tokens
stderr.txt     verbatim
<everything else>   every file the reference wrote into its working directory
```

The four groups are `core` (one analysis per dataset), `apps` (the report and clustering
passes), `ibdseg` (`--ibdseg` and its parameters) and `params` (flag plumbing and error
probes, dataset-independent).

**Re-capturing existing cases** — the normal path. `--update` replays each case with the
*reference* and overwrites the capture:

```bash
python3 tests/parity/run_parity.py --impl target/release/king \
        --ref "/path/to/reference/king" --update --filter ibdseg/
```

`--update` requires `--ref` and honours `--filter`, so you can re-capture one group. It
rewrites `stdout.txt`, `stderr.txt`, `exitcode.txt` and every produced file, and deletes
golden files the reference no longer writes.

**Adding new cases.** Create the directory, write `cmd.txt` with the placeholders, then
run `--update --filter <that case>` to fill in everything else. `capture_params.py` is the
worked example: it builds a `params` case per (dataset × flag combination), runs the
reference with `cwd` set to the case directory, and rewrites absolute input paths to
`{DATA}` / `{ALT}` so the capture is replayable on another machine.

**Three things to get right when capturing:**

1. **Run with `--cpus 1` whenever the case can write `<prefix>X.kin0`.** That writer races
   in the reference and produces a different file every time (`docs/PARITY.md` §5.2).
   Captures made without it cannot be diffed and the harness excludes them.
2. **Never hand-edit a capture.** If it looks wrong, it is telling you something about the
   reference. The `kingX.seg` whose header names 11 columns while its rows carry 9 values
   is real, and "fixing" it would have hidden a genuine finding.
3. **A capture is evidence, not a target to satisfy.** If a new case fails, the first
   question is what the reference is doing, not how to make the diff go away.

---

## 6. Running the parity suite

```bash
# everything
python3 tests/parity/run_parity.py --impl target/release/king

# one group or one case, with diffs
python3 tests/parity/run_parity.py --impl target/release/king --filter ibdseg/ -v
python3 tests/parity/run_parity.py --impl target/release/king \
        --filter core/bigish__related_degree3 -v

# prove the harness itself is sound: must be 480/480
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"

# how big the differences are, not just that they exist
python3 tests/parity/measure_gaps.py --impl target/release/king -q
python3 tests/parity/measure_gaps.py --impl target/release/king --by-dataset king.seg
```

Both scripts are Python 3 standard library only. Exit status of `run_parity.py`: 0 all
passed, 1 at least one failed, 2 harness error. Useful flags: `--jobs N`, `--timeout SEC`,
`--max-diff-lines N`, `--keep` (keep the per-case temp directories), `--json OUT`,
`--include-analysis` (also run the 10 `core/_analysis/` captures, which are kept for
analysis rather than as targets).

**The self-check is not optional.** `--impl <reference>` must be 480/480. If it is not,
the normalization in `run_parity.py` is either incomplete (a genuinely non-reproducible
line is being diffed) or too aggressive (it is hiding a real difference), and every other
number the harness prints is unreliable until that is fixed.

**Interpreting a failure.** The note after `FAIL` names the file and how it differs:
`!=(bytes)` is a raw byte difference, `!=(num)` means the harness found the difference to
be numeric and will show you a column-wise report with `-v`, `missing:<file>` means we did
not write a file the reference wrote. `measure_gaps.py` turns any of these into rows and
columns and errors.

---

## 7. Adding an analysis

Working backwards from the output, which is the order that keeps you honest:

1. **Capture first.** Add golden cases for the new flag across all 13 datasets (and
   whatever parameter combinations matter) *before* writing any Rust. You cannot fit to
   evidence you have not collected, and the captures will contradict your assumptions —
   that is what they are for.
2. **Write the module.** `crates/king-cli/src/analysis/<flag>.rs`, declared in
   `analysis/mod.rs`. It owns its output files and its console body. Domain logic —
   anything that computes rather than prints — belongs in `king-core`.
3. **Wire the flag.** Add it to the option table in `cli.rs` if it is not already there
   (many out-of-scope flags are already parsed so the banner stays exact), and to the
   `ANALYSES` list so `Options in effect:` includes it.
4. **Wire the dispatch.** `main.rs` has three groups: passes that print the loader
   preamble, passes that do not, and `--ibdseg`, which has its own small-sample downgrade.
   Put yours in the right one — the console layout differs between them, and the blank
   line before `Options in effect:` comes from a different place in each.
5. **Match the console exactly.** Line order, the indent of the `… ends at` timestamp,
   the wording of every announced filename. `console.rs` has the C-compatible `%f` / `%g`
   / `%e` formatting; use it rather than Rust's `{}`, which differs on edge values.
6. **Run the suite.** Then `measure_gaps.py`, which will tell you *how* wrong you are on
   anything still failing.
7. **Update `docs/PARITY.md`** — the matrix, the per-file table, and the measured size of
   whatever gap you left. It is the project's public claim; if it overstates, that is worse
   than the gap it is hiding.

---

## 8. The fixture rigs, and the canvas technique

This section is the one a newcomer cannot reconstruct from the code. The IBD-segment
algorithm is unpublished, the corpus is 13 noisy datasets, and **every constant in
`king-core::ibdseg` was measured rather than derived**. These are the instruments that made
that possible, in increasing order of power.

### 8.1 Construct the answer — `fixlab.py`

`docs/research/fixtures/fixlab.py` builds a PLINK1 fileset in which one designated pair's IBD
state is exact **by construction**: explicit shared haplotypes, forced IBS0 / IBS1 markers at
named positions, and if you want, every genotype of every sample written by hand. A threshold
measured in "number of markers of kind X" then becomes a threshold in the statistic itself,
with no sampling noise to fit through. This is how `MIN_INFORMATIVE = 10` was pinned to the
unit, and how the `--degree 1` IBD2 clause was found — the corpus can see neither.

`Fixture(name, chroms, nsample, maf, seed)` then `.build(dir)` writes the `.bed`/`.bim`/`.fam`;
`probe()` runs the reference and parses back `.seg` / `allsegs.txt`. Samples 0 and 1 are the
pair; the rest are padding, present only because `--ibdseg` downgrades below 5 samples.

### 8.2 Move one thing

Deleting the first *m* markers of a fileset shifts the global 64-marker word grid by *m* and
changes nothing else, so genotypes and segment lengths stay fixed while the reference's
verdict moves. 512 such invocations validated the informativeness gate against data that had
no part in choosing it.

### 8.3 The canvas — the project's key instrument

**The problem it solves.** The reference prints `IBD1Seg` / `IBD2Seg` as a *proportion of the
genome, to four decimals*. That is far too coarse to see a rule that moves a segment boundary
by one marker — which is exactly the resolution at which the caller's remaining disagreements
live. Worse, the printed number confounds three unknowns at once: which runs were called, how
many calls each run produced, and where each call's endpoints landed.

**The trick.** Make the printed column *a ruler*. A canvas is a two-chromosome fileset:

* **chr1** carries a small fixed carrier segment, present only so the pair earns a `.seg` row
  at all (a pair with nothing called is not printed, so there would be nothing to read).
* **chr2** is the canvas proper. It is painted **one complete 64-marker word at a time** —
  every word is a single named composition (`n` HetHet, `n` opposite homozygotes, `n`
  het-vs-hom mismatches, …) — and the region is bounded by **all-IBS0 walls**, words in which
  every marker is an opposite homozygote, which no caller will ever cross.
* **the spacing is chosen** so that the usable total `D` sits just above the reference's
  100 Mb floor and **one unit in the last printed place of the column is a known small
  fraction of one marker gap** (about a fifth for `segcanvas.py`, a ninth for `ibd1canvas.py`).

`segcanvas.mk(cv, res, what)` then multiplies the printed proportion back out by `D`,
subtracts the carrier's fixed contribution and divides by the spacing, and what comes out is
**the number of marker intervals called on chr2** — an integer, recoverable exactly despite
four-decimal printing.

**Decoding calls from that integer.** Calls inside one usable segment are emitted adjacent:
each new call starts one marker past the previous call's end. So `c` calls covering `w` whole
words measure `64w − c` markers, and

```
c = (−M) mod 64          words = (M + c) / 64
```

recovers **both** the number of calls and the number of words from the single printed number.
`segcanvas.wc(m)` does this. That is the whole instrument: a four-decimal proportion turned
into an exact readout of what the caller did, word by word.

**Using it.** Write a word sequence as a string over the rig's alphabet and read the answer:

```bash
cd docs/research/fixtures
python3 segcanvas.py            # the whole 17-seg-caller.md document, from cache, ~1 s
python3 segcanvas.py 5 8        # only sections 5 and 8
python3 ibd1canvas.py           # the same for 18-ibd1-caller.md
```

`ibd1canvas.py` is the same rig built the other way up. To read `IBD1Seg` you must silence
IBD2, and it does that **by paint, not by walls**: a word carrying 34 het-vs-hom mismatches is
perfectly usable to the IBD1 pass and refused outright by the IBD2 one. Every fixture there is
checked to report `IBD2Seg 0.0000` before its `IBD1Seg` is read — an isolation assertion, not
an assumption.

**Two hazards.**

1. **The caches are reference-only.** `segcanvas_measured.json` (6 416 answers) and
   `ibd1canvas_measured.json` (1 013) are what make these rigs re-run in a second. Both rigs
   write whatever they measure back into them and take the binary from `$KING`. Pointing
   `$KING` at our build in the tree silently replaces the reference's answers with our own,
   and nothing will tell you afterwards.
2. **The reference aborts at random on small constructed filesets.** It has a major-allele QC
   check seeded from the clock which raises
   `FATAL ERROR - Too many first alleles as the major allele` on maybe one run in three for
   skewed fixtures. `segcanvas.py` retries up to 24 times with a sleep for exactly this
   reason. **A one-shot probe of a fixture is not evidence** — run it a dozen times and count.
   (`docs/PARITY.md` §5.10's two divergences were each measured 12 times per condition for
   this reason.)

### 8.4 Grading *our* binary on the canvases — `gradebinary.py`

```bash
python3 docs/research/fixtures/gradebinary.py target/release/king          # IBD2Seg, 6000 canvases
python3 docs/research/fixtures/gradebinary.py target/release/king --ibd1   # IBD1Seg, 600
```

This replays exactly the cached canvases with an open-king build and compares marker-interval
readings, so what is graded is **the Rust engine against KING 2.3.2 — no Python model on
either side**. It reads the caches and never writes them; its own answers go to `$TMPDIR`
(`$GRADE_CACHE` to relocate). Exit status is 0 iff every canvas in a closed family matches.

Current: **6 000/6 000** on the IBD2 families, **540/540** on the closed IBD1 families, and
**43/60** on the one family that is deliberately open (`18-ibd1-caller.md` §9, the
`--seglength`-triggered run merge). Note that the binary scores *better* here than
`ibd1canvas.py`'s own model does on the mixed families, because that model pairs `predict1()`
with the pre-§14 IBD2 rule while the binary carries the corrected one — another reason to
grade the thing that ships.

### 8.5 The rule for landing a change

When a fitted rule and a constructed fixture disagree, **the fixture wins**. The corpus is 13
datasets; the fixture is an experiment.

That cuts both ways, and `17-seg-caller.md` §14.10 is the case to remember: a rule change that
moves **no** corpus row can still be a real correction, if a fixture family separates it from
the rule it replaces. The bar for landing such a change is:

* the `.seg` scorecard must not move **at all** — exact rows, `IBD1Seg`, `IBD2Seg`, mean and
  worst `PropIBD`, extra/missing, at 3, 5 and 10 Mb (`seg17.py`, `seg18.py`);
* the canvas count must go **up**, on the binary (`gradebinary.py`), with each clause of the
  change shown independently necessary by ablation;
* `check_mirror.py` must still print `MIRROR OK`, which means `fit/engine.py` was updated to
  match.

Anything that trades one of those against another is a different, worse change. And when the
scorecard *does* move, say by how much in both directions: the `IBD1Seg` overlap rule improved
every headline figure at 3 Mb and made the *worst row* at 5 and 10 Mb slightly worse, and
`docs/PARITY.md` §4.4 says so.
