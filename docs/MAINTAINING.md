# Maintaining open-king

Everything you need to change this project safely: the one rule that cannot be broken, how
the tree is laid out, how the test corpus and its golden output are produced, and how to
add an analysis.

Read `docs/PARITY.md` first if you want to know what currently works.

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
  research/              the investigation log, numbered in the order it happened
  research/fixtures/     the fixture rigs: filesets whose answer is forced by
                         construction, used to pin constants the corpus cannot see

tests/parity/
  generate_corpus.py     builds the 13 input datasets from a fixed seed
  make_alt_inputs.py     builds the alternate --fam / --bim inputs
  capture_params.py      captures the `params` golden group from the reference
  run_parity.py          the acceptance test: replay every capture, diff the bytes
  measure_gaps.py        how *big* each remaining difference is
  golden/<group>/<case>/ one captured reference invocation per directory
  work/                  generated inputs (gitignored) — data/ and alt/
  fit/                   analysis scripts used while fitting rules
  probes/                small one-off reference probes worth keeping
```

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

CI (`.github/workflows/`) runs exactly these four on Linux, macOS and Windows. All four
must be clean before the parity suite is worth looking at.

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

### If you cannot work a rule out

Two techniques from `docs/research/` that repeatedly worked where staring at the corpus
did not:

* **Construct the answer.** `docs/research/fixtures/fixlab.py` builds a fileset in which
  one pair's IBD state is exact by construction — explicit shared haplotypes, and if you
  want, every genotype of every sample written by hand. A threshold measured in "number of
  markers of kind X" then becomes a threshold in the statistic itself, with no sampling
  noise to fit through. This is how `MIN_INFORMATIVE = 10` was pinned to the unit, and how
  the `--degree 1` IBD2 clause was found — the corpus cannot see either.
* **Move one thing.** Deleting the first *m* markers of a fileset shifts the global
  64-marker word grid by *m* and changes nothing else, so genotypes and segment lengths
  stay fixed while the reference's verdict moves. 512 such invocations validated the
  informativeness gate against data that had no part in choosing it.

And when a fitted rule and a constructed fixture disagree, **the fixture wins**. The corpus
is 13 datasets; the fixture is an experiment.
