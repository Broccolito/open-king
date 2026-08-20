# Maintaining open-king

Everything you need to change this project safely: the one rule that cannot be broken, how
the tree is laid out, how the test corpus and its golden output are produced, how to run and
extend the parity suite and its regression baseline, **the fixture technique the whole
segment engine rests on (§8)**, and how to add an analysis.

Read `docs/PARITY.md` first if you want to know what currently works.

If you are picking this up cold and intend to work on the IBD-segment engine, read in this
order: `docs/PARITY.md` §5.0 (what is solved, what is not, what to run next), then §4.6
(the out-of-sample measurement, because it is the only grader that still discriminates),
then §8 below (the instruments), then `docs/research/17-seg-caller.md` and
`18-ibd1-caller.md` for the caller, `19-` and `20-` for what closed it, `21-push-merge.md`
because it corrects `17-…` §6 and three clauses of `20-…`, and **`23-gap-bound.md` last,
because it corrects `21-…`'s standing diagnosis and is the current word on the floor test
and the merge**.

The one-line summary of where things stand: **all 480 cases and all 982 segment rows at all
three captured floors are byte-exact**, so the corpus is now only a regression guard.
Held-out graders still distinguish the exact-64 safety divergence, one segment acceptance
counterexample, rare numeric ties/QC behavior, and unusual reconstruction pedigrees
(`PARITY.md` §4.6, §5.11, §6.2).

---

## 1. The clean-room constraint

**Never read, download, decompile or copy KING's source code.** Not to check a formula,
not to settle an argument, not "just to look".

open-king exists so that a permissively licensed implementation of KING's file formats and
estimators can be embedded in other software. That is only defensible if the implementation
is genuinely independent. The moment a contributor has read the original source, everything
they touch afterwards is derivative, and the fact is not recoverable after the fact — a
reviewer cannot tell by looking at a diff. So the rule is absolute rather than
best-effort, and it is stated at the top of `crates/open-king-core/src/ibdseg.rs` as well as
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
Cargo.toml         the 3-crate workspace; Cargo.lock has 15 packages in total
LICENSE            MIT, covering open-king's own code only
CITATION.cff       machine-readable credit: Manichaikul et al. 2010 for the method,
                   Wei-Min Chen et al. for KING itself, this project third
README.md          the cold-reader entry point; must claim nothing PARITY.md lacks
.gitattributes     -text on every byte-compared tree; see §3, it is load-bearing
.github/workflows/ CI: fmt, clippy, test, build on 3 OSes + the parity baseline gate

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
  research/              the investigation log, numbered in the order it happened.
                         17-seg-caller.md + 18-ibd1-caller.md derive the .seg
                         caller constant by constant; 19-ibd2seg-residual.md adds
                         the segment fringe; 20-seg-writer.md is the two WRITER
                         rules (PropIBD from the printed columns, 16-sample block
                         row order) that finished the default floor;
                         20-seglength-floor.md is the --seglength RUN MERGE, the
                         only caller rule that is dormant at the default floor;
                         21-push-merge.md corrects that merge's IBD2 half and makes
                         the one-word push conditional -- --seglength 5 becomes as
                         byte-exact as the default floor; 22-screen.md is the
                         --related screening count -- it proves the screen is NOT
                         the kinship over any marker subset (exact algebra, covers
                         every subset AND every weighting), then closes marker
                         merging and proves the count is not a function of the
                         markers the budget keeps at all. NOTHING LANDED, on
                         purpose: its §13 lists what is closed and the three leads
                         still worth an experiment. Read it before touching
                         --related's screen;
                         23-gap-bound.md is the LAST segment document -- the gate
                         window's own length bound and the IBD1 merge's budget word
                         set, which make --seglength 10 byte-exact too and refute
                         both diagnoses 21-... left standing. Read it last.
  research/fixtures/     the fixture rigs (see §8 — read it before touching the
                         segment caller): filesets whose answer is forced by
                         construction, used to pin constants the corpus cannot see.
                         fixlab.py builds a fileset and drives the reference;
                         gate8.py brackets the --degree 1 IBD2 clause;
                         segcanvas.py is the .seg canvas (+ its measured cache);
                         push1.py sweeps --seglength as a continuous instrument
                         (21-push-merge.md; keep every sweep inside 1 <= L <= 10);
                         ibd1canvas.py is the same canvas built IBD1-side up;
                         fringecanvas.py builds a segment that does NOT start on a
                         word boundary, which the other two cannot (19-…);
                         mergelab.py bisects the run merge's five conditions at a
                         RAISED floor, where alone it can fire (20-seglength-…);
                         screencanvas.py is the --related two-stage screening rig
                         (single-pair probe + clone canvas), PARITY.md §5.7,
                         screenweight.py is its differential MAF-band probe, which
                         refuted the frequency-standardised lead, and
                         screendeflate.py is the dilution bisection that closed the
                         subset search for good (22-screen.md) -- read it before
                         proposing any subset rule, and note its warning that KING
                         rejects a fileset whose A1 is not the minor allele (§8.3
                         hazard 2);
                         screenfold.py is the FOURTH screening rig and the one that
                         closed the space: it refutes marker MERGING and proves the
                         count is not a function of the markers the budget KEEPS
                         (separation), finds a second necessary condition with no
                         budget involved (gate), and shows the threshold is sharp
                         (sharp). Its ladder fileset -- 48 pairs climbing through
                         the cutoff -- reads the effective threshold in ONE run,
                         17x cheaper than a bisection. `facts` re-measures all of
                         it. Read 22-screen.md §13 before attempting the screen;
                         window1.py is 23-...'s window-bound canvas and carries
                         its held-out draws;
                         chrprobe.py reads the reference ONE CHROMOSOME AT A TIME on
                         the corpus's own data, by MUTING the others for the probe
                         pair (never subset the .bim -- it re-phases the word grid);
                         it is the instrument that localised the 10 Mb residual
                         after two campaigns guessed wrong;
                         oosseg.py is the OUT-OF-SAMPLE differential (PARITY.md
                         §4.6): whole fresh filesets on unused seeds, byte-diffed
                         against the reference -- the grader to use now that the
                         corpus scorecard is saturated;
                         gradebinary.py replays those canvases with OUR binary and
                         grades it against the cached reference answers;
                         segwriter.py proves the two writer rules from the CAPTURES
                         alone — no binary, no fileset, no engine (§8.5);
                         avfs.py / avfs_score.py are the --build AV.FS rig, and
                         buildlog.py scores <prefix>build.log itself (rules |
                         blanks | cut | order | pairs) with build_shapes.py
                         generating held-out shapes. `rules` is the --build
                         scorecard: 53 of 59 held-out shapes byte-identical on the
                         lines we write. Three more --build rigs, all black-box,
                         all scored on shapes the corpus does not contain:
                         clusternum.py pins the STAGED MERGE QUEUE that numbers
                         KING1, KING2, ... and the clustering GATE (19 shapes; its
                         `seeds` mode kills the largest-kinship hypothesis on 4 of
                         8 fresh seeds); dupkeep.py pins which copy of a duplicate
                         is removed and proves the line is rule-half, not
                         inference-half (10 shapes x 3 seeds); siborder.py bounds
                         the ONE open --build question, the sibship member order --
                         it is a hash-table iteration order over a family-scoped,
                         id-keyed container, NOT any sort by f(id) (13 subsets
                         contradict each other 91 times). Reproducing it means
                         identifying the hash; it is worth 3 of 59 shapes and no
                         corpus case

tests/parity/
  generate_corpus.py     builds the 13 input datasets from a fixed seed
  make_alt_inputs.py     builds the alternate --fam / --bim inputs
  capture_params.py      captures the `params` golden group from the reference
  run_parity.py          the acceptance test: replay every capture, diff the bytes
  measure_gaps.py        how *big* each remaining difference is
  verify_formulas.py     checks docs/VERIFIED_FORMULAS.md against captured output
  verify_row_order.py    checks the .kin0 / .ibs0 block-tiled row ordering
  golden/<group>/<case>/ one captured reference invocation per directory
  BASELINE.txt           the recorded outcome of all 480 cases, per case AND per
                         output file; run_parity.py --baseline gates on it (§3)
  work/                  generated inputs (gitignored) — data/ and alt/
  fit/                   analysis scripts used while deriving rules. engine.py is
                         the mirror (below); check_mirror.py asserts it at all
                         three captured floors; scorecard.py is the row-level .seg
                         scorecard measured from the BINARY against the goldens at
                         3 / 5 / 10 Mb, which is what PARITY.md §4.4 quotes
  probes/                small one-off reference probes worth keeping —
                         degree_filter.py (38 298 cases) and xseg_probe.py, which
                         builds 1 040 reference-vs-open-king X.seg runs from fresh
                         seeds using xgen.py's pedigree/map generator
```

**One file in `fit/` is not a scratch script.** `fit/engine.py` is a line-for-line Python
mirror of `crates/open-king-core/src/ibdseg.rs` with every disputed rule exposed as a `Params`
field, so a candidate rule can be scored over the whole corpus in a second instead of a
rebuild-and-replay. It is a *mirror, not a second source of truth*: `fit/check_mirror.py`
asserts that with default `Params` it reproduces the built binary's own `.seg` columns and
`MaxIBD2` on every corpus row. **If you change a rule in `ibdseg.rs`, either update
`engine.py` to match or expect `check_mirror.py` to fail** — and when the two disagree,
the Rust is right and the mirror has the bug. `fit/seg17.py` … `fit/seg21.py` are the
scorecards built on it: each prints the committed rule and the one it replaced side by side
over all 982 primary rows at 3, 5 and 10 Mb — `seg17.py` for the `.seg` IBD2 caller,
`seg18.py` for the `IBD1Seg` overlap rule, `seg19.py` for the IBD2 fringe, `seg20.py` for the
run merge, `seg21.py` for the conditional push and the corrected IBD2 merge, `seg23.py` for
the gate window's length bound and the IBD1 merge's budget word set — and
`R17(...)` … `R23(...)` expose every knob for a candidate. `seg17.py grid19`, `seg18.py grid`,
`seg19.py grid`, `seg20.py grid`, `seg21.py grid` and `seg23.py grid` sweep them all.
`seg21.py` also carries `predict()`, which grades a candidate rule against fresh reference
canvases rather than against the corpus — that is the function to reach for first (§8.7).

**`check_mirror.py` runs at 3, 5 and 10 Mb, and that is load-bearing rather than thorough.**
The run merge cannot fire at the default floor on the corpus's marker spacings. While the
check ran at the default floor only, it passed with the merge committed to `Scan` and absent
from `engine.py` — green light, wrong mirror. The general rule: **a rule whose predicate
reads a CLI parameter must be exercised at a value where it is live**, or the check only
proves the rule is dormant. The same reasoning is why `fit/scorecard.py` exists alongside the
`ibdseg_parity` Rust test, which covers the default floor alone.

`engine.py` also pins four **named parameter bundles**, so every scorecard quoted anywhere in
`docs/research/17-` … `23-` re-runs from that one file: `RETIRED` (the word-aligned geometry,
705 exact rows at 3 Mb), `FRINGE18` (before the IBD2 fringe, 747), `PROP19` (before the writer
rules, 806) and the committed `BASE` (982). All three retired bundles pin `merge=False`,
because each is "the engine as it stood at write-up *N*" and the merge did not land until
`20-`; they additionally pin `merge21=False, push_fraction=None, window_fraction=None,
merge_span="unusable"`, the knobs that step `BASE` back past `21-` and `23-`. Setting the
last two alone gives the tree as `21-` shipped (982 / 970 exact at 5 and 10 Mb against
`BASE`'s 982 / 982); setting all four gives the tree as `20-` shipped (947 / 943). `seg18.py` and `seg19.py` pin the
same thing for the IBD1 pass they borrow from `engine.py`, which is what keeps their
raised-floor rows the "before" the merge is measured against. If you retire a rule, add a
bundle or a knob for it and pin everything that postdates it — that is what keeps the
historical numbers in the research log reproducible instead of merely recorded, and
`check_mirror.py` will not catch you if you don't, because it only checks `BASE`.

**Five committed files are measurement caches, and a non-reference binary will silently
corrupt them.** `docs/research/fixtures/segcanvas_measured.json` (6 416 answers),
`ibd1canvas_measured.json` (1 013), `fringecanvas_measured.json` (576),
`mergelab_measured.json` (shared by `mergelab.py`, `push1.py` and `window1.py`) and
`segfit_measured.json` hold readings *measured from the reference*, which is what lets those
rigs re-run in under a second without it. Each writes whatever it measures back into its
cache and takes the binary from `$KING`, so **never point `$KING` at our build while running
`segcanvas.py`, `ibd1canvas.py`, `fringecanvas.py`, `mergelab.py`, `push1.py`, `window1.py`
or `segfit.py` in the tree.** To grade our own build, use `gradebinary.py` (§8), which reads
those caches, never writes them, and keeps its own answers in `$TMPDIR`, or `oosseg.py`,
which caches nothing at all. `git diff --stat docs/research/fixtures/` before committing;
those five files should only ever change when the reference was the thing being run.

A second hazard in that rig, worth knowing before you trust a single probe: the reference
has a **major-allele QC check**, which aborts a run with
`FATAL ERROR - Too many first alleles as the major allele`. It fires for two unrelated
reasons — a synthetic map that codes `A1` as the *major* allele (deterministic, fatal to
every run, and the cause of many a phantom "no bracket"), and a clock seed (random, roughly
one run in three on small skewed fixtures). `segcanvas.py` retries up to 24 times with a
sleep between attempts for the second. A one-shot probe of a fixture is not evidence — run it
a dozen times and count. **§8.3 hazards 2 and 3 spell both out; read them before you build a
fileset or write a rig.**

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
cargo build --release                                   # -> target/release/open-king
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --check
```

CI (`.github/workflows/ci.yml`) runs exactly these four on Linux, macOS and Windows, and
**then replays all 480 captured reference invocations on Linux** against
`tests/parity/BASELINE.txt`:

```bash
python3 tests/parity/run_parity.py --impl ./target/release/open-king --baseline
```

That is a two-sided gate. `BASELINE.txt` records the outcome of every case *with its
per-file notes*, and the step fails on any difference in either direction:

* a case that now fails, or fails for a *different reason* than recorded — a case that
  failed on `king.seg` numerics before a change and on a missing file after it has the same
  status and the same totals, and a summary count would call that no change;
* a case that now **passes** and is still recorded as failing. This is deliberate. It means
  the committed record can never drift away from what the tree actually does, and it forces
  the person who fixed something to say so in the same commit.

Regenerate with `--baseline --write-baseline` and commit the diff **alongside the change
that earned it**, never on its own.

The suite needs no reference binary — the goldens are committed and the input corpus
regenerates from a seed in about 20 seconds — which is why it can run in CI at all. Only the
480/480 self-check needs the reference, and that stays a local step. The parity job is
restricted to Linux: the goldens were captured on macOS and the suite is run there
constantly, but the full replay has never been checked on Windows (the binary is `king.exe`
there).

**`.gitattributes` is load-bearing, not housekeeping.** 486 of the goldens contain bare `CR`
characters — the reference's own `\r`-separated progress tokens, which are *mid-line* rather
than line terminators. Git for Windows defaults to `core.autocrlf=true`, which would rewrite
those files on checkout and fail cases on a platform where nothing is wrong. `.gitattributes`
marks `tests/parity/golden/**`, `BASELINE.txt` and the `*_measured.json` caches `-text`, so
no end-of-line translation happens anywhere. Do not relax it, and if you add a new tree of
captured output, add it there too.

### Before tagging a release

```bash
cargo fmt --all --check
cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace                                             # 355 passed, 0 failed
cargo clean && cargo build --release            # must work from a clean checkout
python3 tests/parity/run_parity.py --impl target/release/open-king      # record the exact count
python3 tests/parity/run_parity.py --impl target/release/open-king --baseline   # must MATCH
python3 tests/parity/run_parity.py --impl "<reference>"            # must be 480/480
python3 tests/parity/measure_gaps.py --impl target/release/open-king -q # the numbers PARITY.md §4 quotes
KING_GOLDEN=tests/parity/golden \
  cargo test -p open-king-core --test ibdseg_parity -- --nocapture      # .seg scorecard, DEFAULT floor
python3 tests/parity/fit/scorecard.py                              # .seg scorecard, ALL THREE floors
cd tests/parity/fit && python3 check_mirror.py                     # must print MIRROR OK (3/5/10 Mb)
python3 seg17.py && python3 seg18.py && python3 seg19.py && python3 seg20.py  # historical bundles
python3 seg21.py && python3 seg23.py                               # the committed rule, 3/5/10 Mb
python3 docs/research/fixtures/gradebinary.py target/release/open-king          # 6000/6000
python3 docs/research/fixtures/gradebinary.py target/release/open-king --ibd1   # 540/540 + 60/60
python3 docs/research/fixtures/segwriter.py                        # the two .seg writer rules
python3 tests/parity/probes/xseg_probe.py --impl target/release/open-king  # 1040/1040, 625/625
python3 tests/parity/probes/degree_filter.py --ref "<reference>"   # 0 false-keep, 0 false-drop
python3 docs/research/fixtures/oosseg.py --ref "<reference>" \
  --expect-known-safe-divergence                                  # OUT OF SAMPLE: pinned 68/72
python3 tests/parity/probes/segment_residuals.py --ref "<reference>" --impl target/release/open-king
python3 docs/research/fixtures/avscore.py 1 work/*                 # reported intervals: 296/297 exact
# the --build rigs: all four are out-of-sample, none is visible to the 480 captures
cd docs/research/fixtures
python3 buildlog.py rules                                          # 53 match, 6 differ
python3 build_shapes.py                                            # 18 OK, 0 MISMATCH, 2 skipped
python3 clusternum.py score                                        # type 19/19, open-king 19/19
python3 dupkeep.py                                                 # connectivity_then_later 27/0
git ls-files | grep -E '\.(bed|bim|fam|vcf|bcf)$'                  # must print nothing
find . -path ./.git -prune -o -size +95M -print                    # must print nothing
# the caches: 0 keys dropped and 0 VALUES changed -- see §8.3, `git diff --stat` cannot
# tell you this, because inserting into a sorted-key JSON reflows the whole file
```

The current tree measures **480 PASS / 0 FAIL / 480**, self-check **480/480**, all workspace
tests passing (1 timing probe ignored), the row scorecard **982/982/982** at 3 / 5 / 10 Mb with MAE 0.000000 at
each, the out-of-sample differential **68/72** runs and 4 of 6 713 rows (0 extra/missing),
the four `--build`
rigs at **53/59**, **18/18**, **19/19** and **27/27**, a clean build in **9.5 s** from a
pristine copy of the tree, and that clean-tree binary re-measured at 480/480
with `baseline: MATCH` — from a cold tree with no `target/` and no pre-generated corpus,
which is the configuration CI runs in. Do not publish a count you have not just re-run: the
parity number is the project's entire claim, and it is cheap to check (the suite takes about
two seconds warm, eight cold).

**Publish both counts, not just the headline.** The 480 is a *whole-file* number — a case
turns `PASS` only when every row of every file it writes is byte-exact — and the two counts
move independently in **both** directions. `docs/PARITY.md` §4.4 is the row-level scoreboard
and §3 the file-level one; a release note that quotes one without the other misleads one way
or the other. Four cases from this project's history make the point: the `17-seg-caller.md`
§14 correction moved **zero** corpus rows and zero cases and was still right (§8 below is how
that was shown); the `20-seg-writer.md` rules moved **zero estimates** and were worth 28
cases; the `20-seglength-floor.md` run merge moved both, 6 cases and 158 rows;
`21-push-merge.md` moved both again, 3 cases and 62 rows; `23-gap-bound.md` moved both a
third time, 2 cases and 24 rows; and **`<prefix>build.log`'s `RULE` lines moved neither**
while being byte-exact as far as they go, because a case is all-or-nothing. **If the headline
does not move, say so plainly** rather than quoting the row gain as though it were a case
gain — and say so equally plainly when a landing is byte-exact and still costs zero cases.

**And publish the row scorecard at every floor, not just the default.** All three currently
read 982 of 982 with MAE 0.000000. Quoting the case count alone would hide that 99.4 % of
cases rests on 100 % of rows; quoting the row count alone would hide that two cases fail on a
console line the rows cannot see. `tests/parity/fit/scorecard.py` prints all three floors.

**And now that both are saturated, publish the out-of-sample number too.** 982/982 at every
floor is a statement about 480 captures of 13 simulated datasets, not about the caller.
`docs/research/fixtures/oosseg.py` grades the same binary on 24 filesets it has never seen,
at all three floors, byte for byte: **68 of 72 runs**, 4 value-differing rows of 6 713 and
0 extra/missing. The four differences are the deliberate exact-64 safety divergence. A
release note that omits it overstates the claim. `PARITY.md` §4.6.

**Three scales, and they are not interchangeable.** `scorecard.py` compares printed column
against printed column, which is what a user diffing two files sees, so an exact floor reads
MAE 0.000000. `fit/engine.py` compares our *unrounded* `PropIBD` to the reference's printed
one, so the same exact floor reads 0.000017 — half of the exact rows sit half a printed ulp
from the value they round to, a property of the ruler and not an error. And `seg17.py` …
`seg23.py` grade `PropIBD` with the retired **`.kin`** rule, so their `exact` column is much
lower again (`seg23.py` prints 806 / 817 / 820 where `scorecard.py` prints 982 / 982 / 982).
Never mix two of them in one table, and name the grader whenever you quote a number.

The `git ls-files` line is not decoration. `.gitignore` excludes
`/docs/research/fixtures/work/` and `/tests/parity/work/`, but **`.gitignore` does not
untrack files already committed** — a fixture fileset added before the rule stays in the
index forever until someone runs `git rm --cached`. Fixture filesets are regenerated by the
rigs in `docs/research/fixtures/`; none of them is evidence that needs preserving, and none
should ever be committed. This project ships **no genomic data of any kind**: the 13 corpus
datasets are built by a seeded generator at test time, and the committed goldens are the
reference's *output* only.

Keep `docs/PARITY.md` and `README.md` in step with the recorded count: the README must claim
nothing `PARITY.md` does not measure.

Two test suites want the goldens and skip themselves without them:

```bash
KING_GOLDEN=tests/parity/golden cargo test -p open-king-core --test ibdseg_parity -- --nocapture
```

prints the per-dataset `.seg` scorecard rather than just passing. Its two totals are not the
same thing — `row=` counts rows byte-exact on all four printed fields, `est=` only on the two
estimate columns — and quoting `est=` alone overstates the engine. Both are 982 at the
default floor.

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
python3 tests/parity/run_parity.py --impl target/release/open-king \
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
python3 tests/parity/run_parity.py --impl target/release/open-king

# one group or one case, with diffs
python3 tests/parity/run_parity.py --impl target/release/open-king --filter ibdseg/ -v
python3 tests/parity/run_parity.py --impl target/release/open-king \
        --filter core/bigish__related_degree3 -v

# prove the harness itself is sound: must be 480/480
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"

# how big the differences are, not just that they exist
python3 tests/parity/measure_gaps.py --impl target/release/open-king -q
python3 tests/parity/measure_gaps.py --impl target/release/open-king --by-dataset king.seg

# the row-level .seg scorecard at all three captured floors (3 / 5 / 10 Mb), measured
# from the binary against the goldens -- PARITY.md §4.4 quotes this table
python3 tests/parity/fit/scorecard.py
python3 tests/parity/fit/scorecard.py --per-dataset   # split by dataset
python3 tests/parity/fit/scorecard.py --residual      # print every non-exact row

# the regression gate: per case AND per output file, against the committed record
python3 tests/parity/run_parity.py --impl target/release/open-king --baseline
python3 tests/parity/run_parity.py --impl target/release/open-king --baseline --write-baseline

# AND THE ONE THE SUITE CANNOT DO: fresh filesets on unused seeds, byte-diffed against
# the reference. The suite above is saturated on the segment columns (982/982 at every
# floor), so it is a regression guard; this is the grader. PARITY.md §4.6.
python3 docs/research/fixtures/oosseg.py --ref "/path/to/reference/king"
```

Both scripts are Python 3 standard library only. Exit status of `run_parity.py`: 0 all
passed, 1 at least one failed, 2 harness error. Useful flags: `--jobs N`, `--timeout SEC`,
`--max-diff-lines N`, `--keep` (keep the per-case temp directories), `--json OUT`,
`--include-analysis` (also run the 10 `core/_analysis/` captures, which are kept for
analysis rather than as targets).

**The suite is now a regression guard, not a measurement of the caller.** 480 of 480 cases
and 982 of 982 `.seg` rows at every captured floor: a change to the segment engine can only
make those numbers worse, never better. Run it to prove you broke nothing, and grade the
change itself on `oosseg.py` and the canvases (§8.6).

**The self-check is not optional.** `--impl <reference>` must be 480/480. If it is not,
the normalization in `run_parity.py` is either incomplete (a genuinely non-reproducible
line is being diffed) or too aggressive (it is hiding a real difference), and every other
number the harness prints is unreliable until that is fixed.

**Interpreting a failure.** The note after `FAIL` names the file and how it differs:
`!=(bytes)` is a raw byte difference, `!=(num)` means the harness found the difference to
be numeric and will show you a column-wise report with `-v`, `missing:<file>` means we did
not write a file the reference wrote. `measure_gaps.py` turns any of these into rows and
columns and errors.

**`!=(num)` is a hint, not a diagnosis, and `-v` is not optional.** The classifier looks at
the *aligned* rows; if the rows are aligned differently on the two sides it says
`STRUCTURAL` and shows you why. `docs/research/20-seg-writer.md` §4 is the case to remember:
`multifam`'s `king.seg` was reported as a numeric difference for the entire life of the
project, and when the numbers were finally exact the `-v` output showed all 104 rows correct
and in the wrong order. Every row-level tool in this tree matches rows on their identifier
columns *before* comparing — that is what makes `0 extra, 0 missing` meaningful — and it is
also what makes them blind to order. Only the byte diff sees it.

---

## 7. Adding an analysis

Working backwards from the output, which is the order that keeps you honest:

1. **Capture first.** Add golden cases for the new flag across all 13 datasets (and
   whatever parameter combinations matter) *before* writing any Rust. You cannot fit to
   evidence you have not collected, and the captures will contradict your assumptions —
   that is what they are for.
2. **Write the module.** `crates/open-king-cli/src/analysis/<flag>.rs`, declared in
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
that possible, in increasing order of power — §8.3 is the one to read if you read only one,
and §8.5 is the one that found the last two rules, which the canvas structurally could not.

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

**A third rig, `fringecanvas.py`, exists because both of the others are word-aligned.**
`segcanvas.py` and `ibd1canvas.py` lay chr2 out so its first marker starts a 64-marker word,
which means neither can see what a call does in the partial word *beyond* the segment's word
grid — and that fringe turned out to carry its own rule (an opposite homozygote there does
not stop an IBD2 call, though one inside the grid disqualifies its whole word). `fringecanvas`
shortens chr1 by `f` markers so chr2's segment opens `f` markers before a boundary while the
painted words stay aligned, keeping the read-back arithmetic intact. `19-ibd2seg-residual.md`
is the write-up; the clause it found was worth 86 exact `IBD2Seg` rows. The general lesson:
**check what your rig holds fixed by construction**, because that is exactly where the
unmeasured clause will be.

`ibd1canvas.py` is the same rig built the other way up. To read `IBD1Seg` you must silence
IBD2, and it does that **by paint, not by walls**: a word carrying 34 het-vs-hom mismatches is
perfectly usable to the IBD1 pass and refused outright by the IBD2 one. Every fixture there is
checked to report `IBD2Seg 0.0000` before its `IBD1Seg` is read — an isolation assertion, not
an assumption.

**A fourth rig sweeps the CLI instead of the paint — `mergelab.py` and `push1.py`.** The
three above vary the *canvas* and read one number. The `--seglength` rules could not be
reached that way, because the thing under test is a function of a flag: a rule whose
predicate reads a CLI parameter has to be exercised at a value where it is live (see the
lesson at the end of `docs/PARITY.md` §5.0 — `check_mirror.py` once ran only at the default
floor and passed while the mirror was wrong). So `push1.py` holds the canvas fixed — two
words is enough, `CLEAN` = 64 HetHet and `WALL` = 64 opposite homozygotes — and sweeps
`--seglength` as a **continuous instrument**. The jumps of `IBD2Seg(L)` *are* the individual
call lengths, so one sweep reads a canvas's entire behaviour, and a bisection over `L`
locates a threshold to the base pair (the push was pinned at `5.080000` armed / `5.080100`
not, reproduced at five spacings). `python3 push1.py` reproduces all four bisections of
`21-push-merge.md`; `python3 push1.py 2` just the push.

**The read-back arithmetic, restated because a newcomer will not invent it.** On these rigs
the spacing is chosen so **one unit in the last printed place of the column is about 0.11 of
a marker gap**. That is the whole point: the ratio is far below 1, so the printed 4-decimal
value determines the **exact integer** number of marker intervals called. A four-decimal
proportion is not an approximation here — it is a lossless encoding of a small integer, and
every rule in the segment engine was read out of one. If you build a new rig and do not check
that ratio, you are reading noise.

**Four hazards.**

1. **The caches are reference-only.** `segcanvas_measured.json` (6 416 answers),
   `ibd1canvas_measured.json` (1 013), `fringecanvas_measured.json` (576) and
   `mergelab_measured.json` (shared by `mergelab.py` and `push1.py`) are what make these rigs
   re-run in a second. All of them write whatever they measure back into their cache and take
   the binary from `$KING`. Pointing
   `$KING` at our build in the tree silently replaces the reference's answers with our own,
   and nothing will tell you afterwards. **`git diff --stat` cannot check this** — the caches
   are JSON with sorted keys, so inserting new entries reflows the file and the diff reports
   thousands of "deletions" that are nothing of the kind. Check the *values*:

   ```bash
   python3 - <<'EOF'
   import json, subprocess
   p = 'docs/research/fixtures/mergelab_measured.json'      # one per cache
   new = json.load(open(p))
   old = json.loads(subprocess.run(['git','show',f'HEAD:{p}'],
                                   capture_output=True, text=True).stdout)
   print(len(old), '->', len(new),
         '| dropped', sum(k not in new for k in old),
         '| CHANGED', sum(k in new and old[k] != new[k] for k in old))
   EOF
   ```

   `dropped 0` and `CHANGED 0` is the only acceptable result. Any changed value means a
   non-reference binary wrote into the cache and the file must be restored from git.
   (This release: `mergelab_measured.json` went 2 357 → 35 088 entries, 0 dropped,
   0 changed.)
2. **THE A1-MINOR-ALLELE TRAP — code `A1` as the *minor* allele in every synthetic fileset
   you build.** This is the single most expensive footgun in the repository, because it fails
   *silently* rather than loudly. The reference runs a QC check on the ratio of first alleles
   that are the major allele across the map; if too many `.bim` rows name the major allele in
   the `A1` column, it aborts the whole run with

   ```
   FATAL ERROR - Too many first alleles as the major allele
   ```

   A rig that treats a non-zero exit as "the reference said no" then reads that abort as a
   *data point*. In a bisection every probe on one side of the bracket dies the same way, the
   rig concludes there is **no bracket**, and you spend a day debugging a search that never
   ran. The tell is a bisection that reports no bracket on a fixture where you can see by
   construction that the answer must lie inside it, or a sweep whose "rejected" side is
   suspiciously uniform. **Any rig that shells out to the reference must treat this string as
   an error, not as an answer** — the rigs in `docs/research/fixtures/` all do, and a new one
   must too.

   The fix at construction time is one line: when you emit a `.bim`, put the allele you gave
   the *lower* frequency in the `A1` (5th) column, and the higher-frequency one in `A2`.
   `generate_corpus.py` and every rig here already do; copy one of them rather than writing a
   `.bim` writer from scratch. Note this is a property of the *whole map*, not of any single
   marker — a handful of major-`A1` rows is fine, a systematically flipped map is fatal.

3. **The same error also fires at random, and that is a different problem.** Independently of
   how you coded `A1`, the check is seeded from the clock, so on small or skewed constructed
   filesets it raises the same string on maybe one run in three even when the map is coded
   correctly. `segcanvas.py` retries up to 24 times with a sleep for exactly this reason.
   **A one-shot probe of a fixture is not evidence** — run it a dozen times and count.
   (`docs/PARITY.md` §5.10's two divergences were each measured 12 times per condition for
   this reason.) Distinguish the two cases by retrying: a miscoded map fails *every* attempt,
   the clock-seeded abort fails *some*.

4. **`--seglength` only behaves like a floor inside `1 ≤ L ≤ 10` Mb** (`21-push-merge.md`
   §8.2). Above ~10 Mb a 14.06 Mb call reports as a constant 8.93 Mb at every larger floor, on
   canvases of quite different content; below 1.0 the flag behaves as though absent. Any sweep
   that wanders outside that window is measuring a different regime. Relatedly, the `.seg`
   floor test is strictly **`>`** — a call of exactly 5 100 000 bp is reported at
   `--seglength 5.099999` and dropped at `5.100000` — and the shipped engine compares `>=`.
   Nothing on the corpus lands on an exact tie, so it never bites there; every fixture in the
   merge rig does, so a new canvas battery that ignores it will misgrade its own boundary
   rows.

### 8.4 Grading *our* binary on the canvases — `gradebinary.py`

```bash
python3 docs/research/fixtures/gradebinary.py target/release/open-king          # IBD2Seg, 6000 canvases
python3 docs/research/fixtures/gradebinary.py target/release/open-king --ibd1   # IBD1Seg, 600
```

This replays exactly the cached canvases with an open-king build and compares marker-interval
readings, so what is graded is **the Rust engine against KING 2.3.2 — no Python model on
either side**. It reads the caches and never writes them; its own answers go to `$TMPDIR`
(`$GRADE_CACHE` to relocate). Exit status is 0 iff every canvas in a closed family matches.

Current: **6 000/6 000** on the IBD2 families, **540/540** on the closed IBD1 families, and
**60/60** on the family that used to be deliberately open. That last one is the clearest
single confirmation of the run merge: `18-ibd1-caller.md` §9 left `mixed seed 61803, L=8` at
**43/60** for two write-ups because the `--seglength`-triggered merge was measured but not
modelled; with `20-seglength-floor.md`'s five conditions committed it closes completely,
without any canvas in it having been used to derive them. Note that the binary scores
*better* here than
`ibd1canvas.py`'s own model does on the mixed families, because that model pairs `predict1()`
with the pre-§14 IBD2 rule while the binary carries the corrected one — another reason to
grade the thing that ships.

### 8.5 When the canvas cannot help: reading the reference's own output

The canvases answer "what did the caller do?". They cannot answer "what did the *writer*
do?", because each one reads a single printed column back as a marker count and never
compares two files. Two of this project's rules live entirely in that blind spot, and the
technique that found them is worth stating separately —
`docs/research/fixtures/segwriter.py` is the whole rig, and it drives nothing at all:

**The captures are themselves an experiment.** `tests/parity/golden/` holds 4 172 `.seg`
rows, 4 248 `.kin` rows and so on, written by the reference. Any hypothesis that is a
*function of printed values* can be tested against all of them without a binary, a fileset
or an engine — and tested to a standard the corpus scorecards cannot reach, because there is
no caller error in the way.

That is how `PropIBD` was settled. The hypothesis "`.seg` computes it from its own printed
columns" is checkable row by row in exact integers; it survived 4 172 rows with 0
refutations, was refuted on `.kin` (42 rows) and `cluster.kin` (3), and the 1 313 rows where
it lands on an exact decimal half then discriminated nine candidate expressions that are all
*mathematically identical*, leaving exactly one. None of that needed the reference to be run
even once.

**Three habits this rewards:**

1. **Diff the reference against itself.** Two files from one invocation that carry the same
   pair are two measurements of the same quantity. When they disagree — `.kin` and `.seg`
   print different `PropIBD` for 54 of the 201 pairs they share — that is not noise to be
   averaged away, it is the statement that there are two rules and you have been looking for
   one.
2. **Count refutations, not agreements, and know what the test could have found.** "0
   inconsistencies in 4 172 rows" is worth nothing until you can say how often the competing
   hypothesis *would* have been caught. Here: 1.7 % of rows put the full-precision value more
   than half a printed ulp from the printed-column combination, so the rival predicts ≈ 71
   refutations. Zero were seen. Do that arithmetic before believing the result.
3. **When a hypothesis has a free tie-break, keep looking.** "Round half up" fitted 84 % of
   the ties and would have improved every scorecard — and it was wrong, with 214
   counterexamples in the reference's own output. The right expression has no free parameter
   and fits all 1 313. A rule that needs a fitted constant to cover its exceptions is usually
   the wrong rule, not a rule with exceptions.

**And the ordering lesson.** The 16-sample block order was visible in principle from the
first day and invisible in practice for the whole project, because every grader in the tree
matches rows on their identifier columns before comparing anything. `measure_gaps.py`
reported `0 extra, 0 missing` on those files throughout. If a file's numbers are exact and
the case still fails, the answer is in the bytes: run `run_parity.py … -v` and read the
`STRUCTURAL` report.

**A fifth rig, for when a corpus row is wrong and you do not know where —
`chrprobe.py`.** The canvases are constructed data; the corpus is real data whose answer
nobody knows. Between them sits the question "which segment of which pair produced this
wrong number?", and for two campaigns it was answered by guessing. `chrprobe.py` answers it
by measurement: pick the pair, **mute** every chromosome but one — set one sample to `A1A1`
and the other to `A2A2` across it, so no caller can find anything there — and the reference
then prints that chromosome's own called length in the same four-decimal column. Twenty-two
runs localise the fault. Two traps it documents, both of which gave a wrong answer first:

* **Mute, never subset the `.bim`.** Deleting rows re-phases the 64-marker word grid and
  changes calls on every later chromosome, so a subsetted probe measures a different
  program.
* **`--seglength` is in Mb and clamped to `1 ≤ L ≤ 10`**, silently falling back to 3 Mb
  outside. A sweep that wanders outside that window is measuring the default floor and
  will look like a rule that "stops behaving like a floor".

**A sixth, for the release claim itself — `oosseg.py`.** Everything above grades a *clause*.
`oosseg.py` grades the *program*: 24 whole filesets built by `generate_corpus.py`'s own
simulator on seeds used nowhere else in the tree, run through both binaries at three floors,
and diffed byte for byte. It is the only instrument here that answers the question a user
actually asks — "will this give me KING's answer on my data?" — and its answer is currently
68 of 72 runs, not 72. It reports extra, missing and value-differing rows separately; the
current result has 0 extra, 0 missing and four value rows caused by KING's exact-64 unsafe
tail read.

### 8.6 The rule for landing a change

When a fitted rule and a constructed fixture disagree, **the fixture wins**. The corpus is 13
datasets; the fixture is an experiment.

**Read this section knowing the corpus scorecard is now saturated.** 982 of 982 rows at 3, 5
and 10 Mb: it can only move **down**. It is a regression guard, not a grader. Every "must go
up" below now means one of the three graders that can still go up — `oosseg.py` (68/72 whole
filesets, `PARITY.md` §4.6), `gradebinary.py` (6 000/6 000 and 600/600 canvases), and the
held-out draws in `window1.py`/`mergelab.py`/`push1.py`.

That cuts both ways, and `17-seg-caller.md` §14.10 is the case to remember: a rule change that
moves **no** corpus row can still be a real correction, if a fixture family separates it from
the rule it replaces. The bar for landing such a change is:

* the `.seg` scorecard must not move **at all** — exact rows, `IBD1Seg`, `IBD2Seg`, mean and
  worst `PropIBD`, extra/missing, at 3, 5 and 10 Mb (`tests/parity/fit/scorecard.py` for the
  binary's own numbers, `engine.py` for a candidate rule, or the `seg17.py` … `seg23.py`
  scorecards for the historical comparisons). Since all three floors read 982, this now means
  literally "no row regresses";
* the out-of-sample count must not go down, and should go up: `oosseg.py --ref <reference>`,
  currently 68 of 72 runs. This is the grader that separated the window bound from the rule it
  replaced (66/72 against 60/72 at that landing, all six gains at `--seglength 10`, none lost)
  and later pinned the merged-call pair filter (68/72, 0 missing);
* the canvas count must go **up**, on the binary (`gradebinary.py`), with each clause of the
  change shown independently necessary by ablation (`seg21.py grid` and `seg23.py grid` are
  the worked examples: in `23-…`'s grid the window bound is worth 970/972 → 982/982 at 10 Mb,
  the budget word set two further `IBD1Seg` rows, the IBD1 side of the window bound **zero on
  the corpus** and 32 of 360 held-out canvases — which is why it is in — and the `pre_merge`
  variant is worse and is out);
* `check_mirror.py` must still print `MIRROR OK`, which means `fit/engine.py` was updated to
  match, **and** the retired parameter bundles it pins (`RETIRED`, `FRINGE18`, `PROP19`, plus
  the `merge21=False, push_fraction=None, window_fraction=None, merge_span="unusable"`
  step-backs to the trees `20-…` and `21-…` shipped) must still reproduce the numbers quoted
  for them in `docs/research/17-` through `23-`.

For a change that *does* move the scorecard, the bar is: **exact rows up and mean error not
worse**, at every floor. Anything that trades one against the other is a different, worse
change — report the trade numerically and leave it out. Two examples, both refused: cutting
`IBD1Seg` by all `IBD2` calls rather than the surviving ones (982 → 950 at 3 Mb), and the
**unconditioned** run merge — joining any two calls separated by less than `--seglength` —
which improves the worst row at 10 Mb while losing exact rows, and which invents 251 pairs at
5 Mb once the merged calls feed the >10 Mb pair filter.

That last one is worth dwelling on, because the conditioned version of the same idea *did*
land (`20-seglength-floor.md`, +6 cases, `IBD1Seg` 844 → 960 at 10 Mb). The difference was
not persistence: it was that the second attempt bisected each of the merge's five conditions
against the reference on constructed canvases, and then validated the whole rule on 360
held-out canvases with unused seeds, instead of turning knobs until the corpus improved. A
rule that improves the corpus is not thereby correct — see §8.7.

**And a landed rule is not thereby finished — twice over.** `23-gap-bound.md` re-opened the
merge a second time and found that **both** diagnoses `21-…` §8.1 published were wrong: it
was neither a second bound on the merge's gap nor an invented merge. The residual was the
floor being asked a *second* time, of the gate window rather than of the reported call, plus
the IBD1 merge's budget reading a wider set of words than its own cap. Two campaigns had
theorised about which segment of which pair was at fault; `fixtures/chrprobe.py` answered it
directly by muting every other chromosome for the probe pair and reading the reference's own
per-chromosome call. **Localise a wrong row before theorising about it** — and mute, never
subset the `.bim`, which re-phases the 64-marker word grid and changes later calls.

`21-push-merge.md` re-opened three clauses of
that same merge — the word cap, what the interruption is measured between, and what `X`
counts — and found all three wrong on the IBD2 pass, plus the one-word push of `17-…` §6
wrong as a fourth. Each was found by building a fixture that *isolates* the clause, not by
re-scoring the corpus: the cap, for instance, is invisible on the corpus at any value, because
at ~50 kb spacing a 10 Mb gap holds at most three words. The lesson to carry: when a rule
lands with a residual, suspect the clauses the fixture that derived it could not vary. `20-…`
§3 measured the two-word cap on the IBD1 pass and the IBD2 pass simply inherited it; that
inheritance was the bug, and it took a fixture built for the IBD2 pass alone to see it.

And when the scorecard does move, say by how much in **both** directions: the `IBD1Seg`
overlap rule improved every headline figure at 3 Mb and made the *worst row* at 5 and 10 Mb
slightly worse, and `docs/PARITY.md` §4.4 says so.

**Do not re-sweep the caller's constants.** Forty single-knob perturbations of
`fit/engine.py`'s `Params` and all 32 combinations of the two IBD1 endpoint rules crossed
with the two IBD1 fringe rules were scored in the final pass: the committed values are the
unique maximum on both exact rows and mean error. Five knobs — `bridge_rule="17"`,
`gate_end="right"`, `inf2_ibs1b=True`, `ibd1_clip_ibd2=True`, `clip_before_len=False` — score
*identically* to the committed engine on every corpus row, so if you change one of those, the
canvases are the only thing that can tell you whether you were right.

---

### 8.7 Never fit to the corpus

**This is the most important rule in this file, and it is the one that has actually been
broken.** A previous session derived several segment-caller constants by searching for the
values that maximised the corpus scorecard. The numbers looked excellent. They were wrong,
they had to be thrown away, and the whole rule had to be re-derived from scratch on
constructed fixtures — weeks of work redone. The tell was that the constants had no
provenance: nobody could say what experiment fixed them, only what score they produced.

**Why the corpus cannot settle a rule.** It is 13 datasets and 982 graded rows. A rule with
three or four free parameters has enough freedom to absorb a good deal of that, so "the
scorecard improved" is weak evidence about the reference's behaviour and strong evidence only
about your search. The failure is silent and it is *anti*-correlated with the signal you want:
the better the fit, the more confident the wrong rule looks. And because the corpus is also
the acceptance test, a fitted rule cannot be caught by running the tests.

**What to do instead.** Establish the rule where you can vary one thing at a time and read
the answer exactly — a constructed fixture, driven against the reference (§8.1–§8.3):

1. **Bisect each constant separately** on a canvas built to isolate it. "9 merges and 10 does
   not, against A1A1/A1A1 loads of 16, 24, 30 and 40" is a measurement. "4 scored best" is
   not.
2. **Validate out of sample, on seeds the derivation never saw**, and say how many. The run
   merge was landed on 360/360 held-out canvases at 5 and 10 Mb with three unused seeds, plus
   600/600 independently drawn interruptions. If a rule cannot survive fresh canvases, it is a
   description of the corpus.
3. **Only then look at the corpus**, as a *check*, under §8.6's bar. It is the last step, not
   the search space.

**A sharp negative is a good outcome; a fitted fiction is not.** If the honest answer is "the
rule is not identified, here is the boundary I bracketed and the models the data refutes",
write that down and leave the code alone. `docs/PARITY.md` §5.7 is exactly this, and it has
now been through all three stages of the discipline.

1. **A promising lead, recorded as a lead.** A frequency-standardised estimate over a
   MAF-selected subset reproduced `bigish`'s own count — and only when it was given the
   corpus's own frequencies. It was written down and *not* committed.
2. **A differential rig refuted it** (`fixtures/screenweight.py`), and refuted the whole
   family it belonged to.
3. **`22-screen.md` closed the search space by proof and by measurement, and declined to fit
   the one free constant that would have passed.** The algebra is exact: at a marker of
   frequency `p`, `E[N_l] = 4pq(1 − 2φ)` and `E[het_l] = 2pq`, so numerator and denominator
   are both proportional to `Σ pq` over *whatever* index set they are summed on — **every**
   subset and every non-negative weighting is unbiased for the same φ. Measurement agrees:
   nine permutations of the marker order print the same count; replicating the map leaves
   every kinship bit-identical and still moves the count. What the screen *does* obey was
   then pinned to 0.2 % — `k_screen = 0.5 + R(k − 0.5)`, with `R ≡ 1` whenever the map holds
   at most 32 768 markers. And `R` was **not** landed, because it swings from 0.998 to 1.085
   with the MAF spectrum: fitting it would have reproduced `bigish` and nothing else. Two
   `--related` cases still fail as a result.

That is the correct trade at every stage — a closed-off search space is worth more than a
rule that scores well and is wrong, and a measured law with an unfitted constant is worth
more than a fitted constant with no law.

**How to tell you are doing it.** If you cannot name the experiment that pins a constant, you
fitted it. If your justification is a scorecard delta rather than a bisection, you fitted it.
If the rule has a knob you set by trying values, you fitted it.
