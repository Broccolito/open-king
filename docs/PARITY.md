# Differential Parity Harness

`tests/parity/run_parity.py` is the acceptance test for open-king: it replays every
captured KING 2.3.2 invocation with **our** binary and diffs the result against the
reference's recorded output. Parity is reached when every case passes.

```
usage: run_parity.py --impl PATH_TO_OUR_KING [--ref PATH_TO_REFERENCE_KING]
                     [--data DIR] [--alt DIR] [--golden DIR] [--filter SUBSTRING]
                     [--update] [--include-analysis] [--jobs N] [--timeout SEC]
                     [--max-diff-lines N] [--keep] [--json OUT] [-q] [-v]
```

Python 3 standard library only — no build step, no dependencies.

---

## 1. Quick start

```bash
# full suite against our binary
python3 tests/parity/run_parity.py --impl target/release/king

# one flag group, with diffs for whatever fails
python3 tests/parity/run_parity.py --impl target/release/king --filter core/ -v

# one case
python3 tests/parity/run_parity.py --impl target/release/king \
        --filter core/bigish__related_degree3 -v

# prove the harness itself is sound (must be 100% PASS)
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"
```

Exit status: **0** = all cases passed, **1** = at least one case failed, **2** = harness
error (missing binary, no cases matched, `--update` without `--ref`). Suitable for CI as-is.

The datasets are **not committed** — only the seeded generator is. The harness regenerates
any missing dataset automatically before the first case runs:

```
[parity] regenerating 13 dataset(s) into tests/parity/work/data: trio nuclear ...
```

That takes ~20 s once. `tests/parity/work/alt/` (the alternate `--fam` / `--bim` inputs the
`params` group needs) is rebuilt the same way from `make_alt_inputs.py`.

---

## 2. What a case is

Each directory under `tests/parity/golden/<group>/<case>/` is one recorded reference run:

| file | meaning |
| --- | --- |
| `cmd.txt` | the argv, with placeholders `{KING}`, `{DATA}`, `{ALT}` |
| `exitcode.txt` | the reference's exit status |
| `stdout.txt` | reference stdout, verbatim (contains `\r` progress tokens) |
| `stderr.txt` | reference stderr, verbatim — **empty in all 480 cases**, fatal errors included: KING writes even `FATAL ERROR` to stdout and exits 1. Anything we write to stderr is a parity failure. |
| *everything else* | the files the reference wrote into its working directory |

`cmd.txt` exists in two captured shapes and the harness accepts both: a single shell line
(`core`, `ibdseg`, `params`) or one token per line (`apps`).

| group | cases | flags exercised |
| --- | ---: | --- |
| `core` | 104 | `--kinship`, `--related [--degree N]`, `--duplicate`, `--ibs` |
| `apps` | 91 | `--unrelated`, `--build`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster` |
| `ibdseg` | 65 | `--ibdseg [--degree N] [--seglength N]`, `--related --ibdseg` |
| `params` | 220 | `--prefix`, `--cpus`, `--minConc`, `--sexchr`, `--fam`, `--bim`, `--degree`, error probes |
| **total** | **480** | 876 output files byte-compared |

`golden/core/_analysis/` holds ten supplementary runs whose output files were deliberately
pruned at capture time (only their stdout is evidence). They are **skipped by default**;
`--include-analysis` replays them and will report six `extra:` failures by design.

## 3. What is compared

For every case the harness runs our binary in a **fresh temp directory** and checks:

1. **exit status** against `exitcode.txt`;
2. **stdout** and **stderr**, after normalization (§4), byte for byte;
3. **the set of files produced** — both files we failed to write and files the reference
   never wrote are failures (`missing:king.kin0`, `extra:king.unexpected`);
4. **the bytes of every shared file**.

Nothing is compared loosely. Aside from the normalization rules below and the one
diff-excluded file class in §6, a passing case means byte-identical output.

## 4. Normalization rules

KING's stdout carries information that cannot be reproduced by any implementation on any
other machine or run. These rules are applied to **both sides** — the golden capture and our
output — before diffing. They are implemented as `LINE_RULES` / `BLOCK_RULES` in
`run_parity.py`; add a rule there and document it here.

| # | Target | Rule | Why |
| --- | --- | --- | --- |
| **R1** | wall-clock timestamps | `Thu Aug 13 18:01:49 2026` → `<TS>`, anywhere in the line | `KING starts/ends at`, `… Inference ends at`, etc. Matched **anywhere**, not line-anchored: KING emits progress with no newline, so the terminator arrives as `0%1%2%<41 spaces>ends at <ts>` on one physical line. |
| **R2** | progress percentages | delete lines that are entirely `N%` tokens; strip a leading `(N%\r?)+` run from every other line | Written as `"%d%%\r"` with no newline, so they glue onto the head of the next message. **How many appear is a function of thread count** (8 threads → `0%…10%`, 1 thread → `0%…90%`). Only *leading* runs are stripped, so a real value like `rate > 100%` survives. |
| **R3** | `N CPU cores are used` | leading integer → `<NCPU>` | Host-dependent when `--cpus` is absent (the 16-core capture machine printed `8`). Indentation and the `.` vs `...` spelling are preserved — those are real signal. |
| **R4** | absolute input paths | any absolute path ending in `.bed`/`.bim`/`.fam` → `{DIR}/<basename>` | Echoed 4–6× per run (`Binary File :`, `Read in PLINK … file …`, `--fam [<p>]`). Depends on where the repo is checked out and on `--data`/`--alt`. |
| **R5** | `--noscreen [-1717986816]` | bracketed integer → `<NOSCREEN>` | KING prints an uninitialized default here. It is constant on the capture host, but it is nobody's documented value and a reimplementation should not be held to it. |
| **R6** | the `parameters are in effect:` banner | fold continuation lines back into their entry, collapse interior whitespace | The block is **word-wrapped at a fixed column**, so the *length of the input path* decides where lines break and which entries share a line. Compared structurally instead of by column. |

Deliberately **not** normalized, because they are real signal:

* the version banner `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`;
* every numeric column of every output file;
* the whole `Relationship summary` table;
* `Total length of N chromosomal segments …` and `Autosome genotypes stored in N words …`;
* the `Options in effect:` block, including KING lower-casing `--bySNP` to `--bysnp`;
* KING's literal typos `3nd-degree` / `4nd-degree` — **reproduce them**.

R4 and R6 were validated by re-running the whole `params` group with `--data` and `--alt`
pointed at directories of very different path lengths from the capture: 220/220 still passed,
so path length does not leak into the comparison.

## 5. Reading a failure

```
FAIL  core/bigish__kinship  king.kin!=(num)
      $ king -b {DATA}/bigish.bed --kinship
  --- king.kin (35035 B golden vs 35033 B ours) ---
    ROUNDING - every delta within 1 ulp of the printed precision [max|d| 1.000e-04, max rel 1.000e+00]
      col  9 Kinship    573 row(s) differ | max|d| 1.000e-04 @row 7 (0.2512 -> 0.2513) | = 1.0 ulp | max rel 1.000e+00
    --- golden/king.kin
    +++ ours/king.kin
    @@ line 2 @@
    -BF01	B01_C1	B01_C2	50000	0.250	0.2500	0.2061	0.0158	0.2505	0
    +BF01	B01_C1	B01_C2	50000	0.250	0.2500	0.2061	0.0158	0.2506	0
```

For any tabular file the harness reports, **per column**, how many rows differ and the worst
absolute and relative delta with the row that produced it. The one-line verdict tells you
which kind of bug you have:

| verdict | meaning |
| --- | --- |
| `ROUNDING` | every delta ≤ 1 ulp of KING's printed precision — a formatting/tie-break difference, not a wrong formula |
| `NEAR-ROUNDING` | ≤ 2 ulp — suspect accumulation order or a rounding mode |
| `ALGORITHMIC` | larger — a wrong formula, wrong denominator, or wrong SNP set |
| `STRUCTURAL` | row counts, field counts, header, or text columns differ — usually **row order** or a row-selection predicate, not arithmetic |

The verdict is decided on **ULPs of the printed precision**, not on the relative delta:
`0.0001` vs `0.0002` is a relative error of 1.0 but only one unit in the last printed place.
Relative delta is still reported, because it is the right measure away from zero.

`-v` adds a unified diff, capped at `--max-diff-lines` (default 20) per file. When both
sides have the same line count the lines are paired positionally (`@@ line N @@`), which
reads far better than one giant hunk on a 19 000-row `.kin0`.

Useful when a case is hard to reproduce by hand: `--keep` leaves each case's temp directory
in place and prints the path; `--json out.json` writes machine-readable results.

## 6. Diff-excluded files (the one exception)

`kingX.kin0` — the **X-chromosome between-family** file — is corrupted by a genuine data race
in KING 2.3.2 whenever it runs with more than one thread: several threads append to one
unlocked `FILE*`, so records tear mid-field. Over 20 runs, `--cpus 1` gave 20/20 identical
output while `--cpus 4` gave 19 distinct files, none correct, with the size wandering between
138 and 664 bytes.

The harness therefore skips the byte-diff of any file whose name ends in **`X.kin0`** in a
case whose command line does not pin `--cpus 1` (match on the *suffix*: `--prefix` renames the
file). Eight files across the corpus are excluded this way; the summary line reports the
count. The authoritative X goldens are the `params/sexchr__kinship_cpus1*` cases, which **are**
compared. Our implementation should match the `--cpus 1` output — do not reproduce the race.

Files listed under `diff_exclude` in `golden/params/runs.json` are honored as well.

## 7. Re-capturing goldens (`--update`)

```bash
python3 tests/parity/run_parity.py --impl target/release/king \
        --ref "/path/to/reference/king" --filter core/multifam --update
```

`--update` re-runs the **reference** binary and overwrites `stdout.txt`, `stderr.txt`,
`exitcode.txt` and every output file in the selected case directories. It refuses to run
without `--ref`, so it can never be triggered by a routine test invocation. Honor `--filter`
to re-capture a subset, and review the result with `git diff` before committing — a golden
that changes without a stated reason means the corpus, not the implementation, moved.

`--update` does not maintain `golden/params/runs.json`; that file is capture-time metadata.

## 8. Harness self-check

Pointing `--impl` at the reference binary itself must be a clean sweep. That is the proof
that the normalization rules are complete and that the harness reports no false positives:

```
$ python3 tests/parity/run_parity.py --impl "<reference king 2.3.2>" -q
parity: 480 PASS, 0 FAIL, 480 total (2.0s wall, 876 output file(s) byte-compared, 8 diff-excluded)
```

Reproduced across three parallel runs, a serial (`--jobs 1`) run, a run against a
freshly-regenerated corpus, and a run with relocated `--data`/`--alt` directories.

The harness was also mutation-tested — each of these deliberately broken wrappers around the
reference binary is caught:

| mutation | detected as |
| --- | --- |
| `Kinship` column +0.0001 (1 ulp) | `king.kin!=(num)` → `ROUNDING` |
| `Kinship` column ×2 | `king.kin0!=(num)` → `ALGORITHMIC`, 2959 ulp |
| `king.kin0` rows re-sorted (breaks the 32×32 tiling) | `king.kin0!=(num)` → `STRUCTURAL` |
| one output file deleted | `missing:king.kin0` |
| one extra output file written | `extra:king.unexpected` |
| wrong exit status | `exit 3 != 0` |
| one real stdout line dropped | `stdout!=` |

Re-run the self-check whenever a normalization rule changes.

---

## 9. Parity matrix

Regenerated from `run_parity.py --json`, not edited by hand. `PASS` means **byte-identical**
output — every file, every column, plus stdout, stderr and exit status under the §4 rules —
and nothing weaker. `FAIL` is any difference at all, however small; §11 says how large each
one actually is.

**Overall: 365 / 480 PASS** in the default suite, and **367 / 490** with `--include-analysis`.

```
$ python3 tests/parity/run_parity.py --impl target/release/king -q
parity: 365 PASS, 115 FAIL, 480 total (1.3s wall, 670 output file(s) byte-compared, 6 diff-excluded)
```

| group | PASS | cases | what is missing |
| --- | ---: | ---: | --- |
| `core` | 57 | 104 | `--related` above the downgrade, `MaxIBD2`/`Pr_IBD2` in `--ibs` |
| `apps` | 73 | 91 | `--autoQC` entirely; the merged-cluster tail of `--build`/`--cluster`/`--unrelated` |
| `ibdseg` | 15 | 65 | the `.seg` estimates and `splitped.txt` |
| `params` | **220** | **220** | — |
| **total** | **365** | **480** | |

The ten supplementary `core/_analysis/` runs are skipped by default because their output files
were pruned at capture time; **the reference binary itself scores only 4/10 on them**, so
484/490 is the ceiling for anyone. We score 2/10: the two extra failures are `--related` cases,
the same gap §11.1 describes.

### 9.1 `core` — `--kinship`, `--related`, `--duplicate`, `--ibs`

| dataset | `duplicate` | `ibs` | `kinship` | `related` | `related_degree1` | `related_degree2` | `related_degree3` | `related_degree4` |
|---|---|---|---|---|---|---|---|---|
| `trio` | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |
| `nuclear` | PASS | FAIL | PASS | PASS | PASS | PASS | PASS | PASS |
| `threegen` | PASS | FAIL | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `multifam` | PASS | FAIL | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `dups` | PASS | FAIL | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `missing` | PASS | FAIL | PASS | PASS | PASS | PASS | PASS | PASS |
| `monomorphic` | PASS | PASS | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `sexchr` | PASS | PASS | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `unrelated` | PASS | PASS | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `admixed` | PASS | FAIL | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |
| `singleton` | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |
| `pair` | PASS | PASS | PASS | PASS | PASS | PASS | PASS | PASS |
| `bigish` | PASS | FAIL | PASS | FAIL | FAIL | FAIL | FAIL | FAIL |

The five datasets whose `--related` passes are exactly the five under ten samples, where the
reference replaces the pass with `--kinship`. `--ibs` passes on the six datasets whose pairs
never reach the `MaxIBD2` gate.

### 9.2 `apps` — `--unrelated`, `--build`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster`

| dataset | `unrelated` | `unrelated_degree2` | `build` | `bysample` | `bySNP` | `autoQC` | `cluster` |
|---|---|---|---|---|---|---|---|
| `trio` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `nuclear` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `threegen` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `multifam` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `dups` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `missing` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `monomorphic` | PASS | PASS | FAIL | PASS | PASS | FAIL | PASS |
| `sexchr` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `unrelated` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `admixed` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `singleton` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `pair` | PASS | PASS | PASS | PASS | PASS | FAIL | PASS |
| `bigish` | FAIL | FAIL | FAIL | PASS | PASS | FAIL | FAIL |

`bigish` is the only fileset with the hundred samples the reference requires before it will
join two families, so it is the only one that reaches the merged-cluster code path in all
three clustering analyses.

### 9.3 `ibdseg` — `--ibdseg` and its parameters

| dataset | `ibdseg` | `ibdseg_degree2` | `ibdseg_seglength5` | `ibdseg_seglength10` | `related_degree2_ibdseg` |
|---|---|---|---|---|---|
| `trio` | PASS | PASS | PASS | PASS | PASS |
| `nuclear` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `threegen` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `multifam` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `dups` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `missing` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `monomorphic` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `sexchr` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `unrelated` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `admixed` | FAIL | FAIL | FAIL | FAIL | FAIL |
| `singleton` | PASS | PASS | PASS | PASS | PASS |
| `pair` | PASS | PASS | PASS | PASS | PASS |
| `bigish` | FAIL | FAIL | FAIL | FAIL | FAIL |

The three that pass are the three under five samples, where `--ibdseg` becomes `--kinship`
and no `.seg` is written at all.

### 9.4 `params` — flag plumbing and error probes

**220 / 220.** Every `--prefix`, `--cpus`, `--minConc`, `--sexchr`, `--degree`, `--fam`,
`--bim` and error-probe case is byte-identical, on all thirteen datasets.

| dataset | PASS/cases | | dataset | PASS/cases |
|---|---|---|---|---|
| `trio` | 27/27 | | `unrelated` | 14/14 |
| `nuclear` | 13/13 | | `admixed` | 13/13 |
| `threegen` | 13/13 | | `singleton` | 16/16 |
| `multifam` | 21/21 | | `pair` | 15/15 |
| `dups` | 19/19 | | `bigish` | 15/15 |
| `missing` | 13/13 | | `monomorphic` | 13/13 |
| `sexchr` | 28/28 | | **total** | **220/220** |

That includes all eight cases that expect a non-zero exit:

| case | expected |
| --- | --- |
| `params/trio__kinship_famnotfound` | exit 1, `Pedigree file … cannot be opened` |
| `params/trio__kinship_bimnotfound` | exit 1, `Map file … cannot be opened` |
| `params/trio__kinship_prefix_subdir` | exit 1, `Cannot open sub/pre$TMP$.ped to write` — the reference converts the pedigree through a temporary `.ped` named off `--prefix` and opens it **while reading the `.fam`**, before the `.bim` and before the duplicate-sample check |
| `params/trio__kinship_bigbim`, `params/singleton__kinship_bigbim`, `params/multifam__kinship_bigbim`, `params/multifam__kinship_bigfam` | exit 1, `Not enough genotypes at the Nth marker` |
| `params/sexchr__kinship_sexchr1` | exit 1, `Sex chromosome 1 out of range.` |

---

## 10. What is byte-identical, and what is not

The short version, for anyone deciding whether to trust an output file:

| analysis | verdict |
| --- | --- |
| `--kinship` | **byte-identical**, autosomes and X, on all 13 datasets and all 220 `params` variants |
| `--duplicate` | **byte-identical** on all 13 datasets |
| `--bysample` / `--bySNP` | **byte-identical** on all 13 datasets |
| `--cluster` | **byte-identical** on 12 of 13; the merged-cluster tail is unimplemented |
| `--build` | **byte-identical** on 11 of 13 |
| `--unrelated` | **byte-identical** on 12 of 13; `bigish` differs by 3 of 84 kept individuals |
| `--related` | **byte-identical** only on the five datasets the reference downgrades to `--kinship`; the sixteen-column pass is unimplemented |
| `--ibs` | `.ibs`/`.ibs0` are byte-identical in **every column except `MaxIBD2` and `Pr_IBD2`** |
| `--ibdseg` | `allsegs.txt` byte-identical on all 10; `.seg` **is not** — see §11.2 |
| `--autoQC` | **unimplemented** |

## 11. The gaps, measured

Nothing below is a rounding difference. Each entry states what is wrong, how big it is, and
what is known about the missing rule.

### 11.1 The IBD-segment engine — 97 of the 115 failures

One unsolved problem accounts for `--ibdseg` (50 cases), the full `--related` pass (40) and
the two trailing `--ibs` columns (7). The parts that *are* solved are solved exactly:

* **`allsegs.txt` is byte-identical** on all ten datasets that emit one, under every
  `--degree`/`--seglength`/`--prefix` variant. The rule — cut the retained autosomal map at
  each chromosome change and each gap over 1 000 000 bp, then between complete 64-marker
  words of the *global* grid whose 64-gap span exceeds 10 000 000 bp, and keep a piece iff it
  holds ≥ 5 complete words *and* its word-aligned span exceeds 10 000 000 bp — is in
  `king_core::ibdseg::usable_segments`.
* Every **file-existence and console decision** around it is right: the `Segments too short.`
  notice, the extra X-chromosome segment line, the `--seglength` clamp, and the silent
  downgrade of `--ibdseg` to `--kinship` below five samples.

What is wrong is the **IBD1 run-acceptance rule**. Measured on `king.seg`, default flags,
all ten datasets that emit one:

| | count |
| --- | ---: |
| reference rows | 982 |
| rows we emit | 1 170 |
| reference rows we **miss** | **0** |
| rows we emit that the reference does not | 188 (182 of them in `bigish`) |
| reference rows we reproduce **byte for byte** | 558 (56.8 %) |
| rows differing in `IBD1Seg` | 357 |
| rows differing in `IBD2Seg` | 160 |
| rows differing in `PropIBD` | 424 |
| rows differing in `InfType` | 9 |

`FID1`/`ID1`/`FID2`/`ID2` and the row order are exact on every row. The error is
**systematic over-calling of IBD1**: a word is treated as IBD1-eligible iff it contains no
IBS0 at all, which is demonstrably the reference's word test (forcing a single opposite
homozygote anywhere in a word splits the reported segment, swept over 2 000 marker
positions), but the reference then accepts far fewer *runs* of such words than any run-length
or physical-length filter explains. `nuclear`'s `N_C3`/`N_C4` is the clearest case: the pair's
genome is 73 % IBS0-free by word, the reference reports `IBD1Seg 0.0939`, and we report
`0.4548`. `king_core::ibdseg::MIN_RUN1 = 2` is a constant **fitted to the corpus**, not
derived, and it is labelled as such in the source; a grid over (IBS0 tolerance, open, close,
head, tail, minimum words, minimum length) never got a third dataset entirely right, so the
missing rule is structurally different rather than mis-tuned.

Two hypotheses have been tested and **refuted**; do not spend time on them again.

* *A longer fixed run.* Raising `MIN_RUN1` to 3 starts losing reference rows outright
  (−257), so no uniform run-length threshold can be right: the reference reports pairs
  whose only runs are short while rejecting most of a sibling pair's long ones.
* *A per-pair run length scaled by the pair's own IBS0 rate.* The natural "positive
  evidence" test — require a run of `L` IBS0-free words where `L` is the smallest length
  whose chance probability `((1-r)^64)^L` falls below `α`, with `r` the pair's genome-wide
  IBS0 rate — was implemented and swept over `α ∈ {0.5, 0.2, 0.05, 0.01, 0.001, 0.0001}`.
  Its best score is **524** exact rows at `α = 0.01`, worse than the fitted constant's 558
  at every value, and by `α = 0.0001` it starts dropping 233 reference rows. The
  criterion is not a function of the pair's aggregate IBS0 rate in this form.

The same engine drives `--ibs`'s last two columns. There the damage is small and precisely
bounded — over the 21 561 rows of every `.ibs`/`.ibs0` the corpus compares:

| column | rows differing | worst delta |
| --- | ---: | --- |
| every column except the two below | **0** | — |
| `MaxIBD2` | 57 (0.26 %) | 3.5e7 bp, on `bigish` row 129 |
| `Pr_IBD2` | 137 (0.64 %) | 0.160, on `missing` row 6 |

Only rows above the `kinship ≥ 2^-3.5` gate can differ; below it both columns are the
literal `-9`/`0.0000` fillers, and those are exact.

Two other `--ibdseg` items are open and independent of the estimator:

* **`<prefix>splitped.txt` is announced but never written.** It is a pedigree-splitting
  artefact — renaming disconnected families `POOL_S1`…, importing a founder's genotyped
  parents into the referencing family, dropping uninformative families — and it is byte-
  identical under every `--degree`/`--seglength`/`--related` variant, so it can be added
  without touching a single number above. It is a `missing:` note on all 50 failing
  `ibdseg` cases and on none of them is it the only note.
* **`<prefix>X.seg` and its console line are absent.** One case (`sexchr__ibdseg_degree2`)
  is one line short because of it; the X segments themselves are already in `allsegs.txt`
  and that part is exact.

Finally, one `--ibdseg` behaviour is **not reproducible by construction**: the fatal
`Too many first alleles as the major allele (~%.1lf%%)` fires or not on identical input
across repeat runs of the reference (5 runs of one fileset: 2 fatals, 3 successes, with the
percentage varying), because it samples markers with an unseeded RNG. No corpus dataset
triggers it and we do not implement it.

### 11.2 `--related` above the downgrade — 40 cases

Under ten samples the reference replaces `--related` with `--kinship`, and that path is
byte-identical (`trio`, `nuclear`, `missing`, `pair`, `singleton`, on all five `--degree`
variants and under `--ibdseg`). At ten samples and up it writes a **sixteen-column** `.kin`,
a **fourteen-column** `.kin0` and `allsegs.txt`; six of those columns (`HetConc`, `HomIBS0`,
`IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`) come from the engine of §11.1. That pass is
**unimplemented**: the run prints its preamble and writes no files, so the harness reports
`missing:` rather than a numeric difference. Implementing it before the IBD1 rule is settled
would replace one honest failure with a subtler one.

### 11.3 `--autoQC` — 13 cases

**Unimplemented.** The pass is a fixed filter pipeline (SNP call rate < 80 %, monomorphic
SNPs, sample call rate < 95 %, SNP call rate < 95 %, then a gender-QC step that only
`sexchr` reaches) writing `<prefix>_autoQC_Summary.txt`,
`<prefix>_autoQC_snptoberemoved.txt`, `<prefix>_autoQC_sampletoberemoved.txt` and, when
X and Y markers are present, `<prefix>_autoQC_updatesex.txt`. It shares nothing with the
relatedness engines and is the largest self-contained piece of work left.

### 11.4 The merged-cluster tail — 5 cases

Family merging needs a hundred samples before the reference will look at a cross-family pair,
so `bigish` is the only fileset that reaches it and the only one that fails here.

* `apps/bigish__cluster` — `<prefix>updateids.txt` and `<prefix>cluster.kin` are not written
  and the four closing console lines are missing. Everything above them matches.
* `apps/bigish__build` — the same, plus the pedigree-reconstruction log. On the nine datasets
  where nothing merges the reference reconstructs nothing, writes a **zero-byte**
  `<prefix>build.log` and a **zero-byte** `<prefix>updateparents.txt`, writes no
  `updateids.txt` at all despite announcing one, and closes with `No pedigrees can be
  reconstructed.` — all of which we reproduce.
* `apps/bigish__unrelated` and `apps/bigish__unrelated_degree2` — the lists are the right
  **size** (84 kept, 116 removed) and differ by exactly **three individuals**: in each merged
  cluster we keep `B01_F`/`B13_F`/`B25_F` where the reference keeps `B02_F`/`B14_F`/`B26_F`.
  The greedy selection visits family members in ascending count of within-family relatives,
  and ties are broken by a permutation that is a pure function of family size — measured
  from the reference for n = 2…70 and shipped as `TIE_ORDER` — but probe fixtures show the
  scramble is not confined to each tied run, so applying the table per run is an
  approximation. These three rows are where it shows.

### 11.5 One console line, one case

`apps/monomorphic__build` fails on two stdout lines and nothing else:

```
Warning: (P_C3 P_C4) does not look like 1st-degree relatives.
please fix within-family errors first before pedigree recontruction.
```

(the misspelling is the reference's). It is emitted for a within-family pair the pedigree
declares 1st-degree whose genotypes disagree — but **not** on the printed kinship: bisecting
`monomorphic`'s `P_C3`/`P_C4` from `Kinship 0.1477` up to `0.2169` keeps the warning at every
step, while `multifam`'s declared sib pair `B_C1`/`B_C2` at `0.1708` never triggers it. The
predicate is therefore something other than the estimate, and it is not implemented rather
than guessed at.

---

## 12. Related documents

* `docs/SPEC.md` — the implementation specification.
* `tests/parity/generate_corpus.py` — the seeded corpus generator (the datasets themselves
  are gitignored; the generator is the committed artifact).
* `tests/parity/golden/<group>/INDEX.md` — per-group capture notes: every quirk, threshold
  and row-ordering rule the capture agents bisected out of the reference. Read the relevant
  one before chasing a `STRUCTURAL` failure — the answer for row order, silent flag
  downgrades at small *n*, and empty-file conventions is usually already there.
* `tests/parity/golden/<group>/NONDETERMINISTIC.txt` — the raw observations behind §4.
