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

Fill a cell with `PASS` as our implementation reaches byte-parity for that case; `n/a` means
no case was captured for that combination. Regenerate the truth with
`run_parity.py --json` rather than editing cells by hand once the suite is mostly green.

**Overall:** _not yet run against open-king_ — 0 / 480.

### 9.1 `core` — `--kinship`, `--related`, `--duplicate`, `--ibs`

| dataset | duplicate | ibs | kinship | related | related_degree1 | related_degree2 | related_degree3 | related_degree4 |
|---|---|---|---|---|---|---|---|---|
| `trio` |  |  |  |  |  |  |  |  |
| `nuclear` |  |  |  |  |  |  |  |  |
| `threegen` |  |  |  |  |  |  |  |  |
| `multifam` |  |  |  |  |  |  |  |  |
| `dups` |  |  |  |  |  |  |  |  |
| `missing` |  |  |  |  |  |  |  |  |
| `monomorphic` |  |  |  |  |  |  |  |  |
| `sexchr` |  |  |  |  |  |  |  |  |
| `unrelated` |  |  |  |  |  |  |  |  |
| `admixed` |  |  |  |  |  |  |  |  |
| `singleton` |  |  |  |  |  |  |  |  |
| `pair` |  |  |  |  |  |  |  |  |
| `bigish` |  |  |  |  |  |  |  |  |

### 9.2 `apps` — `--unrelated`, `--build`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster`

| dataset | unrelated | unrelated_degree2 | build | bysample | bySNP | autoQC | cluster |
|---|---|---|---|---|---|---|---|
| `trio` |  |  |  |  |  |  |  |
| `nuclear` |  |  |  |  |  |  |  |
| `threegen` |  |  |  |  |  |  |  |
| `multifam` |  |  |  |  |  |  |  |
| `dups` |  |  |  |  |  |  |  |
| `missing` |  |  |  |  |  |  |  |
| `monomorphic` |  |  |  |  |  |  |  |
| `sexchr` |  |  |  |  |  |  |  |
| `unrelated` |  |  |  |  |  |  |  |
| `admixed` |  |  |  |  |  |  |  |
| `singleton` |  |  |  |  |  |  |  |
| `pair` |  |  |  |  |  |  |  |
| `bigish` |  |  |  |  |  |  |  |

### 9.3 `ibdseg` — `--ibdseg` and its parameters

| dataset | ibdseg | ibdseg_degree2 | ibdseg_seglength5 | ibdseg_seglength10 | related_degree2_ibdseg |
|---|---|---|---|---|---|
| `trio` |  |  |  |  |  |
| `nuclear` |  |  |  |  |  |
| `threegen` |  |  |  |  |  |
| `multifam` |  |  |  |  |  |
| `dups` |  |  |  |  |  |
| `missing` |  |  |  |  |  |
| `monomorphic` |  |  |  |  |  |
| `sexchr` |  |  |  |  |  |
| `unrelated` |  |  |  |  |  |
| `admixed` |  |  |  |  |  |
| `singleton` |  |  |  |  |  |
| `pair` |  |  |  |  |  |
| `bigish` |  |  |  |  |  |

### 9.4 `params` — flag-plumbing and error probes

220 cases over 49 combinations, sparse by design. Cells give the **case count** for that
probe family; record `k/n` as they turn green.

| dataset | baseline | cpus | degree | minConc | prefix | sexchr | alt-input | total |
|---|---|---|---|---|---|---|---|---|
| `trio` | (3) | (3) | (3) | (3) | (5) | (1) | (9) | 27 |
| `nuclear` | (2) | (2) | (3) | (1) | (1) | (1) | (3) | 13 |
| `threegen` | (2) | (2) | (3) | (1) | (1) | (1) | (3) | 13 |
| `multifam` | (3) | (3) | (4) | (1) | (2) | (1) | (7) | 21 |
| `dups` | (2) | (4) | (3) | (5) | (1) | (1) | (3) | 19 |
| `missing` | (2) | (2) | (3) | (1) | (1) | (1) | (3) | 13 |
| `monomorphic` | (2) | (2) | (3) | (1) | (1) | (1) | (3) | 13 |
| `sexchr` | (2) | (3) | (6) | (2) | (2) | (8) | (5) | 28 |
| `unrelated` | (2) | (3) | (3) | (1) | (1) | (1) | (3) | 14 |
| `admixed` | (2) | (2) | (3) | (1) | (1) | (1) | (3) | 13 |
| `singleton` | (2) | (2) | (3) | (1) | (1) | (1) | (6) | 16 |
| `pair` | (2) | (2) | (3) | (1) | (1) | (1) | (5) | 15 |
| `bigish` | (2) | (4) | (3) | (1) | (1) | (1) | (3) | 15 |
| **total** | **28** | **34** | **43** | **20** | **19** | **20** | **56** | **220** |

Eight of these cases expect a **non-zero exit** and are as much a part of parity as the
numeric ones — they pin KING's input-validation behavior:

| case | expected |
| --- | --- |
| `params/trio__kinship_famnotfound` | exit 1, `Pedigree file … cannot be opened` |
| `params/trio__kinship_bimnotfound` | exit 1, `Map file … cannot be opened` |
| `params/trio__kinship_prefix_subdir` | exit 1, `Cannot open sub/pre$TMP$.ped to write` — KING probes writability **before** reading the `.bim` |
| `params/trio__kinship_bigbim`, `params/singleton__kinship_bigbim`, `params/multifam__kinship_bigbim`, `params/multifam__kinship_bigfam` | exit 1, `Not enough genotypes at the Nth marker` |
| `params/sexchr__kinship_sexchr1` | exit 1, `Sex chromosome 1 out of range.` |

Note the asymmetry these pin down: KING validates only the `.bed` **byte length**. A short
`.fam` or short `.bim` is accepted *silently* (fewer samples/SNPs, exit 0); only an
over-long one that grows the required byte width is an error.

---

## 10. Related documents

* `docs/SPEC.md` — the implementation specification.
* `tests/parity/generate_corpus.py` — the seeded corpus generator (the datasets themselves
  are gitignored; the generator is the committed artifact).
* `tests/parity/golden/<group>/INDEX.md` — per-group capture notes: every quirk, threshold
  and row-ordering rule the capture agents bisected out of the reference. Read the relevant
  one before chasing a `STRUCTURAL` failure — the answer for row order, silent flag
  downgrades at small *n*, and empty-file conventions is usually already there.
* `tests/parity/golden/<group>/NONDETERMINISTIC.txt` — the raw observations behind §4.
