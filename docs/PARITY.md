# Parity with KING 2.3.2 — the measured claim

This is open-king's authoritative statement of what it reproduces and what it does not.
Every number in it was produced by the commands quoted below, against the reference binary
`king` 2.3.2. Nothing here is an estimate.

> **Headline: 392 of the 480 captured reference invocations are byte-identical.**
> **85 of the 88 that are not** trace to a single cause: the IBD-segment caller places a
> called segment's endpoints within ±1 scan word of the reference, so the four segment
> columns (`IBD1Seg`, `IBD2Seg`, `PropIBD`, and the `InfType` derived from them) are close
> but not equal on 30 % of the rows that carry them. The *set* of pairs reported is exactly
> right everywhere: 0 extra rows and 0 missing rows across all 4 169 captured `.seg` rows.
> **The remaining 3** are structural and named in §6: `<prefix>X.seg` is not written at all
> (2 cases), and `--build`'s pedigree reconstruction produces nothing on the one dataset
> that exercises it (1 case).

---

## 1. Reproducing every number below

```bash
cd /path/to/open-king
cargo build --release

# pass/fail for all 480 cases  -> "392 PASS, 88 FAIL, 480 total"
python3 tests/parity/run_parity.py --impl ./target/release/king

# how big each remaining gap is: rows, columns, mean and worst absolute error
python3 tests/parity/measure_gaps.py --impl ./target/release/king -q

# per-dataset roll-up for one output file
python3 tests/parity/measure_gaps.py --impl ./target/release/king -q --by-dataset king.seg

# harness self-check: the reference against its own captures must be 480/480
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"

# the two differential probes that are not part of the capture corpus
python3 tests/parity/probes/degree_filter.py --ref "/path/to/reference/king"
cd docs/research/fixtures && python3 gate8.py
```

All four are Python 3 standard library only. `run_parity.py` and `measure_gaps.py`
regenerate the input corpus automatically on first run (~20 s) and need no reference
binary; the two probes drive the reference directly, and `gate8.py` reads its path from
`KING` at the top of `docs/research/fixtures/fixlab.py`. `run_parity.py` exits 0 when every
case passed, 1 when at least one failed, 2 on a harness error.

Measured on the tree this document describes:

| command | result |
| --- | --- |
| `run_parity.py --impl target/release/king` | **392 PASS, 88 FAIL, 480 total**, 874 output files byte-compared, 8 diff-excluded |
| `run_parity.py --impl <reference>` | **480 PASS, 0 FAIL** — the normalization is complete and the goldens are self-consistent |
| `probes/degree_filter.py --ref <reference>` | 38 298 cases, **0 false-keep, 0 false-drop** |
| `docs/research/fixtures/gate8.py` | brackets the `--degree 1` IBD2 clause to (0.0789, 0.0812] |
| `cargo test --workspace` | 269 passed, 0 failed, 1 ignored |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --check` | clean |

---

## 2. The matrix: every analysis against every dataset

One cell is every captured invocation of that analysis on that dataset — for `--related`
and `--ibdseg` that is five cases each (bare, four `--degree`/`--seglength` variants).
**Bold** means every captured invocation is byte-identical: every output file, every
column, plus stdout, stderr and exit status.

| analysis | trio | nuclear | threegen | multifam | dups | missing | monomorphic | sexchr | unrelated | admixed | singleton | pair | bigish | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `--kinship` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--duplicate` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--bysample` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--bySNP` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--autoQC` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--cluster` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 12/13 |
| `--build` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 12/13 |
| `--unrelated` | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | 0/2 | 24/26 |
| `--related` | **5/5** | **5/5** | **5/5** | 0/5 | 0/5 | **5/5** | 0/5 | 0/5 | **5/5** | 0/5 | **5/5** | **5/5** | 0/5 | 35/65 |
| `--ibdseg` | **5/5** | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | 0/5 | **5/5** | 0/5 | **5/5** | **5/5** | 0/5 | 20/65 |
| `--ibs` | **1/1** | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | **1/1** | 0/1 | **1/1** | **1/1** | 0/1 | 4/13 |
| flag plumbing + error probes (`params`) | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | **220/220** |
| | | | | | | | | | | | | | | **392/480** |

The `params` group is 220 invocations that exercise the command-line surface rather than
one dataset: every `--prefix` shape, `--cpus`, `--sexchr`, `--degree`, `--minConc`,
`--seglength`, alternate `--fam`/`--bim` inputs, malformed and missing files, and the
banner in each case. All 220 are byte-identical.

Where `--related` and `--ibdseg` pass, it is for one of two reasons, both real:

* **The reference declines to run the segment pass.** Below 10 samples `--related` emits
  the 10-column `--kinship` form (`trio`, `nuclear`, `missing`, `singleton`, `pair`), and
  `trio`/`nuclear`/`threegen`/`missing` additionally get a **zero-byte** `king.kin` under
  `--related` even though stdout announces the file. `--ibdseg` takes the same downgrade
  below 5 samples, which is why `trio` (3), `singleton` (1) and `pair` (2) pass it.
* **The caller is exactly right on that dataset.** `unrelated` reports one pair and gets
  it exact; `threegen` passes `--related` (empty `.kin`) but fails `--ibdseg`.

---

## 3. What is byte-identical, per output file

Counted over every case that produces the file, across all four groups.
`python3 tests/parity/measure_gaps.py --impl target/release/king -q` prints this table.

| output file | cases producing it | byte-identical |
| --- | ---: | --- |
| `<prefix>allsegs.txt` | 163 (+2 custom prefix) | **all** |
| `<prefix>splitped.txt` | 50 | **all** |
| `<prefix>.con` (`--duplicate`) | 46 | **all** |
| `<prefix>bySample.txt` | 15 (+2 custom prefix) | **all** |
| `<prefix>bySNP.txt` | 13 | **all** |
| `<prefix>_autoQC_Summary.txt` | 13 | **all** |
| `<prefix>_autoQC_snptoberemoved.txt` | 13 | **all** |
| `<prefix>_autoQC_sampletoberemoved.txt` | 13 | **all** |
| `<prefix>_autoQC_updatesex.txt` | 1 | **all** |
| `<prefix>updateids.txt` | 2 | **all** |
| `<prefix>X.kin0` | 5 diffable of 13 (see §5.2) | **all 5** |
| `<prefix>.kin` | 187 | 151 — every `--kinship` case; the 36 that differ are `--related` |
| `<prefix>.kin0` | 168 | 160 — same split |
| `<prefix>X.kin` | 15 | 9 — the 6 that differ are `--related`'s three segment columns |
| `<prefix>.ibs` | 13 | 4 |
| `<prefix>.ibs0` | 8 | 6 |
| `<prefix>.seg` | 50 | 5 |
| `<prefix>cluster.kin` | 1 | 0 |
| `<prefix>unrelated.txt` | 26 | 24 |
| `<prefix>unrelated_toberemoved.txt` | 26 | 24 |
| `<prefix>updateparents.txt` | 8 | 7 |
| `<prefix>build.log` | 8 | 7 |
| `<prefix>X.seg` | 2 | 0 — **never written**, see §6.1 |

stdout, stderr and exit status are byte-identical (after the normalization of §7) on all
480 cases: no case in the suite fails on console output alone.

---

## 4. The gaps, measured

Rows are matched on their identifier columns before any comparison, so "extra" and
"missing" below mean a pair only one side reports, never a numeric disagreement. Errors
are measured only over rows both sides report.

### 4.1 The segment columns — 85 of the 88 failures

Every one of these is `IBD1Seg` / `IBD2Seg` / `PropIBD` (and the `InfType` that follows
from them) being close but not equal.

| file | rows differing | of | +extra | −missing | worst column | mean abs err | worst |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| `king.seg` | 1 271 | 4 169 | **0** | **0** | `PropIBD` 1 246 rows | 0.00562 | 0.2109 |
| | | | | | `IBD1Seg` 878 rows | 0.02150 | 0.4218 |
| | | | | | `IBD2Seg` 802 rows | 0.01794 | 0.4218 |
| | | | | | `InfType` 7 rows | — | — |
| `king.kin` (`--related`) | 834 | 3 978 | 0 | 0 | `IBD1Seg` 828 rows | 0.01963 | 0.4218 |
| `king.kin0` (`--related`) | 28 | 296 | 0 | 0 | `IBD1Seg` 28 rows | 0.01461 | 0.0537 |
| `kingX.kin` (`--related`) | 12 | 90 | 0 | 0 | `IBD2Seg` 12 rows | 0.07905 | 0.1460 |
| `kingcluster.kin` | 30 | 165 | 0 | 0 | `IBD2Seg` 30 rows | 0.00600 | 0.0195 |
| `king.ibs` | 148 | 759 | 0 | 0 | `Pr_IBD2` 148 rows | 0.04340 | 0.5116 |
| | 18 | 759 | | | `MaxIBD2` (bp) | 1.42 × 10⁷ | 5.75 × 10⁷ |
| `king.ibs0` | 8 | 19 477 | 0 | 0 | `Pr_IBD2` 8 rows | 0.03648 | 0.0553 |
| `kingunrelated.txt` | 6 lines | 170 | | | membership | — | — |
| `kingunrelated_toberemoved.txt` | 30 lines | 234 | | | membership | — | — |

`MaxIBD2` — the length in base pairs of the single longest IBD2 segment, and the sharpest
per-segment grader available — is exact on **21 542 of 21 560** `.ibs`/`.ibs0` rows; of
the 158 rows where the reference reports a non-zero `MaxIBD2`, **145 are exact**.

`kingunrelated.txt` and `kingcluster.kin` are downstream: both consume `PropIBD`, so a
row that lands on the wrong side of a threshold moves a sample between sets. Both fail
only on `bigish`.

**Scorecard for the primary capture** (`<dataset>__ibdseg`, the default 3 Mb reporting
floor, 982 rows over 10 datasets):

* **820 of 982 rows** have both `IBD1Seg` and `IBD2Seg` exact at the printed four
  decimals. This is the figure `cargo test -p king-core --test ibdseg_parity` prints; run
  it with `KING_GOLDEN=tests/parity/golden` to see the per-dataset breakdown.
* **705 of 982 rows** are exact on all four printed columns. `PropIBD` is computed from
  the unrounded estimates, so it can differ by one ulp of the fourth decimal even where
  both rounded inputs agree; that accounts for the 115-row gap between the two counts.
* `IBD1Seg` alone exact on 822, `IBD2Seg` alone on 822, `InfType` on 981.
* Mean absolute `PropIBD` error **0.00137**, worst 0.2109.
* **0 extra and 0 missing pairs**, on every dataset.

Held out at other reporting floors, rules unchanged: `--seglength 5` gives 697/982 exact
on all four columns (MAE 0.00142), `--seglength 10` gives 665/982 (MAE 0.00160).

Per dataset, on that capture:

| dataset | exact | of | mean abs `PropIBD` err | worst |
| --- | ---: | ---: | ---: | ---: |
| `bigish` | 557 | 763 | 0.00045 | 0.0115 |
| `multifam` | 73 | 104 | 0.00097 | 0.0180 |
| `threegen` | 28 | 39 | 0.00052 | 0.0124 |
| `admixed` | 12 | 16 | 0.00136 | 0.0125 |
| `dups` | 2 | 3 | 0.01260 | 0.0378 |
| `unrelated` | 1 | 1 | 0.00000 | 0.0000 |
| `sexchr` | 8 | 14 | 0.00644 | 0.0341 |
| `nuclear` | 8 | 14 | 0.00371 | 0.0168 |
| `missing` | 8 | 14 | 0.00844 | 0.0548 |
| `monomorphic` | 8 | 14 | 0.04026 | 0.2109 |

The last four are the small ones — each reports the 14 within-family pairs of a single
six-person nuclear family, over 5 000 to 10 000 markers against `bigish`'s 50 000. Read
them with §5 in hand: on those four the reference's own numbers are not the biological
answer either, so an exact-row count there grades very little.

**The residual is two-sided**, which is why it reads as a boundary problem rather than a
missing rule: of the 277 inexact rows, `IBD1Seg` is too high on 139 and too low on 21,
`IBD2Seg` too low on 121 and too high on 39. What is *not* wrong any more is which pairs
get reported, which was the previous release's largest error (188 spurious rows).

### 4.2 `--build` pedigree reconstruction on `bigish` — 1 failure

`apps/bigish__build` writes an **empty** `kingbuild.log` and `kingupdateparents.txt` where
the reference writes 19 and 34 lines. This is not a numeric gap: the reconstruction rules
the reference applies here (`RULE FS0`, `INFERENCE AV.FS` with a `Join3/Join2` statistic,
`INFERENCE HS.UN2`) are unimplemented, so nothing is emitted at all. The other 12 `--build`
datasets are byte-identical because they need none of those rules.

### 4.3 `<prefix>X.seg` — 2 failures

Not written. See §6.1 for everything measured about it.

---

## 5. Known limitations

### 5.1 The four small filesets do not grade the caller

`nuclear` (6 samples, 10 000 markers), `missing` (6, 10 000), `sexchr` (10, 6 000) and
`monomorphic` (12, 5 000) each report the same 14 `.seg` rows: the within-family pairs of
one six-person nuclear family. They are an order of magnitude smaller than `bigish`, and
`monomorphic` additionally cycles every tenth marker through *monomorphic*, *monomorphic
with A1 written as PLINK's missing allele*, *MAF 0.001* and *MAF 0.001 with one forced
carrier* — half its markers carry no usable information. On these the reference's own
segment numbers are nowhere near the pedigree the generator built, so agreeing or
disagreeing with the reference there says little about the caller.

The clearest case is `monomorphic P_C1/P_C2`. `tests/parity/generate_corpus.py` builds
`P_C1`…`P_C4` as the four children of `P_F` and `P_M` by Mendelian transmission, so
`P_C1`/`P_C2` are **full siblings**, expected (π1, π2) ≈ (0.5, 0.25). The reference prints

```
FID1  ID1    FID2  ID2    IBD1Seg  IBD2Seg  PropIBD  InfType
MONO  P_C1   MONO  P_C2   0.9800   0.0000   0.4900   PO
MONO  P_C1   MONO  P_C3   0.9007   0.0000   0.4504   PO
```

— one shared haplotype over 98 % of the genome, never two, and the pair labelled
parent–offspring. Its own `--kinship` on the same fileset gives `P_C1`/`P_C2` a kinship of
**0.3384**, while φ = π2/2 + π1/4 applied to its own segment numbers gives **0.2450**.
Both cannot be right. open-king prints 0.5582 / 0.4218 for that pair, which is not the
truth either (0.42 IBD2 for sibs is as wrong in the other direction). **On this fileset
neither implementation recovers the underlying IBD**, and the 0.2109 worst-case `PropIBD`
error in §4.1 is this row.

Two corrections to earlier project notes, both from
`docs/research/13-informativeness-gate.md`:

* `nuclear` and `missing` were previously described as equally unusable, on the strength
  of `nuclear N_C1`/`N_C3` — the reference prints `IBD1Seg 0.1057`, `IBD2Seg 0.3144`
  where the pedigree implies ≈ 0.5 / 0.25. That gap was **mostly the informativeness gate,
  not a reference error**: the reference was declining to call runs with fewer than ten
  informative markers, and once open-king applies the same rule it prints 0.1240 / 0.2975
  for that pair. `nuclear`'s mean absolute `PropIBD` error is now 0.0037 and `missing`'s
  0.0084, against `monomorphic`'s 0.0403.
* The corresponding claim that these datasets are "poisoned for fitting" therefore holds
  only for `monomorphic`. The rules in the current engine were nonetheless established
  without conditioning on any dataset: no branch in `crates/*/src/` tests a dataset name.
  Dataset names appear in the crates only in `crates/*/tests/`, as the list of datasets a
  scorecard iterates over. (`unrelated` and `sexchr` also appear as *flag* names in the
  CLI table, which is unrelated.)

### 5.2 The reference races on `<prefix>X.kin0`

The between-family X-chromosome writer is not serialised. Six identical invocations of
`king -b sexchr.bed --kinship` produced **six different files** — sizes 665, 662, 662,
662, 187 and 662 bytes — with records torn mid-field and identifier columns from different
pairs interleaved:

```
SEX  S_DAU2  SU3  S_UM  FM  1500  0.323  0.17  0  -0.0351      <- one run
SU3  S_UM    SU4  S_UF  MF  150000  0.331  0.1707  -0.0151     <- another
```

Adding `--cpus 1` makes it deterministic (3 runs, 1 distinct file). No capture made
without `--cpus 1` can be a golden, so the harness excludes `<prefix>X.kin0` from those
cases — 8 of the 13 in the suite, which is the whole of the "8 diff-excluded" line
`run_parity.py` prints. The other 5, all captured with `--cpus 1`, **are** diffed and
**are** byte-identical.

### 5.3 The reference prints an uninitialised value in its own banner

Every run prints `--noscreen [<int>]` where the integer is uninitialised memory. Across
the 484 captured stdouts that print it, it takes three values: `-1717986816` (469 times),
`-858993408` (14) and `-515396096` (1). It is stable within one build and environment, so
the harness normalizes it on both sides rather than trying to reproduce it.

### 5.4 `<prefix>.segments.gz` is never produced

The manual documents it; the 2.3.2 build ships without zlib in its segment writer and
writes no such file on any invocation in the corpus. It is not a parity target.

### 5.5 Scope

Out of scope for v1, and rejected at dispatch rather than silently ignored: `--pca`,
`--mds`, `--roh`, `--lmm`, `--tdt`, `--gdt`, `--risk`, `--makeGRM`, `--plink`, the R
plotting flags (`--rplot`, `--pngplot`, `--rpath`), and multi-dataset input. They are still
*accepted* by the parser so the banner stays byte-exact.

---

## 6. The two structural gaps, in detail

### 6.1 `<prefix>X.seg`

**Not implemented.** What is measured about it, from the reference:

* **When it is written.** `--ibdseg` writes `<prefix>X.seg` exactly when `--degree` is
  non-zero *and* the fileset has usable X-chromosome segments. Bare `--ibdseg`,
  `--ibdseg --degree 0` and `--ibdseg --seglength 5` write none; `--degree` (bare, = 1),
  `--degree 1`, `--degree 2`, `--degree 3` and `--degree -1` all do. Of the 13 corpus
  datasets only `sexchr` has usable X segments, so only 2 of the 480 cases are affected.
* **Which rows.** Exactly the pairs the autosomal `<prefix>.seg` carries, in the same
  order — verified at `--degree -1` (0 rows in both), 1, 2 and 3 (14 rows in both).
* **Its shape is malformed in the reference.** The header names 11 columns
  (`FID1 ID1 FID2 ID2 Sex1 Sex2 MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD`); every data row
  carries 10 tab-separated fields, the last empty. The three numbers written are
  `IBD1Seg`, `IBD2Seg`, `PropIBD` — checked by arithmetic: row `S_SON1`/`S_SON2` reads
  `0.1462  0.6393  0.7124`, and 0.6393 + 0.1462/2 = 0.7124. They land in the `MaxIBD1`,
  `MaxIBD2` and `IBD1Seg` column positions, and the last two columns are never written.
* **Announced on stdout** as `Additional summary statistics of X-Chr IBD segments saved in
  file <prefix>X.seg`, after the autosomal line.

Implementing it needs an X-chromosome segment caller, which would inherit the ±1-word
residual of §4.1; the two affected cases would very likely still not be byte-identical.

### 6.2 `--build` on `bigish`

See §4.2. The missing piece is the pedigree-reconstruction rule set, not any numeric
estimator.

---

## 7. How a case is judged

Each directory under `tests/parity/golden/<group>/<case>/` is one captured reference
invocation: `cmd.txt` (argv, with `{KING}` / `{DATA}` / `{ALT}` placeholders),
`exitcode.txt`, `stdout.txt`, `stderr.txt` verbatim, and every output file the reference
wrote into its working directory. The harness replays `cmd.txt` with our binary in a fresh
temporary directory and compares **exit status, stdout, stderr, the set of files produced,
and the bytes of every file**.

490 invocations are captured; 480 are replayed. The 10 under `core/_analysis/` are
alternate-parameter captures kept for analysis (`--cpus` variants, prefix probes) and are
included only with `--include-analysis`.

Five kinds of line are genuinely non-reproducible and are normalized on **both** sides
before comparison — never on ours alone:

1. **Timestamps** — `ctime()` output on every `… starts at` / `… ends at` line.
2. **Thread counts** — `N CPU cores are used`, which follows the host.
3. **Progress tokens** — the `\r`-separated `10%20%…` sequences, whose granularity
   follows timing.
4. **Absolute input paths** — the reference echoes the `.bed`/`.bim`/`.fam` it was given.
5. **`--noscreen [<int>]`** — §5.3.

The one file excluded from diffing entirely is `<prefix>X.kin0` on captures taken without
`--cpus 1` (§5.2). Running the harness with `--impl <reference>` proves the normalization
is not hiding anything: reference against its own captures is 480/480.

---

## 8. Related documents

* `docs/MAINTAINING.md` — the clean-room rule, repo layout, regenerating the corpus,
  re-capturing goldens, adding an analysis.
* `docs/SPEC.md` — the reference's observable behaviour, flag by flag.
* `docs/BEHAVIOR.md` — raw sweeps behind the rules.
* `docs/VERIFIED_FORMULAS.md` — the estimators, with the experiment that fixed each.
* `docs/research/` — the investigation log. `13-informativeness-gate.md` is the one that
  removed the 188 spurious `.seg` rows; `11-segment-rule-fit.md` §9 describes the ±1-word
  boundary ambiguity that is now the whole residual.
