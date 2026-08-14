# KING 2.3.2 golden capture — applications & QC flag group

Reference binary: `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`, macOS arm64, at
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`.

Datasets: regenerated from the committed seed with
`python3 tests/parity/generate_corpus.py --outdir tests/parity/work/data`
(verified byte-identical to `tests/parity/golden/*.bed` — e.g. `trio.bed`
sha256 `18cafb1a…`, `bigish.bed` sha256 `66748a7d…`).

**91 runs = 13 datasets x 7 flag combinations. Every run exited 0. Every
`stderr.txt` is 0 bytes.** Nothing was skipped for time — the slowest run
(`bigish --autoQC`) took ~1 s.

Each directory `<dataset>__<flagslug>/` contains `cmd.txt` (one argv token per
line, `-b` path templated as `{DATA}/<dataset>.bed`), `stdout.txt`,
`stderr.txt`, `exitcode.txt`, plus every file the run produced.

Line-ending normalization for `stdout.txt` is specified in
[`NONDETERMINISTIC.txt`](NONDETERMINISTIC.txt). **Output files need no
normalization** — all of them are byte-stable across repeat runs and across
`--cpus 1` vs `--cpus 4`.

| flagslug | argv after `-b` |
|---|---|
| `unrelated` | `--unrelated` |
| `unrelated_degree2` | `--unrelated --degree 2` |
| `build` | `--build` |
| `bysample` | `--bysample` |
| `bySNP` | `--bySNP` |
| `autoQC` | `--autoQC` |
| `cluster` | `--cluster` |

---

## 1. Which output columns are affected by `--degree`?

**None. `--degree` does not change a single byte of any output file in this
flag group.**

Swept `--degree 1|2|3|4` against `--unrelated`, and `--degree 1|2|3` against
`--cluster` and `--build`, on `bigish`:

| file | md5 at deg 1 / 2 / 3 / 4 |
|---|---|
| `kingunrelated.txt` | `bc1b4ced…` identical at all four |
| `kingunrelated_toberemoved.txt` | identical |
| `kingcluster.kin` | `307190e5…` identical (166 rows at every degree) |
| `kingupdateids.txt` | `8a48e698…` identical |
| `kingupdateparents.txt` | `aa2e1110…` identical |

Cross-checked across all 13 datasets: `<ds>__unrelated` vs
`<ds>__unrelated_degree2` produce byte-identical `kingunrelated.txt`,
`kingunrelated_toberemoved.txt` and `kingallsegs.txt` everywhere.

`--degree` **does** change three things, all outside the output files:

1. **stdout banner** — `--degree, --noscreen …` becomes `--degree [N], …`, and
   `--degree N` is appended to the `Options in effect:` block.
2. **stdout clustering line + relationship-summary table** —
   `Clustering up to 1st-degree relatives in families...` becomes `2nd`, `3nd`
   (sic — literal typo in the reference), and the
   `Relationship summary (total relatives: 0 by pedigree, N by inference)`
   counts grow dramatically. bigish `--cluster`: `0 0 3 0 0 0` at degree 1 →
   `0 0 3 23 0 0` at degree 2 → `0 226 114 102 52 0` at degree 3.
3. **`--degree >= 3` switches the inference engine and emits an extra file.**
   At degree 1–2 stdout says `N CPU cores are used to compute the pairwise
   kinship coefficients...`; at degree >= 3 it says `N CPU cores are used for
   autosome inference...`, drops the `Autosome genotypes stored in … words` /
   `Sorting autosomes...` lines, and additionally writes **`king.seg`** (IBD
   segment output). The application's own outputs remain identical.

**Implication for open-king:** the unrelated-set extraction, family clustering
result, and pedigree reconstruction all use a *fixed internal* relatedness
threshold; `--degree` only steers the reporting/inference path. Do not wire
`--degree` into the selection logic.

---

## 2. Does `--prefix` change ALL filenames, or only some?

**In this flag group: all of them.** `--prefix` is applied as *bare string
concatenation* with no separator inserted, and the value may contain a path
component. Default prefix is `king`.

`--prefix ZZZ_` on `bigish`:

| flag | default names | with `--prefix ZZZ_` |
|---|---|---|
| `--unrelated` | `kingunrelated.txt`, `kingunrelated_toberemoved.txt`, `kingallsegs.txt` | `ZZZ_unrelated.txt`, `ZZZ_unrelated_toberemoved.txt`, `ZZZ_allsegs.txt` |
| `--build` | `kingbuild.log`, `kingupdateids.txt`, `kingupdateparents.txt`, `kingallsegs.txt` | `ZZZ_build.log`, `ZZZ_updateids.txt`, `ZZZ_updateparents.txt`, `ZZZ_allsegs.txt` |
| `--bysample` | `kingbySample.txt` | `ZZZ_bySample.txt` |
| `--bySNP` | `kingbySNP.txt` | `ZZZ_bySNP.txt` |
| `--autoQC` | `king_autoQC_Summary.txt`, `king_autoQC_sampletoberemoved.txt`, `king_autoQC_snptoberemoved.txt`, `king_autoQC_updatesex.txt` | `ZZZ__autoQC_Summary.txt`, … (**double underscore** — proof of raw concatenation) |
| `--cluster` | `kingcluster.kin`, `kingupdateids.txt`, `kingallsegs.txt` | `ZZZ_cluster.kin`, `ZZZ_updateids.txt`, `ZZZ_allsegs.txt` |

`--prefix sub/pfx` writes `sub/pfxbySample.txt` and `sub/pfxallsegs.txt` — the
prefix is a path prefix, not just a basename prefix, and KING does **not**
create the directory (it must already exist).

### The naming asymmetry

There are two distinct suffix families, and the difference is *only* whether the
suffix begins with a dot:

* **Suffix-concatenated (no dot):** `<prefix>bySample.txt`, `<prefix>bySNP.txt`,
  `<prefix>unrelated.txt`, `<prefix>unrelated_toberemoved.txt`,
  `<prefix>allsegs.txt`, `<prefix>build.log`, `<prefix>updateids.txt`,
  `<prefix>updateparents.txt`, `<prefix>cluster.kin`, `<prefix>_autoQC_*.txt`
* **Dot-extension (prefix behaves as a basename):** `<prefix>.kin`,
  `<prefix>.kin0`, `<prefix>.seg`

Verified directly: `--kinship --prefix ZZZ_` yields `ZZZ_.kin` + `ZZZ_.kin0`;
`--related --prefix ZZZ_` yields `ZZZ_.kin` + `ZZZ_allsegs.txt` **in the same
run** — so `king.kin` and `kingallsegs.txt` are not two conventions across
subcommands but two conventions inside one subcommand. Note also
`<prefix>cluster.kin` is concatenated, *not* `<prefix>.cluster.kin`.

Implement this as a per-output-file literal suffix table, never as
`format!("{prefix}.{ext}")`.

---

## 3. Does `--cpus` change numeric output?

**No.** 24 head-to-head comparisons (`{bigish, admixed, unrelated, multifam}` x
`{--unrelated, --build, --bysample, --bySNP, --autoQC, --cluster}`, `--cpus 1`
vs `--cpus 4`): **every produced file was byte-identical, and the produced file
lists matched.**

Summation order therefore does **not** leak into the results for this flag
group, and open-king is free to choose its own threading decomposition.

`--cpus` does change `stdout.txt`:

* `Computing Parameter : --cpus [N]` and the `--cpus N` line in `Options in effect:`
* `N CPU cores are used to compute the pairwise kinship coefficients...` /
  `Scanning autosomes for QC-by-SNP with N CPU cores...`
* the interleaved `N%` progress tokens — `--cpus 1` on `bigish --cluster` emits
  `7%`…`27%` lines that `--cpus 4` never prints. **Progress output is
  thread-scheduling dependent and must be stripped before diffing.**

Repeat-run determinism of the reference itself was also confirmed: two
back-to-back runs of each of the six applications on `bigish` produced
byte-identical outputs.

---

## 4. Row ORDER of each output file

| file | row order |
|---|---|
| `kingbySample.txt` | **exactly `.fam` order.** Verified on all 13 datasets. 1 header line. |
| `kingbySNP.txt` | **`.bim` order for autosome-only filesets** (verified on 12/13). With sex chromosomes present the file is regrouped: **autosomes (1,2,…) → chr25/XY → X → Y → MT**, with `.bim` order preserved *within* each group. See `sexchr__bySNP`: run lengths `2000×1, 2000×2, 150×25, 1500×X, 300×Y, 50×MT` against `.bim` `2000×1, 2000×2, 1500×23, 300×24, 150×25, 50×26`. Chromosome 23/24/26 are relabelled `X`/`Y`/`MT` in the `Chr` column; **25 stays numeric** and is counted as an autosome (`4150 autosome SNPs (including 150 XY SNPs)`). |
| `kingunrelated.txt` / `kingunrelated_toberemoved.txt` | **neither `.fam` order nor sorted.** Grouped by family. Unmerged families appear in `.fam` first-appearance order; **families that clustering merged into a `KING<n>` family are moved to the END** (still printed with their *original* FIDs) — in `bigish` the `BF01/BF02/BF13/BF14/BF25/BF26` members occupy rows 76–84. Within a family the members are in **greedy-selection order**, e.g. `unrelated`'s POOL family emits `P04 P03 P10 P09 P02 P05 P01 P08 P06 P07`. Single-member families keep `.fam` order. No header line. |
| `kingcluster.kin` | 1 header line, then **all C(n,2) pairs of each new cluster family**, sorted **lexicographically by (FID, ID1, ID2)** — verified exactly against `sort -k1,1 -k2,2 -k3,3`. `bigish` gives 3 clusters × 55 rows (each merges an 5-member + 6-member family = C(11,2)=55). |
| `kingupdateids.txt` | **sorted by (original FID, original IID)** ascending — verified against `sort -k1,1 -k2,2`. No header. Columns: `origFID origIID newFID newIID`. |
| `kingupdateparents.txt` | **same row set and same order as `kingupdateids.txt`**, keyed on the *new* IDs. No header. Columns: `newFID newIID FA MO`. |
| `kingallsegs.txt` | 1 header line, then one row per usable chromosome segment in **ascending chromosome order**, `Segment` column = 1-based row index. |
| `king_autoQC_snptoberemoved.txt` | 1 header (`SNP\tREASON`), then **`.bim` order** among removed SNPs. |
| `king_autoQC_sampletoberemoved.txt` | 1 header (`FID\tIID\tREASON`), then `.fam` order. |
| `king_autoQC_Summary.txt` | fixed step order; the step-2.x gender block is present **only** when X/Y SNPs exist. |
| `king_autoQC_updatesex.txt` | no header; `.fam` order among sex-inferred samples. |
| `kingbuild.log` | family-by-family in new-`KING<n>` order; free text. |

---

## 5. Surprising behavior (all reproducible)

**a. `--build` and `--cluster` are hard-disabled below 10 samples.** stdout
prints `This function is currently disabled for tiny dataset with sample size
< 10.`, exit code is still **0**, and **no files at all are written** — not even
`kingallsegs.txt`. Affects `trio`(3), `nuclear`(6), `missing`(6),
`singleton`(1), `pair`(2) → 10 of the 91 runs produce zero output files.

**b. `--build` announces a file it does not write.** On `threegen`, `multifam`,
`dups`, `monomorphic`, `sexchr`, `unrelated`, `admixed`, stdout says
`Update-ID information is saved in file kingupdateids.txt` but **no
`kingupdateids.txt` exists**. Only `bigish` (where families actually merge)
writes it. Same runs write **0-byte `kingbuild.log` and 0-byte
`kingupdateparents.txt`** while stdout says `No pedigrees can be reconstructed.`
So: file created-and-empty vs never-created is *not* uniform.

**c. `kingbySample.txt` / `kingbySNP.txt` have data-dependent column sets.**
Four `bySample` header variants and three `bySNP` variants, selected by the
**pedigree-declared** PO/trio counts (stdout: `There are N parent-offspring
pairs and M trios, and K full-sibling pairs according to the pedigree.`):

| condition | bySample columns added | bySNP columns added |
|---|---|---|
| 0 PO, 0 trios (`unrelated`, `singleton`, `pair`) | — (8 cols) | — (11 cols) |
| >=1 PO, 0 trios (`dups`) | `N_pair N_MIp Err_MIp MI_Removal` | `N_PO N_HomPO N_errPO Err_InPO Err_InHomPO` |
| >=1 PO, >=1 trio (8 datasets) | + `N_trio N_MIt Err_MIt` | + `N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio` |
| X/Y/MT present (`sexchr`) | inserts `N_xSNP xHeterozygosity N_ySNP N_yHetero N_mtSNP N_mtHetero` **after `Heterozygosity`** | (no change) |

**d. `--autoQC` summary arithmetic does not close, by exactly 16 SNPs, on 7 of
13 datasets.** `Raw − Σ(removed) ≠ Final`:

| dataset | raw | Σ removed | raw−removed | final | delta |
|---|---|---|---|---|---|
| nuclear | 10000 | 3726 | 6274 | 6258 | **−16** |
| threegen | 20000 | 2556 | 17444 | 17428 | **−16** |
| dups | 10000 | 707 | 9293 | 9277 | **−16** |
| missing | 10000 | 4900 | 5100 | 5086 | **−14** |
| unrelated | 20000 | 40 | 19960 | 19944 | **−16** |
| admixed | 20000 | 361 | 19639 | 19623 | **−16** |
| bigish | 50000 | 0 | 50000 | 49984 | **−16** |

(`trio`, `multifam`, `monomorphic`, `sexchr`, `singleton`, `pair` close exactly.)
The most degenerate case is `bigish`: every step reports `(0)` removed,
`king_autoQC_snptoberemoved.txt` is header-only, yet
`Auto-QC step 7: Final check` reports `200 samples, 49984 autosome SNPs`. So 16
SNPs are dropped silently, never listed and never counted. I could not derive
the rule from black-box observation (it does not track sample count,
chromosome count, or divisibility of the surviving SNP count); it is captured
verbatim and must be reproduced bug-compatibly.

**e. `--autoQC` never writes `kingallsegs.txt`**, unlike every other flag in
this group (which writes it whenever n >= 10). Conversely `--bysample` and
`--bySNP` write `kingallsegs.txt` even for **n = 1** (`singleton`), where the
other applications refuse to run at all.

**f. `--autoQC` says "Generate Final Study Files" but produces no PLINK
fileset** — no `.bed/.bim/.fam` is written without `--plink`.

**g. `sexchr` is the only dataset producing `king_autoQC_updatesex.txt`** (2
rows: the two sex-0 samples, both inferred female/`2`), and the only one with
the step-2.1…2.9 gender block in the summary.

**h. Degenerate columns worth asserting on:** `monomorphic__bySNP` has 1639 rows
with `Freq_A` exactly `0.0000`/`1.0000`; `missing__bySNP` has 70 rows that are
entirely zero (`N=0, CallRate=0.0000`); `missing__bysample` shows
`Missing 0.5199` for `M_C3` — those are the intended stress rows.

**i. `--noscreen [-1717986816]` is printed in every banner.** It is an
uninitialized default in the reference, but it is *stable*. Reproduce verbatim;
do not "fix" it.

**j. KING lower-cases flags in the `Options in effect:` block** — `--bySNP` is
echoed as `--bysnp`. The block is ordered by KING's internal check order, not
by command-line order.

---

## 6. Full run inventory (91 runs, all exit 0)

`(NL)` = line count.

| run | output files |
|---|---|
| `trio__unrelated` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(1L)  |
| `trio__unrelated_degree2` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(1L)  |
| `trio__build` | **(no output files)** |
| `trio__bysample` | `kingallsegs.txt`(4L) `kingbySample.txt`(4L)  |
| `trio__bySNP` | `kingallsegs.txt`(4L) `kingbySNP.txt`(5001L)  |
| `trio__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(1830L)  |
| `trio__cluster` | **(no output files)** |
| `nuclear__unrelated` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(4L)  |
| `nuclear__unrelated_degree2` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(4L)  |
| `nuclear__build` | **(no output files)** |
| `nuclear__bysample` | `kingallsegs.txt`(6L) `kingbySample.txt`(7L)  |
| `nuclear__bySNP` | `kingallsegs.txt`(6L) `kingbySNP.txt`(10001L)  |
| `nuclear__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(3727L)  |
| `nuclear__cluster` | **(no output files)** |
| `threegen__unrelated` | `kingallsegs.txt`(22L) `kingunrelated.txt`(5L) `kingunrelated_toberemoved.txt`(7L)  |
| `threegen__unrelated_degree2` | `kingallsegs.txt`(22L) `kingunrelated.txt`(5L) `kingunrelated_toberemoved.txt`(7L)  |
| `threegen__build` | `kingallsegs.txt`(22L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `threegen__bysample` | `kingallsegs.txt`(22L) `kingbySample.txt`(13L)  |
| `threegen__bySNP` | `kingallsegs.txt`(22L) `kingbySNP.txt`(20001L)  |
| `threegen__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(2557L)  |
| `threegen__cluster` | `kingallsegs.txt`(22L)  |
| `multifam__unrelated` | `kingallsegs.txt`(19L) `kingunrelated.txt`(8L) `kingunrelated_toberemoved.txt`(12L)  |
| `multifam__unrelated_degree2` | `kingallsegs.txt`(19L) `kingunrelated.txt`(8L) `kingunrelated_toberemoved.txt`(12L)  |
| `multifam__build` | `kingallsegs.txt`(19L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `multifam__bysample` | `kingallsegs.txt`(19L) `kingbySample.txt`(21L)  |
| `multifam__bySNP` | `kingallsegs.txt`(19L) `kingbySNP.txt`(15001L)  |
| `multifam__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(1437L)  |
| `multifam__cluster` | `kingallsegs.txt`(19L)  |
| `dups__unrelated` | `kingallsegs.txt`(15L) `kingunrelated.txt`(8L) `kingunrelated_toberemoved.txt`(2L)  |
| `dups__unrelated_degree2` | `kingallsegs.txt`(15L) `kingunrelated.txt`(8L) `kingunrelated_toberemoved.txt`(2L)  |
| `dups__build` | `kingallsegs.txt`(15L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `dups__bysample` | `kingallsegs.txt`(15L) `kingbySample.txt`(11L)  |
| `dups__bySNP` | `kingallsegs.txt`(15L) `kingbySNP.txt`(10001L)  |
| `dups__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(708L)  |
| `dups__cluster` | `kingallsegs.txt`(15L)  |
| `missing__unrelated` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(4L)  |
| `missing__unrelated_degree2` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(4L)  |
| `missing__build` | **(no output files)** |
| `missing__bysample` | `kingallsegs.txt`(6L) `kingbySample.txt`(7L)  |
| `missing__bySNP` | `kingallsegs.txt`(6L) `kingbySNP.txt`(10001L)  |
| `missing__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(3L) `king_autoQC_snptoberemoved.txt`(4899L)  |
| `missing__cluster` | **(no output files)** |
| `monomorphic__unrelated` | `kingallsegs.txt`(3L) `kingunrelated.txt`(8L) `kingunrelated_toberemoved.txt`(4L)  |
| `monomorphic__unrelated_degree2` | `kingallsegs.txt`(3L) `kingunrelated.txt`(8L) `kingunrelated_toberemoved.txt`(4L)  |
| `monomorphic__build` | `kingallsegs.txt`(3L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `monomorphic__bysample` | `kingallsegs.txt`(3L) `kingbySample.txt`(13L)  |
| `monomorphic__bySNP` | `kingallsegs.txt`(3L) `kingbySNP.txt`(5001L)  |
| `monomorphic__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(1640L)  |
| `monomorphic__cluster` | `kingallsegs.txt`(3L)  |
| `sexchr__unrelated` | `kingallsegs.txt`(4L) `kingunrelated.txt`(6L) `kingunrelated_toberemoved.txt`(4L)  |
| `sexchr__unrelated_degree2` | `kingallsegs.txt`(4L) `kingunrelated.txt`(6L) `kingunrelated_toberemoved.txt`(4L)  |
| `sexchr__build` | `kingallsegs.txt`(4L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `sexchr__bysample` | `kingallsegs.txt`(4L) `kingbySample.txt`(11L)  |
| `sexchr__bySNP` | `kingallsegs.txt`(4L) `kingbySNP.txt`(6001L)  |
| `sexchr__autoQC` | `king_autoQC_Summary.txt`(18L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(802L) `king_autoQC_updatesex.txt`(2L)  |
| `sexchr__cluster` | `kingallsegs.txt`(4L)  |
| `unrelated__unrelated` | `kingallsegs.txt`(22L) `kingunrelated.txt`(30L) `kingunrelated_toberemoved.txt`(0L)  |
| `unrelated__unrelated_degree2` | `kingallsegs.txt`(22L) `kingunrelated.txt`(30L) `kingunrelated_toberemoved.txt`(0L)  |
| `unrelated__build` | `kingallsegs.txt`(22L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `unrelated__bysample` | `kingallsegs.txt`(22L) `kingbySample.txt`(31L)  |
| `unrelated__bySNP` | `kingallsegs.txt`(22L) `kingbySNP.txt`(20001L)  |
| `unrelated__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(41L)  |
| `unrelated__cluster` | `kingallsegs.txt`(22L)  |
| `admixed__unrelated` | `kingallsegs.txt`(22L) `kingunrelated.txt`(34L) `kingunrelated_toberemoved.txt`(6L)  |
| `admixed__unrelated_degree2` | `kingallsegs.txt`(22L) `kingunrelated.txt`(34L) `kingunrelated_toberemoved.txt`(6L)  |
| `admixed__build` | `kingallsegs.txt`(22L) `kingbuild.log`(0L) `kingupdateparents.txt`(0L)  |
| `admixed__bysample` | `kingallsegs.txt`(22L) `kingbySample.txt`(41L)  |
| `admixed__bySNP` | `kingallsegs.txt`(22L) `kingbySNP.txt`(20001L)  |
| `admixed__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(362L)  |
| `admixed__cluster` | `kingallsegs.txt`(22L)  |
| `singleton__unrelated` | `kingunrelated.txt`(1L) `kingunrelated_toberemoved.txt`(0L)  |
| `singleton__unrelated_degree2` | `kingunrelated.txt`(1L) `kingunrelated_toberemoved.txt`(0L)  |
| `singleton__build` | **(no output files)** |
| `singleton__bysample` | `kingallsegs.txt`(3L) `kingbySample.txt`(2L)  |
| `singleton__bySNP` | `kingallsegs.txt`(3L) `kingbySNP.txt`(5001L)  |
| `singleton__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(3246L)  |
| `singleton__cluster` | **(no output files)** |
| `pair__unrelated` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(0L)  |
| `pair__unrelated_degree2` | `kingunrelated.txt`(2L) `kingunrelated_toberemoved.txt`(0L)  |
| `pair__build` | **(no output files)** |
| `pair__bysample` | `kingallsegs.txt`(3L) `kingbySample.txt`(3L)  |
| `pair__bySNP` | `kingallsegs.txt`(3L) `kingbySNP.txt`(5001L)  |
| `pair__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(2488L)  |
| `pair__cluster` | **(no output files)** |
| `bigish__unrelated` | `kingallsegs.txt`(23L) `kingunrelated.txt`(84L) `kingunrelated_toberemoved.txt`(116L)  |
| `bigish__unrelated_degree2` | `kingallsegs.txt`(23L) `kingunrelated.txt`(84L) `kingunrelated_toberemoved.txt`(116L)  |
| `bigish__build` | `kingallsegs.txt`(23L) `kingbuild.log`(18L) `kingupdateids.txt`(33L) `kingupdateparents.txt`(33L)  |
| `bigish__bysample` | `kingallsegs.txt`(23L) `kingbySample.txt`(201L)  |
| `bigish__bySNP` | `kingallsegs.txt`(23L) `kingbySNP.txt`(50001L)  |
| `bigish__autoQC` | `king_autoQC_Summary.txt`(8L) `king_autoQC_sampletoberemoved.txt`(1L) `king_autoQC_snptoberemoved.txt`(1L)  |
| `bigish__cluster` | `kingallsegs.txt`(23L) `kingcluster.kin`(166L) `kingupdateids.txt`(33L)  |
