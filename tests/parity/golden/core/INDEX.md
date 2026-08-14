# Golden capture — core relatedness flag group (KING 2.3.2)

Reference binary: KING 2.3.2, macOS arm64
(`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`).
Captured 2026-08-13. Host: 16 logical CPUs (KING's OpenMP default reports 8).

**Matrix:** 13 datasets x 8 flag combinations = **104 runs, all exit 0, all stderr empty.**
Nothing was skipped: the whole matrix, including `bigish` (200 samples x 50 000 SNPs),
finished in **under 30 s total**; the slowest single run was ~0.3 s, so the >3 min
skip rule never triggered.

Layout — one directory per run, `<dataset>__<flagslug>/`:

```
cmd.txt        exact argv, with the -b path written as {DATA}/<dataset>.bed
stdout.txt     stderr.txt     exitcode.txt
<every output file the run produced>
```

Replay: `{DATA}` = `tests/parity/work/data`, regenerate with
`python3 tests/parity/generate_corpus.py --outdir tests/parity/work/data`.
Run each `cmd.txt` from inside its own directory (KING writes into `$PWD`).

`_analysis/` holds 10 extra runs backing the analysis below (`--cpus 1/4`,
`--prefix ZZZ_`, `--noscreen`); same capture format. Its `bigish` `--cpus` runs keep
`stdout.txt` plus an `MD5SUMS.txt` but not the output files themselves — those were
verified byte-identical to the core-matrix run and deleted to avoid storing a third copy
of a 2.3 MB `king.ibs0` (see `_analysis/README.txt`). Total capture: 6.0 MB, 681 files.
Nondeterministic stdout lines: **`NONDETERMINISTIC.txt`** in this directory.

---

## 1. Which output columns are affected by `--degree`?

**None. `--degree` changes only which ROWS appear in `king.kin0`, never a column.**

* `king.kin` (within-family) is **byte-identical** across `--related`, `--degree 1/2/3/4`
  (verified by md5 on multifam: `c4f9529a…` for all five runs). Within-family pairs are
  never degree-filtered.
* `king.kin0` (between-family) keeps the same 14 columns
  (`FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType`)
  at every degree; only the row set grows. The cutoff is printed verbatim in stdout and is
  **2^-(degree+1.5)**:

  | degree | printed threshold | bigish kin0 rows |
  |---|---|---|
  | 1 (= default `--related`) | `kinship >= 0.17678` | 3 |
  | 2 | `kinship >= 0.08839` | 26 |
  | 3 | `kinship >= 0.04419` | 59 |
  | 4 | `kinship >= 0.02210` | 60 |

* **The default for bare `--related` is degree 1**, not 2: `bigish__related` and
  `bigish__related_degree1` print the same 0.17678 threshold and are byte-identical.
* Degree also drives three *stdout* things: the `Options in effect` echo, the summary-table
  header (`MZ PO FS 2nd 3rd 4th` for the between-family table vs `… 3rd OTHER` for the
  pedigree one), and a wording bug — KING prints `1st-degree`, `2nd-degree`, then
  **`3nd-degree`, `4nd-degree`**. Reproduce the typo.
* The row set is *not* a pure threshold on the `Kinship` column: multifam degree 3 emits
  54 rows where only 52 pairs have Kinship >= 0.04419, and the file contains two `3rd`-typed
  rows with Kinship 0.0388 / 0.0355 (below cutoff) plus two `4th`-typed rows. Selection is on
  the inferred type / PropIBD, and the count in stdout ("52 pairs … are identified") is the
  *typed* count, not the row count.
* `--degree <= 2` and `--degree >= 3` take **different code paths**: <=2 runs a screening
  stage (stdout: `A subset of informative SNPs will be used to screen close relatives.` +
  `Sorting autosomes...`), >=3 goes straight to full inference. See quirk Q3 — the screening
  path needs 32 768 SNPs and silently finds nothing below that.

## 2. Does `--prefix` change ALL filenames?

**Yes — all of them, and it is a literal string concatenation, not a `king` -> prefix
substitution.** `--prefix ZZZ_` on four different runs (`_analysis/*prefixZZZ*`):

| default name | with `--prefix ZZZ_` | template |
|---|---|---|
| `king.kin` | `ZZZ_.kin` | prefix + `.kin` |
| `king.kin0` | `ZZZ_.kin0` | prefix + `.kin0` |
| `kingX.kin` | `ZZZ_X.kin` | prefix + `X.kin` |
| `kingX.kin0` | `ZZZ_X.kin0` | prefix + `X.kin0` |
| `king.ibs` / `king.ibs0` | `ZZZ_.ibs` / `ZZZ_.ibs0` | prefix + `.ibs`/`.ibs0` |
| `king.con` | `ZZZ_.con` | prefix + `.con` |
| `kingallsegs.txt` | `ZZZ_allsegs.txt` | prefix + `allsegs.txt` |
| `kingbySample.txt` | `ZZZ_bySample.txt` | prefix + `bySample.txt` |

So the asymmetry the task flagged is **only in the suffix table, not in the prefix logic**:
`king` is just the default value of `--prefix`, and some suffixes start with `.`
(`.kin`, `.kin0`, `.ibs`, `.ibs0`, `.con`) while others do not (`X.kin`, `allsegs.txt`,
`bySample.txt`). Implement as `format!("{prefix}{suffix}")` with those exact suffixes; do
**not** special-case the string "king". The names are also echoed inside stdout
("saved in file ZZZ_.kin0"), so the prefix must reach the log strings too.

## 3. Does `--cpus` change numeric output?

**No.** `--cpus 1, 2, 4, 8, 16` on `bigish` x {`--kinship`, `--ibs`, `--related --degree 4`}
produce **byte-identical** `king.kin`, `king.kin0`, `king.ibs`, `king.ibs0` — 15 runs, one md5
per file. They are also identical to the default (8-thread) golden capture. Summation order
does not leak into any printed value, so we are free to choose our own decomposition.

Only three stdout lines move: the `--cpus [N]` parameter echo, the `--cpus N` line in
`Options in effect`, and `N CPU cores are used.` — that last one is **host-dependent** when
`--cpus` is omitted (this 16-core Mac prints 8), so it must be normalized. Note the two
spellings: `--kinship`/`--ibs` print `used.`, `--related` prints `used...`, `--duplicate`
prints two leading spaces then `used...`.

## 4. Row order of each output file

Let `f` = position of the sample in the `.fam` file (0-based). Two different orders are in use.

**a) Within-family files — `king.kin`, `king.ibs`, `kingX.kin`.**
Family blocks appear in order of first appearance of the FID in `.fam`. **Within a family the
individuals are re-sorted lexicographically by IID string** (NOT `.fam` order), and pairs are
the upper triangle of that sorted list. E.g. multifam FAM1 is `A_F A_M A_C1 A_C2 A_C3` in
`.fam` but emits `A_C1-A_C2, A_C1-A_C3, A_C1-A_F, A_C1-A_M, A_C2-A_C3, …, A_F-A_M`.
Verified with zero violations over every `.kin`/`.ibs`/`kingX.kin` in the capture.

**b) Between-family files — `king.kin0`, `king.ibs0`.** Pairs are `(i, j)` with `i < j` in raw
`.fam` index, emitted by a **square-tiled** loop, tile-row major:
`sort key = (i / B, j / B, i, j)`.

* `king.kin0` (from `--kinship`): **B = 32**. Rules out B=8/16 on admixed (40 samples) and
  bigish (200); with n <= 32 the tiling is invisible and the order is plain ascending `(i,j)`,
  which is why multifam/unrelated/dups look unsorted-free.
* `king.ibs0` (from `--ibs`): **B = 8**. Confirmed on all 7 datasets big enough to
  distinguish; plain `(i,j)` order is *wrong* for `.ibs0` even at n = 10.
* `king.kin0` from `--related` never contradicts B = 32 (its row sets are small).
* `king.con` (from `--duplicate`) had at most 2 rows anywhere in the corpus — order is
  consistent with ascending `(i,j)`; not strongly determined by this corpus.
* `kingallsegs.txt` is ordered by segment number, i.e. chromosome then start position, and is
  identical for a given dataset across every flag that writes it.

---

## Run inventory

Row counts exclude the header. `0B` = the file exists but is completely empty (no header).

| dataset | flags | files (rows, excl. header) |
|---|---|---|
| trio | --kinship | king.kin(0B) |
| trio | --related | king.kin(0B) |
| trio | --related --degree 1 | king.kin(0B) |
| trio | --related --degree 2 | king.kin(0B) |
| trio | --related --degree 3 | king.kin(0B) |
| trio | --related --degree 4 | king.kin(0B) |
| trio | --duplicate | king.con(0) |
| trio | --ibs | king.ibs(3) kingallsegs.txt(3) |
| nuclear | --kinship | king.kin(0B) |
| nuclear | --related | king.kin(0B) |
| nuclear | --related --degree 1 | king.kin(0B) |
| nuclear | --related --degree 2 | king.kin(0B) |
| nuclear | --related --degree 3 | king.kin(0B) |
| nuclear | --related --degree 4 | king.kin(0B) |
| nuclear | --duplicate | king.con(0) |
| nuclear | --ibs | king.ibs(15) kingallsegs.txt(5) |
| threegen | --kinship | king.kin(0B) |
| threegen | --related | king.kin(0B) kingallsegs.txt(21) |
| threegen | --related --degree 1 | king.kin(0B) kingallsegs.txt(21) |
| threegen | --related --degree 2 | king.kin(0B) kingallsegs.txt(21) |
| threegen | --related --degree 3 | king.kin(0B) kingallsegs.txt(21) |
| threegen | --related --degree 4 | king.kin(0B) kingallsegs.txt(21) |
| threegen | --duplicate | king.con(0) |
| threegen | --ibs | king.ibs(66) kingallsegs.txt(21) |
| multifam | --kinship | king.kin(40) king.kin0(150) |
| multifam | --related | king.kin(40) kingallsegs.txt(18) |
| multifam | --related --degree 1 | king.kin(40) kingallsegs.txt(18) |
| multifam | --related --degree 2 | king.kin(40) kingallsegs.txt(18) |
| multifam | --related --degree 3 | king.kin(40) king.kin0(54) kingallsegs.txt(18) |
| multifam | --related --degree 4 | king.kin(40) king.kin0(65) kingallsegs.txt(18) |
| multifam | --duplicate | king.con(0) |
| multifam | --ibs | king.ibs(40) king.ibs0(150) kingallsegs.txt(18) |
| dups | --kinship | king.kin(2) king.kin0(43) |
| dups | --related | king.kin(2) kingallsegs.txt(14) |
| dups | --related --degree 1 | king.kin(2) kingallsegs.txt(14) |
| dups | --related --degree 2 | king.kin(2) kingallsegs.txt(14) |
| dups | --related --degree 3 | king.kin(2) king.kin0(1) kingallsegs.txt(14) |
| dups | --related --degree 4 | king.kin(2) king.kin0(1) kingallsegs.txt(14) |
| dups | --duplicate | king.con(2) |
| dups | --ibs | king.ibs(2) king.ibs0(43) kingallsegs.txt(14) |
| missing | --kinship | king.kin(0B) |
| missing | --related | king.kin(0B) |
| missing | --related --degree 1 | king.kin(0B) |
| missing | --related --degree 2 | king.kin(0B) |
| missing | --related --degree 3 | king.kin(0B) |
| missing | --related --degree 4 | king.kin(0B) |
| missing | --duplicate | king.con(0) |
| missing | --ibs | king.ibs(15) kingallsegs.txt(5) |
| monomorphic | --kinship | king.kin(15) king.kin0(51) |
| monomorphic | --related | king.kin(15) kingallsegs.txt(2) |
| monomorphic | --related --degree 1 | king.kin(15) kingallsegs.txt(2) |
| monomorphic | --related --degree 2 | king.kin(15) kingallsegs.txt(2) |
| monomorphic | --related --degree 3 | king.kin(15) king.kin0(0) kingallsegs.txt(2) |
| monomorphic | --related --degree 4 | king.kin(15) king.kin0(1) kingallsegs.txt(2) |
| monomorphic | --duplicate | king.con(0) |
| monomorphic | --ibs | king.ibs(15) king.ibs0(51) kingallsegs.txt(2) |
| sexchr | --kinship | king.kin(15) king.kin0(30) kingX.kin(15) kingX.kin0(2) |
| sexchr | --related | king.kin(15) kingX.kin(15) kingallsegs.txt(3) |
| sexchr | --related --degree 1 | king.kin(15) kingX.kin(15) kingallsegs.txt(3) |
| sexchr | --related --degree 2 | king.kin(15) kingX.kin(15) kingallsegs.txt(3) |
| sexchr | --related --degree 3 | king.kin(15) king.kin0(0) kingX.kin(15) kingX.kin0(0) kingallsegs.txt(3) |
| sexchr | --related --degree 4 | king.kin(15) king.kin0(1) kingX.kin(15) kingX.kin0(1) kingallsegs.txt(3) |
| sexchr | --duplicate | king.con(0) |
| sexchr | --ibs | king.ibs(15) king.ibs0(30) kingallsegs.txt(3) |
| unrelated | --kinship | king.kin(45) king.kin0(390) |
| unrelated | --related | king.kin(45) kingallsegs.txt(21) |
| unrelated | --related --degree 1 | king.kin(45) kingallsegs.txt(21) |
| unrelated | --related --degree 2 | king.kin(45) kingallsegs.txt(21) |
| unrelated | --related --degree 3 | king.kin(45) kingallsegs.txt(21) |
| unrelated | --related --degree 4 | king.kin(45) king.kin0(0) kingallsegs.txt(21) |
| unrelated | --duplicate | king.con(0) |
| unrelated | --ibs | king.ibs(45) king.ibs0(390) kingallsegs.txt(21) |
| admixed | --kinship | king.kin(18) king.kin0(762) |
| admixed | --related | king.kin(18) kingallsegs.txt(21) |
| admixed | --related --degree 1 | king.kin(18) kingallsegs.txt(21) |
| admixed | --related --degree 2 | king.kin(18) kingallsegs.txt(21) |
| admixed | --related --degree 3 | king.kin(18) king.kin0(0) kingallsegs.txt(21) |
| admixed | --related --degree 4 | king.kin(18) king.kin0(2) kingallsegs.txt(21) |
| admixed | --duplicate | king.con(0) |
| admixed | --ibs | king.ibs(18) king.ibs0(762) kingallsegs.txt(21) |
| singleton | --kinship | king.kin0(0) |
| singleton | --related | king.kin0(0) |
| singleton | --related --degree 1 | king.kin0(0) |
| singleton | --related --degree 2 | king.kin0(0) |
| singleton | --related --degree 3 | king.kin0(0) |
| singleton | --related --degree 4 | king.kin0(0) |
| singleton | --duplicate | king.con(0) |
| singleton | --ibs | king.ibs(0) kingallsegs.txt(2) |
| pair | --kinship | king.kin0(1) |
| pair | --related | king.kin0(1) |
| pair | --related --degree 1 | king.kin0(1) |
| pair | --related --degree 2 | king.kin0(1) |
| pair | --related --degree 3 | king.kin0(1) |
| pair | --related --degree 4 | king.kin0(1) |
| pair | --duplicate | king.con(0) |
| pair | --ibs | king.ibs(0) king.ibs0(1) kingallsegs.txt(2) |
| bigish | --kinship | king.kin(573) king.kin0(19327) |
| bigish | --related | king.kin(573) king.kin0(3) kingallsegs.txt(22) |
| bigish | --related --degree 1 | king.kin(573) king.kin0(3) kingallsegs.txt(22) |
| bigish | --related --degree 2 | king.kin(573) king.kin0(26) kingallsegs.txt(22) |
| bigish | --related --degree 3 | king.kin(573) king.kin0(59) kingallsegs.txt(22) |
| bigish | --related --degree 4 | king.kin(573) king.kin0(60) kingallsegs.txt(22) |
| bigish | --duplicate | **(none)** |
| bigish | --ibs | king.ibs(573) king.ibs0(19327) kingallsegs.txt(22) |

---

## Surprising behavior (parity landmines)

**Q1 — `--related` is silently rewritten to `--kinship` when the fileset has < 10 samples.**
`trio` (3), `nuclear` (6), `missing` (6), `pair` (2), `singleton` (1) all print
`Options in effect:` / `--kinship` even though `--related` was passed. Bisected with PLINK
subsets of multifam: n=7,8,9 downgrade; n=10,11 do not. **The threshold is exactly n >= 10.**
Under the downgrade the run produces the kinship schema (8-column `.kin0`, `.kin` with
`Error`, no `InfType`/`IBD1Seg`/`IBD2Seg`/`PropIBD`), writes no `kingallsegs.txt`, and
**`--degree` is ignored entirely** (sub9: `--kinship`, `--related`, `--related --degree 1`
and `--degree 4` all emit the same 20 `.kin0` rows).

**Q2 — a single-family fileset writes a 0-byte `king.kin`** (no header, no rows) for both
`--kinship` and `--related`: `trio`, `nuclear`, `threegen`, `missing`. stdout still says
`Within-family kinship data saved in file king.kin` and then `There is only one family.`.
`--ibs` on the same input writes a full `king.ibs` with header and rows, so this is specific
to the `.kin` writer. Multi-FID datasets write `.kin` normally.

**Q3 — the `--degree <= 2` cross-family path needs >= 32 768 SNPs or it finds nothing.**
stdout on bigish reads `Stages 1&2 (with 32768 SNPs): 18 pairs of relatives are detected`.
Every other dataset in the corpus has <= 20 000 SNPs, so the whole screening block is skipped
and KING reports `No close relatives are inferred.` — even for `dups`, whose cross-family
exact duplicate has Kinship 0.5000 in the `--kinship` capture, and for `multifam`, where
8 cross-family pairs exceed the degree-1 cutoff. **`--noscreen` does not fix it**
(`_analysis/multifam__related_degree1__noscreen` still emits nothing), so this is not the
screening filter but the stage itself being skipped. Consequence: `--related`/`--degree 1`/
`--degree 2` produce **no `king.kin0` at all** on 12 of 13 datasets. `--degree >= 3` takes the
non-screening path and finds them normally.

**Q4 — `kingX.kin0` from `--kinship` is malformed.** `sexchr__kinship/kingX.kin0` has a
9-column header (`FID1 ID1 FID2 ID2 Sex N_SNP Het IBS0 KinshipX`) but **10-field data rows**
(`SU3 S_UM SU4 S_UF M FM 15 0.365 0.1627 0.0547`), and **the single pair is printed twice**,
byte-identically. `N_SNP` shows `15` where `kingX.kin` shows `1500` for the same data. The
`--related --degree 3/4` variant of the same file is well formed and uses a different schema
(`FID1 ID1 FID2 ID2 Sex1 Sex2 IBD1Seg IBD2Seg PropIBD`). Byte parity here means reproducing
the duplication and the column misalignment.

**Q5 — `--duplicate` writes no file at all when nothing is found on a large set, but a
header-only file on a small one.** `bigish__duplicate` produced **zero output files**;
`trio`/`nuclear`/`multifam`/`admixed`/`unrelated`/... produced a header-only `king.con`.
Both print `No duplicates are found with heterozygote concordance rate > 80%.`. The small-set
runs additionally print `N additional pairs from screening stage not confirmed in the final
stage` (trio: 3), so `--duplicate` also has a two-stage screen with its own threshold.

**Q6 — `--ibs` emits the sentinel `-9` in `MaxIBD2` and `Pr_IBD2`** for unrelated pairs in
`.ibs0` (19 303 of 19 327 rows on bigish; 390 of 390 on `unrelated`). `.ibs` (within-family)
prints `0.000` / real values instead. Do not treat `-9` as a parse error.

**Q7 — `--kinship` writes `kingX.kin`/`kingX.kin0` only when chromosome 23 is present**
(`sexchr` only); `--ibs` never writes X output at all, even on `sexchr`.
`kingallsegs.txt` is written by `--ibs` on every dataset but by `--related` only when the
run was not downgraded per Q1.

**Q8 — `missing` (6 samples, one family) is the null case for this whole group:** `.kin` is
0 bytes, no `.kin0` exists, and only `--ibs`/`--duplicate` produce readable output. Its
`.ibs` `N_SNP` column does vary per pair (3 816 … 9 732), which is the useful signal there —
pairwise-complete SNP counts, not a global count.

**Q9 — `monomorphic` and `admixed` are numerically well behaved**, no NaN/Inf anywhere in the
capture; the only sentinel is `-9` (Q6). `admixed` cross-population unrelated pairs carry
strongly negative Kinship, which is expected KING-robust behavior, not a defect.

## Cheap sanity gate for the Rust port

1. Load + write the `--kinship` capture for `multifam` (has both `.kin` and `.kin0`) — this
   pins the within-family lexicographic sort, the 32-tile order, and all numeric formatting.
2. `dups --duplicate` — pins `.con`, exact-duplicate `1.00000` and the 0.2 %-error MZ row.
3. `multifam --ibs` — pins the 8-tile order and the 21-column `.ibs` schema.
4. `bigish --related --degree 3` — the only dataset that exercises the real between-family
   inference path end to end.
5. `trio --related` — pins Q1 + Q2 (downgrade to `--kinship`, 0-byte `.kin`).
