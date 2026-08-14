# Golden capture: `--ibdseg` flag group (KING 2.3.2)

Reference binary: `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2, macOS arm64)
Datasets: `/Users/wgu/Desktop/open-king/tests/parity/golden/*.{bed,bim,fam}` (regenerate with
`python3 /Users/wgu/Desktop/open-king/tests/parity/generate_corpus.py --outdir <dir>`)

65 runs = 13 datasets x 5 flag combinations. **All exited 0; every `stderr.txt` is empty (0 bytes).**
Nothing was skipped: the slowest run (`bigish`, 200 samples x 50k SNPs) took 0.09 s.

Each run directory holds `cmd.txt`, `stdout.txt`, `stderr.txt`, `exitcode.txt` plus every file the
run produced. `cmd.txt` uses two placeholders: `{KING}` = the reference binary, `{DATA}` = the
directory holding the `.bed/.bim/.fam`. Runs were executed with the run directory as CWD.

Timestamp / core-count / progress-bar lines in `stdout.txt` are not reproducible; the patterns to
normalize are listed in `NONDETERMINISTIC.txt` (same directory).

---

## The four analysis questions

### 1. Which output columns are affected by `--degree`, and how?

**None.** `--degree` never changes a column *value*; it changes which **rows** are emitted, and (on
X-bearing data) whether an extra **file** is written.

* `king.seg` — `--degree N` keeps only rows with **`PropIBD > 2^-(N+0.5)`**. Verified exactly
  (predicted row set == actual row set, no mismatches) on all 10 datasets that produce a `.seg`.
  Row counts on `threegen` (39 rows unfiltered): `--degree 1` -> 18, `2` -> 29, `3` -> 37, `4` -> 39,
  `5` -> 39. The surviving numbers are bit-identical to the unfiltered run.
* The cutoff is the same one `--related` reports for `king.kin0` as a **kinship** threshold
  `2^-(N+1.5)` (KING prints `kinship >= 0.04419` for `--degree 3`); since `PropIBD ~= 2 * kinship`
  the two are the same rule.
* The filter is on `PropIBD`, **not** on the `InfType` label. They disagree at the boundary: on
  `threegen --degree 1`, row `TG_GM1 TG_C4` (`PropIBD 0.3576 > 2^-1.5 = 0.35355`) survives while
  still being **labelled `2nd`**. `InfType` is assigned from the `IBD1Seg`/`IBD2Seg` pattern
  (`PO` needs `IBD1Seg ~ 1`, `FS` needs `IBD2Seg > 0`), not from `PropIBD` alone.
* `unrelated` degenerates to a header-only `king.seg` under `--degree 2` (its one reported pair has
  `PropIBD 0.0061`).
* **Surprise:** on `sexchr`, `--degree <anything>` additionally emits **`kingX.seg`**, which plain
  `--ibdseg` never writes (confirmed: 3 repeat plain runs, no `kingX.seg`; `--degree 1,2,3,4,5` all
  write it). stdout gains the line `Additional summary statistics of X-Chr IBD segments saved in
  file kingX.seg`. The *value* of `--degree` is irrelevant to this - only its presence.
* `kingallsegs.txt` and `kingsplitped.txt` are **byte-identical** under every `--degree` /
  `--seglength` / `--related` variant.

For contrast, `--seglength` is the flag that *does* move numbers: it rewrites `IBD1Seg`, `IBD2Seg`
and `PropIBD` (and can flip `InfType`) without changing the row set, and adds
`Minimum segment length is set as <N>000000 bp` to stdout. Trap: the stdout note block still says
`Short IBD segments (<3Mb) are not reported/utilized` even at `--seglength 10` - that text is
hard-coded and must be reproduced verbatim.

### 2. Does `--prefix` change ALL filenames, or only some?

**All of them** - `--prefix` replaces the literal stem `king`, it is not a separate token. With
`--prefix ZZ_`:

| default | with `--prefix ZZ_` | join |
|---|---|---|
| `king.seg` | `ZZ_.seg` | `<prefix>` + `.seg` |
| `kingX.seg` | `ZZ_X.seg` | `<prefix>` + `X.seg` |
| `king.kin` | `ZZ_.kin` | `<prefix>` + `.kin` |
| `kingX.kin` | `ZZ_X.kin` | `<prefix>` + `X.kin` |
| `kingallsegs.txt` | `ZZ_allsegs.txt` | `<prefix>` + `allsegs.txt` |
| `kingsplitped.txt` | `ZZ_splitped.txt` | `<prefix>` + `splitped.txt` |

So the naming asymmetry is real but it is **not** a per-file rule: the concatenation is always
`prefix + suffix`, and it only *looks* asymmetric because some suffixes begin with `.` (`.kin`,
`.seg`) and others do not (`allsegs.txt`, `splitped.txt`, `bySample.txt`). A reimplementation must
model this as string concatenation with a default prefix of `"king"`, never as
`prefix + "." + ext`. stdout echoes the prefixed names too. (Aside: `--prefix out/` with a
non-existent directory silently produces **no output files at all** and still exits 0.)

### 3. Does `--cpus` change numeric output?

**No - not one byte.** `--cpus 1 / 2 / 4 / 8 / 16` on `admixed` (40 samples) and `--cpus 1 / 4` on
`bigish` (200 samples) produced **byte-identical** `king.seg`, `kingallsegs.txt` and
`kingsplitped.txt`. Five repeat runs of `bigish --ibdseg` at the default thread count are likewise
identical to each other and to the committed golden file. **We do not have to match KING's threading
decomposition to reach byte parity on output files.**

`--cpus` *does* change `stdout.txt` in two places, both listed in `NONDETERMINISTIC.txt`:
the `N CPU cores are used...` line, and the length of the progress bar
(1 thread -> `0%...90%`, 4 -> `0%...21%`, 8 -> `0%...10%`), which is emitted without a trailing
newline and therefore ends up glued to the following `ends at <timestamp>` line.

### 4. What is the row ORDER of each output file?

* **`king.seg` - a 16x16 tiled pair sweep, not a plain nested loop.** With `i`, `j` the 0-based
  `.fam` row indices, rows are ordered by **`(i/16, j/16, i, j)`** (integer division), i.e. KING
  walks 16-sample blocks (block-row ascending, block-column ascending, upper triangle only) and
  within each tile emits `i < j` lexicographically. Verified as an exact sort key on **all 48
  non-empty `king.seg` files** in this capture (every dataset x every flag combination), including
  `bigish` (763 rows) and `multifam` (104 rows); tile sizes 8 and 32 both fail on those two, so the
  constant is exactly **16**. This is invisible on any dataset with <= 16 samples
  (where it collapses to plain lexicographic order) - `multifam` (20) and `bigish` (200) are the
  only datasets in the corpus that expose it. A naive `for i { for j>i }` implementation will emit
  the correct *set* of rows in the wrong *order*.
* **`kingX.seg`** - the same pairs, in the same order, as the `king.seg` of the same run
  (`sexchr` has 10 samples, so the tiling is not exercised there).
* **`kingallsegs.txt`** - one row per usable chromosomal segment, ordered by the `Segment` column
  (1..N), which is chromosome-ascending. Independent of every flag in this group.
* **`kingsplitped.txt`** - `.fam` order of the *rewritten* pedigree; KING appends synthesised
  founders (e.g. `POFAM KING1`) and renames FIDs when it splits a family
  (`MZFAM MZ_1` -> `MZFAM_S1`), so it is not simply the input `.fam` order. Space-separated,
  9 columns, **no header**.
* **`king.kin`** (from `--related`, and from the small-sample fallback) - FID blocks are contiguous
  and appear in `.fam` FID order, but **within a family the members are sorted alphabetically by
  individual ID string**, *not* in `.fam` order, and pairs are then emitted `i < j` over that
  alphabetical order. E.g. `bigish` family `BF01` is `B01_F B01_M B01_C1 B01_C2 B01_C3` in the
  `.fam` but its `king.kin` rows start `B01_C1 x B01_C2`, `B01_C1 x B01_C3`, `B01_C1 x B01_F`, ...
  Verified on all 4 `multifam` families and all 29 `bigish` families that emit rows; 298 of the 573
  `bigish` rows have `idx(ID1) > idx(ID2)` in `.fam` terms, so this is not a subtle effect.
* **`king.kin0`** - `bigish` (the only dataset here with `.kin0` rows from a real run: 26 rows)
  is consistent with both plain `(i, j)` lexicographic order and the 16-tile order; the two are
  indistinguishable on this row set, so the tiling is unconfirmed for `.kin0`.

---

## Behaviours a reimplementation must reproduce

1. **`--ibdseg` silently degrades to `--kinship` when N <= 4.** Verified by bisection with PLINK
   `--keep` subsets of `nuclear`: N=2,3,4 print
   `--kinship analysis carried out instead for such a small sample size.` and write `king.kin` /
   `king.kin0`; N=5,6,7 run the real analysis and write `king.seg` + `kingallsegs.txt` +
   `kingsplitped.txt`. In this corpus that hits **`trio` (3), `pair` (2), `singleton` (1)** - none of
   them produces a `.seg` file for any of the five flag combinations. With `--related --degree 2
   --ibdseg` the message instead reads `--related is replaced with --kinship for a small sample
   size.` and the whole kinship block is printed **twice** (once for the `--related` stage, once for
   the `--ibdseg` stage).
2. **`trio__*/king.kin` is 0 bytes** - no header at all, the known single-family quirk, here reached
   through the small-sample fallback path.
3. **`--related --degree 2 --ibdseg` runs two analyses back to back**, printing two
   `Options in effect:` blocks (`--related` + `--degree 2`, then `--ibdseg` + `--degree 2`). Its
   `king.seg` is **byte-identical to the corresponding `ibdseg_degree2` run** on all 10 datasets,
   and its `king.kin`/`king.kin0`/`kingallsegs.txt` are byte-identical to `--related --degree 2`
   run alone - the two stages do not influence each other. `--ibdseg` adds `kingsplitped.txt` and
   `king.seg`; it does **not** add columns to `.kin`/`.kin0` (KING 2.3.2's `--related` already
   emits `IBD1Seg`/`IBD2Seg`/`PropIBD`/`InfType` there on its own).
4. **`--seglength` prints a stray `.` prefix**: the line after the banner is
   `.Loading genotype data in PLINK binary format...` (leading dot) instead of
   `Loading genotype data...`.
5. **No `.segments.gz` is ever written by `--ibdseg` in 2.3.2.** Per-segment detail is simply not
   produced; adding `--rplot` yields `king_ibd1vsibd2.{R,Rout,pdf,ps}` but still no segment dump.

## Degenerate / suspicious output worth knowing about

* **`kingX.seg` is malformed.** Header has **11** tab-separated fields
  (`FID1 ID1 FID2 ID2 Sex1 Sex2 MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD`) but every data row has
  only **9 values plus a trailing tab** (10 fields). Two columns of data are missing and everything
  after `Sex2` is misaligned against the header. `kingX.kin` (written by
  `--related --degree 2 --ibdseg`) is well formed by contrast: 9 header fields, 9 data fields.
* **`multifam --related --degree 2` reports `No close relatives are inferred.` and writes no
  `king.kin0`** - in the *same run* whose `king.seg` lists 25 cross-family pairs at 2nd degree or
  closer (including `FAM1 A_F` x `FAM2 B_F` as `FS`, `PropIBD 0.4887`). The across-family screen
  ("A subset of informative SNPs will be used to screen close relatives") misses them at degree 1
  and 2 and only finds them at `--degree 3`. The `.kin0` screen and the `.seg` engine disagree.
* **`unrelated`**: `king.seg` has exactly **one** row out of 435 possible pairs
  (`POOL P04` x `SNG19 S19`, `PropIBD 0.0061`, `InfType UN`) - a pure false positive of the
  ">10Mb long segment" pair filter, and it vanishes under `--degree 2`.
* **`monomorphic`**: full sibs are called **`PO`** (`P_C1` x `P_C2`, `IBD1Seg 0.9800`,
  `IBD2Seg 0.0000`) and `P_C3` x `P_C4` is called `2nd`. The 500 monomorphic / 500 `A1='0'` /
  500 MAF-0.001 SNPs collapse the IBD2 signal. This is the dataset most likely to expose a
  divide-by-zero or allele-frequency edge case.
* **`sexchr`**: `S_SON2` x `S_DAU1` and `S_SON2` x `S_DAU2` (full sibs) are called **`PO`** with
  `PropIBD` 0.7274 / 0.6436 - i.e. `PropIBD > 0.5` for a non-duplicate pair, driven by the haploid
  X/Y/MT markers leaking into the autosomal estimate.
* **`dups`**: clean and useful as an anchor - `DUP_A` x `DUP_A_COPY` = `IBD2Seg 1.0000`,
  `PropIBD 1.0000`, `Dup/MZ`; the 0.2%-error MZ pair = `IBD2Seg 0.9223`, `PropIBD 0.9441`. Only 3
  rows: the 4 unrelated samples never clear the long-segment filter.
* **`admixed`**: only 16 rows for 40 samples; the cross-population unrelated pairs are filtered out
  entirely rather than reported with negative values, which is how `--ibdseg` differs from
  `--kinship` on structured data.

---

## Every run

| run directory | exit | output files (rows = data rows, excl. header) |
|---|---|---|
| `trio__ibdseg` | 0 | `king.kin`(0 bytes) |
| `trio__ibdseg_degree2` | 0 | `king.kin`(0 bytes) |
| `trio__ibdseg_seglength5` | 0 | `king.kin`(0 bytes) |
| `trio__ibdseg_seglength10` | 0 | `king.kin`(0 bytes) |
| `trio__related_degree2_ibdseg` | 0 | `king.kin`(0 bytes) |
| `nuclear__ibdseg` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `nuclear__ibdseg_degree2` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `nuclear__ibdseg_seglength5` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `nuclear__ibdseg_seglength10` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `nuclear__related_degree2_ibdseg` | 0 | `king.kin`(0 bytes), `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `threegen__ibdseg` | 0 | `king.seg`(39), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `threegen__ibdseg_degree2` | 0 | `king.seg`(29), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `threegen__ibdseg_seglength5` | 0 | `king.seg`(39), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `threegen__ibdseg_seglength10` | 0 | `king.seg`(39), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `threegen__related_degree2_ibdseg` | 0 | `king.kin`(0 bytes), `king.seg`(29), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `multifam__ibdseg` | 0 | `king.seg`(104), `kingallsegs.txt`(18), `kingsplitped.txt`(21) |
| `multifam__ibdseg_degree2` | 0 | `king.seg`(68), `kingallsegs.txt`(18), `kingsplitped.txt`(21) |
| `multifam__ibdseg_seglength5` | 0 | `king.seg`(104), `kingallsegs.txt`(18), `kingsplitped.txt`(21) |
| `multifam__ibdseg_seglength10` | 0 | `king.seg`(104), `kingallsegs.txt`(18), `kingsplitped.txt`(21) |
| `multifam__related_degree2_ibdseg` | 0 | `king.kin`(40), `king.seg`(68), `kingallsegs.txt`(18), `kingsplitped.txt`(21) |
| `dups__ibdseg` | 0 | `king.seg`(3), `kingallsegs.txt`(14), `kingsplitped.txt`(4) |
| `dups__ibdseg_degree2` | 0 | `king.seg`(3), `kingallsegs.txt`(14), `kingsplitped.txt`(4) |
| `dups__ibdseg_seglength5` | 0 | `king.seg`(3), `kingallsegs.txt`(14), `kingsplitped.txt`(4) |
| `dups__ibdseg_seglength10` | 0 | `king.seg`(3), `kingallsegs.txt`(14), `kingsplitped.txt`(4) |
| `dups__related_degree2_ibdseg` | 0 | `king.kin`(2), `king.seg`(3), `kingallsegs.txt`(14), `kingsplitped.txt`(4) |
| `missing__ibdseg` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `missing__ibdseg_degree2` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `missing__ibdseg_seglength5` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `missing__ibdseg_seglength10` | 0 | `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `missing__related_degree2_ibdseg` | 0 | `king.kin`(0 bytes), `king.seg`(14), `kingallsegs.txt`(5), `kingsplitped.txt`(5) |
| `monomorphic__ibdseg` | 0 | `king.seg`(14), `kingallsegs.txt`(2), `kingsplitped.txt`(5) |
| `monomorphic__ibdseg_degree2` | 0 | `king.seg`(14), `kingallsegs.txt`(2), `kingsplitped.txt`(5) |
| `monomorphic__ibdseg_seglength5` | 0 | `king.seg`(14), `kingallsegs.txt`(2), `kingsplitped.txt`(5) |
| `monomorphic__ibdseg_seglength10` | 0 | `king.seg`(14), `kingallsegs.txt`(2), `kingsplitped.txt`(5) |
| `monomorphic__related_degree2_ibdseg` | 0 | `king.kin`(15), `king.seg`(14), `kingallsegs.txt`(2), `kingsplitped.txt`(5) |
| `sexchr__ibdseg` | 0 | `king.seg`(14), `kingallsegs.txt`(3), `kingsplitped.txt`(5) |
| `sexchr__ibdseg_degree2` | 0 | `king.seg`(14), `kingX.seg`(14), `kingallsegs.txt`(3), `kingsplitped.txt`(5) |
| `sexchr__ibdseg_seglength5` | 0 | `king.seg`(14), `kingallsegs.txt`(3), `kingsplitped.txt`(5) |
| `sexchr__ibdseg_seglength10` | 0 | `king.seg`(14), `kingallsegs.txt`(3), `kingsplitped.txt`(5) |
| `sexchr__related_degree2_ibdseg` | 0 | `king.kin`(15), `king.seg`(14), `kingX.kin`(15), `kingX.seg`(14), `kingallsegs.txt`(3), `kingsplitped.txt`(5) |
| `unrelated__ibdseg` | 0 | `king.seg`(1), `kingallsegs.txt`(21), `kingsplitped.txt`(9) |
| `unrelated__ibdseg_degree2` | 0 | `king.seg`(0), `kingallsegs.txt`(21), `kingsplitped.txt`(9) |
| `unrelated__ibdseg_seglength5` | 0 | `king.seg`(1), `kingallsegs.txt`(21), `kingsplitped.txt`(9) |
| `unrelated__ibdseg_seglength10` | 0 | `king.seg`(1), `kingallsegs.txt`(21), `kingsplitped.txt`(9) |
| `unrelated__related_degree2_ibdseg` | 0 | `king.kin`(45), `king.seg`(0), `kingallsegs.txt`(21), `kingsplitped.txt`(9) |
| `admixed__ibdseg` | 0 | `king.seg`(16), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `admixed__ibdseg_degree2` | 0 | `king.seg`(15), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `admixed__ibdseg_seglength5` | 0 | `king.seg`(16), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `admixed__ibdseg_seglength10` | 0 | `king.seg`(16), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `admixed__related_degree2_ibdseg` | 0 | `king.kin`(18), `king.seg`(15), `kingallsegs.txt`(21), `kingsplitped.txt`(11) |
| `singleton__ibdseg` | 0 | `king.kin0`(0) |
| `singleton__ibdseg_degree2` | 0 | `king.kin0`(0) |
| `singleton__ibdseg_seglength5` | 0 | `king.kin0`(0) |
| `singleton__ibdseg_seglength10` | 0 | `king.kin0`(0) |
| `singleton__related_degree2_ibdseg` | 0 | `king.kin0`(0) |
| `pair__ibdseg` | 0 | `king.kin0`(1) |
| `pair__ibdseg_degree2` | 0 | `king.kin0`(1) |
| `pair__ibdseg_seglength5` | 0 | `king.kin0`(1) |
| `pair__ibdseg_seglength10` | 0 | `king.kin0`(1) |
| `pair__related_degree2_ibdseg` | 0 | `king.kin0`(1) |
| `bigish__ibdseg` | 0 | `king.seg`(763), `kingallsegs.txt`(22), `kingsplitped.txt`(188) |
| `bigish__ibdseg_degree2` | 0 | `king.seg`(442), `kingallsegs.txt`(22), `kingsplitped.txt`(188) |
| `bigish__ibdseg_seglength5` | 0 | `king.seg`(763), `kingallsegs.txt`(22), `kingsplitped.txt`(188) |
| `bigish__ibdseg_seglength10` | 0 | `king.seg`(763), `kingallsegs.txt`(22), `kingsplitped.txt`(188) |
| `bigish__related_degree2_ibdseg` | 0 | `king.kin`(573), `king.kin0`(26), `king.seg`(442), `kingallsegs.txt`(22), `kingsplitped.txt`(188) |
