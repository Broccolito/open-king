# Golden corpus — parameter flag group

Reference binary: **KING 2.3.2** (macOS arm64), black-box observation only.
Capture harness: `tests/parity/capture_params.py`
Alternate `--fam` / `--bim` inputs: `tests/parity/make_alt_inputs.py` (deterministic,
pure function of the corpus files; regenerate into `tests/parity/work/alt/`).

Flag group covered: `--prefix`, `--cpus`, `--minConc`, `--sexchr`, `--fam`, `--bim`,
plus `--degree` and bare `--kinship` / `--duplicate` baselines so each parameter's
effect can be isolated by diffing against its baseline.

**220 runs**, 13 datasets. 212 exited 0, 8 exited 1. Every run completed in under
0.15 s, so **nothing was skipped for time** — `bigish` (200 samples × 50 000 SNPs,
19 327 cross-family pairs) runs `--kinship` in ~60 ms.

Reproduce:

```sh
python3 tests/parity/generate_corpus.py  --outdir  tests/parity/work/data
python3 tests/parity/make_alt_inputs.py  --datadir tests/parity/work/data \
                                         --outdir  tests/parity/work/alt
python3 tests/parity/capture_params.py   --datadir tests/parity/work/data \
                                         --altdir  tests/parity/work/alt \
                                         --outroot tests/parity/golden/params
```

Each `<dataset>__<flagslug>/` holds `cmd.txt` (replayable argv, real paths rewritten
to `{DATA}` / `{ALT}`), `stdout.txt`, `stderr.txt`, `exitcode.txt`, and every file the
run produced. **Read `NONDETERMINISTIC.txt` before diffing anything** — one data file
in this corpus is race-corrupted by the reference binary itself.

---

## 1. Does `--cpus` change numeric output?

**For autosomes: no. For the X chromosome: yes, catastrophically — but it is a bug,
not a summation-order difference.**

`king.kin`, `king.kin0` and `kingX.kin` are byte-identical across `--cpus 1/2/4/8/16`
and the unflagged default, on all 13 datasets. Hashing `bigish`'s `king.kin0`
(19 327 pairs) over 10 runs at default threading gave 10/10 identical. **We do not
need to match KING's threading decomposition for autosomal results**; summation order
does not leak into the output.

`kingX.kin0` is different. KING appends to it from several threads with no lock:

| setting | 20 runs, distinct outputs | matches the `--cpus 1` result |
|---|---|---|
| `--cpus 1` | 1 | 20/20 |
| `--cpus 2` | 20 | 1/20 |
| `--cpus 3` | 20 | 0/20 |
| `--cpus 4` | 19 | 0/20 |
| `--cpus 8` | 20 | 0/20 |
| `--cpus 16` | 16 | 4/20 |

Records interleave mid-field; file size wandered 138–664 bytes against a correct 662.
Real captured rows from `sexchr__kinship_cpus4/kingX.kin0`:

```
SEX	S_DA	SU4	S_UF	MF	150000	0.348	0.0633	0.0105     <- ID truncated, N_SNP fused
SEX	S_DAU2	SU3	S_UM	FM	15	0.365	0	0.1727	-0.0351   <- field spliced in
SEX	S_DAU2	SU4	S_UF	FF	1500.365	0.1627	0.0547     <- fields merged
```

**Action:** pin every X-chromosome golden to `--cpus 1`. The `sexchr__kinship_cpus1*`
dirs are authoritative; the multi-thread `sexchr` dirs are retained only as evidence
and must be excluded from byte-diff. Do not reproduce the race.

`--cpus` does change stdout: the echoed value, the `<N> CPU cores are used.` line, and
the progress-percentage run (which scales *inversely* with thread count — 28 tokens at
`--cpus 1`, 3 at `--cpus 4`). `--cpus 0` is accepted silently and falls back to the
host default.

---

## 2. Which output columns are affected by `--degree`?

**None. `--degree` changes no column, no header, and no value — it is a pure row
filter on `king.kin0` only, plus one file-level side effect.**

* `king.kin` (within-family) is byte-identical at every degree — the filter never
  touches it.
* Retained `king.kin0` rows are **verbatim byte-identical** to the corresponding
  baseline rows (0 rows differ), so nothing is recomputed or re-rounded.
* Row order is preserved (still 32×32 tiled, see §4).

The predicate is exactly **`Kinship > 2^-(degree + 1.5)`**, strict inequality:

| degree | threshold | |
|---|---|---|
| 1 | 0.176777 | |
| 2 | 0.088388 | |
| 3 | 0.044194 | |
| 4 | 0.022097 | |

Verified across 8 datasets × 4 degrees: **0 missing, 0 extra** vs. filtering the
baseline `.kin0` by that predicate. No prescreen loss on this corpus (relevant
because KING advertises a `--noscreen` escape hatch).

Row counts, `bigish`: 19 327 → 3 (deg 1) / 24 (deg 2) / 59 (deg 3).
`multifam`: 150 → 8 / 32 / 52 / 63 (deg 4).

**File-level side effect worth pinning:** `--degree` **suppresses X-chromosome
analysis entirely.** `sexchr --kinship` emits `kingX.kin` + `kingX.kin0`;
`sexchr --kinship --degree 2` emits neither, and the `X-chromosome analysis...`
section vanishes from stdout.

---

## 3. Does `--prefix` change all filenames?

**Yes — all of them, with no exceptions.** The apparent `kingbySample.txt` vs
`king.kin` asymmetry is not some files ignoring the prefix; the default prefix is
literally the string `king`, and every output name is **`prefix + fixed_suffix`** by
raw string concatenation. The suffixes just belong to three shapes:

| shape | suffix begins with | default | `--prefix ZZZ` |
|---|---|---|---|
| dot-extension | `.` | `king.kin`, `king.kin0`, `king.con`, `king.ibs`, `king.ibs0` | `ZZZ.kin`, `ZZZ.kin0`, `ZZZ.con`, `ZZZ.ibs`, `ZZZ.ibs0` |
| bare concatenation | *(letter)* | `kingbySample.txt`, `kingbySNP.txt`, `kingallsegs.txt`, `kingunrelated.txt`, `kingunrelated_toberemoved.txt`, `kingbuild.log`, `kingupdateparents.txt` | `ZZZbySample.txt`, `ZZZbySNP.txt`, `ZZZallsegs.txt`, … |
| infix | `X.` | `kingX.kin`, `kingX.kin0` | `ZZZX.kin`, `ZZZX.kin0` |

Concatenation is verbatim, with no normalization — proven by the degenerate cases:

* `--prefix custom.` → `custom..kin` (double dot)
* `--prefix cus.tom` → `cus.tom.kin`
* `--prefix sub/pre` where `sub/` does not exist → **exit 1**, and the failure happens
  *before the .bim is even read*, with `FATAL ERROR - Cannot open sub/pre$TMP$.ped to
  write.` KING probes writability at startup by creating `<prefix>$TMP$.ped`. That
  temp file is cleaned up on success (no run dir contains one), but our implementation
  must reproduce the early-exit ordering to match error-path stdout.

---

## 4. What is the row ORDER of each output file?

Two different rules — and neither is `.fam` order.

**`king.kin`, `kingX.kin`, `king.con` — global order, all pairs `i<j`:**
families in `.fam` first-appearance order; within a family, members sorted
**ASCII-ascending by IID**. This is *not* `.fam` line order. `sexchr`'s `SEX` family
is `S_F, S_M, S_SON1, S_SON2, S_DAU1, S_DAU2` in the file but emits as
`S_DAU1, S_DAU2, S_F, S_M, S_SON1, S_SON2`; `dups`' `POFAM` emits `PO_C, PO_P` though
`PO_P` is the earlier line. Verified on all datasets that produce these files.

**`king.kin0` — 32×32 TILED, not sorted:**
the sort key is **`(i/32, j/32, i, j)`** over `.fam` index, ascending, `i<j`. KING
walks 32-sample blocks: outer loop over block-row `bi`, inner over block-column
`bj >= bi`, emitting `(i,j)` ascending within each tile.

This is invisible below 33 samples (one tile ⇒ looks plainly sorted), which is why it
only shows on `bigish` and `admixed`. Proof: `bigish` (200 samples ⇒ 7 blocks ⇒ 28
tiles) has **exactly 21** order inversions against plain `(i,j)` sort — and
`sum(6-bi) for bi in 0..6 = 21` is precisely the tile-transition count. `admixed`
(40 samples ⇒ 2 blocks ⇒ 3 tiles) has exactly 1. Verified `(i/32, j/32, i, j)`
ascending on **every** dataset, and on every `--degree`-filtered file.

**Any implementation that emits `.kin0` in naive `i<j` order will be byte-wrong on any
input with more than 32 samples.** The tiling is CPU-count independent — it is the
serial blocking of the algorithm, not a thread decomposition.

---

## 5. Surprising behavior

**KING validates only the `.bed` byte length — never the `.fam`/`.bim` record counts.**
The check is effectively `filesize >= 3 + M*ceil(N/4)`; excess bytes are ignored.

| probe | result |
|---|---|
| `--fam` with **one fewer** sample | **exit 0, silent.** `trio` "loaded: 2 samples", analysed 2. No warning. |
| `--bim` with **one fewer** SNP | **exit 0, silent.** `trio` "consist of 4999 autosome SNPs". |
| `--fam` with **one more** sample | exit 0 **when `ceil(N/4)` is unchanged** (`trio` 3→4, `pair` 2→3, `singleton` 1→2 all silently analyse a sample made of `.bed` pad bits); exit 1 only when the extra sample pushes the byte width up (`multifam` 20→21): `FATAL ERROR - Not enough genotypes at the 12500th marker` |
| `--bim` with **one more** SNP | exit 1 always: `FATAL ERROR - Not enough genotypes at the 5000th marker` |

Silent truncation and silent pad-bit-as-genotype are the dangerous cases; decide
deliberately whether to replicate them or diverge.

**`--sexchr N` shifts a four-chromosome window, and gates on X-SNP count.**
`N` designates X; then Y=`N+1`, XY=`N+2`, MT=`N+3`; everything else becomes autosome.
Confirmed on `sexchr` (chr 1,2,23,24,25,26 = 4000/1500/300/150/50 SNPs):
`--sexchr 23` ⇒ 4150 autosome (incl. 150 XY) / 1500 X / 300 Y / 50 MT;
`--sexchr 24` ⇒ 5550 autosome (incl. 50 XY) / 300 X / 150 Y;
`--sexchr 25` ⇒ 5800 autosome / 150 X / 50 Y.

* `--sexchr 0` and `--sexchr 1` ⇒ exit 1, `FATAL ERROR - Sex chromosome N out of range.`
* Any `N != 23` prints `Non-human samples are analyzed, with N pairs of chromosomes`.
* `--sexchr 23` is byte-identical to omitting the flag (23 is the default).
* **X analysis runs iff the designated X chromosome carries ≥ 512 SNPs.** Bisected
  exactly: 511 SNPs ⇒ no `kingX.*`; 512 ⇒ both files. (So `--sexchr 24`/`25` above
  produce no X output despite having 300/150 X SNPs — it is the count, not
  "non-human" mode: `--sexchr 2` puts 2000 SNPs on X and *does* emit `kingX.*`.)

**`--duplicate` has two code paths, split at N = 100, differing only when the result
is empty.** With no duplicates found, N < 100 writes a **header-only `king.con`**;
N ≥ 100 writes **no file at all** — same stdout (`No duplicates are found with
heterozygote concordance rate > 80%.`), same exit 0. Bisected exactly on synthetic
filesets: N=99 ⇒ file, N=100 ⇒ no file. `bigish` (200) is the only corpus dataset on
the far side, which is why `bigish__duplicate` is one of two exit-0 runs producing
nothing. When a duplicate *does* exist both paths emit identical rows (verified at
N = 50/99/100/120/200), so this is purely an empty-result artifact.

**`--minConc X` thresholds the `HetConc` column strictly (`> X`)**, and is echoed into
the message as a percentage. `dups` (45 pairs): `--minConc 0` ⇒ all 45 rows;
`0.5`/`0.9`/`0.99` ⇒ 2 rows (DUP HetConc 1.00000, MZ 0.99512); `--minConc 1` ⇒ 0 rows
and `No duplicates are found with heterozygote concordance rate > 100%.`

**Empty-file conventions are inconsistent between the autosomal and X writers.**
`king.kin` is created **0 bytes — no header at all** when the input has a single
family (`trio`, `nuclear`, `threegen`), but is **not created** when there are no
within-family pairs at all (`--fam` giving every sample its own FID). `kingX.kin`
under the same no-within-family-pairs condition **is** created, with a header.

**`--bim` proves chromosome/SNP-ID columns are inert for autosomal kinship.** The
`altbim` variant (every chromosome forced to 1, every SNP renamed) is byte-identical
to baseline on all 12 autosome-only datasets — only `sexchr` changes, where it demotes
X/Y/MT to autosomes and `kingX.*` disappears. Autosomal kinship depends solely on the
genotype matrix and the chromosome-class partition.

**`--fam` overrides are fully honoured** (FID, parents, and sex all take effect):
`altfam` moves every pair out of `king.kin` into `king.kin0`, and the file `king.kin`
is not created.

**Open question for the sexchr/X owner** (outside this flag group, flagged not
resolved): the `.fam` **sex** column materially changes X output, and one case looks
wrong. Isolating the two edits in `altfam` on `sexchr`, all at `--cpus 1`:

| `--fam` variant | `kingX.kin` rows | `kingX.kin0` rows |
|---|---|---|
| baseline `.fam` | 15 | 13 |
| FIDs changed only (sex preserved) | 0 | 22 |
| **sex flipped only** (FIDs preserved) | 0 | **0** |
| parents zeroed only | 15 | 13 |

Flipping sex silently empties all X output — consistent with KING dropping samples
whose X heterozygosity contradicts their declared sex (the corpus generates true males
haploid on X, so flipping makes them heterozygous "males"). Separately, the
FID-only variant retains 22 of 28 possible pairs with a `Sex` tally of FF=6, FM/MF=16,
**MM=0**, i.e. male–male pairs absent — yet the pedigree-grouped baseline *does* emit
`MM` rows. Those two facts are not obviously reconcilable and want a dedicated look.

---

## 6. Run inventory

`exit` is the process exit code; row counts exclude the header line.
Files marked `(0 B, EMPTY)` have no header either.

| dataset | run dir | flags | exit | files produced (data rows, header excluded) |
|---|---|---|---|---|
| admixed | `admixed__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| admixed | `admixed__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| admixed | `admixed__kinship` | `--kinship` | 0 | king.kin (18 rows), king.kin0 (762 rows) |
| admixed | `admixed__kinship_altbim` | `--kinship --bim {ALT}/admixed.altbim.bim` | 0 | king.kin (18 rows), king.kin0 (762 rows) |
| admixed | `admixed__kinship_altfam` | `--kinship --fam {ALT}/admixed.altfam.fam` | 0 | king.kin0 (780 rows) |
| admixed | `admixed__kinship_altfam_altbim` | `--kinship --fam {ALT}/admixed.altfam.fam --bim {ALT}/admixed.altbim.bim` | 0 | king.kin0 (780 rows) |
| admixed | `admixed__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (18 rows), king.kin0 (762 rows) |
| admixed | `admixed__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (18 rows), king.kin0 (762 rows) |
| admixed | `admixed__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (18 rows), king.kin0 (0 rows) |
| admixed | `admixed__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (18 rows), king.kin0 (0 rows) |
| admixed | `admixed__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (18 rows), king.kin0 (0 rows) |
| admixed | `admixed__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (18 rows), custom.kin0 (762 rows) |
| admixed | `admixed__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (18 rows), king.kin0 (762 rows) |
| bigish | `bigish__duplicate` | `--duplicate` | 0 | **none** |
| bigish | `bigish__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | **none** |
| bigish | `bigish__kinship` | `--kinship` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| bigish | `bigish__kinship_altbim` | `--kinship --bim {ALT}/bigish.altbim.bim` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| bigish | `bigish__kinship_altfam` | `--kinship --fam {ALT}/bigish.altfam.fam` | 0 | king.kin0 (19900 rows) |
| bigish | `bigish__kinship_altfam_altbim` | `--kinship --fam {ALT}/bigish.altfam.fam --bim {ALT}/bigish.altbim.bim` | 0 | king.kin0 (19900 rows) |
| bigish | `bigish__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| bigish | `bigish__kinship_cpus16` | `--kinship --cpus 16` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| bigish | `bigish__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| bigish | `bigish__kinship_cpus8` | `--kinship --cpus 8` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| bigish | `bigish__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (573 rows), king.kin0 (3 rows) |
| bigish | `bigish__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (573 rows), king.kin0 (24 rows) |
| bigish | `bigish__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (573 rows), king.kin0 (59 rows) |
| bigish | `bigish__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (573 rows), custom.kin0 (19327 rows) |
| bigish | `bigish__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (573 rows), king.kin0 (19327 rows) |
| dups | `dups__duplicate` | `--duplicate` | 0 | king.con (2 rows) |
| dups | `dups__duplicate_cpus1` | `--duplicate --cpus 1` | 0 | king.con (2 rows) |
| dups | `dups__duplicate_cpus4` | `--duplicate --cpus 4` | 0 | king.con (2 rows) |
| dups | `dups__duplicate_minConc0` | `--duplicate --minConc 0` | 0 | king.con (45 rows) |
| dups | `dups__duplicate_minConc0.5` | `--duplicate --minConc 0.5` | 0 | king.con (2 rows) |
| dups | `dups__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (2 rows) |
| dups | `dups__duplicate_minConc0.99` | `--duplicate --minConc 0.99` | 0 | king.con (2 rows) |
| dups | `dups__duplicate_minConc1` | `--duplicate --minConc 1` | 0 | king.con (0 rows) |
| dups | `dups__kinship` | `--kinship` | 0 | king.kin (2 rows), king.kin0 (43 rows) |
| dups | `dups__kinship_altbim` | `--kinship --bim {ALT}/dups.altbim.bim` | 0 | king.kin (2 rows), king.kin0 (43 rows) |
| dups | `dups__kinship_altfam` | `--kinship --fam {ALT}/dups.altfam.fam` | 0 | king.kin0 (45 rows) |
| dups | `dups__kinship_altfam_altbim` | `--kinship --fam {ALT}/dups.altfam.fam --bim {ALT}/dups.altbim.bim` | 0 | king.kin0 (45 rows) |
| dups | `dups__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (2 rows), king.kin0 (43 rows) |
| dups | `dups__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (2 rows), king.kin0 (43 rows) |
| dups | `dups__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (2 rows), king.kin0 (1 row) |
| dups | `dups__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (2 rows), king.kin0 (1 row) |
| dups | `dups__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (2 rows), king.kin0 (1 row) |
| dups | `dups__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (2 rows), custom.kin0 (43 rows) |
| dups | `dups__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (2 rows), king.kin0 (43 rows) |
| missing | `missing__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| missing | `missing__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| missing | `missing__kinship` | `--kinship` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_altbim` | `--kinship --bim {ALT}/missing.altbim.bim` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_altfam` | `--kinship --fam {ALT}/missing.altfam.fam` | 0 | king.kin0 (15 rows) |
| missing | `missing__kinship_altfam_altbim` | `--kinship --fam {ALT}/missing.altfam.fam --bim {ALT}/missing.altbim.bim` | 0 | king.kin0 (15 rows) |
| missing | `missing__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (0 B, EMPTY) |
| missing | `missing__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (0 B, EMPTY) |
| missing | `missing__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (0 B, EMPTY) |
| monomorphic | `monomorphic__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| monomorphic | `monomorphic__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| monomorphic | `monomorphic__kinship` | `--kinship` | 0 | king.kin (15 rows), king.kin0 (51 rows) |
| monomorphic | `monomorphic__kinship_altbim` | `--kinship --bim {ALT}/monomorphic.altbim.bim` | 0 | king.kin (15 rows), king.kin0 (51 rows) |
| monomorphic | `monomorphic__kinship_altfam` | `--kinship --fam {ALT}/monomorphic.altfam.fam` | 0 | king.kin0 (66 rows) |
| monomorphic | `monomorphic__kinship_altfam_altbim` | `--kinship --fam {ALT}/monomorphic.altfam.fam --bim {ALT}/monomorphic.altbim.bim` | 0 | king.kin0 (66 rows) |
| monomorphic | `monomorphic__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (15 rows), king.kin0 (51 rows) |
| monomorphic | `monomorphic__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (15 rows), king.kin0 (51 rows) |
| monomorphic | `monomorphic__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| monomorphic | `monomorphic__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| monomorphic | `monomorphic__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| monomorphic | `monomorphic__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (15 rows), custom.kin0 (51 rows) |
| monomorphic | `monomorphic__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (15 rows), king.kin0 (51 rows) |
| multifam | `multifam__bysample` | `--bysample` | 0 | kingallsegs.txt (18 rows), kingbySample.txt (20 rows) |
| multifam | `multifam__bysample_prefix_custom` | `--bysample --prefix custom` | 0 | customallsegs.txt (18 rows), custombySample.txt (20 rows) |
| multifam | `multifam__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| multifam | `multifam__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| multifam | `multifam__kinship` | `--kinship` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| multifam | `multifam__kinship_altbim` | `--kinship --bim {ALT}/multifam.altbim.bim` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| multifam | `multifam__kinship_altfam` | `--kinship --fam {ALT}/multifam.altfam.fam` | 0 | king.kin0 (190 rows) |
| multifam | `multifam__kinship_altfam_altbim` | `--kinship --fam {ALT}/multifam.altfam.fam --bim {ALT}/multifam.altbim.bim` | 0 | king.kin0 (190 rows) |
| multifam | `multifam__kinship_badbim` | `--kinship --bim {ALT}/multifam.badbim.bim` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| multifam | `multifam__kinship_badfam` | `--kinship --fam {ALT}/multifam.badfam.fam` | 0 | king.kin0 (171 rows) |
| multifam | `multifam__kinship_bigbim` | `--kinship --bim {ALT}/multifam.bigbim.bim` | 1 | **none** |
| multifam | `multifam__kinship_bigfam` | `--kinship --fam {ALT}/multifam.bigfam.fam` | 1 | **none** |
| multifam | `multifam__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| multifam | `multifam__kinship_cpus2` | `--kinship --cpus 2` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| multifam | `multifam__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| multifam | `multifam__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (40 rows), king.kin0 (8 rows) |
| multifam | `multifam__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (40 rows), king.kin0 (32 rows) |
| multifam | `multifam__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (40 rows), king.kin0 (52 rows) |
| multifam | `multifam__kinship_degree4` | `--kinship --degree 4` | 0 | king.kin (40 rows), king.kin0 (63 rows) |
| multifam | `multifam__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (40 rows), custom.kin0 (150 rows) |
| multifam | `multifam__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (40 rows), king.kin0 (150 rows) |
| nuclear | `nuclear__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| nuclear | `nuclear__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| nuclear | `nuclear__kinship` | `--kinship` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_altbim` | `--kinship --bim {ALT}/nuclear.altbim.bim` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_altfam` | `--kinship --fam {ALT}/nuclear.altfam.fam` | 0 | king.kin0 (15 rows) |
| nuclear | `nuclear__kinship_altfam_altbim` | `--kinship --fam {ALT}/nuclear.altfam.fam --bim {ALT}/nuclear.altbim.bim` | 0 | king.kin0 (15 rows) |
| nuclear | `nuclear__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (0 B, EMPTY) |
| nuclear | `nuclear__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (0 B, EMPTY) |
| pair | `pair__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| pair | `pair__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| pair | `pair__kinship` | `--kinship` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_altbim` | `--kinship --bim {ALT}/pair.altbim.bim` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_altfam` | `--kinship --fam {ALT}/pair.altfam.fam` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_altfam_altbim` | `--kinship --fam {ALT}/pair.altfam.fam --bim {ALT}/pair.altbim.bim` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_badfam` | `--kinship --fam {ALT}/pair.badfam.fam` | 0 | king.kin0 (0 rows) |
| pair | `pair__kinship_bigfam` | `--kinship --fam {ALT}/pair.bigfam.fam` | 0 | king.kin0 (3 rows) |
| pair | `pair__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin0 (1 row) |
| pair | `pair__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin0 (1 row) |
| pair | `pair__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin0 (1 row) |
| sexchr | `sexchr__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| sexchr | `sexchr__duplicate_cpus1` | `--duplicate --cpus 1` | 0 | king.con (0 rows) |
| sexchr | `sexchr__duplicate_cpus1_minConc0.9` | `--duplicate --cpus 1 --minConc 0.9` | 0 | king.con (0 rows) |
| sexchr | `sexchr__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| sexchr | `sexchr__kinship` | `--kinship` | 0 | king.kin (15 rows), king.kin0 (30 rows), kingX.kin (15 rows), kingX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_altbim` | `--kinship --bim {ALT}/sexchr.altbim.bim` | 0 | king.kin (15 rows), king.kin0 (30 rows) |
| sexchr | `sexchr__kinship_altfam` | `--kinship --fam {ALT}/sexchr.altfam.fam` | 0 | king.kin0 (45 rows), kingX.kin (0 rows), kingX.kin0 (0 rows) |
| sexchr | `sexchr__kinship_altfam_altbim` | `--kinship --fam {ALT}/sexchr.altfam.fam --bim {ALT}/sexchr.altbim.bim` | 0 | king.kin0 (45 rows) |
| sexchr | `sexchr__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (15 rows), king.kin0 (30 rows), kingX.kin (15 rows), kingX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_cpus1_altbim` | `--kinship --cpus 1 --bim {ALT}/sexchr.altbim.bim` | 0 | king.kin (15 rows), king.kin0 (30 rows) |
| sexchr | `sexchr__kinship_cpus1_altfam` | `--kinship --cpus 1 --fam {ALT}/sexchr.altfam.fam` | 0 | king.kin0 (45 rows), kingX.kin (0 rows), kingX.kin0 (0 rows) |
| sexchr | `sexchr__kinship_cpus1_degree1` | `--kinship --cpus 1 --degree 1` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| sexchr | `sexchr__kinship_cpus1_degree2` | `--kinship --cpus 1 --degree 2` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| sexchr | `sexchr__kinship_cpus1_degree3` | `--kinship --cpus 1 --degree 3` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| sexchr | `sexchr__kinship_cpus1_prefix_custom` | `--kinship --cpus 1 --prefix custom` | 0 | custom.kin (15 rows), custom.kin0 (30 rows), customX.kin (15 rows), customX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_cpus1_sexchr2` | `--kinship --cpus 1 --sexchr 2` | 0 | king.kin (15 rows), king.kin0 (30 rows), kingX.kin (15 rows), kingX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_cpus1_sexchr23` | `--kinship --cpus 1 --sexchr 23` | 0 | king.kin (15 rows), king.kin0 (30 rows), kingX.kin (15 rows), kingX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_cpus1_sexchr24` | `--kinship --cpus 1 --sexchr 24` | 0 | king.kin (15 rows), king.kin0 (30 rows) |
| sexchr | `sexchr__kinship_cpus1_sexchr25` | `--kinship --cpus 1 --sexchr 25` | 0 | king.kin (15 rows), king.kin0 (30 rows) |
| sexchr | `sexchr__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (15 rows), king.kin0 (30 rows), kingX.kin (15 rows), kingX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| sexchr | `sexchr__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| sexchr | `sexchr__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (15 rows), king.kin0 (0 rows) |
| sexchr | `sexchr__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (15 rows), custom.kin0 (30 rows), customX.kin (15 rows), customX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_sexchr1` | `--kinship --sexchr 1` | 1 | **none** |
| sexchr | `sexchr__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (15 rows), king.kin0 (30 rows), kingX.kin (15 rows), kingX.kin0 (13 rows) |
| sexchr | `sexchr__kinship_sexchr24` | `--kinship --sexchr 24` | 0 | king.kin (15 rows), king.kin0 (30 rows) |
| sexchr | `sexchr__kinship_sexchr25` | `--kinship --sexchr 25` | 0 | king.kin (15 rows), king.kin0 (30 rows) |
| singleton | `singleton__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| singleton | `singleton__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| singleton | `singleton__kinship` | `--kinship` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_altbim` | `--kinship --bim {ALT}/singleton.altbim.bim` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_altfam` | `--kinship --fam {ALT}/singleton.altfam.fam` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_altfam_altbim` | `--kinship --fam {ALT}/singleton.altfam.fam --bim {ALT}/singleton.altbim.bim` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_badbim` | `--kinship --bim {ALT}/singleton.badbim.bim` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_bigbim` | `--kinship --bim {ALT}/singleton.bigbim.bim` | 1 | **none** |
| singleton | `singleton__kinship_bigfam` | `--kinship --fam {ALT}/singleton.bigfam.fam` | 0 | king.kin0 (1 row) |
| singleton | `singleton__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin0 (0 rows) |
| singleton | `singleton__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin0 (0 rows) |
| singleton | `singleton__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin0 (0 rows) |
| threegen | `threegen__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| threegen | `threegen__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| threegen | `threegen__kinship` | `--kinship` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_altbim` | `--kinship --bim {ALT}/threegen.altbim.bim` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_altfam` | `--kinship --fam {ALT}/threegen.altfam.fam` | 0 | king.kin0 (66 rows) |
| threegen | `threegen__kinship_altfam_altbim` | `--kinship --fam {ALT}/threegen.altfam.fam --bim {ALT}/threegen.altbim.bim` | 0 | king.kin0 (66 rows) |
| threegen | `threegen__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (0 B, EMPTY) |
| threegen | `threegen__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__bysample` | `--bysample` | 0 | kingallsegs.txt (3 rows), kingbySample.txt (3 rows) |
| trio | `trio__bysample_prefix_custom` | `--bysample --prefix custom` | 0 | customallsegs.txt (3 rows), custombySample.txt (3 rows) |
| trio | `trio__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| trio | `trio__duplicate_minConc0` | `--duplicate --minConc 0` | 0 | king.con (3 rows) |
| trio | `trio__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| trio | `trio__duplicate_minConc1` | `--duplicate --minConc 1` | 0 | king.con (0 rows) |
| trio | `trio__kinship` | `--kinship` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_altbim` | `--kinship --bim {ALT}/trio.altbim.bim` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_altfam` | `--kinship --fam {ALT}/trio.altfam.fam` | 0 | king.kin0 (3 rows) |
| trio | `trio__kinship_altfam_altbim` | `--kinship --fam {ALT}/trio.altfam.fam --bim {ALT}/trio.altbim.bim` | 0 | king.kin0 (3 rows) |
| trio | `trio__kinship_badbim` | `--kinship --bim {ALT}/trio.badbim.bim` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_badfam` | `--kinship --fam {ALT}/trio.badfam.fam` | 0 | king.kin0 (1 row) |
| trio | `trio__kinship_bigbim` | `--kinship --bim {ALT}/trio.bigbim.bim` | 1 | **none** |
| trio | `trio__kinship_bigfam` | `--kinship --fam {ALT}/trio.bigfam.fam` | 0 | king.kin0 (6 rows) |
| trio | `trio__kinship_bimnotfound` | `--kinship --bim {ALT}/trio.no_such_file.bim` | 1 | **none** |
| trio | `trio__kinship_cpus0` | `--kinship --cpus 0` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (0 B, EMPTY) |
| trio | `trio__kinship_famnotfound` | `--kinship --fam {ALT}/trio.no_such_file.fam` | 1 | **none** |
| trio | `trio__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (0 B, EMPTY) |
| trio | `trio__kinship_prefix_dotted` | `--kinship --prefix cus.tom` | 0 | cus.tom.kin (0 B, EMPTY) |
| trio | `trio__kinship_prefix_subdir` | `--kinship --prefix sub/pre` | 1 | **none** |
| trio | `trio__kinship_prefix_traildot` | `--kinship --prefix custom.` | 0 | custom..kin (0 B, EMPTY) |
| trio | `trio__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (0 B, EMPTY) |
| unrelated | `unrelated__duplicate` | `--duplicate` | 0 | king.con (0 rows) |
| unrelated | `unrelated__duplicate_minConc0.9` | `--duplicate --minConc 0.9` | 0 | king.con (0 rows) |
| unrelated | `unrelated__kinship` | `--kinship` | 0 | king.kin (45 rows), king.kin0 (390 rows) |
| unrelated | `unrelated__kinship_altbim` | `--kinship --bim {ALT}/unrelated.altbim.bim` | 0 | king.kin (45 rows), king.kin0 (390 rows) |
| unrelated | `unrelated__kinship_altfam` | `--kinship --fam {ALT}/unrelated.altfam.fam` | 0 | king.kin0 (435 rows) |
| unrelated | `unrelated__kinship_altfam_altbim` | `--kinship --fam {ALT}/unrelated.altfam.fam --bim {ALT}/unrelated.altbim.bim` | 0 | king.kin0 (435 rows) |
| unrelated | `unrelated__kinship_cpus1` | `--kinship --cpus 1` | 0 | king.kin (45 rows), king.kin0 (390 rows) |
| unrelated | `unrelated__kinship_cpus2` | `--kinship --cpus 2` | 0 | king.kin (45 rows), king.kin0 (390 rows) |
| unrelated | `unrelated__kinship_cpus4` | `--kinship --cpus 4` | 0 | king.kin (45 rows), king.kin0 (390 rows) |
| unrelated | `unrelated__kinship_degree1` | `--kinship --degree 1` | 0 | king.kin (45 rows), king.kin0 (0 rows) |
| unrelated | `unrelated__kinship_degree2` | `--kinship --degree 2` | 0 | king.kin (45 rows), king.kin0 (0 rows) |
| unrelated | `unrelated__kinship_degree3` | `--kinship --degree 3` | 0 | king.kin (45 rows), king.kin0 (0 rows) |
| unrelated | `unrelated__kinship_prefix_custom` | `--kinship --prefix custom` | 0 | custom.kin (45 rows), custom.kin0 (390 rows) |
| unrelated | `unrelated__kinship_sexchr23` | `--kinship --sexchr 23` | 0 | king.kin (45 rows), king.kin0 (390 rows) |

Total runs: 220
