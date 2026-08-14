# Parity with KING 2.3.2 — the measured claim

This is open-king's authoritative statement of what it reproduces and what it does not.
Every number in it was produced by the commands in §1, against the reference binary `king`
2.3.2. Nothing here is an estimate, and nothing is rounded in our favour.

> **Headline: 403 of the 480 captured reference invocations are byte-identical (84.0 %).**
>
> **This number did not move in the final pass.** The `.seg` caller was replaced during it
> and every measure of *closeness* improved sharply (§4.4), but the residual is spread thin
> enough that no whole file crossed into byte-identity. 403/480 is the same count the
> previous caller scored, and the same 77 cases fail. That is stated here rather than buried
> because the interesting numbers below are all in §4, and none of them is a case count.
>
> **76 of the 77 that are not** trace to one cause: the `.seg` IBD-segment caller places a
> called IBD2 segment's endpoints a few markers from the reference, so the segment columns
> (`IBD1Seg`, `IBD2Seg`, `PropIBD`) are close but not equal on a minority of the rows that
> carry them — 157 of the 982 primary rows, by 0.00036 of the genome on average and never
> more than 0.0089. The *set* of pairs reported is exactly right everywhere: **0 extra and
> 0 missing rows** on every output file in the corpus, and so are all 982 `InfType` labels
> and every `Error`. **The remaining 1** is structural: `--build`'s pedigree reconstruction
> is unimplemented (§6.2) — and it is *also* segment-blocked, now that its one missing
> statistic has been identified and measured, so it is not an independent second problem.
> `<prefix>X.seg`, also unwritten (§6.1), costs 2 further cases that are already inside the
> 76 — both fail on the segment columns as well.
>
> Seven analyses are byte-identical on every dataset that runs them: `--kinship`,
> `--duplicate`, `--bysample`, `--bySNP`, `--autoQC` and **`--ibs`** at 13/13 each and
> `--unrelated` at 26/26, plus the whole 220-case `params` group. `<prefix>allsegs.txt` —
> the per-segment listing that underlies everything above — is byte-identical in all 163
> cases that produce it.
>
> `--ibs` joined that list when `Scan::ibd2_words` was replaced by the **chunk scan** of
> `docs/research/16-segment-extension.md`: its two IBD2 columns, `MaxIBD2` and `Pr_IBD2`,
> are now exact on all 21 561 `.ibs`/`.ibs0` rows. That closed the `--ibs` IBD2 caller and
> nothing else — `--ibdseg`'s `.seg` caller is a *different* caller (§5.8), and it took its
> own campaign and its own canvas (`docs/research/17-seg-caller.md`) to pin down: every
> constant of `Scan::ibd2` is now a bisection off the reference, `IBD2Seg` is exact on
> **896 of 982** rows against 822, and — scoring both callers the same way, with
> `tests/parity/fit/seg17.py` — the mean `PropIBD` error is **0.00037** against 0.00138 and
> the worst row **0.0089** against 0.2109. **No parity case flipped**, which is exactly
> why §4 and not §3 is where that work is graded.
>
> **Two divergences live outside the corpus entirely** and therefore cost no case here, but
> would be visible to a user: `--ibdseg` does not apply the reference's 100 Mb usable-total
> floor, and `<prefix>splitped.txt` is written unconditionally where the reference sometimes
> writes none. Both were measured against the reference in this pass and are new; §5.10.

---

## 1. Reproducing every number below

```bash
cd /path/to/open-king
cargo build --release

# pass/fail for all 480 cases  -> "403 PASS, 77 FAIL, 480 total"
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

# the Python engine mirror must still reproduce the binary's own output
cd tests/parity/fit && python3 check_mirror.py     # -> "MIRROR OK"
```

`run_parity.py` and `measure_gaps.py` are Python 3 standard library only, regenerate the
input corpus automatically on first run (~20 s) and need no reference binary. The probes
drive the reference directly. `run_parity.py` exits 0 when every case passed, 1 when at
least one failed, 2 on a harness error.

Measured on the tree this document describes:

| command | result |
| --- | --- |
| `run_parity.py --impl target/release/king` | **403 PASS, 77 FAIL, 480 total**, 874 output files byte-compared, 8 diff-excluded |
| `run_parity.py --impl <reference>` | **480 PASS, 0 FAIL**, 876 files byte-compared — the normalization is complete and the goldens are self-consistent |
| `probes/degree_filter.py --ref <reference>` | 38 298 cases, **0 false-keep, 0 false-drop** |
| `docs/research/fixtures/gate8.py` | brackets the `--degree 1` IBD2 clause to (0.0789, 0.0829] — its ladder refuses at `PropIBD` 0.0789 and accepts at 0.0829 |
| `tests/parity/fit/check_mirror.py` | **MIRROR OK** — the mirror reproduces the binary's `.seg` columns and every printed `MaxIBD2`, 13 datasets |
| `tests/parity/fit/seg17.py` | committed rule `exact 709  both 825  ibd2 896  MAE 0.00037  worst 0.0089`; retired geometry `705 / 820 / 822 / 0.00138 / 0.2109` |
| `cargo test --workspace` | **285 passed, 0 failed, 1 ignored** |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| `cargo build --release --offline` from a clean tree | succeeds in **8.1 s**; `Cargo.lock` has 15 packages — the 3 workspace crates and 12 external |
| that clean-tree binary, re-run through `run_parity.py` | **403 PASS, 77 FAIL** — the published count does not depend on a warm `target/` |

By capture group: `apps` **89/91**, `core` **74/104**, `ibdseg` **20/65**, `params`
**220/220**. By analysis: `--kinship`, `--duplicate`, `--bysample`, `--bySNP`, `--autoQC`,
`--ibs` **13/13** each; `--unrelated` **26/26**; `--cluster` and `--build` **12/13**;
`--related` **35/65**; `--ibdseg` **16/52**; `--related --ibdseg` **4/13**; `params`
**220/220**.

---

## 2. The matrix: every analysis against every dataset

One cell is every captured invocation of that analysis on that dataset. `--related` is five
cases per dataset (bare plus four `--degree`), `--ibdseg` four (bare, `--degree 2`,
`--seglength 5`, `--seglength 10`), `--related --ibdseg` one, `--unrelated` two.
**Bold** means every captured invocation is byte-identical: every output file, every column,
plus stdout, stderr and exit status.

| analysis | trio | nuclear | threegen | multifam | dups | missing | monomorphic | sexchr | unrelated | admixed | singleton | pair | bigish | total |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `--kinship` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--duplicate` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--bysample` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--bySNP` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--autoQC` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--unrelated` | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **2/2** | **26/26** |
| `--cluster` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 12/13 |
| `--build` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 12/13 |
| `--related` | **5/5** | **5/5** | **5/5** | 0/5 | 0/5 | **5/5** | 0/5 | 0/5 | **5/5** | 0/5 | **5/5** | **5/5** | 0/5 | 35/65 |
| `--ibdseg` | **4/4** | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | **4/4** | 0/4 | **4/4** | **4/4** | 0/4 | 16/52 |
| `--related --ibdseg` | **1/1** | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | **1/1** | 0/1 | **1/1** | **1/1** | 0/1 | 4/13 |
| `--ibs` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| flag plumbing + error probes (`params`) | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | **220/220** |
| | | | | | | | | | | | | | | **403/480** |

By capture group: `apps` **89/91**, `core` **74/104**, `ibdseg` **20/65**, `params`
**220/220**.

The `params` group is 220 invocations that exercise the command-line surface rather than one
dataset: every `--prefix` shape, `--cpus`, `--sexchr`, `--degree`, `--minConc`,
`--seglength`, alternate `--fam`/`--bim` inputs, malformed and missing files, and the banner
in each case. All 220 are byte-identical.

Where `--related` and `--ibdseg` pass, it is for one of two reasons, both real:

* **The reference declines to run the segment pass.** Below 10 samples `--related` emits the
  10-column `--kinship` form (`trio`, `nuclear`, `missing`, `singleton`, `pair`), and
  `trio`/`nuclear`/`threegen`/`missing` additionally get a **zero-byte** `king.kin` under
  `--related` even though stdout announces the file. `--ibdseg` takes the same downgrade
  below 5 samples, which is why `trio` (3), `singleton` (1) and `pair` (2) pass it.
* **The caller is exactly right on that dataset.** `unrelated` reports one pair and gets it
  exact, and `threegen` passes `--related` (empty `.kin`) while failing `--ibdseg`. `--ibs`
  is now in this second category on **all thirteen** datasets.

---

## 3. What is byte-identical, per output file

Counted over **every** case that produces the file, across all four groups. "Rows" is every
data row in every such case, so these percentages are corpus-wide row accuracy — not the
narrower denominator §4 uses.

| output file | cases | byte-identical cases | data rows | rows differing | rows exact |
| --- | ---: | ---: | ---: | ---: | ---: |
| `<prefix>allsegs.txt` | 163 (+2 custom) | **all** | 2 196 | 0 | **100 %** |
| `<prefix>splitped.txt` | 50 | **all** | — | 0 | **100 %** |
| `<prefix>.con` (`--duplicate`) | 46 | **all** | 62 | 0 | **100 %** |
| `<prefix>bySample.txt` | 15 (+2 custom) | **all** | — | 0 | **100 %** |
| `<prefix>bySNP.txt` | 13 | **all** | — | 0 | **100 %** |
| `<prefix>_autoQC_Summary.txt` | 13 | **all** | — | 0 | **100 %** |
| `<prefix>_autoQC_snptoberemoved.txt` | 13 | **all** | 23 725 | 0 | **100 %** |
| `<prefix>_autoQC_sampletoberemoved.txt` | 13 | **all** | 2 | 0 | **100 %** |
| `<prefix>_autoQC_updatesex.txt` | 1 | **all** | — | 0 | **100 %** |
| `<prefix>updateids.txt` | 2 | **all** | — | 0 | **100 %** |
| `<prefix>unrelated.txt` | 26 | **all** | — | 0 | **100 %** |
| `<prefix>unrelated_toberemoved.txt` | 26 | **all** | — | 0 | **100 %** |
| `<prefix>X.kin0` | 5 diffable of 13 (§5.2) | **all 5** | 52 | 0 | **100 %** |
| `<prefix>.ibs` | 13 | **all** | 807 | 0 | **100 %** |
| `<prefix>.ibs0` | 8 | **all** | 20 754 | 0 | **100 %** |
| `<prefix>.kin0` | 168 | 160 | 207 986 | 28 | **99.99 %** |
| `<prefix>.kin` | 187 | 151 | 12 081 | 810 | 93.30 % |
| `<prefix>X.kin` | 15 | 9 | 195 | 12 | 93.85 % |
| `<prefix>cluster.kin` | 1 | 0 | 165 | 30 | 81.82 % |
| `<prefix>.seg` | 50 | 5 | 4 172 | 1 252 | 69.99 % |
| `<prefix>build.log` | 8 | 7 | 19 | 19 | see §6.2 |
| `<prefix>updateparents.txt` | 8 | 7 | 34 | 34 | see §6.2 |
| `<prefix>X.seg` | 2 | 0 | — | — | **never written**, §6.1 |

Every `--kinship` case is byte-identical; the 36 differing `<prefix>.kin` cases are all
`--related`. Row identity is matched on the identifier columns before any comparison, so
across the whole corpus there are **0 extra and 0 missing rows** on every file.

**stdout, stderr and exit status.** 475 of the 480 cases match stdout byte-for-byte after
the normalization of §7. **5 cases differ on stdout**, and every one of them also differs
on a file — no case in the suite fails on console output alone:

| cases | stdout line that differs | cause |
| ---: | --- | --- |
| 2 | `bigish --related --degree 2` — `Stages 1&2 (with 32768 SNPs): 36 pairs` vs `50 pairs` | the two-stage screening bound, §5.7 |
| 2 | `sexchr --ibdseg --degree 2` — the `…X-Chr IBD segments saved in file kingX.seg` line is absent | §6.1 |
| 1 | `apps/bigish__build` | §6.2 |

The six `monomorphic --related*` cases used to sit at the head of that table, differing on
the `Inference` relationship-count row (`0 10 2 2 0 1` against `0 9 3 2 0 1`) because a
segment `InfType` downstream of the old IBD2 caller was wrong. They now match: those six
cases still fail, but on `king.kin` / `king.seg` numbers alone (§4.4).

---

## 4. The gaps, measured

Row counts in this section use `measure_gaps.py`'s denominator, which is **rows inside the
cases that differ**, not rows corpus-wide — it is the tighter, less flattering number. §3
gives the corpus-wide view of the same data.

### 4.1 The segment columns — 76 of the 77 failures

Every one of these is `IBD1Seg` / `IBD2Seg` / `PropIBD` (and the `InfType` / `Error` that
follow from them) being close but not equal.

| file | rows differing | of | +extra | −missing | column | mean abs err | worst |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| `king.seg` | 1 252 | 4 169 | **0** | **0** | `PropIBD` 1 229 rows | 0.00264 | 0.0679 |
| | | | | | `IBD1Seg` 856 rows | 0.00717 | 0.0793 |
| | | | | | `IBD2Seg` 465 rows | 0.00153 | 0.1031 |
| | | | | | `InfType` 1 row | — | — |
| `king.kin` (`--related`) | 810 | 3 978 | **0** | **0** | `IBD1Seg` 810 rows | 0.00399 | 0.0125 |
| | | | | | `PropIBD` 810 rows | 0.00217 | 0.0089 |
| | | | | | `IBD2Seg` 420 rows | 0.00050 | 0.0043 |
| `king.kin0` (`--related`) | 28 | 296 | **0** | **0** | `IBD1Seg` 28 rows | 0.00545 | 0.0130 |
| | | | | | `PropIBD` 28 rows | 0.00296 | 0.0074 |
| | | | | | `IBD2Seg` 28 rows | 0.00041 | 0.0018 |
| `kingX.kin` (`--related`) | 12 | 90 | **0** | **0** | `IBD1Seg` 12 rows | 0.00120 | 0.0018 |
| | | | | | `PropIBD` 12 rows | 0.00060 | 0.0009 |
| `kingcluster.kin` | 30 | 165 | **0** | **0** | `IBD1Seg` 30 rows | 0.00424 | 0.0094 |
| | | | | | `PropIBD` 30 rows | 0.00234 | 0.0059 |
| | | | | | `IBD2Seg` 15 rows | 0.00055 | 0.0021 |

Everything in that table moved when `Scan::ibd2` was replaced by the `.seg`-native caller of
`docs/research/17-seg-caller.md` (§4.4). Two columns left it altogether: **`Error` and
`InfType` now differ on no `.kin`, `.kin0` or `X.kin` row anywhere**, where the caller they
replaced missed 12 and 6, and `kingX.kin`'s `IBD2Seg` is exact. The `.seg` worst row fell
from 0.2109 to 0.0679 and mean `PropIBD` error from 0.00562 to 0.00264, over all 50 captured
`.seg` cases and all three `--seglength` floors.

`king.ibs` and `king.ibs0` are **no longer in this table**: since the chunk scan of §5.8
they are byte-identical in all 13 and all 8 cases. That removes the sharpest per-segment
grader from the toolbox as well as from the gap — `MaxIBD2`, the length in base pairs of a
*single* IBD2 segment, is now exact on all **21 561** `.ibs`/`.ibs0` rows and on all **149**
`.ibs` rows where the reference reports a non-zero one, and so has no residual left to point
at. Anything that grades the remaining gap has to be read off `.seg`, whose caller is a
different one (§5.8).

`kingcluster.kin` fails only on `bigish`, and only on its three segment columns: the pair
set, the ordering and the other eleven columns are exact. `--cluster` is `--related` re-run
inside the merged families, so this is the segment residual and nothing else — there is no
clustering rule left to find.

### 4.2 The 16-column `--related` layer is complete

Measured directly, over **4 805** comparable `.kin` / `.kin0` / `X.kin` rows (**3 925** of
which have both `IBD1Seg` and `IBD2Seg` byte-exact):

| column | rows differing where the segments are **exact** | rows differing where the segments already differ |
| --- | ---: | ---: |
| `IBD1Seg` | — | 880 |
| `IBD2Seg` | — | 463 |
| `PropIBD` | **0** | 880 |
| `InfType` | **0** | **0** |
| `Error` | **0** | **0** |

`InfType` and `Error` are now exact on **every one of the 4 805 rows**, differing nowhere at
all — under the retired caller they missed 6 and 12. And no row anywhere differs on `N_SNP`,
`Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0` or `Kinship`. So `HetConc`,
`HomIBS0`-as-union, the `InfType` table, the `Error` grader, row order, the `.kin0`
`N >= 100` gate, the `< 10` sample downgrade and the
`Kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)` inclusion disjunct are all exact, and
**100 % of the `--related` residual is the segment caller**. Every row on which `PropIBD`
differs is a row whose `IBD1Seg` already differs — the 880 are the same 880.

Two rules in this layer are not what an earlier draft of this document claimed, and both are
now pinned against purpose-built pedigree probes (`tests/parity/probes/pederr.py`):

* **`Error` is not a label-degree comparison.** It is `--kinship`'s own multiplicative
  grader fed the *segment* kinship, but only for the middle degrees:
  `InfType ∈ {2nd, 3rd, 4th} → kinship::error_flag(PropIBD / 2, Phi)`, otherwise `0` if
  `InfType` equals the pedigree's label and `1` if not. 0 mismatches over 6 069 rows
  (4 248 corpus + 1 821 probe) against 11 for the old rule.
* **`.kin`'s `InfType` is not `.seg`'s.** The `Dup/MZ` clause additionally requires
  `HetConc > 0.8`; the same pair at `IBD1Seg 0.0182 / IBD2Seg 0.8128` prints `Dup/MZ` in
  `.seg` and `FS` in `.kin`. Bracketed to (0.7986, 0.8004] on a 72-pair ladder; `--minConc`
  does not move it. `king_core::ibdseg::inf_type` stays ungated because `.seg` genuinely
  uses the plain `IBD2Seg > 0.7`; the gate lives in the `.kin` writers.

### 4.3 `.seg`'s `PropIBD` is not a formula we can fix

`PropIBD` differs on more `.seg` rows (1 229) than either of its own inputs (856, 465), which
looks like an arithmetic bug and is not one.

On the primary capture, **116 of the 825 rows whose `IBD1Seg` and `IBD2Seg` are both exact
still print a different `PropIBD`** — and every one of the 116 differs by **exactly ±1 in the
fourth decimal**: ours reads one ulp *low* on 102 of them and one ulp *high* on 14. This
count is essentially untouched by the §4.4 caller change (it was 115 of 820), which is the
point: it is not a segment problem. Fifteen distinct
formulations were scored over those rows — `ibd2 + ibd1/2`, `(ibd1 + 2·ibd2)/2`, `(b2 + b1/2)/d`, `(2·b2 + b1)/(2·d)`,
integer halving of the base pairs, `f32` intermediates, and others — and **all fifteen scored
exactly 705/820** at the time (against the caller §5.0 now lists as retired). No
reassociation and no precision change moved a single row.

What is actually happening is visible in the reference's own output. Run the reference once:

```bash
king -b bigish.bed --related --degree 2 --ibdseg --cpus 1 --prefix r
```

**147** pairs appear in both `r.kin` and `r.seg`; **all 147** carry identical `IBD1Seg` and
`IBD2Seg` in the two files, and **43 of them carry a different `PropIBD`** — for example
`IBD1Seg 0.3852 / IBD2Seg 0.3123` prints `PropIBD` **0.5048** in `.kin` and **0.5049** in
`.seg`. Same invocation, same pair, same printed inputs, two answers: the reference's two
writers do not agree with each other. open-king computes `PropIBD` once, from the unrounded
estimates, and that choice matches `.kin` exactly (§4.2, 0 of 3 925) at the cost of these
116 `.seg` rows.

This is recorded as a **known limitation, not a to-do**: an earlier handoff note claimed
~386 `.seg` rows were recoverable "for free" by reassociating the expression. That claim is
refuted by the fifteen-formulation sweep above.

### 4.4 The primary `--ibdseg` scorecard

`<dataset>__ibdseg`, the default 3 Mb reporting floor, 982 rows over 10 datasets:

* **825 of 982 rows** have both `IBD1Seg` and `IBD2Seg` exact at the printed four decimals;
  the `IBD2Seg` column alone is exact on **896** and `IBD1Seg` alone on **826**. This is
  what `cargo test -p king-core --test ibdseg_parity` prints; run it with
  `KING_GOLDEN=tests/parity/golden` for the per-dataset breakdown.
* **709 of 982 rows** are exact on all four printed columns; the 116-row gap is §4.3.
* Mean absolute `PropIBD` error **0.00036**, worst **0.0089**.
* **All 982 `InfType` labels** are right.
* **0 extra and 0 missing pairs**, on every dataset.

> The mean is printed column against printed column, which is what a user diffing two files
> sees. `tests/parity/fit/seg17.py` scores the *unrounded* estimate against the reference's
> printed value and so reports **0.00037** on the same 982 rows; the difference is our own
> rounding to four decimals, not a different rule. Comparisons against the retired geometry
> below are quoted on `seg17.py`'s scale, where both sides are measured the same way.

Held out at other reporting floors, rules unchanged: `--seglength 5` gives **701/982** exact
on all four columns (MAE 0.00066, worst 0.0552), `--seglength 10` gives **668/982**
(MAE 0.00158, worst 0.0679). Neither floor was used to fit anything.

Those numbers are the `.seg`-native caller of `docs/research/17-seg-caller.md`, committed to
`Scan::ibd2`. Against the geometry it replaced, both scored by `seg17.py` in the same run:
exact rows 705 → **709**, both-columns 820 → **825**, `IBD2Seg` 822 → **896**, mean error
0.00138 → **0.00037** (÷3.7), worst row 0.2109 → **0.0089** (÷24), and 0 extra / 0 missing
either way. `InfType` went 981 → **982** on this capture and, corpus-wide, from 6 wrong rows
to none (§4.2).

Per dataset, on that capture. "all four" is `IBD1Seg`, `IBD2Seg`, `PropIBD` and `InfType`
all byte-exact; "both est." is the two estimate columns only, which is the figure
`cargo test -p king-core --test ibdseg_parity` prints. Retired caller in brackets where it
differs:

| dataset | all four | both est. | of | mean abs `PropIBD` err | worst |
| --- | ---: | ---: | ---: | ---: | ---: |
| `bigish` | 557 | 649 | 763 | 0.00034 | 0.0059 |
| `multifam` | 73 | 87 | 104 | 0.00046 | 0.0073 |
| `threegen` | 28 | 36 | 39 | 0.00024 | 0.0058 |
| `admixed` | 12 | 13 | 16 | 0.00043 | 0.0046 |
| `dups` | 2 | 2 | 3 | 0.00297 | 0.0089 |
| `unrelated` | 1 | 1 | 1 | 0.00000 | 0.0000 |
| `sexchr` | 8 | 8 | 14 | 0.00071 | 0.0047 |
| `nuclear` | 8 | 8 | 14 | 0.00060 | 0.0045 |
| `missing` | 8 | 9 | 14 | 0.00042 | 0.0030 |
| `monomorphic` | 12 [8] | 12 [8] | 14 | 0.00003 | 0.0003 |
| **total** | **709** | **825** | **982** | **0.00036** | **0.0089** |

The last four are the small ones — each reports the 14 within-family pairs of a single
six-person nuclear family, over 5 000 to 10 000 markers against `bigish`'s 50 000. Read them
with §5.1 in hand.

**The residual is one-sided on `IBD1Seg` and two-sided on `IBD2Seg`**: of the **157** rows
whose `IBD1Seg` or `IBD2Seg` differs (982 − 825), `IBD1Seg` is too high on 156 and too low on
none, `IBD2Seg` too high on 52 and too low on 34. The retired caller was two-sided on both
(139/21 and 39/121), and its `IBD2Seg` bias — too low on 121 rows — is gone.

**Detection is finished; what is left is length.** Splitting the same 982 rows by whether the
reference reports any IBD2 at all (`docs/research/14-ibd2-geometry.md` §2 drew this table
against the retired caller; the second column is what moved):

| reference row | rows | both estimate columns exact |
| --- | ---: | ---: |
| `IBD2Seg == 0.0000` | 823 | **823** [819] |
| `IBD2Seg > 0` | 159 | **2** [1] |

The first line is now clean. Its four former exceptions were all `monomorphic`, where the
reference prints `IBD2Seg 0.0000` and the retired rule invented 0.42, 0.27, 0.03 and 0.08 of
the genome out of words carrying two to four mismatches (§5.1). So the IBD1 caller, its
boundary refinement and the "is there any IBD2 here at all" question are all finished, and
every failing `ibdseg/*__ibdseg` case now fails only on the *length* of the 159 rows that do
carry IBD2 — by 0.00036 of the genome on average and never by more than 0.0089.

**How to grade further work on it.** Not with the parity case count, and not with the
exact-row count, which the 823 IBD2-free rows dominate: a rule that cut mean error by a
factor of 3.7 and the worst row by a factor of 24 moved the case count by **zero** and the
exact-row count only 705 → 709. Not with `--ibs` either — `Pr_IBD2` and `MaxIBD2` have been
exact under every candidate since the chunk scan of §5.8. Grade with the **`IBD2Seg` column
and the mean** (`tests/parity/fit/seg17.py` prints both, plus the retired caller for
reference), and out of sample with the `.seg`-native canvas of
`docs/research/fixtures/segcanvas.py`, which inverts a printed `IBD2Seg` back to the number
of calls and the number of words and carries 872 cached reference answers so it re-runs in
seconds without the binary.

---

## 5. Known limitations

### 5.0 The segment residual: what is solved, what is not

Everything in §4 says the same thing from different angles, so here it is once, as a ledger.
Numbered 5.0 rather than 5.1 on purpose: §5.1…§5.9 are cross-referenced from the crates and
from `docs/research/`, so their numbers must not move.

**Solved, and not worth re-deriving** (each measured, with the experiment named):

| piece | status | where |
| --- | --- | --- |
| **The acceptance gate** — a run is called iff popcount over its own *complete* 64-marker words of `inf1 = p1_i & p1_j & (p0_i \| p0_j)` (IBD1) or `inf2 = p1_i & p1_j` (IBD2) is **≥ 10** | exact **and unique** | `13-informativeness-gate.md`, `tests/parity/fit/gate_*.py` |
| **Which pairs are reported** (`--degree` inclusion, the `.kin0` `N ≥ 100` gate, the `< 10` and `< 5` sample downgrades) | **0 extra, 0 missing rows on every output file in the corpus**; the degree filter itself **0 false-keep, 0 false-drop over 38 298 cases** | §3, §4.2, `probes/degree_filter.py`, `fixtures/gate8.py` |
| **The per-segment listing itself** — `allsegs.txt` | byte-identical in **all 163** cases (+2 custom) | §3 |
| **Denominators and thresholds** — `D` = sum over autosomal `allsegs.txt` rows; `--seglength` inclusive; nothing merges across a gap | exact at 3, 5 and 10 Mb | §4.4 |
| **The aggregate** — `PropIBD = IBD2Seg + IBD1Seg/2` in `f64` | 0 differences on all 3 925 rows whose segments are exact | §4.2 |
| **The entire non-segment surface** — `N_SNP`, `Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship` | **no row anywhere differs**, over 4 805 rows | §4.2 |
| **`InfType` and `Error`** | **no row anywhere differs**, over all 4 805 rows — not merely where the segments are exact | §4.2 |
| **The IBD1 caller and its boundary refinement** | **all 823** IBD2-free rows exact | §4.4 |
| **The `--ibs` IBD2 caller** (`Scan::ibd2_words`, the chunk scan) | exact on all **21 561** rows | §5.8 |
| **The `.seg` IBD2 caller's word predicate, gate, reach and push** (`Scan::ibd2`) | every constant bisected on a `.seg`-native canvas; `IBD2Seg` exact on **896 of 982** rows, all 982 `InfType` | `17-seg-caller.md`, `fixtures/segcanvas.py` |

**Not solved — two clauses, in one function.** `Scan::ibd2` is where the whole remaining
`.seg` residual lives, and two of its parts are fitted rather than bisected:

1. **The bridging condition** — "a lone unusable word carrying no opposite homozygote is
   absorbed iff the next word is mismatch-free and the usable words from there carry
   `inf2 >= 10`". Its *lookahead* is guessed: it stops at the next unusable word, and the
   canvases `Cyzy` and `zyzy` say it should not.
2. **The left end of a run that opens on a mismatch-carrying word** — a two-marker effect
   the corpus cannot see, and the larger of the two out of sample (9 of the 11 exhaustive
   misses, against the bridging clause's 2).

Everything else — the word predicate, the gate and where it starts counting, the 63-marker
reach, the whole-word IBS0 block, the one-word push, the segment fringes — is a bisection
off the reference (`docs/research/17-seg-caller.md` §3–§7).

**The precise shape of what is left** (so the next person can recognise it):

* **Detection is done; length is not.** Split the 982 primary rows by whether the reference
  reports any IBD2: `IBD2Seg == 0` → **all 823** exact; `IBD2Seg > 0` → **2 of 159**, off by
  0.00036 of the genome on average and never more than 0.0089.
* **`IBD1Seg` is now one-sided**: too high on 156 of the 157 imperfect rows, never too low.
  `IBD2Seg` is two-sided, 52 high against 34 low (§4.4).
* **It is measurable out of sample, and the residual splits two ways.** Against the
  reference the rule reproduces **329 of 340** exhaustive length-≤4 word-sequence canvases
  and **223 of 240** held-out random ones (three seeds, 93 %), where the geometry it
  replaced scores **4 of 340** and **42 of 240**. The 11 exhaustive misses are not all one
  clause, and the split is the useful part:

  | miss | what it is | count |
  | --- | --- | ---: |
  | `xyC`, `xxyC`, `xyCC`, `xyCz`, `xyCx`, `xyCy`, `xyzC`, `xyzx`, `yxyC` | the run **opens on a mismatch-carrying word** — the two-marker left-end effect | 9 |
  | `Cyzy`, `zyzy` | the **bridging lookahead**, which stops at the next unusable word | 2 |

  So the left-end effect, previously filed as the smaller of the two open items, is
  responsible for more of the exhaustive residual than the bridging clause is. Corpus-wise
  the bridging clause is worth `IBD2Seg` 896 against **874** with no bridging at all, but
  costs a little mean (0.00037 against 0.00033) and tail (0.0089 against 0.0071) — the
  signature of a clause that is *nearly* right. Both figures come from
  `python3 seg17.py` with `bridge="clean"` and `bridge="none"`.
* **The Rust is an exact port of the fitted model, so the residual is the model's.** Driving
  `segcanvas.py` with our own release binary instead of the reference, the binary agrees
  with `segcanvas.predict()` on **340 of 340** exhaustive canvases and **240 of 240** random
  ones. Everything above is therefore model-versus-reference, not a porting defect. (Do this
  in a *copy* of `docs/research/fixtures/` with an empty cache — see the warning below.)
* **The `--ibs` solution does not port, and that is measured, not assumed.** `.seg` is **not
  a quantised confirmation scan**: no chunk quantum, no confirmation count, HetHet and
  A1A1/A1A1 interchangeable where `--ibs` ignores the latter, and nothing ever cut
  (`17-seg-caller.md` §8). Porting the chunk geometry unchanged moved exact rows 705 → 709
  and the worst row 0.2109 → 0.1490 but nearly tripled mean `PropIBD` error,
  0.00138 → 0.00356 — a *partly* right rule, so it was not committed
  (§5.8, §5.9, `tests/parity/fit/segtry.py`).

**The instruments exist — start with them, not with the corpus.**
`docs/research/fixtures/segcanvas.py` paints chromosome 2 one complete word at a time between
all-IBS0 walls and picks the marker spacing so that one ulp of the printed `IBD2Seg` is a
fifth of a marker gap — so the column reads back the number of marker intervals called, and
`c = (−M) mod 64` recovers the number of calls and the number of words exactly. 872 reference
answers are cached, so it re-runs in seconds without the binary; `tests/parity/fit/seg17.py`
scores a candidate against the whole corpus in a second and prints the retired caller beside
it. **Do not drive `segcanvas.py` with a non-reference binary** — it writes whatever it
measures into that cache. Grade a copy instead. (The committed cache is 872 entries; check
`git diff` on `segcanvas_measured.json` before committing anything from that directory.)

**The next experiment worth running**, in the order that gets the most for the least:

1. **The left end, first — it is now the bigger of the two clauses and the cheaper to
   isolate.** `python3 segcanvas.py 5` reproduces the nine `xy…` misses in about a minute.
   The alphabet is already there (`ALPHA` in `segcanvas.py`): sweep the *bit position* of
   the opening word's mismatches the way §4 of `17-seg-caller.md` swept the closing word's,
   and read the left endpoint back off `IBD2Seg` directly. §5 of that write-up fixed the
   right end at "63 markers past the nearest mismatch" by exactly this sweep; the left end
   was assumed symmetric and never swept. If it is not symmetric, nine of the eleven misses
   close at once and the model becomes fully bisected on both ends.
2. **Then the bridging lookahead**, against `Cyzy` and `zyzy`. `predict()` takes `bridge` as
   a parameter, so a candidate can be scored on the exhaustive battery and the corpus
   without touching Rust. Note the trade-off already measured: `bridge="none"` gives a
   better mean and tail but 22 fewer exact `IBD2Seg` rows, so a correct lookahead should
   beat *both* on `IBD2Seg` **and** on the mean before it is believed.
3. **Grade any candidate on the union, not just `IBD2Seg`.** `IBD1Seg + IBD2Seg` is the
   figure `--build`'s `AV.FS` statistic actually reads (§6.2), and the §17 rule bought 74
   rows of `IBD2Seg` while leaving the union at 826 of 982. A rule that fixes lengths
   without fixing the union will not move `apps/bigish__build`.

And one negative to save the next person a week: **do not expect a case count to move.**
The whole `.seg` residual is worth 76 cases, but they only flip when a *file* becomes exact,
and the errors are spread across nearly every dataset. §4.4 is the scoreboard, not §2.

### 5.1 The four small filesets do not grade the caller

`nuclear` (6 samples, 10 000 markers), `missing` (6, 10 000), `sexchr` (10, 6 000) and
`monomorphic` (12, 5 000) each report the same 14 `.seg` rows: the within-family pairs of one
six-person nuclear family. They are an order of magnitude smaller than `bigish`, and
`monomorphic` additionally cycles every tenth marker through *monomorphic*, *monomorphic with
A1 written as PLINK's missing allele*, *MAF 0.001* and *MAF 0.001 with one forced carrier* —
half its markers carry no usable information. On these the reference's own segment numbers
are nowhere near the pedigree the generator built, so agreeing or disagreeing with the
reference there says little about the caller.

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
**0.3384**, while φ = π2/2 + π1/4 applied to its own segment numbers gives **0.2450**. Both
cannot be right. **The reference does not recover the underlying IBD on this fileset**, so it
still grades nothing — but open-king now reproduces it: since the caller change of §4.4 that
row is byte-identical (`0.9800 / 0.0000 / 0.4900 / PO`), where the retired caller printed
0.5582 / 0.4218 and owned the project's worst `PropIBD` error. `monomorphic` is now 12 of 14
rows exact with a mean error of 0.00003, the *lowest* of the ten datasets.

`nuclear` and `missing` were once described as equally unusable, on the strength of
`nuclear N_C1`/`N_C3`. That gap was **mostly the informativeness gate, not a reference
error**: once open-king applies the same rule it prints 0.1240 / 0.2975 for that pair, and
`nuclear`'s mean absolute `PropIBD` error is now 0.0006 and `missing`'s 0.0004. The
"poisoned for fitting" warning is about what the *reference* prints on these four filesets,
and it still holds — agreement there says little about the caller in either direction. No branch in `crates/*/src/` tests a dataset name; dataset names appear in the
crates only in `crates/*/tests/`, as the list a scorecard iterates over.

### 5.2 The reference races on `<prefix>X.kin0`

The between-family X-chromosome writer is not serialised. Six identical invocations of
`king -b sexchr.bed --kinship` produced **six different files** — sizes 665, 662, 662, 662,
187 and 662 bytes — with records torn mid-field and identifier columns from different pairs
interleaved:

```
SEX  S_DAU2  SU3  S_UM  FM  1500  0.323  0.17  0  -0.0351      <- one run
SU3  S_UM    SU4  S_UF  MF  150000  0.331  0.1707  -0.0151     <- another
```

Adding `--cpus 1` makes it deterministic (3 runs, 1 distinct file). No capture made without
`--cpus 1` can be a golden, so the harness excludes `<prefix>X.kin0` from those cases — 8 of
the 13 in the suite, which is the whole of the "8 diff-excluded" line `run_parity.py` prints.
The other 5, all captured with `--cpus 1`, **are** diffed and **are** byte-identical.

### 5.3 The reference prints an uninitialised value in its own banner

Every run prints `--noscreen [<int>]` where the integer is uninitialised memory. Across the
484 captured stdouts that print it, it takes three values: `-1717986816` (469 times),
`-858993408` (14) and `-515396096` (1). It is stable within one build and environment, so the
harness normalizes it on both sides rather than trying to reproduce it. open-king prints
`-1717986816`, which matches the reference on this host but is not meaningfully "correct".

### 5.4 `<prefix>.segments.gz` is never produced

The manual documents it; the 2.3.2 build ships without zlib in its segment writer and writes
no such file on any invocation in the corpus. It is not a parity target.

### 5.5 Scope

Out of scope for v1, and rejected at dispatch rather than silently ignored: `--pca`, `--mds`,
`--roh`, `--lmm`, `--tdt`, `--gdt`, `--risk`, `--makeGRM`, `--plink`, the R plotting flags
(`--rplot`, `--pngplot`, `--rpath`), and multi-dataset input. They are still *accepted* by the
parser so the banner stays byte-exact.

### 5.6 `--unrelated`'s greedy is derived, not tabulated

The order in which `--unrelated` offers a family's members to its greedy used to be a 29-entry
table of measured permutations. The sort behind it is now identified: the reference sorts
members by relative count **descending** with **Sedgewick's median-of-three quicksort** and
walks the sorted array **backwards**; the median-of-three is taken only while `r - l > 7`, and
below that the same partition runs with the plain last-element pivot. Equal keys make every
comparison false, so the algorithm collapses to one fixed permutation per size — which is what
the old table was measuring.

What fixed each piece, all against the reference binary: the **partition**, by inverting the
recursion out of the measured permutations for *n* = 2…30 (twelve independent agreements for
*n* = 9…20); the **cut-off**, uniquely — of 3…11 only `r - l > 7` reproduces all 29
permutations; the **small-subfile partition**, by searching 17 280 parameterised two-pointer
partitions for one matching all seven permutations at *n* = 2…8, which the plain last-element
pivot uniquely does; and the **median-of-three itself**, which equal keys cannot see, on 94
filesets whose counts differ (dropping it breaks 30 of 49).

It then **extrapolates**: *n* = 31, 35, 40, 55 and 70 were predicted before being measured and
all five matched. End to end, 55 fresh filesets the corpus does not contain produce
byte-identical `unrelated.txt` and `unrelated_toberemoved.txt`. Two rules moved with it: the
**screening summary** now uses `--related`'s own inclusion test (reported iff
kinship ≥ `2^-(d+1.5)` **or** `PropIBD` > `2^-(d+0.5)`, bucketed by segment `InfType`), and
**absent parents are materialised** so two rows naming the same unlisted parents are declared
full sibs.

### 5.7 The two-stage screening count is not reproduced

`bigish --related --degree 2` prints `Stages 1&2 (with 32768 SNPs): 36 pairs of relatives are
detected`; open-king prints `50`. This is the only stdout line `--related` gets wrong, and it
is **not** a marker-subset problem. A held-out ladder (*m* = 32 768 / 33 280 / 36 864 /
40 960 / 45 056 / 50 000, both degrees) shows every candidate subset — first, last, evenly
spaced, and highest-MAF markers, and the same four over whole 64-marker words — **overshoots
at degree 2 on every map past 33 280 while matching at degree 1**, and none undershoots. The
stage is lossy; reproducing the number means reproducing the two-stage bound itself. No rule
was fitted to it. It costs 2 cases.

### 5.8 `--ibs` and `.seg` do not share a segment rule

`--ibs`'s `Pr_IBD2` / `MaxIBD2` and `.seg`'s `IBD2Seg` share the word grid, the usable
segments and the per-word masks, but **not the calls**. A controlled fixture settles it: take
a *W*-word all-HetHet block bounded by words in which every marker is an opposite homozygote.
`--ibdseg` reports an `IBD2Seg` worth the block; `--ibs` reports a `MaxIBD2` worth the *whole
usable segment*, running straight through the IBS0 words. So `Scan::ibd2_words` (`--ibs`) and
`Scan::ibd2` (`.seg`) are separate functions over the same masks.

The `--ibs` rules, each fixture-measured (`docs/research/15-ibs-ibd2-rules.md` and
`…/16-segment-extension.md`): a word breaks a run **iff it carries ≥ 5 het-vs-hom
mismatches** — opposite homozygotes and missing calls are irrelevant at any density; one
dirty word is absorbed, two are not; the run is then **confirmed in chunks of five
mismatches**, each needing **≥ 95 HetHet over ≥ 3 words**, and the reported interval stops at
the last confirmed chunk plus a one-mismatch overhang, the whole test being waived where the
run reaches the segment's last two words; and `Pr_IBD2`'s 10 Mb rule gates the **pair**, not
the call, while `MaxIBD2` is never gated. The `.seg` engine's word rule is "dirty iff
IBS0 > 0 or IBS1 ≥ 2" and it must **not** adopt the IBS0-irrelevance finding: that is settled
by §3 of `15-…` and again, from the `.seg` side, by `17-seg-caller.md` §3, where one opposite
homozygote anywhere in a word disqualifies it at any HetHet density and an IBS0 word is never
absorbed. The two callers also disagree on the mismatch threshold (2 against 5) and on
whether a call can be cut short of its run (it cannot, in `.seg`).

**The `--ibs` residual is now zero**, on both IBD2 columns and every dataset. What used to
sit here — a boundary that was "not linear and order-dependent" (smallest reported *h* 19 at
*m*=1, 32 at *m*=2; `m=1,h=16` refused while `m=2,h=32` at the same ratio reported; eight
clean words then eight at *m*=4 reported whole, the reverse order not) — was never an
acceptance boundary at all. It was the chunk quantum seen edge-on: those are `5×19 = 95` and
`3×32 = 96` against `5×16 = 80`, counted *within a chunk* from wherever that chunk starts,
which is exactly why reversing the order changes the answer.

Out of sample the chunk rule is not perfect and is not claimed to be: on 200 fresh random
word canvases the reference agrees with it on **186 (93 %)**, the misses being sequences that
alternate 20- and 64-mismatch words with near-empty ones. Nothing in the corpus resembles
those, which is why the corpus is exact and this is not. `docs/research/16-segment-extension.md`
§10.3 is the standing description of that residual.

### 5.9 The IBD2 geometry: what was unexplained, and what closed it

Recorded so the next person does not re-derive it. The first two entries were the standing
open measurements of `docs/research/14-ibd2-geometry.md`; both are now reproduced by the
caller of `docs/research/17-seg-caller.md`, and the rest is marked where it was superseded:

* **`ibd2gap.py` = 510 — now reproduced.** chr2 entirely IBD2 with one forced IBS0: an
  interior word costs exactly 2 words + 1 marker, a first/last word exactly 1 word,
  **independent of bit position**. This was "the sharpest open measurement in the project"
  against the retired rule, which predicted 595/575/638/607/580 for the five bit positions
  swept. The caller of §4.4 prints **510 for all five**, matching the reference exactly, and
  it falls out of the rule rather than being fitted to it: an IBS0 blocks the reach
  whole-word, and the run after the break is pushed one word.
* **`ibd2end.py` — now reproduced.** A *W*-word IBD2 run bounded by all-IBS0 words reports
  `64W−1`; the retired rule reported `64(W+1)−1`. Our build now prints
  63 / 127 / 191 / 255 / 383 / 511 for *W* = 1, 2, 3, 4, 6, 8, which is the reference's own
  answer on every one.
* **Superseded by `17-seg-caller.md`:** the "4 declined right extensions and 2 denied bridges
  with no separating statistic" and the `{last,first} × {−1,0,+1}` local-optimum sweep were
  both properties of the retired geometry, whose ends were word-aligned. The ends are
  marker-level (63 markers past the nearest mismatch, §5 of that write-up) and the bridge is
  conditional, so neither observation describes the committed caller.
* **The `--ibs` chunk scan does not carry over, and this is measured, not assumed**
  (`docs/research/16-segment-extension.md` §9, `17-seg-caller.md` §8). `.seg` is not a
  quantised confirmation scan at all: no chunk quantum, no confirmation count, HetHet and
  A1A1/A1A1 interchangeable, and a run called whole or refused whole. Porting the chunk
  geometry to `.seg` unchanged (`tests/parity/fit/segtry.py`) moved exact rows 705 → **709**
  and the worst row 0.2109 → **0.1490** but nearly tripled mean `PropIBD` error,
  0.00138 → **0.00356**, so it was **not committed**. The rule that was committed instead
  scores 709 with mean error **0.00037** and worst row **0.0089**.

### 5.10 Two divergences the corpus cannot see

Both were found in the final pass, by driving the reference over a fixture whose sample and
family structure the 13 corpus datasets do not cover. **Neither costs a parity case**, and
neither is a segment-caller problem. They are recorded here because they are the only known
differences that a *user* could hit while the suite stays at 403/480, and because both are
small, well-localised pieces of work rather than open research.

**1. `--ibdseg` does not apply the 100 Mb usable-total floor.** The reference refuses a
fileset whose usable total `D` is under **100 000 000 bp**: it prints

```text
Segments too short.
  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.
```

and writes `<prefix>allsegs.txt` and **nothing else** — no `.seg`, no `splitped.txt`, and
none of the "IBD segment analysis starts at" block. Bisected to the base pair and
deterministic: on a canvas at `D = 99 999 999` the reference took the refusal on **12 of 12**
runs and at `D = 100 000 000` produced a row on **12 of 12**. open-king produces a `.seg` at
both. (Twelve runs rather than one because the reference has a clock-seeded major-allele QC
check that aborts a run at random; a single probe of this fileset is not evidence.)

The constant is **already in the tree and already correct** —
`crates/king-cli/src/analysis/segments.rs`, `INFORMATIVE_BP = 100_000_000`, with
`Segments::informative()` and `console::SEGMENTS_TOO_SHORT` — and `--ibs` already consults
it: on the same below-floor fileset both binaries print `Segments too short.` and write the
same three files (`.ibs`, `.ibs0`, `allsegs.txt`). Only the `--ibdseg` path skips the
check. The site is `analysis::ibdseg::run`, between the `allsegs.txt` write and the
`IBD segment analysis starts at` line. This was *not* fixed in
the final pass on purpose: the same fixture exposes divergence 2 on the same code path, so a
blind one-line patch would have half-fixed a path nobody had measured properly. Measure the
`--related` and `--build` behaviour below the floor first — `--related` was checked and
agrees, but only on a fileset small enough to take the `< 10` sample downgrade, which is not
a real test.

**2. `<prefix>splitped.txt` is written unconditionally.** On a 6-sample fileset whose
families are all singletons, the reference writes **no** `splitped.txt` and prints no
`… is generated for certain pedigree plot applications` line, at `D` above the floor as well
as below; open-king writes and announces it. The corpus cannot distinguish the two rules:
every dataset that reaches the segment pass has a family of at least 4 members
(`unrelated` 10, `bigish` 9, `nuclear` 6, `admixed` 4), and the only datasets with no
multi-member family at all — `singleton` and `pair`, largest family 1 — sit below the `< 5`
sample downgrade with `trio` and never run the pass. So `kingsplitped.txt` is byte-identical in all 50 corpus cases and
still wrong off-corpus.

The obvious hypothesis — *the reference writes it only when some family has more than one
member* — is consistent with every observation above but rests on **one** fixture, so it is
a hypothesis and not a rule. Bisect it with `fixlab.py`: hold the sample count fixed and
vary only the largest family size across 1, 2 and 3.

---

## 6. The two structural gaps, in detail

### 6.1 `<prefix>X.seg`

**Not implemented.** What is measured about it, from the reference:

* **When it is written.** `--ibdseg` writes `<prefix>X.seg` exactly when `--degree` is
  non-zero *and* the fileset has usable X-chromosome segments. Bare `--ibdseg`,
  `--ibdseg --degree 0` and `--ibdseg --seglength 5` write none; `--degree` (bare, = 1),
  `--degree 1`, `--degree 2`, `--degree 3` and `--degree -1` all do. Of the 13 corpus datasets
  only `sexchr` has usable X segments, so only 2 of the 480 cases are affected.
* **Which rows.** Exactly the pairs the autosomal `<prefix>.seg` carries, in the same order —
  verified at `--degree -1` (0 rows in both), 1, 2 and 3 (14 rows in both).
* **Its shape is malformed in the reference.** The header names 11 columns
  (`FID1 ID1 FID2 ID2 Sex1 Sex2 MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD`); every data row
  carries 10 tab-separated fields, the last empty. The three numbers written are `IBD1Seg`,
  `IBD2Seg`, `PropIBD` — checked by arithmetic: row `S_SON1`/`S_SON2` reads
  `0.1462  0.6393  0.7124`, and 0.6393 + 0.1462/2 = 0.7124. They land in the `MaxIBD1`,
  `MaxIBD2` and `IBD1Seg` column positions, and the last two columns are never written.
* **Announced on stdout** as `Additional summary statistics of X-Chr IBD segments saved in
  file <prefix>X.seg`, after the autosomal line.

Implementing it needs an X-chromosome segment caller, which would inherit the residual of
§4.1; the two affected cases would very likely still not be byte-identical.

### 6.2 `--build` on `bigish`

`apps/bigish__build` writes an **empty** `kingbuild.log` and `kingupdateparents.txt` where the
reference writes 19 and 34 lines. This is not a numeric gap: the reconstruction rules the
reference applies here are unimplemented, so nothing is emitted. The other 12 `--build`
datasets are byte-identical because they need none of those rules.

The case's third file, `kingupdateids.txt`, **already matches byte for byte** — the family
numbering, which original families each `KING<n>` absorbs, and the row order are all
correct, so none of the three obvious structural suspects (family numbering, parent
tie-breaking, sex assignment) is what fails here.

Of the two files, `updateparents.txt` alone is structurally derivable — every clustered
member on a tab-separated `FID IID FATHER MOTHER` row in `updateids.txt` order, keeping the
`.fam`'s parents, except that each mutually-full-sib group declaring no parents takes the
next synthetic pair (1,2), (3,4), (5,6) — **one pair per group regardless of its size**,
verified on a three-father sibship. `kingbuild.log` additionally needs `INFERENCE AV.FS`
with a `Join3/Join2` statistic printed to three decimals (`bigish`: 0.778, 0.801, 0.779,
0.827, 0.803) and `INFERENCE HS.UN2`.

**That statistic is no longer unknown, and it confirms the case is blocked on the segment
caller.** Writing `IBD(x, y)` for the union of a pair's called IBD1 and IBD2 segments as a
set of base pairs, for a triple `(R; N1, N2)`:

```text
Join2 = | IBD(R,N1) ∩ IBD(R,N2) |
Join3 = | IBD(R,N1) ∩ IBD(R,N2) ∩ IBD(N1,N2) |
```

Where `R` is IBD to both sibs, a grandparent forces them to have inherited the same
parental haplotype (ratio → 1); an avuncular does not (→ 2/3), which is exactly what the
reference's two message variants discriminate. Scored against **all 53 `AV.FS` values the
reference emitted over 19 filesets** — `bigish` plus 18 held-out pedigrees — the formula is
one-sided high by a mean of **+0.0035** (range −0.0001 … +0.0118), the signature of the
known IBD1 over-call. **Only 5 of the 53 round to the printed three decimals, and none of
`bigish`'s five do.** So even a complete reconstruction implementation leaves all five log
lines wrong: `apps/bigish__build` is blocked on §4.1 exactly as `apps/bigish__cluster` is,
and it is **not** an independent 78th problem.

**The blame is now allocated exactly, and the whole residual is accounted for.** The two
inputs to the statistic fall on opposite sides of §4.1:

| input | reference `IBD2Seg` | our reported union `IBD1Seg+IBD2Seg` |
| --- | --- | --- |
| `IBD(R,N1)`, `IBD(R,N2)` — avuncular, so **`Join2` only** | `0.0000` | exact on **823 of 823** corpus rows |
| `IBD(N1,N2)` — full sibs, entering **`Join3` only** | `> 0` | exact on **3 of 159**, always ours too big |

Measured per triple, `dU` for both `R` pairs is `±0.0000` every time, so the denominator is
right and only the numerator's third set is wrong. An over-call `ΔS` in that set can raise
`Join3` by at most `ΔS` and the ratio by at most `ΔS / Join2`; over **39 triples** (the
corpus five plus fixtures across eight two- and three-family shapes) every residual is
positive and lies inside `[0, ΔS / Join2]` — **39 of 39, nothing left over for a second
cause**. Three alternative readings of the statistic were tried and are worse, so none is a
missing correction: SNP counts instead of base pairs (identical to 4 dp), word-aligned
instead of refined endpoints (mean −0.025), and re-calling at 0…10 Mb minimum length (no
change below 5 Mb, worse at 10). Rig: `docs/research/fixtures/avfs_score.py`.

Two further rules were measured the same way, and one is a sharp negative:

* **The named sib pair belongs to the sibship, not to `R`** — every `AV.FS` line raised
  against one sibship names the same pair whatever `R` is, and where the sibship is the
  `RULE FS0`/`FS1` one it is that sibship's first two members in the order the rule line
  prints them.

  For a **declared** sibship — one `.fam` couple's children, which is what `bigish` names —
  the order is not the `.fam` order, and the previous reading of it as *data-dependent* was
  **wrong**. It is invariant under complete genotype reseeding (4:4 gives `(A_C3 A_C4)` and
  `(B_C1 B_C2)` on all **nine** seeds tried; seven further shapes agree across three seeds
  each) even though every printed `Join3/Join2` moves; invariant under each child's sex
  (five patterns, all-male and all-female included); and invariant under sliding the
  pedigree behind 0…8 extra leading singletons, so it is a *position* in the sibship rather
  than a sample index. It is **not** a function of the sibship's size alone, and not of any
  pairwise statistic: over 19 triples the named pair's rank on `Join2`, `Join3`, the ratio,
  the sibs' mutual `PropIBD` and `Kinship`, and `PropIBD` to `R` each runs from first to
  last. A four-child second family names positions `(1,2)`, `(1,3)`, `(3,4)` or `(3,2)`
  depending on the *first* family's size and on how many unrelated singletons pad the
  cohort — the only input found that moves it while the genotypes do not.
  `avfs_score.py --pairs` re-measures the map. **The generating rule remains
  unidentified**, but it is now known to be searchable without genotypes.
* **The verdict is a cut on the ratio**, `uncle|aunt` below against
  `grandfather|grandmother, HS, or nephew|niece` above, the words following `R`'s sex.
  Bracketed to **(0.848, 0.901)** over the 53 values — which does not separate 0.85, 0.875
  and 0.9.

No partial implementation was attempted, deliberately: writing `updateparents.txt` alone
cannot flip the case, the neighbouring `RULE PO.*` family and the phantom-materialisation
path a non-sibship merge triggers are unrecovered, and the corpus exercises exactly one
shape. The derivation lives in `crates/king-cli/src/analysis/build.rs`'s module doc;
`docs/research/fixtures/avfs.py` regenerates the held-out pedigree shapes and
`docs/research/fixtures/avfs_score.py` — the scorer that was lost, now rebuilt and
committed — drives the reference over them and prints the scorecard and the accounting
bound above in about twenty seconds.

---

## 7. How a case is judged

Each directory under `tests/parity/golden/<group>/<case>/` is one captured reference
invocation: `cmd.txt` (argv, with `{KING}` / `{DATA}` / `{ALT}` placeholders),
`exitcode.txt`, `stdout.txt`, `stderr.txt` verbatim, and every output file the reference wrote
into its working directory. The harness replays `cmd.txt` with our binary in a fresh temporary
directory and compares **exit status, stdout, stderr, the set of files produced, and the bytes
of every file**.

490 invocations are captured; 480 are replayed. The 10 under `core/_analysis/` are
alternate-parameter captures kept for analysis and are included only with
`--include-analysis`.

Five kinds of line are genuinely non-reproducible and are normalized on **both** sides before
comparison — never on ours alone:

1. **Timestamps** — `ctime()` output on every `… starts at` / `… ends at` line.
2. **Thread counts** — `N CPU cores are used`, which follows the host.
3. **Progress tokens** — the `\r`-separated `10%20%…` sequences, whose granularity follows
   timing.
4. **Absolute input paths** — the reference echoes the `.bed`/`.bim`/`.fam` it was given.
5. **`--noscreen [<int>]`** — §5.3.

The one file excluded from diffing entirely is `<prefix>X.kin0` on captures taken without
`--cpus 1` (§5.2). Running the harness with `--impl <reference>` proves the normalization is
not hiding anything: reference against its own captures is **480/480**.

---

## 8. Related documents

* `docs/MAINTAINING.md` — the clean-room rule, repo layout, regenerating the corpus,
  re-capturing goldens, adding an analysis.
* `docs/SPEC.md` — the reference's observable behaviour, flag by flag.
* `docs/BEHAVIOR.md` — raw sweeps behind the rules.
* `docs/VERIFIED_FORMULAS.md` — the estimators, with the experiment that fixed each.
* `docs/research/` — the investigation log, numbered in the order it happened.
  `13-informativeness-gate.md` removed the 188 spurious `.seg` rows; `14-ibd2-geometry.md`
  localises everything left to the IBD2 caller; `15-ibs-ibd2-rules.md` and
  `16-segment-extension.md` are §5.8 — the latter derives the chunk scan, records its 93 %
  out-of-sample accuracy (§10.3) and, in §9, the measured negative that it does not port to
  `.seg`. **`17-seg-caller.md` is the most recent and the one to read first** for anything
  touching `Scan::ibd2`: it derives the committed `.seg` caller constant by constant and
  §8 is the sharp negative that `.seg` is not a quantised confirmation scan at all.
* `docs/research/fixtures/` — the rigs. `fixlab.py` builds a fileset and drives the
  reference (`$KING` repoints it at our build); `gate8.py` brackets the `--degree 1` clause;
  `segfit.py` is the chunk-scan canvas; **`segcanvas.py` is the `.seg`-native canvas of §5.0**
  and carries 872 cached reference answers in `segcanvas_measured.json`; `avfs.py`
  regenerates the `--build` pedigree shapes of §6.2 and `avfs_score.py` scores the `AV.FS`
  statistic over them. Their `work/` output is gitignored and disposable — the JSON cache is
  not, and must only ever be written by the reference binary.
* `tests/parity/fit/` — Python mirrors of the committed engine, kept honest by
  `check_mirror.py`. `chunk.py` keeps the superseded `--ibs` rule alive beside the committed
  one so the before/after scorecard reproduces; `segtry.py` is the `.seg` port trial of §5.9;
  **`seg17.py` scores the committed `.seg` caller and any variant of it over the whole
  corpus in about a second**, printing the retired geometry beside it.
