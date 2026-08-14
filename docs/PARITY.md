# Parity with KING 2.3.2 — the measured claim

This is open-king's authoritative statement of what it reproduces and what it does not.
Every number in it was measured against the reference binary `king` 2.3.2 — nothing here is
an estimate, and nothing is rounded in our favour. §1 lists the commands that reproduce the
headline, the per-file table and both segment scorecards on the tree this document ships
with; every figure outside those comes from a rig that is **named where the figure is
quoted**, so any claim here can be re-run rather than taken on trust.

> **Headline: 408 of the 480 captured reference invocations are byte-identical (85.0 %).**
> The harness self-check — the reference replayed against its own captures — is **480/480**.
>
> **71 of the 72 that are not** trace to one cause: the `.seg` IBD-segment caller places a
> called **IBD2** segment's endpoints a few markers from the reference, so `IBD2Seg` and the
> `PropIBD` computed from it are close but not equal on a minority of the rows that carry
> them. On the primary `--ibdseg` capture that is **86 of 982 rows**, by 0.00006 of the
> genome on average and never more than 0.0042. The **remaining 1** is structural:
> `--build`'s pedigree reconstruction is unimplemented (§6.2) — and it is *also*
> segment-blocked (its one missing statistic has been identified and measured), so it is not
> an independent second problem. `<prefix>X.seg`, also unwritten (§6.1), costs 2 further
> cases that are already inside the 71.
>
> **What the set of reported pairs does:** it is exactly right everywhere. **0 extra and 0
> missing rows on every output file in the corpus**, all 982 `InfType` labels, and every
> `Error`.
>
> **Seven analyses are byte-identical on every dataset that runs them:** `--kinship`,
> `--duplicate`, `--bysample`, `--bySNP`, `--autoQC` and `--ibs` at 13/13 each and
> `--unrelated` at 26/26, plus the whole 220-case `params` group. `<prefix>allsegs.txt` — the
> per-segment listing that underlies everything above — is byte-identical in all 165 cases
> that produce it.
>
> **The headline moved by five in this project's final campaign, and that understates the
> work.** §3 grades whole files; §4.4 grades rows. Two rule corrections landed since the
> count last stood at 403: the `IBD1Seg` overlap rule (`docs/research/18-ibd1-caller.md`),
> which took `IBD1Seg` from 826 exact corpus rows to **all 982** and flipped **five** cases,
> and the `.seg` IBD2 bridge/gate correction (`17-seg-caller.md` §14), which took the binary
> from 5 723 to **6 000 of 6 000** on constructed canvases and flipped **zero** cases while
> changing **no corpus row at all**. A case turns `PASS` only when *every* row of *every*
> file it writes is byte-exact, and the residual is spread thinly across nearly every
> dataset, so large row-level gains can move the count by nothing. §5.0 says which grader to
> use for what, and why.
>
> **Two divergences live outside the corpus entirely** and therefore cost no case here, but
> a user could hit them: `--ibdseg` does not apply the reference's 100 Mb usable-total floor,
> and `<prefix>splitped.txt` is written unconditionally where the reference sometimes writes
> none. Both were re-measured against the reference for this release; §5.10.

---

## 1. Reproducing every number below

```bash
cd /path/to/open-king
cargo build --release

# pass/fail for all 480 cases  -> "408 PASS, 72 FAIL, 480 total"
python3 tests/parity/run_parity.py --impl ./target/release/king

# how big each remaining gap is: rows, columns, mean and worst absolute error
python3 tests/parity/measure_gaps.py --impl ./target/release/king -q

# per-dataset roll-up for one output file
python3 tests/parity/measure_gaps.py --impl ./target/release/king -q --by-dataset king.seg

# harness self-check: the reference against its own captures must be 480/480
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"

# the row-level .seg scorecard, per dataset
KING_GOLDEN=tests/parity/golden cargo test -p king-core --test ibdseg_parity -- --nocapture

# our binary against the reference on the constructed canvases (§5.0)
GRADE_CACHE=$TMPDIR/g2.json python3 docs/research/fixtures/gradebinary.py target/release/king
GRADE_CACHE=$TMPDIR/g1.json python3 docs/research/fixtures/gradebinary.py target/release/king --ibd1

# the two differential probes that are not part of the capture corpus
python3 tests/parity/probes/degree_filter.py --ref "/path/to/reference/king"
cd docs/research/fixtures && python3 gate8.py

# the Python engine mirror must still reproduce the binary's own output
cd tests/parity/fit && python3 check_mirror.py     # -> "MIRROR OK"
python3 seg17.py && python3 seg18.py               # the corpus scorecards
```

`run_parity.py` and `measure_gaps.py` are Python 3 standard library only, regenerate the
input corpus automatically on first run (~20 s) and need no reference binary. The probes and
the canvas rigs drive the reference directly. `run_parity.py` exits 0 when every case passed,
1 when at least one failed, 2 on a harness error.

Measured on the tree this document describes:

| command | result |
| --- | --- |
| `run_parity.py --impl target/release/king` | **408 PASS, 72 FAIL, 480 total**, 874 output files byte-compared, 8 diff-excluded |
| `run_parity.py --impl <reference>` | **480 PASS, 0 FAIL**, 876 files byte-compared — the normalization is complete and the goldens are self-consistent |
| `probes/degree_filter.py --ref <reference>` | 38 298 cases, **0 false-keep, 0 false-drop** |
| `docs/research/fixtures/gate8.py` | brackets the `--degree 1` IBD2 clause to (0.0789, 0.0829] — its ladder refuses at `PropIBD` 0.0789 and accepts at 0.0829 |
| `gradebinary.py target/release/king` | **6 000 / 6 000** canvases — the release binary against the reference's own readings, `IBD2Seg` |
| `gradebinary.py target/release/king --ibd1` | **540 / 540** on the closed families; **43 / 60** on the one family that is deliberately open (§5.0 item 2) |
| `cargo test -p king-core --test ibdseg_parity` | `TOTAL gold=982 exact=896 infType=982 missing=0 extra=0 meandPropIBD=0.0001 worst=0.0042` |
| `tests/parity/fit/check_mirror.py` | **MIRROR OK** — the mirror reproduces the binary's `.seg` columns on all 982 rows and 861 `MaxIBD2` values, 13 datasets |
| `tests/parity/fit/seg17.py` | committed `exact 747  both 896  ibd2 896  MAE 0.00007  worst 0.0042`; retired geometry `705 / 820 / 822 / 0.00138 / 0.2109` |
| `tests/parity/fit/seg18.py` | committed `exact 747  ibd1 982  ibd2 896  MAE 0.000067  worst 0.0042`; retired overlap rule `709 / 826 / 896 / 0.000365 / 0.0089` |
| `cargo test --workspace` | **290 passed, 0 failed, 1 ignored** |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| a pristine copy of the tree, `cargo build --release --offline` | succeeds in **8.07 s**; `Cargo.lock` has 15 packages — the 3 workspace crates and 12 external |
| that clean-tree binary, re-run through `run_parity.py` and `cargo test` | **408 PASS, 72 FAIL** and **290 passed** — the published counts do not depend on a warm `target/` |

By capture group: `apps` **89/91**, `core` **79/104**, `ibdseg` **20/65**, `params`
**220/220**. By analysis: `--kinship`, `--duplicate`, `--bysample`, `--bySNP`, `--autoQC`,
`--ibs` **13/13** each; `--unrelated` **26/26**; `--cluster` and `--build` **12/13**;
`--related` **40/65**; `--ibdseg` **16/52**; `--related --ibdseg` **4/13**; `params`
**220/220**.

The 72 failures, by shape: **45** in the `ibdseg` group (nine datasets × four `--ibdseg`
captures plus their `--related --ibdseg` case), **25** `core/*__related*` (five datasets ×
five `--degree` variants), **1** `apps/bigish__cluster`, **1** `apps/bigish__build`.

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
| `--related` | **5/5** | **5/5** | **5/5** | 0/5 | 0/5 | **5/5** | **5/5** | 0/5 | **5/5** | 0/5 | **5/5** | **5/5** | 0/5 | 40/65 |
| `--ibdseg` | **4/4** | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | 0/4 | **4/4** | 0/4 | **4/4** | **4/4** | 0/4 | 16/52 |
| `--related --ibdseg` | **1/1** | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | 0/1 | **1/1** | 0/1 | **1/1** | **1/1** | 0/1 | 4/13 |
| `--ibs` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| flag plumbing + error probes (`params`) | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | **220/220** |
| | | | | | | | | | | | | | | **408/480** |

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
  exact; `monomorphic` passes all five `--related` cases; `threegen` passes `--related`
  (empty `.kin`) while failing `--ibdseg`. `--ibs` is in this second category on **all
  thirteen** datasets.

---

## 3. What is byte-identical, per output file

Counted over **every** case that produces the file under any `--prefix`, across all four
groups. "Rows" is every data row in every such case, so these percentages are corpus-wide row
accuracy — not the narrower denominator §4 uses. `<prefix>` stands for `king`, `custom`,
`cus.tom` and `custom.` alike.

| output file | cases | byte-identical cases | data rows | rows differing | rows exact |
| --- | ---: | ---: | ---: | ---: | ---: |
| `<prefix>allsegs.txt` | 165 | **all** | 2 217 | 0 | **100 %** |
| `<prefix>splitped.txt` | 50 | **all** | 1 320 | 0 | **100 %** |
| `<prefix>.con` (`--duplicate`) | 46 | **all** | 62 | 0 | **100 %** |
| `<prefix>bySample.txt` | 17 | **all** | 398 | 0 | **100 %** |
| `<prefix>bySNP.txt` | 13 | **all** | 181 000 | 0 | **100 %** |
| `<prefix>_autoQC_Summary.txt` | 13 | **all** | 101 | 0 | **100 %** |
| `<prefix>_autoQC_snptoberemoved.txt` | 13 | **all** | 23 725 | 0 | **100 %** |
| `<prefix>_autoQC_sampletoberemoved.txt` | 13 | **all** | 2 | 0 | **100 %** |
| `<prefix>_autoQC_updatesex.txt` | 1 | **all** | 1 | 0 | **100 %** |
| `<prefix>updateids.txt` | 2 | **all** | 64 | 0 | **100 %** |
| `<prefix>unrelated.txt` | 26 | **all** | 358 | 0 | **100 %** |
| `<prefix>unrelated_toberemoved.txt` | 26 | **all** | 300 | 0 | **100 %** |
| `<prefix>X.kin0` | 5 diffable of 13 (§5.2) | **all 5** | 52 | 0 | **100 %** |
| `<prefix>.ibs` | 13 | **all** | 807 | 0 | **100 %** |
| `<prefix>.ibs0` | 8 | **all** | 20 754 | 0 | **100 %** |
| `<prefix>X.kin` | 17 | **all** | 225 | 0 | **100 %** |
| `<prefix>.kin0` | 178 | 170 | 228 770 | 28 | **99.99 %** |
| `<prefix>.kin` | 201 | 171 | 12 804 | 450 | **96.49 %** |
| `<prefix>cluster.kin` | 1 | 0 | 165 | 16 | 90.30 % |
| `<prefix>.seg` | 50 | 5 | 4 172 | 1 086 | 73.97 % |
| `<prefix>build.log` | 8 | 7 | 18 | 18 | see §6.2 |
| `<prefix>updateparents.txt` | 8 | 7 | 33 | 33 | see §6.2 |
| `<prefix>X.seg` | 2 | 0 | 28 | 28 | **never written**, §6.1 |

Every `--kinship` case is byte-identical; the 30 differing `<prefix>.kin` cases are all
`--related` (25 in `core`, 5 in `ibdseg`), and all 8 differing `<prefix>.kin0` cases are too.
Row identity is matched on the identifier columns before any comparison, so across the whole
corpus there are **0 extra and 0 missing rows** on every file — the only rows we fail to
produce anywhere are the 28 in the two unwritten `X.seg` files.

**stdout, stderr and exit status.** 475 of the 480 cases match stdout byte-for-byte after
the normalization of §7. **5 cases differ on stdout**, and every one of them also differs
on a file — no case in the suite fails on console output alone:

| cases | stdout line that differs | cause |
| ---: | --- | --- |
| 2 | `bigish --related --degree 2` — `Stages 1&2 (with 32768 SNPs): 36 pairs` vs `50 pairs` | the two-stage screening bound, §5.7 |
| 2 | `sexchr --ibdseg --degree 2` — the `…X-Chr IBD segments saved in file kingX.seg` line is absent | §6.1 |
| 1 | `apps/bigish__build` | §6.2 |

The five `monomorphic --related*` cases used to sit at the head of that table, differing on
the `Inference` relationship-count row (`0 10 2 2 0 1` against `0 9 3 2 0 1`) because a
segment `InfType` downstream of the old IBD2 caller was wrong. They are now **byte-identical
outright** — they are the five cases the `IBD1Seg` work of §4.4 flipped.

---

## 4. The gaps, measured

Row counts in this section use `measure_gaps.py`'s denominator, which is **rows inside the
cases that differ**, not rows corpus-wide — it is the tighter, less flattering number. §3
gives the corpus-wide view of the same data.

### 4.1 The segment columns — 71 of the 72 failures

| file | rows differing | of | +extra | −missing | column | mean abs err | worst |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| `king.seg` | 1 086 | 4 169 | **0** | **0** | `PropIBD` 1 050 rows | 0.000689 | 0.0874 |
| | | | | | `IBD2Seg` 465 rows | 0.001531 | 0.1031 |
| | | | | | `IBD1Seg` 214 rows | 0.003113 | 0.0436 |
| `king.kin` (`--related`) | 450 | 3 888 | **0** | **0** | `IBD2Seg` 420 rows | 0.000503 | 0.0043 |
| | | | | | `PropIBD` 414 rows | 0.000506 | 0.0042 |
| `king.kin0` (`--related`) | 28 | 296 | **0** | **0** | `IBD2Seg` 28 rows | 0.000407 | 0.0018 |
| | | | | | `PropIBD` 28 rows | 0.000400 | 0.0017 |
| `kingcluster.kin` | 16 | 165 | **0** | **0** | `IBD2Seg` 15 rows | 0.000547 | 0.0021 |
| | | | | | `PropIBD` 14 rows | 0.000557 | 0.0020 |

`kingX.kin` is **not** in this table: it is byte-identical in all 17 cases.

**`IBD1Seg` appears on one line only.** It differed on 856 `.seg` rows, 810 `.kin` rows, 28
`.kin0` rows, 12 `X.kin` rows and 30 `cluster.kin` rows before
`docs/research/18-ibd1-caller.md`; it now differs on **214 `.seg` rows and nothing else**,
and all 214 are in the `--seglength 5` (73 rows) and `--seglength 10` (141 rows) cases. At
the default 3 Mb the column is exact on **all 982** primary rows. §4.4 and
`…/18-ibd1-caller.md` §9 say why those two floors are different: above the default floor the
reference merges two IBD1 runs across a short interruption, and the gap it tolerates is 65
marker intervals — which at real marker spacings cannot be under 3 Mb.

`king.ibs` and `king.ibs0` are **not** in this table either: since the chunk scan of §5.8
they are byte-identical in all 13 and all 8 cases, which is `MaxIBD2` and `Pr_IBD2` exact on
all **21 561** rows. That removes the sharpest per-segment grader from the toolbox as well as
from the gap: `MaxIBD2`, the length in base pairs of a *single* IBD2 segment, has no residual
left to point at. Anything that grades the remaining gap has to be read off `.seg`, whose
caller is a different one (§5.8).

`kingcluster.kin` fails only on `bigish`, and only on its segment columns: the pair set, the
ordering and the other columns are exact. `--cluster` is `--related` re-run inside the merged
families, so this is the segment residual and nothing else — there is no clustering rule left
to find.

### 4.2 The 16-column `--related` layer is complete

Measured directly over **4 805** corpus rows — every `.kin`, `.kin0`, `X.kin` and
`cluster.kin` row that carries the 16-column form, in 42 + 17 + 6 + 1 cases. **4 342** of
them have both `IBD1Seg` and `IBD2Seg` byte-exact.

| column | rows differing where the segments are **exact** | rows differing where the segments already differ |
| --- | ---: | ---: |
| `N_SNP`, `Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship` | **0** | **0** |
| `IBD1Seg` | — | **0** |
| `IBD2Seg` | — | 463 |
| `InfType` | **0** | **0** |
| `Error` | **0** | **0** |
| `PropIBD` | 31 | 425 |

`InfType` and `Error` are exact on **every one of the 4 805 rows**, differing nowhere at all
— under the retired caller they missed 6 and 12. No row anywhere differs on `N_SNP`, `Z0`,
`Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0` or `Kinship`. So `HetConc`,
`HomIBS0`-as-union, the `InfType` table, the `Error` grader, row order, the `.kin0`
`N ≥ 100` gate, the `< 10` sample downgrade and the
`Kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)` inclusion disjunct are all exact, and
**100 % of the `--related` residual is the segment caller** — now the IBD2 half of it alone.

**The 31 `PropIBD` rows in the top-right cell are §4.3, not a fourteenth rule.** They are all
on `bigish` (five per `--related` case in six cases, plus one in `cluster.kin`), all
**exactly +1 in the fourth decimal** with the reference higher, and every one of them is an
exact half-way tie: `IBD2Seg + IBD1Seg/2` computed from the *printed* four-decimal columns
lands precisely on a rounding boundary (e.g. `0.5157 / 0.2587` → 0.51655, printed 0.5166
against our 0.5165). They appeared only when `IBD1Seg` became exact — the rows had to enter
the "segments exact" column before the disagreement could show there — and their count grows
with it, which is the signature §4.3 describes. **It is not "compute from the printed
values":** rounding `IBD2Seg + IBD1Seg/2` off the printed columns reproduces the reference's
own `PropIBD` on 4 219 of 4 805 `.kin` rows half-up and 4 280 half-even, and on 3 958 and
3 494 of 4 172 `.seg` rows — neither is the rule, and both are far worse than what is
committed.

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

`PropIBD` differs on more `.seg` rows (1 050) than either of its own inputs (214, 465), which
looks like an arithmetic bug and is not one.

On the primary capture, **149 of the 896 rows whose `IBD1Seg` and `IBD2Seg` are both exact
still print a different `PropIBD`**. Across all 50 `.seg` cases that is 519 rows: 431 where
the reference is one unit of the fourth decimal higher and 88 where it is lower, and 496 of
the 519 are exact half-way ties on the printed inputs. That count has now survived three
successive segment callers (115 of 820, then 116 of 825, now 149 of 896) and it *grows* with
the number of rows whose segments are exact, which is the point: it is not a segment problem.
Fifteen distinct formulations were scored over those rows — `ibd2 + ibd1/2`,
`(ibd1 + 2·ibd2)/2`, `(b2 + b1/2)/d`, `(2·b2 + b1)/(2·d)`, integer halving of the base pairs,
`f32` intermediates, and others — and **all fifteen scored exactly 705/820** at the time.
No reassociation and no precision change moved a single row, and §4.2 adds a sixteenth
negative: rounding off the printed columns is worse than what is committed, in both
directions and in both files.

What is actually happening is visible in the reference's own output. Run the reference once:

```bash
king -b bigish.bed --related --degree 2 --ibdseg --cpus 1 --prefix r
```

**147** pairs appear in both `r.kin` and `r.seg`; **all 147** carry identical `IBD1Seg` and
`IBD2Seg` in the two files, and **43 of them carry a different `PropIBD`** — for example
`IBD1Seg 0.3852 / IBD2Seg 0.3123` prints `PropIBD` **0.5048** in `.kin` and **0.5049** in
`.seg`. Same invocation, same pair, same printed inputs, two answers: the reference's two
writers do not agree with each other, so no single rule can match both. open-king computes
`PropIBD` once, from the unrounded estimates, and that choice matches `.kin` on **4 311 of
the 4 342** rows whose two estimate columns are exact, at the cost of these `.seg` rows.

This is recorded as a **known limitation, not a to-do**: an earlier handoff note claimed
~386 `.seg` rows were recoverable "for free" by reassociating the expression. That claim is
refuted by the sweeps above.

### 4.4 The primary `--ibdseg` scorecard

`<dataset>__ibdseg`, the default 3 Mb reporting floor, 982 rows over 10 datasets. "all four"
is `IBD1Seg`, `IBD2Seg`, `PropIBD` and `InfType` all byte-exact; "both est." is the two
estimate columns only, which is the figure `cargo test -p king-core --test ibdseg_parity`
prints:

| dataset | all four | both est. | `IBD1Seg` | `IBD2Seg` | `InfType` | of | mean abs `PropIBD` err | worst |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `bigish` | 582 | 697 | 763 | 697 | 763 | 763 | 0.00005 | 0.0020 |
| `multifam` | 77 | 95 | 104 | 95 | 104 | 104 | 0.00008 | 0.0017 |
| `threegen` | 30 | 38 | 39 | 38 | 39 | 39 | 0.00011 | 0.0035 |
| `admixed` | 13 | 15 | 16 | 15 | 16 | 16 | 0.00002 | 0.0002 |
| `missing` | 10 | 12 | 14 | 12 | 14 | 14 | 0.00022 | 0.0028 |
| `monomorphic` | 13 | 14 | 14 | 14 | 14 | 14 | 0.00001 | 0.0001 |
| `nuclear` | 9 | 9 | 14 | 9 | 14 | 14 | 0.00009 | 0.0005 |
| `sexchr` | 10 | 13 | 14 | 13 | 14 | 14 | 0.00004 | 0.0003 |
| `dups` | 2 | 2 | 3 | 2 | 3 | 3 | 0.00140 | 0.0042 |
| `unrelated` | 1 | 1 | 1 | 1 | 1 | 1 | 0.00000 | 0.0000 |
| **total** | **747** | **896** | **982** | **896** | **982** | **982** | **0.00006** | **0.0042** |

* **`IBD1Seg` is exact on all 982 rows.** `IBD2Seg` on 896, `InfType` on 982, `PropIBD` on
  752 — the 896 − 747 = **149**-row gap between "both est." and "all four" is §4.3.
* **0 extra and 0 missing pairs**, on every dataset.
* The last four datasets are the small ones — each reports the 14 within-family pairs of a
  single six-person nuclear family, over 5 000 to 10 000 markers against `bigish`'s 50 000.
  Read them with §5.1 in hand.

> The mean above is printed column against printed column, which is what a user diffing two
> files sees. `tests/parity/fit/seg17.py` and `seg18.py` score the *unrounded* estimate
> against the reference's printed value and so report **0.000067** on the same 982 rows; the
> difference is our own rounding to four decimals, not a different rule. Comparisons against
> the retired rules below are quoted on that scale, where both sides are measured the same
> way.

Against the geometry it replaced, both scored by `seg17.py` in the same run: exact rows
705 → **747**, both-columns 820 → **896**, `IBD2Seg` 822 → **896**, mean error
0.00138 → **0.00007** (÷20), worst row 0.2109 → **0.0042** (÷50), 0 extra / 0 missing either
way, and `InfType` from 6 wrong rows corpus-wide to none.

**Held out at other reporting floors, rules unchanged** — neither floor was used to fit
anything:

| floor | all four | `IBD1Seg` | `IBD2Seg` | of | MAE | worst | extra | missing |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 Mb (default) | **747** | **982** | **896** | 982 | 0.00006 | 0.0042 | 0 | 0 |
| `--seglength 5` | 729 | 909 | 880 | 982 | 0.00017 | 0.0598 | 0 | 0 |
| `--seglength 10` | 692 | 841 | 877 | 982 | 0.00039 | 0.0874 | 0 | 0 |

The `IBD1Seg` shortfall at 5 and 10 Mb — and it alone — is the measured-but-unmodelled run
merge of `docs/research/18-ibd1-caller.md` §9 (§5.0 item 2). The worst row at those two
floors got *worse* when the overlap rule landed (0.0552 → 0.0598 and 0.0679 → 0.0874) while
everything else improved; that is recorded rather than hidden, and at the default floor worst
improved 0.0089 → 0.0042.

**The residual is entirely `IBD2Seg`, and it is two-sided**: of the **86** rows whose
`IBD1Seg` or `IBD2Seg` differs (982 − 896), `IBD1Seg` is wrong on **none** — too high on 0,
too low on 0 — and `IBD2Seg` is too high on 52 and too low on 34. `IBD1Seg`'s old one-sided
residual (too high on 156 rows, too low on none) was the signature of the overlap rule of
`…/18-ibd1-caller.md` §6, and it is gone.

**Detection is finished; what is left is length.** Splitting the same 982 rows by whether the
reference reports any IBD2 at all:

| reference row | rows | both estimate columns exact | union `IBD1Seg+IBD2Seg` exact |
| --- | ---: | ---: | ---: |
| `IBD2Seg == 0.0000` | 823 | **823** | **823** |
| `IBD2Seg > 0` | 159 | **73** | **73** |

The first line is clean. The second moved from **2** to **73** when the overlap rule landed.
So the IBD1 caller, its boundary refinement, the way `IBD1Seg` is assembled from it, and the
"is there any IBD2 here at all" question are all finished, and every failing
`ibdseg/*__ibdseg` case now fails only on the *length* of the 86 rows whose IBD2 calls are
still a few markers out.

**How to grade further work on it.** Not with the parity case count, and not with the
exact-row count, which the 823 IBD2-free rows dominate. Not with `--ibs` either — `Pr_IBD2`
and `MaxIBD2` have been exact under every candidate since §5.8. Grade with the **`IBD2Seg`
column and the mean** (`tests/parity/fit/seg17.py` prints both, plus the retired caller), and
out of sample with `gradebinary.py`, which replays the constructed canvases with the shipping
binary and compares against the reference's own cached readings.

---

## 5. Known limitations

### 5.0 The segment residual: what is solved, what is not

Everything in §4 says the same thing from different angles, so here it is once, as a ledger.
Numbered 5.0 rather than 5.1 on purpose: §5.1…§5.10 are cross-referenced from the crates and
from `docs/research/`, so their numbers must not move.

**Solved, and not worth re-deriving** (each measured, with the experiment named):

| piece | status | where |
| --- | --- | --- |
| **The entire non-segment surface** — `N_SNP`, `Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`, and the whole command line | **no row anywhere differs**, over 4 805 rows; 220/220 `params` cases | §3, §4.2 |
| **The acceptance gate** — a run is called iff popcount over its own *complete* 64-marker words of `inf1 = p1_i & p1_j & (p0_i \| p0_j)` (IBD1) or `inf2` (IBD2) is **≥ 10**. `inf2` is `p1_i & p1_j & ~ibs1` — HetHet plus A1A1/A1A1, a het-vs-A1A1 marker being uninformative (`17-…` §14.3, bisected at 10 against 20) | exact **and unique** | `13-informativeness-gate.md`, `17-seg-caller.md` §14.3, `tests/parity/fit/gate_*.py` |
| **Which pairs are reported** (`--degree` inclusion, the `.kin0` `N ≥ 100` gate, the `< 10` and `< 5` sample downgrades) | **0 extra, 0 missing rows on every output file in the corpus**; the degree filter itself **0 false-keep, 0 false-drop over 38 298 cases** | §3, §4.2, `probes/degree_filter.py`, `fixtures/gate8.py` |
| **The per-segment listing itself** — `allsegs.txt` | byte-identical in **all 165** cases | §3 |
| **Denominators and thresholds** — `D` = sum over autosomal `allsegs.txt` rows; `--seglength` inclusive, and applied to each surviving `IBD1Seg` piece on its own | exact at 3, 5 and 10 Mb | §4.4 |
| **The aggregate** — `PropIBD = IBD2Seg + IBD1Seg/2` in `f64`, computed once from the unrounded estimates | 0 differences on 4 311 of the 4 342 rows whose segments are exact; the other 31 are §4.3's rounding ties | §4.2, §4.3 |
| **`InfType` and `Error`** | **no row anywhere differs**, over all 4 805 rows — not merely where the segments are exact | §4.2 |
| **The IBD1 caller, its boundary refinement, its gate and its `IBD1Seg` overlap rule** (`Scan::ibd1`, `ibd1_pieces`) | every clause bisected on an IBD1-native canvas; `IBD1Seg` exact on **all 982** primary rows and on every `.kin`/`.kin0`/`X.kin`/`cluster.kin` row; the binary matches the reference on **540 of 540** canvases in the closed families | `18-ibd1-caller.md`, `fixtures/ibd1canvas.py`, `gradebinary.py --ibd1` |
| **The `--ibs` IBD2 caller** (`Scan::ibd2_words`, the chunk scan) | exact on all **21 561** rows | §5.8 |
| **The `.seg` IBD2 caller's word predicate, gate, reach, push and bridge** (`Scan::ibd2`) | every constant bisected on a `.seg`-native canvas; the binary reproduces the reference on **6 000 of 6 000** canvases; `IBD2Seg` exact on **896 of 982** corpus rows | `17-seg-caller.md` §3–§7 and §14, `fixtures/segcanvas.py`, `gradebinary.py` |

**Not solved — two clauses, and neither is a whole subsystem.**

1. **The length of a called IBD2 segment, on 86 of the 982 primary rows.** This is the whole
   `.seg` residual at the default floor and it lives in `Scan::ibd2`. The one clause with any
   evidence against it is **the left end of a run that opens on a mismatch-carrying word** —
   a two-marker effect. Everything else in that function — the word predicate, the gate, its
   statistic and where it starts counting, the 63-marker reach, the whole-word IBS0 block,
   the one-word push, the bridge and the segment fringes — is a bisection off the reference
   (`docs/research/17-seg-caller.md` §3–§7, §14).

2. **The `--seglength`-triggered IBD1 run merge** (`docs/research/18-ibd1-caller.md` §9).
   Measured and bisected twice — the reference merges two IBD1 runs across a short
   interruption when the interrupting words carry at most **5** opposite homozygotes and the
   gap they leave, `pos[first marker of the later run] − pos[last marker of the earlier
   run]`, is strictly under `--seglength`. That gap is 65 marker intervals, so at real
   spacings it cannot fire at 3 Mb and the default capture never sees it; it is the whole
   `IBD1Seg` residual at `--seglength 5` and `10` (`IBD1Seg` 909 and 841 of 982), and the
   whole of `gradebinary.py --ibd1`'s one open family (43 of 60). It is **not** implemented,
   because the obvious generalisation of it makes those two floors much worse (`IBD1Seg` 795
   against 909 at 5 Mb) — unrelated pairs are full of one-IBS0 interruptions, so whatever
   else the reference requires is still missing.

**The precise shape of what is left** (so the next person can recognise it):

* **Detection is done; length is not.** Split the 982 primary rows by whether the reference
  reports any IBD2: `IBD2Seg == 0` → **all 823** exact; `IBD2Seg > 0` → **73 of 159**.
* **`IBD1Seg` is exact on all 982 rows** at the default floor; the 86 imperfect rows are
  `IBD2Seg` alone, and it is two-sided, 52 high against 34 low (§4.4).
* **The canvas residual is zero, which is itself the problem.** Grading the *binary* against
  the reference's own cached readings:

  | battery | canvases exact | canvases lost |
  | --- | ---: | ---: |
  | `gradebinary.py` — six `IBD2Seg` families | **6 000 / 6 000** | — |
  | ...with §7's fitted bridging lookahead | 5 754 | 246 |
  | ...with §7's gate window (following the reach into the word after next) | 5 979 | 21 |
  | ...with §7's `inf2 = p1 & p1` (het-vs-A1A1 informative) | 5 988 | 12 |
  | ...§7 throughout — the rule before this one | 5 723 | 277 |
  | `gradebinary.py --ibd1` — closed families | **540 / 540** | — |
  | `gradebinary.py --ibd1` — the `--seglength 8` family (item 2 above) | 43 / 60 | 17 |

  Each of the three clauses `17-…` §14 corrected is independently necessary, and the losses
  compose: 246 + 21 + 12 = 279 against 277 when all three are reverted together, so they
  overlap on two canvases and nothing else. Each row is a **separately built binary** —
  `Scan::ibd2` reverted one clause at a time in a scratch copy of the tree — graded against
  the reference's own cached readings.

  **And all four of those binaries are indistinguishable on the corpus.** The §7-throughout
  build scores 408 PASS / 72 FAIL and a `.seg` scorecard identical to the committed one on
  every figure in §4.4 — 747 / 896 / 982 / 896, MAE 0.00006, worst 0.0042, 823 of 823 and 73
  of 159. That is the point of this table: §14 was landed on canvas evidence alone, with the
  corpus held constant to the digit, and only an instrument like this could tell the two
  rules apart.

  So the discriminating cases for the remaining 86 corpus rows are **not in the current
  canvas alphabet**. Both instruments are word-uniform: every canvas word is painted from one
  composition and the usable segment is word-aligned. The corpus's failing rows are not —
  they carry ragged fringes, per-word mismatch counts that vary along a run, and usable
  segments cut by `MAX_MARKER_GAP`.
* **The Rust is the thing being graded, not a model of it.** Those numbers are our release
  binary against the reference binary. `tests/parity/fit/engine.py` is separately asserted
  equal to the binary by `check_mirror.py` on all 982 `.seg` rows and 861 `MaxIBD2` values
  across the 13 datasets.
  Note that the binary scores *better* than the Python model on the mixed IBD1 canvases —
  60/60 against `ibd1canvas.py`'s 57–58/60 — because the model there pairs `predict1()` with
  the §7-era IBD2 rule while the binary carries the §14 one.
* **The `--ibs` solution does not port, and that is measured, not assumed.** `.seg` is **not
  a quantised confirmation scan**: no chunk quantum, no confirmation count, HetHet and
  A1A1/A1A1 interchangeable where `--ibs` ignores the latter, and nothing ever cut
  (`17-seg-caller.md` §8). Porting the chunk geometry unchanged moved exact rows 705 → 709
  and the worst row 0.2109 → 0.1490 but nearly tripled mean `PropIBD` error,
  0.00138 → 0.00356 — a *partly* right rule, so it was not committed
  (§5.8, §5.9, `tests/parity/fit/segtry.py`).

**The instruments exist — start with them, not with the corpus.** The canvas technique is
described for a newcomer in `docs/MAINTAINING.md` §8; in one line, `segcanvas.py` paints
chromosome 2 one complete 64-marker word at a time between all-IBS0 walls and picks the
marker spacing so that one ulp of the printed `IBD2Seg` is a fraction of a marker gap — so
the column reads back the number of marker intervals called, and `c = (−M) mod 64` recovers
the number of calls and the number of words exactly. 6 416 reference answers are cached in
`segcanvas_measured.json` and 1 013 in `ibd1canvas_measured.json`, so both rigs re-run in
under a second without the reference binary. **Never drive either rig with a non-reference
binary** — they write whatever they measure into those caches. `gradebinary.py` exists
precisely so that our build can be graded on the same canvases without touching them.

**The next experiment worth running**, in the order that gets the most for the least:

1. **Build a canvas the residual is actually visible on.** This is the blocker for
   everything else: the closed families are 6 540/6 540 and the corpus is still 86 rows
   short, so no existing fixture separates the committed rule from the right one. Extend
   `segcanvas.Canvas` to paint a *fringe* — a partial word at each end of the usable segment
   — and to vary composition word by word inside one run, then re-run the batteries. Until a
   fixture separates two candidate rules, the corpus cannot choose between them: at 3 Mb it
   did not move for `17-…` §14 at all.
2. **The left end of a run that opens on a mismatch-carrying word.** `xyC` and its eight
   relatives were §7 misses and are exact under §14, so the direct *evidence* for this clause
   is gone and only the corpus residual remains. Sweep the *bit position* of the opening
   word's mismatches the way `17-seg-caller.md` §5 swept the closing word's, before assuming
   the left end is still wrong.
3. **Find what else the §9 run merge requires.** The gap threshold and the 5-IBS0 tolerance
   are both bisected; the missing condition is whatever stops it firing on the one-IBS0
   interruptions that fill unrelated pairs. This is the only lever on the 5 and 10 Mb floors,
   and `gradebinary.py --ibd1`'s 43/60 family is a ready-made scoreboard for it.
4. **Grade any candidate on the union, not just `IBD2Seg`.** `IBD1Seg + IBD2Seg` is the
   figure `--build`'s `AV.FS` statistic actually reads (§6.2). That union was the diagnostic
   that led to `docs/research/18-ibd1-caller.md`: the `17-…` IBD2 rule bought 74 rows of
   `IBD2Seg` while leaving the union at 826 of 982, and the campaign that followed found the
   cause was not the IBD2 geometry at all but the way `IBD1Seg` is assembled from it.

And one negative worth keeping in view: **most of the residual is not worth a case count.**
The `.seg` residual is worth 71 cases, but they only flip when a *file* becomes exact, and
the errors are spread across nearly every dataset — the IBD1 fix flipped five, the most any
single rule change in this project has moved the count, and the §14 correction flipped none
while being demonstrably right. §4.4 is the scoreboard, not §2.

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
still grades nothing — but open-king now reproduces it exactly, where the retired caller
printed 0.5582 / 0.4218 and owned the project's worst `PropIBD` error. `monomorphic` is now
14 of 14 rows exact on both estimate columns with a mean error of 0.00001, the *lowest* of
the ten datasets, and all five of its `--related` cases are byte-identical.

`nuclear` and `missing` were once described as equally unusable, on the strength of
`nuclear N_C1`/`N_C3`. That gap was **mostly the informativeness gate, not a reference
error**: once open-king applies the same rule it prints 0.1240 / 0.2975 for that pair, and
`nuclear`'s mean absolute `PropIBD` error is now 0.00009 and `missing`'s 0.00022. The
"poisoned for fitting" warning is about what the *reference* prints on these four filesets,
and it still holds — agreement there says little about the caller in either direction. No
branch in `crates/*/src/` tests a dataset name; dataset names appear in the crates only in
`crates/*/tests/`, as the list a scorecard iterates over.

### 5.2 The reference races on `<prefix>X.kin0`

The between-family X-chromosome writer is not serialised. Six identical invocations of
`king -b sexchr.bed --kinship` produce **six different files**. Re-measured for this release:
6 runs, 6 distinct files, against 3 runs with `--cpus 1` giving 1. An earlier sampling caught
the writer mid-truncation as well — sizes 665, 662, 662, 662, **187** and 662 bytes — with
records torn mid-field and identifier columns from different pairs interleaved:

```
SEX  S_DAU2  SU3  S_UM  FM  1500  0.323  0.17  0  -0.0351      <- one run
SU3  S_UM    SU4  S_UF  MF  150000  0.331  0.1707  -0.0151     <- another
```

Adding `--cpus 1` makes it deterministic. No capture made without
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
  swept. The committed caller prints **510 for all five**, matching the reference exactly, and
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
  (`docs/research/16-segment-extension.md` §9, `17-seg-caller.md` §8). Porting the chunk
  geometry to `.seg` unchanged (`tests/parity/fit/segtry.py`) moved exact rows 705 → **709**
  and the worst row 0.2109 → **0.1490** but nearly tripled mean `PropIBD` error,
  0.00138 → **0.00356**, so it was **not committed**.

### 5.10 Two divergences the corpus cannot see

Both were found by driving the reference over a fixture whose sample and family structure the
13 corpus datasets do not cover, and both were **re-measured for this release** on a fresh
6-sample fixture, 12 runs per binary per condition, deterministic every time. **Neither costs
a parity case**, and neither is a segment-caller problem. They are recorded here because they
are the only known differences a *user* could hit while the suite stays at 408/480.

The one fixture shows both at once. One autosome, 10 000 markers at 10 kb spacing, six
samples in singleton families — usable total `D` = 99 990 000 bp, just under the floor; and
the same fixture with two more markers, `D` = 100 010 000 bp, just over:

| `D` | binary | prints `Segments too short.` | files written |
| --- | --- | --- | --- |
| 99 990 000 | reference | **yes**, 12/12 | `allsegs.txt` |
| 99 990 000 | open-king | no | `.seg`, `allsegs.txt`, `splitped.txt` |
| 100 010 000 | reference | no | `.seg`, `allsegs.txt` |
| 100 010 000 | open-king | no | `.seg`, `allsegs.txt`, **`splitped.txt`** |

**1. `--ibdseg` does not apply the 100 Mb usable-total floor.** Below it the reference prints

```text
Segments too short.
  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.
```

and writes `<prefix>allsegs.txt` and **nothing else** — no `.seg`, no `splitped.txt`, and
none of the "IBD segment analysis starts at" block. open-king produces a `.seg` at both.

The constant is **already in the tree and already correct** —
`crates/king-cli/src/analysis/segments.rs`, `INFORMATIVE_BP = 100_000_000`, with
`Segments::informative()` and `console::SEGMENTS_TOO_SHORT` — and `--ibs` already consults
it (`analysis/ibs.rs` calls `segs.informative()` in four places; `analysis/ibdseg.rs` calls it
nowhere). On the same below-floor fileset both binaries print `Segments too short.` under
`--ibs` and write the same three files. The site of the missing check is
`analysis::ibdseg::run`, between the `allsegs.txt` write and the `IBD segment analysis starts
at` line. It was left unfixed on purpose: the same fixture exposes divergence 2 on the same
code path, so a blind one-line patch would half-fix a path nobody had measured properly.
Measure the `--related` and `--build` behaviour below the floor first — `--related` was
checked and agrees, but only on a fileset small enough to take the `< 10` sample downgrade,
which is not a real test.

**2. `<prefix>splitped.txt` is written unconditionally.** On the *above-floor* run of the same
6-sample fixture, whose families are all singletons, the reference writes **no**
`splitped.txt` and prints no `… is generated for certain pedigree plot applications` line;
open-king writes and announces it. The corpus cannot distinguish the two rules: every dataset
that reaches the segment pass has a family of at least 4 members (`unrelated` 10, `bigish` 9,
`nuclear` 6, `admixed` 4), and the only datasets with no multi-member family at all —
`singleton` and `pair` — sit below the `< 5` sample downgrade and never run the pass. So
`kingsplitped.txt` is byte-identical in all 50 corpus cases and still wrong off-corpus.

The obvious hypothesis — *the reference writes it only when some family has more than one
member* — is consistent with every observation above but rests on **one** fixture shape, so
it is a hypothesis and not a rule. Bisect it with `fixlab.py`: hold the sample count fixed and
vary only the largest family size across 1, 2 and 3.

---

## 6. The two structural gaps, in detail

### 6.1 `<prefix>X.seg`

**Not implemented.** What is measured about it, from the reference:

* **When it is written.** `--ibdseg` writes `<prefix>X.seg` exactly when `--degree` is
  non-zero *and* the fileset has usable X-chromosome segments. Bare `--ibdseg`,
  `--ibdseg --degree 0` and `--ibdseg --seglength 5` write none; `--degree` (bare, = 1),
  `--degree 1`, `--degree 2`, `--degree 3` and `--degree -1` all do. Of the 13 corpus datasets
  only `sexchr` has usable X segments, so only 2 of the 480 cases are affected, and both fail
  on their autosomal `.seg` numbers as well.
* **Which rows.** Exactly the pairs the autosomal `<prefix>.seg` carries, in the same order —
  verified at `--degree -1` (0 rows in both), 1, 2 and 3 (14 rows in both; the golden file is
  15 lines, header included).
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
reference writes 18 and 33 lines. This is not a numeric gap: the reconstruction rules the
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
known IBD1 over-call the overlap rule of §4.4 has since corrected. **Only 5 of the 53 round
to the printed three decimals, and none of `bigish`'s five do.** So even a complete
reconstruction implementation leaves all five log lines wrong: `apps/bigish__build` is
blocked on §4.1 exactly as `apps/bigish__cluster` is, and it is **not** an independent second
cause on top of the segment residual.

**The blame is allocated exactly, and the whole residual is accounted for.** The two
inputs to the statistic fall on opposite sides of §4.1:

| input | reference `IBD2Seg` | our reported union `IBD1Seg+IBD2Seg` |
| --- | --- | --- |
| `IBD(R,N1)`, `IBD(R,N2)` — avuncular, so **`Join2` only** | `0.0000` | exact on **823 of 823** corpus rows |
| `IBD(N1,N2)` — full sibs, entering **`Join3` only** | `> 0` | exact on **73 of 159** (was 2 before §4.4's `IBD1Seg` rule), the rest still ours too big |

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
`docs/research/fixtures/avfs_score.py` drives the reference over them and prints the
scorecard and the accounting bound above in about twenty seconds.

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
not hiding anything: reference against its own captures is **480 PASS, 0 FAIL**, 876 files
byte-compared.

---

## 8. Related documents

* `docs/MAINTAINING.md` — the clean-room rule, repo layout, regenerating the corpus,
  re-capturing goldens, running the suite, **the canvas technique (§8)**, and adding an
  analysis.
* `docs/SPEC.md` — the reference's observable behaviour, flag by flag.
* `docs/BEHAVIOR.md` — raw sweeps behind the rules.
* `docs/VERIFIED_FORMULAS.md` — the estimators, with the experiment that fixed each.
* `docs/research/` — the investigation log, numbered in the order it happened. It is a
  **log, not a spec**: each document records what was true when it was written, so its
  conclusions may have been superseded by a later-numbered one and a few of its
  `docs/PARITY.md §…` cross-references point at a section numbering this file no longer uses
  (`§11.1` was the IBD1 run-acceptance problem, now §5.0). Where a research doc and this file
  disagree, this file is current. Read the log for *how* a rule was established.
  `13-informativeness-gate.md` removed the 188 spurious `.seg` rows; `14-ibd2-geometry.md`
  localises everything left to the IBD2 caller; `15-ibs-ibd2-rules.md` and
  `16-segment-extension.md` are §5.8 — the latter derives the chunk scan, records its 93 %
  out-of-sample accuracy (§10.3) and, in §9, the measured negative that it does not port to
  `.seg`. **`17-seg-caller.md` and `18-ibd1-caller.md` are the two to read first** for
  anything touching `Scan::ibd2` / `Scan::ibd1`: they derive the committed `.seg` caller
  constant by constant. `17-…` §8 is the sharp negative that `.seg` is not a quantised
  confirmation scan at all; `17-…` §14 is the corrected bridge and gate; `18-…` §6 is the
  `IBD1Seg` overlap rule and §9 the one clause deliberately left out.
* `docs/research/fixtures/` — the rigs. `fixlab.py` builds a fileset and drives the
  reference (`$KING` repoints it at our build); `gate8.py` brackets the `--degree 1` clause;
  `segfit.py` is the chunk-scan canvas; **`segcanvas.py` is the `.seg`-native canvas of §5.0**
  (6 416 cached reference answers) and **`ibd1canvas.py` the same canvas built IBD1-side up**
  (1 013); **`gradebinary.py` grades our build on both** without touching either cache;
  `avfs.py` regenerates the `--build` pedigree shapes of §6.2 and `avfs_score.py` scores the
  `AV.FS` statistic over them. Their `work/` output is gitignored and disposable — the JSON
  caches are not, and must only ever be written by the reference binary.
* `tests/parity/fit/` — Python mirrors of the committed engine, kept honest by
  `check_mirror.py`. `chunk.py` keeps the superseded `--ibs` rule alive beside the committed
  one so the before/after scorecard reproduces; `segtry.py` is the `.seg` port trial of §5.9;
  **`seg17.py` scores the `.seg` IBD2 caller and `seg18.py` the `IBD1Seg` overlap rule** over
  the whole corpus in about a second, each printing the retired rule beside the committed one.
