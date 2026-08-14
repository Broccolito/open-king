# Parity with KING 2.3.2 — the measured claim

This is open-king's authoritative statement of what it reproduces and what it does not.
Every number in it was produced by the commands in §1, against the reference binary `king`
2.3.2. Nothing here is an estimate, and nothing is rounded in our favour.

> **Headline: 397 of the 480 captured reference invocations are byte-identical (82.7 %).**
>
> **82 of the 83 that are not** trace to one cause: the IBD-segment caller places a called
> segment's endpoints within about one 64-marker scan word of the reference, so the segment
> columns (`IBD1Seg`, `IBD2Seg`, `PropIBD`, and the `InfType`/`Error` derived from them) are
> close but not equal on a minority of the rows that carry them. The *set* of pairs reported
> is exactly right everywhere: **0 extra and 0 missing rows** on every output file in the
> corpus. **The remaining 1** is structural: `--build`'s pedigree reconstruction is
> unimplemented (§6.2). `<prefix>X.seg`, also unwritten (§6.1), costs 2 further cases that
> are already inside the 82 — both fail on the segment columns as well.
>
> Nine analyses are byte-identical on all 13 datasets: `--kinship`, `--duplicate`,
> `--bysample`, `--bySNP`, `--autoQC`, `--unrelated`, plus the whole 220-case `params`
> group. `<prefix>allsegs.txt` — the per-segment listing that underlies everything above —
> is byte-identical in all 163 cases that produce it.

---

## 1. Reproducing every number below

```bash
cd /path/to/open-king
cargo build --release

# pass/fail for all 480 cases  -> "397 PASS, 83 FAIL, 480 total"
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

`run_parity.py` and `measure_gaps.py` are Python 3 standard library only, regenerate the
input corpus automatically on first run (~20 s) and need no reference binary. The probes
drive the reference directly. `run_parity.py` exits 0 when every case passed, 1 when at
least one failed, 2 on a harness error.

Measured on the tree this document describes:

| command | result |
| --- | --- |
| `run_parity.py --impl target/release/king` | **397 PASS, 83 FAIL, 480 total**, 874 output files byte-compared, 8 diff-excluded |
| `run_parity.py --impl <reference>` | **480 PASS, 0 FAIL** — the normalization is complete and the goldens are self-consistent |
| `probes/degree_filter.py --ref <reference>` | 38 298 cases, **0 false-keep, 0 false-drop** |
| `docs/research/fixtures/gate8.py` | brackets the `--degree 1` IBD2 clause to (0.0789, 0.0812] |
| `cargo test --workspace` | **277 passed, 0 failed, 1 ignored** |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| `cargo build --release` from a clean tree, offline | succeeds (~8 s, 12 dependencies) |

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
| `--ibs` | **1/1** | 0/1 | **1/1** | 0/1 | **1/1** | 0/1 | 0/1 | 0/1 | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 7/13 |
| flag plumbing + error probes (`params`) | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | **220/220** |
| | | | | | | | | | | | | | | **397/480** |

By capture group: `apps` **89/91**, `core` **68/104**, `ibdseg` **20/65**, `params`
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
  exact; `threegen` and `dups` now pass `--ibs`, and `threegen` passes `--related` (empty
  `.kin`) while failing `--ibdseg`.

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
| `<prefix>X.kin0` | 5 diffable of 13 (§5.2) | **all 5** | 39 | 0 | **100 %** |
| `<prefix>.kin0` | 168 | 160 | 207 986 | 28 | **99.99 %** |
| `<prefix>.ibs0` | 8 | 7 | 20 754 | 2 | **99.99 %** |
| `<prefix>.kin` | 187 | 151 | 12 081 | 834 | 93.10 % |
| `<prefix>X.kin` | 15 | 9 | 195 | 12 | 93.85 % |
| `<prefix>.ibs` | 13 | 7 | 807 | 59 | 92.69 % |
| `<prefix>cluster.kin` | 1 | 0 | 165 | 30 | 81.82 % |
| `<prefix>.seg` | 50 | 5 | 4 172 | 1 271 | 69.53 % |
| `<prefix>build.log` | 8 | 7 | 19 | 19 | see §6.2 |
| `<prefix>updateparents.txt` | 8 | 7 | 34 | 34 | see §6.2 |
| `<prefix>X.seg` | 2 | 0 | — | — | **never written**, §6.1 |

Every `--kinship` case is byte-identical; the 36 differing `<prefix>.kin` cases are all
`--related`. Row identity is matched on the identifier columns before any comparison, so
across the whole corpus there are **0 extra and 0 missing rows** on every file.

**stdout, stderr and exit status.** 469 of the 480 cases match stdout byte-for-byte after
the normalization of §7. **11 cases differ on stdout**, and every one of them also differs
on a file — no case in the suite fails on console output alone:

| cases | stdout line that differs | cause |
| ---: | --- | --- |
| 5 | `core/monomorphic__related*` — the `Inference` relationship-count row (`0 10 2 2 0 1` vs `0 9 3 2 0 1`) | downstream of the segment `InfType` |
| 1 | `ibdseg/monomorphic__related_degree2_ibdseg` — same line | same |
| 2 | `bigish --related --degree 2` — `Stages 1&2 (with 32768 SNPs): 36 pairs` vs `50 pairs` | the two-stage screening bound, §5.7 |
| 2 | `sexchr --ibdseg --degree 2` — the `…X-Chr IBD segments saved in file kingX.seg` line is absent | §6.1 |
| 1 | `apps/bigish__build` | §6.2 |

---

## 4. The gaps, measured

Row counts in this section use `measure_gaps.py`'s denominator, which is **rows inside the
cases that differ**, not rows corpus-wide — it is the tighter, less flattering number. §3
gives the corpus-wide view of the same data.

### 4.1 The segment columns — 82 of the 83 failures

Every one of these is `IBD1Seg` / `IBD2Seg` / `PropIBD` (and the `InfType` / `Error` that
follow from them) being close but not equal.

| file | rows differing | of | +extra | −missing | column | mean abs err | worst |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| `king.seg` | 1 271 | 4 169 | **0** | **0** | `PropIBD` 1 246 rows | 0.00562 | 0.2109 |
| | | | | | `IBD1Seg` 878 rows | 0.02150 | 0.4218 |
| | | | | | `IBD2Seg` 802 rows | 0.01794 | 0.4218 |
| | | | | | `InfType` 7 rows | — | — |
| `king.kin` (`--related`) | 834 | 3 978 | **0** | **0** | `IBD1Seg` 828 rows | 0.01963 | 0.4218 |
| | | | | | `IBD2Seg` 822 rows | 0.01702 | 0.4218 |
| | | | | | `PropIBD` 822 rows | 0.00806 | 0.2109 |
| | | | | | `Error` 12 rows | 0.75 | 1.0 |
| | | | | | `InfType` 6 rows | — | — |
| `king.kin0` (`--related`) | 28 | 296 | **0** | **0** | `IBD1Seg` 28 rows | 0.01461 | 0.0537 |
| | | | | | `IBD2Seg` 28 rows | 0.00943 | 0.0448 |
| | | | | | `PropIBD` 28 rows | 0.00331 | 0.0179 |
| `kingX.kin` (`--related`) | 12 | 90 | **0** | **0** | `IBD2Seg` 12 rows | 0.07905 | 0.1460 |
| | | | | | `IBD1Seg` 12 rows | 0.06815 | 0.1236 |
| `kingcluster.kin` | 30 | 165 | **0** | **0** | `IBD2Seg` 30 rows | 0.00600 | 0.0195 |
| | | | | | `IBD1Seg` 29 rows | 0.00880 | 0.0180 |
| `king.ibs` | 59 | 673 | **0** | **0** | `Pr_IBD2` 59 rows | 0.02447 | 0.5116 |
| | 13 | 673 | | | `MaxIBD2` (bp) | 1.60 × 10⁷ | 5.75 × 10⁷ |
| `king.ibs0` | 2 | 19 327 | **0** | **0** | `Pr_IBD2` 2 rows | 0.00895 | 0.0115 |

`MaxIBD2` — the length in base pairs of the single longest IBD2 segment, and the sharpest
per-segment grader available — is exact on **21 547 of 21 560** `.ibs`/`.ibs0` rows. All 13
differing values are in `.ibs`; `.ibs0`'s 20 753 rows are exact on `MaxIBD2` throughout.
Restricted to the rows that actually grade the IBD2 caller — the **149** `.ibs` rows where
the reference reports a non-zero `MaxIBD2` — **139 are exact** and 10 differ.

`kingcluster.kin` fails only on `bigish`, and only on its three segment columns: the pair
set, the ordering and the other eleven columns are exact. `--cluster` is `--related` re-run
inside the merged families, so this is the segment residual and nothing else — there is no
clustering rule left to find.

### 4.2 The 16-column `--related` layer is complete

Measured directly, over **4 805** comparable `.kin` / `.kin0` rows (3 901 of which have both
`IBD1Seg` and `IBD2Seg` byte-exact):

| column | rows differing where the segments are **exact** | rows differing where the segments already differ |
| --- | ---: | ---: |
| `PropIBD` | **0** | 891 |
| `InfType` | **0** | 6 |
| `Error` | **0** | 12 |

No row anywhere differs on `N_SNP`, `Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0` or
`Kinship`. So `HetConc`, `HomIBS0`-as-union, the `InfType` table, the `Error` grader, row
order, the `.kin0` `N >= 100` gate, the `< 10` sample downgrade and the
`Kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)` inclusion disjunct are all exact, and
**100 % of the `--related` residual is the segment caller**.

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

`PropIBD` differs on more `.seg` rows (1 246) than either of its own inputs (878, 802), which
looks like an arithmetic bug and is not one.

On the primary capture, **115 of the 820 rows whose `IBD1Seg` and `IBD2Seg` are both exact
still print a different `PropIBD`** — and every one of the 115 differs by **exactly ±1 in the
fourth decimal** (101 high, 14 low). Fifteen distinct formulations were scored over those
820 rows — `ibd2 + ibd1/2`, `(ibd1 + 2·ibd2)/2`, `(b2 + b1/2)/d`, `(2·b2 + b1)/(2·d)`,
integer halving of the base pairs, `f32` intermediates, and others — and **all fifteen score
exactly 705/820**. No reassociation and no precision change moves a single row.

What is actually happening is visible in the reference's own output. In a *single*
`--related --degree 2 --ibdseg` invocation, 175 pairs appear in both `king.kin` and
`king.seg`; all 175 carry identical `IBD1Seg` and `IBD2Seg` in the two files, and **50 of
them carry a different `PropIBD`** (e.g. `IBD1Seg 0.5298 / IBD2Seg 0.2884` → `.kin` 0.5532,
`.seg` 0.5533). Same invocation, same pair, same printed inputs, two answers: the
reference's two writers do not agree with each other. open-king computes `PropIBD` once,
from the unrounded estimates, and that choice matches `.kin` exactly (§4.2, 0 of 3 901) at
the cost of these 115 `.seg` rows.

This is recorded as a **known limitation, not a to-do**: an earlier handoff note claimed
~386 `.seg` rows were recoverable "for free" by reassociating the expression. That claim is
refuted by the fifteen-formulation sweep above.

### 4.4 The primary `--ibdseg` scorecard

`<dataset>__ibdseg`, the default 3 Mb reporting floor, 982 rows over 10 datasets:

* **820 of 982 rows** have both `IBD1Seg` and `IBD2Seg` exact at the printed four decimals.
  This is what `cargo test -p king-core --test ibdseg_parity` prints; run it with
  `KING_GOLDEN=tests/parity/golden` for the per-dataset breakdown.
* **705 of 982 rows** are exact on all four printed columns; the 115-row gap is §4.3.
* Mean absolute `PropIBD` error **0.00137**, worst 0.2109.
* **0 extra and 0 missing pairs**, on every dataset.

Held out at other reporting floors, rules unchanged: `--seglength 5` gives 697/982 exact on
all four columns (MAE 0.00142), `--seglength 10` gives 665/982 (MAE 0.00160).

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
six-person nuclear family, over 5 000 to 10 000 markers against `bigish`'s 50 000. Read them
with §5.1 in hand.

**The residual is two-sided**, which is why it reads as a boundary problem rather than a
missing rule: of the 277 inexact rows, `IBD1Seg` is too high on 139 and too low on 21,
`IBD2Seg` too low on 121 and too high on 39.

**And it is entirely the IBD2 caller.** Splitting the same 982 rows by whether the reference
reports any IBD2 at all (`docs/research/14-ibd2-geometry.md` §2):

| reference row | rows | both estimate columns exact |
| --- | ---: | ---: |
| `IBD2Seg == 0.0000` | 823 | **819** |
| `IBD2Seg > 0` | 159 | **1** |

The four exceptions in the first line are all `monomorphic` (§5.1). The IBD1 caller and its
boundary refinement are finished; every failing `ibdseg/*__ibdseg` case fails on IBD2. Two
measurements bound what is wrong, and neither is an offset: `--ibs`'s `Pr_IBD2` — the *same*
usable segments and masks under a word-aligned ruler — is biased **too high** where `.seg`'s
`IBD2Seg` is too **low** (`docs/research/14-ibd2-geometry.md` §3); and on `nuclear`
`N_C1`/`N_C2` the reference's own `.seg` total exceeds its own `Pr_IBD2` total by 22.64 Mb
where the widest the current "word-aligned plus usable-segment fringe" model can reach on
that fileset is 13.58 Mb. Grade IBD2 work with `Pr_IBD2` and `MaxIBD2`, never with the
exact-row count above, which the 823 IBD2-free rows pin to 705 ± 1 under every rule variant
tried.

---

## 5. Known limitations

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
cannot be right. open-king prints 0.5582 / 0.4218 for that pair, which is not the truth
either (0.42 IBD2 for sibs is as wrong in the other direction). **On this fileset neither
implementation recovers the underlying IBD**, and the 0.2109 worst-case `PropIBD` error in
§4.1 is this row.

`nuclear` and `missing` were once described as equally unusable, on the strength of
`nuclear N_C1`/`N_C3`. That gap was **mostly the informativeness gate, not a reference
error**: once open-king applies the same rule it prints 0.1240 / 0.2975 for that pair, and
`nuclear`'s mean absolute `PropIBD` error is now 0.0037, `missing`'s 0.0084, against
`monomorphic`'s 0.0403. The "poisoned for fitting" warning therefore holds only for
`monomorphic`. No branch in `crates/*/src/` tests a dataset name; dataset names appear in the
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

The `--ibs` rules, each fixture-measured (`docs/research/15-ibs-ibd2-rules.md`): a word breaks
a run **iff it carries ≥ 5 het-vs-hom mismatches** — opposite homozygotes and missing calls are
irrelevant at any density; one dirty word is absorbed, two are not; a call is kept iff its
measured words hold **≥ 95 HetHet markers**, waived when the run reaches the segment's last two
words; and `Pr_IBD2`'s 10 Mb rule gates the **pair**, not the call, while `MaxIBD2` is never
gated. This contradicts the `.seg` engine's "dirty iff IBS0 > 0 or IBS1 ≥ 5" word rule;
whether `.seg` should adopt the IBS0-irrelevance finding is **open**, and §3 of that document
is the evidence that the two currently differ.

The residual is a sustained low-grade mismatch in one region: uniform 8-word blocks with *m*
mismatches and *h* HetHet per word are reported only above a boundary that is **not linear and
is order-dependent** (smallest reported *h* is 19 at *m*=1, 32 at *m*=2, 63 at *m*=3;
`m=1,h=16` is refused while `m=2,h=32` — same ratio — is reported; eight clean words followed
by eight at *m*=4 form one interval covering all sixteen, the reverse order reports only the
clean eight). Every remaining `--ibs` error lives there. Per dataset: `bigish` 52 rows (mean
0.0042, worst 0.0141), `monomorphic` 3, and `missing` / `multifam` / `nuclear` / `sexchr` 1
each.

### 5.9 What is unexplained about the IBD2 geometry, precisely

Recorded so the next person does not re-derive it (`docs/research/14-ibd2-geometry.md`):

* **`ibd2gap.py` = 510.** chr2 entirely IBD2 with one forced IBS0: an interior word costs
  exactly 2 words + 1 marker, a first/last word exactly 1 word, **independent of bit
  position** — so interior IBD2 boundaries are word-quantised, not marker-refined. The
  committed rule predicts 595. Read-back resolution is ±6 markers; the gap is 64. This is the
  sharpest open measurement in the project.
* **`ibd2end.py`**: a *W*-word IBD2 run bounded by all-IBS0 words reports `64W−1`; the
  committed rule reports `64(W+1)−1`. The fix costs 2 corpus IBD2 columns and does not satisfy
  the point above, so it was not taken.
* **4 declined right extensions** and **2 denied bridges** that no per-word statistic
  separates from the 91 and 14 that do extend/absorb. Deliberately not fitted — the only
  correlate is left-run length over 16 observations.
* Refinement is at a strict local optimum: sweeping `{last,first} × {−1,0,+1}` for all five
  marker-level rules, the committed choice is best or tied-best on every axis and the nearest
  alternatives cost 20–400 rows. The residual is *not* an off-by-one.

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

Of the two files, `updateparents.txt` alone is structurally derivable — `RULE FS0` gives the
two full-sib fathers synthetic parents numbered globally (KING1 → (1,2), KING2 → (3,4),
KING3 → (5,6)) and everyone else keeps real parents. `kingbuild.log` is not: it needs
`INFERENCE AV.FS` with a `Join3/Join2` statistic reproduced to three decimals (0.778, 0.801,
0.779, 0.827, 0.803) and `INFERENCE HS.UN2`. That statistic is almost certainly
segment-derived, so the case would still fail even if the rules were recovered. No partial
implementation was attempted.

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
* `docs/research/` — the investigation log. `13-informativeness-gate.md` removed the 188
  spurious `.seg` rows; `14-ibd2-geometry.md` localises everything left to the IBD2 caller and
  lists what no feature separates; `15-ibs-ibd2-rules.md` is §5.8.
