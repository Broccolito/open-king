# Parity with KING 2.3.2 — the measured claim

This is open-king's authoritative statement of what it reproduces and what it does not.
Every number in it was measured against the reference binary `king` 2.3.2 — nothing here is
an estimate, and nothing is rounded in our favour. §1 lists the commands that reproduce the
headline, the per-file table and both segment scorecards on the tree this document ships
with; every figure outside those comes from a rig that is **named where the figure is
quoted**, so any claim here can be re-run rather than taken on trust.

> **What "the reference" means here, and the caveat that applies to every number below.**
> One binary: `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`, Mach-O 64-bit **arm64**, run on
> macOS (Darwin 25.5.0, Apple silicon). Every rule in `open-king-core::ibdseg` and every figure in
> this document was measured against that one build on that one host. Two consequences,
> both of which a reader should carry:
>
> * **KING's segment numerics are not stable across its own releases**, on KING's own
>   account. The version history transcribed in `docs/research/03-website-manual.md` §14
>   records, verbatim: 2.1.2 *"IBD segment algorithm improved"*; 2.1.3 *"`--ibdseg`,
>   `--related`, `--roh` algorithms improved"*; 2.1.4 up-to-4th-degree inference for
>   `--related`/`--ibdseg`; 2.2.1 a `maxIBD1`/`maxIBD2` fix; 2.2.5 *"`--ibdseg` is
>   substantially improved"*; and 2.2.7 the 2.2.5 `--ibdseg` bug *"completly fixed"* [sic].
>   The algorithm behind those changes has never been published. So "byte-identical to
>   KING" here means **to 2.3.2**, and open-king is **not** claimed to reproduce 2.1.x,
>   2.2.x, or whatever follows 2.3.2 — the `.seg` columns in particular should be expected
>   to differ against those builds. If you are comparing against a different KING,
>   re-capture the goldens (`MAINTAINING.md` §5) before believing any number in this file.
> * **Everything except the segment engine is far less exposed to this.** The kinship and
>   IBS estimators are published (Manichaikul *et al.* 2010) and the file formats are
>   stable; it is the unpublished segment caller, and only it, whose parity is pinned to a
>   single build. The one host-dependent output — §5.3's uninitialised banner integer —
>   is normalized on both sides rather than reproduced.
>
> No cross-build or cross-platform differential has been run. That is a gap in the
> evidence, stated as one, not a claim that none exists.

> **Headline: all 480 captured reference invocations are byte-identical (100 %).**
> The harness self-check — the reference replayed against its own captures — is also
> **480/480**. The implementation run compares **876 output files**; eight documented
> host-unstable reference files are excluded symmetrically.
>
> `<prefix>.seg` is byte-identical in all 50 captures at every captured reporting floor.
> The row-level scorecard is saturated too:
>
> | `--seglength` | all four columns | `IBD1Seg` | `IBD2Seg` | of | extra | missing | MAE | worst |
> | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
> | 3 Mb (default) | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
> | 5 Mb | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
> | 10 Mb | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
>
> Saturating the corpus does not prove universal behavioral identity. §4.6 and §5.10–§5.12
> record held-out testing. The 24-fileset segment battery is 68/72 whole-run exact with
> four value differences among 6,713 rows and no row-set differences; all four require an
> exact 40,000-marker array and reproduce KING's uninitialized multiple-of-64 tail read.
> Safe Rust deliberately does not emulate that undefined behavior. Rare held-out residuals
> in a segment acceptance gate, the sparse PO/FS cutoff, `HomIBS0` ties, `MI_Removal`, and
> unusual reconstruction pedigrees remain quantified in §5/§6 and tracked by issue #11.

---

## 1. Reproducing every number below

```bash
cd /path/to/open-king
cargo build --release

# pass/fail for all 480 cases  -> "480 PASS, 0 FAIL, 480 total"
python3 tests/parity/run_parity.py --impl ./target/release/open-king

# the same, as a regression gate: compare per case AND per output file against the
# recorded outcome, and fail on any difference in either direction
python3 tests/parity/run_parity.py --impl ./target/release/open-king --baseline

# per-file row/column scorecards; every captured file is byte-identical in all its cases
python3 tests/parity/measure_gaps.py --impl ./target/release/open-king -q

# per-dataset roll-up for one output file
python3 tests/parity/measure_gaps.py --impl ./target/release/open-king -q --by-dataset king.seg

# harness self-check: the reference against its own captures must be 480/480
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"

# the row-level .seg scorecard at the default floor, per dataset
KING_GOLDEN=tests/parity/golden cargo test -p open-king-core --test ibdseg_parity -- --nocapture

# the row-level .seg scorecard at ALL THREE captured floors (3 / 5 / 10 Mb) -- §4.4.
# Measured from the binary against the goldens; the table above covers only the default.
python3 tests/parity/fit/scorecard.py                 # add --per-dataset or --residual

# our binary against the reference on the constructed canvases (§5.0)
GRADE_CACHE=$TMPDIR/g2.json python3 docs/research/fixtures/gradebinary.py target/release/open-king
GRADE_CACHE=$TMPDIR/g1.json python3 docs/research/fixtures/gradebinary.py target/release/open-king --ibd1

# the two writer rules, from the captures alone -- no binary, no engine (§4.3, §4.5)
python3 docs/research/fixtures/segwriter.py

# the two differential probes that are not part of the capture corpus
python3 tests/parity/probes/degree_filter.py --ref "/path/to/reference/king"
cd docs/research/fixtures && python3 gate8.py

# THE GRADER THAT STILL DISCRIMINATES (§4.6). The corpus is saturated -- 982/982 rows at
# every floor -- so grade the segment caller on filesets it has never seen, byte for byte.
python3 docs/research/fixtures/oosseg.py --ref "/path/to/reference/king"

# the Python engine mirror must still reproduce the binary's own output, at all three
# floors -- the merge of 20-... is dormant at the default one, so a default-only check
# cannot see it (and for a while did not)
cd tests/parity/fit && python3 check_mirror.py     # -> "MIRROR OK"
python3 seg17.py && python3 seg18.py && python3 seg19.py && python3 seg20.py  # historical
python3 seg21.py && python3 seg23.py               # the committed rule, floor by floor
```

`run_parity.py` and `measure_gaps.py` are Python 3 standard library only, regenerate the
input corpus automatically on first run (~20 s) and need no reference binary. The probes and
the canvas rigs drive the reference directly. `run_parity.py` exits 0 when every case passed,
1 when at least one failed, 2 on a harness error; with `--baseline` it exits 0 only when the
outcome matches `tests/parity/BASELINE.txt` exactly.

Measured on the tree this document describes:

| command | result |
| --- | --- |
| `run_parity.py --impl target/release/open-king` | **480 PASS, 0 FAIL, 480 total**, 876 output files byte-compared, 8 diff-excluded |
| `run_parity.py --impl target/release/open-king --baseline` | `baseline: MATCH (480 case(s))` |
| `run_parity.py --impl <reference>` | **480 PASS, 0 FAIL**, 876 files byte-compared — the normalization is complete and the goldens are self-consistent |
| `probes/degree_filter.py --ref <reference>` | 38 298 cases, **0 false-keep, 0 false-drop** |
| `probes/xseg_probe.py --impl target/release/open-king` | `<prefix>X.seg` out of sample on 1 040 built runs: emission gate **1 040/1 040**, bytes at the default floor **625/625** of the runs whose autosomal `.seg` also matches, raised floors **160/160** (§6.1) |
| `docs/research/fixtures/gate8.py` | brackets the `--degree 1` IBD2 clause to (0.0789, 0.0829] — its ladder refuses at `PropIBD` 0.0789 and accepts at 0.0829 |
| `gradebinary.py target/release/open-king` | **6 000 / 6 000** canvases — the release binary against the reference's own readings, `IBD2Seg` |
| `gradebinary.py target/release/open-king --ibd1` | **540 / 540** on the closed families **and 60 / 60** on the one that was open before the run merge — 600 / 600 in total |
| `segwriter.py` | `.seg`'s `PropIBD` rule consistent on **4 172 / 4 172** reference rows with **0** refutations, refuted on `.kin` (42 rows) and `cluster.kin` (3); row-order block size **uniquely 16** over 2..80 across all 50 `.seg` captures |
| `cargo test -p open-king-core --test ibdseg_parity` | `TOTAL gold=982 row=982 est=982 infType=982 missing=0 extra=0 meandPropIBD=0.000017 worst=0.0001` |
| `tests/parity/fit/scorecard.py` | the three-floor row scorecard of §4.4: `982/982/982` at 3 Mb, at 5 **and at 10**; 0 extra and 0 missing at all three; MAE and worst row a true `0.000000 / 0.0000` at each |
| `tests/parity/fit/check_mirror.py` | **MIRROR OK** — the mirror reproduces the binary's `.seg` columns on all 982 rows **at each of 3 / 5 / 10 Mb** (2 946 rows) and 861 `MaxIBD2` values, 10 datasets |
| `tests/parity/fit/seg19.py` | the `19-…` scorecard, unchanged: `806 / 982 / 982` at 3 Mb, `755 / 910 / 946` at 5, `713 / 844 / 937` at 10 — the "before" the merge is measured against |
| `tests/parity/fit/seg20.py` | the merge's own before/after, both under the retired `.kin` `PropIBD` rule: `19` (no merge) `806 / 755 / 713` exact at 3 / 5 / 10 Mb, `20` (merge) `806 / 795 / 793` |
| `tests/parity/fit/seg21.py` | the push + IBD2-merge corrections, same scale: `20` (previous commit) `806 / 795 / 793`, `21` `806 / 817 / 811` exact, `IBD1Seg` `982 / 982 / 970` and `IBD2Seg` `982 / 982 / 972` |
| `tests/parity/fit/seg23.py` | the window bound and the budget word set, same scale: `21` (previous commit) `806 / 817 / 811` exact, **`23` (committed) `806 / 817 / 820`**, `IBD1Seg` and `IBD2Seg` both `982 / 982 / 982`, MAE at 10 Mb 0.000067 → **0.000022** |
| `docs/research/fixtures/oosseg.py --ref <reference>` | **out of sample**, 24 fresh filesets × 3 floors on 8 unused seeds: **68 / 72** runs byte-identical, **4 of 6 713 rows** value-differing — 0 extra, 0 missing; all four are the deliberate exact-64 safety divergence (§4.6) |
| `tests/parity/probes/segment_residuals.py --ref <reference> --impl target/release/open-king` | merged IBD1 and IBD2 calls both feed the >10 Mb pair filter; 39 999/40 001-marker controls exact; exactly four expected value divergences at 40 000 markers (§4.6) |
| `tests/parity/fit/seg18.py` | `18-…`'s own numbers, unchanged: committed `exact 747  ibd1 982  ibd2 896  MAE 0.000067`; retired overlap rule `709 / 826 / 896` |
| full cached `--build` research replay | **277 / 347** logs byte-identical; another **52** have the same distinct lines and differ only in repetition/count residue; of 18 semantic-or-stale residuals, 2 are truncated debugger artifacts and 2 are held-out `<FID>-><IID>` renaming shapes (§6.2) |
| `docs/research/fixtures/build_shapes.py` | **out of sample**, `updateparents.txt` + the console tail over 20 held-out merge shapes: **18 OK, 0 MISMATCH, 2 skipped** (the still-open renaming shapes) |
| `docs/research/fixtures/clusternum.py score` | **out of sample**, the merge queue over 19 discriminating shapes: queue **19/19** and our binary **19/19**, against family order 7, size 7, largest-kinship 11 (§6.2) |
| `docs/research/fixtures/avscore.py 1 work/*` | **out of sample**, `Join3/Join2` over **297** captured `AV.FS` lines from ~120 filesets: **296 exact** at `%.3lf` on the *reported*-segment reading, against **13** on the raw-call reading it replaces; the miss is the reference's own `2.555` (§6.2) |
| `docs/research/fixtures/bandcut.py sweep` | genotype surgery walking one triple through each verdict edge: `uncle`/silent bisected to **(0.848718, 0.851164]**, silent/ambiguous to **(0.896895, 0.900106]** — two cuts where the doc had one, superseding the 0.056-wide `(0.846, 0.902)` |
| `docs/research/fixtures/battery.py band` | **out of sample**, 14 fresh-seed filesets of forced two-child sibships: **27 / 27** verdicts agree with the three-branch band |
| `docs/research/fixtures/battery.py hs` | **out of sample**, the half-sib candidate gate over 26 candidate pairs: `PropIBD` bracketed to **(0.1868, 0.1878)**, 0 refutations; `Kinship` **refuted**, its silent and firing ranges overlapping |
| `docs/research/fixtures/dupkeep.py` | **out of sample**, which duplicate copy is removed, 10 shapes × 3 seeds: **27/27**, against 21 for "keep the later id" and 6 for "keep the earlier"; the line appears with no `INFERENCE` line in 23 of 27 runs, which is what makes it rule-half |
| `docs/research/fixtures/screenfold.py separation` | **the screen impossibility result**: the kept 32 768 markers' genotypes held bit-identical while the printed count falls **46 → 37** (§5.7) |
| `cargo test --workspace` | all workspace tests pass (one timing probe ignored) |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| a pristine copy of the tree, `cargo build --release` | succeeds in **9.5 s** from cold; `Cargo.lock` has 15 packages — the 3 workspace crates and 12 external |
| that clean-tree binary, re-run through `run_parity.py --baseline` | **480 PASS, 0 FAIL**, `baseline: MATCH` — the published counts do not depend on a warm `target/` or a pre-generated corpus |

By capture group: `apps` **91/91**, `core` **104/104**, `ibdseg` **65/65**, `params`
**220/220**. Every analysis row in the matrix below is complete.

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
| `--cluster` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--build` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--related` | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **65/65** |
| `--ibdseg` | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **52/52** |
| `--related --ibdseg` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| `--ibs` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| flag plumbing + error probes (`params`) | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | **220/220** |
| | | | | | | | | | | | | | | **480/480** |

The `params` group is 220 invocations that exercise the command-line surface rather than one
dataset: every `--prefix` shape, `--cpus`, `--sexchr`, `--degree`, `--minConc`,
`--seglength`, alternate `--fam`/`--bim` inputs, malformed and missing files, and the banner
in each case. All 220 are byte-identical.

**Read the `--ibdseg` row carefully.** It is now **52 of 52** — every dataset, at the bare
invocation, at `--degree 2`, at `--seglength 5` and at `--seglength 10`, `sexchr` included
now that `<prefix>X.seg` is written (§6.1). The two `3/4`s an earlier revision of this table
carried were the `--seglength 10` captures of `multifam` and `bigish`, and they closed with
`23-gap-bound.md`.

The two cells that are not bold are both **`bigish`**, and neither is a segment estimate.
`--related` and `--related --ibdseg` lose one capture each on a single stdout line — the
screening count of §5.7 — while every output file in both cases, `.kin`, `.kin0` and `.seg`
alike, is byte-identical. `--build` loses `bigish` on `<prefix>build.log` (§6.2).

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
| `<prefix>updateparents.txt` | 8 | **all** | 33 | 0 | **100 %** |
| `<prefix>unrelated.txt` | 26 | **all** | 358 | 0 | **100 %** |
| `<prefix>unrelated_toberemoved.txt` | 26 | **all** | 300 | 0 | **100 %** |
| `<prefix>X.kin0` | 5 diffable of 13 (§5.2) | **all 5** | 52 | 0 | **100 %** |
| `<prefix>.ibs` | 13 | **all** | 807 | 0 | **100 %** |
| `<prefix>.ibs0` | 8 | **all** | 20 754 | 0 | **100 %** |
| `<prefix>X.kin` | 17 | **all** | 225 | 0 | **100 %** |
| `<prefix>.kin0` | 178 | **all** | 228 770 | 0 | **100 %** |
| `<prefix>.kin` | 201 | **all** | 12 804 | 0 | **100 %** |
| `<prefix>cluster.kin` | 1 | **all** | 165 | 0 | **100 %** |
| `<prefix>.seg` | 50 | **all** | 4 172 | 0 | **100 %** |
| `<prefix>build.log` | 8 | **all 8** | 18 | 0 | **100 %**, §6.2 |
| `<prefix>X.seg` | 2 | **all 2** | 28 | 0 | **100 %**, §6.1 |

**No compared output file differs anywhere in the 480-case corpus.** Everything above —
including `<prefix>.seg`, which carried 12 differing
rows in the previous revision of this table, and the entire 16-column `--related` layer,
which had 450 differing `.kin` rows and 28 `.kin0` rows two campaigns ago and 16
`cluster.kin` rows one campaign ago — is now byte-identical in **every** case that produces
it. Row identity is matched on the identifier columns before any comparison, so across the
whole corpus there are **0 extra and 0 missing rows** on every file, `X.seg` included now
that it is written (§6.1).

`<prefix>.seg` is byte-identical in **all 50** captures, on all 4 172 rows, at the default
floor and at `--seglength 5` and `--seglength 10` alike. `measure_gaps.py` prints
`king.seg  byte-identical in all 50 case(s)`; there is no numeric residual anywhere in this
table. What that does **not** say is that the caller is exactly right in general — the
corpus can no longer tell. §4.6 is the out-of-sample measurement that can, and it is not
clean.

**stdout, stderr and exit status.** All 480 cases match stdout byte-for-byte after the
normalization of §7, and every case matches stderr and exit status. The former `bigish`
screen and `build.log` failures are resolved in §5.7 and §6.2.

---

## 4. The gaps, measured

Row counts in this section use `measure_gaps.py`'s denominator, which is **rows inside the
cases that differ**, not rows corpus-wide — it is the tighter, less flattering number. §3
gives the corpus-wide view of the same data.

### 4.1 The segment columns — closed on this corpus

**There is no numeric gap left in the corpus.** `measure_gaps.py` reports every output file
this project writes as `byte-identical in all N case(s)` except `<prefix>build.log`, which
carries no numbers at all in the half we write. The table this section used to hold — 12
differing `king.seg` rows of 867, `IBD1Seg` mean absolute error 0.009742 and `IBD2Seg`
0.009080 — is empty.

| file | rows differing | of | +extra | −missing | mean abs err | worst |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `king.seg` | **0** | 4 172 | **0** | **0** | 0.000000 | 0.0000 |
| every other data file | **0** | — | **0** | **0** | 0.000000 | 0.0000 |

The history, on `measure_gaps.py`'s own denominator (rows inside the cases that still
differ): 232 of 1 862 rows before the run merge, 74 of 1 658 after it, 12 of 867 after the
push and IBD2-merge corrections, and **0** now. `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`,
`Error`, `MaxIBD2` and `Pr_IBD2` are exact on all 4 172 `.seg` rows, all 4 805 `--related`
rows and all 21 561 `--ibs` rows.

**What closed it** (`docs/research/23-gap-bound.md`). Both diagnoses the previous revision of
this section recorded were wrong: it was neither the merge's gap acquiring a second bound nor
an invented merge. Two independent faults, both floor-dependent, found by reading the
reference one chromosome at a time (`fixtures/chrprobe.py` mutes every other chromosome for
the probe pair rather than subsetting the `.bim`, which would re-phase the 64-marker word
grid):

* **The floor is asked twice, and the second question is about the gate *window*, not the
  reported call.** A run is emitted only if the span of the words its informativeness gate
  counts over reaches `--seglength / 2` (integer division; IBD2 keeps equality, IBD1 does
  not). It is asked at emit, **after** the merge, so a run the bound refuses still merges and
  the merged window is measured whole. Bisected to the base pair on two independent corpus
  calls — an 11.2066 Mb IBD2 call on `multifam` kept at `--seglength 6.290751` and dropped at
  6.290752, twice its 3 145 375 bp one-word window; the same at `2w+1` on `bigish` — and on
  purpose-built canvases at four marker spacings, 4 of 4 each for IBD1 and IBD2.
* **The IBD1 merge's budget is summed over every word between the two runs**, a gate-refused
  run's words included, while the `MERGE_MAX_WORDS` cap still counts only the *unusable*
  ones. That is `20-seglength-floor.md` §11 item 4, left undecided there. Bisected by
  sweeping a refused run's own het-vs-A1A1 load: the merge turns on at 2, which is exactly
  where those markers take the budget's `V` from 8 to 10.

**What it was worth, on both scales.** Case count 475 → **477**; row scorecard at 10 Mb
`IBD1Seg` 970 → **982** and `IBD2Seg` 972 → **982**, printed-column MAE 0.000046 → **0**. On
`fit/seg23.py`'s retired-`.kin`-rule scale, exact rows 811 → **820** at 10 Mb with MAE
0.000067 → **0.000022** — exact rows up **and** MAE down, which is the landing gate of
`MAINTAINING.md` §8.6. Nothing at 3 or 5 Mb, where neither clause can fire.

**Read this section together with §4.6.** "No numeric gap in the corpus" is a statement about
1 862 rows of simulated data captured in 2026, not a proof of the caller. Out of sample the
same binary gets 6 rows of 6 713 wrong.

### 4.2 The 16-column `--related` layer is complete

Measured directly over **4 805** corpus rows — every `.kin`, `.kin0`, `X.kin` and
`cluster.kin` row that carries the 16-column form, in 42 + 17 + 6 + 1 cases. **All 4 805**
are byte-identical, on every column:

| column | rows differing |
| --- | ---: |
| `N_SNP`, `Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship` | **0** |
| `IBD1Seg`, `IBD2Seg`, `PropIBD` | **0** |
| `InfType`, `Error` | **0** |

So `HetConc`, `HomIBS0`-as-union, the `InfType` table, the `Error` grader, row order, the
`.kin0` `N ≥ 100` gate, the `< 10` sample downgrade and the
`Kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)` inclusion disjunct are all exact, and the
layer is finished. Three campaigns ago it had 450 differing `.kin` rows; the last of them
were `PropIBD` rows that turned out to belong to §4.3.

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
  does not move it. `open_king_core::ibdseg::inf_type` stays ungated because `.seg` genuinely
  uses the plain `IBD2Seg > 0.7`; the gate lives in the `.kin` writers.

### 4.3 `.seg`'s `PropIBD` is a different number from `.kin`'s — and now an exact one

This section used to be titled "not a formula we can fix" and it was wrong. Full write-up:
`docs/research/20-seg-writer.md`; reproduce it with
`python3 docs/research/fixtures/segwriter.py`.

**The reference contradicts itself.** Run it once —
`king -b bigish.bed --related --degree 2 --ibdseg --cpus 1 --prefix r` — and **147** pairs
land in both `r.kin` and `r.seg`. All 147 carry identical `IBD1Seg` and `IBD2Seg` in the two
files. **43** carry a different `PropIBD`, in both directions:

| `IBD1Seg` | `IBD2Seg` | `.kin` | `.seg` |
| --- | --- | --- | --- |
| 0.4885 | 0.2974 | 0.5417 | **0.5416** |
| 0.3852 | 0.3123 | 0.5048 | **0.5049** |

Corpus-wide, 6 captures write both files, 201 pairs appear in both with identical estimates,
and the two files disagree on **54** of them (26.9 %). No single expression matches both, so
each writer gets its own.

**`.seg`'s rule.** It computes `PropIBD` from the two columns it is about to print, after
rounding them to four decimals — not from the underlying totals. With `i1`, `i2` the printed
columns scaled by 10 000:

```
PropIBD = i2 * 1e-4 + i1 * 5e-5          , printed "%.4lf"
```

Over **all 4 172** `.seg` rows in the corpus this is consistent on **4 172** — **0
refutations**. 2 859 rows determine it outright and 1 313 land on an exact decimal half,
where the reference rounds up 1 099 times and down 214; the expression above agrees on every
one of those ties, and `(i1+2·i2)/20000`, `i2/10000+i1/20000`, `(i1/2+i2)/10000`,
`(i1+2·i2)·5e-5`, the printed values as doubles, and integer round-half-up score between
3 804 and 4 086. A 1 313-way coin flip does not come out right by luck.

**That the inputs are the printed columns and not the totals** is decided by 1.7 % of rows:
on those the full-precision value sits at least half a printed ulp from the printed-column
combination, so the full-precision hypothesis predicts ≈ 71 refutations across the 4 172.
Zero were observed. And the rule is `.seg`'s alone — the same test refutes it on `king.kin`
(42 rows) and `kingcluster.kin` (3), both of which open-king reproduces byte for byte using
the full-precision value.

Committed as `open_king_core::ibdseg::seg_prop_ibd`, reaching **one column of one file**:
`InfType`, the `--degree` filter, `--unrelated`'s greedy and `--related`'s `Error` grader
all still read `Segments::prop_ibd`. Worth **806 → 982** byte-exact rows at 3 Mb,
**755 → 900** at 5 Mb and **713 → 832** at 10 Mb.

The sixteen negatives an earlier draft of this section listed — fifteen reassociations of
`ibd2 + ibd1/2` and "round off the printed columns half-up" — were all scored *globally*, on
both files at once, which is why none of them moved: the `.kin` half of every one of them was
already right and the change broke it. The finding was only visible once the two writers were
allowed to differ.

### 4.4 The primary `--ibdseg` scorecard

`<dataset>__ibdseg`, the default 3 Mb reporting floor, 982 rows over 10 datasets. "all four"
is `IBD1Seg`, `IBD2Seg`, `PropIBD` and `InfType` all byte-exact — the `row=` column of
`cargo test -p open-king-core --test ibdseg_parity`; "both est." is the two estimate columns only,
its `est=`.

| dataset | all four | both est. | `IBD1Seg` | `IBD2Seg` | `InfType` | of |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `bigish` | 763 | 763 | 763 | 763 | 763 | 763 |
| `multifam` | 104 | 104 | 104 | 104 | 104 | 104 |
| `threegen` | 39 | 39 | 39 | 39 | 39 | 39 |
| `admixed` | 16 | 16 | 16 | 16 | 16 | 16 |
| `missing` | 14 | 14 | 14 | 14 | 14 | 14 |
| `monomorphic` | 14 | 14 | 14 | 14 | 14 | 14 |
| `nuclear` | 14 | 14 | 14 | 14 | 14 | 14 |
| `sexchr` | 14 | 14 | 14 | 14 | 14 | 14 |
| `dups` | 3 | 3 | 3 | 3 | 3 | 3 |
| `unrelated` | 1 | 1 | 1 | 1 | 1 | 1 |
| **total** | **982** | **982** | **982** | **982** | **982** | **982** |

**Exact, on every row of every dataset**, with 0 extra and 0 missing pairs. Mean and worst
absolute `PropIBD` error are both 0.0000 measured column against column, which is what a user
diffing two files sees. (`ibdseg_parity` prints `meandPropIBD=0.000017`: that statistic
compares our *unrounded* value to the reference's printed one, and since `.seg`'s `PropIBD`
is by construction a multiple of half an ulp, half of the exact rows sit 0.5 ulp from the
printed value they round to. It is a property of the ruler, not an error.)

The three generations of caller behind that row, all on the same 982 rows and the same
unrounded scale, scored by `tests/parity/fit/engine.py`'s pinned parameter bundles:

| | all four | `IBD1Seg` | `IBD2Seg` | MAE | worst |
| --- | ---: | ---: | ---: | ---: | ---: |
| `RETIRED` — word-aligned geometry (`17-…`) | 705 | 822 | 822 | 0.001376 | 0.2109 |
| `FRINGE18` — before the IBD2 fringe (`18-…`) | 747 | 982 | 896 | 0.000067 | 0.0042 |
| `PROP19` — before the writer rule (`19-…`) | 806 | 982 | 982 | 0.000023 | 0.0001 |
| **committed** (`20-…`, `21-…`) | **982** | **982** | **982** | 0.000017 | 0.0001 |

**The three-floor scorecard.** The raised floors were never used to fit the caller's
geometry, and had no part in finding the two writer rules; they *were* the evidence that
something was wrong with the merge and with the floor test, which is why each of those rules
was bisected on constructed canvases against the reference and validated on unused seeds
rather than tuned here (§5.0, "never fit to the corpus"). Measured from the binary against
the goldens by `tests/parity/fit/scorecard.py`:

| floor | all four | `IBD1Seg` | `IBD2Seg` | of | extra | missing | MAE | worst |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 Mb (default) | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
| `--seglength 5` | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
| `--seglength 10` | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |

`MAE`/`worst` here are printed-column against printed-column — what a user diffing two files
sees — so every floor reads a true 0. "All four" equals the number of rows whose two estimate
columns are right at every floor: `PropIBD` adds no error of its own anywhere.

**This table is saturated and can no longer grade a change.** That is the point of §4.6.

**What the last three campaigns moved**, on these same 982 rows and on one consistent scale.
Measured from `fit/engine.py`'s pinned bundles — `replace(BASE, merge=False, merge21=False,
push_fraction=None, window_fraction=None, merge_span="unusable")` is the tree as `19-…` left
it, dropping `merge=False` is `20-…`, dropping `push_fraction`/`window_fraction` in turn is
`21-…`, and `BASE` is what ships:

| floor | | all four | `IBD1Seg` | `IBD2Seg` | MAE | worst |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3 Mb | `19-` no merge | 982 | 982 | 982 | 0.000017 | 0.00005 |
| | `20-` run merge | 982 | 982 | 982 | 0.000017 | 0.00005 |
| | `21-` push + IBD2 merge | 982 | 982 | 982 | 0.000017 | 0.00005 |
| | **`23-` committed** | **982** | **982** | **982** | **0.000017** | **0.00005** |
| 5 Mb | `19-` no merge | 900 | 910 | 946 | 0.000156 | 0.0641 |
| | `20-` run merge | 947 | 959 | 947 | 0.000080 | 0.0111 |
| | `21-` push + IBD2 merge | 982 | 982 | 982 | 0.000017 | 0.00005 |
| | **`23-` committed** | **982** | **982** | **982** | **0.000017** | **0.00005** |
| 10 Mb | `19-` no merge | 832 | 844 | 937 | 0.000384 | 0.0916 |
| | `20-` run merge | 943 | 960 | 945 | 0.000129 | 0.0111 |
| | `21-` push + IBD2 merge | 970 | 970 | 972 | 0.000062 | 0.0081 |
| | **`23-` committed** | **982** | **982** | **982** | **0.000016** | **0.00005** |

**Exact rows up and mean error down at every step, at the raised floors, and nothing
whatever at the default** — where none of these rules can fire on real marker spacings. That
is the bar a change has to clear to land here at all (`MAINTAINING.md` §8.6), and it is why
each of them did. **Nothing has been landed in any campaign that traded exact rows against
MAE.** The clauses of `21-…` and of `23-…` were each graded together and separately
(`fit/seg21.py grid`, `fit/seg23.py grid`); in `23-…`'s grid the window bound is worth
970/972 → 982/982 at 10 Mb, the budget word set two further `IBD1Seg` rows, the IBD1 side of
the window bound **zero on this corpus** (it was landed on canvas evidence — 360/360 held-out
against 328 without it — and it costs nothing here), and the `pre_merge` variant is *worse*
(982/980, worst 0.0536) and was rejected.

**Two scales, and the graders that use them.** `engine.py`'s MAE compares our *unrounded*
`PropIBD` to the reference's printed one, so a fully exact floor reads 0.000017 rather than
0 — half of the exact rows sit half a printed ulp from the value they round to.
`fit/scorecard.py` compares printed against printed, which is what a user diffing two files
sees, so all three floors read **0.000000**. And `fit/seg19.py`, `seg20.py`, `seg21.py` and
`seg23.py` grade `PropIBD` with the retired **`.kin`** rule, so they report the same two
estimate columns but a different `exact`: `seg23.py` prints 806 / 817 / 820 where the tables
here print 982 / 982 / 982. All of them are correct on their own scale; none is
interchangeable with another. When quoting a number, name the grader.

**Detection is finished, and so is length on this corpus.** Splitting the 982 primary rows by
whether the reference reports any IBD2:

| reference row | rows | both estimate columns exact |
| --- | ---: | ---: |
| `IBD2Seg == 0.0000` | 823 | **823** |
| `IBD2Seg > 0` | 159 | **159** |

**How to grade further work on it.** Not with the case count, and not with the exact-row
count at any floor — all four are saturated. Three graders still discriminate, in this order:

1. **`docs/research/fixtures/oosseg.py`** — whole fresh filesets, unused seeds, byte diff
   against the reference. Currently **68 of 72** runs, with no extra or missing rows; §4.6.
2. **The canvases** — `gradebinary.py` (6 000/6 000 `IBD2Seg`, 600/600 `IBD1Seg`),
   `mergelab.py`, `push1.py`, `window1.py`. They grade constructed word sequences, one clause
   at a time, and they are where every landed rule was bisected.
3. **`chrprobe.py`** — when a real row does go wrong, this localises it to one chromosome of
   one pair on the corpus's own data before anyone theorises. Both diagnoses `21-…` §8.1
   recorded were wrong, and this is the instrument that showed it.

A candidate rule lands only if the exact-row counts **and** the MAE both improve (or MAE is
unchanged); anything that trades one for the other is reported with its numbers and rejected.
With the corpus saturated, "improve" now has to be read on graders 1 and 2.

### 4.5 `.seg`'s row order is not `.kin`'s either

`<prefix>.seg` lists its pairs **by 16-sample block**: for each block `b1`, for each block
`b2 ≥ b1`, every reported pair with `i` in `b1` and `j` in `b2`, `i` then `j` ascending.
Every other pairwise file this project writes uses plain index order, and `.seg` does too on
any fileset of 16 samples or fewer — which is nine of the thirteen datasets, and why this
survived every earlier campaign behind the `PropIBD` residual.

Sweeping the block size over 2..80 against the row order of **all 50** captured `.seg`
files, exactly one value reproduces every one: **16**. `threegen` (12 samples) rules out
everything below 12, `multifam` (20) rules out 20 and above including plain index order, and
only 16 survives `bigish` (200 samples, 13 blocks). It is not a threading artifact — the
reference gives the identical order at `--cpus` 1, 2, 4 and 8.

It is a pure permutation: same rows, same values. `measure_gaps.py` matches rows on their
identifier columns before comparing and reported **0 extra, 0 missing** both before and
after, so only a byte diff could see it. `docs/research/20-seg-writer.md` §4;
`analysis::ibdseg::seg_pair_order`.

### 4.6 Out of sample — what the corpus can no longer see

Every number above is measured on 480 captures of 13 simulated datasets. All of them are
now exact on the segment columns, which means the corpus has stopped being a test of the
caller and become a regression guard. So this release also measures the same binary on
**filesets it has never seen**: `docs/research/fixtures/oosseg.py` builds 24 of them with
`generate_corpus.py`'s own simulator on **8 seeds used nowhere else in this repository**, in
three pedigree shapes, runs both binaries at 3, 5 and 10 Mb, and diffs `<prefix>.seg` byte
for byte.

```
72 runs, 6 713 reference rows   ->   68 byte-identical
rows: 0 extra, 0 missing, 4 value-differing
```

**94.4 % of whole-file runs and 99.94 % of rows** — and note how far apart those two are,
which is the same whole-file-versus-row effect the headline warns about, seen off-corpus.
The four differing rows have one shape, recorded here rather than smoothed away:

| shape | rows | what differs |
| --- | ---: | --- |
| one full-sib pair on 2 exact-40 000-marker filesets, at 3 and 5 Mb | 4 | `IBD1Seg` **exact**; safe Rust `IBD2Seg` is lower by 0.0181–0.0182 because it does not reproduce KING's uninitialised tail read; `PropIBD` follows it |

Two controls make the cause precise.

* **Marker-count control.** The same two seeds, genotypes and target pair are exact at
  39 999 markers and after appending one all-missing marker (40 001); only 40 000 differs.
  This is the exact-multiple-of-64 uninitialised read already independently measured in
  §5.11. It is a memory-safety bug in KING, not a caller law, and is deliberately not
  emulated. `tests/parity/probes/segment_residuals.py` runs all twelve controls.
* **The former two missing rows are fixed.** On the held-out distant pair, two IBD1 calls
  of 7.71 and 6.84 Mb remain separate at 3 Mb, then merge to 14.60 Mb at 5/10 Mb. KING
  reports the pair only after that merge. An independent canvas makes the same distinction
  for IBD2. `pair_segments` now applies the >10 Mb filter to the conditioned merged calls;
  the 480-case corpus admits 0 extra/0 missing under both hypotheses, while the held-out
  pair and canvas reject the former unmerged-call implementation.

**This rig also validated the window bound out of sample.** At the time that clause landed,
disabling `WINDOW_FRACTION` scored **60 of 72** against 66: six further filesets wrong,
every one at `--seglength 10`, and none right that the candidate got wrong. The later
pair-filter correction raises the current result to 68/72 without changing that one-way
ablation result.

**Take §4.4 and this section together.** "982 of 982 at every floor" is true, and is a
statement about this corpus. "68 of 72, with no row-set errors" is also true. The remaining
four differences are a deliberate safe-product divergence from undefined C++ behavior, not
an unresolved Rust segment rule.

---

## 5. Known limitations

### 5.0 The ledger: what is solved, what is not

Everything in §4 says the same thing from different angles, so here it is once, as a ledger.
Numbered 5.0 rather than 5.1 on purpose: §5.1…§5.10 are cross-referenced from the crates and
from `docs/research/`, so their numbers must not move.

**Solved, and not worth re-deriving** (each measured, with the experiment named):

| piece | status | where |
| --- | --- | --- |
| **The entire non-segment surface** — `N_SNP`, `Z0`, `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`, and the whole command line | **no row anywhere differs**, over 4 805 rows; 220/220 `params` cases | §3, §4.2 |
| **The acceptance gate** — a run is called iff popcount over its own *complete* 64-marker words of `inf1 = p1_i & p1_j & (p0_i \| p0_j)` (IBD1) or `inf2` (IBD2) is **≥ the gate that run reads off a marker-count table** — 10 below 400 000 markers, **20** from there to 2 000 000, 100 above (§5.14); every corpus dataset is in the first row, which is why 10 looked universal. `inf2` is `p1_i & p1_j & ~ibs1` — HetHet plus A1A1/A1A1, a het-vs-A1A1 marker being uninformative (`17-…` §14.3, bisected at 10 against 20) | exact **and unique** within a tier; the tier boundaries bisected to the marker | `13-informativeness-gate.md`, `17-seg-caller.md` §14.3, `tests/parity/fit/gate_*.py`, §5.14 |
| **Which pairs are reported** (`--degree` inclusion, the `.kin0` `N ≥ 100` gate, the `< 10` and `< 5` sample downgrades) | **0 extra, 0 missing rows on every output file in the corpus**; the degree filter itself **0 false-keep, 0 false-drop over 38 298 cases** | §3, §4.2, `probes/degree_filter.py`, `fixtures/gate8.py` |
| **The per-segment listing itself** — `allsegs.txt` | byte-identical in **all 165** cases | §3 |
| **Denominators and thresholds** — `D` = sum over autosomal `allsegs.txt` rows; `--seglength` inclusive, and applied to each surviving `IBD1Seg` piece on its own | exact at 3, 5 and 10 Mb | §4.4 |
| **The two `PropIBD` rules.** `.kin`/`.kin0`/`X.kin`/`cluster.kin` print `IBD2Seg + IBD1Seg/2` at full precision; `.seg` prints `i2*1e-4 + i1*5e-5` off its own four-decimal columns. The reference disagrees with itself on 54 of 201 pairs it writes into both, so this is two rules and not one | `.kin` family byte-identical on **all 4 805** rows; `.seg` rule consistent on **all 4 172** captured rows with **0** refutations | §4.3, `20-seg-writer.md`, `fixtures/segwriter.py` |
| **`.seg`'s row order** — by 16-sample block, then by index | block size **uniquely 16** over 2..80 across all 50 `.seg` captures | §4.5, `20-seg-writer.md` §4 |
| **`InfType` and `Error`** | **no row anywhere differs**, over all 4 805 rows — not merely where the segments are exact | §4.2 |
| **The IBD1 caller, its boundary refinement, its gate and its `IBD1Seg` overlap rule** (`Scan::ibd1`, `ibd1_pieces`) | every clause bisected on an IBD1-native canvas; `IBD1Seg` exact on **all 982** primary rows and on every `.kin`/`.kin0`/`X.kin`/`cluster.kin` row; the binary matches the reference on **600 of 600** IBD1 canvases — 540 closed plus the 60 of the family that was open until the run merge landed | `18-ibd1-caller.md`, `20-seglength-floor.md`, `fixtures/ibd1canvas.py`, `gradebinary.py --ibd1` |
| **The `--ibs` IBD2 caller** (`Scan::ibd2_words`, the chunk scan) | exact on all **21 561** rows | §5.8 |
| **The `.seg` IBD2 caller** (`Scan::ibd2`) — word predicate, gate, reach, push, bridge and fringe | every constant bisected on a `.seg`-native canvas; the binary reproduces the reference on **6 000 of 6 000** word-aligned canvases and **504 of 504** fringe canvases; `IBD2Seg` exact on **all 982** primary rows at **all three** floors | `17-seg-caller.md` §3–§7 and §14, `19-ibd2seg-residual.md`, `fixtures/segcanvas.py`, `fixtures/fringecanvas.py`, `gradebinary.py` |
| **The gate window's own length bound** (`WINDOW_FRACTION`) — the floor is asked a *second* time, of the span of the gate window rather than of the reported call, at emit and after the merge; IBD2 keeps `>= L/2`, IBD1 is one unit tighter | bisected **to the base pair** on two independent corpus calls and on canvases at four spacings (4/4 each pass); out of sample **353/360** IBD2 and **360/360** IBD1 held-out canvases against 344 and 328 without it, and at landing time **66/72** whole filesets against 60 | `23-gap-bound.md` §1–§4, `fixtures/chrprobe.py`, `fixtures/window1.py`, `fixtures/oosseg.py` |
| **The IBD1 merge's budget word set** — summed over *every* word between the two runs, a gate-refused run's included, while the word cap still counts only the unusable ones | bisected on a refused run's own het-vs-A1A1 load: merge off at 0–1, on at 2, which is where `V` crosses 8 → 10 | `23-gap-bound.md` §5, `fixtures/mergelab.py` |

**The `--seglength` run merge and the one-word push** (`20-seglength-floor.md`, corrected by
`21-push-merge.md` and `23-gap-bound.md`) — `Scan::merge_ok`, `Scan::join_runs`,
`Scan::join_runs2`, the `armed` flag. Both are committed and both are exact at 3, 5 and
10 Mb. **The two passes are not the same rule**, which is what `20-…` got wrong:

* **Shared.** After the gate has refused what it refuses, two runs are joined iff the gap
  between them is **strictly** under `--seglength` and a budget `cost·(bad − 2) ≤ X` passes
  over the interrupting words. A run the gate refused lies *inside* an interruption rather
  than ending one. The conditioned merged calls from both passes **do** feed the >10 Mb
  pair filter: one held-out IBD1 pair and a separate IBD2 canvas distinguish that rule from
  the former unmerged-call implementation (§4.6).
  The gap rule itself is exactly right and was bisected to the base pair on real data — a
  `multifam` IBD1 pair splits at a run-to-run gap of 9 652 629 bp and merges at 9 652 630
  (`23-…` §5).
* **IBD1.** At most **2** unusable words; the gap measured run-to-run; `bad` = opposite
  homozygotes; `cost` 4; `X` = A1A1/A1A1 unless the het-vs-A1A1 markers alone reach 10.
* **The cap and the budget do not read the same words** on the IBD1 pass. The cap counts
  only the *unusable* words, so a gate-refused run between two runs is stepped over by it;
  the budget is summed over **every** word of the interruption, that run's included
  (`23-…` §5). `20-…` §11 item 4 left this open and it is now bisected.
* **IBD2.** **No word cap at all** — a purpose-built fixture joins fifteen unusable words
  where the IBD1 pass refuses three at any floor. The interruption runs between the two
  runs' **gate windows**, not between the runs, so the word an earlier run reaches into is
  not part of it (and that holds after any usable word, a gate-refused run's reach word
  included). `bad` = opposite homozygotes **plus** het-vs-hom mismatches; `cost` 3;
  `X` = **HetHet**, switching to A1A1/A1A1 below 10 — the IBD1 pass's own clause, not `inf2`.
* **The push is conditional.** `17-…` §6 read "every call after the first in a usable
  segment starts one word later" as unconditional; it is armed only when a call reaches
  **half** the floor, measured from its own gate-start word, and once armed it stays armed
  for the rest of the segment. Bisected to the base pair on three spacings: at a floor of
  5 080 001 bp a 2 540 000 bp call arms it, at 5 080 100 it does not. At the default floor
  the condition is almost always true, which is why §6 could not see it.

Every constant was bisected on `fixtures/mergelab.py`, `fixtures/push1.py` and
`fixtures/window1.py` against the reference and validated out of sample — 360/360 held-out
canvases at 5 and 10 Mb on three unused seeds for the IBD1 merge, 357/360 on three further
unused seeds for the corrected IBD2 rule (the committed-at-the-time rule scored 343/360 on
the same canvases), 353/360 and 360/360 on three more unused seeds for the two sides of the
window bound (against 344 and 328 for the rule it replaced), plus 600/600 independently drawn
interruptions and, after the pair-filter correction, 68/72 whole filesets (§4.6). **Never
fitted to the corpus.** Between them
they took the headline 464 → 472 → 475 → **477**, `IBD1Seg` at 10 Mb from 844 to 982 of 982,
and both raised floors to byte-exact.

**There is no segment residual on this corpus at any captured floor.** All 982 primary rows
are byte-exact on all four printed fields at 3, 5 and 10 Mb, printed-column MAE a true
0.000000 at each.

**Current held-out ledger.** The exact-multiple-of-64 safety divergence below remains. The
other supported-core residuals are one constructed segment-acceptance counterexample, the
unknown data-derived sparse PO/FS cutoff, last-digit `HomIBS0` ties, the approximate
`MI_Removal` predicate, and rare reconstruction trigger/repetition and cross-family parent
shapes. [`CONTINUATION.md`](CONTINUATION.md#remaining-supported-core-work) records the
current discriminator and evidence for each.

**Historical ledger.** Items 2 and 3 below describe the state before the build-log and
two-stage-screen fixes. They are retained as derivation history and are superseded by the
current 480/480 result and §6.2.

1. **KING's exact-multiple-of-64 tail read is intentionally not reproduced.** 0 corpus
   rows, **4 of 6 713 rows on 24 fresh filesets** (§4.6, `fixtures/oosseg.py`). The four
   value differences occur on two independent 40 000-marker seeds at 3/5 Mb. Holding the
   target pair and its genotypes fixed, 39 999 and 40 001 markers are exact; only 40 000
   differs. This is the uninitialised array-tail read independently pinned in §5.11, not an
   unresolved segment rule. `probes/segment_residuals.py` asserts the safe divergence.

2. **Before the reconstruction fix, `<prefix>build.log` was written only down to its
   `RULE` lines** (§6.2). 1 case. Its
   header, `Duplicate … is removed.` and `RULE FS0`/`FS1` lines now land — 6 of `bigish`'s 18
   lines, 243 of its 806 bytes, every one byte-identical, and byte-identical on **53 of 59**
   held-out shapes.

   **The segment blocker is gone.** `Join3/Join2` intersects three pairs' segment *sets*, so
   it is the only statistic in KING's output that reads segment **placement** rather than a
   total — and this release reproduces it. The sets are the pairs' **reported** segments
   (IBD2 calls plus the IBD1 pieces that clear `--seglength` once the IBD2 calls are cut
   out), not the raw calls; on that reading the ratio is exact at `%.3lf` on **296 of 297**
   captured `AV.FS` lines, `bigish`'s five included, with no triple off by more than 0.0005.
   The single miss is the reference printing `2.555`, a ratio above 1, in the one capture
   where the line sits under a `Reconstruct parent-offspring pair`. The previous entry here
   — that byte-identical totals do not imply identical segments, so `--build` was blocked by
   the caller — measured a formula error, not a caller error; it is **withdrawn**, and the
   caller is exonerated. `fixtures/segprobe`, `fixtures/avscore.py`.

   What remains is one thing: **the sibship's internal member order**, which every inference
   line uses to pick the two people it names. It gates all five `AV.FS` lines, and through
   them the blank lines and the `HS` block. The verdict rule around it is settled — a
   three-branch band, `< 0.85` uncle, `> 0.90` grandparent-or-HS-or-nephew, **nothing at all
   in between** — with both edges bisected to about 0.003 by genotype surgery, and the HS
   candidate gate is `PropIBD > 0.1875`, bracketed to 0.001. See §6.2.

3. **Before the screening fix, `--related`'s two-stage count differed**, 2 cases (§5.7).
   Not a segment problem: one
   stdout line, `36 pairs` against our `50`; every output file in both cases, `.kin0` and
   `.seg` alike, is byte-identical. This is the one gap the project has closed *negatively*
   rather than positively, and the negative is the deliverable:

   * the count is **not the kinship over any subset of markers** — proved by exact algebra
     (`E[N_l] = 4pq(1−2φ)`, `E[het_l] = 2pq`, so numerator and denominator are both
     proportional to `Σ pq` over *whatever* index set they are summed on, making every subset
     and every non-negative weighting unbiased for the same φ), and confirmed three ways by
     measurement;
   * it is **not a merge** of markers into 32 768 slots — refuted by a no-step-at-budget
     scan, by three arrangements of one duplicate multiset giving one count, and by scoring
     every rank grouping under every idempotent operation;
   * it is **not a function of the kept markers at all** — the kept set's genotypes are held
     bit-identical while the printed count falls 46 → 37;
   * it is **sharp**, not two noisy estimates intersected — 0 inversions over 48 ladder pairs.

   What is left is a measured law (`k_screen = 0.5 + R(k − 0.5)`, `R ≡ 1` for `m ≤ 32 768`,
   pinned to 0.2 %) with one constant, `R`, that is **deliberately not fitted** because it
   swings 0.998–1.085 with the MAF spectrum, and a **second necessary condition** that binds
   with no budget involved. §5.7 carries the full proof and an explicit list of what a future
   maintainer should *not* attempt.

**The next experiment worth running, followed by standing warnings.** The two segment
targets formerly at the head of this list are closed: the four value rows are the deliberate
exact-64 safety divergence, and the merged-call pair filter is pinned on both passes.

1. **Chase the screen's *second* condition, not its budget** (§5.7). It is the most
   tractable thing left in `--related`, and the reason is methodological: it is **binary**
   (accept / `No close relatives are inferred.`), its effect is **huge** (a pair at kinship
   0.20006 refused outright), and it fires at `m = 32 768` where no budget is involved and the
   screen is otherwise the exact whole-map kinship to 4e-6 — so one run reads one point and
   nothing has to be bisected. `screenfold.py gate` is the rig. The companion search is for a
   statistic that **degrades with the discarded markers' informativeness while the kept bits
   are held fixed**, which is the one shape `separation` leaves open; a deterministic
   bound-based early exit over the informativeness-sorted map has exactly that shape and
   would supply the second condition for free. **Do not fit `R`** — that is the whole point
   of §5.7, and a fitted `R` reproduces `bigish` and nothing else.

2. **Three hazards any new rig must respect** (`21-…` §8.2, §8.3), or it will misgrade its
   own boundary rows:
   * The reference **stops behaving like a floor** outside `1 ≤ L ≤ 10` Mb. Above ~10 Mb a
     14.06 Mb call reports as a constant 8.93 Mb at every larger floor; below 1.0 the flag
     behaves as though absent. No rule here is measured in either regime and none should be.
   * `.seg`'s floor test is strictly **`>`**, not `>=`: a call measuring exactly 5 100 000 bp
     is reported at `--seglength 5.099999` and dropped at `5.100000`. The engine compares
     `>=` and it never bites on the corpus — nothing there lands on an exact tie — but every
     fixture in the push/merge rig does.
   * Canvas read-back is a *measurement*, not an inspection: see `MAINTAINING.md` §8.3.
3. **Do not re-sweep the caller's constants.** Forty single-knob perturbations and all 32
   combinations of the two IBD1 endpoint rules crossed with the two IBD1 fringe rules were
   scored: none improves exact rows, none beats the committed MAE, and the committed values
   are the unique maximum of that grid (`20-seg-writer.md` §6). Likewise the merge's own
   knobs, swept in `fit/seg20.py grid` and `fit/seg21.py grid` — where dropping `reach` costs
   982 → 959/947 at 5 Mb, `hethet` 982 → 981/980, `push_half` one row on each column at
   10 Mb, and `no_cap` nothing on this corpus but four canvases out of sample.
4. **Five knobs the corpus cannot see at all** — `bridge_rule="17"`, `gate_end="right"`,
   `inf2_ibs1b=True`, `ibd1_clip_ibd2=True`, `clip_before_len=False` all score identically to
   the committed engine on every corpus row. They were settled on the canvases (`17-…` §14)
   and the canvases remain the only evidence for them. If you change one, grade it there.

**And the standing caveat, restated where it bites hardest.** Everything above is measured
against **one** KING build: `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`, Mach-O 64-bit arm64,
on macOS Darwin 25.5.0 — see the note at the top of this file. The segment caller is the part
of KING whose numerics its own release notes say changed most across 2.1.x–2.2.x (2.1.2 "IBD
segment algorithm improved"; 2.1.3 "`--ibdseg`, `--related`, `--roh` algorithms improved";
2.2.1 a `maxIBD1`/`maxIBD2` fix; 2.2.5 "`--ibdseg` is substantially improved"; 2.2.7 that
same change "completly fixed" [sic]), and the algorithm behind those changes is unpublished.
**"982 of 982 at every floor" is a statement about 2.3.2 on this host and nothing wider.**
Against a 2.1.x or 2.2.x build the `.seg` columns should be expected to differ, and no
cross-build or cross-platform differential has been run — that is a gap in the evidence,
stated as one. If you are comparing against a different KING, re-capture the goldens
(`MAINTAINING.md` §5) before believing any segment number in this file. The published
estimators (`--kinship`, `--ibs`, the 16-column `--related` layer) and the file formats are
far less exposed to this: they are documented, and they have not moved across releases.

**Closed since the previous revision of this section**, recorded so it is not re-derived:

* **`dups`' duplicate pair**, previously the largest single error in the corpus (0.0641 at
  5 Mb, 0.0916 at 10), is **exact at all three floors** — it was the run merge, on the IBD2
  pass. The old note asked "what does the reference do with an IBD2 call the floor would
  drop?"; the answer is that it never had one, because it had merged the runs first.
* **The whole `--seglength 5` floor.** It was 35 rows and 3 cases short, one-sided the other
  way, and it is now as byte-exact as the default. `missing` is exact at every floor.
* **The push-counter hypothesis this section used to lead with** — "when a merge joins two
  runs that would each have been counted separately, does the reference increment `emitted`
  once or twice?" — is answered and it was the wrong question. There is no counter: the push
  is armed by a *length* test on a call, not by a count of calls (`21-…` §2). The
  measurement that would have settled it — build it both ways on a canvas above the default
  floor — is the one that did, which is the useful part of the lead.
* **The whole `--seglength 10` floor**, and with it the last two `.seg` parity cases. It was
  12 rows in 2 cases and it is now byte-exact. **Both** of the diagnoses the previous
  revision of this section published were wrong — it was not a second bound on the merge's
  gap (the gap rule is exactly right, bisected to the bp on real data) and it was not an
  invented merge. It was the floor being asked a second time, of the gate window, and the
  IBD1 merge's budget reading a wider set of words than its cap (`23-gap-bound.md`). The
  lesson is item 1 of the experiment list above: **localise a wrong row before theorising
  about it.** Two campaigns guessed which segment of which pair was at fault and both
  guessed wrong; `chrprobe.py` answered it in an afternoon by muting chromosomes.
* **"An absolute cap on the merged span, or on the gap"** — ruled out, measured. `20-…` §2's
  gap rule needs no second bound of any kind.

**And one lesson about graders.** For most of this project the case count was the *wrong*
number to optimise: the `.seg` residual was spread thinly across nearly every dataset, so a
correction could win 277 canvases and flip zero cases (`17-…` §14), while the row-level
scoreboard moved every time. That inverted in the final pass. Once the *numbers* were right,
the remaining differences were a printing rule and a row order — things only a byte diff can
see, which the row-level graders had been silently normalising away for the whole project
(`measure_gaps.py` matches rows on their identifiers before comparing, and reported 0 extra
and 0 missing throughout). Keep both scoreboards, and when one saturates, look at the other.

**And one about where a check is pointed.** `check_mirror.py` asserts that `fit/engine.py`
reproduces the binary exactly, and it is the guard that keeps the research mirror honest.
For an interval it ran **only at the default floor** — and so passed while the binary had
the run merge and the mirror did not, because the merge cannot fire at 3 Mb on the corpus's
spacings. The check was green and the mirror was wrong. It now runs at 3, 5 and 10 Mb
(2 946 rows), and the general rule it teaches is: **a rule whose predicate reads a CLI
parameter must be exercised at a value of that parameter where it is live**, or the check
only proves the rule is dormant. `fit/scorecard.py` exists for the same reason — the
`ibdseg_parity` Rust test covers the default floor alone.

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

Deliberately outside the minimal product scope: `--pca`, `--mds`, `--roh`, `--lmm`, `--tdt`,
`--gdt`, `--risk`, `--makeGRM`, `--plink`, the R plotting flags (`--rplot`, `--pngplot`,
`--rpath`), and multi-dataset input. They remain accepted by the compatibility parser but do
not run an analysis. They are documented product exclusions rather than relatedness-parity
failures; see [SCOPE.md](SCOPE.md).

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

### 5.7 The two-stage screening derivation — resolved

The dense and sparse screening paths now reproduce both captured candidate counts exactly;
the two former failing invocations pass byte for byte. The remainder of this subsection is
the historical negative-search record that constrained the final implementation.

Before the fix, `bigish --related --degree 2` printed
`Stages 1&2 (with 32768 SNPs): 36 pairs of relatives are detected` in KING and `50` in
open-king. That stdout line cost 2 cases — `core/bigish__related_degree2` and
`ibdseg/bigish__related_degree2_ibdseg`. The consequence was contained to that one line:
`.kin0`'s row set comes from the exhaustive re-estimate below it and is byte-correct at every
degree, including in both of those cases, and both of their `.seg` files are byte-identical.

**What the stage is, measured.** Four rigs, each driving the reference over constructed
filesets and reading the count off its stdout: `docs/research/fixtures/screencanvas.py`
(single-pair probe, clone canvas), `…/screenweight.py` (what the statistic weights),
`…/screendeflate.py` (the affine law and the subset proof) and `…/screenfold.py` (the merge
and subset-function refutations, and the second condition). The last two each carry a `facts`
subcommand that re-measures everything below; the full record is
**`docs/research/22-screen.md`**. `screenfold.py`'s core instrument is a **ladder fileset** —
48 pairs climbing through the cutoff, so a single run reads the effective threshold instead of
bisecting for it, 17× cheaper.

The sharp instrument is a **dilution bisection**: replace one member of a real `bigish` pair,
at a growing random marker set, with genotypes drawn from that fileset's own allele
frequencies — synthetic and unrelated to everybody, so no new relative pair appears — and
bisect. The boundary lands to one marker, ~1e-5 in kinship.

* **The stage is per-pair.** Run the candidate pairs one at a time and the count is 0 or 1
  each time; it is 1 on **36** of them at degree 2 — the whole fileset's number. (At degree 1
  the same sum is 17 against the fileset's 18, so per-pair is exact only to ±1.) Not a budget,
  a cap, a ranking or a per-block bound.
* **The law is affine about 0.5, one `R` per fileset.** With
  `k_screen = 0.5 + R*(k − 0.5)`, 36 bisections at each cutoff on `bigish` give
  `R = 1.02257 ± 0.00065` at 0.0625 (boundary kinship 0.07216) and `1.02079 ± 0.00062` at
  0.1250 (0.13264) — 0.2 % apart, where a multiplicative rule needs 0.866 and 0.943 and a
  constant offset 0.0097 and 0.0076. A synthetic flat-MAF fileset, deflating four times harder,
  gives the same verdict with a 25× longer lever arm: `R` = 1.0798 / 1.0838 while `cut/k*`
  moves 0.659 → 0.812.
* **`R` is exactly 1 whenever *m* ≤ 32 768** — 0.99999 ± 0.00001 on `bigish`, 1.00000 on three
  synthetic MAF spectra. This is the constraint that kills most candidate estimators outright.
* **The deflation is systematic, not sampling noise.** Realisation spread of the boundary is
  0.0018 against a deflation of 0.0089, and the per-pair labels are a sharp threshold: every
  `bigish` pair above kinship 0.0731 accepted, every one below 0.0718 rejected, one inversion
  inside a 0.0009-wide window.

**It is not the kinship over any subset of markers, and that is now a proof.** At a marker of
frequency *p*, `E[N_l] = 4pq(1 − 2φ)` and `E[het_l] = 2pq` — the `p²q²` terms cancel exactly —
so the numerator `het_i + het_j + 4·IBS0 − 2·HetHet` and the denominator `min(het_i, het_j)`
are both proportional to `Σ pq` over whatever index set they are summed on. **Every subset and
every non-negative per-marker weighting is unbiased for the same φ.** Three measurements agree:
top-K-by-MAF subsets of `bigish` count 47/45/44/48 pairs over `2^−4` at K = 50 000/32 768/25 000/16 384
and 41/41/41/40 on its first 16 384 markers, flat, where the reference gives 36; **replicating
a map r times** leaves every kinship bit-identical (KING's own `.kin0` confirms it) and still
moves the count, 41 → 36 → 33 → 29 → 27 at r = 2…6, i.e. `R` = 1.000/1.021/1.037/1.055/1.065;
and the one loophole — a subset chosen from data that includes the pair — is simulated and
closed, top-32 768 by in-sample MAF giving `R` = 0.995 ± 0.002 and by in-sample heterozygote
count 0.916 ± 0.003, a bias of the *wrong sign*.

Nine permutations of `bigish`'s marker order print 36/18 every time, which retires prefixes,
strides and word decimations for good; the boundary bisection, forty times finer, moves by
0.0004 — 5 % of the deflation, the size of a tie-break inside an informativeness ranking and no
more.

**What the deflation does track.** It needs the markers overflowing the 32 768 budget to be
*informative*, and grows with how much equally-informative material overflows. Appending 17 232
markers at MAF 0.02 to `bigish`'s first 32 768 leaves the count at exactly its *m* = 32 768
value (50/18, where a `bigish`-sized deflation would read 42); appending the real tail gives
36/18. Two-point MAF maps at *m* = 65 536 put `R`'s minimum (1.0081) exactly where the budget
need not split a tied group, rising to 1.0596 at K = 40 000 and 1.0776 at K = 50 000. Across
spectra at *m* = 50 000: flat 1.080, uniform 1.033, `bigish` 1.022, low-MAF-heavy beta 1.007 —
and one beta point sits *below* one (0.9980 ± 0.0003), so the earlier "never below 1" is
retired. `R` is not a function of `(m, n)` alone: `bigish` at *m* = 50 000 reads 1.0216 while
its first 25 000 markers replicated twice — same *m*, same *n* — read 1.0280.

What survives is a shape, recorded as a shape and not as a rule: when a map holds more
equally-informative markers than 32 768, the reference reaches its budget by something lossy
applied **uniformly across markers** rather than by keeping some and dropping others. Measured
directly — on a flat-MAF map a contiguous clone block grown from marker 0, from 20 000, from
32 768, or backwards from the tail hits the boundary at 0.0957/0.0957/0.0930/0.0940, with no
preference for the head of the file.

#### The two mechanisms that survived the proof are now closed too

A second round (`22-screen.md` §§7–13, instrument `docs/research/fixtures/screenfold.py`,
whose `facts` subcommand re-measures all of it) took the two candidates §5 had left standing
— a merge of markers into 32 768 slots grouped by informativeness rank, and two stages
intersected — and refuted both. **Nothing was landed; the negative got sharper.**

* **Merging into 32 768 slots is dead, three ways.** (i) *No step at one marker over
  budget*: appending `bigish`'s real tail one marker at a time to its 32 768-marker prefix
  prints 50 at every `m` through 33 024 and only then ramps, where any block merge
  (`blockSize = ceil(m/32768)`) flips every slot to a pair at `m = 32 769` and must step
  there. (ii) *Same multiset, three arrangements, identical count*: 32 768 markers plus
  8 192 duplicates print **41 against a true 47** whether the copies are appended in order
  (so `j mod 32768` pairs each marker with its own copy), interleaved directly after their
  originals (so consecutive-pair merging does), or shuffled. Every idempotent merge
  operation — `or`, `and`, `max`, take-the-more-informative — is *lossless* for at least one
  of those three, and all three lose the same six pairs. A stable rank sort makes tied copies
  adjacent, so rank-block grouping dies on the same fixture. (iii) *Scored*: rank-stride and
  rank-block groupings under `or`/`and`/`xor`/saturating-sum, on the sparse encoding where
  merging junk is free, either accept every pair (`or`, `and` — wrong sign) or destroy the
  estimate (`xor`, sum — already wrong where the reference is exact).

* **The statistic is not a function of the markers the budget keeps** — which retires the
  whole "function of a marker subset" family, not just the specific subsets tried. Hold
  `m = 50 000` as 32 768 markers at MAF 0.45 plus 17 232 at MAF `x`, and vary `x`. The
  top-32 768 by allele count is the MAF-0.45 group with **zero index swaps** through
  `x = 0.25`, and — the column the rig prints and the one that matters — the kept markers'
  **genotypes stay bit-identical** through `x = 0.30`, so every pair's kinship over the kept
  set is unchanged to the last bit. The printed count nevertheless falls **46 → 46 → 45 → 43
  → 39 → 37**. Whatever the screen computes, it reads markers it did not keep.

* **And it is deterministic, so "two stages intersected" is out in its noisy-estimates
  form.** Labelling 48 ladder pairs one at a time on the `x = 0.25` map gives **zero
  inversions** in both the whole-map and the kept-subset kinship, with the threshold
  displaced by 0.018 — twelve times the realisation spread of a fixed subset (0.0015).
  A sharp threshold in the wrong place is not two noisy estimates ANDed together.

* **A second necessary condition exists, with no budget in sight.** At `m = 32 768`, where
  the screen is the exact whole-map kinship to 4e-6, a pair cloned across every marker of a
  MAF-0.20 stratum and left untouched on an equally sized MAF-0.45 one is refused outright —
  `No close relatives are inferred.` at kinship **0.20006**, a number KING's own `--kinship`
  prints and agrees with. Same at 0.10/0.13890, at 0.25/0.21731, and at 24 576 @ 0.15 +
  8 192 @ 0.45 / 0.30669. The accept region is flat in the low stratum's clone fraction: past
  the cutoff, extra sharing among uninformative markers buys nothing. It is not a subset
  kinship, not IBS0 under any normalisation (a pair at 0.1047 IBS0/marker and kinship 0.0655
  is **accepted** while one at 0.0602 and 0.20006 is **refused**), not HetHet, not
  contiguity, and not the `Dup/MZ` path. It never binds for uniformly related pairs, which is
  why `bigish` and every Round 1 instrument miss it entirely.

* **In-sample ascertainment is a real component and not the whole story.** It is the one
  loophole in the proof above, and it has the right sign: only a ranking key that leans
  against the pair's own heterozygosity deflates at all. Minor-homozygote count, the simplest
  such key, reproduces the flat-MAF magnitude (model 1.0615 against a measured 1.0654) and
  moves `bigish` from 47 to 41 — against the reference's 36. But it is **4× short** wherever
  the ranking actually resolves, and `R − 1` falls only 2.1× between `n = 110` and `n = 700`
  where a 2-in-`n` selection effect demands 6.4×.

**Two false leads are recorded so nobody re-burns them.** `N = 8·IBS0` and
`φ + 2·IBS0/min_het = 0.5` both look constant across the gate boundary; both are *identities*
for clone-block pairs whose untouched markers sit at MAF 0.45, and neither is a rule.

**What a future maintainer should not attempt.** Do not look for a better marker subset, a
better weighting, a merge, a rank grouping, or a MAF/informativeness threshold — the algebra
rules out the first two for *any* choice, and the measurements above rule out the rest. The
three leads that are still worth an experiment are listed in `22-screen.md` §13: chase the
second condition directly (it is binary, its effect is huge, and one run reads one point);
look for a statistic that degrades with the **discarded** markers' informativeness while the
kept bits are held fixed — a deterministic bound-based early exit over the sorted map has
exactly that shape and would supply the gate for free; and do not fit `R`.

Nothing is landed. Fitting `R` would reproduce `bigish` and nothing else, since `R` swings from
0.998 to 1.085 with the MAF spectrum. The placeholder prefix stays: it reproduces the degree-1
count (18) on every map tried, and swapping it for the whole map would lose that for nothing at
degree 2 (47 against 36). The cost is two cases and one integer on one stdout line; the `.kin0`
row set below that line comes from the exhaustive re-estimate and is byte-correct regardless.

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
a parity case**, and neither is a segment-caller problem. They are retained because the
corpus cannot prove behavior on these shapes. (§4.6 adds the exact-64 safety divergence.)

The one fixture shows both at once. One autosome, 10 000 markers at 10 kb spacing, six
samples in singleton families — usable total `D` = 99 990 000 bp, just under the floor; and
the same fixture with two more markers, `D` = 100 010 000 bp, just over:

| `D` | binary | prints `Segments too short.` | files written |
| --- | --- | --- | --- |
| 99 990 000 | reference | **yes**, 12/12 | `allsegs.txt` |
| 99 990 000 | open-king | **yes** | `allsegs.txt`, **`splitped.txt`** |
| 100 010 000 | reference | no | `.seg`, `allsegs.txt` |
| 100 010 000 | open-king | no | `.seg`, `allsegs.txt`, **`splitped.txt`** |

**1. The `--ibdseg` 100 Mb usable-total floor — fixed.** Below it both binaries print

```text
Segments too short.
  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.
```

and stop before the `IBD segment analysis starts at` block or `.seg` write. The boundary is
closed: a held-out pedigree-bearing fixture at 99,999,999 / 100,000,000 / 100,000,001 bp
matches the reference in normalized console text, exact file set and every output byte. At
the lower point the files are `allsegs.txt` and the independently justified `splitped.txt`;
at and above the boundary `.seg` is added. The comparison is
`tests/parity/probes/segment_floor.py`.

`--ibs` already used the same constant. A new 20-sample cross-check also found that
`--related` treats a non-empty but below-floor map as a kinship-only fallback and that its
small-marker between-family flow differs from the current Rust path; that separate fallback
work remains tracked under §5.12 rather than being hidden inside this `--ibdseg` gate.

**2. Conditional `<prefix>splitped.txt` generation — fixed.** On the *above-floor* run of the same
6-sample fixture, whose families are all singletons, the reference writes **no**
`splitped.txt` and prints no `… is generated for certain pedigree plot applications` line.
Both binaries now do the same. The corpus cannot distinguish the two rules: every dataset
that reaches the segment pass has a family of at least 4 members (`unrelated` 10, `bigish` 9,
`nuclear` 6, `admixed` 4), and the only datasets with no multi-member family at all —
`singleton` and `pair` — sit below the `< 5` sample downgrade and never run the pass. So
`kingsplitped.txt` is byte-identical in all 50 corpus cases and on the focused off-corpus
presence sweep.

The exact rule is now pinned by holding 20 samples and their genotypes fixed while sweeping
maximum family size 1, 2 and 3: size 1 writes and announces nothing; sizes 2 and 3 write and
announce the file, with byte-identical contents. A singleton that names a parent also emits
it, matching the renderer's pre-existing rule. `tests/parity/probes/splitped_presence.py`
checks the three `--ibdseg` shapes and verifies that `--related` never owns this artefact.

### 5.11 Three more behaviors the corpus cannot see

All three were found while validating `<prefix>X.seg` out of sample (§6.1) on built filesets.
**None costs a parity case**, and none is in the X.seg writer — each is in a shared component
that X.seg merely exercises harder than the corpus does.

**1. The `.seg` acceptance gate admits a pair the reference omits.** On a built 8-sample
single-family fileset the reference's `king.seg` has no row for one brother–sister pair while
open-king writes `IBD1Seg 0.0260  IBD2Seg 0.0893  PropIBD 0.1023  3rd`. It is the *acceptance*
gate, not the `--degree` filter: the row is absent from the reference at `--degree 0` as well,
and present in ours at −1, 0, 1, 3 and 4 (at 2 the `PropIBD` cutoff removes it anyway). This
is the first observed counterexample to `inf1/inf2 >= 10` — the corpus has **0 extra and 0
missing** `.seg` rows — and it propagates into `X.seg`, which mirrors `.seg`'s rows.
*Reproduce:* `probes/xseg_probe.py -v`, case `onefam/x1000#0` — an 8-member nuclear family
(3 sons, 3 daughters) with 1 000 X markers. It is the only default-floor case in that probe's
1 040 runs whose `X.seg` differs, and it differs because `king.seg` does.

**This is not §5.14's marker-count table, and saying so is the point.** The `>= 10` this
counterexample contradicts is now known to be the *first row* of a table the reference reads
off the total marker count, and the row below it — 20, from 400 000 markers — closes a
second and much larger counterexample. It does not touch this one. The fileset here carries
5 000 markers in all (two autosomes of 2 000 apiece plus the 1 000 X markers), which is that
same first row, so the gate over it is 10 both before and after that fix and open-king's
answer on this case is unchanged to the byte. Issue #11 item 1 stands, and what it still
needs is a discriminator that lives **inside** the first tier.

**2. `.fam` SEX fields outside {0,1,2} — fixed.** Measured by sweeping 43 spellings through
`X.seg`'s raw `Sex` column:

| reference prints | `.fam` field |
| --- | --- |
| `0` | `0` `00` `0.0` `-0` `-9` `-9.0` `x` `?` `NA` `na` `b2` |
| `1` | `1` `-1` `-2` `3` `9` `10` `12` `02` `002` `0002` `007` `+2` `1.9` `1e0` `M` `m` `male` `MALE` |
| `2` | `2` `20` `21` `22` `2.5` `2.9` `2x` `20x` `2e0` `F` `f` `female` `FEMALE` |

i.e. a leading `2`/`F`/`f` is female, a leading `M`/`m` is male, an otherwise-numeric field is
male unless it evaluates to `0` or `-9`, and anything else is unknown. `parse_sex` now
implements and unit-tests every measured spelling. Every corpus `.fam` uses only `0`/`1`/`2`,
so the focused test is the regression guard for this rule.

**3. The reference reads past the end of a marker array whose length is an exact multiple of
64.** On such an array it adds an *absolute* coordinate — `pos[last]`, not a difference — to a
pair's IBD2 total, so `IBD2Seg` comes out larger than 1. Swept at 320/384/448 markers (anomaly
on some pairs, always exactly `+pos[last]/D`) against 319/321/383/385/449 (no anomaly on any
pair, any seed). It shows identically in `X.kin` and `X.seg`, so it is the shared caller, and
it is unreachable from the corpus, whose arrays all end mid-word. Not emulated: it is an
uninitialised read, not a rule.


### 5.12 Three divergences found while writing the user documentation

Measured on 2026-08-18 against the same reference build, on filesets derived from the corpus
by the recipes given below, and reproduced from a cold `cargo build --release`. **None costs a
parity case** — every one of them needs an input shape the 480 captures do not contain — and
all three are visible to a user, which is why they are recorded rather than left in a
transcript. §5.10 and §5.11 remain the older, independently measured set.

**1. The "no informative IBD segments" fallback and shared screen are fixed.** Any panel
too sparse for the segment caller lands on it, which on a 200-sample fileset means roughly
12,500 markers genome-wide.

*Reproduce* (`reshape.py` is the thinning script published in
[`INTERPRETING.md`](INTERPRETING.md#appendix--reshapepy)):

```text
python3 reshape.py /tmp/kingdocs/bigish thin4 --every 4     # 200 samples x 12 500 SNPs
king -b thin4.bed --related --cpus 1
```

| | reference | open-king |
| --- | --- | --- |
| console | `No informative IBD segments.` + `Relationship inference will be based on kinship estimation only.` | **matches**; 15 candidates (3 FS, 12 second-degree) |
| `.kin` | **12 columns** (`… Kinship Error`) | **byte-identical**, all 574 lines |
| relationship summary | 436 by inference | **matches**, 436 by inference |
| `.kin0` rows | 15 | **byte-identical, all 15 rows** |

The reference detects that the map yields no usable segment and falls back to pure
kinship-based inference — a documented mode of KING, and the short `.kin` layout is its
signature. open-king now takes that path for `--related`, including the short headers,
kinship-based within-family error/summary rules and the screened between-family row shape.
The shared `ConvertLGtoSLG` score, unstable QuickIndex ordering, progressive degree gates,
and packed-tail behavior now reproduce the reference's exact 15-candidate row set.
The same fallback now runs through `--unrelated`, `--cluster` and `--build`. On the held-out
fixture, both unrelated selection files are byte-identical; `--cluster` writes the exact
`updateids.txt` and correctly omits its segment-only `cluster.kin`; and `--build` writes
byte-identical `build.log`, `updateids.txt` and `updateparents.txt`. The console's inferred
relationship split is exact: 3 full-sibling and 12 second-degree pairs. Both binaries print
`Cutoff value for IBS0 between FS and PO is set at 0.0050`. `tests/parity/probes/sparse_fallback.py`
checks all five analysis paths and every artifact. On `--ibdseg`, both binaries print the
PLINK sorting note after `No informative IBD segments.` and write only `splitped.txt`.

A second held-out shape pins the adjacent non-empty case: at 99,999,999 bp of usable map,
both binaries print `Segments too short.` and take the same short kinship-only `--related`
layout.

**2. Unsorted `.bim` validation — fixed.** Two shapes, both derived from `multifam` by the
`fixtures.py` script published in [`CLI.md`](CLI.md#10-the-derived-filesets-used-above):

| map | both binaries |
| --- | --- |
| positions descending inside each chromosome | `Positions unsorted: rs1_1009689 at 65904473, rs1_1055261 at 65851170.` + the PLINK note; writes only `splitped.txt` under `--ibdseg` |
| chromosomes 22 → 1, positions ascending inside each | `Chromosomes unsorted: rs22_14205438 on chr 22, rs21_1002722 on chr 21.` + the PLINK note; writes only `splitped.txt` under `--ibdseg` |

`tests/parity/probes/map_order.py` checks both shapes through `--related`, `--ibs`,
`--unrelated`, `--build`, `--bysample`, `--bySNP`, `--cluster` and `--ibdseg`: all 16
normalized console streams, file sets and file bytes match the reference.

**3. Sample IDs colliding only in case — fixed.** The reference folds ASCII case when checking
`(FID, IID)` uniqueness — established independently in
[`BEHAVIOR.md`](BEHAVIOR.md#q6--the-sample-id-sort-comparator), which records `{A, a}` and
`{ab, aB}` being rejected at load. `open-king-io` now canonicalises only the identity key while
retaining the original spelling for output, so exact duplicates and case-only FID/IID
collisions stop at the same pedigree-validation point.

*Reproduce:*

```text
awk 'BEGIN{OFS=" "} {if ($1=="FAM1" && $2=="A_M") $2="a_f"; print}' \
    /tmp/kingdocs/multifam.fam > case.fam
king -b /tmp/kingdocs/multifam.bed --fam case.fam --kinship    # both: exit 1, same duplicate diagnostic
```

**A1-major inputs are now rejected at the same stable boundary.** Black-box sweeps identified
the first 4,096 retained autosomal markers as the window, a strict ten-percent cutoff
(409 passes, 410 aborts), and the affected analysis/sample-size surface. The fatal console,
exit status and complete pre-fatal artifact set match in `tests/parity/probes/a1_major.py`.
`--kinship`, `--duplicate`, `--autoQC` and KING's disabled/downgraded small-sample paths remain
exempt. For maps shorter than 4,096 markers the reference reads unstable tail state and can
abort valid data nondeterministically; open-king deliberately skips that unsafe check.
[`CLI.md` §3](CLI.md#two-hard-requirements-that-are-easy-to-miss) states the input contract.

### 5.13 A second build of KING 2.3.2 agrees on every output file but one, and that one is `--noscreen`

§5.3 records that the reference prints uninitialised memory in its own banner, and notes the
value "is stable within one build and environment". That sentence understates the
consequence, so here is the cross-build measurement it implies.

KING 2.3.2 was compiled a second time from the published `KINGcode.tar.gz`
(sha256 `b6c636ac…`) with a different compiler from the one that produced the capture
binary: Homebrew GCC 16, `g++-16 -lm -lz -O2 -fopenmp`, arm64. Call it **build B**; the
binary the goldens were captured from is **build A**. The source was compiled and then
deleted without being read, so §1's clean-room rule is intact.

Two replays were run, and they answer different questions.

**The graded 480, through the project's own harness.**
`run_parity.py --impl <build B>` is the harness's reference-vs-reference self-check, and it
reports `5 PASS, 475 FAIL, 480 total (876 output file(s) byte-compared, 8 diff-excluded)`.
Every one of the 475 failures is `stdout!=`. **Not one is an output-file difference**, over
all 876 files. So the entire disagreement between the two builds, on the graded corpus, is
on standard output. Two lines carry it, measured on
`--related -b multifam.bed`: build A prints `--noscreen [-1717986816],` where build B prints
`--noscreen,`, which reflows the two-line `Inference Parameter` block and defeats the
harness's normalization of the integer; and the `N CPU cores are used...` line reports the
host's core count, which is not a build property at all.

**All 490 captures, compared directly.** Widening past the graded set and comparing output
files byte for byte outside the harness: **488 of 490 cases produce byte-identical output
files.** Six of those are `_analysis` captures that store `MD5SUMS.txt` instead of the files
themselves; their checksums were compared directly and all match. The two exceptions are:

| case | difference | what it is |
| --- | --- | --- |
| `core/sexchr__kinship` | `kingX.kin0` differs | the documented threads>1 race on the X between-family writer, §5.2. Not a build difference; the file differs between two runs of the same binary. |
| `core/_analysis/multifam__related_degree1__noscreen` | build B writes a `king.kin0` that build A does not | a real behavioural difference between the two builds. |

The second one is worth stating plainly, because it is the one option whose banner value is
undefined:

* **Build A** on `--related --degree 1 --noscreen`: writes `king.kin` and `kingallsegs.txt`,
  no `king.kin0`, and stdout ends `No close relatives are inferred.`
* **Build B**, same invocation: additionally writes `king.kin0` with eight between-family
  pairs, including `FS` and `PO` calls.
* **open-king** matches build A.

So `--noscreen` bypasses the two-stage screen in build B and does not in build A. `--noscreen`
takes an integer, and §5.3 establishes that the integer it carries is uninitialised memory. An
undefined value reaching a branch is exactly how the same source yields two behaviours, so
this is best read as one symptom rather than two: the option's value is undefined, and
therefore so is its effect.

Two consequences for this project:

1. **open-king treating `--noscreen` as inert is not a gap.** It reproduces build A, which is
   the only build the corpus can speak for. Making it bypass the screen would reproduce
   build B and fail 480 cases. There is no implementation that satisfies both, because the
   reference does not agree with itself.
2. **§5.7's screening derivation must not be pushed further using a rebuilt reference.** A
   probe compiled locally may answer `--noscreen` questions differently from the binary the
   goldens came from, and would look like a discovery rather than a build artifact.

The wider reading is the encouraging one. KING's *computed* output is reproducible across
independent compilations: every kinship, IBS, segment and QC file in the corpus is byte-equal
between build A and build B. What is not reproducible is its banner and one option whose
value was never initialised. That strengthens the goldens as a target and narrows the
standing "one build" caveat at the end of §5.0 to the places where it actually bites.

### 5.14 `--ibdseg`'s segment-acceptance parameters are a table keyed on the marker count

**Found by counterexample on real data, fixed in v0.1.1, and invisible to all 480 captured
cases** — every corpus dataset is at most 50 000 markers, and the smallest fileset that can
see this is 400 000. Reported as
[#13](https://github.com/Broccolito/open-king/issues/13) and
[#14](https://github.com/Broccolito/open-king/issues/14).

The instrument is a 663 197-marker × 157-sample autosomal panel — an LD-pruned common-SNV
grid, markers about 4 300 bp apart, an order of magnitude more markers and an order of
magnitude denser than anything in the corpus or the fixture rigs. Against KING 2.3.2 on it,
`--kinship` was already exact (770 `.kin0` rows and 11 476 `.kin` rows, byte-identical), and
`--ibdseg` was wrong on **every one** of the 7 pairs it reports: `IBD1Seg` and `PropIBD` on
all seven and `InfType` on one, which read `3rd` against the reference's `4th`. Two separate
faults, and they had to be fixed together.

**1. The informativeness gate is a table, not a constant.** The reference chooses it once
per run from the total marker count:

| markers | gate | minimum candidate length |
| ---: | ---: | ---: |
| `< 400 000` | 10 | 400 000 bp |
| `400 000 … 2 000 000` | **20** | 400 000 bp |
| `> 2 000 000` | 100 | 100 000 bp |

`docs/research/13-informativeness-gate.md`'s measurement of 10 is not withdrawn: it is
exact, it is validated out of sample on 1 170 pairs with no overlap, and every fixture and
dataset behind it is in the first row. What was wrong was promoting a first-row measurement
to a universal constant. That the trigger is the **count** and nothing about the markers
took 283 controlled reference runs and a single-marker bisection: appending one all-missing
marker — uncallable for any pair and uncountable by any gate — to take a fileset from
399 999 markers to 400 000 moves 27 Mb of called IBD, and padding the count alone with the
span and the usable-segment set held fixed reproduces the step. Both boundaries are bisected
to the marker, at exactly 400 000 and again at exactly 2 000 001. The second column is
recorded but not implemented: this caller has no separate candidate-length gate to key on,
the first two rows agree on 400 000 bp, and nothing here reaches the third row.

**2. The IBD1 merge cap of two unusable words is a first-tier reading.**
`20-seglength-floor.md` §3 bisected it as absolute — three unusable words merge at no floor,
however little they carry and however short the gap — on a rig whose markers sit 20 000 bp
apart, and `21-push-merge.md` §4 reproduced it. It is a genuine count and not a distance in
disguise; `23-gap-bound.md` §5 finds a corpus merge across two unusable words and
9 652 629 bp. But on this panel the reference joins across **four**, which no cap of two and
no cap of 129 marker intervals permits — so `20-…` §11 item 3, "two words or 129 marker
intervals?", is answered *neither*. The cap now applies only to the first tier.

**Which of the two things that move together — the marker count or the spacing a dense
panel comes with — the merge cap is really keyed on is not separated.** Two observations at
opposite ends of both cannot tell them apart. The tier is used because it is the boundary
that is independently bisected here, and because keying it there leaves every first-tier
measurement reproduced exactly as measured. A fixture holding the marker count inside one
tier and moving only the spacing would settle it.

**Why they had to land together.** On this panel, correcting the merge cap alone makes the
result *worse*: an eighth pair appears that the reference does not report. The cap was
masking the gate.

**Verification, all of it.**

| check | result |
| --- | --- |
| `ibd.seg` vs KING 2.3.2 on the panel | **byte-identical** (`cmp`), 7 pairs, all four columns |
| `ibdallsegs.txt`, `ibdsplitped.txt` | byte-identical |
| `--kinship` on the same panel | `kin.kin` 11 476 rows and `kin.kin0` 770 rows, byte-identical — unregressed |
| captured parity corpus | **480 PASS / 0 FAIL**, 876 files byte-compared, baseline `MATCH` — not one case moves |
| determinism | `--ibdseg` run three times, byte-identical each time and to the reference |

The corpus not moving is a result and not a formality: it is the check that the fix is a
*table* rather than a new constant. Had any of the 480 changed, the first row would have
been wrong.

---

## 6. Structural analyses and held-out residuals

### 6.1 `<prefix>X.seg` — **implemented, both cases pass**

Rules, gate and evidence: `crate::analysis::xseg`. What was measured against the reference,
including the corrections this section previously got wrong:

* **When it is written.** `--ibdseg` writes `<prefix>X.seg` exactly when `--degree` is
  non-zero *and* the fileset has usable X-chromosome segments. Bare `--ibdseg`,
  `--ibdseg --degree 0` and `--ibdseg --seglength 5` write none; `--degree` (bare, = 1),
  `--degree 1`, `--degree 2`, `--degree 3` and `--degree -1` all do. Of the 13 corpus datasets
  only `sexchr` has usable X segments, so only 2 of the 480 cases were affected; **both now
  pass**. Two things the corpus could not decide, settled on built filesets: there is **no
  512-marker threshold** (that belongs to `--kinship`'s X pass alone — 320 X markers over
  30 Mb write `X.seg`, 319 do not, and 640 markers inside 5 Mb do not), and there is **no
  family-count condition** (a single-family fileset writes it).
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
* **Its `PropIBD` follows `.kin`'s rule, not `.seg`'s** — this section previously said the
  opposite, and building it both ways settled it: `seg_prop_ibd` fails both captures
  (`kingX.seg!=(num)`), full-precision `IBD2Seg + IBD1Seg/2` passes both. The integer test of
  §4.3 does not discriminate here because it treats an exact decimal tie as compatible with
  either direction; the actual double arithmetic does. Two rows refute the `.seg` rule
  outright: `S_SON1`/`S_DAU1` at `IBD1Seg 0.4257` prints `0.2128` where `4257 * 5e-5` lands
  just above the half and renders `0.2129`, and `S_SON1`/`S_DAU2` at `0.9067` prints `0.4533`
  against `0.4534`. So `X.seg` and `X.kin` carry the same three numbers for the same pair.

* **There is no X caller.** The autosomal segment caller runs unchanged over the X marker
  array and the X bit planes; there is no male/female branch, and none is needed, because a
  hemizygous male is stored **homozygous** and so is never heterozygous and never IBS0
  against a heterozygote. That is why a father–son pair scores 0, a father–daughter pair
  `IBD1Seg 1.0000`, and brothers who drew the same maternal X score `IBD2Seg`.
* **Samples of unknown sex are not excluded**, unlike in `--kinship`'s `X.kin`; their rows
  appear with the `.fam` code printed raw.

**Held out** — `python3 tests/parity/probes/xseg_probe.py --impl target/release/open-king`, which
prints `presence 1040/1040`, `default 875/880`, `default_given_autosome_ok 625/625`,
`nondefault 140/160` (105/160 before the push and IBD2-merge corrections landed). The two
captures are 28 rows of one 6-sample family, so nothing here was
fitted to them. 1 040 reference-vs-open-king runs were built from fresh pedigrees
(4 to 48 samples, one to six families, three-generation, unknown-sex, all-singleton, and one
below the `< 5` downgrade) crossed with five X maps (one, two and three usable segments; 333
to 1 500 markers; X alone and X + Y + XY + MT) and thirteen flag combinations: the **emission
gate agrees on 1 040 of 1 040**, and at the default floor `X.seg` is byte-identical on
**625 of 625** runs whose autosomal `.seg` is also byte-identical. Every remaining difference
is inherited, and every one of the 20 is a `--seglength 10` run: the residual of §5.0 (none
at `--seglength 5` any more), or the one `.seg` row the autosomal
acceptance gate disagrees about (§5.11).

### 6.2 `--build` on `bigish`

`apps/bigish__build` is now **byte-identical across stdout and all four compared output
files**. Its `kingbuild.log` contains all 18 reference lines (806 bytes), including the
inference block and blank-line layout. The other 12 captured `--build` datasets remain
byte-identical.

The implementation now includes FS2, PO.S orientation and synthetic-parent consumption,
AV.FS/AV.HS/HS.UN2, reported-segment `Join3/Join2`, the verdict dead band, and the exact
internal sibling order. That last order is not a hash table: the KING 2.3.2 binary uses
libStatGen's unstable `QuickSort` first by `(FID, IID)` and then by
`(FID, father, mother)`, leaving a deterministic swap residue among equal-parent siblings.

The cached research replay is **277 / 347** byte-exact. Another **52** logs contain the
same distinct lines and differ only in repetition/count residue from the reference's
unstable inference loop. Of the remaining 18, two reference logs were truncated during a
debugger probe and two require the still-open cross-FID `<FID>-><IID>` rename.
The remaining rare constructed-pedigree residuals stay recorded here; they do not affect
the now-passing 480-case corpus claim.

The remainder of this section records the derivation history. Intermediate scores and
open-question wording below are superseded by the current status above, but are retained
as evidence for the triggers and counterexamples that led to the implementation.

`kingupdateparents.txt` **is written and byte-identical**, and so is `kingupdateids.txt`.
The two are in **different orders**, which is a thing this section previously had wrong
because `bigish` cannot show it: `updateparents.txt` is in cluster order (`KING1`'s rows,
then `KING2`'s), while `updateids.txt` is in original-`(FID, IID)` order, so a fileset
whose clusters are `KING1 = Z1+Z2`, `KING2 = M1+M2`, `KING3 = A1+A2` writes the `A` rows
first. On every corpus fileset the two orders coincide.

**The cluster-numbering bug recorded here is fixed.** Merged clusters are not numbered in
family order and not by the relationship type of the joining pair either, though the two
agree on the three shapes that first exposed it. They are numbered in the order a
**staged merge queue** creates them:

* the qualifying cross-family pairs are worked through **by relationship type** —
  `Dup/MZ`, then `PO`, then `FS`, then anything weaker — the scan order surviving inside
  each type;
* a cluster is created, and takes its `KING<k>`, the first time the queue joins two
  families that are not already together;
* a cluster's `OriginalFamID` list is in **absorption** order, not file order: a cluster
  whose `Dup/MZ` edge is `QBB–QBC` and whose `FS` edge is `QBA–QBB` prints
  `QBB,QBC,QBA`, and that list order is the clearest evidence that the queue is staged
  rather than sorted afterwards.

`docs/research/fixtures/clusternum.py` builds nineteen shapes that discriminate this from
family order, cluster size, and the largest kinship of the joining pair: the queue rule is
**19 of 19**, family order 7, size 7, kinship 11. `seeds` re-runs one two-type shape over
eight fresh seeds and **4 of 8** contradict the kinship ordering outright — the `FS` pair
scores φ = 0.30 against the `PO` pair's 0.25 and the `PO` cluster is still `KING1`. The
corpus is indifferent to all of it (every merge in `bigish` is `FS`, and no other corpus
fileset merges), and the harness stayed at 477 across the change.

**The same rig found that our clustering *gate* was wrong**, and that is fixed too. It was
`kinship > 2^-2.5`; the reference admits a pair on the disjunction `--related` uses for a
`.kin0` row at `--degree d`, `kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)`. A 3/4-sib
pair at `kinship 0.1749` — *under* `2^-2.5` — merges its two families on `PropIBD 0.3646`
alone, and `clusternum.py gate` scores the disjunction 19 of 19 against the kinship rule's
18. The cut also follows `--degree`: a fileset whose only cross-family link is a half-sib
pair reports `No families were found to be connected.` at `--degree 1` and merges at
`--degree 2`. No corpus capture moves, because no corpus fileset has a cross-family
2nd-degree pair outside a cluster that merges anyway.

Note the two gates are **not** the same one. Reconstruction keeps the plain 1st-degree band
edge: the 3/4-sib cluster above merges and then reconstructs nothing at all, raising no
header and no `RULE FS0`, which is why `unrelated::InfTypes` exposes `merging` and
`first_degree` separately.

**Which log lines are in which half.** This is the other thing the round corrected. The log
splits into a *rule* half (written) and an *inference* half (not), and the split is not the
one the indentation suggests:

| template | half | trigger |
| --- | --- | --- |
| `Family KING<k>:` | — | once, before the cluster's first line |
| `Duplicate <a> (of <b>) is removed.` | **rule** | an inferred `Dup/MZ` pair, if the cluster raises something else |
| `RULE FS0` | **rule** | a component of *inferred FS* ∪ *declares the same couple* that the inference created |
| `RULE FS1` | **rule** | one more member joining a component that already had a sibship |
| `RULE FS2` | rule | two declared sibships in one component — never observed |
| `Reconstruct parent-offspring pair (X, Y)...` | **inference** | an inferred `PO` pair, in a cluster whose inference block also speaks |
| `…'s sibship is used to determine…`, `RULE PO.S`, `<n> is created as …'s mother.` | inference | that `PO` pair, when a sibship orients it |
| `INFERENCE AV.FS` | inference | an `R` inferred 2nd-degree to **both** named members of a sibship |
| `INFERENCE AV.HS`, `HS <a> unrelated to <b>`, `INFERENCE HS.UN2` | inference | a half-sib pair the avuncular pass turned up |

`Reconstruct parent-offspring pair` was assumed to be a rule line and is not: **42 of 42**
clusters that print it also print an `INFERENCE` line, and a `PO` merge between two
families with no sibship anywhere — two one-person families and two childless couples,
three seeds each — prints nothing at all, not even a header. Writing it would have been a
guess that happens to fit two shapes and fails at least nine.

`Duplicate … is removed.` is the opposite and is now written. `dupkeep.py` scores it over
ten shapes × three seeds: **23 of 27** runs print it in a file with no `INFERENCE` line, so
it is rule-half. Which copy goes is measured the same way — the reference keeps the copy
with more **declared 1st-degree relatives that the fileset carries** (named parents, full
sibs naming the same couple, children naming it) and breaks ties on the ID comparator,
keeping the later id: **27 of 27**, against 21 for "keep the later id" and 6 for "keep the
earlier". That also corrects the old clause "a cluster holding an inferred duplicate
contributes no rows at all" — it contributes them whenever removing the duplicate leaves an
`FS` or `PO` pair behind.

**At that stage, the rule half scored 53 of 59 held-out shapes**, byte-identical, up from 23 of 30
(`buildlog.py rules`, over `build_shapes.py`'s twenty, `avfs.py`'s ten, `clusternum.py`'s
nineteen and `dupkeep.py`'s ten). The six that differ are two `<FID>-><IID>` renaming shapes
that were not yet reproduced, three that differ only in the *order* a sibship's members are listed,
and one where the then-unimplemented `PO.S` branch consumes a synthetic id so the next sibship
takes `(4 5)` where we write `(3 4)`. `build_shapes.py` is **18 of 18** in scope on
`updateparents.txt` and the console tail, up from 15.

**At that stage, the blank lines were not written**, because their count is a function of the
inference half. Two rules fit; the one this section carried is the weaker of them:

* **block** — one blank before each sibship's block until the family prints its first
  inference, and one per block if it never does;
* **reject** — one blank opens the section, and one more per candidate `R` *examined and
  turned down* before the first line prints.

`buildlog.py blanks` scores both at **107 of 113** clusters, on different failure sets, with
a scorer that has to guess the block order and the candidate order. What separates them by
hand are the two clusters `block` provably misses: `three_fs`, whose first sibship faces two
candidate uncles and prints **three** blanks where `block` says two, and `ord3`, whose two
sibships face no candidate at all and prints **one** where `block` says two. `reject` gets
both, and reproduces 1, 1, 2 for `bigish`, 3, 2, 1 for `three_clusters` and 1, 2 for
`mixed_po_fs`.

**At that stage, the sibship member order was the open question, narrowed as follows.** The
`RULE FS1` line prints a sibship's members in an internal order — `(A_C2 A_C3 A_C1)` for
`A_C1..A_C3` — and the same order picks the pair an `AV.FS` line names. Earlier rounds ruled
out genotypes, `.fam` row order, absolute sample index, sibship size and position, and every
pairwise statistic. `docs/research/fixtures/siborder.py` closes the remaining space:

* it **is** a function of the members' **id strings** — three distinct position-orders over
  eight id sets in one fixed pedigree, `A_C1..A_C3` giving `(2,3,1)` where `K1..K3` and
  `1001..1003` give the identity and `zeta,alpha,mu` gives `(3,2,1)`;
* it is **not a per-id ranking**: over thirteen subsets of one eight-id pool the pairwise
  precedences contradict each other **91** times, so no `sort by f(id)` can reproduce them —
  the order moves when the *set* changes, which is what a hash table's capacity does;
* the container is scoped to the **individual's own family**: renaming the sibship's parents
  changes the kids' order, while the joiner's id, the other family's ids, the padding
  vocabulary (four) and the total sample count (102 … 142, nine values) all leave it
  byte-identical;
* four `.fam` permutations of the same three ids give one order, confirming it is keyed by
  the strings and not by position.

So it is an iteration order over a family-scoped, id-keyed, capacity-sensitive container —
a hash table — and reproducing it means identifying the hash. That is the next piece of work
on this case, and it is worth exactly three of the 59 shapes: `bigish`'s log has no `FS1`
line at all, so the corpus case does not turn on it.

**The verdict cut is tightened** from (0.846, 0.902) to **(0.8495, 0.9005)**. `buildlog.py
cut` re-reads every `AV` line every rig has produced — **259** of them, 133 `uncle|aunt`
against 126 ambiguous, up from 53 — and the largest `uncle` prints `0.850` while the
smallest ambiguous prints `0.900`. At `%.3lf` those stand for `[0.8495, 0.8505)` and
`[0.8995, 0.9005)`, so all three of 0.85, 0.875 and 0.9 still survive; the bracket is now
symmetric about 0.875 and one printed step from closing.

**That statistic is no longer unknown, and it confirms the case is blocked on the segment
caller.** Writing `IBD(x, y)` for the union of a pair's called IBD1 and IBD2 segments as a
set of base pairs, for a triple `(R; N1, N2)`:

```text
Join2 = | IBD(R,N1) ∩ IBD(R,N2) |
Join3 = | IBD(R,N1) ∩ IBD(R,N2) ∩ IBD(N1,N2) |
```

Where `R` is IBD to both sibs, a grandparent forces them to have inherited the same
parental haplotype (ratio → 1); an avuncular does not (→ 2/3), which is exactly what the
reference's two message variants discriminate. Re-scored against **34 `AV.FS` values the
reference emitted over 16 filesets** — `bigish` plus 15 held-out pedigrees with sibships of
2…6 — the formula is one-sided high on every one: mean **+0.0039**, range
**+0.0003 … +0.0102**, and **only 1 of 34 rounds to the printed three decimals**, none of
`bigish`'s five among them. So even a complete reconstruction implementation leaves all five
of `bigish`'s log lines wrong.

**And with the shipped binary that is now a hard result rather than a deferral.** Re-run for
this release on the default sweep — 14 triples, `bigish`'s five among them — the residual is
mean **+0.0040**, range **+0.0008 … +0.0084**, and **0 of 14** round to the printed three
decimals. On **10 of the 14**, every one of the three pair totals behind the triple matches
the reference's own printed `IBD1Seg + IBD2Seg` **exactly** (`dU = 0.0000`), `bigish`'s five
included — and the residual on those ten is unchanged. Byte-identical totals do not imply
identical segments. `Join3/Join2` intersects three segment *sets*, so it is the only
statistic in KING's entire output that reads **placement**, and it says our placement is not
the reference's even where every printed number agrees. (The other four triples, on two
held-out `4:4` shapes, are also the only out-of-corpus rows in that sweep where our
`.seg` totals themselves differ — the same phenomenon §4.6 measures.)

**A retracted argument, and what replaces it.** An earlier revision of this section claimed
the residual was *entirely* accounted for by our sib-pair union over-call `ΔS`: since `ΔS`
can raise `Join3` by at most `ΔS`, and the two `R`-to-nephew sets are avuncular (reference
`IBD2Seg 0.0000`, our reported union exact on **823 of 823** such corpus rows), the ratio
should be high by at most `ΔS / Join2` — quoted as *39 of 39 triples inside*
`[0, ΔS / Join2]`.

**Re-run on fresh shapes and unused seeds, that is 11 of 34.** And on the **8 triples where
all three pair unions match the reference's reported `IBD1Seg + IBD2Seg` to four decimals**
the residual is still `+0.0003 … +0.0049`, mean `+0.0031`. The bound was never evidence:
`dU` is a *rounded total*, whereas `Join3/Join2` intersects three sets and therefore reads
segment **placement**, about which a matching total says nothing. The argument is withdrawn.

The conclusion survives on its own footing: **no variant of the statistic removes the
residual**, so the fault is in the inputs, not the arithmetic. Swept over the same 34
triples —

| variant | exact / 34 | mean residual |
| --- | ---: | ---: |
| as defined above (base pairs, refined endpoints, `IBD1 ∪ IBD2`) | 1 | +0.0039 |
| marker counts instead of base pairs | 1 | +0.0039 |
| minimum piece length on `Join3` (0.1…3 Mb) | ≤ 4 | +0.0015 … +0.0039 |
| minimum piece length on both (0.25…5 Mb) | ≤ 4 | +0.0014 … +0.0038 |
| eroding every set by 1…63 markers | ≤ 10 | +0.0034 … −0.0287 |
| word-aligned instead of refined endpoints | — | −0.025 |
| re-calling at 0…10 Mb minimum segment length | — | no change below 5 Mb, worse at 10 |

— nothing approaches exactness, and the residual is *heteroscedastic*, spanning thirtyfold
across triples, which is the signature of data-dependent caller error rather than of a
constant the formula is missing. (Eroding by 6 markers zeroes the mean at 10 of 34 exact.
It is an unprincipled knob and is not landed.)

**What that means for this case, stated plainly.** The claim an earlier revision made —
"`apps/bigish__build` closes when §4.1 closes" — is **refuted**, not merely retracted. §4.1
has closed: every `.seg` row in the corpus is byte-exact at every floor. The `AV.FS`
residual did not move. So `--build`'s `INFERENCE` half is blocked on something the corpus
cannot show and no output file exposes: the base-pair **placement** of called segments,
plus the unidentified sib-pair ordering below. Anyone resuming it should start there and
not with the `.seg` columns.

Two further rules were measured the same way, and one is a sharp negative:

* **Which `R` a line can be raised for** — `R` must be an inferred **2nd-degree** relative
  of *both* named members, and nothing weaker reproduces the candidate set: the
  three-father shape names exactly the third family's two children against the father
  sibship, the four-father one exactly the six children of the third and fourth families,
  and each family's own children are excluded for being 1st-degree to their own father. One
  `R` may print the same line two to four times; that repeat count is per `(R, sibship)`,
  is not the number of sib pairs, and is **not identified**.
* **The named sib pair belongs to the sibship, not to `R`** — every `AV.FS` line raised
  against one sibship names the same pair whatever `R` is (four sibships against three
  distinct `R` each, two more against two, including cases whose verdicts differ), and
  where the sibship is the `RULE FS0`/`FS1` one it is that sibship's first two members in
  the order the rule line prints them. So the `AV.FS` pair and the `FS1` member list are
  one ordering, and it is **unidentified**. Its candidate space, however, is now closed by
  measurement rather than open:
  - **not genotype-derived** — four fresh seeds at each of three sibship sizes give
    byte-identical `FS1` orders (`C2 C3 C1`, `C3 C4 C2 C1`, `C4 C1 C5 C3 C2`) while the
    sibship's own kinships move over a 0.10 range (`buildlog.py order`);
  - **not the `.fam` row order** — permuting a sibship's rows, genotypes with them, leaves
    the named pair on the same two individuals at new positions;
  - **not the absolute sample index** — moving all 80 padding singletons of a four-family
    fixture to the front leaves the whole log byte-identical;
  - **not the sibship's size or position** — four three-child sibships in one cluster print
    four different orders, and `bigish`'s structurally identical `B01` and `B13` name
    `(C2, C3)` and `(C2, C1)`;
  - **not any pairwise statistic** — over the 27 measured sibships of ≥ 3 children, no
    `argmin`/`argmax` of `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`, `IBD1Seg`,
    `IBD2Seg`, `PropIBD`, `N_SNP`, `Z0` or `Error` picks the named pair more than **11 of
    27** times against a 1-in-3-or-worse baseline, and nor does any of ten segment-level
    statistics computed here (best 6 of 20). The pair's rank on the ratio itself runs from
    first to last.

  The earlier reading — "a function of the pedigree shape alone", with a positional map —
  is **withdrawn**: it was measured only on the first family of two-family fixtures, where
  the answer happens to be constant, and the four-sibship fixture refutes it outright.
  `buildlog.py order` and `buildlog.py pairs` re-measure all of it.
* **The verdict is a cut on the ratio**, `uncle|aunt` below against
  `grandfather|grandmother, HS, or nephew|niece` above, all three word pairs following
  `R`'s sex. Bracketed at the top of this section to **(0.8495, 0.9005)** over 259 values,
  which still does not separate 0.85, 0.875 and 0.9.

`updateparents.txt` and `updateids.txt` are both written and byte-identical, and the rule
half of the log with them — none of which can flip the case, but all of which generalise.
The `build.log` derivation lives in `crates/open-king-cli/src/analysis/build.rs`'s module doc;
the rigs are `docs/research/fixtures/avfs.py` (held-out pedigree shapes),
`avfs_score.py` (the `Join3/Join2` scorecard, about twenty seconds),
`build_shapes.py` (twenty merge shapes and the `updateparents.txt` scorecard),
`clusternum.py` (nineteen shapes for the merge queue and the merge gate: `score`, `seeds`,
`gate`, `dump`), `dupkeep.py` (ten shapes for the duplicate rule), `siborder.py` (the
sibship-order container: `names`, `popn`, `setsize`, `perm`, `subsets`, `family`, `sizes`)
and `buildlog.py` (`rules`, `blanks`, `cut`, `order`, `pairs`).

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

**Nothing else is normalized, and that matters more than it sounds.** Row order is compared,
not sorted; trailing whitespace is compared; column widths are compared. Both rules of §4.3
and §4.5 were found only because this comparison is a byte diff — every row-level grader in
the project matches rows on their identifier columns first, and had been reporting `0 extra,
0 missing` on a file whose rows were in the wrong order for the entire life of the project.

**The regression baseline.** `tests/parity/BASELINE.txt` records the outcome of all 480
cases *with their per-file notes*, and `run_parity.py --baseline` fails on any difference in
either direction — a new failure, a changed reason for an existing failure, or an
**unrecorded improvement**. That last one is deliberate: a change that fixes a case must
also update the file that says it was broken. CI runs it on every push. Regenerate with
`--baseline --write-baseline` and commit the diff alongside the change that earned it.

For this release the baseline diff is exactly three lines, all in the same direction:

```
-FAIL ibdseg/bigish__ibdseg_seglength5   | king.seg!=(num)   ->  PASS
-FAIL ibdseg/missing__ibdseg_seglength5  | king.seg!=(num)   ->  PASS
-FAIL ibdseg/missing__ibdseg_seglength10 | king.seg!=(num)   ->  PASS
```

The other 477 rows — status *and* per-file note — are byte-identical to the previous
baseline, which is the harness's own statement that nothing regressed.

**A note on checkouts.** The goldens are compared byte for byte and 486 of them contain bare
`CR` characters (the reference's own progress tokens, mid-line rather than line-terminating).
`.gitattributes` marks the whole golden tree, `BASELINE.txt` and the fixture caches `-text`
so that no platform's git rewrites them on checkout. Without it a Windows clone would fail
cases for a reason that has nothing to do with the code.

---

## 8. Related documents

* `docs/MAINTAINING.md` — the clean-room rule and why it is absolute, repo layout,
  regenerating the corpus, re-capturing goldens, running the suite and its regression
  baseline, **the fixture rigs and the canvas technique with its read-back arithmetic
  (§8)**, the never-fit-to-the-corpus rule with the incident that motivates it (§8.7), and
  adding an analysis. Read it before changing anything in `open-king-core::ibdseg`.
* `README.md` — the cold-reader entry point. It quotes this file and claims nothing this
  file does not measure; `CITATION.cff` carries the same credit machine-readably.
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
  `IBD1Seg` overlap rule and §9 the one clause deliberately left out. **`19-…` closes
  `IBD2Seg` with the segment-fringe clause** and **`20-seg-writer.md` is the two writer rules
  of §4.3 and §4.5** — the only document in the log whose evidence is the captured output
  itself rather than a constructed fixture. **`20-seglength-floor.md` derives the run merge
  and `21-push-merge.md` corrects it** — the latter is the last word on the IBD2 merge and on
  the one-word push, and it supersedes `17-…` §6 and `20-…` §3/§5/§7 where they disagree.
  **`22-screen.md` is §5.7** — the two-stage screening count, and the only document in the
  log that ends in a deliberate non-landing: it proves by algebra and by measurement that the
  statistic is not the kinship over any marker subset, pins the law it does obey
  (`k_screen = 0.5 + R(k − 0.5)`, `R ≡ 1` below 32 768 markers), and declines to fit `R`.
  **`23-gap-bound.md` is the last segment document and closes `--seglength 10`** — the gate
  window's own length bound and the IBD1 merge's budget word set. Its §8 supersedes
  `21-…` §8 as the standing description of what is open; read it with §4.6 here, which
  measures the same engine out of sample.
* `docs/research/fixtures/` — the rigs. `fixlab.py` builds a fileset and drives the
  reference (`$KING` repoints it at our build); `gate8.py` brackets the `--degree 1` clause;
  `segfit.py` is the chunk-scan canvas; **`segcanvas.py` is the `.seg`-native canvas of §5.0**
  (6 416 cached reference answers) and **`ibd1canvas.py` the same canvas built IBD1-side up**
  (1 013); **`gradebinary.py` grades our build on both** without touching either cache;
  `avfs.py` regenerates the `--build` pedigree shapes of §6.2 and `avfs_score.py` scores the
  `AV.FS` statistic over them; **`screencanvas.py` is the `--related` screening canvas of
  §5.7** — the single-pair probe and the clone-fraction boundary, with `--facts` re-measuring
  every number that section quotes — and **`screenweight.py` is what that section weights**,
  the differential MAF-band probe whose measurements retired the frequency-standardised lead;
  **`mergelab.py` is the run-merge canvas of `20-…` and `push1.py` the two-word instrument of
  `21-…`**, which reproduces all four of that document's bisections in one run;
  **`window1.py` is `23-…`'s window-bound canvas** and carries its held-out draws (§7–§8);
  **`chrprobe.py` reads the reference one chromosome at a time on the corpus's own data**, by
  muting every other chromosome for the probe pair rather than subsetting the `.bim`, and is
  the instrument that localised the 10 Mb residual after two campaigns had guessed wrong;
  **`oosseg.py` is the out-of-sample differential of §4.6** — whole fresh filesets on unused
  seeds, byte-diffed against the reference, and the grader to use now that the corpus is
  saturated; **`screendeflate.py` is `22-…`'s instrument**, whose `facts` subcommand
  re-measures every number in §5.7 and which documents the trap that cost that campaign
  hours (KING dies with `FATAL ERROR - Too many first alleles as the major allele` unless a
  synthetic fileset codes A1 as the **minor** allele, which silently turns every bisection
  into "no bracket"); and **`buildlog.py` scores `<prefix>build.log` itself** (`rules`,
  `order`, `pairs`). Their `work/` output is gitignored and disposable — the JSON caches are
  not, and must only ever be written by the reference binary.
* `tests/parity/fit/` — Python mirrors of the committed engine, kept honest by
  `check_mirror.py`. `chunk.py` keeps the superseded `--ibs` rule alive beside the committed
  one so the before/after scorecard reproduces; `segtry.py` is the `.seg` port trial of §5.9;
  **`seg17.py` scores the `.seg` IBD2 caller, `seg18.py` the `IBD1Seg` overlap rule,
  `seg19.py` the IBD2 fringe, `seg20.py` the run merge, `seg21.py` the push and IBD2-merge
  corrections and `seg23.py` the window bound and budget word set** over the whole corpus in
  about a second each, every one printing the retired rule beside the one that replaced it.
  `engine.py` itself pins four named parameter bundles — `RETIRED`, `FRINGE18`, `PROP19` and
  the committed `BASE` — plus the knobs `merge21`, `push_fraction`, `window_fraction` and
  `merge_span` that step `BASE` back to the trees `20-…` and `21-…` shipped, so every
  scorecard quoted in `17-…` through `23-…` re-runs from one file (§4.4). `check_mirror.py`
  asserts the mirror still reproduces the binary at **all three** floors; a default-only
  check once passed while the mirror was wrong, and that is why.
