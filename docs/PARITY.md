# Parity with KING 2.3.2 — the measured claim

This is open-king's authoritative statement of what it reproduces and what it does not.
Every number in it was measured against the reference binary `king` 2.3.2 — nothing here is
an estimate, and nothing is rounded in our favour. §1 lists the commands that reproduce the
headline, the per-file table and both segment scorecards on the tree this document ships
with; every figure outside those comes from a rig that is **named where the figure is
quoted**, so any claim here can be re-run rather than taken on trust.

> **What "the reference" means here, and the caveat that applies to every number below.**
> One binary: `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`, Mach-O 64-bit **arm64**, run on
> macOS (Darwin 25.5.0, Apple silicon). Every rule in `king-core::ibdseg` and every figure in
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

> **Headline: 475 of the 480 captured reference invocations are byte-identical (99.0 %).**
> The harness self-check — the reference replayed against its own captures — is **480/480**.
>
> **Every one of the 5 that are not** is one of three named causes, and none of them is new:
>
> | cases | cause | where |
> | ---: | --- | --- |
> | 2 | `IBD1Seg`/`IBD2Seg` under `--seglength 10` — the residual the push and IBD2-merge corrections did not close; `--seglength 5` is now exact | §4.4, §5.0, `21-push-merge.md` |
> | 2 | `--related`'s two-stage screening count on stdout | §5.7 |
> | 1 | `--build`'s `<prefix>build.log` is unimplemented | §6.2 |
>
> **And the row-level scorecard, at all three captured reporting floors** — the same 982
> pairs, graded on the printed columns (`tests/parity/fit/scorecard.py`). The case count
> above and this table are different measurements of the same tree, and **both belong in any
> summary of it**; see "Two numbers" below:
>
> | `--seglength` | all four columns | `IBD1Seg` | `IBD2Seg` | of | extra | missing | MAE | worst |
> | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
> | 3 Mb (default) | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
> | 5 Mb | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
> | 10 Mb | 970 | 970 | 972 | 982 | 0 | 0 | 0.000046 | 0.0081 |
>
> **At the default floor and at `--seglength 5` the segment engine is exact** — all 982 rows
> byte-exact on all four printed fields, mean `PropIBD` error a true 0, 0 extra and 0 missing
> pairs. Every remaining `.seg` failure is a `--seglength 10` case. Reporting only the first
> row of that table would be true and misleading; reporting only the case count would hide
> that 99.0 % of cases rests on 98.8–100 % of rows depending on the floor.
>
> **Both numbers moved this pass, and by different amounts.** The case headline went
> 472 → **475** (three cases). On the same tree and on the same scale as the table above,
> the row scorecard went 947/959/947 → **982/982/982** at 5 Mb and 943/960/945 →
> **970/970/972** at 10 — 62 exact rows, and the last 35 rows at the 5 Mb floor. Neither
> figure is a restatement of the other.
>
> **What the set of reported pairs does:** it is exactly right everywhere. **0 extra and 0
> missing rows on every output file in the corpus**, all 982 `InfType` labels, and every
> `Error` — at every floor.
>
> **Twenty-nine of the thirty-one output files this project writes are byte-identical in
> every case that produces them** — `.kin`, `.kin0`, `X.kin`, `X.kin0`, `X.seg`,
> `cluster.kin`, `.ibs`, `.ibs0`, `.con`, `allsegs.txt`, `splitped.txt`, `unrelated.txt`,
> `updateparents.txt`, the `autoQC` set and the rest. Only `<prefix>.seg` (48 of 50 cases)
> and `<prefix>build.log` (7 of 8) differ anywhere. **Eleven analyses are byte-identical on
> every dataset that runs them**, `--related` is 64/65 and `--build`/`--ibdseg` are the rest.
>
> **Two numbers, and neither one alone is honest.** §3 and the case count grade **whole
> files**: a case turns `PASS` only when *every* row of *every* file it writes is byte-exact.
> §4.4 grades **rows**. They move independently, in both directions, and quoting either by
> itself misleads:
>
> * A change can win thousands of rows and **no** cases. The `.seg` IBD2 bridge/gate
>   correction (`17-seg-caller.md` §14) took the binary from 5 723 to **6 000 of 6 000** on
>   constructed canvases while changing no corpus row and no case at all.
> * A change can win 28 cases without touching **a single estimate**. `<prefix>.seg` computes
>   `PropIBD` from its own printed columns rather than from the underlying totals, and lists
>   its rows in 16-sample blocks rather than by sample index (`20-seg-writer.md`). Both are
>   writer rules; between them they flipped 28 cases and moved no number.
> * And a change can move both, which is what the run merge did: 464 → 472 cases, and
>   `IBD1Seg` from 844 to 960 of 982 rows at the 10 Mb floor.
>
> So the headline is quoted here **with** the row-level scorecard of §4.4, never instead of
> it. The case count is the stricter of the two and is the one the CI baseline asserts.
>
> **How the headline got here:** 436 → 464 on the two writer rules, → 466 when
> `<prefix>X.seg` landed (§6.1), → 472 when the `--seglength` run merge landed
> (`20-seglength-floor.md`), → **475** when the one-word push and the IBD2 side of that
> merge were re-measured (`21-push-merge.md`). §5.0 says which grader to use for what.
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

# pass/fail for all 480 cases  -> "475 PASS, 5 FAIL, 480 total"
python3 tests/parity/run_parity.py --impl ./target/release/king

# the same, as a regression gate: compare per case AND per output file against the
# recorded outcome, and fail on any difference in either direction
python3 tests/parity/run_parity.py --impl ./target/release/king --baseline

# how big each remaining gap is: rows, columns, mean and worst absolute error
python3 tests/parity/measure_gaps.py --impl ./target/release/king -q

# per-dataset roll-up for one output file
python3 tests/parity/measure_gaps.py --impl ./target/release/king -q --by-dataset king.seg

# harness self-check: the reference against its own captures must be 480/480
python3 tests/parity/run_parity.py --impl "/path/to/reference/king"

# the row-level .seg scorecard at the default floor, per dataset
KING_GOLDEN=tests/parity/golden cargo test -p king-core --test ibdseg_parity -- --nocapture

# the row-level .seg scorecard at ALL THREE captured floors (3 / 5 / 10 Mb) -- §4.4.
# Measured from the binary against the goldens; the table above covers only the default.
python3 tests/parity/fit/scorecard.py                 # add --per-dataset or --residual

# our binary against the reference on the constructed canvases (§5.0)
GRADE_CACHE=$TMPDIR/g2.json python3 docs/research/fixtures/gradebinary.py target/release/king
GRADE_CACHE=$TMPDIR/g1.json python3 docs/research/fixtures/gradebinary.py target/release/king --ibd1

# the two writer rules, from the captures alone -- no binary, no engine (§4.3, §4.5)
python3 docs/research/fixtures/segwriter.py

# the two differential probes that are not part of the capture corpus
python3 tests/parity/probes/degree_filter.py --ref "/path/to/reference/king"
cd docs/research/fixtures && python3 gate8.py

# the Python engine mirror must still reproduce the binary's own output, at all three
# floors -- the merge of 20-... is dormant at the default one, so a default-only check
# cannot see it (and for a while did not)
cd tests/parity/fit && python3 check_mirror.py     # -> "MIRROR OK"
python3 seg17.py && python3 seg18.py && python3 seg19.py && python3 seg20.py  # historical
python3 seg21.py                                   # the committed rule, floor by floor
```

`run_parity.py` and `measure_gaps.py` are Python 3 standard library only, regenerate the
input corpus automatically on first run (~20 s) and need no reference binary. The probes and
the canvas rigs drive the reference directly. `run_parity.py` exits 0 when every case passed,
1 when at least one failed, 2 on a harness error; with `--baseline` it exits 0 only when the
outcome matches `tests/parity/BASELINE.txt` exactly.

Measured on the tree this document describes:

| command | result |
| --- | --- |
| `run_parity.py --impl target/release/king` | **475 PASS, 5 FAIL, 480 total**, 876 output files byte-compared, 8 diff-excluded |
| `run_parity.py --impl target/release/king --baseline` | `baseline: MATCH (480 case(s))` |
| `run_parity.py --impl <reference>` | **480 PASS, 0 FAIL**, 876 files byte-compared — the normalization is complete and the goldens are self-consistent |
| `probes/degree_filter.py --ref <reference>` | 38 298 cases, **0 false-keep, 0 false-drop** |
| `probes/xseg_probe.py --impl target/release/king` | `<prefix>X.seg` out of sample on 1 040 built runs: emission gate **1 040/1 040**, bytes at the default floor **625/625** of the runs whose autosomal `.seg` also matches, raised floors **140/160** (§6.1) |
| `docs/research/fixtures/gate8.py` | brackets the `--degree 1` IBD2 clause to (0.0789, 0.0829] — its ladder refuses at `PropIBD` 0.0789 and accepts at 0.0829 |
| `gradebinary.py target/release/king` | **6 000 / 6 000** canvases — the release binary against the reference's own readings, `IBD2Seg` |
| `gradebinary.py target/release/king --ibd1` | **540 / 540** on the closed families **and 60 / 60** on the one that was open before the run merge — 600 / 600 in total |
| `segwriter.py` | `.seg`'s `PropIBD` rule consistent on **4 172 / 4 172** reference rows with **0** refutations, refuted on `.kin` (42 rows) and `cluster.kin` (3); row-order block size **uniquely 16** over 2..80 across all 50 `.seg` captures |
| `cargo test -p king-core --test ibdseg_parity` | `TOTAL gold=982 row=982 est=982 infType=982 missing=0 extra=0 meandPropIBD=0.000017 worst=0.0001` |
| `tests/parity/fit/scorecard.py` | the three-floor row scorecard of §4.4: `982/982/982` at 3 Mb, `982/982/982` at 5, `970/970/972` at 10; 0 extra and 0 missing at all three |
| `tests/parity/fit/check_mirror.py` | **MIRROR OK** — the mirror reproduces the binary's `.seg` columns on all 982 rows **at each of 3 / 5 / 10 Mb** (2 946 rows) and 861 `MaxIBD2` values, 10 datasets |
| `tests/parity/fit/seg19.py` | the `19-…` scorecard, unchanged: `806 / 982 / 982` at 3 Mb, `755 / 910 / 946` at 5, `713 / 844 / 937` at 10 — the "before" the merge is measured against |
| `tests/parity/fit/seg20.py` | the merge's own before/after, both under the retired `.kin` `PropIBD` rule: `19` (no merge) `806 / 755 / 713` exact at 3 / 5 / 10 Mb, `20` (merge) `806 / 795 / 793` |
| `tests/parity/fit/seg21.py` | the push + IBD2-merge corrections, same scale: `20` (previous commit) `806 / 795 / 793`, `21` (committed) `806 / 817 / 811` exact, `IBD1Seg` `982 / 982 / 970` and `IBD2Seg` `982 / 982 / 972` |
| `tests/parity/fit/seg18.py` | `18-…`'s own numbers, unchanged: committed `exact 747  ibd1 982  ibd2 896  MAE 0.000067`; retired overlap rule `709 / 826 / 896` |
| `cargo test --workspace` | **314 passed, 0 failed, 1 ignored** |
| `cargo clippy --workspace --all-targets -- -D warnings` | clean |
| `cargo fmt --all --check` | clean |
| a pristine copy of the tree, `cargo build --release` | succeeds in **8.3 s** (two cold builds: 8.52 and 8.28); `Cargo.lock` has 15 packages — the 3 workspace crates and 12 external |
| that clean-tree binary, re-run through `run_parity.py --baseline` | **475 PASS, 5 FAIL**, `baseline: MATCH` — the published counts do not depend on a warm `target/` or a pre-generated corpus, and the corpus regenerates itself from Python 3 stdlib in the same run |

By capture group: `apps` **90/91**, `core` **103/104**, `ibdseg` **62/65**, `params`
**220/220**. By analysis: `--kinship`, `--duplicate`, `--bysample`, `--bySNP`, `--autoQC`,
`--ibs`, `--cluster` **13/13** each; `--unrelated` **26/26**; `--build` **12/13**;
`--related` **64/65**; `--ibdseg` **50/52**; `--related --ibdseg` **12/13**; `params`
**220/220**.

The 5 failures, by shape: **2** in the `ibdseg` group (`bigish` and `multifam` at
`--seglength 10` only), **2** `bigish --related --degree 2` cases differing on one stdout
line, **1** `apps/bigish__build`.

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
| `--build` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 12/13 |
| `--related` | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | **5/5** | 4/5 | 64/65 |
| `--ibdseg` | **4/4** | **4/4** | **4/4** | 3/4 | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | **4/4** | 3/4 | 50/52 |
| `--related --ibdseg` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | 0/1 | 12/13 |
| `--ibs` | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **1/1** | **13/13** |
| flag plumbing + error probes (`params`) | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | \- | **220/220** |
| | | | | | | | | | | | | | | **475/480** |

The `params` group is 220 invocations that exercise the command-line surface rather than one
dataset: every `--prefix` shape, `--cpus`, `--sexchr`, `--degree`, `--minConc`,
`--seglength`, alternate `--fam`/`--bim` inputs, malformed and missing files, and the banner
in each case. All 220 are byte-identical.

**Read the `--ibdseg` row carefully.** Both `3/4`s in it are a **`--seglength 10`** case and
nothing else — `multifam` and `bigish` lose that one capture each. The bare invocation,
`--degree 2` and `--seglength 5` pass on all thirteen datasets, `sexchr` included now that
`<prefix>X.seg` is written (§6.1). `bigish`'s `--related --degree 2` cases lose on one stdout
line (§5.7).

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
| `<prefix>.seg` | 50 | 48 | 4 172 | 12 | **99.71 %** |
| `<prefix>build.log` | 8 | 7 | 18 | 18 | see §6.2 |
| `<prefix>X.seg` | 2 | **all 2** | 28 | 0 | **100 %**, §6.1 |

**Two files in the whole project differ anywhere**, and one of those is a file we do not
write at all. Everything above the `.seg` line — including the entire 16-column `--related`
layer, which had 450 differing `.kin` rows and 28 `.kin0` rows two campaigns ago and 16
`cluster.kin` rows one campaign ago — is now byte-identical in **every** case that produces
it. Row identity is matched on the identifier columns before any comparison, so across the
whole corpus there are **0 extra and 0 missing rows** on every file, `X.seg` included now
that it is written (§6.1).

Of `<prefix>.seg`'s 12 differing rows, **every one is in a `--seglength 10` capture**, and
they sit in two of the fifty captures. At the default floor and at `--seglength 5` the file
is byte-identical on all 13 datasets.

**stdout, stderr and exit status.** 477 of the 480 cases match stdout byte-for-byte after
the normalization of §7. **3 cases differ on stdout**, and two of them differ on *nothing
else* — the only cases in the suite that fail on console output alone:

| cases | stdout line that differs | cause |
| ---: | --- | --- |
| 2 | `bigish --related --degree 2` — `Stages 1&2 (with 32768 SNPs): 36 pairs` vs `50 pairs` | the two-stage screening bound, §5.7 |
| 1 | `apps/bigish__build` | §6.2 |

Those three cases are also three of the **five** remaining failures. The other two fail on
`king.seg` bytes alone, and both are `--seglength 10` captures:
`bigish__ibdseg_seglength10` and `multifam__ibdseg_seglength10`. At the default floor and at
`--seglength 5`, `<prefix>.seg` is byte-identical on every dataset.

---

## 4. The gaps, measured

Row counts in this section use `measure_gaps.py`'s denominator, which is **rows inside the
cases that differ**, not rows corpus-wide — it is the tighter, less flattering number. §3
gives the corpus-wide view of the same data.

### 4.1 The segment columns — 2 of the 5 failures

Re-measured after the one-word push and the IBD2 merge were corrected
(`docs/research/21-push-merge.md`). The history of this table, on `measure_gaps.py`'s own
denominator: 232 of 1 862 rows before the run merge, 74 of 1 658 after it, **12 of 867**
now — the denominator shrinks with the numerator because it counts rows inside the cases
that still differ, and two cases now do rather than five.

| file | rows differing | of | +extra | −missing | column | mean abs err | worst |
| --- | ---: | ---: | ---: | ---: | --- | ---: | ---: |
| `king.seg` | 12 | 867 | **0** | **0** | `IBD1Seg` 12 rows | 0.009742 | 0.0277 |
| | | | | | `IBD2Seg` 10 rows | 0.009080 | 0.0168 |
| | | | | | `PropIBD` 12 rows | 0.003733 | 0.0081 |

That is the whole numeric gap in the project. `king.kin`, `king.kin0`, `kingX.kin`,
`kingX.seg`, `kingcluster.kin`, `king.ibs` and `king.ibs0` are **not** in this table: every
one of them is byte-identical in every case, which is `IBD1Seg`, `IBD2Seg`, `PropIBD`,
`InfType`, `Error`, `MaxIBD2` and `Pr_IBD2` exact on all 4 805 `--related` rows and all
21 561 `--ibs` rows.

**All 12 rows are in a `--seglength 10` capture.** `PropIBD` appears here only because it
is derived from the two columns beside it (§4.3): every one of its rows also has a wrong
`IBD1Seg` or `IBD2Seg`, and it contributes no error of its own at any floor.

Where the 12 sit, per case (`fit/scorecard.py --residual` lists them individually):

| case | rows | differing | `IBD1Seg` | `IBD2Seg` |
| --- | ---: | ---: | --- | --- |
| `bigish` 10 Mb | 763 | 7 | 7, **all too low** | 6, **all too high** |
| `multifam` 10 Mb | 104 | 5 | 5, **all too low** | 4, **all too high** |

The one-sidedness is the shape of the residual and is discussed in §5.0: it is one fault
seen twice, not two, and after `21-push-merge.md` it points the other way from the residual
that document opened with. This caller now merges an IBD2 pair of runs where the reference
does not; the IBD2 territory it invents *is* then subtracted from `IBD1Seg`, so the same
event pushes one column up and the other down on the same rows.

**The arithmetic says so directly.** Writing `d1`, `d2` for ours minus the reference's
printed column, `d1` is negative on all 12 rows and `d2` positive on all 10 that differ, and
`|d1| = |d2|` to within one printed ulp on **8** of the 12. The other four are the two rows
where the reference reports no IBD2 at all (`multifam C_M/D_C2`, `bigish B21_G_M/B21_C2`:
`d2 = 0`, only `IBD1Seg` moves) and two where `|d1| > |d2|` (`multifam C_M/D_M` −0.0277
against +0.0168; `bigish B26_C1/B26_C2` −0.0070 against +0.0046) — so a second, smaller
effect is present on top of the invented merge, and those two rows are the sharpest handle
on it.

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
  does not move it. `king_core::ibdseg::inf_type` stays ungated because `.seg` genuinely
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

Committed as `king_core::ibdseg::seg_prop_ibd`, reaching **one column of one file**:
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
`cargo test -p king-core --test ibdseg_parity`; "both est." is the two estimate columns only,
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

**The three-floor scorecard, rules unchanged.** The raised floors were never used to fit the
caller's geometry, and had no part in finding the two writer rules; they *were* the evidence
for the run merge, which is why that rule was bisected on constructed canvases against the
reference and validated on unused seeds rather than tuned here (§5.0, "never fit to the
corpus"). Measured from the binary against the goldens by `tests/parity/fit/scorecard.py`:

| floor | all four | `IBD1Seg` | `IBD2Seg` | of | extra | missing | MAE | worst |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 Mb (default) | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
| `--seglength 5` | **982** | **982** | **982** | 982 | 0 | 0 | 0.000000 | 0.0000 |
| `--seglength 10` | 970 | 970 | 972 | 982 | 0 | 0 | 0.000046 | 0.0081 |

`MAE`/`worst` here are printed-column against printed-column — what a user diffing two files
sees — so the default floor reads a true 0.

`--seglength 5` now reads a true 0 as well: after `21-push-merge.md` it is exactly as
byte-exact as the default floor. At 10 Mb **"all four" equals the number of rows whose two
estimate columns are right** (970): `PropIBD` adds no error of its own anywhere. Everything
left is `IBD1Seg` and `IBD2Seg` at the 10 Mb floor — §4.1 and §5.0.

**What the last two campaigns moved**, on these same 982 rows and on one consistent scale.
Measured from `fit/engine.py`'s pinned bundles — `replace(BASE, merge=False, merge21=False,
push_fraction=None)` is the tree as `19-…` left it, `replace(BASE, merge21=False,
push_fraction=None)` is the tree as `20-…` left it, and `BASE` is what ships:

| floor | | all four | `IBD1Seg` | `IBD2Seg` | MAE | worst |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3 Mb | `19-` no merge | 982 | 982 | 982 | 0.000017 | 0.00005 |
| | `20-` run merge | 982 | 982 | 982 | 0.000017 | 0.00005 |
| | **`21-` committed** | **982** | **982** | **982** | **0.000017** | **0.00005** |
| 5 Mb | `19-` no merge | 900 | 910 | 946 | 0.000156 | 0.0641 |
| | `20-` run merge | 947 | 959 | 947 | 0.000080 | 0.0111 |
| | **`21-` committed** | **982** | **982** | **982** | **0.000017** | **0.00005** |
| 10 Mb | `19-` no merge | 832 | 844 | 937 | 0.000384 | 0.0916 |
| | `20-` run merge | 943 | 960 | 945 | 0.000129 | 0.0111 |
| | **`21-` committed** | **970** | **970** | **972** | **0.000062** | **0.0081** |

**Exact rows up and mean error down at both raised floors at every step, and nothing
whatever at the default** — where neither rule can fire on real marker spacings. That is the
bar a change has to clear to land here at all (`MAINTAINING.md` §8.6), and it is why these
two did. Nothing was landed this campaign that traded exact rows against MAE; the four
clauses of `21-…` were graded together and each one separately (`fit/seg21.py grid`), and
each is either neutral or an improvement on both.

**Two scales, and three graders that use them.** `engine.py`'s MAE compares our *unrounded*
`PropIBD` to the reference's printed one, so a fully exact floor reads 0.000017 rather than
0 — half of the exact rows sit half a printed ulp from the value they round to (see the note
above the ladder). `fit/scorecard.py` compares printed against printed, which is what a user
diffing two files sees, so the same three floors read **0.000000 / 0.000000 / 0.000046**.
And `fit/seg19.py`, `seg20.py` and `seg21.py` grade `PropIBD` with the retired **`.kin`**
rule, so they report the same two estimate columns but a different `exact`: `seg21.py` prints
806 / 817 / 811 where the tables here print 982 / 982 / 970. All three are correct on their
own scale; none of them is interchangeable with another.

**Detection is finished; what is left is length, at raised floors only.** Splitting the 982
primary rows by whether the reference reports any IBD2:

| reference row | rows | both estimate columns exact |
| --- | ---: | ---: |
| `IBD2Seg == 0.0000` | 823 | **823** |
| `IBD2Seg > 0` | 159 | **159** |

**How to grade further work on it.** Not with the case count and not with the exact-row
count at 3 or 5 Mb — all three are saturated. Grade at `--seglength 10`, on `IBD1Seg` and
`IBD2Seg` (`tests/parity/fit/scorecard.py`, or `engine.py`'s `score_seg(suffix=…,
min_bp=…)`), and out of sample with `gradebinary.py --ibd1`, which is now **600/600** — its
one previously open family closed when the run merge landed, so it no longer localises the
residual and a *new* canvas will be needed to. A candidate rule lands only if the
exact-row counts **and** the MAE both improve (or MAE is unchanged); anything that trades
one for the other is reported with its numbers and rejected.

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
| **The two `PropIBD` rules.** `.kin`/`.kin0`/`X.kin`/`cluster.kin` print `IBD2Seg + IBD1Seg/2` at full precision; `.seg` prints `i2*1e-4 + i1*5e-5` off its own four-decimal columns. The reference disagrees with itself on 54 of 201 pairs it writes into both, so this is two rules and not one | `.kin` family byte-identical on **all 4 805** rows; `.seg` rule consistent on **all 4 172** captured rows with **0** refutations | §4.3, `20-seg-writer.md`, `fixtures/segwriter.py` |
| **`.seg`'s row order** — by 16-sample block, then by index | block size **uniquely 16** over 2..80 across all 50 `.seg` captures | §4.5, `20-seg-writer.md` §4 |
| **`InfType` and `Error`** | **no row anywhere differs**, over all 4 805 rows — not merely where the segments are exact | §4.2 |
| **The IBD1 caller, its boundary refinement, its gate and its `IBD1Seg` overlap rule** (`Scan::ibd1`, `ibd1_pieces`) | every clause bisected on an IBD1-native canvas; `IBD1Seg` exact on **all 982** primary rows and on every `.kin`/`.kin0`/`X.kin`/`cluster.kin` row; the binary matches the reference on **600 of 600** IBD1 canvases — 540 closed plus the 60 of the family that was open until the run merge landed | `18-ibd1-caller.md`, `20-seglength-floor.md`, `fixtures/ibd1canvas.py`, `gradebinary.py --ibd1` |
| **The `--ibs` IBD2 caller** (`Scan::ibd2_words`, the chunk scan) | exact on all **21 561** rows | §5.8 |
| **The `.seg` IBD2 caller** (`Scan::ibd2`) — word predicate, gate, reach, push, bridge and fringe | every constant bisected on a `.seg`-native canvas; the binary reproduces the reference on **6 000 of 6 000** word-aligned canvases and **504 of 504** fringe canvases; `IBD2Seg` exact on **all 982** primary rows | `17-seg-caller.md` §3–§7 and §14, `19-ibd2seg-residual.md`, `fixtures/segcanvas.py`, `fixtures/fringecanvas.py`, `gradebinary.py` |

**The `--seglength` run merge and the one-word push** (`20-seglength-floor.md`, corrected by
`21-push-merge.md`) — `Scan::merge_ok`, `Scan::join_runs`, `Scan::join_runs2`, the `armed`
flag. Both are committed and both are exact at 3 and 5 Mb. **The two passes are not the same
rule**, which is what `20-…` got wrong:

* **Shared.** After the gate has refused what it refuses, two runs are joined iff the gap
  between them is **strictly** under `--seglength` and a budget `cost·(bad − 2) ≤ X` passes
  over the interrupting words. A run the gate refused lies *inside* an interruption rather
  than ending one, and a merged call may **not** satisfy the >10 Mb pair filter — which is
  why `pair_segments` computes the longest segment from the *unmerged* call sets.
* **IBD1.** At most **2** unusable words; the gap measured run-to-run; `bad` = opposite
  homozygotes; `cost` 4; `X` = A1A1/A1A1 unless the het-vs-A1A1 markers alone reach 10.
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

Every constant was bisected on `fixtures/mergelab.py` and `fixtures/push1.py` against the
reference and validated out of sample — 360/360 held-out canvases at 5 and 10 Mb on three
unused seeds for the IBD1 merge, 357/360 on three further unused seeds for the corrected
IBD2 rule (the committed-at-the-time rule scored 343/360 on the same canvases), plus 600/600
independently drawn interruptions. **Never fitted to the corpus.** Between them they took the
headline 464 → 472 → **475**, `IBD1Seg` at 10 Mb from 844 to 970 of 982, and `--seglength 5`
from 900 to 982 of 982 on both estimate columns.

**At the default 3 Mb floor and at `--seglength 5` there is no segment residual at all.** All
982 primary rows are byte-exact on all four printed fields, printed-column MAE a true
0.000000 at both. What follows is entirely about `--seglength 10`.

**Not solved — one residual, plus two unrelated gaps.**

1. **An invented IBD2 merge at `--seglength 10`, 2 of the 5 cases.** 12 of 867 rows in the
   two captures that still differ, in two datasets (`bigish` 7, `multifam` 5); the whole
   corpus is exact at 3 and 5 Mb. The shape is sharper than it has ever been, and it is
   **inverted** from the residual `20-…` left behind:

   * **One-sided, and one fault seen twice.** Every wrong `IBD2Seg` is too **high** (10 of
     10) and every wrong `IBD1Seg` too **low** (12 of 12). This caller merges an IBD2 pair
     of runs where the reference does not; the IBD2 it invents is then subtracted from
     `IBD1Seg` by `ibd1_pieces`, so one event moves both columns. `|d1| = |d2|` to within a
     printed ulp on 8 of the 12 (§4.1); on the two rows where the reference reports
     `IBD2Seg 0.0000` the invented merge is the only IBD2 the pair has.
   * **Each wrong row is worth about one merged segment**: `d1` runs −5.5 to −19.2 Mb and
     `d2` +10.0 to +11.6 Mb across the twelve.
   * **Ruled out, measured** (`21-…` §8.1): the **word cap** — capping the IBD2 merge at 2,
     3 or 4 words changes nothing on this corpus, because at ~50 kb spacing a 10 Mb gap
     holds at most three words; the **budget constants**, bisected at equality on both
     `bad` and `X`; and the **push**, worth one row at 10 Mb and none at 5.
   * **Also ruled out earlier, and still ruled out:** allele frequency (padding-only changes
     flip nothing); rate-over-merged-segment readings; per-word rather than per-interruption
     budgets; and any linear function of the nine marker-kind counts (a perceptron over all
     nine leaves 36/600 misclassified — the switch at 10 makes the region non-convex).

2. **`<prefix>build.log` is never written** (§6.2). 1 case, and still segment-blocked: it
   needs exact segment *placement* plus the still-unidentified named-sib-pair rule. All 34
   measured `AV.FS` triples miss the printed 3 dp, mean residual +0.0039, one-sided high on
   every one — and `bigish`'s own five are +0.0008…+0.0052 — so the log's five
   `Join3/Join2` values would print wrong even with the rule scaffolding written.

3. **`--related`'s two-stage screening count**, 2 cases (§5.7). Not a segment problem: one
   stdout line, `36 pairs` against our `50`. §5.7 now closes off the whole "find the right
   marker subset" search space by measurement.

**The next experiments worth running, in order.**

1. **Sweep the gap, at both floors.** The gap is the only floor-dependent term left in the
   merge, the residual appears only at the larger floor, and every wrong row is worth about
   one merged segment. Build a canvas whose interruption gap sits **between 5 and 10 Mb** and
   sweep `--seglength` across both, on `fixtures/push1.py`'s instrument (the jumps of
   `IBD2Seg(L)` are the individual call lengths, so one sweep reads the whole behaviour).
   The question it answers: does `gap < seglength` acquire a **second, absolute** bound above
   5 Mb? This is the single most likely explanation of the 10 Mb column.
2. **Explain the four rows where `|d1| ≠ |d2|`.** Two have `IBD2Seg 0.0000` in the reference
   and only `IBD1Seg` moves; two (`multifam C_M/D_M`, `bigish B26_C1/B26_C2`) have `|d1|`
   strictly larger, so a second and smaller effect rides on top of the invented merge.
   `multifam` at 10 Mb is the smallest reproduction of the whole residual — 5 rows;
   `fit/scorecard.py --residual` prints them.
3. **Three hazards any new rig must respect** (`21-…` §8.2, §8.3), or it will misgrade its
   own boundary rows:
   * The reference **stops behaving like a floor** outside `1 ≤ L ≤ 10` Mb. Above ~10 Mb a
     14.06 Mb call reports as a constant 8.93 Mb at every larger floor; below 1.0 the flag
     behaves as though absent. No rule here is measured in either regime and none should be.
   * `.seg`'s floor test is strictly **`>`**, not `>=`: a call measuring exactly 5 100 000 bp
     is reported at `--seglength 5.099999` and dropped at `5.100000`. The engine compares
     `>=` and it never bites on the corpus — nothing there lands on an exact tie — but every
     fixture in the push/merge rig does.
   * Canvas read-back is a *measurement*, not an inspection: see `MAINTAINING.md` §8.3.
4. **Do not re-sweep the caller's constants.** Forty single-knob perturbations and all 32
   combinations of the two IBD1 endpoint rules crossed with the two IBD1 fringe rules were
   scored: none improves exact rows, none beats the committed MAE, and the committed values
   are the unique maximum of that grid (`20-seg-writer.md` §6). Likewise the merge's own
   knobs, swept in `fit/seg20.py grid` and `fit/seg21.py grid` — where dropping `reach` costs
   982 → 959/947 at 5 Mb, `hethet` 982 → 981/980, `push_half` one row on each column at
   10 Mb, and `no_cap` nothing on this corpus but four canvases out of sample.
5. **Five knobs the corpus cannot see at all** — `bridge_rule="17"`, `gate_end="right"`,
   `inf2_ibs1b=True`, `ibd1_clip_ibd2=True`, `clip_before_len=False` all score identically to
   the committed engine on every corpus row. They were settled on the canvases (`17-…` §14)
   and the canvases remain the only evidence for them. If you change one, grade it there.

**And the standing caveat, restated where it bites hardest.** Everything above is measured
against **one** KING build (2.3.2, macOS arm64 — see the note at the top of this file). The
segment caller is the part of KING whose numerics its own release notes say changed most
across 2.1.x–2.2.x, and the algorithm is unpublished. A residual of 12 rows at one raised
floor is a statement about 2.3.2 on this host and nothing wider.

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
detected`; open-king prints `50`. This is the only stdout line `--related` gets wrong. It
costs 2 cases — `core/bigish__related_degree2` and `ibdseg/bigish__related_degree2_ibdseg` —
and it is the only `--related` failure left. The consequence is contained to that one line:
`.kin0`'s row set comes from the exhaustive re-estimate below it and is byte-correct at every
degree, including in both of those cases, and both of their `.seg` files are byte-identical.

**What the stage is, measured.** `docs/research/fixtures/screencanvas.py` turns the aggregate
count into per-pair and per-marker measurements by driving the reference over two constructed
canvases: `bigish`'s 167 unrelated fillers (which alone print `No close relatives are
inferred.`) plus exactly one candidate pair, and a pair whose kinship one cloned marker set
tunes continuously. `--facts` re-measures all of the following in about two minutes; none of
it is read off `bigish`'s own answer.

* **The stage is per-pair.** Run the 47 candidate pairs one at a time and the count is 0 or 1
  each time; it is 1 on **36** of them — the whole fileset's number. So it is not a budget, a
  cap, a ranking or a per-block bound, and six sample permutations do not move it.
* **It is not a marker prefix, stride or decimation**, which is what this repo implements. A
  10 000-marker clone window at `[40000, 50000)`, where the first-32 768 estimate reads
  0.0020, is accepted at the same true kinship as one at `[0, 10000)`; stride-2/3/4 clone
  sets at every offset agree.
* **It is a threshold on kinship that sits above the printed cutoff** — 0.0700 against 0.0625
  at *n* = 167, *m* = 50 000, lossy in exactly the direction that turns 50 into 36. Written
  `k_screen = 0.5 + R*(k − 0.5)`, one `R` fits both degrees (1.0186 at cutoff 0.1250, 1.0176
  at 0.0625); an additive offset disagrees by 18 % between them, a plain multiplicative one by
  a factor of two.
* **`R` is exactly 1 while *m* ≤ 32 768** (0.99995 at 32 768, 0.99993 at 33 280), grows with
  the map (1.0106 / 1.0128 / 1.0204 / 1.0176 at 36 864 / 40 000 / 45 000 / 50 000) and *falls*
  with the sample count (1.033 at *n* = 100 to 1.018 at *n* = 167, smoothly, no 16/32-sample
  step), and varies pair to pair (1.018…1.026 over six filler pairs). It was never below 1.
* **The screening statistic reads the other samples' genotypes.** Drive the fillers to one
  homozygote at 17 232 random markers and a pair related *only* there is rejected at kinship
  0.154, while the same pair related on the rest is accepted at 0.062; replace the whole
  background with HWE-consistent random genotypes and the threshold passes 0.25. No marker
  subset can do that — a uniformly spread clone set meets any subset in its own proportion —
  so the stage is estimating with **sample-level allele frequencies**, which is also what the
  *n* dependence is.

**What it weights** (`docs/research/fixtures/screenweight.py`, all of it on constructed
canvases). The `m ≤ 32768` reading is exact and is the binding constraint: 12 fresh bisections
at m = 20 000 / 30 000 / 32 768 give `R = 1.00010 ± 0.00005` at a boundary kinship of 0.06252,
so below the threshold the screen *is* the robust kinship on the whole map under a strict `>`,
and no change of estimator, denominator or standardisation can hide there. Above it the
deflation is **affine, by prediction rather than by fit**: `R = 1.01779 ± 0.00062` at degree 2
predicts a degree-1 boundary of 0.13161 against 0.13979 for a multiplicative deflation, and the
measured value is 0.1316. The rule is **additive** over markers — mixing half of one band's
boundary set needs 3 906 of another band against a parameter-free prediction of 4 065, and a
quarter needs 6 001 against 6 098. And the per-marker weight is **`2p(1−p)` below 32 768
markers** (flat to `1.03 ± 0.05` over ten MAF bands, as any consistent kinship estimator must
be) but **not above it**: at m = 50 000 it ramps, `w/2pq` = 0.48…1.00 over MAF 0.16…0.35, and
markers below MAF ≈ 0.15 carry none at all — four bands stay rejected with all 2 133 / 3 278 /
3 772 / 3 739 of their markers cloned.

**The lead above is now refuted, not merely unvalidated.** A frequency-standardised estimate
weights every marker equally, so it predicts one boundary *marker count* whatever the band; the
measured counts are 3 978, 4 552 and 8 130, a factor of two. Het-weighting predicts one
boundary het *mass* and the top two bands give 1 910 and 1 898, agreeing to 0.6 %. This does not
depend on which frequencies the estimate uses, which is what its earlier out-of-sample failure
turned on.

**Nor is it this estimator on any subset of these markers.** Eight measured brackets — five
random clone families and three confined to a MAF band — each pin `k_S(size−1) ≤ 0.0625 <
k_S(size)` for the true screening set. Over 3 selectors × 11 sizes from 500 to 50 000, none of
the 33 satisfies all eight, structurally: a uniformly drawn clone family meets any subset in
its own proportion, so the five random boundaries read 0.069…0.072 for *every* subset and never
the cutoff. A subset bends the band boundaries but cannot supply the uniform deflation, so
**swapping the placeholder prefix for a better subset cannot close this line**. No rule is
fitted to it and the prefix stays: it reproduces the degree-1 count (18) on every map tried,
and swapping it for the whole map would lose that for nothing at degree 2 (47 against 36).

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
are the only known differences a *user* could hit while the suite stays at 475/480.

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
vary only the largest family size across 1, 2 and 3. It gained independent support from the
X.seg work: `--related`'s `<prefix>X.kin` obeys exactly that rule (twelve singleton families
write none; add one two-member family and it appears), which is now implemented.

### 5.11 Three more divergences the corpus cannot see

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

**2. `.fam` SEX fields outside {0,1,2} are read differently.** `king-io::fam::parse_sex` maps
everything but `1` and `2` to `0`; the reference does not. Measured by sweeping 43 spellings
through `X.seg`'s raw `Sex` column:

| reference prints | `.fam` field |
| --- | --- |
| `0` | `0` `00` `0.0` `-0` `-9` `-9.0` `x` `?` `NA` `na` `b2` |
| `1` | `1` `-1` `-2` `3` `9` `10` `12` `02` `002` `0002` `007` `+2` `1.9` `1e0` `M` `m` `male` `MALE` |
| `2` | `2` `20` `21` `22` `2.5` `2.9` `2x` `20x` `2e0` `F` `f` `female` `FEMALE` |

i.e. a leading `2`/`F`/`f` is female, a leading `M`/`m` is male, an otherwise-numeric field is
male unless it evaluates to `0` or `-9`, and anything else is unknown. Deliberately **not
implemented**: `parse_sex` is shared with `--autoQC`'s `updatesex`, `--bysample`, `splitped`
and `--kinship`'s X classification, none of which the corpus exercises with a field outside
`{0,1,2}`, so the change could not be regression-tested where it matters. Every corpus `.fam`
uses only `0`/`1`/`2`, where the two rules agree.

**3. The reference reads past the end of a marker array whose length is an exact multiple of
64.** On such an array it adds an *absolute* coordinate — `pos[last]`, not a difference — to a
pair's IBD2 total, so `IBD2Seg` comes out larger than 1. Swept at 320/384/448 markers (anomaly
on some pairs, always exactly `+pos[last]/D`) against 319/321/383/385/449 (no anomaly on any
pair, any seed). It shows identically in `X.kin` and `X.seg`, so it is the shared caller, and
it is unreachable from the corpus, whose arrays all end mid-word. Not emulated: it is an
uninitialised read, not a rule.

---

## 6. The remaining structural gap, in detail

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

**Held out** — `python3 tests/parity/probes/xseg_probe.py --impl target/release/king`, which
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

`apps/bigish__build` writes an **empty** `kingbuild.log` where the reference writes 18
lines. The other 12 `--build` datasets are byte-identical because they need no
reconstruction rules at all.

`kingupdateparents.txt` **is now written and byte-identical.** It is the half of the
reconstruction that reads no segment statistic — it carries only what the `RULE FS*` lines
decided — and it was pinned on **twenty held-out merge shapes**, not on `bigish`:
`docs/research/fixtures/build_shapes.py` builds them and re-runs the scorecard. Eighteen
are in scope (two are excluded because the reference renames every individual to
`<FID>-><IID>` when a `.fam` names a parent living in another family, a feature this binary
does not implement at all), **fifteen of them byte-identical on the file and on the console
tail**, and the remaining three carry identical `(IID, FATHER, MOTHER)` rows and differ only
in which cluster is named `KING1` — see the numbering bug below. The rules the shapes
pinned, none of which `bigish` alone shows:

* a sibship is a connected component of *inferred FS* ∪ *declares the same named couple*,
  and only a component the inference touched is reassigned;
* it takes the couple one of its members already declares if there is one, else the next
  synthetic pair `1 2`, `3 4`, … — one pair per sibship whatever its size, counted across
  the whole run;
* two sibships inside one cluster take consecutive pairs ordered by their first member
  under the ID comparator, **not** by `.fam` row order;
* a cluster holding an inferred **duplicate** contributes no rows at all, while one merged
  by **PO** alone contributes identity rows with nobody's parents changed — without ages
  the reference will not orient a parent-offspring pair;
* **nothing is written unless some sibship got parents**, so a run whose only merges are PO
  or duplicate leaves a zero-byte file, no `Update-parent information is saved…` line and
  the closing `No pedigrees can be reconstructed.` — even though its clusters *are* in
  `updateids.txt`. The old code keyed that tail on "did clustering merge anything", which
  is wrong on any such fileset; that is fixed.

**A held-out clustering bug the same rig found, and which is *not* fixed.** Merged clusters
are not numbered in family order: they are numbered by the *relationship type* of the pair
that joined them — every `Dup/MZ`-joined cluster first, then the `PO`-joined ones, then the
`FS`-joined ones, ties broken by family order. A fixture whose three clusters are joined, in
family order, by `FS`, `PO` and a duplicate is numbered `KING3`, `KING2`, `KING1`.
`unrelated::clusters` uses family order alone. Every merge in `bigish` is `FS`, and no other
corpus fileset merges at all, so all 480 captures are indifferent to it and nothing moves
either way; the fix belongs in `unrelated.rs`.

The case's third file, `kingupdateids.txt`, **already matches byte for byte** — the family
numbering, which original families each `KING<n>` absorbs, and the row order are all
correct, so none of the three obvious structural suspects (family numbering, parent
tie-breaking, sex assignment) is what fails here.

`kingbuild.log` needs more: `INFERENCE AV.FS` with a `Join3/Join2` statistic printed to
three decimals (`bigish`: 0.778, 0.801, 0.779, 0.827, 0.803) and `INFERENCE HS.UN2`. It is
left **empty** rather than half written, since a log carrying its `RULE` lines and not its
`INFERENCE` ones is no closer to the capture than no log at all.

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
log lines wrong: `apps/bigish__build` is blocked on §4.1 exactly as `apps/bigish__cluster`
is, and it is **not** an independent second cause on top of the segment residual.

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
It is an unprincipled knob and is not landed.) So §4.1 closing is **necessary** for
`apps/bigish__build`, and no longer demonstrably **sufficient**.
Rig: `docs/research/fixtures/avfs_score.py`.

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
  Bracketed to **(0.848, 0.901)** over the values measured — which does not separate 0.85, 0.875
  and 0.9.

`updateparents.txt` has since been written anyway — it cannot flip the case, but it is a
rule that generalises (see the top of this section), and shipping it removed one of the
three diffs and fixed the console tail on PO-only merges. The `build.log` derivation lives in `crates/king-cli/src/analysis/build.rs`'s module doc;
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
  adding an analysis. Read it before changing anything in `king-core::ibdseg`.
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
  the one-word push, it supersedes `17-…` §6 and `20-…` §3/§5/§7 where they disagree, and its
  §8 is the standing description of everything still open at `--seglength 10`.
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
  `21-…`**, which reproduces all four of that document's bisections in one run. Their `work/`
  output is gitignored and disposable — the JSON caches are not, and must only ever be written
  by the reference binary.
* `tests/parity/fit/` — Python mirrors of the committed engine, kept honest by
  `check_mirror.py`. `chunk.py` keeps the superseded `--ibs` rule alive beside the committed
  one so the before/after scorecard reproduces; `segtry.py` is the `.seg` port trial of §5.9;
  **`seg17.py` scores the `.seg` IBD2 caller, `seg18.py` the `IBD1Seg` overlap rule,
  `seg19.py` the IBD2 fringe, `seg20.py` the run merge and `seg21.py` the push and IBD2-merge
  corrections** over the whole corpus in about a second each, every one printing the retired
  rule beside the one that replaced it. `engine.py` itself pins four named parameter bundles —
  `RETIRED`, `FRINGE18`, `PROP19` and the committed `BASE` — plus the two knobs `merge21` and
  `push_fraction` that step `BASE` back to the tree `20-…` shipped, so every scorecard quoted
  in `17-…` through `21-…` re-runs from one file (§4.4).
