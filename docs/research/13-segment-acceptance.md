# Which runs the reference accepts — the last gap in the IBD-segment caller

**Status:** measurement. No KING source was read. Everything below is either a reading
taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) or a
score computed against the captured corpus. Where a rule is not pinned, it says so.

This document supersedes the open questions in `10-segment-rule-fixtures.md` §6 and
`11-segment-rule-fit.md` §7/§11 by **removing several hypotheses from the board** and
handing over three new instruments. It does not close the gap.

New scripts, all Python 3 + numpy, all under `tests/parity/fit/`:

| script | what it does |
| --- | --- |
| `score_seg.py` | runs the *implementation* on every dataset and scores `.seg` per column |
| `informative.py` | measures marker informativeness under each called segment |
| `matched.py` | brute-forces a decision stump over 25 features of the borderline runs |
| `subset.py` | re-asks the reference about one pair inside a 30-sample subset |
| `probe_seg.py` | **stretches the map** so the reference prints lengths it normally hides |
| `gridshift.py`, `sweep.py`, `sweepdata.py` | sweep the 64-marker word grid under one pair |

---

## 1. Reconciling the three recon reports

Three prior investigations disagreed in four places. Each was settled by running the
reference or by scoring against all 982 captured rows, never by averaging.

| question | `10-…fixtures` says | `11-…fit` says | settled |
| --- | --- | --- | --- |
| IBD2 minimum run | 1 word (measured on a constructed fixture) | 1 word | **1**. The corpus cannot tell 1 from 2 (626 exact rows either way), so the fixture decides. |
| IBD2 right boundary | clipped to the run's own words: `64W−1` intervals | reaches into the ending word | **reaches in.** Scored both: 626 exact / MAE 0.00327 extending, 625 / 0.00534 clipping. 982 real rows outrank one constructed fixture, but the disagreement is real and unexplained — see §6. |
| IBD2 word predicate | (not measured) | IBS1 ≤ 0 by `.seg`, ≤ 4 by `MaxIBD2` | **≤ 0**, i.e. no het mismatch at all: 626 exact at 0, 623 at 1, 622 at 4. |
| `InfType` | — | one FS clause | **two** FS clauses, per `12-segment-aggregates.md` §2. Re-verified here: computing `InfType` from each captured row's *own printed columns* reproduces **8 722 / 8 722** values across every `.seg`, `.kin` and `.kin0` in the corpus, against 8 715 for the one-clause rule. |

Two further corrections from `10-…fixtures` §7 were adopted: the pair filter is `>=`
10 000 000 bp, not `>` (bisected on a fixture: 9 990 000 absent, 10 000 000 present), and
the `--seglength` floor is likewise inclusive (the committed code already had the latter).

The engine in `crates/king-core/src/ibdseg.rs` now implements exactly the rule of
`11-segment-rule-fit.md` §9 with those corrections, and reproduces its score to the row:

```
rows 982   exact 626   IBD1 626   IBD2 824   InfType 973   extra 188   missing 0
MAE(PropIBD) 0.00327   worst 0.1753
--seglength 5 -> 626 exact,  --seglength 10 -> 731 exact,  same 188 extras, 0 missing
```

`<prefix>splitped.txt` is now written and is **byte-identical on all ten datasets**
(`crates/king-cli/src/analysis/splitped.rs`, ported from
`tests/parity/probes/splitped.py`). Every remaining `ibdseg` failure is now the numbers in
`king.seg` alone — plus `kingX.seg` on the single `sexchr__ibdseg_degree2` case.

---

## 2. The gap, stated precisely

Cross-tabulating the length of the IBS0-free **word run** underlying each pair's longest
called segment against whether the reference reports that pair, on `bigish`:

| clean words in the run | reference reports | reference refuses |
| ---: | ---: | ---: |
| 1 | — | (never called by us either) |
| **2** | **255** | **182** |
| 3 | 12 | 0 |
| ≥ 4 | 496 | 0 |

Every run of three or more IBS0-free words is accepted. Runs of exactly two are accepted
58 % of the time. **That single split is the whole of the 188 extra rows**, and it is the
whole of the remaining `--ibdseg` gap.

---

## 3. Instrument 1 — the 30-sample subset (`subset.py`)

Rebuild the fileset with the probe pair plus 28 fixed padding samples. The `.bim` is
copied verbatim, so the word grid, the usable segments and the denominator are identical
to the full run; only the allele frequencies change.

```
reported pairs reproduced in the subset: 10 / 10
extra    pairs reproduced in the subset:  0 / 10
```

**The verdict is carried entirely by the pair's own genotypes.** This kills the
"informativeness comes from the other samples" reading of
`10-segment-rule-fixtures.md` §5.3 as an explanation *for the corpus*: cutting the sample
from 200 to 30 changes every allele-frequency estimate and moves nothing.

It also makes the reference a **cheap oracle**: one verdict costs about 150 ms.

## 4. Instrument 2 — stretching the map (`probe_seg.py`)

A pair whose longest segment is under 10 Mb prints nothing at all, so a refusal and a
short call look identical. To tell them apart without touching a genotype:

* keep every marker, so the 64-marker grid — indexed by position in the retained-autosome
  array, not by base pair — is unchanged;
* multiply the target chromosome's positions by `K`, so a segment of length `L` measures
  `K·L` and clears the fixed 10 Mb pair filter;
* squash every other chromosome to 1 kb spacing so it drops out of the usable-segment list,
  leaving the denominator equal to the target chromosome alone — which buys back the print
  resolution the stretch costs, since one ulp of `%.4lf` is `D/10000`.

`K = 2` is the largest safe stretch for this corpus: a usable segment is cut wherever one
word spans over 10 Mb, and 64 × 2 × 50 kb = 6.4 Mb still clears it.

Result, lengths divided by `K` to return to the original scale:

```
group     pair              our core  our ext    reference
reported  B01_F/B10_C2         6.350   10.903       10.909
reported  B01_M/B04_C1         6.352   10.851       10.855
reported  B02_M/B18_P1         6.350   11.147       11.147
extra     B01_F/B22_F          6.346   10.203       absent
extra     B01_F/B25_F          6.349   10.852       absent
extra     B01_M/B08_M          6.353   10.054       absent
```

Two facts, both new:

1. **Where the reference calls a two-word run, our boundary rule is right to ~4 kb** —
   under a tenth of a marker. The boundary convention of `11-…fit` §4 is confirmed on real
   data, not just on `MaxIBD2` inversions.
2. **Where it does not, it calls nothing at all.** At `K = 2` any call over 5 Mb would have
   been printed. So the refusals are not short calls, and not a length-threshold effect.

## 5. Instrument 3 — sweeping the word grid (`sweepdata.py`)

Deleting the first `m` markers of the fileset shifts the global word grid by `m` and
changes nothing else. Sweeping `m = 0…63` under eight pairs (512 reference invocations,
~45 s) gives 512 labelled observations in which the genotypes are constant.

**The verdict moves with the alignment, in both directions.** Excerpt for one refused pair
(`B01_F/B22_F`), whose IBS0-free stretch is 158 markers:

```
 m  words  head  tail   our lo   our hi   our Mb   reference
 0    2     14    16     21426    21630   10.203      no
13    2     27     3     21426    21630   10.203      no
14    2     28     2     21426    21630   10.203      YES 10.207
16    2     30     0     21426    21630   10.203      YES 10.207
50    2      0    30     21426    21613    9.356      YES  9.363
```

At `m = 13` and `m = 14` the run is the same two words, the refined start and end are the
*same two marker indices*, the length is the same to the base pair — and the reference
answers differently. The only thing that changed is where the word boundaries fall inside
the IBS0-free stretch (`head`/`tail` are the clean markers spilling into the flanking
words).

So the missing rule is a function of the **bit positions of the IBS0 markers within their
words**, not of any count, length or density. That is precisely why no aggregate separates
the two groups.

The converse also happens: a pair the reference reports at `m = 0…41` is refused at
`m = 42…63` with a *longer* computed segment (9.599 Mb accepted, 9.795 Mb refused).

## 6. What is now ruled out

Each of these was tested here and fails; do not re-derive them.

* **Marker informativeness / allele frequency** (`informative.py`). Over the longest called
  segment of every reported and every extra pair in `bigish`, the two groups are
  indistinguishable in mean MAF (0.2514 vs 0.2450 median), in the count of markers above
  MAF 0.05/0.10/0.20, and in the expected IBS0 count under IBD0 (`Σ 2p²q²`: reported
  12.6–296, extra 12.4–16.8, fully overlapping the reported pairs of the same size). The
  fixture effect of `10-…fixtures` §5 is real but does not fire on this corpus.
* **Any local summary of the run** (`matched.py`). A decision stump over 25 features —
  extended length, core length, SNP count, left/right extension, IBS0 count in each
  flanking word, IBS1 and HetHet in the core, missing calls, marker spacing, distance to
  each end of the usable segment, chromosome, start-word index mod 2 and mod 4 — reaches
  **0.648** accuracy against a 0.583 majority baseline on the 369 borderline pairs. There
  is no separating feature among them.
* **The length of the IBS0-free marker stretch** and its physical span: reported 142–222
  markers, extra 143–218; 7.6–11.1 Mb against 7.4–10.9 Mb.
* **A pair-level screen using the wider sample.** §3.
* **A shorter call the 10 Mb filter then drops.** §4.
* **A clean-word count.** §5: verdicts flip while the count is fixed at two.

## 7. Where to look next

1. The reference's acceptance is alignment-sensitive but the *reported boundaries* are not
   (§5: identical `lo`/`hi` on both sides of a flip). Any candidate rule therefore has to
   be a second, independent test applied to the run — not a different way of placing its
   ends. A scan that walks markers inside the flanking words, or a scoring/DP formulation
   whose optimum happens to coincide with the last-IBS0 convention, both fit that shape.
2. The sweep is a far better fitting target than the corpus: 64 labelled observations per
   pair with genotypes held constant, at ~150 ms each. `sweepdata.py` writes it as CSV;
   eight pairs is 512 rows and the rule must explain every one.
3. The IBD2 boundary conflict of §1 is still open, and the fixture that produced it
   (`docs/research/fixtures/rig2.py`) is cheap to re-run. Deciding it properly would also
   settle whether IBD1 and IBD2 really do share one boundary rule.
4. `<prefix>X.seg` is emitted only under `--degree` and is visibly a half-finished writer
   in the reference — an 11-column header over 9 written fields, every row ending `\t\n`.
   It blocks exactly one case (`sexchr__ibdseg_degree2`), which also fails on `king.seg`,
   so it is worth nothing until §2 is closed.
