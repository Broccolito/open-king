# Fitting the IBD-segment rule from the golden corpus

**Status:** measurement + fit. No KING source was read; everything below comes from the
captured corpus, from re-running the reference binary on that same corpus with different
*flags*, and from the corpus generator's own simulated haplotypes. Scripts:
`tests/parity/fit/` (Python 3 + numpy only).

**Scope.** This document is about the *caller*: which stretches become segments and how long
they are. The aggregates downstream of it — `allsegs.txt`, the `InfType` table, the
`.kin`/`.kin0` layouts — are `docs/research/12-segment-aggregates.md`, which supersedes the
`InfType` rule used for the "InfType ok" counts below (they are computed with the older
one-clause FS test from `02-ibdseg.md`, so treat that one column as indicative only).

**Headline.** The denominator and the three reported columns are now closed questions. The
boundary convention is measured, not guessed: segments are **word-aligned**, refined by the
flanking word's **last IBS0**. The best rule found reproduces **626 / 982** `.seg` rows at
all four printed decimals with **0 missing** and **188 extra** pairs, mean |ΔPropIBD|
**0.00327** (committed engine: 625 / 982, same 188 extras, MAE 0.00488). On `bigish`, which
is 78 % of the corpus, MAE drops 0.00329 → **0.00145** and the median row is exact. The
residual is not noise: it is **±1 scan word** of ambiguity at segment ends, and it is
*anti-correlated* between the two columns — total span is unbiased (mean −0.05 words) while
IBD1 runs **+2.0 words** long and IBD2 **−2.0 words** short per pair.

---

## 1. Reproducing

```bash
# datasets (once)
python3 tests/parity/generate_corpus.py --outdir tests/parity/work/data

cd tests/parity/fit
python3 truth.py             # replay the simulator; must print 10/10 datasets reproduced
python3 fit.py baseline      # mirrors crates/king-core/src/ibdseg.rs exactly
python3 report.py            # the fitted rule, per dataset and per --seglength variant
```

Two probe sets are produced by re-running the reference on the **existing** corpus (no new
fixtures — only different flags on the same `.bed`):

```bash
KING="/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
D=tests/parity/work/data
for ds in nuclear threegen multifam dups missing monomorphic sexchr unrelated admixed bigish; do
  "$KING" -b $D/$ds.bed --ibs --prefix ibs_$ds            # MaxIBD2 per pair
  for L in 1 2 3 4 5 6 7 8 9 10; do
    "$KING" -b $D/$ds.bed --ibdseg --seglength $L --prefix ${ds}_L$L
  done
done
python3 boundary.py <that-dir>      # fit the word predicate / boundary convention
python3 fit_ibd2.py <that-dir>      # fit the IBD2 rule against MaxIBD2
```

## 2. Four measuring instruments

| instrument | what it gives | size |
| --- | --- | --- |
| `golden/ibdseg/*/king.seg` | two aggregates per pair | 982 rows |
| **`--ibs`'s `MaxIBD2` column** | **the length in bp of one single segment**, exact | 158 pairs |
| `--seglength 1..10` sweep | the length distribution of the called segments | 10 points/pair |
| `truth.py` | the **true** IBD state at every marker | every pair |

`MaxIBD2` is the important one and it was previously unused: `--ibs` prints
`MaxIBD2 <bp>.000` for every pair with kinship ≥ 2^−3.5 on a map with ≥ 100 Mb usable. It is
a single segment's length, so a candidate rule either reproduces it exactly or does not —
no averaging, no cancellation. Inverting it (search every marker interval whose span equals
the printed value) is what pinned the boundary convention in §4.

`truth.py` re-runs `generate_corpus.py`'s meioses while tracking which founder haplotype
each transmitted allele came from, adding **no** RNG draws, and asserts that the replayed
genotypes reproduce the committed `.bed` — 10/10 datasets, every marker. Truth is not a
rule (the reference only sees genotypes); it is used to tell "the reference is conservative
here" apart from "our caller is wrong here".

## 3. SOLVED — the denominator and the three columns

```
D        = sum over the AUTOSOMAL rows of <prefix>allsegs.txt of (pos[StopSNP] - pos[StartSNP])
IBD1Seg  = (bp called IBD1 and not IBD2) / D
IBD2Seg  = (bp called IBD2) / D
PropIBD  = IBD2Seg + IBD1Seg/2, in f64, formatted once
```

* `D` equals the console's `Total length of N chromosomal segments … is X Mb` on **all 10**
  datasets that emit a `.seg` (`report.py` re-derives both).
* It is `last − first` marker position per segment, not `+1` and not a marker count: a
  parent–offspring pair then reads **exactly** `1.0000 / 0.0000 / 0.5000`, which 314 of the
  316 `PO`-labelled rows do. (The two exceptions are sib pairs the `InfType` rule *labels*
  `PO`; their IBD1Seg is 0.9800 and 0.9007.)
* IBD1 excludes IBD2: a duplicate pair reads `0.0000 / 1.0000`.
* `PropIBD` is computed at full precision — **87 of 982** rows disagree with
  `round(printed IBD2 + printed IBD1/2, 4)`, so it cannot be derived from the printed columns.
* X-chromosome segments are excluded from `D` (they are listed in the same `allsegs.txt`).

The `--seglength` sweep additionally proves the two column definitions **interact**: on
`nuclear N_C1×N_C2`, going from `--seglength 8` to `9` *raises* IBD1 (167.3 → 179.0 Mb) while
lowering IBD2 (131.2 → 123.1 Mb). Dropping a short IBD2 segment hands its territory back to
IBD1, so IBD1 must be computed as a set difference and not as an independent scan.

## 4. MEASURED — segments are word-aligned, and where the boundary lands

Inverting `MaxIBD2` over **all** marker pairs (not just word-aligned ones) on `nuclear`: all
**6 of 6** longest-IBD2 segments resolve to an interval `[64u, 64v+63]` — start ≡ 0 (mod 64),
end ≡ 63 (mod 64); one pair additionally matched a second, non-aligned interval whose
interior is full of IBS0, i.e. a coincidence of spans. Across the whole corpus, 154 of the
158 `MaxIBD2` values localise to exactly one plausible (IBS0-free) such interval. The
boundary convention that fits them, and that also makes PO read 1.0000, is one rule for both
segment types:

```
run of "good" words [u..v] inside usable segment [lo..hi], complete words [w0..w1]

lo_seg = lo                                   if u == w0   (with the fringe rule below)
       = 64(u-1) + lastIBS0(word u-1) + 1     if word u-1 carries an IBS0
       = 64u                                  otherwise
hi_seg = hi                                   if v == w1   (with the fringe rule below)
       = 64(v+1) + lastIBS0(word v+1)         if word v+1 carries an IBS0
       = 64(v+1) + 63                         otherwise
length = pos[hi_seg] - pos[lo_seg]
```

Read that carefully: the segment **reaches into the word that ended it**, stopping on that
word's *last* IBS0 — and when the word that ended it holds no IBS0 at all (the IBD2 case,
where runs end on het-mismatch) it swallows that word whole. Evidence, over the 154 located
segments, given the observed start word:

| candidate | start reproduced | end reproduced |
| --- | ---: | ---: |
| good word = IBS0-free ∧ IBS1 ≤ 4, end = last good + 1 word | **154 / 154** | **136 / 154** |
| … IBS1 ≤ 3 | 148 | 128 |
| … IBS1 ≤ 5 | 147 | 136 |
| break only when two consecutive words are busy | 62 | 147 |

The **fringe rule** matters and is separately decisive. Where a run touches the first/last
*complete* word of a usable segment, the segment takes the usable segment's own end, pulled
in only by IBS0 markers among the segment's own markers in that incomplete word. Scoring
the whole corpus: fringe-aware **626** exact rows vs **581** if the run simply snaps to
`lo`/`hi`, vs **211** if the ordinary flanking-word formula is applied blindly and clamped
(that variant destroys every PO row).

## 5. FITTED — the IBD1 rule, on rows where IBD2 cannot interfere

823 of the 982 rows have `IBD2Seg == 0.0000`, so `IBD1Seg` alone determines the rule there —
no set-difference, no cross-talk. Fitting on those (`fit.py`-style scoring, IBD1 column only):

| IBD1 rule | exact / 823 | MAE |
| --- | ---: | ---: |
| **IBS0-free words, ≥ 2 complete words, boundary as §4** | **580** | **0.00196** |
| same, boundary = last-IBS0 left / *first*-IBS0 right | 312 | 0.00355 |
| same, boundary = first-IBS0 left / last-IBS0 right | 330 | 0.00344 |
| same, boundary = first-IBS0 both ends | 323 | 0.00250 |
| ≥ 1 word | 315 | 0.01361 |
| ≥ 3 words | 361 | 0.00510 |
| ≥ 1 word but run span ≥ 6 Mb | 573 | 0.00199 |
| tolerate 1 IBS0 per word | 325 | 0.00960 |
| absorb an isolated bad word (bridge) | 537 | 0.00240 |

So: **zero** IBS0 tolerance, **no** bridging, minimum run **2 complete words**, and both
boundaries anchored on the flanking word's *last* IBS0. The 313 PO rows are exact under any
of these, so the discriminating power is the other 510 rows.

**Caveat, stated rather than hidden:** every corpus dataset has ~50 kb marker spacing, so
"2 words" (≈ 6.4 Mb) and "run span ≥ 6 Mb" are nearly indistinguishable here (580 vs 573).
Separating them needs a map with markedly different density, which this corpus does not
contain.

## 6. FITTED — the IBD2 rule

A word joins an IBD2 run if it is IBS0-free **and** its IBS1 (het-vs-hom) count is ≤ t₂. The
two instruments disagree mildly on t₂:

* `MaxIBD2` start-boundary fit: t₂ = 4 (154/154 starts; 148 at t₂ = 3, 147 at t₂ = 5);
* whole-corpus `.seg` totals: t₂ = 0 → 626 exact, t₂ = 1 → 623, t₂ = 2 → 623, t₂ = 3 → 622,
  t₂ = 4 → 622. The spread is 4 rows: the corpus barely cares, `MaxIBD2` does.

Reproducing `MaxIBD2` *exactly* (length of the longest IBD2 segment, all 158 pairs) reaches
**51 / 158** with t₂ = 3–4 plus bridging, mean relative error 3.9 %; every alternative
(`end = last good word`, `start = one word back`) is worse. The residual is again ≈ one word.

Swapping in the committed engine's two-word het-break contingency rule instead of a per-word
count changes the corpus score by ≤ 2 rows (579–580 vs 581 at the same settings) — the
corpus cannot separate them, so the simpler per-word count is what the doc recommends.

## 7. The 188 extra rows — characterised, not fixed

Under both the committed engine and the fitted rule the corpus produces the **same 188**
extra pairs (182 in `bigish`, 3 `monomorphic`, 2 `unrelated`, 1 `multifam`). They are not a
separate inclusion threshold:

* the pair filter really is **`max segment length > 10 Mb`** on the *reported* (extended)
  length. Every alternative statistic is much worse — core-only span: 0 extra but **270
  missing**; core + right extension: 2 extra / 257 missing; core + left extension: 0 / 260;
  total IBD length or segment count: no separation at all;
* the filter is **one-sided and exact in one direction**: of the 763 pairs the reference
  reports in `bigish`, not one has a computed longest segment ≤ 10.00 Mb, and 862 pairs sit
  in the 9–10 Mb band with none reported. Whatever is wrong, it is never that we under-state;
* our lengths are right where they can be checked. For the 243 reported `bigish` pairs where
  the rule calls exactly one segment, `IBD1Seg × D` is that segment's length to ±0.13 Mb:
  **76.5 %** agree within the printing resolution and the largest *excess* over the whole set
  is **+0.11 Mb** (the tail runs the other way, down to −16 Mb, i.e. segments we miss entirely);
* and yet the borderline pairs cannot be told apart. Restricting to pairs whose longest
  segment has a **2-word core** and measures > 10 Mb — 255 reported, 182 extra — the two
  groups match quartile for quartile on every quantity computed: length (median 10.45 vs
  10.41 Mb), left extension (1.26 vs 1.25 Mb), right extension (2.95 vs 3.00 Mb), SNPs in the
  segment (210 vs 209), HetHet (30 vs 32), IBS1 (96 vs 92), pair kinship from `--ibs`
  (−0.0007 vs −0.0014), number of called segments, and the parity of the start word mod 2 and
  mod 4. There is no feature here that separates them;
* PropIBD is not the discriminator either: the reference reports pairs all the way down to
  PropIBD 0.0020, which is exactly where the extras sit.

So the honest statement is: **for ~42 % of two-word runs the reference's own length falls the
other side of 10 Mb for a reason invisible in the pair's genotypes as we summarise them.**
The earlier guess — that the extras are the ones whose flanking word is IBS0-sparse, so that
we extend too far — is *wrong*, and was an artefact of comparing the largest extras against
the smallest reported pairs; matched on core size the extensions are identical.

**This remains the single highest-value open problem: it gates ~97 parity cases, and it is
now known to live in the *calling* of two-word runs, not in the length measurement or in the
pair-level threshold.**

## 8. Ground truth: where the *reference* is the odd one

`truth.py` gives the true π₁/π₂ per pair. Comparing reference against truth:

| dataset | samples | mean |ref − truth| |
| --- | ---: | ---: |
| bigish, admixed, multifam, dups, sexchr, threegen | 10–200 | 0.01–0.03 |
| monomorphic | 12 | 0.06 |
| **nuclear, missing** | **6** | **0.15–0.40** |

On `nuclear N_C1×N_C3` truth is π₁ = 0.5326 / π₂ = 0.2663, the reference prints
0.1057 / 0.3144, and *both* the committed engine (0.5647 / 0.2528) and the fitted rule land
near truth. The reference simply loses ~190 Mb of genuinely IBD1 genome on that pair, and
its own `--ibs` kinship for it is 0.2669 — inconsistent with its own PropIBD of 0.4201. On
`monomorphic` it reports `IBD2Seg 0.0000` for two pairs whose truth is 0.41 and 0.26.

Consequence for anyone continuing this work: **do not tune against `nuclear`, `missing` or
`monomorphic`.** They are 42 of 982 rows, they are where every rule scores worst, and the
reference's answers there are not reconstructions of the underlying IBD — fitting them would
be fitting a defect. All the fits above were checked to be driven by the other 940 rows.

## 9. The best-fitting rule, stated

```
WORD = 64 markers of the global retained-autosome grid
for each usable segment [lo..hi] with complete words [w0..w1]:
    good1(w) = IBS0count(w) == 0
    good2(w) = good1(w) and IBS1count(w) <= t2          # t2 = 0 by .seg totals, 4 by MaxIBD2
    IBD1 runs = maximal runs of good1, length >= 2 words
    IBD2 runs = maximal runs of good2, length >= 1 word
    boundaries: §4 (flanking word's last IBS0; usable-segment ends via the fringe rule)
    drop any segment shorter than --seglength (default 3 Mb)
IBD2 bp = sum of IBD2 segment lengths
IBD1 bp = sum of IBD1 segment lengths minus their overlap with IBD2 segments
report the pair iff its longest segment (either kind) > 10 Mb
IBD1Seg = IBD1bp/D, IBD2Seg = IBD2bp/D, PropIBD = IBD2Seg + IBD1Seg/2
```

Scorecard (`report.py`), against all 982 captured rows:

| dataset | rows | exact | IBD1 ok | IBD2 ok | extra | missing | InfType | MAE | worst |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| nuclear | 14 | 8 | 8 | 8 | 0 | 0 | 14 | 0.04275 | 0.1753 |
| threegen | 39 | 23 | 23 | 36 | 0 | 0 | 37 | 0.00737 | 0.0571 |
| multifam | 104 | 58 | 58 | 87 | 1 | 0 | 100 | 0.00366 | 0.0299 |
| dups | 3 | 2 | 2 | 2 | 0 | 0 | 3 | 0.01146 | 0.0344 |
| missing | 14 | 8 | 8 | 8 | 0 | 0 | 12 | 0.03993 | 0.1550 |
| monomorphic | 14 | 11 | 11 | 12 | 3 | 0 | 13 | 0.00959 | 0.0971 |
| sexchr | 14 | 8 | 8 | 8 | 0 | 0 | 14 | 0.00501 | 0.0180 |
| unrelated | 1 | 1 | 1 | 1 | 2 | 0 | 1 | 0.00001 | 0.0000 |
| admixed | 16 | 12 | 12 | 13 | 0 | 0 | 16 | 0.00274 | 0.0142 |
| bigish | 763 | 495 | 495 | 649 | 182 | 0 | 763 | 0.00145 | 0.0160 |
| **ALL** | **982** | **626** | **626** | **824** | **188** | **0** | **973** | **0.00327** | **0.1753** |

Same rule against the `--seglength` captures: **626 / 982** exact at `--seglength 5`,
**731 / 982** at `--seglength 10`, 188 extras and 0 missing in both — the reporting floor
and the fixed 10 Mb pair filter are modelled correctly.

Committed engine for comparison: 625 exact, IBD1 625, IBD2 822, 188 extra, 0 missing,
InfType 973, MAE 0.00488, worst 0.1679. The fitted rule is worth **+2 IBD2 columns and a
33 % MAE reduction** (55 % on `bigish`) for two changes: word-aligned/`+1`-word IBD2 ends
instead of the het-break geometry, and the per-word IBS1 count instead of the two-word
contingency test. It does **not** move exact-row parity, and it does not touch the 188
extras.

## 10. Ruled out

* IBD2 boundaries at marker resolution — 154/158 `MaxIBD2` values are word-aligned.
* A separate pair-inclusion threshold on PropIBD, on total IBD length, or on segment count.
* An IBS0 tolerance per word (t₁ ≥ 1) for IBD1 — halves the exact-row count.
* Bridging an isolated bad word in IBD1 runs (580 → 537).
* Whole-run classification (a run is *either* IBD1 *or* IBD2): the assignments consistent
  with the reference's totals put runs with IBS1 rate 0.35 in IBD2 and runs with rate 0.02
  outside it. Reference segments are strict sub-intervals of the IBS0-free runs.
* Any rule reading the `.bim` cM column, or the pair's kinship, into the segment scan.

## 11. Next experiments, in the order they are worth doing

1. **Two-word runs** (§7). The 437 borderline `bigish` pairs are a ready-made labelled set
   (255 called, 182 not) that no per-segment feature separates. The next thing to try is a
   feature the scan can see but this analysis has not summarised: the *positions* of the IBS0
   markers inside the two core words and their immediate neighbours (not just counts), and
   whether the run is adjacent to a second, sub-threshold run. `analyze.py`/`report.py`
   already assemble the two groups; a decision-tree over per-marker patterns would find it if
   it exists.
2. **`MaxIBD2` as a regression target.** 158 exact segment lengths, currently 51 reproduced.
   Getting that to 158 almost certainly gets the IBD2 column to parity, and it is a far
   sharper gradient than the `.seg` totals.
3. **`--roh`'s `MaxROH`/`F_ROH`** is the same machinery on one sample instead of two, with a
   documented 5 Mb floor; an inbred fixture would expose the run rule with no pair geometry
   at all. (This corpus is outbred: every `MaxROH` is 0.0.)
4. **Density.** One fileset at 10 kb and one at 200 kb spacing would separate "2 words" from
   "6 Mb" (§5) — the one place this corpus is structurally blind.
