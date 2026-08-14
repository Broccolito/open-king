# Attack 4 — the IBD2 caller is the whole remaining `.seg` residual

**Status:** measurement. No KING source was read; every number is either a score against
the captured corpus or a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) run on
the parity corpus or on constructed filesets.

**Outcome, stated first: no rule changed.** Every variant tried scored the same as or
worse than the committed engine on all four graders, so `crates/king-core/src/ibdseg.rs`
carries only a documentation update. What this attack did produce is a much sharper
localisation of the fault, three new instruments, and four measurements that rule out the
model the current geometry is built on. §8 lists what is still unexplained.

---

## 1. Headline

| claim | evidence |
| --- | --- |
| **The IBD1 caller is finished.** Of the 823 corpus rows where the reference prints `IBD2Seg 0.0000`, open-king reproduces **819** on both estimate columns. All 4 exceptions are `monomorphic`, the fileset `docs/PARITY.md` §5.1 already documents as ungradeable. | §2 |
| **The IBD2 caller is wrong on essentially every row that carries IBD2.** Of the 159 rows where the reference prints `IBD2Seg > 0`, **3** have both columns right. 154 have both wrong. | §2 |
| The `.seg` and `--ibs` columns disagree about the *direction* of the error: `.seg`'s `IBD2Seg` is **too low** on 121 of 159 rows, while `--ibs`'s `Pr_IBD2` — the same calls under the word-aligned ruler — is **too high on 150 of 158** (too low on 1), mean `+0.0342` of the genome. | §4 |
| Those two facts are **not simultaneously satisfiable** by the committed geometry: on `nuclear N_C1/N_C2` the reference's `.seg` total exceeds its own `Pr_IBD2` total by 22.64 Mb, and the largest gain the "word-aligned plus usable-segment fringe" model can produce on that fileset is **13.58 Mb**, over *all five* usable segments and *all* pairs. | §5 |
| A constructed fixture confirms `10-…fixtures` §3 independently: an IBD2 run of `W` clean words bounded by all-IBS0 words reports exactly `64W - 1` marker intervals for `W` = 1,2,3,4,6,8. The committed rule reports `64(W+1) - 1`. | §6 |
| A second fixture — chr2 entirely IBD2, one IBS0 forced into a chosen word — is reproduced by **no** endpoint rule in the parameter space, including the one that fixes the first fixture. | §6.2 |

---

## 2. Where the residual lives

`tests/parity/fit/engine.py` is a Python mirror of the committed engine, asserted equal to
the Rust binary on all 982 `.seg` rows and all 861 diffable `MaxIBD2` values by
`check_mirror.py`. Splitting its scorecard by whether the *reference* reports any IBD2:

| reference row | rows | both columns exact | IBD1 only | IBD2 only | both wrong |
| --- | ---: | ---: | ---: | ---: | ---: |
| `IBD2Seg == 0.0000` | 823 | **819** | 0 | 0 | 4 |
| `IBD2Seg > 0` | 159 | **1** | 2 | 2 | 154 |

The four failures in the first row are `monomorphic` `P_C1/P_C2`, `P_C1/P_C3`,
`P_C1/P_C4`, `P_C3/P_C4` — the fileset whose own `--kinship` (0.3384 for `P_C1/P_C2`)
contradicts its own `.seg` (φ = 0.2450 from the printed columns) and whose reference
answers are nowhere near the pedigree the generator built.

So: **every one of the 36 failing `ibdseg/*__ibdseg` parity cases is failing because of
the IBD2 caller**, and the IBD1 boundary work is done. The corollary matters for
prioritising: `.seg` exact-row count is a *bad* gradient for IBD2 work, because it is
dominated by the 823 rows IBD2 cannot touch — every variant in §7 scores 705 ± 1.

## 3. The instruments

Four graders, in increasing sharpness, all in `tests/parity/fit/`:

| instrument | what it grades | resolution | script |
| --- | --- | --- | --- |
| `.seg` `IBD1Seg`/`IBD2Seg`/`PropIBD` | two totals per pair | 4 dp | `engine.score_seg` |
| ~~`--ibs` `Pr_IBD2`~~ | the word-aligned IBD2 **total**, no length filter | 4 dp | `engine.score_pr` |
| ~~`--ibs` `MaxIBD2`~~ | one exact segment length per pair, word-aligned | 1 bp | `invert.py` |
| **`--seglength` bisection** | the **length of every individual segment under 10 Mb** | 1 bp | `seglen_probe.py` |

> **Two of these four are spent.** `docs/research/16-segment-extension.md` solved the
> `--ibs` caller, so `Pr_IBD2` and `MaxIBD2` are now exact — and a grader that every
> candidate passes grades nothing. Worse, they were never grading *this* function: they are
> produced by `Scan::ibd2_words`, a different caller over the same masks (§3 of `15-…`).
> What is left for `.seg` work is a bad total-based gradient and one sharp instrument that
> costs a reference invocation per bisection step. §10 of `16-…` argues the next campaign's
> first job is building a `.seg`-native canvas out of opposite homozygotes.

### 3.1 `MaxIBD2` inversion — now 157 of 158 localised

`invert.py` searches every `(u, e)` word pair inside every usable segment for a span equal
to the printed `MaxIBD2`, keeping only intervals an IBD2 caller could produce (opening on
an IBS0-free word, IBS0 confined to the last two). **157 of 158 targets localise to
exactly one interval** — up from the 154 of `11-segment-rule-fit.md` §4, because the
plausibility filter is now stated in terms of the extension words rather than requiring the
whole interval to be IBS0-free. The one that does not localise is `dups MZ_1/MZ_2`.

A marker-level search (all marker pairs, not just word-aligned ones) over the six
disagreeing pairs finds no non-word-aligned alternative that is IBS0-free, so the
word-aligned reading of `MaxIBD2` survives the wider search.

### 3.2 `--seglength` bisection — the new instrument

`--seglength` accepts a float in Mb and is honoured to the base pair (`--seglength
6.250001` prints `Minimum segment length is set as 6250001 bp`), and a segment is kept iff
its length clears it. `IBD2Seg(L)` is therefore a step function of `L` whose jumps sit at
individual segment lengths, so an adaptive bisection over the *vector* of printed values
recovers segment lengths one at a time — 20 ms per reference invocation, a few hundred
invocations per dataset.

Two limits found while building it, both worth recording on their own:

* **`--seglength` above 10 Mb is silently ignored.** At 10.5, 12, 20, 300 the
  "Minimum segment length is set as …" line disappears and the output is identical to the
  default 3 Mb. open-king already reproduces this (it prints
  `KING supports minimum segment length from 1 to 10 Mb at the moment.`), so it is not a
  parity gap — but it caps this instrument at 10 Mb.
* **Dropping a short segment can make `IBD2Seg` go *up*.** On `dups MZ_1/MZ_2`,
  `IBD2Seg` climbs 0.9223 → 0.9385 → 0.9547 → 0.9715 → 0.9877 across four jumps while
  `IBD1Seg` falls 0.0436 → 0.0300 → 0.0145 → 0.0145 → 0.0000. Both columns jump at the
  *same* four values of `L`. So the reference's `--seglength` filter drops a segment
  **before** it can clip its neighbour — which the committed engine already models — and,
  less comfortably, dropping a segment that shows up in the IBD1 column releases territory
  that shows up in the IBD2 column. The two passes are not independent in the reference.

## 4. The two rulers disagree about the sign of the error

Over the 158 pairs the reference grades with a non-zero `MaxIBD2`:

| grader | exact | bias |
| --- | ---: | --- |
| `MaxIBD2` (longest IBD2 segment, word-aligned) | **145 / 158** | — |
| `Pr_IBD2` (total IBD2, word-aligned, unfiltered) | **7 / 158** | **+0.0342** of the genome; too high on 150, too low on 1 |
| `.seg` `IBD2Seg` (total IBD2, segment-end ruler, 3 Mb floor) | 3 / 159 | too low on 121, too high on 35 |

The longest segment is right and the total is not, in both directions at once. Whatever is
wrong is in the *other* segments, and it changes sign between the two rulers.

## 5. The `.seg` ruler is not "word-aligned plus fringe" — a counting proof

The committed model says a `.seg` IBD2 call is the word-aligned interval `[64u, 64e+63]`,
widened to the usable segment's own first/last marker when it touches `w0`/`w1`. The
widening is bounded: per usable segment it is at most
`pos[64·w0] - pos[lo]` plus `pos[hi] - pos[64·w1+63]`.

Summed over **every** usable segment of `nuclear`, that ceiling is **13 583 989 bp**
(0.0272 of `D`). The reference's own two columns for `N_C1/N_C2` are `IBD2Seg 0.2626` and
`Pr_IBD2 0.2173`, a difference of **0.0453 of `D` = 22.64 Mb** — 1.7× the ceiling, on one
pair. The model cannot produce it.

`threegen` `TG_P1/TG_P2` is the same story with room to spare in the ceiling (62.1 Mb) but
the same shape: the reference gains 73.8 Mb from `.seg`'s ruler where open-king gains 21.5.

So either the `.seg` IBD2 interval is materially wider than the word-aligned one on
*interior* calls, or `Pr_IBD2` is measured over a narrower interval than `MaxIBD2` is.
Both were tried in §7; neither scores.

## 6. Two fixtures, and what they say

Run from `docs/research/fixtures/` (they drive the reference directly):

```bash
python3 ibd2end.py      # §6.1
python3 ibd2gap.py      # §6.2
python3 seglen_edge.py  # §6.3
```

### 6.1 An IBD2 run is not extended — reproduced

`ibd2end.py` carves a `W`-word IBD2 block out of the solid-IBS0 canvas of `rig2.py`:

| `W` | reported IBD2 marker intervals | `64W-1` | `64(W+1)-1` |
| ---: | ---: | ---: | ---: |
| 1 | **63** | 63 | 127 |
| 2 | **127** | 127 | 191 |
| 3 | **191** | 191 | 255 |
| 4 | **255** | 255 | 319 |
| 6 | **383** | 383 | 447 |
| 8 | **511** | 511 | 575 |

Six for six on `64W - 1`, independently reproducing `10-segment-rule-fixtures.md` §3. The
committed rule ends an IBD2 call on the flanking word's **last** IBS0, which on this canvas
is bit 63, and so reports `64(W+1) - 1` — one whole word too many. The rule that fits is
"stop one marker before the flanking word's **first** IBS0, and take the whole word when it
has none", which is `Params(ibd2_right=("first", -1))` in the harness.

That variant was scored on the corpus: `.seg` exact 705 (unchanged), `IBD2Seg` columns
822 → **820**, `MaxIBD2` 145 (unchanged), MAE 0.00138 (unchanged). It is fixture-correct
and two corpus rows worse, and §6.2 shows it is not the missing rule either, so it was
**not** committed — swapping one wrong endpoint rule for another wrong one is churn.

### 6.2 …but a single IBS0 in the flanking word is reproduced by nothing

`ibd2gap.py` makes chr2 (640 markers = 10 complete words, no fringe) **entirely IBD2**, so
the pair has no natural opposite homozygote anywhere on it, then forces exactly one IBS0
at a chosen word and bit. The reported IBD2 length is then a pure statement about the two
endpoints either side of that one marker.

| poke | reported IBD2 marker intervals |
| --- | ---: |
| no poke | 639 (the whole chromosome) |
| word 0, bit 20 | 575 |
| word 1, bit 20 | 510 |
| word 5, bit 20 | 510 |
| word 5, bit 0 | 510 |
| word 5, bit 63 | 510 |
| word 9, bit 20 | 575 |

Three things fall out. **The bit position does not matter at all** — bits 0, 20 and 63 of
word 5 give the identical answer, so the endpoints either side of an interior IBS0 are
word-quantised, not marker-refined. **An IBS0 in the first or last word costs exactly one
word** (639 − 64 = 575). **An IBS0 in an interior word costs two words and one marker**
(639 − 129 = 510), whatever word it is in.

No rule in the harness produces 510. "Own words only" predicts 574; the committed rule
predicts 595; `("first", -1)` predicts 594. The read-back resolution here is ±6 marker
intervals (one ulp of the printed `IBD2Seg` is 12.8 markers on this fixture), which is
nowhere near enough to explain a 64-marker gap. **This is the sharpest single unexplained
measurement in the project**, and it is cheap to iterate on: one reference invocation per
candidate.

### 6.3 The `--seglength` comparison is not simply inclusive

`10-segment-rule-fixtures.md` §4 records the floor as `length >= seglength`. Bisecting
`--seglength` to the base pair against blocks of known exact length says the truth is
finer:

| block | reported intervals | first `--seglength` that drops it | that, in intervals |
| --- | ---: | ---: | ---: |
| IBD1, 1 clean word | 127 | 6 300 000 | 126.00 |
| IBD1, 2 clean words | 191 | 9 550 001 | 191.00 + 1 bp |
| IBD2, 2 clean words | 127 | 6 350 000 | 127.00 |
| IBD2, 3 clean words | 191 | 9 550 000 | 191.00 |

The 2-word IBD1 block is inclusive on its reported length exactly as documented. The two
IBD2 blocks flip one base pair earlier than that, and the 1-word IBD1 block flips a whole
marker earlier. So the length the filter compares is not always the length the totals
report.

**Corpus impact at the captured floors: none at 3 Mb.** No call open-king makes on any
corpus dataset has a length in `[3 000 000, 3 000 000 + one marker gap)`; one IBD1 call
lands in that band above 5 Mb, and 63 calls land in it above 10 Mb. Since the primary
captures use 3 Mb this cannot be part of the current residual, and it was left alone.

## 7. What was tried and did not work

All scored with `python3 tests/parity/fit/sweep2.py`, which prints all four graders at
once. `seg` is exact rows of 982, `Max` of 158, `Pr` of 158.

| variant | seg | both | Max | Pr | Pr bias |
| --- | ---: | ---: | ---: | ---: | ---: |
| **committed** | **705** | **820** | **145** | **7** | **+0.0342** |
| no bridging of a lone dirty word | 705 | 820 | 139 | 7 | +0.0340 |
| IBD2 not extended (`ext = 0`) | 705 | 820 | 63 | 3 | +0.0122 |
| `ext = 0` with `--ibs` measuring one word wider | 705 | 820 | 145 | 7 | +0.0342 |
| IBD2 borrows IBD1's endpoint refinement | 705 | 820 | 142 | 7 | +0.0340 |
| tail-snap window 0 / 1 / 3 (committed: 2) | 705 | 820 | 142 | 7 | ±0.001 |
| dirty threshold IBS1 ≥ 3 / 4 / 6 / 7 / 8 / 10 | 659–705 | 771–820 | 106–142 | 3–7 | — |
| IBD2 minimum run 2 words / 3 words | 705 / 706 | 820 / 821 | 145 | 12 / 24 | +0.027 / +0.016 |
| IBD2 start refined by the previous word's last IBS0 | 705 | 820 | 145 | 7 | +0.0349 |
| clip against the previous call after the length test | 705 | 820 | 145 | 7 | +0.0342 |
| every ±1-marker variant of all five refinement rules | ≤ 705 | ≤ 820 | 145 | 7 | — |

Two of these deserve a note.

* **Raising the IBD2 minimum run length is the only thing that moves `Pr_IBD2`** (7 → 24 at
  three words, bias +0.0342 → +0.0162). It is also flatly contradicted by the fixture — a
  single clean IBD2 word reports its 63 intervals (§6.1, `W = 1`) — so it is a fit, not a
  rule, and was not taken.
* **The five marker-level refinement rules are at a strict local optimum.** Sweeping
  `{last, first} × {−1, 0, +1}` independently for the IBD1 right end, the IBD1 left end,
  the IBD2 right end and both usable-segment fringe rules, the committed choice is the best
  or tied-best on every axis; the nearest alternatives cost 20 to 400 exact rows. The
  ±1-marker residual is therefore **not** an off-by-one in the refinement.

`PropIBD`'s formula was also re-checked, since 115 rows print both estimates correctly and
`PropIBD` wrong. Six arrangements of `IBD2Seg + IBD1Seg/2` — including one-division forms
like `(2·ibd2 + ibd1) / (2·D)` — give **identical** results on all 982 rows. Those 115 rows
are a genuine sub-printing-resolution difference in the called base pairs, of order one
marker interval, and 92 of them are `bigish`.

## 8. What is still unexplained, precisely

> **Items 3, 4 and 5 below were `--ibs` observations, and all three are now explained** —
> they are chunk refusals under the confirmation scan of
> `docs/research/16-segment-extension.md` §8.1, which reproduces `MaxIBD2` 158/158 and
> `Pr_IBD2` 158/158. They are kept here as written because the *reasoning* about them (that
> no per-word statistic separates the exceptions from the rest) was correct and is the
> reason the answer turned out to be a rule about chunks rather than about words. Items 1
> and 2 are `.seg` observations and remain open.

1. **§6.2's 510.** An IBS0 in an interior word of an otherwise fully-IBD2 chromosome costs
   two words and one marker of IBD2 coverage, and costs one word when it is in the first or
   last word of the usable segment. Nothing in the current model produces that.
2. **The `.seg`/`Pr_IBD2` sign inversion** (§4, §5). The reference's `.seg` IBD2 total
   exceeds its word-aligned total by more than the fringe model allows.
3. **The four declined right extensions.** `bigish` `B06_C1/B06_C2`, `B09_C1/B09_C2`,
   `B26_C1/B26_C2`, `B26_C1/B26_C4`: the reference's interval ends on the run's own last
   clean word where 91 other interior cases end one word later. `endfit.py` tabulates
   `b − v` against IBS1 at `v+1` and `v+2`, IBS0 at `v+1`, informative count at `v+1`, run
   length, distance to `w1`, whether the run was bridged, and the gap to and length of the
   next run. **None separates them**: the four sit inside the interquartile range of the 91
   on every column. (The three `db = 2` cases and all 48 `db = 0` cases at a usable-segment
   end *are* explained — by the tail-snap rule, 48/48 and 3/3.)
4. **Two denied bridges.** `missing M_C2/M_C3` word 133 (IBS1 = 8) and `bigish
   B11_C2/B11_C3` word 40 (IBS1 = 29) are *not* absorbed by the reference, while 14 other
   lone dirty words with clean neighbours are — including IBS1 = 27, 23, 22, 22, 20. No
   per-word statistic separates the 2 from the 14. The one visible correlate is the length
   of the run to their **left** (3 and 2 words, against 8–14 for the absorbed ones), on 16
   observations — far too few to fit a threshold to, and deliberately not fitted.
5. **`dups MZ_1/MZ_2`'s `MaxIBD2`** (41 542 807) localises to no plausible interval at all.

## 9. Reproducing everything here

```bash
cd tests/parity/fit
python3 check_mirror.py            # the mirror is the Rust engine, to the byte
python3 engine.py                  # the baseline scorecard, all graders
python3 sweep2.py                  # §7's table
python3 invert.py                  # §3.1: localisation census + every MaxIBD2 miss
python3 endfit.py                  # §8.3: the contingency tables
python3 segdiff.py dups -v         # row-level .seg diff with the calls behind it
python3 seglen_probe.py nuclear    # §3.2: per-segment lengths under 10 Mb
python3 seglen_invert.py nuclear   # those lengths, back to marker endpoints

cd ../../../docs/research/fixtures
python3 ibd2end.py                 # §6.1
python3 ibd2gap.py                 # §6.2
python3 seglen_edge.py             # §6.3
```

`engine.py`, `sweep2.py`, `invert.py`, `endfit.py` and `segdiff.py` need no reference
binary; the rest drive it (path at the top of `fixlab.py`, or `$KING`).
