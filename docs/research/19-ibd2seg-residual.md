# The `.seg` residual: the fringe clause, measured

**Status: measured, validated out of sample, committed to Rust.** This closes the 86 wrong
`IBD2Seg` rows `docs/research/18-ibd1-caller.md` left behind, and — more importantly —
converts the clause that closes them from a **corpus fit into a bisection**. It also
relocates what is left: at the default floor the `.seg` residual is no longer in either
estimate column, it is in `PropIBD`.

> **§1-§8 are current; §9 and §10 are superseded by `20-seg-writer.md`.** This document's
> diagnosis of the leftover `PropIBD` rows — that they are a sub-ulp shortfall in the IBD1
> pass — is wrong. They are a *writer* rule: `.seg` computes `PropIBD` from its own printed
> columns. The fringe clause derived in §1-§7 is unaffected and remains committed. See the
> banner above §9.

No KING source was read. Every rule below is a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) run on
filesets built for the purpose, or a score against the captured parity corpus.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| A **new instrument**: a canvas whose usable segment does *not* start on a word boundary. `17-…` §5 and `18-…` both say in as many words that their rigs cannot see a segment edge; this one can, by spending markers instead of words — chromosome 1 is shortened by `f` markers so chromosome 2 opens `f` markers before a word boundary, while the painted words stay word-aligned. | §0 |
| The partial word beyond a segment's word grid **is** scanned, marker by marker, and an IBD2 call stops **one marker short of the nearest het-vs-hom mismatch** in it. Bisected at 16 positions on each side; `off=0`, `extend` and `none` are all excluded at every one. | §1, §3 |
| **An opposite homozygote in a fringe does not stop an IBD2 call.** One anywhere in a whole word *inside* the grid disqualifies that word outright, so this is a real asymmetry, and it is the clause `17-…` §5 would never have guessed. Missing calls and A1A1/A1A1 do not stop it either; het-vs-A1A1 does. | §2 |
| The reach **snaps**: once an end lands on the grid's own edge the marker scan takes over, whether that moves the end out or pulls it back in. Measured on a run that opens at the *second* complete word, where §5-as-written predicts a stop 24 markers short of what the reference prints. | §4 |
| The **IBD1** fringe is the exact mirror and the committed rule is already right: an IBD1 call stops one marker past the **last opposite homozygote** in a head fringe and one marker before the **first** in a tail fringe, and a het-vs-hom mismatch does *not* stop it. Each pass stops at its own breaking marker. | §5 |
| Out of sample on **504** fresh canvases — two unused random seeds and an exhaustive 3-word composition sweep crossed with six fringe shapes — the rule reproduces the reference on **504**. The Rust port reproduces it on **416**. | §6, §7 |
| Corpus, default floor: `IBD2Seg` exact **982 of 982** (was 896), whole rows **806** (was 747), mean `PropIBD` error **0.000023** (was 0.000067), worst row **0.0001** (was 0.0042). Parity **436 of 480** (was 408). | §8 |
| What is left is **not** `IBD2Seg`. 176 rows still differ, every one of them on `PropIBD` alone with both estimate columns exact; **115 of those carry no IBD2 at all on either side**, so they are the IBD1 pass's, sub-ulp, and short by about one marker gap. Two candidate explanations are measured dead. | §9, §10 |

---

## 0. The instrument: a segment that does not start on a word boundary

`docs/research/fixtures/segcanvas.py` builds chromosome 1 with `64·nw1` markers and
chromosome 2 with `64·nw2`. Both chromosomes therefore begin at a marker index divisible
by 64, every usable segment it makes is word-aligned at both ends, and the partial word
beyond a segment's grid — which the engine calls the *fringe* — never exists.
`ibd1canvas.py` imports that class unchanged, so the same is true there. Both campaigns
recorded the gap and took the clause from the corpus instead (`17-…` §5: *"The canvas
cannot see this (its chromosome 2 is word-aligned); the corpus can"*). A clause chosen to
make 982 corpus rows come out right is a fit, and this project does not keep fits.

`docs/research/fixtures/fringecanvas.py` closes it. The trick is to spend **markers**, not
words:

    chr1   64*nw1 - f markers      the carrier, so the pair still earns a `.seg` row
    chr2   f + 64*nw2 + t markers  f head-fringe markers, nw2 complete words, t tail

Because chromosome 1 is short by exactly `f`, chromosome 2 opens `f` markers before a word
boundary. Its usable segment is then

    lo = 64*nw1 - f          first_word = nw1            head fringe = chr2 [0, f)
    hi = 64*(nw1+nw2) + t-1  last_word  = nw1+nw2-1      tail fringe = the last t markers

and the `nw2` painted words stay exactly word-aligned to the global grid, so every rule
`17-…` and `18-…` measured still applies to them untouched. Only the two fringes vary, and
they are painted one marker at a time.

`nw1 = 6` keeps chromosome 1 over `MIN_WORDS` once it has lost its `f` markers; chromosome
2's uniform spacing is chosen by the same rule as the other rigs, so one ulp of the printed
`%.4lf` is about an eighth of a marker gap and the column reads back the number of chr2
marker intervals called, exactly.

### 0.1 The controls

```
  no fringe, 8 clean words then walls
    nseg 2 (want 2)   D 100617000   Dref 100617000   one ulp = 0.117 markers
  IBD2Seg                         511.041   word-aligned=511   [word-aligned]
```

Two usable segments, `D` identical to the reference's own `allsegs.txt` total to the base
pair, and a word-aligned call over 8 words reading `64·8 − 1 = 511` to four hundredths of a
marker. The ruler is sound.

### 0.2 Is the fringe reached at all?

```
  head= 1 clean                   512.012   extend=512  none=511   [extend]
  head= 8 clean                   518.993   extend=519  none=511   [extend]
  head=32 clean                   543.005   extend=543  none=511   [extend]
  head=63 clean                   573.979   extend=574  none=511   [extend]

  tail= 1 clean                   511.946   extend=512  none=511   [extend]
  tail= 8 clean                   518.947   extend=519  none=511   [extend]
  tail=32 clean                   543.033   extend=543  none=511   [extend]
  tail=63 clean                   573.945   extend=574  none=511   [extend]
```

Yes, on both sides, and all the way to the segment's own first and last marker when the
fringe carries nothing. `none` — the rule that never leaves the word grid — is dead at
eight independent widths. This much `17-…` §5 had right.

### 0.3 What the rig cannot see

The head-fringe word is **shared**: its low bits are chromosome 1's last markers, its high
bits are chromosome 2's first. Each segment owns only its own side. Because the carrier's
markers all sit at bits *below* the fringe, a stop computed from them can never land past
`seg.lo`, so whether KING masks the other chromosome's bits out of that word is
**indistinguishable here** — both readings give the same answer. That is a limitation of
the rig, stated rather than papered over. The engine masks.

---

## 1. The head fringe, bisected

A clean fringe of 32 markers with a single het-vs-hom mismatch at position `q`, swept:

| mismatch at | reference | `off=1` | `off=0` | `extend` | `none` |
| --- | ---: | ---: | ---: | ---: | ---: |
| head[ 0] | **542.051** | 542 | 543 | 543 | 511 |
| head[ 2] | **540.023** | 540 | 541 | 543 | 511 |
| head[ 4] | **537.995** | 538 | 539 | 543 | 511 |
| head[ 6] | **535.967** | 536 | 537 | 543 | 511 |
| head[ 8] | **534.058** | 534 | 535 | 543 | 511 |
| head[10] | **532.030** | 532 | 533 | 543 | 511 |
| head[12] | **530.002** | 530 | 531 | 543 | 511 |
| head[14] | **527.974** | 528 | 529 | 543 | 511 |
| head[16] | **525.946** | 526 | 527 | 543 | 511 |
| head[18] | **524.038** | 524 | 525 | 543 | 511 |
| head[20] | **522.010** | 522 | 523 | 543 | 511 |
| head[22] | **519.982** | 520 | 521 | 543 | 511 |
| head[24] | **517.954** | 518 | 519 | 543 | 511 |
| head[26] | **516.045** | 516 | 517 | 543 | 511 |
| head[28] | **514.018** | 514 | 515 | 543 | 511 |
| head[30] | **511.990** | 512 | 513 | 543 | 511 |

Sixteen for sixteen on `off=1`, and `off=0` is excluded at every one of them. The call
starts **one marker past** the mismatch:

    left = (last breaking marker in [lo, 64*w0 - 1]) + 1,   else lo

---

## 2. Which markers break a fringe — and which do not

Same rig, 32-marker head fringe, one or two markers changed:

| head fringe | reference | reading |
| --- | ---: | --- |
| mismatches at 5 **and** 20 | **522.010** | stops after 20 — the **last** one wins |
| **IBS0 at 20 only** | **543.005** | `extend` — an opposite homozygote **does not break** |
| IBS0 at 20, mismatch at 5 | **537.040** | stops after 5 — the IBS0 is ignored |
| mismatch at 20, IBS0 at 5 | **522.010** | stops after 20 |
| het-vs-A1A1 at 20 | **522.010** | breaks, exactly like het-vs-A2A2 |
| both missing at 20 | **543.005** | does not break |
| A1A1/A1A1 at 20 | **543.005** | does not break |

The middle row is the result. **One opposite homozygote anywhere in a complete word inside
the grid disqualifies that word outright at any HetHet density** (`17-…` §3) — and the
same marker, one position further left in the partial word, is invisible. The fringe scan
tests the *mismatch* predicate only.

That asymmetry is why the fitted `fringe="mis0"` variant — stop at a mismatch *or* an
opposite homozygote, the reading a symmetry argument would pick — is wrong, and it is why
this clause had to be measured rather than reasoned about. It also matches the breaking set
of the word predicate: the kinds that break a fringe are exactly the kinds counted by the
`ibs1` mask, `ibs1b` included.

---

## 3. The tail fringe: the mirror

| mismatch at | reference | `off=1` | `off=0` | `extend` | `none` |
| --- | ---: | ---: | ---: | ---: | ---: |
| tail[ 0] | **511.040** | 511 | 512 | 543 | 511 |
| tail[ 2] | **512.972** | 513 | 514 | 543 | 511 |
| tail[ 4] | **515.024** | 515 | 516 | 543 | 511 |
| tail[ 6] | **516.956** | 517 | 518 | 543 | 511 |
| tail[ 8] | **519.008** | 519 | 520 | 543 | 511 |
| tail[10] | **520.940** | 521 | 522 | 543 | 511 |
| tail[12] | **522.992** | 523 | 524 | 543 | 511 |
| tail[14] | **525.045** | 525 | 526 | 543 | 511 |
| tail[16] | **526.976** | 527 | 528 | 543 | 511 |
| tail[18] | **529.029** | 529 | 530 | 543 | 511 |
| tail[20] | **530.960** | 531 | 532 | 543 | 511 |
| tail[22] | **533.013** | 533 | 534 | 543 | 511 |
| tail[24] | **534.944** | 535 | 536 | 543 | 511 |
| tail[26] | **536.997** | 537 | 538 | 543 | 511 |
| tail[28] | **539.049** | 539 | 540 | 543 | 511 |
| tail[30] | **540.981** | 541 | 542 | 543 | 511 |

and the kinds behave the same way, with **first** replacing last:

| tail fringe | reference | reading |
| --- | ---: | --- |
| mismatches at 8 and 24 | **519.008** | stops before 8 — the **first** one wins |
| IBS0 at 8 only | **543.033** | `extend` — does not break |
| IBS0 at 24, mismatch at 8 | **519.008** | stops before 8 |
| het-vs-A1A1 at 8 | **519.008** | breaks |

    right = (first breaking marker in [64*(w1+1), hi]) - 1,   else hi

---

## 4. Both ends at once, and the snap

A call spanning the whole segment reaches into both fringes independently:

```
  clean head AND tail, full span 1071.000   both=1071  head only=1047  tail only=1047  neither=1023
  mismatch head[10] and tail[9]  1045.026   both stop=1045  extend both=1071  neither=1023
```

The open question `17-…` §5 could not phrase was what happens when the *reach* — the rule
that carries a call up to 63 markers past the mismatch in the word that bounds it — lands
on the grid's own edge. Make block word 0 unusable but IBS0-free (two mismatches), so the
run is words 1..7 and the reach carries its left end back over word 0 to exactly `64·w0`:

```
  run from word 1, clean head     535.056   snap to fringe=535  stop at grid=511  no reach=447
  run from word 1, mismatch at head[10]  524.022   snap to fringe=524  stop at grid=511  no reach=447
```

Both readings are `snap`, 24 markers away from the alternative. **Once an end reaches the
grid's own edge the word scan is over and the marker scan takes over** — and it takes over
whether that moves the end further out (a clean fringe) or pulls it back in (a fringe that
breaks). The clause is therefore not "a call *touching* `w0` runs to `seg.lo`", which is
how `17-…` §5 stated it and which tests the *run's* first word; it is a property of the
computed end. That distinction is worth 24 markers here and is invisible to a word-aligned
rig.

---

## 5. The IBD1 fringe is the mirror image, and was already right

The same rig with the canvas painted in `18-…`'s letter `K` — 34 het-vs-hom mismatches per
word, which the IBD2 pass refuses outright — holds `IBD2Seg` at 0.0000 and lets `IBD1Seg`
be read. `IBD1Seg` also carries the carrier chromosome, so it is read differentially
against a canvas with the same chromosome 1 and no chr2 IBD1 at all.

Head (an IBD1 run reaches into the word that ended it out to that word's last IBS0, and
here that word is a full wall, so the word-aligned top is 607):

| head fringe | reference | reading |
| --- | ---: | --- |
| clean | **606.944** | `extend` — runs to `seg.lo` |
| IBS0 at head[ 0] | **605.990** | stops one past marker 0 |
| IBS0 at head[ 8] | **597.997** | stops one past marker 8 |
| IBS0 at head[20] | **585.949** | stops one past marker 20 |
| IBS0 at head[31] | **574.974** | stops one past marker 31 |
| IBS0 at 5 **and** 20 | **585.949** | the **last** one wins |
| **mismatch at head[20]** | **606.944** | `extend` — a mismatch **does not break** an IBD1 call |

Tail:

| tail fringe | reference | reading |
| --- | ---: | --- |
| clean | **543.033** | `extend` — runs to `seg.hi` |
| IBS0 at tail[ 0] | **511.040** | stops one before marker 0 |
| IBS0 at tail[12] | **522.992** | stops one before marker 12 |
| IBS0 at tail[24] | **535.065** | stops one before marker 24 |
| IBS0 at tail[31] | **542.067** | stops one before marker 31 |
| IBS0 at 12 **and** 24 | **522.992** | the **first** one wins |
| **mismatch at tail[12]** | **543.033** | `extend` — does not break |

This is exactly what `Scan::left_end` and `Scan::right_end` already did, so nothing changed
in the IBD1 pass — but it was a corpus fit before this section and is a bisection now. Put
beside §1–§3 it states the symmetry cleanly:

> Past the segment's word grid the scan is marker by marker, and **each pass stops at its
> own breaking marker** — an opposite homozygote for IBD1, a het-vs-hom mismatch for IBD2.
> Neither pass can see the other's.

---

## 6. Out of sample

`seg19.predict` is the rule written over word *descriptions*, the exact mirror of
`Scan::ibd2`. Graded against the reference on canvases that had no part in choosing
anything:

```
  random canvases with random head/tail fringes (unused seeds):
  seed 4919    60 canvases:    60 agree,   0 differ
  seed 6271    60 canvases:    60 agree,   0 differ

  exhaustive: every 3-word block over {clean, quiet, 1-mis, 2-mis}
  crossed with 6 fringe shapes on each side
  384 compositions x fringes: 384 agree, 0 differ

  TOTAL 504 canvases: 504 agree, 0 differ
```

The random canvases draw fringe widths uniformly in `[0, 40]` on each side and fill them
from `{mismatch, IBS0, HetHet, A2A2}`, so IBS0-in-a-fringe and mismatch-in-a-fringe both
occur throughout; the exhaustive sweep crosses all 64 three-word compositions with six
fringe shapes including a leading IBS0.

## 7. The port

`target/release/open-king` put through the same rig — the two 16-point sweeps of §1 and §3 plus
§6's exhaustive sweep, 416 filesets — reproduces the reference on **416 of 416**
(`fringecanvas.py` §7, which caches our binary's answers separately so it can never
contaminate the reference's).

---

## 8. The corpus, before and after

`tests/parity/fit/seg19.py`, all ten datasets, 982 graded rows:

| floor | rule | exact rows | `IBD1Seg` | `IBD2Seg` | mean `PropIBD` err | worst |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| 3 Mb | `17-…`/`18-…` | 747 | 982 | 896 | 0.000067 | 0.0042 |
| 3 Mb | **§1–§4** | **806** | **982** | **982** | **0.000023** | **0.0001** |
| 5 Mb | `17-…`/`18-…` | 729 | 909 | 880 | 0.000177 | 0.0598 |
| 5 Mb | **§1–§4** | **755** | **910** | **946** | **0.000161** | 0.0641 |
| 10 Mb | `17-…`/`18-…` | 692 | 841 | 877 | 0.000399 | 0.0874 |
| 10 Mb | **§1–§4** | **713** | **844** | **937** | **0.000389** | 0.0917 |

The "before" rows are `tests/parity/fit/seg18.py`, which scores the committed rule.
`seg19.py` also prints a `fringe=extend` baseline of its own; that one reads 729/982/862 at
3 Mb and is **not** the same thing — turning this file's fringe knob off also moves where
the retired extend clause sits relative to the push (§4), so it is a third rule and not the
one that shipped. Compare against `seg18.py`.

Exact rows and mean error improve at all three floors, which is the bar this project
requires before a rule lands. The reported pair set is unchanged and still exactly right
(0 extra, 0 missing) everywhere. The worst single row rises slightly at 5 and 10 Mb; §10
says what that is.

Whole-invocation parity went **408 → 436 of 480**. `core/*__related` fell from 25 failures
to 1 and `ibdseg/*` from 45 to 42. (One further case, `apps/*__cluster`, was fixed by
concurrent work in `analysis/cluster.rs` and is not this campaign's.)

`tests/parity/fit/port19.py` checks the shipped engine against `seg19.py` row by row on the
corpus, so the scorecard above is a statement about the binary and not only about a Python
model.

---

## 9. What is left, and where it is

> **SUPERSEDED by `20-seg-writer.md`, and the error is instructive.** Everything from here
> to the end of this document rests on one assumption stated in the first line of §9 —
> that `.seg`'s `PropIBD` is a *third rounding of the same two base-pair totals*, so a
> disagreement there must mean our totals are slightly wrong. It is not. `.seg` computes
> `PropIBD` from the **printed** columns (`i2*1e-4 + i1*5e-5`), a fact checkable on the
> captured output alone with no caller in the way, and once each of the reference's two
> writers is given its own rule the 176 rows below go to **zero** — at 3 Mb and at both
> raised floors — with no change to any estimate. The measurements in §9 and §10 are real;
> the conclusion drawn from them ("the IBD1 pass is systematically about a marker short on
> a tenth of the corpus") is false. §10's two negatives still stand on their own evidence
> and were re-confirmed later. Kept unedited as a record of how a well-built instrument
> aimed at the wrong question produces confident, consistent, wrong answers for a whole
> campaign.


At the default floor the residual has **left both estimate columns**:

```
=== --seglength 3 Mb: 982 rows, 176 not byte-identical
  wrong IBD1Seg 0   wrong IBD2Seg 0   wrong PropIBD 176   wrong InfType 0
  PropIBD sign: over 33  under 143
  |PropIBD delta| in ulps: 0.5:26  0.6:57  0.7:66  0.8:8  0.9:7  1.0:7  1.1:4  1.2:1
  by InfType: 2nd 35/144  3rd 24/80  4th 7/16  Dup/MZ 0/2  FS 60/152  PO 0/315  UN 50/273
```

`PropIBD` is `IBD2Seg + IBD1Seg/2`, a **third** rounding of the same two base-pair totals,
and it therefore reads about one more bit of them than either column does alone. Every one
of the 176 rows has both columns exact; what they disagree about is a sub-ulp difference in
the totals, and the sign is biased — we are short on 143 of 176.

`tests/parity/fit/prop19.py` turns that into an instrument. The three printed columns bound
the reference's `(ibd1, ibd2)` to the intersection of a square and a diagonal band, so each
row yields a *directed, bounded* correction rather than a bare disagreement:

```
  rows with IBD2 == 0 on both sides: 115 of 823 wrong   <- PropIBD is IBD1Seg/2 alone
  rows carrying IBD2:                 61 of 159 wrong
  IBD2 == 0 rows   short 101  long  14   |correction| bp: min 608  median 22668  max 61592
  IBD2 > 0 rows    short  42  long  19   |correction| bp: min  13  median 38719  max 173138
```

**115 of the 176 carry no IBD2 at all on either side.** On those rows `PropIBD` is
`IBD1Seg/2` and nothing else, so the missing base pairs are unambiguously the **IBD1**
pass's — and 101 of the 115 need us to report *more*, by a median of 22 668 bp on the
`PropIBD` axis, which is 45 kb of `IBD1Seg`, about **one marker gap**.

That is the sharpest statement the corpus can make: `18-…`'s "`IBD1Seg` is exact on all 982
rows, the IBD1 pass is finished" is true **to four decimal places only**. One digit further
down, the IBD1 pass is systematically about a marker short on a tenth of the corpus.

## 10. What is ruled out

**Not a merge.** Two IBD1 calls separated by a single bad word come out adjacent — the
first ends on that word's last opposite homozygote, the second starts one marker later — so
reporting them as one segment instead of two adds exactly `pos[lo2] - pos[hi1]`, one marker
interval. The corpus has 287 such touching pairs, and one marker gap is exactly the deficit
§9 measures. It is still wrong, twice over: `18-…` §5 already measured the reference at
510.0 against a merged 511 at four different IBS0 placements, and the corpus screen
(`tests/parity/fit/merge19.py`) prices it out —

```
  split    exact  806  ibd1  982  ibd2  982  prop  806  MAE 0.000023
  touch    exact  768  ibd1  917  ibd2  982  prop  815  MAE 0.000023
```

nine more `PropIBD` rows bought with sixty-five `IBD1Seg` rows. That is the shape of a
fitted fiction and it is refused.

**Not a single one-marker endpoint rule.** `tests/parity/fit/where19.py` enumerates every
single-marker perturbation of our own call set on the 115 rows and asks which would land
`PropIBD` inside the allowed band:

```
  fixable by moving ONE endpoint by ONE marker: 65 of 115
  rows by the set of perturbations that work:
    extend hi,extend lo   56        none   50        shrink hi,shrink lo    8
```

**50 of 115 rows cannot be fixed by any one-marker move at all**, so they need two or more;
and on the 56 that can, extending *either* end works, so the probe does not discriminate
between them and cannot name a clause. Whatever is left is either diffuse across several
calls per row or lives somewhere other than an IBD1 endpoint.

**Not the IBD1 fringe.** §5 measured all four of its cases against the reference and the
committed rule reproduces every one.

**Not a rounding convention.** The hypothesis that KING derives `PropIBD` from the two
already-rounded columns fits all six of `nuclear`'s sibling rows perfectly — three that
pass and three that fail — and then dies on the corpus:
reconstructing `PropIBD` from the reference's own printed `IBD1Seg` and `IBD2Seg` in exact
rational arithmetic reproduces 3 958 of 4 172 golden rows under round-half-up and 3 494
under round-half-even, and the failures fall on **both** sides of the tie. `PropIBD` is
computed from the raw totals, and the raw totals are what differ.

## 11. The other floor

At 5 and 10 Mb `IBD2Seg` is 946 and 937 of 982, not 982, and those failures are *large*
(worst row 0.0641 and 0.0917 against 0.0001 at 3 Mb). A call whose length is right to
within a fraction of an ulp can still fall on the wrong side of a 5 Mb threshold, and then
the whole call is kept or dropped — so a residual that is invisible at 3 Mb becomes a
multi-megabyte error at 5. This is consistent with §9's finding and is probably the same
defect seen through a threshold rather than a separate one. `18-…` §9's deliberately
omitted `--seglength`-triggered IBD1 run merge lives at these floors too and is still out.

## 12. Where a continuation should start

1. The IBD1 pass is short by about a marker on ~10 % of rows and the corpus cannot say
   where. **`fringecanvas.py` can now build a segment with an arbitrary marker offset**, so
   the natural next instrument is an IBD1 canvas that varies the *interior* endpoint rule
   under a ruler fine enough to see one marker — `18-…` §2's endpoints were measured at a
   coarser one.
2. `tests/parity/fit/seglocal.py` (from this campaign's predecessor) localises a pair's
   IBD2 total to a single usable segment by neutralising the rest of the fileset. The same
   trick applied to `IBD1Seg` would turn §9's per-row deficit into a per-segment one, which
   is what `where19.py` needed and did not have.
3. The 100 Mb usable-total floor (`17-…` §2) is still not modelled. No corpus dataset is
   near it.
