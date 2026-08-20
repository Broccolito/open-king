# The gate window has a length, and the floor asks about it

**Status: measured, validated out of sample, committed to Rust.** This closes the last two
`.seg` parity cases (**475 → 477 of 480**) and takes `--seglength 10` from 970/972 to
**982 of 982 on both estimate columns**. All three captured floors — 3, 5 and 10 Mb — are
now byte-exact on every one of the 982 corpus rows, with MAE 0.000000 and worst row 0.0000.

`docs/research/21-push-merge.md` §8.1 handed on a sharp profile and a named suspect: at
10 Mb every wrong `IBD2Seg` was too high and every wrong `IBD1Seg` too low, with
`d1 = -d2` on ten of the twelve rows, and "the gap is the only floor-dependent term left …
a canvas whose gap sits between 5 and 10 Mb, swept at both floors, would say whether
`gap < seglength` acquires a second, absolute bound."

**The gap is not it, and neither is the merge.** The profile was right about the arithmetic
and wrong about the cause, and it was wrong because it was inferred from two genome-wide
numbers. The first thing this campaign built was an instrument that does not have to infer.

No KING source was read. Every rule below is a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2), run
either on the corpus's own filesets or on filesets built for the purpose.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| **The reference's answer can be read one chromosome at a time**, on the corpus's own data, by *muting* every other chromosome for the probe pair. `IBD?Seg × denom` is then that chromosome's own called length in bp, to ±35 kb. Every campaign before this one had to guess which segment of which pair a wrong row came from; both guesses `21-…` §8.1 recorded were wrong. | §1 |
| **The floor is asked twice, and the second question is not about the reported call.** A run is emitted only if its **gate window** — `gs` through `ge_of(e)`, the very span the informativeness gate counts over — spans at least `seglength / 2`. On the corpus an 11.2066 Mb IBD2 call is reported at `--seglength 6.290751` and gone at 6.290752, which is twice its one-word 3.145375 Mb window plus one. | §2 |
| Bisected to the base pair on **two independent corpus calls** in two datasets, and then on purpose-built canvases at **four spacings**. The division is integer, the same one the push uses (`21-…` §2.5). | §2, §3 |
| **Both ends of the window are the gate's own**, each separated from its alternative by a canvas: the left end is the gate-start word `gs`, not the run's first word; the right end is `ge_of(e)`, the word the call reaches into, not the run's last word. | §3 |
| **The IBD1 pass has the same bound**, over the run's own complete words, with the comparison one unit of `seglength / 2` **tighter**: IBD2 keeps a window of exactly `min_bp / 2` and IBD1 does not. Bisected at four spacings; the corpus cannot see this clause at all and 32 of 360 held-out canvases can. | §4, §7 |
| **A gate-refused run's own markers do enter the merge budget.** `20-…` §11 item 4 left it undecided. They do: the IBD1 budget is summed over *every* word between the two runs while [`MERGE_MAX_WORDS`] still counts only the unusable ones, so a refused run is stepped over by the cap and paid for by the budget. | §5 |
| **The bound is asked at emit, not with the gate.** A run it refuses still merges with its neighbour, and the merged window is then measured whole. Asking it with the gate scores 980/980 at 10 Mb against 982/982. | §6 |
| Corpus: 3 and 5 Mb **do not move**. At 10 Mb `IBD1Seg` 970 → **982**, `IBD2Seg` 972 → **982**, MAE 0.000067 → 0.000022, worst row 0.0081 → **0.0001**. Measured from the shipped binary against the goldens, all three floors are 982/982/982/982 with MAE 0.000000. | §7 |
| Held out: IBD2 **353 of 360** against 344 for `21-…` on three unused seeds, winning all nine canvases where the two rules differ; IBD1 **360 of 360** against 328. | §7 |
| Harness: **475 → 477 of 480**. Closed: `ibdseg/{bigish,multifam}__ibdseg_seglength10`. Self-check 480/480; 320 tests, clippy and fmt clean; `check_mirror.py` green on all 2 946 corpus rows and 861 `MaxIBD2` values. `--ibs` is untouched and all 13 of its cases stay byte-identical. | §7 |

---

## 1. The instrument: read one chromosome at a time

`--ibdseg` prints two genome-wide proportions. On `multifam` one wrong 11 Mb call in 690 Mb
of usable genome moves `IBD2Seg` by 0.0162 and says nothing at all about *where*. Every
previous campaign closed that gap by pattern-matching the residual against candidate
mechanisms, and `21-…` §8's two candidates — an invented merge, then the merge's gap — were
both wrong.

`docs/research/fixtures/chrprobe.py` measures instead. For a probe pair `(i, j)` it rewrites
the fileset with every marker outside the chromosomes of interest **muted for that pair
only**: sample `i` set to `A1A1` and sample `j` to `A2A2`, so every muted marker is an
opposite homozygote. A chromosome muted that way carries no usable word on either pass and
contributes exactly zero, and the difference between two runs that differ by one chromosome
is that chromosome's own called length.

Three details make it exact rather than suggestive, and each was a wrong answer first.

* **Mute, do not subset.** Deleting the other chromosomes' `.bim` rows changes the answer:
  KING packs the whole retained marker list into 64-marker words, so removing markers
  re-phases every later usable segment and moves its calls. The first version of this rig
  did that and produced a table in which *fifteen* chromosomes disagreed, including one
  whose IBD2 went *negative*.
* **Recompute the denominator.** The printed total is `%.1lf Mb`. `Pallsegs.txt` names each
  usable segment's first and last SNP, so the exact denominator comes back out of the
  `.bim`; what is left is the 4-dp rounding of `IBD?Seg`, about ±35 kb on `multifam`.
* **`--seglength` is in Mb and KING clamps it to 1..10.** Outside that range it falls back
  to its 3 Mb default with no diagnostic — `Minimum segment length is set as 3000000 bp`
  and nothing else. Passing base pairs by mistake reads as "the same answer at every
  floor", which is exactly what a floor-independent rule looks like. `chrprobe.run()` takes
  base pairs and asserts the range. (This is also the whole of `21-…` §8 item 2: the
  reference "stops behaving like a floor above ~10 Mb and below 1 Mb" because outside
  1..10 Mb the flag is discarded.)

Two controls: muting nothing reproduces the golden `.seg` row exactly, and at
`--seglength 5`, where the corpus is already exact, every chromosome of every probe pair
agrees.

With it, the twelve wrong rows at 10 Mb resolve immediately. `multifam 2/4`:

```
--- floor  5.0 Mb (carrier chr1) ---   (every chromosome agrees)
--- floor 10.0 Mb (carrier chr1) ---   chr   ref_ibd1   our_ibd1 |   ref_ibd2   our_ibd2
    13    29.1833    17.9505 <<<|     0.0000    11.2066 <<<
```

One chromosome, one call, one direction: **we report an 11.2066 Mb IBD2 call at
`--seglength 10` and the reference does not.** No merge is involved — this caller performs
no IBD2 merge anywhere on that pair, at any floor.

## 2. The floor is asked twice

Sweeping `--seglength` against that one chromosome turns the question into a bisection.
`chrprobe.flip('multifam', 2, 4, 13, 1, 1_000_000, 10_000_000)`:

```
  kept at 6 290 751 bp        dropped at 6 290 752 bp
```

A call that measures 11.2066 Mb is dropped at a floor of 6.29 Mb. So the reported length is
not what the floor is compared against — or rather, not the *only* thing.

The call is one word wide: scan word 7 of `multifam`'s chr13 segment, markers 11 328–11 391,
whose span is **3 145 375 bp**. Its reported length comes from the endpoint reach, which
runs 119 markers back on the left and 42 forward on the right. And

    2 x 3 145 375 + 1 = 6 290 751

is exactly the last floor at which it survives. With integer division the rule reads

    span(gate window) >= seglength / 2

— `6 290 751 / 2 = 3 145 375` still passes and `6 290 752 / 2 = 3 145 376` does not.

**A second, independent corpus call agrees to the base pair.** `bigish 66/69` chr6, a call
of 10.1511 Mb over a one-word window of 3 148 672 bp:

```
  kept at 6 297 345 bp        dropped at 6 297 346 bp        2 x 3 148 672 + 1 = 6 297 345
```

Different dataset, different marker grid, different window; same law, no free parameter.

### 2.1 On a canvas built for it

`docs/research/fixtures/window1.py` §2 reproduces the geometry away from the corpus.
The block is `[M0, M0, CLEAN]` at 50 kb spacing, where

```
CLEAN   64 HetHet                       usable, inf2 = 64, mismatch-free
M0      2 het-vs-hom mismatches at bits 0..1, 62 HetHet   unusable; ml = 1
```

`M0`'s last mismatch sits at bit 1, so the left end of the call — `pos[64(u-1) + ml - 63]` —
lands inside the word *before* it, and the reported call is 189 marker intervals (9.450 Mb)
over a run, and therefore a window, of 63 (3.150 Mb):

```
  --seglength sweep      1.000:189   6.400:0
  bisected: kept at 6 300 001 bp, dropped at 6 300 002 bp        2 x 63 x 50 000 + 1 = 6 300 001
```

Bisected again at 30 kb, 45 kb and 70 kb spacing (§4's table): the IBD2 flip is
`(2w + 1, 2w + 2)` every time.

## 3. Both ends of the window are the gate's own

The window is `[gs .. ge_of(e)]` — the same one `13-informativeness-gate.md` counts `inf2`
over, and the same one `21-…` §3 measures the merge's interruption against. Each end is
separated from its alternative by one canvas.

**The right end is `ge_of(e)`, not `e`.** Put a mismatch-only word after the run:

```
  [M0, M0, CLEAN, MISW]     window = 2 words = 6.350 Mb, call = 253 markers = 12.65 Mb
  --seglength sweep         1.000:253                    (constant, 1 .. 10 Mb)
```

Nothing under 10 Mb can reach `2 × 6.350`, and nothing kills it. Had the window ended at
the run's last word it would have died at 6.300 Mb, exactly as §2's canvas does.

**The left end is `gs`, not the run's first word.** Lead the run with a word that is usable
but not mismatch-free, so its gate-start word is the second:

```
  [M0, M0, U1, CLEAN]       run = 2 words, window = 1 word, call = 253 markers
  --seglength sweep         1.000:253   6.400:0
  bisected: kept at 6 300 001 bp, dropped at 6 300 002 bp
```

The run is two words wide and the bound still cuts at one word's span. Had it measured the
run, this canvas would have behaved like the previous one.

## 4. The IBD1 pass has it too, one unit tighter

The IBD1 window is the run's own complete words — the span `Scan::informative` counts
`inf1` over; there is no `ge_of` on that pass. `window1.py` §4 builds the mirror fixture:
a one-word IBD1 run flanked by words whose single opposite homozygote sits at bit 0 on the
left and bit 63 on the right, so `left_end`/`right_end` reach almost a whole word each way
and the reported call is again 189 markers over a 63-marker run.

```
  spacing  window     IBD2 flip                    IBD1 flip
   30 000  1 890 000  (3 780 001, 3 780 002)       (3 779 999, 3 780 000)
   45 000  2 835 000  (5 670 001, 5 670 002)       (5 669 999, 5 670 000)
   50 000  3 150 000  (6 300 001, 6 300 002)       (6 299 999, 6 300 000)
   70 000  4 410 000  (8 820 001, 8 820 002)       (8 819 999, 8 820 000)
```

Both passes have the bound; the IBD1 flip is two base pairs of `--seglength` — one unit of
`seglength / 2` — earlier at every spacing. The IBD2 pass keeps a window of exactly
`min_bp / 2` and the IBD1 pass does not:

    IBD1:  span >  seglength / 2            IBD2:  span >= seglength / 2

A two-word IBD1 run (window 127 markers = 6.350 Mb) survives every floor under 10 Mb, as it
must.

**One reading this rig cannot separate.** "IBD1 compares strictly" and "the IBD1 span is
measured one base pair shorter" make identical predictions at every spacing and for every
window, including windows of odd length: `Q > L/2` and `Q - 1 >= L/2` are the same predicate
over the integers. The committed form is the strict comparison.

## 5. A gate-refused run is stepped over by the cap and paid for by the budget

`20-…` §11 item 4: "Whether a gate-refused run's own markers also enter `X` … A canvas with
a gate-refused run carrying close to 9 A1A1/A1A1 markers between two interruptions whose
budget is exactly on the boundary would settle it." The corpus produced exactly that canvas
by itself. `multifam 11/18`, chr2, the second of the two faults §1's table localises:

```
  runs        [10..12] inf1 26 GATE-OK        [14..14] inf1 4 refused        [16..18] inf1 21 GATE-OK
  gap 9 652 629 bp        unusable words 13 and 15        word 14 is the refused run
  word 13:  z 2  U 0  V 4          word 14:  z 0  U 0  V 4          word 15:  z 2  U 3  V 4
```

`chrprobe.flip` bisects the merge at **`min_bp` = 9 652 630**, and the run-to-run gap is
9 652 629 — `gap < seglength`, exactly `20-…` §2's rule, no new gap clause anywhere. What
differs is the budget. Over the unusable words alone, `bad = 4`, `V = 8 < MIN_INFORMATIVE`
so `X = U = 3`, and `4 × (4 − 2) = 8 ≤ 3` fails. Over every word between the runs, word 14's
own four het-vs-A1A1 markers take `V` to 12, the switch of `20-…` §5 hands `X` those, and
`8 ≤ 12` passes.

`window1.py` §5 turns that into a sweep. Four-word runs either side, interruption
`[INTW, refused(vr), INTW]` where the two `INTW` carry two opposite homozygotes and four
het-vs-A1A1 each and `refused(vr)` is a run of one word with `vr` het-vs-A1A1 markers and
nothing else, at `--seglength 5`:

```
   vr   markers  (words,calls)   span=unusable  span=all
    0       638  (10, 2)        638            638
    1       638  (10, 2)        638            638
    2       767  (12, 1)        638            767
    ...
    9       767  (12, 1)        638            767
   10       767  (12, 1)        767            767
   11       767  (12, 1)        767            767
```

One call is a merge, two are a split. The flip is at `vr = 2`, which is precisely where the
refused run's markers take `V` from 8 to 10; `span=unusable` is wrong on eight consecutive
rows and `span=all` on none. (At `vr ≥ 10` the middle run passes the gate itself and is an
endpoint rather than an interruption, so both readings merge for a different reason.)

The **cap** is unchanged and still counts unusable words only — that is `20-…` §6, which a
canvas of its own forced and which this campaign does not touch. So the two halves of the
merge test read different word sets, deliberately: the cap steps over a refused run, the
budget pays for it.

## 6. The bound is asked at emit, not with the gate

The gate is asked before the merge (`20-…` §6). The length bound is not.
`[CLEAN, Z1, CLEAN × 4]`: the first run's window is one word (3.150 Mb), the second's is
four (12.750 Mb), and the one interrupting word is cheap enough to merge across.

```
  L = 1.0 Mb -> 254 markers      (gap 3.25 Mb >= floor: no merge, two calls)
  L = 5.0 .. 9.0 Mb -> 383       (one merged call over six words)
```

At 7, 8 and 9 Mb the bound would refuse the first run on its own, and the answer is still
the merged 383. Asked with the gate, the refused run could not be an endpoint and the answer
would be the second run alone, 255. On the corpus the same knob is worth two rows on each
column at 10 Mb: emit-only 982/982, pre-merge 982/980, and pre-merge's worst row is 0.0536
against 0.0001.

## 7. The corpus, and out of sample

`python3 tests/parity/fit/seg23.py grid`:

```
--seglength 3 Mb
  21 (committed)        exact 806  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
  23 (window bound)     exact 806  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
--seglength 5 Mb
  21 (committed)        exact 817  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
  23 (window bound)     exact 817  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
--seglength 10 Mb
  21 (committed)        exact 811  ibd1 970  ibd2 972  extra 0  miss 0  MAE 0.000067  worst 0.0081
  23 (window bound)     exact 820  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000022  worst 0.0001
  23 without window     exact 813  ibd1 972  ibd2 972  extra 0  miss 0  MAE 0.000060  worst 0.0081
  23 without window1    exact 820  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000022  worst 0.0001
  23 without span_all   exact 818  ibd1 980  ibd2 982  extra 0  miss 0  MAE 0.000028  worst 0.0039
  23 window pre-merge   exact 818  ibd1 982  ibd2 980  extra 0  miss 0  MAE 0.000079  worst 0.0536
```

The IBD2 bound is worth ten rows on each column, the budget two more, and **`window1` — the
IBD1 bound — is worth nothing here**: no corpus run is in the regime where it can fire. It
is committed anyway, on the same footing as `21-…`'s missing word cap: it is bisected
directly against the reference at four spacings (§4) and it is worth 32 of 360 canvases out
of sample. This is the second clause the corpus plainly cannot see and the reference plainly
has.

The same measured from the shipped binary against the goldens, with `.seg`'s own `PropIBD`
rule (`python3 tests/parity/fit/scorecard.py`):

```
  floor   rows  exact   ibd1   ibd2  extra  missing        MAE    worst
    3 Mb    982    982    982    982      0        0   0.000000   0.0000
    5 Mb    982    982    982    982      0        0   0.000000   0.0000
   10 Mb    982    982    982    982      0        0   0.000000   0.0000
```

**Every row at every captured floor is byte-exact.** Before this campaign 10 Mb stood at
970 exact rows.

`python3 tests/parity/run_parity.py --impl ./target/release/open-king`: **475 → 477 of 480.**
Closed: `ibdseg/bigish__ibdseg_seglength10` and `ibdseg/multifam__ibdseg_seglength10` — the
last two `.seg` cases in the suite. Still open: `apps/bigish__build`,
`core/bigish__related_degree2` and `ibdseg/bigish__related_degree2_ibdseg`, none of them a
`.seg` difference and all three unchanged by this work (verified by rebuilding with the
clause reverted). Self-check 480/480, 320 tests, clippy and fmt clean, and `check_mirror.py`
green on all 2 946 corpus rows at the three floors plus the 861 `MaxIBD2` values.

`--ibs` is untouched — nothing here is in `Scan::ibd2_words` — and all 13 of its parity
cases stay byte-identical.

### Held out

`window1.py` §7 and §8. Three unused seeds, 60 random canvases each, at the two spacings
where the bound *can* decide a call. That regime is worth stating, because it is why the
corpus sees this clause only at 10 Mb: a one-word window spans `63s` and an endpoint can
reach at most two more words, so the bound decides a call only when `126s < seglength ≤ 189s`.
At the corpus's 50 kb spacing that interval is 6.3–9.5 Mb, which straddles the 10 Mb floor
and nothing else. Here it is reproduced twice, at 30 kb against `--seglength 5` and at 60 kb
against `--seglength 10`.

```
  IBD2      21: 344/360      23: 353/360      (the two rules differ on 9; 23 is right on all 9)
  IBD1      20: 328/360      +span: 328/360   +window: 360/360
```

The IBD1 result is the sharp one: a clause the corpus scores at exactly zero fixes 32 of 360
held-out canvases and leaves none wrong. The 15 canvases both IBD2 rules miss are a
pre-existing residual of the shared geometry in this harsher draw — the two rules agree on
every one of them — and are recorded in §8 below.

Nothing in either battery chose a constant. The samplers were tuned (that is a property of
the *draw*, not of the rule): `rblock4` draws canvases as run/interruption units rather than
word by word, and `_mk` places a word's mismatches in a contiguous block two times in five,
because scattering `k` bits over 64 puts the last one near bit 63 almost always and a
uniform draw therefore almost never builds the long reach the bound decides.

## 8. What is still open

1. **Both `.seg` estimate columns are saturated at all three captured floors**, so the
   corpus can no longer grade this caller at all. Anything further has to be graded on
   canvases, and the two batteries above are where the remaining residual is visible:
   15 of 360 IBD2 canvases in the §7 draw are wrong under both `21-…` and `23-…`. That draw
   is harsher than `mergelab`'s — many adjacent unusable words, mismatch bits in contiguous
   blocks, and endpoints reaching two words out — which is exactly the region `17-…` §14's
   bridge and `19-…`'s fringe were measured away from. It is the next thing to look at, and
   `chrprobe.py` cannot help with it because the corpus has no wrong row left to localise.
2. **Why a half.** `seglength / 2` now appears twice: arming the push (`21-…` §2) and
   bounding the gate window (§2 here), both with the reference's own integer division. Two
   independent bisections landing on the same constant is suggestive of one shared quantity
   in the reference and is no more than that.
3. **Why the two passes differ in strictness** (§4), and whether the difference is a
   comparison or a one-base-pair difference in the measured span. Not separable by any
   fixture — the two predicates are equal over the integers.
4. **The 100 Mb usable-total floor** and the 1..10 Mb `--seglength` clamp (§1) are both
   bisected and neither is modelled by this crate. No corpus dataset approaches either.

## 9. Reproducing everything here

```bash
cd docs/research/fixtures
python3 chrprobe.py multifam 2 4 5 10          # §1 — the fault, localised to chr13
python3 -c "import chrprobe as C; print(C.flip('multifam', 2, 4, 13, 1, 1000000, 10000000))"
python3 -c "import chrprobe as C; print(C.flip('bigish', 66, 69, 6, 2, 1000000, 10000000))"
python3 -c "import chrprobe as C; print(C.flip('multifam', 11, 18, 2, 1, 1000000, 10000000, what=0))"
python3 window1.py                             # §2 §3 §4 §5 §6, then both batteries
python3 window1.py 2 4                         # just the two bisections

cd ../../../tests/parity/fit
python3 seg23.py                               # the corpus scorecard at 3 / 5 / 10 Mb
python3 seg23.py grid                          # each clause dropped in turn
python3 scorecard.py                           # the same from the shipped binary
python3 check_mirror.py                        # engine.py still equals that binary
cd ../../..
python3 tests/parity/run_parity.py --impl ./target/release/open-king
```

`window1.py` shares `mergelab_measured.json`; `chrprobe.py` caches nothing and drives the
reference directly, so it needs `$KING` pointing at KING 2.3.2 (its default). As with the
other rigs: do not point `$KING` at a non-reference build in this directory, because
`window1.py` writes whatever it measures into that cache.
