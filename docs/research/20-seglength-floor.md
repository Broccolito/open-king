# The `--seglength` run merge, measured on both passes

**Status: measured, validated out of sample, committed to Rust.** This closes eight of the
sixteen remaining parity cases (**464 → 472 of 480**), including the whole `sexchr` pair
and six of the eleven `--seglength 5/10` `.seg` cases, and it does not touch a single byte
at the default floor.

`docs/research/18-ibd1-caller.md` §9 measured this clause on two of its five conditions,
found that implementing it that way made the corpus much worse, and left it out with a
note. This document measures the other three, which is what turns it from a liability into
the largest single correction the `.seg` caller has had since §6 of that campaign.

No KING source was read. Every rule below is a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) run on
filesets built for the purpose, or a score against the captured parity corpus.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| The diagnosis handed to this campaign — "the geometry is right, something about *how the floor is applied* is wrong" — is **half right**. Nothing about the floor's application changed. What is floor-*dependent* is a **merge**, and the reference's own captured output proves one exists before any canvas is built: `IBD1Seg` **rises** with `--seglength` on 80 corpus rows and `IBD2Seg` on 27. An estimate that grows when the minimum segment length grows cannot come from a filter. | §1 |
| **Both passes merge**, with one geometry. Two runs are joined iff at most **two unusable words** lie between them, the run-to-run gap is **strictly** under `--seglength`, and a budget test passes. The three-word cap is absolute — a three-word interruption never merges, at any floor, however little it carries and however short the gap. | §2, §3 |
| The budget is `cost · (bad − 2) ≤ X`. **Two bad markers are free and each further one costs `cost` informative markers** — `cost = 4` on the IBD1 pass, `3` on the IBD2 one. Bisected on ten values of `bad` (IBD1) and on the HetHet sweep (IBD2); the "2" is bisected independently off the all-A2A2 interruption, which carries no informative marker at all and merges at two opposite homozygotes but not at three, **on both passes**. | §4, §7 |
| What counts as *bad* and as *informative* is each pass's own business, and the split is not the obvious one. IBD1: `bad` is opposite homozygotes only — 64 het-vs-hom mismatches in the interruption change nothing — and `X` is **A1A1/A1A1** markers, *not* `inf1`. IBD2: `bad` is opposite homozygotes **plus** het-vs-hom mismatches, and `X` is `inf2`. | §4, §7 |
| One clause no symmetry argument produces: on the IBD1 pass, if the het-vs-A1A1 markers in the interruption alone reach **10** — the informativeness gate's own constant — then `X` becomes *those* instead of the A1A1/A1A1 count. Bisected at 9 / 10 against A1A1/A1A1 loads from 16 to 40, and it is why the merge region in the `(A1A1/A1A1, het-vs-A1A1)` plane is **not convex**. | §5 |
| **The gate runs before the merge.** A run under `inf1 ≥ 10` is refused outright; it then lies *inside* a later interruption rather than ending one, so it can never be a merged segment's endpoint but does not stop a merge either. Both halves are forced: without the first, four of 60 random canvases are wrong; without the second, two others are. | §6 |
| **A merged call may not satisfy the ">10 Mb" pair-reporting filter** — the reference reports the *same pair set* at 3, 5 and 10 Mb on all ten corpus datasets, while the merge is floor-dependent by construction. Letting merged calls feed it invents 251 pairs at 5 Mb and 297 at 10. This closes `18-…` §11's third open item. | §8 |
| Corpus: at 3 Mb **nothing moves** (982 / 982 / 806, MAE 0.000023). At 5 Mb `IBD1Seg` **910 → 959**, `IBD2Seg` 946 → 947, whole rows 755 → 795, worst row 0.0641 → **0.0111**. At 10 Mb `IBD1Seg` **844 → 960**, `IBD2Seg` 937 → 945, whole rows 713 → 793, worst 0.0917 → **0.0111**. 0 extra / 0 missing everywhere. | §9 |
| Held out: **360 of 360** random canvases on three unused seeds, graded at **5 and 10 Mb specifically**, with the merge firing on most of them. | §10 |

---

## 0. The instruments

`docs/research/fixtures/mergelab.py`, which reuses `segcanvas.Canvas` and `ibd1canvas`'s
alphabet unchanged and keeps its own answer cache (`mergelab_measured.json`) so the two
earlier campaigns' caches are untouched.

Two things had to change about the rig before the merge was visible at all.

* **Spacing.** `ibd1canvas` runs at `nw2 = 16`, which forces a chr2 spacing of 88 000 bp to
  clear the 100 Mb usable-total floor, and there a one-word interruption leaves a
  5.72 Mb gap — over the legal `--seglength` range for most of it. Widening to `nw2 = 70`
  drops the spacing to 20 000 bp, so a one-word gap is 1.30 Mb and a two-word gap 2.58 Mb,
  and the merge can be swept from below 3 Mb to 10 Mb. **This is why `18-…` §9 saw only a
  sliver of the clause.**
* **A readout that is not one marker wide.** With a single interrupting word the merged and
  split answers differ by exactly one marker, because consecutive IBD1 calls come out
  adjacent (`18-…` §2.2). Choosing run lengths so that *both* split calls fall under the
  floor and the merged one clears it turns the reading into a clean binary — `0` or the
  merged length — worth hundreds of markers instead of one.

## 1. The merge exists, and the captured corpus already said so

Before any fixture: for each of the 982 corpus rows compare the reference's own
`IBD1Seg`/`IBD2Seg` at 3, 5 and 10 Mb.

```
rows 982   IBD1Seg rises with --seglength on 80 rows   IBD2Seg on 27
  dups  (6,7)   IBD2Seg  0.9223 -> 0.9877 -> 0.9877
  nuclear (2,3) IBD1Seg  0.3484 -> 0.3348 -> 0.3582
```

Raising the minimum segment length can only *drop* segments, so a column that rises is
proof that something joins them. `dups`' duplicate pair is the sharpest: 0.9223 at 3 Mb and
0.9877 at 5 — the reference recovers nearly the whole genome as one IBD2 segment once the
floor is high enough to swallow the breaks a duplicate's genotyping errors leave. That row
is the largest single error in the corpus and it is a merge.

## 2. The gap, and that the merge is a merge

`mergelab` E1. Two runs of `p` and `q` words separated by a two-word interruption
(`Z{0}`, `Z{63}`), swept over `--seglength` at 20 000 bp spacing, where the run-to-run gap
is 129 marker intervals = 2.580 Mb:

```
            L=1    L=2.5   L=2.57   L=2.59    L=2.6      L=3      L=4
 p1 q1    190.9    126.8      0.2    319.1    319.1    319.1    319.1
 p1 q3    319.1    255.0    255.0    447.2    447.2    447.2    447.2
 p3 q1    319.1    319.1    191.9    447.2    447.2    447.2    447.2
 p3 q3    447.2    447.2    447.2    574.8    574.8    574.8    574.8
```

The threshold is the gap, to the base pair, and the comparison is **strict**: `2.57` splits,
`2.59` merges. That reproduces `18-…` §9.1 at a third spacing.

**And the floor plays no part in *which* runs may merge.** The `p1 q1` row is the one that
matters: at `L = 2.59` both calls are individually under the floor (1.28 Mb and 2.54 Mb) —
at `L = 2.57` the row reads 0, so they really are dropped — and they merge anyway, into
6.38 Mb. So the merge is not a join of *kept* segments; it happens on the runs, before the
floor.

It is also genuinely a merge and not a floor-dependent relaxation of the word predicate,
which would look identical on the canvas above. E7 separates them:

```
  lead z=1     L=1    319.1   L=5    319.1        (a Z word next to the wall)
  lone z=5     L=1      0.2   L=5      0.2        (a Z word walled on both sides)
  K4 only      L=1    319.1   L=5    319.1
  K5 only      L=1    383.1   L=5    383.1
```

A word carrying an opposite homozygote is never absorbed into a run beside it (`lead`
reads 319, the four-word answer, not 383) and never forms a run of its own (`lone` reads
zero) at any floor. **The clause needs a run on both sides.** The word predicate of
`18-…` §2.1 stands exactly as measured.

## 3. At most two unusable words — and that cap is absolute

E3, at 20 000 bp spacing so every gap up to seven words stays under 10 Mb:

| interruption | run gap | L = 3 | L = 5 | L = 10 |
| --- | ---: | ---: | ---: | ---: |
| 1 word | 65 int = 1.30 Mb | **639** | **639** | **639** |
| 2 words | 129 int = 2.58 Mb | **703** | **703** | **703** |
| 3 words | 193 int = 3.86 Mb | 575 | 575 | 0 |
| 4 words | 257 int = 5.14 Mb | 575 | 575 | 0 |
| 5 words | 321 int = 6.42 Mb | 575 | 575 | 0 |
| 6 words | 385 int = 7.70 Mb | 575 | 575 | 0 |
| 7 words | 449 int = 8.98 Mb | 575 | 575 | 0 |

(merged = the bold values, split = 575 or, once both calls fall under the floor, 0.)

One and two merge at every floor above their own gap. **Three never merges** — not at
`L = 10`, where its 3.86 Mb gap is comfortably under the floor, and not with one opposite
homozygote per word (E2 sweeps the middle word from 1 to 64 opposite homozygotes at
`L = 5` and `L = 8`: all split). The gap condition of §2 is necessary and not sufficient,
and the second condition is a hard cap on the interruption's width.

The rig cannot separate "at most two words" from "at most 129 marker intervals", because a
run gap is always `64j + 1` markers. They make identical predictions on any word-aligned
scan, which is all of them; §11 records it as open.

## 4. The budget: two bad markers free, then a fixed price

Hold the interruption at one word, `z` opposite homozygotes at fixed bits, and vary what
else the word carries. E9 first, because it is the surprise:

```
  z=4 z=5 z=6      filler for the other 64-z markers
   .   .   .       A2A2/A2A2         (uninformative)
   .   .   .       HetHet
   .   .   .       het-vs-A2A2       (a het-vs-hom mismatch)
   .   .   .       missing
   M   M   M       A1A1/A1A1
   M   M   M       het-vs-A1A1
```

Sixty-four het-vs-hom mismatches in the interruption do not block a merge and sixty-four
HetHet do not enable one. Only the two `inf1` kinds move it. So on the IBD1 pass *bad* is
opposite homozygotes and nothing else.

EA and EB then sweep the plane. With the word carrying `z` opposite homozygotes and `h`
A1A1/A1A1 markers:

```
  ibs0= 3  merges from h = 4      ibs0= 8  merges from h = 24
  ibs0= 4  merges from h = 8      ibs0= 9  merges from h = 28
  ibs0= 5  merges from h = 12     ibs0=10  merges from h = 32
  ibs0= 6  merges from h = 16     ibs0=11  merges from h = 36
  ibs0= 7  merges from h = 20     ibs0=12  merges from h = 40
```

Ten bisections on one line: `h ≥ 4(z − 2)`, i.e.

    4 * (bad - 2) <= X                          bad = opposite homozygotes
                                                X   = A1A1/A1A1 markers

The "2" is free and independent: with `h = 0` the word merges at `z = 1` and `z = 2` and
splits at `z = 3`, which is the same statement read at the origin.

**Two independent confirmations, neither of which chose a constant.**

* *Two-word interruptions.* Words carrying 12 A1A1/A1A1 markers each, IBS0 counts `a` and
  `b`, at `L = 5`: the merge region is the **triangle** `a + b ≤ 8`, not a rectangle. The
  law predicts `4(a + b − 2) ≤ 24`, i.e. `a + b ≤ 8`, exactly — so `bad` and `X` are summed
  over the interruption, not tested per word.
* *Run lengths and floors are irrelevant.* E8 fixes the interrupting word and sweeps the
  two runs over `(1,1) (2,2) (3,3) (4,4) (1,4) (4,1) (1,6) (6,1) (2,5)` with the floor
  moved to match: the threshold is `z ≤ 5` in all nine. Nothing outside the interruption
  enters the test — which also kills every "rate over the merged segment" reading.

## 5. The switch at 10, and why the merge region is not convex

`X` is the A1A1/A1A1 count, but not always. Take an interruption with 6 opposite
homozygotes and vary A1A1/A1A1 (`U`) against het-vs-A1A1 (`V`) — EE:

```
   U\V   0  2  4  6  8 10 12 14 16 18 20
    0    .  .  .  .  .  .  .  .  M  M  M
   ...
   14    .  .  .  .  .  .  .  .  M  M  M
   16    M  M  M  M  M  .  .  .  M  M  M
   18    M  M  M  M  M  .  .  .  M  M  M
   20    M  M  M  M  M  .  .  .  M  M  M
```

Adding het-vs-A1A1 markers to an interruption that merges **stops** it merging, and adding
more starts it again. No linear function of marker counts does that — a perceptron over all
nine marker kinds leaves 36 of 600 random canvases misclassified — so a variable is missing,
and the shape says it is a switch rather than a weight.

EH bisects it. At 6 opposite homozygotes, against A1A1/A1A1 loads of 16, 24, 30 and 40:

```
  U\V     6   7   8   9  10  11  12
   16     M   M   M   M   .   .   .
   24     M   M   M   M   .   .   .
   30     M   M   M   M   .   .   .
   40     M   M   M   M   .   .   .
```

**9 merges, 10 does not, whatever `U` is.** The reading:

    X = V  if V >= 10  else  U

and 10 is [`MIN_INFORMATIVE`], the informativeness gate's own constant. Every earlier family
falls out of it: the pure-`V` sweep merges from `max(10, 4(z−2))` — 10 at `z = 3` and `z = 4`
where `4(z−2)` is only 4 and 8, and 12, 16 at `z = 5, 6` where it binds — which is what made
that family look inconsistent with the pure-`U` one. EG confirms the switch is symmetric in
the pair (`het-i vs A1A1-j` and `A1A1-i vs het-j` behave identically) and that HetHet,
het-vs-A2A2, missing and polymorphic-A2A2 markers are inert.

Scored on 600 random one- and two-word interruptions drawn independently (`/tmp/rndword.py`,
seeds 7001 and 7002): **600 of 600**.

## 6. The gate runs first — and a refused run is stepped over, not stopped at

This is the condition that decides whether the clause helps or hurts, and only the random
battery found it. Of 60 random 12-word canvases at `L = 5`, four were wrong; all four had a
run that fails `inf1 ≥ 10` at one end of the merge the model made. Refusing such runs
*before* merging fixed all four — and broke two others, in which a gate-failing run sits in
the *middle* of a span the reference merges across.

Both halves are needed, and together they are one sentence: **the gate is asked first; a run
it refuses is not an endpoint and not a barrier.** Concretely, seed 424242 canvas 16 —

```
  w0 IBS0-wall | w1 w2 w3 (inf1 64) | w4 IBS0 | w5 (inf1 8) | w6 IBS0 | w7 w8 (inf1 35)
```

— reads 573.8 marker intervals, which is `[w1 … w8]` as one call. `w5` passes the word
predicate but carries 8 informative markers, so it is refused; the merge then sees **two**
unusable words between the survivors (`w4` and `w6`), inside the cap, and joins them. Count
`w5` toward the cap and the model refuses; treat `w5` as a run and the model merges `w1..w3`
with it and reports a different, shorter call.

So the two-word cap counts **unusable words only**, and the budget of §4 is summed over
those same words. (Whether a stepped-over run's own markers also enter `X` is not decided
here — see §11 — but it can only make the test easier and no canvas has yet separated the
two.)

## 7. The IBD2 pass: same shape, its own two numbers

The `dups` row of §1 says the IBD2 pass merges too. G3–G5 measure it on `segcanvas`'s
IBD2-native rig, with the interruption made unusable by an opposite homozygote so the
`17-…` §14 bridge (which needs a mismatch-only word) cannot fire instead.

The gap rule is *identical*, bisected on the same two widths:

```
  1-word interruption, run gap 65 int = 1.300 Mb : 1.30 -> 127 (split)   1.31 -> 319 (merged)
  2-word interruption, run gap 129 int = 2.580 Mb: 2.58 -> 0   (split)   2.59 -> 383 (merged)
```

The budget has the same shape and different contents. G5, with the interruption holding 10
opposite homozygotes and 30 HetHet and `b` markers of one further kind added:

```
  kind       b=0  b=1  b=2  b=3  b=4  b=6  b=8 b=12
  het-vs-A2A2  M    M    M    .    .    .    .    .
  het-vs-A1A1  M    M    M    .    .    .    .    .
  missing      M    M    M    M    M    M    M    M
  A1A1/A1A1    M    M    M    M    M    M    M    M
```

Every **het-vs-hom mismatch** costs exactly what an opposite homozygote costs — three
markers of margin buys three of them — while missing calls and A1A1/A1A1 do not. So

    3 * (bad - 2) <= X       bad = opposite homozygotes + het-vs-hom mismatches
                             X   = inf2 = HetHet + A1A1/A1A1  (`share & ~ibs1`)

Three checks that did not choose a constant: with a HetHet filler the threshold is `z ≤ 17`
(`3z − 6 ≤ 64 − z` gives 17.5); with an all-A2A2 filler it is `z ≤ 2`, the same free
allowance as the IBD1 pass; and with 10 opposite homozygotes the HetHet count must reach
**24**, against `3(10 − 2) = 24` predicted. A het-vs-A1A1 filler never merges at any `z`,
which the law gets right for a reason worth stating: on this pass that marker is counted on
the *left*, as bad, and excluded from `inf2` on the right — so it is charged twice, and
`18-…`'s IBD1 switch has no analogue here.

## 8. What the ">10 Mb" pair filter reads

Merged calls make new long segments, so the clause immediately threatens the reported pair
set. The reference settles it without a fixture:

```
  nuclear threegen multifam dups missing monomorphic sexchr unrelated admixed bigish
  ... identical pair set at all three floors, all ten datasets (982 rows)
```

The merge is floor-dependent by construction; the pair set is not; therefore the filter does
not read merged calls. Feeding them to it invents **251** pairs at 5 Mb and **297** at 10 —
the same failure `18-…` §9 predicted at 1 054 for its cruder rule. `pair_segments` therefore
computes the longest segment from the unmerged call sets and the two estimate columns from
the merged ones. This closes `18-…` §11's third open item.

## 9. The corpus, before and after

`python3 tests/parity/fit/seg20.py`:

```
--seglength 3 Mb
  19 (no merge)              exact  806  ibd1  982  ibd2  982  extra 0  miss 0  MAE 0.000023  worst 0.0001
  20 (merge, both passes)    exact  806  ibd1  982  ibd2  982  extra 0  miss 0  MAE 0.000023  worst 0.0001
--seglength 5 Mb
  19 (no merge)              exact  755  ibd1  910  ibd2  946  extra 0  miss 0  MAE 0.000161  worst 0.0641
  20 IBD1 merge only         exact  794  ibd1  958  ibd2  946  extra 0  miss 0  MAE 0.000150  worst 0.0629
  20 IBD2 merge only         exact  756  ibd1  911  ibd2  947  extra 0  miss 0  MAE 0.000095  worst 0.0111
  20 (merge, both passes)    exact  795  ibd1  959  ibd2  947  extra 0  miss 0  MAE 0.000086  worst 0.0111
--seglength 10 Mb
  19 (no merge)              exact  713  ibd1  844  ibd2  937  extra 0  miss 0  MAE 0.000389  worst 0.0917
  20 IBD1 merge only         exact  787  ibd1  957  ibd2  937  extra 0  miss 0  MAE 0.000242  worst 0.0905
  20 IBD2 merge only         exact  716  ibd1  846  ibd2  945  extra 0  miss 0  MAE 0.000271  worst 0.0111
  20 (merge, both passes)    exact  793  ibd1  960  ibd2  945  extra 0  miss 0  MAE 0.000134  worst 0.0111
```

**The default floor does not move at all** — the merge cannot fire there on real spacings
often enough to matter, which is the same fact that made `IBD1Seg` exact on all 982 rows at
3 Mb in the first place. At the two floors that were wrong, `IBD1Seg` gains 49 and 116 rows
and `IBD2Seg` 1 and 8, the worst row falls by a factor of six at both, and the pair set stays
exact. The two passes' merges are close to independent: each helps on its own and the two
together are the sum.

`python3 tests/parity/run_parity.py --impl ./target/release/king`: **464 → 472 of 480.**
Closed: `ibdseg/sexchr__{ibdseg_degree2, related_degree2_ibdseg}` (both were also missing
`kingX.seg` entirely) and six of the eleven `.seg` floor cases —
`{admixed, dups, threegen}__ibdseg_seglength{5,10}`. Self-check stays 480/480; 307 tests,
clippy and fmt clean.

## 10. Out of sample

Nothing in this section chose a constant.

`mergelab.battery(seed, n, seglen)` draws 12-word canvases with IBS0 density and placement,
`inf1` content, mismatch load and missingness all randomised, holds them IBD2-free so the
printed `IBD1Seg` is the chromosome-2 IBD1 call alone, drives the reference at one floor and
grades the model. The spacing is the §0 spacing, so most canvases contain at least one
merge-eligible interruption.

```
  seed 424242   L=5      60 / 60      seed 424242   L=10     60 / 60
  seed 777001   L=5      60 / 60      seed 777001   L=10     60 / 60
  seed 31415    L=5      60 / 60      seed 31415    L=10     60 / 60
                                                     -> 360 / 360
```

Three unused seeds, **at 5 and 10 Mb specifically** — the previous campaigns' held-out
batteries were almost all at the default floor, which is exactly why this clause survived
them. The `(bad, X)` law of §4–§5 was additionally scored on 600 independently drawn one-
and two-word interruptions (seeds 7001, 7002): 600 of 600.

The IBD2 mirror was validated the same way in a regime where no merge can fire (88 000 bp
spacing, `L = 3`, so a one-word gap is 5.72 Mb): **60 of 60**, which is what says the
transcription of `19-…`'s caller into `mergelab.ibd2` is faithful before the merge is put on
trial. With the merge live at 20 000 bp spacing it is 56–58 of 60 — see §11.

## 11. What is still open

**Superseded by `docs/research/21-push-merge.md`.** Items 1 and 2 were the whole of that
campaign and are closed: the push is conditional (`21-…` §2), the IBD2 merge has no word
cap (§4), its interruption runs between gate windows rather than between runs (§3), and
its `X` is the HetHet count rather than `inf2` (§5) — so §7's last clause above is the one
reading in this document that `21-…` retracts. Item 3 is closed the other way: §3's
two-word cap is confirmed on the **IBD1** pass on a fresh fixture, and it is the IBD2 pass
that has none. Items 4 and 5 stand.

1. **The IBD2 merge is not closed.** Its gap rule, its cap and its budget are bisected
   (§7) and it is worth 1 row at 5 Mb and 8 at 10 on the corpus, but on random IBD2-native
   canvases it is 56/60 at `L = 3` and 58/60 at 5 and 10, in **both** directions. The
   residual is not the budget: two of the failures are canvases where the reference's
   `IBD2Seg` *appears* as the floor rises past 2.58 Mb in a configuration whose two-word
   gap should need 2.59, and two where it vanishes at a floor its call clears by a
   megabase. Both smell like an interaction with the one-word **push** (`17-…` §6) — which
   is counted over gate-passing calls whether or not they survive the floor — rather than
   with the merge. That is the first thing to look at.
2. **`bigish`, `missing` and `multifam` at 5/10 Mb** are the five parity cases still open.
   The residual is small, concentrated, and at 5 Mb **one-sided**, which names its cause:

   ```
   === 3 Mb ===   982 rows   ibd1 0            ibd2 0
   === 5 Mb ===   982 rows   ibd1 23 (+23/-0)  ibd2 35 (+0/-35)
        missing   14 rows    ibd1  0           ibd2  1 (-1, max 0.0111)
        bigish   763 rows    ibd1 23 (+23/-0, max 0.0042)   ibd2 34 (-34, max 0.0053)
   === 10 Mb ==   982 rows   ibd1 22 (+8/-14)  ibd2 37 (+14/-23)
        multifam 104 rows    ibd1  5 (-5, max 0.0277)   ibd2  4 (+4, max 0.0168)
        missing   14 rows    ibd1  0                    ibd2  1 (-1, max 0.0111)
        bigish   763 rows    ibd1 17 (+8/-9, max 0.0096)   ibd2 32 (+10/-22, max 0.0077)
   ```

   At 5 Mb **every** wrong `IBD2Seg` is too low and **every** wrong `IBD1Seg` is too high,
   on the same dataset — which is one fault, not two: the reference merges IBD2 where this
   caller does not, and the extra IBD2 is then subtracted from `IBD1Seg` (`18-…` §6). So
   the 5 Mb residual is entirely item 1 above, seen from the corpus. At 10 Mb both columns
   go both ways and `multifam`'s five `IBD1Seg` rows (all too low, worst 0.0277) are the
   largest single block left; `missing` is the cheapest place to start, being 14 rows with
   one wrong one at both floors.
3. **"Two words" or "129 marker intervals"?** Indistinguishable on a word-aligned scan
   (§3). A fringe canvas whose usable segment starts mid-word could in principle make a run
   gap that is not `64j + 1`, and would separate them.
4. **Whether a gate-refused run's own markers enter `X`.** §6 leaves it undecided; the
   committed rule sums only the unusable words. A canvas with a gate-refused run carrying
   close to 9 A1A1/A1A1 markers between two interruptions whose budget is exactly on the
   boundary would settle it.
5. **Why 4 and 3, and why the switch at 10.** The constants are bisected but not explained.
   The shape `cost · (bad − 2) ≤ informative` reads like an error-tolerance budget, and 10 is
   the informativeness gate; that is as far as measurement goes.

## 12. Reproducing everything here

```bash
cd docs/research/fixtures
python3 -c "import mergelab as M; print(M.battery(424242, 60, 5.0)[:2])"   # §10
python3 -c "import mergelab as M; print(M.battery2(90210, 60, 3.0)[:2])"   # §11.1

cd ../../../tests/parity/fit
python3 seg20.py            # §9 — the corpus scorecard at 3 / 5 / 10 Mb
python3 seg20.py grid       # every knob of the merge swept
python3 where20.py rows     # the residual, row by row
```

`mergelab.py` drives the KING 2.3.2 reference by default, or whatever `$KING` names, and
caches every answer in `mergelab_measured.json`. As with the other two rigs: do not point
`$KING` at a non-reference build in this directory, because it writes whatever it measures
into that cache. Grade a copy.
