# The `.seg` IBD2 caller, measured on a canvas built out of opposite homozygotes

**Status: measurement, fitted, validated out of sample — not committed to Rust.**
`docs/research/16-segment-extension.md` §10 said the next campaign's first job was to
build a `.seg`-native canvas, because the `--ibs` rig drives `IBD2Seg` to zero. That
canvas is `docs/research/fixtures/segcanvas.py`, and this is what it says.

No KING source was read. Every rule below is a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) run on
filesets built for the purpose, or a score against the captured parity corpus.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| The canvas works. Uniform chr2 spacing chosen so `D` lands just over the reference's 100 Mb floor makes one ulp of the printed `IBD2Seg` a fifth of a marker gap, so the printed column reads back **the number of marker intervals called, to a tenth of a marker** — and, when the calls are word-aligned, the number of words *and* the number of calls. | §1 |
| A **new hard constant**: the reference refuses a fileset whose usable total is under **100 000 000 bp** ("Segments too short."), bisected to the base pair. | §2 |
| The `.seg` IBD2 word predicate is **`IBS0 == 0` and `het-vs-hom mismatches ≤ 1`** — not `≤ 4` (the committed rule) and not `mismatches == 0` (`13-…` §1). Bisected two ways and confirmed by the corpus. A **lone** unusable word is absorbed, but only when the run picks up cleanly on the other side. | §3, §7 |
| The gate is `inf2 ≥ 10`, **counted from the run's first mismatch-free word**. The `≥ 10` is bisected at 9/10 on HetHet and on A1A1/A1A1 separately; the starting point is what makes `zx` a call and `xz` not. | §4 |
| The endpoints are **not** word-aligned. A call reaches **63 markers past the nearest het-vs-hom mismatch** in the word that bounds it — right past the *first*, left before the *last* — and a word carrying any opposite homozygote blocks the reach **whole-word**. Swept over all 64 bit positions. | §5 |
| Every call after the first **in the same usable segment starts one word late**, and it is the emitted call that causes it, not the break: a run refused by the gate pushes nothing, and the push survives its cause being dropped by `--seglength`. | §6 |
| **The sharp negative: `.seg` is not a quantised confirmation scan.** There is no chunk quantum, no confirmation count and no HetHet threshold. Uniform blocks at one mismatch per word are called in full at every width; the staircase that defines `--ibs` does not exist here. | §8 |
| Out of sample on the corpus: **709 of 982 exact rows** against 705, **`IBD2Seg` exact on 896 rows against 822**, mean `PropIBD` error **0.00037 against 0.001376**, worst row **0.0089 against 0.2109**, and the reported pair set still exactly right (0 extra, 0 missing). | §10.1 |
| Held out: on 160 random canvases from two seeds never used in the fit the rule reproduces the reference on **147 (92 %)**, where the committed geometry reproduces **18 %**. | §10.2 |

---

## 1. The instrument

`docs/research/fixtures/segcanvas.py`.

* **chr1 — the carrier.** Five complete words that are IBD1-clean (no opposite
  homozygote) and IBD2-dirty (34 het-vs-hom mismatches per word), carrying just enough
  `inf1` to clear the gate. It gives the pair one IBD1 segment of 10.527 Mb, which is what
  earns it a `.seg` row at all, and it contributes nothing to `IBD2Seg`.
* **chr2 — the canvas.** Complete words painted marker by marker from an explicit
  composition, walled at both ends by all-IBS0 words and padded out to a fixed word count
  so `D` is constant across a family.
* **the ruler.** chr2's spacing is uniform and chosen so `D` clears §2's floor by as
  little as possible — 100.551 Mb at the default 16-word canvas. One ulp of the printed
  `%.4lf` is then `D/10000 = 10 055 bp`, against an 88 000 bp marker gap, so

  ```
  markers called = IBD2Seg × D / spacing        to ±0.06 markers
  ```

  A word-aligned call over `n` words measures `64n − 1`, so a total `M` from `c` calls over
  `w` words is `64w − c`: **`c = (−M) mod 64` recovers the number of calls and the number
  of words exactly.** That is the whole advance over every previous `.seg` rig, which could
  only see an aggregate.
* **a graded ruler**, used where a count is not enough: per-word gaps `88 000 + 3 000 k`,
  so two candidate callings of the same width have different lengths and the printed total
  picks one (`segcanvas.pick`).
* **`--seglength`** separates calls where the count is ambiguous.

Two checks before anything else: the reference's own `allsegs.txt` agrees with the
computed `D` to the printed precision on every canvas, and an all-wall canvas reports
`IBD2Seg 0.0000` with `IBD1Seg` equal to the carrier alone — so chr2 is silent unless
painted.

**One trap worth recording.** KING's "Too many first alleles as the major allele" QC
samples with an **unseeded** RNG and fatals at random, and repeated invocations inside one
second reproduce the same draw. `segcanvas` retries with a one-second gap, and the carrier
keeps A1A1/A1A1 markers scarce because those are what the check trips on.

## 2. A new hard constant — the 100 Mb usable-total floor

Under the floor the reference prints `Segments too short.` and writes no `.seg` at all.
Solving `D = 319·sp1 + 895·spacing` for exact totals puts it on the base pair:

| `D` | result |
| ---: | --- |
| 99 999 999 | no `.seg` |
| **100 000 000** | analysed |
| 100 000 001 | analysed |

open-king does not model this today; it is one line and it is not part of the current
parity gap, because every corpus dataset is far above it.

## 3. The word predicate

A block of 8 pure-HetHet words with `j` consecutive words replaced by `(m mismatches,
64 − m HetHet)`, read back as (words called, calls):

```
 j\m      0       1       2       3       4       5       8      16      64
  1     8/1     8/1     8/1     8/1     8/1     8/1     8/1     8/1     8/1
  2     8/1     8/1     6/2     6/2     6/2     6/2     6/2     6/2     6/2
  3     8/1     8/1     5/2     5/2     5/2     5/2     5/2     5/2     5/2
```

**A word with one het-vs-hom mismatch is usable; two make it unusable.** The bisection is
unique — `m = 1` never splits at any `j`, `m = 2` always does. The `j = 1` row says the
second half of the predicate: a **lone** unusable word is absorbed whatever it carries (64
mismatches and the block is still one call over all eight words), while two in a row are
not, and both are then dropped from the call rather than shortening it. The absorption is
conditional and §7 states the condition; here every word after the lone one is a full
HetHet word, which is the easy case.

Opposite homozygotes are absolute: one IBS0 anywhere in a word
makes it unusable at any HetHet, and an IBS0 word is never absorbed. Missing calls,
A2A2/A2A2 and A1A1/A1A1 markers are all irrelevant to the predicate; a het-vs-A1A1
mismatch counts exactly like a het-vs-A2A2 one.

The corpus agrees, and its gradient is sharp — the `IBD2Seg` column is exact on

```
dirty threshold   1       2       3       5
IBD2Seg exact    867     896     868     837
```

**The IBD1 caller is confirmed unchanged**, on the same canvas with `inf1`-carrying words
in place of the HetHet ones (`segcanvas.py 3`). One opposite homozygote breaks an IBD1 run;
34 or even 64 het-vs-hom mismatches in a word do not, and neither do 64 missing calls. Its
asymmetric refinement reads back to the marker: with a single IBS0 at bit 0 of the word
after the run the call ends at `64·(v+1) + 0` and with two at bits 0,1 it ends at
`64·(v+1) + 1` — the flanking word's **last** IBS0 — while a run that opens after such a
word starts one marker past it. That is exactly `Scan::right_end`/`Scan::left_end`, and it
is *not* the `.seg` IBD2 rule of §5: the two passes really do have different geometry.

## 4. The gate

One word carrying exactly `k` informative markers, walled:

| `k` HetHet | 0 | 8 | 9 | **10** | 11 | 20 |
| --- | --- | --- | --- | --- | --- | --- |
| called | no | no | no | **yes** | yes | yes |

and the identical bisection on A1A1/A1A1 markers (9 no, 10 yes), and 5 HetHet + 5 A1A1 in
one word accepted where 4 + 5 is refused, and two words at 5 each accepted where two at 4
are refused. So the statistic is `inf2` — both samples carry A1 — exactly as
`13-informativeness-gate.md` measured for IBD1's `inf1`, and the threshold is the same 10.

**Where the count starts is the new part.** Writing `C` for a pure HetHet word, `z` for a
pure A2A2 word (clean, `inf2 = 0`) and `x` for one mismatch + 63 HetHet:

```
z   refused        zx   called (2 words)       zzx   called (3 words)
x   refused        xz   refused                xzz   refused
C   called         xC   called (2 words)       xzx   called (3 words)
```

`zx` and `xz` are the same two words with the same `inf2` total and the same run, and one
is a call and the other is not. **The count starts at the run's first mismatch-free word**
and runs to the end of the reach (§5): `xz` starts counting at the `z` and finds 0, `zx`
starts at the `z` and picks up the `x`'s 63. A run with no mismatch-free word at all is
refused outright, which is why *every* uniform block at one mismatch per word reports
nothing however wide it is (§8).

On the corpus this clause alone is worth 4 exact rows and more than halves the mean error
(`gate_from_clean=False`: 705 exact, MAE 0.00091, worst 0.2000).

## 5. The endpoints

A block of 6 HetHet words (383 markers on its own) with one boundary word beside it, whose
mismatches sit at named bits:

```
flanking word      left of the block      right of the block
bits 0,1                   +126                 +64
bits 30,31                  +96                 +94
bits 62,63                  +64                +126
bits 0..19                 +108                 +64
bits 44..63                 +64                +108
all 64                      +64                 +64
```

```
extension left  = 127 − (last mismatch bit of the word before the run)
extension right =  64 + (first mismatch bit of the word after the run)
```

which is the same statement twice: **the call reaches 63 markers past the nearest
het-vs-hom mismatch on that side.** Left end = `(marker of the last mismatch before the
run) − 63`; right end = `(marker of the first mismatch after the run) + 63`. It is a
marker-level rule, and 63 is pinned by the bit sweep to a unique integer (bit 31 gives
exactly +96, so 62 would give 95 and 64 would give 97).

**Opposite homozygotes block it whole-word, not marker by marker.** A flanking word
carrying one IBS0 stops the call dead on the run's own last marker whatever bit that IBS0
sits at (0, 32 and 63 all give +0), and an IBS0 in the *second* word out caps the reach at
that word's boundary (+64). This asymmetry — marker-level for mismatches, word-level for
IBS0 — is exactly why `14-ibd2-geometry.md` §6.2's "two words and one marker" looked
unexplainable: it is one word plus the 63-marker reach.

A call that reaches the usable segment's own first or last complete word runs on to the
segment's first or last **marker**, which is the fringe rule the committed engine already
had; the duplicate pair in `dups` reads `IBD2Seg 1.0000`, not the word-aligned 0.8984.
The canvas cannot see this (chr2 is word-aligned), the corpus can: it is worth 37 exact
`IBD2Seg` values and takes the worst row from 0.0508 to 0.0071.

## 6. The push

Two runs separated by dirty words, read back as markers:

| canvas | reference | model, no push | model, push 1 word |
| --- | ---: | ---: | ---: |
| `C C y y C C` | **254** | 383 | 254 |
| `C C y y y C C` | **254** | 444 | 254 |
| `C C W W C C` | **190** | 254 | 190 |
| `C W C C C C C C` | **382** | 446 | 382 |
| `z0 W C C C C C C` | **383** | 383 | 383 |

(`y` = 2 mismatches, `W` = an all-IBS0 wall, `z0` = a clean word with `inf2 = 0`.)

**Every call after the first in a usable segment starts one word later** — whatever
separates them, and however wide the separation. The graded ruler says which word is lost:
the *first* of the later run, never the last of the earlier one (`C W C6` decodes to
`[0,0] + [3,7]`, not `[0,0] + [2,6]`). The word it counts from is the run's **gate-start**
word — the first word with no mismatch (§4) — not the run's first word; `CyxC`, whose
second run opens on a mismatch word, is the case that separates the two.

Two controls make the clause precise. It is the **emitted call** that pushes, not the
break: in `z0 W C6` the first run fails the gate, nothing is emitted, and the second run
is not pushed. And the push **survives its cause being dropped**: at `--seglength 6` the
1-word first call (5.54 Mb) disappears from the total and the second call is still short by
one word — so the clip is applied before the length filter, as `13-…` §1 already had it.

The counter is **per usable segment**, not per pair: making the carrier IBD2-callable puts
a call on chr1 and chr2's first call is not pushed.

## 7. The rule

Per pair, per usable segment `[w0, w1]` covering markers `[lo, hi]`, over the retained
autosome 64-marker word grid, from three per-word quantities — the IBS0 mask, the
het-vs-hom mismatch mask, and the `inf2` popcount.

```
a word is USABLE iff it carries no IBS0 and at most one het-vs-hom mismatch
a lone unusable word carrying no IBS0, between two usable words, is ABSORBED iff
    the word after it carries no mismatch at all, and the usable words from there on
    carry >= 10 inf2 between them
runs are the maximal stretches of usable (or absorbed) words

for each run [a, b], in order:
    left  = 64a ; right = 64b + 63
    if word a-1 is inside the segment and carries no IBS0:
        left = (marker of the last mismatch in word a-1) - 63
        if a-2 is outside the segment or carries an IBS0: left = max(left, 64(a-1))
    if word b+1 is inside the segment and carries no IBS0:
        right = (marker of the first mismatch in word b+1) + 63
        if b+2 is outside the segment or carries an IBS0: right = min(right, 64(b+2) - 1)

    g = the first word of [a, b] with no mismatch      # no such word -> refuse the run
    if inf2 over words [g .. right/64] < 10:  refuse the run

    if a call has already been emitted in this segment:  left = max(left, 64(g+1))
    if a == w0: left  = min(left, lo)                  # the segment's own fringes
    if b == w1: right = max(right, hi)
    clamp to [lo, hi], then to the previous call:  left = max(left, previous right)
    keep it if pos[right] - pos[left] >= seglength
```

`predict()` in `docs/research/fixtures/segcanvas.py` is this, over per-word
`wordinfo` tuples; `ibd2_17()` in `tests/parity/fit/seg17.py` is the same rule over the
corpus's bit planes, and the two agree on every canvas.

## 8. The sharp negative — `.seg` is not a confirmation scan

The obvious hypothesis after `16-…` was that `.seg` runs the same **quantised confirmation
scan** as `--ibs` with different constants: a chunk of `N` mismatches, confirmed by `H`
markers over `Y` words, with the interval cut at the last confirmed chunk. It does not.

* **There is no chunk quantum.** `--ibs`'s signature is the staircase of §4 of `16-…`:
  sweeping the width `W` of a uniform block moves the reported end in steps of `⌈5/m⌉`
  words. Here, sweeping `W` over 1…10 on a uniform block of clean words gives exactly
  `64W − 1` markers at every `W`, with no staircase at any `(m, h)` that is called at all.
* **There is no confirmation count.** A block of `W` words at one mismatch each is
  refused at *every* `W` from 2 to 10 and at every HetHet from 0 to 63 — and the same
  words with a single clean word in front of them are called **in full**, all `W + 1` of
  them, at every `W`. A confirmation scan cannot do that: the credit a chunk needs would
  have to grow with the block, and here one clean word is worth an unbounded tail.
* **HetHet is not the currency.** `--ibs` counts HetHet and ignores A1A1/A1A1 entirely
  (`15-…` §5). `.seg` counts them interchangeably: 10 of either passes the gate, 9 of
  either fails (§4).
* **Nothing is ever cut.** Across 340 exhaustive canvases and 240 random ones, a run that
  is called is called over its whole extent plus the §5 reach. The `--ibs` behaviour of
  reporting *part* of a uniform block has no analogue here: a call is the whole run or
  the run is refused.

So the two callers are not the same machine under different constants, and the shape
`16-…` §9 warned about — porting the chunk geometry — was not merely mis-tuned, it was the
wrong kind of rule. That port scored 709 exact rows with MAE 0.00356 and a worst row of
0.1490; this rule scores 709 with MAE 0.00037 and a worst row of 0.0089.

## 9. The constants, and how each is pinned

| constant | value | how |
| --- | ---: | --- |
| usable-total floor | 100 000 000 bp | exact bisection, §2 |
| IBD2 word predicate | IBS0 = 0 **and** mismatches ≤ 1 | §3: `m = 1` never splits, `m = 2` always does, at every run width; corpus `IBD2Seg` 896 at 2 against 867/868/837 at 1/3/5 |
| gate | `inf2 ≥ 10` | §4: 9 refused / 10 accepted, four independent bisections (HetHet, A1A1/A1A1, mixed, split across two words) |
| gate start | the run's first mismatch-free word | §4: `zx` called / `xz` refused at equal totals; corpus 709 exact vs 705, MAE 0.00037 vs 0.00091 |
| reach | 63 markers | §5: extension = `127 − last bit` / `64 + first bit`, swept over bits 0…63; corpus `IBD2Seg` 896 at 63 against 824/825/830/826 at 0/32/62/64 |
| push | 1 word, from the gate-start word | §6, direct; corpus `IBD2Seg` 896 against 835 without |
| clip between calls | calls may touch | corpus cannot see it (identical either way); `CyC` needs it on the canvas |
| bridging condition | next word mismatch-free **and** `inf2 ≥ 10` after it | the least-supported clause — §11.1 |

`segcanvas_measured.json` caches every reference answer used above (872 invocations),
so all of it re-runs in seconds without the binary.

## 10. Out of sample

Nothing in this section had any part in choosing a constant.

### 10.1 The corpus — 982 `.seg` rows over ten datasets

```
$ python3 tests/parity/fit/seg17.py
committed (engine.py)   exact  705  both  820  ibd2  822  MAE 0.00138  bias +0.00000  worst 0.2109
17 fitted               exact  709  both  825  ibd2  896  MAE 0.00037  bias +0.00034  worst 0.0089
```

`extra 0, missing 0` in both: the set of pairs reported is still exactly right. Per
dataset the exact-row counts are unchanged everywhere except `monomorphic`, 8 → 12 — and
that is not a shrug, it is the point: the reference reports `IBD2Seg 0.0000` for
`P_C1/P_C2`, `P_C1/P_C3`, `P_C1/P_C4` and `P_C3/P_C4`, where the committed rule invents
0.42, 0.27, 0.03 and 0.08 of the genome out of words carrying two to four mismatches. The
`dups` MZ pair moves from 0.9718 to 0.9207 against a reference 0.9223.

The headline is not the exact-row count, which `14-…` §2 already showed is a bad gradient.
It is that **`IBD2Seg` is now exact on 896 of 982 rows instead of 822**, that the mean
`PropIBD` error falls by a factor of 3.7, and that the worst row falls by a factor of 24 —
0.2109 → 0.0089. There is no row left where the caller is grossly wrong.

### 10.2 Held-out random canvases

`segcanvas.battery(seed, n)` draws word sequences with IBS0 density, mismatch count,
mismatch placement and informative content all randomised, drives the reference, and
compares the printed total to `predict()`:

```
seed 101 (used to read misses during the fit)   76 / 80
seed 777   (never used)                         71 / 80
seed 8081  (never used)                         76 / 80
------------------------------------------------------------
held out (777 + 8081)                          147 / 160   92 %
all three seeds, this rule                     223 / 240   93 %
all three seeds, the committed geometry         42 / 240   18 %
```

and on the 340-canvas exhaustive battery of `segcanvas.py 5` (every word sequence of
length ≤ 4 over `{C, z, x, y}`), this rule reproduces **329**, the committed geometry **4**.
Thirty-eight further named families — the whole predicate table of §3, both endpoint
families of §5, every push canvas of §6 and every interior-IBS0 canvas — agree **38 / 38**.

### 10.3 What it is not

**It was not committed to Rust.** The rule is right about the corpus's numbers but wrong
about one canvas in eight, and `crates/king-core/src/ibdseg.rs` is a byte-parity engine;
`tests/parity/run_parity.py` grades whole files, and at 709 of 982 rows no `*__ibdseg`
case becomes byte-identical, so landing it moves the harness by zero. That is the only
reason it is not in `crates/king-core/src/ibdseg.rs`: on every gradient that can see it —
mean, tail, and the `IBD2Seg` column — it beats the committed rule by a wide margin, and
`seg17.py` makes the port mechanical once §11.1 closes.

## 11. What is still open

1. **The bridging condition is the least-supported clause.** "The next word carries no
   mismatch and the usable words from there carry `inf2 ≥ 10`" is fitted from about thirty
   canvases and it is the only clause with a free shape rather than a bisected integer.
   The corpus likes it — `IBD2Seg` exact on 896 rows against 874 with no bridging at all —
   but it costs a little mean (MAE 0.00037 against 0.00033) and a little tail (worst 0.0089
   against 0.0071), which is what a clause that is *nearly* right looks like. The lookahead
   in particular is guessed: it stops at the next unusable word, and `Cyzy` (§11.2) says it
   should not.
2. **11 of the 340 exhaustive canvases**, all of them a lone 2-mismatch word with a run
   continuing past it. Two shapes: `Cyzy` and `zyzy` are reported over all four words where
   the model splits them (the lookahead above), and `xyC`, `xyCC`, `xyzC` … are reported
   two markers *shorter* than the model — a left end at `64·g − 62` rather than `64·a`,
   i.e. the §5 reach applied at the gate-start word rather than the run start. That second
   shape is a two-marker effect on four-word canvases and is invisible on the corpus.
3. **The 8 % random residual**, concentrated on canvases with several short runs separated
   by single dirty words — again item 1. Nothing in the corpus looks like those sequences,
   which is why the corpus is nearly exact and this is not, and why 92 % must not be read
   as "solved".
4. **The push's mechanism.** Measured to the word and reproduced everywhere, but it is
   stated as "the second and later calls start one word late" because no clipping
   formulation reproduces it: the gap between the runs is irrelevant, so it is not
   `max(a, previous_end + k)` for any `k`.
5. **The 100 Mb floor** (§2) is not implemented in open-king. One line, no corpus impact,
   worth adding for a user who runs a small fileset.
6. **The residual left after this rule is in the *IBD1* caller, not the IBD2 one** —
   measured from the `--build` side, §13.

## 12. Reproducing everything here

```bash
cd docs/research/fixtures
python3 segcanvas.py 0     # §1, §2 — the rig, the floor, the marker ruler
python3 segcanvas.py 1     # a pure-HetHet block of W words is exactly 64W-1 markers
python3 segcanvas.py 2     # §3 — the IBD2 word predicate
python3 segcanvas.py 3     # §3 — the IBD1 word predicate
python3 segcanvas.py 4     # §5 — the endpoints, and what blocks the reach
python3 segcanvas.py 5     # §10.2 — the exhaustive length-<=4 battery
python3 segcanvas.py 6     # §4 — the gate
python3 segcanvas.py 7     # §6 — the push
python3 segcanvas.py 8     # §10.2 — the random batteries, this rule against the old one

cd ../../../tests/parity/fit
python3 seg17.py           # §10.1 — the corpus scorecard
python3 seg17.py grid      # §9 — every constant swept against the corpus
```

`segcanvas.py` drives the KING 2.3.2 reference by default, or whatever `$KING` names, and
caches every answer in `segcanvas_measured.json`; `seg17.py` reads the captured corpus and
needs no binary at all.

---

## 13. Where the residual now lives: the *union*, and therefore the IBD1 caller

Measured while chasing `--build`'s `INFERENCE AV.FS` (`docs/PARITY.md` §6.2), because the
statistic that case needs reads the **union** `IBD1 ∪ IBD2` rather than either column.
That union is observable in the reference's own output: since a call's `IBD1Seg` is the
part of it not already IBD2, `IBD1Seg + IBD2Seg` **is** the union as a fraction of `D`.

Scored against every corpus row, with the rule of §7 committed:

| quantity | exact rows |
| --- | ---: |
| `IBD2Seg` alone | 896 / 982 |
| both columns | 825 / 982 |
| **the union `IBD1Seg + IBD2Seg`** | **826 / 982** |
| the union, on rows whose reference `IBD2Seg` **is 0** | **823 / 823** |
| the union, on rows whose reference `IBD2Seg` is **> 0** | **3 / 159** |

So §7 bought 74 rows of `IBD2Seg` and left the union where it was. Splitting the 156
union-wrong rows by whether our IBD2 call sticks out past our IBD1 call:

* **135 of 156** have IBD2 lying entirely inside IBD1, so their union is our IBD1 call's
  extent and nothing else — the IBD2 geometry cannot be blamed for them;
* on **0 of 156** would using the IBD1 total alone have matched, so the 21 rows where IBD2
  does stick out are not the explanation either;
* and every union-wrong row has reference `IBD2Seg > 0`, while **every** row with
  `IBD2Seg == 0` is exact.

That is the shape of an *interaction*: our IBD1 caller never looks at IBD2, is exact on all
823 IBD2-free rows, and is wrong on 156 of the 159 rows where the pair also carries IBD2.
Whatever the reference does to an IBD1 run in the presence of an IBD2 one, it is not the
obvious candidate — clipping IBD1 against IBD2 and re-applying the 3 Mb minimum to the
fragments is a **no-op** on this corpus at 1, 2 and 3 Mb (no fragment is ever short enough
to drop) and costs 20 rows at 5 Mb.

Reproduce: `python3 tests/parity/fit/seg17.py` for the column scorecard, and
`docs/research/fixtures/avfs_score.py` for the union columns and the `--build` consequence.
