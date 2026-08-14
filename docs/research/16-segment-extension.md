# How an IBD2 run is started, extended and ended — solved for `--ibs`

**Status: COMMITTED.** `Scan::ibd2_words` in `crates/king-core/src/ibdseg.rs` is the rule
below; the harness went 397 → **403 of 480** with zero regressions (§8.3). Solved for the
`--ibs` IBD2 caller — exactly on the corpus, about 93 % out of sample on adversarial
constructed sequences (§10.3). **Not** solved for the `.seg` caller, which is a different
caller and is untouched (§9).

No KING source was read. Every rule below is either a reading taken off the reference
binary `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2)
run on filesets built for the purpose, or a score against the captured parity corpus.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| The residual of `15-ibs-ibd2-rules.md` §6 is **not** a running score, not an X-drop, and not a Viterbi path. It is a **quantised confirmation scan**: the run is confirmed in *chunks of five het-vs-hom mismatches*, and a chunk is confirmed iff it spans **≥ 3 words** and carries **≥ 95 HetHet** markers. | §3–§6 |
| Every constant is pinned to a unique integer by bisection: dirty-word threshold **5**, chunk **5** mismatches, chunk **95** HetHet, chunk **3** words, measured-interval overhang **1** mismatch. | §7 |
| Out of sample on the corpus the rule reproduces **`MaxIBD2` 158/158** and **`Pr_IBD2` 158/158** — both `--ibs` IBD2 columns, exactly, on every pair the reference grades. The rule it replaced scored 148/158 and 100/158. | §8 |
| On 658 constructed filesets, including 200 *random* word sequences the rule was never shown, it reproduces the reference's called interval on **614 (93.3 %)**; a fresh 200-canvas battery scores **186 (93.0 %)**. That residual ships — §10.3. | §8.2, §10.3 |
| The `.seg` IBD2 caller is **a different caller**, not just a different ruler: its word predicate refuses a word with even one het-vs-hom mismatch, so none of these fixtures produce a `.seg` call at all. Porting the chunk geometry to `.seg` naively gains 4 exact rows and loses on MAE. | §9 |

---

## 1. What the previous agent's lead asked, and the answer

`15-ibs-ibd2-rules.md` §6 left a boundary that was non-linear in `(mismatches, HetHet)`
and order-dependent, and proposed three shapes: **(a)** a running score with an X-drop,
**(b)** a two-state HMM/Viterbi path, **(c)** a greedy extension from a seed.

**(a) and (b) are ruled out by measurement, not by preference.**

* Both are *aggregate-optimal* rules over the words they see. A uniform block of `W`
  identical words bounded by walls admits no interior optimum under either: an X-drop
  score is monotone in the block, and a Viterbi path with fixed transition costs is
  all-IBD2 or no-IBD2. The reference reports **partial** intervals of a uniform block
  (§3), and the length of the part depends on `W` (§4), which no forward-causal score
  can produce.
* A running score also predicts that credit accumulates. It does not: a 20-word prefix
  of all-HetHet words (1 280 HetHet markers) buys the tail **exactly one** chunk of five
  mismatches, never two (§5). The counters are reset, not carried.

**(c) is the closest of the three and still not right.** The scan *is* one forward pass,
and a run does absorb words that could not have started it — but the unit of extension is
not a word, it is a *chunk of five mismatches*, and a chunk is accepted or refused as a
whole. That single change explains every measurement the previous agent recorded.

## 2. The instrument

`docs/research/fixtures/segfit.py`. chr1 is a 60-word IBD1 carrier that keeps the pair
above `--ibs`'s gates; chr2 is a canvas painted one complete word at a time from an
explicit composition. Two rulers:

* **per-word spacing** — word `w` of chr2 gets marker gap `40 000 + 137 w` bp, so the
  reported `MaxIBD2` inverts to exactly one word interval `[u, e]`. This is the whole
  advance over the previous rig, which could only see *whether* a pair was reported.
* **per-marker spacing** — gap `40 000 + 17 i` after marker `i`, so a length inverts to a
  unique *marker* interval. Used once, to confirm that both endpoints are exactly word
  aligned (`[64u, 64e+63]`) and that mismatch offsets inside a word never move them.

Two facts checked before anything else: KING keeps monomorphic markers (a 60 + 14 word
canvas reports "stored in 74 words"), so the word grid is exactly as constructed; and
`Pr_IBD2 × D` equals `MaxIBD2` to inside one printing ulp on every uniform-block fixture,
so each of these canvases produces exactly **one** call and the decoded interval is that
call.

```bash
cd docs/research/fixtures
python3 segfit.py 0 1 2 3 4 5 6      # ~1 500 reference invocations, a few minutes
```

## 3. The reported interval is a *part* of the run

Eight words each carrying `m = 1` mismatch and `h = 19` HetHet, bounded by walls, are
reported as the interval `[0, 5]` — six of the ten words, not all eight. Raising `h`
to 32 makes it `[0, 8]`. So the previous agent's "smallest reported `h`" table was not
measuring an acceptance threshold at all: it was measuring where the interval is **cut**.

| `m` | `h` | `W` | reported interval |
| ---: | ---: | ---: | --- |
| 0 | 12 | 8 | `[0, 8]` |
| 0 | 11 | 8 | refused |
| 1 | 19 | 8 | `[0, 5]` |
| 1 | 18 | 8 | refused |
| 1 | 32 | 8 | `[0, 8]` |
| 2 | 32 | 8 | `[0, 5]` |
| 2 | 62 | 8 | `[0, 8]` |

## 4. The cut position depends on where the *wall* is — so the scan is not word-local

Sweeping the block width `W` under a fixed composition gives a staircase, not a constant:

```
m=1 h=19    W:  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
            e:  5  5  5  5  9 10 10 10 10 14 15 15 15 15 19
```

`e = min(W - 1, 5⌊W/5⌋)` fits it exactly. The same nine words that are cut at `e = 5`
when `W = 9` are *not* cut when `W = 10`: adding one further identical word past them
moves the reported end from word 5 to word 9. A left-to-right rule that decides word by
word cannot do that, and this is what made the residual look like a running score.

The staircase's period is `⌈5/m⌉` words — the number of words needed to accumulate five
mismatches — for every `m` tried (5 words at `m = 1`, 3 at `m = 2`, 2 at `m = 3, 4`,
1 at `m = 5`). The reported interval is invariant to the block's absolute position, to
the carrier length, and to the number of walls on either side; and mismatch offsets
*inside* a word never move it (verified at offsets 0, 16, 30–33, 62–63).

## 5. Only the *trailing* part matters, and 95 HetHet buys one chunk of five mismatches

Put `j` words carrying one mismatch each into a 24-word all-HetHet block. At the front or
in the middle they cost nothing at all; only at the **end** do they cut the interval. So
the mechanism is a trailing trim, and it can be measured cleanly against an established
run: 20 all-HetHet words, then `j` words at `(m, h)`, then a wall. The figure below is
how many of the `j` words the reported interval covers.

```
        j:  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
m=1 h= 0:   0  1  2  3  4  5  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6
m=1 h=18:   0  1  2  3  4  5  6  6  6  6  6  6  6  6  6  6  6  6  6  6  6
m=1 h=19:   0  1  2  3  4  5  6  6  6  6 10 11 11 11 11 15 16 16 16 16 20
m=1 h=23:   0  1  2  3  4  5  6  6  6  6 10 11 11 11 11 15 16 16 16 16 20
m=1 h=24:   0  1  2  3  4  5  6  6  6  9 10 11 11 11 14 15 16 16 16 19 20
m=1 h=31:   0  1  2  3  4  5  6  6  6  9 10 11 11 11 14 15 16 16 16 19 20
m=1 h=32:   0  1  2  3  4  5  6  6  8  9 10 11 11 13 14 15 16 16 18 19 20
m=1 h=47:   0  1  2  3  4  5  6  6  8  9 10 11 11 13 14 15 16 16 18 19 20
m=1 h=48:   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
m=1 h=63:   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
```

Three readings, and they are the whole rule.

1. **A 20-word all-HetHet prefix buys exactly five trailing mismatches, never ten.** At
   `h = 0` the interval stops after the fifth mismatch whatever `j` is. The 1 280 HetHet
   markers of the prefix pay for one chunk of five and are then gone: the counters reset.
2. **The thresholds are 19, 24, 32, 48** — and `5 × 19 = 95`, `4 × 24 = 96`,
   `3 × 32 = 96`, `2 × 48 = 96`, against `5 × 18 = 90`, `4 × 23 = 92`, `3 × 31 = 93`,
   `2 × 47 = 94` refused. Four independent bisections; their intersection is
   **`HetHet ≥ 95`**, a unique integer, and it is the same 95 as `13-…`/`15-…` §5.
3. **The unit is the chunk of five mismatches.** At `h = 19` the covered count steps
   1…6, then 10, 11, then 15, 16 — every fifth mismatch, plus a one-mismatch overhang.

The **start** trims by the same machinery, mirrored: `j` words at `(m, h)` followed by 20
all-HetHet words report a start of `5⌊j/5⌋ − 1` at `m = 1, h ≤ 18`, `3⌊j/3⌋` at
`m = 2, h ≤ 31`, `2⌊j/2⌋` at `m = 3, 4`, and `0` at `m = 1, h ≥ 19` / `m = 2, h ≥ 32`.
A refused chunk closes the segment *and* starts the next one.

`m = 3` and `m = 4` are refused at the start **at every `h`**, including `h = 61`
(122 HetHet). A 2-word chunk with 122 HetHet is refused while a 3-word chunk with 96 is
accepted, so the chunk carries a **word-count** requirement as well as a HetHet one; the
grid search of §7 puts it at **≥ 3 words**, the same three-word floor `15-…` §5 measured
on the interval.

## 6. The rule

Per pair, per usable segment, over the retained-autosome 64-marker word grid, using two
per-word popcounts: `m_w` = het-vs-hom mismatches (`(het_i & hom_j) | (hom_i & het_j)`)
and `h_w` = HetHet (`het_i & het_j`). Opposite homozygotes and missing calls are
irrelevant to this pass at any density (`15-…` §2).

```
a word is DIRTY iff m_w >= 5
a run is a maximal stretch of non-dirty words, with a LONE dirty word bridged
   (two dirty words in a row end it)

for each run [a, b]:
    scan k from a to min(b+1, w1)        # one word past the run: that word fires the
    u    = a                             # counter, and its HetHet counts for the chunk
    mis  = het = 0 ; conf = none ; cstart = a
    at each k:  mis += m_k ; het += h_k
        if mis >= 5:                              # the chunk closes here
            if (het >= 95 and k - cstart + 1 >= 3)  # ...and is confirmed
               or (the run reaches the segment's last two words):
                    conf = k ; mis = het = 0 ; cstart = k+1
            else:                                  # ...and is refused
                    emit [u, extend(conf)] if conf exists
                    u = restart(k) ; mis = het = 0 ; conf = none ; cstart = u ; rescan
    when the scan runs out:
        if the run reaches the segment's last two words:  emit [u, w1]
        elif conf exists:                                 emit [u, extend(conf)]

extend(conf) = the furthest word past `conf`, within the run, reachable while picking up
               at most ONE further mismatch
restart(k)   = the word after the 4th mismatch's word, when the word that closed the
               chunk holds only the 5th mismatch; otherwise the word after k

then, in order: clip each interval against the previous one, drop any interval spanning
fewer than 3 words, and measure it word-aligned: pos[64e+63] − pos[64u].
```

Everything the previous agent recorded falls out.

* **The staircase** is the chunk quantum: the interval can only end where a chunk was
  confirmed (plus the overhang).
* **The `W`-dependence** is not non-causality. The run's end is confirmed by a chunk, and
  the chunk that would confirm word 8 is only *closed* by the fifth mismatch, which at
  `W = 9` lies in the wall (HetHet 0 → refused) and at `W = 10` lies in a normal word
  (HetHet 19 × 5 = 95 → confirmed).
* **The order-dependence** (`8 clean + 8 at m = 4` → one interval of sixteen; reversed →
  only the clean eight) is that HetHet is counted *within a chunk*, from wherever that
  chunk starts. Reversed, the leading `(4, 60)` words form 2-word chunks, which fail the
  three-word floor, and each failure restarts the segment.
* **`m = 1, h = 16` refused while `m = 2, h = 32` reported, at the same ratio**: the
  chunks are 5 × 16 = 80 < 95 and 3 × 32 = 96 ≥ 95. Nothing about the ratio matters.

## 7. The constants, and how each is pinned

| constant | value | how |
| --- | ---: | --- |
| dirty word | `m_w ≥ 5` | `15-…` §2, re-confirmed here: two words at `q = 5` split a run, at `q = 4` merge, at any HetHet |
| chunk closes at | 5 mismatches | the staircase period is `⌈5/m⌉` words for `m = 1…5` |
| chunk needs | 95 HetHet | four bisections, `(5,19) (4,24) (3,32) (2,48)` accept vs `(5,18) (4,23) (3,31) (2,47)` refuse → the intersection `(94, 95]` |
| chunk needs | 3 words | 2-word chunks refused at 122 HetHet; grid search over `{1,2,3,4,5}` scores 614/658 at 3, against 599 at 1–2 and 516 at 4 |
| interval overhang | 1 mismatch | grid search over `{0,1,2}` further mismatches scores 614/658 at 1, against 577 at 0 and 563 at 2 |
| interval spans | ≥ 3 words | `15-…` §5, unchanged |

Both grids are in the fit log of `docs/research/fixtures/segfit_measured.json` (658 cached
reference answers) and can be re-run without touching the binary.

The **restart** rule is the one clause the corpus cannot see: `restart = "after the
refusing word"` scores identically to the fitted rule on all 316 corpus targets. It is
fitted only from constructed filesets, where fourteen irregular patterns separate the two
(`[2,2,1] → 2` against `[2,2,2] → 3`, `[4,1] → 1` against `[4,4] → 2`). Treat it as the
least-supported part of the rule.

## 8. Out-of-sample accuracy

### 8.1 The corpus — both `--ibs` IBD2 columns, exactly

```
$ python3 tests/parity/fit/chunk.py
rule              MaxIBD2      Pr_IBD2    Pr bias
committed       148/158      100/158      +0.0029
chunk fit       158/158      158/158      +0.0000
chunk after     158/158      158/158      +0.0000
chunk at        158/158      157/158      +0.0000
```

`MaxIBD2` is one exact segment length per pair to the base pair; `Pr_IBD2` is the total
over all of that pair's calls at 4 dp. **158/158 on both.** These 316 targets span
`nuclear`, `threegen`, `multifam`, `dups`, `missing`, `monomorphic`, `sexchr`,
`unrelated`, `admixed` and `bigish`; none of them was used to fit any constant. One
detail *was* corpus-informed — that a chunk closed by the usable segment's own last word
is exempt from the HetHet test, found on `missing M_C2/M_C4` — and it is the same tail
exemption `15-…` §5 had already measured on a fixture, applied in the right place.

This closes `14-ibd2-geometry.md` §8.5 (`dups MZ_1/MZ_2`'s unlocalisable `MaxIBD2`),
§8.4 (the two denied bridges) and §8.3 (the four declined right extensions) — all of them
are chunk refusals.

### 8.2 Constructed filesets, including sequences the rule never saw

658 filesets: 120 uniform blocks, 200 **random** word sequences (each word's mismatch
count drawn from `{0,0,0,1,1,2,3,4,5,8,20,64}` and its HetHet uniform on what is left),
and 338 prefix/tail families. Agreement with the reference on the exact called interval:

```
agree 614 (93.3 %)   disagree 44
   of the 44:  end short by 1–10 words 30,  start or both wrong 12,  we call / ref refuses 2
```

A fresh random battery (`python3 segfit.py 6`, seed 5, 60 sequences the fit never saw)
scores **59/60**. The residual is concentrated on adversarial sequences that alternate
20- and 64-mismatch words with near-empty ones; nothing in the corpus looks like that,
which is why the corpus is exact and this is not.

### 8.3 What this was worth to the implementation — **committed**

`Scan::ibd2_words` in `crates/king-core/src/ibdseg.rs` now *is* §6. The harness moved
**397 → 403 of 480**, and the six cases are exactly `core/{nuclear,multifam,missing,
monomorphic,sexchr,bigish}__ibs` — every remaining `--ibs` failure, and nothing else.
Diffing the before and after FAIL lists gives **zero regressions**: the other 77 lines are
identical down to their per-file failure detail. As predicted, it closes no `*__ibdseg`
case, because those are decided by `king.seg` (§9).

`.seg` is untouched and unmoved: 705 of 982 exact rows at the default floor, mean absolute
`PropIBD` error 0.00138. `--ibs` is now byte-identical on all 13 datasets.

### 8.4 The port was checked against the model, not just against the corpus

Three independent checks, because "the corpus passes" only proves the port on the corpus:

* `tests/parity/fit/check_mirror.py` — the Python mirror `engine.py` (updated to §6 in the
  same change) reproduces the Rust binary's `.seg` columns **and** every `MaxIBD2` value it
  prints, across all 13 datasets: `MIRROR OK`. Before the mirror was updated it diverged on
  13 values, which is what a stale mirror is supposed to do.
* 200 fresh random word canvases driven through **our binary** and compared to `predict()`:
  **200 agree, 0 disagree.** (`fixlab.py` now honours `$KING`, so the same fixture rig can
  drive either binary; the docs already claimed it did.)
* `python3 segfit.py 6` against our binary scores **60/60** where the reference scores
  59/60 — i.e. our binary is the model, and the model is 93 % of the reference (§8.2).

Four unit tests in `ibdseg.rs` pin the chunk behaviour in-process, on a word-composition
canvas that mirrors this rig: the §4 staircase (nine words at `(1, 19)` cut at block word
5, ten words not cut at all), the 18-vs-19 HetHet bisection, the three-word chunk floor
`(3, 61)` refused against `(2, 48)` confirmed, and the counters-reset measurement of §5.

## 9. The `.seg` caller is a different caller — a sharp negative

The obvious hope is that `.seg` runs the same scan under a different ruler. It does not.

* Every fixture in §3–§5 reports `IBD2Seg 0.0000` under `--ibdseg` while `--ibs` reports a
  segment, because the `.seg` IBD2 word predicate refuses a word with **any** het-vs-hom
  mismatch (`13-segment-acceptance.md` §1: "≤ 0"), where `--ibs` tolerates four. So this
  whole fixture family is invisible to `--ibdseg`, and the chunk constants cannot be
  measured through it. A `.seg` campaign needs its own canvas, built out of opposite
  homozygotes rather than het-vs-hom mismatches.
* Porting the chunk geometry to `.seg` unchanged (`tests/parity/fit/segtry.py`) scores:

  | | exact rows | both columns | IBD2 column | MAE(PropIBD) | worst |
  | --- | ---: | ---: | ---: | ---: | ---: |
  | committed | 705 | 820 | 822 | 0.001376 | 0.2109 |
  | chunk geometry | **709** | **824** | **824** | 0.003561 | **0.1490** |

  Four rows better and the worst row much better, but the mean error nearly triples. That
  is what a *partly* right rule looks like: the chunk mechanism is doing real work on the
  rows that were badly wrong, and the `.seg` pass's own word predicate, informativeness
  gate and marker-level refinement are interacting with it in ways this port does not
  model. **It was not committed.**

## 10. What is still open

1. **The `.seg` IBD2 caller**, which is the whole remaining `ibdseg/*__ibdseg` gap. The
   method that worked here transfers directly: build a canvas out of opposite
   homozygotes (the mask `.seg` actually splits on), give each word its own spacing,
   drive `--ibdseg`, and read `IBD2Seg` back through the `--seglength` bisection of
   `14-…` §3.2 to recover individual segment lengths. Expect its own chunk constants;
   do not assume 5/95/3 carry over.

   **And note what §8.3 cost.** `--ibs`'s `Pr_IBD2`/`MaxIBD2` were the sharp graders for
   IBD2 work; they are now exact, so every candidate `.seg` rule scores an identical
   158/158 on them and they discriminate nothing. The `.seg` exact-row count was always a
   bad gradient (705 ± 1 under every variant, dominated by 823 IBD2-free rows). So the
   canvas above is not an optimisation of the next campaign — it is its **precondition**.
2. **The restart clause** (§7). Fitted from fourteen irregular patterns; the corpus
   cannot see it. It *is* now written into Rust — flagged in place in `ibd2_words` as the
   least-supported clause in the function, with `restart="after"` (which scores identically
   on all 316 corpus targets) reachable from `engine.py` for anyone re-testing it. Worth
   twenty more patterns before it is believed.
3. **44 constructed filesets** (§8.2) where the model's end is short — and this residual
   ships. Out of sample the committed rule reproduces the reference on 93 % of random word
   canvases, not 100 %: a fresh battery of 200 (seed 11, never used in the fit) scores
   **186 agree / 14 disagree**, our binary matching the model on all 200. Every miss is a
   sequence alternating 20- and 64-mismatch words with near-empty ones; the misses split
   into "our end is short", "our start is early" and "the reference refuses what we call",
   so it is not one sign of error. A fifth constant may be hiding there, or the bridging
   rule may interact with the chunk scan in a way §6 does not capture. Nothing in the
   corpus resembles those sequences, which is exactly why the corpus is exact and this is
   not — and why that must not be read as the rule being finished.

## 11. Reproducing everything here

```bash
cd docs/research/fixtures
python3 segfit.py 0        # the word grid is not disturbed by monomorphic markers
python3 segfit.py 1        # §3, the reported interval of a uniform block
python3 segfit.py 2        # §4, the staircase against block width
python3 segfit.py 3        # the (m, h) plane
python3 segfit.py 4        # §5, the trailing trim and the 19/24/32/48 thresholds
python3 segfit.py 5        # §5, the start trim
python3 segfit.py 6        # §8.2, the model against the reference on random sequences

cd ../../../tests/parity/fit
python3 chunk.py           # §8.1, the corpus scorecard (no reference binary needed)
python3 chunk.py -v        # ...with every miss
python3 segtry.py          # §9, the .seg port and why it was not committed
```

`segfit.py` and `segtry.py` drive a binary — the KING 2.3.2 reference by default, or
whatever `$KING` names, which is how §8.4 points the same rig at our own build:

```bash
cd docs/research/fixtures
KING=../../../target/release/king python3 segfit.py 6    # 60/60 — our binary is the model

cd ../../../tests/parity/fit
python3 check_mirror.py    # `engine.py` is the committed Rust engine: MIRROR OK
```

`chunk.py` reads the captured corpus only and needs no binary at all; its `committed` row
is now this rule and its `pre-chunk` row is the one it replaced. `segfit_measured.json` is
658 cached reference answers, so the two grid searches of §7 re-run in a second.
