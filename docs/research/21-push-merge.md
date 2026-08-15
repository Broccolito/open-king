# The push, and the IBD2 merge re-measured

**Status: measured, validated out of sample, committed to Rust.** This closes three of the
five remaining parity cases (**472 → 475 of 480**) and takes the `--seglength 5` corpus
from 959/947 to **982/982** on both estimate columns — the same score the default floor
has. `--seglength 10` goes 960/945 → **970/972**, and its residual is now one-sided the
other way (§6), which is the sharpened negative this write-up hands on.

`docs/research/20-seglength-floor.md` §11 named the first thing to look at: the IBD2 merge
scored 56-58 of 60 on random IBD2-native canvases while the same mirror scored 60/60 where
no merge can fire, and the failures "smell like an interaction with the one-word push
(`17-…` §6) rather than with the merge". Half right. The push **is** wrong — clause 1
below — but so are three separate clauses of the IBD2 merge, and the reason `20-…` did not
see them is that every one of its IBD2 fixtures happened to hold the confounding variable
constant.

No KING source was read. Every rule below is a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) run on
filesets built for the purpose, or a score against the captured parity corpus.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| **The one-word push is conditional.** `17-…` §6 has *every* call after the first in a usable segment starting one word later. A call arms the push only when it reaches **half** the floor — and the length that counts is measured from its own **gate-start word**, not from its left end. `armed \|= pos[hi] - pos[64*gs] >= seglength/2`, sticky for the rest of the segment. Bisected to the base pair on five spacings; the integer division is the reference's own. | §2 |
| **The IBD2 merge has no word cap.** `20-…` §3 bisected "at most two unusable words, absolute" on the **IBD1** pass and the IBD2 pass simply inherited it. Fifteen unusable words merge on the IBD2 pass when the gap and the budget allow. The same fixture on the IBD1 pass still refuses three, so the two passes really do differ. | §4 |
| **The interruption runs between the two runs' *gate windows*, not between the runs.** A run's window is `[gs .. ge_of(b)]` — from its gate-start word through the one word its right end reaches into. The gap is measured window-to-window, and a word a window covers is not part of the interruption. The exclusion holds after *any* usable word, so a gate-refused run's own reach word is skipped too. | §3 |
| **`X` is the HetHet count, not `inf2`.** `20-…` §7 read the IBD2 budget's informative count as `inf2` = HetHet + A1A1/A1A1. It is HetHet, with the IBD1 pass's own switch: `X = HetHet if HetHet >= 10 else A1A1/A1A1`, the `10` bisected at 9/10 on this pass as well. Every fixture in `20-…` §7 used a HetHet filler, which is exactly the case the two readings agree on. | §5 |
| Corpus: **3 Mb does not move** (982 / 982, MAE 0.000023). At 5 Mb `IBD1Seg` 959 → **982**, `IBD2Seg` 947 → **982**, MAE 0.000086 → **0.000023**, worst row 0.0111 → **0.0001** — 5 Mb is now as exact as the default floor. At 10 Mb 960 → **970** and 945 → **972**, MAE 0.000134 → **0.000067**. 0 extra / 0 missing everywhere. | §6 |
| Held out: **357 of 360** random IBD2 canvases on three unused seeds at 5 and 10 Mb, against **343 of 360** for the committed rule. | §7 |
| Harness: **472 → 475 of 480**. Closed: `ibdseg/{bigish,missing}__ibdseg_seglength5` and `ibdseg/missing__ibdseg_seglength10`. Self-check stays 480/480; 123 core tests, clippy and fmt clean; `check_mirror.py` green on all 2 946 corpus rows. | §6 |

---

## 0. The instrument

`docs/research/fixtures/push1.py`, which reuses `segcanvas.Canvas` and `mergelab`'s answer
cache. Everything in it is built out of two words —

```
CLEAN   64 HetHet         usable, inf2 = 64, no mismatch
WALL    64 opposite homs  unusable, unbridgeable, never absorbed by anything
```

— so a canvas is a sentence in runs and walls and the printed `IBD2Seg` decodes to the
number of marker intervals called. A word-aligned call over `k` words measures `64k - 1`
markers; the same call pushed one word measures `64k - 65`. That 64-marker difference is
the whole readout for §2.

The second instrument is `push1.steps(canvas, lo, hi, step)`: sweep `--seglength` and
record every value at which the printed column changes. The jumps of that step function
**are** the individual call lengths, so one sweep reads a canvas's entire behaviour
instead of one point of it. Every section below is a step function.

**One hard boundary on the rig.** Above roughly 10 Mb the reference stops behaving like a
floor: a canvas whose only call is 14.06 Mb long reports it at `--seglength 10.0` and
reports something else entirely at 10.5. Below 1.0 Mb it behaves as though the flag were
absent. Both regimes are outside the corpus's three floors and no rule here is measured in
them; every sweep is run over `1.0 <= L <= 10.0`.

## 1. Why `20-…` could not see any of this

`20-…` §7 measured the IBD2 merge on three families of canvas, and each one pins the
variable that this campaign had to vary:

* its interruptions were **one and two words wide**, so a cap of two never bound;
* its interruption was **made unusable by an opposite homozygote**, so the earlier run
  never had a word to reach into and runs and gate windows coincided;
* its filler was **HetHet**, so `X = HetHet` and `X = inf2` are the same number.

None of that was careless — each choice was made to stop a *different* clause (the `17-…`
§14 bridge) from firing instead. It just means the IBD2 merge was measured on a
one-dimensional slice of a four-dimensional rule, and the other three dimensions were
filled in by symmetry with the IBD1 pass. Three of the four guesses were wrong.

## 2. The push is armed at half the floor, from the gate-start word

### 2.1 The fixture

Two runs of `k1` and `k2` clean words with one wall between them. The wall blocks the
endpoint reach in both directions and can never be merged across, so the geometry is
exactly

```
call 1               = 64*k1 - 1 markers
call 2, not pushed   = 64*k2 - 1
call 2, pushed       = 64*k2 - 65
```

and sweeping `--seglength` past `call 1`'s own length asks what happens to the push when
its cause is dropped.

### 2.2 The reading

```
runs [1, 4]     1.00:254  1.26:191  2.54:255  5.10:0
runs [2, 5]     1.00:382  2.54:255  5.10:319  6.38:0
```

Take `[2, 5]`: `c1 = 127` markers (2.54 Mb), `c2` is 319 unpushed and 255 pushed.

* `L < 2.54` — both calls, `c2` pushed: `127 + 255 = 382`.
* `2.54 <= L < 5.10` — `c1` is gone and `c2` is **still pushed**: `255`. So the push does
  survive its cause being dropped, exactly as `17-…` §6 said.
* `5.10 <= L < 6.38` — `c2` is **no longer pushed**: `319`. The push has switched off, at
  a floor of `2 x 2.54` Mb.

That last step is the whole finding. A call arms the push while `2 x len >= seglength`,
and stops arming it above. `[1, 4]` says the same at `2 x 1.26 = 2.54`.

### 2.3 It is sticky, and it is not a chain

With three runs the flag is set by *any* earlier call, not only the immediately preceding
one:

```
runs [6, 2, 6]   1.00:765  1.26:702  6.38:383   7.66:0
runs [6, 1, 6]   1.00:702  6.38:383  7.66:0
```

In `[6, 1, 6]` the middle call is empty once pushed and can arm nothing, yet the third
call is pushed all the way to 6.38 — armed by the first, whose `2 x 7.66` covers the whole
range. A per-step chain predicts an unpushed third call there and is wrong by 64 markers.
Because the flag is sticky, the question "does a *pushed* call arm the next one on its
pushed or its unpushed length?" cannot arise: a call is pushed only when the flag is
already set, and once set it never clears.

### 2.4 The length is measured from the gate-start word

Two canvases with the *same* call length disagree, which is what forces the last clause:

```
[MISW, CLEAN, WALL, CLEAN x6]    1.00:446  2.54:383            un-push at 2.54
[CLEAN x2,   WALL, CLEAN x6]     1.00:446  2.54:319  5.10:383  un-push at 5.10
```

Both first calls measure 127 markers (2.54 Mb): the first because a mismatch-only word
before the run pulls its left end back a whole word, the second because the run is two
words wide. The thresholds are `2 x 1.26` and `2 x 2.54` — that is, the left extension does
**not** count. Pushing the gate-start word further in confirms which point the length is
measured from:

```
run1 = [CLEAN x4]                 gs = word 0   un-push above 2 x 255 markers
run1 = [U1, CLEAN x3]             gs = word 1   un-push at 7.66 = 2 x 191 markers
run1 = [U1, U1, CLEAN x2]         gs = word 2   un-push at 5.10 = 2 x 127 markers
```

(`U1` carries one mismatch: usable, but not a gate-start word.) So

    armed |= pos[hi] - pos[64 * gs] >= seglength / 2

### 2.5 The comparison, to the base pair

Bisected on `runs [1, 2, 6]`, whose second call measures 2 540 000 bp:

```
--seglength 5.080000  pushed      5.080001  pushed      5.080100  NOT pushed
```

`2 x 2 540 000 = 5 080 000` is not `>= 5 080 001`, so the test is not `2*len >= L`. It is
`len >= L / 2` with **integer division**: `5 080 001 / 2 = 2 540 000` still passes and
`5 080 100 / 2 = 2 540 050` does not. Reproduced at 7.640000 / 7.640100 on a canvas with a
3.82 Mb call, and the factor of two is confirmed at five spacings (15/20/25/30/35/50 kb),
each time landing on the first grid point above `2 x len`.

### 2.6 Why `17-…` §6 read it as unconditional

§6's control was one canvas at one floor: a 5.54 Mb first call at `--seglength 6`. Half of
6 is 3, and 5.54 >= 3, so the flag is armed and the second call is short by one word —
which is what §6 saw and correctly recorded. The clause only becomes visible when the
first call is *between* half the floor and the floor, and no fixture before this one put a
call in that window.

## 3. The interruption is between gate windows, not between runs

A run's gate window is `[gs .. ge_of(b)]`: from its first mismatch-free word through the
one word its right end reaches into (`ge_of(b) = b+1` when word `b+1` carries no opposite
homozygote, which at a run boundary means it carries mismatches). Both ends of the merge
test read that window rather than the run.

**The earlier end.** Same interruption, two orders:

```
4 CLEAN, [MISW, Z1], 4 CLEAN      merges at 1.31       gap 1.30 Mb = from MISW's end
4 CLEAN, [Z1, MISW], 4 CLEAN      merges at 2.59       gap 2.58 Mb = from the run's end
```

`MISW` is mismatch-only, so in the first canvas the earlier run reaches into it and the
gap is measured from *its* end, one word further along; `Z1` carries an opposite
homozygote, so in the second there is nothing to reach into. The threshold moves by
exactly one word. Varying `MISW`'s mismatch count from 2 to 10 and their positions from
bits 0-3 to bits 60-63 does not move it at all.

**The later end.** Give the later run some leading words that are usable but not
mismatch-free, so its gate-start is not its first word:

```
run2 = [CLEAN x4]              merges at 1.31    gap 1.30 Mb   (gs = first word)
run2 = [U1, CLEAN x3]          merges at 2.59    gap 2.58 Mb   (gs = second word)
run2 = [U1, U1, CLEAN x2]      merges at 3.87    gap 3.86 Mb   (gs = third word)
```

One word of threshold per word of `gs`. So

    gap = pos[64 * gs2] - pos[64 * (ge_of(b1) + 1) - 1]   <  seglength   (strict)

**And the same words are dropped from the budget.** With a two-word interruption whose
first word is mismatch-only, that word is free: sweeping its mismatch load from 2 to 64
never blocks the merge, while the *second* word's load has a threshold at 17, which is the
single-word budget `3*(m - 2) <= 64 - m` exactly. The exclusion is stated on the word, not
on the pair: an unusable, IBS0-free word immediately after **any** usable word is covered
by that word's window. That "any" is load-bearing — a run the gate refused has a window
too, and the corpus rows `bigish 106/108` and `bigish 120/122` at 10 Mb turn on it.

## 4. The IBD2 merge has no word cap

`20-…` §3's "three never merges, at any floor" was measured on the **IBD1** pass. Run the
same sweep on both passes, four clean words either side and `j` interrupting words the
budget allows, at a floor above the widest gap:

```
  j     gap        IBD2 pass          IBD1 pass
  1   1.30 Mb      MERGED  575        MERGED  639
  2   2.58 Mb      MERGED  639        MERGED  703
  3   3.86 Mb      MERGED  703        split     0
  4   5.14 Mb      MERGED  767        split     0
  5   6.42 Mb      MERGED  831        split     0
  6   7.70 Mb      MERGED  895        split     0
  7   8.98 Mb      MERGED  959        split     0
```

At 10 kb spacing, where fifteen words still fit under a 10 Mb gap, the IBD2 pass merges
all fifteen. There is no cap on that pass and the IBD1 cap of two is reproduced on a fresh
fixture in the same run, so this is a real difference between the passes and not a
retraction of `20-…` §3.

The budget itself is unchanged and is still summed over the whole interruption, not tested
per word: six words carrying two opposite homozygotes and five HetHet each merge at
`3*(12-2) = 30 <= 30` and split at 24, and a per-word test would merge both.

## 5. `X` is the HetHet count, with the switch at 10

One interruption word carrying 10 opposite homozygotes, and the 54 remaining markers split
between HetHet (`h`) and A1A1/A1A1 (`u`). `3 * (10 - 2) = 24`:

```
   h    u   h+u    reference   X=HetHet   X=inf2
   0   54    54    MERGED      M          M
   5   49    54    MERGED      M          M
   9   45    54    MERGED      M          M
  10   44    54    split       .          M
  12   12    24    split       .          M
  20   34    54    split       .          M
  23   31    54    split       .          M
  24   30    54    MERGED      M          M
  24    0    24    MERGED      M          M
  23    0    23    split       .          .
  30   20    50    MERGED      M          M
```

`X = inf2` gets six of these wrong. The non-convexity — merge at `h = 9`, split at
`h = 10..23`, merge again at `h >= 24` — is `20-…` §5's clause on the other pass:

    X = HetHet          if HetHet >= MIN_INFORMATIVE
      = A1A1/A1A1       otherwise

bisected here at 9/10 against a constant `h + u = 54`, cleanly, one row per value. The
symmetry with the IBD1 pass is exact once each pass's *primary* informative kind is named:
IBD1 prefers its het-vs-A1A1 markers over its A1A1/A1A1 ones, IBD2 prefers HetHet over
A1A1/A1A1, and both switch at the informativeness gate's own constant.

Confirmed off the corpus's own words as well. `bigish 40/44`'s interruption at 5 Mb holds
20 mismatches against 24 + 30 `inf2` and 4 + 10 HetHet; swapping in synthetic words with
the same `inf2` but all-HetHet content makes the reference merge, and raising the real
words' HetHet content to 44 and 50 finds the boundary at exactly `3*(20-2) = 54`.

## 6. The corpus, before and after

`python3 tests/parity/fit/seg21.py`:

```
--seglength 3 Mb
  20 (committed)        exact 806  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
  21 (push + merge)     exact 806  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
--seglength 5 Mb
  20 (committed)        exact 795  ibd1 959  ibd2 947  extra 0  miss 0  MAE 0.000086  worst 0.0111
  21 (push + merge)     exact 817  ibd1 982  ibd2 982  extra 0  miss 0  MAE 0.000023  worst 0.0001
--seglength 10 Mb
  20 (committed)        exact 793  ibd1 960  ibd2 945  extra 0  miss 0  MAE 0.000134  worst 0.0111
  21 (push + merge)     exact 811  ibd1 970  ibd2 972  extra 0  miss 0  MAE 0.000067  worst 0.0081
```

**5 Mb is now exactly as good as the default floor**: 982 of 982 on both estimate columns,
the same MAE, the same worst row. `python3 seg21.py grid` drops each clause in turn
(`ibd1` / `ibd2` / MAE):

```
                      --seglength 5 Mb                --seglength 10 Mb
  21 (all four)       982 / 982 / 0.000023            970 / 972 / 0.000067
  without reach       959 / 947 / 0.000086            963 / 947 / 0.000133
  without hethet      981 / 980 / 0.000026            970 / 970 / 0.000070
  without push_half   982 / 982 / 0.000023            969 / 971 / 0.000068
  without no_cap      982 / 982 / 0.000023            970 / 972 / 0.000067
```

`reach` and `hethet` carry the corpus; `push_half` is worth one row on each column at
10 Mb; and `no_cap` is worth nothing here, because at the corpus's ~50 kb spacing a 10 Mb
gap holds at most three words. Out of sample it is worth four canvases of 240 (§7), and it
is measured directly in §4 — this is a case where the corpus cannot see a clause that the
reference plainly has.

The same thing measured from the shipped binary against the goldens, with `.seg`'s own
`PropIBD` rule rather than the retired `.kin` one `seg20.py` grades with
(`python3 tests/parity/fit/scorecard.py`):

```
  floor   rows  exact   ibd1   ibd2  extra  missing        MAE    worst
    3 Mb    982    982    982    982      0        0   0.000000   0.0000
    5 Mb    982    982    982    982      0        0   0.000000   0.0000
   10 Mb    982    970    970    972      0        0   0.000046   0.0081
```

**Every one of the 982 rows at 5 Mb is now byte-exact**, against 947 before
(`20-…` §9's note gives the same runs as 947 and 943 exact under this rule).

`python3 tests/parity/run_parity.py --impl ./target/release/king`: **472 → 475 of 480.**
Closed: `ibdseg/bigish__ibdseg_seglength5`, `ibdseg/missing__ibdseg_seglength5` and
`ibdseg/missing__ibdseg_seglength10`. Still open: `ibdseg/bigish__ibdseg_seglength10`,
`ibdseg/multifam__ibdseg_seglength10`, and the three non-`.seg` cases that predate this
campaign. Self-check 480/480, 123 core tests, clippy and fmt clean, and `check_mirror.py`
green on all 2 946 corpus rows at the three floors plus the 158 `MaxIBD2` values, so the
Python mirror and the shipped engine are still the same rule.

`--ibs`'s IBD2 columns are untouched — nothing here is in `Scan::ibd2_words` — and they are
still exact on all 21 560 rows.

## 7. Out of sample

`push1.py` chose every constant. The grading is on canvases none of it saw.

```
  seed 20260814   L=5     60 / 60      L=10    59 / 60
  seed 777333     L=5     59 / 60      L=10    59 / 60
  seed 41414141   L=5     60 / 60      L=10    60 / 60
                                             -> 357 / 360
```

against **343 / 360** for the committed `20-…` rule on the same 360 canvases. Three unused
seeds, graded at 5 and 10 Mb specifically, with the merge live on most of them. The grader
is `seg21.predict` — the function the write-up quotes — not the scratch lab.

The IBD1 battery of `20-…` §10 was re-run on two of those seeds as a control: 238 of 240,
unchanged by this campaign (the IBD1 pass is not touched), and the two misses are the
IBD1-side residual `20-…` already had.

## 8. What is still open

1. **`--seglength 10` is now wrong in the other direction.** The residual is 12 rows, and
   it is as one-sided as 5 Mb's was before this campaign — but inverted:

   ```
   === 5 Mb ===   982 rows   ibd1 0             ibd2 0
   === 10 Mb ==   982 rows   ibd1 12 (+0/-12)   ibd2 10 (+10/-0)
        bigish     7 rows    d1 -10.0 .. -17.5 Mb   d2 +10.0 .. +11.5 Mb
        multifam   5 rows    d1  -5.5 .. -19.2 Mb   d2 +11.2 .. +11.6 Mb
   ```

   Every wrong `IBD2Seg` is too **high** and every wrong `IBD1Seg` too **low**, and on ten
   of the twelve rows `d1 = -d2` to the fourth decimal. That is one fault seen twice, the
   mirror image of `20-…` §11.2: **this caller now merges IBD2 where the reference does
   not**, and the IBD2 it invents is then subtracted from `IBD1Seg`. The two rows where
   `d2 = 0` (`bigish 130/136`, `multifam 11/18`) have `IBD2Seg 0.0000` in the reference, so
   there the invented merge is the *only* IBD2 the pair has.

   What it is not: the word cap (capping the IBD2 merge at 2, 3 or 4 words changes nothing
   on this corpus — at the corpus's ~50 kb spacing a 10 Mb gap holds at most three words),
   the budget's constants (they are bisected in §4 and §5 at equality), or the push (which
   is worth one row at 10 Mb and none at 5). The natural next suspect is the *gap*: it is
   the only floor-dependent term in the merge, the residual appears only at the larger
   floor, and each wrong row is worth about one merged segment. A canvas whose gap sits
   between 5 and 10 Mb, swept at both floors, would say whether `gap < seglength` acquires
   a second, absolute bound above 5 Mb.

2. **The reference stops behaving like a floor above ~10 Mb.** A canvas whose single call
   is 14.06 Mb reports it at `--seglength 10.0` and reports a constant 8.93 Mb at 10.5 and
   at every larger floor, on canvases of quite different content. Below `--seglength 1.0`
   the same canvases behave as though the flag were absent (a 2.54 Mb call is dropped at
   0.4 and kept at 1.0). Neither regime touches the corpus, and no rule here is measured in
   them — but any future rig must stay inside `1 <= L <= 10`, and the two boundaries are
   probably the same fact seen twice.

3. **The `.seg` floor test is strictly greater, not `>=`.** A call measuring exactly
   5 100 000 bp is reported at `--seglength 5.099999` and dropped at `5.100000`. The engine
   compares `>=` and scores 982/982 at 3 and 5 Mb, so nothing on the corpus lands on an
   exact tie; it is recorded here because every fixture in this rig does, and a future
   canvas battery that ignores it will misgrade its own boundary rows.

4. **Why half.** `seglength / 2` arming the push, `3 * (bad - 2)` costing a merge, `10`
   switching `X` — the constants are bisected but not explained. That the push's threshold
   and the merge's `X`-switch both reuse a number the caller already has (the floor, and
   `MIN_INFORMATIVE`) is suggestive and no more.

## 9. Reproducing everything here

```bash
cd docs/research/fixtures
python3 push1.py                # §2 §3 §4 §5, all four bisections
python3 push1.py 2              # just the push

cd ../../../tests/parity/fit
python3 seg21.py                # the corpus scorecard at 3 / 5 / 10 Mb, 20 vs 21
python3 seg21.py grid           # each of the four clauses dropped in turn
python3 check_mirror.py         # engine.py still equals the shipped binary
cd ../../..
python3 tests/parity/run_parity.py --impl ./target/release/king
```

`push1.py` drives the KING 2.3.2 reference by default and shares `mergelab_measured.json`.
As with the other rigs: do not point `$KING` at a non-reference build in this directory,
because it writes whatever it measures into that cache.
