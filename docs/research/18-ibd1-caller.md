# The `.seg` IBD1 caller, measured on its own canvas

**Status: measured, committed to Rust, validated out of sample. One clause measured and
deliberately left out — §9.**

`docs/research/17-seg-caller.md` built a canvas whose read-back column is `IBD2Seg` and
closed the IBD2 caller with it; its §13 then said the residual had moved to the *union*
`IBD1Seg + IBD2Seg`, which is the IBD1 side. This is that campaign. The instrument is
`docs/research/fixtures/ibd1canvas.py`, the same canvas built the other way up.

No KING source was read. Every rule below is a reading taken off the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) run on
filesets built for the purpose, or a score against the captured parity corpus.

**Headline, stated first.**

| claim | evidence |
| --- | --- |
| The canvas works, and the isolation holds: every fixture below reports **`IBD2Seg 0.0000`**, so the printed `IBD1Seg` is the chromosome-2 IBD1 call and nothing else. A block of `W` callable words walled by all-IBS0 words reads back exactly `64(W+1) − 1` marker intervals. | §1 |
| **The geometry was already right.** The word predicate, both endpoints, the gate, the absence of a push and the absence of bridging are all exactly what `Scan::ibd1` already did — but they were inferred from `--ibs` inversions and early fixtures, and they are now **bisections**. | §2–§5 |
| **What was wrong is one line outside the caller**: how an IBD1 call and the IBD2 calls inside it combine into `IBD1Seg`. The cut is at **marker** granularity and **excludes the IBD2 call's own end markers**, and **each surviving piece faces the `--seglength` floor on its own**. The floor runs *before* the subtraction: a dropped IBD2 call is not subtracted at all. | §6 |
| On the corpus, at the default 3 Mb floor, that takes **`IBD1Seg` from 826 exact rows to all 982**, exact-on-all-four-columns 709 → **747**, mean `PropIBD` error 0.00037 → **0.000067**, worst row 0.0089 → **0.0042**, with 0 extra / 0 missing either way. **Five parity cases flipped** — the first time any single `.seg` rule change moved the harness. | §8.1 |
| Held out: **240 of 240** random canvases from three unused seeds, **120 of 120** on the exhaustive length-≤4 battery, **48 of 48** canvas × `--seglength` combinations of the overlap family, and on mixed canvases *every* miss is one where the `IBD2Seg` column is also wrong. | §8.2 |
| **One clause measured and not implemented.** Above the default floor the reference *does* merge two IBD1 runs across a short interruption — bisected twice, on the IBS0 count (5 merges, 6 does not) and on the gap (to the base pair, at two spacings). It cannot fire at 3 Mb, it is the whole `IBD1Seg` residual at 5 and 10 Mb, and the obvious generalisation of it makes those floors much worse. | §9 |

---

## 1. The instrument

`docs/research/fixtures/ibd1canvas.py`, which reuses `segcanvas.Canvas` unchanged.

* **chromosome 1 — the carrier.** The same five `inf1`-carrying, IBD2-dirty words as
  `17-…`. It gives the pair one 10.527 Mb IBD1 segment, which is what earns it a `.seg`
  row, and `segcanvas.mk(cv, res, 1)` subtracts its fixed contribution.
* **chromosome 2 — the canvas.** Complete words painted marker by marker, walled at both
  ends by all-IBS0 words.
* **the isolation.** `17-…` §3 established that the IBD2 pass refuses any word carrying two
  or more het-vs-hom mismatches. Every word painted here carries **thirty-four**, so the
  IBD2 pass has nothing to call and `IBD2Seg` reads `0.0000` — printed and checked on every
  fixture. That is the whole trick, and it is the mirror image of `17-…`, where IBS0 was
  the paint and the callable words had to be mismatch-free.
* **`inf1`.** IBD1 informativeness is `p1_i & p1_j & (p0_i | p0_j)` — both carry A1 and at
  least one is homozygous — so it comes from A1A1/A1A1 markers and from het-vs-A1A1 ones.
  The latter are also het-vs-hom mismatches, which is convenient: one marker kind supplies
  `inf1` *and* keeps the IBD2 pass out.
* **the ruler.** Unchanged: chromosome 2's uniform spacing puts `D` just over the 100 Mb
  floor of `17-…` §2, so one ulp of the printed `%.4lf` is about a ninth of a marker gap
  and the column reads back the number of marker intervals called.
* **the decode.** IBD1 calls inside one usable segment come out adjacent — each starts one
  marker past the previous call's end (§2) — so `c` calls covering `w` whole words measure
  `64w − c`, and `c = (−M) mod 64` recovers both counts, exactly as on the IBD2 canvas.

Sanity, before anything else:

```
    all wall                       0.01 mk   ibd2 0.0000
    K x 1                        126.95 mk = 2 words / 1 calls    ibd2 0.0000
    K x 2                        191.05 mk = 3 words / 1 calls    ibd2 0.0000
    K x 3                        255.04 mk = 4 words / 1 calls    ibd2 0.0000
    K x 4                        319.03 mk = 5 words / 1 calls    ibd2 0.0000
    K x 8                        574.98 mk = 9 words / 1 calls    ibd2 0.0000
```

`64(W+1) − 1`, not `64W − 1`: the call swallows the whole trailing wall word, because a
wall's last opposite homozygote sits at bit 63 and that is where the right end goes (§2).
The one-word row, 127 markers, is the same number `13-informativeness-gate.md` §6 reported
from a different rig.

**The alphabet.** `K` = 34 mismatches + 12 A1A1/A1A1 (IBD1-callable, IBD2-dead, `inf1` 12);
`k` = the same with no `inf1`; `W` = 64 opposite homozygotes; `Z(bits)` = `K` with opposite
homozygotes forced at named bits; `B` = 12 A1A1/A1A1 + 52 HetHet (**the only word both
passes can use**); `b` = 64 HetHet (IBD2 only — no `inf1`).

## 2. The word predicate and the endpoints

### 2.1 One opposite homozygote, no tolerance, and nothing else

Eight callable words with `j` consecutive words replaced by one carrying `z` opposite
homozygotes, read back as (words called, calls):

```
 j\z    0      1      2      3      5     64
  1    9/1    9/2    9/2    9/2    9/2    9/2
  2    9/1    8/2    8/2    8/2    8/2    8/2
  3    9/1    7/2    7/2    7/2    7/2    7/2
```

**One IBS0 breaks the word; the split is complete at `z = 1` and nothing about it changes
up to 64.** The same block with words replaced by other content, at `j = 1` and `j = 2`:

```
    ibs1=64            9/1   9/1        hethet=62 + 2 mis   9/1   9/1
    ibs1b=64           9/1   9/1        hom1=62 + 2 mis     9/1   9/1
    miss=64            9/1   9/1        zero=62 + 2 mis     9/1   9/1
    missA=64           9/1   9/1        ibs1=63 + ibs0=1    9/2   8/2
```

64 het-vs-hom mismatches, 64 missing calls, 64 HetHet, 64 A1A1/A1A1, 64 het-vs-A1A1: none
of them breaks an IBD1 run. Add **one** opposite homozygote to the same word and it does.
So the IBD1 word predicate is `ibs0 == 0`, and it has no other term.

### 2.2 Both ends read the *last* opposite homozygote of the immediately flanking word

Six callable words with one boundary word beside them, whose opposite homozygotes sit at
named bits. Six words walled measure 447 markers; the run's own word-aligned span is 383.

| IBS0 bits | left of the block | right of the block |
| --- | ---: | ---: |
| `{0}` | 510.0 (+63) | 384.1 (+1) |
| `{1}` | 509.0 (+62) | 385.1 (+2) |
| `{5}` | 505.0 (+58) | 389.1 (+6) |
| `{31}` | 479.0 (+32) | 414.9 (+32) |
| `{32}` | 478.0 (+31) | 415.9 (+33) |
| `{62}` | 448.0 (+1) | 446.0 (+63) |
| `{63}` | 447.0 (+0) | 447.0 (+64) |
| `{0,1}` | 509.0 (+62) | 385.1 (+2) |
| `{62,63}` | 447.0 (+0) | 447.0 (+64) |
| `{0,63}` | 447.0 (+0) | 447.0 (+64) |
| `{0..19}` | 491.0 (+44) | 402.9 (+20) |
| `{44..63}` | 447.0 (+0) | 447.0 (+64) |
| all 64 | 447.0 (+0) | 447.0 (+64) |

```
left  = 64(a-1) + (last IBS0 bit of word a-1) + 1     -> extension 63 - lastbit
right = 64(b+1) + (last IBS0 bit of word b+1)         -> extension  1 + lastbit
```

Both ends read the **last** opposite homozygote — `{0,63}` reads exactly like `{63}` and
not at all like `{0}` — and both read the **immediately flanking word and no further**:

```
    R: K6 Z{5}           389.06 mk          L: Z{10} Z{5} K6     504.99 mk
    R: K6 Z{5} Z{10}     389.06 mk          L: WALL Z{5} K6      504.99 mk
    R: K6 Z{5} WALL      389.06 mk          L: Z{5} WALL K6      446.95 mk
    R: K6 WALL Z{5}      446.95 mk
```

An opposite homozygote two words out moves nothing. This is a **different geometry from
the IBD2 pass**, which reaches 63 markers past the nearest *het-vs-hom mismatch* and is
blocked whole-word by IBS0 (`17-…` §5). The two passes really do have different endpoint
rules, and each is now bisected on its own canvas.

The asymmetry is only apparent: the call runs from *one past* the last IBS0 before it to
*exactly* the last IBS0 after it. That is what makes consecutive calls come out adjacent,
which is what makes `c = (−M) mod 64` decode.

## 3. There is no push

`17-…` §6 measured that every `.seg` IBD2 call after the first in a usable segment starts
one word late. **The IBD1 pass has no such clause**, at any number of calls:

| canvas | reference | model, no push | model, push 1 word |
| --- | ---: | ---: | ---: |
| `K2 W K2` | **382.1** | 382 | 318 |
| `K2 W K2 W K2` | **573.0** | 573 | 445 |
| `K2 W2 K2 W2 K2` | **573.0** | 573 | 445 |
| `K W K4` | **446.0** | 446 | 382 |
| `K4 W K` | **446.0** | 446 | 382 |
| `k W K4` | **319.0** | 319 | 319 |

Three calls in a row and the third is not pushed either. (The `k` row is the control: with
the first run's `inf1` at zero it is refused by the gate, so only one call is emitted.)

## 4. The gate: `inf1 ≥ 10` over the run's own complete words

```
    1 word, 0 A1A1/A1A1          refused        2 words, 4 each        refused
    1 word, 8 A1A1/A1A1          refused        2 words, 5 each        called (3 words)
    1 word, 9 A1A1/A1A1          refused        1 word, 5+5 mixed      called (2 words)
    1 word, 10 A1A1/A1A1         called         1 word, 62 HetHet      refused
    1 word, 9 het-vs-A1A1        refused        9 + rich right flank   refused
    1 word, 10 het-vs-A1A1       called         9 + rich left flank    refused
                                                10 + rich right flank  called

    4 A1A1/A1A1 + 5 het-vs-A1A1  refused        5 A1A1/A1A1 + 4 het-vs-A1A1  refused
    5 A1A1/A1A1 + 5 het-vs-A1A1  called
```

Four independent bisections at 9/10 — A1A1/A1A1, het-vs-A1A1, the two mixed in one word
(9 either way round refused, 10 called), and five each across two words. HetHet is worth nothing to this pass, which is the
`inf1`/`inf2` distinction `13-informativeness-gate.md` §5 established from the other side.

The last three rows are the clause that matters and had not been re-tested on a canvas:
**the markers the call reaches into do not count.** A run carrying 9 is still refused when
the flanking word its call swallows whole carries 40. So the window is the run's own
complete words, exactly as `Scan::informative` has it — and unlike the IBD2 gate, which
counts from the run's first mismatch-free word *through* the words the right end reaches
into (`17-…` §4).

## 5. Bridging: never

A lone bad word between two long runs, at `--seglength 1`:

| interruption | reference | no bridge | bridge |
| --- | ---: | ---: | ---: |
| 1 IBS0 at bit 0 | **510.0** | 510 | 511 |
| 1 IBS0 at bit 31 | **510.0** | 510 | 511 |
| 1 IBS0 at bit 63 | **510.0** | 510 | 511 |
| all 64 IBS0 | **510.0** | 510 | 511 |
| no bad word | **511.0** | 511 | 511 |

Two calls, always, wherever the IBS0 sits and however few there are. The `j = 1` row of
§2.1's table says the same thing at six different IBS0 counts. **On the IBD1 pass a lone
unusable word is never absorbed** — which is the opposite of the IBD2 pass, where it is
absorbed conditionally (`17-…` §3, §7).

That is the answer at the default floor. It is not the whole answer: §9.

## 6. The overlap — what `IBD1Seg` subtracts, and in what order

This is where the residual was. `IBD1Seg` is the part of an IBD1 call that is *not*
already IBD2, and the canvas can put a real IBD2 call inside a real IBD1 call: `B` (no
IBS0, no mismatch, `inf1` 12) is usable to both passes, `K` bounds it and is usable only to
IBD1, so the IBD1 call runs out to the walls while the IBD2 call spills one word out of the
`B` block.

### 6.1 The cut is at marker granularity and excludes the IBD2 call's end markers

```
      K B3 K       ref ibd1    63.0 mk   exclusive   63.0   inclusive   64.0
      K B6 K       ref ibd1    63.0 mk   exclusive   63.0   inclusive   64.0
      K4 B4 K      ref ibd1   224.0 mk   exclusive  224.0   inclusive  226.0
      K B4 K B4 K  ref ibd1    63.0 mk   exclusive   63.0   inclusive   64.0
      k b4 k       ref ibd1     0.0 mk   exclusive    0.0   inclusive    0.0
      K b4 K       ref ibd1    63.0 mk   exclusive   63.0   inclusive   64.0
```

"inclusive" is the naive `length − overlap` the engine computed before this campaign;
"exclusive" removes the IBD2 call's own markers, so an IBD1 call `[lo, hi]` cut by an IBD2
call `[a, b]` leaves `[lo, a−1]` and `[b+1, hi]`. `K4 B4 K` separates them twice over — the
model's two pieces are 162 and 64 markers under the naive rule and 161 and 63 under the
measured one, and the reference reports 224.

The canvas ruler is uniform, so "one marker short" alone does not say *which* marker. A
**graded ruler** (per-word gaps `70 000 + 7 000·k`) does. On `K4 B4 K` the printed
numerator is 23 532 958 bp, against candidates:

| pieces | bp | |
| --- | ---: | --- |
| `[0,161] + [576,639]` (exclusive) | 23 527 000 | **match**, 0.4 ulp |
| `[0,162] + [575,639]` (inclusive) | 23 723 000 | 14.1 ulp out |
| `[0,161] + [575,639]` | 23 604 000 | 5.3 ulp out |
| `[0,162] + [576,639]` | 23 646 000 | 8.4 ulp out |
| `[1,162] + [576,640]` (shifted) | 23 646 000 | 8.4 ulp out |

`k b4 k` is the control for the *shape* of the subtraction: `k` carries no `inf1`, so the
IBD1 gate refuses the run while the IBD2 pass calls the `b` block anyway. `IBD1Seg` reads
`0.0000`, not a negative number and not a reduced carrier — so the subtraction is **per
IBD1 call**, never an aggregate over the pair or over the usable segment.

### 6.2 Every piece faces the `--seglength` floor on its own

`K B3 K` leaves a single piece of exactly 4 410 000 bp inside an IBD1 call of 26.8 Mb:

```
      --seglength 4.409000 Mb -> piece counted        --seglength 4.410001 Mb -> dropped
      --seglength 4.409999 Mb -> piece counted        --seglength 4.411000 Mb -> dropped
      --seglength 4.410000 Mb -> piece counted
```

Bisected to the base pair, and the comparison is `>=`. The floor is applied to the
**piece**, not to the call it came from: the call is five times the floor and `IBD1Seg`
still reads zero. `K4 B4 K` shows the same thing with two pieces of different sizes — at
`--seglength 5` the 14.2 Mb piece is counted and the 5.5 Mb one is not (161 markers, not
224).

### 6.3 The floor runs before the subtraction

A canvas with a short IBD2 call — `[K B K]` at 35 000 bp spacing, IBD1 call 8.921 Mb, IBD2
call 6.686 Mb, piece 2.205 Mb — swept over `--seglength`:

| `--seglength` | ref `IBD1` | ref `IBD2` | model |
| --- | ---: | ---: | --- |
| 1, 2.2 | 62.99 mk | 191.03 mk | 63 / 191 |
| 2.3 … 6.6 | 0 | 191.03 mk | 0 / 191 |
| **6.7, 8** | **254.88 mk** | **0** | **255 / 0** |
| 9, 10 | 0 | 0 | 0 / 0 |

At 6.7 Mb the IBD2 call falls under the floor, disappears from `IBD2Seg` — **and the whole
IBD1 call comes back**. So the calls subtracted are the ones that survived `--seglength`,
not every call that cleared the gate. Scoring the other way round on the corpus confirms
it: subtracting the unfiltered IBD2 calls gives `IBD1Seg` 949 of 982 against 982.

The corpus settles the ordering against the IBD2 **push** the same way: the calls
subtracted are the *pushed, clipped, filtered* ones — the ones printed in `IBD2Seg`.
Substituting the un-pushed IBD2 calls costs 141 rows (`IBD1Seg` 841 of 982) and multiplies
the mean `PropIBD` error by 27 (0.00178 against 0.000067).

### 6.4 The rule

```
for each IBD1 call [lo, hi], in order:
    pieces = [lo, hi] with the closed marker ranges of the surviving IBD2 calls removed
    IBD1Seg numerator += sum of pos[p.hi] - pos[p.lo] over pieces with that >= seglength
```

`ibd1_pieces` in `crates/open-king-core/src/ibdseg.rs` is this; `pieces()` in
`docs/research/fixtures/ibd1canvas.py` and in `tests/parity/fit/seg18.py` are the same
function, and the three agree everywhere. Over the whole overlap family × eight
`--seglength` floors the model reproduces the reference on **48 of 48** combinations.

### 6.5 A cross-check the campaign fell over, already implemented

Sweeping `--seglength` past 10 Mb makes the floor *stop applying*, non-monotonically. It is
not a rule, it is the reference's own input validation, printed verbatim:

```
KING supports minimum segment length from 1 to 10 Mb at the moment.
Default seglength of 3Mb is used.
```

Accepted at `--seglength 10.001` ("Minimum segment length is set as 10001000 bp"), rejected
at `10.01` — the stored value is rounded to two decimals before the range test. `open-king-cli`
already models this (`console::SEGLENGTH_MIN/MAX`, `analysis::ibdseg::seglength_bp`), which
is why the canvas reproduces it; recorded here because it looks exactly like a caller rule
until you read the console.

## 7. The rule, end to end

```
per pair, per usable segment [w0, w1] covering markers [lo, hi]:

  a word is USABLE iff it carries no opposite homozygote          (no tolerance, §2.1)
  a lone unusable word is NEVER absorbed                          (§5, but see §9)
  runs are the maximal stretches of usable words

  for each run [a, b], in order:
      if inf1 over words a..b < 10:  refuse the run               (§4)
      left  = 64a       if a == w0 or word a-1 carries no IBS0    (cannot happen inside)
            = 64(a-1) + (last IBS0 bit of word a-1) + 1           (§2.2)
            = the segment's own `lo` when a is the segment's first complete word
      right = 64(b+1) + (last IBS0 bit of word b+1)               (§2.2)
            = 64(b+1) + 63 when word b+1 carries no IBS0
            = the segment's own `hi` when b is the segment's last complete word
                (or, in the trailing fringe, one marker before its first IBS0)
      clip to the previous emitted call:  left = max(left, previous right + 1)
      keep it if pos[right] - pos[left] >= seglength

  IBD1Seg numerator = sum over kept IBD1 calls of                 (§6)
      sum over the pieces of that call not covered by a surviving IBD2 call,
      each piece measured over its own markers and kept only if it is >= seglength
```

`predict1()` in `docs/research/fixtures/ibd1canvas.py` is the caller; `subtract()` there and
`pieces()` in `tests/parity/fit/seg18.py` are §6; `Scan::ibd1` and `ibd1_pieces` in
`crates/open-king-core/src/ibdseg.rs` are the port.

## 8. Out of sample

Nothing in this section had any part in choosing a constant.

### 8.1 The corpus — 982 `.seg` rows over ten datasets

```
$ python3 tests/parity/fit/seg18.py
--seglength 3 Mb
  retired (length minus overlap)   exact  709  ibd1  826  ibd2  896  extra 0  miss 0  MAE 0.000365  worst 0.0089
  18 measured                      exact  747  ibd1  982  ibd2  896  extra 0  miss 0  MAE 0.000067  worst 0.0042
--seglength 5 Mb
  retired                          exact  701  ibd1  810  ibd2  880  MAE 0.000661  worst 0.0552
  18 measured                      exact  729  ibd1  909  ibd2  880  MAE 0.000177  worst 0.0598
--seglength 10 Mb
  retired                          exact  668  ibd1  766  ibd2  877  MAE 0.001580  worst 0.0679
  18 measured                      exact  692  ibd1  841  ibd2  877  MAE 0.000399  worst 0.0874
```

**`IBD1Seg` is exact on all 982 rows at the default floor** — and on every `.kin`, `.kin0`,
`X.kin` and `cluster.kin` row in the whole corpus, which took `<prefix>X.kin` to
byte-identical in all 15 cases. `tests/parity/run_parity.py` moved **403 → 408 of 480**: the
five `monomorphic --related*` cases. The knob grid says both clauses are needed and neither
is sufficient:

| | `frag drop` | `frag whole` | `frag keep` |
| --- | ---: | ---: | ---: |
| **cut exclusive** | **982** | 859 | 847 |
| **cut inclusive** | 827 | 827 | 826 |

(`IBD1Seg` exact rows; `frag keep` + `cut inclusive` is the retired rule.)

The split that made this findable is the one `17-…` §13 drew. Before: `IBD1Seg` exact on
**all 823** rows whose reference `IBD2Seg` is 0 and on **3 of 159** where it is not — an
IBD1 caller that is perfect except in the presence of IBD2 is not an IBD1 caller problem.
After: 823 and **159**. And the residual's *sign* was the other clue — before the fix
`IBD1Seg` was too high on 156 rows and too low on none, which is what a missing subtraction
looks like.

### 8.2 Held-out canvases

`ibd1canvas.battery(seed, n)` draws word sequences with IBS0 density, IBS0 placement,
`inf1` content and mismatch load all randomised, drives the reference, and compares.

```
IBD2-free canvases (IBD1Seg graded):
  seed 201     80 / 80      seed 5150   80 / 80      seed 99991  80 / 80     -> 240 / 240
exhaustive length-<=4 over {K, k, W}:                                          120 / 120
the overlap family x 8 --seglength floors:                                      48 /  48
mixed canvases, BOTH columns graded (so §6 is on trial):
  seed 3300  L=1   58/60      seed 61803 L=1   57/60      seed 3300 L=4   58/60
```

On the mixed canvases at `L = 1` and `L = 4`, **every miss is a canvas whose `IBD2Seg` is
also wrong** — 0 of 180 have `IBD2Seg` right and `IBD1Seg` wrong. That is the same
statement as the corpus's 982/982: at the default floor, wherever the IBD2 caller is right,
this rule is right. At `L = 8` it is not (17 of 60), and that is §9.

## 9. The one clause measured and deliberately not implemented

§5 says a lone bad word is never absorbed. That is true at the default 3 Mb and it is not
the whole story. Raise `--seglength` past the gap the bad word leaves and the two runs
**merge into one**. Two bisections:

**The gap**, to the base pair, at two different marker spacings. The interruption is one
word, so the earlier run ends at marker `64a−1` and the later starts at `64(a+1)` — 65
marker intervals apart:

```
  s = 88 000 bp   gap = 5 720 000   --seglength 5.719999 -> split   5.720001 -> MERGED
  s = 47 000 bp   gap = 3 055 000   --seglength 3.054999 -> split   3.055001 -> MERGED
```

so the condition is `pos[first marker of the later run] − pos[last marker of the earlier
run] < --seglength`, strictly, and the threshold is a marker count and not a fixed bp.

**The IBS0 the interruption may carry**, at `--seglength 8`: 4 merges, **5 merges, 6 does
not**, 7, 8, 16 and 64 do not — and it does not matter where the bits sit (five at bits
`0..4` and five at `0,16,32,48,63` both merge; six scattered do not). It is the total over
the interrupting words, which may be more than one: two words at one IBS0 each merge when
the resulting 129-interval gap is under the floor (`L = 8`, `s = 47 000`) and split when it
is not (`L = 6`); three words at one each, 193 intervals, split at `L = 8`.

**It cannot fire at the default floor**, because 65 marker intervals at any real spacing is
over 3 Mb — which is exactly why `IBD1Seg` is exact on all 982 corpus rows at 3 Mb, on 909
at 5 Mb and on 841 at 10 Mb, and why the merge shows up on the canvas only above ~5.7 Mb
(its spacing is 88 000 bp).

**And the obvious generalisation is wrong.** Implemented as "merge two runs when the words
between them carry at most 5 opposite homozygotes and the gap is under `--seglength`", with
the merged run then taking the gate, the endpoints and the floor:

```
  --seglength  5 Mb   IBD1Seg 795 of 982   (against 909 without the clause)
  --seglength 10 Mb   IBD1Seg 732 of 982   (against 841)
```

and, if the merged calls are also allowed to feed the ">10 Mb" pair filter, **1 054 extra
pairs**. Real unrelated pairs are full of one-IBS0 interruptions, so the reference must
require something further that the canvas has not yet asked for — a bound on how many
merges, or on the runs' own lengths, or on their `inf1`. Until that is measured the clause
stays out of `crates/open-king-core`, where on net it would lose 114 `IBD1Seg` rows at 5 Mb and
109 at 10 Mb.

## 10. Reproducing everything here

```bash
cd docs/research/fixtures
python3 ibd1canvas.py 0     # §1 — the rig, the isolation check, the marker ruler
python3 ibd1canvas.py 1     # §2.1 — the word predicate
python3 ibd1canvas.py 2     # §2.2 — the endpoints and how far out they read
python3 ibd1canvas.py 3     # §3 — no push
python3 ibd1canvas.py 4     # §4 — the gate
python3 ibd1canvas.py 5     # §5 — no bridging
python3 ibd1canvas.py 6     # §6 — the overlap, its two bisections, and the ordering
python3 ibd1canvas.py 7     # §8.2 — the exhaustive length-<=4 battery
python3 ibd1canvas.py 8     # §8.2 — the random batteries
python3 ibd1canvas.py 9     # §9 — the --seglength run merge

cd ../../../tests/parity/fit
python3 seg18.py            # §8.1 — the corpus scorecard at 3 / 5 / 10 Mb
python3 seg18.py grid       # §8.1 — the knob grid
```

`ibd1canvas.py` drives the KING 2.3.2 reference by default, or whatever `$KING` names, and
caches every answer in `ibd1canvas_measured.json` — **1 013 invocations**, so the whole
document re-runs in 0.3 s without the binary. As with `segcanvas.py`: do not point `$KING`
at a non-reference build in this directory, because it writes whatever it measures into
that cache. Grade a copy.

## 11. What is still open

1. **§9's merge.** Bisected on both of its constants and not modelled. It is the entire
   `IBD1Seg` residual at `--seglength 5` and `10` and cannot appear at 3 Mb.
2. **The fringes.** A call touching the usable segment's first or last complete word runs
   to the segment's own first or last marker, and in the trailing fringe stops one marker
   before the first IBS0 there. The canvas cannot see this — its chromosome 2 is
   word-aligned — so it rests on the corpus, where it is exact on all 982 rows. A canvas
   with a chromosome whose marker count is not a multiple of 64 would bisect it.
3. **Whether the ">10 Mb" pair filter reads the IBD1 call or its pieces.** The corpus
   cannot tell (0 extra, 0 missing either way) and the canvas as built cannot either,
   because its carrier always supplies a 10.5 Mb segment. A canvas with a sub-10 Mb carrier
   would.
4. **The residual is now entirely `Scan::ibd2`** — `IBD2Seg` exact on 896 of 982, 86 rows
   two-sided by at most 0.0042. `17-…` §11 is the list; nothing in it moved.
