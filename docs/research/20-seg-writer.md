# 20 — `<prefix>.seg` is not the `.kin` table with different columns

Two rules, both in the **writer** rather than the caller, both found from the reference's
own captured output with no genotypes involved, and together worth 28 of the 480 parity
cases. Reproduce every number here with:

```bash
python3 docs/research/fixtures/segwriter.py
```

It reads `tests/parity/golden/` and nothing else — no reference binary, no engine, no
fixtures. That is the point of this document.

---

## 0. Why these were invisible until now

`17-seg-caller.md` and `18-ibd1-caller.md` closed `IBD2Seg` and `IBD1Seg`; `19-…` closed
the last of `IBD2Seg` with the fringe clause. That left, at the default 3 Mb floor:

```
982 rows   IBD1Seg exact on 982   IBD2Seg exact on 982   byte-exact rows 806
```

**176 rows had both estimate columns exact and still printed a different `PropIBD`.** No
change to a segment caller can produce that shape: `PropIBD = IBD2Seg + IBD1Seg/2` reads
the same two totals the columns read. Either the third column is computed from something
else, or the residual is a sub-half-ulp difference in the totals that only the third
rounding can see.

`tests/parity/fit/prop19.py` and `where19.py` pursued the second reading for a whole
campaign and got a long way with it — every needed correction is under half an ulp of
`IBD1Seg` (max 0.497, median 0.246), which is exactly why `IBD1Seg` stayed at 982 — and
they were chasing the wrong thing. The first reading is correct, and it is checkable
without a caller at all.

**The instrument for this class of bug is a byte diff of a file whose numbers are already
right.** The canvases (§8 of `docs/MAINTAINING.md`) cannot see a writer: they read one
printed column back as a marker count and never compare two files.

---

## 1. The reference contradicts itself, and that is the finding

Run the reference once:

```bash
king -b bigish.bed --related --degree 2 --ibdseg --cpus 1 --prefix r
```

**147** pairs land in both `r.kin` and `r.seg`. **All 147** carry identical `IBD1Seg` and
`IBD2Seg` in the two files. **43** carry a different `PropIBD`:

| `IBD1Seg` | `IBD2Seg` | `.kin` `PropIBD` | `.seg` `PropIBD` |
| --- | --- | --- | --- |
| 0.4885 | 0.2974 | 0.5417 | **0.5416** |
| 0.3852 | 0.3123 | 0.5048 | **0.5049** |
| 0.5207 | 0.1808 | 0.4411 | **0.4412** |

Same invocation, same pair, same printed inputs, two answers — and the disagreement goes
both ways, so it is not "one of them is rounded up". Over the whole corpus (§1 of the rig):
**6** captures write both files, **201** pairs appear in both with identical estimates, and
the two files disagree on **54** of them — 26.9 %.

So there is no single expression that reproduces both. Each writer needs its own, and
open-king now gives them one each: `Segments::prop_ibd` (full precision) for the `.kin`
family, `ibdseg::seg_prop_ibd` for `.seg`.

---

## 2. `.seg` computes `PropIBD` from its own printed columns

Test the hypothesis in exact integers, so no floating point enters the test. Let `i1`, `i2`
and `m` be the printed `IBD1Seg`, `IBD2Seg` and `PropIBD` scaled by 10 000. The
printed-column combination in units of half an ulp is

```
n = 2*i2 + i1
```

and it is consistent with `m` iff `n == 2m` (the combination lands mid-cell, so any
rounding agrees) or `|n − 2m| == 1` (it lands exactly on a cell boundary — an exact decimal
tie, where the rounding decides). Over **all 4 172** `.seg` rows the corpus captures:

| | rows |
| --- | ---: |
| unambiguous, and correct | 2 859 |
| exact tie, the reference rounded **up** | 1 099 |
| exact tie, the reference rounded **down** | 214 |
| **inconsistent — would refute the rule** | **0** |

Zero refutations. And the ties go both ways in a 5:1 ratio, so the tie-break is arithmetic
and not a convention: any "round half up" or "round half even" rule is wrong on one of
those two groups by construction.

### 2.1 How strong is "0 refutations"?

Not vacuous, and this is the number that decides it. The competing hypothesis is that
`.seg` computes `PropIBD` from the **unrounded** totals, as `.kin` does. Under that
hypothesis the printed value can sit further from the printed-column combination than half
an ulp, and the test above would flag it. How often? Measure it on our own engine's
unrounded values over the 982 primary rows:

```
|unrounded − printed-column combination|,  in printed ulps:
   median 0.090   p90 0.234   max 0.693
   ≥ 0.25 ulp : 76 rows (7.7 %)
   ≥ 0.50 ulp : 17 rows (1.7 %)
```

1.7 % of rows put the unrounded value at least half an ulp away. Over 4 172 reference rows
that predicts **≈ 71** refutations under the unrounded hypothesis. **0** were observed.

### 2.2 It is `.seg`'s rule and no other file's

The same integer test on every other captured file carrying the three columns:

| file | rows | verdict |
| --- | ---: | --- |
| `king.seg` | 4 172 | **consistent, 0 refutations** |
| `king.kin` | 4 248 | **refuted on 42** |
| `kingcluster.kin` | 165 | **refuted on 3** |
| `king.kin0` | 302 | consistent (too few near-tie rows to separate) |
| `kingX.kin` | 90 | consistent (likewise) |

The two files with enough rows to discriminate both refute it. open-king reproduces all
four of them byte for byte using the full-precision value, which is the direct evidence
that they use it.

---

## 3. Which expression, exactly

The 1 313 ties are the discriminator: they are decided by which side of the exact decimal
half the reference's double lands on, so every candidate that is *mathematically* the same
expression gives a different answer. Scored over all 4 172 rows:

| expression on the printed integers | exact |
| --- | ---: |
| **`i2 * 1e-4 + i1 * 5e-5`** | **4 172 / 4 172** |
| `i2 * 0.0001 + i1 * 0.00005` (the same doubles) | 4 172 / 4 172 |
| `(i1 + 2*i2) * 5e-5` | 4 086 |
| integer round-half-up | 3 958 |
| `i2/10000 + i1/20000` | 3 825 |
| the printed values as doubles, `b + a/2` | 3 825 |
| the printed values as doubles, half-up | 3 825 |
| `(i1 + 2*i2)/20000` | 3 804 |
| `(i1/2 + i2)/10000` | 3 804 |

One expression reproduces every row including all 1 313 ties, in both directions. A
1 313-way coin flip does not come out right by luck.

**What this says about the reference.** It holds the two proportions in units of 1e-4 —
integers, or values already quantised to the fourth decimal — and forms `PropIBD` by
scaling them, rather than by halving a full-precision proportion. This document does not
claim to know which; the rule is what is measured, and it is exact.

Committed as `open_king_core::ibdseg::seg_prop_ibd`, used by the `.seg` writer only.
`InfType`, the `--degree` filter, `--unrelated`'s greedy and `--related`'s `Error` grader
all still read the full-precision value: this reaches one column of one file.

---

## 4. The row order: 16-sample blocks

With `PropIBD` fixed, `multifam` and `bigish` still failed — with **every value correct**
and the rows in a different order. All 104 `multifam` rows matched on key, none on
position.

The reference does not list pairs by sample index. It walks **blocks of 16 samples**: for
each block `b1`, for each block `b2 ≥ b1`, every reported pair with `i` in `b1` and `j` in
`b2`, `i` then `j` ascending. `multifam` (20 samples, 2 blocks) shows it plainly — the
first block is finished at `(13, 14)` before the first pair reaching into the second,
`(11, 16)`, is written:

```
 ... (11,12) (11,13) (11,14) (12,13) (12,14) (13,14)   <- block 0 x block 0 ends
     (11,16) (11,17) (11,18) (11,19) (12,16) ...       <- block 0 x block 1 begins
     ... (14,19) (15,17) (15,18) (15,19)
     (16,17) (16,18) (16,19) (17,18) (17,19) (18,19)   <- block 1 x block 1
```

**The block size is 16 and nothing else.** Sweeping it over 2..80 against the row order of
**all 50** captured `.seg` files, exactly one value reproduces every one. `threegen`
(12 samples) rules out everything below 12; `multifam` rules out 20 and up, plain index
order included; only 16 survives `bigish` (200 samples, 13 blocks) inside that window.

It is not a thread artifact: the reference gives the identical order at `--cpus` 1, 2, 4
and 8, and matches the multi-threaded golden in every case.

Committed as `analysis::ibdseg::seg_pair_order`. Nine of the thirteen datasets have 16
samples or fewer, where one block makes the two orders identical — which is why this
survived every earlier campaign.

---

## 5. What the two rules are worth

`.seg` row-exactness, all four printed fields, on the 982-row primary capture and at the
two floors that had no part in finding either rule:

| floor | before | after | rows whose two estimate columns are exact |
| --- | ---: | ---: | ---: |
| 3 Mb (default) | 806 | **982** | 982 |
| `--seglength 5` | 755 | **900** | 900 |
| `--seglength 10` | 713 | **832** | 832 |

At every floor row-exactness now equals the number of rows whose estimates are right:
**`PropIBD` contributes no error at all**. Mean `PropIBD` error at 3 Mb went 0.000023 → and
the column is exact. Whatever `.seg` still gets wrong is `IBD1Seg` or `IBD2Seg`, at
`--seglength 5` and `10` only, and that is `18-ibd1-caller.md` §9.

Parity: **436 → 464 of 480**, 28 cases flipped, none regressed. `--ibdseg` went 18/52 to
40/52 and `--related --ibdseg` 5/13 to 11/13.

---

## 6. Negatives measured in the same pass

Recorded so nobody re-runs them. All are corpus scores at the default floor unless stated.

**The caller is at a sharp local optimum.** Forty single-knob perturbations of
`tests/parity/fit/engine.py`'s `Params` — every endpoint offset, `MIN_INFORMATIVE` at
8/9/11/12, the word-dirtiness threshold at 1/3/4/5, both minimum run lengths, both fringe
rules, the reach — and all 32 combinations of the two `IBD1` endpoint rules crossed with
the two IBD1 fringe rules. **None improves exact rows. None beats the committed MAE.** The
committed values are the unique maximum of that 32-cell grid.

**Five knobs the corpus cannot see at all** — identical 982/982/982 and identical MAE:
`bridge_rule="17"`, `gate_end="right"`, `inf2_ibs1b=True`, `ibd1_clip_ibd2=True`,
`clip_before_len=False`. These were settled on the canvases (`17-…` §14), not here, and
that remains the only evidence for them.

**Cutting `IBD1` by all `IBD2` calls rather than the surviving ones** (`ibd1_cut="all"`) —
the natural explanation for `dups`' duplicate pair at raised floors, where the reference
reports `IBD1Seg 0.0000` and we report 0.0436. It is worse everywhere: 982 → 950 at 3 Mb,
900 → 861 at 5 Mb, 832 → 783 at 10 Mb.

**Merging calls separated by less than `--seglength`**, the obvious reading of
`18-…` §9's measured run merge, applied as a post-pass:

| variant | 3 Mb | 5 Mb | 10 Mb |
| --- | ---: | ---: | ---: |
| committed | **982** | **900** | **832** |
| merge IBD1 calls | 781 | 748 | 712 |
| merge IBD2 calls | 982 | 898 | 807 |
| merge both | 781 | 746 | 707 |

Merging IBD2 does improve `IBD1Seg` at 5 Mb (910 → 914) and the worst row at 10 Mb
(0.0916 → 0.0424) while losing exact rows and MAE — a trade, not a win, and refused under
the landing rule. `18-…` §9's conclusion stands: the reference requires a further condition
the canvas has not yet been asked for.

---

## 7. What to try next

The whole remaining `.seg` residual is **`IBD1Seg` and `IBD2Seg` at `--seglength 5` and
`10`** — 210 and 81 rows respectively across the 11 failing captures, 0 extra and 0 missing
pairs, and nothing at all at the default floor. Two concrete handles:

1. **`dups`' duplicate pair is one row and the largest single error in the corpus**
   (0.0641 at 5 Mb, 0.0916 at 10 Mb). The reference reports the *same* `IBD2Seg 0.9877` at
   both raised floors while ours falls 0.9018 → 0.8804; its answer is floor-independent
   above 3 Mb and ours is not. One pair, one obvious question: what does the reference do
   with an IBD2 call that the floor would drop?
2. **The run merge of `18-…` §9**, whose gap condition and 5-opposite-homozygote bound are
   both bisected against the reference but whose missing side condition is not. §6 above
   prices the unconditioned version; `docs/research/fixtures/ibd1canvas.py` is the rig, and
   the 17 open canvases it already fails (`gradebinary.py --ibd1`, 43/60 on that family)
   are the smallest reproductions in the tree.
