# The informativeness gate — solved

**Status:** measured, then validated out-of-sample, and now **implemented** —
`MIN_INFORMATIVE` in `crates/open-king-core/src/ibdseg.rs`, gating both `Scan::runs` (with
`inf1`) and `Scan::ibd2` (with `inf2`), with `MIN_RUN1` dropped to 1. §10 below described
the patch before it landed; it landed as written. The corpus scorecard that resulted is in
`docs/PARITY.md` §4.1 and differs slightly from §1's table here, because `Scan::ibd2` was
rewritten concurrently — that trade is described in `docs/PARITY.md`.

Every statement is a reading taken off
the reference binary `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(KING 2.3.2), or a score computed against the captured corpus. **KING's C++ source was not
downloaded, opened or read.**

This document closes the open problem carried by `10-segment-rule-fixtures.md` §5/§6,
`11-segment-rule-fit.md` §7/§11 and `13-segment-acceptance.md` §2/§7 — the 188 extra
`.seg` rows and, with them, the "which two-word runs are accepted" question.

---

## 1. The rule

```
WORD = 64 markers of the global retained-autosome grid
carriesA1(s, m) = sample s is A1A1 or heterozygous at marker m   (plane1)
hom(s, m)       = sample s is homozygous at marker m             (plane0)

inf1(i, j, m) = carriesA1(i,m) & carriesA1(j,m) & (hom(i,m) | hom(j,m))
inf2(i, j, m) = carriesA1(i,m) & carriesA1(j,m)

A run of good words [u..v] inside a usable segment becomes a segment only if

        popcount over markers 64u .. 64(v+1)-1  of  inf  >=  10

with inf = inf1 for an IBD1 run and inf2 for an IBD2 run.
Runs that fail are dropped outright — not shortened, not merged, not re-scored.
```

* **The constant is 10**, and the test is `>=`.
* **The aggregation window is the run's own complete words** — `[64u, 64(v+1))`. Not a
  per-word rate, not a sliding window, not the reported (extended) interval. Markers in
  the flanking words lengthen the segment but contribute nothing to the count.
* **`A1` is the first allele column of the `.bim`**, taken literally. Not the minor
  allele, not a frequency estimated from the cohort.
* **There is no minimum run length.** A single clean word carrying 10 informative markers
  is a valid IBD1 segment; `MIN_RUN1 = 2` was this gate in disguise.

Read `inf1` as *"markers at which an IBS0 could have been observed"*: an IBS0 needs one
sample A1A1 and the other A2A2, so inside a clean word "at least one of the pair is A1A1,
both called" and `inf1` are the same set. The gate is therefore the natural statistical
requirement that the absence of an IBS0 only counts as evidence where an IBS0 had enough
opportunities to appear.

### Scorecard, all 982 captured rows

| | IBD1 **and** IBD2 exact | all three columns exact | IBD1 ok | IBD2 ok | **extra** | **missing** | MAE(PropIBD) | worst |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| committed engine (no gate, min run 2) | 626 | 559 | 626 | 824 | **188** | 0 | 0.00327 | 0.1753 |
| **gate, C = 10, min run 1** | **824** | **709** | **825** | 824 | **0** | **0** | **0.00148** | 0.0971 |

The first column is `11-segment-rule-fit.md` §9's metric — its headline 626 — so the
comparable number is **626 → 824 of 982**. The second additionally requires the printed
`PropIBD` to match, which 67 of those 626 rows already failed before this change and 115
of the 824 still do; that residual is the ±1 word boundary ambiguity of §9, not the gate.

**The set of reported pairs is now exactly right on all ten datasets** — no extras, no
omissions — which is what the 188 rows were.

---

## 2. Method: a rig where every genotype is written by hand

`rig2.py` draws the block's genotypes from an allele frequency, which confounds four
things at once: the per-variant frequency across the cohort, the pair's heterozygosity,
its HetHet count and its het-vs-hom count. `gatelab.py` keeps `rig2.py`'s layout — chr1 a
fully-IBD1 carrier so the pair always clears the 10 Mb filter, chr2 a canvas forced to
opposite homozygotes at every marker with a block of `W` complete words carved out — but
sets **every genotype of every sample explicitly** inside the block. A threshold measured
in "number of markers of kind X" is then a threshold in the statistic itself, with no
sampling noise to fit through.

The read-out is quantitative, not just yes/no: a run of `W` clean words on the solid
background reports exactly `64(W+1) - 1` marker intervals for IBD1 and `64W - 1` for IBD2
(`10-segment-rule-fixtures.md` §2.2/§3), so the reported length says *which* words were
called.

```bash
cd docs/research/fixtures
python3 gate2.py vocab     # which pair genotypes support a call
python3 gate3.py perword   # threshold vs run width, deterministic counts
python3 gate4.py exact     # the constant, three layouts x five widths
python3 gate4.py window    # do flanking-word markers count?
python3 gate4.py weight    # per-genotype weights
python3 gate5.py allele    # A1 column vs minor allele
python3 gate5.py ibd2 ; python3 gate6.py fine ; python3 gate6.py mix
python3 gate7.py sweep     # the "one clean word" open item of 10-...fixtures §2.4
```

One change to `fixlab.py`: a `Fixture.noflip` set exempting chosen markers from the
A1-minor re-orientation, so a fixture can put the **major** allele in the A1 column. It
defaults to empty and changes nothing for existing callers.

---

## 3. The per-marker predicate

A `W = 4` block filled with one repeated genotype vector, all six samples written
explicitly. "chr2" is the marker intervals reported on the canvas; a fully-called block is
319 for IBD1, 255 for IBD2.

| pair genotypes | cohort | chr2 | reading |
| --- | --- | ---: | --- |
| A1A1 / het | A2A2 ×4 | 319 | called, IBD1 |
| het / A1A1 | A2A2 ×4 | 319 | called, IBD1 |
| A1A1 / A1A1 | het, A2A2 ×3 | 63 + 255 | called (IBD2, plus its IBD1 extension word) |
| het / het | anything | 255 | called, IBD2 |
| **het / A2A2** | het ×4 | **0** | refused |
| **A2A2 / A2A2** | A1A1 ×2, het ×2 | **0** | refused |
| **A2A2 / A2A2** | A2A2 ×4 (monomorphic) | **0** | refused |
| A1A1 / missing | A2A2 ×4 | 0 | refused — a missing call never counts |

Sweeping the *number* of each kind among filler markers (`gate4.py weight`, threshold at a
weight-1 marker is 10):

| candidate marker | first n that passes | weight |
| --- | ---: | ---: |
| A1A1 / het, het / A1A1, A1A1 / A1A1 | 10 | **1** |
| het / het | never (tested to 12) | **0** for IBD1 |
| het / A2A2, A2A2 / A2A2, A1A1 / missing, het / missing, missing / missing | never | **0** |

So for an IBD1 run: **both members carry the A1 allele and at least one of them is
homozygous.** Note `A1A1 / A1A1` weighs 1, not 2 — it is a count of markers, not of shared
alleles.

### 3.1 It is the `.bim`'s A1 column, not the minor allele, and not the cohort

Two controls with the A1-minor re-orientation switched off (`gate5.py allele`), so the
block's markers reach the reference with A1 as the **major** allele:

| block markers | first n that passes |
| --- | --- |
| pair A2A2/A2A2, cohort A1A1 — pair is hom for the **minor** allele | never (to 12) |
| pair A2A2/het, cohort A1A1 | never (to 12) |
| pair A1A1/het, cohort A2A2 — control | 10 |
| pair A1A1/A1A1, cohort A2A2 — control | 10 |

A pair homozygous for the *minor* allele counts for nothing when that allele sits in the
A2 column. The gate reads the allele **columns**, and KING's insistence that A1 be the
minor allele (`Too many first alleles as the major allele…`) is what makes that a
frequency filter in practice.

And the cohort itself is irrelevant: holding the pair at het/het and sliding the padding
samples from all-A2A2 (marker MAF 2/12) to all-het (MAF 6/12) leaves the call identical at
every point (`gate2.py freq`).

---

## 4. The aggregation window: a total over the run, and only over its own words

Deterministic counts, `k` informative markers per word evenly spread, filler het/A2A2
(`gate3.py perword`). `Y` = the block is called:

```
  W  k=0  1    2    3    4    5    6    7    8    9   10   11 ...
  1  .    .    .    .    .    .    .    .    .    .    Y    Y
  2  .    .    .    .    .    Y    Y    Y    Y    Y    Y    Y
  3  .    .    .    .    Y    Y    Y    Y    Y    Y    Y    Y
  4  .    .    .    Y    Y    Y    Y    Y    Y    Y    Y    Y
  6  .    .    Y    Y    Y    Y    Y    Y    Y    Y    Y    Y
  8  .    .    Y    Y    Y    Y    Y    Y    Y    Y    Y    Y
 10  .    Y    Y    Y    Y    Y    Y    Y    Y    Y    Y    Y
 14  .    Y    Y    Y    Y    Y    Y    Y    Y    Y    Y    Y
```

The smallest passing **total** `k·W` is 10 at `W` = 1, 2 and 10; the largest failing total
is 9 (`W` = 3, `k` = 3). Placing the markers by explicit total instead of per word
(`gate4.py exact`) makes it exact — the flip is at 10 for every width and every layout:

| W | 6 | 7 | 8 | 9 | **10** | 11 | 12 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 2 / 3 / 4 / 8 / 14, spread | . | . | . | . | **Y** | Y | Y |
| 2 / 3 / 4 / 8 / 14, packed into the first word | . | . | . | . | **Y** | Y | Y |
| 2 / 3 / 4 / 8 / 14, split between the two ends | . | . | . | . | **Y** | Y | Y |

Ten markers packed into the *first word* of a fourteen-word run pass. So there is **no
per-word component at all** — hypotheses (a), (c) and (d) of the brief are dead, and (b) is
right in shape: the window is larger than a word. It is exactly the run.

### 4.1 The window stops at the run's own words

The reported segment reaches into the flanking words (`10-…fixtures` §2.2). Those markers
do **not** count. Left flank IBS0 except for its last `j` markers, right flank informative
in its first `j` (`gate4.py window`), so all `2j` sit inside the reported interval:

| informative in the core | j per flank | total inside [lo,hi] | chr2 |
| ---: | ---: | ---: | ---: |
| 8 | 0 / 1 / 2 / 4 | 8 / 10 / 12 / 16 | **0 / 0 / 0 / 0** |
| 9 | 0 / 1 / 2 / 4 | 9 / 11 / 13 / 17 | **0 / 0 / 0 / 0** |
| 10 | 0 / 1 / 2 / 4 | 10 / 12 / 14 / 18 | 191 / 192 / 193 / 195 |

A sub-threshold core is never rescued, and the length still grows by `j` — the two are
independent. The corpus agrees: counting over every word the reported segment touches
never reaches parity (186 extra at its best constant), and counting the usable segment's
own markers in incomplete edge words leaves 11 extra. Only the run's complete words give
0/0.

Missing calls inside a run neither break it nor count (`gate4.py miss`): 128 markers of
missing/missing in a 4-word run leave the call at its full 319.

---

## 5. The IBD2 side

An IBD2-eligible run (no IBS0 and no het-vs-hom anywhere, so only IBS2 genotypes occur)
with `k` evidence markers among A2A2/A2A2 filler (`gate6.py fine`, `mix`):

| a = HetHet | b = A1A1/A1A1 | a + b | called |
| ---: | ---: | ---: | --- |
| 0 | 9 | 9 | no |
| 9 | 0 | 9 | no |
| 5 | 4 | 9 | no |
| 4 | 5 | 9 | no |
| 0 | 10 | 10 | **yes** |
| 10 | 0 | 10 | **yes** |
| 5 | 5 | 10 | **yes** |
| 9 | 1 | 10 | **yes** |

**Same constant, one pooled count, broader mask.** HetHet and A1A1/A1A1 both weigh 1;
A2A2/A2A2 weighs 0. Within an IBD2-eligible word that set is exactly "both carry A1", so
`inf2` drops the `(hom(i) | hom(j))` clause that `inf1` carries. HetHet is worth 0 to an
IBD1 run and 1 to an IBD2 run — the two masks are genuinely different, verified in both
directions.

The corpus cannot see this difference (every IBD2 variant scores 709/0/0), so the fixture
decides — but the corpus does *prefer* it: substituting `inf1` for `inf2` on the IBD2 pass
triples the mean PropIBD error, 0.00192 → 0.00610.

---

## 6. There is no minimum run length

`10-segment-rule-fixtures.md` §2.4 left "exactly one clean word" open: a reproducible but
irregular kept/dropped pattern, a function of the word index, not explained by physical
length or `--seglength`. It is the gate. Re-running that sweep and computing the
informative count per word (`gate7.py sweep`):

```
maf 0.5, one clean word swept across the canvas       (a 1-word run reports 127 intervals)
 word   informative  predicted   reference
    1             6       drop           0
    2             8       drop           0
    3            10       call         127
    4             9       drop           0
    5            12       call         127
    6             9       drop           0
    7             7       drop           0
    8            12       call         127
```

24 of 24 across MAF 0.5 / 0.3 / 0.2, with the knife edge visible at 9 against 10. A
deterministic 1-word fixture confirms it directly: `k = 9` reports nothing, `k = 10`
reports its 127 intervals.

On the corpus, dropping `MIN_RUN1` from 2 to 1 **with** the gate is worth 112 exact rows
(597 → 709) and 151 IBD1 columns (674 → 825), still at 0 extra and 0 missing. Without the
gate it is a disaster (316 exact). That is the whole reason the two-word floor ever fitted.

---

## 7. Out-of-sample validation

The constant was measured on constructed fixtures. Everything below is held out from it.

### 7.1 The corpus separates perfectly, and densely, at 10

For every pair, `M` = the largest informative count over a run whose reported segment
clears 10 Mb. The pair is reported iff `M >= 10`:

| M | reference refuses | reference reports |
| ---: | ---: | ---: |
| 3–8 | 126 | 0 |
| **9** | **62** | **0** |
| **10** | **0** | **60** |
| 11–19 | 0 | 210 |
| ≥ 20 | 0 | 712 |

1 170 pairs, **zero overlap**, and both sides of the edge are heavily populated — 62
refusals at exactly 9, 60 acceptances at exactly 10. The constant is pinned to the unit by
data that had no part in choosing it. Sweeping it: C = 9 → 62 extra, **C = 10 → 0 extra,
0 missing**, C = 11 → 60 missing.

### 7.2 The mask shape is confirmed too

Eight candidate masks, each given **its own** best constant (`gate_masks.py`), scored on
pair inclusion over all 982 rows:

| mask | best C | extra | missing | max refused / min accepted |
| --- | ---: | ---: | ---: | ---: |
| **both carry A1 & ≥1 hom (fitted)** | **10** | **0** | **0** | **9 / 10** |
| both carry A1 | 25 | 129 | 10 | 39 / 20 |
| both A1A1 | 1 | 110 | 73 | 4 / 0 |
| HetHet | 12 | 182 | 4 | 31 / 7 |
| ≥1 hom, either allele | 107 | 125 | 55 | 175 / 94 |
| both non-missing (N_SNP) | 0 | 188 | 0 | 192 / 128 |
| both carry A1 & both hom | 1 | 110 | 73 | 4 / 0 |
| ≥1 carries A1 | 76 | 119 | 49 | 93 / 62 |

The fitted mask is the only one that separates the corpus at all, by a wide margin.

### 7.3 The word-grid sweep — the sharpest test

`13-segment-acceptance.md` §5's instrument: delete the first `m` markers of the fileset,
which shifts the global word grid by `m` and changes nothing else, so genotypes, refined
endpoints and segment length in bp are all constant along a sweep while the verdict moves
in both directions. 8 borderline pairs × 64 shifts = **512 labelled reference
invocations** (`gate_sweep.py`).

**Agreement: 511 / 512 = 99.8 %.** Confusion (reference, predicted): 358 both-refuse, 153
both-accept, 0 false accepts, 1 false refusal.

Excerpt for `B01_F/B22_F`, the pair `13-…acceptance` §5 tabulates, where the reference
flips between `m = 13` and `m = 14` with identical `lo`, `hi` and length:

```
 m      0..13   14 15 16   17..49   50..63
 ref    no      YES        no       YES
 pred   no      YES        no       YES
```

The one disagreement is at a shift where the run has a single complete word and the
reference calls it — i.e. the pre-existing `MIN_RUN1` question of §6, in the version of
the model still carrying `min1 = 2`; it is not a gate failure.

### 7.4 Held-out `--seglength` captures

The same rule, unchanged, against the corpus captured at other reporting floors:

| capture | rows | exact | extra | missing |
| --- | ---: | ---: | ---: | ---: |
| `__ibdseg` (3 Mb) | 982 | 709 | 0 | 0 |
| `__ibdseg_seglength5` | 982 | 701 | 0 | 0 |
| `__ibdseg_seglength10` | 982 | 668 | 0 | 0 |

---

## 8. What this overturns

1. **`10-segment-rule-fixtures.md` §5.3 — "informativeness comes from the other samples"
   is an artefact of the fixture builder, not a property of KING.** `fixlab.py`
   re-orients A1 to the observed minor allele *using all six samples*, so changing the
   padding cohort changes which markers get flipped, which changes the A1 column, which
   changes the gate — with the pair's dosages byte-identical. Reproduced directly
   (`gate2.py vocab`): pair dosage (0,0) with cohort A1A1 flips and is called; the same
   (0,0) with a cohort that does not trigger the flip is refused. KING never estimates a
   frequency.
2. **`13-segment-acceptance.md` §3 is therefore not in conflict with it.** `subset.py`
   copies the `.bim` verbatim, so the A1 column — and the verdict — cannot move. Both
   observations are consequences of one rule.
3. **`13-segment-acceptance.md` §6's "marker informativeness ruled out" was measuring the
   wrong statistic.** `informative.py` compares MAF-derived quantities over the *reported*
   segment. The gate is a pair statistic over the *core words*, keyed on allele columns.
4. **`MIN_RUN1 = 2` is wrong**; the floor is one word (§6).
5. **`11-segment-rule-fit.md` §8's "poisoned" datasets are largely rehabilitated.** The
   reference was not losing IBD — it was applying this gate, and it bites hardest exactly
   where evidence is thinnest. On `nuclear` `N_C1/N_C3` the reference prints `IBD1Seg
   0.1057`, the old engine 0.5499, and the gated caller **0.1240**; `nuclear` MAE goes
   0.04275 → 0.00793 and `missing` 0.03993 → 0.00650, with nothing tuned on either.
   `nuclear` and `missing` can be used for fitting again. `monomorphic` still has one bad
   row (`P_C2/P_C3`) and that one is IBD2 geometry, not this.
6. **The 188 extras are not a "two-word run" phenomenon at all.** Two-word runs were
   simply the width at which a typical corpus run sits near ten informative markers.

---

## 9. What is still open

1. **±1 word at segment ends.** 709 of 982 rows are exact; the residual is the boundary
   ambiguity `11-…fit` §9 already describes, now the only thing between this caller and
   byte parity. `IBD1 ok` 825 and `IBD2 ok` 824 against 709 both-exact says the two
   columns fail on disjoint rows.
2. **Whether `inf1` is `carriesA1_i & carriesA1_j & (hom_i | hom_j)` or "at least one is
   A1A1 and both are called".** The two masks are identical on every word a run can
   contain (they differ only where an IBS0 is present, which disqualifies the word), so no
   experiment inside a run can separate them. Either is safe to implement.
3. **Where the gate sits relative to the other filters.** Implemented as: gate the run,
   then clip against the previous segment, then apply `--seglength`. That ordering scores
   0 extra / 0 missing; the alternatives have not been separated.
4. **Whether `--roh` uses the same constant** on one sample instead of two.

---

## 10. The patch

`crates/open-king-core/src/ibdseg.rs` was **not** changed by this investigation beyond its
documentation: at the time of writing another session was mid-rewrite of `Scan::ibd2`, and
two agents editing one file loses work. The change is small and additive:

* add to `WordDiff` / `Scan` two per-word masks alongside `ibs0`/`ibs1`, from the same
  planes `word_diff` already loads —
  `inf1 = p1i & p1j & (p0i | p0j)` and `inf2 = p1i & p1j`;
* in `Scan::runs`, after the `min_run` check and **before** the overlap clip, reject a run
  `[k0..k1]` whose `inf` popcount over words `k0..=k1` is under `MIN_INFORMATIVE = 10`,
  with `inf1` for the IBD1 pass and `inf2` for the IBD2 pass;
* set the IBD1 run floor to 1 word.

The reference implementation, scored end to end, is `tests/parity/fit/rules3.py`.

```bash
cd tests/parity/fit
python3 gate_corpus.py            # the scorecard of §1, in three stages
python3 gate_corpus.py sweep      # the constant x scope sweep
python3 gate_masks.py             # §7.2, and the histogram of §7.1
python3 gate_minrun.py            # §6
python3 gate_final.py sweeplen    # §7.4
python3 gate_sweep.py 4 4 64      # §7.3 — 512 reference invocations, ~45 s
```
