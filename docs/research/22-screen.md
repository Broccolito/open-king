# 22 — the `--related` screening count (resolved)

> **Resolution, 2026-08-18.** The black-box campaign below accurately falsified many
> plausible screens but drew the wrong final conclusion because it did not have KING's
> marker score or its unstable tie order. A later arm64 disassembly established the
> complete mechanism: `ConvertLGtoSLG` ranks markers by genotype-class informativeness
> with libStatGen `QuickIndex::Sort`; degree-specific outlined loops apply progressive
> homozygote/IBS0 gates and robust-kinship comparisons. Degree 1 additionally counts the
> unused bits in its last packed word while reconstructing IBS1. The implementation now
> matches all **480/480** captures and the held-out 12,500-marker sparse fixture exactly.
> Sections 0–21 are retained as a historical falsification record; §22 supersedes their
> open-mechanism conclusions.

The one stdout line the corpus still gets wrong, on two cases (`core/bigish__related_degree2`
and `ibdseg/bigish__related_degree2_ibdseg`):

```
-  Stages 1&2 (with 32768 SNPs): 36 pairs of relatives are detected (with kinship > 0.0625)
+  Stages 1&2 (with 32768 SNPs): 50 pairs of relatives are detected (with kinship > 0.0625)
```

The task set for this round was: find which 32768 SNPs the screen uses when `m > 32768`,
and how the kinship is computed on them.

**The answer is that there is no such subset.** This round closes the "find the subset"
line of enquiry with four independent refutations, one of them algebraic, and replaces it
with a law that is measured to 0.1 % but whose mechanism is still open. Nothing was
landed: no rule survived that would generalise off `bigish`, and fitting one to `bigish`
would be a fitted fiction. `related.rs` is unchanged apart from its documentation.

Instrument: `docs/research/fixtures/screendeflate.py` (`facts` re-measures everything
below; each subcommand re-measures one section). It builds on `screencanvas.py` and
`screenweight.py` and supersedes several of their conclusions where noted.

---

## 0. The instrument

`screencanvas.py`'s clone canvas tunes a pair's kinship by cloning marker sets. This round
adds a **dilution bisection**, which is both finer and usable on `bigish`'s *real*
between-family pairs rather than only on constructed ones:

> Take one `bigish` between-family pair in a 169-sample fileset — the 167 fillers that on
> their own print `No close relatives are inferred.`, plus the pair. Replace one member's
> genotypes, at a growing random marker set, with genotypes drawn from that fileset's own
> allele frequencies. The replacement is **synthetic and unrelated to everybody**, not a
> copy of another sample, so no new relative pair appears and the printed count stays a
> clean read-out of this one pair. Kinship falls linearly in the replaced-marker count;
> bisecting it locates the acceptance boundary to one marker, i.e. to ~1e-5 in kinship.

The boundary is quoted as the deflation `R` of `screencanvas.py` §2:

```
k_screen = 0.5 + R*(k - 0.5)      R = (cut - 0.5) / (k* - 0.5)
```

One bisection is ~17 reference runs and a reference run on `bigish` is ~90 ms, so a
36-bisection mean costs about a minute. That is what makes the numbers below tight enough
to separate hypotheses that the earlier instruments could not.

**A trap.** KING rejects a fileset whose A1 is the major allele —
`FATAL ERROR - Too many first alleles as the major allele (~100.0%)` — and also rejects a
50/50 mix (`~49.2%`). A synthetic generator must follow PLINK and code A1 as the *minor*
allele. Get this wrong and every run dies, every bisection reports "no bracket", and the
few filesets that squeak past the check produce garbage `R`. Two hours of this round were
spent chasing a "deflation at m = 32768" that was only this.

---

## 1. The law: affine in `0.5 - k`, one `R` per fileset

`bigish`, m = 50000, n = 169, 36 dilution bisections at each cutoff:

| cutoff | `R` | boundary kinship `k*` |
|---|---|---|
| 0.0625 (`--degree 2`) | **1.02257 ± 0.00065** | 0.07216 |
| 0.1250 (`--degree 1`) | **1.02079 ± 0.00062** | 0.13264 |

The two agree to 0.2 %. The competing readings do not survive: a multiplicative
`k_screen = c*k` needs `c` = 0.866 and 0.943 at the two cutoffs; a constant offset needs
0.0097 and 0.0076.

The lever arm is 25× bigger on a synthetic flat-MAF fileset, where the deflation is four
times larger (`screendeflate.py boundary`):

| fileset | cutoff 0.0625 | cutoff 0.1250 |
|---|---|---|
| flat MAF 0.25, m = 50000, n = 140 | `R` = 1.0798, `k*` = 0.0948 | `R` = 1.0838, `k*` = 0.1540 |

`R` agrees to 0.4 % across the two cutoffs while `cut/k*` — what a multiplicative rule
would have to hold fixed — moves from 0.659 to 0.812. **The law is affine about 0.5.**

Equivalently, with `D = 1 - 2*phi` (the non-IBD allele fraction, which is what the KING
estimator actually accumulates: `D = N / (2*min(het_i, het_j))` with
`N = het_i + het_j + 4*IBS0 - 2*HetHet`), the screen computes `R * D`. The screen inflates
the pair's *discordance* by a constant factor, and everything below is about that factor.

**`R` is exactly 1 whenever `m <= 32768`**: 0.99999 ± 0.00001 on `bigish`, and 1.00000,
1.00001, 0.99997 on three synthetic MAF spectra. This is the constraint that kills most
candidate estimators outright, and it is now confirmed on four independent filesets rather
than one.

---

## 2. The deflation is systematic, not sampling noise

Two independent readings, because the whole question is whether a subset's sampling error
could masquerade as a bias.

* **Realisation spread.** Eight dilution seeds on each of four pairs: the boundary kinship
  has sd **0.0018** while the deflation is **0.0089**. Five to one.
* **Per-pair labels, no bisection at all.** Rebuild `bigish` 59 times as its 167 fillers
  plus exactly one candidate pair and read the printed count (`screendeflate.py labels`).
  The accept/reject split is a *sharp threshold* at whole-map kinship ≈ 0.0719:

  ```
  ... 0.07456 ACCEPT  0.07339 ACCEPT  0.07301 ACCEPT  0.07274 reject
      0.07193 ACCEPT  0.07192 ACCEPT  0.07181 reject   0.07044 reject ...
  ```

  One inversion, inside a 0.0009-wide window. A subset-sampling model with sd 0.0018 would
  scatter accepts and rejects over a ~0.007-wide band and would flip pairs at 0.08 and at
  0.063; it does not happen.

The degree-2 column sums to **36**, the whole fileset's number — the stage is per-pair, as
`screencanvas.py` found. At degree 1 the same sum is **17** against the fileset's 18, so
"per-pair" is exact to ±1 and there is a weak dependence on the rest of the sample.

---

## 3. It is not the kinship over any marker subset — four refutations

### 3.1 Algebraically, every subset and every weighting is unbiased

For a pair with IBD probabilities `(k0, k1, k2)` and `phi = k2/2 + k1/4`, at a marker of
frequency `p`:

```
E[N_l]   = 4pq + 8*k0*p²q² - 2*(k2*2pq + k1*pq + k0*4p²q²) = 4pq(1 - 2phi)
E[het_l] = 2pq
```

The `p²q²` terms cancel **exactly** — the same identity `4pq(p²+q²) + 8p²q² = 4pq` that
`screenweight.py` §4 used. So `N` and `min(het_i, het_j)` are both proportional to `Σ pq`
over whatever index set they are summed on, their ratio is free of that set, and

> the KING robust kinship over **any** marker subset, and under **any** non-negative
> per-marker weighting, is an unbiased estimate of the same `phi`.

There is no subset and no weighting that deflates. This is not a search that came up
empty; it is a proof that the search cannot succeed.

### 3.2 Measured on `bigish`

`kinship` over top-K-by-MAF subsets, counting between-family pairs over the cutoffs:

```
m=50000  top-50000: 47 / 21     m=16384  top-16384: 41 / 13
         top-32768: 45 / 23              top-10922: 41 / 13
         top-25000: 44 / 23              top- 8192: 41 / 13
         top-16384: 48 / 22              top- 5461: 40 / 13
```

Flat, as §3.1 says it must be. The reference prints 36 / 18.

### 3.3 Replicating a map holds every kinship *bit-identical* and the screen still moves

`r` exact copies of one map (distinct positions, identical genotypes) multiplies every
count in the estimator by `r` and leaves every kinship unchanged to the last bit — KING's
own `.kin0` confirms it. So any change in the count is the screen, with the data held
fixed. From `bigish`'s first 16384 markers, 41 pairs over 0.0625:

| copies | m | printed count | implied threshold | `R` |
|---|---|---|---|---|
| ×2 | 32768 | 41 | — | 1.000 |
| ×3 | 49152 | 36 | 0.0716 | 1.021 |
| ×4 | 65536 | 33 | 0.0780 | 1.037 |
| ×5 | 81920 | 29 | 0.0854 | 1.055 |
| ×6 | 98304 | 27 | 0.0892 | 1.065 |

Any sub-multiset of a replicated map is a weighting of the base map, hence unbiased by
§3.1, hence 41. The reference gives 36 at ×3.

### 3.4 Selecting on an in-sample statistic does not rescue it either

The one loophole in §3.1 is a subset chosen from data that *includes the pair*, which is
no longer a fixed index set. Simulated (`screendeflate.py subset`), on a flat-MAF map where
the selection is pure sampling noise and the ascertainment is therefore maximal:

```
top-32768 by in-sample MAF              R = 0.995 ± 0.002    (no bias)
top-32768 by in-sample heterozygote count  R = 0.916 ± 0.003 (bias of the wrong sign)
```

So `screencanvas.py`'s "the estimator reads sample-level allele frequencies" is right that
the screen is sample-dependent, but the dependence is not this.

---

## 4. What the deflation *does* track

### 4.1 It needs the discarded markers to be informative

Take `bigish`'s first 32768 markers (`m = 32768`, screen exact, prints 50 / 18) and append
17232 more (`screendeflate.py append`):

```
base 32768 alone                        deg2 = 50  deg1 = 18
+ 17232 markers at MAF 0.02             deg2 = 50  deg1 = 18     <- nothing at all
+ 17232 markers at MAF 0.05             deg2 = 45  deg1 = 17
+ 17232 markers at MAF 0.10             deg2 = 33  deg1 = 13
+ the real bigish tail                  deg2 = 36  deg1 = 18
```

The MAF-0.02 row is exact, and the base map has 42 pairs over 0.0719 — so a deflation of
the `bigish` size would have shown as 42, not 50. Appending junk costs **zero**.

### 4.2 It scales with how much equally-informative material overflows the budget

Synthetic two-point maps, `m = 65536`, MAF 0.45 on `K` markers and 0.12 on the rest, so the
32768 budget lands cleanly between the groups exactly at `K = 32768`:

```
K = 20000   R = 1.0283       K = 33536   R = 1.0136
K = 28000   R = 1.0153       K = 40000   R = 1.0596
K = 32768   R = 1.0081  <-   K = 50000   R = 1.0776
                             K = 65536   R = 1.0612
```

`R` has its minimum exactly where the budget does not have to split a tied group, and
climbs steeply once it does.

The same fact from the other side — `R` against the MAF spectrum at fixed `m`
(`screendeflate.py spectrum`, n = 140):

| spectrum | m=32768 | m=40960 | m=50000 | m=65536 |
|---|---|---|---|---|
| flat 0.25 | 1.00001 | 1.0586 | 1.0798 | 1.0809 |
| flat 0.45 | 1.00000 | 1.0606 | 1.0845 | 1.0612 |
| uniform(0.05, 0.5) | 1.00002 | 1.0150 | 1.0333 | 1.0439 |
| beta(0.6, 2.2) | 0.99997 | 0.9980 | 1.0069 | 1.0143 |
| `bigish` (n = 169) | 1.00000 | 1.0115 | 1.0216 | 1.033 (×4 replicate) |

A flat spectrum — every marker equally informative, so the budget can only be met by
discarding useful markers — deflates four times harder than `bigish`. A spectrum with a
long low-MAF tail, where the budget is met by discarding junk, barely deflates at all, and
one point sits *below* 1 (0.9980 ± 0.0003), which retires `screencanvas.py`'s "`R` was
never measured below 1".

### 4.3 `R` against `m` and `n`, on `bigish`

Dilution bisections, 8 pairs each:

```
m     32768   34816   36864   38912   40960   43008   45056   47104   49152   50000
R     1.0000  1.0030  1.0087  1.0126  1.0128  1.0156  1.0166  1.0202  1.0203  1.0216
```

and against the sample count at m = 50000: `R` = 1.0314 (n = 112), 1.0231 (n = 142),
1.0216 (n = 169). Below n = 100 the screening path does not run at all
(`SCREEN_MIN_SAMPLES`). The n-dependence is real but not clean: `(R-1)*n` is 10.2, 9.5,
10.6 across those three, and the full 200-sample `bigish` needs `R` ≈ 1.022, i.e. no fall
at all from n = 169. Two configurations with the same `m` and `n` differ by 4.5 sd
(`bigish` 50000 → 1.0216, `bigish`-25000-replicated-×2 → 1.0280), so `R` is not a function
of `(m, n)` alone.

### 4.4 Marker file order matters, but only through ties

Nine permutations of `bigish`'s marker order — same marker set, same kinships — print
36 / 18 every time, which retires prefixes, strides and word decimations for good. The
boundary bisection is forty times finer and does move: 0.07100, 0.07070, 0.07107, 0.07100
over the identity and three permutations, a spread of 0.0004, about 5 % of the deflation.
So the screen has an order-dependent component of exactly the size you would expect from
ties in an informativeness ranking (`bigish`'s 200-sample MAF lives on a 1/400 grid, so
tie groups are large), and no more.

---

## 5. Where that leaves the mechanism

Everything above is compatible with exactly one shape, and it is not a shape any subset
rule can take:

> When the map holds more equally-informative markers than 32768, KING reaches its budget
> by something **lossy applied uniformly across markers** — every marker keeps a reduced
> share of its evidence — rather than by keeping some markers and dropping others. On a
> flat-MAF map the loss inflates the pair's discordance `D = 1 - 2*phi` by 8 %; when the
> overflow is junk the loss is nil. The uniformity is measured directly
> (`screendeflate.py blocks`): on a flat-MAF map a contiguous clone block grown from
> marker 0, from marker 20000, from marker 32768, or backwards from the tail hits the
> boundary at kinship 0.09566 / 0.09567 / 0.09300 / 0.09401 — the same to ~2 %, with no
> preference at all for the head of the file.

Candidates that survive this and are worth the next round's time:

1. **Merging rather than selecting.** 32768 slots each carrying an aggregate of `m/32768`
   markers, grouped by informativeness rank, would be uniform, order-invariant up to ties,
   lossless when the partner is junk and maximally lossy on a flat spectrum. It has to be
   reconciled with §3.3, where a stable rank sort separates the `r` copies of one marker by
   `m0` positions so that groups do mix distinct markers — which is consistent, but it has
   not been shown that any concrete merge operation on KING's two bit planes produces
   `R = 1.021` there and `1.080` on a flat map.
2. **Two stages with different budgets, intersected.** "Stages 1&2" is two stages. An
   intersection of two thresholds is biased upward and vanishes when the two stages see the
   same markers — which would give `R = 1` at `m <= 32768` for free. The obstacle is §2:
   the bias would then be of the order of the realisation spread, and it is five times it.
   A *stage 1 with a deliberately tightened cutoff* escapes that, but then it has to explain
   why the tightening is exactly zero at `m = 32768`.

Candidates now closed, so that nobody spends another round on them:

* the map's first 32768 markers, any stride, any word decimation — §4.4;
* top-32768 by in-sample MAF, by heterozygote count, or by any other key, with the KING
  estimator on the selected set — §3.1 (proof), §3.2, §3.3 (measurement);
* any per-marker weighting, including frequency-standardised ones — §3.1;
* ascertainment from selecting on an in-sample statistic — §3.4;
* a multiplicative deflation, a constant offset — §1;
* a noise or intersection-of-noisy-estimates account of the *size* of the effect — §2.

---

## 6. Consequence for the implementation

`related.rs` keeps the whole-map estimate on the map's first `min(m, 32768)` markers. It is
right whenever `m <= 32768` — which §1 now proves is not a coincidence but the reference's
own behaviour — and it is right at degree 1 on `bigish` (18) by luck. It is wrong at degree
2 on `bigish` (50 against 36).

The consequence is contained: the count reaches stdout and nothing else. `.kin0`'s rows
come from the exhaustive re-estimate that follows the screen and are byte-correct at every
degree, including on the two cases whose stdout this line spoils. The reference's screen is
*lossy in the safe direction* here — the 14 pairs it drops all sit below the 0.08839
reporting threshold, so no reported row depends on it.

Landing the affine law with a fitted `R` would reproduce `bigish` and nothing else: §4.2
and §4.3 show `R` swinging from 0.998 to 1.085 with the MAF spectrum, and not being a
function of `(m, n)`. That is the definition of the fitted fiction this project refuses.

---
---

# Round 2 — merging is dead, and the screen reads markers it does not keep

The task set for this round was the hypothesis at the top of §5.1: that the budget is met
by **merging** markers into 32768 slots rather than by selecting them, which would dissolve
§3.1's impossibility proof because a merged slot is not a subset.

It is refuted, three independent ways, one of them quantitative. In its place this round
puts two facts that are harder than anything in Round 1, and that between them rule out
every mechanism so far proposed, including both of §5's survivors:

* the screen carries a **second necessary condition** that has nothing to do with the
  budget — at `m = 32768`, where §1 proves the screen is the exact whole-map kinship to
  4e-6, it refuses pairs whose sharing avoids the informative markers, at kinships up to
  0.31 that KING's own `--kinship` confirms;
* above the budget the screen's statistic **depends on markers it did not keep** — with the
  kept 32768 markers held *bit-identical* and every pair's kinship over them identical to
  the last bit, changing only the discarded markers moves the printed count from 46 to 37.

Still nothing landed. `related.rs` is unchanged apart from its documentation.

Instrument: `docs/research/fixtures/screenfold.py` (`facts` re-measures everything below).

## 7. The instrument: a ladder, not a bisection

A dilution bisection costs ~17 reference runs per number, which is what kept Round 1 to a
handful of configurations. A **ladder fileset** costs one run: build 48 pairs whose
kinships climb through the cutoff in even steps, add unrelated fillers, and read the
printed count. The count *is* the effective threshold, to the ladder's step — ~0.0015 in
kinship, ~0.004 in `R`. Twenty configurations take five seconds. Every scan below uses it;
the bisection is kept only where a single pair has to be tracked.

## 8. Merging is dead

### 8.1 One marker over budget is not a step

`bigish`'s first 32768 markers, then its real tail one marker at a time
(`screenfold.py step`), printed against the whole-map truth:

```
m = 32768  +0     50 / 50      m = 33792  +1024    51 / 51
m = 32769  +1     50 / 50      m = 34816  +2048    47 / 49
m = 32770  +2     50 / 50      m = 36864  +4096    43 / 50
m = 32772  +4     50 / 50      m = 40960  +8192    42 / 47
m = 32832  +64    50 / 50      m = 49152 +16384    38 / 47
m = 33024  +256   50 / 50      m = 50000 +17232    36 / 47
```

Exact through +256 and then a ramp. A block merge — `blockSize = ceil(m/32768)`, the
shape any fixed-slot compression takes — flips every slot to a pair at the first
overflowing marker and so predicts a step at `m = 32769`. There is none.

### 8.2 The same multiset, three arrangements, the same count

`bigish`'s first 32768 markers plus 8192 duplicate markers, arranged three ways
(`screenfold.py arrange`). All three hold the *same marker multiset*, so all three have
bit-identical kinships; only the file order differs.

```
self-aligned  (copies appended in order)        printed 41   true 47
shuffled      (copies appended in random order) printed 41   true 47
adjacent      (each copy after its original)    printed 41   true 47
```

A merge needs a grouping, and every grouping is *lossless* for at least one of these
arrangements under any idempotent operation — or, and, max, "take the more informative":

* `slot j = j mod 32768` merges each of the first 8192 markers with **its own copy** in the
  self-aligned arrangement;
* consecutive-pair grouping does the same in the adjacent arrangement;
* a stable rank sort puts a marker and its copy adjacent (they are tied), so rank-block
  grouping does too.

All three lose the same six pairs. The only grouping left standing after this is a
rank-*stride* (slot `i` takes ranks `i`, `i+32768`, …), which separates the copies.

### 8.3 The surviving grouping, scored

Fold `m` markers into 32768 slots and score the KING estimator on the folded data against
the reference, on the separation scan of §10 (`screenfold.py merge`). The encoding is the
sparse one — hom-major = `00`, het = `10`, hom-minor = `11` — which is the only encoding in
which merging junk is free, as §4.1 requires.

```
   x   printed  true | stride/or  stride/and  stride/xor  stride/sum  block/or
0.05      46      46 |    48          45          40          37         48
0.15      45      46 |    48          48          33          23         48
0.25      39      46 |    48          48          28          15         48
0.35      37      46 |    48          48          20          13         48
0.45      35      46 |    48          48          15          14         48
```

`or` and `and` accept every pair — merging makes pairs look *more* related, the wrong sign
— while `xor` and a saturating dosage sum destroy far too much, and both are already wrong
at `x = 0.05` where the reference is exact. No operation is anywhere near the reference's
46 → 35.

**Merging is closed.** With it goes the whole "the budget is met by something lossy applied
uniformly across markers" reading of §5.

## 9. A second necessary condition, and no budget in sight

At `m = 32768` the screen is the exact whole-map kinship — §1 measures the boundary at
0.0625 to 4e-6, and §7's ladder confirms it on four spectra. It is nevertheless possible to
build a pair the screen simply refuses (`screenfold.py gate`), at `m = 32768` exactly:

```
16384 markers @ MAF 0.20 cloned + 16384 @ MAF 0.45 untouched   kinship 0.20006   refused
16384 @ MAF 0.10 cloned          + 16384 @ MAF 0.45 untouched   kinship 0.13890   refused
16384 @ MAF 0.25 cloned          + 16384 @ MAF 0.45 untouched   kinship 0.21731   refused
24576 @ MAF 0.15 cloned          +  8192 @ MAF 0.45 untouched   kinship 0.30669   refused
```

"Refused" is `No close relatives are inferred.` — not a low count, no `Stages 1&2` line at
all. KING's own `--kinship` on the very same fileset prints the same kinship to four places
(`0.1519`, `0.1530` on the three-stratum variants), so the reference agrees the pairs are
relatives and its screen throws them away regardless.

The whole accept region in the `(fA, fB)` plane — clone fraction inside the MAF-0.20
stratum against the MAF-0.45 one (`screenfold.py region`) — shows the shape:

```
fA\fB  0.00  0.10  0.20  0.30 ... 1.00        whole-map kinship along the fB=0 column
 0.00   .     .     Y     Y        Y                     0.001
 0.20   .     Y     Y     Y        Y                     0.042
 0.50   .     Y     Y     Y        Y                     0.101
 1.00   .     Y     Y     Y        Y                     0.200
```

The boundary is **flat in `fA`** at `fB ≈ 0.045`: past the cutoff, extra sharing among the
MAF-0.20 markers buys exactly nothing, and a pair at kinship 0.20 is refused for want of
~700 shared markers among the informative ones.

What it is not:

* a kinship threshold on a top-K subset — the required `k_T` at the boundary runs from
  0.006 to 0.026 across designs, and no `K` from 512 to 32 768 makes it constant;
* an IBS0 threshold, in any normalisation — the refused pair above carries IBS0 at 0.0602
  per marker (0.1463 per het), while on a homogeneous MAF-0.45 map of the same size a pair
  at 0.1047 per marker (0.2112 per het) and a third of the kinship, 0.06554, is **accepted**;
* a HetHet threshold, or any monotone function of the whole-map kinship — the refused pairs
  score higher than accepted ones on every one of them;
* contiguity or segments — §5's block test already shows a scattered clone set and a
  contiguous one hit the boundary at the same kinship;
* a duplicate/MZ path — the refused pair's genotype concordance is 0.6872, and a real
  duplicate on the very same map is detected (and reported as MZ).

It does not bind for uniformly-related pairs, which is why `bigish` never shows it and why
Round 1's instrument, which only ever diluted uniformly, could not see it. On `bigish` it
is invisible directly: diluting a real pair over the top three quarters by MAF still brackets
at 0.06250 exactly (`screenfold.py strata`).

## 10. Above the budget, the screen reads markers it did not keep

Take `m = 50000` as 32768 markers at MAF 0.45 plus 17232 at MAF `x`, and walk `x`
(`screenfold.py separation`). Through `x = 0.30` this holds the kept set *fixed in every
sense*: the top-32768 by allele count is exactly the MAF-0.45 group with **zero** markers
swapped in, it is the same index set at every `x`, those markers' genotypes are
**bit-identical** across the scan, and therefore every ladder pair's kinship over them is
identical to the last bit. Only the discarded 17232 markers change.

```
   x   printed  true      R   swaps  kept bits identical
0.05      46      46  1.0007      0        -
0.10      46      46  1.0053      0      True
0.15      45      46  1.0074      0      True
0.20      43      46  1.0153      0      True
0.25      39      46  1.0371      0      True
0.30      37      46  1.0535     33      True
0.41      35      46  1.0643   5145      True
0.45      35      46  1.0654      -      True (one tied group)
```

The printed count falls from 46 to 37 with the kept data held bit-identical. **The screen's
statistic is not a function of the markers it keeps.** That is a stronger statement than
§3.1's: not only is it not the kinship over a subset, it is not *any* function of any fixed
subset, because the deflation grows while the candidate subset's bits do not move.

And the loss is **deterministic, not the intersection of two noisy tests**. Labelling all
48 ladder pairs one at a time on the `x = 0.25` map (`screenfold.py sharp`) gives a perfectly
sharp threshold:

```
Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y Y . . . . . . . . . . .
lowest accepted  kinship 0.08319 (over the kept 32768: 0.08791)
highest rejected kinship 0.07848 (over the kept 32768: 0.07902)
inversions: 0 in the whole-map kinship, 0 in the kept-subset kinship
```

Zero inversions over 48 pairs, in a design where a fixed-subset estimate deviates from the
whole map by sd 0.0015 and the observed displacement is 0.018 — twelve times that. §5's
second survivor, "two stages intersected", is closed with it: an intersection of noisy
estimates buys a bias of the order of its own scatter, and here the scatter is zero.

## 11. In-sample ascertainment: real, right sign, and not enough

§3.1's proof has exactly one loophole, and §3.4 probed it with two keys and moved on. It is
worth reopening, because on a **flat-MAF map every marker is exchangeable**, so a subset can
differ from the whole map *only* through a selection that leans on the pair's own genotypes
— and flat maps deflate hardest of all (`R` = 1.08).

Which key matters, and only keys that lean *against* the pair's heterozygosity deflate at
all. Predicted counts on `bigish` for `kinship over top-32768 by <key>`
(`screenfold.py keys`; the reference prints 36 / 18):

```
whole map          47 / 21        hom-minor count    41 / 18
allele count       45 / 23        carriers           47 / 23
het count          53 / 23
```

Ranking by the count of minor homozygotes — which at a fixed allele count is a ranking
*against* heterozygosity — is the only candidate that moves in the right direction, and
offline it reproduces the flat-MAF numbers well (model 1.0615 against a measured 1.0654 at
MAF 0.45, m = 50000; 1.0852 against 1.0798 at MAF 0.25). It also explains §4.2's tie-group
result exactly: `R` is at its minimum where the budget need not choose among equals, and
climbs once it must.

Two measurements stop it being the answer:

* **It is 4× short in the middle of §10's scan.** At `x = 0.25` the ranking resolves with
  zero swaps, so no ascertainment of any kind is available, and the model gives 46 against
  the reference's 39.
* **The `n` scaling is wrong.** A 2-in-`n` ascertainment must fall as `1/n`. On flat MAF
  0.25 at m = 50000, `R - 1` = 0.085, 0.064, 0.051, 0.040 at n = 110, 200, 400, 700 — a
  fall of 2.1× where `1/n` demands 6.4×.

So ascertainment is a real component of the effect and cannot be the whole of it.

## 12. Where round 2 leaves the mechanism

Everything measured is consistent with one shape, and it is not a shape any of the
mechanisms proposed so far can take. The screen's statistic is a **deterministic** function
of the pair which

1. equals the whole-map robust kinship exactly when `m <= 32768` (§1, and the ladder on four
   spectra),
2. equals it still when the overflow is junk (§4.1; §10 at `x = 0.05`, printed = true = 46),
3. is deflated smoothly as the discarded markers approach the kept ones in informativeness,
   with the kept markers' bits held fixed (§10) — so it is computed from more than 32768
   markers, or from something other than genotypes at 32768 markers,
4. is deflated further when the ranking at rank 32768 is unresolvable, in the manner and
   very nearly the size of an in-sample ascertainment (§11, §4.2), and
5. carries a second necessary condition, active at every `m`, that asks where in the MAF
   spectrum the pair's sharing sits (§9).

Closed this round, on top of §5's list:

* **merging or compressing markers into 32768 slots**, under any idempotent operation with
  an index-, rank-block- or rank-stride grouping, and under `or`/`and`/`xor`/saturating-sum
  with either grouping — §8.1, §8.2, §8.3;
* **any function of any fixed marker subset**, ascertained or not — §10 holds the candidate
  subset bit-identical and still moves the count;
* **two stages intersected**, in the noisy-estimates form — §10's labels have zero
  inversions where the model needs scatter twelve times the observed displacement;
* **"the screen is a threshold on the whole-map kinship"** as a complete description — §9.

Worth the next round's time, in order:

1. **Chase §9's second condition to a rule.** It is the only phenomenon here that is
   binary, enormous (a factor of three in kinship) and reproducible in one run, which makes
   it far cheaper to pin than the 2 % deflation — and it is plausibly the same machinery.
   The `(fA, fB)` region is the right canvas; vary the informative stratum's size and MAF
   and find what is conserved along the boundary. Beware two false leads already burned:
   `N = 8*IBS0` and `φ + 2*IBS0/min_het = 0.5` are *identities* for clone-block pairs whose
   untouched markers sit at MAF 0.45, not rules.
2. **Ask what statistic can degrade with the discarded markers' informativeness while the
   kept bits are fixed.** A running estimate over the sorted map with an early exit has that
   shape — the exit test sees partial sums whose *bound* depends on what remains — and it
   would also give §9 for free, since a pair sharing only at the bottom of the ranking looks
   unrelated in every prefix. Round 1's §5.2 dismissed staging on noise grounds; §10's
   sharpness says the staging must be *deterministic*, which a bound-based early exit is.
3. **Do not fit `R`.** §11 shows how easy it is to reproduce a spectrum or two with a key
   chosen after the fact, and §10 shows any such fit is wrong by 4× one column over.

## 13. Consequence for the implementation, unchanged

`related.rs` still estimates on the map's first `min(m, 32768)` markers. That remains right
whenever `m <= 32768` and right at degree 1 on `bigish` (18), and wrong at degree 2 (50
against 36). Nothing about that changed this round, and nothing should until a rule survives
§10. The blast radius is still one stdout line: `.kin0`'s rows come from the exhaustive
re-estimate that follows the screen and are byte-correct at every degree, and the 14 pairs
the reference's screen drops all sit below the 0.08839 reporting threshold.

Harness at the end of that historical round: **477 / 480 byte-identical, self-check
480 / 480**. See §22 for the resolved result.

## 22. Resolution from the executable

The missing pieces were recovered from KING 2.3.2's arm64 symbols
`KingEngine::ConvertLGtoSLG`, `QuickIndex::Sort`, and
`Engine::ScreenCloseRelativesInSubset64Bit`:

1. For every marker, let `A`, `H`, and `C` be the counts of A1 homozygotes,
   heterozygotes, and A2 homozygotes; `p = (A + H/2) / (A+H+C)`. Monomorphic markers
   score `-4`; otherwise the score is
   `max(0.096N, min(H, 4AC/H, 2p(1-p)N_called))`.
2. libStatGen's unstable quicksort orders markers by that score. The screen packs the
   highest-scoring `min(m,32768)` markers in descending order. Reproducing its equal-key
   instability is necessary; a stable sort changes the selected boundary.
3. Degree 2 applies homozygote/IBS0 gates after 1/16, 1/8, 1/4, and all of the packed
   panel, with factors `0.5`, `0.475`, `0.45`, and `0.375`, then requires robust
   between-family kinship above `0.0625`.
4. Degree 1 applies corresponding gates after 1,024, 2,048, 3,072, and 4,096 markers,
   with factors `0.5`, `0.4`, `0.375`, and `0.35`, then requires robust kinship above
   `0.125`. Its final IBS1 reconstruction treats the unused bits in the last 64-bit word
   as missing, lowering the comparison by
   `padding / min(N_Het1,N_Het2)`. On the held-out 12,500-marker panel, 44 unused bits
   remove exactly the two candidates KING removes (17 to 15); on the 32,768-marker
   corpus panel the adjustment is zero.

Differential evidence after landing the transcription:

```text
parity: 480 PASS, 0 FAIL, 480 total
held-out sparse fallback: console, .kin, .kin0, unrelated, cluster, and build exact
```

This is a mechanism-derived implementation, not the rejected fitted affine law.

---
---

# Round 3 — the early exit is dead: rank does not matter, a MAF knee does

The task set for this round was §12.2: the one hypothesis two rounds of negatives had
converged on — a **deterministic bound-based early exit** over the informativeness-sorted
map. It predicts every surviving observation. It reads the discarded markers (they enter
the bound); it is deterministic, so per-pair labels are sharp; `R = 1` below the budget
because no exit is possible; the deflation tracks the MAF spectrum because that is what
the bound's slack depends on; and a conservative bound can only reject, which matches the
sign.

**It is refuted, twice over, and the refutations are cheap and reproducible.**

* An early *accept* on a lower bound can only ever **under**-detect, so `R >= 1` always.
  Measured directly: `R = 0.980` — the screen accepts pairs whose whole-map kinship is
  *below* the printed cutoff (§16). That also settles the question Round 1 left open about
  the beta spectrum's `R = 0.998`: a real acceptance, not measurement slack.
* An early *reject* on a prefix is a statement about **rank**. Rank is not what matters.
  A block of unshared markers ranked *above* the pair's entire sharing costs the pair
  **nothing at all** at MAF 0.40 — at any block size up to 20 480, five eighths of the map
  — and is **lethal** at MAF 0.45 (§17, §18). No scan over a ranking can tell those two
  blocks apart: both are unshared, both outrank the sharing, and they differ by 3 % in
  informativeness.

In their place this round puts the two hardest facts in the file: the screen reads the
pair's **genotypes** at markers the budget discards, and the second condition of §9 is a
**knee in MAF**, not a position in a ranking.

Still nothing landed. `related.rs` is unchanged apart from its documentation.

Instrument: `docs/research/fixtures/screenexit.py` (`facts` re-measures everything below).

## 14. The instrument: two canvases and one correction

**The kept/discarded canvas.** Build `m` markers as `nA` at MAF 0.45 plus the rest at MAF
`x < 0.42`, so the ranking separates the groups with zero swaps, and clone each ladder pair
at fraction `fA` inside the high group and `fB` inside the low one, drawing the two clone
sets from **separate generators keyed on the rung**. Every arm of a scan over `fB` then
holds the high group's genotypes bit-identical — §10's condition — while the *pair's own
genotypes* at the discarded markers change. §10 could only vary the discarded markers'
frequencies; this varies the pair.

**The block canvas.** `u` markers at MAF `a` that the pair never shares, above a window of
`w` markers at a low MAF that it clones at fraction `f`, with junk at MAF 0.03 filling out
to `m`. A 48-rung ladder over `f` reads the accept threshold in one run. Quoting the
threshold **in whole-map kinship** is what makes the scan fair: it already charges the pair
for whatever information the block adds to the map.

**A correction to Round 2's reading.** The accept set is the top `det` pairs by *realised*
kinship, not by ladder index. Near `f = 1` successive rungs differ by less in kinship than
realisation noise, so an index-ordered reading inverts. Every number below reports the
`det`-th largest kinship, which is what `screenfold.py`'s `read_threshold` already did;
single-pair labels away from the boundary agree with it, and a 24-pair labelling run at
`fB = 0.30` puts the flip between kinship 0.0696 and 0.0715 with no inversions.

## 15. The screen reads the pair's genotypes at markers it discards

`m = 50 000` = 32 768 @ MAF 0.45 + 17 232 @ MAF 0.25, `n = 200`, ladder in `fA` inside the
kept group (`screenexit.py tail`):

```
    fB  printed thr(whole) thr(kept)      R   kept bits identical
  0.00       30     0.0589    0.0854  0.9919        -
  0.05       30     0.0665    0.0854  1.0093      True
  0.10       30     0.0734    0.0854  1.0256      True
  0.20       31     0.0811    0.0747  1.0444      True
  0.40       32     0.1012    0.0640  1.0970      True
  0.70       37     0.1358    0.0528  1.2013      True
  1.00       38     0.1732    0.0428  1.3389      True
```

The kept group's genotypes are bit-identical down the column, so every pair's kinship over
them is too. The screen's demand on the kept group nevertheless falls by a factor of two.
Out of sample on a fresh canvas (`m = 45 056` = 30 000 @ 0.45 + 15 056 @ 0.20, `n = 240`,
seed 404): printed 32 → 33 → 35 → 38 → 40 as `fB` runs 0 → 1, kept bits identical
throughout.

This is stronger than §10 and it cuts the other way from what one might hope: sharing in
the discarded markers **does** count, but at roughly a quarter of the weight its
information deserves — `thr(whole)` climbing from 0.059 to 0.173 is the screen charging the
pair three times over for evidence it is only partly willing to spend.

## 16. `R` below 1 is real, and it kills the bound

Sharing confined to the *most informative* markers, `m = 50 000` (`screenexit.py below`):

```
     nA  MAF A  MAF B  thr(whole)      R
  32768   0.45   0.25      0.0589   0.9919
  32768   0.45   0.15      0.0571   0.9878
  24576   0.45   0.20      0.0536   0.9801
```

The screen accepts pairs at whole-map kinship 0.0536 against a printed cutoff of 0.0625.
A conservative bound — the whole content of the early-exit hypothesis — can only ever
reject pairs a full computation would accept. It cannot do this. Round 1's lone
sub-unity point (beta spectrum, 0.9980 ± 0.0003) was the same effect seen through a
spectrum that happens to put a pair's sharing slightly above the map's centre of
informativeness, and it is reproducible at four times the size.

## 17. The second condition is a knee in MAF, not a place in the ranking

`m = 32 768`; `u = 4096` unshared markers at MAF `a`, above 8192 @ MAF 0.15 that the pair
clones at the rung's `f`, junk @ 0.03 filling the rest. The block outranks the window at
every `a` below, by at least seven in-sample standard deviations
(`screenexit.py band`, three seeds):

```
   a    2pq      n = 150              n = 300              n = 600
0.34  0.4488  0.0666 0.0629 0.0712  0.0637 0.0627 0.0626  0.0628 0.0690 0.0648
0.38  0.4712  0.0663 0.0628 0.0688  0.0639 0.0646 0.0651  0.0628 0.0633 0.0637
0.40  0.4800  0.0660 0.0631 0.0691  0.0653 0.0634 0.0650  0.0639 0.0631 0.0626
0.41  0.4838  0.0641 0.0660 0.0687  0.0641 0.0655 0.0651  0.0631 0.0643 0.0631
0.42  0.4872  0.0655 0.0650 0.0754  0.0638 0.0655 0.0638  0.0654 0.0693 0.0786
0.43  0.4902  0.0768 0.0646 0.0815  0.0843 0.0801 0.0758  0.0866 0.0834 0.0928
0.44  0.4928  0.0791 0.0766 0.0823  0.1149 0.1101 0.0970  0.1040 0.1180 0.1115
0.45  0.4950  0.1032 0.0900 0.0894  0.1467 0.1370 0.1435  0.1512 0.1568 0.1417
0.46  0.4968  0.1062 0.0999 0.1025  0.1694 0.1664 0.1631  0.1591 0.1827 0.1555
0.48  0.4992  0.1494 0.1301 0.1302  0.1837 0.1760 0.1744  0.1805 0.1911 0.1946
```

Flat at the printed cutoff to `a = 0.42`, then a knee: the same 4096 markers, in the same
place in the ranking, differing by 3 % in `2pq`, have more than **doubled** what the pair
must score by `a = 0.45` and nearly tripled it by 0.48. The knee does not move between
`n = 150` and `n = 600`, so it is not a binomially smeared cut on the in-sample frequency
sitting somewhere else. Out of sample it reproduces at `m = 45 056` (above the budget) and
`m = 20 480` (well below it) with a window at MAF 0.22, junk at 0.06, `n = 260` and unused
seeds — flat at 0.36 and 0.40, and 0.11–0.26 from 0.43 up — so it is **not a budget
effect**.

**The knee moves with the cutoff.** In one geometry (window 16 384 @ MAF 0.20, `n = 300`)
the degree-2 knee sits at `a ≈ 0.435` while at degree 1 a block at `a = 0.34` already costs
1.25× and everything from 0.38 up is out of reach — with the no-block control landing on
0.1278/0.1390/0.1375 against the 0.1250 cutoff, as it must. So "MAF ≳ 0.42" is the knee's
position at degree 2, not a constant of the algorithm.

## 18. Rank, size and coding, each held apart

Same canvas, degree 2, window 8192 @ MAF 0.15 (`screenexit.py size`, `rank`, `coding`):

```
    u   |  a = 0.40 (below the knee)   |  a = 0.45 (above it)
    0   |  0.0626 0.0660 0.0646        |  0.0626 0.0660 0.0646
 1024   |  0.0682 0.0675 0.0637        |  0.0643 0.0679 0.0654
 2048   |  0.0673 0.0687 0.0705        |  0.0858 0.1083 0.0935
 4096   |  0.0653 0.0634 0.0650        |  0.1467 0.1370 0.1435
 8192   |  0.0643 0.0640 0.0636        |  0.1416   n/a    n/a
12288   |  0.0641 0.0635 0.0650        |    n/a    n/a    n/a
16384   |  0.0670 0.0651 0.0684        |    n/a    n/a    n/a
20480   |  0.0642 0.0655 0.0628        |    n/a    n/a    n/a
```

* **Rank is free.** Below the knee the block is free at *every* size: 20 480 unshared
  markers, five eighths of the map, every one of them ranked above the pair's entire
  sharing, leave the threshold on the cutoff. A prefix scan cannot produce that column.
* **Above the knee it is not.** 2048 markers already cost 1.5×; 8192 put the pair beyond
  reach at any kinship the canvas can build (0.15).
* **What rescues a pair is entering the band, not climbing the ranking.** With `u = 12 288`
  at MAF 0.45 fixed and the window's own MAF raised: 0.15 and 0.30 refused outright, 0.38
  refused-then-0.087, 0.42 and 0.45 back on the cutoff (0.0659–0.0702). MAF 0.38 outranks
  two thirds of the map and is still refused.
* **It is a band in MAF, not in allele coding.** Swapping hom-A1 and hom-A2 at a marker
  leaves het, HetHet and IBS0 — and so every kinship — bit-identical. Recoding the 45
  markers whose in-sample A1 is the major allele, or the whole 8192-marker block, leaves
  the printed count and the threshold unchanged to the last digit.

Also checked and negative: an IBS0 reweighting `phi - lambda*IBS0/min_het` needs
`lambda` = 0.005, 0.059, 0.120 and 0.160 at four points of the §17 sweep; a MAF-thresholded
subset `{MAF >= t}` on `bigish` counts 44–49 pairs at every `t` from 0 to 0.46 (§3.1's proof
again, measured); and `--noscreen`, new in 2.3.2 and advertised for `--related`, prints the
identical `Stages 1&2 (with 32768 SNPs): 36 pairs` line — it is a no-op on this path.

## 19. The band is not the deflation

Delete the band from `bigish` and the deflation does not move (`screenexit.py bigish`):

```
     c   m left  printed2  true2  printed1  true1  R(deg2)
  1.01    50000        36     47        18     21   1.0220
  0.45    47963        36     47        18     21   1.0221
  0.44    47303        37     47        18     21   1.0193
  0.42    45384        38     47        18     20   1.0154
  0.40    43345        37     47        18     22   1.0150
  0.38    41021        36     47        21     23   1.0152
  0.35    37342        45     47        18     22   1.0053
  0.30    30955        50     50        18     21   1.0009
```

`R` tracks `m` against the 32 768 budget and nothing else — the band-deleted maps land on
§4.3's prefix-truncation curve at the same `m` (1.0087 at 36 864, 1.0128 at 40 960) — and
the deflation vanishes only when `m` itself falls under the budget. So Round 1's deflation
and Round 3's band are **two different phenomena**. `bigish`'s uniformly related pairs
share inside the band in proportion to it, which is why the band never binds there and why
no instrument before this one saw it.

## 20. Where round 3 leaves the mechanism

Closed this round, on top of §5's and §12's lists:

* **a bound-based early exit**, in the early-accept form — §16 measures `R = 0.980`, and a
  conservative bound cannot accept what the full computation rejects;
* **a bound-based early exit**, in the early-reject form, and with it *every* mechanism
  whose state is a prefix of a sorted scan — §18's `a = 0.40` column is flat at the cutoff
  with 20 480 unshared markers ranked above the sharing;
* **any account of §9's second condition in terms of rank, contiguity, block size or
  informativeness** — the knee is 3 % wide in `2pq` and 60 % of the map wide in rank;
* **allele coding**, **an IBS0 reweighting of the estimator**, and **a MAF-thresholded
  subset** — §18;
* **`--noscreen` as a lever** — it does not change this line.

What a mechanism now has to do, all at once:

1. equal the whole-map kinship exactly when `m <= 32768` and the pair's sharing is spread
   (§1, §17's flat columns, §19's `c = 0.30` row);
2. read the pair's genotypes at discarded markers, and weight them at roughly a quarter of
   their information (§15);
3. accept below the printed cutoff when the sharing sits high in the spectrum (§16);
4. deflate above the budget as a function of `m/32768` and the spectrum, *independently* of
   the band (§19);
5. carry a condition that is nearly a step function of a marker's MAF, with the step at
   ≈ 0.42 at degree 2 and below 0.34 at degree 1, insensitive to rank, size below the step,
   and sample count (§17, §18).

Worth the next round's time, in order:

1. **Chase the knee, not the deflation.** It is a factor of three in kinship, it is one run
   per point, it reproduces out of sample at three map sizes and two degrees, and item 5 is
   a far tighter constraint on any candidate mechanism than the 2 % deflation ever was. The
   right next scan is the knee's position against the *cutoff* at fixed geometry — degree 1
   and 2 differ, so a third point (a fileset engineered so degree 3's exhaustive path still
   prints a screen line, or a `--degree 2` map whose ladder reaches far higher) would say
   whether the knee is a function of the cutoff alone.
2. **Ask what a marker at MAF 0.45 has that one at MAF 0.40 does not**, when the pair does
   not share it. `2pq` differs by 3 %, `2p²q²` by 6 %, `p²` by 27 %. Nothing smooth in `p`
   has a knee there, so the answer is likely a comparison — a marker being *the most*
   informative available, or crossing a fixed integer count — rather than a function.
3. **Keep the two phenomena apart.** §19 shows the deflation survives deleting the band and
   the band survives dropping below the budget; a mechanism that explains one and not the
   other is still progress, and conflating them cost this round two hours.
4. **Do not fit `R`.** Unchanged from §12.3.

## 21. Consequence for the implementation, unchanged

`related.rs` still estimates on the map's first `min(m, 32768)` markers. Right whenever
`m <= 32768`, right at degree 1 on `bigish` (18), wrong at degree 2 (50 against 36). The
blast radius is still one stdout line: `.kin0`'s rows come from the exhaustive re-estimate
that follows the screen and are byte-correct at every degree, and the 14 pairs the
reference's screen drops all sit below the 0.08839 reporting threshold.

Harness: **477 / 480 byte-identical, self-check 480 / 480** — unchanged.
