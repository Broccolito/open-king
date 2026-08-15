# 22 — the `--related` screening count

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
