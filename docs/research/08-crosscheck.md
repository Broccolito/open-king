# 08 — Cross-check of the KING-robust estimator against independent implementations

**Purpose.** Independent confirmation of the KING-robust kinship math, so that our clean-room
Rust reimplementation of KING 2.3.2 rests on the *published* estimator plus corroboration from
four unrelated codebases — never on KING's own source.

**Method.** Read the Manichaikul et al. (2010) paper equations directly (equation images
transcribed by eye), then check each third-party implementation's arithmetic against them.
Third-party code was read to extract *formulas only*; no code was transcribed, and no
third-party code will appear in our implementation.

**Verdict up front: there is no genuine disagreement about the math.** Every source that
implements KING-robust computes bit-for-bit the same estimator. All the real parity risk is in
*which estimator is selected*, *which variants enter the counts*, and *degenerate-input
behaviour* — not in the formula.

---

## 0. Licence posture of each source (read this before reusing anything)

| Source | Licence | Safe to read for math? | Safe to copy code? |
|---|---|---|---|
| Manichaikul et al. 2010, *Bioinformatics* 26(22):2867-73 (PMC3025716) | Journal article; math is uncopyrightable fact | **Yes — this is our authority** | n/a |
| PLINK 2.0 / plink-ng (`plink2_matrix_calc.cc`) | **GPL-3.0-or-later** | Yes (formulas are not copyrightable) | **No** |
| SNPRelate (`src/genKING.cpp`) | **GPL-3.0** | Yes | **No** |
| Hail (`hail/python/hail/methods/relatedness/king.py` + docs) | **MIT** | Yes | Would be attributable, but we don't need it |
| Illumina **akt** (`kin.cpp`) | **PolyForm Strict License 1.0.0** | ⚠️ see below | **No — and worse than GPL** |

> ⚠️ **akt licence flag (correction to the task brief).** The brief assumed akt is
> "permissively/GPL licensed". It is not. Illumina akt ships under **PolyForm Strict License
> 1.0.0**, which is *not* an open-source licence. Its copyright grant explicitly excludes
> *"distributing the software or making changes or new works based on the software"*, and
> permitted purposes are limited to noncommercial/personal/nonprofit use. Reading it as an
> independent check is defensible under fair use, and its contribution here turned out to be
> **exactly zero new information** — its expression is a literal transcription of the paper's
> equation (11) right-hand side. **Recommendation: strike akt from our reference set.** It is
> cited below only to record that it was checked and that it agrees; nothing in our
> implementation should trace to it.

MIT-licensing our own reimplementation is unaffected: the estimator is a published formula,
and formulas carry no copyright.

---

## 1. The authority: exact paper equations

Transcribed from the equation images in PMC3025716 (`btq559m2/m5/m7/m9/m11.jpg`, plus the
unnumbered identity `btq559um1.jpg`). These are the paper's own equations, verbatim in
structure.

### 1.1 Notation

| Symbol | Meaning |
|---|---|
| `X_m^(i)` | genotype dosage of individual *i* at SNP *m*, ∈ {0, 1, 2} |
| `M_ij` | number of SNPs with **non-missing genotypes in both** *i* and *j* |
| `N_Aa^(i)` | number of SNPs at which *i* is heterozygous — **counted over the M_ij pairwise-non-missing SNPs only** |
| `N_Aa^(j)` | same, for *j* |
| `N_Aa,Aa` | number of SNPs at which **both** are heterozygous |
| `N_AA,aa` | number of SNPs at which the two are **opposite homozygotes** (this is IBS0) |
| `Ĥ_ij` | the pair's heterozygosity estimate — the *only* thing that differs between the estimators |

Paper prose, §2.2, defining these counts (verbatim):

> "…where N Aa , Aa , N Aa ( i ) and N Aa ( j ) are the total numbers of SNPs at which both
> individuals of the pair are heterozygous, and the total number of heterozygotes for the i -th
> and j -th individual, respectively, **excluding those SNPs with missing genotypes in either
> individual of the pair**."

That sentence settles the pairwise-missingness question at the source. `N_Aa^(i)` is **not**
individual *i*'s global heterozygote count; it is recomputed per pair.

### 1.2 The genetic-distance identity (unnumbered, §2.2)

```
(X^(i) − X^(j))²  =  4·I_AA,aa  −  2·I_Aa,Aa  +  I_Aa^(i)  +  I_Aa^(j)
```

Summed over the pairwise-non-missing SNPs:

```
Σ_m (X_m^(i) − X_m^(j))²  =  4·N_AA,aa  −  2·N_Aa,Aa  +  N_Aa^(i)  +  N_Aa^(j)
```

This is the bridge between the "sum of squared dosage differences" form and the "genotype
counts" form. Every implementation below uses one side or the other.

### 1.3 Equation (5) — the master form

```
φ̂_ij = 1/2  −  [ Σ_m (X_m^(i) − X_m^(j))² ] / ( 4·Ĥ_ij )
```

### 1.4 Equation (7) — the master form in genotype counts

```
φ̂_ij = (N_Aa,Aa − 2·N_AA,aa) / (2·Ĥ_ij)  +  1/2  −  (N_Aa^(i) + N_Aa^(j)) / (4·Ĥ_ij)
```

### 1.5 Equation (9) — KING-robust **WITHIN**-family

Plug `Ĥ_ij = (N_Aa^(i) + N_Aa^(j)) / 2` (the *average* het count) into (5)/(7):

```
φ̂_ij = 1/2 − (1/2)·[ Σ_m (X_m^(i) − X_m^(j))² ] / ( N_Aa^(i) + N_Aa^(j) )

     = ( N_Aa,Aa − 2·N_AA,aa ) / ( N_Aa^(i) + N_Aa^(j) )        ← denominator is the SUM
```

The right-hand equality is printed in the paper itself; the `+1/2` and `−(N_i+N_j)/(4Ĥ)` terms
of (7) cancel exactly when `Ĥ_ij` is the average. **The within-family form has no additive
1/2.**

### 1.6 Equation (11) — KING-robust **BETWEEN**-family

Paper prose immediately preceding: *"we consider the smaller of the observed heterozygosity
rates, min( N Aa ( i ) / M ij , N Aa ( j ) / M ij ), as an alternative to E (2 P (1 − P )).
Without loss of generality, suppose the i -th individual has lower heterozygosity than the j
-th individual."*

So `Ĥ_ij = N_Aa^(i)` where **i is the individual with the SMALLER het count**:

```
φ̂_ij = 1/2 − (1/4)·[ Σ_m (X_m^(i) − X_m^(j))² ] / N_Aa^(i)

     = ( N_Aa,Aa − 2·N_AA,aa ) / ( 2·N_Aa^(i) )
       + 1/2
       − ( N_Aa^(i) + N_Aa^(j) ) / ( 4·N_Aa^(i) )                ← denominator is the MIN
```

with `N_Aa^(i) = min(N_Aa^(i), N_Aa^(j))`.

**The extra `+1/2 − (N_i+N_j)/(4·min)` correction term is essential and is the single most
commonly dropped piece of this formula.** It vanishes only when `N_Aa^(i) = N_Aa^(j)`, so a
naive implementation that writes `(N_Aa,Aa − 2·N_AA,aa)/(2·min)` will pass an MZ-twin test and
an equal-heterozygosity test and then be silently wrong on every real pair.

### 1.7 Which is used where — the paper's own rule

> "We use estimator (Equation 9) for **within-family** relationship checking and estimator
> (Equation 11) for **between-family** relationship checking, naming this combined approach
> KING-robust."

And on bounds:

> "The estimator above is no larger than the estimator in Equation (9), and **both estimators
> are bounded above by 0.5**."

Note that 0.5 is an *analytic* upper bound under the model, not a clamp — see §5.4.

### 1.8 Supporting equations (needed only if we also implement KING-homo / `--homog`)

Equation (2), IBD0 proportion:
```
π̂_0,ij = N_AA,aa / Σ_m 2·p̂_m²·(1 − p̂_m)²
```

Equation (5) with the homogeneous-population plug-in `Ĥ_ij = 2·Σ_m p̂_m(1 − p̂_m)`:
```
φ̂_ij = 1/2 − [ Σ_m (X_m^(i) − X_m^(j))² ] / ( 8·Σ_m p̂_m(1 − p̂_m) )
```
Allele frequency `p̂_m` from equation (3) is estimated from genotype frequencies across the
**entire sample**, not the pair.

---

## 2. PLINK 2.0

### 2.1 What the documentation says (verbatim)

From <https://www.cog-genomics.org/plink/2.0/distance#make_king>:

- "**--make-king** writes KING-robust coefficients in matrix form to `plink2.king[.zst]` or `plink2.king.bin`"
- "**--make-king-table** writes them in table form to `plink2.kin0[.zst]`"
- "KING kinship coefficients are scaled such that duplicate samples have kinship 0.5, not 1. First-degree relations (parent-child, full siblings) correspond to ~0.25, second-degree relations correspond to ~0.125, etc."
- "Only autosomes are included in this computation."
- "**Pedigree information is currently ignored; the between-family estimator is used for all pairs.**"
- "For multiallelic variants, REF allele counts are used."
- "KING-robust … doesn't require MAFs at all" but "KING-robust underestimates kinship when the parents are very different populations."
- "**--king-table-filter** causes only kinship coefficients ≥ the given threshold to be reported."
- Modifiers: `counts` (counts instead of frequencies), `rel-check` (same-FID pairs only), `concordance-check` (same-FID-and-IID pairs only), `cols=`.
- On the relationship to KING itself, the docs only say: "See also the original KING software package, which has some useful two-step workflows directly built in…". **No claim of numeric or byte identity with KING's `--kinship` is made anywhere.**

### 2.2 What the source actually computes

`plink2_matrix_calc.cc` (GPL-3+), `ComputeKinship()` around lines 1565–1573, and the same
expression inlined in the `.kin0` writer at lines 2292–2298. Restated in our notation (formula
only — no code reproduced):

The pair accumulators are, over SNPs where **neither** call is missing:
`ibs0` (opposite homs), `hethet`, `het1hom2` (sample 1 het, sample 2 hom), `het2hom1`, `homhom`.

```
smaller_het  =  hethet + min(het1hom2, het2hom1)
KINSHIP      =  0.5  −  ( 4·ibs0 + het1hom2 + het2hom1 ) / ( 4·smaller_het )
```

Two identities make this exactly paper equation (11):

- `N_Aa^(1) = hethet + het1hom2` and `N_Aa^(2) = hethet + het2hom1`, therefore
  `smaller_het = min(N_Aa^(1), N_Aa^(2))` — **min denominator, confirmed**;
- `4·ibs0 + het1hom2 + het2hom1 = Σ_m (X^(1) − X^(2))²` by the §1.2 identity — so this is
  literally `0.5 − Σ(ΔX)²/(4·min)`, i.e. equation (11) left-hand side.

**Missing data:** pairwise. `het1hom2` requires sample 2 to be a called hom, so the het counts
can only accumulate on doubly-called sites — structurally identical to the paper's definition.

**Within-family:** not implemented. Equation (11) is applied to every pair including same-FID
pairs. There is no equation (9) code path.

### 2.3 `.kin0` output shape

Header emitted by `AppendKingTableHeader()` (lines ~1611–1653), tab-separated, `#`-prefixed:

```
#FID1  IID1  [SID1]  FID2  IID2  [SID2]  NSNP  HETHET  IBS0  HET1_HOM2  HET2_HOM1  IBS  KINSHIP
```

- FID columns present only when FIDs are in play; SID columns only when SIDs are.
- **The docs page lists `ID1`/`ID2`; the current source emits `IID1`/`IID2`.** The source
  carries an explicit comment that the name was `ID1` before alpha 3 and that "the header line
  still doesn't perfectly match KING due to e.g. capitalization" — an outright statement that
  plink2 does **not** target byte parity with KING.
- `NSNP = hethet + het1hom2 + het2hom1 + homhom` (line ~2316) — i.e. `M_ij`. Note `homhom`
  includes the opposite-hom sites, so `ibs0 ⊆ homhom`.
- `IBS` = Hamming distance = `2·ibs0 + het1hom2 + het2hom1` (line ~2352).
- HETHET/IBS0/HET1_HOM2/HET2_HOM1 are **proportions of NSNP by default**; the `counts`
  modifier switches them to raw counts.
- Numbers are written with plink2's `dtoa_g` (shortest-round-trip `%g`-style), **not** a fixed
  `%.4f`. Formatting is therefore not comparable to KING.

### 2.4 Matrix output and ordering

- `.king` is lower-triangular text excluding the diagonal by default; `square`, `square0`,
  `triangle`, `bin`, `bin4` modifiers as for `--make-rel`. `.king.id` carries the sample IDs.
- Source comment at lines ~2011–2013: results "are always reported in **lower-triangular
  order, rather than KING's upper-triangular order**". A direct, citable statement about
  KING's own pair emission order.

### 2.5 Verdict

PLINK 2.0 computes **exactly paper equation (11)**, pairwise-missing, min-denominator,
between-family-only. It **claims no identity with KING's `--kinship`**, and documents two
deliberate departures: pedigree is ignored, and output ordering/naming differ.

---

## 3. Hail — `hl.king`

Docs (<https://hail.is/docs/0.2/methods/relatedness.html>) and the MIT-licensed docstring in
`hail/python/hail/methods/relatedness/king.py` both state:

```
φ̂_ij^between = 1/2 + ( 2·N^(Aa,Aa)_ij − 4·N^(AA,aa)_ij − N^(Aa)_i − N^(Aa)_j )
                     / ( 4 · min( N^(Aa)_i , N^(Aa)_j ) )
```

- **Denominator: min.** Docstring: *"The estimator replaces the average count of heterozygous
  genotypes with the minimum count of heterozygous genotypes: (N^Aa_i + N^Aa_j)/2 ⇝
  min(N^Aa_i, N^Aa_j)"* — the clearest available statement of the *why*, and it matches the
  paper's `Ĥ_ij` substitution exactly.
- **Missing data: pairwise.** *"The three counts above, N^Aa, N^Aa,Aa, and N^AA,aa, exclude
  variants where one or both individuals have missing genotypes."*
- **Within-family: not implemented.** *"This function, king(), only implements the
  'between-family' estimator, φ̂_ij^between."*
- Genotype score `X_i,s` from `n_alt_alleles()` ∈ {0,1,2}.

This is the same estimator as PLINK 2.0, written in the equation-(7) expanded form rather than
the equation-(5) squared-difference form. Expanding `N^(Aa)_i + N^(Aa)_j = 2·hethet +
het1hom2 + het2hom1` turns Hail's numerator into `−(4·ibs0 + het1hom2 + het2hom1)`, i.e.
plink2's numerator with the sign folded into the leading `1/2`. Algebraically identical
(verified numerically in §6).

No claim of parity with the KING binary.

---

## 4. SNPRelate — `snpgdsIBDKING`

The most informative of the four, because it is the only third-party implementation that
implements **both** paper estimators and the within/between *selection rule*.

### 4.1 Counting (`src/genKING.cpp`, GPL-3, scalar path ~lines 404–424)

Per pair, over a `mask` defined as "called in individual 1 AND called in individual 2":

| accumulator | meaning (documented in the struct at lines 276–280) |
|---|---|
| `nLoci` | popcount of the pairwise mask = `M_ij` |
| `IBS0` | opposite homozygotes = `N_AA,aa` |
| `SumSq` | struct comment: `\sum_m (g_m^{(i)} - g_m^{(j)})^2` = `Σ(ΔX)²` |
| `N1_Aa` | "the number of hetet loci for the first individual" |
| `N2_Aa` | second individual |

`N1_Aa` and `N2_Aa` are each ANDed with the pairwise `mask` before popcount → **pairwise
missingness, explicitly**.

`SumSq` is accumulated as `popcount(het-xor-term) + 4·popcount(ibs0)` — the §1.2 identity again
(`het1hom2 + het2hom1 + 4·ibs0`).

### 4.2 The estimator selection (lines ~630–638 matrix path, ~655–664 vector path)

```
kinship = (same non-NA family id)
            ?  0.5 − SumSq / ( 2 · (N1_Aa + N2_Aa) )        ← paper eq (9), SUM
            :  0.5 − SumSq / ( 4 · min(N1_Aa, N2_Aa) )      ← paper eq (11), MIN
```

These are *literally* the equation-(5) master form with `Ĥ_ij` set to the average and to the
min respectively — the cleanest independent confirmation available that the two KING-robust
estimators differ **only** in `Ĥ_ij`.

The family-id switch comes from the R wrapper's `family.id=` argument (`R/IBD.R`):
`family.id=NULL` (the default) fills all ids with `NA`, and the C++ requires
`f1 == f2 && f1 != NA` for the within-family branch. So **SNPRelate's default behaviour equals
PLINK 2.0's** (between-family everywhere), but supplying `family.id` reproduces KING's own
within/between split.

### 4.3 Other SNPRelate details worth mirroring in tests

- Diagonal is written as `kinship = 0.5`, `IBS0 = 0` — self-kinship is *asserted*, not computed.
- The `IBS0` output is a **proportion**: `IBS0 / nLoci`, or `NaN` when `nLoci == 0`.
- Any non-finite kinship is coerced to `NaN` (`R_FINITE` guard) — so a zero denominator surfaces
  as `NaN`, not `-Inf`. **This differs from PLINK 2.0** (see §5.3).
- `type="KING-homo"` computes `theta = 0.5 − SumSq/(8·Σp(1−p))` and `k0 = IBS0/(2·Σp²(1−p)²)`,
  `k1 = 2 − 2·k0 − 4·theta` — exactly paper equations (5)-with-homo-plug-in and (2). Independent
  confirmation of §1.8.
- Docs describe the two types as: "KING-robust" = "robust relationship inference within or
  across families in the presence of population substructure"; "KING-homo" = "relationship
  inference in a homogeneous population".

---

## 5. Illumina akt — `kin.cpp` (recorded for completeness; **do not use**, see §0)

`Kinship::estimateKinship`, `method == 1`, lines ~195–207. The kinship line is a direct
transcription of the right-hand side of paper equation (11):

```
minhet = min(Nhet_1, Nhet_2)
ks     = (Nhet_12 − 2·ibd0) / (2·minhet)  +  0.5  −  0.25·(Nhet_1 + Nhet_2) / minhet
```

with `ibd0 = N_AA,aa`, `Nhet_12 = N_Aa,Aa`. The source's own inline comments annotate the
variables as `//NAa^i`, `//NAa^j`, `//NAa,Aa`, confirming the mapping.

Missing data is pairwise: the het counts are masked by `(missing_1 | missing_2).flip()`, and the
code notes that the hethet count needs no mask because het-het already implies both called.

`ibd0` is a `float&`, so the divisions are floating point — no integer-division trap despite the
integer-looking counts.

**Contributes no independent information** (it is the paper's equation verbatim) and is under a
non-open-source licence. Excluded from the reference set.

---

## 6. Numerical verification that all stated forms are one formula

Implemented all nine published/observed expressions independently in Python and evaluated them
on 20,000 random genotype pairs (20–60 SNPs each, alleles drawn uniformly from {0, 1, 2,
missing}), skipping only pairs with a zero denominator.

**Between-family group** — all agreed to < 1e-12 on every trial:
1. paper eq (11) LHS: `0.5 − 0.25·Σ(ΔX)²/min`
2. paper eq (11) RHS: `(N_AaAa − 2N_AAaa)/(2·min) + 0.5 − 0.25(N_i+N_j)/min`
3. PLINK 2.0: `0.5 − (4·ibs0 + het1hom2 + het2hom1)/(4·(hethet + min(het1hom2, het2hom1)))`
4. Hail: `0.5 + (2·hethet − 4·ibs0 − N_i − N_j)/(4·min)`
5. akt: `(hethet − 2·ibs0)/(2·min) + 0.5 − 0.25(N_i+N_j)/min`
6. SNPRelate: `0.5 − SumSq/(4·min)`

**Within-family group** — all agreed to < 1e-12:
1. paper eq (9) LHS: `0.5 − 0.5·Σ(ΔX)²/(N_i+N_j)`
2. paper eq (9) RHS: `(N_AaAa − 2N_AAaa)/(N_i+N_j)`
3. SNPRelate: `0.5 − SumSq/(2·(N_i+N_j))`

Both bridging identities held on every trial:
- `Σ(ΔX)² == 4·ibs0 + het1hom2 + het2hom1`
- `min(N_i, N_j) == hethet + min(het1hom2, het2hom1)`

**Mismatches: 0 / 20,000.**

### Golden vectors for our unit tests

Genotypes as 0/1/2, `.` = missing. Counts are over pairwise-called sites.

| case | g1 | g2 | M | hethet | het1hom2 | het2hom1 | ibs0 | N_i | N_j | **between (min)** | **within (sum)** |
|---|---|---|---|---|---|---|---|---|---|---|---|
| MZ / duplicate | `1111002020` | `1111002020` | 10 | 4 | 0 | 0 | 0 | 4 | 4 | **0.5** | **0.5** |
| unequal het, no IBS0 | `1111102020` | `1102002020` | 10 | 2 | 3 | 0 | 0 | 5 | 2 | **0.125** | **0.285714285714…** |
| IBS0 + missing | `021120.120` | `2011201.00` | 8 | 2 | 1 | 0 | 3 | 3 | 2 | **−1.125** | **−0.8** |

Row 2 is the discriminating test: it is the case where dropping the `+1/2 − (N_i+N_j)/(4·min)`
correction, or swapping min for sum, produces a *different but plausible-looking* answer
(0.2 for the naive `(hethet−2·ibs0)/(2·min)` form vs the correct 0.125). Any implementation
that only tests MZ twins will not catch that.

---

## 7. CONSENSUS FORMULA SET

Unanimous across the paper and all four implementations. This is what we implement.

### 7.1 Per-pair accumulation (single pass, integers only)

Over autosomal variants that pass the variant filter, restricted to sites where **both**
individuals have a called genotype:

```
hethet    = # sites both heterozygous                            (= N_Aa,Aa)
het1hom2  = # sites i heterozygous, j homozygous
het2hom1  = # sites i homozygous,   j heterozygous
ibs0      = # sites opposite homozygotes                         (= N_AA,aa)
homhom    = # sites both homozygous                              (⊇ ibs0)

M_ij      = hethet + het1hom2 + het2hom1 + homhom                (= N_SNP)
N_i       = hethet + het1hom2                                    (= N_Aa^(i), pairwise)
N_j       = hethet + het2hom1                                    (= N_Aa^(j), pairwise)
```

### 7.2 Between-family estimator — paper eq (11), used for cross-FID pairs → `.kin0`

```
min_het = hethet + min(het1hom2, het2hom1)          ( = min(N_i, N_j) )

φ̂ = 0.5 − ( 4·ibs0 + het1hom2 + het2hom1 ) / ( 4 · min_het )
```

### 7.3 Within-family estimator — paper eq (9), used for same-FID pairs → `.kin`

```
φ̂ = ( hethet − 2·ibs0 ) / ( N_i + N_j )
```

### 7.4 Implementation notes that follow from the consensus

- **Accumulate integers, divide once.** Both forms are an exact-integer numerator over an
  exact-integer denominator with a single floating-point division. This eliminates
  summation-order sensitivity entirely and makes results bit-reproducible regardless of
  SIMD/threading. Do **not** accumulate `Σ(ΔX)²` as a float.
- **The `min` picks the individual, not the term.** `min_het` is `min(N_i, N_j)`; computing it
  as `hethet + min(het1hom2, het2hom1)` is equivalent and avoids materialising `N_i`/`N_j`.
- **Never use a global per-sample heterozygote count.** `N_i` must be recomputed per pair under
  the pairwise-called mask. This is the single subtlest requirement in the whole estimator and
  the most likely source of a near-miss that only shows up on samples with differing call rates.
- **KING prints `%.4lf`.** With four decimal places, last-ULP differences between algebraically
  equivalent orderings are invisible except at exact `.00005` rounding ties, so we do not need
  to guess KING's exact operation order — but we should keep the single-division form so that
  ties are at least deterministic on our side.

### 7.5 Corroboration from the KING 2.3.2 binary itself

A light `strings` pass over `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(observing embedded constants — facts about behaviour, not source) corroborates the above:

- Column tokens present as distinct strings: `N_SNP`, `HetHet`, `IBS0`, `Kinship`, `HomIBS0`,
  `PropIBD`, `InfType`, `NHetHet`, `N_IBS0`, `KinshipX` — the row format is assembled at
  runtime, so no single row-format string exists to read off.
- The embedded R plotting script exposes the exact degree thresholds KING uses internally:
  `0.04419`, `0.08839`, `0.17678`, `0.35355`. These are `2^-4.5`, `2^-3.5`, `2^-2.5`, `2^-1.5`
  — i.e. the cutoffs are `2^-(d + 3/2)`, matching Table 1 of the paper (`>0.354`, `[0.177,
  0.354]`, `[0.0884, 0.177]`, `[0.0442, 0.0884]` per the KING manual) at full precision.
  **Use the exact powers of two, not the rounded manual values.**
- `"1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf"` — confirms IBS0 is
  the PO-vs-FS discriminator and that KING formats it with 4 decimals.
- KING's `.kin` columns: `FID ID1 ID2 N_SNP Z0 Phi HetHet IBS0 Kinship Error`;
  `.kin0` columns: `FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship`. Note KING uses `ID1`/`ID2`,
  **not** plink2's `IID1`/`IID2`, and `N_SNP`, **not** plink2's `NSNP`.

---

## 8. Where the sources DISAGREE — the parity-test targets

None of these are disagreements about the estimator's algebra. They are disagreements about
*application*, and they are precisely where our KING parity work must concentrate.

### 8.1 ⚠️ HIGH — Which estimator is applied to which pair

| Source | same-FID pairs | cross-FID pairs |
|---|---|---|
| **KING (paper §2.3 + `.kin`/`.kin0` split)** | **eq (9), sum denominator** | eq (11), min denominator |
| SNPRelate with `family.id=` | eq (9), sum | eq (11), min |
| SNPRelate default (`family.id=NULL`) | eq (11), min | eq (11), min |
| **PLINK 2.0** | **eq (11), min** — "Pedigree information is currently ignored" | eq (11), min |
| **Hail** | **eq (11), min** — "only implements the 'between-family' estimator" | eq (11), min |
| akt | eq (11), min | eq (11), min |

**Consequence for us:** if we validate our `--kinship` only against PLINK 2.0, every
within-family pair will appear to disagree, and we will "fix" the wrong thing. Our `.kin`
writer must use eq (9); our `.kin0` writer must use eq (11). Validate `.kin0` against plink2,
and `.kin` against SNPRelate-with-`family.id` (or against the KING binary directly).

Sanity check on magnitude: the two forms diverge whenever `N_i ≠ N_j` (golden-vector row 2:
0.125 vs 0.2857 — a factor of 2.3). Divergence is zero when het counts are equal, so
same-ancestry, same-call-rate pairs will look deceptively fine.

### 8.2 ⚠️ HIGH — Which variants enter the counts

No two sources agree, and **the paper does not specify this at all** — it is a program-level
policy, not part of the estimator:

- **PLINK 2.0:** "Only autosomes are included in this computation." No MAF filter, no call-rate
  filter, no LD pruning by default. "For multiallelic variants, REF allele counts are used."
- **SNPRelate / Hail:** whatever variant set the caller passes in; no implicit filtering.
- **KING 2.3.2:** the binary contains
  `"%d autosome SNPs with MAF>%.3lf and call rate>%d%% are used."`, plus separate call-rate
  filter messages for autosomes, chrX and chrY. **KING pre-filters variants by MAF and call
  rate by default.** The exact default thresholds must be read off a real run — they are not in
  the paper and are the single largest threat to `N_SNP` parity (and therefore to every count
  column, though *not* to the kinship value's formula).

**Action:** treat "which variants KING kept" as a separate, independently-verified stage. Get
`N_SNP` matching before comparing `Kinship`. A `Kinship` mismatch with a matching `N_SNP` is a
formula bug; a `Kinship` mismatch with a mismatched `N_SNP` is a filtering bug, and they need
completely different fixes.

### 8.3 ⚠️ MEDIUM — Zero-denominator behaviour

Reached when `min(N_i, N_j) = 0` (a pair sharing no called heterozygous site in the smaller
sample) — realistic on tiny variant sets, on inbred samples, or on badly-overlapping call sets.

| Source | Result |
|---|---|
| **PLINK 2.0** | `−Inf`. The source carries an explicit dated note (18 Nov 2017) that "kinship_coeff can be -inf when smaller_het_ct is zero", and deliberately **does not** filter those rows out via `--king-table-filter`. |
| **SNPRelate** | `NaN` — a `R_FINITE` guard rewrites any non-finite value. |
| **Hail** | unspecified in docs. |
| **KING 2.3.2** | **UNKNOWN — must be observed.** |

**Action:** construct a fixture that forces `min_het = 0` and record what KING actually prints
(`-inf`? `nan`? `-1.#IND`? a skipped row?). This is a genuine open question, not resolvable
from any documentation.

### 8.4 ⚠️ MEDIUM — Clamping

**No implementation clamps.** The paper says both estimators "are bounded above by 0.5", but
that is an asymptotic property of the estimator, not an enforced range — and the *lower* bound
is unbounded (our golden vector row 3 yields −1.125, and the paper's equation (10) explicitly
discusses large negative values as a *signal* of population heterogeneity: "an extreme negative
value indicates the pair of individuals is drawn from two distinct populations").

**Action:** do not clamp. Verify against KING that large negative kinships pass through
unmodified rather than being floored at 0 or −1.

### 8.5 ⚠️ MEDIUM — Pair ordering and output layout

- PLINK 2.0 source, verbatim: results "are always reported in **lower-triangular order, rather
  than KING's upper-triangular order**". So KING emits pairs in upper-triangular order —
  a citable third-party statement about KING's behaviour, but confirm it against the binary.
- PLINK 2.0's `.king` matrix omits the diagonal by default; SNPRelate writes `0.5` on the
  diagonal. KING's matrix conventions must be observed separately.

### 8.6 LOW — Column naming and number formatting

Not a disagreement about math, but fatal to byte parity:

| | KING 2.3.2 | PLINK 2.0 |
|---|---|---|
| within-family file | `.kin` | *(not produced)* |
| between-family file | `.kin0` | `.kin0` |
| ID columns | `FID1 ID1 FID2 ID2` | `#FID1 IID1 FID2 IID2` |
| SNP count column | `N_SNP` | `NSNP` |
| het-het column | `HetHet` | `HETHET` |
| kinship column | `Kinship` | `KINSHIP` |
| count columns | counts | **proportions by default** (`counts` modifier to switch) |
| number format | `%.4lf` style fixed-decimal | `dtoa_g` shortest round-trip |

plink2's own source concedes "the header line still doesn't perfectly match KING due to e.g.
capitalization". **PLINK 2.0 is a correctness oracle for the kinship *value* only; it is
useless as a format oracle.** Format parity must come from the binary.

### 8.7 Settled — no disagreement (do not re-litigate)

- **Denominator, between-family: `min(N_Aa^i, N_Aa^j)`.** Unanimous: paper eq (11), plink2,
  Hail, SNPRelate, akt. The `sum` denominator belongs to eq (9) *only*.
- **Missing data: pairwise, always.** Unanimous, and stated outright in the paper
  ("excluding those SNPs with missing genotypes in either individual of the pair"), in Hail's
  docstring, in plink2's `.kin0` NSNP definition ("neither call missing"), in SNPRelate's mask,
  and in akt's mask.
- **`N_Aa^(i)` is pair-specific, never a global per-sample het count.** Unanimous.
- **The `+1/2 − (N_i+N_j)/(4·min)` correction is part of eq (11).** Unanimous.
- **Scaling: duplicates → 0.5, not 1.** Unanimous.

---

## 9. Recommended parity-test plan (falls out of §8)

1. **Stage the comparison.** Match `N_SNP` first (variant filtering, §8.2), then the four count
   columns, only then `Kinship`. A kinship-only comparison cannot distinguish a filtering bug
   from a formula bug.
2. **Test both estimators separately.** Build a fixture with ≥2 individuals sharing an FID and
   ≥2 in different FIDs, so `.kin` (eq 9) and `.kin0` (eq 11) are exercised in one run.
3. **Use golden-vector row 2 (§6) as the regression guard** against min/sum confusion and
   against dropping the correction term. MZ-twin tests are worthless here — every wrong variant
   returns 0.5.
4. **Force the degenerate cases:** `min_het = 0`, `M_ij = 0`, a pair with zero overlap in called
   sites, and a pair with an extreme negative kinship. Record KING's literal output bytes for
   each (§8.3, §8.4).
5. **Cross-validate the between-family value against `plink2 --make-king-table counts`** on the
   *same* variant set (feed plink2 the exact post-filter variant list so §8.2 cannot confound
   it). Expect agreement to full double precision on `KINSHIP`, and exact agreement on
   `HETHET`/`IBS0`/`HET1_HOM2`/`HET2_HOM1` counts. Expect the header, column names, ordering
   and number formatting to differ (§8.5, §8.6) — do not chase those.
6. **Cross-validate the within-family value against `snpgdsIBDKING(..., family.id=)`**, which is
   the only readable implementation of equation (9).

---

## Appendix A — artefacts on disk

Downloaded reference material (read for math only; none of it is to be copied into the
implementation):

```
…/scratchpad/research/src/plink2_matrix_calc.cc      PLINK 2.0 / plink-ng, GPL-3+
…/scratchpad/research/src/snprelate_genKING.cpp      SNPRelate, GPL-3
…/scratchpad/research/src/snprelate_IBD.R            SNPRelate R wrapper, GPL-3
…/scratchpad/research/src/akt_kin.cpp                Illumina akt, PolyForm Strict 1.0.0 — DO NOT USE
…/scratchpad/research/src/pmc3025716.html            Manichaikul et al. 2010, PMC copy
…/scratchpad/research/src/eq/btq559m{1,2,5,7,9,10,11}.jpg   paper equation images (+ *_big.png upscales)
…/scratchpad/research/src/eq/btq559um1.jpg           the unnumbered §2.2 identity
```

## Appendix B — sources consulted

- Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen W-M. "Robust relationship
  inference in genome-wide association studies." *Bioinformatics* 2010;26(22):2867-73.
  <https://pmc.ncbi.nlm.nih.gov/articles/PMC3025716/>
- PLINK 2.0 documentation — distance/relationship matrices:
  <https://www.cog-genomics.org/plink/2.0/distance#make_king>
- PLINK 2.0 documentation — file formats (`.king`, `.king.id`, `.kin0`):
  <https://www.cog-genomics.org/plink/2.0/formats>
- plink-ng source: <https://github.com/chrchang/plink-ng> (`2.0/plink2_matrix_calc.cc`)
- Hail documentation — relatedness methods: <https://hail.is/docs/0.2/methods/relatedness.html>
- Hail source: `hail/python/hail/methods/relatedness/king.py`
- SNPRelate documentation — `snpgdsIBDKING`: <https://rdrr.io/bioc/SNPRelate/man/snpgdsIBDKING.html>
- SNPRelate source: <https://github.com/zhengxwen/SNPRelate> (`src/genKING.cpp`, `R/IBD.R`)
- Illumina akt: <https://github.com/Illumina/akt> (`kin.cpp`, `LICENSE`)
- KING manual: <https://www.kingrelatedness.com/manual.shtml>
- KING 2.3.2 binary (embedded constants only):
  `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
