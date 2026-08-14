# 08 — Cross-check of the KING-robust estimator against independent implementations

**Purpose.** Independent confirmation of the KING-robust kinship math, so that our clean-room
Rust reimplementation of KING 2.3.2 rests on the *published* estimator plus corroboration from
several unrelated codebases — never on KING's own source.

**Method.** Read the Manichaikul et al. (2010) paper prose and equations directly, then check
each third-party implementation's arithmetic against them. Third-party code was read to extract
*formulas only*; no code was transcribed, and no third-party code will appear in our
implementation. Every quotation below was re-verified against the on-disk artefact at the time
of writing (see Appendix A) — nothing here is from memory.

**Verdict up front: there is no genuine disagreement about the math.** Every source that
implements KING-robust computes the same estimator, and a 20 000-trial numeric cross-check of
all nine published/observed algebraic forms found **zero** mismatches. All the real parity risk
is in *which estimator is selected*, *which variants enter the counts*, *what the output columns
mean* (proportions vs counts), and *degenerate-input behaviour* — not in the formula.

> **Revision note.** This document supersedes an earlier draft that contained two substantive
> errors, both now corrected and re-derived from scratch:
> 1. the §6 golden-vector table was internally inconsistent (row 3's stated counts did not match
>    its own genotype strings) and its "naive form" counter-example value was wrong;
> 2. §8.6 claimed KING emits **counts** in `HetHet`/`IBS0`. It does not — the KING manual
>    defines both as **proportions**. This inverts the recommended plink2 cross-validation
>    recipe. See §8.6.

---

## 0. Licence posture of each source (read this before reusing anything)

| Source | Licence (verified) | Safe to read for math? | Safe to copy code? |
|---|---|---|---|
| Manichaikul et al. 2010, *Bioinformatics* 26(22):2867-73 (PMC3025716) | Journal article; math is uncopyrightable fact | **Yes — this is our authority** | n/a |
| PLINK 2.0 / plink-ng (`plink2_matrix_calc.cc`) | **GPL-3.0-or-later** — file header: "either version 3 of the License, or (at your option) any later version" | Yes (formulas are not copyrightable) | **No** |
| SNPRelate (`src/genKING.cpp`) | **GPL-3.0** — file header: "under the terms of the GNU General Public License Version 3" | Yes | **No** |
| Hail (docs + `hail/python/hail/methods/relatedness/king.py`) | **MIT** | Yes | Would be attributable, but we don't need it |
| Illumina **akt** (`kin.cpp`) | **PolyForm Strict License 1.0.0** | ⚠️ see below | **No — and worse than GPL** |

> ⚠️ **akt licence flag (correction to the task brief).** The brief assumed akt is
> "permissively/GPL licensed". It is not. Illumina akt ships under **PolyForm Strict License
> 1.0.0**, verified by fetching `Illumina/akt/LICENSE`. Its grant covers everything
> *"other than distributing the software or making changes or new works based on the software"*,
> and permitted purposes are restricted to noncommercial use. Reading it as an independent check
> is defensible, and its contribution here turned out to be **exactly zero new information** —
> its expression is a literal transcription of the paper's equation (11) right-hand side.
> **Recommendation: strike akt from our reference set.** It is recorded below only to note that
> it was checked and that it agrees; nothing in our implementation should trace to it.

MIT-licensing our own reimplementation is unaffected: the estimator is a published formula, and
formulas carry no copyright.

---

## 1. The authority: the paper

### 1.1 Notation

| Symbol | Meaning |
|---|---|
| `X_m^(i)` | genotype dosage of individual *i* at SNP *m*, ∈ {0, 1, 2} |
| `M_ij` | number of SNPs with **non-missing genotypes in both** *i* and *j* |
| `N_Aa^(i)` | number of SNPs at which *i* is heterozygous — **counted over the M_ij pairwise-non-missing SNPs only** |
| `N_Aa^(j)` | same, for *j* |
| `N_Aa,Aa` | number of SNPs at which **both** are heterozygous |
| `N_AA,aa` | number of SNPs at which the two are **opposite homozygotes** (this is IBS0) |
| `Ĥ_ij` | the pair's heterozygosity estimate — the *only* thing that differs between the two estimators |

**Paper prose defining the counts (§2.2, quoted verbatim from PMC3025716, re-verified):**

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
                          =  4·N_AA,aa  +  N_Aa-hom  +  N_hom-Aa
```

(the second line because `N_Aa^(i) + N_Aa^(j) − 2·N_Aa,Aa` is exactly the count of sites where
precisely one of the pair is heterozygous). This is the bridge between the "sum of squared
dosage differences" form and the "genotype counts" form. Every implementation below uses one
side or the other. **Verified numerically on 20 000 random pairs: held every time.**

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

The `+1/2` and `−(N_i+N_j)/(4Ĥ)` terms of (7) cancel exactly when `Ĥ_ij` is the average, so
**the within-family form has no additive 1/2 and no correction term.** Hail's documentation
states this closed form independently (§3), which is our second witness for it.

### 1.6 Equation (11) — KING-robust **BETWEEN**-family

Paper prose immediately preceding (quoted verbatim, re-verified):

> "In order to guard against potential estimation inflation due to departure from
> individual-level HWE, we consider **the smaller of the observed heterozygosity rates, min ( N
> Aa ( i ) / M ij , N Aa ( j ) / M ij )**, as an alternative to E (2 P (1 − P )). **Without loss
> of generality, suppose the i -th individual has lower heterozygosity than the j -th
> individual.** Then, the robust estimator is (11)"

So `Ĥ_ij = N_Aa^(i)` where **i is the individual with the SMALLER het count** (`M_ij` cancels
out of the ratio, so the *rate* min and the *count* min pick the same individual):

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
every equal-heterozygosity test and then be silently wrong on real pairs — see golden vectors
B and D′ in §6, where the naive form returns 0.5000 and −0.5000 against correct values of
0.1250 and −0.8750.

### 1.7 Which is used where — the paper's own rule (verbatim, re-verified)

> "The estimator above is no larger than the estimator in Equation ( 9 ), and **both estimators
> are bounded above by 0.5**. **We use estimator (Equation 9 ) for within-family relationship
> checking and estimator (Equation 11 ) for between-family relationship checking**, naming this
> combined approach KING-robust."

Note that 0.5 is an *analytic* upper bound under the model, not a clamp — see §8.4.

On negative values (verbatim, re-verified) — this is why clamping is wrong:

> "…is a consistent estimator of a parameter with a negative value (10) Thus, the robust
> estimator also can be used to determine the extent of population heterogeneity between the
> pair of individuals; **an extreme negative value indicates the pair of individuals is drawn
> from two distinct populations**."

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
**entire sample**, not the pair. SNPRelate implements exactly these two (§4.3) — independent
confirmation.

---

## 2. PLINK 2.0

### 2.1 What the documentation says (verbatim, fetched)

From <https://www.cog-genomics.org/plink/2.0/distance#make_king>:

- "**--make-king** writes KING-robust coefficients in matrix form to `plink2.king[.zst]` or `plink2.king.bin`"; "**--make-king-table**" writes table form to `plink2.kin0[.zst]`.
- "KING kinship coefficients are scaled such that **duplicate samples have kinship 0.5, not 1**. First-degree relations (parent-child, full siblings) correspond to ~0.25, second-degree relations correspond to ~0.125, etc."
- "**Only autosomes are included in this computation.**"
- "**Pedigree information is currently ignored; the between-family estimator is used for all pairs.**"
- "For multiallelic variants, REF allele counts are used."
- KING-robust "doesn't require MAFs at all" and "can also be mostly trusted on mixed-population datasets", but "KING-robust **underestimates kinship when the parents are from very different populations**."
- `--king-table-filter` reports only coefficients ≥ threshold; `--king-table-subset` restricts to pairs listed in an existing `.kin0` (the documented two-step workflow).
- Modifiers: `counts`, `rel-check`, `concordance-check`, `cols=`.
- **No claim of numeric or byte identity with KING's `--kinship` is made anywhere on the page.**
  The only reference is a pointer to "the original KING software package, which has some useful
  two-step workflows directly built in".

### 2.2 What the source actually computes

`plink2_matrix_calc.cc` (GPL-3+), `ComputeKinship()` at line 1567 and the same expression
inlined in the `.kin0` writer at line 2298. Restated in our notation (formula only — no code
reproduced):

Per-pair accumulators, over SNPs where **neither** call is missing: `ibs0` (opposite homs),
`hethet`, `het1hom2` (sample 1 het, sample 2 hom), `het2hom1`, `homhom`.

```
smaller_het  =  hethet + min(het1hom2, het2hom1)
KINSHIP      =  0.5  −  ( 4·ibs0 + het1hom2 + het2hom1 ) / ( 4·smaller_het )
```

Two identities make this exactly paper equation (11):

- `N_Aa^(1) = hethet + het1hom2` and `N_Aa^(2) = hethet + het2hom1`, therefore
  `smaller_het = min(N_Aa^(1), N_Aa^(2))` — **min denominator, confirmed**;
- `4·ibs0 + het1hom2 + het2hom1 = Σ_m (X^(1) − X^(2))²` by the §1.2 identity — so this is
  literally `0.5 − Σ(ΔX)²/(4·min)`, i.e. equation (11).

**Missing data:** pairwise. `het1hom2` requires sample 2 to be a called hom, so the het counts
can only accumulate on doubly-called sites — structurally identical to the paper's definition.

**Within-family:** not implemented. Equation (11) is applied to every pair including same-FID
pairs. There is no equation (9) code path.

**One thing that will confuse a reader of the source.** The *live* `ComputeKinship` takes four
extra `singleton_het/hom` arguments and folds them into `ibs0_ct`, `het1hom2_ct`, `het2hom1_ct`
before applying the formula; the plain textbook version is present in the file only as a
commented-out block immediately above (line 1555). This is plink2's sparse/singleton-variant
optimisation — very rare variants are accumulated per-sample and re-added per-pair instead of
being scanned pairwise. **It changes nothing about the estimator**; it is a different route to
the same integer counts. Do not mistake it for a variant of the formula. (Note also the source
comment "'2' here refers to the larger index, so this is swapped" — an index-convention wrinkle
local to plink2, not a semantic difference.)

### 2.3 `.kin0` output shape

Header emitted by `AppendKingTableHeader()` (line 1611), tab-separated, `#`-prefixed:

```
#FID1  IID1  [SID1]  FID2  IID2  [SID2]  NSNP  HETHET  IBS0  HET1_HOM2  HET2_HOM1  IBS  KINSHIP
```

- FID columns present only when FIDs are in play; SID columns only when SIDs are.
- **The docs page says `ID1`/`ID2`; the current source emits `IID1`/`IID2`.** The source carries
  an explicit comment: it "Was 'ID1' before alpha 3, but that's inconsistent with other plink2
  commands, and in the meantime **the header line still doesn't perfectly match KING due to e.g.
  capitalization**." That is an outright statement that plink2 does **not** target byte parity
  with KING.
- `NSNP = het1hom2 + het2hom1 + homhom + hethet` (writer, `nonmiss_ct`) — i.e. `M_ij`. Note
  `homhom` includes the opposite-hom sites, so `ibs0 ⊆ homhom`.
- `IBS` = Hamming distance = `2·ibs0 + het1hom2 + het2hom1` in `counts` mode; in the default mode
  it is `0.5 · hamming / NSNP` — **half** the Hamming distance as a proportion, not the raw
  proportion. Easy to trip over.
- `HETHET`/`IBS0`/`HET1_HOM2`/`HET2_HOM1` are **proportions of NSNP by default**; the `counts`
  modifier switches them to raw counts.
- Numbers are written with plink2's `dtoa_g` (shortest-round-trip `%g`-style), **not** a fixed
  `%.4f`. Formatting is therefore not comparable to KING.

### 2.4 plink2's *reader* is third-party evidence about KING's own format

`KingCutoffBatchTable()` (line ~685) and `CalcKingTableSubset()` (line ~3391) both parse a
`.kin0`, with the comment **"Make this work with both KING- and plink2-generated .kin0 files."**
What they accept is a statement about KING's column naming:

| position | accepted spellings |
|---|---|
| first token | `#FID1` **or** `FID` |
| ID column 1 | `ID1` **or** `IID1` |
| optional | `SID1` |
| FID column 2 | `FID2` |
| ID column 2 | `ID2` **or** `IID2` |
| kinship column (searched for by name) | `KINSHIP` **or** `Kinship` |

The `ID1`/`ID2`/`Kinship` alternatives exist solely to accommodate KING, and they match the KING
manual exactly (§8.6). ⚠️ **Open discrepancy worth a live check:** the manual shows KING's
`.kin0` header starting with `FID1` (no `#`), which matches *neither* accepted spelling
(`#FID1` nor bare `FID`) under plink2's exact-match comparison. Either recent KING emits a
`#`-prefixed header, or plink2's compatibility path is aimed at a different KING file. Resolve
by running the binary; it does not affect our estimator, only our header bytes.

### 2.5 Matrix output and ordering

- `.king` is lower-triangular text excluding the diagonal by default; `square`, `square0`,
  `triangle`, `bin`, `bin4` modifiers as for `--make-rel`. `.king.id` carries the sample IDs.
- Source comment at line 2011: results "are always reported in **lower-triangular order, rather
  than KING's upper-triangular order**, since the former plays more nicely with…". A direct,
  citable third-party statement about KING's own pair emission order.

### 2.6 Verdict

PLINK 2.0 computes **exactly paper equation (11)**, pairwise-missing, min-denominator,
between-family-only. It **claims no identity with KING's `--kinship`**, and documents three
deliberate departures: pedigree is ignored, output ordering differs, header naming differs.

---

## 3. Hail — `hl.king`

Docs (<https://hail.is/docs/0.2/methods/relatedness.html>, fetched) state **both** estimators —
making Hail the only third-party *documentation* that writes out equation (9) in closed form:

**Within-family (documented, not implemented):**
```
φ̂_ij^within  = ( N^(Aa,Aa)_ij − 2·N^(AA,aa)_ij ) / ( N^(Aa)_i + N^(Aa)_j )
```

**Between-family (documented and implemented):**
```
φ̂_ij^between = 1/2 + ( 2·N^(Aa,Aa)_ij − 4·N^(AA,aa)_ij − N^(Aa)_i − N^(Aa)_j )
                     / ( 4 · min( N^(Aa)_i , N^(Aa)_j ) )
```

- **The within form is character-for-character our §1.5 equation (9).** Independent confirmation
  that the within-family denominator is the **sum** and that there is no additive 1/2.
- **Denominator, between: min.** The docs describe the substitution explicitly — the between
  estimator "replaces the average heterozygote count with the minimum".
- **Missing data: pairwise.** Verbatim: *"The three counts above, N^{Aa}, N^{Aa,Aa}, and
  N^{AA,aa}, exclude variants where one or both individuals have missing genotypes."*
- **Within-family: not implemented.** Verbatim: *"This function, king(), only implements the
  'between-family' estimator."*
- Genotype score from `n_alt_alleles()` ∈ {0,1,2}.

Expanding `N^(Aa)_i + N^(Aa)_j = 2·hethet + het1hom2 + het2hom1` turns Hail's between numerator
into `−(4·ibs0 + het1hom2 + het2hom1)`, i.e. plink2's numerator with the sign folded into the
leading `1/2`. Algebraically identical (verified numerically in §6).

No claim of parity with the KING binary.

---

## 4. SNPRelate — `snpgdsIBDKING`

The most informative source, because it is the only third-party implementation that implements
**both** paper estimators **and** the within/between selection rule.

### 4.1 Counting (`src/genKING.cpp`, GPL-3)

Per pair, over a `mask` defined as "called in individual 1 AND called in individual 2":

| accumulator | meaning (from the struct's own doc comments, lines 276-280) |
|---|---|
| `nLoci` | popcount of the pairwise mask = `M_ij` |
| `IBS0` | "the number of loci sharing no allele" = `N_AA,aa` |
| `SumSq` | `\sum_m (g_m^{(i)} - g_m^{(j)})^2` = `Σ(ΔX)²` |
| `N1_Aa` | "the number of hetet loci for the first individual" |
| `N2_Aa` | second individual |

`N1_Aa` / `N2_Aa` are each ANDed with the pairwise `mask` before popcount → **pairwise
missingness, explicitly**. `SumSq` is accumulated as `popcount(het-term) + 4·popcount(ibs0)` —
the §1.2 identity again.

### 4.2 The estimator selection (matrix path line ~636, vector path line ~661)

```
kinship = (f1 == f2 && f1 != NA)                         // same, non-missing family id
            ?  0.5 − SumSq / ( 2 · (N1_Aa + N2_Aa) )     ← paper eq (9),  SUM denominator
            :  0.5 − SumSq / ( 4 · min(N1_Aa, N2_Aa) )   ← paper eq (11), MIN denominator
if (!R_FINITE(kinship)) kinship = NaN
```

These are *literally* the equation-(5) master form with `Ĥ_ij` set to the average and to the
min respectively — the cleanest independent confirmation available that the two KING-robust
estimators differ **only** in `Ĥ_ij`.

The family-id switch comes from the R wrapper's `family.id=` argument (`R/IBD.R`):
`family.id=NULL` (the default) fills all ids with `NA`, and the C++ requires `f1 == f2 && f1 !=
NA` for the within-family branch. So **SNPRelate's default behaviour equals PLINK 2.0's**
(between-family everywhere), but supplying `family.id` reproduces KING's own within/between
split. This makes it our only readable oracle for equation (9).

### 4.3 Other SNPRelate details worth mirroring in tests

- Diagonal is written as `kinship = 0.5`, `IBS0 = 0` — self-kinship is *asserted*, not computed.
- The `IBS0` output is a **proportion**: `IBS0 / nLoci`, or `NaN` when `nLoci == 0`.
- Any non-finite kinship is coerced to `NaN` (`R_FINITE` guard) — so a zero denominator surfaces
  as `NaN`, not `-Inf`. **This differs from PLINK 2.0** (see §8.3).
- `type="KING-homo"` computes `theta = 0.5 − SumSq/(8·SumAFreq)` and `k0 = IBS0/(2·SumAFreq2)`,
  `k1 = 2 − 2·k0 − 4·theta` — exactly paper equation (5)-with-homo-plug-in and equation (2).
  Independent confirmation of §1.8.

---

## 5. Illumina akt — `kin.cpp` (recorded for completeness; **do not use**, see §0)

`Kinship::estimateKinship`, `method == 1` (`-M 1`, "king-robust"), lines 195-207. The kinship
line is a direct transcription of the right-hand side of paper equation (11):

```
minhet = min(Nhet_1, Nhet_2)
ks     = (Nhet_12 − 2·ibd0) / (2·minhet)  +  0.5  −  0.25·(Nhet_1 + Nhet_2) / minhet
```

with `ibd0 = N_AA,aa`, `Nhet_12 = N_Aa,Aa`. The source's own inline comments annotate the
variables `//NAa^i`, `//NAa^j`, `//NAa,Aa`, confirming the mapping to the paper.

- **Missing data: pairwise** — the het counts are masked by `(missing_1 | missing_2).flip()`,
  and the source notes the hethet count "no mask needed here" because het∧het already implies
  both called. Same reasoning as plink2's.
- `ibd0` is a `float&`, so the first division is float — no integer-division trap despite the
  integer-looking counts. But note the accumulators are `float`, i.e. **single** precision:
  akt's own output is therefore not a full-precision oracle.
- `estimateIBD()` is called *after* `ks` is computed, so `ibd0` is still a raw count in the
  kinship expression.

**Contributes no independent information** (it is the paper's equation verbatim) and is under a
non-open-source licence. Excluded from the reference set.

---

## 6. Numerical verification that all stated forms are one formula

Implemented all nine published/observed expressions independently in Python and evaluated them
on 20 000 random genotype pairs (20-60 SNPs each, alleles drawn uniformly from {0, 1, 2,
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
2. paper eq (9) RHS / Hail within: `(N_AaAa − 2N_AAaa)/(N_i+N_j)`
3. SNPRelate: `0.5 − SumSq/(2·(N_i+N_j))`

Both bridging identities were asserted on every trial and held every time:
- `Σ(ΔX)² == 4·ibs0 + het1hom2 + het2hom1`
- `min(N_i, N_j) == hethet + min(het1hom2, het2hom1)`
- `M_ij == hethet + het1hom2 + het2hom1 + homhom`

**Mismatches: 0 / 20 000.**

### 6.1 Golden vectors for our unit tests (recomputed from scratch, exact rationals)

Genotypes as dosages `0/1/2`, `.` = missing, read left to right as one SNP per character.
Counts are over pairwise-called sites only. Values are given as exact fractions and as KING's
`%.4f` rendering.

| # | g1 | g2 | M | hethet | het1hom2 | het2hom1 | ibs0 | homhom | N_i | N_j | Σ(ΔX)² | **between eq(11)** | **within eq(9)** | naive (no correction) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A** duplicate/MZ | `1111002020` | `1111002020` | 10 | 4 | 0 | 0 | 0 | 6 | 4 | 4 | 0 | `1/2` → **0.5000** | `1/2` → **0.5000** | 0.5000 |
| **B** unequal het, no IBS0 | `1111102020` | `1102002020` | 10 | 2 | 3 | 0 | 0 | 5 | 5 | 2 | 3 | `1/8` → **0.1250** | `2/7` → **0.2857** | **0.5000** ✗ |
| **C** IBS0 + missing, equal het | `021120.120` | `2011201.00` | 8 | 2 | 0 | 0 | 3 | 6 | 2 | 2 | 12 | `−1` → **−1.0000** | `−1` → **−1.0000** | −1.0000 |
| **D′** IBS0 + missing + unequal het | `02.2000110` | `11.001211.` | 8 | 2 | 0 | 3 | 2 | 3 | 2 | 5 | 11 | `−7/8` → **−0.8750** | `−2/7` → **−0.2857** | **−0.5000** ✗ |
| **E** unequal het + IBS0, mild | `1210112010` | `1200110210` | 10 | 4 | 1 | 0 | 2 | 5 | 5 | 4 | 9 | `−1/16` → **−0.0625** | `0` → **0.0000** | 0.0000 ✗ |
| **F** zero denominator | `0022002020` | `1212112121` | 10 | 0 | 0 | 6 | 1 | 4 | 0 | 6 | 10 | **min = 0 → undefined** | `−1/3` → **−0.3333** | undefined |
| **G** no overlap | `01....2...` | `....12..01` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **undefined** | **undefined** | undefined |

Proportion columns for the same vectors (KING emits proportions — see §8.6):

| # | HetHet = hethet/M (`%.3f`) | IBS0 = ibs0/M (`%.4f`) |
|---|---|---|
| A | `0.400` | `0.0000` |
| B | `0.200` | `0.0000` |
| C | `0.250` | `0.3750` |
| D′ | `0.250` | `0.2500` |
| E | `0.400` | `0.2000` |
| F | `0.000` | `0.1000` |
| G | undefined (M = 0) | undefined (M = 0) |

**How to use these:**

- **B and D′ are the discriminating tests.** They are the only rows where all three candidate
  formulas (correct-between, correct-within, naive-without-correction) give three *different*
  plausible-looking answers. B catches "dropped the correction term" spectacularly: the naive
  form returns **0.5000**, i.e. it would report an unrelated pair as a duplicate.
- **A is worthless as a test on its own** — every wrong variant returns 0.5.
- **C is a trap**: it has IBS0 *and* missingness but equal het counts, so between and within
  coincide. A test suite containing only A and C cannot distinguish min from sum.
- **E** is the gentlest realistic case where between ≠ within (−0.0625 vs 0.0000) — good for
  catching a sign or off-by-one in the correction term without extreme values.
- **F and G** are the degenerate-input fixtures for §8.3.

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

  ≡ 0.5 − ( 4·ibs0 + het1hom2 + het2hom1 ) / ( 2 · (N_i + N_j) )
```

### 7.4 Reported proportions

```
HetHet = hethet / M_ij          (KING prints %.3f)
IBS0   = ibs0   / M_ij          (KING prints %.4f)
N_SNP  = M_ij                   (integer)
```

### 7.5 Implementation notes that follow from the consensus

- **Accumulate integers, divide once.** Both estimators are an exact-integer numerator over an
  exact-integer denominator with a single floating-point division. This eliminates
  summation-order sensitivity entirely and makes results bit-reproducible regardless of
  SIMD/threading. Do **not** accumulate `Σ(ΔX)²` as a float (akt's `float` accumulators are a
  cautionary example).
- **The `min` picks the individual, not the term.** `min_het` is `min(N_i, N_j)`; computing it as
  `hethet + min(het1hom2, het2hom1)` is equivalent (proved in §6) and avoids materialising
  `N_i`/`N_j`.
- **Never use a global per-sample heterozygote count.** `N_i` must be recomputed per pair under
  the pairwise-called mask. This is the subtlest requirement in the whole estimator and the most
  likely source of a near-miss that only shows up on samples with differing call rates.
- **A two-bit genotype encoding is the paper's own suggestion**, not just an optimisation: "When
  each genotype is stored in two bits, N_Aa^(i), N_Aa^(j), N_Aa,Aa and N_AA,aa can be computed
  using only bit operations (i.e. AND, OR, XOR and NOT), eliminating multiplication and division
  during the process of scanning the genotypes." All three readable implementations do exactly
  this. Our popcount-based accumulator is on the published path, not derived from anyone's code.
- **KING prints `%.4lf` for Kinship.** With four decimals, last-ULP differences between
  algebraically equivalent orderings are invisible except at exact `.00005` rounding ties, so we
  need not guess KING's exact operation order — but keep the single-division form so ties are at
  least deterministic on our side.

### 7.6 Corroboration from the KING 2.3.2 binary and manual

- KING manual, `.kin` (within-family, `--kinship`) header, verbatim:
  `FID  ID1  ID2  N_SNP  Z0  Phi  HetHet  IBS0  Kinship  Error`
- KING manual, `.kin0` (between-family) header, verbatim:
  `FID1  ID1  FID2  ID2  N_SNP  HetHet  IBS0  Kinship`
- Manual's own column definitions: "N_SNP: The number of SNPs that do not have missing genotypes
  in either of the individual" (= `M_ij`, pairwise — third confirmation); "HetHet: **Proportion**
  of SNPs with double heterozygotes"; "IBS0: **Porportion** [sic] of SNPs with zero IBS".
- Manual on negatives: "The reason that a negative kinship coefficient is not set to zero is a
  very negative value may indicate the population structure between the two individuals." —
  **explicit confirmation that KING does not clamp.**
- Manual degree thresholds: ">0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884]" for
  duplicate/MZ, 1st, 2nd, 3rd degree. The binary's embedded R plotting script carries the same
  cutoffs at full precision — `0.04419`, `0.08839`, `0.17678`, `0.35355` — which are `2^-4.5`,
  `2^-3.5`, `2^-2.5`, `2^-1.5`, i.e. `2^-(d + 3/2)`. **Use the exact powers of two, not the
  rounded manual values.**
- Binary string `"1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf"` —
  confirms IBS0 is the PO-vs-FS discriminator and that KING formats it with 4 decimals.

---

## 8. Where the sources DISAGREE — the parity-test targets

None of these are disagreements about the estimator's algebra. They are disagreements about
*application*, and they are precisely where our KING parity work must concentrate.

### 8.1 ⚠️ HIGH — Which estimator is applied to which pair

| Source | same-FID pairs | cross-FID pairs |
|---|---|---|
| **KING (paper §2.3 + the `.kin`/`.kin0` split)** | **eq (9), sum denominator** | eq (11), min denominator |
| SNPRelate with `family.id=` | eq (9), sum | eq (11), min |
| SNPRelate default (`family.id=NULL`) | eq (11), min | eq (11), min |
| **PLINK 2.0** | **eq (11), min** — "Pedigree information is currently ignored" | eq (11), min |
| **Hail** | **eq (11), min** — "only implements the 'between-family' estimator" | eq (11), min |
| akt | eq (11), min | eq (11), min |

**Consequence for us:** if we validate `--kinship` only against PLINK 2.0, every within-family
pair will appear to disagree and we will "fix" the wrong thing. Our `.kin` writer must use
eq (9); our `.kin0` writer must use eq (11). Validate `.kin0` against plink2, and `.kin` against
SNPRelate-with-`family.id` (or against the KING binary directly).

Magnitude of the divergence: golden vector B gives 0.1250 vs 0.2857 (factor 2.3); D′ gives
−0.8750 vs −0.2857 (factor 3.1). Divergence is **exactly zero** when `N_i = N_j`, so
same-ancestry, same-call-rate pairs will look deceptively fine.

### 8.2 ⚠️ HIGH — Which variants enter the counts

No two sources agree, and **the paper does not specify this at all** — it is program policy, not
part of the estimator:

- **PLINK 2.0:** "Only autosomes are included in this computation." No MAF filter, no call-rate
  filter, no LD pruning by default. "For multiallelic variants, REF allele counts are used."
- **SNPRelate / Hail:** whatever variant set the caller passes in; no implicit filtering.
- **KING 2.3.2:** the binary contains
  `"%d autosome SNPs with MAF>%.3lf and call rate>%d%% are used."`, plus separate call-rate
  filter messages for autosomes, chrX and chrY. **KING pre-filters variants by MAF and call rate
  by default.** The exact default thresholds must be read off a real run — they are not in the
  paper, and they are the largest single threat to `N_SNP` parity (and therefore to every
  proportion column, though *not* to the kinship formula itself).

**Action:** treat "which variants KING kept" as a separate, independently-verified stage. Get
`N_SNP` matching before comparing `Kinship`. A `Kinship` mismatch with a matching `N_SNP` is a
formula bug; a `Kinship` mismatch with a mismatched `N_SNP` is a filtering bug, and the two need
completely different fixes.

### 8.3 ⚠️ MEDIUM — Zero-denominator behaviour

Reached when `min(N_i, N_j) = 0` (golden vector F) or `M_ij = 0` (vector G).

| Source | Result |
|---|---|
| **PLINK 2.0** | `−Inf`. The source carries a dated note (18 Nov 2017): "kinship_coeff can be -inf when smaller_het_ct is zero. Don't filter those lines out when --king-table-filter wasn't specified." So the row is **emitted**, with `-inf`. |
| **SNPRelate** | `NaN` — the `R_FINITE` guard rewrites any non-finite value; `IBS0` becomes `NaN` when `nLoci == 0`. |
| **Hail** | unspecified in the docs. |
| **KING 2.3.2** | **UNKNOWN — must be observed.** |

**Action:** build fixtures F and G, run the real binary, and record its literal output bytes
(`-inf`? `nan`? `-1.#IND`? `0.0000`? a skipped row?). This is a genuine open question, not
resolvable from any documentation, and it is cheap to answer once.

### 8.4 ⚠️ MEDIUM — Clamping

**No implementation clamps, and KING's manual says it deliberately does not.** The paper says
both estimators "are bounded above by 0.5", but that is an asymptotic property of the estimator,
not an enforced range; the *lower* bound is unbounded (golden vector C yields −1.0, and the
paper's equation (10) discusses large negative values as a *signal* of population heterogeneity).

**Action:** do not clamp, in either direction. Verify against KING that large negative kinships
pass through unmodified rather than being floored at 0 or −1.

### 8.5 ⚠️ MEDIUM — Pair ordering and output layout

- PLINK 2.0 source, verbatim: results "are always reported in **lower-triangular order, rather
  than KING's upper-triangular order**". So KING emits pairs in upper-triangular order — a
  citable third-party statement, but confirm it against the binary.
- PLINK 2.0's `.king` matrix omits the diagonal by default; SNPRelate writes `0.5` on the
  diagonal. KING's matrix conventions must be observed separately.
- KING's `.kin` groups by family (the manual's example is sorted by FID, then ID1, then ID2);
  `.kin0` iterates family-pairs. Order must come from the binary, not from us.

### 8.6 ⚠️ HIGH (corrected) — Column semantics: **proportions, not counts**

This is the item the earlier draft got backwards, and it changes the cross-validation recipe.

| | KING 2.3.2 | PLINK 2.0 |
|---|---|---|
| within-family file | `.kin` | *(not produced)* |
| between-family file | `.kin0` | `.kin0` |
| `.kin0` ID columns | `FID1  ID1  FID2  ID2` | `#FID1  IID1  FID2  IID2` |
| SNP count column | `N_SNP` | `NSNP` |
| het-het column | `HetHet` | `HETHET` |
| kinship column | `Kinship` | `KINSHIP` |
| `HetHet` semantics | **proportion** `hethet/N_SNP` | **proportion by default**, counts with the `counts` modifier |
| `IBS0` semantics | **proportion** `ibs0/N_SNP` | **proportion by default**, counts with `counts` |
| extra plink2 columns | — | `HET1_HOM2`, `HET2_HOM1`, `IBS` (half-Hamming) |
| number format | fixed-decimal: `HetHet` `%.3f`, `IBS0` `%.4f`, `Kinship` `%.4f`, `Z0` `%.3f`, `Phi` `%.4f` | `dtoa_g` shortest round-trip |

Evidence for the KING side: the manual's own definitions ("HetHet: Proportion of SNPs with
double heterozygotes"; "IBS0: Porportion of SNPs with zero IBS") plus its worked example rows,
e.g. `28  3  117  1  2360618  0.143  0.0267  0.1356` — three decimals for HetHet, four for
IBS0 and Kinship, on a 2.36 M-SNP run where counts would be seven-digit integers.

**Consequence for parity testing (revised):**
- plink2's **default** output (proportions) matches KING's *semantics* — use it to check
  `HetHet` and `IBS0` values.
- plink2's **`counts`** modifier is what you want to verify our *integer accumulators*, since it
  exposes `hethet`, `ibs0`, `het1hom2`, `het2hom1` as exact integers that we can reconcile
  against our own before any division happens.
- Run both. Neither alone covers both concerns.
- plink2's own source concedes "the header line still doesn't perfectly match KING due to e.g.
  capitalization". **PLINK 2.0 is a correctness oracle for the kinship *value* and the counts;
  it is useless as a byte-format oracle.** Format parity must come from the binary.

### 8.7 Settled — no disagreement (do not re-litigate)

- **Denominator, between-family: `min(N_Aa^i, N_Aa^j)`.** Unanimous: paper eq (11), plink2, Hail,
  SNPRelate, akt. The `sum` denominator belongs to eq (9) *only*. **The task brief's question
  "min or sum?" has a clean answer: both, in different estimators — and nobody disagrees about
  which goes where.**
- **Missing data: pairwise, always.** Unanimous, and stated outright in the paper ("excluding
  those SNPs with missing genotypes in either individual of the pair"), in Hail's docs, in
  plink2's `NSNP` definition, in SNPRelate's mask, in akt's mask, and in KING's own manual
  ("N_SNP: The number of SNPs that do not have missing genotypes in either of the individual").
- **`N_Aa^(i)` is pair-specific, never a global per-sample het count.** Unanimous.
- **The `+1/2 − (N_i+N_j)/(4·min)` correction is part of eq (11).** Unanimous.
- **Within-family form has no additive 1/2** — `(hethet − 2·ibs0)/(N_i + N_j)`. Paper and Hail
  state the identical closed form; SNPRelate computes the equivalent `0.5 − SumSq/(2·(N_i+N_j))`.
- **Scaling: duplicates → 0.5, not 1.** Unanimous.
- **No clamping.** Unanimous, and explicit in KING's manual.

---

## 9. Recommended parity-test plan (falls out of §8)

1. **Stage the comparison.** Match `N_SNP` first (variant filtering, §8.2), then the count-level
   quantities, only then `Kinship`. A kinship-only comparison cannot distinguish a filtering bug
   from a formula bug.
2. **Test both estimators separately.** Build a fixture with ≥2 individuals sharing an FID and
   ≥2 in different FIDs, so `.kin` (eq 9) and `.kin0` (eq 11) are exercised in one run.
3. **Use golden vectors B and D′ (§6.1) as the regression guard** against min/sum confusion and
   against dropping the correction term. MZ-twin tests are worthless here — every wrong variant
   returns 0.5. Assert on the exact rationals (`1/8`, `2/7`, `−7/8`, `−2/7`), not on the `%.4f`
   strings, in the unit tests; assert on the strings only in the format tests.
4. **Force the degenerate cases:** vector F (`min_het = 0`), vector G (`M_ij = 0`), and a pair
   with an extreme negative kinship. Record KING's literal output bytes for each (§8.3, §8.4).
5. **Cross-validate the between-family value against plink2 twice** on the *same* variant set
   (feed plink2 the exact post-filter variant list so §8.2 cannot confound it):
   `--make-king-table counts` to reconcile our integer accumulators, and the default
   (proportion) form to reconcile `HetHet`/`IBS0` semantics against KING's. Expect agreement to
   full double precision on `KINSHIP`. Expect the header, column names, ordering and number
   formatting to differ (§8.5, §8.6) — do not chase those.
6. **Cross-validate the within-family value against `snpgdsIBDKING(..., family.id=)`**, the only
   readable implementation of equation (9).
7. **Resolve the `#FID1` vs `FID1` header question** (§2.4) from a real KING run before writing
   the `.kin0` writer's header bytes.

---

## Appendix A — artefacts on disk (all re-verified for this document)

Read for math only; none of it is to be copied into the implementation.

```
…/scratchpad/research/src/plink2_matrix_calc.cc      PLINK 2.0 / plink-ng, GPL-3-or-later
…/scratchpad/research/src/snprelate_genKING.cpp      SNPRelate, GPL-3
…/scratchpad/research/src/snprelate_IBD.R            SNPRelate R wrapper, GPL-3
…/scratchpad/research/src/akt_kin.cpp                Illumina akt, PolyForm Strict 1.0.0 — DO NOT USE
…/scratchpad/research/src/pmc3025716.html            Manichaikul et al. 2010, PMC copy
…/scratchpad/research/src/eq/btq559*.jpg             paper equation images
…/scratchpad/research/txt/manual.txt                 kingrelatedness.com manual, text extraction
…/scratchpad/gv.py                                   golden-vector generator (exact rationals)
…/scratchpad/verify_king.py                          20 000-trial equivalence harness
```

Key line anchors in the third-party sources (for a reviewer who wants to re-check the claims —
line numbers, not content):

- `plink2_matrix_calc.cc`: 1555 (commented plain `ComputeKinship`), 1567 (live formula), 1611
  (`AppendKingTableHeader`), 2011 (triangular-order comment), 2298 (inlined formula), 2303
  (`-inf` note), 686 / 3391 (KING-compatible `.kin0` readers).
- `snprelate_genKING.cpp`: 276-280 (accumulator doc comments), 418-422 (scalar counting),
  636 / 661 (the within/between branch).
- `akt_kin.cpp`: 195-207 (`method == 1`).

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
- SNPRelate source: <https://github.com/zhengxwen/SNPRelate> (`src/genKING.cpp`, `R/IBD.R`)
- Illumina akt: <https://github.com/Illumina/akt> (`kin.cpp`, `LICENSE`)
- KING manual: <https://www.kingrelatedness.com/manual.shtml>
- KING 2.3.2 binary (embedded string constants only — facts about output format, not source):
  `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
