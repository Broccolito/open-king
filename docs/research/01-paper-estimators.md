# KING 2.3.2 — Estimator Recon from the Source Paper

**Target paper:** Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen W-M.
"Robust relationship inference in genome-wide association studies."
*Bioinformatics* 2010;26(22):2867–2873. PMC3025716.

**Purpose:** clean-room MIT reimplementation of KING's relatedness inference.

---

## 0. Provenance and confidence labelling

Every claim below is tagged so the implementer knows how much to trust it.

| Tag | Meaning |
|---|---|
| **[PAPER]** | Stated in the 2010 paper; extracted via PMC full text. |
| **[DERIVED]** | Not readable directly (PMC renders the display equations as images), but **algebraically forced** by moment identities the paper does state, and cross-checked to reproduce an independently-documented form exactly. Treat as high confidence — the derivation is reproduced so you can verify it yourself. |
| **[HAIL]** | Independently documented in Hail's `hl.king` docs, which explicitly summarise this paper's Methods. Used as a second witness. |
| **[MANUAL]** | kingrelatedness.com/manual.shtml — the program documentation (post-2010 behaviour, KING 2.x). |
| **[BINARY]** | Observed from the reference binary: embedded string constants + `--help` output. Facts about the output format, not source code. |
| **[OPEN]** | Not resolved here; needs a live run of the binary. |

**Legal note.** No KING C++ source was read or consulted. Mathematical identities are
facts, not copyrightable expression. Binary string constants are facts about the output
format. Nothing in this document is transcribed from an implementation.

Reference binary: `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(Mach-O 64-bit arm64, banner `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`).

---

## 1. Notation

Fix a pair of individuals `i`, `j` and a set of autosomal biallelic SNPs indexed `m = 1..M`.
Let `A` be the reference allele and `a` the alternate; `p_m` = frequency of `A` at SNP `m`.

**Genotype score** **[PAPER]** — `x_i(m)` is the count of the reference allele in
individual `i` at SNP `m`:

```
x = 0  for genotype AA
x = 1  for genotype Aa
x = 2  for genotype aa
```

(The paper phrases this as the number of the reference allele. The orientation is
irrelevant to every estimator below — see §8.1.)

**Pairwise-complete marker set.** `M_ij` = the number of SNPs at which **both** `i` and `j`
have a non-missing genotype. **All counts below are taken over exactly this set.**

**Genotype-pair counts over `M_ij`:**

| Symbol | Definition |
|---|---|
| `N_Aa^i` | # SNPs where `i` is heterozygous (`Aa`) |
| `N_Aa^j` | # SNPs where `j` is heterozygous (`Aa`) |
| `N_Aa,Aa^ij` (a.k.a. `N_HetHet`) | # SNPs where **both** are heterozygous |
| `N_AA,aa^ij` (a.k.a. **IBS0**) | # SNPs where the two are **opposite homozygotes** (`AA`/`aa` or `aa`/`AA`) |
| `N_IBS0`, `N_IBS1`, `N_IBS2` | # SNPs sharing 0, 1, 2 alleles IBS (§3) |

**CRITICAL PARITY DETAIL [PAPER, verified]:** `N_Aa^i` and `N_Aa^j` are the heterozygote
counts for `i` and `j` **excluding SNPs missing in either member of the pair**. They are
*not* each individual's own marginal heterozygote count. This is stated explicitly in the
paper's definition of Equation (9). Getting this wrong produces subtly wrong kinship for
every pair with differential missingness.

---

## 2. Equation inventory of the paper **[PAPER]**

The Methods section runs Equations (1)–(12):

| Eq | Subject |
|---|---|
| (1) | Expected zero-IBS (opposite-homozygote) count/proportion under HWE, as a function of `π₀` |
| (2) | Estimator of `π₀` (Pr[IBD=0]) from the observed IBS0 count |
| (3) | Sample allele-frequency estimator `p̂_m` |
| (4) | Moment relation: expected squared genotype difference as a function of `p_m` and `φ_ij` |
| (5) | KING-homo kinship estimator, general form |
| (6) | KING-homo plug-in estimator using `p̂_m` |
| (7) | The same relation re-expressed in genotype counts (`N_Aa,Aa`, `N_AA,aa`, …) |
| (8) | Individual heterozygosity `E(2P(1−P))` estimated by `N_Aa / M_ij` |
| (9) | **KING-robust, within-family** (average heterozygosity denominator) |
| (10) | Expected value of the robust estimator for unrelated pairs from *different* populations (negative) |
| (11) | **KING-robust, between/across-family** (minimum heterozygosity denominator) |
| (12) | Variance of allele frequencies across populations, `Var(P)` |

---

## 3. IBS definitions and counting **[PAPER]**

IBS = number of alleles shared identical-by-state between the two genotypes at a SNP.
For biallelic SNPs the full table is:

| geno `i` | geno `j` | IBS | `(x_i − x_j)²` | contributes to |
|---|---|---|---|---|
| AA | AA | 2 | 0 | `N_IBS2`, `N_HomHom` |
| aa | aa | 2 | 0 | `N_IBS2`, `N_HomHom` |
| Aa | Aa | 2 | 0 | `N_IBS2`, `N_Aa,Aa` (HetHet) |
| AA | Aa | 1 | 1 | `N_IBS1` |
| Aa | AA | 1 | 1 | `N_IBS1` |
| Aa | aa | 1 | 1 | `N_IBS1` |
| aa | Aa | 1 | 1 | `N_IBS1` |
| **AA** | **aa** | **0** | **4** | **`N_IBS0` = `N_AA,aa`** |
| **aa** | **AA** | **0** | **4** | **`N_IBS0` = `N_AA,aa`** |

Consequences that the implementation should rely on (all exact, not approximate):

```
N_IBS0 = N_AA,aa                                  (IBS0 <=> opposite homozygotes)
N_IBS1 = N_Aa^i + N_Aa^j - 2 * N_Aa,Aa            (IBS1 <=> exactly one of the pair is het)
N_IBS2 = M_ij - N_IBS0 - N_IBS1
```

**Master identity [DERIVED, exact]** — links the squared-difference form to the count form:

```
SUM_m (x_i(m) - x_j(m))^2  =  N_IBS1 + 4 * N_IBS0
                           =  N_Aa^i + N_Aa^j - 2*N_Aa,Aa + 4*N_AA,aa
```

Denote this quantity `D_ij`. Every estimator below is `1/2 − D_ij / (2·denominator)`.

**Note on IBS1 [PAPER]:** the paper's parent–offspring criterion rests on the fact that a
PO pair shares at least one allele IBD at every autosomal locus, so IBS is 1 or 2 at every
SNP — i.e. `N_IBS0 = 0` for a true PO pair apart from genotyping error. This is the basis
of the PO-vs-FS split in §6.3.

---

## 4. KING-homo (homogeneous-population estimator)

### 4.1 Allele frequency **[PAPER]**, Eq (3)

```
p̂_m = (2 * #AA_m + #Aa_m) / (2 * N_m)
```
where `N_m` is the number of individuals with a non-missing genotype at SNP `m`, and
`#AA_m`, `#Aa_m` are genotype counts in the whole sample. (Sample-wide, not pair-wise.)

### 4.2 The moment relation **[DERIVED]**, Eq (4)

Under HWE at SNP `m`, with kinship `φ_ij`:

```
Var(x_i)          = 2 * p_m * (1 - p_m)
Cov(x_i, x_j)     = 4 * φ_ij * p_m * (1 - p_m)
E[(x_i - x_j)^2]  = Var(x_i) + Var(x_j) - 2*Cov(x_i, x_j)
                  = 4 * p_m * (1 - p_m) * (1 - 2 * φ_ij)
```

Sanity: `φ=0.5` (MZ) ⇒ 0. `φ=0` (unrelated) ⇒ `4p(1−p)` = twice the per-SNP heterozygosity.

> One PMC summarisation pass returned `2p(1−p)(1−4φ)` for this. That is **wrong** — it
> fails at `φ=0` (gives `2p(1−p)`, half the correct independent-genotype value) and it does
> not reproduce the count-form estimator in §5, which is independently confirmed by [HAIL].
> Use `4p(1−p)(1−2φ)`.

### 4.3 KING-homo kinship estimator **[DERIVED from Eq (4)+(6)]**

Solving the moment relation summed over SNPs:

```
                    SUM_m (x_i(m) - x_j(m))^2                      D_ij
  φ̂_ij^homo  =  1/2  -  ------------------------------  =  1/2  -  --------------------
                    4 * SUM_m 2*p̂_m*(1 - p̂_m)                  4 * Ĥ_exp
```

where `Ĥ_exp = SUM_m 2*p̂_m*(1−p̂_m)` is the expected heterozygote count over `M_ij`.
Equivalently, expanded fully:

```
  φ̂_ij^homo  =  1/2  -  ---------------------------------------------------
                          8 * SUM_m p̂_m * (1 - p̂_m)
```
(numerator `D_ij` as above).

Equivalent count form (this is Eq (7)'s shape):

```
                        2*N_Aa,Aa - 4*N_AA,aa - N_Aa^i - N_Aa^j
  φ̂_ij^homo  =  1/2  +  ---------------------------------------
                          4 * SUM_m 2*p̂_m*(1 - p̂_m)
```

Verification: MZ twins ⇒ `D_ij = 0` ⇒ `φ̂ = 0.5`. Unrelated ⇒ `E[D_ij] = SUM 4p(1−p) = 2·Ĥ_exp`
⇒ `φ̂ = 1/2 − 1/2 = 0`. ✓

**Property:** KING-homo requires allele frequencies and is therefore *not* robust to
population structure — mis-specified `p̂_m` under admixture/stratification inflates `φ̂`.
That is the entire motivation for KING-robust.

### 4.4 `π₀` (IBD0) estimator **[PAPER Eq (2), denominator DERIVED]**

Expected opposite-homozygote probability at SNP `m` for a pair sharing **0** alleles IBD:

```
P(i=AA, j=aa) + P(i=aa, j=AA) = p_m^2 * (1-p_m)^2 + (1-p_m)^2 * p_m^2 = 2 * p_m^2 * (1-p_m)^2
```

Pairs sharing 1 or 2 alleles IBD can never be opposite homozygotes, so Eq (1) is:

```
E[N_AA,aa] = π₀ * SUM_m 2 * p̂_m^2 * (1 - p̂_m)^2
```

and Eq (2) inverts it:

```
                     N_AA,aa
  π̂₀  =  ------------------------------------
          SUM_m 2 * p̂_m^2 * (1 - p̂_m)^2
```

> **Exponent warning.** A PMC summarisation pass rendered the denominator as
> `2p(1−p)²`. That is a transcription artifact — it is dimensionally wrong (it is not a
> probability of two independent homozygous genotypes). The correct term is
> `2 * p² * (1−p)²`, derived above from HWE genotype frequencies. Verify with a
> simulation before shipping.

---

## 5. KING-robust — the main estimator

KING-robust replaces the allele-frequency-based denominator with the **observed
heterozygote counts of the pair itself**, which makes it allele-frequency-free and
therefore insensitive to population structure. Eq (8) is the justification:
`E(2P(1−P)) ≈ N_Aa / M_ij`, i.e. an individual's observed heterozygote rate estimates
their own expected heterozygosity, whatever population they come from.

### 5.1 Within-family estimator **[PAPER Eq (9)] [HAIL confirmed]**

Count form (this is the canonical form; **implement this one**):

```
                     N_Aa,Aa^ij  -  2 * N_AA,aa^ij
  φ̂_ij^within  =  -------------------------------
                        N_Aa^i  +  N_Aa^j
```

Equivalent squared-difference form:

```
                              D_ij
  φ̂_ij^within  =  1/2  -  ----------------------
                          2 * (N_Aa^i + N_Aa^j)
```

Equivalent "expanded" form (same value; useful for cross-checking against other write-ups):

```
                            2*N_Aa,Aa - 4*N_AA,aa - N_Aa^i - N_Aa^j
  φ̂_ij^within  =  1/2  +  ---------------------------------------
                                2 * (N_Aa^i + N_Aa^j)
```

**Proof of equivalence** (do this check in a unit test):
`D_ij = N_Aa^i + N_Aa^j − 2·N_Aa,Aa + 4·N_AA,aa`, so
`1/2 − D_ij / (2(N_i+N_j)) = [(N_i+N_j) − N_i − N_j + 2N_Aa,Aa − 4N_AA,aa] / (2(N_i+N_j))
= (2N_Aa,Aa − 4N_AA,aa) / (2(N_i+N_j)) = (N_Aa,Aa − 2N_AA,aa)/(N_i+N_j)`. ∎

Sanity: MZ twins — `N_Aa,Aa = N_Aa^i = N_Aa^j = H`, `N_AA,aa = 0` ⇒ `H / 2H = 0.5`. ✓

The denominator `N_Aa^i + N_Aa^j` is (twice) the **average** heterozygosity of the pair —
appropriate when `i` and `j` are known to be from the same family, hence the same ancestry.

### 5.2 Between-family / across-family estimator **[PAPER Eq (11)] [HAIL confirmed]**

```
                              2*N_Aa,Aa^ij - 4*N_AA,aa^ij - N_Aa^i - N_Aa^j
  φ̂_ij^between  =  1/2  +  --------------------------------------------
                                 4 * min(N_Aa^i, N_Aa^j)
```

Equivalent squared-difference form:

```
                               D_ij
  φ̂_ij^between  =  1/2  -  --------------------------
                           4 * min(N_Aa^i, N_Aa^j)
```

Sanity: MZ twins ⇒ `D_ij = 0` ⇒ `0.5`. ✓
Unrelated, same population, `N_Aa^i = N_Aa^j = H` ⇒ `E[D_ij] = 2H` ⇒ `1/2 − 2H/4H = 0`. ✓

**Unified implementation.** Both variants are the same expression with a different
denominator `D`:

```
  φ̂ = 1/2 - D_ij / (2 * Denom)
      Denom = N_Aa^i + N_Aa^j          -> within-family
      Denom = 2 * min(N_Aa^i, N_Aa^j)  -> between-family
```
They coincide exactly when `N_Aa^i == N_Aa^j`.

### 5.3 Why `min()` **[PAPER Eq (10), (12)]**

- The average-heterozygosity denominator is inflated when the two individuals have
  different ancestry (their heterozygosities differ), so Eq (9) applied across families can
  be biased. Using the **smaller** of the two heterozygote counts is the conservative
  choice and guards against inflation from departures from individual-level HWE.
- Eq (10): for two **unrelated** individuals drawn from **different** populations, the
  robust estimator converges to a **negative** value, roughly of the shape
  `E[1 − 2P_i(1−P_i) − 2P_j(1−P_j)] / (2·E[2P_i(1−P_i)])` in the paper's notation. A
  strongly negative `φ̂` is therefore a *signal of divergent ancestry*, not of an error.
  **Do not clamp negative kinship to zero** — KING reports it and downstream users read it.
- Eq (12): variance of allele frequencies per individual, estimated as
  `Var(P) = [ (N_Aa^i/M_ij) − (N_Aa^i/M_ij)^2 ] / 2`, quantifying the ancestry divergence
  that drives the Eq (10) bias.
- **Residual bias [PAPER]:** when a pair is *both* related *and* from different
  populations, KING-robust is no longer consistent; the paper reports the effect is small
  for relationships out to 3rd degree.

### 5.4 KING's own usage rule **[PAPER]**

- Eq (9) (within-family, average) → pairs **within** a declared family (`.kin` output).
- Eq (11) (between-family, min) → pairs **across** families (`.kin0` output).

This matches the binary's two output files exactly (§7).

---

## 6. Relationship classification

### 6.1 Kinship cutoffs **[MANUAL] [PAPER Table 1]**

The boundaries are **geometric midpoints** between adjacent degrees. A degree-`d` relative
has expected `φ = 2^-(d+1)`; the boundary between degree `d` and `d+1` is
`sqrt(2^-(d+1) · 2^-(d+2)) = 2^-(d + 3/2)`.

| Boundary | Exact value | Decimal (10 dp) | Rounded as published |
|---|---|---|---|
| Dup/MZ vs 1st-degree | `2^(-3/2)` | 0.3535533906 | **0.354** |
| 1st vs 2nd degree | `2^(-5/2)` | 0.1767766953 | **0.177** |
| 2nd vs 3rd degree | `2^(-7/2)` | 0.0883883476 | **0.0884** |
| 3rd vs 4th degree | `2^(-9/2)` | 0.0441941738 | **0.0442** |
| 4th degree vs unrelated | `2^(-11/2)` | 0.0220970869 | **0.0221** |

Classification table:

| Inferred relationship | `InfType` label **[BINARY]** | Kinship interval | Expected `φ` |
|---|---|---|---|
| Duplicate / MZ twin | `Dup/MZ` | `φ̂ > 0.354` | 0.5 |
| 1st-degree (PO or FS) | `PO` / `FS` | `0.177 < φ̂ ≤ 0.354` | 0.25 |
| 2nd-degree | `2nd` | `0.0884 < φ̂ ≤ 0.177` | 0.125 |
| 3rd-degree | `3rd` | `0.0442 < φ̂ ≤ 0.0884` | 0.0625 |
| 4th-degree | `4th` | `0.0221 < φ̂ ≤ 0.0442` | 0.03125 |
| Unrelated | `UN` | `φ̂ ≤ 0.0442` (or `≤ 0.0221` with `--degree 4`) | 0 |

**[BINARY]** confirms the label string `Dup/MZ` is embedded verbatim, and `InfType` is a
column in `.kin`, `.kin0` and the `--ibdseg` outputs. `--degree D` controls how deep the
reporting goes; `--degree 4` extends the reported range down to the `0.0221` boundary.

> **Boundary-inclusivity is [OPEN].** Whether KING uses `>=` or `>` at each boundary must
> be pinned by running the binary on constructed edge cases. Ties are measure-zero on real
> data but matter for byte-parity tests.

### 6.2 `π₀` cutoffs **[PAPER Table 1]**

The paper's Table 1 gives a second inference axis using `π̂₀`:

| Relationship | `φ` | `φ` criterion | `π₀` | `π₀` criterion |
|---|---|---|---|---|
| MZ twin / duplicate | 0.5 | `> 0.354` | 0 | `< 0.1` |
| Parent–offspring | 0.25 | `(0.177, 0.354)` | 0 | `< 0.1` |
| Full sibling | 0.25 | `(0.177, 0.354)` | 0.25 | `(0.1, 0.365)` |
| 2nd-degree | 0.125 | `(0.0884, 0.177)` | 0.75 | `(0.365, 0.9)` |
| 3rd-degree | 0.0625 | `(0.0442, 0.0884)` | 0.9375 | `(0.9, 0.9912)` |
| Unrelated | 0 | `< 0.0442` | 1 | `> 0.9912` |

Expected IBD triples for the standard relationships (`φ = π₁/4 + π₂/2`):

| Relationship | `π₀` | `π₁` | `π₂` | `φ` |
|---|---|---|---|---|
| MZ / duplicate | 0 | 0 | 1 | 0.5 |
| Parent–offspring | 0 | 1 | 0 | 0.25 |
| Full sibling | 0.25 | 0.5 | 0.25 | 0.25 |
| Half sib / avuncular / grandparent (2nd) | 0.5 | 0.5 | 0 | 0.125 |
| First cousin (3rd) | 0.75 | 0.25 | 0 | 0.0625 |
| 4th degree | 0.875 | 0.125 | 0 | 0.03125 |
| Unrelated | 1 | 0 | 0 | 0 |

### 6.3 Parent–offspring vs full sibling **[PAPER] + [BINARY]**

Both have `φ = 0.25`, so kinship alone cannot separate them. The discriminator is **IBS0**:
a true PO pair shares ≥1 allele IBD at every autosomal locus, so `N_AA,aa = 0` up to
genotyping error; a FS pair has `π₀ = 0.25` and therefore a clearly non-zero IBS0 rate.

The binary implements this as a **thresholded, data-derived** cutoff on the **IBS0
proportion**. Two embedded format strings confirm it:

```
1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf
Cutoff value for IBS0 between FS and PO is set at %.4f
```

Both print with **4 decimal places**. The fact that the cutoff is *printed at runtime*
rather than being a documented constant means it is **computed from the data** (it must
scale with the MAF spectrum of the SNP set, since the unrelated-pair IBS0 rate does).

> **[OPEN] — highest-priority follow-up.** The exact formula for this cutoff is the single
> biggest parity risk in the whole reimplementation. Determine it by running the reference
> binary on several synthetic PLINK filesets with deliberately different MAF spectra and
> reading back the printed `%.4f` value. Hypothesis to test first: it is a fixed fraction
> of the mean observed IBS0 rate among inferred-unrelated pairs (FS expectation is
> `0.25 ×` the unrelated rate; a natural cutoff sits well below that, e.g. `0.1×`).

Alternatively, in `--ibdseg` mode PO vs FS is trivially separable from the segment
statistics: PO has `IBD1Seg ≈ 1.0` and `IBD2Seg ≈ 0`, whereas FS has `IBD2Seg ≈ 0.25`.

---

## 7. IBD estimators, PropIBD, and the output columns

### 7.1 `π̂₁`, `π̂₂` from `π̂₀` and `φ̂` **[PAPER]**

The paper derives `π₁` and `π₂` from the two fundamental relations:

```
π₀ + π₁ + π₂ = 1
φ           = π₁/4 + π₂/2
```

Solving:

```
  π̂₂  =  4*φ̂  +  π̂₀  -  1
  π̂₁  =  1  -  π̂₀  -  π̂₂   =   2  -  2*π̂₀  -  4*φ̂
```

(Derivation: `4φ = π₁ + 2π₂ = (1 − π₀ − π₂) + 2π₂ = 1 − π₀ + π₂`.)
These require `π̂₀`, hence allele frequencies — they belong to the **KING-homo** branch,
not to KING-robust.

### 7.2 `IBD1Seg`, `IBD2Seg`, `PropIBD` **[MANUAL] [BINARY]**

These come from KING 2.x's `--ibdseg` mode (segment-based IBD inference), **not** from the
2010 paper — the 2010 paper predates them. They are the proportions of the genome covered
by inferred IBD1 and IBD2 **segments**.

```
  PropIBD  =  IBD2Seg  +  IBD1Seg / 2
```

Note this is the segment analogue of `φ = π₁/4 + π₂/2` scaled by 2: `PropIBD = 2φ`.
So `PropIBD` for MZ = 1.0, PO = 0.5, FS = 0.5, 2nd = 0.25, 3rd = 0.125.

### 7.3 Output columns observed in the binary **[BINARY]**

`.kin` (within-family), in embedded order:

```
FID  ID1  ID2  N_SNP  Z0  Phi  HetHet  IBS0  HetConc  HomIBS0  Kinship  IBD1Seg  IBD2Seg  PropIBD  InfType  Error
```

`.kin0` (across-family), in embedded order:

```
FID1  ID1  FID2  ID2  N_SNP  HetHet  IBS0  HetConc  HomIBS0  Kinship  IBD1Seg  IBD2Seg  PropIBD  InfType
```

A reduced `.kin0` variant also appears (`FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship`) and a
reduced `.kin` variant (`... N_SNP HetHet IBS0 Kinship Error`) — KING emits different column
sets depending on which analysis produced the file.

X-chromosome variants **[BINARY]**: `X.kin` uses `Sex1 Sex2 N_SNP PhiX IBD1Seg IBD2Seg PropIBD`;
`X.kin0` uses `FID1 ID1 FID2 ID2 N_SNP IBS0 KinshipX`.

`.con` (duplicate detection, `--duplicate`) **[BINARY]**:
```
FID1  ID1  FID2  ID2  N_IBS0  N_IBS1  N_IBS2  Concord  HomConc  HetConc
```

`.ibs` (`--ibs`) **[BINARY]**:
```
... N_SNP  N_IBS0  N_IBS1  N_IBS2  NHetHet  NHomHom  N_Het1  N_Het2  Dist  HetConc  Het2|1  Het1|2  HomConc  Kinship  MaxIBD2  Pr_IBD2
```

**Number formatting [BINARY]:** the dominant format specifier throughout all relatedness
writers is **`%.4lf`** (4 decimals). A `%.3lf` appears at the head of each row group —
consistent with the pedigree-derived `Z0` / `Phi` columns being printed at **3 decimals**
while all data-derived columns (`HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`,
`IBD1Seg`, `IBD2Seg`, `PropIBD`) use 4. `N_SNP` is an integer (`%d`).

> **[OPEN]** The exact column↔format binding and the field separator / padding must be
> confirmed by running the binary and diffing bytes. Do not guess.

`Error` column semantics **[MANUAL]**: flags disagreement between the inferred and the
pedigree-reported relationship — **`1` = error, `0.5` = warning**, `0` otherwise.
`Z0` and `Phi` are the **pedigree-expected** `π₀` and `φ`, not estimates.

---

## 8. Heterozygote concordance and error rates

**The 2010 paper does not define a heterozygote-concordance rate or a genotyping-error-rate
formula** — verified against the full text. These are KING 2.x implementation features.
What can be established:

### 8.1 Duplicate detection threshold **[BINARY]**

`--help` reports the default:

```
--minConc [0.80]
```

so duplicate/MZ pairs are declared when the **heterozygote concordance rate exceeds 0.80**
by default. Runtime messages print it as an integer percentage:
`%d pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s`.
Note the strict `>`.

### 8.2 `HetConc` and `HomIBS0` definitions — **[OPEN, inferred]**

The `.ibs` column list (`NHetHet`, `NHomHom`, `N_Het1`, `N_Het2`, `HetConc`, `Het2|1`,
`Het1|2`, `HomConc`) makes the intended algebra fairly clear. The natural reading, which
gives exactly 1.0 for a perfect duplicate:

```
  HetConc  =  N_HetHet / (N_Het1 + N_Het2 - N_HetHet)     [union / Jaccard form]
  Het2|1   =  N_HetHet / N_Het1                            [conditional]
  Het1|2   =  N_HetHet / N_Het2                            [conditional]
  HomConc  =  (N_HomHom - N_IBS0) / N_HomHom               [concordance among both-hom sites]
  HomIBS0  =  N_IBS0 / N_HomHom                            [IBS0 rate among both-hom sites]
  Concord  =  (N_IBS1 + N_IBS2) / N_SNP  or  N_IBS2 / N_SNP  -- ambiguous
```

The presence of both `HetConc` **and** the two conditionals `Het2|1`/`Het1|2` in `.ibs`
strongly implies `HetConc` is the symmetric (union) form, since otherwise it would
duplicate one of the conditionals.

**These are inferences from column names, not confirmed.** Confirm every one by running
the binary against a synthetic fileset with hand-computable counts.

### 8.3 Error-rate handling **[BINARY]**

An embedded message shows an error rate is used when building Mendelian-inconsistency
removal lists, printed at 4 decimals:
`Error rate is set at %.4lf to determine the MI removal list.`
Scope is `--build`/pedigree reconstruction, not the core kinship estimator.

---

## 9. Missing genotypes, QC, and adjustments

### 9.1 Missing data **[PAPER, verified]**

- **Pairwise-complete only.** Only markers non-missing in **both** members of the pair enter
  any count; `M_ij` is that count, and it varies per pair. This is the `N_SNP` column.
- `N_Aa^i`, `N_Aa^j` are recomputed **per pair** over `M_ij` (see §1). This is the most
  commonly botched detail in reimplementations.
- No imputation, no mean-filling. Missing genotypes are simply dropped from the pair.

### 9.2 Sample exclusion **[BINARY]**

Samples with too few non-missing SNPs are dropped before kinship, with messages:
```
The following %d samples are excluded from the kinship analysis (M<%d):
Sample FID=%s,IID=%s is removed with only %d non-missing SNPs.
Fam %s ID %s excluded for having %d non-missing SNPs
```
**[OPEN]** the numeric `M<%d` default threshold.

### 9.3 Inflation / deflation adjustments **[PAPER]**

KING applies **no post-hoc inflation correction** to the kinship estimate. Robustness is
structural, achieved by three design choices:

1. **Allele-frequency-free.** KING-robust never uses `p̂_m`, so stratification cannot bias it
   through mis-specified frequencies. (This is the whole point versus KING-homo and versus
   PLINK's `--genome`.)
2. **`min()` denominator across families** (§5.3) — a deliberate *deflation* guard against
   heterozygosity excess / individual-level HWE departure.
3. **Negative estimates are meaningful, not clipped** — an extreme negative `φ̂` indicates the
   pair is drawn from two distinct populations (Eq 10).

The only accuracy caveat the paper states is the residual, small bias for pairs that are
simultaneously related and cross-population.

### 9.4 Two-stage screening **[BINARY]**

For scalability across families KING screens on a SNP subset then confirms on all SNPs:
```
Stage 1 (with %d SNPs) screening ends at %s
Stage 2 (with all SNPs) inference ends at %s
Stages 1&2 (with %d SNPs): %lli pairs of relatives are detected (with kinship > %.4lf)
```
`--noscreen` disables it. **This is a performance optimisation and must not change
results** — a from-scratch implementation may compute all pairs directly and still match,
provided screening is disabled or is genuinely lossless in the reference. Verify by
comparing a `--noscreen` run against a default run on the same data.

Also **[PAPER]**: KING computes counts with word-parallel bit operations (AND/OR/XOR/NOT)
over packed 2-bit genotypes, reporting >60× speedup versus PLINK. **[BINARY]** confirms
`Autosome genotypes stored in %d words for each of %d individuals` and 64-bit words for the
X chromosome. This is an implementation strategy the reimplementation is free to adopt
independently — it is described in the paper, and popcount-based genotype counting is
standard practice.

---

## 10. Implementation checklist for byte-parity

1. Pack genotypes 2 bits/sample; derive per-pair `N_Aa^i`, `N_Aa^j`, `N_Aa,Aa`, `N_AA,aa`,
   `M_ij` with popcount over the pairwise-non-missing mask. Nothing else is needed for
   KING-robust.
2. `Kinship` (across families) `= 1/2 + (2·N_HetHet − 4·N_IBS0 − N_Aa^i − N_Aa^j) / (4·min(N_Aa^i, N_Aa^j))`.
3. `Kinship` (within family) `= (N_HetHet − 2·N_IBS0) / (N_Aa^i + N_Aa^j)`.
4. Reported `HetHet` and `IBS0` columns are **proportions** — divide the counts by `N_SNP`
   (`= M_ij`). **[MANUAL]** describes both as proportions.
5. Print data columns with `%.4f`; pedigree columns `Z0`/`Phi` with `%.3f`. Confirm bindings
   against a real run.
6. Do not clamp negative kinship.
7. Orientation invariance: `(x_i − x_j)²` and all four counts are invariant under
   `x → 2 − x`, so REF/ALT assignment does not affect any result. Use this as a property test.
8. Property tests to write first:
   - `SUM (x_i − x_j)² == N_IBS1 + 4·N_IBS0` on random genotype vectors.
   - within-family count form == squared-difference form, exactly, on random inputs.
   - within == between whenever `N_Aa^i == N_Aa^j`.
   - simulated MZ ⇒ 0.5; simulated PO ⇒ 0.25 with `IBS0 == 0`; simulated FS ⇒ 0.25 with
     `IBS0 ≈ 0.25 ×` unrelated rate; simulated unrelated ⇒ ~0.

---

## 11. Open items ranked by parity risk

| # | Item | Risk | How to close |
|---|---|---|---|
| 1 | Exact FS/PO IBS0 cutoff formula (§6.3) | **High** — changes `InfType` on every 1st-degree pair | Run binary on synthetic sets with varied MAF spectra; read the printed `%.4f` |
| 2 | Exact `.kin`/`.kin0` field separator, padding, column↔format binding (§7.3) | **High** — pure byte-parity | Run binary, hexdump output |
| 3 | `HetConc`, `HomIBS0`, `Concord` exact definitions (§8.2) | Medium | Synthetic fileset with hand-computable counts |
| 4 | Boundary inclusivity (`>` vs `>=`) at each kinship cutoff (§6.1) | Medium | Construct edge-case pairs |
| 5 | Sample-exclusion `M<%d` default (§9.2) | Low | `--help` / run with sparse data |
| 6 | Whether default two-stage screening is lossless (§9.4) | Low | Diff default vs `--noscreen` |
| 7 | `π̂₀` denominator exponent `2p²(1−p)²` (§4.4) | Low — derivation is solid | Simulation |

---

## 12. Sources

- [Manichaikul et al. 2010, PMC3025716](https://pmc.ncbi.nlm.nih.gov/articles/PMC3025716/)
- [Hail — Relatedness methods (`hl.king`)](https://hail.is/docs/0.2/methods/relatedness.html) — second witness for the count-form estimators
- [KING manual](https://www.kingrelatedness.com/manual.shtml) — output columns, cutoffs, `PropIBD`
- [SNPRelate `snpgdsIBDKING`](https://zhengxwen.github.io/SNPRelate/release/help/snpgdsIBDKING.html)
- Reference binary string constants + `--help`, KING 2.3.2 arm64
