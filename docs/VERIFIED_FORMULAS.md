# Verified formulas

Every formula on this page was **checked numerically against output produced by the
reference KING 2.3.2 binary**, not merely read out of the paper. Each one reproduces the
reference's printed value to the full 4 decimal places it prints, on both within-family
and between-family pairs, including a pair with genotyping error and a pair with
substantial missingness.

Reproduce the check with `tests/parity/verify_formulas.py`.

Anything **not** on this page is unverified and must be treated as a hypothesis until a
parity case confirms it. Unverified items are listed in [Open questions](#open-questions).

## Per-pair counts

All counts are over the set of SNPs where **both** members of the pair are non-missing.
Call that set size `M_ij` — this is the `N_SNP` column. Missing data is handled
**pairwise**: `N_SNP` differs from pair to pair.

| Symbol | Meaning |
| --- | --- |
| `M_ij` | SNPs non-missing in both `i` and `j` (printed as `N_SNP`) |
| `Het_i` | SNPs where `i` is heterozygous, counted over `M_ij` (printed as `N_Het1`) |
| `Het_j` | SNPs where `j` is heterozygous, counted over `M_ij` (printed as `N_Het2`) |
| `HetHet` | SNPs where **both** are heterozygous (printed as `NHetHet`) |
| `IBS0` | SNPs where the two are **opposite homozygotes** (printed as `N_IBS0`) |
| `HomHom` | SNPs where **both** are homozygous (printed as `NHomHom`) |

> **Critical parity detail.** `Het_i` and `Het_j` are counted over the *pairwise*
> non-missing set `M_ij`, **not** over each sample's own non-missing set. Using a
> per-sample heterozygote count computed once and reused across pairs gives subtly wrong
> kinship whenever missingness differs between samples. This was confirmed by the
> `missing` parity dataset.

Derived identities (exact, verified):

```
N_IBS1 = Het_i + Het_j - 2*HetHet
N_IBS2 = M_ij - IBS0 - N_IBS1
```

## Kinship

### Within-family (`.kin`, the `Kinship` column)

```
                HetHet - 2*IBS0
    phi_ij  =  -----------------
                  Het_i + Het_j
```

Sanity: for an MZ pair, `HetHet = Het_i = Het_j = H` and `IBS0 = 0`, giving `H/2H = 0.5`. ✓

### Between-family (`.kin0`, the `Kinship` column)

```
                     2*HetHet - 4*IBS0 - Het_i - Het_j
    phi_ij  =  0.5 + ---------------------------------
                          4 * min(Het_i, Het_j)
```

Sanity: for unrelated samples from one population, `Het_i = Het_j = H` and
`E[2*HetHet] = H`, giving `0.5 + (H - 2H)/(4H)`... which tends to 0. ✓

The two forms **coincide exactly when `Het_i == Het_j`**; they differ otherwise. The
`min()` denominator is what makes the between-family estimator robust to population
structure — it is the whole point of "KING-robust".

**Do not** use the within-family form for cross-family pairs or vice versa. Which form
applies is decided purely by whether the two samples share an `FID`.

### Verified against the reference

| Pair | Kind | Computed | Reference |
| --- | --- | --- | --- |
| `f1dad`/`f1kid1` | within | 0.252297 | `0.2523` |
| `f1dad`/`f1kid2` | within | 0.237633 | `0.2376` |
| `f1dad`/`f2dad` | between | 0.235215 | `0.2352` |
| `f1dad`/`f2mom` | between | −0.006831 | `-0.0068` |

## The `.ibs` / `.ibs0` derived columns

All verified to 4 decimals:

```
IBS      = (N_IBS1 + 2*N_IBS2) / M_ij          mean IBS allele sharing
Dist     = (N_IBS1 + 4*IBS0)   / M_ij          mean squared genotype distance
HetConc  = HetHet / (Het_i + Het_j - HetHet)   heterozygote Jaccard concordance
Het2|1   = HetHet / Het_i
Het1|2   = HetHet / Het_j
HomConc  = (HomHom - IBS0) / HomHom
```

Note `Dist` is **not** `2 - IBS`; that identity only holds when `IBS0 = 0`.

## The `.kin` / `.kin0` proportion columns

In `.kin` and `.kin0`, `HetHet` and `IBS0` are printed as **proportions of `N_SNP`**,
not as raw counts (the raw counts appear in `.ibs`):

```
HetHet_column = HetHet / M_ij
IBS0_column   = IBS0   / M_ij
```

Verified: `HetHet = 302`, `N_SNP = 1970` → `0.1533`; `IBS0 = 28`, `N_SNP = 1960` → `0.0143`.

## `--related`'s extra columns

`--related` widens `.kin` to sixteen columns and `.kin0` to fourteen. Four of the extra
ones (`IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`) come from the IBD-segment engine; the
other two are plain pairwise statistics:

```
HetConc_column = HetHet / (HetHet + IBS1)                       -- same as `.ibs`'s HetConc
HomIBS0_column = IBS0 / |{i hom-A1} u {j hom-A1}|               -- over M_ij
```

`HomIBS0` is undocumented and is **not** `IBS0 / HomHom` and **not** `1 - HomConc`: its
denominator is the number of variants of `M_ij` at which *either* sample is homozygous for
A1 — a union, not an intersection, and A1 specifically rather than "any homozygote".
Re-derived from raw `.bed` genotypes on the 727 within-family rows of `dups`, `multifam`,
`monomorphic` and `admixed` with zero mismatches. A pair with no A1 homozygote on either
side has a zero denominator and the reference prints `nan`.

`PropIBD = IBD2Seg + IBD1Seg / 2` in full `double` precision, not recomputed from the
printed 4-decimal columns: 87 of the corpus's 982 `.seg` rows disagree with the recomputed
value in the last digit.

## Pedigree-expected columns (`.kin` only)

`Z0` and `Phi` are **expected values derived from the declared pedigree**, not estimates:
`Phi` is the pedigree kinship coefficient and `Z0` the pedigree Pr[IBD = 0]. Observed:

| Pedigree relationship | `Z0` | `Phi` |
| --- | --- | --- |
| Parent–offspring | `0.000` | `0.2500` |
| Full siblings | `0.250` | `0.2500` |
| Unrelated within family | `1.000` | `0.0000` |

Formats: `Z0` is `%.3f`, `Phi` is `%.4f`.

### The `Error` column — graded, not a flag

`Error` measures disagreement between the **inferred** class and the
**pedigree-declared** class, and it is **not an integer**. It is graded:

| Value | Meaning |
| --- | --- |
| `0` | inferred class matches the pedigree |
| `0.5` | off by exactly one degree |
| `1` | off by more than one degree |

Confirmed by scanning every `.kin` in the golden corpus: the value set is exactly
`{0, 0.5, 1}`. An implementation that writes this column with `%d` is wrong — it will
print `0` for every half-step disagreement.

Worked example from the `tiny` capture: `f2dad`/`f2dup` (pedigree unrelated, inferred 2nd
degree) and `f2dup`/`f2kid1` (pedigree unrelated, inferred 3rd degree) are flagged, while
`f2dad`/`f2mom` and `f2dup`/`f2mom` both infer unrelated and are flagged `0`.

**`--kinship` and `--related` grade it differently, and their `Error` columns disagree on
nine corpus rows.**

* `--kinship` compares the kinship *estimate* against `Phi` on a multiplicative scale:
  within a factor of `sqrt(2)` is `0`, within a factor of 2 is `0.5`, beyond that `1`.
* `--related` compares the pedigree's relationship *label* — read off `(Phi, Z0)` on the
  usual `2^-(k+1/2)` bands, with `Z0 == 0` splitting the first degree into `PO` and `FS` —
  against `InfType`. An exact label match is `0`; otherwise `0.5` when the two degrees
  differ by exactly one **and both are 2nd degree or more distant**, and `1` in every other
  case. So a pedigree `PO` inferred `FS` scores `1` although both are first degree, while
  `2nd` against `3rd` scores `0.5`.

The `InfType` `--related` compares against is **not the one it prints**: the first `FS`
clause's `IBD2Seg >= 0.08` guard is dropped for the comparison, so a pair printed `2nd` can
be graded as `FS`. Both rules were validated over the 4 550 `InfType`-carrying rows of the
golden `.kin`/`.kin0` corpus with zero mismatches.

## Relationship inference cutoffs

Kinship thresholds, from the paper and the KING manual:

| Class | Kinship range |
| --- | --- |
| Duplicate / MZ twin | `> 0.354` |
| 1st degree | `0.177 – 0.354` |
| 2nd degree | `0.0884 – 0.177` |
| 3rd degree | `0.0442 – 0.0884` |
| 4th degree | `0.0221 – 0.0442` |
| Unrelated | `< 0.0221` |

The boundaries are `2^(-3/2)`, `2^(-5/2)`, `2^(-7/2)`, … i.e. successive halvings of the
kinship coefficient on a `2^(-k/2)` grid.

Parent–offspring is separated from full siblings **by IBS0, not by kinship** — both have
`phi = 0.25`. A true PO pair has `IBS0 ≈ 0` (no opposite-homozygote sites are possible
without genotyping error); a FS pair has `Pr[IBD=0] = 0.25` and hence a clearly non-zero
IBS0 rate. The reference binary uses a **data-derived threshold** on the IBS0 proportion
rather than a hard-coded constant — see [Open questions](#open-questions).

## Output-file field formats

| Column | Format |
| --- | --- |
| `N_SNP`, all `N_*` counts | `%d` |
| `Z0` | `%.3f` |
| `Phi`, `HetHet`, `IBS0`, `Kinship` | `%.4f` |
| `IBS`, `Dist`, `HetConc`, `Het2\|1`, `Het1\|2`, `HomConc` | `%.4f` |
| `Concord`, `HomConc`, `HetConc` in `.con` | `%.5f` |
| `Error` | **not** `%d` — see below |

Separators: `.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`, `unrelated.txt` are **tab**
separated. `bySample.txt` and `bySNP.txt` are **space** separated. This asymmetry is
real and verified.

## Row ordering

There are **two different rules**, and the between-family one is not what it looks like
on small inputs.

### Between-family (`.kin0`, `.ibs0`) — square-tiled, not row-major

Pairs are `i < j` over **`.fam` index order**, but they are emitted by a *block-tiled*
loop, sorted by the key:

```
(i / B, j / B, i, j)        integer division
```

with a **different block size per file**: `B = 32` for `.kin0`, `B = 8` for `.ibs0`.

Plain ascending `(i, j)` is **wrong**, and dangerously so: it coincides with the tiled
order whenever `n <= B`, so it looks correct on small fixtures and silently diverges on
real data. Verified against golden output by generating both candidate orderings and
comparing to the emitted rows:

| File | n | Rows | Uniquely matches |
| --- | --- | --- | --- |
| `dups` `.ibs0` | 10 | 43 | `B = 8` |
| `unrelated` `.ibs0` | 30 | 390 | `B = 8` |
| `bigish` `.ibs0` | 200 | 19,327 | `B = 8` |
| `bigish` `.kin0` | 200 | 19,327 | `B = 32` |
| `unrelated` `.kin0` | 30 | 390 | `B = 32`, 64 and plain all agree — `n <= 32`, so undiscriminating |

That last row is exactly the trap: at `n = 30` every candidate agrees.

### Within-family (`.kin`, `.ibs`) — natural sort on IID

Grouped by `FID` in first-appearance order; **within a family the members are re-sorted by
a natural (numeric-aware) sort of the IID, independent of `.fam` order**, and pairs are
the `i < j` upper triangle of that sorted list.

The `.kin` ordering was established with a discriminating experiment: the same `.bed`
was analysed with the family's members listed in several different `.fam` orders, and the
emitted row order never changed.

| IDs in family | `.fam` orders tried | `.kin` order emitted |
| --- | --- | --- |
| `alpha` `mike` `zeta` | `zeta,mike,alpha` and `alpha,zeta,mike` | `alpha` `mike` `zeta` |
| `a2` `a9` `a10` | `a10,a9,a2` | `a2` `a9` `a10` |
| `2` `9` `10` | `2,10,9` | `2` `9` `10` |
| `Dad` `kid` `mom` | `kid,Dad,mom` | `Dad` `kid` `mom` |

So the order is **natural sort** (numeric-aware): `a2 < a9 < a10` and `2 < 9 < 10`, which
plain lexicographic ordering would get wrong (`a10 < a2 < a9`, `10 < 2 < 9`). Using
`.fam` order — the obvious guess — is wrong, and using lexicographic order is also wrong.

> **Contested and re-settled.** A second independent pass over the corpus concluded this
> sort was *lexicographic*. It is not. That pass was under-powered: every sample ID in the
> corpus (`T_F`, `MZ_1`, `U1`, …) sorts identically under both rules, so the corpus cannot
> distinguish them. The discriminating case is a family whose IIDs are `2`, `9`, `10`:
> lexicographic predicts `10, 2, 9`, numeric predicts `2, 9, 10`, and the reference emits
> `(2,9) (2,10) (9,10)` — numeric. Reproduced under three different `.fam` orderings.

Pairs are then emitted as the `i < j` upper triangle over that sorted order.

### The exact comparator

"Natural sort" is a good approximation but not the rule. The exact comparator, fitted
against 44 probe families in `docs/BEHAVIOR.md`, walks the two IDs together:

* a **run of digits** compares against another digit run by **length first, then bytes**
  — so `7 < 70 < 007` (lengths 1, 2, 3), which is what the reference emits and what plain
  natural sort gets wrong (it would tie `007` with `7`);
* a **non-digit** character compares ASCII with **uppercase folding**;
* a non-digit sorts **before** a digit;
* a shorter ID that is a prefix of a longer one sorts first.

Because digit runs compare by length first, this coincides with numeric ordering for
IDs without leading zeros — which is why `2 < 9 < 10` and `a2 < a9 < a10` — while also
explaining the zero-padded case. The same comparator orders the `FID` blocks.

## Output files are not unconditional

An implementation can compute every number correctly and still fail the diff. Full detail
in `docs/BEHAVIOR.md`; the essentials:

### Truncation, not emptiness

The single-family `.kin` is **not "empty" as a rule — it is truncated to whole flushed
64 KiB chunks**, because the reference never closes the file. A zero-byte `.kin` is just
the small-data case of that: fewer than 64 KiB of rows were buffered, so nothing reached
disk. Holding the data fixed and padding the `FID` to push past the boundary changes the
row count while the byte count stays pinned near 65,536.

Every dataset in this corpus is small enough that the effect always presents as
zero-byte, so the corpus alone would have led to the wrong rule.

### Existence rules

| File | Created when |
| --- | --- |
| `.kin` | some family has ≥ 2 members (content then subject to the truncation above) |
| `.kin0` | ≥ 2 distinct `FID`s; for `--related`, additionally **N ≥ 100** |
| `.ibs` | always — header-only (139 bytes) when no within-family pair exists; never truncated |
| `.con` | always for N < 100; for N ≥ 100 only if a duplicate is found |

> An earlier note here claimed `--ibs` on `trio`/`nuclear` created neither `.ibs` nor
> `.ibs0`. **That did not reproduce** in a clean directory — both write a `.ibs`. The
> original observation came from a run whose fatal-error exit was masked by a redirect
> (see the A1/major-allele gate below). It has been removed rather than corrected.

### The `.ibs0` column set depends on the marker map

`MaxIBD2` and `Pr_IBD2` are appended to **both** `.ibs` and `.ibs0` iff the total
*usable* segment length is **≥ 100 Mb**, where a usable segment is a maximal run of SNPs
with a base-pair gap **≤ 156,250** (= 10 Mb / 64; 156,250 is usable, 156,251 is not). The
`.bim` cM column is ignored entirely. This is why the same 8-family fileset emits the long
header genome-wide and the short header when subset to chromosomes 1–2.

```
short: ... HetConc Het2|1 Het1|2 HomConc Kinship
long:  ... HetConc Het2|1 Het1|2 HomConc Kinship MaxIBD2 Pr_IBD2
```

Both are the literal `-9` unless `Kinship >= 2^-3.5`. `MaxIBD2` is the longest IBD2
segment in bp at `%.3f`; `Pr_IBD2` is a genome proportion at `%.4f`. `--related`'s `.kin`
and `.kin0` gain `IBD1Seg`/`IBD2Seg`/`PropIBD`/`InfType` on the same trigger.

## SNP inclusion — nothing is filtered

Every `.bim` record on chromosome `1`–`22`, `25` or `XY` enters the computation, **in file
order**. There is **no** monomorphic filter, **no** call-rate threshold, **no** MAF
threshold; `0`-coded alleles, duplicate IDs and duplicate positions are all kept. A SNP
called in only 2 of 8 samples still counts for that pair.

Note **chromosome 25 / `XY` is pooled with the autosomes**. Codes `23`/`X`, `24`/`Y` and
`26`/`MT` are held aside (X gets its own `<prefix>X.kin`), and **any other code — `0`,
`27`, `M`, `-1`, `chr1` — is dropped at map load** and is not counted in
`PLINK maps loaded`.

## `--autoQC` — the one pass that does filter

Everything above describes the relatedness path, which filters nothing. `--autoQC` is the
exception, and its rules are all verified against the reference. Full derivation, with the
experiment that fixes each rule, in `crates/open-king-cli/src/analysis/autoqc.rs`.

| step | test | applies to |
| --- | --- | --- |
| 1 | `missing / n > 1 - t1`, `t1 = min(0.8, 0.1 * (trunc(callrateM * 10) - 1))` | autosomes, X |
| 1 | monomorphic — every call the same homozygote (an all-het marker is **not**) | autosomes, X, Y |
| 2 | `missing / m > 1 - callrateN`, over the autosomes step 1 left | samples |
| 3 | `missing / n > 1 - callrateM`, over the samples step 2 left | autosomes, X |
| 4a / 6a | Y call rate in males, at `t1` then at `callrateM` | Y |
| 4b / 6b | X heterozygosity in males `> 5%` then `> 1%`, over **called** males | X |
| 4c / 6c | Y call rate in females `> 10%` then `> 2%` | Y |

`--callrateN` and `--callrateM` both default to **0.95**, and the default is applied inside
the pass, not in the option table: the banner prints `--callrateN` bare until the flag is
given, and an explicit `--callrateM 0` is honoured literally.

The four report files, all named by **concatenation** (`--prefix ZZ_` gives
`ZZ__autoQC_Summary.txt`, double underscore). The three `*toberemoved*` / `updatesex` files
are tab separated; `_autoQC_Summary.txt` is fixed-width and space-padded and contains no tab
at all (measured: `grep -c $'\t'` returns 0):

| file | header | rows |
| --- | --- | --- |
| `<p>_autoQC_snptoberemoved.txt` | `SNP\tREASON` | one per removed marker, **grouped by the step that removed it**, autosomes then X then Y within a group, `.bim` order within a class. Reasons: `CallRateLessThan<pct>`, `Monomorphic`, `xHeterozygosityInMale`, `YSNPInFemales` — `<pct>` is `round(100 * callrateM)` for *both* Y call-rate filters, even step 4a, which applies `t1` |
| `<p>_autoQC_sampletoberemoved.txt` | `FID\tIID\tREASON` | one per removed sample, grouped by check: `MissingMoreThan<pct>` (`round(100 * (1 - callrateN))`), then `MislabeledAsMale`, `MislabeledAsFemale`, `GenderQC`, `.fam` order within a group |
| `<p>_autoQC_updatesex.txt` | none | `FID\tIID\tsex` for every sample whose `.fam` sex was 0, `.fam` order. **Not created** when there is none |
| `<p>_autoQC_Summary.txt` | `Step Description Subjects SNPs` | the step table, byte-identical to the block on stdout. `%-5s%-55s%-10s%-10s` for the two count rows, `%-5s%-65s(%d)` for a SNP counter and `%-5s%-55s(%d)` for a subject counter — the step-2.x block only when the map has both X and Y markers |

The summary's `1.1`/`1.3`/`1.4` labels are **fixed text**: they say `< 80%` and `< 95%`
whatever `--callrateN`/`--callrateM` were, while the console lines above them print the
thresholds actually applied.

Three traps, none of them derivable from the manual:

* **The tests are on the missing rate, not the call rate** — `missing / n > 1 - t`. At
  exact boundaries the two differ: a marker missing in 2 of 20 samples is removed under
  `--callrateM 0.9` (because `1 - 0.9` is `0.09999999999999998`) while a marker at exactly
  95 % survives the default.
* **The dropped word.** The pass packs 16 markers to a word — it reports `ceil(m/16)` words
  where every other analysis reports `ceil(m/64)` — and scans `16 * (words - 1) + m % 16`
  of them, so a class whose marker count is a multiple of 16 never looks at its last 16
  markers and never counts them in the final tally. Autosomes, X and Y each drop their own;
  MT is never scanned at all, and is also never filtered. The step-5 per-sample statistics
  are the exception that proves it: they walk the dropped word and count its 16 slots as
  called and homozygous.
* **The console headers do not agree with the data.** They print the *untruncated* count
  (`step 2 … (with 5232 SNPs)` on `missing` counts 16 markers the pass cannot see), and
  step 4's header prints the Y count from before step 1's removals while step 6's prints
  the current one.

## The A1/major-allele fatal gate

`--related`, `--ibs` and `--build` (but **not** `--kinship`) abort with

```
FATAL ERROR - Too many first alleles as the major allele (~X%)
```

when A1 is not the observed minor allele often enough. Any synthetic fixture must
re-orient alleles per SNP after drawing genotypes, exactly as `--make-bed` does.
`generate_corpus.py` does this; ad-hoc `p = 0.5` fixtures do not, and the resulting fatal
exit is easy to mistake for "the file was not created".

## Open questions

Questions 2–8 were attacked empirically; the experiments, raw output and resulting rules
live in **[BEHAVIOR.md](BEHAVIOR.md)**. Status summary below.

1. ~~**`.kin` within-family row order.**~~ **RESOLVED** — natural sort on sample ID,
   independent of `.fam` order. See [Row ordering](#row-ordering).
2. **PO vs FS IBS0 threshold.** **PARTIALLY RESOLVED** —
   [BEHAVIOR.md § Q2](BEHAVIOR.md#q2--parentoffspring-vs-full-sibling-discrimination).
   The message is `Cutoff value for IBS0 between FS and PO is set at %.4f`, printed by
   `--build`/`--cluster`/`--unrelated` and only when IBD-segment analysis is unavailable;
   when segments exist, PO/FS is decided by IBD2 sharing instead. The *application* is
   established (1st-degree pair is PO iff its IBS0 proportion is below the cutoff, a pair
   exactly on the printed value counting as PO). The *value* is still open: it is
   deterministic in the data but was shown **not** to be fitted to the observed
   1st-degree IBS0 distribution, nor a function of `N`, SNP count, sample order, relative
   count, or any single allele-frequency summary tried. It sits in `[0.0035, 0.0060]`,
   is `0.0055` on ordinary data, and never reaches an output file.
3. **SNP inclusion rules.** **RESOLVED** —
   [BEHAVIOR.md § Q3](BEHAVIOR.md#q3--snp-inclusion-rules). Every `.bim` record on
   chromosome `1`–`22`, `25` or `XY` is used, in file order. Monomorphic SNPs, `0`
   alleles, low call rate, low MAF, duplicate IDs and duplicate positions are **all**
   retained; `23`/`X`, `24`/`Y`, `26`/`MT` are held aside; any other chromosome code is
   dropped at map load. Missingness is pairwise only.
4. **`--cpus` determinism.** **RESOLVED** —
   [BEHAVIOR.md § Q4](BEHAVIOR.md#q4----cpus-determinism). Every output file is
   byte-identical across `--cpus 1/2/4/8` on a 200 × 50 000 dataset; only stdout progress
   percentages differ.
5. **`--degree` filtering semantics.** **RESOLVED** —
   [BEHAVIOR.md § Q5](BEHAVIOR.md#q5----degree-semantics). `--degree d` filters `.kin0`
   only (never `.kin`, never `.ibs`/`.ibs0`), on the kinship **estimate**, keeping pairs
   with `kinship >= 2^-(d+1.5)` compared against the exact double rather than the printed
   `%.5lf`. `--degree 0` means unset. Negative degrees behave inconsistently between
   `--kinship` and `--related` and remain unspecified.
6. **Zero-padded numeric ID ordering.** **RESOLVED** —
   [BEHAVIOR.md § Q6](BEHAVIOR.md#q6--the-sample-id-sort-comparator). No integer parse is
   involved. Chunked comparison: digit runs by (length, then bytes), non-digits one
   character at a time ASCII-uppercase-folded, a non-digit sorting before a digit, and a
   shorter prefix first. The same comparator orders the FID blocks. IDs are
   case-insensitively unique and an IID of `0` is rejected. Rust implementation given in
   BEHAVIOR.md.
7. **Output-file existence.** **RESOLVED** —
   [BEHAVIOR.md § Q7](BEHAVIOR.md#q7--output-file-existence). Includes the correction that
   the single-family `.kin` is not "empty" but *truncated to flushed 64 KiB chunks* — a
   zero-byte file is just the small-data case of that bug.
8. **`.ibs`/`.ibs0` column-set variation.** **RESOLVED** —
   [BEHAVIOR.md § Q8](BEHAVIOR.md#q8--ibs--ibs0-column-set-variation). `MaxIBD2` and
   `Pr_IBD2` appear on both files iff the total usable IBD-segment length is ≥ 100 Mb,
   where a usable segment is a maximal run of SNPs with base-pair gaps ≤ 156 250.
   `--related`'s `.kin0` gains its `IBD1Seg/IBD2Seg/PropIBD/InfType` block on the same
   trigger.

Newly opened by that work, and recorded in
[BEHAVIOR.md § Side findings](BEHAVIOR.md#side-findings): the `.kin` `Error` column is
**not** `%d` (it takes the value `0.5`), and the reference aborts with
`Too many first alleles as the major allele` when too many `.bim` A1 alleles are the major
allele.
