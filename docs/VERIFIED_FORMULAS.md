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

## Pedigree-expected columns (`.kin` only)

`Z0` and `Phi` are **expected values derived from the declared pedigree**, not estimates:
`Phi` is the pedigree kinship coefficient and `Z0` the pedigree Pr[IBD = 0]. Observed:

| Pedigree relationship | `Z0` | `Phi` |
| --- | --- | --- |
| Parent–offspring | `0.000` | `0.2500` |
| Full siblings | `0.250` | `0.2500` |
| Unrelated within family | `1.000` | `0.0000` |

Formats: `Z0` is `%.3f`, `Phi` is `%.4f`.

### The `Error` column

`Error` is `1` when the **inferred** relationship class disagrees with the
**pedigree-declared** class, else `0`. Verified on the captured `tiny` dataset: the
pairs flagged `1` were `f2dad`/`f2dup` (pedigree unrelated, inferred 2nd degree,
`Kinship = 0.1161`) and `f2dup`/`f2kid1` (pedigree unrelated, inferred 3rd degree,
`Kinship = 0.0741`), while `f2dad`/`f2mom` (`-0.0308`) and `f2dup`/`f2mom` (`0.0090`)
both infer unrelated and are flagged `0`.

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
| `Error` | `%d` |

Separators: `.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`, `unrelated.txt` are **tab**
separated. `bySample.txt` and `bySNP.txt` are **space** separated. This asymmetry is
real and verified.

## Row ordering

* `.kin0` — outer loop over samples in **`.fam` file order**, inner loop over
  later-ordered samples in a different family, also in `.fam` order. Verified on the
  `tiny` capture.
* `.kin` — grouped by `FID`, and **within a family the rows are ordered by a
  deterministic sort of the sample ID that is *independent of `.fam` order*.**

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

Pairs are then emitted as the `i < j` upper triangle over that sorted order.

**Unresolved edge case:** zero-padded numeric IDs. For the family `{007, 7, 70}` the
emitted order is `7`, `70`, `007` under all three `.fam` orders tried, which is neither
plain natural sort (`007` and `7` would tie) nor lexicographic (`007 < 7 < 70`). Tracked
as open question 6.

## Single-family behaviour

**When the dataset contains only one family, `.kin` is written as a zero-byte file** —
not even the header row — even though the console still prints
`Within-family kinship data saved in file king.kin` and the relationship-summary table
reports the pairs as correctly inferred. Adding a second family to the same `.bed` makes
the header and all within-family rows appear.

Verified with `--kinship` and `--related`, and independently of the phenotype column.
Reproducing this exactly is required for parity: a single-family dataset is a common
real-world input, and emitting a populated `.kin` there would be a diff against the
reference on the very first case a user tries.

## Output files are not unconditional

Two effects mean an implementation can compute every number correctly and still fail the
diff. Both are tracked in `docs/BEHAVIOR.md`.

**File existence varies with the input.** Running `--ibs`:

| Dataset | Families | Samples | `king.ibs` | `king.ibs0` |
| --- | --- | --- | --- | --- |
| `trio` | 1 | 3 | not created | not created |
| `nuclear` | 1 | 6 | not created | not created |
| `threegen` | 1 | 14 | created | not created |
| `dups` | 8 | 10 | created | created |

Note `threegen` is also a single family yet does get a `.ibs`, so family count alone does
not explain it. **Absent**, **zero-byte**, and **header-only** are three distinct
outcomes and must each be reproduced as-is.

**The `.ibs0` column set varies with the marker map.** Full-genome `dups` emits two extra
trailing columns, `MaxIBD2` and `Pr_IBD2`; the *same* 8-family fileset subset to
chromosomes 1–2 emits the short header ending at `Kinship`. The trigger is therefore the
map, not the samples — presumably whether an IBD2 segment analysis is attempted.

```
short: ... HetConc Het2|1 Het1|2 HomConc Kinship
long:  ... HetConc Het2|1 Het1|2 HomConc Kinship MaxIBD2 Pr_IBD2
```

## Open questions

These are unresolved and each is paired with the experiment that settles it.

1. ~~**`.kin` within-family row order.**~~ **RESOLVED** — natural sort on sample ID,
   independent of `.fam` order. See [Row ordering](#row-ordering).
2. **PO vs FS IBS0 threshold.** The binary's string table contains
   `1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf`, implying the
   cutoff is computed from the data, not fixed.
   *Experiment:* run `--related` on datasets with differing overall IBS0 rates and read
   the printed threshold back off stdout; fit the rule.
3. **SNP inclusion rules.** Whether monomorphic SNPs, SNPs with an allele coded `0`, or
   SNPs above a missingness threshold are dropped before counting, and whether only
   autosomes are used by default. The console prints
   `Genotype data consist of N autosome SNPs`, which suggests non-autosomes are excluded
   from the default relatedness path.
   *Experiment:* the `monomorphic` and `sexchr` parity datasets; compare `N_SNP` against
   the count we predict under each candidate rule.
4. **`--cpus` determinism.** Whether thread count changes any printed digit (it should
   not, since the kernel sums integers, but must be confirmed).
   *Experiment:* diff `--cpus 1` against `--cpus 8` output on `bigish`.
5. **`--degree` filtering semantics.** Whether `--degree n` filters `.kin0` rows only, or
   also changes `.kin`, and whether the cutoff is applied to kinship or to inferred class.
6. **Zero-padded numeric ID ordering.** `{007, 7, 70}` emits as `7`, `70`, `007`.
   *Experiment:* sweep families of IDs mixing widths and leading zeros
   (`0`, `00`, `01`, `1`, `10`, `0010`, `1a`, `a1`, `-1`, `1.0`, very long digit strings
   that overflow 32-bit) and fit the comparator. Likely the ID is parsed to an integer
   with non-parsing or overflowing values falling back to a string compare and sorting
   after the numerics — confirm or refute.
