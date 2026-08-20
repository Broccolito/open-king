# Interpreting the output

How to read the numbers `open-king` writes, and how to avoid reading them wrongly.

This page is for the person holding a `.kin`, `.kin0`, `.seg` or `.ibs` file who has to
decide what it means. It assumes you know what a kinship coefficient is. It does **not**
assume you know this codebase.

* Command surface, flags and file formats → [`SPEC.md`](SPEC.md)
* Exact formulas and how each was verified → [`VERIFIED_FORMULAS.md`](VERIFIED_FORMULAS.md)
* Undocumented behaviours of KING 2.3.2 that this program reproduces →
  [`BEHAVIOR.md`](BEHAVIOR.md)
* Byte-level agreement with the reference binary → [`PARITY.md`](PARITY.md)

**Every number on this page comes from a real run.** File excerpts are pasted verbatim; where
a block is an aggregate over many rows, or has had columns trimmed to fit, it says so. No
output here is invented or predicted. Every command is reproducible from a clean checkout:

```
cd /path/to/open-king
cargo build --release                                    # -> target/release/open-king
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
```

That writes 13 synthetic PLINK filesets whose true pedigrees are known by construction —
which is what makes it possible to say, below, that a given call is *wrong*. A handful of
examples need a fileset reshaped (thinned, MAF-filtered, error-injected); the short script
that does it is in the [appendix](#appendix--reshapepy).

---

## Contents

* [1. What the kinship coefficient is](#1-what-the-kinship-coefficient-is)
* [2. The two estimators](#2-the-two-estimators)
* [3. The relationship cutoffs](#3-the-relationship-cutoffs)
* [4. Parent–offspring vs full siblings](#4-parentoffspring-vs-full-siblings)
* [5. IBD1Seg, IBD2Seg and PropIBD](#5-ibd1seg-ibd2seg-and-propibd)
* [6. The Error column](#6-the-error-column)
* [7. Pitfalls](#7-pitfalls)
* [8. What relatedness inference cannot tell you](#8-what-relatedness-inference-cannot-tell-you)
* [9. A reading checklist](#9-a-reading-checklist)
* [Appendix — reshape.py](#appendix--reshapepy)

---

## 1. What the kinship coefficient is

The kinship coefficient `φ_ij` is the probability that an allele drawn at random from
individual `i` and an allele drawn at random from individual `j`, at the same autosomal
locus, are identical by descent.

| Relationship | `φ` | Degree |
| --- | --- | --- |
| self / MZ twin / duplicate sample | 0.5 | — |
| parent–offspring, full siblings | 0.25 | 1st |
| half sibs, grandparent–grandchild, avuncular | 0.125 | 2nd |
| first cousins, half-avuncular | 0.0625 | 3rd |
| first cousins once removed | 0.03125 | 4th |
| unrelated | 0 | — |

The `Kinship` column is the **KING-robust** estimator of Manichaikul A, Mychaleckyj JC,
Rich SS, Daly K, Sale M, Chen W-M. *Robust relationship inference in genome-wide
association studies.* **Bioinformatics 2010;26(22):2867–2873** (PMC3025716), equations (9)
and (11).

Three properties follow from how it is defined, and all three matter when reading output:

1. **It uses no allele frequencies.** Not sample frequencies, not a reference panel.
   Everything comes from per-pair genotype counts. That is what "robust" refers to: the
   estimator does not inherit the bias that a mis-specified allele-frequency vector injects
   into frequency-based estimators.
2. **It is computed pairwise.** Every count is taken over the SNPs where *both* members of
   the pair are non-missing. There is no per-sample precomputation you can reuse; two rows
   in the same file can rest on very different marker sets. The `N_SNP` column tells you
   which.
3. **It is a ratio, and it is unbounded below.** Negative values are ordinary output and
   carry information — see [§7.3](#73-population-structure).

### The counts

Fix a pair `i`, `j`. Let `M_ij` be the SNPs called in both. Over exactly that set:

| Symbol | Column in `.ibs` / `.ibs0` | Meaning |
| --- | --- | --- |
| `M_ij` | `N_SNP` | SNPs non-missing in both |
| `Het_i`, `Het_j` | `N_Het1`, `N_Het2` | SNPs at which `i` (resp. `j`) is heterozygous |
| `HetHet` | `NHetHet` | SNPs at which **both** are heterozygous |
| `IBS0` | `N_IBS0` | SNPs at which the two are **opposite homozygotes** |

In `.kin` and `.kin0` the `HetHet` and `IBS0` columns are printed as **proportions of
`N_SNP`**, not raw counts. The raw counts are in `.ibs` / `.ibs0`. This trips people up
when they try to reconcile the two files by hand.

### Before you read anything: the file may not be there

Three existence rules surprise people, and all three are faithful reproductions of KING
2.3.2 rather than bugs in this program:

* **A single-family dataset produces a zero-byte `.kin`.** The run reports the numbers on
  stdout and then writes nothing:

  ```
  $ open-king -b /tmp/kingdocs/threegen.bed --kinship --prefix tg
  Within-family kinship data saved in file tg.kin

  Relationship summary (total relatives: 39 by pedigree, 38 by inference)
    Source	MZ	PO	FS	2nd	3rd	OTHER
    ===========================================================
    Pedigree	0	14	3	14	8	27
    Inference	0	14	3	12	9	28

  There is only one family.

  $ wc -c tg.kin
         0 tg.kin
  ```

  The general rule is that `.kin` is truncated to whole flushed 64 KiB chunks; on small
  single-family inputs that rounds down to nothing. If you need the `.kin` numbers for a
  one-family dataset, split the FIDs (see [§2](#2-the-two-estimators)) and read `.kin0`, or
  use `--ibs`, which is never truncated.
* **`--related` writes `.kin0` only when N ≥ 100.** `--kinship` has no such gate. On a
  20-sample fileset, `--related` gives you a `.kin` and no `.kin0` at all.
* **`.ibs` is always created**, header-only (139 bytes) if no within-family pair exists.

Full rules and the experiments behind them: [`BEHAVIOR.md` § Q7](BEHAVIOR.md#q7--output-file-existence).

---

## 2. The two estimators

There are two, they give different answers, and which one you get is decided by one thing:
**whether the two samples share a FID in the `.fam` file.**

### Within-family — used for pairs sharing a FID, written to `.kin` and `.ibs`

```
                HetHet - 2*IBS0
    phi_ij  =  -----------------
                  Het_i + Het_j
```

The denominator is (twice) the **average** heterozygosity of the pair. That is the right
choice when you already know the two individuals come from the same ancestry — which is
what a shared FID asserts.

### Between-family — used for pairs in different FIDs, written to `.kin0` and `.ibs0`

```
                     2*HetHet - 4*IBS0 - Het_i - Het_j
    phi_ij  =  0.5 + ---------------------------------
                          4 * min(Het_i, Het_j)
```

The denominator uses the **minimum** heterozygosity of the pair.

### The arithmetic, on real rows

```
$ open-king -b /tmp/kingdocs/admixed.bed --ibs --prefix fam
$ awk '{$1=$2; print}' /tmp/kingdocs/admixed.fam > allsplit.fam   # every sample its own FID
$ open-king -b /tmp/kingdocs/admixed.bed --fam allsplit.fam --ibs --prefix nofam
```

```
within-family  X_C1 x X_C2: Het_i=6568 Het_j=6590 HetHet=4114 IBS0=263
   (HetHet - 2*IBS0)/(Het_i+Het_j) = (4114 - 2*263)/(6568+6590) = 0.272686   file prints 0.2727
between-family Y_F x Y_C1 : Het_i=6640 Het_j=7364 HetHet=3293 IBS0=0  min=6640
   0.5 + (2*3293 - 4*0 - 6640 - 7364)/(4*6640) = 0.220708   file prints 0.2207
   the same counts under the within-family form: 0.235147
```

### Why `min(Het_i, Het_j)`

The average-heterozygosity denominator is only meaningful if the two heterozygosities are
measuring the same thing. When `i` and `j` have different ancestry their heterozygosities
differ, the average is inflated relative to either one, and the estimate is biased.
Taking the smaller of the two is the conservative choice: it can only pull the estimate
down, never invent relatedness. It also guards against individual-level departures from
Hardy–Weinberg.

The two forms are **algebraically identical when `Het_i == Het_j`**, and they diverge in
proportion to the gap. That is directly measurable. The run above analyses one fileset
twice — once with its declared families, once with every sample in its own FID — so the
counts are byte-for-byte the same and the *only* difference is which formula was applied:

```
pair            Het_i  Het_j  gap    HetHet  IBS0   .kin/.ibs   .kin0/.ibs0   delta
X_C1   X_C2     6568   6590   0.3%    4114   263    0.2727        0.2723  -0.0004
X_C1   X_F      6568   6607   0.6%    3252     0    0.2468        0.2461  -0.0007
X_C1   X_M      6568   6635   1.0%    3222     0    0.2440        0.2427  -0.0013
X_C2   X_F      6590   6607   0.3%    3231     0    0.2448        0.2445  -0.0003
X_C2   X_M      6590   6635   0.7%    3303     0    0.2498        0.2489  -0.0009
X_F    X_M      6607   6635   0.4%    2700  1396   -0.0069       -0.0080  -0.0011
Y_C1   Y_C2     7364   7464   1.4%    4389   448    0.2356        0.2338  -0.0018
Y_C1   Y_F      7364   6640  10.9%    3293     0    0.2351        0.2207  -0.0144
Y_C1   Y_M      7364   6666  10.5%    3331     0    0.2374        0.2237  -0.0137
Y_C2   Y_F      7464   6640  12.4%    3369     0    0.2389        0.2227  -0.0162
Y_C2   Y_M      7464   6666  12.0%    3368     0    0.2384        0.2227  -0.0157
Y_F    Y_M      6640   6666   0.4%    2423  1959   -0.1124       -0.1136  -0.0012
Z_C1   Z_C2     6568   6619   0.8%    3830   510    0.2131        0.2120  -0.0011
Z_C1   Z_F      6568   6520   0.7%    3237     0    0.2473        0.2464  -0.0009
Z_C1   Z_M      6568   6556   0.2%    3322     0    0.2531        0.2529  -0.0002
Z_C2   Z_F      6619   6520   1.5%    3312     0    0.2521        0.2502  -0.0019
Z_C2   Z_M      6619   6556   1.0%    3295     0    0.2501        0.2489  -0.0012
Z_F    Z_M      6520   6556   0.6%    2561  1312   -0.0048       -0.0062  -0.0014
```

`admixed` holds two populations at F<sub>ST</sub> 0.10. Families X and Z are within one
population; family Y is a cross-population mating, so its children (`Y_C1`, `Y_C2`) are
~11 % more heterozygous than either parent. Where the gap is under 2 % the two estimators
agree to 0.002; on the four cross-ancestry parent–offspring pairs, where the gap is 11–12 %,
they differ by up to 0.016 and the between-family form is always the lower.

### What this means for you

* **A `.kin` number and a `.kin0` number are not the same statistic.** Do not pool them
  into one distribution, one histogram, or one threshold sweep without checking that
  heterozygosity is homogeneous across your samples.
* **Your `.fam` FIDs are an analysis choice, not metadata.** Declaring an ancestrally
  heterogeneous set of people as one "family" silently switches them to the less
  conservative estimator. If you do not have real pedigree families, give every sample its
  own FID and read `.kin0`.
* Never apply the within-family form across populations by hand. See
  [§7.6](#76-using-the-within-family-estimator-across-populations).

---

## 3. The relationship cutoffs

| Class | Kinship range |
| --- | --- |
| Duplicate / MZ twin | `> 0.354` |
| 1st degree | `0.177 – 0.354` |
| 2nd degree | `0.0884 – 0.177` |
| 3rd degree | `0.0442 – 0.0884` |
| 4th degree | `0.0221 – 0.0442` |
| Unrelated | `< 0.0221` |

These are not round numbers chosen for convenience. A degree-`k` relative has expected
`φ = 2^-(k+1)`, and each boundary is the **geometric midpoint** of two adjacent expectations,
`sqrt(2^-(k+1) · 2^-(k+2)) = 2^-(k+3/2)`:

```
$ python3 -c "
for k in range(1,6):
    e=-(k+0.5); print('2^%-5s = %.10f  ->  %s' % (e, 2**e, ('0.354','0.177','0.0884','0.0442','0.0221')[k-1]))"
2^-1.5  = 0.3535533906  ->  0.354
2^-2.5  = 0.1767766953  ->  0.177
2^-3.5  = 0.0883883476  ->  0.0884
2^-4.5  = 0.0441941738  ->  0.0442
2^-5.5  = 0.0220970869  ->  0.0221
```

Three consequences worth internalising:

* **The bands halve.** The 1st-degree band is 0.177 wide; the 4th-degree band is 0.022 wide.
  A fixed absolute error of, say, ±0.02 in the estimate is negligible for a sibling and
  fatal for a first-cousin-once-removed. Most classification mistakes happen at 3rd and 4th
  degree, and that is arithmetic, not a defect.
* **`--degree d` filters on the same grid.** `--degree 2` keeps `.kin0` pairs with
  kinship `>= 2^-3.5`, compared against the full-precision double, not the printed value.
* **Negative kinship is real output.** It is not an error and must not be clamped to zero;
  it means the pair is *less* allele-sharing than the estimator's zero point, which is
  usually ancestry divergence ([§7.3](#73-population-structure)).

`InfType`, where it appears, uses the labels `Dup/MZ`, `PO`, `FS`, `2nd`, `3rd`, `4th`,
`UN`.

---

## 4. Parent–offspring vs full siblings

Both have `φ = 0.25`. **Kinship cannot separate them, and does not try.** IBS0 does.

A true parent–offspring pair shares one allele IBD at every locus, so an opposite-homozygote
site (`AA` vs `aa`) is genetically impossible: `IBS0 = 0` exactly. Full siblings have
`Pr[IBD = 0] = 0.25`, so a quarter of their genome can produce opposite homozygotes and the
IBS0 rate is clearly non-zero.

```
$ open-king -b /tmp/kingdocs/multifam.bed --kinship --prefix mfk
```

`multifam` declares four families and contains undeclared cross-family relatives. Within
FAM1, `A_F`/`A_M` are the parents and `A_C1..A_C3` the children:

```
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	Kinship	Error
FAM1	A_C1	A_C2	15000	0.250	0.2500	0.2219	0.0153	0.2721	0
FAM1	A_C1	A_C3	15000	0.250	0.2500	0.2171	0.0214	0.2476	0
FAM1	A_C1	A_F	15000	0.000	0.2500	0.1735	0.0000	0.2475	0
FAM1	A_C1	A_M	15000	0.000	0.2500	0.1717	0.0000	0.2444	0
FAM1	A_C2	A_C3	15000	0.250	0.2500	0.2251	0.0124	0.2837	0
FAM1	A_C2	A_F	15000	0.000	0.2500	0.1773	0.0000	0.2523	0
FAM1	A_C2	A_M	15000	0.000	0.2500	0.1751	0.0000	0.2485	0
FAM1	A_C3	A_F	15000	0.000	0.2500	0.1765	0.0000	0.2508	0
FAM1	A_C3	A_M	15000	0.000	0.2500	0.1753	0.0000	0.2484	0
FAM1	A_F	A_M	15000	1.000	0.0000	0.1368	0.0693	-0.0027	0
```

Every sib pair and every parent–child pair sits in the 1st-degree band. The `IBS0` column
splits them cleanly: `0.0000` for the six parent–child pairs, `0.0124`–`0.0214` for the three
sib pairs. The unrelated couple `A_F`/`A_M` has `IBS0 = 0.0693`.

The same separation appears across families. `A_F` is genetically a full sib of `B_F` and a
parent of `C_F`, neither declared:

```
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	Kinship
FAM1	A_F	FAM2	B_F	15000	0.2169	0.0204	0.2501
FAM1	A_F	FAM2	B_M	15000	0.1386	0.0685	0.0001
FAM1	A_F	FAM2	B_C1	15000	0.1536	0.0301	0.1319
FAM1	A_F	FAM2	B_C2	15000	0.1533	0.0402	0.1039
FAM1	A_F	FAM2	B_C3	15000	0.1589	0.0275	0.1451
FAM1	A_F	FAM3	C_F	15000	0.1772	0.0000	0.2526
FAM1	A_F	FAM3	C_M	15000	0.1377	0.0719	-0.0118
```

Kinship 0.2501 with IBS0 0.0204 → siblings. Kinship 0.2526 with IBS0 0.0000 → parent and
child. Same kinship, different relationship, and only IBS0 says which.

**The threshold.** There is no fixed IBS0 cutoff, and none is written to any output file.
When IBD-segment analysis is available — which it is for every corpus example on this page —
PO and FS are decided by IBD2 sharing instead ([§5](#5-ibd1seg-ibd2seg-and-propibd)), and
that is the path you will normally be on. When it is not, the family-clustering commands
(`--build`, `--cluster`, `--unrelated`) fall back to a cutoff on the IBS0 proportion,
computed from the data, reported on stdout as
`Cutoff value for IBS0 between FS and PO is set at %.4f`, and applied so that a 1st-degree
pair is PO iff its IBS0 proportion is at or below it. On ordinary genome-wide array data it
lands near 0.0055. How it is applied is established; how the value itself is derived is not
— see
[`BEHAVIOR.md` § Q2](BEHAVIOR.md#q2--parentoffspring-vs-full-sibling-discrimination).

The number 0.0055 is used as a reference point in the pitfalls below because it is the value
the reference binary produces on ordinary data. Treat it as a scale, not a constant.

Remember what that cutoff assumes: that a real PO pair produces essentially no opposite
homozygotes. Genotyping error ([§7.4](#74-genotyping-error)) and rare-variant panels
([§7.5](#75-ascertained-or-filtered-snp-panels)) both break that assumption, in opposite
directions.

---

## 5. IBD1Seg, IBD2Seg and PropIBD

`--ibdseg`, and `--related` on data with sufficient marker density, add a second, independent
line of evidence: the genome is scanned for **contiguous IBD segments**, and the results are
reported as genome proportions.

| Column | Meaning |
| --- | --- |
| `IBD1Seg` | proportion of the analysable genome shared IBD on **exactly one** of the two chromosomes |
| `IBD2Seg` | proportion shared IBD on **both** chromosomes |
| `PropIBD` | `IBD2Seg + IBD1Seg / 2` — a segment-based relatedness proportion, ≈ `2 × φ` |
| `InfType` | the relationship called from those three |

`PropIBD` is computed in full double precision from the underlying values, not from the
printed 4-decimal columns; 87 of the corpus's 982 `.seg` rows disagree with the naive
recomputation in the last digit.

```
$ open-king -b /tmp/kingdocs/threegen.bed --ibdseg --prefix tgs
Total length of 21 chromosomal segments usable for IBD segment analysis is 982.7 Mb.
```

`threegen` is one three-generation pedigree. Excerpt from `tgs.seg` (`TG_GF`/`TG_GM1` are
grandparents, `TG_P*` their children, `TG_C*` their grandchildren):

```
FID1	ID1	FID2	ID2	IBD1Seg	IBD2Seg	PropIBD	InfType
TG	TG_GF	TG	TG_P1	1.0000	0.0000	0.5000	PO
TG	TG_GF	TG	TG_C1	0.5543	0.0000	0.2772	2nd
TG	TG_GF	TG	TG_C2	0.3433	0.0000	0.1716	3rd
TG	TG_P1	TG	TG_P2	0.3336	0.3479	0.5147	FS
TG	TG_P1	TG	TG_P3	0.4571	0.0000	0.2286	2nd
TG	TG_P1	TG	TG_C1	1.0000	0.0000	0.5000	PO
TG	TG_P2	TG	TG_C1	0.4339	0.0000	0.2170	2nd
TG	TG_P3	TG	TG_C4	0.1040	0.0000	0.0520	4th
TG	TG_C1	TG	TG_C2	0.4298	0.3389	0.5538	FS
TG	TG_C1	TG	TG_C4	0.1559	0.0000	0.0780	4th
TG	TG_C3	TG	TG_C4	0.6017	0.2444	0.5453	FS
```

Read the shape, not just the number:

* **Parent–offspring is a signature, not an estimate.** `IBD1Seg 1.0000`, `IBD2Seg 0.0000`,
  `PropIBD 0.5000` — exact, because a child shares exactly one chromosome IBD with a parent
  everywhere. Any deviation from those three values on a pair you believe is PO is a data
  problem, not a biology problem.
* **`IBD2Seg` is what separates FS from PO.** Both have `PropIBD ≈ 0.5`; full sibs put a
  quarter of that into `IBD2Seg` and parent–child put none.
* **Everything below 1st degree has `IBD2Seg = 0`** and is distinguished only by how much
  IBD1 is left.

### What segments see that kinship does not

`multifam` family FAM2 contains a sib pair that realised unusually little sharing:

```
$ open-king -b /tmp/kingdocs/multifam.bed --kinship --prefix mfk
FAM2	B_C1	B_C2	15000	0.250	0.2500	0.1791	0.0295	0.1708	0.5

$ open-king -b /tmp/kingdocs/multifam.bed --related --prefix mfr
FAM2	B_C1	B_C2	15000	0.250	0.2500	0.1791	0.0295	0.3420	0.2642	0.1708	0.4582	0.1273	0.3564	FS	0
```

Kinship is 0.1708 in both runs — just under the 0.177 boundary, so kinship alone calls this
pair **2nd degree** and flags a pedigree mismatch. The segment scan sees `IBD1Seg 0.4582` and
`IBD2Seg 0.1273`: two-chromosome sharing over an eighth of the genome is not something a
2nd-degree relationship produces at all. `InfType` is `FS`, and the `Error` flag clears.
The dataset summaries say the same thing:

```
--kinship
Relationship summary (total relatives: 36 by pedigree, 36 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	24	12	0	0	4
  Inference	0	24	11	1	0	4

--related
Relationship summary (total relatives: 36 by pedigree, 36 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	24	12	0	0	4
  Inference	0	24	12	0	0	4
```

The general point: kinship is a genome-wide average and a borderline average is ambiguous;
the *distribution* of sharing along the chromosome is not. Prefer the segment columns when
they are available.

### When the segment columns exist at all

The IBD-segment engine needs physical marker density, and it reports whether it has it. Look
for this line on stdout:

```
Total length of NN chromosomal segments usable for IBD segment analysis is XXX Mb.
```

A "usable" segment is a maximal run of SNPs whose consecutive base-pair gaps are all
**≤ 156,250 bp**, and the segment columns are produced only when the usable total reaches
**100 Mb**. The `.bim` centimorgan column is ignored entirely — only base-pair positions
matter. Two of the corpus datasets (`trio`, `pair` — 5,000 SNPs spread over 22 chromosomes)
fall below the bar, as does any genome-wide fileset thinned much past ~20,000 markers:

```
$ python3 reshape.py /tmp/kingdocs/bigish thin4 --every 4       # 200 samples x 12500 SNPs
$ open-king -b thin4.bed --ibdseg --prefix oks
No informative IBD segments.
```

**If the `usable for IBD segment analysis` line is absent, IBD-segment analysis did not run,
and any `IBD1Seg` / `IBD2Seg` / `PropIBD` / `InfType` values in your output are not
measurements.** On that input KING 2.3.2 falls back to a kinship-only inference and says so;
open-king does not, and reports `UN` for every pair instead — a known divergence, measured in
[PARITY.md §5.12](PARITY.md#512-three-divergences-found-while-writing-the-user-documentation).
So the check below is not optional here. Fall back to reading kinship and IBS0 as in
[§3](#3-the-relationship-cutoffs) and [§4](#4-parentoffspring-vs-full-siblings), and treat
the relationship calls with corresponding caution.

This is one more reason not to LD-prune or MAF-filter before running `open-king`: thinning
markers destroys the density the segment caller depends on, and it does so silently.

---

## 6. The Error column

`.kin` carries three columns that are **not estimates**. `Phi` and `Z0` are read straight off
the declared pedigree — `Phi` is the pedigree kinship coefficient, `Z0` the pedigree
`Pr[IBD = 0]` — and `Error` grades the disagreement between what the pedigree says and what
the genotypes say.

| Pedigree relationship | `Z0` | `Phi` |
| --- | --- | --- |
| Parent–offspring | `0.000` | `0.2500` |
| Full siblings | `0.250` | `0.2500` |
| Unrelated within a family | `1.000` | `0.0000` |

`Error` is **graded and not an integer**. Its value set is exactly `{0, 0.5, 1}`:

| Value | Meaning |
| --- | --- |
| `0` | the data agree with the pedigree |
| `0.5` | off by exactly one degree |
| `1` | off by more than one degree |

Parsing this column as an integer silently turns every half-step disagreement into `0`.

### The two commands grade it differently

This is real and it is not a rounding artifact.

* **`--kinship`** compares the kinship *estimate* against `Phi` multiplicatively: within a
  factor of `sqrt(2)` is `0`, within a factor of 2 is `0.5`, beyond that `1`.
* **`--related`** compares the pedigree's relationship *label* against `InfType`. Exact match
  is `0`; `0.5` when the two degrees differ by one **and both are 2nd degree or more
  distant**; `1` otherwise. A pedigree `PO` inferred `FS` therefore scores `1`, even though
  both are 1st degree.

Across the 573 within-family pairs of `bigish`, four rows carry a non-zero `Error` under one
command or the other, and three of the four disagree between the two:

```
FID   ID1     ID2       Phi     Z0    Kinship  PropIBD  InfType  Err(--kinship)  Err(--related)
BF03  B03_C1  B03_C3   0.0625  0.750 0.0442   0.0822   4th      0.5             0.5
BF15  B15_C2  B15_G_F  0.1250  0.500 0.0924   0.1756   3rd      0               0.5
BF21  B21_C2  B21_G_F  0.1250  0.500 0.0886   0.1747   3rd      0               0.5
BF27  B27_C3  B27_G_M  0.1250  0.500 0.0879   0.1778   2nd      0.5             0
```

(A join of `bgk.kin` from `--kinship` with `bgr.kin` from `--related`, keeping the rows where
either `Error` is non-zero. Those four are the only such rows in the file.)

`BF15`/`BF21`: the kinship estimate is still inside the 2nd-degree band and within
`sqrt(2)` of `Phi`, so `--kinship` is content; the segment scan calls them 3rd degree, so
`--related` flags them. `BF27` is the mirror image — kinship falls 0.0005 below the
2nd/3rd boundary while the segments still say 2nd.

### How to use it

`Error` is a **pedigree-versus-data flag**, not a data-quality metric. A non-zero value has
three possible causes, and telling them apart is your job:

1. **A pedigree error or a sample swap.** The usual reason, and the reason the column exists.
   Look for a *consistent* pattern: one sample mismatching all of its declared relatives at
   once is a swap; one pair mismatching is a pedigree error.
2. **Genuine variance at 3rd/4th degree.** `BF03_C1`/`BF03_C3` above are true first cousins
   (`Phi 0.0625`) that realised kinship 0.0442 — exactly on the 3rd/4th boundary. Nothing is
   wrong with the data. See [§3](#3-the-relationship-cutoffs).
3. **A relationship the pedigree simply does not declare.** In `dups`, `MZ_1` and `MZ_2` are
   listed in one family with no stated relationship, so `Phi` is `0.0000`:

```
$ open-king -b /tmp/kingdocs/dups.bed --related --prefix dup
$ cat dup.kin
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType	Error
MZFAM	MZ_1	MZ_2	10000	1.000	0.0000	0.3468	0.0009	0.9951	0.0164	0.4962	0.0436	0.9223	0.9441	Dup/MZ	1
POFAM	PO_C	PO_P	10000	0.000	0.2500	0.1680	0.0000	0.3236	0.0000	0.2445	1.0000	0.0000	0.5000	PO	0
```

   `Error 1` here means "the pedigree says nothing and the genotypes say identical twins",
   which is a finding, not a fault.

`Error` appears only in `.kin`. Between-family pairs have no declared relationship to compare
against, so `.kin0` has no `Phi`, no `Z0` and no `Error`.

---

## 7. Pitfalls

Each of the following is demonstrated on corpus data whose truth is known by construction.

### 7.1 Too few SNPs

The estimator's sampling error scales like `1/sqrt(M)`, and the relationship bands are
narrow at the bottom. `unrelated` holds 30 mutually unrelated founders — every one of the
435 pairs has true kinship 0 — so every non-zero value below is pure noise. Thinning it and
re-running:

```
$ for k in 1 4 20 40 100; do
    python3 reshape.py /tmp/kingdocs/unrelated u$k --every $k
    open-king -b u$k.bed --kinship --prefix ku$k
  done
```

```
SNPs    file   pairs     mean       sd       min      max  | pairs called 4th/3rd/2nd degree
 20000 .kin      45   0.0001   0.0043  -0.0095   0.0086  |    0    0    0
 20000 .kin0    390  -0.0030   0.0064  -0.0224   0.0144  |    0    0    0
  5000 .kin      45   0.0025   0.0106  -0.0163   0.0251  |    1    0    0
  5000 .kin0    390  -0.0080   0.0128  -0.0357   0.0490  |    7    1    0
  1000 .kin      45  -0.0025   0.0349  -0.0669   0.0686  |    8    4    0
  1000 .kin0    390  -0.0109   0.0272  -0.0850   0.0984  |   38    7    1
   500 .kin      45  -0.0139   0.0449  -0.1284   0.0776  |    8    2    0
   500 .kin0      0   (header only: every sample dropped by the M<512 screen)
   200 .kin      45   0.0056   0.0585  -0.1360   0.1346  |   10    7    4
   200 .kin0      0   (header only: every sample dropped by the M<512 screen)
```

* At 20,000 SNPs the standard deviation is 0.006 and not one unrelated pair is misclassified.
* At 5,000 it is 0.013 and 8 unrelated pairs cross into the 4th-degree band; one reaches 3rd.
* At 1,000 it is 0.027–0.035 and 46 of 435 unrelated pairs look like 4th-degree relatives,
  11 like 3rd, and one like a **2nd-degree relative** at kinship 0.0984.
* At 200 SNPs, four of 45 pairs land in the 2nd-degree band and the maximum is 0.1346.

The `.kin0` spread tracks `0.9/sqrt(M)` almost exactly here (0.0064 / 0.0128 / 0.0272 against
predictions 0.0064 / 0.0127 / 0.0285). The constant depends on your allele-frequency
spectrum, so refit it on your own data rather than importing this one — but the `1/sqrt(M)`
shape holds, and it is the thing to reason with: **quartering your marker count doubles your
noise.**

Practical floors. Reliable calls down to 3rd degree want tens of thousands of markers.
Detecting only 1st-degree relatives is far more forgiving. And there is a hard gate: below
513 called autosomal markers a sample is dropped from the **between-family** stage
altogether, with

```
The following 30 samples are excluded from the kinship analysis (M<512):
```

`.kin0` then contains a header and nothing else. Within-family pairs are unaffected. (The
*names* printed in that message are a known KING display bug faithfully reproduced here —
the count is right, the names are just the first `count` samples in `.fam` order. See
[`SPEC.md` §5.2](SPEC.md).)

### 7.2 Missingness — and why it is worse when it is asymmetric

**Every count is pairwise.** `N_SNP` is not a property of the dataset, it is a property of
the pair, and `Het_i` is recounted for each pair over that pair's own marker set.

The `missing` dataset is one nuclear family with per-sample missingness of 0 %, 1 %, 5 %,
20 %, 50 % and 0 %, over 10,000 SNPs:

```
$ open-king -b /tmp/kingdocs/missing.bed --ibs --prefix misi
FID	ID1	ID2	Z0	Phi	N_SNP	N_IBS0	N_IBS1	N_IBS2	NHetHet	NHomHom	N_Het1	N_Het2	IBS	Dist	HetConc	Het2|1	Het1|2	HomConc	Kinship	MaxIBD2	Pr_IBD2
MIS	M_C1	M_C2	0.250	0.2500	7417	172	2267	4978	1420	3730	2532	2575	1.6480	0.3984	0.3851	0.5608	0.5515	0.9539	0.2107	76753302.000	0.1727
MIS	M_C1	M_C3	0.250	0.2500	4546	118	1627	2801	742	2177	1541	1570	1.5902	0.4617	0.3132	0.4815	0.4726	0.9458	0.1626	9548806.000	0.0000
MIS	M_C1	M_C4	0.250	0.2500	9266	147	2699	6420	1835	4732	3196	3173	1.6770	0.3547	0.4047	0.5742	0.5783	0.9689	0.2420	31948285.000	0.2109
MIS	M_C1	M_F	0.000	0.2500	9268	0	3156	6112	1641	4471	3191	3247	1.6595	0.3405	0.3421	0.5143	0.5054	1.0000	0.2549	0.000	0.0000
MIS	M_C2	M_C3	0.250	0.2500	3816	106	1028	2682	797	1991	1316	1306	1.6751	0.3805	0.4367	0.6056	0.6103	0.9468	0.2231	54351654.000	0.1407
MIS	M_C2	M_C4	0.250	0.2500	7779	179	2623	4977	1374	3782	2721	2650	1.6168	0.4292	0.3438	0.5050	0.5185	0.9527	0.1892	12748785.000	0.0255
MIS	M_C3	M_C4	0.250	0.2500	4768	65	1319	3384	988	2461	1649	1646	1.6961	0.3312	0.4283	0.5992	0.6002	0.9736	0.2604	35152194.000	0.0703
MIS	M_C3	M_F	0.000	0.2500	4764	0	1662	3102	814	2288	1649	1641	1.6511	0.3489	0.3288	0.4936	0.4960	1.0000	0.2474	0.000	0.0000
MIS	M_C4	M_F	0.000	0.2500	9732	0	3415	6317	1670	4647	3342	3413	1.6491	0.3509	0.3284	0.4997	0.4893	1.0000	0.2472	0.000	0.0000
MIS	M_F	M_M	1.000	0.0000	9650	658	4026	4966	1332	4292	3384	3306	1.4464	0.6899	0.2486	0.3936	0.4029	0.8467	0.0024	0.000	0.0000
```

(Ten of the file's fifteen rows.)

* `N_SNP` runs from **3,816 to 9,732** inside a single family. The pair that lost the most
  markers — `M_C2` (20 % missing) × `M_C3` (50 % missing) — is the one whose estimate rests
  on the least data, and its effective marker count is *multiplicative* in the two rates,
  which is why asymmetric missingness bites harder than it looks.
* **Parent–offspring survives it.** `M_C3`/`M_F`, on 4,764 markers, gives 0.2474 with
  `IBS0 = 0`. The PO signature is structural and does not need many markers.
* **Full sibs do not.** The six sib pairs estimate 0.2107, 0.1626, 0.2420, 0.2231, 0.1892 and
  0.2604 against a true 0.25, and `M_C1`/`M_C3` at 0.1626 falls into the 2nd-degree band.
  The summary line reports `FS 5  2nd 1` where the pedigree says `FS 6`.

Read `N_SNP` on every row you are about to act on, and get the per-sample picture first:

```
$ open-king -b /tmp/kingdocs/missing.bed --bysample --prefix q
$ cat qbySample.txt
FID IID FA MO SEX N_SNP Missing Heterozygosity N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal
MIS M_F 0 0 1 9784 0.0216 0.3505 31542 0 0.0000 31219 0 0.0000 0
MIS M_M 0 0 2 9711 0.0289 0.3420 31305 0 0.0000 31219 0 0.0000 0
MIS M_C1 M_F M_M 1 9337 0.0663 0.3449 18465 0 0.0000 9167 0 0.0000 0
MIS M_C2 M_F M_M 2 7831 0.2169 0.3504 15497 0 0.0000 7698 0 0.0000 0
MIS M_C3 M_F M_M 1 4801 0.5199 0.3456 9497 0 0.0000 4720 0 0.0000 0
MIS M_C4 M_F M_M 2 9788 0.0212 0.3441 19388 0 0.0000 9634 0 0.0000 0
```

(`bySample.txt` and `bySNP.txt` are **space** separated; `.kin`, `.kin0`, `.ibs`, `.ibs0`
and `.con` are **tab** separated. The asymmetry is real, and it catches parsers.)

A caveat on the mechanism: differential missingness also biases `Het_i` and `Het_j`
themselves whenever the missing markers are not a random sample of the genome (they rarely
are — call rate correlates with MAF and with assay difficulty). The heterozygosity difference
that results then feeds straight into the `min()` versus average denominator of
[§2](#2-the-two-estimators). Heavy, uneven missingness is therefore a *structure-like*
problem, not only a sample-size problem.

### 7.3 Population structure

KING-robust is **robust to** population structure. It is not **immune to** it, and the
difference is measurable.

`admixed` contains two populations at F<sub>ST</sub> 0.10 (11 unrelated founders each), plus
six unrelated founders with ancestry fractions α = 0.10, 0.25, 0.40, 0.50, 0.60, 0.90. Every
one of these people is unrelated to every other; true kinship is 0 throughout.

```
$ open-king -b /tmp/kingdocs/admixed.bed --ibs --prefix adm
```

The unrelated founders of `adm.ibs0` (between-family estimator), grouped by which population
each member came from:

```
class    n    kinship: mean       sd        min
P1-P1    55         -0.0045   0.0048    -0.0165
P2-P2    55         +0.0017   0.0055    -0.0128
P1-P2   121         -0.1139   0.0064    -0.1321
```

Within either population the estimator is centred on 0, as advertised. Across the two it is
centred on **−0.114**. That is not noise: its standard deviation is 0.006, so the bias is
eighteen standard deviations wide.

The six admixed founders make the mechanism plain — their pairwise kinship is a function of
the *difference in ancestry* and nothing else:

```
d_ancestry  pair                       kinship
   0.10     ADM_4 (0.50) x ADM_5 (0.60)   +0.0024
   0.15     ADM_1 (0.10) x ADM_2 (0.25)   -0.0078
   0.25     ADM_2 (0.25) x ADM_4 (0.50)   -0.0192
   0.30     ADM_1 (0.10) x ADM_3 (0.40)   -0.0222
   0.40     ADM_1 (0.10) x ADM_4 (0.50)   -0.0290
   0.50     ADM_3 (0.40) x ADM_6 (0.90)   -0.0487
   0.65     ADM_2 (0.25) x ADM_6 (0.90)   -0.0586
   0.80     ADM_1 (0.10) x ADM_6 (0.90)   -0.0767
```

(Eight of the fifteen pairs, sorted by ancestry gap.) The estimate slides over 0.08 kinship
units — almost two full relationship degrees' worth — with no relatedness anywhere in the
picture.

**Why this is still "robust".** The bias is *downward*. Structure makes unrelated people look
less related, never more. So a positive KING-robust kinship remains trustworthy evidence of
relatedness in a structured sample, which is exactly the guarantee frequency-based estimators
cannot give. What you lose is the other direction.

**What breaks, concretely.**

* **Related pairs whose members differ in ancestry are pulled down.** In `admixed`, families
  X and Z are within-population and family Y is cross-population. Parent–offspring rows of
  `adm.ibs` (within-family estimator):

  ```
  same-population PO (FAMX, FAMZ)   n=8  mean 0.2485  range 0.2440-0.2531
  cross-population PO (FAMY)        n=4  mean 0.2374  range 0.2351-0.2389
  ```

  A 4 % downward shift for a parent–child pair where one member is 50 % admixed — under the
  between-family estimator, which is what you would actually use on undeclared samples, the
  same four pairs read 0.2207–0.2237, a 10 % shift ([§2](#2-the-two-estimators)). Scale that
  to a fully cross-population 2nd- or 3rd-degree pair and the relationship can drop a whole
  band, or below the reporting threshold entirely. The 2010 paper notes the same residual
  bias for pairs that are both related and cross-population, and reports it as small out to
  3rd degree.
* **A single global threshold is wrong in a structured sample.** Any cut you place on
  kinship implicitly means something different for a within-ancestry pair than for a
  cross-ancestry pair.
* **Negative kinship is a diagnostic.** Do not clamp it, and do not discard it. A block of
  strongly negative values between two groups of samples is telling you those groups have
  different ancestry, before you have run any PCA.

**Segment columns do not carry this bias.** Segments are called from local haplotype
sharing, not from a genome-wide allele-sharing count, so ancestry does not enter. The same
three families under `--related` (columns trimmed to `FID ID1 ID2 Phi Kinship IBD1Seg
IBD2Seg PropIBD InfType`):

```
$ open-king -b /tmp/kingdocs/admixed.bed --related --prefix admr
$ cut -f1,2,3,6,11,12,13,14,15 admr.kin
FID	ID1	ID2	Phi	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
FAMX	X_C1	X_C2	0.2500	0.2727	0.5298	0.2884	0.5532	FS
FAMX	X_C1	X_F	0.2500	0.2468	1.0000	0.0000	0.5000	PO
FAMX	X_C1	X_M	0.2500	0.2440	1.0000	0.0000	0.5000	PO
FAMX	X_C2	X_F	0.2500	0.2448	1.0000	0.0000	0.5000	PO
FAMX	X_C2	X_M	0.2500	0.2498	1.0000	0.0000	0.5000	PO
FAMX	X_F	X_M	0.0000	-0.0069	0.0000	0.0000	0.0000	UN
FAMY	Y_C1	Y_C2	0.2500	0.2356	0.4375	0.1895	0.4083	FS
FAMY	Y_C1	Y_F	0.2500	0.2351	1.0000	0.0000	0.5000	PO
FAMY	Y_C1	Y_M	0.2500	0.2374	1.0000	0.0000	0.5000	PO
FAMY	Y_C2	Y_F	0.2500	0.2389	1.0000	0.0000	0.5000	PO
FAMY	Y_C2	Y_M	0.2500	0.2384	1.0000	0.0000	0.5000	PO
FAMY	Y_F	Y_M	0.0000	-0.1124	0.0000	0.0000	0.0000	UN
FAMZ	Z_C1	Z_C2	0.2500	0.2131	0.3728	0.2467	0.4331	FS
FAMZ	Z_C1	Z_F	0.2500	0.2473	1.0000	0.0000	0.5000	PO
FAMZ	Z_C1	Z_M	0.2500	0.2531	1.0000	0.0000	0.5000	PO
FAMZ	Z_C2	Z_F	0.2500	0.2521	1.0000	0.0000	0.5000	PO
FAMZ	Z_C2	Z_M	0.2500	0.2501	1.0000	0.0000	0.5000	PO
FAMZ	Z_F	Z_M	0.0000	-0.0048	0.0000	0.0000	0.0000	UN
```

Kinship for the cross-population `FAMY` pairs is visibly depressed — 0.2351 against 0.2473
for the equivalent `FAMZ` pair, and −0.1124 against −0.0048 for the unrelated couple. The
segment columns are identical across all three families: `1.0000 / 0.0000 / 0.5000 / PO` for
every parent–child pair regardless of ancestry, and `0.0000 / 0.0000 / 0.0000 / UN` for every
unrelated couple. In a structured sample the segment columns are the more trustworthy of the
two lines of evidence — when marker density permits them at all
([§5](#5-ibd1seg-ibd2seg-and-propibd)).

### 7.4 Genotyping error

Error creates opposite homozygotes out of nothing. Because the PO signature *is*
`IBS0 = 0`, parent–offspring pairs are the most fragile thing in the output.

`dups` contains an exactly duplicated sample and an MZ-twin pair simulated at a 0.2 %
per-genotype error rate:

```
$ open-king -b /tmp/kingdocs/dups.bed --duplicate --prefix dupd
FID1	ID1	FID2	ID2	N	N_IBS0	N_IBS1	N_IBS2	Concord	HomConc	HetConc
DUPA	DUP_A	DUPB	DUP_A_COPY	10000	0	0	10000	1.00000	1.00000	1.00000
MZFAM	MZ_1	MZFAM	MZ_2	10000	9	17	9974	0.99740	0.99862	0.99512
```

0.2 % error is enough to produce 9 opposite-homozygote sites and drop kinship from 0.5000 to
0.4962. Both pairs still call `Dup/MZ` — that band is wide.

Now the same experiment on a parent–offspring pair, injecting error into the child only:

```
$ open-king -b /tmp/kingdocs/dups.bed --related --prefix R0            # the rate-0 baseline
$ for r in 0.001 0.002 0.005 0.01 0.02 0.05; do
    python3 reshape.py /tmp/kingdocs/dups g$r --error $r --samples PO_C --seed 7
    open-king -b g$r.bed --related --prefix R$r
  done
```

```
rate    N_IBS0  IBS0_prop   Kinship   IBD1Seg   IBD2Seg  InfType
0            0     0.0000    0.2445    1.0000    0.0000       PO
0.001        1     0.0001    0.2438    0.9883    0.0000       PO
0.002        4     0.0004    0.2435    0.9995    0.0000       PO
0.005       11     0.0011    0.2401    0.9136    0.0000       PO
0.01        38     0.0038    0.2312    0.7451    0.0000      2nd
0.02        79     0.0079    0.2202    0.5379    0.0000      2nd
0.05       177     0.0177    0.1855    0.2503    0.0000      3rd
```

Three things to take from this.

* **IBS0 climbs roughly linearly in the error rate.** At 2 % it reaches 0.0079, past the
  ~0.0055 mark an IBS0-based PO/FS rule would sit at
  ([§4](#4-parentoffspring-vs-full-siblings)); at 5 % its IBS0 of 0.0177 lands squarely
  inside the FS range measured on clean data in the same section (0.0124–0.0214). Wherever
  PO and FS are told apart by IBS0 rather than by IBD2, **genotyping error turns parents into
  siblings**, and the resulting call is confidently wrong rather than uncertain.
* **The segment caller breaks first.** At 1 % error the kinship estimate is still 0.2312 —
  comfortably 1st degree — but `IBD1Seg` has already fallen to 0.745, because scattered
  errors chop long IBD tracts into pieces too short to survive the 3 Mb reporting floor.
  `InfType` says `2nd`. The segment columns are more informative *and* more error-sensitive
  than kinship; they are not a free upgrade.
* **Kinship decays gently.** 0.2445 → 0.1855 across a 50-fold increase in error rate. The
  ratio form absorbs a lot. If your kinship values look sane but your `InfType` calls do not,
  suspect genotyping quality before suspecting the pedigree.

`--bysample`'s `Err_MIp` / `Err_MIt` columns (Mendelian-inconsistency rates in declared
parent–offspring pairs and trios) are the direct measurement to make here.

### 7.5 Ascertained or filtered SNP panels

The kinship estimator is a ratio of genotype counts, and it survives a shifted allele-
frequency spectrum surprisingly well. **IBS0 does not**, and neither does the segment caller.

Three 12,000-SNP panels drawn from the same 200-sample fileset — rare-only, common-only, and
unfiltered:

```
$ python3 reshape.py /tmp/kingdocs/bigish rare   --maf-max 0.15 --first 12000
$ python3 reshape.py /tmp/kingdocs/bigish common --maf-min 0.35 --first 12000
$ python3 reshape.py /tmp/kingdocs/bigish allmaf                --first 12000
$ for t in rare common allmaf; do open-king -b $t.bed --kinship --prefix k_$t; done
```

```
panel   PO pairs             FS pairs                        unrelated (.kin0, n=19327)
        n   kinship  IBS0    n   kinship  IBS0 median range        mean      sd
rare    226  0.2480 0.00000  111  0.2478  0.0037 0.0013-0.0070    -0.0076   0.0118
        FS pairs with IBS0 below the 0.0055 PO/FS cutoff: 106 of 111
common  226  0.2511 0.00000  111  0.2508  0.0293 0.0137-0.0464    -0.0020   0.0108
        FS pairs with IBS0 below the 0.0055 PO/FS cutoff: 0 of 111
allmaf  226  0.2503 0.00000  111  0.2513  0.0158 0.0020-0.0413    -0.0017   0.0105
        FS pairs with IBS0 below the 0.0055 PO/FS cutoff: 5 of 111
```

* **Kinship is essentially unaffected**: 0.2480 / 0.2511 / 0.2503 for parent–offspring and
  0.2478 / 0.2508 / 0.2513 for full sibs, all against a true 0.25, and the unrelated spread
  barely moves (sd 0.0118 vs 0.0108).
* **IBS0 collapses on the rare panel.** Full-sib IBS0 falls eightfold, from a median of
  0.0293 on the common panel to 0.0037 on the rare one — because an opposite-homozygote site
  needs *both* alleles to be reasonably frequent. **106 of 111 true sib pairs now sit below
  the ~0.0055 PO/FS cutoff and would be called parent–offspring.** This is the exact mirror
  image of the genotyping-error failure in §7.4, and it is caused by nothing more exotic than
  a MAF filter.
* The same shift shows up at the other end of the scale: over the 137 within-family
  *unrelated* pairs, mean IBS0 is 0.1160 on the common panel, 0.0676 unfiltered and 0.0167 on
  the rare panel. Nothing about these people changed; only the markers did.

Two corollaries.

**`N_SNP` counts markers, not information.** Nothing is filtered out: every `.bim` record on
chromosomes 1–22 and 25/`XY` enters the computation, in file order, monomorphic markers
included — no MAF threshold, no call-rate threshold, no LD pruning, no de-duplication. On
`monomorphic` — 5,000 SNPs of which 1,639 are fixed in the sample and 745 more have MAF
below 5 % — `N_SNP` reads a confident `5000` for every pair, while the six true sib pairs
estimate:

```
$ open-king -b /tmp/kingdocs/monomorphic.bed --related --prefix mono
$ head -1 mono.kin; awk -F'\t' '$3 ~ /^P_C/' mono.kin
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType	Error
MONO	P_C1	P_C2	5000	0.250	0.2500	0.1576	0.0008	0.5194	0.0284	0.3384	0.9800	0.0000	0.4900	PO	1
MONO	P_C1	P_C3	5000	0.250	0.2500	0.1330	0.0044	0.3951	0.1146	0.2645	0.9007	0.0000	0.4504	PO	1
MONO	P_C1	P_C4	5000	0.250	0.2500	0.1084	0.0058	0.3203	0.1374	0.2167	0.8973	0.0000	0.4487	2nd	0
MONO	P_C2	P_C3	5000	0.250	0.2500	0.1608	0.0046	0.5407	0.1565	0.3309	0.6661	0.2457	0.5788	FS	0
MONO	P_C2	P_C4	5000	0.250	0.2500	0.1296	0.0146	0.4238	0.3883	0.2306	0.3506	0.2653	0.4406	FS	0
MONO	P_C3	P_C4	5000	0.250	0.2500	0.1064	0.0204	0.3152	0.4513	0.1477	0.4812	0.0000	0.2406	2nd	1
```

0.3384 (near-duplicate territory) down to 0.1477 (2nd degree) for six pairs of the same
relationship, with two called `PO` and two called `2nd`. Meanwhile all eight parent–offspring
pairs in the same file (`P_C1..P_C4` × `P_F`, `P_M`, estimates 0.2274–0.2693) are called
correctly — the PO signature survives what the FS signature does not. If you want to know how
much information your marker set really carries, get it from `--bySNP`, whose `Freq_A` and
`CallRate` columns give the per-marker picture, not from `N_SNP`.

**Do not prune or filter before running `open-king`.** KING's own guidance is not to LD-prune or
MAF-filter inputs, and the reasons are on this page: filtering shifts the IBS0 scale that
PO/FS discrimination depends on, and thinning destroys the marker density the segment engine
requires ([§5](#5-ibd1seg-ibd2seg-and-propibd)). Give it your QC-passing SNP set as-is.

Ascertainment that is *differential across your samples* — an array designed on one
population, applied to several — is the harder version of this problem, because it makes
heterozygosity depend on ancestry by construction and therefore behaves like
[§7.3](#73-population-structure) even when the samples are genuinely one population. There
is no synthetic demonstration of that here; treat it as the composition of the two effects
above.

### 7.6 Using the within-family estimator across populations

This is the mistake that follows from the previous two, and it is easy to make by accident.

The within-family form divides by the **average** heterozygosity. Its correctness rests on
the two heterozygosities estimating the same quantity — which is only true if the two people
share ancestry. The between-family form divides by the **minimum** precisely because that
assumption is unavailable across families.

You choose between them without meaning to, every time you write a `.fam` file:

* **Putting unrelated, ancestrally heterogeneous samples under a shared FID** — a cohort
  label, a batch id, a placeholder like `FAM1` for everyone — routes every pair through the
  within-family estimator, into `.kin`. On the cross-ancestry pairs of
  [§2](#2-the-two-estimators) that is a systematic 0.014–0.016 upward shift relative to the
  conservative form, in the direction that manufactures relatedness.
* **Splitting a genuine family across FIDs** does the opposite, and also costs you the `Phi`,
  `Z0` and `Error` columns, since those are read off the declared pedigree.

The rule: **a shared FID is an assertion of shared ancestry, not a grouping convenience.**
If you do not have real pedigree families, give every sample its own FID (`awk '{$1=$2; print}'`)
and read `.kin0`. If you do, make sure the FIDs are the actual families.

And when you compare numbers across files: a `.kin` value and a `.kin0` value computed from
identical counts differ by up to 0.016 in the table in §2. That is small next to a
relationship band, and large next to the 3rd/4th-degree boundary.

---

## 8. What relatedness inference cannot tell you

The estimator sees allele sharing. Several genuinely distinct situations produce the same
sharing, and no amount of data resolves them from genotypes alone.

**Direction.** Parent and child are symmetric. In `threegen`, the grandfather–father pair and
the father–child pair print the same three numbers:

```
TG	TG_GF	TG	TG_P1	1.0000	0.0000	0.5000	PO
TG	TG_P1	TG	TG_C1	1.0000	0.0000	0.5000	PO
```

`PO` means "one of these two is the parent of the other". Which one requires information
that is not in the genotypes. `--build`, which reconstructs pedigrees from the inferred
relationships, says so out loud:

```
$ open-king -b /tmp/kingdocs/multifam.bed --build --prefix bl
Reconstructing pedigree...
Age information not provided.
```

**Which 2nd-degree relationship.** Half siblings, grandparent–grandchild and avuncular pairs
all have `φ = 0.125`, `Pr[IBD=1] = 0.5`, `Pr[IBD=2] = 0`. All three occur in `threegen`, and
all three get the same label:

```
TG	TG_GF	TG	TG_C1	0.5543	0.0000	0.2772	2nd     <- grandparent / grandchild
TG	TG_P1	TG	TG_P3	0.4571	0.0000	0.2286	2nd     <- half siblings
TG	TG_P2	TG	TG_C1	0.4339	0.0000	0.2170	2nd     <- aunt / nephew
```

(The three differ subtly in the *length distribution* of their IBD segments, which some other
methods exploit. KING does not report that distinction, and neither does this program.)

**Twins versus a duplicated sample.** In `dups`, `DUP_A`/`DUP_A_COPY` is one person entered
twice and `MZ_1`/`MZ_2` is a twin pair; both land in the same `.con` block in
[§7.4](#74-genotyping-error) and both are labelled `Dup/MZ`. MZ twins are genetically
identical, so no genotype-based statistic can separate the two cases. The concordance
difference in that block is the simulated genotyping error, not the twinship — a real
re-genotyped duplicate would show the same kind of discordance. Deciding between them is a
sample-tracking question, not a genetics one.

**Generation, and therefore pedigree shape.** Kinship is symmetric and generation-blind. A
set of pairwise coefficients constrains a pedigree but rarely determines one.

**Anything about a single individual.** These are pair statistics. Individual inbreeding,
ancestry proportion and runs of homozygosity are separate analyses, and none of them is
implemented here — `--roh`, `--pca` and `--mds` are accepted on the command line but produce
no output. Use another tool for those.

**Whether a *specific* pair is related, at 3rd degree and beyond.** [§7.1](#71-too-few-snps)
is the honest version of this: at 3rd/4th degree the sampling distributions of "distant
relative" and "unrelated" overlap. Aggregate statements ("this cohort contains ~40
3rd-degree pairs") survive; per-pair statements at that range need corroboration.

---

## 9. A reading checklist

Before acting on a relatedness table:

1. **Which file is the row from?** `.kin` and `.kin0` use different estimators
   ([§2](#2-the-two-estimators)). Never merge them blind.
2. **What is `N_SNP` on that row?** Not the dataset's SNP count — the row's
   ([§7.2](#72-missingness--and-why-it-is-worse-when-it-is-asymmetric)). And check that those
   markers are informative, not monomorphic
   ([§7.5](#75-ascertained-or-filtered-snp-panels)).
3. **Did IBD-segment analysis run?** Look for the `usable for IBD segment analysis` line
   ([§5](#5-ibd1seg-ibd2seg-and-propibd)). If it did, prefer `IBD2Seg`/`PropIBD` over kinship
   for 1st-degree calls; if it did not, the segment columns are not measurements.
4. **Are there negative kinship values?** They indicate ancestry divergence
   ([§7.3](#73-population-structure)), which also means any positive value in that same
   comparison is an underestimate.
5. **Is the pair 3rd degree or more distant?** Then the band is narrower than your noise
   unless you have tens of thousands of markers ([§7.1](#71-too-few-snps)).
6. **Does a PO call rest on `IBS0 = 0`?** Check genotyping error
   ([§7.4](#74-genotyping-error)) and the MAF spectrum
   ([§7.5](#75-ascertained-or-filtered-snp-panels)) before believing it — both can flip PO
   and FS, in opposite directions.
7. **Is `Error` non-zero?** Decide which of the three causes it is before touching the
   pedigree ([§6](#6-the-error-column)). And parse the column as a float.

---

## Appendix — reshape.py

Three examples on this page ([§7.1](#71-too-few-snps), [§7.4](#74-genotyping-error),
[§7.5](#75-ascertained-or-filtered-snp-panels)) need a reshaped PLINK fileset. This is the
exact script used to produce them — standard library only, no PLINK required. The FID
rewrites in [§2](#2-the-two-estimators) need only `awk` and `--fam`.

```python
#!/usr/bin/env python3
"""reshape.py IN OUT [--every K] [--maf-min X] [--maf-max X] [--first N]
                     [--error R --samples IID,... --seed S]
Rewrites a PLINK1 SNP-major fileset: thin SNPs, filter on sample MAF, and/or
inject genotyping error into named samples. Standard library only."""
import argparse, random, sys
MAGIC = bytes((0x6C, 0x1B, 0x01))
TO_DOS = {0b00: 2, 0b10: 1, 0b11: 0, 0b01: None}
TO_BED = {2: 0b00, 1: 0b10, 0: 0b11, None: 0b01}

a = argparse.ArgumentParser(); a.add_argument("inp"); a.add_argument("out")
a.add_argument("--every", type=int); a.add_argument("--first", type=int)
a.add_argument("--maf-min", type=float); a.add_argument("--maf-max", type=float)
a.add_argument("--error", type=float); a.add_argument("--samples")
a.add_argument("--seed", type=int, default=1)
o = a.parse_args()

fam = [l.split() for l in open(o.inp + ".fam") if l.strip()]
bim = [l.rstrip("\n").split() for l in open(o.inp + ".bim") if l.strip()]
n, bpr = len(fam), (len(fam) + 3) // 4
body = open(o.inp + ".bed", "rb").read()[3:]
geno = [[TO_DOS[(body[k * bpr + (i >> 2)] >> (2 * (i & 3))) & 3] for i in range(n)]
        for k in range(len(bim))]

def maf(r):
    c = [d for d in r if d is not None]
    f = sum(c) / (2.0 * len(c)) if c else 0.0
    return min(f, 1.0 - f)

keep = range(len(bim))
if o.maf_min is not None: keep = [k for k in keep if maf(geno[k]) >= o.maf_min]
if o.maf_max is not None: keep = [k for k in keep if maf(geno[k]) <= o.maf_max]
if o.every: keep = list(keep)[::o.every]
if o.first: keep = list(keep)[:o.first]
bim, geno = [bim[k] for k in keep], [geno[k] for k in keep]

if o.error:
    rng = random.Random(o.seed)
    want = set(o.samples.split(",")) if o.samples else {r[1] for r in fam}
    idx = [i for i, r in enumerate(fam) if r[1] in want]
    for row in geno:
        for i in idx:
            if row[i] is not None and rng.random() < o.error:
                row[i] = rng.choice([d for d in (0, 1, 2) if d != row[i]])

open(o.out + ".fam", "w").writelines(" ".join(r) + "\n" for r in fam)
open(o.out + ".bim", "w").writelines("\t".join(r) + "\n" for r in bim)
with open(o.out + ".bed", "wb") as fh:
    fh.write(MAGIC)
    for row in geno:
        blk = bytearray(bpr)
        for i, d in enumerate(row): blk[i >> 2] |= TO_BED[d] << (2 * (i & 3))
        fh.write(bytes(blk))
sys.stderr.write("%s: %d samples x %d SNPs\n" % (o.out, n, len(bim)))
```

Reshaped filesets keep the original `.bim` positions for the markers they retain, which is
why thinning changes IBD-segment usability: the base-pair gaps grow even though the
coordinates do not move.

---

## See also

* [`README.md`](README.md) — the documentation index
* [`CLI.md`](CLI.md) — every command-line option, and which analyses each one affects
* [`OUTPUTS.md`](OUTPUTS.md) — every output file: columns, formats, row order, existence rules
* [`COOKBOOK.md`](COOKBOOK.md) — task-oriented recipes, from finding duplicates to diffing
  against KING
* [`PARITY.md`](PARITY.md) — the authoritative statement of what is byte-identical to
  KING 2.3.2, measured per file and per row

---

## References

Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen W-M. Robust relationship
inference in genome-wide association studies. *Bioinformatics* 2010;26(22):2867–2873.
doi:10.1093/bioinformatics/btq559. PMC3025716.

Formula-by-formula provenance, including which equations were verified numerically against
the reference binary and which are unverified, is in
[`VERIFIED_FORMULAS.md`](VERIFIED_FORMULAS.md). No KING source code was read in the writing
of this program or this page.
