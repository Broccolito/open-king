# KING 2.3.2 — Website & Manual Recon (clean-room reference)

**Recon date:** 2026-08-13
**Source:** <https://www.kingrelatedness.com/> (public documentation only)
**Reference binary cross-check:** `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` — banner + embedded string constants only (facts about output format). **No C++ source was read or transcribed.**

---

## 0. LEGAL / PROVENANCE NOTE

Everything in this document comes from one of two permitted sources:

1. **Public HTML documentation** on kingrelatedness.com (fetched with `curl`, converted to text locally).
2. **Observed binary facts**: the no-argument banner printed by the 2.3.2 executable, and `strings(1)` output of filename/format constants and the *embedded R plotting code* (which KING writes verbatim to `*.R` files at runtime — i.e. it is program output, not program source).

The website links a `KINGcode.tar.gz` C++ tarball. **It was not downloaded and not read.** Its stated distribution terms, quoted verbatim from `Download.shtml`, are:

> "Feel free to use KING for your research, but please do not redistribute AND make profits."

That is **not** an OSI license and is **not** MIT-compatible. It reinforces that the reimplementation must be clean-room: derive behavior from the paper (Manichaikul et al. 2010) + the documented I/O formats below + observed binary behavior.

---

## 1. PAGES FETCHED (complete site map)

All pages under the root. Every internal `href` on the site was enumerated; the list below is exhaustive for HTML pages.

| Page | Bytes | "Last updated" footer |
|---|---|---|
| `index.shtml` (home) | 6051 | — |
| `manual.shtml` (Relationship Inference — **the main manual**) | 39479 | July 28, 2023 |
| `KINGvisualization.shtml` (Visualization of Families) | 13268 | May 21, 2019 |
| `QC.shtml` (Quality Control) | 12960 | February 21, 2018 |
| `kingpopulation.shtml` (Population Structure) | 9705 | Oct 11, 2019 |
| `ancestry/` (Ancestry Inference) | 13616 | May 28, 2021 |
| `genemapping.shtml` (Association Mapping) | 8647 | October 10, 2022 |
| `riskprediction.shtml` (Risk Prediction) | 10077 | August 24, 2018 |
| `Download.shtml` | 7538 | September 8, 2023 |
| `history.shtml` (Download History / changelog) | 34091 | — |
| `flagindex.shtml` (Complete Flag Index) | 9772 | — |
| `forum/` | (phpBB, not scraped) | — |

Non-HTML assets referenced: `ex.tar.gz` (1.35 MB example dataset), `KINGcode.tar.gz`, `Linux-king.tar.gz`, `Windows-king.zip`, `KGref.{bed,bim,fam}.xz`, `hapmapkin.{pdf,R}`, `hapmapkin0.{pdf,R}`, `king_segments_plot.R`, `ex_relplot.pdf`, `executables/*` (every historical release).

**Example dataset (used throughout the tutorial):** `ex.tar.gz` — 332 HapMap samples (165 CEU + 167 YRI), each genotyped at **18,290 SNPs**.

---

## 2. ⚠️ CRITICAL CAVEATS BEFORE USING ANY OF THIS FOR BYTE PARITY

These three points matter more than anything else in this document.

### 2.1 The website's whitespace is NOT the file's whitespace

I dumped the raw HTML: **`manual.shtml` contains ZERO tab characters (`grep -c $'\t'` → 0).** Every "column" in every `<PRE>` example on the site is space-padded to 8-column stops — the author expanded tabs when writing the HTML.

Real KING output files are **tab-delimited** (`\t`). The site's alignment is a rendering artifact. **Do not derive field widths from the website.** Column *names* and *order* from the site are authoritative; column *separators* and *padding* must come from running the real binary.

### 2.2 The website documents ~2.2.7; the shipped binary is 2.3.2

The manual page still says "The latest version is KING 2.3.1" in prose while `Download.shtml` announces 2.3.2, and the option lists on the site are stale. Diffing the site's flag index against the **actual 2.3.2 banner**:

| Change | Detail |
|---|---|
| **`--homog` REMOVED** | Retired in 2.3.0: *"--homog is retired since it does not do a good job on either relationship inference or GRM"*. Still documented on `manual.shtml`. **Not in the 2.3.2 banner.** |
| **`--mtscore` REMOVED** | Superseded by `--lmm`. Still in flag index. Not in 2.3.2 banner. |
| **`--makeGRM` ADDED** (2.3.0) | Not documented anywhere on the website. |
| **`--lmm` ADDED** (2.3.0) | "work in progress" per 2.3.1 notes. |
| **`--gdt` ADDED** (2.3.1) | Family-based association (Chen et al. 2009). |
| **`--minConc` ADDED** (2.3.1) | Heterozygote concordance cutoff for `--duplicate`; **default 0.80**. |
| **`--noscreen` ADDED** (2.3.2) | For `--related` and `--duplicate`; skips the SNP-subset screening pass. |
| **`--phefile`, `--covfile`, `--prunedsnp` ADDED** | Optional input files, undocumented on site. |
| **`--lessmem` RETIRED** (2.1.6) | Flag index still lists it as retired. |
| `--maxP` moved | From "Association Parameter" to "Association Model" group. |

### 2.3 `--noscreen` prints an uninitialized default

The 2.3.2 banner prints `--noscreen [-1717986816]`. That is `0xCCCCCCC0`-ish garbage — an uninitialized integer being formatted as a default. It is a **bug in KING's banner**, not a meaningful default. Do not replicate the number; replicate the flag.

---

## 3. THE 2.3.2 OPTION BANNER (verbatim, observed)

Printed by the binary with no arguments. **This is the authoritative option list for 2.3.2.** Bracketed values are defaults.

```
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen

The following parameters are in effect:
                   Binary File :                 (-bname)

Additional Options
         Close Relative Inference : --related, --duplicate
   Pairwise Relatedness Inference : --kinship, --ibdseg, --ibs, --makeGRM
              Inference Parameter : --degree, --noscreen [-1717986816],
                                    --seglength, --minConc [0.80]
         Relationship Application : --unrelated, --cluster, --build
                        QC Report : --bysample, --bySNP, --roh, --autoQC
                     QC Parameter : --callrateN, --callrateM
             Population Structure : --pca, --mds
              Structure Parameter : --projection, --pcs
          Quantitative Trait GWAS : --lmm
                Binary Trait GWAS : --tdt, --gdt
                Association Model : --trait [], --covariate [], --maxP
     Association Method Parameter : --invnorm
               Genetic Risk Score : --risk, --model [], --prevalence, --noflip
              Computing Parameter : --cpus
                   Optional Input : --fam [], --bim [], --phefile [],
                                    --covfile [], --prunedsnp [],
                                    --sexchr [23]
                           Output : --rplot, --pngplot, --plink
                 Output Parameter : --prefix [king], --rpath []


FATAL ERROR - 
Genotype files are required. e.g.,
  king -b ex.bed --related

Please check the reference paper Manichaikul et al. 2010 Bioinformatics,
					Chen et al. 2024,
          or the KING website at kingrelatedness.com
```

Note: the fatal-error footer contains literal tab characters before "Chen et al. 2024".

### 3.1 Historical banner shape (KING 2.2.7, verbatim from ancestry page)

Useful as a diff baseline; shows the ON-flag notation `--pca [ON]`, `--projection [1]`, `--rplot [ON]`, `--prefix [ex]`:

```
KING 2.2.7 - (c) 2010-2021 Wei-Min Chen

The following parameters are in effect:
                   Binary File : ../data/KGref,../data/ex (-bname)

Additional Options
         Close Relative Inference : --related, --duplicate
   Pairwise Relatedness Inference : --kinship, --ibdseg, --ibs, --homog
              Inference Parameter : --degree, --seglength
         Relationship Application : --unrelated, --cluster, --build
                        QC Report : --bysample, --bySNP, --roh, --autoQC
                     QC Parameter : --callrateN, --callrateM
             Population Structure : --pca [ON], --mds
              Structure Parameter : --projection [1], --pcs
              Disease Association : --tdt
   Quantitative Trait Association : --mtscore
                Association Model : --trait [], --covariate []
            Association Parameter : --invnorm, --maxP
               Genetic Risk Score : --risk, --model [], --prevalence, --noflip
              Computing Parameter : --cpus
                   Optional Input : --fam [], --bim [], --sexchr [23]
                           Output : --rplot [ON], --pngplot, --plink
                 Output Parameter : --prefix [ex], --rpath []
```

---

## 4. COMPLETE FLAG INDEX (verbatim from `flagindex.shtml`)

Two-column table, OPTION → FUNCTION. Reproduced exactly (this page is stale — see §2.2).

| OPTION | FUNCTION |
|---|---|
| `--autoQC` | Quality control (QC) including call rate and gender checking |
| `--bim` | Specify .bim file as alternative input |
| `--build` | Reconstruct pedigrees using information from SNP data |
| `--bysample` | Sample-level QC |
| `--bySNP` | SNP-level QC |
| `--callrateM` | Specify SNP-level call rate for QC |
| `--callrateN` | Specify sample-level call rate for QC |
| `--cluster` | Cluster individuals into families according to inferred relatedness |
| `--covariate` | Specify covariate names to be adjusted in association analysis |
| `--cpus` | Specify number of CPU cores for parallel computing |
| `--degree` | Specify degree of relatedness for all relatives to be inferred |
| `--duplicate` | Identify duplicate pairs (including MZ twins) using autosome SNP data |
| `--fam` | Specify .fam file as alternative input |
| `--homog` | Infer relatedness assuming a homogeneous population |
| `--ibdseg` | Infer IBD segments shared between any two samples using SNP data |
| `--ibs` | Provide IBS statistics between any two samples using autosome SNP data |
| `--invnorm` | Inverse normal transformation for quantitative traits prior to association scan |
| `--kinship` | Estimate kinship coefficients between any two samples using SNP data |
| `--lessmem` | Request less memory usage, now retired in KING 2.1.6 and later |
| `--maxP` | Specify maximum p-values for being included in the output files |
| `--mds` | SNP-based multi-demensional scaling (MDS) for ancestry inference |
| `--model` | Specify a model template file for risk prediction |
| `--mtscore` | GWAS scan with a many traits version of score test |
| `--noflip` | Specify no-flip flag for risk prediction |
| `--pca` | Compute principal components of ancestry using autosome SNP data |
| `--pcs` | Specify the number of PCs used for PCA/MDS, e.g., 10 as default |
| `--prefix` | Specify prefix of files for inference results, e.g., "king" as default |
| `--prevalence` | Specify disease prevalence in the general population for risk prediction |
| `--projection` | Project samples onto the principal component space of reference samples |
| `--projection N` | Relatedness inference between two subsets where the first subset includes the first N samples |
| `--related` | Fast and integrated relationship inference to identify close relatives |
| `--risk` | Predict disease risk using genetic risk scores (GRS) |
| `--roh` | Scan for runs of homozygosity |
| `--rpath` | Full path of the R program. e.g., --rpath R |
| `--rplot` | Plot inference results using R code |
| `--sexchr` | Specify pair number of the sex chromosome for non-human species |
| `--tdt` | GWAS scan with transmission-disequilibrium test |
| `--trait` | Specify trait names for association scan |
| `--unrelated` | Extract a subset of unrelated individuals |

**Documented defaults, consolidated:**

| Option | Default |
|---|---|
| `--prefix` | `king` |
| `--pcs` | `10` (flag index says "e.g., 10 as default"; population page says "The default pcs is 10"; MDS prose says "20 by default" — **contradictory, see §12**) |
| `--cpus` | half the total number of *logical* cores |
| `--sexchr` | `23` |
| `--callrateN` | 0.95 (95%) |
| `--callrateM` | 0.95 (95%) |
| `--minConc` | 0.80 (2.3.1+) |
| `--degree` | unset; `--related` alone ⇒ 1st degree only |

---

## 5. INPUT FILE CONVENTIONS

Verbatim from `manual.shtml`, section **GENERAL INPUT FILES**:

> The input files for KING need to be in PLINK binary format, which include a binary genotype file, a family file, and a map file, e.g., `ex.bed`, `ex.fam`, and `ex.bim`.

Memory model, verbatim:

> The amount of computer memory required by KING analysis is modest, at ~N ✕ M / 4 (where N is the number of samples and M is the number of SNPs) plus a small percentage of overhead cost. E.g., for a dataset consisting of 100,000 samples each genotyped at 1,000,000 SNPs, the required memory size is ~25GB.

### 5.1 Documented input example command lines (verbatim)

```
  prompt> king -b ex.bed --related
  prompt> king -b ex.bed --fam ex.fam --bim ex.bim --related
```

> In the first example, although only `ex.bed` is specified, the other two input files are pre-assumed to be `ex.fam` and `ex.bim`.

### 5.2 Multiple datasets

```
  prompt> king -b ex.bed,mystudy.bed --duplicate
```

> KING reads in two sets of data (`ex.bed`, `ex.fam`, `ex.bim`, `mystudy.bed`, `mystudy.fam`, `mystudy.bim`) and then identifies all duplicate pairs, within and across datasets.

Rules, verbatim:

- "users do not need to worry about allele strands, which are well taken care in KING by either autoflip (at unambiguous SNPs) or removal (at ambiguous SNPs)"
- "Note KING is only looking at the SNP names when merging multiple datasets and other information such as chr and positions is not utilized."
- "all IDs (combinations of FID and IID) need to be unique within and across datasets."
- "One exception is majority of samples can be overlapped across multiple datasets, in which case ID names are modified with either `"REF_"` or `"QRY_"`"

```
 prompt> king -b ex.bed,ex.bed --duplicate
```

The ancestry page restates the merge rules verbatim as a numbered list:

```
1. SNPs with unambiguous allele labels can be auto-flipped before merging
2. SNPs with ambiguous allele labels (i.e., A/T, or C/G) are excluded
3. SNPs with inconsistent allele labels (e.g., >3 alleles) are excluded
```

`--fam` and `--bim` also accept comma-separated lists (from `KINGvisualization.shtml`):

```
  prompt> king -b ex.bed,ex.bed --fam ex.fam,ex2.fam --duplicate --rplot
```

### 5.3 VCF conversion (documented)

```
 prompt> plink1.9 --vcf example.vcf.gz --make-bed --out ex
```

### 5.4 Gene-mapping extra inputs

> Besides the standard PLINK binary format (`ex.fam`, `ex.bim`, `ex.bed`), two other files can be specified, including `ex.phe` for phenotypes and `ex.cov` for covariates. KING searches for all 5 files automatically even though only one file (`ex.bed`) needs to be specified in the command line.

(2.3.2 also exposes explicit `--phefile` / `--covfile`.)

### 5.5 Pre-processing guidance (verbatim, matters for parity of SNP counts)

> Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference, unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in a WGS study, in which case rare variants can be filtered out. **LD pruning is not recommended in KING.**

But for PCA specifically: "Please run LD-pruning prior to PCA analysis."

---

## 6. RELATIONSHIP INFERENCE — `--kinship`

**What it does:** estimates pairwise kinship coefficients using the **KING-robust** algorithm from the 2010 paper. Robust to population structure. Accurate "up to 2nd-degree".

**Pedigree handling, verbatim:**
> If pedigrees are documented in the .fam file, kinship coefficients can be estimated within families. Note if each FID is unique and no pedigrees are provided, then the within-family inference will be skipped.
> The output files are separated for relationships that are within or between families. Note an unrelated individual is treated as a family of size one. If the datasets only consist of unrelated individuals as reported, then all results are saved in the between-family output.

### 6.1 Output files

| File | Content |
|---|---|
| `{prefix}.kin` | within-family pairs |
| `{prefix}.kin0` | between-family pairs |
| `{prefix}X.kin`, `{prefix}X.kin0` | X-chromosome counterparts (observed in binary strings) |

Default prefix `king` ⇒ `king.kin`, `king.kin0`.

### 6.2 `{prefix}.kin` — columns, in order

```
FID  ID1  ID2  N_SNP  Z0  Phi  HetHet  IBS0  Kinship  Error
```

Verbatim column definitions:

```
FID: Family ID for the pair
ID1: Individual ID for the first individual of the pair
ID2: Individual ID for the second individual of the pair
N_SNP: The number of SNPs that do not have missing genotypes in either of the individual
Z0: Pr(IBD=0) as specified by the provided pedigree data
Phi: Kinship coefficient as specified by the provided pedigree data
HetHet: Proportion of SNPs with double heterozygotes (e.g., AG and AG)
IBS0: Porportion of SNPs with zero IBS (identical-by-state) (e.g., AA and GG)
Kinship: Estimated kinship coefficient from the SNP data
Error: Flag indicating differences between the estimated and specified kinship coefficients (1 for error, 0.5 for warning)
```

(Note: "Porportion" is a typo **in the official docs**, not mine.)

**Documented example block (verbatim; site whitespace, not real tabs):**

```
FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    Kinship Error
28      1       2       2359853 0.000   0.2500  0.162   0.0008  0.2459  0
28      1       3       2351257 0.000   0.2500  0.161   0.0008  0.2466  0
28      2       3       2368538 1.000   0.0000  0.120   0.0634  -0.0108 0
117     1       2       2354279 0.000   0.2500  0.163   0.0006  0.2477  0
117     1       3       2358957 0.000   0.2500  0.164   0.0006  0.2490  0
117     2       3       2348875 1.000   0.0000  0.122   0.0616  -0.0017 0
1344    1       12      2372286 0.000   0.2500  0.149   0.0003  0.2480  0
1344    1       13      2370435 0.000   0.2500  0.148   0.0003  0.2465  0
1344    12      13      2374888 1.000   0.0000  0.117   0.0582  0.0003  0
```

**Inferred number formatting from this block (to be confirmed against the binary):**

| Column | Apparent format | Evidence |
|---|---|---|
| `N_SNP` | `%d` | `2359853` |
| `Z0` | `%.3f` | `0.000`, `1.000` |
| `Phi` | `%.4f` | `0.2500`, `0.0000` |
| `HetHet` | **`%.3f`** | `0.162`, `0.120`, `0.117` — **3 dp here, but 4 dp in `.kin` from `--related` (see §8.2). This is a real inconsistency between the two writers.** |
| `IBS0` | `%.4f` | `0.0008`, `0.0634` |
| `Kinship` | `%.4f` | `0.2459`, `-0.0108`, `0.0003` |
| `Error` | `%d` or `%.2g` | `0`; docs say 0.5 is possible for warning |

### 6.3 `{prefix}.kin0` — columns, in order

```
FID1  ID1  FID2  ID2  N_SNP  HetHet  IBS0  Kinship
```

**Documented example block (verbatim):**

```
FID1    ID1     FID2    ID2     N_SNP   HetHet  IBS0    Kinship
28      3       117     1       2360618 0.143   0.0267  0.1356
28      3       117     2       2352628 0.161   0.0009  0.2441
28      3       117     3       2354540 0.120   0.0624  -0.0119
28      3       1344    1       2361807 0.093   0.1095  -0.2295
28      3       1344    12      2367180 0.094   0.1091  -0.2225
28      3       1344    13      2364816 0.093   0.1082  -0.2224
117     1       1344    1       2362787 0.094   0.1093  -0.2312
117     1       1344    12      2368467 0.095   0.1088  -0.2230
117     1       1344    13      2365036 0.094   0.1084  -0.2253
117     2       1344    1       2354855 0.094   0.1084  -0.2281
117     2       1344    12      2361351 0.095   0.1078  -0.2206
117     2       1344    13      2357936 0.095   0.1067  -0.2190
117     3       1344    1       2357771 0.094   0.1102  -0.2348
117     3       1344    12      2364365 0.095   0.1086  -0.2232
117     3       1344    13      2361061 0.094   0.1096  -0.2301
```

Here `HetHet` is again `%.3f`, `IBS0` and `Kinship` are `%.4f`.

> This analysis shows the "unrelated" families 28 and 117 are actually connected through an unreported parent-offspring pair (FID 28 IID 3, and FID 117 IID 2).

### 6.4 Semantics of negative kinship (verbatim — important, do not clamp to 0)

> A negative kinship coefficient estimation indicates an unrelated relationship. The reason that a negative kinship coefficient is not set to zero is a very negative value may indicate the population structure between the two individuals.

### 6.5 Degree cutoffs (verbatim)

> an estimated kinship coefficient range >0.354, [0.177, 0.354], [0.0884, 0.177] and [0.0442, 0.0884] corresponds to duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree relationships respectively.

```
  prompt> king -b ex.bed --kinship --degree 2
```
> In this example, only pairs with kinship coefficient > 0.0884 are saved in the king.kin0 output file.

### 6.6 `--kinship --projection N` (KING 2.2.2+)

> `--kinship --projection N` estimates the kinship coefficients between any two samples each from a different subset where the first subset includes the first N samples.

Verbatim example (note the `--proj` abbreviation is used):

```
  prompt> king -b subset1,subset2 --kinship --proj 100000 --prefix subset12
  prompt> king -b subset1,subset3 --kinship --proj 100000 --prefix subset13
  prompt> king -b subset1,subset4 --kinship --proj 100000 --prefix subset14
  prompt> king -b subset2,subset3 --kinship --proj 100000 --prefix subset23
  prompt> king -b subset2,subset4 --kinship --proj 100000 --prefix subset24
  prompt> king -b subset3,subset4 --kinship --proj 100000 --prefix subset34
  prompt> king -b subset1.bed --kinship --prefix subset1
  prompt> king -b subset2.bed --kinship --prefix subset2
  prompt> king -b subset3.bed --kinship --prefix subset3
  prompt> king -b subset4.bed --kinship --prefix subset4
```

> the kinship estimates from the `--kinship --projection N` inference should be idential (with no numerical differences) to the standard `--kinship` inference without splitting.

**Parity note:** this is a hard guarantee — the estimator must be exactly pairwise (no dataset-wide allele frequencies), so splitting cannot change a single digit. Good invariant to test.

---

## 7. IBD SEGMENT INFERENCE — `--ibdseg`

Introduced in **KING 2.1**. Accurate "up to 3rd- or 4th-degree (depending on array or WGS)".

**Documented example command lines:**

```
  prompt> king -b ex.bed --ibdseg
  prompt> king -b ex.bed --ibdseg --degree 3 --rplot --prefix ex
```

> The second command specifies only pairs with IBD proportion > 0.0884 will be saved in the output.

`--projection N` is also available for `--ibdseg`.

Since 2.2.3: "IBD segment analysis and Run of Homozygosity anaysis now apply to chromosome X as well."

### 7.1 Output files

| File | Content |
|---|---|
| `{prefix}.seg` | per-pair IBD summary |
| `{prefix}.segments.gz` | gzipped per-segment detail |
| `{prefix}allsegs.txt` | genome-wide usable-segment map (observed in binary strings + embedded R) |
| `{prefix}X.seg` | X-chromosome summary (observed) |

Observed runtime messages (format strings):
- `Summary statistics of IBD segments for individual pairs saved in file %s`
- `Additional summary statistics of X-Chr IBD segments saved in file %s`
- `Total length of %d chromosomal segments usable for IBD segment analysis is %.1lf Mb.`
- `In addition to autosomes, %d segments of length %.1lf Mb on X-chr can be further used.`
- `Information of these chromosomal segments can be found in file %s`
- `Short IBD segments (<3Mb) are not reported/utilized.`
- `Sample pairs without any long IBD segments (>10Mb) are excluded.`

**Those last two are load-bearing algorithm parameters: min reported segment length 3 Mb; pair inclusion requires ≥1 segment >10 Mb.**

### 7.2 `{prefix}.seg` — columns, in order

```
FID1  ID1  FID2  ID2  IBD1Seg  IBD2Seg  PropIBD  InfType
```

Verbatim column definitions:

```
FID1: Family ID for the first individual of the pair
ID1: Individual ID for the first individual of the pair
FID2: Family ID for the second individual of the pair
ID2: Individual ID for the second individual of the pair
IBD1Seg: Total length of IBD1 segments divided by total length of all segments, estimate of π1=Pr(IBD=1)
IBD2Seg: Total length of IBD2 segments divided by total length of all segments, estimate of π2=Pr(IBD=2)
PropIBD: Proportion of genomes shared identical-by-descent, estimated by IBD2Seg + IBD1Seg/2, estimate of π=π2+π1/2
InfType: Inferred relationship type, such as Dup/MZTwin, PO, FS, 2nd, 3rd, 4th, UN
```

**Documented example block (verbatim):**

```
FID1    ID1     FID2    ID2     IBD1Seg IBD2Seg PropIBD InfType
1330    NA12335 1330    NA12340 0.9976  0.0000  0.4988  PO
1330    NA12335 1330    NA12341 1.0000  0.0000  0.5000  PO
1330    NA12336 1330    NA12342 0.9969  0.0000  0.4985  PO
1330    NA12336 1330    NA12343 0.9987  0.0000  0.4994  PO
1334    NA10846 1334    NA12144 1.0000  0.0000  0.5000  PO
1334    NA10846 1334    NA12145 1.0000  0.0000  0.5000  PO
1334    NA10847 1334    NA12146 1.0000  0.0000  0.5000  PO
1334    NA10847 1334    NA12239 0.9999  0.0000  0.5000  PO
1340    NA06994 1340    NA07029 1.0000  0.0000  0.5000  PO
```

Formats: `IBD1Seg`, `IBD2Seg`, `PropIBD` are all **`%.4f`**. `InfType` is a bare string.

**Exact `InfType` vocabulary (from the docs):** `Dup/MZTwin`, `PO`, `FS`, `2nd`, `3rd`, `4th`, `UN`.

### 7.3 `{prefix}.segments.gz` — columns, in order

```
FID1  ID1  FID2  ID2  IBDType  Chr  StartMB  StopMB  StartSNP  StopSNP  N_SNP  Length
```

Verbatim column definitions:

```
FID1: Family ID for the first individual of the pair
ID1: Individual ID for the first individual of the pair
FID2: Family ID for the second individual of the pair
ID2: Individual ID for the second individual of the pair
IBDType: Type of IBD segments: IBD1 or IBD2
Chr: Chromosome number.
StartMB: Start position of the IBD segment (in Mb)
StopMB: Stop position of the IBD segment (in Mb)
StartSNP: Start SNP of the IBD segment
StopSNP: Stop SNP of the IBD segment
N_SNP: The number of SNPs in the IBD segment
Length: Total Length of the IBD segment (in Mb)
```

**Documented example block (verbatim):**

```
FID1    ID1     FID2    ID2     IBDType Chr     StartMB StopMB  StartSNP        StopSNP         N_SNP   Length
1330    NA12335 1330    NA12340 IBD1    1       51.799  95.862  rs7534689       rs1858111       294     44.1
1330    NA12335 1330    NA12340 IBD1    1       148.175 247.083 rs1868992       rs12058711      692     98.9
1330    NA12335 1330    NA12340 IBD1    2       0.143   88.714  rs408209        rs7581608       619     88.6
1330    NA12335 1330    NA12340 IBD1    2       165.994 242.590 rs1835889       rs10186231      484     76.6
1330    NA12335 1330    NA12340 IBD1    3       0.080   90.221  rs990284        rs9877833       643     90.1
1330    NA12335 1330    NA12340 IBD1    3       113.243 165.061 rs1844925       rs4519708       320     51.8
1330    NA12335 1330    NA12340 IBD1    5       70.869  180.626 AFFX-SNP_7697354__rs276593      rs876154        738     109.8
1330    NA12335 1330    NA12340 IBD1    6       0.131   58.178  rs736864        rs3863230       482     58.0
1330    NA12335 1330    NA12340 IBD1    6       88.258  170.736 rs3778671       rs734249        560     82.5
```

**Formats:** `StartMB`/`StopMB` are **`%.3f`** (`51.799`, `0.080`, `242.590` — note the trailing zero is preserved, so it is fixed-precision not `%g`). `Length` is **`%.1f`** (`44.1`, `58.0`). `Chr` and `N_SNP` are `%d`. SNP names are bare strings and **can contain no whitespace but do contain `-` and `__`** (`AFFX-SNP_7697354__rs276593`).

Header note: the manual says "The header of `zcat ex.segments.gz` looks like this" — so the gz file **does** carry a header line.

### 7.4 IBD segment plotting helper

```
  prompt> Rscript king_segments_plot.R ex ibdseg
```
> All pairs of close relatives are plotted in their own files/plots, and all plots/files are gzipped together in a single file `ex_ibdseg_rplots.tar.gz`.

Requires R packages `ggplot2` and `parallel`.

---

## 8. INTEGRATED INFERENCE — `--related`

**Recommended default.** "The largest dataset we have successfully analyzed on a single server using the `--related` option consists of ~10 million samples (i.e., ~50,000,000,000,000 pairs!)."

**Documented example command lines:**

```
  prompt> king -b ex.bed --related
  prompt> king -b ex.bed --related --degree 2 --rplot --prefix ex
  prompt> king -b ex.bed --related --degree 2
```

> The first command identifies close relatives up to the first degree, and the second command specifies close relatives up to the second degree.
> `--related --degree 2` specifies that only related pairs (up to the 2nd-degree in this case) between families are included in the output. Specifically all pairs across families with a kinship coefficient less than **0.0884** will be excluded from the output.
> `--related` without the `--degree` option is highly recommended. Although distant relatedness that is higher than 2 is allowed, no fast algorithm is available at the moment and computation is substantially slower than `--related --degree 2`.

### 8.1 Output files

`{prefix}.kin` and `{prefix}.kin0` — **same file names as `--kinship` but a WIDER schema.** This is a trap for reimplementation: the `.kin` written by `--related` has 16 columns; the `.kin` written by `--kinship` has 10.

Observed runtime message: `Between-family relatives (kinship >= %.5lf) saved in file %s` — the cutoff is printed with **5 decimal places** (e.g. `0.08839`).

### 8.2 `{prefix}.kin` from `--related` — columns, in order (16)

```
FID  ID1  ID2  N_SNP  Z0  Phi  HetHet  IBS0  HetConc  HomIBS0  Kinship  IBD1Seg  IBD2Seg  PropIBD  InfType  Error
```

Verbatim column definitions given on the site (note: the site's definition list **omits `HetConc` and `HomIBS0`** even though they are in the header — a documentation gap):

```
FID: Family ID for the pair
ID1: Individual ID for the first individual of the pair
ID2: Individual ID for the second individual of the pair
N_SNP: The number of SNPs that do not have missing genotypes in either of the individual
Z0: Pr(IBD=0) as specified by the provided pedigree data
Phi: Kinship coefficient as specified by the provided pedigree data
HetHet: Proportion of SNPs with double heterozygotes (e.g., AG and AG)
IBS0: Porportion of SNPs with zero IBS (identical-by-state) (e.g., AA and GG)
Kinship: Estimated kinship coefficient (φ) from the SNP data
IBD1Seg: Total length of IBD1 segments divided by total length of all segments, estimate of π1=Pr(IBD=1)
IBD2Seg: Total length of IBD2 segments divided by total length of all segments, estimate of π2=Pr(IBD=2)
PropIBD: Proportion of genomes shared identical-by-descent, estimated by IBD2Seg + IBD1Seg/2, estimate of π=π2+π1/2
InfType: Inferred relationship type, such as Dup/MZTwin, PO, FS, 2nd, 3rd, 4th, UN
Error: Flag Indicating differences between inferred and reported relationship (1 for error, 0.5 for warning)
```

**Undocumented columns, inferred from the embedded R plot code (which plots them by name):**
- `HetConc` — **heterozygote concordance rate** (x-axis of the "Kinship vs Heterozygote Concordance" plot). Same statistic `--minConc` thresholds for `--duplicate`.
- `HomIBS0` — **IBS0 restricted to homozygote-informative SNPs** (x-axis of a second diagnostic plot). Related to the QC page's "informative" PO definition ("at least one carries the minor homozygote").

**Documented example block (verbatim):**

```
FID     ID1     ID2     N_SNP   Z0      Phi     HetHet  IBS0    HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType Error
Y001    NA18484 NA18486 18250   0.000   0.2500  0.2324  0.0002  0.3368  0.0006  0.2515  0.9750  0.0000  0.4875  PO      0
Y001    NA18484 NA18488 18249   0.000   0.2500  0.2332  0.0002  0.3379  0.0004  0.2522  1.0000  0.0000  0.5000  PO      0
Y001    NA18486 NA18488 18270   1.000   0.0000  0.2141  0.1053  0.3036  0.2460  0.0039  0.0000  0.0000  0.0000  UN      0
Y002    NA18485 NA18487 18276   0.000   0.2500  0.2349  0.0002  0.3413  0.0005  0.2541  1.0000  0.0000  0.5000  PO      0
Y002    NA18485 NA18489 18275   0.000   0.2500  0.2274  0.0003  0.3298  0.0007  0.2474  1.0000  0.0000  0.5000  PO      0
Y002    NA18487 NA18489 18269   1.000   0.0000  0.2098  0.1120  0.2999  0.2601  -0.0157 0.0000  0.0000  0.0000  UN      0
Y003    NA18497 NA18498 18262   0.000   0.2500  0.2232  0.0003  0.3254  0.0009  0.2448  0.9897  0.0000  0.4949  PO      0
Y003    NA18497 NA18499 18258   0.000   0.2500  0.2310  0.0002  0.3387  0.0006  0.2525  0.9612  0.0000  0.4806  PO      0
Y003    NA18498 NA18499 18280   1.000   0.0000  0.2100  0.1043  0.2997  0.2452  0.0015  0.0000  0.0000  0.0000  UN      0
```

**Formats here:** `Z0` = `%.3f`; **`HetHet` = `%.4f`** (`0.2324` — *differs from the `%.3f` in `--kinship`'s `.kin`!*); `Phi`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`, `IBD1Seg`, `IBD2Seg`, `PropIBD` all `%.4f`. `N_SNP` `%d`. `Error` integer-looking.

### 8.3 `{prefix}.kin0` from `--related`

Not shown verbatim on the site, but `--rplot` prose says plots use "output files ex.kin and ex.kin0", and the embedded R reads both with the same column names (`Kinship`, `PropIBD`, `IBD1Seg`, `IBD2Seg`, `InfType`, `HetConc`, `HomIBS0`). Expected schema by analogy (drop `Z0`, `Phi`, `Error`; split FID):

```
FID1  ID1  FID2  ID2  N_SNP  HetHet  IBS0  HetConc  HomIBS0  Kinship  IBD1Seg  IBD2Seg  PropIBD  InfType
```

**⇒ MUST be confirmed by running the binary. Flagged as a gap.**

### 8.4 Classification thresholds embedded in KING's own R output

KING writes these expressions verbatim into its generated `*_relplot.R`. They are the *plotting* classifier and reflect the documented degree boundaries (powers of 2: 2^-1.5 = 0.35355, 2^-2.5 = 0.17678, 2^-3.5 = 0.08839, 2^-4.5 = 0.04419):

```
allpair <- data$PropIBD>0 | data$Kinship>0.04419
d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08
d2 <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)
d3 <- data$PropIBD>0.08839 & data$PropIBD<=0.17678
d4 <- data$PropIBD>0.04419 & data$PropIBD<=0.08839
dU <- data$PropIBD<=0.04419
dU <- data$PropIBD>0 & data$PropIBD<=0.04419
```

Also observed as C-level messages:
- `1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf`
- `Cutoff value for IBS0 between FS and PO is set at %.4f`

**⇒ The PO-vs-FS split is an IBS0 threshold that KING computes and prints at `%.4f`. Capture the actual value at runtime.**

Plot titles observed (for `--rplot` parity):
- `"Kinship vs IBS0 in %s Families"`
- `"Kinship vs Heterozygote Concordance In %s Families"`
- `"Kinship vs Proportion IBD (Corr=", round(cor(...),digit=3),") in %s Families"`
- `"IBD Segments In %s Families"` / `"IBD Segments In Inferred %s Relatives"`

---

## 9. `--duplicate`

> implements the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time is in seconds, unless the number of samples is > 1,000,000 in which case a few minutes may be needed.

**Documented example command lines:**

```
  prompt> king -b ex.bed --duplicate
  prompt> king -b ex.bed,mystudy.bed --duplicate
  prompt> king -b ex.bed,ex.bed --fam ex.fam,ex2.fam --duplicate --rplot
```

**Output file:** `{prefix}.con` (from binary string constant `.con`; the website never names it).

Observed runtime messages (both a 32-bit and 64-bit count variant exist):
```
%d pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s
%lli pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s
```

⇒ the concordance threshold is printed as an **integer percent** (`> 80%`), driven by `--minConc` (default 0.80, 2.3.1+).

**Column headers for `.con` are NOT documented on the website — GAP, must come from the binary.**

`--rplot` for `--duplicate` emits `{prefix}_duplicateplot.R` (default `king_duplicateplot.R`) and needs R package `igraph`.

Toy-dataset recipe from the visualization page (verbatim):
```
  prompt> head -232 ex.fam > ex2.fam
  prompt> awk 'NR>232 && NR%2==1' ex.fam >> ex2.fam
  prompt> awk 'NR>232 && NR%2==0' ex.fam >> ex2.fam
```

---

## 10. OTHER RELATEDNESS OPTIONS

### 10.1 `--ibs`

> provides summary statistics such as the counts of IBS0, IBS1, IBS2, the average of IBS, in additional to the kinship estimates.

```
  prompt> king -b ex.bed --ibs
```

Observed message: `Between-family IBS data saved in file %s`. **Column headers NOT documented — GAP.** Historically `--ibs` widens `.kin`/`.kin0` with `N_IBS0 N_IBS1 N_IBS2 Dist`-style columns; must be confirmed against the binary.

### 10.2 `--homog` (RETIRED in 2.3.0 — not in the 2.3.2 binary)

> estimates pair-wise kinship coefficients assuming a homogeneous population. The best application of `--homog` may be for the linear mixed models (LMM) ... Although `--homog` is not recommended as a good method to infer relatedness in general populations, it provides inference results comparable to multiple alternative methods.

```
  prompt> king -b ex.bed --homog
```

Retirement note (2.3.0): "--homog is retired since it does not do a good job on either relationship inference or GRM". **Do not implement.** Binary still contains the string `Please do not run --ibdseg together with --homog` (dead code path).

### 10.3 `--unrelated`

```
  prompt> king -b ex.bed --unrelated
  prompt> king -b ex.bed --unrelated --degree 2
```

> This example estimates relatedness in the data first, followed by extracting a list of individuals that contains no pairs of individuals with a 1st- or 2nd-degree relationship. The detailed algorithm is described in this reference: **Manichaikul et al. 2012** [PDF → `https://www.chen.kingrelatedness.com/publications/pdf/PLoS8e1002640.pdf`]

**Output files** (from binary strings):
- `{prefix}unrelated.txt` — observed message: `A list of %d unrelated individuals saved in file %s`
- `{prefix}unrelated_toberemoved.txt` — observed message: `An alternative list of %d to-be-removed individuals saved in file %s`

Default prefix ⇒ `kingunrelated.txt`, `kingunrelated_toberemoved.txt`.
2.2.9: "`--cluster` and `--unrelated` now allow `--degree 3` and higher."
**Column headers NOT documented — GAP** (likely bare `FID IID`, no header).

### 10.4 `--build`

> reconstructs pedigrees using SNP data without the need of specifying pedigrees (although the pedigree information can still be incorporated)

```
  prompt> king -b ex.bed --build
  prompt> king -b ex.bed --build --degree 2
```

> The output includes two files: **`kingupdateids.txt`** and **`kingupdateparents.txt`**.

Downstream usage documented verbatim:
```
  prompt> plink1.9 --bfile ex --update-ids kingupdateids.txt --make-bed --out ex2
  prompt> plink1.9 --bfile ex2 --update-parents kingupdateparents.txt --make-bed --out ex3
```

⇒ These files must match **PLINK 1.9's** `--update-ids` and `--update-parents` formats exactly:
- `--update-ids`: `oldFID oldIID newFID newIID` (4 cols, no header)
- `--update-parents`: `FID IID newPatID newMatID` (4 cols, no header)

Observed message: `Update-ID information is saved in file %s`.

> The current `--build` algorithm connects all 1st-degree relatives with high accuracy. Known scenarios that `--build` does well are families that consist of at least a pair of full siblings, and/or a parent-child trio, etc.

2.1.8: "`--build` has improved and KING pedigree reconstruction can now incorporate 2nd-degree relatedness inference."
2.1.6: "`--build` can now take a rather small sample size (which crashes in previous versions when N<100)".

`--rplot` with `--build` emits `{prefix}_buildplot.R`; needs R package `kinship2`.
Also observed: `{prefix}splitped.txt`, `{prefix}_pedplot.R`.

### 10.5 `--cluster`

> is both a standalone parameter and a parameter to go with other options. As a standalone option, it clusters relatives into families by generating an updateid file which can then be used to update the pedigrees (e.g., using PLINK `--update-ids`). `--cluster` can also be used to group cyptic relatives together prior to association analysis

```
  prompt> king -b ex.bed --cluster
  prompt> king -b ex.bed --cluster --tdt
  prompt> king -b ex.bed --cluster --bySNP
  prompt> king -b ex.bed --cluster --bysample
  prompt> king -b ex.bed --prefix ex --cluster --degree 2 --rplot
```

Outputs: `{prefix}updateids.txt`; observed constant `cluster.kin`. `--rplot` emits `{prefix}_clusterplot.R` (needs `igraph`).

---

## 11. QUALITY CONTROL (`QC.shtml`)

### 11.1 `--bySNP` → `{prefix}bySNP.txt` (default `kingbySNP.txt`)

```
  prompt> king -b ex.bed --bySNP
  prompt> king -b ex.bed --cluster --bySNP
```

Columns, in order (verbatim definitions):

```
SNP: SNP name
Chr: Chromosome number of the SNP
Pos: Position of the SNP
Label_A: Label of the reference allele
Label_a: Label of the alternative allele
Freq_A: Frequency of the reference allele
N: Total number of samples with non-missing genotypes
N_AA: Total number of samples with genotype AA
N_Aa: Total number of samples with genotype Aa
N_aa: Total number of samples with genotype aa
CallRate: Proportion of samples with non-missing genotypes
N_MZ: Total number of MZ twins or duplicates
N_errMZ: Total number of inconsistencies between duplicates
Err_InMZ: Error rate in duplicates
N_PO: Total number of parent-offspring (PO) pairs
N_HomPO: Total number "informative" PO pairs (at least one carries the minor homozygote)
N_errPO: Total number of Medelian inconsistencies (MI) (AA->aa or aa->AA) in PO pairs 
Err_InPO: Error rate in PO pairs (N_errPO / N_PO)
Err_InHomPO: N_errPO / N_HomPO
N_trio: Total number of parent-offspring (PO) trios
N_HetOff: Total number of heterozygote offspring
N_errTrio: Total number of Medelian inconsistencies (MI) (AA x aa -> Aa) in PO trios
Err_InTrio: Error rate in PO trios (N_MIt / N_trio)
Err_InHetTrio: N_Mit / N_Het
```

Header line, in order (24 columns):
```
SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate N_MZ N_errMZ Err_InMZ N_PO N_HomPO N_errPO Err_InPO Err_InHomPO N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio
```
(Docs typos preserved: "Medelian", and the inconsistent `N_MIt`/`N_Mit`/`N_Het` references in the last two definitions.)
**No example block is given — number formats are a GAP.**
Observed message: `QC statistics by SNPs saved in file %s`.

### 11.2 `--bysample` → `{prefix}bySample.txt` (default `kingbySample.txt`)

Note the **capital `S` in the filename** vs the lowercase `s` in the flag.

```
  prompt> king -b ex.bed --bysample
  prompt> king -b ex.bed --cluster --bysample
```

Columns, in order (verbatim definitions):

```
FID: Family ID
IID: Individual ID
FA: Father ID
MO: Mother ID
SEX: Sex
N_SNP: Total number of non-missing SNPs on autosomes
Missing: SNP missing rate on autosomes
Heterozygosity: Heterozygosity on autosomes
N_Pair: Total number of SNPs that are not missing for the parent-offspring (PO) pair that the individual is involved
N_MIp: Total number of Mendelian inconsistencies (MI) (AA -> aa or aa -> AA) in the PO pair
Err_MIp: Error rate in the PO pair
N_trio: Total number of SNPs that are not missing for the PO trio
N_MIt: Total number of MIs (AA x aa -> Aa) in the PO trio
Err_MIt: Error rate in the PO trio
MI_Removal: Flag for removal 
```

Header line, in order (15 columns):
```
FID IID FA MO SEX N_SNP Missing Heterozygosity N_Pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal
```

**Corroborating binary string:** `FID IID FA MO SEX N_SNP Missing Heterozygosity` — a **space-separated** prefix of exactly this header. Interesting: this suggests `bySample.txt` may be *space*-delimited, unlike the `.kin`/`.seg` family. Must confirm.

Observed message: `QC statistics by samples saved in file %s`.
Bug history: 2.2.9 "The bug that crashed `--bySample` before is now fixed"; 2.3.1 "A bug in `--bySNP` and `--bysample` is fixed, for the scenario of presence of families".

### 11.3 `--autoQC`

> performs a straightforward QC pipeline, including sample-level QC (at call rate 95% by default, or a different call rate set by `--callrateN`), SNP-level QC (at call rate 95% by default, or a different call rate set by `--callrateM`), and gender QC. This analysis generates a list of SNPs to be removed, and a list of samples to be removed.

```
  prompt> king -b ex.bed --autoQC
```

Output files (names from binary string constants; website names none):
- `{prefix}_autoQC_Summary.txt`
- `{prefix}_autoQC_sampletoberemoved.txt`
- `{prefix}_autoQC_snptoberemoved.txt`
- `{prefix}_autoQC_updatesex.txt`
- `{prefix}_gender_autodata.txt`, `{prefix}_gender_autoplot.R`, `{prefix}_gender_autoplot.png`

Constraint observed: `Please do not run --ibdseg together with --autoQC`.
**Column headers — GAP.**

### 11.4 `--roh`

```
  prompt> king -b ex.bed --roh
  prompt> king -b ex.bed --roh --rplot
```

> Inbreeding coefficient for each sample is generated in file **`king.roh`**, and the exact ROH segments are saved in a gzipped file **`king.rohseg.gz`**.

#### `{prefix}.roh` — columns, in order

```
FID  ID  FA  MO  SEX  MaxROH  FInbred
```

**Documented example block (verbatim):**

```
FID     ID      FA      MO      SEX     MaxROH  FInbred
1328    NA06984 0       0       1       0.0     0.0000
1328    NA06989 0       0       2       0.0     0.0000
1330    NA12335 NA12340 NA12341 1       0.0     0.0000
1330    NA12336 NA12342 NA12343 2       0.0     0.0000
1330    NA12340 0       0       1       0.0     0.0000
1330    NA12341 0       0       2       0.0     0.0000
1330    NA12342 0       0       1       31.3    0.0449
1330    NA12343 0       0       2       0.0     0.0000
1334    NA10846 NA12144 NA12145 1       0.0     0.0000
```

Formats: `MaxROH` = **`%.1f`** (Mb), `FInbred` = **`%.4f`**, `SEX` = `%d`, missing parent = literal `0`.
Note the ID column here is named **`ID`**, not `IID`.
Observed message: `Run of homozygosity summary saved in file %s`.

#### `{prefix}.rohseg.gz` — columns, in order

```
FID  ID  Chr  StartMB  StopMB  StartSNP  StopSNP  N_SNP  Length
```

**Documented example block (verbatim):**

```
FID     ID      Chr     StartMB StopMB  StartSNP        StopSNP         N_SNP   Length
1330    NA12342 5       70.869  97.455  AFFX-SNP_7697354__rs276593      rs10866786      156     26.6
1330    NA12342 5       136.510 167.849 rs11745163      rs582906        224     31.3
1330    NA12342 6       25.472  31.787  rs13215347      rs805286        69      6.3
1346    NA10852 2       30.803  51.636  rs2681682       rs2698026       162     20.8
1459    NA12874 1       148.175 247.083 rs1868992       rs12058711      692     98.9
Y045    NA19201 5       70.869  85.089  AFFX-SNP_7697354__rs276593      rs10063186      103     14.2
Y057    NA19224 10      90.301  105.270 rs7901991       rs12268628      96      15.0
Y079    NA19113 17      0.116   10.842  rs4247500       rs4792080       114     10.7
```

Same numeric conventions as `.segments.gz`: `StartMB`/`StopMB` `%.3f`, `Length` `%.1f`.
This is exactly `.segments.gz` minus `FID2/ID2/IBDType`, confirming a shared segment writer.

Embedded R confirms the filename construction and the column subset used:
```
segments_name <- paste0(prefix, ".rohseg.gz")
all_seg_name <- paste0(prefix, "allsegs.txt")
segments <- subset(segments, select = c(FID, ID, Chr, StartMB, StopMB))
if( !(file.exists(seg_name) & file.exists(segments_name) & file.exists(all_seg_name)) ) stop("Missing RoH files")
```

ROH plot threshold (visualization page): plots individuals "with proportion of their genomes being ROH > **4.4%**, which corresponds to being offspring of parents that are 2nd-degree or closer." Requires R package `ggplot2`; emits `{prefix}_rohplot.R`.

---

## 12. POPULATION STRUCTURE (`kingpopulation.shtml`)

### 12.1 `--mds`

```
 prompt> king -b ex.bed --mds
 prompt> king -b ex.bed --mds --rplot
```

> Top principal components / ancestry coordinates (**20 by default**) are saved in files `kingpc.txt`.

**⚠️ CONTRADICTION:** this prose says 20; `--pcs` is documented as defaulting to 10 in both the flag index and the same page's "OTHER PARAMETERS" section ("The default pcs is 10"); the ancestry page's real screen printout says "**10** principal components saved in file expc.txt". The example block below nonetheless shows **PC1..PC20**. Treat **10** as the default and the 20-column example as pre-2.2.4 output (`--pcs` was only added in 2.2.4). **Confirm at runtime.**

### 12.2 `{prefix}pc.txt` — columns, in order

```
FID IID FA MO SEX AFF PC1 PC2 ... PCn
```

> Each row provides summary information for a sample. The top 10 principal components / ancestry coordinates are in the 7th to the 16th columns.

**Documented example block (verbatim — note this one is SPACE-delimited in the docs, single spaces, unlike the tab-style blocks above):**

```
FID IID FA MO SEX AFF PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20
1328 NA06984 0 0 1 1 -0.0545 0.0117 -0.0179 0.0081 -0.0293 0.0126 -0.0077 0.0143 -0.0061 0.0159 -0.0055 0.0260 -0.0184 0.0079 -0.0121 0.0143 0.0024 -0.0112 0.0204 -0.0265
1328 NA06989 0 0 2 1 -0.0542 0.0031 -0.0030 0.0115 0.0070 -0.0110 -0.0242 -0.0006 0.0078 0.0079 -0.0094 0.0137 0.0087 0.0036 -0.0299 -0.0031 -0.0149 -0.0054 0.0348 -0.0082
1330 NA12335 NA12340 NA12341 1 1 -0.0550 0.0063 -0.0353 -0.0021 0.1184 0.0747 -0.0337 -0.1091 -0.0734 0.0203 0.0146 0.0174 -0.1601 -0.0513 -0.0819 -0.0141 -0.0115 -0.0557 0.0547 0.0286
1330 NA12336 NA12342 NA12343 2 1 -0.0548 0.0380 -0.0058 0.0276 -0.0665 -0.0796 0.0319 0.0224 -0.1627 -0.0613 -0.1429 -0.1600 0.0735 -0.0596 0.0093 -0.0936 -0.1194 -0.1304 -0.0086 -0.0362
1330 NA12340 0 0 1 1 -0.0549 0.0095 -0.0274 -0.0035 0.0664 0.0375 -0.0347 -0.0736 -0.0648 0.0177 0.0151 -0.0015 -0.1308 -0.0444 -0.0518 -0.0305 -0.0007 -0.0363 0.0370 0.0228
1330 NA12341 0 0 2 1 -0.0534 0.0020 -0.0204 0.0067 0.0910 0.0755 -0.0168 -0.0788 -0.0445 0.0119 0.0090 0.0288 -0.0859 -0.0281 -0.0617 0.0051 -0.0188 -0.0488 0.0261 0.0324
1330 NA12342 0 0 1 1 -0.0547 0.0295 -0.0072 0.0189 -0.0535 -0.0451 0.0273 0.0192 -0.1190 -0.0273 -0.1062 -0.1218 0.0484 -0.0391 0.0045 -0.0634 -0.0821 -0.0998 0.0171 -0.0216
1330 NA12343 0 0 2 1 -0.0546 0.0219 -0.0140 0.0202 -0.0423 -0.0595 0.0185 0.0123 -0.1038 -0.0466 -0.0910 -0.1054 0.0649 -0.0589 0.0123 -0.0600 -0.0764 -0.0885 -0.0220 -0.0264
1334 NA10846 NA12144 NA12145 1 1 -0.0561 0.0068 -0.0126 -0.0318 -0.0476 0.0463 -0.0614 0.0760 0.0276 0.0580 0.0319 -0.0216 0.0215 -0.0286 0.0122 -0.0671 0.0552 0.1801 -0.1843 0.0054
```

All PC values are **`%.4f`**. Corroborating binary string: `FID IID FA MO SEX AFF`.

### 12.3 `--pca`

```
 prompt> king -b ex.bed --pca
 prompt> king -b ex.bed --pca --rplot
```
> Please run LD-pruning prior to PCA analysis.
> The top 10 pincipal components / ancestry coordinates are saved in files `kingpc.txt`, which has the same format as `kingpc.txt` from the `--mds` analysis.

Also observed: `{prefix}_eigenvalue.txt`, `{prefix}_pcplot.R`, `{prefix}_Dist.txt`, `{prefix}_popref.txt`.
Observed messages: `%d principal components saved in file %s`, `Distances to %d reference samples saved in file %s`, `Population distances saved in file %s`.

Warning verbatim: "Precompile KING binaries with versions lower than 2.2.3 are not suitable for population structure analysis in larger datasets for lacking LAPACK libraries."
Constraint observed: `Please do not run --ibdseg together with --pca`.

### 12.4 `--makeGRM` (2.3.0+; UNDOCUMENTED on the website)

Banner group: "Pairwise Relatedness Inference". Observed strings:
```
--makeGRM
--grm %s
--grm is misspecified: king is default, lmm is alternative.
--grm-lmm
--ibdGRM cannot run without ZLIB
GRM saved in file %s
GRM starts at %s
GRM ends at %s
Preparing GRM-%s at %s
PCA of GRM starts at %s
R plot for --makegrm is not available.
  %d 1st-degree relatives found in GRM.
  No 1st-degree relatives found in GRM.
GRM_Raw
diagGRM = %.3lf
```
Output files: `{prefix}_grm.txt`, and for the LMM variant `{prefix}_grm_grm-lmm.txt`, `{prefix}_eigenvalue_grm-lmm.txt`.
⇒ `--grm` takes a value: **`king`** (default) or **`lmm`**. Requires ZLIB for the IBD-GRM path.
**Everything else about this option is a GAP.**

---

## 13. ANCESTRY INFERENCE (`ancestry/`)

Requires R with the **`e1071`** package (SVM).

Reference files: `KGref.bed.xz` (489 MB), `KGref.fam.xz` (3 KB), `KGref.bim.xz` (37 MB).
Superpopulation labels: **`AFR`, `AMR`, `EAS`, `EUR`, `SAS`**.

**Documented example command lines (verbatim):**

```
 prompt> king -b KGref.bed,ex.bed
 prompt> king -b KGref.bed,ex.bed --pca --projection --rplot
 prompt> king -b KGref.bed,ex.bed --pca --projection --rplot --prefix ex
 prompt> king -b KGref.bed,ex.bed --pca --projection --pngplot
```

### 13.1 Full screen printout (verbatim, KING 2.2.7) — the best available model of stdout

```
KING starts at Thu May 27 18:29:56 2021
Read in PLINK fam files
        ../data/KGref.fam...
        ../data/ex.fam...
  PLINK pedigrees loaded: 2741 samples
Read in PLINK bim files
        ../data/KGref.bim...
        ../data/ex.bim...
  Genotype data consist of 16824 autosome SNPs
  PLINK maps loaded: 16824 SNPs
Read in PLINK bed files
        ../data/KGref.bed...
        ../data/ex.bed...
  PLINK binary genotypes loaded: 2741 samples
  KING format genotype data successfully converted

Options in effect:
        --pca
        --projection
        --rplot
        --prefix ex

PCA projection starts at Thu May 27 18:30:43 2021
2409 1000 Genomes samples are detected and used as reference.
Preparing matrix (2409 x 2409) for PCA...
  16824 SNPs are used in PCA.
SVD starts at Thu May 27 18:30:44 2021
  LAPACK is being used...
Largest 10 eigenvalues: 2059.73 1366.88 700.55 618.08 276.33 265.47 238.72 217.65 198.71 189.81
Projecting 332 samples starts at Thu May 27 18:30:44 2021
PCA projection ends at Thu May 27 18:30:44 2021
10 principal components saved in file expc.txt
Ancestry populations are inferred as in ex_InferredAncestry.txt
Ancestry plots are generated in ex_ancestryplot.pdf
KING ends at Thu May 27 18:30:56 2021
```

**Stdout parity notes:** 8-space indent for input filenames; 2-space indent for summary lines; "Options in effect:" block uses 8-space indent; timestamps are C `ctime()`-style `Thu May 27 18:29:56 2021`; eigenvalues printed at `%.2f`, space-separated.

### 13.2 `{prefix}_InferredAncestry.txt` — columns, in order

```
FID  IID  PC1  PC2  Anc_1st  Pr_1st  Anc_2nd  Pr_2nd  Ancestry
```

Confirmed by the embedded R that writes it:
```
colnames(pred.out) <- c("FID", "IID", "PC1", "PC2", "Anc_1st", "Pr_1st", "Anc_2nd",  "Pr_2nd", "Ancestry")
write.table(pred.out, paste0(prefix, "_InferredAncestry.txt"), sep = "...
print(paste("Results are saved to", paste0(prefix, "_InferredAncestry.txt"), date()))
```

**Documented example block (verbatim):**

```
FID     IID     PC1     PC2     Anc_1st Pr_1st  Anc_2nd Pr_2nd  Ancestry
1328    NA06984 -0.011  0.0268  EUR     0.9934  AFR     0.0032  EUR
1328    NA06989 -0.0104 0.0276  EUR     0.9962  AFR     0.0019  EUR
1330    NA12335 -0.0109 0.0267  EUR     0.9948  AFR     0.0024  EUR
1330    NA12336 -0.0101 0.0277  EUR     0.9965  AFR     0.0019  EUR
1330    NA12340 -0.0105 0.0288  EUR     0.9958  AFR     0.0023  EUR
1330    NA12341 -0.0102 0.0265  EUR     0.9924  AFR     0.0036  EUR
1330    NA12342 -0.0106 0.0279  EUR     0.9958  AFR     0.002   EUR
1330    NA12343 -0.0104 0.0273  EUR     0.9953  AFR     0.0023  EUR
1334    NA10846 -0.0115 0.0271  EUR     0.9959  AFR     0.0019  EUR
```

**⚠️ This file is written by R's `write.table`, NOT by C++ printf.** Evidence: `-0.011` (3 dp) and `0.002` (3 dp) alongside `0.9934` (4 dp) — trailing zeros are **stripped**. That is R's `signif`/default numeric formatting, not `%.4f`. A byte-parity reimplementation must reproduce R-style significant-digit trimming for this file only.

Other ancestry outputs: `{prefix}_ancestryplot.R`, `{prefix}_ancestryplot.pdf`, `{prefix}_popdistplot.R`, `{prefix}_popref.txt`.
Observed messages: `Ancestry results saved in file %s`, `Format: FID IID Population`, `FID IID Population`, `FID IID HetProj HetRef MinDist Kinship Closest RefID`.
Alternative-reference standalone R code: <https://github.com/chenlab-uva/AncestryInference_KING>

---

## 14. ASSOCIATION MAPPING (`genemapping.shtml`)

```
  prompt> king -b ex.bed --tdt    
  prompt> king -b ex.bed --cov, --mtscore --maxP 5E-8 --invnorm    
```
(The `--cov,` in the second line is a typo on the site, reproduced verbatim.)

- `--tdt`: "implements the well-known Transmission/Disequilibrium Test for family data that consist of parent-affected child trios."
- `--lmm`: "implements a linear mixed models for association between a SNP and a quantitative trait ... especially for a lot of traits, e.g., in eQTL/pQTL/meQTL/mQTL analysis"
- `--gdt` (2.3.1+): family-based association, Chen et al. 2009.
- `--invnorm`, `--trait`, `--covariate`, `--maxP`, `--prefix`, `--cpus`.

Output files (binary strings): `{prefix}tdt.txt`, `{prefix}TDTuninfo.txt`, `{prefix}_gdt.txt`, `{prefix}_gdt_kinship.txt`, `{prefix}_gdt_ibdseg.txt`, `{prefix}_gdt_lmm.txt`, `{prefix}_gdt_ped.txt`, `{prefix}_linear.txt`, `{prefix}_lmm_disease.txt`, `{prefix}_lmmking_disease.txt`, `{prefix}_lmmpc.txt`, `{prefix}_novclmm.txt`, `{prefix}_poodt.txt`.
Observed message: `LMM scan results (lambda_GC=%.4lf) saved in file %s`.
**Column headers — GAP. Out of scope for relatedness.**

---

## 15. RISK PREDICTION (`riskprediction.shtml`)

```
  prompt> king -b ex.bed --risk --model model.txt --prevalence 0.004 --noflip
```

### 15.1 Model input file columns

```
SNP: SNP name
EA: effect allele
AF: allele frequency of the effect allele
WT: weight at the effect allele
CHR: chromosome of the SNP
POS: position of the SNP
OA: other allele
```

Example (verbatim):
```
SNP             EA      AF      WT      CHR     POS             OA
rs9273363       A       0.131   1.702   6       32626272        C
rs9271594       G       0.095   1.801   6       32591213        A
rs2187668       T       0.076   1.367   6       32605884        C
rs34850435      T       0.345   0.839   6       32583299        C
rs34303755      C       0.216   1.079   6       32450613        A
rs689           T       0.265   0.403   11      2182224         A
rs2290400       C       0.459   0.295   17      38066240        T
```

### 15.2 `{prefix}grs.txt` columns

```
FID: family ID
IID: individual ID
InfoSNP: call rate
InfoVar: proportion of GRS variance at non-missing SNPs
GRS: genetic risk score, as the original form of weighted sum
Zscore: GRS divided by the total GRS variance (a function of the model that is independent of the test data)  
Percent: Estimated percentage of GRS in the general population
ScaledGRS: transformed GRS, in the range (0, 1)
Status: given disease status, if the 6th column of the .fam file is available
```

### 15.3 Stdout summary table (verbatim)

```
Risk Cutoff     0.1     0.2     0.3     0.4     0.5     0.6     0.7     0.8     0.9
TruePositives   1015    987     965     944     923     909     889     857     751
FalsePositives  2616    1959    1524    1306    1126    1030    868     716     444
TrueNegatives   312     969     1404    1622    1802    1898    2060    2212    2484
FalseNegatives  6       34      56      77      98      112     132     164     270
Sensitivity     0.9941  0.9667  0.9452  0.9246  0.9040  0.8903  0.8707  0.8394  0.7356
Specificity     0.1066  0.3309  0.4795  0.5540  0.6154  0.6482  0.7036  0.7555  0.8484
Positive PV     0.0056  0.0072  0.0090  0.0103  0.0117  0.0126  0.0145  0.0170  0.0238
Negative PV     0.9997  0.9995  0.9994  0.9993  0.9992  0.9992  0.9991  0.9989  0.9984
AUC (Area under the ROC curve) = 0.8708
AUC among 1234 males is 0.8686
AUC among 2715 females is 0.8693
Genetic risk scores are saved in file exgrs.txt
```
Rates at `%.4f`; AUC at `%.4f`.

---

## 16. `--rplot` / `--pngplot` / `--rpath`

> `--rplot` generates R code first and then calls R program to make plots in a PDF file.
> `--pngplot` generates R code first and then calls R program to make plots in a PNG file for certain applications.
> `--rpath` specifies the full path of the R program in case "R" command without a full path cannot properly run.

**Key behavior for reimplementation:** KING **always writes the `.R` file**, and only *additionally* runs R if the required package is present. Verbatim: "otherwise only the R code `ex_buildplot.R` is generated without the actual plots in PDF."

Required R packages by analysis:

| Command | R package(s) | Generated R file |
|---|---|---|
| `--rplot` (bare, plots .fam pedigrees) | `kinship2` | `{prefix}_pedplot.R` |
| `--build --rplot` | `kinship2` | `{prefix}_buildplot.R` |
| `--ibdseg --degree 2 --rplot` | `igraph` | `{prefix}_uniqfamplot.R` |
| `--related --degree 2 --rplot` | `igraph` + `kinship2` | `{prefix}_relplot.R`, `{prefix}_uniqfamplot.R` |
| `--cluster --degree 2 --rplot` | `igraph` | `{prefix}_clusterplot.R` |
| `--duplicate --rplot` | `igraph` | `{prefix}_duplicateplot.R` |
| `--roh --rplot` | `ggplot2` | `{prefix}_rohplot.R` |
| `king_segments_plot.R` (standalone) | `ggplot2`, `parallel` | → `{prefix}_ibdseg_rplots.tar.gz` |
| `--pca --projection --rplot` (ancestry) | `e1071` | `{prefix}_ancestryplot.R` |

Documented plot commands (verbatim):
```
  prompt> king -b ex.bed --prefix ex --rplot
  prompt> king -b ex.bed --prefix ex --build --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --ibdseg --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --related --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --cluster --degree 2 --rplot
  prompt> king -b ex.bed --roh --rplot
```

Differences between `--ibdseg --rplot` and `--related --rplot` (verbatim):
> 1) `--related --rplot` only visualizes unique family configurations that are cryptic (between families);
> 2) `--related` is expected to be orders of magnitude faster; and
> 3) `--related --rplot` also visualizes (within-family) pedigree errors.

Also observed R-generator names: `_MIerrorplot.R`, `_ibd1vsibd2.R`, `_ibdmapplot.R`, `_aucmapplot.R`, `_herplot.R`, `_mthomoplot.R`, `_nplplot.R`, `_poprohplot.R`, `_popdistplot.R`, `_pcplot.R`.

---

## 17. MASTER OUTPUT-FILE CATALOG

`{p}` = value of `--prefix` (default `king`). Files marked ✅ have their columns documented on the website.

| File | Produced by | Columns documented? |
|---|---|---|
| `{p}.kin` | `--kinship` (10 col) / `--related` (16 col) | ✅ both |
| `{p}.kin0` | `--kinship` (8 col) / `--related` (?) | ✅ kinship only |
| `{p}X.kin`, `{p}X.kin0` | X-chr variants | ❌ |
| `{p}.con` | `--duplicate` | ❌ |
| `{p}.seg` | `--ibdseg`, `--related` | ✅ |
| `{p}.segments.gz` | `--ibdseg` | ✅ |
| `{p}X.seg` | `--ibdseg` (X-chr) | ❌ |
| `{p}allsegs.txt` | `--ibdseg`, `--roh` | ❌ |
| `{p}.roh` | `--roh` | ✅ |
| `{p}.rohseg.gz` | `--roh` | ✅ |
| `{p}bySNP.txt` | `--bySNP` | ✅ (names only, no example) |
| `{p}bySample.txt` | `--bysample` | ✅ (names only, no example) |
| `{p}_autoQC_Summary.txt` | `--autoQC` | ❌ |
| `{p}_autoQC_sampletoberemoved.txt` | `--autoQC` | ❌ |
| `{p}_autoQC_snptoberemoved.txt` | `--autoQC` | ❌ |
| `{p}_autoQC_updatesex.txt` | `--autoQC` | ❌ |
| `{p}unrelated.txt` | `--unrelated` | ❌ |
| `{p}unrelated_toberemoved.txt` | `--unrelated` | ❌ |
| `{p}updateids.txt` | `--build`, `--cluster` | ✅ implicitly (PLINK `--update-ids`) |
| `{p}updateparents.txt` | `--build` | ✅ implicitly (PLINK `--update-parents`) |
| `{p}splitped.txt` | `--build` | ❌ |
| `{p}pc.txt` | `--pca`, `--mds` | ✅ |
| `{p}_eigenvalue.txt` | `--pca`, `--mds` | ❌ |
| `{p}_Dist.txt`, `{p}_popref.txt` | `--pca --projection` | ❌ |
| `{p}_InferredAncestry.txt` | `--pca --projection --rplot` | ✅ |
| `{p}_grm.txt` | `--makeGRM` | ❌ |
| `{p}grs.txt` | `--risk` | ✅ |
| `{p}tdt.txt`, `{p}TDTuninfo.txt` | `--tdt` | ❌ |
| `{p}_gdt*.txt` | `--gdt` | ❌ |
| `{p}_lmm*.txt`, `{p}_linear.txt` | `--lmm` | ❌ |
| `{p}het.txt`, `{p}flip.txt`, `{p}_relatives.txt`, `{p}_relative_removed.txt` | misc | ❌ |
| `{p}_*.R`, `*.pdf`, `*.png` | `--rplot` / `--pngplot` | n/a |
| `{p}.ped/.dat/.map/.bed/.bim/.fam`, `{p}_pat.*`, `{p}_mat.*`, `{p}af.ped`, `{p}hap.ped` | `--plink` / format conversion | ❌ |

---

## 18. VERSION HISTORY — ENTRIES THAT AFFECT OUTPUT SEMANTICS

Only changes that could alter numbers or file layout are listed; the full changelog runs from 1.0.0 (Oct 5, 2010) to 2.3.2.

| Version | Date | Relevant change (verbatim) |
|---|---|---|
| 2.3.2 | Sept 8, 2023 | "Crashing when running `--related --degree 3` (or higher), a bug in 2.3.1 only is now fixed"; "`--noscreen` is newly available for `--related` and `--duplicate` analyses, for potentially more precise inference by skipping screening relatives with a subset of SNPs" |
| 2.3.1 | July 28, 2023 | `--gdt` added; "A bug in `--bySNP` and `--bysample` is fixed, for the scenario of presence of families"; `--minConc` added; "`--lmm` is work in progress" |
| 2.3.0 | Oct 10, 2022 | `--lmm` added; **`--makeGRM` added**; **`--homog` retired** |
| 2.2.9 | Sept 20, 2022 | "`--cluster` and `--unrelated` now allow `--degree 3` and higher"; `--bySample` crash fixed; `--risk` bug fixed |
| 2.2.8 | May 10, 2022 | minor `--pca --projection` fix |
| 2.2.7 | May 18, 2021 | `--ibdseg` bug from 2.2.5 "completly fixed" |
| 2.2.6 | Mar 23, 2021 | `--ibdseg` hang fixed; "`--pca --projection --rplot` can be used for ancestry inference" |
| 2.2.5 | June 5, 2020 | "A bug in `--unrelated` is fixed"; "**`--ibdseg` is substantially improved**"; "`--build` is improved" |
| 2.2.4 | Oct 11, 2019 | population structure much faster; **`--pcs` added, default 10** |
| 2.2.3 | Aug 9, 2019 | static LAPACK for Linux; Windows binaries; "**IBD segment analysis and Run of Homozygosity anaysis now apply to chromosome X as well**"; "A bug is fixed for `--related` analysis when the number of SNPs is less than 4096" |
| 2.2.2 | May 29, 2019 | fix for reading multiple datasets with identical SNP sets; **`--projection N` for `--kinship`/`--ibdseg` introduced** |
| 2.2.1 | May 14, 2019 | "A minor bug in `--ibdseg` and `--roh` is now fixed (regarding maxIBD1 and maxIBD2, not affecting the main inference)"; `--build` bug fixed |
| 2.2 | Mar 28, 2019 | visualization of families; duplicate IDs allowed between first (Reference) and other datasets |
| 2.1.8 | Feb 27, 2019 | `--build` incorporates 2nd-degree; `kinship2` pedigree drawing |
| 2.1.6 | Nov 28, 2018 | `--lessmem` retired; N>10,000,000 supported; multi-dataset input; **progress percentage printed to screen** |
| 2.1.5 | Aug 24, 2018 | `--sexchr` added (default 23); `--risk` fixed; `--ibdseg` up to 700,000 samples (was 256,000) |
| 2.1.4 | June 6, 2018 | "KING can now accurately infer up to the 4th-degree relatedness (`--related`, `--ibdseg`), while the original KING method (`--kinship`) remains accurate up to the second-degree relatedness (even in the presence of admixture)" |
| 2.1.3 | Feb 13, 2018 | `--ibdseg`, `--related`, `--roh` algorithms improved |
| 2.1.2 | Dec 14, 2017 | IBD segment algorithm improved |
| 2.1 | Oct 24, 2017 | **`--ibdseg` introduced**; `--related` made integrative |
| 2.0 | Oct 17, 2016 | multi-core; `--tdt`, `--mtscore` |
| 1.9 | Oct 10, 2015 | fast `--duplicate` / `--related` algorithm |
| 1.4.2 | 2013 | `--pca --projection` released; "reference samples need to be indicated as 1 in the phenotype (6th) column, and the study samples need to 2" |
| 1.4 | Dec 14, 2011 | `--unrelated` introduced |

**Implication:** IBD-segment numeric output changed materially at 2.1.2, 2.1.3, 2.2.1, 2.2.5, 2.2.6 and 2.2.7. Byte-parity is only meaningful against **2.3.2 specifically** — never against published example output from the website, which predates several of these.

---

## 19. REFERENCES CITED BY THE SITE

- **Primary (cite this):** Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) *Robust relationship inference in genome-wide association studies.* **Bioinformatics 26(22):2867-2873.**
  Abstract: <http://bioinformatics.oxfordjournals.org/content/26/22/2867.abstract>
  PDF: <https://www.chen.kingrelatedness.com/publications/pdf/BI26_2867.pdf>
- **`--unrelated` algorithm:** Manichaikul et al. 2012 — <https://www.chen.kingrelatedness.com/publications/pdf/PLoS8e1002640.pdf> (PLoS Genet 8:e1002640)
- **`--gdt`:** Chen et al. 2009
- **Cited by the 2.3.2 binary itself:** "Chen et al. 2024" — a newer paper not linked anywhere on the website.
- **`--risk` application:** Onengut-Gumuscu S, … Rich SS (2019) *Type 1 Diabetes Risk in African-Ancestry Participants and Utility of an Ancestry-Specific Genetic Risk Score.* Diabetes Care 42(3):406-415.

---

## 20. GAPS — WHAT THE WEBSITE DOES NOT ANSWER (must come from running the binary)

Ordered by importance for a relatedness-focused reimplementation.

1. **Field delimiter.** Website has zero tabs. Confirm `\t` vs space, per file. `bySample.txt`'s binary string is space-separated, suggesting it may differ from `.kin`/`.seg`.
2. **`.kin0` schema for `--related`** — not shown anywhere.
3. **`--duplicate` → `{p}.con` column headers** — completely undocumented.
4. **`--ibs` output columns** — the flag index promises "counts of IBS0, IBS1, IBS2, the average of IBS" but never names the columns.
5. **`HetConc` and `HomIBS0` exact definitions** — present in `--related`'s `.kin` header but omitted from the site's own column glossary.
6. **The PO-vs-FS IBS0 cutoff value** — KING prints it (`Cutoff value for IBS0 between FS and PO is set at %.4f`) but never documents how it is derived.
7. **Row ordering** within every output file (pair enumeration order, family order, sort keys).
8. **Header row presence** for `unrelated.txt`, `updateids.txt`, `updateparents.txt`, `.con`, autoQC files.
9. **Negative-zero and rounding behavior** — e.g. does `-0.0000` ever print?
10. **`%.3f` vs `%.4f` for `HetHet`** — genuinely differs between the `--kinship` writer and the `--related` writer per the site's own examples. Confirm.
11. **`--makeGRM` output format** entirely.
12. **`--pcs` true default** (10 vs the 20-column example).
13. **Exact stdout text for 2.3.2** — the only full printout on the site is 2.2.7.
14. **`--noscreen` and `--minConc` semantics** beyond the one-line changelog entries.
15. **X-chromosome outputs** (`X.kin`, `X.kin0`, `X.seg`) — undocumented.

---

## 21. CONSOLIDATED COMMAND-LINE CORPUS (every example on the site, verbatim)

```
  prompt> king -b ex.bed --related
  prompt> king -b ex.bed --fam ex.fam --bim ex.bim --related
  prompt> king -b ex.bed,mystudy.bed --duplicate
  prompt> king -b ex.bed,ex.bed --duplicate
  prompt> plink1.9 --vcf example.vcf.gz --make-bed --out ex
  prompt> king -b ex.bed --kinship
  prompt> king -b ex.bed --ibdseg
  prompt> king -b ex.bed --ibs
  prompt> king -b ex.bed --homog
  prompt> king -b ex.bed --duplicate
  prompt> king -b ex.bed --related --degree 2
  prompt> king -b ex.bed --unrelated
  prompt> king -b ex.bed --build
  prompt> king -b ex.bed --cluster
  prompt> king -b ex.bed --unrelated --degree 2
  prompt> king -b ex.bed --build --degree 2
  prompt> plink1.9 --bfile ex --update-ids kingupdateids.txt --make-bed --out ex2
  prompt> plink1.9 --bfile ex2 --update-parents kingupdateparents.txt --make-bed --out ex3
  prompt> king -b ex.bed --cluster --tdt
  prompt> king -b ex.bed --kinship --degree 2
  prompt> king -b subset1,subset2 --kinship --proj 100000 --prefix subset12
  prompt> king -b subset1.bed --kinship --prefix subset1
  prompt> king -b ex.bed --ibdseg --degree 3 --rplot --prefix ex
  prompt> Rscript king_segments_plot.R ex ibdseg
  prompt> king -b ex.bed --related --degree 2 --rplot --prefix ex
  prompt> king -b ex.bed --bySNP
  prompt> king -b ex.bed --cluster --bySNP
  prompt> king -b ex.bed --bysample
  prompt> king -b ex.bed --cluster --bysample
  prompt> king -b ex.bed --autoQC
  prompt> king -b ex.bed --roh
  prompt> king -b ex.bed --mds
  prompt> king -b ex.bed --mds --rplot
  prompt> king -b ex.bed --pca
  prompt> king -b ex.bed --pca --rplot
  prompt> king -b KGref.bed,ex.bed
  prompt> king -b KGref.bed,ex.bed --pca --projection --rplot
  prompt> king -b KGref.bed,ex.bed --pca --projection --rplot --prefix ex
  prompt> king -b KGref.bed,ex.bed --pca --projection --pngplot
  prompt> king -b ex.bed --tdt
  prompt> king -b ex.bed --cov, --mtscore --maxP 5E-8 --invnorm
  prompt> king -b ex.bed --risk --model model.txt --prevalence 0.004 --noflip
  prompt> king -b ex.bed --prefix ex --rplot
  prompt> king -b ex.bed --prefix ex --build --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --ibdseg --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --related --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --cluster --degree 2 --rplot
  prompt> king -b ex.bed,ex.bed --fam ex.fam,ex2.fam --duplicate --rplot
  prompt> king -b ex.bed --roh --rplot
```

Build instructions from `Download.shtml` (verbatim — useful only to know the toolchain, source not read):
```
wget https://www.kingrelatedness.com/KINGcode.tar.gz
tar -xzvf KINGcode.tar.gz
c++ -lm -lz -O2 -fopenmp -o king *.cpp
```
⇒ KING links **libm**, **libz** (zlib, for the `.gz` outputs), and uses **OpenMP**. Confirms `.segments.gz` / `.rohseg.gz` are standard gzip.

---

*Raw fetched HTML and converted text are preserved at `<scratchpad>/web03/`.*
