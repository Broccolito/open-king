# KING 2.3.2 — Website & Manual Recon (clean-room reference)

**Source:** https://www.kingrelatedness.com/ — every sub-page fetched as raw HTML on 2026-08-13.
**Purpose:** exact output formats, column headers, file names, option list for a clean-room MIT
reimplementation. Author: Wei-Min Chen (U. Virginia).

**Provenance of every fact below is marked:**
- `[WEB]` — published on kingrelatedness.com (documentation = fact about I/O contract).
- `[BIN]` — observed from the reference binary's *usage screen* or *embedded string constants*
  (facts about output format; **no source code was read**).
- `[INF]` — inference/derivation made by me from the above.

Raw HTML kept at `research/raw/*.shtml`; text renderings at `research/txt/*.txt`.

---

## 0. LEGAL NOTE (read this first)

`[WEB]` Download.shtml states verbatim:

> Feel free to use KING for your research, but please do not redistribute AND make profits.

KING is **not** open-source and carries **no OSI license**. Source is distributed as
`KINGcode.tar.gz` / `KING2.3.2code.tar.gz` under that restriction only. This confirms the
clean-room approach is mandatory: **do not fetch, open, or transcribe `KINGcode.tar.gz`.**
Everything in this document derives from published documentation and observed binary behavior.

The build line the site publishes (a documented fact about the toolchain, not source):
```
c++ -lm -lz -O2 -fopenmp -o king *.cpp
```
→ single-TU-per-file C++, OpenMP, zlib. Implies gzip output written via zlib.

**Citation** `[WEB]` (on every page):
> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) Robust relationship
> inference in genome-wide association studies. *Bioinformatics* 26(22):2867-2873

`[BIN]` The 2.3.2 binary's error footer additionally cites **"Chen et al. 2024"** (unpublished at
the time of the website's last update — the `--ibdseg` manuscript).

---

## 1. SITE MAP (all pages, all fetched)

| Page | URL | Last updated `[WEB]` | Covers |
|---|---|---|---|
| Home | `/index.shtml` | — | overview only |
| Relationship Inference (main manual) | `/manual.shtml` | July 28, 2023 | `--kinship --ibdseg --related --duplicate --ibs --homog --unrelated --build --cluster` |
| Visualization of Families | `/KINGvisualization.shtml` | May 21, 2019 | `--rplot` variants |
| Quality Control | `/QC.shtml` | February 21, 2018 | `--bySNP --bysample --autoQC --roh` |
| Population Structure | `/kingpopulation.shtml` | Oct 11, 2019 | `--mds --pca --pcs --projection` |
| Ancestry Inference | `/ancestry/` | May 28, 2021 | `--pca --projection --rplot` + SVM |
| Association Mapping | `/genemapping.shtml` | October 10, 2022 | `--tdt --mtscore --lmm` |
| Risk Prediction | `/riskprediction.shtml` | August 24, 2018 | `--risk --model --prevalence --noflip` |
| Complete Flag Index | `/flagindex.shtml` | — | one-line description per option |
| Download | `/Download.shtml` | September 8, 2023 | binaries, source, `ex.tar.gz` |
| Download History | `/history.shtml` | March 23, 2021 (stale) | per-version changelog |
| User Forum | `/forum/` | — | not scraped |

**Documentation is stale relative to 2.3.2.** manual.shtml still says "The latest version is KING
2.3.1". flagindex.shtml predates 2.3.0 (still lists `--homog`, `--mtscore`, `--lessmem`; missing
`--lmm --makeGRM --gdt --noscreen --minConc --plink --phefile --covfile --prunedsnp`). Trust the
binary over the site where they disagree — see §3.2.

---

## 2. GLOBAL CONVENTIONS

### 2.1 Input files `[WEB]` manual.shtml §GENERAL INPUT FILES

PLINK 1 binary format only: `.bed` + `.fam` + `.bim`.

```
  prompt> king -b ex.bed --related
  prompt> king -b ex.bed --fam ex.fam --bim ex.bim --related
```
> In the first example, although only ex.bed is specified, the other two input files are
> pre-assumed to be ex.fam and ex.bim.

- `-b` / `-bname` — the .bed file(s). **Also accepts a bare prefix** (see the `--projection`
  example `king -b subset1,subset2` with no `.bed` suffix). `[INF]` strip a trailing `.bed` if
  present, else treat the token as the prefix.
- `--fam`, `--bim` — override the derived paths; comma-separated in the same order as `-b`
  (`--fam ex.fam,ex2.fam` in the duplicate-visualization example).
- **Multiple datasets:** comma-separated, *no spaces*: `king -b ex.bed,mystudy.bed --duplicate`.
  Merging rules `[WEB]` ancestry page, verbatim:
  ```
  1. SNPs with unambiguous allele labels can be auto-flipped before merging
  2. SNPs with ambiguous allele labels (i.e., A/T, or C/G) are excluded
  3. SNPs with inconsistent allele labels (e.g., >3 alleles) are excluded
  ```
  > Note KING is only looking at the SNP names when merging multiple datasets and other
  > information such as chr and positions is not utilized.
  > all IDs (combinations of FID and IID) need to be unique within and across datasets.
  One exception: when most samples overlap across datasets, IDs are prefixed **`REF_`** or
  **`QRY_`** (e.g. `king -b ex.bed,ex.bed --duplicate`).
- `[BIN]` extra optional inputs not on the site: `--phefile`, `--covfile`, `--prunedsnp`.
- `[WEB]` genemapping.shtml: `ex.phe` (phenotypes) and `ex.cov` (covariates) are **auto-discovered**
  from the `-b` prefix — "KING searches for all 5 files automatically".
- No dosage input (`-d`) is documented anywhere on the site. `[INF]` `-d` is not a KING 2.3.2 flag.

**Memory model** `[WEB]`: `~N × M / 4` bytes plus small overhead (N samples, M SNPs);
100,000 × 1,000,000 ⇒ ~25 GB.

**Pre-processing guidance** `[WEB]`, important because it constrains expected accuracy:
> Please do not prune or filter any "good" SNPs that pass QC prior to any KING inference,
> unless the number of variants is too many to fit the computer memory, e.g., > 100,000,000 as in
> a WGS study, in which case rare variants can be filtered out. **LD pruning is not recommended in
> KING.**
Exception `[WEB]` kingpopulation.shtml: "Please run LD-pruning prior to PCA analysis."

VCF conversion the site recommends:
```
 prompt> plink1.9 --vcf example.vcf.gz --make-bed --out ex
```

### 2.2 `--prefix` and output-name composition `[WEB]` + `[BIN]`

Default prefix is **`king`** (`[BIN]` usage screen prints `--prefix [king]`).

Two distinct composition rules — **this matters for byte parity**:

| Rule | Pattern | Examples |
|---|---|---|
| **A. dotted suffix** | `{prefix}` + `.ext` | `king.kin`, `ex.seg`, `king.roh`, `king.rohseg.gz`, `ex.segments.gz` |
| **B. glued suffix** | `{prefix}` + `word` (no separator) | `kingpc.txt`/`expc.txt`, `kingbySNP.txt`, `kingbySample.txt`, `kingupdateids.txt`, `exgrs.txt` |
| **C. underscored** | `{prefix}` + `_word` | `ex_InferredAncestry.txt`, `ex_relplot.pdf`, `ex_buildplot.R` |

Rule B is confirmed by the ancestry page's screen printout, where `--prefix ex` yields
**`expc.txt`** (not `ex.pc.txt`, not `ex_pc.txt`) `[WEB]`, and by the binary's string table which
contains the bare fragments `pc.txt`, `bySNP.txt`, `bySample.txt`, `updateids.txt`,
`updateparents.txt`, `grs.txt` `[BIN]`.

Note the capitalization asymmetry: option `--bysample` (lower s) → file `kingbySample.txt`
(**capital S**); option `--bySNP` → file `kingbySNP.txt` `[WEB]` QC.shtml.

### 2.3 Delimiters and number formatting — **critical finding**

`[INF]` **The raw website HTML contains ZERO tab characters** (verified:
`perl -ne '$c += tr/\t//' raw/manual.shtml` → `0`). Every example block was pasted with tabs
expanded to spaces at **8-column tab stops**. Verify on the `king.kin` header:

```
FID·····ID1·····ID2·····N_SNP···Z0······Phi·····HetHet··IBS0····Kinship·Error
```
`FID`(3)+5sp=8, `ID1`(3)+5=8, `N_SNP`(5)+3=8, `IBS0`(4)+4=8, `Kinship`(7)+1=8 — a perfect
8-column grid. Data rows agree: `2359853`(7)+1=8, `-0.0108`(7)+1=8.

⇒ **`.kin`, `.kin0`, `.seg`, `.segments.gz`, `.roh`, `.rohseg.gz`, `*_InferredAncestry.txt`, the GRS
model file and the risk-prediction table are TAB-delimited (`\t`), one tab between fields, no
padding.** The website's visual alignment is an artifact — do **not** emit space padding.

⇒ **Counter-example: `{prefix}pc.txt` is genuinely SPACE-delimited.** Its rows read
`1328 NA06984 0 0 1 1 -0.0545 0.0117 ...` — single spaces, no 8-column grid. If it were
tab-expanded, `1328` would be followed by 4 spaces. `[INF]` MDS/PCA output uses `' '`, relatedness
output uses `'\t'`. Do not unify these.

**Numeric field widths observed** `[WEB]` (see per-file sections for the authoritative list):

| Field | Format | Evidence |
|---|---|---|
| `Z0` | `%.3f` | `0.000`, `1.000` |
| `Phi` | `%.4f` | `0.2500`, `0.0000` |
| `HetHet` (in `.kin`/`.kin0`) | `%.3f` | `0.162`, `0.120` |
| `HetHet` (in `.kin` from `--related`) | `%.4f` | `0.2324`, `0.2141` |
| `IBS0` (in `.kin`/`.kin0`) | `%.4f` | `0.0008`, `0.0634` |
| `HetConc`, `HomIBS0` | `%.4f` | `0.3368`, `0.0006` |
| `Kinship` | `%.4f` | `0.2459`, `-0.0108`, `0.1356` |
| `IBD1Seg`, `IBD2Seg`, `PropIBD` | `%.4f` | `0.9976`, `0.0000`, `0.4988` |
| `Error` | `%g`-like | `0` (and `1`, `0.5` per the docs) |
| `StartMB`, `StopMB` | `%.3f` | `51.799`, `247.083`, `0.080`, `0.116` |
| `Length` | `%.1f` | `44.1`, `98.9`, `109.8`, `6.3` |
| `MaxROH` | `%.1f` | `0.0`, `31.3` |
| `FInbred` | `%.4f` | `0.0000`, `0.0449` |
| PC columns in `{prefix}pc.txt` | `%.4f` | `-0.0545`, `0.0117` |
| `PC1`/`PC2` in `_InferredAncestry.txt` | **`%g` (trailing zeros stripped)** | `-0.011`, `0.0268`, `-0.0104` |
| `Pr_1st`/`Pr_2nd` | **`%g` (trailing zeros stripped)** | `0.9934`, `0.002` (not `0.0020`) |

⚠️ The `_InferredAncestry.txt` columns are **not** fixed-precision: row 7 shows `0.002` where a
`%.4f` would print `0.0020`, and `-0.011` where `%.4f` would print `-0.0110`. `[INF]` These are
written by the R helper via `write.table`, not by C++ `printf`.

`[BIN]` Confirmed C-style format specifiers seen in the string table include `%.4lf` and `%.5lf`
(e.g. `"Cutoff value for IBS0 between FS and PO is set at %.4f"`,
`"Between-family relatives (kinship >= %.5lf) saved in file %s"`).

### 2.4 Example dataset `[WEB]`

`ex.tar.gz` [1.35 MB] — **332 HapMap samples (165 CEU + 167 YRI), 18,290 SNPs**. Used throughout
the tutorial. Reference panel for ancestry: `KGref.bed.xz` [489 MB], `KGref.bim.xz` [37 MB],
`KGref.fam.xz` [3 KB].

---

## 3. COMPLETE OPTION LIST

### 3.1 Published flag index `[WEB]` `/flagindex.shtml` — verbatim, in the site's order

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

Note `--projection` is listed **twice** with two different meanings — see §3.4.

### 3.2 Actual 2.3.2 option menu `[BIN]` — verbatim from `king` with no arguments

This is the authoritative list for 2.3.2 and supersedes flagindex.shtml.

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

Notes on this block, for parity of our own `--help`:
- Label column is right-aligned to width **34**, then `" : "`, then the option list.
- Continuation lines indent to the same column-37 start.
- Defaults are printed in `[...]`; empty string options print `[]`.
- **`--noscreen [-1717986816]`** — that is `0x9999...`-style garbage; the flag is a *boolean* whose
  integer slot is printed uninitialized. A KING bug, not a meaningful default. Our reimplementation
  should print `--noscreen` with no default.
- Note the literal **tab characters** in the "Chen et al. 2024," line of the footer, and the
  trailing space after `FATAL ERROR - `.

### 3.3 Delta: website vs 2.3.2 binary `[INF]`

| Option | flagindex.shtml | 2.3.2 binary | Notes |
|---|---|---|---|
| `--homog` | present | **gone** | retired in 2.3.0 `[WEB]` history: "does not do a good job on either relationship inference or GRM". manual.shtml still documents it. |
| `--mtscore` | present | **gone** | replaced by `--lmm` in 2.3.0 |
| `--lessmem` | present (marked retired) | gone | retired in 2.1.6 |
| `--makeGRM` | absent | **present** | new in 2.3.0, undocumented on the site |
| `--lmm` | absent | **present** | new in 2.3.0; mentioned only in genemapping.shtml prose |
| `--gdt` | absent | **present** | new in 2.3.1 (Chen et al. 2009) |
| `--noscreen` | absent | **present** | new in 2.3.2 — "for potentially more precise inference by skipping screening relatives with a subset of SNPs" `[WEB]` history |
| `--minConc [0.80]` | absent | **present** | new in 2.3.1 — "to set the heterozygote concordance rate" for `--duplicate` |
| `--plink` | absent | **present** | output format flag; also in the 2.2.7 menu on the ancestry page |
| `--phefile`, `--covfile`, `--prunedsnp` | absent | **present** | optional inputs |
| `--pngplot` | absent from index | present | documented on ancestry page + manual |
| `--rpath` | present | present | |

`[WEB]` The ancestry page preserves the **2.2.7** menu, useful for dating: it shows
`--homog`, `--mtscore`, `--projection [1]`, `--pca [ON]`, `--rplot [ON]`, `--prefix [ex]` — i.e.
boolean flags print `[ON]` when enabled and `--projection` printed its integer argument.

### 3.4 Documented defaults and semantics

| Option | Arg | Default | Semantics `[WEB]` |
|---|---|---|---|
| `--prefix` | string | `king` | "specifies the name of the output files that store various inference results" |
| `--cpus` | int | **half of the total number of (logical) cores** | stated identically on manual, QC, genemapping pages |
| `--sexchr` | int | **23** | "the pair number of the sex chromosome … useful for non-human species" |
| `--pcs` | int | **10** | flagindex says "e.g., 10 as default"; kingpopulation.shtml says "The default pcs is 10" but *also* says MDS saves "Top principal components / ancestry coordinates (**20 by default**)" and the `kingpc.txt` example has **PC1..PC20**. ⚠️ contradiction — see §9. |
| `--degree` | int | none (unfiltered) | 1st/2nd/3rd/…; filters output |
| `--seglength` | float (Mb) | undocumented | "minimum length of IBD segments that are considered towards the relationship inference" |
| `--callrateN` | float | **0.95** | sample-level call rate for `--autoQC` |
| `--callrateM` | float | **0.95** | SNP-level call rate for `--autoQC` |
| `--minConc` | float | **0.80** `[BIN]` | heterozygote concordance cutoff for `--duplicate` |
| `--projection` | int N *or* flag | `1` when bare `[BIN 2.2.7]` | dual meaning, see below |
| `--maxP` | float | none | max p-value for inclusion in output |
| `--prevalence` | float | none | e.g. `0.004` |
| `--rpath` | path | none | "full path of the R program in case 'R' command without a full path cannot properly run" |

**`--degree` → kinship cutoff mapping** `[WEB]` manual.shtml, stated explicitly:
> `--related --degree 2` … Specifically all pairs across families with a kinship coefficient less
> than **0.0884** will be excluded from the output.
and
> In this example [`--kinship --degree 2`], only pairs with kinship coefficient > **0.0884** are
> saved in the king.kin0 output file.
and for `--ibdseg --degree 3`:
> only pairs with IBD proportion > **0.0884** will be saved in the output

⚠️ Note the site uses `0.0884` for *both* `--degree 2` (kinship) and `--degree 3` (PropIBD). `[INF]`
This is a documentation slip; the real ladder is the halving series
`2^-(d+1.5)`: d=1 → 0.17678, d=2 → 0.08839, d=3 → 0.04419. The binary's own R code `[BIN]`
uses exactly `0.35355 / 0.17678 / 0.08839 / 0.04419`, and one plot line uses the full
`0.08838835`. Use `1/(2^(d+1.5))` computed in double precision, not the rounded literals.

**`--projection` has two distinct meanings** `[WEB]`:
1. With `--pca` / `--mds`: *no argument* — "project affected samples to the reference samples' PC
   space". Reference samples are marked `1` and study samples `2` in the **6th (phenotype) column**
   of the `.fam` `[WEB]` history 1.4.2. With multiple `-b` datasets the first dataset is the
   reference `[WEB]` ancestry page ("2409 1000 Genomes samples are detected and used as reference").
2. With `--kinship` / `--ibdseg`: *takes an integer N* — "estimates the kinship coefficients between
   any two samples each from a different subset where the first subset includes the first N
   samples". Abbreviatable as `--proj`. Available in KING 2.2.2+.
   > the kinship estimates from the --kinship --projection N inference should be idential (with no
   > numerical differences) to the standard --kinship inference without splitting

### 3.5 Every example command line published on the site `[WEB]`

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
  prompt> king -b subset1,subset3 --kinship --proj 100000 --prefix subset13
  prompt> king -b subset1,subset4 --kinship --proj 100000 --prefix subset14
  prompt> king -b subset2,subset3 --kinship --proj 100000 --prefix subset23
  prompt> king -b subset2,subset4 --kinship --proj 100000 --prefix subset24
  prompt> king -b subset3,subset4 --kinship --proj 100000 --prefix subset34
  prompt> king -b subset1.bed --kinship --prefix subset1
  prompt> king -b subset2.bed --kinship --prefix subset2
  prompt> king -b subset3.bed --kinship --prefix subset3
  prompt> king -b subset4.bed --kinship --prefix subset4
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
  prompt> king -b ex.bed --prefix ex --rplot
  prompt> king -b ex.bed --prefix ex --build --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --ibdseg --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --related --degree 2 --rplot
  prompt> king -b ex.bed --prefix ex --cluster --degree 2 --rplot
  prompt> king -b ex.bed,ex.bed --fam ex.fam,ex2.fam --duplicate --rplot
  prompt> king -b ex.bed --roh --rplot
  prompt> king -b ex.bed --tdt
  prompt> king -b ex.bed --cov, --mtscore --maxP 5E-8 --invnorm
  prompt> king -b ex.bed --risk --model model.txt --prevalence 0.004 --noflip
```
(`--cov,` in the mtscore line is a typo in the source page — reproduced verbatim.)

---

## 4. OUTPUT FILE INVENTORY

Consolidated. `{p}` = `--prefix` value, default `king`.

| File | Produced by | Delim | Documented? |
|---|---|---|---|
| `{p}.kin` | `--kinship`, `--related` (within-family) | TAB | ✅ full columns |
| `{p}.kin0` | `--kinship`, `--related` (between-family) | TAB | ✅ columns for `--kinship` form |
| `{p}.seg` | `--ibdseg` | TAB | ✅ full columns |
| `{p}.segments.gz` | `--ibdseg` | TAB, gzip | ✅ full columns |
| `{p}.roh` | `--roh` | TAB | ✅ full columns |
| `{p}.rohseg.gz` | `--roh` | TAB, gzip | ✅ full columns |
| `{p}bySNP.txt` | `--bySNP` | TAB `[INF]` | ✅ columns, no example rows |
| `{p}bySample.txt` | `--bysample` | TAB `[INF]` | ✅ columns, no example rows |
| `{p}pc.txt` | `--mds`, `--pca` | **SPACE** | ✅ with example |
| `{p}_InferredAncestry.txt` | `--pca --projection --rplot` | TAB | ✅ with example |
| `{p}grs.txt` | `--risk` | TAB `[INF]` | ✅ columns |
| `{p}updateids.txt` | `--build` | — | ⚠️ name only |
| `{p}updateparents.txt` | `--build` | — | ⚠️ name only |
| `{p}.con` `[BIN]` | `--duplicate` | TAB `[INF]` | ❌ undocumented |
| `{p}.ibs`, `{p}.ibs0` `[BIN]` | `--ibs` | TAB `[INF]` | ❌ undocumented |
| `{p}unrelated.txt` `[INF]` | `--unrelated` | — | ❌ undocumented |
| `{p}_grm.txt` `[BIN]` | `--makeGRM` | — | ❌ undocumented |
| `{p}_eigenvalue.txt` `[BIN]` | `--pca`/`--mds` | — | ❌ undocumented |
| `{p}_relatives.txt`, `{p}_popref.txt`, `{p}_Dist.txt` `[BIN]` | ancestry/cluster helpers | — | ❌ undocumented |
| `{p}_relplot.pdf`, `{p}_ancestryplot.pdf` | `--rplot` | PDF | ✅ mentioned |
| `{p}_relplot.R`, `{p}_buildplot.R`, `{p}_uniqfamplot.R`, `{p}_clusterplot.R`, `{p}_duplicateplot.R`, `{p}_rohplot.R`, `{p}_pcplot.R`, `{p}_ancestryplot.R`, `{p}_pedplot.R`, `{p}_MIerrorplot.R` | `--rplot` | R source | partially |
| `{p}_ibdseg_rplots.tar.gz` | external `king_segments_plot.R` | tar.gz | ✅ |

`[BIN]` Additional dotted extensions present in the string table (analyses outside our scope):
`.anc .aucmap .con .cov .dat .dis .dst .her .homomap .ibdgdt .ibdmap .ibs .ibs0 .ih2 .king .map
.mthomo .npl .popibd .poproh .por .rohdiff` and `_gdt.txt _linear.txt _lmmpc.txt _novclmm.txt
_poodt.txt`.

---

## 5. PER-ANALYSIS SPECIFICATIONS

### 5.1 `--kinship` — KING-robust pairwise kinship

`[WEB]` "estimates pair-wise kinship coefficients using the KING-Robust algorithm described in the
original KING paper." Accurate **up to 2nd-degree**. Robust to population structure.

> If pedigrees are documented in the .fam file … kinship coefficients can be estimated within
> families. **Note if each FID is unique and no pedigrees are provided, then the within-family
> inference will be skipped.** … If the datasets only consist of unrelated individuals as reported,
> then all results are saved in the between-family output.
> Note an unrelated individual is treated as a family of size one.

#### Output `{p}.kin` (within-family) — header, verbatim, then 9 example rows verbatim

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

Columns in order, with the site's verbatim definitions:
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
(`Porportion` is the site's typo — do not replicate in our docs, but it confirms the source text.)

Formats: `Z0`=`%.3f`, `Phi`=`%.4f`, `HetHet`=`%.3f`, `IBS0`=`%.4f`, `Kinship`=`%.4f`, `Error` integer-ish.

#### Output `{p}.kin0` (between-family) — header + 15 example rows verbatim

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
Column order: `FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship` — **no `Z0`, `Phi`, or `Error`**
(there is no reported relationship to compare against across families).

⚠️ `[BIN]` The binary's string table also contains `IID1`, `IID2`, `N_IBS0`, `KinshipX` — so some
`.kin0` variants (X-chromosome, `--related`) use `IID1/IID2` rather than `ID1/ID2`. The website only
ever shows `ID1`/`ID2`. Verify empirically before fixing our header.

#### Interpretation ladder `[WEB]` verbatim — the canonical KING kinship cutoffs

> an estimated kinship coefficient range **>0.354, [0.177, 0.354], [0.0884, 0.177] and
> [0.0442, 0.0884]** corresponds to **duplicate/MZ twin, 1st-degree, 2nd-degree, and 3rd-degree**
> relationships respectively.

> A negative kinship coefficient estimation indicates an unrelated relationship. The reason that a
> negative kinship coefficient is not set to zero is a very negative value may indicate the
> population structure between the two individuals.

`[INF]` These are rounded `2^-(k+1.5)`: 0.35355, 0.17678, 0.08839, 0.04419.

#### Scaling
- `--kinship --degree` tested to **1,000,000 samples**.
- `--kinship --projection N` — see §3.4; identical numeric results to unsplit run.
- Companion plots: `hapmapkin.pdf` / `hapmapkin.R` (within), `hapmapkin0.pdf` / `hapmapkin0.R`
  (between) — linked from manual.shtml.

---

### 5.2 `--ibdseg` — IBD segment inference

`[WEB]` New in 2.1. Accurate **up to 3rd–4th degree** (array vs WGS). "as fast as estimating kinship
coefficients, e.g., seconds in 1000s of samples". Tested to 1,000,000 samples with `--degree`.
Applies to chromosome X as well since 2.2.3.

#### Output `{p}.seg` — header + 9 example rows verbatim

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

Column definitions verbatim:
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

**`InfType` value set** `[WEB]`: `Dup/MZTwin`, `PO`, `FS`, `2nd`, `3rd`, `4th`, `UN`.
`[BIN]` The R plotting code's `Inf.type` vector has 5 named levels plus a fallback `6`, and adds an
"other" bucket `dO` — so a 6th label exists in practice (likely `Other`/`Unknown`). Note `.seg`
example shows `PropIBD = IBD2Seg + IBD1Seg/2` exactly: `0.9976/2 = 0.4988` ✓, `0.9969/2 = 0.49845`
→ printed `0.4985` (round-half-up or round-half-even on the 5) ⚠️ and `0.9987/2 = 0.49935` →
`0.4994`, `0.9999/2 = 0.49995` → `0.5000`. All three ties round **up** ⇒ `%.4f` with round-half-away
-from-zero, i.e. plain C `printf` behaviour on a value slightly above the tie. Compute in double and
let `printf("%.4f")` do it.

#### Output `{p}.segments.gz` — gzip; header + 9 example rows verbatim

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

Column definitions verbatim:
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

`IBDType` ∈ {`IBD1`, `IBD2`}. `StartMB`/`StopMB` = `%.3f`; `Length` = `%.1f`.
⚠️ `Length` ≠ `StopMB − StartMB` exactly: row 1 gives `95.862 − 51.799 = 44.063` → printed `44.1`;
row 3 `88.714 − 0.143 = 88.571` → `88.6`; row 5 `90.221 − 0.080 = 90.141` → `90.1`. Consistent with
`%.1f` of the exact difference computed at full precision. Row 7: `180.626 − 70.869 = 109.757` →
`109.8` ✓. So `Length` **is** `StopMB − StartMB` at `%.1f`.
Ordering: grouped by pair, then by `IBDType`, then ascending `Chr`, then ascending `StartMB`.

The **verbatim** header line was also verified against the raw HTML at `raw/manual.shtml` line 560
region and is tab-separated in the real file (see §2.3).

#### Segment plotting helper `[WEB]`
```
  prompt> Rscript king_segments_plot.R ex ibdseg
```
Requires `ggplot2` and `parallel`. Output: `ex_ibdseg_rplots.tar.gz` (one plot file per close pair).

---

### 5.3 `--related` — integrated fast inference

`[WEB]` "integrative, fast, and accurate inference for close relationships … highly recommended,
especially when dealing with biobank-level datasets." Tested to **~10 million samples**
(~5×10¹³ pairs). Without `--degree` it identifies **1st-degree** relatives; `--degree 2` for 2nd.
> Although distant relatedness that is higher than 2 is allowed, no fast algorithm is available at
> the moment and computation is substantially slower than --related --degree 2.

Writes **`{p}.kin` and `{p}.kin0`** — same names as `--kinship` but with a **wider schema**.

#### Output `{p}.kin` from `--related` — header + 9 example rows verbatim

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

**16 columns, in order:**
`FID ID1 ID2 N_SNP Z0 Phi HetHet IBS0 HetConc HomIBS0 Kinship IBD1Seg IBD2Seg PropIBD InfType Error`

Definitions verbatim (note the site's list **omits `HetConc` and `HomIBS0`** even though they are in
the header — a documentation gap):
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
`[INF]` Missing definitions, from the plot axis labels `[BIN]`
(`xlab="Heterozygote Concordance Rate"`):
- `HetConc` — heterozygote concordance rate = Pr(both het | at least one het); `%.4f`
- `HomIBS0` — IBS0 proportion restricted to homozygote-informative pairs; `%.4f`

⚠️ **`HetHet` is `%.4f` here** (`0.2324`) but `%.3f` in the `--kinship`-only `.kin` (`0.162`).
Two different writers. Confirm empirically.

`Error` semantics differ between the two `.kin` flavors:
- `--kinship`: "differences between the **estimated and specified kinship coefficients**"
- `--related`: "differences between **inferred and reported relationship**"

#### `--related --degree 2` filter `[WEB]`
> only related pairs (up to the 2nd-degree in this case) between families are included in the
> output. Specifically all pairs across families with a kinship coefficient less than 0.0884 will be
> excluded from the output.

`[BIN]` The runtime message is `Between-family relatives (kinship >= %.5lf) saved in file %s`
→ the cutoff is printed at **5 decimals** (e.g. `0.08839`), and the comparison is **`>=`**, not `>`.

#### Plots
`--related --rplot` → `{p}_relplot.pdf` (example: `ex_relplot.pdf`) from `{p}.kin` + `{p}.kin0`.
Since 2.2: pedigree errors visualized (documented pedigree left, inferred relatedness right);
cryptic relatedness summarized as unique family configurations with counts. Requires R `igraph`
+ `kinship2`.

#### `[BIN]` Inference decision rules — from the *R plotting code* embedded in the binary

These are literal R expressions in the emitted plot script (a published-artifact string, not C++
source), and they encode the `InfType` classification boundaries exactly:

```r
d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)
d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08
d2    <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)
d3    <- data$PropIBD>0.08839 & data$PropIBD<=0.17678
d4    <- data$PropIBD>0.04419 & data$PropIBD<=0.08839
dU    <- data$PropIBD<=0.04419
allpair <- data$PropIBD>0 | data$Kinship>0.04419
```
Plot colors, which reveal the category ordering: `d0`=purple (Dup/MZ), `d1.PO`=red, `d1.FS`=green,
`d2`=blue, `d3`=magenta, `d4`=gold, `dO`(other)=gold, `dU`=black.

`[BIN]` Two related runtime messages that expose tunable cutoffs:
```
1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf
Cutoff value for IBS0 between FS and PO is set at %.4f
```

---

### 5.4 `--duplicate` — duplicates / MZ twins

`[WEB]` "the fastest (and accurate) algorithm to identify duplicates or MZ twins. The running time
is in seconds, unless the number of samples is > 1,000,000 in which case a few minutes may be
needed." Tested to ~10 million samples. Primary use: cross-study duplicate detection via multiple
`-b` datasets.

⚠️ **The website never documents `--duplicate`'s output file name or columns.** From `[BIN]`:
- Output extension `.con` ⇒ `{p}.con` `[INF]`
- Runtime messages:
  ```
  %d pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s
  %lli pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s
  ```
  (two variants for `int` vs `long long` pair counts; the concordance threshold is printed as an
  **integer percent**, i.e. `--minConc 0.80` prints as `80`)
- Header field strings present: `FID1 ID1 FID2 ID2 N_SNP N_IBS0 Concord`, plus `N_Het1 N_Het2
  N_IBS1 N_IBS2 HetConc` for the `--ibs` family. `[INF]` Likely `.con` header is
  `FID1 ID1 FID2 ID2 N_SNP N_IBS0 Concord`. **Must be confirmed by running the binary.**
- `--minConc` default `0.80` `[BIN]`; new in 2.3.1 `[WEB]`.
- `--duplicate --rplot` → `{p}_duplicateplot.R` (+ PDF if `igraph` installed); directed igraph
  showing sample-swap / shift patterns `[WEB]` KINGvisualization.

---

### 5.5 `--ibs` — IBS summary statistics

`[WEB]` one sentence only:
> --ibs provides summary statistics such as the counts of IBS0, IBS1, IBS2, the average of IBS, in
> additional to the kinship estimates.

⚠️ No file name, no columns documented. `[BIN]`: extensions `.ibs` and `.ibs0` exist; field-name
strings `N_IBS0 N_IBS1 N_IBS2 N_Het1 N_Het2 HetConc Concord` exist; runtime message
`Between-family IBS data saved in file %s`. ⇒ `{p}.ibs` (within-family) and `{p}.ibs0`
(between-family), mirroring the `.kin`/`.kin0` split `[INF]`. **Confirm by running.**

---

### 5.6 `--homog` — homogeneous-population kinship (RETIRED)

`[WEB]` manual.shtml still documents it:
> --homog estimates pair-wise kinship coefficients assuming a homogeneous population. The best
> application of --homog may be for the linear mixed models (LMM) … Although --homog is not
> recommended as a good method to infer relatedness in general populations, it provides inference
> results comparable to multiple alternative methods.

`[WEB]` history.shtml 2.3.0: **"--homog is retired since it does not do a good job on either
relationship inference or GRM"**. `[BIN]` Absent from the 2.3.2 menu. ⇒ **Do not implement.**
Its role was taken over by `--makeGRM`.

---

### 5.7 `--unrelated` — extract a maximal unrelated subset

`[WEB]`
```
  prompt> king -b ex.bed --unrelated --degree 2
```
> This example estimates relatedness in the data first, followed by extracting a list of individuals
> that contains no pairs of individuals with a 1st- or 2nd-degree relationship.

Algorithm reference `[WEB]`: **Manichaikul et al. 2012**, PDF at
`https://www.chen.kingrelatedness.com/publications/pdf/PLoS8e1002640.pdf` (PLoS Genetics 8:e1002640).
`--degree 3` and higher allowed since 2.2.9 `[WEB]`.

⚠️ Output file name/columns undocumented. `[BIN]` runtime messages:
```
A list of %d unrelated individuals saved in file %s
An alternative list of %d to-be-removed individuals saved in file %s
```
⇒ **two** files: the keep-list and the remove-list. `[INF]` likely `{p}unrelated.txt` and
`{p}unrelated_toberemoved.txt`; format almost certainly two columns `FID IID` (the string
`FID IID` exists standalone in the binary). **Confirm by running.**

---

### 5.8 `--build` — pedigree reconstruction

`[WEB]`
> --build reconstructs pedigrees using SNP data without the need of specifying pedigrees (although
> the pedigree information can still be incorporated)
> The output includes two files: **kingupdateids.txt** and **kingupdateparents.txt**.

Consumed by PLINK:
```
  prompt> plink1.9 --bfile ex --update-ids kingupdateids.txt --make-bed --out ex2
  prompt> plink1.9 --bfile ex2 --update-parents kingupdateparents.txt --make-bed --out ex3
```
`[INF]` Column layout is therefore fixed by PLINK 1.9's contract, not by KING:
- `--update-ids`: `oldFID oldIID newFID newIID` (4 cols)
- `--update-parents`: `FID IID newPatID newMatID` (3–4 cols: FID IID father mother)

> The current --build algorithm connects all 1st-degree relatives with high accuracy. Known
> scenarios that --build does well are families that consist of at least a pair of full siblings,
> and/or a parent-child trio.

2.1.8: "--build has improved and KING pedigree reconstruction can now incorporate 2nd-degree
relatedness inference." 2.1.6: works for N<100. 2.2.5 / 2.2.1: further fixes.
`--build --rplot` → `{p}_buildplot.R` + pedigree PDF (needs R `kinship2`).

---

### 5.9 `--cluster` — cluster relatives into families

`[WEB]`
> --cluster is both a standalone parameter and a parameter to go with other options. As a standalone
> option, it clusters relatives into families by generating an **updateid file** which can then be
> used to update the pedigrees (e.g., using PLINK --update-ids). --cluster can also be used to group
> cyptic relatives together prior to association analysis, e.g.,
```
  prompt> king -b ex.bed --cluster --tdt
```
Also composes with QC: `--cluster --bySNP`, `--cluster --bysample` → "family-based QC without using
any reported pedigrees".
`[BIN]` field name `KING_FID` / `S KING_FID` ⇒ the cluster output introduces a synthetic family ID
column named **`KING_FID`**. `--degree 3`+ allowed since 2.2.9.
`--cluster --rplot` → `{p}_clusterplot.R` + igraph PDF per clustered family.

---

### 5.10 `--bySNP` — SNP-level QC

`[WEB]` Output: **`{p}bySNP.txt`** (default `kingbySNP.txt`). Columns, verbatim and in order:
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

⚠️ **The binary's actual header differs from the docs.** `[BIN]` exact header fragments:
```
 Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate
 N_MZ N_HetMZ N_errMZ Err_InMZ Err_InHetMZ
 N_PO N_HomPO N_errPO Err_InPO Err_InHomPO
 N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio
```
The MZ group in 2.3.2 has **5** fields (`N_MZ N_HetMZ N_errMZ Err_InMZ Err_InHetMZ`) where the docs
list **3** (`N_MZ N_errMZ Err_InMZ`) — `N_HetMZ` and `Err_InHetMZ` were added after Feb 2018.
The PO and trio groups match the docs exactly. Trust the binary.
Note the header fragments begin with a **leading space** in the string table, suggesting the header
is assembled with a leading separator per group.
No example rows are published. `[WEB]` A bug in `--bySNP` and `--bysample` "for the scenario of
presence of families" was fixed in **2.3.1** — so pre-2.3.1 output is not a valid parity target.

---

### 5.11 `--bysample` — sample-level QC

`[WEB]` Output: **`{p}bySample.txt`** (default `kingbySample.txt` — capital `S`). Columns verbatim:
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
`[BIN]` corroborating strings: `FID IID FA MO SEX N_SNP Missing Heterozygosity` (one contiguous
literal — confirms these 8 are the leading fields and are **space-or-tab joined in that exact
order**) and a standalone ` MI_Removal` (leading space) confirming it is the final column.
No example rows published. `--bySample` crash fixed in 2.2.9; families bug fixed in 2.3.1 `[WEB]`.

---

### 5.12 `--autoQC` — automated QC pipeline

`[WEB]`
> --autoQC option performs a straightforward QC pipeline, including sample-level QC (at call rate
> **95% by default**, or a different call rate set by --callrateN), SNP-level QC (at call rate
> **95% by default**, or a different call rate set by --callrateM), and **gender QC**. This analysis
> generates a **list of SNPs to be removed, and a list of samples to be removed.**
```
  prompt> king -b ex.bed --autoQC
```
⚠️ File names undocumented. `[BIN]` string `CallRateLessThan%d` suggests generated
names/labels embed the integer percent (e.g. `CallRateLessThan95`). **Confirm by running.**

---

### 5.13 `--roh` — runs of homozygosity

`[WEB]` Two outputs: **`{p}.roh`** (per-sample inbreeding) and **`{p}.rohseg.gz`** (segments).
Applies to chromosome X since 2.2.3.

#### `{p}.roh` — header + 9 rows verbatim
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
Column order: `FID ID FA MO SEX MaxROH FInbred` — **note `ID`, not `IID`**, unlike
`{p}bySample.txt` which uses `IID`. Missing parents are `0`. `MaxROH`=`%.1f` (Mb),
`FInbred`=`%.4f`. No per-column prose definitions are given on the site; `MaxROH` is the longest
single ROH in Mb and `FInbred` the genomic inbreeding coefficient `[INF]`.

#### `{p}.rohseg.gz` — header + 8 rows verbatim
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
Same tail schema as `{p}.segments.gz` minus the pair/IBDType columns:
`FID ID Chr StartMB StopMB StartSNP StopSNP N_SNP Length`.
`Length` = `StopMB − StartMB` at `%.1f` (`97.455−70.869=26.586`→`26.6` ✓;
`167.849−136.510=31.339`→`31.3` ✓; `31.787−25.472=6.315`→`6.3` ✓).
Consistency check: `NA12342`'s `MaxROH` in `{p}.roh` is `31.3`, matching its longest segment ✓.
Sort order: by `FID`, then `ID`, then `Chr` ascending, then `StartMB` ascending.

`--roh --rplot` `[WEB]` KINGvisualization: plots ROH for "all individuals with proportion of their
genomes being ROH > **4.4%**, which corresponds to being offspring of parents that are 2nd-degree or
closer". Requires R `ggplot2`. Emits `{p}_rohplot.R`. `[BIN]` extension `.poproh`, `.rohdiff` also
exist.

---

### 5.14 `--mds` / `--pca` — population structure

`[WEB]` kingpopulation.shtml.
> The Multidimensional Scaling (MDS) with the Euclidean distance is **highly recommended** for the
> identification of population substructure.
> Principal Component Analysis (PCA) … **Please run LD-pruning prior to PCA analysis.**

Warning verbatim:
> Precompile KING binaries with versions lower than 2.2.3 are not suitable for population structure
> analysis in larger datasets for lacking LAPACK libraries.

#### Output `{p}pc.txt` — SPACE-delimited; header + 9 rows verbatim
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
Fixed leading 6 columns `FID IID FA MO SEX AFF`, then `PC1..PCk`. `[WEB]` "The top 10 principal
components / ancestry coordinates are in the 7th to the 16th columns." PC values `%.4f`.
`[BIN]` `FID IID FA MO SEX AFF` exists as one contiguous string literal ⇒ header prefix confirmed.
`--pca` "has the same format as kingpc.txt from the --mds analysis" `[WEB]`.
`[BIN]` companion `{p}_eigenvalue.txt` and `{p}_pcplot.R`; runtime message
`%d principal components saved in file %s`.

`[BIN]` The 2.2.7 ancestry screen printout shows eigenvalues echoed to stdout:
```
Largest 10 eigenvalues: 2059.73 1366.88 700.55 618.08 276.33 265.47 238.72 217.65 198.71 189.81
```
⇒ `%.2f`, space-separated, count = `--pcs`.

---

### 5.15 Ancestry inference (`--pca --projection --rplot`)

`[WEB]` ancestry page. Requires R + the **`e1071`** package (SVM). Superpopulations inferred:
**AFR, AMR, EAS, EUR, SAS**. Reference = 1000 Genomes (`KGref.*`), passed as the **first** `-b`
dataset.

Full screen printout, verbatim (KING 2.2.7 — useful as a template for our own stdout):
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
Timestamps use C `ctime()` format (`Thu May 27 18:29:56 2021`). Note `-bname` echoes the prefixes
with `.bed` **stripped**.

#### Output `{p}_InferredAncestry.txt` — header + 9 rows verbatim
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
9 columns: `FID IID PC1 PC2 Anc_1st Pr_1st Anc_2nd Pr_2nd Ancestry`.
`[BIN]` exactly corroborated by the R literal:
```r
colnames(pred.out) <- c("FID", "IID", "PC1", "PC2", "Anc_1st", "Pr_1st", "Anc_2nd",  "Pr_2nd", "Ancestry")
```
(note the double space before `"Pr_2nd"` in the source literal — cosmetic only).
Numeric columns are `%g`-style, trailing zeros stripped (see §2.3 warning).
`[BIN]` also: `FID IID Population` and `Format: FID IID Population` ⇒ an optional user-supplied
population-label input file with that 3-column format; and
`FID IID HetProj HetRef MinDist Kinship Closest RefID` ⇒ a projection-diagnostics output.
Alternative-reference workflow: stand-alone R at
`https://github.com/chenlab-uva/AncestryInference_KING` `[WEB]`.

---

### 5.16 `--makeGRM` `[BIN]` only

New in 2.3.0, replacing `--homog`. **Completely undocumented on the website.**
`[BIN]` runtime message `GRM saved in file %s`; file suffix `_grm.txt`. ⇒ `{p}_grm.txt` `[INF]`.

---

### 5.17 Association: `--tdt`, `--gdt`, `--lmm`, `--mtscore`

`[WEB]` genemapping.shtml (thin).
- `--tdt` — "the well-known Transmission/Disequilibrium Test for family data that consist of
  parent-affected child trios."
- `--lmm` — linear mixed model, SNP × quantitative trait; "quite efficient, especially for a lot of
  traits, e.g., in eQTL/pQTL/meQTL/mQTL analysis". New in 2.3.0.
- `--gdt` — new in 2.3.1, family-based association, **Chen et al. 2009**.
- `--mtscore` — many-traits score test; **gone in 2.3.2**.
- Modifiers: `--trait`, `--covariate`, `--invnorm`, `--maxP`, `--prefix`, `--cpus`.
```
  prompt> king -b ex.bed --tdt    
  prompt> king -b ex.bed --cov, --mtscore --maxP 5E-8 --invnorm    
```
⚠️ No output columns documented at all. `[BIN]`: `{p}_gdt.txt`, `{p}_linear.txt`, `{p}_lmmpc.txt`,
`{p}_novclmm.txt`, `{p}_poodt.txt`; messages `LMM scan results (lambda_GC=%.4lf) saved in file %s`,
`GDT GWAS summary statistics saved in file %s`, `%s GWAS summary statistics saved in file %s`,
`IDs of %d TDT uninformative pedigrees saved in file %s.`
**Out of scope for a relatedness reimplementation.**

---

### 5.18 `--risk` — genetic risk scores

`[WEB]` riskprediction.shtml.

#### Input model file (`--model model.txt`) — verbatim example
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
```
SNP: SNP name
EA: effect allele
AF: allele frequency of the effect allele
WT: weight at the effect allele
CHR: chromosome of the SNP
POS: position of the SNP
OA: other allele
```
(Tab-delimited; `SNP` header occupies 2 tab stops.)

#### Command
```
  prompt> king -b ex.bed --risk --model model.txt --prevalence 0.004 --noflip
```
`--prevalence` and `--noflip` are optional. `--prevalence` only affects PPV/NPV. `--noflip` is for
when strands already match the model.

#### Screen output — verbatim
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
Rates `%.4f`; risk cutoffs `0.1`–`0.9` step `0.1`; AUC `%.4f`.

#### Output `{p}grs.txt` — columns verbatim
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
Second reference `[WEB]`: Onengut-Gumuscu S, … Rich SS (2019) *Diabetes Care* 42(3):406-415.
`[BIN]` `.aucmap`, `{p}_aucmapplot.R` relate to AUC mapping.

---

## 6. VISUALIZATION (`--rplot` / `--pngplot`) `[WEB]` KINGvisualization.shtml

`--rplot` "generates R code first and then calls R program to make plots in a PDF file";
`--pngplot` the PNG variant. If the required R package is missing, **only the `.R` file is written**
— a useful degradation contract.

| Command | R script emitted | R packages required |
|---|---|---|
| `king -b ex.bed --prefix ex --rplot` (no analysis) | `ex_pedplot.R` `[BIN]` | `kinship2` |
| `--build --degree 2 --rplot` | `ex_buildplot.R` | `kinship2` |
| `--ibdseg --degree 2 --rplot` | `ex_uniqfamplot.R` | `igraph` |
| `--related --degree 2 --rplot` | `ex_uniqfamplot.R` (+ pedigree-error plots) | `igraph` **and** `kinship2` |
| `--cluster --degree 2 --rplot` | `ex_clusterplot.R` | `igraph` |
| `--duplicate --rplot` | `king_duplicateplot.R` | `igraph` |
| `--roh --rplot` | `{p}_rohplot.R` | `ggplot2` |
| `--pca --projection --rplot` | `{p}_ancestryplot.R` | `e1071` |

Differences between `--ibdseg --rplot` and `--related --rplot`, verbatim:
> 1) --related --rplot only visualizes unique family configurations that are cryptic (between
> families); 2) --related is expected to be orders of magnitude faster; and 3) --related --rplot
> also visualizes (within-family) pedigree errors.

Toy-data recipe for the duplicate plot (verbatim):
```
  prompt> head -232 ex.fam > ex2.fam
  prompt> awk 'NR>232 && NR%2==1' ex.fam >> ex2.fam
  prompt> awk 'NR>232 && NR%2==0' ex.fam >> ex2.fam
  prompt> king -b ex.bed,ex.bed --fam ex.fam,ex2.fam --duplicate --rplot
```

`[BIN]` R plot titles worth mirroring:
```
main = "Kinship vs IBS0 in %s Families"
main = "Kinship vs Heterozygote Concordance In %s Families"
main = paste("Kinship vs Proportion IBD (Corr=", round(cor(...),digit=3),") in %s Families",sep="")
xlab="Heterozygote Concordance Rate", ylab = "Estimated Kinship Coefficient"
```

---

## 7. VERSION HISTORY — entries that affect output parity `[WEB]` history.shtml

Reproduced only where behavior/format changed. Full text in `research/txt/history.txt`.

| Ver | Date | Parity-relevant change |
|---|---|---|
| **2.3.2** | Sept 8, 2023 | (1) crash on `--related --degree 3+` (2.3.1-only bug) fixed; (2) `--noscreen` added for `--related`/`--duplicate` — "skipping screening relatives with a subset of SNPs" |
| **2.3.1** | July 28, 2023 | `--gdt` added; **bug in `--bySNP`/`--bysample` fixed for presence of families**; `--minConc` added for `--duplicate` |
| **2.3.0** | Oct 10, 2022 | `--lmm` added; **`--makeGRM` added**; **`--homog` retired** |
| **2.2.9** | Sept 20, 2022 | `--cluster`/`--unrelated` now allow `--degree 3`+; `--bySample` crash fixed; `--risk` bug fixed |
| **2.2.8** | May 10, 2022 | `--pca --projection` minor fix |
| **2.2.7** | May 18, 2021 | `--ibdseg` bug from 2.2.5 fully fixed |
| **2.2.6** | Mar 23, 2021 | `--ibdseg` hang fixed; `--pca --projection --rplot` ancestry inference |
| **2.2.5** | June 5, 2020 | `--unrelated` bug fixed; **`--ibdseg` substantially improved**; `--build` improved |
| **2.2.4** | Oct 11, 2019 | population structure much faster; **`--pcs` added, default 10** |
| **2.2.3** | Aug 9, 2019 | static LAPACK; **IBD segment + ROH now apply to chromosome X**; `--related` bug when #SNPs < 4096 fixed |
| **2.2.2** | May 29, 2019 | multi-dataset identical-SNP-set bug fixed; **`--projection N` introduced** |
| **2.2.1** | May 14, 2019 | `--rplot` improved; `--ibdseg`/`--roh` maxIBD1/maxIBD2 minor bug ("not affecting the main inference"); `--build` bug fixed |
| **2.2** | Mar 28, 2019 | family visualization; duplicate IDs allowed between first (reference) dataset and others; `--ibdseg --degree 2` bug (2.1.8-only) fixed |
| **2.1.8** | Feb 27, 2019 | `--build` incorporates 2nd-degree; pedigree drawing via `kinship2` |
| **2.1.6** | Nov 28, 2018 | `--lessmem` retired; sample-size cap removed (N>10M); multi-dataset input; % progress printed |
| **2.1.5** | Aug 24, 2018 | **`--sexchr` added, default 23**; `--risk` fixed; `--ibdseg` up to 700,000 samples (was 256,000) |
| **2.1.4** | Jun 6, 2018 | **accuracy statement: `--related`/`--ibdseg` accurate to 4th degree; `--kinship` to 2nd degree** |
| **2.1.3** | Feb 13, 2018 | `--ibdseg`, `--related`, `--roh` improved; `--cluster --bySNP`, `--build` bugs fixed |
| **2.1.2** | Dec 14, 2017 | IBD segment algorithm improved |
| **2.1** | Oct 24, 2017 | **`--ibdseg` introduced**; `--related` integrative |
| **2.0** | Oct 17, 2016 | multi-core; `--tdt`, `--mtscore` added |
| **1.9** | Oct 10, 2015 | fast `--duplicate` / `--related` |
| **1.4.2** | 2013 | `--pca --projection`; **reference samples = phenotype col 6 value 1, study = 2** |
| **1.4** | Dec 14, 2011 | `--unrelated` introduced |

⚠️ **Parity target must be 2.3.2 specifically.** `--ibdseg` changed materially in 2.1.2, 2.1.3,
2.2.5, 2.2.6, 2.2.7; `--bySNP`/`--bysample` changed in 2.3.1. Published example blocks in
manual.shtml date from various eras and may not byte-match a 2.3.2 run.

---

## 8. CONSOLIDATED THRESHOLD CONSTANTS

| Name | Value `[WEB]` | Exact `[INF]` | Use |
|---|---|---|---|
| Dup/MZ | > 0.354 | 2^-1.5 = 0.353553390593… | kinship |
| 1st degree | [0.177, 0.354] | 2^-2.5 = 0.176776695297… | kinship |
| 2nd degree | [0.0884, 0.177] | 2^-3.5 = 0.088388347648… | kinship |
| 3rd degree | [0.0442, 0.0884] | 2^-4.5 = 0.044194173824… | kinship |
| `--degree d` cutoff | — | `2^-(d+1.5)` | output filter, compared with **`>=`** `[BIN]` |
| PO vs FS split | — | `IBD2Seg` 0.08; `IBD1Seg+IBD2Seg` 0.9 / 0.96 `[BIN]` | InfType |
| ROH plot threshold | 4.4% of genome | 0.044194… | `--roh --rplot` |
| `--minConc` | 0.80 `[BIN]` | | `--duplicate` |
| `--callrateN` / `--callrateM` | 0.95 | | `--autoQC` |

---

## 9. CONTRADICTIONS, GAPS, AND OPEN QUESTIONS

Things the website does **not** answer that must be resolved by running the binary:

1. **`--pcs` default: 10 or 20?** flagindex + kingpopulation prose say 10; the same page's prose
   says "20 by default" and its `kingpc.txt` example has PC1..PC20. `[BIN]` message
   `%d principal components saved in file %s` printed `10` in the 2.2.7 log with default settings.
   → Likely **10** since 2.2.4; the 20-column example predates it.
2. **`--duplicate` output** — file name and header entirely undocumented. `[BIN]` points to
   `{p}.con` with `FID1 ID1 FID2 ID2 N_SNP N_IBS0 Concord`. **Must verify.**
3. **`--ibs` output** — undocumented. `[BIN]` points to `{p}.ibs` / `{p}.ibs0`. **Must verify.**
4. **`--unrelated` output** — two files (keep-list + remove-list), names undocumented.
5. **`--autoQC` output** — two lists, names undocumented.
6. **`--makeGRM`** — undocumented entirely; `{p}_grm.txt`.
7. **`ID1/ID2` vs `IID1/IID2`** in `.kin0` — the binary contains both spellings.
8. **`HetHet` precision** — `%.3f` in `--kinship`'s `.kin`, `%.4f` in `--related`'s `.kin`.
9. **`HetConc` / `HomIBS0`** appear in `--related`'s `.kin` header but have **no published
   definition**.
10. **`--bySNP` MZ block** — binary has 5 fields (`N_MZ N_HetMZ N_errMZ Err_InMZ Err_InHetMZ`),
    docs list 3. Docs are pre-2018.
11. **`--seglength` default** — never stated.
12. **Sort order** of `.kin0` / `.seg` rows — inferable from examples (FID1, then ID1, then FID2,
    then ID2, ascending as read from the .fam) but never stated.
13. **Whether a trailing newline / trailing tab** terminates each row — not determinable from HTML.
14. **`--noscreen`** — new in 2.3.2, no documentation beyond the one-line changelog entry.
15. **X-chromosome outputs** — `[BIN]` `X-Chr IBD-sharing inference saved in file %s`,
    `Additional summary statistics of X-Chr IBD segments saved in file %s`, field `KinshipX`.
    Entirely undocumented on the site.
16. **`ex.tar.gz`** (332 samples / 18,290 SNPs) is the only published fixture; the `.kin`/`.kin0`
    examples in manual.shtml use ~2.36M SNPs and FIDs `28/117/1344` — i.e. a **different, unpublished
    dataset**. Only the `--related` `.kin` example (18,250–18,280 SNPs, FIDs `Y001`–`Y003`) and the
    `--roh`/`.seg` examples (`NA*` IDs) come from `ex`. Byte-comparison against published blocks is
    only possible for those.

---

## 10. RECOMMENDED NEXT STEPS FOR PARITY

1. Download `ex.tar.gz` (public, 1.35 MB) and run 2.3.2 for every option; diff against §5 blocks.
2. Capture real output files and `xxd` the first rows to settle tab-vs-space and trailing-newline.
3. Resolve every item in §9 empirically.
4. Treat §5's verbatim blocks as regression fixtures for `--related`, `--ibdseg`, `--roh`.
