# open-king command-line reference

`open-king` reads one PLINK1 fileset and writes relatedness results as plain text. This page is
the complete reference for its command line: every option the parser accepts, what each one
does, which analyses it affects, and the handful of parser behaviours that surprise people.

It is written for someone who has `.bed`/`.bim`/`.fam` files and wants kinship coefficients
or IBD segments out of them — and who may be replacing KING 2.3.2 with this binary, or
running both and diffing. It assumes you know what a kinship coefficient is; it assumes
nothing about this codebase.

Every command below was run against `target/release/open-king` and every output block is pasted
from that run. Reproduce them with the corpus the parity suite ships:

```
cargo build --release
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
```

The examples run from a directory holding those filesets (`multifam`, `dups`, `bigish`,
`sexchr`, `trio`, `missing`) with `target/release` on `PATH`. Timestamps and the CPU count
differ on your machine; nothing else should.

A dozen examples need a **deliberately malformed or relabelled** fileset that the corpus does
not ship — an A1-major recode, an unsorted map, alphabetic chromosome codes. Every one of them
is built by the short script in [§10](#10-the-derived-filesets-used-above), which derives them
all from the corpus with the Python standard library. Their names are used verbatim below.

**Related pages.** [PARITY.md](PARITY.md) is the authoritative statement of what does and
does not match KING 2.3.2 — read it before trusting a number in a cross-check.
[BEHAVIOR.md](BEHAVIOR.md) records the measured reference behaviour behind most rules here.
[SPEC.md](SPEC.md) is the implementation spec. This page does not repeat their evidence; it
links to it.

---

## Contents

* [1. Synopsis](#1-synopsis)
* [2. Quick start](#2-quick-start)
* [3. Input](#3-input)
* [4. Output](#4-output)
* [5. Option reference](#5-option-reference)
* [6. How the parser behaves](#6-how-the-parser-behaves)
* [7. Exit status and fatal errors](#7-exit-status-and-fatal-errors)
* [8. Accepted compatibility spellings outside product scope](#8-accepted-compatibility-spellings-outside-product-scope)
* [9. Differences from the reference](#9-differences-from-the-reference)
* [10. The derived filesets used above](#10-the-derived-filesets-used-above)

---

## 1. Synopsis

```
open-king -b <file>.bed [analysis ...] [parameter ...]
```

`-b` names the input. At least one *analysis* option must be given or nothing is computed.
Parameters modify whatever analyses are running. Options may appear in any order, before or
after `-b`.

There is **no `--help` and no `--version`.** Running `open-king` with no arguments prints the
banner, the full option table, and a fatal error:

```
$ open-king
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

That table is the option list. It is printed on **every** run, with the value of each option
in brackets, which makes it the quickest way to check that a command line was understood the
way you meant it. `[ON]` marks a switch that is on; an option with no bracket is off or
unset. The `--noscreen [-1717986816]` you see there is not a typo — see
[`--noscreen`](#--noscreen-n).

`open-king --help` is not an error you can act on either; it parses as an undefined option and the
run continues:

```
$ open-king -b multifam.bed --kinship --help
WARNING - 
Problems encountered parsing command line:

Command line parameter --help is undefined
```

---

## 2. Quick start

Kinship for every pair in a fileset:

```
$ open-king -b multifam.bed --kinship
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen

The following parameters are in effect:
                   Binary File :    multifam.bed (-bname)

Additional Options
         Close Relative Inference : --related, --duplicate
   Pairwise Relatedness Inference : --kinship [ON], --ibdseg, --ibs, --makeGRM
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

KING starts at Fri Aug 14 22:47:11 2026
Loading genotype data in PLINK binary format...
Read in PLINK fam file multifam.fam...
  PLINK pedigrees loaded: 20 samples
Read in PLINK bim file multifam.bim...
  Genotype data consist of 15000 autosome SNPs
  PLINK maps loaded: 15000 SNPs
Read in PLINK bed file multifam.bed...
0%6%13%19%25%31%38%44%50%56%63%69%75%81%88%  PLINK binary genotypes loaded.
94%  KING format genotype data successfully converted.
Autosome genotypes stored in 235 words for each of 20 individuals.

Options in effect:
	--kinship

Within-family kinship data saved in file king.kin

Relationship summary (total relatives: 36 by pedigree, 36 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	24	12	0	0	4
  Inference	0	24	11	1	0	4

Relationship inference across families starts at Fri Aug 14 22:47:11 2026
16 CPU cores are used.
                                         ends at Fri Aug 14 22:47:11 2026
Between-family kinship data saved in file king.kin0
Note --kinship --degree <n> can filter & speed up the kinship computing.
KING ends at Fri Aug 14 22:47:11 2026
```

Two tab-separated files, both with a header line. `.kin` holds pairs **inside** a family,
`.kin0` pairs **across** families:

```
$ head -3 king.kin
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	Kinship	Error
FAM1	A_C1	A_C2	15000	0.250	0.2500	0.2219	0.0153	0.2721	0
FAM1	A_C1	A_C3	15000	0.250	0.2500	0.2171	0.0214	0.2476	0

$ head -3 king.kin0
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	Kinship
FAM1	A_F	FAM2	B_F	15000	0.2169	0.0204	0.2501
FAM1	A_F	FAM2	B_M	15000	0.1386	0.0685	0.0001
```

Everything goes to **stdout**; nothing is written to stderr. Output files land in the
current directory unless `--prefix` says otherwise.

The three commands most people want:

| goal | command |
| --- | --- |
| kinship for all pairs | `open-king -b study.bed --kinship` |
| relatives up to 2nd degree, with IBD segments | `open-king -b study.bed --related --degree 2` |
| find duplicate / MZ samples | `open-king -b study.bed --duplicate` |

---

## 3. Input

### `-b <file>.bed`, `-B <file>.bed`

The one short option, and the only way to name input. Case-insensitive, and the value may be
attached or separate — `-b study.bed` and `-bstudy.bed` are the same. Only **one** fileset
per run.

The argument names the `.bed` **file itself**, not the fileset stem, and it must end in a
literal lower-case `.bed`. A stem is not accepted:

```
$ open-king -b multifam --kinship

FATAL ERROR - 
Genotype file multifam cannot be opened
```

A readable file whose name does not end in `.bed` is rejected even when its contents are a
valid `.bed` (the same file copied to `plain`, with `plain.bim`/`plain.fam` alongside):

```
$ open-king -b plain --kinship

FATAL ERROR - 
Please use PLINK binary format as input.
```

`UPPER.BED` fails the same way — the suffix check is case-**sensitive**, unlike the option
names. Both behaviours match the reference binary exactly.

`.bim` and `.fam` are derived by replacing the `.bed` suffix, so `sub/study.bed` implies
`sub/study.bim` and `sub/study.fam`. Override either with [`--bim`](#--bim-file) /
[`--fam`](#--fam-file).

### File requirements

| file | requirement |
| --- | --- |
| `.bed` | PLINK 1 binary, magic bytes `0x6c 0x1b`, **SNP-major** (mode byte `0x01`), long enough for `ceil(N/4) × M` genotype bytes |
| `.bim` | exactly 6 whitespace-separated columns `CHR ID CM BP A1 A2`; extra trailing columns are ignored, fewer than 6 is fatal |
| `.fam` | exactly 6 whitespace-separated columns `FID IID FA MO SEX PHENO`; `(FID, IID)` must be unique under ASCII case-folding |

Each has its own fatal error. A corrupted magic:

```
FATAL ERROR - 
Please use either PLINK or KING binary format as input.
```

An individual-major `.bed` (mode byte `0x00`):

```
FATAL ERROR - 
Currently only SNP-major mode can be analyzed.
```

A `.bed` shorter than the map and sample count require (here truncated to 5 000 bytes):

```
FATAL ERROR - 
Not enough genotypes at the 999th marker
```

A short `.fam` line:

```
FATAL ERROR - 
t.fam: line 1 has 5 fields, expected 6
```

A repeated `(FID, IID)`:

```
FATAL ERROR - 
Please correct problems with pedigree structure
```

Nothing else in the map is filtered. Monomorphic SNPs, `0` alleles, low call rate, low MAF,
and duplicate SNP IDs or positions are all **kept** — there is no MAF or call-rate filter on
the relatedness path, and missingness is handled pairwise, per pair, per SNP. That is
measured in [BEHAVIOR.md Q3](BEHAVIOR.md#q3--snp-inclusion-rules).

### Two hard requirements that are easy to miss

**1. A1 must be the minor allele.** PLINK 1.9's `--make-bed` writes it that way, so filesets
that came through PLINK are fine; hand-built or hand-edited ones may not be. The reference
binary refuses such a fileset outright. Recoding a corpus fileset so that the A1 homozygote
is the *common* genotype and running KING 2.3.2 on it:

```
$ king -b major.bed --ibs          # KING 2.3.2

FATAL ERROR -
Too many first alleles as the major allele (~77.9%). Please use plink1.9 --make-bed to regenerate the genotype data again.

$ open-king -b major.bed --ibs          # open-king
FATAL ERROR -
Too many first alleles as the major allele (~77.9%). Please use plink1.9 --make-bed to regenerate the genotype data again.
```

Both binaries now enforce the same stable rule: inspect the first 4,096 retained autosomal
markers, count a marker when its observed A1/A1 call count is greater than its A2/A2 count,
and abort when more than 10% qualify. Thus 409/4,096 passes and 410/4,096 aborts; markers
after the window do not enter the percentage. `--kinship`, `--duplicate` and `--autoQC` are
exempt, as are the small-sample `--related`/`--ibdseg`/clustering paths that KING replaces
or disables. The fatal console, exit status and files written before the fatal are compared
against KING by `tests/parity/probes/a1_major.py`.

KING's checker reads unstable tail state when fewer than 4,096 autosomal markers exist—it
has produced random false fatals on a valid eight-sample fixture. open-king deliberately
skips the gate on such short maps instead of reproducing that unsafe behavior. Most
relatedness columns are allele-orientation invariant, but `HomIBS0` and segment-derived
classification are not, so **regenerate questionable inputs with `plink1.9 --make-bed`.**

**2. The `.bim` must be sorted by `(chromosome, position)`, ascending,** for anything that
calls IBD segments — which is every analysis except `--kinship` and `--duplicate`. Both
binaries now refuse segment work on an unsorted map at the same point and with the same
diagnostic. On a map whose positions run backwards within each chromosome:

```
$ open-king -b unsortedpos.bed --ibdseg     # both -> kingsplitped.txt only
Positions unsorted: rs1_1009689 at 65904473, rs1_1055261 at 65851170.
  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.
```

On a map whose chromosomes run 22 → 1 with positions ascending inside each:

```
$ open-king -b unsortedchr.bed --ibdseg     # both -> kingsplitped.txt only
Chromosomes unsorted: rs22_14205438 on chr 22, rs21_1002722 on chr 21.
  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.
```

Pairwise *counting* is order-independent, so `--kinship`, `--duplicate` and the IBS columns
of `--ibs` are unharmed by map order — but every segment caller reads the map as a sequence,
and an unsorted one is meaningless to it. Sort the map.

### Chromosome codes

Codes are read from `.bim` column 1 and classified before anything else runs. At the default
`--sexchr 23`:

| code | class |
| --- | --- |
| `1`–`22` | autosome |
| `25`, `XY` | autosome — **pooled with the autosomes**, not held aside |
| `23`, `X` | X chromosome (its own analysis) |
| `24`, `Y` | Y chromosome (excluded) |
| `26`, `MT` | mitochondrial (excluded) |
| anything else (`0`, `27`, `chr1`, contig names) | dropped at map load; not even counted |

The `sexchr` fileset carries 2 000 SNPs each on chr 1 and 2, 1 500 on 23, 300 on 24, 150 on
25 and 50 on 26:

```
$ open-king -b sexchr.bed --kinship
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs, 300 Y-chromosome SNPs, 50 mitochondrial SNPs
  PLINK maps loaded: 6000 SNPs
```

`4150 = 2000 + 2000 + 150`. The XY count is inside the autosome total, not beside it.
Alphabetic spellings give the identical partition:

```
$ open-king -b alphachr.bed --kinship   # same data, chromosomes written X / Y / XY / MT
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs, 300 Y-chromosome SNPs, 50 mitochondrial SNPs
  PLINK maps loaded: 6000 SNPs
```

Codes not in any class are removed and reported. Rewriting chr 24 as `0` and chr 26 as
`chr1`:

```
$ open-king -b unknownchr.bed --kinship
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs
  350 other SNPs are removed.
  PLINK maps loaded: 5650 SNPs
```

`--sexchr` moves the whole partition; see [`--sexchr`](#--sexchr-n).

---

## 4. Output

### Where files go

Every output path is `<prefix>` followed by a fixed suffix, and `<prefix>` defaults to
`king`, so files land in the **current working directory** as `king.kin`, `king.kin0`,
`kingallsegs.txt` and so on. `--prefix` may contain a directory; the directory must already
exist.

### `--prefix` is a concatenation, not a stem

There is **no separator inserted**. `--prefix ZZ_` gives `ZZ_.kin`:

```
$ open-king -b multifam.bed --kinship --prefix ZZ_
Within-family kinship data saved in file ZZ_.kin
Between-family kinship data saved in file ZZ_.kin0
$ ls ZZ_*
ZZ_.kin
ZZ_.kin0
```

and the QC and segment files, whose suffixes have no leading dot, come out as
`<prefix>bySample.txt`, not `<prefix>.bySample.txt`:

```
$ open-king -b multifam.bed --bysample --prefix ZZ_
QC statistics by samples saved in file ZZ_bySample.txt
$ ls ZZ_*
ZZ_allsegs.txt
ZZ_bySample.txt
```

So `--prefix study` gives `study.kin` *and* `studyallsegs.txt`. If you want
`study_allsegs.txt`, pass `--prefix study_` and accept `study_.kin`. A prefix ending in `/`
puts everything in a directory; `--prefix out/study` gives `out/study.kin`.

### What each analysis writes

Measured by running each analysis alone in an empty directory. `<p>` is the prefix.

| analysis | files |
| --- | --- |
| `--kinship` | `<p>.kin`, `<p>.kin0`, plus `<p>X.kin` and `<p>X.kin0` when the X gate opens |
| `--related` | `<p>.kin`, `<p>.kin0`, `<p>allsegs.txt`, plus `<p>X.kin` and — when a `.kin0` is written on a map with usable X segments — `<p>X.kin0` ([OUTPUTS.md](OUTPUTS.md#prefixxkin0)) |
| `--duplicate` | `<p>.con` |
| `--ibs` | `<p>.ibs`, `<p>.ibs0`, `<p>allsegs.txt` |
| `--ibdseg` | `<p>.seg`, `<p>allsegs.txt`, `<p>splitped.txt`, plus `<p>X.seg` |
| `--unrelated` | `<p>unrelated.txt`, `<p>unrelated_toberemoved.txt`, `<p>allsegs.txt` |
| `--cluster` | `<p>cluster.kin`, `<p>allsegs.txt`, plus `<p>updateids.txt` |
| `--build` | `<p>updateparents.txt`, `<p>build.log`, `<p>allsegs.txt`, plus `<p>updateids.txt` |
| `--bysample` | `<p>bySample.txt`, `<p>allsegs.txt` |
| `--bySNP` | `<p>bySNP.txt`, `<p>allsegs.txt` |
| `--autoQC` | `<p>_autoQC_Summary.txt`, `<p>_autoQC_snptoberemoved.txt`, `<p>_autoQC_sampletoberemoved.txt`, and `<p>_autoQC_updatesex.txt` when a sex was inferred from the X data (of the corpus, only `sexchr` produces it) |

"plus" marks a conditional file. The X files need X data and their own gate (below);
`updateids.txt` appears only when families actually merge — on `bigish` both `--cluster` and
`--build` write it, on `multifam` neither does. `--autoQC` is the one analysis whose names
carry an underscore of their own, so the default prefix gives `king_autoQC_Summary.txt`.

Several files are **conditional**, which trips people up when a file they expected is
absent. The rules, all measured
([BEHAVIOR.md Q7](BEHAVIOR.md#q7--output-file-existence)):

* `.kin0` needs at least two distinct FIDs; `.kin` needs at least one family with two or
  more members.
* `--related`'s `.kin0` additionally needs the between-family stage to confirm relatives,
  which needs **N ≥ 100 samples**. Below that the run prints `No close relatives are
  inferred.` and writes no `.kin0` even when duplicates exist. This is a reference bug that
  open-king reproduces deliberately.
* `.con` is written for every run with N < 100 (header-only if there are no duplicates); at
  N ≥ 100 it appears **only if a duplicate pair was found**. On `bigish` (200 samples, no
  duplicates) `--duplicate` writes nothing at all and prints
  `No duplicates are found with heterozygote concordance rate > 80%.`
* When the dataset has exactly one distinct FID, `.kin` is **truncated**: rows are flushed
  every 65 536 bytes and the final partial buffer is never written, so a small single-family
  dataset yields a zero-byte `.kin`. Reproduced from the reference on purpose.
* The X files have their own gates, and they are not the same gate. `--kinship` writes
  `X.kin`/`X.kin0` only with **≥ 512 X markers, more than one family, and no `--degree`**;
  `--ibdseg` writes `X.seg` only **with** a `--degree`.

All three conditions are visible on the same fileset (`sexchr`, 1 500 X markers, 5
families):

```
$ open-king -b sexchr.bed --kinship            -> king.kin king.kin0 kingX.kin kingX.kin0
$ open-king -b sexchr.bed --kinship --degree 2 -> king.kin king.kin0
$ open-king -b sexchr.bed --ibdseg             -> king.seg kingallsegs.txt kingsplitped.txt
$ open-king -b sexchr.bed --ibdseg --degree 2  -> king.seg kingX.seg kingallsegs.txt kingsplitped.txt
```

and cutting the same fileset down to 400 X markers removes the X files from the first
command:

```
$ open-king -b x400.bed --kinship              -> king.kin king.kin0
```

### Column headers

```
.kin    (--kinship)  FID ID1 ID2 N_SNP Z0 Phi HetHet IBS0 Kinship Error
.kin    (--related)  FID ID1 ID2 N_SNP Z0 Phi HetHet IBS0 HetConc HomIBS0 Kinship
                     IBD1Seg IBD2Seg PropIBD InfType Error
.kin0   (--kinship)  FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 Kinship
.kin0   (--related)  FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 HetConc HomIBS0 Kinship
                     [IBD1Seg IBD2Seg PropIBD InfType]
.ibs                 FID ID1 ID2 Z0 Phi N_SNP N_IBS0 N_IBS1 N_IBS2 NHetHet NHomHom
                     N_Het1 N_Het2 IBS Dist HetConc Het2|1 Het1|2 HomConc Kinship
                     [MaxIBD2 Pr_IBD2]
.ibs0                same, without Z0 Phi
.con                 FID1 ID1 FID2 ID2 N N_IBS0 N_IBS1 N_IBS2 Concord HomConc HetConc
.seg                 FID1 ID1 FID2 ID2 IBD1Seg IBD2Seg PropIBD InfType
cluster.kin          FID ID1 ID2 Sex1 Sex2 N_SNP HetHet IBS0 HetConc HomIBS0 Kinship
                     IBD1Seg IBD2Seg PropIBD InfType
allsegs.txt          Segment Chr StartMB StopMB Length N_SNP StartSNP StopSNP
bySample.txt         FID IID FA MO SEX N_SNP Missing Heterozygosity N_pair N_MIp
                     Err_MIp N_trio N_MIt Err_MIt MI_Removal
bySNP.txt            SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate
                     N_PO N_HomPO N_errPO Err_InPO Err_InHomPO N_trio N_HetOff
                     N_errTrio Err_InTrio Err_InHetTrio
```

Bracketed columns are present only when the map yields ≥ 100 Mb of usable IBD segment
([BEHAVIOR.md Q8](BEHAVIOR.md#q8--ibs--ibs0-column-set-variation)).

Files are tab-separated, **except** `bySample.txt`, `bySNP.txt` and `splitped.txt`, which
are space-separated, and `_autoQC_Summary.txt`, which is a fixed-width space-padded report.
(The two `_autoQC_*toberemoved.txt` files are tab-separated like the rest.)

---

## 5. Option reference

Sections are the binary's own banner groups, in banner order. 46 long options plus `-b`.

Types: **switch** takes no value; **int** and **double** take a numeric token; **string**
takes the next token whatever it is. What "takes" means precisely is
[§6](#6-how-the-parser-behaves) — the rules are unusual and worth reading once.

### Close Relative Inference

#### `--related`

*Switch, default off.* The headline analysis: kinship for every pair, IBD segments for the
close ones, and an inferred relationship label. Writes `<p>.kin` (16 columns) and
`<p>.kin0` (14 columns) with `IBD1Seg`, `IBD2Seg`, `PropIBD` and `InfType`
(`Dup/MZ`, `PO`, `FS`, `2nd`, `3rd`, `4th`), plus `<p>allsegs.txt`.

Reports pairs at **degree 1 by default** — the only analysis with a non-zero default degree.
Widen it with [`--degree`](#--degree-d):

```
$ open-king -b bigish.bed --related
  Final Stage (with 50000 SNPs): 3 pairs of relatives (up to 1st-degree) are confirmed
Between-family relatives (kinship >= 0.17678) saved in file king.kin0
Note only duplicates and 1st-degree relatives are included in the inference.
  Specifying '--degree 2' if a higher degree relationship inference is needed.
```

A between-family pair reaches `.kin0` if `Kinship >= 2^-(d+1.5)` **or**
`PropIBD > 2^-(d+0.5)` — a disjunction, unlike `--kinship`'s plain kinship cut.

Two size gates. Under **10 samples** the run silently becomes a `--kinship` run — same files,
same columns, and `--degree` is discarded:

```
$ open-king -b trio.bed --related
--related is replaced with --kinship for a small sample size.
Autosome genotypes stored in 79 words for each of 3 individuals.

Options in effect:
	--kinship
```

Under **100 samples** the between-family stage reports `No close relatives are inferred.`
and writes no `.kin0` regardless of what the data holds. Both are reference behaviour,
reproduced deliberately.

#### `--duplicate`

*Switch, default off.* Finds duplicate and MZ pairs by heterozygote concordance and writes
`<p>.con`. The concordance floor is [`--minConc`](#--minconc-x) (default 0.80):

```
$ open-king -b dups.bed --duplicate
2 pairs of duplicates with heterozygote concordance rate > 80% are saved in file king.con
```

See [§4](#what-each-analysis-writes) for when `.con` is absent rather than empty — at
N ≥ 100 with no duplicates in the data, no file is written at all.

### Pairwise Relatedness Inference

#### `--kinship`

*Switch, default off.* KING-robust kinship for every pair, and nothing else: no segments, no
`InfType`. Writes the 10-column `<p>.kin` and 8-column `<p>.kin0`, and — with ≥ 512 X
markers, > 1 family and no `--degree` — `<p>X.kin` and `<p>X.kin0`.

No sample-size gate: it works on 3 samples. This is the analysis to reach for when you want
numbers for every pair rather than a filtered relative list.

#### `--ibdseg`

*Switch, default off.* Pairwise IBD-segment inference. Writes `<p>.seg`
(`IBD1Seg IBD2Seg PropIBD InfType`) and `<p>allsegs.txt`; it also writes
`<p>splitped.txt` when at least one family has two members or a singleton names a parent.
With a `--degree`, it may also write `<p>X.seg`.

Affected by [`--seglength`](#--seglength-mb) (the minimum length a called segment must
reach) and [`--degree`](#--degree-d) (which pairs are reported). Under **5 samples** it
becomes a `--kinship` run:

```
$ open-king -b trio.bed --ibdseg
--kinship analysis carried out instead for such a small sample size.
```

Requires a `(chr, bp)`-sorted `.bim` — see [§3](#two-hard-requirements-that-are-easy-to-miss).

#### `--ibs`

*Switch, default off.* Full IBS and concordance statistics for every pair: `<p>.ibs`
(20 or 22 columns), `<p>.ibs0` (19 or 21), plus `<p>allsegs.txt`. The two trailing columns
`MaxIBD2` and `Pr_IBD2` appear only when the usable segment total reaches 100 Mb, and carry
the literal `-9` for pairs below `Kinship 2^-3.5`.

`--degree` does **not** filter `.ibs`/`.ibs0`. On `bigish` the row counts are 573 and 19 327
at degrees 1, 2 and 3 alike.

Note that `--ibs`'s `Pr_IBD2` and `--ibdseg`'s `IBD2Seg` are computed by two different
callers and legitimately disagree — that is faithful to the reference, not a bug here
([BEHAVIOR.md](BEHAVIOR.md), Q10).

#### `--makeGRM`

*Switch, default off.* **Not implemented** — see [§8](#8-accepted-compatibility-spellings-outside-product-scope).

### Inference Parameter

#### `--degree <d>`

*Int, default 0 (= unset).* The relationship degree to report. It is a **reporting filter**,
not a computation switch: everything is still computed, fewer rows are written.

* **`--kinship`** — filters `<p>.kin0` and nothing else. `.kin` is never filtered, not its
  rows and not its columns:

  ```
  --degree 1: .kin0 = 8 rows, .kin = 40 rows
  --degree 2: .kin0 = 32 rows, .kin = 40 rows
  --degree 3: .kin0 = 52 rows, .kin = 40 rows
  --degree 4: .kin0 = 63 rows, .kin = 40 rows
  (none)    : .kin0 = 150 rows, .kin = 40 rows
  ```

  The test is `Kinship >= 2^-(d + 1.5)` against the **exact IEEE double**, not against the
  five-decimal value the console prints. At `d = 2` the cut is 0.088388…; on `multifam` all
  32 kept rows are ≥ 0.0883883 and the largest excluded pair sits at 0.0843. A pair at
  0.1767775 — above `2^-2.5` but below the printed `0.17678` — is kept at `--degree 1`
  ([BEHAVIOR.md Q5](BEHAVIOR.md#q5----degree-semantics)). The console echoes the rounded form:

  ```
  $ open-king -b multifam.bed --kinship --degree 2
  Between-family kinship data (up to degree 2, 32 pairs in total) saved in file king.kin0
  ```

* **`--related`** — same `.kin0`-only filtering, but the test is the disjunction
  `Kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)`, and the **default is degree 1** rather
  than unset.

* **`--ibdseg`** — filters `<p>.seg` rows on the segment estimates:
  `PropIBD > 2^-(d+0.5)`, plus, at `d = 1` only, `IBD2Seg >= 0.08`. A **negative** degree
  inverts the comparison and reports the complement, which is exactly what it sounds like:

  ```
  --degree  2 -> 442 of 763 rows
  --degree -2 -> 321 of 763 rows      (442 + 321 = 763)
  ```

* **`--unrelated`, `--cluster`, `--build`** — sets the clustering degree
  (`Clustering up to 2nd-degree relatives in families...`). `--degree 0`, a negative degree
  and an absent one all cluster at 1st degree. On the shipped corpus the selected sets do
  not change with degree.

* **`--ibs`** — no effect on row counts.

`--degree 0` is treated as *unset*: it is not echoed in the banner and filters nothing.
Bare `--degree` with no value means **1** (see [§6](#6-how-the-parser-behaves)). Degrees
above the useful range work as the formula says. Note the reference's typo, reproduced:

```
$ open-king -b bigish.bed --related --degree 3
  59 pairs of relatives (up to 3nd-degree) are identified
```

#### `--noscreen [n]`

*Int.* **Accepted and echoed; it changes nothing.** In the reference it was meant to disable
the two-stage screen. Running `bigish --related --degree 2` with and without
`--noscreen 1` gives byte-identical output files and an identical screening line on stdout;
the only difference anywhere is the banner echo. It has no reader in the codebase outside
the parser: no analysis pass consults the value, so accepting it is a parse-surface
obligation and nothing else.

Its default is the notorious `-1717986816`. That is not a value — in the reference it is
uninitialised memory that overlaps `--minConc`'s storage, and open-king reproduces the
overlap byte for byte, including the way `--minConc 0.9` changes it to `-858993408`.
Ignore it.

**Being inert here is a deliberate match to the capture binary, not an omission.** A second
compilation of KING 2.3.2 from the published source behaves differently: on
`multifam --related --degree 1 --noscreen` it bypasses the screen and writes a `king.kin0`
with eight between-family pairs, where the binary the goldens were captured from writes no
such file and reports `No close relatives are inferred.` open-king matches the capture
binary. Since the option carries an undefined value, its effect is undefined too, and no
single behaviour satisfies both builds. The full measurement is in
[`PARITY.md` §5.13](PARITY.md#513-a-second-build-of-king-232-agrees-on-every-output-file-but-one-and-that-one-is---noscreen).

#### `--seglength <Mb>`

*Double, default 3 Mb.* Minimum length, in megabases, for a called IBD segment. Affects
`--ibdseg`, `--related`, and every analysis that runs the segment engine.

**Accepted only strictly inside (0.99, 10.01) Mb.** Anything else silently reverts to 3 Mb,
with a notice:

```
$ open-king -b multifam.bed --ibdseg --seglength 5
KING starts at Fri Aug 14 22:35:54 2026
Minimum segment length is set as 5000000 bp
.Loading genotype data in PLINK binary format...

$ open-king -b multifam.bed --ibdseg --seglength 0.99
KING supports minimum segment length from 1 to 10 Mb at the moment.
Default seglength of 3Mb is used.

$ open-king -b multifam.bed --ibdseg --seglength 10.01
KING supports minimum segment length from 1 to 10 Mb at the moment.
Default seglength of 3Mb is used.
```

`1`, `1.0001`, `10` and `10.009` are all accepted; `0.5`, `0`, `11` and `12` are not. The
revert is genuine, not cosmetic — on `bigish` the `.seg` md5 for `--seglength 0.5`,
`--seglength 0`, `--seglength 11` and no flag at all are the same
`4fdddc6b00be91bf9a29bd5df51b2a15`. (That stray `.` on the next line is the reference's own
missing newline, reproduced.)

**`--seglength` is never echoed under `Options in effect:`**, by any analysis, even on the
runs where it changes the output bytes. The `Minimum segment length is set as <n> bp` line
above the block is the only report of the value:

```
$ open-king -b bigish.bed --ibdseg --seglength 5
Minimum segment length is set as 5000000 bp
...
Options in effect:
	--ibdseg
```

That is reference behaviour, not an omission here. See [`--minConc`](#--minconc-x) for the
rest of the echo rules.

#### `--minConc <x>`

*Double, default 0.80.* Heterozygote-concordance floor for `--duplicate`. Echoed to two
decimals and printed as a percentage:

```
$ open-king -b dups.bed --duplicate --minConc 0.99
2 pairs of duplicates with heterozygote concordance rate > 99% are saved in file king.con
```

Out of `[0, 1]` it warns and is used anyway, percentage and all:

```
$ open-king -b dups.bed --duplicate --minConc 1.5
minConc value is out of range and not specified.
No duplicates are found with heterozygote concordance rate > 150%.
```

**It is echoed under `Options in effect:` for eight analyses and dropped by three.** The
eight that echo it are `--duplicate`, `--ibs`, `--autoQC`, `--unrelated`, `--build`,
`--bysample`, `--bySNP` and `--cluster`, which share one echo list. `--kinship`, `--related`
and `--ibdseg` build their own lists, which carry only `--degree`, `--cpus` and `--prefix`,
so the value is dropped there without a word:

```
$ open-king -b bigish.bed --duplicate --minConc 0.9
Options in effect:
	--duplicate
	--minConc 0.9

$ open-king -b bigish.bed --kinship --minConc 0.9
Options in effect:
	--kinship
```

The value is in effect either way; only the echo differs. Both halves are reference
behaviour, and both matter to anyone diffing stdout.

### Relationship Application

#### `--unrelated`

*Switch, default off.* Greedy maximal unrelated subset. Writes `<p>unrelated.txt` (the
subset to keep) and `<p>unrelated_toberemoved.txt` (the complement), plus
`<p>allsegs.txt`. Both are two-column `FID<TAB>IID`, no header.

```
$ open-king -b bigish.bed --unrelated --degree 2
Clustering up to 2nd-degree relatives in families...
A list of 84 unrelated individuals saved in file kingunrelated.txt
An alternative list of 116 to-be-removed individuals saved in file kingunrelated_toberemoved.txt
```

#### `--cluster`

*Switch, default off.* Merges families connected by inferred relatedness and reports the
pairs inside the merged clusters. Writes `<p>allsegs.txt` when the map has usable segments
and `<p>updateids.txt` when families actually merge. `<p>cluster.kin` additionally requires
usable segments, because its relationship columns are segment-derived:

```
$ open-king -b bigish.bed --cluster   -> kingallsegs.txt kingcluster.kin kingupdateids.txt
$ open-king -b multifam.bed --cluster -> kingallsegs.txt
```

Unlike every other analysis, selecting `--cluster` suppresses the
`The following analyses will run separately:` line entirely, whatever else was asked for.

#### `--build`

*Switch, default off.* Pedigree reconstruction from the genotypes. Writes
`<p>updateparents.txt`, `<p>build.log`, `<p>allsegs.txt`, and `<p>updateids.txt` when
families merge. `build.log` is the one output file in the project that is not byte-identical
to the reference everywhere — see [PARITY.md §6.2](PARITY.md).

**`--build` announces `updateids.txt` on runs where it writes no such file.** The file is
written once, at the point family clustering merges two or more families, and announced there
with the `The following families are found to be connected` table above it. The line in the
reconstruction tail is printed unconditionally, whether or not that write happened, so a run
with no merges names a file that is not on disk:

```
$ open-king -b multifam.bed --build     # no families merge
Details of pedigree reconstruction are available in log file kingbuild.log
Update-ID information is saved in file kingupdateids.txt
No pedigrees can be reconstructed.
$ ls king*
kingallsegs.txt  kingbuild.log  kingupdateparents.txt

$ open-king -b bigish.bed --build      # three pairs of families merge
$ ls king*
kingallsegs.txt  kingbuild.log  kingupdateids.txt  kingupdateparents.txt
```

On `bigish` the line therefore appears twice, once per site. This is reference behaviour,
reproduced deliberately. [`--cluster`](#--cluster) does not share it: it writes and announces
the file at the same point, so its line and the file always agree.

### QC Report

#### `--bysample`

*Switch, default off.* Per-sample QC: call rate, heterozygosity, Mendelian-error counts.
Writes `<p>bySample.txt` and `<p>allsegs.txt`. Six header variants depending on what the
pedigree supports.

#### `--bySNP`

*Switch, default off.* Per-SNP QC: allele frequency, genotype counts, call rate,
parent-offspring and trio error rates. Writes `<p>bySNP.txt` and `<p>allsegs.txt`. Three
header variants.

Note the spelling: the option is `--bySNP`, matched case-insensitively, but it is echoed as
`--bysnp` in the "will run separately" line. (So are `--makeGRM` → `--grm` and
`--lmm` → `--mtscore`.)

#### `--roh`

*Switch, default off.* **Not implemented** — see [§8](#8-accepted-compatibility-spellings-outside-product-scope).

#### `--autoQC`

*Switch, default off.* The packaged call-rate and sex QC pipeline. Writes
`<p>_autoQC_Summary.txt`, `<p>_autoQC_snptoberemoved.txt`,
`<p>_autoQC_sampletoberemoved.txt`, and `<p>_autoQC_updatesex.txt` when a sex was inferred
from the X data. It writes no `.kin` and does not need `--kinship`.

```
$ open-king -b missing.bed --autoQC
Auto-QC step 1: Apply SNP call rate filter 80.0% on 10000 SNPs (in 6 samples)
  1569 autosome SNPs have call rate < 80.0%
  0 X-chr SNPs have call rate < 80.0%
Auto-QC step 2: Apply sample call rate filter 95.0% on 6 samples (with 5232 SNPs)
  2 samples have call rate < 95.0%
Auto-QC step 3: Apply SNP call rate filter 95.0% on 5232 SNPs (in 4 samples)
  130 SNPs have call rate < 95.0%
  0 chr-X SNPs have call rate < 95.0%

$ cat king_autoQC_Summary.txt
Step Description                                            Subjects  SNPs      
1    Raw data counts                                        6         10000     
1.1  SNPs with very low call rate < 80% (removed)                     (1569)
1.2  Monomorphic SNPs (removed)                                       (3199)
1.3  Sample call rate < 95% (removed)                       (2)
1.4  SNPs with call rate < 95% (removed)                              (130)
3    Generate Final Study Files                             
     Final QC'ed data                                       4         5086
```

### QC Parameter

#### `--callrateN <x>`

*Double, default 0.95.* The **sample** call-rate threshold used by `--autoQC` step 2. Only
`--autoQC` reads it.

#### `--callrateM <x>`

*Double, default 0.95.* The **SNP** call-rate threshold used by `--autoQC` step 3. Only
`--autoQC` reads it.

**It also sets step 1's pre-filter, which is not the fixed 80 % it looks like.** The rule is
`min(0.8, 0.1 * (trunc(callrateM * 10) - 1))`: one decimal digit of `--callrateM`, truncated
toward zero, less one tenth, capped at 80 %. It reads as fixed only because it saturates from
`--callrateM 0.9` upward, which is where the 0.95 default and every other example on this
page sit. Swept on `missing`:

| `--callrateM` | step 1 filter | step 3 filter |
| --- | ---: | ---: |
| 0.1 | 0.0% | 10.0% |
| 0.3 | 20.0% | 30.0% |
| 0.5 | 40.0% | 50.0% |
| 0.7 | 60.0% | 70.0% |
| 0.75 | 60.0% | 75.0% |
| 0.8 | 70.0% | 80.0% |
| 0.85 | 70.0% | 85.0% |
| 0.9 | 80.0% | 90.0% |
| 0.95 | 80.0% | 95.0% |
| 1.0 | 80.0% | 100.0% |

```
$ open-king -b missing.bed --autoQC --callrateM 0.5
Auto-QC step 1: Apply SNP call rate filter 40.0% on 10000 SNPs (in 6 samples)
  220 autosome SNPs have call rate < 40.0%
  0 X-chr SNPs have call rate < 40.0%
```

There is no floor, only the 80 % ceiling: `--callrateM 0` gives a step 1 filter of `-10.0%`
and `--callrateM -0.55` gives `-60.0%`, which is also what pins the truncation as toward zero
rather than downward. `--callrateN` does not move step 1; only `--callrateM` does.

Both together:

```
$ open-king -b missing.bed --autoQC --callrateN 0.99 --callrateM 0.9
Auto-QC step 1: Apply SNP call rate filter 80.0% on 10000 SNPs (in 6 samples)
  1569 autosome SNPs have call rate < 80.0%
Auto-QC step 2: Apply sample call rate filter 99.0% on 6 samples (with 5232 SNPs)
  3 samples have call rate < 99.0%
Auto-QC step 3: Apply SNP call rate filter 90.0% on 5232 SNPs (in 3 samples)
  22 SNPs have call rate < 90.0%
```

The `Summary.txt` rows still read `< 95%` whatever you pass — that text is a fixed label in
the reference, and reproduced:

```
1.3  Sample call rate < 95% (removed)                       (3)
1.4  SNPs with call rate < 95% (removed)                              (22)
```

### Population Structure

#### `--pca`
#### `--mds`

*Switches, default off.* **Not implemented** — see
[§8](#8-accepted-compatibility-spellings-outside-product-scope).

### Structure Parameter

#### `--projection [n]`
#### `--pcs [n]`

*Ints, default 0.* Parameters of `--pca`/`--mds`. Parsed and echoed in the banner, then
**rejected**: naming either one is fatal, exit 1, before the input is opened.
See [§8](#8-accepted-compatibility-spellings-outside-product-scope).

### Quantitative Trait GWAS

#### `--lmm`

*Switch, default off.* **Not implemented.** Echoed as `--mtscore` in the "will run
separately" line.

### Binary Trait GWAS

#### `--tdt`
#### `--gdt`

*Switches, default off.* **Not implemented.**

### Association Model

#### `--trait <name>`
#### `--covariate <name>`

*Strings, default empty.* Parameters of the association analyses. Parsed and echoed in the
banner, then **rejected**: naming either one is fatal, exit 1, before the input is opened.
See [§8](#8-accepted-compatibility-spellings-outside-product-scope).

The string-swallowing rule of [§6](#6-how-the-parser-behaves) still decides *which* name the
fatal reports. A string option takes the next token unconditionally, option or not, so
`--trait --related` sets the trait to the literal `--related`, leaves `--related` off, and
the run is rejected for `--trait` alone:

```
$ open-king -b multifam.bed --trait --related --kinship

FATAL ERROR - 
open-king's minimal relatedness product does not implement: --trait.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
```

#### `--maxP <p>`

*Double, unset by default.* Parameter of the association analyses. Parsed and echoed in the
banner, then **rejected**: naming it is fatal, exit 1, before the input is opened, whatever
value follows.

```
$ open-king -b multifam.bed --kinship --maxP 0.05

FATAL ERROR - 
open-king's minimal relatedness product does not implement: --maxP.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
```

Measured, `0`, `0.05`, `1`, `2`, `3`, `-1` and a bare `--maxP` all produce that same block and
exit 1. The gate fires on the option's presence, not on its value.

A range check on `p` does exist in `main.rs`, and it is **unreachable**: it sits after the
product-scope gate, so `p-value [x] outside range in ninv()` and the `0 < p < 2` window it
enforces cannot be observed from the command line. See
[§7](#7-exit-status-and-fatal-errors).

### Association Method Parameter

#### `--invnorm`

*Switch, default off.* Parameter of the association analyses. Parsed and echoed in the
banner, then **rejected**: naming it is fatal, exit 1, before the input is opened.
See [§8](#8-accepted-compatibility-spellings-outside-product-scope).

### Genetic Risk Score

#### `--risk`

*Switch, default off.* **Outside the minimal product scope.** Naming it is fatal, exit 1,
before the input is opened:

```
$ open-king -b multifam.bed --risk

FATAL ERROR - 
open-king's minimal relatedness product does not implement: --risk.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
```

`--risk --model m.txt` names both options in one fatal. Two reference behaviours around
`--risk` sit behind that gate and are therefore **unreachable**: the
`Please use --model <file> to specify a risk model.` fatal for `--risk` without `--model`,
and the `--risk --model` path in the loader that skips the early `.bed` probe, which in the
reference reports a missing fileset as `Pedigree file <name>.fam cannot be opened` rather
than `Genotype file <name>.bed cannot be opened`. Both are implemented; no command line
reaches either. See [§7](#7-exit-status-and-fatal-errors).

#### `--model <file>`

*String, default empty.* The risk model file. Parsed and echoed in the banner, then
**rejected**: naming it is fatal, exit 1, before the input is opened. It is never read.

#### `--prevalence <x>`
#### `--noflip`

*Double / switch.* Parameters of `--risk`. Parsed and echoed in the banner, then
**rejected**: naming either one is fatal, exit 1, before the input is opened.
See [§8](#8-accepted-compatibility-spellings-outside-product-scope).

### Computing Parameter

#### `--cpus <n>`

*Int, default = the machine's available parallelism* (16 on the host these examples ran on).
Worker threads. Echoed, and reported:

```
$ open-king -b multifam.bed --kinship --cpus 3
3 CPU cores are used.
```

**`--cpus` changes no printed digit in any output file.** It is a pure performance knob.
Verified two ways — `bigish --kinship` at 1, 4, 16 and default threads gives the same md5 for
both files, and a combined `--related --ibs --ibdseg` run on `multifam` at 1 and 16 threads
gives six byte-identical files:

```
$ diff <(cd c1 && md5 -q *) <(cd c16 && md5 -q *) \
      && echo 'ALL OUTPUT FILES IDENTICAL across --cpus 1 and --cpus 16'
ALL OUTPUT FILES IDENTICAL across --cpus 1 and --cpus 16
```

The only stdout difference is the echoed value, the `N CPU cores are used` line, and how far
the `0%1%2%…` progress counter gets before the work finishes. Normalise stdout before
diffing runs at different thread counts.

### Optional Input

#### `--fam <file>`

*String, default empty.* Replaces the derived `.fam` path. Used **verbatim** — `--fam alt`
reads a file literally named `alt`, with no suffix appended.

#### `--bim <file>`

*String, default empty.* Same, for the `.bim`.

```
$ open-king -b geno.bed --fam alt.ped --bim alt.map --kinship
Read in PLINK fam file alt.ped...
  PLINK pedigrees loaded: 20 samples
Read in PLINK bim file alt.map...
  PLINK maps loaded: 15000 SNPs
Read in PLINK bed file geno.bed...
```

The `.bed` path is always exactly what `-b` said; there is no `--bed`.

#### `--phefile <file>`
#### `--covfile <file>`
#### `--prunedsnp <file>`

*Strings, default empty.* Inputs to the association and structure analyses. Parsed and
echoed in the banner, then **rejected**: naming any of the three is fatal, exit 1, before the
input is opened. The path itself is never looked at, so a nonexistent one gives the same
fatal as a real one. See [§8](#8-accepted-compatibility-spellings-outside-product-scope).

#### `--sexchr <n>`

*Int, default 23.* Which chromosome code is X. It also moves Y, XY and MT, and therefore
**changes which codes count as autosomal** — which changes every relatedness number, because
the autosome set is the analysis set.

```
autosomes = 1 .. n-1,  plus n+2       (the XY/PAR code, pooled with the autosomes)
X  = n      Y  = n+1      XY = n+2      MT = n+3
every other numeric code is dropped
```

On the `sexchr` fileset (chr 1, 2, 23, 24, 25, 26):

```
$ open-king -b sexchr.bed --kinship                 # default 23
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs, 300 Y-chromosome SNPs, 50 mitochondrial SNPs
  PLINK maps loaded: 6000 SNPs

$ open-king -b sexchr.bed --kinship --sexchr 26
Non-human samples are analyzed, with 26 pairs of chromosomes
  Genotype data consist of 5950 autosome SNPs, 50 X-chromosome SNPs
  PLINK maps loaded: 6000 SNPs

$ open-king -b sexchr.bed --kinship --sexchr 2
Non-human samples are analyzed, with 2 pairs of chromosomes
  Genotype data consist of 2000 autosome SNPs, 2000 X-chromosome SNPs
  2000 other SNPs are removed.
  PLINK maps loaded: 4000 SNPs
```

Any value other than 23 prints the `Non-human samples are analyzed` line. Values below 2 are
fatal:

```
$ open-king -b sexchr.bed --kinship --sexchr 1

FATAL ERROR - 
Sex chromosome 1 out of range.
```

**`--sexchr` only re-partitions numeric codes.** The alphabetic spellings `X`, `Y`, `XY`,
`MT` are classified by name and are unaffected — the same fileset written with letters
instead of numbers gives the default partition at `--sexchr 22`, while the numeric version
is re-partitioned:

```
$ open-king -b alphachr.bed --kinship --sexchr 22    # chromosomes written X / Y / XY / MT
Non-human samples are analyzed, with 22 pairs of chromosomes
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs, 300 Y-chromosome SNPs, 50 mitochondrial SNPs
  PLINK maps loaded: 6000 SNPs

$ open-king -b sexchr.bed --kinship --sexchr 22      # same data, numeric codes
Non-human samples are analyzed, with 22 pairs of chromosomes
  Genotype data consist of 4300 autosome SNPs (including 300 XY SNPs), 1500 Y-chromosome SNPs, 150 mitochondrial SNPs
  50 other SNPs are removed.
  PLINK maps loaded: 5950 SNPs
```

Also note that bare `--sexchr` with no value sets it to **0**, which is fatal — see
[§6](#6-how-the-parser-behaves).

### Output

#### `--rplot`
#### `--pngplot`

*Switches, default off.* In the reference these shell out to `R CMD BATCH` with embedded
scripts. **Outside the minimal product scope**: open-king has no R dependency. A recognized
plotting request exits 1 before opening the input:

```
$ open-king -b multifam.bed --rplot

FATAL ERROR - 
open-king's minimal relatedness product does not implement: --rplot.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
```

`--pngplot` gives the same block with its own name in place of `--rplot`.

#### `--plink`

*Switch, default off.* **Outside the minimal product scope.** It is rejected through the same
pre-I/O scope gate:

```
$ open-king -b multifam.bed --plink

FATAL ERROR - 
open-king's minimal relatedness product does not implement: --plink.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
```

### Output Parameter

#### `--prefix <string>`

*String, default `king`.* Prepended to every output filename by plain concatenation. See
[§4](#--prefix-is-a-concatenation-not-a-stem).

The prefix must be writable — the loader probes it by creating and removing
`<prefix>$TMP$.ped` before it parses a single row, so a nonexistent directory fails early:

```
$ open-king -b multifam.bed --kinship --prefix nodir/x
Read in PLINK fam file multifam.fam...

FATAL ERROR - 
Cannot open nodir/x$TMP$.ped to write.
```

#### `--rpath <path>`

*String, default empty.* Path to the R installation for the plotting flags. Parsed and
echoed in the banner, then **rejected**: naming it is fatal, exit 1, before the input is
opened. See [§8](#8-accepted-compatibility-spellings-outside-product-scope).

---

## 6. How the parser behaves

The command line is not parsed by a normal argument library. It reproduces KING 2.3.2's
hand-written parser, quirks included, because parity is measured on the banner as well as on
the output files. Five rules cover everything.

### 1. Names match case-insensitively, by unique prefix

`--RELATED`, `--related`, `--Related` and `--rel` are the same option. Any unambiguous
prefix works:

```
$ open-king -b multifam.bed --rel
Options in effect:
	--related
```

An **ambiguous** prefix is an error, not a resolution — it is reported and then ignored:

```
$ open-king -b multifam.bed --r
WARNING - 
Problems encountered parsing command line:

Command line parameter --r is ambiguous
```

Ambiguous two-character prefixes include `by`, `ca`, `co`, `ib`, `ma`, `no`, `pc`, `pr`,
`rp` and `se`. A bare `--` is ambiguous too (it is the empty prefix). A name *longer* than
the option is undefined, so `--related=1` is not a way to pass a value:

```
$ open-king -b multifam.bed --related=1 --kinship
WARNING - 
Problems encountered parsing command line:

Command line parameter --related=1 is undefined
```

Neither warning aborts the run. Parsing always continues and the analysis proceeds — which
means a typo can silently change what you ran. **Check the banner echo.**

### 2. Switches TOGGLE, they do not set

Repeating a switch turns it back off. `--related --related` leaves `--related` **off**,
which is why the run below falls through to the "no analysis selected" notice:

```
$ open-king -b multifam.bed --related --related
Please specify one of the following 24 options: --related --kinship --autoQC --mtscore --risk --ibs --homog --ibdseg --mds --pca --cluster --build --bysample --bysnp --tdt --unrelated --duplicate --roh --grm --gdt --pc -- --pcgdt --
```

Three occurrences turn it back on. This matters when a wrapper script appends a flag that
your config already set.

### 3. A value is consumed only if it *looks* like one

| type | with a value token | with no value token |
| --- | --- | --- |
| switch | never consumes one | toggles |
| int | consumes iff the token matches `^[+-]?[0-9]+$` | **toggles between 0 and 1** |
| double | consumes iff the token starts with `[0-9.+-]` and any exponent is complete | left unchanged |
| string | consumes **whatever comes next**, option or not | left unchanged |

A rejected token is not silently dropped — it comes back as an ignored argument, and the
option keeps its no-value behaviour:

```
$ open-king -b multifam.bed --kinship --degree 3.0
              Inference Parameter : --degree [1], --noscreen [-1717986816],
WARNING - 
Problems encountered parsing command line:

Command line parameter 3.0 (#5) ignored
```

`--degree 3.0` therefore means `--degree 1`. So do `--degree 3x`, `--degree 0x10` and
`--degree .5`. Integers are decimal only; doubles accept C's `strtod` forms including hex
(`--callrateN 0x10` is 16.0).

The int toggle is why bare `--degree` (default 0) becomes **1** while bare `--sexchr`
(default 23) becomes **0** — and 0 is out of range:

```
$ open-king -b multifam.bed --kinship --sexchr

FATAL ERROR - 
Sex chromosome 0 out of range.
```

Strings are the dangerous case, because they take the next token unconditionally and issue
no warning:

```
$ open-king -b multifam.bed --trait --related --kinship
                Association Model : --trait [--related], --covariate [],
	--kinship
```

`--related` was eaten as the trait name and never turned on.

### 4. Everything else is "ignored"

A bare token, or a single-dash token other than `-b`/`-B`, is reported with its 1-based
argument position and skipped:

```
$ open-king -b multifam.bed --kinship -h foo
WARNING - 
Problems encountered parsing command line:

Command line parameter -h (#4) ignored
Command line parameter foo (#5) ignored
```

`-b` followed by a dash-leading token consumes nothing and **clears** the fileset, so a
later bare `-b` throws away an earlier one:

```
$ open-king -bmultifam.bed -b --kinship

FATAL ERROR - 
Genotype files are required. e.g.,
  king -b ex.bed --related
```

A lone `--` borrows the following non-dash token as an option name, and the borrowed option
never takes a value of its own:

```
$ open-king -b multifam.bed -- related
Options in effect:
	--related
```

### 5. The warning block starts with a BEL

The `WARNING - ` line is preceded by a `0x07` byte, faithfully to the reference. If you are
grepping stdout for `^WARNING` it will not match; strip control characters first
(`tr -d '\a'`).

---

## 7. Exit status and fatal errors

Two statuses. **0** for a completed run; **1** for a fatal error. Nothing is written to
stderr — fatal errors go to stdout like everything else, in a block of the form

```

FATAL ERROR - 
<message>

```

Measured, on the corpus:

| command | exit |
| --- | ---: |
| `open-king` | 1 |
| `open-king -b multifam.bed --kinship` | 0 |
| `open-king -b multifam.bed` (no analysis) | 0 |
| `open-king -b multifam.bed --pca` (excluded analysis) | 1 |
| `open-king --xyz -b multifam.bed --kinship` (bad option) | 0 |
| `open-king -b nosuch.bed --kinship` | 1 |
| `open-king -b multifam.bed --kinship --sexchr 1` | 1 |
| `open-king -b multifam.bed --kinship --maxP 0` (excluded parameter) | 1 |
| `open-king -b multifam.bed --risk` (excluded analysis) | 1 |
| `open-king -b multifam.bed --kinship --prefix no/such/` | 1 |

Two things a shell script must not infer from the exit status:

* **A bad option is not an error.** `--xyz` warns and the run proceeds, exit 0. Grep the
  warning block if you need to catch typos in CI.
* **An excluded product-scope request is an error.** `--pca`, its associated parameters,
  and every other entry in §8 exit 1 before the input is opened. This is a deliberate
  open-king safety contract, not reference-console parity.

The fatal messages, in the order they can fire:

| message | cause |
| --- | --- |
| `Genotype files are required. e.g., …` | no `-b` |
| `open-king's minimal relatedness product does not implement: …` | any option or input form in [§8](#8-accepted-compatibility-spellings-outside-product-scope) |
| `Sex chromosome n out of range.` | `--sexchr` below 2 |
| `Genotype file <path> cannot be opened` | `.bed` missing or unreadable |
| `Please use PLINK binary format as input.` | `-b` argument does not end in `.bed` |
| `Please use either PLINK or KING binary format as input.` | bad `.bed` magic |
| `Cannot open <prefix>$TMP$.ped to write.` | prefix directory missing or read-only |
| `Pedigree file <path> cannot be opened` | `.fam` missing |
| `Map file <path> cannot be opened` | `.bim` missing |
| `<path>: line N has K fields, expected 6` | short `.fam` or `.bim` line |
| `Please correct problems with pedigree structure` | duplicate `(FID, IID)` |
| `No autosome SNPs are available. Please check your map file.` | nothing survived the chromosome filter |
| `Currently only SNP-major mode can be analyzed.` | `.bed` mode byte is 0 |
| `Not enough genotypes at the Nth marker` | `.bed` shorter than the map requires |

**Two of the reference's fatal messages cannot fire here, and one of its quirks cannot be
observed.** All three sit after the product-scope gate, which rejects an excluded option on
its presence alone, whatever value follows it, so nothing reaches them:

* `p-value [x] outside range in ninv()`, the `0 < p < 2` check on `--maxP`.
* `Please use --model <file> to specify a risk model.`, for `--risk` without `--model`.
* the `--risk --model` path in the loader that skips the early `.bed` probe, which reports a
  missing fileset as `Pedigree file <name>.fam cannot be opened` rather than
  `Genotype file <name>.bed cannot be opened`.

The code is kept for reference fidelity. It is listed here so nobody hunts for a command line
that reaches it; there is none.

---

## 8. Accepted compatibility spellings outside product scope

Twenty-four options and one input form are outside open-king's deliberately minimal
relatedness and QC scope. They are parsed and echoed only so the banner and parse surface stay
byte-exact against the reference. They are not planned functionality and are not limitations
of the supported core; see [SCOPE.md](SCOPE.md).

Eleven analyses and output modes, plus the input form:

`--pca` · `--mds` · `--roh` · `--lmm` · `--tdt` · `--gdt` · `--risk` · `--makeGRM` ·
`--plink` · `--rplot` · `--pngplot` · multi-dataset input

Their associated parameters are rejected too, thirteen of them: `--projection`, `--pcs`,
`--trait`, `--covariate`, `--maxP`, `--invnorm`, `--model`, `--prevalence`, `--noflip`,
`--phefile`, `--covfile`, `--prunedsnp` and `--rpath`.

The list is `Options::unsupported_requests()` in `crates/open-king-cli/src/cli.rs`, and the two
groups above are its 24 entries. Nothing in [§5](#5-option-reference) is exempt from it: an
option described there as parsed and echoed is still rejected if it appears here.

The parser continues to recognize and echo these names for a familiar KING surface. After
the banner, open-king emits one fatal block naming every excluded request, exits 1, and does
not print `KING starts at`, probe the input, or write an output file. For example:

```text
FATAL ERROR -
open-king's minimal relatedness product does not implement: --pca.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
```

**Multi-dataset input** — the reference accepts `-b a.bed,b.bed` (and comma-separated
`--fam`/`--bim`) and merges the filesets. open-king rejects this form before file lookup:

```
$ open-king -b multifam.bed,dups.bed --kinship          # open-king, exit 1

FATAL ERROR -
open-king's minimal relatedness product does not implement: comma-separated multi-fileset input.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.

$ king -b multifam.bed,dups.bed --kinship          # KING 2.3.2
Read in PLINK bim files
	multifam.bim...
	dups.bim...
```

Merge the filesets with PLINK first.

---

## 9. Differences from the reference

[**PARITY.md**](PARITY.md) is the authoritative statement, with the measurement behind every
claim: all 480 captured reference invocations are byte-identical. Read its held-out section
before treating that finite corpus as universal proof. This section only lists what a *user*
driving the command line can hit.

**Byte-identical on every capture that produces them:** `--kinship` (including the X pass),
`--duplicate`, `--ibs`, `--unrelated`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster`,
`--ibdseg` and `--related` at all three captured `--seglength` floors, plus `X.kin`,
`X.kin0` and `X.seg`.

**There are no known gaps in the captured suite.** The former two-stage-screen console
differences and incomplete primary `build.log` are fixed; stdout, stderr, exit status, and
all compared files match in every case.

**Divergences outside the captured suite**, which the corpus cannot see — they need an input
shape the 480 captures do not contain. PARITY.md §4.6, §5.10, §5.11 and §5.12 enumerate the
known set; these are the ones a command line can reach:

* **The sparse-map fallback and screen are implemented.** On a panel too thin for the segment
  caller (roughly under 12 500 markers at 200 samples), `--related` now prints
  `No informative IBD segments.`, switches to the reference's 12-column `.kin` and infers
  relationships from kinship alone. Its held-out comparison is byte-exact within families
  and on all 15 between-family candidates. `--unrelated`, `--cluster` and `--build` take the same
  kinship-only clustering path; both unrelated lists, `updateids.txt`, `updateparents.txt`
  and `build.log` are byte-identical, and `--cluster` correctly omits its segment-only
  `cluster.kin`. Check for the `usable for IBD segment analysis` line before trusting any
  segment column
  ([PARITY.md §5.12](PARITY.md#512-three-divergences-found-while-writing-the-user-documentation)).
  The fallback's PO/FS split still uses a fixed `0.0050`; KING's printed cutoff is
  data-derived and its derivation remains a held-out open question.
* **Unsorted `.bim` maps are rejected for segment work.** Both binaries emit the same
  `Positions unsorted: …` / `Chromosomes unsorted: …` diagnostic and suppress the same
  segment outputs. See [§3](#two-hard-requirements-that-are-easy-to-miss).
* **Sample IDs colliding only in case are rejected.** Both binaries treat `(FID, IID)` as
  unique under ASCII case-folding, name the second spelling in the duplicate diagnostic and
  abort with `Please correct problems with pedigree structure`.
* **The A1-major QC check is implemented for the stable 4,096-marker window.** Its boundary,
  analysis surface, fatal placement and pre-fatal artifacts match KING. Shorter maps skip
  the check because KING's own tail read is nondeterministic
  ([§3](#two-hard-requirements-that-are-easy-to-miss)).
* **`--ibdseg` applies the closed 100 Mb usable-segment floor.** Below it both binaries
  print `Segments too short.` and suppress `.seg`; exactly 100,000,000 bp proceeds
  ([PARITY.md §5.10](PARITY.md)).

**Reference bugs that are reproduced on purpose**, because parity requires it: the
`--noscreen` uninitialised default, the `3nd-degree` typo, the truncated single-family
`.kin`, `--related`'s 100-sample floor on `.kin0`, and `X.seg`'s eleven-name header over
nine-value rows. None is a defect in this implementation, and none will be fixed while
byte-parity is the target.

Finally, one thing that is *not* a difference: KING's own segment numerics changed at
2.1.2, 2.1.3, 2.2.1, 2.2.5, 2.2.6 and 2.2.7. Parity here means parity with **2.3.2**. If you
are diffing against another KING build, expect the `.seg` columns to differ, and re-capture
your goldens before believing any number in PARITY.md.

---

## 10. The derived filesets used above

Twelve examples on this page need a fileset the corpus does not ship: an A1-major recode, two
unsorted maps, two relabelled chromosome sets, a thinned X array, and a `.bed` under a name
that does not end in `.bed`. This script builds every one of them from the corpus, using only
the Python standard library, so each of those commands is reproducible exactly as written.

```python
#!/usr/bin/env python3
# fixtures.py -- build docs/CLI.md's derived filesets from the corpus.
#   python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
#   python3 fixtures.py /tmp/kingdocs .
import shutil, sys
from pathlib import Path
C, OUT = Path(sys.argv[1]), Path(sys.argv[2]); OUT.mkdir(parents=True, exist_ok=True)

def read(stem):
    bim = [l.split() for l in (C / (stem + ".bim")).read_text().split("\n") if l.strip()]
    return bim, (C / (stem + ".fam")).read_text(), (C / (stem + ".bed")).read_bytes()

def write(stem, bim, fam, bed):
    (OUT / (stem + ".bim")).write_text("".join("\t".join(r) + "\n" for r in bim))
    (OUT / (stem + ".fam")).write_text(fam)
    (OUT / (stem + ".bed")).write_bytes(bed)

def nsamp(fam): return len([l for l in fam.split("\n") if l.strip()])
def groups(bim):
    g, i = [], 0
    while i < len(bim):
        j = i
        while j < len(bim) and bim[j][0] == bim[i][0]: j += 1
        g.append((i, j)); i = j
    return g

# plain -- a valid .bed whose filename does not end in .bed
bim, fam, bed = read("multifam")
write("plain_src", bim, fam, bed)
for a, b in (("plain_src.bed", "plain"), ("plain_src.bim", "plain.bim"),
             ("plain_src.fam", "plain.fam")): shutil.copy(OUT / a, OUT / b)

# minor / major -- the same genotypes with A1 as the minor, then the major, allele
write("minor", bim, fam, bed)
FLIP = {0b00: 0b11, 0b11: 0b00, 0b10: 0b10, 0b01: 0b01}
tbl = bytes(sum(FLIP[(b >> (2 * k)) & 3] << (2 * k) for k in range(4)) for b in range(256))
write("major", [r[:4] + [r[5], r[4]] for r in bim], fam, bed[:3] + bed[3:].translate(tbl))

# unsortedpos -- positions descending inside each chromosome
rows = [list(r) for r in bim]
for i, j in groups(rows):
    for k, pos in enumerate([r[3] for r in rows[i:j]][::-1]): rows[i + k][3] = pos
write("unsortedpos", rows, fam, bed)

# unsortedchr -- chromosomes 22 -> 1, positions ascending inside each
nb = (nsamp(fam) + 3) // 4
g = groups(bim)
write("unsortedchr", [r for a, b in reversed(g) for r in bim[a:b]], fam,
      bed[:3] + b"".join(bed[3 + a * nb: 3 + b * nb] for a, b in reversed(g)))

# alphachr / unknownchr -- the sexchr map relabelled two ways
bim, fam, bed = read("sexchr")
write("alphachr", [[{"23": "X", "24": "Y", "25": "XY", "26": "MT"}.get(r[0], r[0])] + r[1:]
                   for r in bim], fam, bed)
write("unknownchr", [[{"24": "0", "26": "chr1"}.get(r[0], r[0])] + r[1:] for r in bim], fam, bed)

# x400 -- sexchr cut from 1500 X markers to 400
nb, keep, seen = (nsamp(fam) + 3) // 4, [], 0
for i, r in enumerate(bim):
    if r[0] == "23":
        seen += 1
        if seen > 400: continue
    keep.append(i)
write("x400", [bim[i] for i in keep], fam,
      bed[:3] + b"".join(bed[3 + i * nb: 3 + (i + 1) * nb] for i in keep))

# geno.bed + alt.ped + alt.map -- multifam under non-default names, for --fam/--bim
bim, fam, bed = read("multifam")
(OUT / "geno.bed").write_bytes(bed)
(OUT / "alt.ped").write_text(fam)
(OUT / "alt.map").write_text("".join("\t".join(r) + "\n" for r in bim))
```

The `c1`/`c16` and `n1`/`n2` directories in [§5](#--cpus-n) and [§8](#8-accepted-compatibility-spellings-outside-product-scope)
are just two empty directories per comparison, each holding one run of the command quoted
beside them.

---

## See also

* [`README.md`](README.md) — the documentation index
* [`OUTPUTS.md`](OUTPUTS.md) — every output file: columns, formats, row order, existence rules
* [`COOKBOOK.md`](COOKBOOK.md) — task-oriented recipes, from finding duplicates to diffing
  against KING
* [`INTERPRETING.md`](INTERPRETING.md) — what the numbers mean, and where they mislead
* [`PARITY.md`](PARITY.md) — the authoritative statement of what is byte-identical to
  KING 2.3.2, measured per file and per row
