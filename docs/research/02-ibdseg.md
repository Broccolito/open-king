# KING 2.3.2 — IBD-segment side (`--ibdseg`, `--related`, `--degree`) — recon for clean-room reimplementation

**Status:** RECON ONLY. Sources used: (a) the KING website (manual / download / version
history pages), (b) the *behaviour* of the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2, Mach-O
arm64), (c) the binary's embedded literal strings (output headers, printf formats, and the
R script it emits — these are facts about the output format, not source code).

**Clean-room statement.** KING's C++ source is publicly downloadable
(`https://www.kingrelatedness.com/KINGcode.tar.gz`, and per-version
`KING2.3.2code.tar.gz` etc.). **It was NOT downloaded, opened, or read.** Nothing below is
transcribed from KING source. Every algorithmic statement is either (i) quoted from public
documentation, or (ii) *inferred from black-box experiments* that are fully reproducible
from the synthetic-data scripts listed in §11. Keep it that way: do not fetch KINGcode.tar.gz.

---

## 1. Executive summary of what must be reproduced

`king -b X.bed --ibdseg --prefix P` writes **exactly three** files in this build:

| file | always? | content |
|---|---|---|
| `P.seg` | yes | one row per reported pair: `IBD1Seg IBD2Seg PropIBD InfType` |
| `Pallsegs.txt` | yes | the "usable chromosomal segments" (the denominator definition) |
| `Psplitped.txt` | only when the pedigree needs splitting | pedigree bookkeeping for plots |

`P.segments.gz` (the per-segment detail file described in the manual) is **NOT produced by
this binary** — see §5.3. Plan for `.seg` byte-parity first; treat `.segments.gz` as a
"documented format we may choose to emit" rather than a parity target.

---

## 2. `--ibdseg` command-line semantics (verified against the binary)

### 2.1 Banner
Running `king` with no genotype file prints (verbatim, relevant lines):

```
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen
...
   Pairwise Relatedness Inference : --kinship, --ibdseg, --ibs, --makeGRM
              Inference Parameter : --degree, --noscreen [-1717986816],
                                    --seglength, --minConc [0.80]
...
Please check the reference paper Manichaikul et al. 2010 Bioinformatics,
					Chen et al. 2024,
          or the KING website at kingrelatedness.com
```

Note `--noscreen [-1717986816]` — that is an uninitialised-looking default printed for a
boolean flag; harmless, but if you mimic the banner, that literal is what 2.3.2 prints.
Note also the literal tab-indentation before `Chen et al. 2024,`.

### 2.2 `--degree <int>`
* **Default: absent.** With no `--degree`, `--ibdseg` reports **every** pair that survives
  the "long segment" filter (§5.4), including pairs inferred `UN`.
* With `--degree d`, only pairs with **`PropIBD` above the d-th-degree cutoff** are written.
  Manual: *"The second command specifies only pairs with IBD proportion > 0.0884 will be
  saved in the output"* (for `--degree 3`).
* Cutoff = `2^-(d+0.5)`  (i.e. twice the kinship cutoff `2^-(d+1.5)`):

  | `--degree` | PropIBD cutoff | exact value |
  |---|---|---|
  | 1 | 2^-1.5 | 0.3535534 |
  | 2 | 2^-2.5 | 0.1767767 |
  | 3 | 2^-3.5 | 0.08838835 |
  | 4 | 2^-4.5 | 0.04419417 |

  Verified empirically: `--degree 1` kept min PropIBD 0.4696 and dropped 0.3023;
  `--degree 2` kept 0.2197 and dropped 0.1048; `--degree 3` kept 0.0888 and dropped 0.0047.
  (The binary's own R output uses the literals `0.3535534 / 0.1767767 / 0.08838835 /
  0.04419` — see §7.)
* `--degree` is *only a filter on what is written*, not on how InfType is computed: a pair
  labelled `2nd` still prints `2nd` under `--degree 1` if it passes the cutoff.
* For `--related`, `--degree` means something different: it is the **depth of relatedness
  searched**, and it changes the algorithm (screening stages). Default for `--related` is
  degree 1 ("Note only duplicates and 1st-degree relatives are included in the inference. /
  Specifying '--degree 2' if a higher degree relationship inference is needed."). The manual
  says *"--related without the --degree option is highly recommended"* and that degrees > 2
  have "no fast algorithm … computation is substantially slower".
* `--degree` also applies to `--kinship`, `--unrelated`, `--build`, `--cluster`.
  Binary strings show `--degree 5` is accepted (used by `--cluster`/`--unrelated`).

### 2.3 `--seglength <int Mb>`
* Manual: *"--seglength specifies the minimum length of IBD segments that are considered
  towards the relationship inference."*
* Binary strings (exact):
  * `Minimum segment length is set as %d bp`
  * `KING supports minimum segment length from 1 to 10 Mb at the moment.`
  * `Default seglength of 3Mb is used.`
* **Default = 3 Mb.** Units of the argument are **Mb**; internally stored as bp.
* **Verified:** `--seglength 5` changes results (e.g. an FS pair moved
  `0.4989/0.2201/0.4696` → `0.5003/0.2187/0.4689`); `--seglength 0` and `--seglength 11`
  produce a file **byte-identical to the default** (md5-identical) ⇒ out-of-range values
  silently fall back to 3 Mb.
* `--seglength` is *not* echoed in the "Options in effect" block (only `--ibdseg`,
  `--degree`, `--prefix`, `--rplot` etc. are).

### 2.4 Other parameters that interact
* `--prefix P` (default `king`), `--cpus N` (default = half the logical cores),
  `--sexchr N` (default 23), `--rplot`, `--pngplot`, `--rpath`.
* `--projection N` — manual: *"--ibdseg --projection N infers IBD segments between two
  samples from two subsets where the first subset consists of N samples."*
* Mutually exclusive (binary strings): `Please do not run --ibdseg together with --autoQC` /
  `--homog` / `--pca`; and `--ibdseg is skipped.`
* `Cannot run --exact analysis without IBD segments.` — `--exact` is an undocumented
  add-on that consumes IBD segments (not present as a working flag in this build).

### 2.5 Console output of a successful `--ibdseg` run (verbatim template)

```
Total length of %d chromosomal segments usable for IBD segment analysis is %.1lf Mb.
  In addition to autosomes, %d segments of length %.1lf Mb on X-chr can be further used.
  Information of these chromosomal segments can be found in file %s

IBD segment analysis starts at %s
%d CPU cores are used for %s inference...            <- %s is "autosome" or "X-chr"
                       ends at %s

Note with relationship inference as the primary goal, the following filters are applied:
  Sample pairs without any long IBD segments (>10Mb) are excluded.
  Short IBD segments (<3Mb) are not reported/utilized.
Summary statistics of IBD segments for individual pairs saved in file %s
Additional summary statistics of X-Chr IBD segments saved in file %s
```
Other messages in the same code path: `No informative IBD segments.`,
`Segments too short.`, `Positions unsorted: %s at %d, %s at %d.`,
`Chromosomes unsorted: %s on chr %d, %s on chr %d.`,
`Too many first alleles as the major allele (~%.1lf%%). Please use plink1.9 --make-bed to
regenerate the genotype data again.`

---

## 3. `<prefix>.seg` — exact format (VERIFIED byte-for-byte)

* Separator: **single TAB** (`\t`). No trailing tab. Unix `\n` line endings. Trailing
  newline present on the last row.
* Header (exact):

```
FID1	ID1	FID2	ID2	IBD1Seg	IBD2Seg	PropIBD	InfType
```

* Row format: `%s\t%s\t%s\t%s\t%.4lf\t%.4lf\t%.4lf\t%s\n`
  (four `%s` IDs, then **three `%.4lf`**, then the InfType string).
* Real observed rows (from a synthetic pedigree run of the reference binary):

```
FAM1	GF	FAM1	FA	1.0000	0.0000	0.5000	PO
FAM1	KID1	FAM1	KID2	0.4002	0.3238	0.5239	FS
FAM1	KID1	FAM4	DUPKID1	0.0000	1.0000	1.0000	Dup/MZ
FAM1	KID1	FAM2	COUSIN	0.2097	0.0000	0.1048	3rd
FAM1	GF	FAM1	GM	0.0040	0.0000	0.0020	UN
```

* **Row ordering:** nested loop over sample index in **.fam / genotype order**,
  `for i in 0..N-1 { for j in i+1..N-1 }`; rows whose pair is filtered out are simply
  skipped. Within-family and between-family pairs are interleaved in that single file
  (unlike `--related`, which splits `.kin` / `.kin0`).
* **`InfType` value set (exact strings, VERIFIED in output):**
  `Dup/MZ`, `PO`, `FS`, `2nd`, `3rd`, `4th`, `UN`.
  ⚠ The **manual says `Dup/MZTwin`; the binary writes `Dup/MZ`**. Match the binary.
* Numbers are always 4 decimals, including `0.0000` and `1.0000`. No sign, no `NA`.

### 3.1 `<prefix>X.seg` (X-chromosome companion)
Not emitted in any of my runs (even with X SNPs present, `--sexchr 23`, `--degree`,
`--rplot`), but the writer exists. From the binary's string table the header is, in order:

```
FID1	ID1	FID2	ID2	Sex1	Sex2	IBD1Seg	IBD2Seg	PropIBD
```
with row format `%s\t%s\t%s\t%s\t%d\t%d\t%.4lf\t%.4lf\t%.4lf` (Sex1/Sex2 are `%d`).
A second, longer header also exists in the same region —
`FID1 ID1 FID2 ID2 Sex1 Sex2 MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD` — matching the
version-history note *"A minor bug in --ibdseg and --roh is now fixed (regarding maxIBD1 and
maxIBD2, not affecting the main inference)"*. Treat X output as **out of scope / low
priority**; there is no `InfType` column on the X file.

---

## 4. `<prefix>allsegs.txt` — the denominator definition (VERIFIED)

This file defines the genome over which IBD proportions are computed. Exact header and a
real row:

```
Segment	Chr	StartMB	StopMB	Length	N_SNP	StartSNP	StopSNP
1	1	0.100	249.000	248.900	2490	rs1_0	rs1_2489
```

* Tab separated. Formats: `Segment %d`, `Chr %d`, `StartMB %.3lf`, `StopMB %.3lf`,
  `Length %.3lf`, `N_SNP %d`, `StartSNP %s`, `StopSNP %s`.
* `Length == StopMB - StartMB` exactly (verified: `249.000 - 0.100 = 248.900`).
  Positions are bp/1e6.
* `Segment` is a 1-based running index across the whole genome (not per chromosome).
* Column order here (`… Length N_SNP StartSNP StopSNP`) is **different** from the
  `.segments.gz` order (`… StartSNP StopSNP N_SNP Length`). Do not mix them up.
* X-chromosome usable segments are appended to the same file (observed: a `23` row).
* The console line `Total length of %d chromosomal segments usable for IBD segment analysis
  is %.1lf Mb.` prints the count of rows and the sum of `Length` (autosomes only).

### 4.1 How chromosomes are cut into "usable segments" (black-box findings)

1. **Cut at large inter-marker gaps.** A chromosome is split wherever the distance between
   consecutive markers is **> 1,000,000 bp**. Verified to the bp:
   gap of exactly `1,000,000` → *no* split; `1,000,001` → split. (This is what makes real
   data yield ~34–40 segments instead of 22: centromeres and assembly gaps.)
2. **Discard pieces that are too small.** A resulting piece is dropped entirely (it does not
   appear in `allsegs.txt` and does not count toward the denominator). Empirical brackets
   from a controlled island sweep:
   * span `10.0` and `10.1` Mb → dropped; span `10.5` Mb → kept (with 2000 markers).
     ⇒ minimum span is in **(10.1, 10.5] Mb**; the "long IBD segment > 10 Mb" language
     suggests the intent is >10 Mb.
   * marker count matters too and behaves in a **quantised, alignment-dependent** way
     (at a fixed 40 Mb span: 256→dropped, 320→dropped, 352→dropped, 353→**kept**,
     360→dropped, 368→dropped, 370→kept, 371→kept, 384→kept, 390→kept, 400→kept).
     The kept/dropped pattern tracks the number of **complete 64-marker words** contained
     in the piece (≥5 complete 64-SNP words ⇒ kept) — consistent with KING packing
     genotypes into 64-marker words and only using whole words. One case (11 Mb / 352
     markers / 5 words) was still dropped, so span and word-count criteria are combined.
   * **Action for the reimplementation:** replicate (1) exactly; for (2) start with
     "span > 10 Mb AND ≥ 5 complete 64-marker words", and re-probe against the binary
     with real `.bim` maps before claiming parity. This directly scales every IBD
     proportion, so it is a first-order parity item.

---

## 5. IBD segment detection — what is publicly known

### 5.1 There is no published algorithm description
The KING manual (last updated **July 28, 2023**) states, verbatim, in the *IBD SEGMENT
INFERENCE* section:

> "IBD (identical by descent) segments can be raplidly and accurately inferred between any
> pair of individuals in KING. **The associated manuscript is yet to be published** but this
> algorithm has been well tested."

So: **the exact IBD-calling algorithm is undocumented.** There is no paper describing the
sliding-window/run logic, the error tolerance, or the IBD1/IBD2 boundary rules. Everything
in §5.2 is inference from the binary's own messages plus the general (unphased,
IBS-based) approach that KING advertises.

### 5.2 What the binary tells us about the calling rules
* Only three numeric knobs are exposed: minimum segment length (`--seglength`, default 3 Mb,
  range 1–10 Mb), the fixed 10 Mb "long segment" pair filter, and `--degree`.
* Both filters are announced verbatim on every run:
  * `Sample pairs without any long IBD segments (>10Mb) are excluded.`
  * `Short IBD segments (<3Mb) are not reported/utilized.` (the `3` here tracks
    `--seglength`, i.e. it is the same parameter).
* KING works from **unphased** genotypes and needs **sorted** maps (it errors with
  `Positions unsorted:` / `Chromosomes unsorted:`), and it wants A1=minor
  (`Too many first alleles as the major allele … use plink1.9 --make-bed`).
* Mechanism (inferred, consistent with all observations): within a usable chromosomal
  segment, IBD2 stretches are runs with **no IBS0 and no "opposite-homozygote/het
  mismatch"** signal, IBD1 stretches are runs with **no IBS0** (an IBD1 region cannot
  produce IBS0 because one haplotype is shared), with some error tolerance and a
  word-granular (64-marker) scan. My synthetic experiments are consistent with a
  no-IBS0-run detector with tolerance: a truly IBD1 block of *L* Mb is recovered as
  ~*L*+1–2 Mb (boundaries extend to the next IBS0), and unrelated pairs in
  linkage-equilibrium data pick up spurious 10–15 Mb "segments" at a rate consistent with
  runs of ~100 markers without IBS0. **Do not treat that as the algorithm — treat it as
  the acceptance test.**
* IBD2 is nested inside IBD1 in the *reporting* sense: the two proportions are reported
  separately and `IBD1Seg + IBD2Seg <= 1` in every observed row (max observed sum = 1.0000).

### 5.3 `.segments.gz` — documented but NOT produced by this binary
The manual documents a per-segment file `ex.segments.gz` ("tar zipped"). **This build never
writes it**: verified across `--ibdseg`, `--ibdseg --degree 2/3/4`, `--ibdseg --seglength
1/5/10`, `--ibdseg --rplot`, `--ibdseg --sexchr 23`, and `--roh` (which likewise produced
`.roh` but no `.rohseg.gz`). The binary *is* linked against `/usr/lib/libz.1.dylib`, and its
string table contains `--ibdall cannot run without ZLIB` / `--ibdGRM cannot run without
ZLIB` (those flags are not accepted by the released build — passing them silently falls
back to `--related`). Conclusion: the per-segment writer is compiled out / gated in 2.3.2.

For completeness, the documented format (manual, verbatim header + rows):

```
FID1    ID1     FID2    ID2     IBDType Chr     StartMB StopMB  StartSNP        StopSNP         N_SNP   Length
1330    NA12335 1330    NA12340 IBD1    1       51.799  95.862  rs7534689       rs1858111       294     44.1
1330    NA12335 1330    NA12340 IBD1    1       148.175 247.083 rs1868992       rs12058711      692     98.9
1330    NA12335 1330    NA12340 IBD1    2       0.143   88.714  rs408209        rs7581608       619     88.6
```

Column semantics (manual, verbatim):
`IBDType`: "Type of IBD segments: IBD1 or IBD2" (values `IBD1` / `IBD2`);
`Chr`: chromosome number; `StartMB`/`StopMB`: start/stop position in Mb;
`StartSNP`/`StopSNP`: start/stop SNP; `N_SNP`: number of SNPs in the segment;
`Length`: "Total Length of the IBD segment (in Mb)".
Observed precision in the manual's example: `StartMB`/`StopMB` **3 decimals**, `Length`
**1 decimal**, `N_SNP` integer. (`Pallsegs.txt` uses 3 decimals for `Length`, so the two
files differ — the `.1lf` for `Length` is the manual's rendering of the gz file.)
There is one row per segment per pair, grouped by pair then chromosome then position.

### 5.4 Pair-level filter
A pair is written to `.seg` **iff at least one detected IBD segment is longer than 10 Mb**
(message is explicit; `--seglength` does not change this 10 Mb constant — runs with
`--seglength 1/5/10` produced identical pair sets). Pairs failing it are absent entirely
(not written with zeros). Verified: unrelated pairs whose total IBD1 was ~11 Mb (0.0040)
*were* written, while others were omitted, exactly tracking "some single segment > 10 Mb".

---

## 6. `IBD1Seg`, `IBD2Seg`, `PropIBD` — definitions (manual, verbatim)

> **IBD1Seg**: Total length of IBD1 segments divided by total length of all segments,
> estimate of π1=Pr(IBD=1)
> **IBD2Seg**: Total length of IBD2 segments divided by total length of all segments,
> estimate of π2=Pr(IBD=2)
> **PropIBD**: Proportion of genomes shared identical-by-descent, estimated by
> IBD2Seg + IBD1Seg/2, estimate of π=π2+π1/2

Operationally, verified:

```
denominator D = sum of `Length` over all rows of <prefix>allsegs.txt   (autosomes)
IBD1Seg = (total Mb called IBD1) / D
IBD2Seg = (total Mb called IBD2) / D
PropIBD = IBD2Seg + IBD1Seg / 2
```

* `PropIBD` is computed from the **full-precision** π1/π2, not from the printed 4-dp values,
  then formatted at 4 dp. Evidence: the row `0.2097 0.0000 0.1048` — half of the *printed*
  0.2097 is 0.10485, which would print as `0.1049` (or `0.1048`) depending on the rounding
  mode, whereas the underlying double was just below 0.10485. Compute in f64 and format
  once with `%.4f`; never round intermediates.
* `PropIBD ≈ 2 × kinship`. KING's own R code plots `Kinship` vs `PropIBD` with a reference
  line at `0.08839` on both axes.
* Denominator is the **autosomal** usable length; X segments are reported separately.
* Sanity anchors from a controlled construction (pair sharing exactly a fraction f of the
  cumulative genome): target π1=1.00 → `1.0000`; target 0.50 → `0.5002`; target 0.20 →
  `0.2021`. Reconstructed lengths agree with `D = 2875.8 Mb` printed by the binary.

---

## 7. `InfType` — the exact decision rule (EMPIRICALLY MAPPED)

KING's own generated R script (`<prefix>_ibd1vsibd2.R`, written by `--ibdseg --rplot`)
contains the classification it uses for plotting. Reproduced verbatim from a run:

```r
d0 <- data$IBD2Seg>0.7
d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)
d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08
d2 <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)
d3 <- data$PropIBD>0.08839 & data$PropIBD<=0.17678
d4 <- data$PropIBD>0.04419 & data$PropIBD<=0.08839
dU <- data$PropIBD>0 & data$PropIBD<=0.04419
...
abline(h = 0.08,       col = "green",   lty = 3, lwd = 2)
abline(a = 0.96,       b = -1,   col = "red",     lty = 3, lwd = 2)
abline(a = 0.3535534,  b = -0.5, col = "green",   lty = 3, lwd = 2)
abline(a = 0.1767767,  b = -0.5, col = "blue",    lty = 3, lwd = 2)
abline(a = 0.08838835, b = -0.5, col = "magenta", lty = 3, lwd = 2)
abline(a = 0.04419,    b = -0.5, col = "gold",    lty = 3, lwd = 2)
```

I then **verified the C++ `InfType` writer matches this rule** by constructing 225 synthetic
pairs with prescribed (π1, π2) and reading back the labels. Boundaries observed:

| transition | last value before | first value after | implied threshold |
|---|---|---|---|
| `UN`→`4th` | PropIBD 0.0411 | 0.0529 | 0.04419417 (2^-4.5) |
| `4th`→`3rd` | 0.0810 | 0.0907 | 0.08838835 (2^-3.5) |
| `3rd`→`2nd` | 0.1701 | 0.1803 | 0.1767767 (2^-2.5) |
| `2nd`→`FS`  | π2 = 0.1022 (PropIBD 0.3521) | π2 = 0.1122 (PropIBD 0.3617) | PropIBD > 0.3535534 (2^-1.5) |
| `PO`→`FS`   | π2 = 0.0723 | π2 = 0.0822 | π2 ≥ 0.08 |
| `2nd`→`PO`  | π1+π2 = 0.8962 | π1+π2 = 0.9004 | sum > 0.9 (with π2 < 0.08) |
| `FS`→`PO`   | sum = 0.9513 | sum = 0.9605 | sum > 0.96 |
| `FS`→`Dup/MZ` | π2 = 0.6924 | π2 = 0.7023 | π2 > 0.7 |

**Canonical rule (first match wins), with π1 = IBD1Seg, π2 = IBD2Seg, π = PropIBD:**

```
if  π2 > 0.7                                        -> "Dup/MZ"
elif (π1 + π2) > 0.96  or ((π1 + π2) > 0.9 and π2 < 0.08)  -> "PO"
elif π > 0.3535534 and π2 >= 0.08                   -> "FS"
elif π > 0.1767767                                  -> "2nd"
elif π > 0.08838835                                 -> "3rd"
elif π > 0.04419417                                 -> "4th"
else                                                -> "UN"
```

Notes that bite:
* **The `2nd` bucket is not bounded above.** A pair with π = 0.4481 but π1+π2 = 0.8962 and
  π2 = 0 is labelled `2nd`, not FS/PO. This is a real, reproducible behaviour (an entire
  arm of the π2=0 axis from π ≈ 0.18 to ≈ 0.45 is `2nd`).
* `FS` requires **π2 ≥ 0.08** (`>=`, not `>`), and `PO` requires **π2 < 0.08** in its second
  clause. Together with the π1+π2 clauses this makes the PO/FS split a two-piece boundary:
  a pair with π2 ≥ 0.08 is still `PO` if π1+π2 > 0.96.
* Thresholds are exact powers of two: 2^-1.5, 2^-2.5, 2^-3.5, 2^-4.5. KING's R writes
  `0.04419` for the last one; the C++ boundary is consistent with 0.04419417 (my sweep can
  only bracket it to (0.0411, 0.0529)). Use `2^-4.5`.
* Operators: use `>` for all PropIBD cuts and `>=` for the π2 ≥ 0.08 cut, as above.

---

## 8. `--related` — the integrated path (IBD columns only)

`--related` fuses kinship + IBD-segment estimates. It writes `<prefix>.kin` (within-family)
and `<prefix>.kin0` (across-family), plus `<prefix>X.kin` / `<prefix>X.kin0` when X data
exist. VERIFIED headers from the reference binary (tab-separated):

`<prefix>.kin`:
```
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType	Error
FAM1	FA	GF	28780	0.000	0.2500	0.1967	0.0000	0.3353	0.0000	0.2511	1.0000	0.0000	0.5000	PO	0
```
Formats: `N_SNP %d`, `Z0 %.3lf`, then `%.4lf` for Phi, HetHet, IBS0, HetConc, HomIBS0,
Kinship, IBD1Seg, IBD2Seg, PropIBD; `InfType %s`; `Error` printed as `%G` (`0`, `0.5`, `1`).
(The manual's column list omits `HomIBS0`; the actual header has it — trust the binary.)

`<prefix>X.kin` (VERIFIED, this one *is* emitted):
```
FID	ID1	ID2	Sex1	Sex2	PhiX	IBD1Seg	IBD2Seg	PropIBD
FAM1	FA	GF	1	1	0.0000	0.0000	1.0000	1.0000
```
(`Sex1`/`Sex2` `%d`; four `%.4lf`; **no InfType** column.)

`<prefix>.kin0` header (from strings; not triggered in my runs):
`FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 HetConc HomIBS0 Kinship` — plus IBD columns and
`InfType` in the integrated path (the manual's `--related` example shows the same IBD1Seg /
IBD2Seg / PropIBD / InfType tail). Console strings for the between-family stage:
```
  Stages 1&2 (with %d SNPs): %lli pairs of relatives are detected (with kinship > %.4lf)
  Final Stage (with %d SNPs): %lli pairs of relatives (up to %d%s-degree) are confirmed
Between-family relatives (kinship >= %.5lf) saved in file %s
```
The `--related` summary table printed to stdout is:
```
Relationship summary (total relatives: %d by pedigree, %d by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	8	1	4	0	6
  Inference	0	9	1	4	0	5
```

---

## 9. "Chen et al. 2024" — unresolved, and probably unfindable

The 2.3.2 binary asks users to cite *"Manichaikul et al. 2010 Bioinformatics, Chen et al.
2024"*. I could not identify that reference:

* Not on kingrelatedness.com. The manual's REFERENCE section (last updated 2023-07-28) lists
  **only** Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) *Robust
  relationship inference in genome-wide association studies*, **Bioinformatics 26(22):
  2867-2873**. The Download page's Reference anchor is the same.
* Not on Wei-Min Chen's own publication list (chen.kingrelatedness.com/publications) — his
  2024 entries are T1D/T2D GWAS papers, no methods paper.
* PubMed (`Chen Wei-Min[Author] AND 2024[dp]`, 19 hits) contains no relatedness/IBD methods
  paper. Crossref author+topic search likewise finds nothing.
* The manual explicitly says the IBD-segment manuscript is *"yet to be published"*.

**Best interpretation:** "Chen et al. 2024" is a forward reference to the (then-in-preparation)
KING IBD-segment / biobank-scale relatedness paper. Note the binary is version 2.3.2 dated
**September 8, 2023** on the website, yet cites a 2024 work — so the citation was
anticipatory, or this Mach-O arm64 build is a later rebuild of the 2.3.2 tree. Either way,
**there is no published algorithm description to implement against.** For our clean-room
purposes this is good news (nothing to copy) and bad news (behavioural parity is the only
spec).

---

## 10. Version history relevant to `--ibdseg` (kingrelatedness.com/history.shtml)

| version | date | ibdseg-relevant change |
|---|---|---|
| 2.1 | 2017-10-24 | `--ibdseg` introduced ("new functions added: --ibdseg for pair-wise IBD segment inference") |
| 2.1.2 | 2017-12-14 | "The IBD segment algorithm is improved" |
| 2.1.3 | 2018-02-13 | `--ibdseg`, `--related`, `--roh` improved |
| 2.1.4 | 2018-06-06 | "KING can now accurately infer up to the **4th-degree** relatedness (--related, --ibdseg), while the original KING method (--kinship) remains accurate up to the second-degree" |
| 2.1.5 | 2018-08-24 | `--ibdseg` up to 700,000 samples |
| 2.1.8/2.2 | 2019-02/03 | bug in `--ibdseg --degree 2` fixed |
| 2.2.1 | 2019-05-14 | "minor bug in --ibdseg and --roh … regarding maxIBD1 and maxIBD2, not affecting the main inference" |
| 2.2.3 | 2019-08-09 | "IBD segment analysis and Run of Homozygosity analysis now apply to **chromosome X** as well" |
| 2.2.5 | 2020-06-05 | "`--ibdseg` is substantially improved" |
| 2.2.6 / 2.2.7 | 2021 | two follow-up `--ibdseg` bug fixes |
| 2.3.2 | 2023-09-08 | `--related --degree 3+` crash fixed; `--noscreen` added |

Accuracy claim to hold ourselves to: **`--ibdseg` is accurate to 3rd degree on arrays and
4th degree on WGS** (manual: "with accuracy up to 3rd- or 4th-degree (depending on array or
WGS) for --related and --ibdseg analyses, and up to 2nd-degree for --kinship analysis").

---

## 11. Reproduction assets (all synthetic, no real genomes)

Generators + KING outputs live in
`/private/tmp/claude-501/.../74b7491e-.../scratchpad/`:

| file | purpose |
|---|---|
| `gen.py` | 13-person synthetic pedigree (PO/FS/2nd/3rd/dup/unrelated) → `run/synth.*` |
| `genx.py` | same + chromosome 23, for X behaviour |
| `sweep.py` | **225 pairs with prescribed (π1, π2)** — the InfType boundary map |
| `seglen.py` | pairs sharing a single block of prescribed length — filter probing |
| `run/*.seg`, `run/*allsegs.txt` | reference-binary outputs, byte-exact |
| `king_strings.txt`, `king_strings2.txt` | `strings` dumps of the binary (headers/formats) |
| `manual.txt`, `pg_history.shtml.txt` | text renderings of the website pages |

Regenerate the boundary map with:
```
python3 sweep.py sw && king -b sw.bed --ibdseg --prefix sw
```

---

## 12. Parity checklist for the Rust implementation

1. `.seg` header string, tab separators, `%.4f` on three columns, `\n` endings. ✅ specified
2. `InfType` = `Dup/MZ` (not `Dup/MZTwin`), `PO`, `FS`, `2nd`, `3rd`, `4th`, `UN`. ✅
3. Decision rule of §7, first-match-wins, with `>=` on the π2 = 0.08 test. ✅
4. `PropIBD = IBD2Seg + IBD1Seg/2` in f64, printed at 4dp. ✅
5. `--degree d` ⇒ keep rows with `PropIBD > 2^-(d+0.5)`; no `--degree` ⇒ keep all. ✅
6. Pair emitted only if some segment > 10 Mb; segments shorter than `--seglength`
   (default 3 Mb, clamp-to-default outside 1..10 Mb) discarded. ✅
7. Denominator = Σ Length over `allsegs.txt`; cut chromosomes at inter-marker gaps
   > 1,000,000 bp; drop tiny pieces (span > 10 Mb + ≥5 complete 64-marker words — **needs
   one more probing round against real maps**). ⚠
8. Row order = `.fam` order, `i<j`, within- and between-family interleaved. ✅
9. `allsegs.txt` header/format (`%.3lf` ×3, `%d` ×3, `%s` ×2), 1-based global `Segment`. ✅
10. Do **not** attempt `.segments.gz` parity — the reference binary never writes it. ✅
11. Open: exact IBD1/IBD2 calling (error tolerance, boundary placement). This is the only
    genuinely unspecified part; drive it by fitting against the reference binary on
    synthetic pairs, never by reading KING source. ⚠
