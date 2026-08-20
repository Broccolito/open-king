# open-king cookbook

Task-oriented recipes for `open-king`. Each one states the goal, gives the exact command, shows
the real output, and explains how to read it.

**Every command on this page was run, and every output block is pasted from that run.**
Console output is long, so blocks are often trimmed; each trim is marked with a
`[... elided ...]` line. Nothing is paraphrased, reformatted, or invented.

**The numbers are from a synthetic demo corpus.** They are reproducible — you can generate
the exact same filesets and get the exact same digits — but they are properties of *that
data*, not universal constants. Each recipe names the dataset it uses. Where a result
depends on the data, the text says so.

Reference material this page deliberately does not duplicate:

| For | Read |
| --- | --- |
| every command-line option, in full | [CLI.md](CLI.md) |
| every output file and every column | [OUTPUTS.md](OUTPUTS.md) |
| what the numbers mean, and how to misread them | [INTERPRETING.md](INTERPRETING.md) |
| the estimator formulas, verified against the reference | [VERIFIED_FORMULAS.md](VERIFIED_FORMULAS.md) |
| which SNPs are used, when each file is written, sort orders | [BEHAVIOR.md](BEHAVIOR.md) |
| the implementation specification | [SPEC.md](SPEC.md) |
| how close to the original KING this is, gap by gap | [PARITY.md](PARITY.md) |
| contributing, re-capturing goldens | [MAINTAINING.md](MAINTAINING.md) |

---

## Contents

* [Before you start](#before-you-start)
* [Which command do I want?](#which-command-do-i-want)
* [1. Find duplicate samples and MZ twins](#1-find-duplicate-samples-and-mz-twins)
* [2. Find cryptic relatives in a supposedly unrelated cohort](#2-find-cryptic-relatives-in-a-supposedly-unrelated-cohort)
* [3. Verify a stated pedigree](#3-verify-a-stated-pedigree)
* [4. Produce a maximal unrelated subset for GWAS](#4-produce-a-maximal-unrelated-subset-for-gwas)
* [5. Relatedness in an admixed or structured cohort](#5-relatedness-in-an-admixed-or-structured-cohort)
* [6. Tell parent–offspring from full siblings](#6-tell-parentoffspring-from-full-siblings)
* [7. IBD segments, and what they can and cannot separate](#7-ibd-segments-and-what-they-can-and-cannot-separate)
* [8. Per-sample QC, per-SNP QC, and auto-QC](#8-per-sample-qc-per-snp-qc-and-auto-qc)
* [9. Reconstruct pedigrees from genotypes](#9-reconstruct-pedigrees-from-genotypes)
* [10. X-chromosome relatedness](#10-x-chromosome-relatedness)
* [11. Large cohorts: `--degree`, `--cpus`, and what actually costs time](#11-large-cohorts---degree---cpus-and-what-actually-costs-time)
* [12. Migrating from the original KING](#12-migrating-from-the-original-king)
* [Traps](#traps)

---

## Before you start

### Build

```
$ cargo build --release
   Compiling open-king-io v0.1.0 (/Users/wgu/Desktop/open-king/crates/open-king-io)
   Compiling rayon v1.12.0
   Compiling open-king-core v0.1.0 (/Users/wgu/Desktop/open-king/crates/open-king-core)
   Compiling open-king-cli v0.1.0 (/Users/wgu/Desktop/open-king/crates/open-king-cli)
    Finished `release` profile [optimized] target(s) in 19.18s
real 19.27
```

A clean build from an empty target directory, on an Apple-silicon Mac. No toolchain beyond
Rust is needed. The binary lands at `target/release/open-king`; the recipes below assume it is on
your `PATH` as `open-king`.

### Demo data

Every recipe uses the project's own test corpus, so you can reproduce all of it:

```
$ python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
  trio            3 samples    5000 SNPs   3 chrom      2 related pairs
  nuclear         6 samples   10000 SNPs   5 chrom     14 related pairs
  threegen       12 samples   20000 SNPs  22 chrom     39 related pairs
  multifam       20 samples   15000 SNPs  22 chrom    104 related pairs
  dups           10 samples   10000 SNPs  22 chrom      3 related pairs
  missing         6 samples   10000 SNPs   5 chrom     14 related pairs
  monomorphic    12 samples    5000 SNPs   2 chrom     14 related pairs
  sexchr         10 samples    6000 SNPs   6 chrom     14 related pairs
  unrelated      30 samples   20000 SNPs  22 chrom      0 related pairs
  admixed        40 samples   20000 SNPs  22 chrom     15 related pairs
  singleton       1 samples    5000 SNPs  22 chrom      0 related pairs
  pair            2 samples    5000 SNPs  22 chrom      1 related pairs
  bigish        200 samples   50000 SNPs  22 chrom    496 related pairs
wrote /tmp/kingdocs/MANIFEST.json (13 datasets)
```

The generator is deterministic (fixed seed), writes ordinary PLINK1 `.bed`/`.bim`/`.fam`
triples, and also writes `MANIFEST.json` — the *true* relationship of every pair, which two
recipes below use as ground truth. The datasets used here:

| Dataset | What it is |
| --- | --- |
| `dups` | an exact duplicate pair across two FIDs, an MZ-like pair with 0.2 % genotype error, one parent–offspring pair, four unrelated founders |
| `multifam` | four declared nuclear families **plus undeclared cross-family relatives** (PO, FS, avuncular, first cousins) |
| `threegen` | one three-generation family: PO, FS, grandparent, avuncular, half-sib, half-avuncular and first-cousin pairs |
| `unrelated` | 30 mutually unrelated founders |
| `admixed` | two populations at F<sub>ST</sub> 0.10, six admixed founders, and three nuclear families — one within each population and one across them |
| `missing` | a nuclear family with per-sample missingness of 0/1/5/20/50 %, 300 high-missingness SNPs, 5 SNPs missing in everyone |
| `sexchr` | autosomes plus X, Y, XY and MT, both sexes, and two samples coded sex 0 |
| `bigish` | 200 people: nuclear families, three-generation units, undeclared cross-family sibships, unrelated singletons |

### Output naming

`--prefix` is a **plain string concatenation**, not a stem plus a separator. The default is
`king`, which is why files below are called `king.kin` and `kingallsegs.txt`:

```
$ open-king -b /tmp/kingdocs/dups.bed --kinship --prefix pfx
Within-family kinship data saved in file pfx.kin
Between-family kinship data saved in file pfx.kin0

$ ls pfx*
pfx.kin
pfx.kin0
```

`--prefix out/run_` gives `out/run_.kin`; `--prefix out/run` gives `out/run.kin`. The
directory must already exist.

---

## Which command do I want?

| You want | Command | Writes |
| --- | --- | --- |
| duplicate / MZ pairs only | `--duplicate` | `.con` |
| kinship for every pair | `--kinship` | `.kin`, `.kin0` (+ `X.kin`, `X.kin0`) |
| kinship **and** IBD segments and a relationship label | `--related` | `.kin`, `.kin0`, `allsegs.txt` (+ `X.kin`) |
| IBD-segment sharing only | `--ibdseg` | `.seg`, `allsegs.txt`, `splitped.txt` (+ `X.seg`) |
| full IBS / concordance statistics | `--ibs` | `.ibs`, `.ibs0`, `allsegs.txt` |
| a maximal unrelated subset | `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt` |
| merge families connected by genotype | `--cluster` | `updateids.txt`, `cluster.kin` |
| reconstruct parents from genotypes | `--build` | `updateids.txt`, `updateparents.txt`, `build.log` |
| sample-level QC | `--bysample` | `bySample.txt` |
| SNP-level QC | `--bySNP` | `bySNP.txt` |
| an automatic QC pass | `--autoQC` | four `_autoQC_*.txt` files |

`--related` is **not** a synonym for `--kinship`: it adds six columns, four of which come
from the IBD-segment engine. It also has sample-size gates that `--kinship` does not — see
[Traps](#traps) before you use it on a small cohort.

Kinship bands, used by every relationship label in every file:

| Class | Kinship |
| --- | --- |
| Duplicate / MZ twin | > 0.354 |
| 1st degree | 0.177 – 0.354 |
| 2nd degree | 0.0884 – 0.177 |
| 3rd degree | 0.0442 – 0.0884 |
| 4th degree | 0.0221 – 0.0442 |
| Unrelated | < 0.0221 |

---

## 1. Find duplicate samples and MZ twins

**Goal:** before anything else, find samples that are the same person twice — a plate
mix-up, a re-hybridised sample, or a genuine MZ twin pair.

```
$ open-king -b /tmp/kingdocs/dups.bed --duplicate
[... 27-line startup banner and the genotype-loading block elided ...]
Options in effect:
	--duplicate

Sorting autosomes...
Computing pairwise genotype concordance starts at Fri Aug 14 22:49:25 2026
  16 CPU cores are used...
        Stage 2 (with all SNPs) inference ends at Fri Aug 14 22:49:25 2026
2 pairs of duplicates with heterozygote concordance rate > 80% are saved in file king.con

  43 additional pairs from screening stage not confirmed in the final stage

KING ends at Fri Aug 14 22:49:25 2026
```

```
$ cat king.con
FID1	ID1	FID2	ID2	N	N_IBS0	N_IBS1	N_IBS2	Concord	HomConc	HetConc
DUPA	DUP_A	DUPB	DUP_A_COPY	10000	0	0	10000	1.00000	1.00000	1.00000
MZFAM	MZ_1	MZFAM	MZ_2	10000	9	17	9974	0.99740	0.99862	0.99512
```

**How to read it.** The test is `HetConc`, the heterozygote concordance rate, and the
default cutoff is 0.80.

* `DUP_A` / `DUP_A_COPY` — `N_IBS0 = 0`, `N_IBS1 = 0`, everything IBS2, all three
  concordances exactly `1.00000`. This is the same genotype vector twice: a duplicated
  sample, and one of the two must go.
* `MZ_1` / `MZ_2` — 9 opposite-homozygote sites and 17 half-mismatches out of 10 000.
  Identical DNA plus genotyping error (0.2 % per call, by construction). You cannot tell an
  MZ twin pair from a sample duplication genetically; the `.con` file will not do it for
  you, and neither will anything else in KING. Decide from the sample manifest.

The "43 additional pairs from screening stage not confirmed" line is the two-stage screen
reporting how many candidates the fast pass raised and the exhaustive pass then rejected.
It is informational.

**Tightening or loosening the cutoff** with `--minConc`:

```
$ open-king -b /tmp/kingdocs/dups.bed --duplicate --minConc 0.30 --prefix loose_
3 pairs of duplicates with heterozygote concordance rate > 30% are saved in file loose_.con

$ cat loose_.con
FID1	ID1	FID2	ID2	N	N_IBS0	N_IBS1	N_IBS2	Concord	HomConc	HetConc
DUPA	DUP_A	DUPB	DUP_A_COPY	10000	0	0	10000	1.00000	1.00000	1.00000
MZFAM	MZ_1	MZFAM	MZ_2	10000	9	17	9974	0.99740	0.99862	0.99512
POFAM	PO_P	POFAM	PO_C	10000	0	3511	6489	0.64890	1.00000	0.32364
```

At 0.30 a parent–offspring pair leaks in at `HetConc 0.32364`. The gap between a true
duplicate (≈ 1.0) and a first-degree relative (≈ 0.33 here) is enormous, so the default 0.80
is not a delicate choice — but do not lower it hoping to catch "near-duplicates".

**A clean cohort** produces a header-only file, not a missing one:

```
$ open-king -b /tmp/kingdocs/unrelated.bed --duplicate --prefix clean_
No duplicates are found with heterozygote concordance rate > 80%.

$ cat clean_.con
FID1	ID1	FID2	ID2	N	N_IBS0	N_IBS1	N_IBS2	Concord	HomConc	HetConc
```

> **Trap.** `--duplicate` scans all pairs regardless of FID, which is what you want. But if
> your cohort has 100 or more samples and contains no duplicates, the `.con` file is not
> created at all — read the console line, not the directory listing.

---

## 2. Find cryptic relatives in a supposedly unrelated cohort

**Goal:** you have a case–control cohort where everyone is supposed to be unrelated. Find
the pairs that are not.

KING splits pairs by FID: pairs *within* one family go to `.kin`, pairs *across* families go
to `.kin0`. In a cohort of nominally unrelated people every sample is its own family, so
`.kin0` is the file you want, and every row in it is by definition a relationship the
pedigree does not declare.

Using `bigish` (200 samples, 50 000 SNPs, with undeclared cross-family sibships built in):

```
$ open-king -b /tmp/kingdocs/bigish.bed --related --degree 2
[... startup banner and loading block elided ...]
Options in effect:
	--related
	--degree 2

Total length of 22 chromosomal segments usable for IBD segment analysis is 2498.9 Mb.
  Information of these chromosomal segments can be found in file kingallsegs.txt

Within-family kinship data saved in file king.kin

Relationship summary (total relatives: 436 by pedigree, 435 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	226	111	81	18	137
  Inference	0	226	111	79	19	138

A subset of informative SNPs will be used to screen close relatives.
Sorting autosomes...
Relationship inference across families starts at Fri Aug 14 22:49:25 2026
16 CPU cores are used...
  Stages 1&2 (with 32768 SNPs): 50 pairs of relatives are detected (with kinship > 0.0625)
                               Screening ends at Fri Aug 14 22:49:26 2026
  Final Stage (with 50000 SNPs): 26 pairs of relatives (up to 2nd-degree) are confirmed
                               Inference ends at Fri Aug 14 22:49:26 2026

Relationship summary (total relatives: 0 by pedigree, 26 by inference)
        	MZ	PO	FS	2nd	3rd	4th
  =========================================================
  Inference	0	0	3	23	0	0


Between-family relatives (kinship >= 0.08839) saved in file king.kin0
KING ends at Fri Aug 14 22:49:26 2026
```

```
$ cat king.kin0
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
BF01	B01_F	BF02	B02_F	50000	0.2274	0.0132	0.4855	0.1175	0.2885	0.4575	0.3676	0.5964	FS
BF01	B01_F	BF02	B02_C1	50000	0.1560	0.0299	0.2873	0.2261	0.1368	0.5684	0.0000	0.2842	2nd
BF01	B01_F	BF02	B02_C2	50000	0.1560	0.0295	0.2878	0.2227	0.1387	0.5856	0.0000	0.2928	2nd
BF01	B01_F	BF02	B02_C3	50000	0.1558	0.0285	0.2877	0.2127	0.1417	0.6069	0.0000	0.3035	2nd
BF01	B01_F	BF02	B02_C4	50000	0.1579	0.0291	0.2925	0.2184	0.1428	0.5971	0.0000	0.2986	2nd
BF01	B01_C1	BF02	B02_F	50000	0.1578	0.0242	0.2945	0.1832	0.1575	0.6610	0.0000	0.3305	2nd
BF01	B01_C2	BF02	B02_F	50000	0.1555	0.0338	0.2872	0.2528	0.1244	0.5244	0.0000	0.2622	2nd
BF01	B01_C2	BF02	B02_C2	50000	0.1506	0.0447	0.2743	0.3318	0.0870	0.3599	0.0000	0.1799	2nd
BF01	B01_C3	BF02	B02_F	50000	0.1570	0.0286	0.2921	0.2146	0.1434	0.6050	0.0000	0.3025	2nd
BF01	B01_C3	BF02	B02_C4	50000	0.1482	0.0431	0.2702	0.3138	0.0882	0.3802	0.0000	0.1901	2nd
BF13	B13_F	BF14	B14_F	50000	0.2304	0.0116	0.4911	0.1045	0.2959	0.4484	0.3767	0.6009	FS
BF13	B13_F	BF14	B14_C1	50000	0.1610	0.0252	0.2988	0.1905	0.1577	0.6232	0.0000	0.3116	2nd
BF13	B13_F	BF14	B14_C2	50000	0.1611	0.0262	0.2977	0.1983	0.1534	0.6211	0.0000	0.3105	2nd
BF13	B13_F	BF14	B14_C3	50000	0.1564	0.0320	0.2886	0.2389	0.1320	0.5460	0.0000	0.2730	2nd
BF13	B13_F	BF14	B14_C4	50000	0.1613	0.0296	0.2983	0.2211	0.1440	0.5755	0.0000	0.2878	2nd
BF13	B13_C1	BF14	B14_F	50000	0.1613	0.0286	0.2979	0.2148	0.1467	0.5998	0.0000	0.2999	2nd
BF13	B13_C2	BF14	B14_F	50000	0.1579	0.0303	0.2917	0.2235	0.1388	0.5735	0.0000	0.2867	2nd
BF13	B13_C3	BF14	B14_F	50000	0.1600	0.0271	0.2965	0.2037	0.1511	0.6159	0.0000	0.3080	2nd
BF25	B25_F	BF26	B26_F	50000	0.2099	0.0169	0.4288	0.1467	0.2504	0.5064	0.2682	0.5214	FS
BF25	B25_F	BF26	B26_C1	50000	0.1579	0.0344	0.2897	0.2568	0.1268	0.5250	0.0000	0.2625	2nd
BF25	B25_F	BF26	B26_C2	50000	0.1590	0.0322	0.2905	0.2456	0.1323	0.5328	0.0000	0.2664	2nd
BF25	B25_F	BF26	B26_C3	50000	0.1574	0.0337	0.2872	0.2553	0.1264	0.5099	0.0000	0.2549	2nd
BF25	B25_F	BF26	B26_C4	50000	0.1589	0.0344	0.2906	0.2628	0.1263	0.4977	0.0000	0.2488	2nd
BF25	B25_C1	BF26	B26_F	50000	0.1530	0.0356	0.2825	0.2654	0.1172	0.4943	0.0000	0.2472	2nd
BF25	B25_C2	BF26	B26_F	50000	0.1559	0.0294	0.2889	0.2209	0.1397	0.5946	0.0000	0.2973	2nd
BF25	B25_C3	BF26	B26_F	50000	0.1576	0.0346	0.2911	0.2626	0.1248	0.5056	0.0000	0.2528	2nd
```

**How to read it.** Twenty-six pairs, in three connected clusters. Three of them are `FS`
(kinship ≈ 0.25–0.30, and non-zero `IBD2Seg`): `B01_F`/`B02_F`, `B13_F`/`B14_F`,
`B25_F`/`B26_F` are full siblings who head two separately declared families. Everything else
is the knock-on: each of those men is a 2nd-degree relative (an uncle) of the other family's
children.

That is the usual shape of a cryptic-relatedness finding — one undeclared link, then a
fan-out of second-degree consequences. Fix the link, not the twenty-three consequences.

`--degree 2` sets how far out to look: it keeps pairs with kinship ≥ 2<sup>-(d+1.5)</sup>,
which the console prints explicitly (`kinship >= 0.08839`). Use `--degree 3` or `4` to reach
further; see [recipe 11](#11-large-cohorts---degree---cpus-and-what-actually-costs-time) for
what that costs.

**A cohort that really is unrelated** gives an empty (header-only) file, and says so:

```
$ open-king -b /tmp/kingdocs/unrelated.bed --kinship --degree 3 --prefix unr_
Between-family kinship data (up to degree 3, 0 pairs in total) saved in file unr_.kin0

$ cat unr_.kin0
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	Kinship
```

### If your cohort has fewer than 100 samples, do not use `--related` for this

`--related`'s between-family stage only writes `.kin0` when N ≥ 100. Below that it prints a
sentence and writes nothing — **even when the data contains an exact duplicate pair**:

```
$ open-king -b /tmp/kingdocs/dups.bed --related --prefix rel_
No close relatives are inferred.

$ ls rel_.kin0
ls: rel_.kin0: No such file or directory
```

`--kinship` has no such gate. On the same 10-sample fileset:

```
$ open-king -b /tmp/kingdocs/dups.bed --kinship --prefix kin_
Between-family kinship data saved in file kin_.kin0

$ grep -E 'DUP_A_COPY' kin_.kin0
UNR1	U1	DUPB	DUP_A_COPY	10000	0.1415	0.0692	0.0007
UNR2	U2	DUPB	DUP_A_COPY	10000	0.1365	0.0693	-0.0093
UNR3	U3	DUPB	DUP_A_COPY	10000	0.1274	0.0685	-0.0143
UNR4	U4	DUPB	DUP_A_COPY	10000	0.1346	0.0712	-0.0150
DUPA	DUP_A	DUPB	DUP_A_COPY	10000	0.3459	0.0000	0.5000
DUPB	DUP_A_COPY	MZFAM	MZ_1	10000	0.1348	0.0684	-0.0044
DUPB	DUP_A_COPY	MZFAM	MZ_2	10000	0.1346	0.0689	-0.0056
DUPB	DUP_A_COPY	POFAM	PO_P	10000	0.1344	0.0686	-0.0060
DUPB	DUP_A_COPY	POFAM	PO_C	10000	0.1289	0.0677	-0.0109
```

There it is: `DUP_A` / `DUP_A_COPY` at kinship `0.5000`. **For cohorts under 100 samples use
`--kinship` (and `--duplicate`), not `--related`.**

---

## 3. Verify a stated pedigree

**Goal:** the `.fam` says who is related to whom. Find the pairs where the genotypes
disagree.

Two independent signals do this, and you want both.

### Signal 1 — the relationship summary and the `Error` column

`--kinship` grades every declared within-family pair against its estimate and puts the result
in the last column of `.kin`. Using `dups`, where `MZ_1` and `MZ_2` share a FID but are
declared as unrelated founders:

```
$ open-king -b /tmp/kingdocs/dups.bed --kinship
[... startup banner and loading block elided ...]
Options in effect:
	--kinship

Within-family kinship data saved in file king.kin

Relationship summary (total relatives: 1 by pedigree, 2 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	1	0	0	0	1
  Inference	1	1	0	0	0	0

Relationship inference across families starts at Fri Aug 14 22:49:26 2026
16 CPU cores are used.
                                         ends at Fri Aug 14 22:49:26 2026
Between-family kinship data saved in file king.kin0
Note --kinship --degree <n> can filter & speed up the kinship computing.
KING ends at Fri Aug 14 22:49:26 2026
```

```
$ cat king.kin
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	Kinship	Error
MZFAM	MZ_1	MZ_2	10000	1.000	0.0000	0.3468	0.0009	0.4962	1
POFAM	PO_C	PO_P	10000	0.000	0.2500	0.1680	0.0000	0.2445	0
```

**How to read it.**

* The summary table is the fastest check in the whole program. Compare the two rows: the
  pedigree declares one PO pair and one "other"; the genotypes say one MZ pair and one PO
  pair. One declared relationship is wrong.
* `Z0` and `Phi` are **pedigree expectations**, not estimates — `Phi` is the pedigree kinship
  coefficient and `Z0` the pedigree Pr[IBD = 0]. `Phi 0.0000` / `Z0 1.000` means "the `.fam`
  says these two are unrelated".
* `Kinship` is the estimate. `0.4962` against a declared `0.0000` is the contradiction.
* `Error` is **not a flag and not an integer**. It takes the values `0`, `0.5` and `1`:
  a match, off by exactly one degree, and off by more than one. Filtering with `$NF == 1`
  and treating the column as an int will silently drop every half-step disagreement.

### A half-step disagreement, and why `--related` grades it differently

On `multifam` (four declared families):

```
$ open-king -b /tmp/kingdocs/multifam.bed --kinship --prefix mf_
Relationship summary (total relatives: 36 by pedigree, 36 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	24	12	0	0	4
  Inference	0	24	11	1	0	4

$ head -1 mf_.kin; awk -F'\t' 'NR>1 && $NF != "0"' mf_.kin
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	Kinship	Error
FAM2	B_C1	B_C2	15000	0.250	0.2500	0.1791	0.0295	0.1708	0.5
```

One declared full-sib pair estimates at `0.1708`, just under the 0.177 first-degree
boundary, so it is graded 2nd degree and flagged `0.5`. Run the same fileset through
`--related`, which decides using IBD segments rather than the kinship point estimate:

```
$ open-king -b /tmp/kingdocs/multifam.bed --related --prefix mfr_
$ awk -F'\t' 'NR==1 || $2=="B_C1" && $3=="B_C2"' mfr_.kin
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType	Error
FAM2	B_C1	B_C2	15000	0.250	0.2500	0.1791	0.0295	0.3420	0.2642	0.1708	0.4582	0.1273	0.3564	FS	0
```

Same kinship (`0.1708`), but `IBD2Seg 0.1273` — these two share both chromosomes over 12.7 %
of the genome, which only full siblings do — so `InfType` is `FS` and `Error` is `0`. The
pedigree was right and the kinship estimate was noisy.

**Rule of thumb: near a boundary, trust the segments over the kinship estimate.** The two
`Error` columns are computed by different rules and legitimately disagree; the details are in
[VERIFIED_FORMULAS.md § `--related`'s extra columns](VERIFIED_FORMULAS.md).

### Signal 2 — relatives the pedigree does not declare

The `.kin` `Error` column can only catch pairs the pedigree already puts in one family.
Relationships between people the pedigree calls strangers show up in `.kin0`:

```
$ head -1 mf_.kin0; awk -F'\t' '$8+0 > 0.177' mf_.kin0
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	Kinship
FAM1	A_F	FAM2	B_F	15000	0.2169	0.0204	0.2501
FAM1	A_F	FAM3	C_F	15000	0.1772	0.0000	0.2526
FAM1	A_M	FAM3	C_F	15000	0.1737	0.0000	0.2466
FAM1	A_C1	FAM3	C_F	15000	0.2011	0.0195	0.2309
FAM1	A_C2	FAM3	C_F	15000	0.2307	0.0151	0.2845
FAM1	A_C3	FAM3	C_F	15000	0.2071	0.0173	0.2440
FAM3	C_M	FAM4	D_M	15000	0.1753	0.0233	0.1843
```

Eight first-degree pairs across family boundaries. Note the `IBS0` column, which tells you
*which kind*: `A_F`/`C_F` and `A_M`/`C_F` have `IBS0 = 0.0000` — that is parent and child
(`C_F` is the son of the `FAM1` couple). `A_F`/`B_F` at `IBS0 0.0204` is a sibling pair, and
`A_C1`…`A_C3` versus `C_F` are `C_F`'s siblings. See
[recipe 6](#6-tell-parentoffspring-from-full-siblings).

---

## 4. Produce a maximal unrelated subset for GWAS

**Goal:** drop the smallest number of samples such that no two survivors are related.

```
$ open-king -b /tmp/kingdocs/bigish.bed --unrelated
[... startup banner and loading block elided ...]
Options in effect:
	--unrelated

Family clustering starts at Fri Aug 14 22:49:26 2026
Autosome genotypes stored in 782 words for each of 200 individuals.
Sorting autosomes...
Total length of 22 chromosomal segments usable for IBD segment analysis is 2498.9 Mb.
  Information of these chromosomal segments can be found in file kingallsegs.txt

16 CPU cores are used to compute the pairwise kinship coefficients...
Clustering up to 1st-degree relatives in families...
Individual IDs are unique across all families.

Relationship summary (total relatives: 0 by pedigree, 3 by inference)
        	MZ	PO	FS	2nd	3rd	4th
  =========================================================
  Inference	0	0	3	0	0	0

The following families are found to be connected
  NewFamID  OriginalFamID
  KING1     BF01,BF02
  KING2     BF13,BF14
  KING3     BF25,BF26


A list of 84 unrelated individuals saved in file kingunrelated.txt
An alternative list of 116 to-be-removed individuals saved in file kingunrelated_toberemoved.txt

Extracting a subset of unrelated individuals ends at Fri Aug 14 22:49:27 2026
KING ends at Fri Aug 14 22:49:27 2026
```

```
$ wc -l kingunrelated.txt kingunrelated_toberemoved.txt
      84 kingunrelated.txt
     116 kingunrelated_toberemoved.txt
     200 total

$ head -5 kingunrelated.txt
BF03	B03_S2
BF03	B03_S1
BF03	B03_G_F
BF03	B03_G_M
BF04	B04_F
```

**How to read it.** Two `FID<TAB>IID` lists that partition the cohort: the ones to keep and
the ones to drop. Feed either straight to PLINK:

```
plink --bfile study --keep kingunrelated.txt --make-bed --out study_unrel
```

Note that four of the first five keepers share FID `BF03` — the subset is maximal, not
one-per-family. Where a family contains people who are not related *to each other* (two
grandparents and two in-marrying spouses, here), they all survive.

**Verify it yourself.** Do not take the file on trust; check every surviving pair against the
kinship run:

```
$ open-king -b /tmp/kingdocs/bigish.bed --kinship --prefix all_

$ python3 check_unrelated.py kingunrelated.txt all_.kin all_.kin0
kept samples           : 84
pairs checked          : 3486 (all 3486 of them)
largest kinship in set : 0.0159  (B18_G_F / BSNG002)
```

All 3 486 pairs among the 84 survivors sit below 0.0159, comfortably under the 0.0221
unrelated cutoff. `check_unrelated.py` is twelve lines and reads both `.kin` (within-family
pairs) and `.kin0` (between-family pairs) — you need both, because `--unrelated` can and does
keep several members of one FID:

```python
import sys
kept = {tuple(l.split()) for l in open(sys.argv[1])}
worst, n = ("", "", -9.0), 0
def look(a, b, k):
    global worst, n
    if a in kept and b in kept:
        n += 1
        if float(k) > worst[2]: globals()["worst"] = (a, b, float(k))
for ln in open(sys.argv[2]):                       # .kin  : within-family pairs
    f = ln.rstrip("\n").split("\t")
    if f[0] != "FID": look((f[0], f[1]), (f[0], f[2]), f[8])
for ln in open(sys.argv[3]):                       # .kin0 : between-family pairs
    f = ln.rstrip("\n").split("\t")
    if f[0] != "FID1": look((f[0], f[1]), (f[2], f[3]), f[7])
print("kept samples           :", len(kept))
print("pairs checked          :", n, "(all %d of them)" % (len(kept)*(len(kept)-1)//2))
print("largest kinship in set : %.4f  (%s / %s)" % (worst[2], worst[0][1], worst[1][1]))
```

**On `--degree`.** `--unrelated --degree 2` changes the family-clustering step (the console
says `Clustering up to 2nd-degree relatives in families`), but on this dataset the selected
subset is byte-identical to the default at degrees 1, 2 and 3 — 84 kept, 116 removed, same
people. That is a property of `bigish`, whose undeclared links are all first degree. Check
your own data rather than assuming.

---

## 5. Relatedness in an admixed or structured cohort

**Goal:** estimate kinship when your samples do not all come from one homogeneous
population — and understand why KING's between-family estimator is the right one there.

### The two estimators

KING computes kinship two different ways and picks between them **purely on whether the two
samples share a FID**:

```
                                HetHet - 2*IBS0
within-family  (.kin, .ibs)  =  ---------------
                                 Het_i + Het_j

                                      2*HetHet - 4*IBS0 - Het_i - Het_j
between-family (.kin0, .ibs0)  = 0.5 + ---------------------------------
                                            4 * min(Het_i, Het_j)
```

The between-family form — "KING-robust" — divides by **four times the smaller of the two
heterozygote counts**, not by their sum. That `min()` is the entire robustness mechanism, and
you can watch it work.

### Watching it work

`admixed` has three nuclear families: `FAMX` inside population 1, `FAMZ` inside population 2,
and `FAMY` **across** the two — a population-1 father, a population-2 mother, and two admixed
children. Run `--ibs` to get the raw counts, then run the same `.bed` again with a `.fam` in
which those three families are split into singleton FIDs, which forces the between-family
form onto exactly the same pairs:

```
$ open-king -b /tmp/kingdocs/admixed.bed --ibs --prefix fam_
Within-family IBS data saved in file fam_.ibs
Between-family IBS data saved in file fam_.ibs0

$ awk '{ if ($1 ~ /^FAM[XYZ]$/) print "S"$2, $2, 0, 0, $5, $6; else print }' /tmp/kingdocs/admixed.fam > split.fam

$ head -2 split.fam; grep -c . split.fam
AP1_01 A1_01 0 0 1 -9
AP1_02 A1_02 0 0 2 -9
40

$ open-king -b /tmp/kingdocs/admixed.bed --fam split.fam --kinship --prefix split_
Between-family kinship data saved in file split_.kin0
```

`--fam` swaps the pedigree without touching the genotypes, so the two runs describe the same
18 pairs. Recomputing both formulas from the counts and asserting them against what KING
wrote:

```
$ python3 two_estimators.py
pair          N_Het1  N_Het2  HetHet   IBS0   within  between     diff
X_C1/X_C2       6568    6590    4114    263   0.2727   0.2723  -0.0004
X_C1/X_F        6568    6607    3252      0   0.2468   0.2461  -0.0008
X_C1/X_M        6568    6635    3222      0   0.2440   0.2427  -0.0013
X_C2/X_F        6590    6607    3231      0   0.2448   0.2445  -0.0003
X_C2/X_M        6590    6635    3303      0   0.2498   0.2489  -0.0009
X_F/X_M         6607    6635    2700   1396  -0.0069  -0.0080  -0.0011
Y_C1/Y_C2       7364    7464    4389    448   0.2356   0.2338  -0.0018
Y_C1/Y_F        7364    6640    3293      0   0.2351   0.2207  -0.0144
Y_C1/Y_M        7364    6666    3331      0   0.2374   0.2237  -0.0137
Y_C2/Y_F        7464    6640    3369      0   0.2389   0.2227  -0.0162
Y_C2/Y_M        7464    6666    3368      0   0.2384   0.2227  -0.0157
Y_F/Y_M         6640    6666    2423   1959  -0.1124  -0.1136  -0.0012
Z_C1/Z_C2       6568    6619    3830    510   0.2131   0.2120  -0.0011
Z_C1/Z_F        6568    6520    3237      0   0.2473   0.2464  -0.0009
Z_C1/Z_M        6568    6556    3322      0   0.2531   0.2529  -0.0002
Z_C2/Z_F        6619    6520    3312      0   0.2521   0.2502  -0.0019
Z_C2/Z_M        6619    6556    3295      0   0.2501   0.2489  -0.0012
```

The script asserts that its `within` column equals the `Kinship` column of `fam_.ibs` and its
`between` column equals the `Kinship` column of `split_.kin0`, to all four printed decimals,
on every row — so these are KING's own numbers, re-derived:

```python
# Both KING-robust forms recomputed from the raw counts in <prefix>.ibs, and
# checked against the two files KING itself wrote.
ibs, rob = {}, {}
for ln in open("fam_.ibs"):            # families as declared -> within-family form
    f = ln.rstrip("\n").split("\t")
    if f[0] != "FID":
        ibs[tuple(sorted((f[1], f[2])))] = (int(f[6]), int(f[9]), int(f[11]), int(f[12]), f[19])
for ln in open("split_.kin0"):         # same .bed, families split -> between-family form
    f = ln.rstrip("\n").split("\t")
    if f[0] != "FID1":
        rob[tuple(sorted((f[1], f[3])))] = f[7]
print("%-12s %7s %7s %7s %6s %8s %8s %8s"
      % ("pair", "N_Het1", "N_Het2", "HetHet", "IBS0", "within", "between", "diff"))
for k in sorted(ibs):
    if k not in rob: continue
    ibs0, hh, h1, h2, kin = ibs[k]
    within  = (hh - 2*ibs0) / (h1 + h2)
    between = 0.5 + (2*hh - 4*ibs0 - h1 - h2) / (4 * min(h1, h2))
    assert "%.4f" % within  == kin,    (k, within,  kin)      # == the .ibs / .kin column
    assert "%.4f" % between == rob[k], (k, between, rob[k])   # == the .kin0 column
    # and the exact difference between the two forms
    assert abs((between - within) - (h1 + h2 - 2*min(h1,h2)) / (4*min(h1,h2)) * (2*within - 1)) < 1e-12
    print("%-12s %7d %7d %7d %6d %8.4f %8.4f %+8.4f"
          % ("/".join(k), h1, h2, hh, ibs0, within, between, between - within))
```

**How to read it.**

* For `FAMX` and `FAMZ` — both parents from one population — the two estimators agree to
  within 0.002. When `N_Het1 ≈ N_Het2`, `min()` and the average are the same thing.
* For `FAMY` — the cross-population family — the parent–child rows diverge by up to
  **0.0162**. The admixed children carry ~7 400 heterozygous sites against their parents'
  ~6 650: admixture raises heterozygosity, `N_Het1 ≠ N_Het2`, and the two formulas part
  company.

The relationship is exact, and holds on all 18 rows (the script asserts this too):

```
between - within  =  (Het_i + Het_j - 2*min) / (4*min)  *  (2*within - 1)
```

Since kinship is at most 0.5, the second factor is never positive, so **the between-family
estimator is always ≤ the within-family one, with equality exactly when the two heterozygote
counts are equal.** Ancestry difference can only push a KING-robust estimate *down*.

### Why that is the property you want

Look at what happens to a pair of *unrelated* people from the two different populations:

```
$ awk -F'\t' '$2=="A1_01" && $4=="A2_01"' split_.kin0
AP1_01	A1_01	AP2_01	A2_01	20000	0.1240	0.0984	-0.1123
```

Kinship `-0.1123`. Population structure drives unrelated cross-population pairs strongly
*negative*, not positive. **KING-robust does not invent relatives out of ancestry.** That is
the guarantee: an estimator that used pooled allele frequencies would read shared ancestry as
shared descent and fill your `.kin0` with spurious "3rd-degree" pairs; this one cannot.

The price is on the other side, and the `FAMY` rows show it: a genuine parent–child pair
whose members differ in ancestry is estimated at `0.2207` instead of the true 0.25. The bias
is downward, so **admixed distant relatives get under-called** — a true 4th-degree admixed
pair can land below the 0.0221 cutoff and vanish. If distant relatedness in an admixed cohort
matters to you, lower your threshold deliberately, or use the IBD-segment columns
([recipe 7](#7-ibd-segments-and-what-they-can-and-cannot-separate)), which do not depend on
the heterozygosity ratio in the same way.

### Practical guidance

* **Do not group people of different ancestry under one FID** in order to get "within-family"
  numbers. The within-family form assumes the pair share a heterozygosity rate; if they do
  not, it is the wrong estimator and KING has no way to know.
* When in doubt, give every sample its own FID and read `.kin0` only. That guarantees the
  robust form everywhere.
* Read a negative kinship as "unrelated, and from a different genetic background", not as an
  error.

---

## 6. Tell parent–offspring from full siblings

**Goal:** both have kinship 0.25. Separate them.

The discriminator is **IBS0** — the fraction of markers where the two are opposite
homozygotes. A parent and child share one allele at every locus by descent, so IBS0 is
structurally impossible except through genotyping error. Full siblings have Pr[IBD = 0] =
0.25, so a quarter of the genome is free to produce opposite homozygotes.

Using `multifam`, whose four declared families contain both kinds:

```
$ open-king -b /tmp/kingdocs/multifam.bed --kinship
$ awk -F'\t' 'NR==1 || $9+0 > 0.177' king.kin | cut -f1-3,6,8,9 | head -20
FID	ID1	ID2	Phi	IBS0	Kinship
FAM1	A_C1	A_C2	0.2500	0.0153	0.2721
FAM1	A_C1	A_C3	0.2500	0.0214	0.2476
FAM1	A_C1	A_F	0.2500	0.0000	0.2475
FAM1	A_C1	A_M	0.2500	0.0000	0.2444
FAM1	A_C2	A_C3	0.2500	0.0124	0.2837
FAM1	A_C2	A_F	0.2500	0.0000	0.2523
FAM1	A_C2	A_M	0.2500	0.0000	0.2485
FAM1	A_C3	A_F	0.2500	0.0000	0.2508
FAM1	A_C3	A_M	0.2500	0.0000	0.2484
FAM2	B_C1	B_C3	0.2500	0.0105	0.3019
FAM2	B_C1	B_F	0.2500	0.0000	0.2480
FAM2	B_C1	B_M	0.2500	0.0000	0.2501
FAM2	B_C2	B_C3	0.2500	0.0254	0.1809
FAM2	B_C2	B_F	0.2500	0.0000	0.2427
FAM2	B_C2	B_M	0.2500	0.0000	0.2512
FAM2	B_C3	B_F	0.2500	0.0000	0.2512
FAM2	B_C3	B_M	0.2500	0.0000	0.2502
FAM3	C_C1	C_C2	0.2500	0.0146	0.2735
FAM3	C_C1	C_C3	0.2500	0.0117	0.2537
```

Splitting every first-degree row by what the pedigree says it is (`Z0 = 0.000` is
parent–offspring, `Z0 = 0.250` is full sibs) and listing the IBS0 values:

```
$ awk -F'\t' 'NR>1 && $9+0>0.177 {t=($5=="0.000")?"PO":"FS"; print t, $8}' king.kin | sort |
      awk '{a[$1]=a[$1]" "$2} END {for (k in a) print k":"a[k]}'
FS: 0.0105 0.0117 0.0124 0.0146 0.0153 0.0163 0.0214 0.0242 0.0244 0.0253 0.0254
PO: 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000 0.0000
```

**Perfect separation.** All 24 parent–offspring pairs are at exactly `0.0000`; all 11 sib
pairs are at 0.0105 or above. Kinship cannot make this call: the PO range here is
0.2427–0.2527 and the sibling range is 0.1809–0.3019, so the parent–offspring interval sits
*entirely inside* the sibling one.

**How to read it in practice.**

* `IBS0 = 0.0000` with kinship in the first-degree band ⇒ parent–offspring.
* `IBS0` of the order of 1–3 % with kinship in the first-degree band ⇒ full siblings.
* This corpus is simulated without genotyping error on those pairs, which is why the PO column
  is *exactly* zero. On real arrays a true PO pair shows a small non-zero IBS0 from error —
  the classic real-data figures are a few tenths of a percent, still an order of magnitude
  below the sibling rate. The reference binary picks its own data-derived cutoff in
  `[0.0035, 0.0060]` rather than using a fixed constant; `0.0055` on ordinary data. It is
  never written to a file. See
  [BEHAVIOR.md § Q2](BEHAVIOR.md#q2--parentoffspring-vs-full-sibling-discrimination).
* `--related` and `--ibdseg` do not use IBS0 for this at all: when IBD segments are available
  they split PO from FS on **IBD2 sharing**, which is cleaner still — see the next recipe.

---

## 7. IBD segments, and what they can and cannot separate

**Goal:** several relationships share a kinship coefficient. Half-siblings, grandparent–
grandchild and avuncular pairs are all φ = 0.125. Can the segment output tell them apart?

`--ibdseg` calls the IBD segments each pair shares and summarises them per pair:

```
$ open-king -b /tmp/kingdocs/threegen.bed --ibdseg
[... startup banner and loading block elided ...]
kingsplitped.txt is generated for certain pedigree plot applications.

Options in effect:
	--ibdseg

Total length of 21 chromosomal segments usable for IBD segment analysis is 982.7 Mb.
  Information of these chromosomal segments can be found in file kingallsegs.txt

IBD segment analysis starts at Fri Aug 14 22:50:08 2026
16 CPU cores are used for autosome inference...
                       ends at Fri Aug 14 22:50:08 2026

Note with relationship inference as the primary goal, the following filters are applied:
  Sample pairs without any long IBD segments (>10Mb) are excluded.
  Short IBD segments (<3Mb) are not reported/utilized.
Summary statistics of IBD segments for individual pairs saved in file king.seg
KING ends at Fri Aug 14 22:50:08 2026
```

`threegen` is one twelve-person family with every 2nd- and 3rd-degree relationship type in
it, and `MANIFEST.json` records the truth for each pair. Joining the truth, the kinship
estimate and the segment summary (kinship comes from `--ibs`, because a single-family `.kin`
comes out empty — see [Traps](#traps)):

```
$ open-king -b /tmp/kingdocs/threegen.bed --ibs --prefix ibs_
$ python3 by_relationship.py
true phi    pair         Kinship   IBD1Seg   IBD2Seg  PropIBD InfType
PO   0.2500 C1/P1         0.2499    1.0000    0.0000   0.5000 PO
PO   0.2500 C1/S1         0.2486    1.0000    0.0000   0.5000 PO
PO   0.2500 C2/P1         0.2497    1.0000    0.0000   0.5000 PO
PO   0.2500 C2/S1         0.2443    1.0000    0.0000   0.5000 PO
PO   0.2500 C3/P2         0.2503    1.0000    0.0000   0.5000 PO
PO   0.2500 C3/S2         0.2475    1.0000    0.0000   0.5000 PO
PO   0.2500 C4/P2         0.2508    1.0000    0.0000   0.5000 PO
PO   0.2500 C4/S2         0.2530    1.0000    0.0000   0.5000 PO
PO   0.2500 GF/P1         0.2541    1.0000    0.0000   0.5000 PO
PO   0.2500 GF/P2         0.2498    1.0000    0.0000   0.5000 PO
PO   0.2500 GF/P3         0.2501    1.0000    0.0000   0.5000 PO
PO   0.2500 GM1/P1        0.2496    1.0000    0.0000   0.5000 PO
PO   0.2500 GM1/P2        0.2486    1.0000    0.0000   0.5000 PO
PO   0.2500 GM2/P3        0.2514    1.0000    0.0000   0.5000 PO
FS   0.2500 C1/C2         0.2792    0.4298    0.3389   0.5538 FS
FS   0.2500 C3/C4         0.2702    0.6017    0.2444   0.5453 FS
FS   0.2500 P1/P2         0.2517    0.3336    0.3479   0.5147 FS
GG   0.1250 C1/GF         0.1390    0.5543    0.0000   0.2772 2nd
GG   0.1250 C1/GM1        0.1135    0.4475    0.0000   0.2238 2nd
GG   0.1250 C2/GF         0.0875    0.3433    0.0000   0.1716 3rd
GG   0.1250 C2/GM1        0.1664    0.6449    0.0000   0.3225 2nd
GG   0.1250 C3/GF         0.1009    0.3548    0.0000   0.1774 2nd
GG   0.1250 C3/GM1        0.1531    0.6079    0.0000   0.3039 2nd
GG   0.1250 C4/GF         0.0783    0.2767    0.0000   0.1384 3rd
GG   0.1250 C4/GM1        0.1736    0.7151    0.0000   0.3576 2nd
AV   0.1250 C1/P2         0.1228    0.4339    0.0000   0.2170 2nd
AV   0.1250 C2/P2         0.1234    0.4737    0.0000   0.2369 2nd
AV   0.1250 C3/P1         0.1329    0.4390    0.0000   0.2195 2nd
AV   0.1250 C4/P1         0.1193    0.4136    0.0000   0.2068 2nd
HS   0.1250 P1/P3         0.1271    0.4571    0.0000   0.2286 2nd
HS   0.1250 P2/P3         0.1163    0.4987    0.0000   0.2494 2nd
HAV  0.0625 C1/P3         0.0666    0.2515    0.0000   0.1258 3rd
HAV  0.0625 C2/P3         0.0510    0.1977    0.0000   0.0989 3rd
HAV  0.0625 C3/P3         0.0519    0.2011    0.0000   0.1006 3rd
HAV  0.0625 C4/P3         0.0258    0.1040    0.0000   0.0520 4th
FC   0.0625 C1/C3         0.0684    0.2240    0.0000   0.1120 3rd
FC   0.0625 C1/C4         0.0448    0.1559    0.0000   0.0780 4th
FC   0.0625 C2/C3         0.0611    0.2118    0.0000   0.1059 3rd
FC   0.0625 C2/C4         0.0730    0.3167    0.0000   0.1584 3rd
```

(`GG` = grandparent–grandchild, `AV` = avuncular, `HS` = half-sib, `HAV` = half-avuncular,
`FC` = first cousins. `PropIBD = IBD2Seg + IBD1Seg/2`, so it is on the 2× kinship scale:
0.5 for PO, ~0.25 for 2nd degree.)

The join, `by_relationship.py`:

```python
import json
truth = {}
for a, b, lab, phi in json.load(open("/tmp/kingdocs/MANIFEST.json"))["datasets"]["threegen"]["pairs"]:
    truth[tuple(sorted((a.split(":")[1], b.split(":")[1])))] = (lab, phi)
kin = {}
for ln in open("ibs_.ibs"):
    f = ln.rstrip("\n").split("\t")
    if f[0] != "FID": kin[tuple(sorted((f[1], f[2])))] = f[19]
print("%-4s %-6s %-11s %8s %9s %9s %8s %-7s"
      % ("true", "phi", "pair", "Kinship", "IBD1Seg", "IBD2Seg", "PropIBD", "InfType"))
rows = []
for ln in open("king.seg"):
    f = ln.rstrip("\n").split("\t")
    if f[0] == "FID1": continue
    k = tuple(sorted((f[1], f[3])))
    rows.append((["PO","FS","GG","AV","HS","HAV","FC"].index(truth[k][0]), k, f))
for _, k, f in sorted(rows):
    lab, phi = truth[k]
    print("%-4s %-6.4f %-11s %8s %9s %9s %8s %-7s"
          % (lab, phi, "/".join(x[3:] for x in k), kin[k], f[4], f[5], f[6], f[7]))
```

**What the segments buy you.**

* **PO becomes exact.** Every one of the 14 parent–offspring pairs reads `IBD1Seg 1.0000,
  IBD2Seg 0.0000, PropIBD 0.5000` — not approximately, exactly, on all fourteen. The kinship
  estimates for the same pairs scatter over 0.2443–0.2541. If you need parent–offspring
  called with certainty, this column is the one.
* **PO versus FS is decided by `IBD2Seg`.** All three full-sib pairs share both chromosomes
  over 24–35 % of the genome; every non-sib pair in the table is at `0.0000`. This is the
  same call [recipe 6](#6-tell-parentoffspring-from-full-siblings) makes with IBS0, but
  without a threshold to tune.

**What the segments do not buy you.** Look at the `GG`, `AV` and `HS` blocks. All three are
φ = 0.125, all three have `IBD2Seg` exactly `0.0000`, and their `PropIBD` ranges overlap
completely (GG 0.1384–0.3576, AV 0.2068–0.2369, HS 0.2286–0.2494). **Nothing in `.seg`
separates half-sib from grandparent from avuncular.** It cannot: all three share, in
expectation, one quarter of the genome on one chromosome only.

What *does* separate them in principle is the **number and length distribution** of the
shared segments — a grandparent–grandchild pair shares fewer, longer segments than an
avuncular pair of the same total. KING's summary files do not carry that: `.seg` reports four
aggregate numbers per pair, and the per-segment file the manual documents (`.segments.gz`) is
**never written by the 2.3.2 reference build** — that build has no zlib in its segment writer
— so open-king does not produce it either. If you need HS/GP/AV resolution, you need a
different tool (or non-genetic information such as age).

**Also worth seeing in that table: 2nd-degree calls are genuinely noisy.** Two true
grandparent–grandchild pairs (`C2/GF`, `C4/GF`) come out `3rd`, and one true half-avuncular
pair (`C4/P3`) comes out `4th`. That is real IBD variance over 22 chromosomes, not an
implementation defect — the same rows are byte-identical to the reference binary. Treat a
single 2nd-versus-3rd-degree call as a coin-flip-ish call, and look at the whole family.

### `--seglength`

`--seglength` sets the reporting floor for a segment, in Mb (default 3):

```
$ for l in 3 10; do open-king -b /tmp/kingdocs/threegen.bed --ibdseg --seglength $l --prefix sl$l\_ >/dev/null; done

$ diff sl3_.seg sl10_.seg | head -14
11c11
< TG	TG_GM1	TG	TG_C1	0.4475	0.0000	0.2238	2nd
---
> TG	TG_GM1	TG	TG_C1	0.4285	0.0000	0.2143	2nd
16,17c16,17
< TG	TG_P1	TG	TG_P2	0.3336	0.3479	0.5147	FS
< TG	TG_P1	TG	TG_P3	0.4571	0.0000	0.2286	2nd
---
> TG	TG_P1	TG	TG_P2	0.3540	0.3090	0.4860	FS
> TG	TG_P1	TG	TG_P3	0.4482	0.0000	0.2241	2nd

$ diff sl3_.seg sl10_.seg | grep -c '^[<>]'
22
```

Raising the floor to 10 Mb changes 11 of the 39 rows' numbers but not the row set and not a
single `InfType`. Short segments contribute little to the totals and are the ones most likely
to be spurious; raising the floor is a reasonable conservatism, not a different analysis.

---

## 8. Per-sample QC, per-SNP QC, and auto-QC

**Goal:** call rate, heterozygosity, Mendelian errors, and a sex check — before you compute
anything else.

### `--bysample`

```
$ open-king -b /tmp/kingdocs/missing.bed --bysample
QC-by-sample starts at Fri Aug 14 22:50:34 2026
There are 8 parent-offspring pairs and 4 trios, and 6 full-sibling pairs according to the pedigree.
QC starts...
  QC-by-sample ends at Fri Aug 14 22:50:34 2026
QC statistics by samples saved in file kingbySample.txt

KING ends at Fri Aug 14 22:50:34 2026

$ cat kingbySample.txt
FID IID FA MO SEX N_SNP Missing Heterozygosity N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal
MIS M_F 0 0 1 9784 0.0216 0.3505 31542 0 0.0000 31219 0 0.0000 0
MIS M_M 0 0 2 9711 0.0289 0.3420 31305 0 0.0000 31219 0 0.0000 0
MIS M_C1 M_F M_M 1 9337 0.0663 0.3449 18465 0 0.0000 9167 0 0.0000 0
MIS M_C2 M_F M_M 2 7831 0.2169 0.3504 15497 0 0.0000 7698 0 0.0000 0
MIS M_C3 M_F M_M 1 4801 0.5199 0.3456 9497 0 0.0000 4720 0 0.0000 0
MIS M_C4 M_F M_M 2 9788 0.0212 0.3441 19388 0 0.0000 9634 0 0.0000 0
```

`bySample.txt` is **space** separated (`.kin`, `.kin0`, `.con`, `.ibs` are tab separated —
this asymmetry is real). `Missing` is the per-sample missing-call rate: `M_C3` at `0.5199` is
half-empty and should go. `Heterozygosity` is stable at ~0.345 across all six — a sample
whose heterozygosity is far from its cohort's mode is a contamination or inbreeding signal.
`N_MIp` / `N_MIt` count Mendelian inconsistencies against declared parents and trios; all
zero here, because this fileset was simulated without transmission errors.

### The sex check hides in `--bysample`

On a fileset that has X and Y markers the header grows, and those extra columns are the sex
check:

```
$ open-king -b /tmp/kingdocs/sexchr.bed --bysample --prefix sex_
$ cat sex_bySample.txt
FID IID FA MO SEX N_SNP Missing Heterozygosity N_xSNP xHeterozygosity N_ySNP N_yHetero N_mtSNP N_mtHetero N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal
SEX S_F 0 0 1 4150 0.0000 0.3489 1500 0.0000 300 0 50 0 16600 0 0.0000 16600 0 0.0000 0
SEX S_M 0 0 2 4150 0.0000 0.3465 1500 0.3393 0 0 50 0 16600 0 0.0000 16600 0 0.0000 0
SEX S_SON1 S_F S_M 1 4150 0.0000 0.3441 1500 0.0000 300 0 50 0 8300 0 0.0000 4150 0 0.0000 0
SEX S_SON2 S_F S_M 1 4150 0.0000 0.3552 1500 0.0000 300 0 50 0 8300 0 0.0000 4150 0 0.0000 0
SEX S_DAU1 S_F S_M 2 4150 0.0000 0.3436 1500 0.3313 0 0 50 0 8300 0 0.0000 4150 0 0.0000 0
SEX S_DAU2 S_F S_M 2 4150 0.0000 0.3518 1500 0.3227 0 0 50 0 8300 0 0.0000 4150 0 0.0000 0
SU1 S_U0A 0 0 0 4150 0.0000 0.3472 1500 0.3733 0 0 50 0 0 0 0.0000 0 0 0.0000 0
SU2 S_U0B 0 0 0 4150 0.0000 0.3518 1500 0.3553 0 0 50 0 0 0 0.0000 0 0 0.0000 0
SU3 S_UM 0 0 1 4150 0.0000 0.3554 1500 0.0000 300 0 50 0 0 0 0.0000 0 0 0.0000 0
SU4 S_UF 0 0 2 4150 0.0000 0.3489 1500 0.3653 0 0 50 0 0 0 0.0000 0 0 0.0000 0
```

Read `SEX` against `xHeterozygosity` and `N_ySNP`: males have `xHeterozygosity 0.0000` and
300 Y markers called; females have `xHeterozygosity ≈ 0.33` and `N_ySNP 0`. Any row where
those disagree with the declared `SEX` is a sample swap or a mislabel. The two samples coded
`SEX 0` (`S_U0A`, `S_U0B`) both look female by this test — that is how you fill in an unknown
sex column.

> The header of `bySample.txt` has six variants depending on which chromosomes and pedigree
> structures your fileset has. Parse by header name, never by column position.

### `--bySNP`

```
$ open-king -b /tmp/kingdocs/missing.bed --bySNP --prefix snp_
QC statistics by SNPs saved in file snp_bySNP.txt

$ head -3 snp_bySNP.txt
SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate N_PO N_HomPO N_errPO Err_InPO Err_InHomPO N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio
rs1_1003166 1 1003166 G C 0.2500 4 0 2 2 0.6667 4 1 0 0.0000 0.0000 2 1 0 0.0000 0.0000
rs1_1050971 1 1050971 C G 0.5000 6 1 4 1 1.0000 8 0 0 0.0000 0.0000 4 2 0 0.0000 0.0000

$ tail -n +2 snp_bySNP.txt | awk '$11 == 0' | wc -l
      70

$ tail -n +2 snp_bySNP.txt | awk '$11 < 0.5' | wc -l
     220
```

Column 11 is `CallRate`. Seventy markers in this fileset were called in nobody, and 220 in
fewer than half the samples. `Freq_A`, the allele frequency, and `N_errPO`/`N_errTrio`, the
per-marker Mendelian error counts, are the other two things people filter on.

### `--autoQC`

`--autoQC` is the one pass in KING that applies filters. Everything else — every kinship,
IBS and segment number on this page — uses **every** marker in the `.bim` with no call-rate,
MAF or monomorphic filter at all.

```
$ open-king -b /tmp/kingdocs/missing.bed --autoQC --prefix qc_
Auto-QC step 1: Apply SNP call rate filter 80.0% on 10000 SNPs (in 6 samples)
  1569 autosome SNPs have call rate < 80.0%
  0 X-chr SNPs have call rate < 80.0%
  3199 SNPs are monomorphic

Auto-QC step 2: Apply sample call rate filter 95.0% on 6 samples (with 5232 SNPs)
  2 samples have call rate < 95.0%

Auto-QC step 3: Apply SNP call rate filter 95.0% on 5232 SNPs (in 4 samples)
  130 SNPs have call rate < 95.0%
  0 chr-X SNPs have call rate < 95.0%
X-Chr SNPs are not available. Gender QC is skipped.
Y-Chr SNPs are not available. Gender QC is skipped.

Auto-QC step 7: Final check
  4 samples, 5086 autosome SNPs

Auto-QC step 8: QC Summary Report

Step Description                                            Subjects  SNPs
1    Raw data counts                                        6         10000
1.1  SNPs with very low call rate < 80% (removed)                     (1569)
1.2  Monomorphic SNPs (removed)                                       (3199)
1.3  Sample call rate < 95% (removed)                       (2)
1.4  SNPs with call rate < 95% (removed)                              (130)
3    Generate Final Study Files
     Final QC'ed data                                       4         5086

QC summary report saved in qc__autoQC_Summary.txt
SNP-removal QC file saved in qc__autoQC_snptoberemoved.txt
Sample-removal QC file saved in qc__autoQC_sampletoberemoved.txt

Auto-QC ends at Fri Aug 14 22:50:34 2026
KING ends at Fri Aug 14 22:50:34 2026
```

Note the file names: `--prefix qc_` produces `qc__autoQC_Summary.txt`, with two underscores,
because the prefix is concatenated.

```
$ ls qc_*
qc__autoQC_Summary.txt
qc__autoQC_sampletoberemoved.txt
qc__autoQC_snptoberemoved.txt

$ cat qc__autoQC_sampletoberemoved.txt
FID	IID	REASON
MIS	M_C2	MissingMoreThan5
MIS	M_C3	MissingMoreThan5

$ head -4 qc__autoQC_snptoberemoved.txt; wc -l qc__autoQC_snptoberemoved.txt
SNP	REASON
rs1_1003166	CallRateLessThan80
rs1_1206939	CallRateLessThan80
rs1_1456164	CallRateLessThan80
    4899 qc__autoQC_snptoberemoved.txt

$ cut -f2 qc__autoQC_snptoberemoved.txt | tail -n +2 | sort | uniq -c
1569 CallRateLessThan80
 130 CallRateLessThan95
3199 Monomorphic
```

**How to use it.** `--autoQC` **does not write a filtered fileset** — it writes removal lists
with a reason per row. Hand them to PLINK:

```
plink --bfile study --exclude <(tail -n +2 study_autoQC_snptoberemoved.txt | cut -f1) \
                    --remove  <(tail -n +2 study_autoQC_sampletoberemoved.txt | cut -f1,2) \
                    --make-bed --out study_qc
```

A fourth file, `<prefix>_autoQC_updatesex.txt`, appears only when some sample's `.fam` sex is
`0`; it was not written here because `missing` declares every sex.

Thresholds are `--callrateN` (samples) and `--callrateM` (markers), both defaulting to 0.95.
Two things to know: the labels in the summary table (`< 80%`, `< 95%`) are **fixed text** and
do not track the thresholds you passed — the console lines above them do; and 3 199
"monomorphic" markers here is an artefact of having only six samples, not a property of the
array. Full rules, including several genuinely surprising boundary behaviours, are in
[VERIFIED_FORMULAS.md § `--autoQC`](VERIFIED_FORMULAS.md).

---

## 9. Reconstruct pedigrees from genotypes

**Goal:** you have genotypes and a `.fam` whose family structure is incomplete. Recover it.

`--build` does two things: it merges declared families that turn out to be genetically
connected, and it infers parent–child links inside the merged families.

```
$ open-king -b /tmp/kingdocs/bigish.bed --build
[... startup banner, loading, and the clustering block elided ...]
Pedigree reconstruction starts at Fri Aug 14 22:50:35 2026
Reconstructing pedigree...
Age information not provided.
Total length of 22 chromosomal segments usable for IBD segment analysis is 2498.9 Mb.
  Information of these chromosomal segments can be found in file kingallsegs.txt

Family KING1:
  Family KING1 RULE FS0: Sibship (B01_F B02_F)'s parents are (1 2)
Family KING2:
  Family KING2 RULE FS0: Sibship (B13_F B14_F)'s parents are (3 4)
Family KING3:
  Family KING3 RULE FS0: Sibship (B25_F B26_F)'s parents are (5 6)

Details of pedigree reconstruction are available in log file kingbuild.log
Update-ID information is saved in file kingupdateids.txt
Update-parent information is saved in file kingupdateparents.txt
Pedigree reconstruction ends at Fri Aug 14 22:50:35 2026
KING ends at Fri Aug 14 22:50:35 2026
```

```
$ head -6 kingupdateids.txt; wc -l kingupdateids.txt
BF01	B01_C1	KING1	B01_C1
BF01	B01_C2	KING1	B01_C2
BF01	B01_C3	KING1	B01_C3
BF01	B01_F	KING1	B01_F
BF01	B01_M	KING1	B01_M
BF02	B02_C1	KING1	B02_C1
      33 kingupdateids.txt

$ head -6 kingupdateparents.txt
KING1	B01_C1	B01_F	B01_M
KING1	B01_C2	B01_F	B01_M
KING1	B01_C3	B01_F	B01_M
KING1	B01_F	1	2
KING1	B01_M	0	0
KING1	B02_C1	B02_F	B02_M
```

**How to read it.**

* `updateids.txt` is `oldFID oldIID newFID newIID`, in PLINK `--update-ids` format. Families
  `BF01` and `BF02` are merged into a new family `KING1` because `B01_F` and `B02_F` turned
  out to be full siblings. 33 of the 200 samples are re-assigned.
* `updateparents.txt` is `FID IID father mother`, in PLINK `--update-parents` format. The
  interesting rows are `B01_F  1  2`: KING has invented two ungenotyped phantom parents
  (`1` and `2`) to hold the newly discovered sibship together. `B01_M`, who married in, keeps
  `0 0`.
* `build.log` records which rule fired for each inference. `RULE FS0` is "these people are
  full sibs, so give them a common pair of parents".

Feed the two files back to PLINK:

```
plink --bfile study --update-ids kingupdateids.txt --make-bed --out study_ids
plink --bfile study_ids --update-parents kingupdateparents.txt --make-bed --out study_ped
```

**Not every fileset reconstructs.** Of the thirteen corpus datasets only `bigish` produces
non-empty `--build` output; on the others KING prints `No families were found to be
connected` and `No pedigrees can be reconstructed`, and writes an empty `build.log` and no
`updateids.txt` — even though `multifam` contains undeclared cross-family parent–offspring
pairs. `--build` reconstructs *within* clusters it can form, and it needs the cluster to
exist first. (The reference binary behaves identically here; this was checked side by side.)

> **This is the one place open-king is known to differ from the reference in an output
> file.** `kingupdateids.txt` and `kingupdateparents.txt` are byte-identical in every case
> that writes them, and `build.log`'s header and `RULE` lines are byte-identical — but
> open-king does not emit the `INFERENCE` lines. Exactly what is missing is shown in
> [recipe 12](#12-migrating-from-the-original-king). If you consume `build.log` rather than
> the two update files, read [PARITY.md](PARITY.md) §6.2 first.

A related, lighter-weight command is `--cluster`, which does the family merge and writes a
`cluster.kin` of every pair inside the newly merged families, but does not infer parents.

---

## 10. X-chromosome relatedness

**Goal:** use the X chromosome — for sex-specific relationship checks, and as an independent
confirmation of an autosomal call.

Nothing extra is needed on the command line. If the `.bim` carries X markers, KING runs a
separate X pass and writes extra files. Load-time reporting tells you what it found:

```
$ open-king -b /tmp/kingdocs/sexchr.bed --kinship --cpus 1
Loading genotype data in PLINK binary format...
Read in PLINK fam file /tmp/kingdocs/sexchr.fam...
  PLINK pedigrees loaded: 10 samples
Read in PLINK bim file /tmp/kingdocs/sexchr.bim...
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs, 300 Y-chromosome SNPs, 50 mitochondrial SNPs
  PLINK maps loaded: 6000 SNPs
```

Note **chromosome 25 / `XY` (the pseudo-autosomal region) is pooled with the autosomes**,
while `23`/`X`, `24`/`Y` and `26`/`MT` are held aside. Any other chromosome code — `0`, `27`,
`chr1` — is dropped at map load and never counted.

```
$ open-king -b /tmp/kingdocs/sexchr.bed --kinship --cpus 1 --prefix x_
X-chromosome analysis...
X-chromosome genotypes stored in 24 64-bit words for each of 10 individuals.
Within-family kinship data saved in file x_X.kin
Relationship inference across families starts at Fri Aug 14 22:50:35 2026
                                         ends at Fri Aug 14 22:50:35 2026
Between-family kinship data saved in file x_X.kin0
KING ends at Fri Aug 14 22:50:35 2026

$ cat x_X.kin
FID	ID1	ID2	Sex	N_SNP	PhiX	Het	IBS0	KinshipX
SEX	S_DAU1	S_DAU2	FF	1500	0.3750	0.327	0.0000	0.3262
SEX	S_DAU1	S_F	FM	1500	0.5000	0.331	0.0000	0.5000
SEX	S_DAU1	S_M	FF	1500	0.2500	0.335	0.0000	0.2435
SEX	S_DAU1	S_SON1	FM	1500	0.2500	0.331	0.0973	0.2062
SEX	S_DAU1	S_SON2	FM	1500	0.2500	0.331	0.0540	0.3370
SEX	S_DAU2	S_F	FM	1500	0.5000	0.323	0.0000	0.5000
SEX	S_DAU2	S_M	FF	1500	0.2500	0.331	0.0000	0.2336
SEX	S_DAU2	S_SON1	FM	1500	0.2500	0.323	0.0247	0.4236
SEX	S_DAU2	S_SON2	FM	1500	0.2500	0.323	0.0493	0.3471
SEX	S_F	S_M	MF	1500	0.0000	0.339	0.1680	0.0049
SEX	S_F	S_SON1	MM	1500	0.0000	0.331	0.3227	-0.2238
SEX	S_F	S_SON2	MM	1500	0.0000	0.331	0.3153	-0.2017
SEX	S_M	S_SON1	FM	1500	0.5000	0.339	0.0000	0.5000
SEX	S_M	S_SON2	FM	1500	0.5000	0.339	0.0000	0.5000
SEX	S_SON1	S_SON2	MM	1500	0.5000	0.331	0.0793	0.5106
```

**How to read it.** The columns are not the autosomal ones. `Sex` is the pair's sex
combination (`FF`, `FM`, `MM`) and **there is a different estimator for each**; `PhiX` is the
pedigree-expected X kinship, `KinshipX` the estimate.

The X-specific facts are all visible in that table:

* **Father–daughter is `PhiX 0.5000`, and the estimate is `0.5000` exactly** — a daughter
  inherits her father's entire X. `S_DAU1`/`S_F` and `S_DAU2`/`S_F` both land on it.
* **Father–son shares nothing on the X**: `S_F`/`S_SON1` has `PhiX 0.0000` and a large IBS0
  (0.3227). A negative `KinshipX` for a declared father–son pair is expected, not an error.
* **Mother–son is `0.5000`**, again exactly (`S_M`/`S_SON1`, `S_M`/`S_SON2`).
* **Sisters are `PhiX 0.3750`**, not 0.25 — they always share their father's X.
* This makes the X a good **sex-assignment cross-check**: a declared father–son pair that
  reads 0.5 on the X is not a father and son.

The cross-family file has the same shape, and shows one further rule:

```
$ cat x_X.kin0
FID1	ID1	FID2	ID2	Sex	N_SNP	Het	IBS0	KinshipX
SEX	S_DAU1	SU3	S_UM	FM	1500	0.331	0.1707	-0.0151
SEX	S_DAU1	SU4	S_UF	FF	1500	0.348	0.0633	0.0105
SEX	S_DAU2	SU3	S_UM	FM	1500	0.323	0.1727	-0.0351
SEX	S_DAU2	SU4	S_UF	FF	1500	0.344	0.0627	0.0039
SEX	S_F	SU3	S_UM	MM	1500	0.331	0.3293	-0.2440
SEX	S_F	SU4	S_UF	MF	1500	0.365	0.1607	0.0602
SEX	S_M	SU3	S_UM	FM	1500	0.339	0.1793	-0.0285
SEX	S_M	SU4	S_UF	FF	1500	0.352	0.0533	0.0416
SEX	S_SON1	SU3	S_UM	MM	1500	0.331	0.3373	-0.2681
SEX	S_SON1	SU4	S_UF	MF	1500	0.365	0.1667	0.0438
SEX	S_SON2	SU3	S_UM	MM	1500	0.331	0.3380	-0.2701
SEX	S_SON2	SU4	S_UF	MF	1500	0.365	0.1640	0.0511
SU3	S_UM	SU4	S_UF	MF	1500	0.365	0.1627	0.0547
```

```
$ grep -c 'S_U0' x_.kin0; grep -c 'S_U0' x_X.kin0
17
0
```

**Samples of unknown sex are dropped from `--kinship`'s X analysis.** The two `SEX 0` samples
appear in 17 rows of the autosomal `.kin0` and in none of `X.kin0`. Fill in the sex column
(recipe 8's X-heterozygosity check tells you what to put there) before running the X pass.

### X IBD segments

```
$ open-king -b /tmp/kingdocs/sexchr.bed --ibdseg --degree 2 --cpus 1 --prefix xs_
  Genotype data consist of 4150 autosome SNPs (including 150 XY SNPs), 1500 X-chromosome SNPs, 300 Y-chromosome SNPs, 50 mitochondrial SNPs
  In addition to autosomes, 1 segments of length 75.0 Mb on X-chr can be further used.
Additional summary statistics of X-Chr IBD segments saved in file xs_X.seg

$ cat xs_X.seg
FID1	ID1	FID2	ID2	Sex1	Sex2	MaxIBD1	MaxIBD2	IBD1Seg	IBD2Seg	PropIBD
SEX	S_F	SEX	S_SON1	1	1	0.0000	0.0000	0.0000
SEX	S_F	SEX	S_SON2	1	1	0.0000	0.0000	0.0000
SEX	S_F	SEX	S_DAU1	1	2	1.0000	0.0000	0.5000
SEX	S_F	SEX	S_DAU2	1	2	1.0000	0.0000	0.5000
SEX	S_M	SEX	S_SON1	2	1	1.0000	0.0000	0.5000
SEX	S_M	SEX	S_SON2	2	1	1.0000	0.0000	0.5000
SEX	S_M	SEX	S_DAU1	2	2	1.0000	0.0000	0.5000
SEX	S_M	SEX	S_DAU2	2	2	1.0000	0.0000	0.5000
SEX	S_SON1	SEX	S_SON2	1	1	0.1462	0.6393	0.7124
SEX	S_SON1	SEX	S_DAU1	1	2	0.4257	0.0000	0.2128
SEX	S_SON1	SEX	S_DAU2	1	2	0.9067	0.0000	0.4533
SEX	S_SON2	SEX	S_DAU1	1	2	0.6397	0.0000	0.3199
SEX	S_SON2	SEX	S_DAU2	1	2	0.7245	0.0000	0.3623
SEX	S_DAU1	SEX	S_DAU2	2	2	0.6464	0.3530	0.6762
```

> **Read this file by position, not by header.** The header names eleven columns; every data
> row has ten tab-separated fields, the last of them empty:
>
> ```
> $ awk -F'\t' 'NR<=2{print NR": "NF" fields"}' xs_X.seg
> 1: 11 fields
> 2: 10 fields
> ```
>
> The three numbers written are **`IBD1Seg`, `IBD2Seg`, `PropIBD`**, sitting in the column
> positions the header labels `MaxIBD1`, `MaxIBD2` and `IBD1Seg`. This is a defect in the
> reference binary that open-king reproduces deliberately for byte parity. Check the
> arithmetic yourself on any row: `S_SON1`/`S_SON2` reads `0.1462  0.6393  0.7124`, and
> 0.6393 + 0.1462/2 = 0.7124 = `IBD2Seg + IBD1Seg/2`, which is what `PropIBD` means.

Two gates to know. `X.kin`/`X.kin0` need **512 or more X markers**, no `--degree`, and more
than one family. `X.seg` needs a usable X segment and a **non-zero `--degree`** — bare
`--ibdseg` writes no `X.seg` at all. And unlike `X.kin`, `X.seg` does *not* exclude samples
of unknown sex; it prints the raw `.fam` code in `Sex1`/`Sex2`.

---

## 11. Large cohorts: `--degree`, `--cpus`, and what actually costs time

The relatedness computation is O(N²) in samples: 200 people is 19 900 pairs, 20 000 people is
199 995 000. Two knobs matter.

### `--degree` — filter the output, and skip work

`--degree d` restricts `.kin0` to pairs with kinship ≥ 2<sup>-(d+1.5)</sup>. On `bigish`
(200 samples: 19 900 pairs, of which 19 327 cross a family boundary and so are `.kin0`'s):

```
$ open-king -b /tmp/kingdocs/bigish.bed --kinship --prefix d0_ >/dev/null
$ for d in 1 2 3 4; do open-king -b /tmp/kingdocs/bigish.bed --kinship --degree $d --prefix d$d\_ >/dev/null; done

$ for f in d0_ d1_ d2_ d3_ d4_; do
      printf '%-4s %6d rows  %8d bytes\n' $f $(($(wc -l < $f.kin0)-1)) $(wc -c < $f.kin0)
  done
d0_   19327 rows   1000337 bytes
d1_       3 rows       191 bytes
d2_      24 rows      1241 bytes
d3_      59 rows      3026 bytes
d4_      60 rows      3077 bytes
```

Unfiltered you get a megabyte of mostly-noise: 19 327 rows of unrelated pairs. `--degree 2`
gives you the 24 rows you actually wanted. On a 20 000-sample cohort the unfiltered file would
be about ten gigabytes of rows whose kinship is indistinguishable from zero. **Always pass
`--degree` on a large cohort**, and pick it from the question: 1 for duplicates and
first-degree, 2 for a GWAS relatedness screen, 3 for cryptic-relatedness discovery.

It also filters `.kin0` only — never `.kin`, never `.ibs`/`.ibs0` — and the comparison is
against the full-precision estimate, not the printed four decimals.

### `--cpus` — free, and does not change the answer

```
$ for c in 1 2 4 8 16; do open-king -b /tmp/kingdocs/bigish.bed --kinship --cpus $c --prefix c$c\_ >/dev/null; done

$ for c in 2 4 8 16; do for f in .kin .kin0; do
      cmp -s c1_$f c$c\_$f && echo "--cpus $c $f identical to --cpus 1" || echo "--cpus $c $f DIFFERS"
  done; done
--cpus 2 .kin identical to --cpus 1
--cpus 2 .kin0 identical to --cpus 1
--cpus 4 .kin identical to --cpus 1
--cpus 4 .kin0 identical to --cpus 1
--cpus 8 .kin identical to --cpus 1
--cpus 8 .kin0 identical to --cpus 1
--cpus 16 .kin identical to --cpus 1
--cpus 16 .kin0 identical to --cpus 1
```

Output is byte-identical across thread counts, so you can raise `--cpus` without worrying
about reproducibility. The default is every core the machine has.

### What things actually cost

```
N = 200 samples, 50 000 SNPs (bigish); each figure is the total of 20 runs

--kinship --cpus 1                      real 1.30 user 1.65 sys 0.15
--kinship --cpus 16                     real 1.32 user 1.64 sys 0.14
--kinship --degree 2 --cpus 16          real 1.17 user 1.44 sys 0.14
--related --degree 2 --ibdseg --cpus 16 real 14.52 user 15.05 sys 0.26
```

**Be careful reading this.** The largest fileset in the corpus is 200 × 50 000, which is small
enough that a whole `--kinship` run is about 65 ms and is dominated by loading the `.bed` and
writing the output — which is exactly why `--cpus 1` and `--cpus 16` are indistinguishable
here. This benchmark cannot tell you how the program scales; it can only tell you the relative
cost of the analyses, which is the useful part:

* **The IBD-segment engine is roughly an order of magnitude more expensive than kinship** —
  0.73 s versus 0.065 s per run here. That is where your time goes on a real cohort. If you
  only need kinship, do not ask for `--related` or `--ibdseg`.
* `--degree` shaves a little even at this size, and saves much more where the output is large.

Two more notes for big runs. The `--related` path automatically runs a two-stage screen (a
fast pass over ~32 768 informative SNPs, then an exhaustive re-estimate of the survivors) —
you can see it in [recipe 2](#2-find-cryptic-relatives-in-a-supposedly-unrelated-cohort)'s
console. And passing `--noscreen` on `bigish` changed neither the console flow nor the `.kin`
and `.kin0` bytes, so do not expect it to be a performance lever.

---

## 12. Migrating from the original KING

### What is identical

open-king is a clean-room reimplementation whose target is **byte-identical output**, not
statistical equivalence. Of 480 captured reference invocations across the 13 corpus datasets,
**all 480 reproduce byte for byte** — every output file, plus stdout, stderr and exit status.

Byte-identical everywhere: `--kinship` (including the X pass), `--duplicate`, `--ibs`,
`--unrelated`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster`, `--ibdseg` and `--related` at
all three captured `--seglength` floors, plus `X.kin`, `X.kin0` and `X.seg`. The command line,
the console banner, the warning and fatal frames, and the reference's own rounding and
row-ordering quirks are all reproduced.

The committed corpus currently has no differing invocation. This is a regression claim about
the recorded inputs, not a claim that every possible KING input is identical; the held-out
counterexamples are listed in [CONTINUATION.md](CONTINUATION.md#remaining-supported-core-work).

### What differs

**Out-of-scope analyses are not implemented.** Their recognized spellings fail before any
input is opened, with exit status 1 and a diagnostic that points to the product scope:

```
$ open-king -b /path/that/does/not/exist.bed --pca
FATAL ERROR - 
open-king's minimal relatedness product does not implement: --pca.
Supported analyses: --related, --duplicate, --kinship, --ibdseg, --ibs, --unrelated, --cluster, --build, --bysample, --bySNP, and --autoQC.
See docs/SCOPE.md for the product-scope contract.
$ echo $?
1
```

The reference binary runs the analysis and writes, for `--pca`, a `<prefix>pc.txt`. If your pipeline calls KING for
PCA, MDS, ROH, GRM, association testing (`--lmm`, `--tdt`, `--gdt`), risk scores, `--plink`
conversion, or the R plotting flags, keep the original binary for those steps.

The former 100 Mb segment-floor and conditional `splitped.txt` differences are fixed. The
remaining supported-core work is deliberately kept separate from the green corpus headline
and summarized in [CONTINUATION.md](CONTINUATION.md#remaining-supported-core-work).

Everything here is measured against **one reference build** — KING 2.3.2, Mach-O arm64, macOS.
KING's release notes record repeated changes to the (unpublished) IBD-segment algorithm across
2.1.x–2.2.x, so "byte-identical" means *to 2.3.2*.

### Run the parity suite yourself

```
$ python3 tests/parity/run_parity.py --impl target/release/open-king -q
[parity] 480 case(s), impl=target/release/open-king, jobs=8
parity: 480 PASS, 0 FAIL, 480 total (876 output file(s) byte-compared, 8 diff-excluded)
```

The wall-clock figure varies. It needs no reference binary — the goldens are committed. Point
`--impl` at the reference binary itself to prove the harness's normalization is sound.

### Diff the two binaries on your own data

The goldens cover the corpus. On your own fileset, run both and compare directly. Output
*files* compare with plain `cmp`; the console needs three things normalized first, because
they are properties of the machine and the moment, not of the analysis: wall-clock
timestamps, the CPU-core count, and the `\r`-terminated progress percentages.

```sh
#!/bin/sh
# scrub.sh - normalize a KING console log for diffing.
tr '\r' '\n' < "$1" \
  | grep -v '^[0-9][0-9]*%$' \
  | sed -E -e 's/(Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z][a-z] +[0-9]+ [0-9:]+ [0-9]{4}/<TS>/g' \
           -e 's/^([[:space:]]*)[0-9]+([[:space:]]*CPU cores are used)/\1<NCPU>\2/'
```

Use `--cpus 1` on both sides: the reference's `X.kin0` writer has a data race and produces a
different file on nearly every multi-threaded run (six identical reference runs gave six
different files, one truncated mid-number). `--cpus 1` makes it deterministic.

```
$ mkdir -p ref new
$ (cd ref && /path/to/king -b /tmp/kingdocs/multifam.bed --related --ibdseg --degree 2 --cpus 1 > stdout.txt 2>&1)
$ (cd new && open-king     -b /tmp/kingdocs/multifam.bed --related --ibdseg --degree 2 --cpus 1 > stdout.txt 2>&1)

$ diff -rq ref new
Files ref/stdout.txt and new/stdout.txt differ

$ ./scrub.sh ref/stdout.txt > ref.norm; ./scrub.sh new/stdout.txt > new.norm
$ diff ref.norm new.norm && echo 'console: identical after normalization'
console: identical after normalization

$ for f in king.kin king.seg kingallsegs.txt kingsplitped.txt; do
      cmp -s ref/$f new/$f && echo "$f: byte-identical" || echo "$f: DIFFERS"
  done
king.kin: byte-identical
king.seg: byte-identical
kingallsegs.txt: byte-identical
kingsplitped.txt: byte-identical
```

That is the normal result. Now the same procedure on the case that does not pass — `bigish`,
`--related --degree 2`:

```
$ mkdir -p g/ref g/new
$ (cd g/ref && /path/to/king -b /tmp/kingdocs/bigish.bed --related --degree 2 --cpus 1 > stdout.txt 2>&1)
$ (cd g/new && open-king     -b /tmp/kingdocs/bigish.bed --related --degree 2 --cpus 1 > stdout.txt 2>&1)
$ ./scrub.sh g/ref/stdout.txt > g/ref.norm; ./scrub.sh g/new/stdout.txt > g/new.norm

$ diff g/ref.norm g/new.norm && echo "console: byte-identical"
console: byte-identical

$ for f in king.kin king.kin0 kingallsegs.txt; do
      cmp -s g/ref/$f g/new/$f && echo "$f: byte-identical" || echo "$f: DIFFERS"
  done
king.kin: byte-identical
king.kin0: byte-identical
kingallsegs.txt: byte-identical
```

The fast screening count and every output byte now match. Check `--build` the same way:

```
$ mkdir -p b/ref b/new
$ (cd b/ref && /path/to/king -b /tmp/kingdocs/bigish.bed --build --cpus 1 >/dev/null 2>&1)
$ (cd b/new && open-king     -b /tmp/kingdocs/bigish.bed --build --cpus 1 >/dev/null 2>&1)

$ for f in kingupdateids.txt kingupdateparents.txt kingbuild.log; do
      cmp -s b/ref/$f b/new/$f && echo "$f: byte-identical" || echo "$f: DIFFERS"
  done
kingupdateids.txt: byte-identical
kingupdateparents.txt: byte-identical
kingbuild.log: byte-identical
```

This primary reconstruction includes the complete `INFERENCE` narration as well as the two
files a pipeline consumes. `PARITY.md` separately records rare held-out pedigree shapes that
the 480 captures do not exercise.

### Running both side by side

Nothing stops you. The binaries take the same arguments and write the same file names, so run
them in separate directories, or with different `--prefix` values, and diff as above. That is
how the whole project is tested.

---

## Traps

Things that cost people an afternoon. All of these are the reference binary's behaviour,
reproduced deliberately; several were verified against it side by side while writing this
page.

**A single-family dataset writes an empty `.kin`.** KING says it saved the file; the file is
zero bytes.

```
$ open-king -b /tmp/kingdocs/threegen.bed --kinship --prefix one_
Within-family kinship data saved in file one_.kin
There is only one family.

$ ls -la one_.kin
-rw-r--r--@ 1 wgu  wheel  0 Aug 14 22:51 one_.kin
```

The reference does exactly the same (checked). The cause is a buffer that is flushed every
64 KiB and never flushed at the end, so a one-family `.kin` is always truncated to whole
64 KiB chunks — under 64 KiB of rows means nothing reaches disk. With two or more distinct
FIDs the file is complete regardless of size. **Workarounds:** use `--ibs`, whose `.ibs` is
complete and carries the same `Kinship` column plus the raw counts —

```
$ open-king -b /tmp/kingdocs/threegen.bed --ibs --prefix ibs_
$ ls -la ibs_.ibs; head -2 ibs_.ibs | cut -f1-9
-rw-r--r--@ 1 wgu  wheel  8849 Aug 14 22:51 ibs_.ibs
FID	ID1	ID2	Z0	Phi	N_SNP	N_IBS0	N_IBS1	N_IBS2
TG	TG_C1	TG_C2	0.250	0.2500	20000	288	5044	14668
```

— or split the family across FIDs with `--fam`, which moves every pair into `.kin0`.

**`--related` silently becomes `--kinship` below 10 samples.**

```
$ open-king -b /tmp/kingdocs/trio.bed --related --prefix t3_
[... startup banner and loading block elided ...]

--related is replaced with --kinship for a small sample size.
Autosome genotypes stored in 79 words for each of 3 individuals.

Options in effect:
	--kinship

Within-family kinship data saved in file t3_.kin
```

Note that `Options in effect` now reads `--kinship`, and the `.kin` has 10 columns rather than
16.

**`--ibdseg` does the same below 5 samples** — observed at N = 2 and N = 3; N = 6 runs the
real segment path. The reference binary behaves identically:

```
$ open-king -b /tmp/kingdocs/pair.bed --ibdseg --prefix p2_
Options in effect:
	--kinship

Relationship inference across families starts at Fri Aug 14 22:51:43 2026
16 CPU cores are used.
                                         ends at Fri Aug 14 22:51:43 2026
Between-family kinship data saved in file p2_.kin0
```

**`--related` writes no `.kin0` below 100 samples** even when duplicates are present. See
[recipe 2](#if-your-cohort-has-fewer-than-100-samples-do-not-use---related-for-this).

**Nothing is filtered.** Outside `--autoQC`, every `.bim` record on chromosomes 1–22, 25 and
`XY` enters the computation in file order: no MAF filter, no call-rate filter, no monomorphic
filter, `0`-coded alleles and duplicate IDs and duplicate positions all kept. Missingness is
handled **pairwise** — `N_SNP` differs from row to row, and a marker called in only two
samples still counts for that pair. QC before you infer relatedness, or accept that a pair
with 4 000 shared calls is reported next to a pair with 50 000.

**The `Error` column is not an integer.** `0`, `0.5`, `1`. Parsing it with `%d` drops every
half-degree disagreement.

**`--prefix` concatenates.** `--prefix ZZ_` gives `ZZ_.kin` and `ZZ__autoQC_Summary.txt`. A
prefix pointing at a non-existent or unwritable directory is a fatal error *while the `.fam`
is being read*, not at output time.

**Missing input is a fatal error with exit 1**, and the message is on stdout:

```
$ open-king -b /tmp/kingdocs/nope.bed --kinship > e.out 2>&1; echo "exit=$?"; tail -4 e.out
exit=1

FATAL ERROR -
Genotype file /tmp/kingdocs/nope.bed cannot be opened
```

**Separators are inconsistent, on purpose.** `.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`,
`unrelated.txt` are **tab** separated; `bySample.txt` and `bySNP.txt` are **space** separated.

**Headers are dynamic.** `bySample.txt` has six header variants and `bySNP.txt` three,
depending on which chromosomes and pedigree structures the fileset has; `.ibs`/`.ibs0` gain
`MaxIBD2` and `Pr_IBD2` only when the map holds at least 100 Mb of usable segments; `.kin`
has 10 columns under `--kinship` and 16 under `--related`. **Always parse by header name.**

---

## Where to read next

* [README.md](README.md) — the documentation index.
* [CLI.md](CLI.md) — the complete command-line reference, including every flag this page
  does not use.
* [OUTPUTS.md](OUTPUTS.md) — every output file, every column, every format and row order.
* [INTERPRETING.md](INTERPRETING.md) — how to read the numbers, and the ways people misread
  them.
* [VERIFIED_FORMULAS.md](VERIFIED_FORMULAS.md) — every estimator, every column, every field
  format, each one checked numerically against the reference.
* [BEHAVIOR.md](BEHAVIOR.md) — which SNPs are used, when each file is created, how rows are
  ordered, and the experiments that established each rule.
* [SPEC.md](SPEC.md) — the implementation specification, for anyone changing the code.
* [PARITY.md](PARITY.md) — the authoritative parity claim, the analysis × dataset matrix, and
  the measured size of every remaining gap.
* [MAINTAINING.md](MAINTAINING.md) — regenerating the corpus, re-capturing goldens,
  contributing.
* `docs/research/` — 26 notes recording how each undocumented rule in the IBD-segment engine
  was recovered by black-box experiment.

**Cite the method, not this implementation.** The estimators are Manichaikul *et al.* 2010,
*Bioinformatics* 26(22):2867–2873, and KING itself is the work of Wei-Min Chen and colleagues.
