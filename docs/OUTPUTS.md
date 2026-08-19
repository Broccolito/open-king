# Output files — the complete reference

Every file `king` can write, what triggers it, what its columns mean, how each column is
formatted, and what order its rows come out in.

This page is for someone who has PLINK filesets and wants to read the results. It does not
explain the estimators — [`VERIFIED_FORMULAS.md`](VERIFIED_FORMULAS.md) does that — and it
is not the parity statement, which is [`PARITY.md`](PARITY.md). It documents open-king's
behaviour, which is byte-identical to KING 2.3.2 on all 480 captured reference
invocations; held-out differences are named in [§ Known divergences](#known-divergences-from-king-232)
and none of them affects a data column of any file below.

**Everything shown here was produced by running the binary.** Reproduce any of it:

```bash
cd /path/to/open-king
cargo build --release                                          # -> target/release/king
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs # 13 test filesets, ~20 s
```

Then run the command quoted above each listing, from an empty directory, with
`target/release/king` on your `PATH` as `king`. Output below is pasted verbatim; the
console excerpts start at `Options in effect:` and drop the load-time preamble (which
echoes your own input paths).

---

## Contents

* [Filenames: `--prefix` is concatenated, not joined](#filenames---prefix-is-concatenated-not-joined)
* [Which analysis writes which file](#which-analysis-writes-which-file)
* [Conventions that cut across every file](#conventions-that-cut-across-every-file)
  * [Separators](#separators)
  * [Numeric formats](#numeric-formats)
  * [Row order — three different rules](#row-order--three-different-rules)
  * [Proportions vs raw counts](#proportions-vs-raw-counts)
  * [`Z0` and `Phi` are pedigree expectations, not estimates](#z0-and-phi-are-pedigree-expectations-not-estimates)
  * [The `Error` column is graded 0 / 0.5 / 1](#the-error-column-is-graded-0--05--1)
  * [`InfType` — the inferred-relationship vocabulary](#inftype--the-inferred-relationship-vocabulary)
* **Relatedness files**
  * [`<prefix>.kin`](#prefixkin)
  * [`<prefix>.kin0`](#prefixkin0)
  * [`<prefix>.con`](#prefixcon)
  * [`<prefix>.ibs`](#prefixibs)
  * [`<prefix>.ibs0`](#prefixibs0)
  * [`<prefix>.seg`](#prefixseg)
  * [`<prefix>allsegs.txt`](#prefixallsegstxt)
  * [`<prefix>splitped.txt`](#prefixsplitpedtxt)
* **Sample-selection and pedigree files**
  * [`<prefix>unrelated.txt` and `<prefix>unrelated_toberemoved.txt`](#prefixunrelatedtxt-and-prefixunrelated_toberemovedtxt)
  * [`<prefix>cluster.kin`](#prefixclusterkin)
  * [`<prefix>updateids.txt`](#prefixupdateidstxt)
  * [`<prefix>updateparents.txt`](#prefixupdateparentstxt)
  * [`<prefix>build.log`](#prefixbuildlog)
* **QC files**
  * [`<prefix>bySample.txt`](#prefixbysampletxt)
  * [`<prefix>bySNP.txt`](#prefixbysnptxt)
  * [The four `<prefix>_autoQC_*` files](#the-four-prefix_autoqc_-files)
* **X-chromosome files**
  * [`<prefix>X.kin`](#prefixxkin)
  * [`<prefix>X.kin0`](#prefixxkin0)
  * [`<prefix>X.seg`](#prefixxseg)
* [Empty, header-only, truncated, absent](#empty-header-only-truncated-absent)
* [Known divergences from KING 2.3.2](#known-divergences-from-king-232)

---

## Filenames: `--prefix` is concatenated, not joined

`--prefix` defaults to `king`. It is **glued directly onto the suffix** — no separator is
inserted, so `--prefix king` gives `king.kin` *and* `kingallsegs.txt` *and*
`king_autoQC_Summary.txt`. A prefix ending in `_` therefore produces a double underscore
in front of the autoQC suffixes, and a prefix containing `/` writes into that directory
(which must already exist).

```bash
king -b /tmp/kingdocs/sexchr.bed --autoQC  --prefix ZZ_
king -b /tmp/kingdocs/sexchr.bed --kinship --prefix ZZ_
king -b /tmp/kingdocs/sexchr.bed --ibdseg --degree 2 --prefix ZZ_
king -b /tmp/kingdocs/sexchr.bed --kinship --prefix run.2026
mkdir -p sub && king -b /tmp/kingdocs/sexchr.bed --kinship --prefix sub/x_
```

```
ZZ_.kin
ZZ_.kin0
ZZ_.seg
ZZ_X.kin
ZZ_X.kin0
ZZ_X.seg
ZZ__autoQC_Summary.txt
ZZ__autoQC_sampletoberemoved.txt
ZZ__autoQC_snptoberemoved.txt
ZZ__autoQC_updatesex.txt
ZZ_allsegs.txt
ZZ_splitped.txt
run.2026.kin
run.2026.kin0
run.2026X.kin
run.2026X.kin0
sub
```

`sub/` holds `x_.kin`, `x_.kin0`, `x_X.kin`, `x_X.kin0`.

Note the shapes this produces: `ZZ_.kin` (prefix + `.kin`), `ZZ_X.kin` (prefix + `X.kin`),
`run.2026X.kin`, and the double underscore in `ZZ__autoQC_Summary.txt`. Throughout this
page `<prefix>` stands for whatever you passed, concatenated exactly like this.

## Which analysis writes which file

Each cell is the **union** of what that analysis can write, and every entry was observed in
a real run's directory listing. No analysis writes all of its files on every input — the
conditions are stated per file below and collected in
[Empty, header-only, truncated, absent](#empty-header-only-truncated-absent).

| flag | files written |
| --- | --- |
| `--kinship` | `.kin`, `.kin0`, and `X.kin`/`X.kin0` when the X pass runs |
| `--related` | `.kin`, `.kin0`, `X.kin`, `X.kin0`, `allsegs.txt` |
| `--duplicate` | `.con` |
| `--ibs` | `.ibs`, `.ibs0`, `allsegs.txt` |
| `--ibdseg` | `.seg`, `X.seg`, `allsegs.txt`, `splitped.txt` |
| `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt`, `allsegs.txt` |
| `--cluster` | `cluster.kin`, `updateids.txt`, `allsegs.txt` |
| `--build` | `build.log`, `updateparents.txt`, `updateids.txt`, `allsegs.txt` |
| `--bysample` | `bySample.txt`, `allsegs.txt` |
| `--bySNP` | `bySNP.txt`, `allsegs.txt` |
| `--autoQC` | `_autoQC_Summary.txt`, `_autoQC_snptoberemoved.txt`, `_autoQC_sampletoberemoved.txt`, `_autoQC_updatesex.txt` |

`--kinship` and `--duplicate` are the only analyses that do **not** run the usable-segment
pre-pass, so they are the only ones that never write `allsegs.txt`. `--autoQC` is a
standalone QC pass and writes only its own four reports.

## Conventions that cut across every file

### Separators

This asymmetry is real, is the reference's, and will bite a naive parser:

| separator | files |
| --- | --- |
| **tab** | `.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`, `.seg`, `X.kin`, `X.kin0`, `X.seg`, `cluster.kin`, `allsegs.txt`, `unrelated.txt`, `unrelated_toberemoved.txt`, `updateids.txt`, `updateparents.txt`, `_autoQC_snptoberemoved.txt`, `_autoQC_sampletoberemoved.txt`, `_autoQC_updatesex.txt` |
| **space** | `bySample.txt`, `bySNP.txt`, `splitped.txt` |
| **fixed-width, space padded** | `_autoQC_Summary.txt` (contains no tab at all) |
| free text | `build.log` |

All files are `\n` terminated and have a header line **except** `unrelated.txt`,
`unrelated_toberemoved.txt`, `updateids.txt`, `updateparents.txt`, `splitped.txt` and
`_autoQC_updatesex.txt`, which have none.

### Numeric formats

| format | columns |
| --- | --- |
| `%d` | every `N_*` / `N` count, `N_SNP`, `Pos`, `Chr` (numeric codes), `Sex1`/`Sex2`/`SEX` (the `Sex` column of `--kinship`'s `X.kin`/`X.kin0` is text, not a number — see below) |
| `%.3f` | `Z0`; `Het` in `X.kin`/`X.kin0`; `MaxIBD2` (a base-pair count, printed with three decimals) |
| `%.4f` | `Phi`, `PhiX`, `HetHet`, `IBS0`, `Kinship`, `KinshipX`, `HetConc`, `HomIBS0`, `IBS`, `Dist`, `Het2\|1`, `Het1\|2`, `HomConc`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `Pr_IBD2`, `Freq_A`, `CallRate`, `Missing`, `Heterozygosity`, the `Err_*` rates |
| `%.5f` | `Concord`, `HomConc`, `HetConc` in `.con` — the only five-decimal fields anywhere |
| shortest round-trip | `Error` (see below) — prints `0`, `0.5` or `1`, **never** `%d` |
| literal `-9` | `MaxIBD2`/`Pr_IBD2` in `.ibs0` for a pair below the analysis gate |
| literal `nan` | a column whose denominator is zero, e.g. `HomIBS0` for a pair with no A1 homozygote on either side |

`MaxIBD2` is a length in **base pairs**, not megabases: `44749371.000` in the example
below is 44.7 Mb. `Pr_IBD2` is a proportion of the usable genome.

### Row order — three different rules

None of the three is `.fam` order, and two of them coincide with the obvious guess on
small inputs, which is exactly what makes them dangerous.

**1. Within-family files (`.kin`, `.ibs`, `X.kin`) — sorted, not `.fam` order.**
Families are ordered by a natural-sort comparator on the FID; inside a family, members are
ordered by the same comparator on the IID; rows are the `i < j` upper triangle of that
sorted member list. The comparator compares runs of digits by **length first, then bytes**
(so `7 < 70 < 007`), folds non-digits to uppercase, and sorts a non-digit before a digit.

Demonstrated by renaming one family's members to numeric IDs whose `.fam` order,
lexicographic order and numeric order all differ:

```bash
python3 -c "
m={'A_F':'10','A_M':'9','A_C1':'2','A_C2':'007','A_C3':'70'}
out=[]
for l in open('/tmp/kingdocs/multifam.fam'):
    f=l.split()
    if f[0]=='FAM1':
        f[1]=m.get(f[1],f[1]); f[2]=m.get(f[2],f[2]); f[3]=m.get(f[3],f[3])
    out.append(' '.join(f))
open('sort.fam','w').write('\n'.join(out)+'\n')"
king -b /tmp/kingdocs/multifam.bed --fam sort.fam --kinship --prefix s
```

`.fam` lists the family as `10, 9, 2, 007, 70`:

```
FAM1 10 0 0 1 -9
FAM1 9 0 0 2 -9
FAM1 2 10 9 1 -9
FAM1 007 10 9 2 -9
FAM1 70 10 9 1 -9
```

`s.kin` emits it as `2, 9, 10, 70, 007` (first three columns of the `FAM1` rows, via
`awk -F'\t' '$1=="FAM1"{print $1"\t"$2"\t"$3}' s.kin`):

```
FAM1	2	9
FAM1	2	10
FAM1	2	70
FAM1	2	007
FAM1	9	10
FAM1	9	70
FAM1	9	007
FAM1	10	70
FAM1	10	007
FAM1	70	007
```

Lexicographic order would have given `007, 10, 2, 70, 9`; plain `.fam` order would have
given `10, 9, 2, 007, 70`. Neither is what comes out.

**2. Between-family files (`.kin0`, `.ibs0`) — block-tiled over `.fam` index order.**
Pairs are `i < j` over `.fam` row index, but sorted by the key `(i/B, j/B, i, j)` with
integer division — a square-tiled walk, not row-major. **`B` differs per file: 32 for
`.kin0`, 8 for `.ibs0`.**

This coincides with plain ascending `(i, j)` whenever `n <= B`, so it looks right on small
files and silently diverges on real data. Checked against a 200-sample run
(`--ibs --cpus 4` and `--kinship --cpus 4` on `bigish`), by regenerating each candidate
order and comparing to the emitted rows:

```
.ibs0   19327 rows   block size that reproduces the order: 8   (rejected: 1, 16, 32, 64)
.kin0   19327 rows   block size that reproduces the order: 32  (rejected: 1, 8, 16, 64)
.seg      442 rows   block size that reproduces the order: 16  (rejected: 1, 8, 32, 64)
```

**3. `.seg` and `X.seg` — 16-sample blocks.** As the third line above shows: for each
block `b1`, for each block `b2 >= b1`, every reported pair with `i` in `b1` and `j` in
`b2`, in index order. Blocks of 16 samples, and 16 uniquely (swept 2..80).

Other files: `.con` is in serial `i < j` `.fam` order; `bySample.txt` is in `.fam` order;
`bySNP.txt` is in `.bim` order within each chromosome class (autosomes and `XY` first, then
X, then Y, then MT); `allsegs.txt` is in map order; `cluster.kin` is by cluster, then the
sorted-member upper triangle; `updateids.txt` is by original `(FID, IID)` under the ID
comparator while `updateparents.txt` is in cluster order.

### Proportions vs raw counts

**In `.kin`, `.kin0`, `cluster.kin` the `HetHet` and `IBS0` columns are proportions of
`N_SNP`. In `.ibs` and `.ibs0` the same quantities appear as raw counts,** under the names
`NHetHet` and `N_IBS0`. Verified across the two files of the same dataset:

```
FAM1 A_C1/A_C2: .ibs NHetHet=3329 N_IBS0=230 N_SNP=15000 -> 3329/15000=0.2219 (.kin HetHet 0.2219), 230/15000=0.0153 (.kin IBS0 0.0153)
FAM1 A_C1/A_C3: .ibs NHetHet=3257 N_IBS0=321 N_SNP=15000 -> 3257/15000=0.2171 (.kin HetHet 0.2171), 321/15000=0.0214 (.kin IBS0 0.0214)
FAM1 A_C1/A_F:  .ibs NHetHet=2602 N_IBS0=0   N_SNP=15000 -> 2602/15000=0.1735 (.kin HetHet 0.1735), 0/15000=0.0000   (.kin IBS0 0.0000)
```

`N_SNP` itself is **pairwise**: the number of markers called in *both* members of the pair.
It differs from row to row whenever your data has missingness.

### `Z0` and `Phi` are pedigree expectations, not estimates

`Z0` (`%.3f`) and `Phi` (`%.4f`) in `.kin` and `.ibs` are computed **from the declared
pedigree in the `.fam`**, not from the genotypes. A parent–offspring pair is
`Z0 0.000  Phi 0.2500`; full siblings are `Z0 0.250  Phi 0.2500`; a within-family pair the
pedigree says is unrelated is `Z0 1.000  Phi 0.0000`. If your `.fam` has no parent columns,
every `Phi` is `0.0000` and every `Z0` is `1.000` — that is a statement about your input,
not about your samples.

`.kin0` has no `Z0`/`Phi`: by construction its pairs are in different families, so the
pedigree says nothing about them.

### The `Error` column is graded 0 / 0.5 / 1

`Error` (last column of `.kin`) measures disagreement between what the genotypes infer and
what the pedigree declares. Its value set is exactly `{0, 0.5, 1}` and it is **not an
integer** — a parser that reads it as `int` will crash or truncate.

* `0` — inference agrees with the pedigree
* `0.5` — off by one degree
* `1` — off by more than one degree, or a first-degree mismatch (see below)

**`--kinship` and `--related` grade the same pair by different rules and can disagree.**
`--kinship` compares the kinship *estimate* against `Phi` multiplicatively (within a factor
of √2 is `0`, within a factor of 2 is `0.5`, beyond that `1`). `--related` compares the
pedigree's relationship *label* against the segment-based `InfType`, which makes an exact
label mismatch inside the first degree score `1`.

A real pair showing both:

```bash
king -b /tmp/kingdocs/monomorphic.bed --kinship --prefix k
king -b /tmp/kingdocs/monomorphic.bed --related --prefix r
```

`k.kin`, filtered to its non-zero `Error` rows with
`awk -F'\t' 'NR==1||$10=="0.5"' k.kin` — the pedigree says full sibs, the estimate 0.1477
is within a factor of 2 of `Phi` 0.25 but not within √2, so `0.5`:

```
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	Kinship	Error
MONO	P_C3	P_C4	5000	0.250	0.2500	0.1064	0.0204	0.1477	0.5
```

`r.kin`, filtered the same way with `awk -F'\t' 'NR==1||$16=="1"' r.kin` — the same pair
now grades `1`, because its `InfType` `2nd` is not the pedigree's label `FS`:

```
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType	Error
MONO	P_C1	P_C2	5000	0.250	0.2500	0.1576	0.0008	0.5194	0.0284	0.3384	0.9800	0.0000	0.4900	PO	1
MONO	P_C1	P_C3	5000	0.250	0.2500	0.1330	0.0044	0.3951	0.1146	0.2645	0.9007	0.0000	0.4504	PO	1
MONO	P_C3	P_C4	5000	0.250	0.2500	0.1064	0.0204	0.3152	0.4513	0.1477	0.4812	0.0000	0.2406	2nd	1
```

(`monomorphic` is a deliberately degenerate fixture — most of its markers carry no
information — which is why its inferences are wrong. That is the point of it.)

`.kin0` never carries an `Error` column: with no shared FID there is no pedigree
expectation to compare against.

### `InfType` — the inferred-relationship vocabulary

Seven values, observed across every `.seg`, `.kin` and `.kin0` this page produced:

```
Dup/MZ   PO   FS   2nd   3rd   4th   UN
```

`Dup/MZ` is a duplicate sample or MZ twin pair; `PO` parent–offspring; `FS` full siblings;
`2nd`/`3rd`/`4th` degree; `UN` unrelated. `PO` and `FS` are separated by IBD2 sharing when
segment analysis is available (a true `PO` pair is `IBD1Seg 1.0000  IBD2Seg 0.0000`), and
by an IBS0 cutoff otherwise.

One subtlety worth knowing if you cross-reference files: **`.kin`'s `InfType` is not
`.seg`'s.** `.kin`'s `Dup/MZ` clause additionally requires `HetConc > 0.8`, so the same
pair can print `Dup/MZ` in `.seg` and `FS` in `.kin`.

---

# Relatedness files

## `<prefix>.kin`

Within-family pairwise relatedness. Two shapes: **10 columns under `--kinship`**, **16
columns under `--related`**.

**Written when** at least one family in the `.fam` has two or more members. Not created at
all when every sample is its own family. **Truncated** — see below — when the dataset
contains exactly one distinct FID.

### `--kinship` — 10 columns

```bash
king -b /tmp/kingdocs/multifam.bed --kinship --prefix mf
```

```
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	Kinship	Error
FAM1	A_C1	A_C2	15000	0.250	0.2500	0.2219	0.0153	0.2721	0
FAM1	A_C1	A_C3	15000	0.250	0.2500	0.2171	0.0214	0.2476	0
FAM1	A_C1	A_F	15000	0.000	0.2500	0.1735	0.0000	0.2475	0
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1 | `FID` | family ID, shared by both members | text |
| 2 | `ID1` | first member's IID | text |
| 3 | `ID2` | second member's IID | text |
| 4 | `N_SNP` | markers called in **both** members | `%d` |
| 5 | `Z0` | pedigree Pr[IBD = 0] | `%.3f` |
| 6 | `Phi` | pedigree kinship coefficient | `%.4f` |
| 7 | `HetHet` | fraction of `N_SNP` at which both are heterozygous | `%.4f` |
| 8 | `IBS0` | fraction of `N_SNP` at which they are opposite homozygotes | `%.4f` |
| 9 | `Kinship` | estimated kinship, within-family estimator `(HetHet − 2·IBS0)/(Het₁ + Het₂)` | `%.4f` |
| 10 | `Error` | pedigree-vs-inference grade, `0` / `0.5` / `1` | not `%d` |

**How to read row 1.** `A_C1` and `A_C2` are declared full sibs (`Z0 0.250`, `Phi 0.2500`).
Over the 15 000 markers both were called at, 22.19 % are het in both and 1.53 % are
opposite homozygotes. The estimate 0.2721 is a first-degree value and agrees with the
pedigree's 0.25, so `Error` is `0`. Row 3 is a parent–offspring pair: `Z0 0.000` because
PO pairs share an allele IBD at every locus, and `IBS0 0.0000` because opposite homozygotes
are impossible for a true PO pair without genotyping error.

Interpretation cutoffs (from the KING paper): `> 0.354` duplicate/MZ, `0.177–0.354` 1st
degree, `0.0884–0.177` 2nd, `0.0442–0.0884` 3rd, `0.0221–0.0442` 4th, `< 0.0221`
unrelated. The boundaries are successive halvings on the `2^(−k/2)` grid.

### `--related` — 16 columns

```bash
king -b /tmp/kingdocs/bigish.bed --related --cpus 4
```

```
FID	ID1	ID2	N_SNP	Z0	Phi	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType	Error
BF01	B01_C1	B01_C2	50000	0.250	0.2500	0.2061	0.0158	0.4202	0.1357	0.2505	0.5328	0.2569	0.5233	FS	0
BF01	B01_C1	B01_C3	50000	0.250	0.2500	0.2056	0.0149	0.4208	0.1288	0.2533	0.5157	0.2587	0.5166	FS	0
BF01	B01_C1	B01_F	50000	0.000	0.2500	0.1724	0.0000	0.3296	0.0000	0.2479	1.0000	0.0000	0.5000	PO	0
```

Columns 1–8 are the `--kinship` set. The six extra ones:

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 9 | `HetConc` | heterozygote concordance, `HetHet / (HetHet + IBS1)` | `%.4f` |
| 10 | `HomIBS0` | `IBS0` over the count of markers at which **either** sample is homozygous for A1 — a union, and A1 specifically. Not `IBS0/HomHom`, not `1 − HomConc`. `nan` when no A1 homozygote exists on either side | `%.4f` |
| 11 | `Kinship` | as above | `%.4f` |
| 12 | `IBD1Seg` | proportion of the usable genome shared IBD1, from called segments | `%.4f` |
| 13 | `IBD2Seg` | proportion shared IBD2 | `%.4f` |
| 14 | `PropIBD` | `IBD2Seg + IBD1Seg/2`, at full precision | `%.4f` |
| 15 | `InfType` | inferred relationship (see the vocabulary above) | text |
| 16 | `Error` | grade — under `--related` this compares *labels*, not kinships | not `%d` |

**How to read row 3.** `B01_C1` and `B01_F` are child and father. `IBD1Seg 1.0000` and
`IBD2Seg 0.0000` — a parent and child share exactly one allele IBD along the whole genome
and never two — so `PropIBD` is exactly `0.5000` and `InfType` is `PO`. Rows 1 and 2 are
sibs: about half the genome IBD1 and a quarter IBD2, `PropIBD ≈ 0.52`, `FS`.

> **`PropIBD` in `.kin` and in `.seg` can differ in the last digit for the same pair.**
> This is the reference's own inconsistency, reproduced deliberately: `.kin` computes
> `IBD2Seg + IBD1Seg/2` at full precision, while `.seg` recombines its own two *printed*
> four-decimal columns. See [`.seg`](#prefixseg).

## `<prefix>.kin0`

Between-family pairwise relatedness — every pair whose two members are in different
families. **8 columns under `--kinship`**, **14 columns under `--related`**.

**Written when:**

* under `--kinship`: at least two distinct FIDs exist. Never truncated.
* under `--related`: only when the between-family stage detects at least one candidate
  pair. At `--degree 1` or `2` (including the default) that stage additionally requires
  **N ≥ 100 samples** — below that the console prints `No close relatives are inferred.`
  and no `.kin0` appears, even if your data contains duplicates. At `--degree 3` and above
  the screening stage is skipped and the sample-count gate does not apply.

Measured:

| run | samples | `.kin0`? |
| --- | ---: | --- |
| `--related` on `multifam` | 20 | no — `No close relatives are inferred.` |
| `--related --degree 3` on `multifam` | 20 | **yes**, 54 rows |
| `--related` on 110 unrelated singletons † | 110 | no — nothing detected |
| `--related` on `bigish` | 200 | **yes**, 3 rows |

Count the rows in the file rather than trusting the console: the `N pairs … are identified`
line reports the relationship-summary tally, which never increments its `4th` column, so
the `--degree 3` run above says `52` while the file holds 54.

† A constructed fileset: no corpus dataset combines ≥ 100 samples with no relatives. Built
by importing `tests/parity/generate_corpus.py` as a module and simulating 110 singleton
families over five autosomes — the same generator that produces the 13 corpus datasets.
[`X.kin0`](#prefixxkin0) below shows the full recipe for a fileset built that way.

### `--kinship` — 8 columns

```bash
king -b /tmp/kingdocs/multifam.bed --kinship --prefix mf
```

```
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	Kinship
FAM1	A_F	FAM2	B_F	15000	0.2169	0.0204	0.2501
FAM1	A_F	FAM2	B_M	15000	0.1386	0.0685	0.0001
FAM1	A_F	FAM2	B_C1	15000	0.1536	0.0301	0.1319
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1–2 | `FID1`, `ID1` | first sample | text |
| 3–4 | `FID2`, `ID2` | second sample, in a different family | text |
| 5 | `N_SNP` | markers called in both | `%d` |
| 6 | `HetHet` | proportion of `N_SNP` het in both | `%.4f` |
| 7 | `IBS0` | proportion of `N_SNP` opposite-homozygous | `%.4f` |
| 8 | `Kinship` | **between-family (KING-robust) estimator**, `0.5 + (2·HetHet − 4·IBS0 − Het₁ − Het₂) / (4·min(Het₁, Het₂))` | `%.4f` |

**The estimator is not the same one `.kin` uses.** The `min()` denominator is what makes
it robust to population structure; the two forms coincide only when the two samples have
equal heterozygosity. Which form applies is decided purely by whether the pair shares an
FID.

**How to read these rows.** `A_F` and `B_F` are in different declared families but come
out at kinship 0.2501 — a cryptic first-degree relationship, which is exactly what `.kin0`
exists to find. `A_F`/`B_M` at 0.0001 with 6.85 % opposite homozygotes is an unrelated
pair. `A_F`/`B_C1` at 0.1319 is second-degree, consistent with being the child of a
first-degree relative.

### `--related` — 14 columns

```bash
king -b /tmp/kingdocs/bigish.bed --related --cpus 4
```

```
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
BF01	B01_F	BF02	B02_F	50000	0.2274	0.0132	0.4855	0.1175	0.2885	0.4575	0.3676	0.5964	FS
BF13	B13_F	BF14	B14_F	50000	0.2304	0.0116	0.4911	0.1045	0.2959	0.4484	0.3767	0.6009	FS
BF25	B25_F	BF26	B26_F	50000	0.2099	0.0169	0.4288	0.1467	0.2504	0.5064	0.2682	0.5214	FS
```

Same columns as the 16-column `.kin` minus `Z0`, `Phi` and `Error` (all three pedigree
quantities, all three meaningless across families), plus the second FID.

**Rows are filtered, unlike `.kin`.** A pair appears only if
`Kinship >= 2^-(d+1.5)` **or** `PropIBD > 2^-(d+0.5)`, where `d` is `--degree` (default 1).
That is why the file above has three rows out of 19 327 possible pairs.

### `--degree` filters `.kin0` and nothing else

```bash
king -b /tmp/kingdocs/multifam.bed --kinship            --prefix a_
king -b /tmp/kingdocs/multifam.bed --kinship --degree 2 --prefix b_
echo "a_ rows $(( $(wc -l < a_.kin0) - 1 ));  b_ rows $(( $(wc -l < b_.kin0) - 1 ))"
```

```
a_ rows 150;  b_ rows 32
```

All 150 cross-family pairs are written unfiltered; `--degree 2` keeps the 32 whose kinship
reaches `2^-3.5 = 0.08839`. The comparison is against the exact double, not the printed
four decimals. `--degree` never filters `.kin`, `.ibs` or `.ibs0`, and `--degree 0` means
"unset".

## `<prefix>.con`

Duplicate/MZ pairs found by `--duplicate`. Eleven tab-separated columns; both within- and
cross-family pairs, in serial `i < j` `.fam` order.

**Written when:** always for N < 100 samples — header-only if nothing passes the threshold.
For **N ≥ 100 only if at least one duplicate is found**; otherwise no file at all.

```bash
king -b /tmp/kingdocs/dups.bed --duplicate
```

```
FID1	ID1	FID2	ID2	N	N_IBS0	N_IBS1	N_IBS2	Concord	HomConc	HetConc
DUPA	DUP_A	DUPB	DUP_A_COPY	10000	0	0	10000	1.00000	1.00000	1.00000
MZFAM	MZ_1	MZFAM	MZ_2	10000	9	17	9974	0.99740	0.99862	0.99512
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1–4 | `FID1`, `ID1`, `FID2`, `ID2` | the pair | text |
| 5 | `N` | markers called in both (note: `N`, not `N_SNP`) | `%d` |
| 6 | `N_IBS0` | opposite homozygotes | `%d` |
| 7 | `N_IBS1` | sharing exactly one allele | `%d` |
| 8 | `N_IBS2` | identical genotypes | `%d` |
| 9 | `Concord` | overall genotype concordance, `N_IBS2 / N` | `%.5f` |
| 10 | `HomConc` | homozygote concordance, `(HomHom − IBS0) / HomHom` | `%.5f` |
| 11 | `HetConc` | heterozygote concordance, `HetHet / (Het₁ + Het₂ − HetHet)` | `%.5f` |

**How to read this.** Row 1 is an exact duplicate — 10 000 of 10 000 genotypes identical.
Row 2 is an MZ twin pair with 26 discordant genotypes, most of them genotyping error:
9 opposite homozygotes and 17 one-allele differences.

`--minConc` (default 0.80) is the threshold on `HetConc`, applied **strictly**:

```bash
king -b /tmp/kingdocs/dups.bed --duplicate --minConc 0.999 --prefix c
```

drops the MZ pair (`HetConc 0.99512`) and keeps only the exact duplicate.

The header-only case — 65 bytes, no data rows:

```bash
king -b /tmp/kingdocs/unrelated.bed --duplicate   # 30 samples, no duplicates
```

```
FID1	ID1	FID2	ID2	N	N_IBS0	N_IBS1	N_IBS2	Concord	HomConc	HetConc
```

and the absent case — `--duplicate` on the 200-sample `bigish`, no duplicates, writes no
`.con` at all.

## `<prefix>.ibs`

Within-family IBS and concordance statistics from `--ibs`. **This is where the raw counts
live**; `.kin` carries the same quantities as proportions.

**Written when:** always — even a fileset with no within-family pair gets a header-only
`.ibs` (139 bytes in the long form). Never truncated.

**Column set varies.** `MaxIBD2` and `Pr_IBD2` are appended **iff the total usable segment
length for IBD analysis is at least 100 Mb** — the same total the console reports as
`Total length of <n> chromosomal segments usable for IBD segment analysis is <x> Mb.`

```bash
king -b /tmp/kingdocs/multifam.bed --ibs      # 691.5 Mb usable -> long form
```

```
FID	ID1	ID2	Z0	Phi	N_SNP	N_IBS0	N_IBS1	N_IBS2	NHetHet	NHomHom	N_Het1	N_Het2	IBS	Dist	HetConc	Het2|1	Het1|2	HomConc	Kinship	MaxIBD2	Pr_IBD2
FAM1	A_C1	A_C2	0.250	0.2500	15000	230	3885	10885	3329	7786	5257	5286	1.7103	0.3203	0.4615	0.6333	0.6298	0.9705	0.2721	44749371.000	0.2401
FAM1	A_C1	A_C3	0.250	0.2500	15000	321	4046	10633	3257	7697	5257	5303	1.6875	0.3553	0.4460	0.6196	0.6142	0.9583	0.2476	31947441.000	0.2124
FAM1	A_C1	A_F	0.000	0.2500	15000	0	5310	9690	2602	7088	5257	5257	1.6460	0.3540	0.3289	0.4950	0.4950	1.0000	0.2475	0.000	0.0000
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1–3 | `FID`, `ID1`, `ID2` | the within-family pair | text |
| 4 | `Z0` | pedigree Pr[IBD = 0] | `%.3f` |
| 5 | `Phi` | pedigree kinship | `%.4f` |
| 6 | `N_SNP` | markers called in both | `%d` |
| 7 | `N_IBS0` | **count** of opposite homozygotes | `%d` |
| 8 | `N_IBS1` | count sharing exactly one allele, `Het₁ + Het₂ − 2·HetHet` | `%d` |
| 9 | `N_IBS2` | count of identical genotypes, `N_SNP − N_IBS0 − N_IBS1` | `%d` |
| 10 | `NHetHet` | **count** of markers het in both | `%d` |
| 11 | `NHomHom` | count of markers homozygous in both | `%d` |
| 12 | `N_Het1` | `ID1`'s heterozygote count **over the pairwise non-missing set** | `%d` |
| 13 | `N_Het2` | `ID2`'s, likewise | `%d` |
| 14 | `IBS` | mean IBS allele sharing, `(N_IBS1 + 2·N_IBS2)/N_SNP`, in `[0, 2]` | `%.4f` |
| 15 | `Dist` | mean squared genotype distance, `(N_IBS1 + 4·N_IBS0)/N_SNP`. **Not** `2 − IBS` unless `IBS0 = 0` | `%.4f` |
| 16 | `HetConc` | `NHetHet / (N_Het1 + N_Het2 − NHetHet)` | `%.4f` |
| 17 | `Het2\|1` | `NHetHet / N_Het1` | `%.4f` |
| 18 | `Het1\|2` | `NHetHet / N_Het2` | `%.4f` |
| 19 | `HomConc` | `(NHomHom − N_IBS0) / NHomHom` | `%.4f` |
| 20 | `Kinship` | within-family estimator | `%.4f` |
| 21 | `MaxIBD2` | longest IBD2 segment, **in base pairs** | `%.3f` |
| 22 | `Pr_IBD2` | proportion of the usable genome shared IBD2 | `%.4f` |

`N_Het1` and `N_Het2` are counted over the *pairwise* non-missing set, not over each
sample's own. With uneven missingness the two differ, and using a per-sample count gives
subtly wrong kinship.

**How to read row 1.** Full sibs: 230 of 15 000 markers opposite-homozygous, 10 885
identical, mean IBS 1.71. Their longest shared IBD2 stretch is 44.7 Mb and 24.0 % of the
usable genome is IBD2 — both hallmarks of a sib pair. Row 3, the PO pair, has
`N_IBS0 = 0`, `HomConc 1.0000` and `MaxIBD2 0.000`: no IBD2 anywhere, which is what
distinguishes PO from FS at the same kinship of 0.25.

**In `.ibs`, a pair below the segment gate prints `0.000 / 0.0000`, not `-9`** — that is
`.ibs0`'s spelling. The gate is `Kinship >= 2^-3.5 = 0.08839`; `A_C1`/`A_F` above is *above*
the gate and genuinely has no IBD2.

The short form, when the map has less than 100 Mb of usable segments:

```bash
king -b /tmp/kingdocs/pair.bed --ibs     # 42.6 Mb usable -> "Segments too short."
```

```
FID	ID1	ID2	Z0	Phi	N_SNP	N_IBS0	N_IBS1	N_IBS2	NHetHet	NHomHom	N_Het1	N_Het2	IBS	Dist	HetConc	Het2|1	Het1|2	HomConc	Kinship
```

(header-only here, 123 bytes, because `pair`'s two samples are in different families).

## `<prefix>.ibs0`

Between-family IBS statistics from `--ibs`. Same statistic columns as `.ibs`, with
`FID1 ID1 FID2 ID2` in place of `FID ID1 ID2` and **no `Z0`/`Phi`**, and the `Kinship`
column carrying the between-family estimator.

**Written when:** at least two distinct FIDs exist. Row order is block-tiled at **B = 8**.
The `MaxIBD2`/`Pr_IBD2` pair is appended on the same ≥ 100 Mb trigger as `.ibs`.

```bash
king -b /tmp/kingdocs/multifam.bed --ibs
```

```
FID1	ID1	FID2	ID2	N_SNP	N_IBS0	N_IBS1	N_IBS2	NHetHet	NHomHom	N_Het1	N_Het2	IBS	Dist	HetConc	Het2|1	Het1|2	HomConc	Kinship	MaxIBD2	Pr_IBD2
FAM1	A_F	FAM2	B_F	15000	306	4031	10663	3253	7716	5257	5280	1.6905	0.3503	0.4466	0.6188	0.6161	0.9603	0.2501	35150873.000	0.2494
FAM1	A_F	FAM2	B_M	15000	1028	6400	7572	2079	6521	5257	5301	1.4363	0.7008	0.2452	0.3955	0.3922	0.8424	0.0001	-9	-9
FAM1	A_F	FAM2	B_C1	15000	451	5937	8612	2304	6759	5257	5288	1.5441	0.5161	0.2796	0.4383	0.4357	0.9333	0.1319	0.000	0.0000
```

**`MaxIBD2` and `Pr_IBD2` are the literal string `-9` for any pair whose `Kinship` is below
`2^-3.5 = 0.08839`** — the pair was never handed to the segment scanner. Verified over
every row of a 200-sample run:

```
.ibs0 -9 rule: 19327 rows, violations of 'MaxIBD2/Pr_IBD2 == -9 iff Kinship < 2^-3.5': 0
   max Kinship on a -9 row: 0.0882   min Kinship on a valued row: 0.1172
```

Row 2 above is that case: kinship 0.0001, so `-9  -9`. Row 3 is above the gate at 0.1319
and was analysed — it simply has no IBD2, hence `0.000  0.0000`. **Do not read `-9` as
zero**: it means "not measured".

## `<prefix>.seg`

Per-pair IBD segment summary from `--ibdseg`. Eight tab-separated columns; both within-
and between-family pairs in one file.

**Written when** `--ibdseg` runs and reaches the segment stage. Rows are **filtered**: a
pair appears only if it has at least one long (> 10 Mb) IBD segment, and only segments
above the `--seglength` floor (default 3 Mb) are counted. Row order is 16-sample blocks.

```bash
king -b /tmp/kingdocs/multifam.bed --ibdseg
```

```
FID1	ID1	FID2	ID2	IBD1Seg	IBD2Seg	PropIBD	InfType
FAM1	A_F	FAM1	A_C1	1.0000	0.0000	0.5000	PO
FAM1	A_F	FAM1	A_C2	1.0000	0.0000	0.5000	PO
FAM1	A_F	FAM1	A_C3	1.0000	0.0000	0.5000	PO
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1–4 | `FID1`, `ID1`, `FID2`, `ID2` | the pair — `FID1 == FID2` for a within-family pair | text |
| 5 | `IBD1Seg` | proportion of the usable genome (the `allsegs.txt` total) covered by IBD1 segments | `%.4f` |
| 6 | `IBD2Seg` | proportion covered by IBD2 segments | `%.4f` |
| 7 | `PropIBD` | `IBD2Seg + IBD1Seg/2`, **computed from the two printed columns above** | `%.4f` |
| 8 | `InfType` | inferred relationship | text |

**How to read this.** `PropIBD` is the single number to sort on: 0.5 for a first-degree
relationship, 0.25 for second, 0.125 for third. The `IBD1Seg`/`IBD2Seg` split is what
separates PO (`1.0000 / 0.0000`) from FS (roughly `0.50 / 0.25`) at the same `PropIBD`.

> **`.seg`'s `PropIBD` is not `.kin`'s.** `.seg` recombines its own four-decimal printed
> columns (`i2·1e-4 + i1·5e-5`); every other file computes the same quantity at full
> precision from the underlying totals. The reference disagrees with itself here and
> open-king reproduces both rules. Measured on a run that writes both files:
>
> ```bash
> king -b /tmp/kingdocs/bigish.bed --related --degree 2 --ibdseg --cpus 4 --prefix r
> ```
>
> ```
> r.kin and r.seg share 147 pairs; PropIBD differs on 43
>     B02_C1 B02_C2 IBD1Seg 0.3852 IBD2Seg 0.3123 -> .kin 0.5048 / .seg 0.5049
>     B02_C3 B02_C4 IBD1Seg 0.4885 IBD2Seg 0.2974 -> .kin 0.5417 / .seg 0.5416
> ```
>
> Both directions of disagreement occur, and it is always the last digit. If you need one
> number, take `.kin`'s.

`--seglength` moves the estimates by changing which segments are reported at all:

```bash
king -b /tmp/kingdocs/multifam.bed --ibdseg                --prefix s3
king -b /tmp/kingdocs/multifam.bed --ibdseg --seglength 5  --prefix s5
king -b /tmp/kingdocs/multifam.bed --ibdseg --seglength 10 --prefix s10
for p in s3 s5 s10; do printf "%-4s " $p; awk -F'\t' '$2=="A_C1" && $4=="A_C2"' $p.seg; done
```

```
s3   FAM1	A_C1	FAM1	A_C2	0.4516	0.3144	0.5402	FS
s5   FAM1	A_C1	FAM1	A_C2	0.4516	0.3144	0.5402	FS
s10  FAM1	A_C1	FAM1	A_C2	0.4379	0.3144	0.5333	FS
```

All three files have the same 104 rows — the floor changes the estimates, not the pair set,
on this dataset.

## `<prefix>allsegs.txt`

The map's **usable segments** — the denominator every IBD proportion is computed against.
Eight tab-separated columns, one row per segment, in map order (autosomal segments first,
then X).

**Written by every analysis that runs the segment pre-pass**: `--ibs`, `--ibdseg`,
`--related`, `--unrelated`, `--cluster`, `--build`, `--bysample`, `--bySNP`. Not written by
`--kinship`, `--duplicate` or `--autoQC`. Not written when the map yields no usable segment
at all.

A segment is a maximal run of markers cut at every chromosome change and every base-pair
gap over 1 Mb, cut again between 64-marker words spanning more than 10 Mb, and then kept
only if it holds at least five complete 64-marker words and spans more than 10 Mb. That
last filter is why a 22-chromosome map can yield 18 segments rather than 22.

```bash
king -b /tmp/kingdocs/multifam.bed --ibdseg
```

```
Segment	Chr	StartMB	StopMB	Length	N_SNP	StartSNP	StopSNP
1	1	1.010	65.904	64.895	1299	rs1_1009689	rs1_65904473
2	2	1.006	64.158	63.152	1264	rs2_1005960	rs2_64158411
3	3	1.010	52.710	51.700	1035	rs3_1009935	rs3_52709533
4	4	1.002	50.604	49.603	993	rs4_1001547	rs4_50604488
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1 | `Segment` | 1-based index, running across the whole file — autosomal segments then X | `%d` |
| 2 | `Chr` | chromosome, numeric (X appears as `23`, or as whatever `--sexchr` was set to) | `%d` |
| 3 | `StartMB` | first marker's position, **in megabases** | `%.3f` |
| 4 | `StopMB` | last marker's position, in megabases | `%.3f` |
| 5 | `Length` | `StopMB − StartMB`, in megabases | `%.3f` |
| 6 | `N_SNP` | markers in the segment | `%d` |
| 7 | `StartSNP` | first marker's ID | text |
| 8 | `StopSNP` | last marker's ID | text |

Summing `Length` over the **autosomal** rows gives the figure the console prints
(`Total length of 18 chromosomal segments usable for IBD segment analysis is 691.5 Mb.`),
which must reach **100 Mb** for `.ibs`/`.ibs0` to gain their two segment columns. X
segments are listed in the same file but totalled separately, on the console's
`In addition to autosomes, <n> segments of length <x> Mb on X-chr…` line — and they are the
denominator for the X files, never for the autosomal ones. A map with X markers therefore
looks like this:

```bash
king -b /tmp/kingdocs/sexchr.bed --ibdseg --degree 2
```

```
Segment	Chr	StartMB	StopMB	Length	N_SNP	StartSNP	StopSNP
1	1	1.007	100.953	99.946	2000	rs1_1007474	rs1_100953159
2	2	1.004	100.952	99.947	2000	rs2_1004059	rs2_100951520
3	23	1.002	75.954	74.952	1500	rs23_1002040	rs23_75953987
```

```
Total length of 2 chromosomal segments usable for IBD segment analysis is 199.9 Mb.
  In addition to autosomes, 1 segments of length 75.0 Mb on X-chr can be further used.
  Information of these chromosomal segments can be found in file kingallsegs.txt
```

## `<prefix>splitped.txt`

A pedigree file `--ibdseg` leaves behind for downstream pedigree-plotting tools. **Nine
space-separated fields, no header.** Nothing in it depends on genotypes, `--degree` or
`--seglength`.

**Written by `--ibdseg`**, before any segment work, when at least one family has two members
or a singleton names a parent — so it survives an early segment exit once that pedigree gate
is met. All-parentless singleton families write and announce nothing.

```
OldFID OldIID NewFID NewIID Father Mother Sex Pheno Dummy
```

| # | field | meaning |
| ---: | --- | --- |
| 1 | `OldFID` | FID as it appears in your `.fam` |
| 2 | `OldIID` | IID as it appears in your `.fam` |
| 3 | `NewFID` | FID after splitting — `<FID>_S<k>` when one declared family turns out to be several disconnected pedigrees, otherwise unchanged |
| 4 | `NewIID` | the IID; never changes |
| 5 | `Father` | father's IID, `0` if none |
| 6 | `Mother` | mother's IID, `0` if none |
| 7 | `Sex` | `1` male, `2` female, `0` unknown |
| 8 | `Pheno` | the `.fam` phenotype, with `-9` rewritten to `0` |
| 9 | `Dummy` | `1` if this person was invented by the program, `0` if genotyped |

Three things it does that your `.fam` does not: it materialises absent parents, invents a
mate for a half-specified parentage, and splits a family that turns out to be several
disconnected pedigrees.

**Materialises absent parents.** `multifam`'s `FAM3` declares `C_F`'s parents as `A_F`/`A_M`
— who live in `FAM1`. They appear inside `FAM3` as dummy founders:

```
FAM3 A_F FAM3 A_F 0 0 1 0 1
FAM3 A_M FAM3 A_M 0 0 2 0 1
FAM3 C_M FAM3 C_M 0 0 2 0 0
FAM3 C_F FAM3 C_F A_F A_M 1 0 0
FAM3 C_C1 FAM3 C_C1 C_F C_M 1 0 0
FAM3 C_C2 FAM3 C_C2 C_F C_M 2 0 0
FAM3 C_C3 FAM3 C_C3 C_F C_M 1 0 0
```

**Invents a mate, and splits a disconnected family.** Both in one file:

```bash
king -b /tmp/kingdocs/dups.bed --ibdseg
```

```
MZFAM MZ_1 MZFAM_S1 MZ_1 0 0 2 0 0
MZFAM MZ_2 MZFAM_S2 MZ_2 0 0 2 0 0
POFAM KING1 POFAM KING1 0 0 2 0 1
POFAM PO_P POFAM PO_P 0 0 1 0 0
POFAM PO_C POFAM PO_C PO_P KING1 2 0 0
```

`MZFAM` holds two people the `.fam` declares no relationship between, so it splits into
`MZFAM_S1` and `MZFAM_S2`. `PO_C` names a father but no mother, so a founder `KING1` is
invented (`Dummy = 1`) and `PO_C`'s row is rewritten to point at her. `KING<n>` numbering
is global, in `.fam` family order.

Rows are ordered by FID under the ID comparator, and within a family by generation depth
then IID — founders first. A family of a single parentless founder is dropped entirely.

---

# Sample-selection and pedigree files

## `<prefix>unrelated.txt` and `<prefix>unrelated_toberemoved.txt`

The output of `--unrelated`: a maximal set of mutually unrelated individuals, and its
complement. **Two tab-separated columns, `FID` and `IID`, no header.** Together they
partition your sample set exactly.

**Written when:** always, by every `--unrelated` run — including on datasets under 10
samples, where the clustering step is disabled and selection falls back to the pedigree.

```bash
king -b /tmp/kingdocs/multifam.bed --unrelated
```

`kingunrelated.txt` has 8 rows — every family's two founders, the children dropped. Its
first four:

```
FAM1	A_F
FAM1	A_M
FAM2	B_F
FAM2	B_M
```

`kingunrelated_toberemoved.txt` has the complementary 12. Its first three:

```
FAM1	A_C1
FAM1	A_C2
FAM1	A_C3
```

The tiny-dataset path still writes both:

```bash
king -b /tmp/kingdocs/trio.bed --unrelated
```

console: `This function is currently disabled for tiny dataset with sample size < 10.`

`kingunrelated.txt`:

```
TRIO	T_F
TRIO	T_M
```

`kingunrelated_toberemoved.txt`:

```
TRIO	T_C1
```

Two things to know before you use these lists:

* **A pair is an edge when either the pedigree or the genotypes make it closer than
  unrelated** — kinship above `2^-5.5`, the 4th-degree band edge, not the 1st-degree one
  the console message mentions. `--degree` does **not** move it.
* **Only pairs inside one cluster are considered.** Cross-family relatives are invisible to
  the selection unless clustering merged their families first, which needs ≥ 100 samples.
  In `dups`, the exact duplicate pair spans two FIDs and **both copies survive**.

The kept list is written in the order the greedy selection picked, which is neither `.fam`
order nor sorted order.

## `<prefix>cluster.kin`

Pairwise relatedness inside families that `--cluster` newly merged. Fifteen tab-separated
columns — the 14-column `--related` `.kin0` set, reorganised: one `FID` (the new cluster
ID), both IIDs, and a `Sex1`/`Sex2` pair no other file carries.

**Written when** `--cluster` actually merges at least two families *and* the map yields an
informative IBD segment. A run that merges nothing writes only `allsegs.txt`. On the
kinship-only sparse fallback, `updateids.txt` still records the merge but `cluster.kin` and
`allsegs.txt` are absent. Merging requires ≥ 100 samples.

```bash
king -b /tmp/kingdocs/bigish.bed --cluster --cpus 4
```

```
FID	ID1	ID2	Sex1	Sex2	N_SNP	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
KING1	B01_C1	B01_C2	1	2	50000	0.2061	0.0158	0.4202	0.1357	0.2505	0.5328	0.2569	0.5233	FS
KING1	B01_C1	B01_C3	1	1	50000	0.2056	0.0149	0.4208	0.1288	0.2533	0.5157	0.2587	0.5166	FS
KING1	B01_C1	B01_F	1	1	50000	0.1724	0.0000	0.3296	0.0000	0.2479	1.0000	0.0000	0.5000	PO
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1 | `FID` | the **new** cluster ID, `KING1`, `KING2`, … | text |
| 2–3 | `ID1`, `ID2` | the pair, both now inside the cluster | text |
| 4–5 | `Sex1`, `Sex2` | `.fam` sex codes | `%d` |
| 6–15 | `N_SNP` … `InfType` | as in the 14-column `.kin0` above; `Kinship` is the **within**-family estimator, since the pair is now in one family | — |

Rows are every pair of every merged cluster: `KING1` above has 11 members and so 55 rows,
and the file has 165 rows for three clusters. Members are ordered by the ID comparator,
pairs are the upper triangle.

Read this together with `updateids.txt`, which tells you which original families each
`KING<k>` absorbed.

## `<prefix>updateids.txt`

The FID rename map produced by clustering. **Four tab-separated columns, no header:**

```
OldFID   OldIID   NewFID   NewIID
```

The IID never changes, so column 4 always repeats column 2 — the file is really a FID
rename table, in the shape PLINK's `--update-ids` wants.

**Written by `--cluster` and `--build`, only when at least one merge happened.** Note that
`--build` prints `Update-ID information is saved in file <p>updateids.txt` **whether or not
it wrote the file** — on an unmerged run the line appears and no file exists. (That is the
reference's behaviour, reproduced.)

```bash
king -b /tmp/kingdocs/bigish.bed --build --cpus 4
```

```
BF01	B01_C1	KING1	B01_C1
BF01	B01_C2	KING1	B01_C2
BF01	B01_C3	KING1	B01_C3
BF01	B01_F	KING1	B01_F
```

Only merged clusters get rows — 33 rows here for three clusters. Rows are ordered by the
**original** `(FID, IID)` under the ID comparator, *not* by cluster.

## `<prefix>updateparents.txt`

The reconstructed parentage produced by `--build`. **Four tab-separated columns, no
header:**

```
FID   IID   Father   Mother
```

in the shape PLINK's `--update-parents` wants.

**Written by every `--build` run that gets past the tiny-dataset gate**, but **empty
(zero bytes) when nothing was reconstructed** — which is also when the console says
`No pedigrees can be reconstructed.`

```bash
king -b /tmp/kingdocs/bigish.bed --build --cpus 4
```

```
KING1	B01_C1	B01_F	B01_M
KING1	B01_C2	B01_F	B01_M
KING1	B01_C3	B01_F	B01_M
KING1	B01_F	1	2
KING1	B01_M	0	0
KING1	B02_C1	B02_F	B02_M
```

| # | column | meaning |
| ---: | --- | --- |
| 1 | `FID` | the cluster ID after merging (`KING1`, …) |
| 2 | `IID` | the person |
| 3 | `Father` | father's IID, `0` if none; a bare integer names an **inferred** phantom parent |
| 4 | `Mother` | likewise |

**How to read row 4.** `B01_F` was a founder in your `.fam`, but reconstruction found he
and `B02_F` are full sibs, so it invented a parent couple for them, numbered `1` and `2`.
Those numbers are the same ones that appear in `build.log`'s `RULE FS0` line. Row 5 shows
`B01_M` still has no parents.

Rows are in **cluster order** (`KING1`'s rows, then `KING2`'s) — a different order from
`updateids.txt`, which is in original-`(FID, IID)` order. The two coincide only when your
FIDs happen to sort the same way.

## `<prefix>build.log`

A free-text narrative of what `--build` reconstructed, echoed to stdout as it is written.

**Written by every `--build` run that gets past the tiny-dataset gate**, and **zero bytes
when there was nothing to reconstruct.**

```bash
king -b /tmp/kingdocs/bigish.bed --build --cpus 4
```

```
Family KING1:
  Family KING1 RULE FS0: Sibship (B01_F B02_F)'s parents are (1 2)
Family KING2:
  Family KING2 RULE FS0: Sibship (B13_F B14_F)'s parents are (3 4)
Family KING3:
  Family KING3 RULE FS0: Sibship (B25_F B26_F)'s parents are (5 6)
```

Two line kinds:

* `Family <FID>:` — a header opening one cluster's block.
* `  Family <FID> RULE <name>: <what was inferred>` — one reconstruction step. `FS0` above
  means "these people are full sibs with no genotyped parents, so give them an invented
  parent couple"; the numbers `(1 2)` are the phantom parent IDs that appear in
  `updateparents.txt`.

```bash
king -b /tmp/kingdocs/multifam.bed --build
```

leaves `kingbuild.log` at **0 bytes** — nothing was reconstructed.

> This file is byte-identical in all captured cases, including the complete `INFERENCE`
> half on `bigish`. See [Known divergences](#known-divergences-from-king-232) for rare
> held-out pedigree shapes.

---

# QC files

## `<prefix>bySample.txt`

Per-sample QC from `--bysample`. **Space separated**, one row per sample in `.fam` order.

**The column list grows and shrinks with the data.** The first eight columns are always
present; X, Y and MT blocks appear only if the map has those markers, and the Mendelian
blocks only if the pedigree has parent–offspring pairs / trios.

Minimal form — a fileset of unrelated singletons with autosomes only:

```bash
king -b /tmp/kingdocs/unrelated.bed --bysample --prefix un
```

```
FID IID FA MO SEX N_SNP Missing Heterozygosity
POOL P01 0 0 1 20000 0.0000 0.3475
POOL P02 0 0 2 20000 0.0000 0.3476
POOL P03 0 0 1 20000 0.0000 0.3486
```

Full form — every optional block present:

```bash
king -b /tmp/kingdocs/sexchr.bed --bysample
```

```
FID IID FA MO SEX N_SNP Missing Heterozygosity N_xSNP xHeterozygosity N_ySNP N_yHetero N_mtSNP N_mtHetero N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal
SEX S_F 0 0 1 4150 0.0000 0.3489 1500 0.0000 300 0 50 0 16600 0 0.0000 16600 0 0.0000 0
SEX S_M 0 0 2 4150 0.0000 0.3465 1500 0.3393 0 0 50 0 16600 0 0.0000 16600 0 0.0000 0
```

| column | present when | meaning | format |
| --- | --- | --- | --- |
| `FID` `IID` | always | the sample | text |
| `FA` `MO` | always | parent IIDs from the `.fam`; a sample naming exactly one parent gets an invented `KING<k>` in the empty slot | text |
| `SEX` | always | `.fam` sex code | `%d` |
| `N_SNP` | always | autosomal (and `XY`) markers called for this sample | `%d` |
| `Missing` | always | autosomal missing rate | `%.4f` |
| `Heterozygosity` | always | autosomal het rate over called markers | `%.4f` |
| `N_xSNP` `xHeterozygosity` | map has X markers | X call count and X het **rate** | `%d`, `%.4f` |
| `N_ySNP` `N_yHetero` | map has Y markers | Y call count and Y het **count** — a count, not a rate; the reference's own asymmetry with the X block | `%d`, `%d` |
| `N_mtSNP` `N_mtHetero` | map has MT markers | MT call count and MT het count | `%d`, `%d` |
| `N_pair` `N_MIp` `Err_MIp` | pedigree has ≥ 1 PO pair | markers examined across this sample's PO pairs, Mendelian inconsistencies among them, and the rate | `%d`, `%d`, `%.4f` |
| `N_trio` `N_MIt` `Err_MIt` | pedigree has ≥ 1 trio | the same for trios | `%d`, `%d`, `%.4f` |
| `MI_Removal` | pedigree has ≥ 1 PO pair | `1` if this sample's Mendelian error rate suggests removing it, else `0` | `%d` |

**How to read `S_M`'s row.** (In this fixture `S_F` is the father, `SEX = 1`, and `S_M` the
mother, `SEX = 2` — the names are the roles, not the sexes.) The mother has 4 150 autosomal
calls, no missingness, 34.65 % heterozygous. Her X heterozygosity is 0.3393, consistent
with two X chromosomes; the father's is `0.0000`, being hemizygous. She has 0 Y calls; he
has 300. `N_pair 16600` is 4 PO pairs × 4 150 markers, with 0 inconsistencies, so nothing
is flagged.

**What counts as a Mendelian inconsistency here:** a PO *pair* is inconsistent when the two
are opposite homozygotes; a *trio* is inconsistent when the offspring is heterozygous and
both parents are homozygous for the same allele — and in no other case. The two checks are
separate, not nested.

`bySample.txt` counts autosomal markers only (`XY`/chromosome 25 pooled in, as everywhere);
X, Y and MT contribute only their own blocks.

## `<prefix>bySNP.txt`

Per-marker QC from `--bySNP`. **Space separated**, one row per retained marker.

**Row order:** autosomes and `XY` first in `.bim` order, then X, then Y, then MT. So
`bySNP.txt` is *not* simply the `.bim` in file order when your map interleaves classes.

```bash
king -b /tmp/kingdocs/sexchr.bed --bySNP
head -1 kingbySNP.txt; awk 'NR>1{if(!seen[$2]++) print}' kingbySNP.txt
```

— the header, and the first real row of each chromosome class (which is also a listing of
the class order):

```
SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate N_PO N_HomPO N_errPO Err_InPO Err_InHomPO N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio
rs1_1007474 1 1007474 T A 0.2000 10 0 4 6 1.0000 8 2 0 0.0000 0.0000 4 2 0 0.0000 0.0000
rs2_1004059 2 1004059 G A 0.1500 10 1 1 8 1.0000 8 8 0 0.0000 0.0000 4 0 0 0.0000 0.0000
rs25_279120 25 279120 C T 0.0000 10 0 0 10 1.0000 8 8 0 0.0000 0.0000 4 0 0 0.0000 0.0000
rs23_1002040 X 1002040 G A 0.5000 10 4 2 4 1.0000 8 4 2 0.2500 0.5000 4 2 0 0.0000 0.0000
rs24_1006052 Y 1006052 T A 0.0000 4 0 0 4 0.4000 2 2 0 0.0000 0.0000 0 0 0 0 0
rs26_1687 MT 1687 T C 0.0000 10 0 0 10 1.0000 0 0 0 0 0 0 0 0 0 0
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1 | `SNP` | marker ID from the `.bim` | text |
| 2 | `Chr` | chromosome — **`X`, `Y` and `MT` are spelled out, everything else is numeric** (note `25`, the pseudoautosomal code, stays a number) | text |
| 3 | `Pos` | base-pair position | `%d` |
| 4 | `Label_A` | the A1 allele | text |
| 5 | `Label_a` | the A2 allele | text |
| 6 | `Freq_A` | A1 frequency among called genotypes | `%.4f` |
| 7 | `N` | samples called at this marker | `%d` |
| 8 | `N_AA` | homozygous A1 | `%d` |
| 9 | `N_Aa` | heterozygous | `%d` |
| 10 | `N_aa` | homozygous A2 | `%d` |
| 11 | `CallRate` | `N` over the total sample count | `%.4f` |
| 12 | `N_PO` | parent–offspring pairs where both are called here | `%d` |
| 13 | `N_HomPO` | of those, pairs where both are homozygous | `%d` |
| 14 | `N_errPO` | of those, Mendelian-inconsistent pairs | `%d` |
| 15 | `Err_InPO` | `N_errPO / N_PO` | `%.4f` |
| 16 | `Err_InHomPO` | `N_errPO / N_HomPO` | `%.4f` |
| 17 | `N_trio` | complete trios called here | `%d` |
| 18 | `N_HetOff` | of those, trios whose offspring is heterozygous | `%d` |
| 19 | `N_errTrio` | Mendelian-inconsistent trios | `%d` |
| 20 | `Err_InTrio` | `N_errTrio / N_trio` | `%.4f` |
| 21 | `Err_InHetTrio` | `N_errTrio / N_HetOff` | `%.4f` |

**How to read the X row.** `rs23_1002040` is called in all 10 samples, A1 frequency 0.5. Of
the 8 PO pairs, 4 are hom–hom and **2 are Mendelian-inconsistent** — `Err_InPO 0.2500`,
`Err_InHomPO 0.5000`. On a real dataset a marker like that is a genotyping failure.

**Two per-class quirks, both real and both visible above.** Y markers get no trio statistics
and MT markers get neither PO nor trio statistics; in those cells the file prints a bare
integer `0` rather than the `0.0000` an autosomal zero-denominator row would print. So
column 20 reads `0.0000` on the autosomal rows and `0` on the Y and MT rows. A strict
fixed-format parser will notice.

## The four `<prefix>_autoQC_*` files

`--autoQC` is the one pass that filters anything. It writes four reports; none of the other
analyses does.

**`<prefix>_autoQC_Summary.txt`** — the step table, byte-identical to the block `--autoQC`
prints on stdout. **Fixed-width, space padded; it contains no tab at all.**

```bash
king -b /tmp/kingdocs/sexchr.bed --autoQC
```

```
Step Description                                            Subjects  SNPs      
1    Raw data counts                                        10        6000      
1.1  SNPs with very low call rate < 80% (removed)                     (0)
1.2  Monomorphic SNPs (removed)                                       (801)
1.3  Sample call rate < 95% (removed)                       (0)
1.4  SNPs with call rate < 95% (removed)                              (0)
2    data counts for gender error checking                  10        5199      
2.1  Y-chr SNPs with call rate < 80% in men (removed)                 (0)
2.2  X-chr SNPs with heterozygosity > 5% in men (removed)             (0)
2.3  Y-chr SNPs with genotypes in >10% women (removed)                (0)
2.4  Mislabeled as male (removed)                           (0)
2.5  Mislabeled as female (removed)                         (0)
2.6  Suspicious gender error (removed)                      (0)
2.7  Y-chr SNPs with call rate < 95% in men (removed)                 (0)
2.8  X-chr SNPs with heterozygosity > 1% in men (removed)             (0)
2.9  Y-chr SNPs with genotypes in >2% women (removed)                 (0)
3    Generate Final Study Files                             
     Final QC'ed data                                       10        5199      
```

Layout: `%-5s%-55s%-10s%-10s` for the two count rows, `%-5s%-65s(%d)` for a SNP counter and
`%-5s%-55s(%d)` for a subject counter. A parenthesised number is a removal count; the two
unbracketed pairs are the surviving `(subjects, SNPs)` before and after gender QC. The
`2.x` block appears only when the map has both X and Y markers.

**Caution:** the `1.1` / `1.3` / `1.4` labels are *fixed text*. They say `< 80%` and `< 95%`
whatever you passed for `--callrateN` / `--callrateM` (both default 0.95); the console lines
above them print the thresholds actually applied.

**`<prefix>_autoQC_snptoberemoved.txt`** — tab separated, header `SNP  REASON`, one row per
removed marker, **grouped by the step that removed it** (autosomes, then X, then Y within a
group; `.bim` order within a class).

```
SNP	REASON
rs1_2353935	Monomorphic
rs1_2506794	Monomorphic
```

Reason strings: `CallRateLessThan<pct>`, `Monomorphic`, `xHeterozygosityInMale`,
`YSNPInFemales`. `<pct>` is `round(100 × callrateM)` for *both* Y call-rate filters, even
the one that actually applies the looser step-1 threshold. On a dataset with real
missingness the file shows the grouping:

```bash
king -b /tmp/kingdocs/missing.bed --autoQC
awk 'NR>1{print $2}' king_autoQC_snptoberemoved.txt | uniq -c
```

```
1569 CallRateLessThan80
3199 Monomorphic
 130 CallRateLessThan95
```

`uniq` without `sort`, so those are *consecutive runs* in file order: step 1's call-rate
filter, then step 1's monomorphic filter, then step 3's call-rate filter. The reasons do
not interleave.

**`<prefix>_autoQC_sampletoberemoved.txt`** — tab separated, header `FID  IID  REASON`, one
row per removed sample, grouped by check, `.fam` order within a group. Header-only when
nothing is removed.

```
FID	IID	REASON
MIS	M_C2	MissingMoreThan5
MIS	M_C3	MissingMoreThan5
```

Reason strings: `MissingMoreThan<pct>` where `<pct>` is `round(100 × (1 − callrateN))`, then
`MislabeledAsMale`, `MislabeledAsFemale`, `GenderQC`.

**`<prefix>_autoQC_updatesex.txt`** — tab separated, **no header**, `FID  IID  sex` for every
sample whose `.fam` sex was `0` and whose sex the pass inferred, in `.fam` order.

```
SU1	S_U0A	2
SU2	S_U0B	2
```

**This file is not created when there is no such sample.** The `missing` run above writes
the other three and not this one.

---

# X-chromosome files

Three different files carry X results and they are **not** variants of one layout — the
same name means different columns under different flags.

## `<prefix>X.kin`

**Under `--kinship`: within-family X kinship, 9 columns.**

The X pass runs only when the map has **≥ 512 X markers**, **no `--degree` was given**, and
the autosomal between-family stage ran (so: at least two families). Both X files are then
written unconditionally.

```bash
king -b /tmp/kingdocs/sexchr.bed --kinship
```

```
FID	ID1	ID2	Sex	N_SNP	PhiX	Het	IBS0	KinshipX
SEX	S_DAU1	S_DAU2	FF	1500	0.3750	0.327	0.0000	0.3262
SEX	S_DAU1	S_F	FM	1500	0.5000	0.331	0.0000	0.5000
SEX	S_DAU1	S_M	FF	1500	0.2500	0.335	0.0000	0.2435
```

| # | column | meaning | format |
| ---: | --- | --- | --- |
| 1–3 | `FID`, `ID1`, `ID2` | the within-family pair | text |
| 4 | `Sex` | the pair's sexes as a two-letter code **in the row's own order** — `FM` is a female `ID1` and a male `ID2` | text |
| 5 | `N_SNP` | X markers called in both | `%d` |
| 6 | `PhiX` | pedigree X-kinship coefficient | `%.4f` |
| 7 | `Het` | the pair's X heterozygosity denominator (see below) | `%.3f` — three decimals, not four |
| 8 | `IBS0` | proportion of `N_SNP` at which they are opposite homozygotes | `%.4f` |
| 9 | `KinshipX` | estimated X kinship | `%.4f` |

Three estimators, chosen by the pair's sexes, since a male carries one X:

```
FF   Het = (Het_i + Het_j) / 2N     KinshipX = (HetHet - 2*IBS0) / (Het_i + Het_j)
FM   Het = Het_female / N           KinshipX = 0.5  - IBS0 / Het_female
MM   Het = H (imputed)              KinshipX = 0.75 - IBS0 / (N * H)
```

`H` for a hemizygous male is imputed as the lower median of the female X heterozygosity
rates in his own family. Unlike the autosomal `.kin0`, the between-family X file uses the
**same three forms** — there is no population-structure-robust variant on X.

**Samples whose `.fam` sex is neither 1 nor 2 are excluded outright** from the X pass, and
appear in no X row.

**How to read row 2.** A daughter and her father: `PhiX 0.5000`, because a father passes
his single X entire to every daughter, and the estimate matches at exactly `0.5000` with
zero opposite homozygotes. Row 1's two sisters get `PhiX 0.3750` — the standard X value for
full sisters, not the autosomal 0.25.

**Under `--related`: within-family X *segment* sharing, 9 different columns.**

```bash
king -b /tmp/kingdocs/sexchr.bed --related
```

console: `Within-family X-chr IBD-sharing inference saved in file kingX.kin`

```
FID	ID1	ID2	Sex1	Sex2	PhiX	IBD1Seg	IBD2Seg	PropIBD
SEX	S_DAU1	S_DAU2	2	2	0.3750	0.6464	0.3530	0.6762
SEX	S_DAU1	S_F	2	1	0.5000	1.0000	0.0000	0.5000
SEX	S_DAU1	S_M	2	2	0.2500	1.0000	0.0000	0.5000
```

Same file name, different table: numeric `Sex1`/`Sex2` instead of the two-letter `Sex`, no
`N_SNP`, no `IBS0`, no `KinshipX`, and the three segment columns instead. `PropIBD` here is
the **full-precision** `IBD2Seg + IBD1Seg/2`. The `--related` X pass has no 512-marker
threshold — its gate is that the X map yields a usable segment — and it is not suppressed
by `--degree`.

## `<prefix>X.kin0`

**Under `--kinship`: between-family X kinship, 9 columns.** Same statistics as `--kinship`'s
`X.kin`, minus `PhiX`, plus the second FID. Written whenever the `--kinship` X pass runs.

```bash
king -b /tmp/kingdocs/sexchr.bed --kinship
```

```
FID1	ID1	FID2	ID2	Sex	N_SNP	Het	IBS0	KinshipX
SEX	S_DAU1	SU3	S_UM	FM	1500	0.331	0.1707	-0.0151
SEX	S_DAU1	SU4	S_UF	FF	1500	0.348	0.0633	0.0105
SEX	S_DAU2	SU3	S_UM	FM	1500	0.323	0.1727	-0.0351
```

Row order here is **family-major** — families in sorted order, members sorted inside them,
plain `i < j` upper triangle — and *not* the autosomal `.kin0`'s 32-tiled `.fam` order.

**How to read row 1.** A daughter from the `SEX` family against an unrelated male: X
kinship −0.0151, i.e. zero within noise. Negative values are normal and mean "no more
sharing than the reference population".

**Under `--related`: between-family X segment sharing, 9 different columns.** Written when
the `--related` between-family stage writes a `.kin0` and the X map has a usable segment —
so, in practice, ≥ 100 samples with detected relatives. Rows are exactly `.kin0`'s pairs,
re-measured on X.

No corpus dataset has both ≥ 100 samples and X markers, so this one is built by importing
the corpus generator as a module — 104 samples in 26 nuclear families, alternate pairs of
which share undeclared grandparents, over three autosomes plus 2 000 X markers:

```python
# bigx.py
import sys
sys.path.insert(0, "/path/to/open-king/tests/parity")
import generate_corpus as g

ped, fam, prev = g.Ped(), 0, None
while ped.n_emitted() < 104:
    fam += 1
    fp = prev if fam % 2 == 0 else g.add_couple(ped, "XPH%02d" % fam, "XPH%02d" % fam, emit=False)
    if fam % 2 == 1:
        prev = fp
    g.add_nuclear(ped, "XF%02d" % fam, "X%02d" % fam, 2, father_parents=fp)

spec = g.Spec("bigx", ped, {1: 4000, 2: 4000, 3: 4000, 23: 2000, 24: 300}, 14300)
g.simulate(spec, g.dataset_seed(20260813, "bigx"), ".")
```

```bash
python3 bigx.py
king -b bigx.bed --related --cpus 4
```

```
FID1	ID1	FID2	ID2	Sex1	Sex2	IBD1Seg	IBD2Seg	PropIBD
XF01	X01_F	XF02	X02_F	1	1	0.0315	0.1916	0.2073
XF03	X03_F	XF04	X04_F	1	1	0.0000	0.1996	0.1996
XF03	X03_F	XF04	X04_C1	1	1	0.0000	0.0000	0.0000
XF03	X03_F	XF04	X04_C2	1	2	0.2221	0.0000	0.1111
```

The columns are `.seg`'s three estimates plus the two sexes; `PropIBD` is full precision.
Rows are exactly `.kin0`'s pairs, in `.kin0`'s order.

## `<prefix>X.seg`

The X half of `--ibdseg` — the same pairs `.seg` reports, measured again over the X array.

**Written when `--degree` is given with a non-zero value *and* the X map yields a usable
segment.** Plain `--ibdseg` on a fileset with 1 500 X markers writes **no** `X.seg`;
`--degree 1`, `2`, `3` … all write it. This is the reference's own oddity, reproduced.
There is no marker-count threshold — the gate is the same usable-segment construction
`allsegs.txt` reports, which is why the file appears exactly when the console prints
`In addition to autosomes, <n> segments of length <x> Mb on X-chr can be further used.`

```bash
king -b /tmp/kingdocs/sexchr.bed --ibdseg --degree 2
cat -e kingX.seg
```

`cat -e` marks the line ends, so the trailing tab on every data row is visible:

```
FID1	ID1	FID2	ID2	Sex1	Sex2	MaxIBD1	MaxIBD2	IBD1Seg	IBD2Seg	PropIBD$
SEX	S_F	SEX	S_SON1	1	1	0.0000	0.0000	0.0000	$
SEX	S_F	SEX	S_SON2	1	1	0.0000	0.0000	0.0000	$
SEX	S_F	SEX	S_DAU1	1	2	1.0000	0.0000	0.5000	$
SEX	S_F	SEX	S_DAU2	1	2	1.0000	0.0000	0.5000	$
SEX	S_M	SEX	S_SON1	2	1	1.0000	0.0000	0.5000	$
SEX	S_M	SEX	S_SON2	2	1	1.0000	0.0000	0.5000	$
SEX	S_M	SEX	S_DAU1	2	2	1.0000	0.0000	0.5000	$
SEX	S_M	SEX	S_DAU2	2	2	1.0000	0.0000	0.5000	$
SEX	S_SON1	SEX	S_SON2	1	1	0.1462	0.6393	0.7124	$
SEX	S_SON1	SEX	S_DAU1	1	2	0.4257	0.0000	0.2128	$
SEX	S_SON1	SEX	S_DAU2	1	2	0.9067	0.0000	0.4533	$
SEX	S_SON2	SEX	S_DAU1	1	2	0.6397	0.0000	0.3199	$
SEX	S_SON2	SEX	S_DAU2	1	2	0.7245	0.0000	0.3623	$
SEX	S_DAU1	SEX	S_DAU2	2	2	0.6464	0.3530	0.6762	$
```

> **The header is wrong, and open-king copies it wrong on purpose.** It names **eleven**
> columns; every data row carries **nine values and a trailing tab**. `MaxIBD1` and
> `MaxIBD2` are never written, so the three numbers a row does carry — `IBD1Seg`,
> `IBD2Seg`, `PropIBD` — sit under the headings `MaxIBD1`, `MaxIBD2`, `IBD1Seg`.
> Reproducing the misalignment is deliberate: it is what the reference emits, and a tool
> built against the reference will already be compensating for it.

So the real layout is:

| # | header says | actually holds | format |
| ---: | --- | --- | --- |
| 1–4 | `FID1` `ID1` `FID2` `ID2` | the pair | text |
| 5–6 | `Sex1` `Sex2` | `.fam` sex codes, **raw** — samples with sex outside `{1,2}` are *not* excluded here, unlike in `--kinship`'s `X.kin` | `%d` |
| 7 | `MaxIBD1` | **`IBD1Seg`** — proportion of the usable X shared IBD1 | `%.4f` |
| 8 | `MaxIBD2` | **`IBD2Seg`** | `%.4f` |
| 9 | `IBD1Seg` | **`PropIBD`** — full precision, unlike `.seg`'s | `%.4f` |
| 10–11 | `IBD2Seg` `PropIBD` | nothing; the row ends with a tab | — |

**Rows are exactly `.seg`'s rows, in exactly its order** — the autosomal pass chose the
pairs, including applying `--degree` to the *autosomal* `PropIBD`. Nothing re-screens a
pair on its X evidence, so a pair with no X sharing at all still gets a row of zeros. An
empty `.seg` gives a header-only `X.seg`.

**How to read this.** `S_F` is the father. Against his sons: all zeros — a father and son
share no X. Against his daughters: `1.0000 / 0.0000 / 0.5000` — he gives each daughter his
entire single X. The two brothers (`S_SON1`/`S_SON2`) come out `0.1462 / 0.6393` — they drew
much of the same maternal X. This is why the X table is worth reading separately: on the
autosomes those three pairs look identical.

---

## Empty, header-only, truncated, absent

Four distinguishable states, and the difference matters when you script around these files.

| file | absent when | zero bytes when | header-only when |
| --- | --- | --- | --- |
| `.kin` | no family has ≥ 2 members | one distinct FID and under 64 KiB of rows | — |
| `.kin0` (`--kinship`) | fewer than 2 distinct FIDs | — | — |
| `.kin0` (`--related`) | nothing detected; and always below 100 samples at `--degree` 1–2 | — | — |
| `.ibs` | never | never | no within-family pair exists |
| `.ibs0` | fewer than 2 distinct FIDs | — | — |
| `.con` | N ≥ 100 and no duplicate found | — | N < 100 and no duplicate found (65 bytes) |
| `build.log` | tiny dataset (< 10 samples) | nothing reconstructed | — |
| `updateparents.txt` | tiny dataset | nothing reconstructed | — |
| `updateids.txt` | no cluster merge (announced anyway by `--build`) | — | — |
| `cluster.kin` | no cluster merge, or no informative IBD segment | — | — |
| `_autoQC_updatesex.txt` | no sample has `.fam` sex `0` | — | — |
| `allsegs.txt` | the map yields no usable segment | — | — |

"Tiny dataset" means **N < 10**, at which point `--build` and `--cluster` write nothing at
all, `allsegs.txt` included — the run stops at
`This function is currently disabled for tiny dataset with sample size < 10.` `--unrelated`
is the exception: it takes the same early exit and still writes both of its lists.

### The single-family `.kin` is truncated, not empty

**If your dataset has exactly one distinct FID, `.kin` loses its tail.** Rows are buffered
and flushed every 65 536 bytes, and the final partial buffer is never written. Under 64 KiB
of content that means a **zero-byte file**:

```bash
king -b /tmp/kingdocs/nuclear.bed --kinship --prefix nuc   # 6 samples, one family
wc -c nuc.kin
```

```
       0 nuc.kin
```

Over 64 KiB it means a file that ends, on a line boundary, at the last flushed chunk.
Relabelling all 200 `bigish` samples into one family:

```bash
awk '{$1="ONE"; print}' /tmp/kingdocs/bigish.fam > one.fam
king -b /tmp/kingdocs/bigish.bed --fam one.fam --kinship --cpus 4 --prefix one
echo "one.kin: $(wc -c < one.kin) bytes, $(( $(wc -l < one.kin) - 1 )) rows of $(( 200*199/2 ))"
```

```
one.kin:  1180057 bytes, 19546 rows of 19900
```

354 pairs are missing, and the byte count is 18 flushes' worth. **With two or more distinct
FIDs the file is complete regardless of size** — the same 200 samples in their real
families give a complete `.kin`.

This is a reference bug that open-king reproduces deliberately, because the alternative is
silently disagreeing with KING. If you have one family and want the whole table, split the
FID, or use `--ibs`, which is never truncated.

### No within-family pairs at all

```bash
awk '{$1="S"NR; print}' /tmp/kingdocs/unrelated.fam > sing.fam
king -b /tmp/kingdocs/unrelated.bed --fam sing.fam --kinship --prefix k
king -b /tmp/kingdocs/unrelated.bed --fam sing.fam --ibs     --prefix i
wc -c i.ibs i.ibs0 iallsegs.txt k.kin0; ls k.kin
```

```
     139 i.ibs
   49313 i.ibs0
    1226 iallsegs.txt
   18779 k.kin0
   69457 total
ls: k.kin: No such file or directory
```

There is no `k.kin` at all, and `i.ibs` is exactly its 139-byte header.

### `--related` on a small dataset is silently `--kinship`

At **N ≤ 9** the console prints

```
--related is replaced with --kinship for a small sample size.
```

and the run is identical to `--kinship`: a 10-column `.kin`, an 8-column `.kin0`, no
segment columns. Check the header before assuming you have the 16-column form.

---

## Known divergences from KING 2.3.2

Complete, current, and quantified in [`PARITY.md`](PARITY.md). The short list, as it
affects output files:

* **`<prefix>build.log` is byte-identical in all captured cases.** Rare constructed
  pedigrees still expose repetition/trigger and cross-family-parent-renaming differences;
  the cached held-out replay is 277/347 whole-log exact (§6.2).
* **`--ibdseg` applies the reference's closed 100 Mb usable-total floor.** Below that
  floor both binaries print `Segments too short.` and suppress `.seg`; exactly
  100,000,000 bp proceeds. `--ibs` uses the same floor (§5.10).
* **`<prefix>splitped.txt` is conditional**: all-parentless singleton families write and
  announce nothing; a family of two or a singleton naming a parent emits it (§5.10).
* **`.fam` `SEX` fields outside `{0,1,2}`** follow the reference's permissive rule:
  `M`/`m` is male, a leading `2` or `F`/`f` is female, and other nonzero numeric prefixes
  except `-9` are male. All 43 measured spellings are pinned (§5.11).
* One deliberate out-of-sample segment safety divergence: 4 value rows in 6 713 over 24
  filesets the corpus has never seen, all caused by KING's uninitialised exact-multiple-of-64
  tail read; open-king has 0 extra and 0 missing rows on that battery (§4.6). One separate
  acceptance-gate counterexample remains (§5.11).
* The segment-unavailable PO/FS path currently prints and uses `0.0050`; held-out reference
  probes show that KING derives a deterministic cutoff from the data, but its rule remains
  unidentified (§5.12 and `VERIFIED_FORMULAS.md`).
* Rare `HomIBS0` exact ties can differ in the last printed digit, and one of seven focused
  `MI_Removal` probes refutes the current approximate predicate. No golden row differs.

Everything else — `.kin`, `.kin0`, `X.kin`, `X.kin0`, `.seg`, `X.seg`, `cluster.kin`,
`.ibs`, `.ibs0`, `.con`, `allsegs.txt`, `splitped.txt`, `unrelated.txt`,
`unrelated_toberemoved.txt`, `updateids.txt`, `updateparents.txt`, `bySample.txt`,
`bySNP.txt` and the four autoQC reports — is byte-identical to KING 2.3.2 in every corpus
case that produces it.

## See also

* [`README.md`](README.md) — the documentation index
* [`CLI.md`](CLI.md) — every command-line option, and which analyses each one affects
* [`INTERPRETING.md`](INTERPRETING.md) — how to read the numbers once you have them
* [`COOKBOOK.md`](COOKBOOK.md) — task-oriented recipes
* [`PARITY.md`](PARITY.md) — what is byte-identical to KING 2.3.2, measured, per file and
  per row
* [`VERIFIED_FORMULAS.md`](VERIFIED_FORMULAS.md) — the estimators behind every column, each
  checked numerically against the reference
* [`BEHAVIOR.md`](BEHAVIOR.md) — the experiments that fixed the file-existence, row-order
  and column-set rules quoted above
* [`SPEC.md`](SPEC.md) — the full specification, including the CLI surface and input
  handling
