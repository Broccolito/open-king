# 09 — Runtime probe of KING 2.3.2 (black-box, real data)

Reference binary: `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(macOS arm64, banner `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`).

**Method.** Everything below is observed I/O from running the shipped binary on
hand-built PLINK filesets. No source was read. All artefacts (inputs, every output
file, every captured stdout) live under
`/private/tmp/claude-501/-Users-wgu-Desktop-GeneQuire-Studio--claude-worktrees-king-relatedness-implementation-ac06bd/74b7491e-2c75-4297-8170-f18e23fe8596/scratchpad/probe/tiny/rp/`.

`cat -A` is GNU-only; macOS equivalent used throughout is **`cat -vet`**
(`^I` = TAB, `$` = LF). All dumps below are `cat -vet` output.

---

## 1. Test fixtures

### 1.1 `tiny` — 8 samples, 2000 SNPs (the required fixture)

Generator: `rp/mk.py` (seed 20260813). 2000 SNPs, 1000 on chr1 + 1000 on chr2,
ACGT alleles, MAF drawn U(0.10, 0.50), 1% missing genotypes, alleles actually
transmitted through meiosis (Haldane recombination at 1 cM/Mb) so the pedigree is
genetically real.

`rp/tiny.fam` (space-delimited, produced by `plink --file tiny --make-bed`):

```
fam1 FA1 0 0 1 1$
fam1 MO1 0 0 2 1$
fam1 CH1 FA1 MO1 1 2$
fam1 CH2 FA1 MO1 2 2$
fam2 FA2 0 0 1 1$
fam2 MO2 0 0 2 1$
fam2 CH3 FA2 MO2 1 2$
fam2 DUP1 0 0 1 2$
```

Truth: fam1 = {FA1×MO1 → CH1, CH2} (4 PO pairs, 1 FS pair, 1 unrelated spouse pair).
fam2 = {FA2×MO2 → CH3} plus **DUP1 = an exact genotype copy of fam1/CH1**, so there
is one cross-family MZ/duplicate pair, and CH2↔DUP1 is a cross-family FS pair.

Conversion: `runtimes/macos-arm64/plink1.9/plink --file tiny --make-bed --out tiny`
(PLINK v1.9.0-b.7.11). Genotyping rate 0.990688.

### 1.2 `big` — 164 samples, 2000 SNPs

Generator: `rp/mkbig.py 40` (seed 777). 40 nuclear families `f01`..`f40` of
(FA, MO, CH1, CH2), plus family `famX` = {DUPCH1 (copy of f01/CH1), GC (grandchild
of f01), GGC (great-grandchild), HSIB (half-sib of f01/CH1)}. **Required** because
KING silently downgrades `--related` on small datasets (§3).

### 1.3 `edge` / `thr` / `thr2` / `ord/*` — targeted probes

`rp/edge/` 300 SNPs (100 poly chr1, 50 monomorphic chr1, 100 poly chr2, 50 chrX)
with a zero-overlap pair, an all-missing sample and an unknown-sex sample.
`rp/thr2/` a 4096-SNP ladder of per-sample call counts. `rp/ord/*.fam` are
`--fam` overrides used to reverse-engineer output row ordering.

---

## 2. File inventory per mode

Exit code is **0** on success, **1** on `FATAL ERROR`. All modes were run in their
own empty directory. Output files (`stdout.txt`/`stderr.txt` are my redirections;
KING writes everything to stdout and nothing to stderr):

| invocation | files written (default prefix `king`) |
|---|---|
| `--related` (N ≥ 10) | `king.kin`, `king.kin0`*, `kingallsegs.txt` |
| `--related` (N < 10) | `king.kin`, `king.kin0` — *silently downgraded to `--kinship`* |
| `--kinship` | `king.kin`, `king.kin0` |
| `--duplicate` | `king.con` |
| `--ibs` | `king.ibs`, `king.ibs0`, `kingallsegs.txt`† |
| `--unrelated` | `kingunrelated.txt`, `kingunrelated_toberemoved.txt`, `kingallsegs.txt`† |
| `--bysample` | `kingbySample.txt`, `kingallsegs.txt`† |
| `--bySNP` | `kingbySNP.txt`, `kingallsegs.txt`† |

\* `king.kin0` is **not created at all** by `--related` when no cross-family
relative clears the degree threshold (stdout: `No close relatives are inferred.`).
† `kingallsegs.txt` is written only when at least one chromosomal segment is long
enough for IBD-segment analysis; otherwise stdout says `No informative IBD segments.`
and the file is absent.

### 2.1 `--prefix` semantics

`--prefix P` (default `king`) is a **literal string concatenation**, not a stem +
separator. Observed with `--prefix ZZ_`:

```
--related    -> ZZ_.kin  ZZ_.kin0  ZZ_allsegs.txt
--kinship    -> ZZ_.kin  ZZ_.kin0
--duplicate  -> ZZ_.con
--ibs        -> ZZ_.ibs  ZZ_.ibs0  ZZ_allsegs.txt
--unrelated  -> ZZ_allsegs.txt  ZZ_unrelated.txt  ZZ_unrelated_toberemoved.txt
--bysample   -> ZZ_allsegs.txt  ZZ_bySample.txt
--bySNP      -> ZZ_allsegs.txt  ZZ_bySNP.txt
```

So the fixed suffixes are exactly:
`.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`, `allsegs.txt`, `unrelated.txt`,
`unrelated_toberemoved.txt`, `bySample.txt`, `bySNP.txt`.
Note the dot belongs to the suffix for the five relatedness files and is **absent**
for the five `*.txt` files. The prefix may contain a directory path
(`--prefix sub/run1.` → `sub/run1..kin`); the directory must already exist.
`--kinship --prefix foo` produced `foo.kin`/`foo.kin0` byte-identical to the
default-prefix run.

### 2.2 Encoding

LF line endings only; every file ends with exactly one `\n` and no trailing blank
line (verified with `tail -c 8 | xxd -p`). No BOM. No trailing field separator.
`.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`, `kingallsegs.txt`,
`kingunrelated*.txt` are **TAB**-delimited. `kingbySample.txt` and
`kingbySNP.txt` are **SPACE**-delimited. `kingunrelated*.txt` have **no header**.

---

## 3. `--related` is not always `--related`

With **N < 10 samples** KING prints

```
--related is replaced with --kinship for a small sample size.
```

and from then on behaves exactly as `--kinship` (10-column `.kin`, 8-column
`.kin0`, `--degree` ignored). Boundary measured by subsetting `big` with
`plink --keep`:

| N | downgraded? | `.kin` cols | `.kin0` cols |
|---|---|---|---|
| 8 | yes | 10 | 8 |
| 9 | yes | 10 | 8 |
| 10 | **no** | 16 | 14 |
| 12…164 | no | 16 | 14 |

`--unrelated` has the same N<10 gate, with a different message:
`This function is currently disabled for tiny dataset with sample size < 10.`
(it still writes both `kingunrelated*.txt` files, just without family clustering).

**Consequence for the task's tiny fixture:** the required `--related` and
`--related --degree 2` runs on 8 samples produce files byte-identical to
`--kinship`. The 16/14-column `--related` formats below were captured on `big`.

---

## 4. `.kin` — within-family pairs

### 4.1 `--kinship` form (10 columns)

Header, then rows, `cat -vet`, from `rp/r_kinship/king.kin` (tiny):

```
FID^IID1^IID2^IN_SNP^IZ0^IPhi^IHetHet^IIBS0^IKinship^IError$
fam1^ICH1^ICH2^I1969^I0.250^I0.2500^I0.2103^I0.0142^I0.2304^I0$
fam1^ICH1^IFA1^I1957^I0.000^I0.2500^I0.1937^I0.0000^I0.2515^I0$
fam1^ICH1^IMO1^I1964^I0.000^I0.2500^I0.1950^I0.0000^I0.2474^I0$
fam1^IFA1^IMO1^I1953^I1.000^I0.0000^I0.1726^I0.0886^I-0.0059^I0$
fam2^ICH3^IDUP1^I1975^I1.000^I0.0000^I0.1691^I0.0927^I-0.0208^I0$
```

Literal header line: `FID<TAB>ID1<TAB>ID2<TAB>N_SNP<TAB>Z0<TAB>Phi<TAB>HetHet<TAB>IBS0<TAB>Kinship<TAB>Error`

### 4.2 `--related` form (16 columns)

From `rp/b_related/king.kin` (big):

```
FID^IID1^IID2^IN_SNP^IZ0^IPhi^IHetHet^IIBS0^IHetConc^IHomIBS0^IKinship^IIBD1Seg^IIBD2Seg^IPropIBD^IInfType^IError$
famX^IDUPCH1^IGC^I1985^I1.000^I0.0000^I0.2005^I0.0000^I0.3350^I0.0000^I0.2509^I1.0000^I0.0000^I0.5000^IPO^I1$
famX^IDUPCH1^IGGC^I1985^I1.000^I0.0000^I0.1844^I0.0443^I0.3045^I0.2708^I0.1212^I0.4917^I0.0000^I0.2458^I2nd^I1$
famX^IDUPCH1^IHSIB^I1988^I1.000^I0.0000^I0.1761^I0.0609^I0.2846^I0.3623^I0.0684^I0.2641^I0.0000^I0.1321^I3rd^I1$
f01^ICH1^ICH2^I1975^I0.250^I0.2500^I0.2309^I0.0041^I0.4007^I0.0268^I0.2760^I0.8263^I0.1085^I0.5217^IFS^I0$
```

`--kinship --degree N` still yields the 10-column form (`--degree` only filters
`.kin0`); `--related --degree N` still yields the 16-column form.

### 4.3 Column semantics / formats

| col | name | format | notes |
|---|---|---|---|
| 1 | `FID` | string | family id |
| 2 | `ID1` | string | first member of the pair in sorted order (§8) |
| 3 | `ID2` | string | |
| 4 | `N_SNP` | `%d` | count of autosomal SNPs non-missing in **both**; monomorphic SNPs are counted |
| 5 | `Z0` | `%.3f` | **pedigree-expected** IBD0 prob, not estimated |
| 6 | `Phi` | `%.4f` | **pedigree-expected** kinship, not estimated |
| 7 | `HetHet` | `%.4f` | |
| 8 | `IBS0` | `%.4f` | |
| 9 | `HetConc` | `%.4f` | `--related` only |
| 10 | `HomIBS0` | `%.4f` | `--related` only |
| 11 (9) | `Kinship` | `%.4f` | estimated; goes negative, e.g. `-0.0059` |
| 12 | `IBD1Seg` | `%.4f` | `--related` only |
| 13 | `IBD2Seg` | `%.4f` | `--related` only |
| 14 | `PropIBD` | `%.4f` | `--related` only |
| 15 | `InfType` | string | `--related` only; one of `Dup/MZ` `PO` `FS` `2nd` `3rd` `4th` `UN` |
| last | `Error` | `%G` | **not** fixed-width: observed values are exactly `0`, `0.5`, `1` |

`Error` is the only column that varies its decimal count — census over 247 rows of
`b_kinship/king.kin`: 217×`0`, 16×`0.5`, 13×`1`. That is `%G`/`%g`-style output
(trailing zeros suppressed), **not** `%.1f`. There is no `Error` column in `.kin0`
or in `.ibs`/`.ibs0`.

Pedigree-derived `Z0`/`Phi` observed values (probe `rp/ord/PED3B`, a 3-generation
pedigree with half-sibs):

| pedigree relationship | `Z0` | `Phi` |
|---|---|---|
| parent–offspring | `0.000` | `0.2500` |
| full sibs | `0.250` | `0.2500` |
| grandparent / avuncular / half-sib | `0.500` | `0.1250` |
| unrelated (incl. spouses) | `1.000` | `0.0000` |

### 4.4 Row selection and quirks

- Pairs with `N_SNP == 0` are **omitted entirely** (verified in `rp/edge`: samples
  A and B have disjoint call sets; the A–B row is absent while A–C and B–C are present).
- A sample with **no** non-missing genotypes never appears in any pair.
- **Single-family datasets:** if the whole fileset is one family, KING prints
  `There is only one family.`, writes **`king.kin` as a 0-byte file** (not even the
  header) and writes no `king.kin0`. Reproduced twice (`rp/ord/A`, `rp/ord/C`).
  Any reimplementation should decide deliberately whether to bug-for-bug match this.

---

## 5. `.kin0` — cross-family pairs

### 5.1 `--kinship` form (8 columns) — all pairs

```
FID1^IID1^IFID2^IID2^IN_SNP^IHetHet^IIBS0^IKinship$
fam1^IFA1^Ifam2^IFA2^I1953^I0.1669^I0.0896^I-0.0274$
fam1^IFA1^Ifam2^IMO2^I1948^I0.1679^I0.0862^I-0.0213$
fam1^ICH1^Ifam2^IDUP1^I1968^I0.3857^I0.0000^I0.5000$
fam1^ICH2^Ifam2^IDUP1^I1969^I0.2092^I0.0147^I0.2218$
```

Literal header: `FID1<TAB>ID1<TAB>FID2<TAB>ID2<TAB>N_SNP<TAB>HetHet<TAB>IBS0<TAB>Kinship`

On `big` this file has 13120 data rows = C(164,2) − 246 within-family pairs, i.e.
**every** cross-family pair, unfiltered.

### 5.2 `--related` form (14 columns) — only pairs above the degree threshold

```
FID1^IID1^IFID2^IID2^IN_SNP^IHetHet^IIBS0^IHetConc^IHomIBS0^IKinship^IIBD1Seg^IIBD2Seg^IPropIBD^IInfType$
f01^IFA^IfamX^IDUPCH1^I1985^I0.2060^I0.0000^I0.3484^I0.0000^I0.2545^I1.0000^I0.0000^I0.5000^IPO$
f01^ICH1^IfamX^IDUPCH1^I1984^I0.4047^I0.0000^I1.0000^I0.0000^I0.5000^I0.0000^I1.0000^I1.0000^IDup/MZ$
f01^ICH2^IfamX^IDUPCH1^I1983^I0.2305^I0.0040^I0.3998^I0.0267^I0.2751^I0.8263^I0.1085^I0.5217^IFS$
```

Same formats as §4.3 (all `%.4f`; `N_SNP` `%d`; `InfType` string; **no `Error`
column**).

### 5.3 `--degree` thresholds

`--degree d` filters `.kin0` at kinship ≥ 2^−(d+1.5), printed in stdout with `%.5f`:

| `--degree` | printed threshold | 2^−(d+1.5) | `.kin0` rows on `big` |
|---|---|---|---|
| (default) / 0 | `0.17678` | 0.1767767 | 7 |
| 1 | `0.17678` | 0.1767767 | 7 |
| 2 | `0.08839` | 0.0883883 | 12 |
| 3 | `0.04419` | 0.0441942 | 59 |
| 4 | `0.02210` | 0.0220971 | 965 |
| 5 | `0.01105` | 0.0110485 | 2455 |

`--degree 0` is treated as 1. `--kinship --degree d` applies the same threshold to
the 8-column `.kin0` (7 / 12 / 59 rows for d = 1/2/3) and leaves `.kin` untouched.

Note the count in stdout (`N pairs of relatives (up to Nnd-degree) are identified`)
counts **classified** relatives, which is smaller than the number of rows written
(13 identified vs 59 rows at degree 3): rows above the raw kinship threshold but
classified `UN` are still written.

### 5.4 The `M ≤ 512` sample exclusion

Before the cross-family stage, KING drops samples with too few called autosomal
SNPs and prints:

```
The following 4 samples are excluded from the kinship analysis (M<512):$
^I(A full1)^I(A n513)^I(A n512)^I(B n511)$
```

i.e. a header line, then one line beginning with TAB holding TAB-separated
`(FID IID)` tuples.

Ladder measurement (`rp/thr2`, 4096 SNPs, per-sample call counts
4096/2048/1024/600/540/520/513/512/511/500/400/300/256/200/128/64):
samples with **M ≥ 513** appear in `.kin0`; samples with **M ≤ 512** do not.
So the real predicate is `M < 513` (`M <= 512` excluded), despite the `M<512` text.

**The printed list is wrong.** The *count* is correct (9 samples excluded in the
ladder run) but the names printed are simply the first 9 samples in `.fam` serial
order, including `m4096a`/`m4096b` which are demonstrably present in `.kin0`. This
is a display bug in KING 2.3.2; only the count is trustworthy.

Excluded samples still appear in `.kin` (the within-family stage does not apply the
filter) — in `rp/edge` all 8 samples were listed as excluded yet `king.kin` has 5
rows and `king.kin0` has a header and nothing else.

---

## 6. `.con` — `--duplicate`

```
FID1^IID1^IFID2^IID2^IN^IN_IBS0^IN_IBS1^IN_IBS2^IConcord^IHomConc^IHetConc$
fam1^ICH1^Ifam2^IDUP1^I1968^I0^I0^I1968^I1.00000^I1.00000^I1.00000$
```

From `big` (three pairs, showing that within- and cross-family pairs share one file
and that near-duplicates are reported too):

```
f01^ICH1^IfamX^IDUPCH1^I1984^I0^I0^I1984^I1.00000^I1.00000^I1.00000$
f03^ICH1^If03^ICH2^I1978^I0^I9^I1969^I0.99545^I1.00000^I0.98810$
f26^ICH1^If26^ICH2^I1981^I0^I137^I1844^I0.93084^I1.00000^I0.83252$
```

Literal header: `FID1<TAB>ID1<TAB>FID2<TAB>ID2<TAB>N<TAB>N_IBS0<TAB>N_IBS1<TAB>N_IBS2<TAB>Concord<TAB>HomConc<TAB>HetConc`

Formats: `N`, `N_IBS0`, `N_IBS1`, `N_IBS2` are `%d`; **`Concord`, `HomConc`,
`HetConc` are `%.5f`** — the only 5-decimal columns KING emits.

Row order is `.fam` **serial** order with an i<j double loop (verified with
`--minConc 0.1` on tiny, which emitted 28 rows starting
`fam1 FA1 fam1 MO1`, `fam1 FA1 fam1 CH1`, …). Selection threshold is `--minConc`
(default 0.80) applied to `HetConc`. If nothing passes, the file is written with the
**header only** and stdout says `No duplicates are found with heterozygote
concordance rate > 80%.`; otherwise `N pairs of duplicates with heterozygote
concordance rate > 80% are saved in file king.con` (the percentage is
`minConc × 100` rendered without decimals: `80`, `10`).

---

## 7. `.ibs` / `.ibs0` — `--ibs`

Full form, 22 / 21 columns (`big`):

```
FID^IID1^IID2^IZ0^IPhi^IN_SNP^IN_IBS0^IN_IBS1^IN_IBS2^INHetHet^INHomHom^IN_Het1^IN_Het2^IIBS^IDist^IHetConc^IHet2|1^IHet1|2^IHomConc^IKinship^IMaxIBD2^IPr_IBD2$
famX^IDUPCH1^IGC^I1.000^I0.0000^I1985^I0^I790^I1195^I398^I797^I807^I779^I1.6020^I0.3980^I0.3350^I0.4932^I0.5109^I1.0000^I0.2509^I0.000^I0.0000$
famX^IDUPCH1^IGGC^I1.000^I0.0000^I1985^I88^I836^I1061^I366^I783^I802^I766^I1.4902^I0.5985^I0.3045^I0.4564^I0.4778^I0.8876^I0.1212^I0.000^I0.0000$
```

```
FID1^IID1^IFID2^IID2^IN_SNP^IN_IBS0^IN_IBS1^IN_IBS2^INHetHet^INHomHom^IN_Het1^IN_Het2^IIBS^IDist^IHetConc^IHet2|1^IHet1|2^IHomConc^IKinship^IMaxIBD2^IPr_IBD2$
f01^IFA^If02^IFA^I1983^I159^I928^I896^I299^I756^I776^I750^I1.3717^I0.7887^I0.2437^I0.3853^I0.3987^I0.7897^I-0.0213^I-9^I-9$
f01^ICH1^IfamX^IDUPCH1^I1968^I0^I0^I1968^I759^I1209^I759^I759^I2.0000^I0.0000^I1.0000^I1.0000^I1.0000^I1.0000^I0.5000^I48941787.000^I0.9627$
```

Note the header contains **pipe characters** in `Het2|1` and `Het1|2`.

`.ibs` header (literal): `FID<TAB>ID1<TAB>ID2<TAB>Z0<TAB>Phi<TAB>N_SNP<TAB>N_IBS0<TAB>N_IBS1<TAB>N_IBS2<TAB>NHetHet<TAB>NHomHom<TAB>N_Het1<TAB>N_Het2<TAB>IBS<TAB>Dist<TAB>HetConc<TAB>Het2|1<TAB>Het1|2<TAB>HomConc<TAB>Kinship<TAB>MaxIBD2<TAB>Pr_IBD2`

`.ibs0` header is the same but with `FID1 ID1 FID2 ID2` in place of `FID ID1 ID2`
and **without `Z0` and `Phi`** (those are pedigree quantities, meaningless across
families).

**Dynamic columns:** when there are no chromosomal segments long enough for IBD
analysis (stdout `No informative IBD segments.`), the trailing `MaxIBD2` and
`Pr_IBD2` columns are **removed from the header and from every row**, giving 20 /
19 columns. Captured on `rp/edge` (chromosomes only ~25 Mb and ~10 Mb long):

```
FID^IID1^IID2^IZ0^IPhi^IN_SNP^IN_IBS0^IN_IBS1^IN_IBS2^INHetHet^INHomHom^IN_Het1^IN_Het2^IIBS^IDist^IHetConc^IHet2|1^IHet1|2^IHomConc^IKinship$
fE1^IA^IC^I1.000^I0.0000^I100^I13^I44^I43^I22^I34^I46^I42^I1.3000^I0.9600^I0.3333^I0.4783^I0.5238^I0.6176^I-0.0455$
```

Formats: `Z0` `%.3f`; `Phi` `%.4f`; `N_SNP`…`N_Het2` `%d`; `IBS`, `Dist`,
`HetConc`, `Het2|1`, `Het1|2`, `HomConc`, `Kinship` `%.4f`;
**`MaxIBD2` `%.3f`, `Pr_IBD2` `%.4f`** — but only for pairs that got an IBD-segment
analysis. For every other pair in `.ibs0` the two fields are the bare literal
strings **`-9`** and **`-9`** (not `-9.000` / `-9.0000`). In `.ibs` (within-family)
unanalysed pairs print `0.000` / `0.0000` instead. `MaxIBD2` is a base-pair
magnitude, so values like `48941787.000` occur — this column is emphatically
**not** `%g`.

---

## 8. Row ordering (this is the hard part)

Two different orders are in play, and getting them wrong breaks byte parity.

- **`.kin0`, `.con`, `kingbySample.txt`** use **`.fam` serial order**, with an i<j
  double loop for pairs.
- **`.kin` and `.ibs`** use a **sorted** order: families sorted, and members sorted
  within each family, then an i<j double loop over the sorted member list.

The comparator was reverse-engineered by feeding KING crafted `--fam` files
(`rp/ord/`). Input `f01,f02,famX,famY` → output `famX,famY,f01,f02`; input
`zz,aa,M1,0` → `aa,M1,zz,0`; input `b10,b9,b2,b1` → `b1,b2,b9,b10`.

A 40-value probe (`rp/ord/D`, family IDs) and an independent 36-value probe
(`rp/ord/E`, person IDs) produced **the same total order**, so one comparator
serves both:

```
input : a1 ab 1a a A0 B 9 09 10 2 z_ z- z. z0 zA z9 z10 z09 _a -a 0a x1y x1z
        x01y x2y abc ABD aBc2 Q q1 007 7 0 00 m~ m} m! m# mZ m0
output: -a a ab abc aBc2 ABD A0 a1 B m! m# mZ m} m~ m0 Q q1 x1y x1z x2y x01y
        z- z. zA zzz41..zzz82 z_ z0 z9 z09 z10 _a 0 0a 1a 2 7 9 00 09 10 007
```

Rules that reproduce it exactly:

1. **Case-insensitive, folding to UPPERCASE.** (`zA` sorts before `z_` because
   `'A'`=65 < `'_'`=95; `zzz…` sorts before `z_` because `'z'`→`'Z'`=90 < 95.
   Confirms uppercase folding rather than lowercase.) Two ids differing only in
   case are *equal*: KING rejects the fileset with
   `Family F: Person a is duplicated` when both `a` and `A` exist in one family.
2. **Digits sort after every non-digit.** `ABD` < `A0`; `zA` < `z0`; `_a` < `0`.
   Within the `m` group: `m!`(0x21) < `m#`(0x23) < `mZ` < `m}`(0x7D) < `m~`(0x7E)
   < `m0` — pure ASCII on the uppercased char, with digits displaced to the top.
3. **Runs of digits compare as a unit: longer run first, then digit-by-digit.**
   `x1y` < `x1z` < `x2y` < `x01y` (the two-digit run `01` beats every one-digit
   run), and `b1 b2 b9 b10`, and `0 < 0a < 1a < 2 < 7 < 9 < 00 < 09 < 10 < 007`.
   This is a natural sort that does **not** normalise leading zeros.
4. **End of string sorts before any character** (`a` < `ab` < `abc`).

Everything else in the comparison is ordinary ASCII on the uppercased byte.

---

## 9. `kingbySample.txt` — `--bysample` (SPACE-delimited)

Full form:

```
FID IID FA MO SEX N_SNP Missing Heterozygosity N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal$
f01 FA 0 0 1 1989 0.0055 0.3917 3954 0 0.0000 3936 0 0.0000 0$
f01 MO 0 0 2 1991 0.0045 0.3867 3957 0 0.0000 3936 0 0.0000 0$
fam2 DUP1 0 0 1 1984 0.0080 0.3841 0 0 0.0000 0 0 0.0000 0$
```

**The header is dynamic.** Measured by holding the genotypes fixed and swapping in
different `--fam` files (`rp/hdr/`):

| X SNPs present? | PO pairs? | trios? | header |
|---|---|---|---|
| no | no | no | `FID IID FA MO SEX N_SNP Missing Heterozygosity` |
| no | yes | no | `… Heterozygosity N_pair N_MIp Err_MIp MI_Removal` |
| no | yes | yes | `… Heterozygosity N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal` |
| yes | no | no | `… Heterozygosity N_xSNP xHeterozygosity` |
| yes | yes | no | `… Heterozygosity N_xSNP xHeterozygosity N_pair N_MIp Err_MIp MI_Removal` |
| yes | yes | yes | `… Heterozygosity N_xSNP xHeterozygosity N_pair N_MIp Err_MIp N_trio N_MIt Err_MIt MI_Removal` |

Block rules: `N_xSNP xHeterozygosity` appear iff the bim contains X-chromosome
SNPs; `N_pair N_MIp Err_MIp` iff the pedigree has ≥1 parent–offspring pair;
`N_trio N_MIt Err_MIt` iff it has ≥1 complete trio; `MI_Removal` iff there is ≥1
PO pair.

Formats: `FID`, `IID`, `FA`, `MO` strings (`0` for a missing parent); `SEX`,
`N_SNP`, `N_xSNP`, `N_pair`, `N_MIp`, `N_trio`, `N_MIt` `%d`;
`Missing`, `Heterozygosity`, `xHeterozygosity`, `Err_MIp`, `Err_MIt` `%.4f`;
**`MI_Removal` is `%G`** — observed values `0` and `0.5` (rp/edge).

Rows are in `.fam` serial order, one per sample, including samples with no
genotypes at all (`fE3 H 0 0 2 0 1.0000 0.0000 …`).

---

## 10. `kingbySNP.txt` — `--bySNP` (SPACE-delimited)

Full form:

```
SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate N_PO N_HomPO N_errPO Err_InPO Err_InHomPO N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio$
rs10001 1 135029 T C 0.4785 163 37 82 44 0.9939 158 42 0 0.0000 0.0000 78 38 0 0.0000 0.0000$
rs10002 1 192745 G C 0.2927 164 13 70 81 1.0000 160 63 0 0.0000 0.0000 80 33 0 0.0000 0.0000$
```

Dynamic header, three variants:

```
no PO pairs : SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate
PO, no trio : … CallRate N_PO N_HomPO N_errPO Err_InPO Err_InHomPO
PO + trio   : … Err_InHomPO N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio
```

Formats: `Chr`, `Pos`, `N`, `N_AA`, `N_Aa`, `N_aa`, `N_PO`, `N_HomPO`, `N_errPO`,
`N_trio`, `N_HetOff`, `N_errTrio` `%d`; `Freq_A`, `CallRate`, `Err_InPO`,
`Err_InHomPO`, `Err_InTrio`, `Err_InHetTrio` `%.4f`.

Other observations:
- **`Chr` is rendered symbolically**: chromosome 23 in the bim prints as `X`
  (`x_000 X 100000 G T …`). Autosomes print as plain integers.
- `--bySNP` covers **all** chromosomes including X (301 lines for 300 SNPs), while
  the relatedness modes use autosomes only.
- A monomorphic SNP whose bim `A1` is `0` prints `Label_A` as `0` and
  `Freq_A 0.0000`.
- Rows are in bim order.

---

## 11. `kingallsegs.txt`

```
Segment^IChr^IStartMB^IStopMB^ILength^IN_SNP^IStartSNP^IStopSNP$
1^I1^I0.135^I51.316^I51.181^I1000^Irs10001^Irs11000$
2^I2^I0.124^I49.628^I49.504^I1000^Irs20001^Irs21000$
```

Literal header: `Segment<TAB>Chr<TAB>StartMB<TAB>StopMB<TAB>Length<TAB>N_SNP<TAB>StartSNP<TAB>StopSNP`
Formats: `Segment`, `Chr`, `N_SNP` `%d`; `StartMB`, `StopMB`, `Length` `%.3f`;
`StartSNP`, `StopSNP` are rs ids.
Segments are numbered from 1. Stdout companion line:
`Total length of 2 chromosomal segments usable for IBD segment analysis is 100.8 Mb.`
(`%.1f` on the Mb total).

---

## 12. `kingunrelated.txt` / `kingunrelated_toberemoved.txt` — `--unrelated`

**No header line.** Two TAB-separated columns, `FID<TAB>IID`:

```
famX^IHSIB$
famX^IGC$
f01^IFA$
```

```
famX^IDUPCH1$
famX^IGGC$
f01^ICH1$
```

The two files partition the sample set: on `tiny`, 5 kept + 3 removed = 8; on
`big`, sizes 116 + 48 = 164. Row order is neither `.fam` order nor the sorted
order — it is the order in which the greedy selection visits individuals, and it is
reproducible run-to-run.

Stdout:
`A list of 5 unrelated individuals saved in file kingunrelated.txt` /
`An alternative list of 3 to-be-removed individuals saved in file kingunrelated_toberemoved.txt`.

---

## 13. stdout structure

Every run prints the banner, the full parameter block (with `[ON]` marking the
selected mode), `KING starts at <ctime>`, the loading trace, an `Options in effect:`
block (mode flags each on their own TAB-indented line), the per-mode body, and
`KING ends at <ctime>`. `--bySNP` echoes its flag lower-cased as `--bysnp`.

The relationship-summary tables are TAB-separated:

```
Relationship summary (total relatives: 7 by pedigree, 7 by inference)$
  Source^IMZ^IPO^IFS^I2nd^I3rd^IOTHER$
  ===========================================================$
  Pedigree^I0^I6^I1^I0^I0^I5$
  Inference^I0^I6^I1^I0^I0^I5$
```

The cross-family variant has a different first cell (8 spaces, no `Source`), a
shorter rule (57 `=` vs 59) and `4th` instead of `OTHER`:

```
Relationship summary (total relatives: 0 by pedigree, 7 by inference)$
        ^IMZ^IPO^IFS^I2nd^I3rd^I4th$
  =========================================================$
  Inference^I1^I5^I1^I0^I0^I0$
```

Typo worth matching or fixing deliberately: degrees ≥ 3 print as `3nd-degree`,
`4nd-degree`, `5nd-degree` (only `1st`/`2nd` are correct).

`Loading genotype data` prints `Genotype data consist of 250 autosome SNPs,
50 X-chromosome SNPs` — the X clause appears only when X SNPs exist.

---

## 14. Failure modes and non-determinism

- Missing/unopenable bed → exit 1, `FATAL ERROR - \nGenotype file <path> cannot be opened`.
- No genotype file at all → exit 1, usage block + `Genotype files are required. e.g.,\n  king -b ex.bed --related`.
- Duplicate person id within a family (case-insensitively) → exit 1,
  `Family F: Person a is duplicated` then `FATAL ERROR - \nPlease correct problems with pedigree structure`.
- A person id colliding with the missing-parent code `0` produces
  `Parental sex codes don't make sense for Person X in Family Y` and a fatal exit.
- KING leaves temp files `king$TMP$.dat` / `king$TMP$.ped` in the cwd when `--fam`
  is used and the run aborts.

**Non-deterministic spurious fatal error.** On the 8-sample fixture, `--bySNP`
aborts at random with

```
FATAL ERROR - 

Too many first alleles as the major allele (~10.8%). Please use plink1.9 --make-bed to regenerate the genotype data again.
```

Measured rate: **9/40 runs with `tiny.fam`, 9/40 with an all-founders fam**, on
byte-identical input. The reported percentage varies run to run (10.4%–11.0%) while
the true fraction of A1-major SNPs in that fileset is **0.00%** (computed directly
from the bed). `--cpus 1` does not fix it (7/40), so this is uninitialised memory
rather than a thread race. It did **not** occur in 40/40 runs on the 164-sample
fileset, nor in 30/30 `--kinship` runs on the tiny fileset. Swapping the bim allele
columns for 0–50% of SNPs on the 164-sample fileset never triggered it, so the check
is not a straightforward A1-major-fraction gate.

Practical consequence: **do not build regression fixtures at ~8 samples for
`--bySNP`** — retry until the run succeeds, or use ≥ 100 samples. When a run does
succeed, its output is fully deterministic (md5-identical across three successful
`--bySNP` runs; `.kin`/`.kin0` identical across three `--kinship` runs).

---

## 15. Number-format summary (one table)

| format | columns |
|---|---|
| `%d` | `N_SNP`, `N`, `N_IBS0/1/2`, `NHetHet`, `NHomHom`, `N_Het1`, `N_Het2`, `Segment`, `Chr`(autosome), `Pos`, `N_AA`, `N_Aa`, `N_aa`, `N_PO`, `N_HomPO`, `N_errPO`, `N_trio`, `N_HetOff`, `N_errTrio`, `SEX`, `N_xSNP`, `N_pair`, `N_MIp`, `N_MIt` |
| `%.3f` | `Z0`, `StartMB`, `StopMB`, `Length`, `MaxIBD2` |
| `%.4f` | `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `IBS`, `Dist`, `Het2\|1`, `Het1\|2`, `HomConc`, `Pr_IBD2`, `Missing`, `Heterozygosity`, `xHeterozygosity`, `Err_MIp`, `Err_MIt`, `Freq_A`, `CallRate`, `Err_InPO`, `Err_InHomPO`, `Err_InTrio`, `Err_InHetTrio` |
| `%.5f` | `Concord`, `HomConc`, `HetConc` (**`.con` only**) |
| `%G` | `Error` (`.kin`), `MI_Removal` (`kingbySample.txt`) — values `0`, `0.5`, `1` |
| literal `-9` | `MaxIBD2`, `Pr_IBD2` in `.ibs0` for pairs without IBD-segment analysis |
| symbolic | `Chr` = `X` for chromosome 23 in `kingbySNP.txt` |

Negative values print with a leading `-` and no space padding (`-0.0059`); there is
no field-width padding anywhere — every field is separated by a single TAB (or
single space in the two QC files) with no alignment.

---

## 16. Checklist for the reimplementation

1. Two `.kin` shapes (10-col kinship, 16-col related) and two `.kin0` shapes
   (8-col, 14-col). Pick by mode, and remember `--related` degrades to `--kinship`
   below 10 samples.
2. `.ibs`/`.ibs0` drop `MaxIBD2`/`Pr_IBD2` when there are no informative segments.
3. `kingbySample.txt` and `kingbySNP.txt` have six and three header variants
   respectively, driven by X-SNP presence and PO/trio counts.
4. `.kin`/`.ibs` use the sorted order (§8); `.kin0`/`.con`/`kingbySample.txt` use
   serial order.
5. Skip pairs with `N_SNP == 0`.
6. Drop samples with ≤ 512 called autosomal SNPs from the cross-family stage only.
7. Space vs tab delimiters, `%G` for `Error`/`MI_Removal`, `%.5f` in `.con`,
   `-9` sentinels in `.ibs0`.
8. Decide explicitly about the three KING bugs: the 0-byte `.kin` on single-family
   input, the mis-indexed exclusion list, and the `3nd-degree` typo.
