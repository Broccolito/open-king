# Reference-binary behaviour

Companion to [VERIFIED_FORMULAS.md](VERIFIED_FORMULAS.md). That page covers the
*arithmetic*; this one covers everything else the reference binary does that we must
reproduce byte for byte: which SNPs enter the computation, which files get created, which
columns those files carry, how rows are ordered, and how the option flags change all of it.

Everything below was established by running the reference binary
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (KING 2.3.2) on
purpose-built synthetic PLINK filesets. **No KING source was read.** Each section states
the hypothesis space, the discriminating experiment, the raw output, and then the rule.

Where a question could not be settled the section says so explicitly and states what would
settle it.

---

## Contents

* [Q2 — parent–offspring vs full-sibling discrimination](#q2--parentoffspring-vs-full-sibling-discrimination)
* [Q3 — SNP inclusion rules](#q3--snp-inclusion-rules)
* [Q4 — `--cpus` determinism](#q4----cpus-determinism)
* [Q5 — `--degree` semantics](#q5----degree-semantics)
* [Q6 — the sample-ID sort comparator](#q6--the-sample-id-sort-comparator)
* [Q7 — output-file existence](#q7--output-file-existence)
* [Q8 — `.ibs` / `.ibs0` column-set variation](#q8--ibs--ibs0-column-set-variation)
* [Side findings](#side-findings)

---

## Q2 — parent–offspring vs full-sibling discrimination

### Hypothesis space

The binary's string table contains two formats mentioning an IBS0 cutoff:

```
1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf
Cutoff value between full siblings and parent-offspring is set at %.4f
```

`%.4lf` / `%.4f` implies the value is computed, not a literal. Candidate rules: a constant;
a function of the mean IBS0 among 1st-degree pairs; a quantile of that distribution; a
gap/valley finder on it; a function of the allele-frequency spectrum; a function of `N` or
of the SNP count.

### Where the message actually comes from

Neither message appears under `--kinship` or under `--related` on a normal dataset. The
one that governs relatedness is emitted by the **family-clustering** path — `--build`,
`--cluster`, `--unrelated` — and **only when IBD-segment analysis is unavailable**:

```
Sorting autosomes...
No informative IBD segments.
  Inference will be based on kinship estimation only.
8 CPU cores are used to compute the pairwise kinship coefficients...
0%Cutoff value for IBS0 between FS and PO is set at 0.0055
Clustering up to 1st-degree relatives in families...
```

(The observed text is `Cutoff value for IBS0 between FS and PO is set at %.4f`.) The
`1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf` string belongs to
the `--tdt` trio-inference path and was never reached in any of these runs.

When usable IBD segments **do** exist (see [Q8](#q8--ibs--ibs0-column-set-variation)),
`--related`/`--build` distinguish PO from FS by IBD-segment sharing instead: PO pairs come
out with `IBD1Seg 1.0000  IBD2Seg 0.0000  PropIBD 0.5000` exactly, FS pairs with
`IBD2Seg ≈ 0.25`. No IBS0 cutoff is printed on that path.

### How the cutoff is applied — established

Datasets with 1st-degree pairs whose IBS0 proportion is known exactly (built by choosing
the per-pair genotype-pattern counts directly — see the construction in
[Q5](#q5----degree-semantics)):

| dataset | 1st-degree pairs' IBS0 | printed cutoff | summary line `MZ PO FS 2nd 3rd 4th` |
| --- | --- | --- | --- |
| `h_lo`   | 20 pairs, all `0.0100` | `0.0050` | `0 0 20 0 0 0` |
| `h_mid`  | 20 pairs, all `0.0350` | `0.0050` | `0 0 20 0 0 0` |
| `h_hi`   | 20 pairs, all `0.0600` | `0.0050` | `0 0 20 0 0 0` |
| `h_bimodal` | 10 @ `0.0200`, 10 @ `0.0600` | `0.0050` | `0 0 20 0 0 0` |
| `h_wide` | 10 @ `0.0050`, 10 @ `0.0900` | `0.0050` | `0 10 10 0 0 0` |

and, from a genotyping-error sweep on pedigree-simulated data (p = 0.5, 122 samples,
20 000 SNPs, error rate ε applied per genotype):

| ε | PO pairs' IBS0 range | printed cutoff | `PO`/`FS` counts |
| --- | --- | --- | --- |
| 0      | `0.0000`          | `0.0055` | 60 / 15 |
| 0.001  | `0.0001–0.0008`   | `0.0055` | 60 / 15 |
| 0.005  | `0.0020–0.0038`   | `0.0050` | 60 / 15 |
| 0.01   | `0.0040–0.0064`   | `0.0050` | 38 / 37 |
| 0.02   | `0.0081–0.0116`   | `0.0040` |  0 / 75 |

Every one of these is exactly reproduced by the rule below (e.g. at ε = 0.01, exactly the
38 PO pairs with IBS0 < 0.0050 stayed PO; at ε = 0.02 all true PO pairs sit above the
cutoff and *all* were misclassified as FS — the reference makes this mistake).

A finer probe pins the comparison side. 33 designed pairs at kinship 0.25 with IBS0
stepped by 1e-4 from `0.0040` to `0.0072` gave `PO = 11`, i.e. `0.0040 … 0.0050` inclusive
are PO and `0.0051` upward are FS, with the cutoff printed as `0.0050`. So the effective
threshold is at or immediately above `0.0050`; a pair whose IBS0 proportion equals the
printed value is classified **PO**.

> **ESTABLISHED RULE (application).** On the kinship-only clustering path, a pair whose
> kinship estimate falls in the 1st-degree band is labelled **PO** when its IBS0
> proportion (`IBS0 / N_SNP`) is below the cutoff and **FS** otherwise. The comparison
> includes a pair sitting exactly on the printed value.

### How the cutoff value is chosen — NOT established

What was ruled out, each by a discriminating run:

* **Not fitted to the 1st-degree pairs.** The five `h_*` datasets above place the *only*
  1st-degree pairs at IBS0 = 0.005 / 0.010 / 0.020 / 0.035 / 0.060 / 0.090 in five
  different arrangements and all five print the identical cutoff `0.0050`. A
  mean/quantile/valley-finder over those values cannot produce a constant.
* **Not the number of 1st-degree pairs** — 5 pairs vs 150 pairs, both `0.0055`.
* **Not the sample count** — N = 100, 105, 110, 122, 150, 200, 300 all `0.0055`.
* **Not the SNP count** — M = 5 000, 10 000, 20 000, 40 000 all `0.0055`.
* **Not the sample order** — four random `.fam` permutations of one genotype matrix, all
  `0.0055`.
* **Not stochastic** — five independent simulation seeds per configuration give the
  identical value every time.
* **Not a monotone function of MAF, heterozygosity, or the unrelated-pair IBS0 rate.**
  Constant-p datasets at p = 0.02, 0.10, 0.20, 0.30, 0.40, 0.50 all give `0.0055`, but
  p = 0.05 gives `0.0060` and a 25 % @ 0.5 / 75 % @ 0.05 mixture gives `0.0040` — a
  non-monotone response that no single summary statistic tried here explains.

Every value observed across ~50 datasets is an exact multiple of `0.0005` and lies in
`[0.0035, 0.0060]`: `0.0035`, `0.0040`, `0.0050`, `0.0055`, `0.0060`.

> **STATUS: OPEN.** The cutoff is deterministic in the genotype data but its functional
> form is unresolved. **Practical guidance:** it never reaches a file — it only affects
> the PO/FS split of the console *Relationship summary* table on the `--build` /
> `--cluster` / `--unrelated` path, and only when IBD segments are unusable. Using the
> constant `0.0055` reproduces the reference on every ordinary (low-error, HWE-like)
> dataset tested. **What would settle it:** an instrumented sweep that holds every
> candidate statistic fixed but one — in particular a designed sweep of the per-SNP MAF
> spectrum at fixed heterozygosity, since p = 0.05 and the 0.5/0.05 mixture are the only
> configurations that moved the value without genotyping error.

---

## Q3 — SNP inclusion rules

### Experiment design

One fileset per candidate class, each with the class present in a **known count**, 8
samples in 2 families, run with `--kinship` (and re-checked with `--ibs`). Compare the
console's `Genotype data consist of N autosome SNPs` and `PLINK maps loaded: M SNPs`
against each pair's `N_SNP`.

### (a) Non-autosomes — EXCLUDED, but chromosome 25 counts as autosomal

`a_chrom`: 100 SNPs on chr 1, 60 on chr 2, 50 on chr 23, 20 on chr 24, 10 on chr 25,
5 on chr 26 (245 total).

```
  Genotype data consist of 170 autosome SNPs (including 10 XY SNPs), 50 X-chromosome SNPs,
      20 Y-chromosome SNPs, 5 mitochondrial SNPs
  PLINK maps loaded: 245 SNPs
F1 A1 A2  N_SNP=170
```

`170 = 100 + 60 + 10`. X, Y and MT are excluded from the autosomal relatedness
computation; **XY (25) is pooled with the autosomes.**

Chromosome-code spelling (`f_chrname`, `k_chrodd`): `1`–`22` autosome; `X` and `23` → X;
`Y` and `24` → Y; `XY` and `25` → autosome; `MT` and `26` → MT. **Every other code — `0`,
`27`, `30`, `M`, `-1`, `chr1`, `Chr2` — is dropped at map-load time**, so it does not even
appear in `PLINK maps loaded`. (`f_chrname` wrote 190 SNPs of which 13 were on chr `0`;
the console reported `PLINK maps loaded: 177 SNPs`.)

When X-chromosome SNPs are present, `--kinship` additionally writes `<prefix>X.kin` and
`<prefix>X.kin0` from an entirely separate X analysis; the autosomal files are unaffected.

### (b) Monomorphic SNPs — KEPT

`b_mono`: 100 polymorphic + 40 all-hom-A1 + 30 all-hom-A2 + 20 all-het = 190.
Console `190 autosome SNPs`, every pair `N_SNP = 190`. **Monomorphic SNPs are not
dropped**; they contribute to `HomHom` and to `M_ij`.

### (c) SNPs with an allele coded `0` — KEPT

`c_zero`: 100 normal + 25 with `A1 = 0` + 15 polymorphic with `A1 = 0` + 10 with
`A2 = 0` = 150. Console `150 autosome SNPs`, every pair `N_SNP = 150`.
`i_allele` adds 20 SNPs with `A1 == A2 == A` and 10 with `0 0`: all 130 kept.
**The allele columns are never used as a filter.**

### (d) SNP call rate — NO GLOBAL THRESHOLD; missingness is purely pairwise

`d_call`: 100 complete SNPs plus 10 each at call rates 87.5 %, 75 %, 62.5 %, 50 %, 37.5 %,
25 %, plus 5 SNPs missing in everyone (165 total). Missingness is placed in the *trailing*
samples so the reference pair `A1/A2` is never affected.

```
  Genotype data consist of 165 autosome SNPs
F1 A1 A2  N_SNP=160     (165 minus the 5 all-missing)
F1 A1 A3  N_SNP=150     (A3 also missing at the ten 25%-call SNPs)
F1 A1 A4  N_SNP=140     (A4 missing at the 25% and 37.5% SNPs)
```

A SNP called in only 2 of 8 samples still counts for the pair that has it. **There is no
SNP call-rate filter at all on the default relatedness path** (the `--autoQC` strings in
the binary belong to a separate opt-in analysis).

### (e) MAF — NO THRESHOLD

`e_maf` (8 samples): 100 polymorphic + 40 singletons + 20 doubletons → all 160 kept.
`h_maf` (60 samples): 200 polymorphic + 50 singletons (MAF = 1/120 ≈ 0.0083) + 30
monomorphic → console `280 autosome SNPs`, every pair `N_SNP = 280`. Together with (b)
(MAF = 0 retained) this is conclusive.

### Duplicate SNP identifiers / positions — KEPT

`g_dupid` (100 unique + 20 SNPs all named `DUP` + 5 IDs used twice) → 130 loaded, 130 used.
`j_duppos` (20 SNPs sharing one bp position) → 120 loaded, 120 used. The binary's
`%d SNPs are removed for appearing more than once` string was never triggered by either.

> **ESTABLISHED RULE.** The SNP set used for autosomal relatedness is:
> *every* record in the `.bim` whose chromosome code is one of `1`–`22`, `25`, or `XY`,
> in `.bim` file order. Nothing else is filtered — not monomorphic sites, not `0` alleles,
> not low call rate, not low MAF, not duplicate IDs or positions. Records on `23`/`X`,
> `24`/`Y`, `26`/`MT` are held aside (X gets its own analysis); records with any other
> chromosome code are discarded during map load and are not counted in
> `PLINK maps loaded`.
> Missing data is handled **pairwise only**: `M_ij` = SNPs in that set called in both
> members.

---

## Q4 — `--cpus` determinism

**Experiment.** `bigish` (200 samples × 50 000 SNPs, 22 chromosomes) run as
`--related --ibs --cpus C` for `C ∈ {1, 2, 4, 8}`, then MD5 every output file.

```
cpus=1  168174eb… 912abcc3… 76b31704… c235983a… bbf17701…   (.ibs .ibs0 .kin .kin0 allsegs.txt)
cpus=2  168174eb… 912abcc3… 76b31704… c235983a… bbf17701…
cpus=4  168174eb… 912abcc3… 76b31704… c235983a… bbf17701…
cpus=8  168174eb… 912abcc3… 76b31704… c235983a… bbf17701…
```

All five files are **byte-identical** across all four thread counts. The only stdout
differences are the echoed flag (`--cpus [1]` vs `--cpus [8]`), the line
`N CPU cores are used`, and how far the `0%1%2%…` progress counter gets before the work
finishes.

> **ESTABLISHED RULE.** `--cpus` changes no printed digit in any output file. Thread count
> is a pure performance knob; a single-threaded implementation is byte-compatible.
> Reproducing the *stdout* progress percentages exactly would require reproducing the
> thread schedule, so stdout must be normalised before diffing.

---

## Q5 — `--degree` semantics

### Experiment design

A dataset whose pairwise kinship values are *placed*, not simulated. For a pair over `M`
SNPs, choose the counts `a` (both het), `b = c` (one het/one hom), `d` (opposite hom),
`e` (same hom); with `Het_i = Het_j = H = a + b` the estimators reduce to

```
between-family  phi = 0.5 - (2d + b) / (2H)
within-family   phi = (a - 2d) / (2H)
```

so any rational `phi` on a `1/(2H)` grid can be dialled in exactly. With `M = 200 000`
(`H = 100 000`) the grid step is `5e-6`. 27 target kinships were emitted twice — once as
cross-family pairs (each member in its own FID) and once inside a shared FID — including
four probes straddling every inference boundary at ±5e-6 and ±2e-5. 108 samples.

### Raw results

```
                       .kin0 rows   .kin rows   header identical to no-degree?
no --degree               4347        1431      (baseline: 8 columns)
--degree 0                4347        1431      yes    option not echoed at all
--degree 1                1330        1431      yes
--degree 2                1779        1431      yes
--degree 3                2284        1431      yes
--degree 4                2772        1431      yes
--degree 5                3338        1431      yes
--degree 9                4091        1431      yes
--degree -1               4132        1431      yes
```

`.kin` is 1431 rows (= C(54,2) + header) in every single run.

Console line, `--kinship`:
`Between-family kinship data (up to degree 2, 1779 pairs in total) saved in file …`
Console line, `--related`:
`Between-family relatives (kinship >= 0.08839) saved in file …` — the printed cutoffs are
`0.17678, 0.08839, 0.04419, 0.02210, 0.01105, 0.00069` for degrees 1, 2, 3, 4, 5, 9,
i.e. `%.5lf` of `2^-(d+1.5)`.

Boundary membership of the designed probes:

| `--degree` | probe kinship | included? |
| --- | --- | --- |
| 1 | 0.176795 / 0.176780 | yes / yes |
| 1 | 0.176770 / 0.176755 | no / no |
| 2 | 0.088410 / 0.088395 | yes / yes |
| 2 | 0.088385 / 0.088370 | no / no |
| 3 | 0.044215 / 0.044200 | yes / yes |
| 3 | 0.044190 / 0.044175 | no / no |
| 4 | 0.022115 / 0.022100 | yes / yes |
| 4 | 0.022090 / 0.022075 | no / no |

A second, higher-resolution fileset (`M = 400 000`, grid step `2.5e-6`) separates the exact
power `2^-2.5 = 0.1767766953` from the *printed* literal `0.17678`:

| exact kinship | `--degree 1` includes it? | consistent with |
| --- | --- | --- |
| 0.1767725 | no  | both |
| 0.1767750 | no  | both |
| **0.1767775** | **yes** | **only `2^-2.5`** (the literal `0.17678` would exclude it) |
| 0.1767800 | yes | both |

`--ibs` line counts are unchanged by `--degree` (`.ibs` 1432 / `.ibs0` 4348 with and
without it).

> **ESTABLISHED RULE.**
> * `--degree d` filters **`.kin0` only**. `.kin` is never filtered — not its row set, not
>   its columns. Neither is `.ibs`/`.ibs0`.
> * The filter is on the **kinship estimate**, not on the inferred class:
>   a between-family pair is written iff `kinship >= 2^-(d + 1.5)`.
> * The threshold is the **exact IEEE double `2^-(d+1.5)`**, not the `%.5lf` value the
>   console prints. Comparison is `>=` (a pair at 0.1767775, above `2^-2.5` but below the
>   printed `0.17678`, is kept).
> * `--degree 0` is treated as *unset*: the flag is not even echoed in `Options in effect`,
>   and `--kinship` writes every between-family pair. Under `--related`, whose own default
>   is degree 1, `--degree 0` behaves as degree 1 (`kinship >= 0.17678`).
> * Degrees above the useful range work as the formula says: `--degree 20` prints
>   `kinship >= 0.00000` and keeps essentially everything.
> * **Negative degrees diverge between the two paths and are not fully characterised.**
>   `--related --degree -1` writes no `.kin0` at all and prints no cutoff line;
>   `--kinship --degree -1` writes 4132 of 4347 pairs with an effective cutoff
>   indistinguishable from `0` (min included kinship `0.0000`, max excluded `-0.0000`),
>   which is *not* `2^-0.5 = 0.7071`. Treat negative degrees as unspecified.

---

## Q6 — the sample-ID sort comparator

### Experiment design

Each probe is its own fileset: one *probe family* whose members are listed in the `.fam`
in a deliberately scrambled order, plus a two-member dummy family (`.kin` needs ≥ 2
families). The `.kin` block for the probe family reveals the comparator's output order.
44 probe families in total, sweeping widths, leading zeros, signs, decimal points,
exponent-looking strings, punctuation, case, and digit strings long enough to overflow 32-
and 64-bit integers.

### Raw results (input order → emitted order)

```
007 7 70                      -> 7 70 007
00 000 1 01 001               -> 1 00 01 000 001
1 01 10 0010 2 02             -> 1 2 01 02 10 0010
5 05 005 0005 00005           -> 5 05 005 0005 00005
9 10 11 100 99                -> 9 10 11 99 100
2147483647 2147483648 1 4294967295 4294967296
                              -> 1 2147483647 2147483648 4294967295 4294967296
9223372036854775806/07/08 18446744073709551616 5
                              -> 5 …806 …807 …808 18446744073709551616
99999999999999999999999999 1 …998   -> 1 …998 …999
1 2 10 20 100 0100 010 0001   -> 1 2 10 20 010 100 0001 0100
a1 a2 a01 a10 b1 1a           -> a1 a2 a01 a10 b1 1a
1 a 1a a1 10 z                -> a a1 z 1 1a 10
-1 1 -2 2 3                   -> -1 -2 1 2 3
+1 1 2 3                      -> +1 1 2 3
1.0 1 1.5 2 10                -> 1 1.0 1.5 2 10
1e3 1000 999 2000             -> 1e3 999 1000 2000
1-2 1-10 1_2 1.2 12           -> 1-2 1-10 1.2 1_2 12
0x10 16 10 0x9                -> 0x9 0x10 10 16
x1y2 x1y10 x2y1 x01y2 x1y02   -> x1y2 x1y02 x1y10 x2y1 x01y2
_1 1 -1 a                     -> -1 a _1 1
a A9 a10 A2                   -> a A2 A9 a10
A B Z aa bb                   -> A aa B bb Z
!1 #1 $1 %1 &1 A1             -> !1 #1 $1 %1 &1 A1
[1 A1 `1 ^1 {1 ~1 z1          -> A1 z1 [1 ^1 `1 {1 ~1
9z z9 a9 9a                   -> a9 z9 9a 9z
ab a1 a2 ac                   -> ab ac a1 a2
a1b a1 a1a a2                 -> a1 a1a a1b a2
aa1 a1a aab a1                -> aab aa1 a1 a1a
a- a1 ab a                    -> a a- ab a1
x x0 x00 x1                   -> x x0 x1 x00
```

Family (FID) blocks in `.kin` are ordered by the **same comparator** applied to the FID,
also independently of `.fam` order:

```
.fam FID order:   10 1 007 70 B a 2 7
.kin block order:  a B 1 2 7 10 70 007
```

### Fitting

* Digit runs cannot be compared character by character: `1 < 01` refutes that
  (`'0' < '1'`), so a **digit run is compared by length first, then lexicographically** —
  which is ordinary numeric comparison when there are no leading zeros, and puts
  zero-padded forms after their unpadded equivalents. This holds for digit strings of
  any length (26-digit strings order correctly), so **no integer parse is involved** and
  there is no overflow behaviour to emulate.
* Non-digit characters *are* compared one at a time, not run at a time: `ab < a1`
  (run-wise would give `a1 < ab`, since the run `"a"` is a prefix of `"ab"`).
* A digit sorts **after** any non-digit: `b1 < 1a`, `z9 < 9a`, `ab < a1`, `_1 < 1`.
  Plain byte order would put `'1'` (0x31) before `'b'` (0x62).
* Case folding is to **upper** case: `z1 < [1` (`'Z'` 0x5A < `'['` 0x5B). Lower-case
  folding would give `[1 < z1`. Confirmed twice (`H02`, `H03`) and by the punctuation
  ladder `!  #  $  %  &  A` which is plain byte order among non-digits.
* Case folding is also applied to **equality**: `{A, a}`, `{ab, aB}` and
  `{abc, ABC, Abc, aBc}` are each rejected at load time with
  `Family G14b: Person a is duplicated` / `FATAL ERROR - Please correct problems with
  pedigree structure`. IDs are therefore case-insensitively unique.
* A string that runs out first sorts first (`a < a1`, `x < x0`, `a1 < a1a`).

> **ESTABLISHED RULE — implementable directly.**
>
> ```rust
> use std::cmp::Ordering;
>
> /// Comparator used for `.kin` FID-block order and for within-family ID order.
> pub fn king_id_cmp(a: &[u8], b: &[u8]) -> Ordering {
>     let (mut i, mut j) = (0usize, 0usize);
>     loop {
>         match (i == a.len(), j == b.len()) {
>             (true, true)  => return Ordering::Equal,
>             (true, false) => return Ordering::Less,     // prefix sorts first
>             (false, true) => return Ordering::Greater,
>             _ => {}
>         }
>         let (ca, cb) = (a[i], b[j]);
>         let (da, db) = (ca.is_ascii_digit(), cb.is_ascii_digit());
>         if da != db {
>             // a non-digit run sorts BEFORE a digit run
>             return if da { Ordering::Greater } else { Ordering::Less };
>         }
>         if da {
>             // maximal digit runs: longer run is larger; equal length -> byte order
>             let ra = &a[i..i + a[i..].iter().take_while(|c| c.is_ascii_digit()).count()];
>             let rb = &b[j..j + b[j..].iter().take_while(|c| c.is_ascii_digit()).count()];
>             match ra.len().cmp(&rb.len()).then_with(|| ra.cmp(rb)) {
>                 Ordering::Equal => { i += ra.len(); j += rb.len(); }
>                 ord => return ord,
>             }
>         } else {
>             // single non-digit characters, ASCII-uppercase folded
>             match ca.to_ascii_uppercase().cmp(&cb.to_ascii_uppercase()) {
>                 Ordering::Equal => { i += 1; j += 1; }
>                 ord => return ord,
>             }
>         }
>     }
> }
> ```
>
> `.kin` emits family blocks in `king_id_cmp` order of FID, and within a block the
> `i < j` upper triangle over the members sorted by `king_id_cmp` of the IID — both
> independent of `.fam` order. This reproduces the `{007, 7, 70} → 7, 70, 007` puzzle
> from VERIFIED_FORMULAS: the digit runs have lengths 1, 2, 3.
>
> **Two hard input constraints discovered alongside it:** an IID of exactly `0` is
> rejected (`Parental sex codes don't make sense for Person 0 in Family …` →
> `FATAL ERROR`, because `0` is the missing-parent sentinel), and two IIDs in one family
> that differ only in case are rejected as duplicates.

---

## Q7 — output-file existence

Three states must be distinguished: **absent**, **zero-byte**, **header-only**.

### Experiment design

One 200-sample master (50 units of father/mother/two full sibs, one unit member replaced
by a duplicate of another sample, 20 000 SNPs over 22 autosomes, A1 the minor allele
throughout). Sample-count subsets `N ∈ {3 … 101}` were emitted in three family layouts —
**one** (all in one FID), **units** (families of 4), **singleton** (every sample its own
FID) — and each subset run with `--kinship`, `--related`, `--ibs`, `--duplicate`.

> A first attempt at this sweep produced spurious "no files" rows: the reference aborts
> with `FATAL ERROR - Too many first alleles as the major allele (~X%)` when too many
> `.bim` A1 alleles are the *major* allele. See [Side findings](#side-findings). Any
> experiment on this binary must use minor-allele-A1 filesets.

### Raw results (excerpt; `rows/bytes`)

```
ONE family                      | --kinship            | --related                          | --ibs                | --duplicate
  N=9    replaced with --kinship| .kin EMPTY  .kin0 -  | .kin EMPTY  .kin0 -                | .ibs 37/4928 .ibs0 - | .con 2/125
  N=25                          | .kin EMPTY  .kin0 -  | .kin EMPTY  .kin0 -                | .ibs 301/40086 …   - | .con 2/125
  N=50                          | .kin 1124/65558      | .kin 681/65593                     | .ibs 1226/163291     | .con 2/125
  N=100                         | .kin 4490/262263     | .kin 4762/459142  .kin0 -          | .ibs 4951/660086     | .con 2/125
families of 4
  N=9    replaced with --kinship| .kin 13/759 .kin0 25 | .kin 13/759  .kin0 25/1231         | .ibs 13   .ibs0 25   | .con 2/127
  N=10                          | .kin 14/818 .kin0 33 | .kin 14/1364 .kin0 -  "No close…"  | .ibs 14   .ibs0 33   | .con 2/127
  N=100                         | .kin 151   .kin0 4801| .kin 151/14664 .kin0 - "No close…" | .ibs 151  .ibs0 4801 | .con 2/127
all singleton families
  N=9    replaced with --kinship| .kin -    .kin0 37   | .kin -      .kin0 37/1964          | .ibs 1/139 .ibs0 37  | .con 2/131
  N=99                          | .kin -    .kin0 4852 | .kin -      .kin0 -   "No close…"  | .ibs 1/139 .ibs0 4852| .con 2/131
  N=100                         | .kin -    .kin0 4951 | .kin -      .kin0 124/11289        | .ibs 1/139 .ibs0 4951| .con 2/131
```

With the duplicate pair removed from the master, `--duplicate` gives `.con` = 1 row / 65 B
(header only) for `N ≤ 99` and **no file at all** for `N ≥ 100`.

### The single-family `.kin` is a truncated buffer, not an empty file

VERIFIED_FORMULAS records "one family ⇒ zero-byte `.kin`". That is the small-data face of
a buffer bug. `.kin` sizes for a one-family dataset:

```
N=50  1124 rows /  65 558 B      (complete would be 1226 rows)
N=80  2245 rows / 131 142 B
N=100 4490 rows / 262 263 B
N=200 19065 rows / 1 114 484 B   (complete would be 19 901 rows)
```

Every size is `k × ~65 545 B`. Padding the FID to change the line length pins it to bytes,
not rows — the same 50 samples, same genotypes, only the FID width changed:

| FID length | rows written | bytes | last line complete? |
| --- | --- | --- | --- |
| 1  | 1162 | 65 589 | yes |
| 3  | 1122 | 65 562 | yes |
| 8  | 1033 | 65 540 | yes |
| 13 |  958 | 65 566 | yes |
| 19 |  881 | 65 572 | yes |
| 33 |  742 | 65 611 | yes |

Every file ends on a line boundary at the *first* whole-line prefix that reaches 65 536
bytes. Re-labelling the same 200 samples as two families gives a **complete** 9901-row /
568 822 B `.kin`, so this affects the single-family case only. `.ibs` for the same
single-family datasets is complete (`N=100` → 4951 rows = C(100,2)+1).

> **ESTABLISHED RULES.**
>
> **`.kin` (`--kinship` and `--related`)**
> * Not created at all when no family has ≥ 2 members (all-singleton `.fam`).
> * Created when at least one family has ≥ 2 members.
> * **When the dataset contains exactly one distinct FID the file is truncated:** rows are
>   accumulated in a buffer that is written out and cleared each time it reaches 65 536
>   bytes, and the final partial buffer is never flushed. Content under 65 536 bytes
>   ⇒ zero-byte file. Model to implement: build the full text (header + rows), then emit
>   only the whole-line prefix produced by repeatedly cutting at the first line boundary
>   at or past each 65 536-byte mark.
> * With ≥ 2 distinct FIDs the file is complete regardless of size.
>
> **`.kin0` (`--kinship`)** — created iff at least one cross-family pair exists (≥ 2
> distinct FIDs). Never truncated.
>
> **`.kin0` (`--related`)** — created only when the between-family screening actually
> confirms relatives. That requires **N ≥ 100**: with 10 ≤ N ≤ 99 the reference prints
> `No close relatives are inferred.` and writes no `.kin0`, even when the data contains
> duplicates and parent–offspring pairs. At N ≥ 100 the screening runs
> (`Stages 1&2 (with … SNPs): … pairs of relatives are detected`) and the file appears.
>
> **`--related` sample-size gate** — with **N ≤ 9** the console prints
> `--related is replaced with --kinship for a small sample size.` and the run is
> byte-identical to `--kinship` (8-column `.kin0`, 10-column `.kin`). From N = 10 the real
> `--related` path runs, with the wider column sets of Q8.
>
> **`.ibs` / `.ibs0` (`--ibs`)** — no sample-size gate (they appear from N = 3).
> `.ibs` is **always** created; when no family has ≥ 2 members it is header-only
> (1 row / 139 B), never absent and never truncated. `.ibs0` is created iff ≥ 2 distinct
> FIDs exist.
>
> **`.con` (`--duplicate`)** — created for every N < 100 (header-only when there are no
> duplicate pairs). For N ≥ 100 it is created **only if at least one duplicate pair is
> found**, and is absent otherwise.
>
> **X-chromosome files** — if the `.bim` carries any X SNPs, `--kinship` also writes
> `<prefix>X.kin` / `<prefix>X.kin0` from a separate X analysis; `--related` writes
> `<prefix>X.kin`.
>
> **`<prefix>allsegs.txt`** — written by every run that performs IBD-segment analysis
> (see Q8), listing the usable segments.
>
> The integrator's observation that `trio` and `nuclear` produced *neither* `.ibs` nor
> `.ibs0` **did not reproduce**: run in a clean directory, `/tmp/kingcorpus/trio.bed`
> (3 samples, 1 family) yields `k.ibs` with 4 rows and `/tmp/kingcorpus/nuclear.bed`
> yields 16 rows; both correctly omit `.ibs0`. `threegen` (1 family, 14 samples) → `.ibs`
> only, `dups` (8 families) → both — matching the rules above.

---

## Q8 — `.ibs` / `.ibs0` column-set variation

### Hypothesis space

Extra trailing columns `MaxIBD2` / `Pr_IBD2` were seen on a whole-genome fileset but not
on a chr 1–2 subset of the same samples. Candidates: number of distinct chromosomes,
total map length, per-chromosome span, or SNP count / density sufficient for segment
analysis.

### Experiment design

A fixed 12-sample / 8-family pedigree (cross-family PO and FS with real recombination),
re-emitted over: 1–22 chromosomes at realistic spans; one chromosome at shrinking span;
one chromosome at falling SNP density; and a chromosome with a dense head and a sparse
tail.

### Raw results

```
set     "Total length … usable for IBD segment analysis"   .ibs0 columns
c01     1 seg / 248.9 Mb                                   21  (MaxIBD2, Pr_IBD2)
c22     14 seg / 1335.6 Mb                                 21
d01000  No informative IBD segments                        19
d02000  1 seg / 248.7 Mb                                   21
s025    1 seg / 62.2 Mb                                    19
s050    1 seg / 124.4 Mb                                   21
```

Span sweep at safe density:

```
usable length 97.9 / 98.9 / 99.9 Mb  -> 19 columns
usable length 100.0 / 100.9 / 104.9  -> 21 columns
```

The `100.0 Mb` case that gets the extra columns has an exact first-to-last-SNP span of
100 000 000 bp; the `99.9 Mb` case that does not has 99 990 000 bp.

Density sweep, uniform inter-SNP gap on one chromosome:

```
gap 156 072 / 156 200 / 156 249 / 156 250 bp  -> segment usable
gap 156 251 / 156 400 / 157 000 / 160 000 bp  -> "No informative IBD segments"
```

The genetic-map column is irrelevant: re-running the 157 kb-gap fileset with cM = 0 for
every SNP and with cM = 10 × Mb gives the same "no segments" result, and the 150 kb fileset
stays usable under all three cM encodings.

Segments are **sub-chromosomal**. A chr 1 with SNPs every 10 kb from 1–121 Mb and every
1 Mb thereafter yields

```
Segment Chr StartMB StopMB  Length  N_SNP StartSNP StopSNP
1       1   1.000   120.990 119.990 12000 rs1      rs12000
```

and the sparse tail is discarded. `<prefix>allsegs.txt` always lists exactly the segments
that were counted.

### What the columns contain

From the 22-chromosome run with real IBD segments (`.ibs0`):

```
FID1 ID1 FID2 ID2 … Kinship  MaxIBD2       Pr_IBD2
A0F  U0F A0C  U0C … 0.2469   0.000         0.0000     (parent-offspring)
A0C  U0C A0S  U0S … 0.2671   82495032.000  0.1592     (full siblings)
…                   0.0112   -9            -9         (unrelated)
```

`MaxIBD2` is the length of the **longest IBD2 segment in base pairs**, `%.3f`; `Pr_IBD2` is
the genome-wide **proportion in IBD2**, `%.4f`. Pairs that were not evaluated carry the
literal token `-9` in both columns (not `-9.000`, not `-9.0000`).

A designed dataset with pairs placed at exact kinships locates the gate:

```
kinship 0.0900 -> MaxIBD2 computed
kinship 0.0880 -> -9            (2^-3.5 = 0.0883883)
```

### `--related`'s `.kin0` varies the same way

```
no usable segments : FID1 ID1 FID2 ID2 N_SNP HetHet IBS0 HetConc HomIBS0 Kinship          (10)
usable segments    : … Kinship IBD1Seg IBD2Seg PropIBD InfType                            (14)
```

> **ESTABLISHED RULES.**
>
> **Usable IBD segments.** Walk each retained autosome in `.bim` order and cut it into
> maximal runs of consecutive SNPs whose base-pair gap is **≤ 156 250 bp** (= 0.15625 Mb =
> 10 Mb / 64; a gap of exactly 156 250 is usable, 156 251 is not). Each surviving run is a
> segment with `Length = last_bp − first_bp`. Only physical positions matter — the
> `.bim` cM column is ignored. `allsegs.txt` lists these; the console prints their count
> and the sum of their lengths as
> `Total length of %d chromosomal segments usable for IBD segment analysis is %.1lf Mb.`
> If no segment survives, the console prints `No informative IBD segments.` instead.
>
> **Column set.** `.ibs` (within-family) and `.ibs0` (between-family) gain the two extra
> trailing columns `MaxIBD2` and `Pr_IBD2` — **both files, always together** — iff the
> total usable segment length is **≥ 100 Mb (100 000 000 bp)**. The boundary is bracketed
> to `(99 990 000, 100 000 000]`. Below that they carry the short 19-/20-column header
> ending at `Kinship`. Neither the chromosome count nor the sample set has any effect.
>
> **Column values.** `MaxIBD2` (`%.3f`, base pairs) and `Pr_IBD2` (`%.4f`) are computed
> only for pairs whose `Kinship >= 2^-3.5` (2nd degree or closer, `0.0883883…`); every
> other pair gets the literal `-9` in both fields.
>
> **`--related`'s `.kin0`/`.kin`** use the same trigger: with usable segments they carry
> the `IBD1Seg / IBD2Seg / PropIBD / InfType` block, without them they stop at `Kinship`.
> `--kinship`'s `.kin0` is always the plain 8-column form regardless of the map.

### Exact headers observed

```
--kinship  .kin   FID  ID1 ID2       N_SNP Z0 Phi HetHet IBS0 Kinship Error
--kinship  .kin0  FID1 ID1 FID2 ID2  N_SNP HetHet IBS0 Kinship
--related  .kin   FID  ID1 ID2       N_SNP Z0 Phi HetHet IBS0 HetConc HomIBS0 Kinship
                  IBD1Seg IBD2Seg PropIBD InfType Error
--related  .kin0  FID1 ID1 FID2 ID2  N_SNP HetHet IBS0 HetConc HomIBS0 Kinship
                  [IBD1Seg IBD2Seg PropIBD InfType]
--ibs      .ibs   FID  ID1 ID2  Z0 Phi N_SNP N_IBS0 N_IBS1 N_IBS2 NHetHet NHomHom
                  N_Het1 N_Het2 IBS Dist HetConc Het2|1 Het1|2 HomConc Kinship
                  [MaxIBD2 Pr_IBD2]
--ibs      .ibs0  FID1 ID1 FID2 ID2  N_SNP N_IBS0 N_IBS1 N_IBS2 NHetHet NHomHom
                  N_Het1 N_Het2 IBS Dist HetConc Het2|1 Het1|2 HomConc Kinship
                  [MaxIBD2 Pr_IBD2]
--duplicate .con  FID1 ID1 FID2 ID2  N N_IBS0 N_IBS1 N_IBS2 Concord HomConc HetConc
allsegs.txt       Segment Chr StartMB StopMB Length N_SNP StartSNP StopSNP
```

Note `.ibs` carries `Z0 Phi` before `N_SNP` while `.ibs0` does not, and `InfType` values
observed are `Dup/MZ`, `PO`, `FS`, `2nd`, `3rd`, `4th`.

---

## Side findings

These fell out of the experiments above and matter for parity.

**1. The `.kin` `Error` column is not `%d`.** VERIFIED_FORMULAS lists `Error` as `%d`;
observed values include `0.5`. The column takes exactly three values, printed as `0`,
`0.5`, `1` (a `%G`-style format that drops trailing zeros). Reading them against the
inferred class:

```
pedigree unrelated, kinship 0.0161 -> inferred unrelated  -> 0
pedigree unrelated, kinship 0.0236 -> inferred 4th degree -> 0.5
pedigree unrelated, kinship 0.0427 -> inferred 4th degree -> 0.5
pedigree unrelated, kinship 0.0720 -> inferred 3rd degree -> 1
```

which is consistent with the tiny capture cited in VERIFIED_FORMULAS. `Error` is a graded
mismatch: `0` when the inferred class equals the pedigree class, `0.5` when it is off by
exactly one degree step, `1` when it is off by more.

**2. A1-major fatal check.** The reference refuses affected analyses when too many retained
autosomal markers encode A1 as the observed *major* allele:

```
FATAL ERROR -
Too many first alleles as the major allele (~11.4%). Please use plink1.9 --make-bed to
regenerate the genotype data again.
```

A controlled 20-sample, 5,000-marker sweep resolves the statistic: only the first **4,096**
retained autosomal markers enter the denominator, and a marker counts when its observed
A1/A1 count is strictly greater than its A2/A2 count. The threshold is strictly greater
than ten percent: 409/4,096 passes, 410/4,096 aborts and prints `~10.0%`; 904 A1-major
markers placed entirely after index 4,095 do not affect the gate.

The stable affected surface is `--related`, `--ibs`, `--unrelated`, `--build`, `--bysample`,
`--bySNP`, `--cluster` and `--ibdseg`. `--kinship`, `--duplicate` and `--autoQC` are exempt.
The `--related` and clustering gates start at ten samples; `--ibdseg` starts at five, matching
their full-pass thresholds. On shorter maps KING reads unstable tail state and can false-fatal;
open-king skips that unsafe case. The differential regression is
`tests/parity/probes/a1_major.py`.

**3. Unrecognised options.** Passing an unknown flag makes the reference print
`Please specify one of the following 24 options: --related --kinship …` and then proceed
with `--related`. (Discovered by accident: `zsh` does not word-split unquoted parameter
expansions, so a whole option string arrived as one argument.)

**4. Chromosome sort order.** If the `.bim` interleaves chromosomes (e.g. `rs22` on chr 22
followed by `rs23` on chr 1), the console prints
`Chromosomes unsorted: rs22 on chr 22, rs23 on chr 1.` and falls back to
`Relationship inference will be based on kinship estimation only.` — i.e. an unsorted map
disables IBD-segment analysis entirely and therefore changes the column sets of Q8.

---

## Q9 — `--related`'s between-family stage: which of three flows runs

### Experiment design

`--related` prints one of three tails after `Relationship inference across families starts
at`. The corpus alone cannot separate them: nine of its thirteen datasets are unrelated
across families, and only `bigish` has 100 samples. So the flows were separated by
constructing filesets — marker-count and sample-count ladders cut from `unrelated`,
`admixed` and `bigish`, plus synthetic 100-sample sets from
`tests/parity/probes/mkpairs.py` with one pair's IBD sharing swept — and reading the
reference's tail on each.

### The three tails

```
(K)  <41sp>ends at <ctime>
     No close relatives are inferred.
                                                 -- and no .kin0 is written

(S)    Stages 1&2 (with <s> SNPs): <d> pairs of relatives are detected (with kinship > <t>)
     <31sp>Screening ends at <ctime>
       Final Stage (with <m> SNPs): <c> pairs of relatives (up to <k>-degree) are confirmed
     <31sp>Inference ends at <ctime>

(E)  <31sp>Inference ends at <ctime>
       <c> pairs of relatives (up to <k>-degree) are identified
```

### Raw results

`--degree 0` behaves as `--degree 1` throughout.

| fileset | n | m | d1 | d2 | d3 | d4 |
| --- | ---: | ---: | --- | --- | --- | --- |
| `unrelated` | 30 | 20 000 | K | K | K | E |
| `unrelated` (m = 8 800) | 30 | 8 800 | K | K | **E** | E |
| `admixed` | 40 | 20 000 | K | K | E | E |
| `admixed` first 34 samples | 34 | 20 000 | K | K | **K** | – |
| `bigish` | 200 | 50 000 | **S** | **S** | E | E |
| `bigish` first 99 samples | 99 | 50 000 | K | – | – | – |
| `bigish` first 100 samples | 100 | 50 000 | **S** | – | – | – |

### Fitting

Three independent rules, each bracketed:

1. **The `A subset of informative SNPs will be used to screen close relatives.` /
   `Sorting autosomes...` header appears iff the effective degree is ≤ 2.** It is printed
   before any decision about what to do, so K and S both carry it and E never does.
2. **S needs 100 samples.** Bisected on `bigish` prefixes: 99 → K, 100 → S, with the same
   genotypes. The gate is unconditional, so below it `--related --degree 1` reports
   `No close relatives are inferred.` on a fileset holding a duplicate pair, an MZ pair and
   a parent–offspring pair across families — every one of which `--kinship --degree 1`
   lists in its `.kin0`. This is a reference bug and we reproduce it.
3. **E runs iff some between-family pair exceeds `2^-(degree + 2.5)`**, and otherwise K.
   Seventeen marker subsets of `unrelated` at degree 3 separate perfectly on the maximum
   between-family kinship: 0.0209 and below → K, 0.0228 and above → E, bracketing
   `2^-5.5 = 0.02210`. The candidate threshold is thus one degree looser than the
   `.kin0` reporting threshold `2^-(degree + 1.5)`; `admixed --degree 3` admits a 0.0254
   pair as a candidate, prints E, and then writes a **header-only** `.kin0` because 0.0254
   is under the 0.04419 reporting threshold.

### The screening counts

`Stages 1&2 (with <s> SNPs)` prints `s = min(m, 32768)`: 5 000 / 10 000 / 20 000 / 30 000
on `bigish` truncated to those maps, then 32 768 at both 40 000 and 50 000. `Final Stage`
always prints `m`.

Where `m <= 32768` there is no subsetting and `<d>` is exactly the number of
between-family pairs whose kinship exceeds `2^-(degree + 2)` — `bigish` truncated to 32 768
markers prints 18 at degree 1 and 50 at degree 2, and its own `.kin0` holds 18 and 50 pairs
over `2^-3` and `2^-4`. Above that the subset matters and **which 32 768 markers is
unresolved**: the reference prints 18 and 36 on the full 50 000-marker `bigish`, while the
first 32 768 markers give 18 and 50, the last give 22 and 46, evenly spaced give 21 and 50,
and the 32 768 highest-MAF give 23 and 45. No candidate reproduces both degrees.

The `Final Stage … <c> … confirmed` and `<c> … identified` counts are the same number: the
between-family summary table's total, which tallies `InfType` over the `.kin0` rows and
**never increments its own `4th` column** — `bigish --degree 4` writes 60 rows and reports
59.

---

## Q10 — `--ibs`'s `Pr_IBD2` is not `--ibdseg`'s `IBD2Seg`

Both are "the share of the usable genome called IBD2", both divide by the same `D`, and the
reference prints different numbers for the same pair in the same fileset. On `nuclear`:

| pair | `king.ibs` `Pr_IBD2` | `king.seg` `IBD2Seg` |
| --- | ---: | ---: |
| `N_C1`/`N_C2` | 0.2173 | 0.2626 |
| `N_C1`/`N_C3` | 0.2749 | 0.3144 |
| `N_C1`/`N_C4` | 0.4669 | 0.5095 |
| `N_C2`/`N_C3` | 0.4604 | 0.5194 |
| `N_C2`/`N_C4` | 0.2942 | 0.3531 |
| `N_C3`/`N_C4` | 0.2812 | 0.3108 |

`.ibs` is the smaller on every pair. So KING runs **two** IBD2 callers, and the tempting
cleanup — pointing `--ibs` at the `--ibdseg` engine — is wrong in a direction that is easy
to miss, because `MaxIBD2` (a maximum) tolerates the difference far better than `Pr_IBD2`
(a sum) does.
