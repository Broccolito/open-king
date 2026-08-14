# Attack 3 — the segment aggregates and the inference table

**Scope.** Everything *downstream* of the IBD-segment caller: `<prefix>allsegs.txt`,
the `InfType` decision table, the sixteen-column `.kin` / fourteen-column `.kin0` that
`--related` writes above its downgrade, and `<prefix>splitped.txt`.

**Clean-room statement.** KING's C++ source was not read. Every rule below is either read
off a golden capture in `tests/parity/golden/` or established by running the reference
binary `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` on synthetic
filesets built for the purpose. Every probe fileset is reproducible from the generators
described in §6.

**Headline results.**

1. `allsegs.txt` is a pure function of the marker map — byte-identical across all eight
   analyses that emit it, under every `--degree` / `--seglength` variant. Our `--ibdseg`
   path already needs no change.
2. `InfType` has a **two-clause FS test**, not the one-clause test in
   `docs/research/02-ibdseg.md`. The missing clause is `PropIBD > 0.32 && IBD2Seg > 0.15`.
   The corrected table reproduces **13 049 of 13 052** labels across the golden corpus and
   4 330 purpose-built probe rows; the three exceptions are rows whose printed 4-dp value
   sits exactly on a threshold.
3. The 16/14-column layouts are fully pinned, including `HomIBS0`, which was undocumented:
   `HomIBS0 = N_IBS0 / |{i is hom-A1} ∪ {j is hom-A1}|`.
4. `splitped.txt` is fully specified and **reproduced byte-for-byte on all ten corpus
   datasets plus four adversarial probe pedigrees**. A reference implementation is in
   §5.6.

---

## 1. `<prefix>allsegs.txt` — confirmed, no `--ibdseg`-specific behaviour

### 1.1 It is a function of the map alone

Over the whole golden corpus there are 163 `kingallsegs.txt` files. Grouped by dataset:

| dataset | files | distinct contents |
| --- | ---: | ---: |
| every one of the 13 | 3 … 18 | **1** |

That is: for a given `.bed`/`.bim`, *every* analysis that writes `allsegs.txt` writes the
same bytes. In particular

* `--seglength 5` and `--seglength 10` do **not** change it (the file describes the
  denominator, not the reported segments);
* `--degree` does not change it;
* `--bysample` / `--bySNP` / `--ibs` / `--ibdseg` / `--related` / `--unrelated` /
  `--build` / `--cluster` all agree.

So the rule already implemented in `king_core::ibdseg::usable_segments` — cut at each
chromosome change and each marker gap > 1 000 000 bp, then between complete 64-marker
words of the **global** grid whose 64-gap span exceeds 10 000 000 bp, keep a piece iff it
holds ≥ 5 complete words *and* its word-aligned span exceeds 10 000 000 bp — is the whole
story for `--ibdseg` too. **No change needed.**

### 1.2 Which analyses emit it

Emission is gated on whether the analysis reaches the segment machinery at all:

| analysis | emits `allsegs.txt` |
| --- | --- |
| `--bysample`, `--bySNP`, `--ibs` | always |
| `--ibdseg` | when `N ≥ 5` (below that it silently becomes `--kinship`) |
| `--related` | when `N ≥ 10` (below that it silently becomes `--kinship`) |
| `--unrelated`, `--build`, `--cluster` | only on `bigish` — i.e. only when the family-merging path is reached (`N ≥ 100`) |

`missing` and `nuclear` (N = 6) therefore have 8 `allsegs.txt` files each while `bigish`
has 17; that difference is entirely the `--related` downgrade.

### 1.3 One maintenance note (not a parity bug)

The usable-segment rule is implemented **twice**:

* `crates/king-core/src/ibdseg.rs::usable_segments` — used by `--ibdseg`
  (`analysis/ibdseg.rs:191`);
* `crates/king-cli/src/analysis/segments.rs::usable` / `usable_x` — used by `--ibs`
  and the QC passes (`analysis/ibs.rs:68`).

Both are byte-correct today. Since §1.1 proves there is only one rule, they should collapse
into one before `--related` adds a third caller.

---

## 2. `InfType` — the exact decision table

### 2.1 The rule

Let `π1 = IBD1Seg`, `π2 = IBD2Seg`, `π = PropIBD` (`= π2 + π1/2`, in f64, before the `%.4lf`
rounding). First match wins:

```
if   π2 > 0.7                                            ->  "Dup/MZ"
elif π1 + π2 > 0.96                                      ->  "PO"
elif π1 + π2 > 0.9   and π2 <  0.08                      ->  "PO"
elif π  > 0.3535534  and π2 >= 0.08                      ->  "FS"     <-- clause A
elif π  > 0.32       and π2 >  0.15                      ->  "FS"     <-- clause B  (NEW)
elif π  > 0.1767767                                      ->  "2nd"
elif π  > 0.08838835                                     ->  "3rd"
elif π  > 0.04419417                                     ->  "4th"
else                                                     ->  "UN"
```

`0.3535534 = 2^-1.5`, `0.1767767 = 2^-2.5`, `0.08838835 = 2^-3.5`, `0.04419417 = 2^-4.5`.
`0.32`, `0.15`, `0.08`, `0.7`, `0.9`, `0.96` are literal decimal constants, not powers of two.

The same table drives `.seg`, the `InfType` column of the 16-column `.kin` and of the
14-column `.kin0` — one implementation serves all three.

### 2.2 Clause B is real and was previously missed

`docs/research/02-ibdseg.md` §7 gives only clause A, transcribed from the R script KING
emits under `--rplot`. That R script is the **plotting** rule; the C++ writer has one more
FS clause. Evidence:

* Corpus counter-examples the one-clause rule gets wrong (all FS in the reference, all
  `2nd` under clause A alone):

  | dataset | π1 | π2 | π | clause B? |
  | --- | ---: | ---: | ---: | --- |
  | `missing` (`M_C1`/`M_C2`, default) | 0.3125 | 0.1888 | 0.3450 | yes |
  | `missing` (same pair, `--seglength 5/10`) | 0.3117 | 0.1805 | 0.3364 | yes |
  | `nuclear` (`N_C3`/`N_C4`, `--seglength 5/10`) | 0.1036 | 0.3012 | 0.3530 | yes |

* Corpus rows that clause B must **not** capture, which pins `π2 > 0.15` from below:
  `multifam` `(0.4310, 0.1310, 0.3465)` and `(0.4546, 0.1176, 0.3449)`, `bigish`
  `(0.4341, 0.1142, 0.3313)` — all `2nd`.

### 2.3 How the constants were bracketed

All brackets are from purpose-built probe filesets: 5 000 markers on a 250 Mb chromosome,
independent two-person families, each pair given a prescribed contiguous IBD2 block and
IBD1 block (§6). The realized `(π1, π2)` is read back from the reference's own `.seg`, so
the bracket is on the *printed* value and needs no modelling.

| boundary | tested on | last value below | first value above | constant |
| --- | --- | ---: | ---: | --- |
| `Dup/MZ` vs `FS` | π2, π1 ≈ 0 | 0.6925 | 0.7039 | **0.7** (π2 is word-quantised here, ±0.006) |
| `FS` → `PO` | π1+π2, π2 = 0.20 | 0.9599 | 0.9601 | **0.96** |
| `2nd` → `PO` | π1+π2, π2 = 0.04 | 0.8974 | 0.9014 | **0.90** |
| clause A's π2 gate | π2, π ≈ 0.42 | 0.0800 (`2nd`) | 0.0802 (`FS`) | **0.08** |
| clause A's π gate | π, π2 ∈ [0.10, 0.145] | 0.3529 | 0.3538 | **2^-1.5 = 0.3535534** |
| clause B's π2 gate | π2, π ≈ 0.33 | 0.1500 (`2nd`) | 0.1508 (`FS`) | **0.15** |
| clause B's π gate | π, π2 ∈ [0.16, 0.33] | 0.3199 | 0.3201 | **0.32** (printed 0.3200 appears with *both* labels) |
| `3rd` → `2nd` | π | 0.1768 (`3rd`) | 0.1773 | **2^-2.5** |
| `4th` → `3rd` | π | 0.0879 | 0.0889 | **2^-3.5** |
| `UN` → `4th` | π | 0.0437 | 0.0447 | **2^-4.5** |

`0.32` is nailed to the fourth decimal: a dense probe (400 pairs at π2 ≈ 0.313, π stepped
by ~2e-4 through 0.317…0.323) gives `2nd` for every printed π ≤ 0.3199, `FS` for every
printed π ≥ 0.3201, and **both labels at printed 0.3200** — exactly what a `> 0.32` test on
the unrounded double produces.

The π2 gates (0.7, 0.15, 0.08) are quantised by the probe's word granularity
(64 markers / 5 000 = 0.0128) except where the realized π2 happened to vary finely; the
0.15 and 0.08 brackets are tight (0.0008 and 0.0002 wide), 0.7 is ±0.006 and is
corroborated by the literal `data$IBD2Seg>0.7` in KING's own R output.

### 2.4 Inclusivity of each comparison

`>` on every PropIBD cut and on `π1+π2`; `>=` on clause A's `π2 >= 0.08`; `<` on the PO
clause's `π2 < 0.08`; `>` on clause B's `π2 > 0.15`. The `>=` in clause A is taken from
KING's R script — the probes cannot separate `>` from `>=` because a segment total that
lands *exactly* on the constant is unreachable (π2 is a ratio of base-pair sums to the
`allsegs` total). Any implementation choice is unobservable there; use the operators above.

### 2.5 Validation

```
13 052 rows   (4 172 golden .seg rows over all --degree/--seglength variants,
               4 248 golden 16-column .kin rows, 302 golden 14-column .kin0 rows,
               4 330 probe rows spanning the whole (π1, π2) simplex)
     3 mismatches, all on a printed threshold value:
       (0.0258, 0.3071, 0.3200) -> FS   [π internally just above 0.32]
       (0.6761, 0.0800, 0.4181) -> 2nd  [π2 internally just below 0.08]
       (0.2623, 0.0456, 0.1768) -> 3rd  [π internally just below 2^-2.5]
```

Nothing else in the (π1, π2) simplex misbehaves: a 231-point 0.05-grid, a 364-point
(π2, π) grid and several thousand random points all agree.

### 2.6 Traps

* **The `2nd` bucket has no upper bound.** `(π1, π2) = (0.8138, 0)` gives π = 0.4069 and
  is labelled `2nd` — π2 = 0 fails clause A and the sum 0.8138 fails both PO clauses.
* **PO can out-rank `Dup/MZ`'s neighbours.** `(0.4871, 0.5127)` → sum 0.9998 > 0.96 → `PO`,
  even though π2 = 0.51.
* **`InfType` follows `--seglength`.** Changing `--seglength` changes the reported π1/π2 and
  therefore the label: 5 `multifam` pairs and 1 `sexchr` pair flip between the default and
  `--seglength 10`. Do **not** compute the label from an unfiltered segment set.
* PO vs FS here is decided by IBD2 sharing, never by the IBS0 cutoff that
  `--build`/`--cluster`/`--unrelated` print.

---

## 3. The 16-column `.kin` and the 14-column `.kin0`

Both are written by `--related` once `N ≥ 10`. Tab separated, `\n`, trailing newline,
one header line.

### 3.1 `<prefix>.kin` — 16 columns

```
FID  ID1  ID2  N_SNP  Z0  Phi  HetHet  IBS0  HetConc  HomIBS0  Kinship  IBD1Seg  IBD2Seg  PropIBD  InfType  Error
```

| # | column | format | definition |
| ---: | --- | --- | --- |
| 1–3 | `FID` `ID1` `ID2` | `%s` | family, then the two IIDs |
| 4 | `N_SNP` | `%d` | SNPs non-missing in **both** (pairwise) |
| 5 | `Z0` | `%.3lf` | pedigree Pr[IBD=0] |
| 6 | `Phi` | `%.4lf` | pedigree kinship |
| 7 | `HetHet` | `%.4lf` | `NHetHet / N_SNP` |
| 8 | `IBS0` | `%.4lf` | `N_IBS0 / N_SNP` |
| 9 | `HetConc` | `%.4lf` | `NHetHet / (N_Het1 + N_Het2 − NHetHet)` |
| 10 | `HomIBS0` | `%.4lf` | `N_IBS0 / #{SNPs where i **or** j is homozygous A1}` |
| 11 | `Kinship` | `%.4lf` | within-family estimator `(HetHet − 2·IBS0)/(Het1 + Het2)` |
| 12–14 | `IBD1Seg` `IBD2Seg` `PropIBD` | `%.4lf` | segment estimates, `PropIBD = IBD2Seg + IBD1Seg/2` in f64 |
| 15 | `InfType` | `%s` | §2 |
| 16 | `Error` | `%G` | `0`, `0.5` or `1` — §3.3 |

**`HomIBS0` is new.** It is *not* `IBS0/HomHom` and *not* `1 − HomConc`. The denominator is
the **union** of the two samples' hom-A1 (first-allele-homozygote) sets over the pairwise
non-missing SNPs. Verified from the raw `.bed` on 727 rows (`bigish` and `multifam`,
`.kin` and `.kin0`) with zero mismatches; `bigish` `B01_C1`/`B01_C2` is `791/5830 = 0.1357`
where `HomHom = 25 472` and `IBS0/HomHom` would be `0.0311`.

`Z0`/`Phi` observed pairs: `(0.000, 0.2500)` PO, `(0.250, 0.2500)` FS,
`(0.500, 0.1250)` 2nd, `(0.750, 0.0625)` 3rd, `(1.000, 0.0000)` unrelated;
`(0.000, 0.5000)` would be MZ.

**Row inclusion:** *every* within-family pair, unfiltered. `--degree` never touches `.kin`
(`bigish` has 573 rows at every degree, identical to the `--kinship` `.kin`).

**Row order:** identical to the `--kinship` `.kin` — families in first-appearance order,
members sorted inside a family by the `docs/BEHAVIOR.md` §Q6 ID comparator, `i < j` upper
triangle. Verified equal, row for row, to `core/bigish__kinship/king.kin`.

**Existence:** written when some family has ≥ 2 members, and subject to the same
never-`fclose` truncation as `--kinship`'s `.kin` — `threegen` (12 samples, one family)
gets a **0-byte** `.kin` from `--related` exactly as it does from `--kinship`.

### 3.2 `<prefix>.kin0` — 14 columns

```
FID1  ID1  FID2  ID2  N_SNP  HetHet  IBS0  HetConc  HomIBS0  Kinship  IBD1Seg  IBD2Seg  PropIBD  InfType
```

Same formats as above. **There is no `Error` column, and no `Z0`/`Phi`.** `Kinship` is the
between-family (KING-robust) estimator `0.5 + (2·HetHet − 4·IBS0 − Het1 − Het2)/(4·min(Het1, Het2))`,
byte-identical to the value `--kinship` prints for the same pair. All six numeric columns
were re-derived from the raw `.bed` on `bigish --related --degree 4` (60 rows, 0 mismatches).

**Row inclusion** — this is the one rule that is *not* the obvious one. With
`d = --degree` (default 1), `kcut = 2^-(d+1.5)`, `pcut = 2^-(d+0.5)`:

```
a cross-family pair is written  iff  Kinship >= kcut   OR   PropIBD > pcut
```

i.e. a pair is kept when **either** estimator puts it within `d` degrees. Verified on all
17 `--related` cases in the corpus that write a 14-column `.kin0`: no expected row is ever
missing, and every "extra" row (2 on `bigish --degree 2`, 2 on `multifam --degree 3`, 2 on
`multifam --degree 4`) is explained by the `PropIBD` disjunct. Using `Kinship >= kcut`
alone loses 6 rows across the corpus.

The announcement line prints only the kinship half: `Between-family relatives (kinship >= %.5lf) saved in file %s`
(`0.17678`, `0.08839`, `0.04419`, `0.02210` for d = 1…4).

**Row order:** the block-tiled order of `docs/VERIFIED_FORMULAS.md` — sort the kept pairs
by `(i/32, j/32, i, j)` on `.fam` index. (`bigish`'s 60-row degree-4 file is consistent with
B ∈ {16, 32, 64}; B = 8 is excluded. B = 32 is already pinned by the `--kinship` `.kin0`.)

### 3.3 The `Error` column

`Error` compares the pedigree-declared class with an inferred class. `--related` and
`--kinship` compute it **differently** — they disagree on 9 rows of the corpus, so this
cannot be shared code.

*`--related`.* Inferred class `C_inf` is §2's table with **clause A's `π2 >= 0.08`
requirement dropped** (so any π > 2^-1.5 that is not PO/MZ counts as `FS`); pedigree class
`C_ped` comes from `(Z0, Phi)`. With `deg(MZ)=0, deg(PO)=deg(FS)=1, deg(2nd)=2, deg(3rd)=3,
deg(4th)=4, deg(UN)=5`:

```
Error = 0     if C_inf == C_ped
      = 0.5   if |deg(C_inf) − deg(C_ped)| == 1 and both degrees >= 2
      = 1     otherwise
```

The "both degrees ≥ 2" guard is what makes a declared-FS pair inferred `2nd` a **full**
error while a declared-2nd pair inferred `3rd` is a half error. Established with a
purpose-built 32-family probe (`perr`) that crosses four declared relationships
(FS / half-sib / PO / unrelated-within-family) with eight prescribed `(π1, π2)` points, and
consistent with all 4 248 golden `.kin` rows, including the two `monomorphic` rows that look
contradictory until the dropped π2 gate is applied (`P_C1`/`P_C4`, `InfType 2nd`,
π = 0.4487 → `C_inf = FS` → `Error 0`; `P_C3`/`P_C4`, `InfType 2nd`, π = 0.2406 → `Error 1`).

*`--kinship`.* Degree only, from the `Kinship` estimate, with no PO/FS distinction and no
"both ≥ 2" guard: 0 if the degrees match, 0.5 if adjacent, 1 otherwise.

*Not established:* the pedigree class for a declared 4th-degree pair (`Phi = 0.03125`) —
no corpus or probe pedigree reaches it.

### 3.4 The two `Relationship summary` tables on stdout

*Within-family* (printed whenever a `.kin` is written):

```
Relationship summary (total relatives: %d by pedigree, %d by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	…
  Inference	…
```

Both rows tally the `.kin` rows into six buckets; `OTHER` collects everything that is not
MZ/PO/FS/2nd/3rd (i.e. 4th and unrelated). The `Pedigree` row classifies by `(Z0, Phi)`, the
`Inference` row by `InfType`. `total relatives` is the sum of the **first five** buckets of
each row. Checked on `bigish` (`0 226 111 81 18 137` / `0 226 111 79 19 138`, totals
436 / 435) and `dups`.

*Between-family* (printed by the exhaustive path):

```
        	MZ	PO	FS	2nd	3rd	4th
  =========================================================
  Inference	…
```

One row, tallying the *written* `.kin0` rows by `InfType` — except that **the `4th` bucket
is never incremented**. `bigish --degree 4` writes a `4th` row and still prints `4th = 0`;
`multifam` prints the identical table at `--degree 3` and `--degree 4` even though the file
grows from 54 to 65 rows. The count in `  %d pairs of relatives (up to %d%s-degree) are identified`
is the sum of that table, so it can be **smaller than the number of rows in the file**
(59 vs 60 on `bigish --degree 4`), and `0 pairs … are identified` / `No cryptic relatedness
(up to the N-degree) is found.` is printed while a non-empty `.kin0` is written
(`admixed --degree 4`: 2 `UN` rows).

### 3.5 Open: which between-family path runs

`--related` has (at least) three between-family control flows, and the choice is **not**
simply `--degree ≥ 3`:

| console signature | seen on |
| --- | --- |
| `Stages 1&2 (with %d SNPs): … / Final Stage (with %d SNPs): …` | `bigish` (N=200) at d ≤ 2 |
| `A subset of informative SNPs will be used to screen close relatives. / Sorting autosomes…` then `No close relatives are inferred.` and **no `.kin0`** | `admixed`, `dups`, `monomorphic`, `multifam`, `sexchr`, `unrelated` at d ≤ 2 |
| `Inference ends at … / %d pairs of relatives (up to %d%s-degree) are identified` | everything at d ≥ 3 |

`unrelated` (N=30) is the exception that breaks the `d ≥ 3` reading: at `--degree 3` it
takes the screening path (no `.kin0` written at all) and only at `--degree 4` the exhaustive
one (header-only `.kin0`). `monomorphic` (N=12) takes the exhaustive path already at
`--degree 3` and writes a header-only `.kin0`. The screening path also *misses* relatives it
should find — `dups --related --degree 2` reports "No close relatives are inferred" on a
fileset whose cross-family duplicate pair `--degree 3` finds immediately. **This is the
`--related` screening algorithm and is a separate research target**; §3.2's inclusion rule
describes what lands in the file once the exhaustive path runs, and matches the
`Stages 1&2` path too on the one dataset (`bigish`) that reaches it.

---

## 4. `--related`'s silent downgrade

`N < 10` ⇒ `--related` runs `--kinship` instead (10-column `.kin`, 8-column `.kin0`, no
`allsegs.txt`, no segment columns). Corpus: `missing` and `nuclear` (6), `trio` (3),
`pair` (2), `singleton` (1). `--ibdseg`'s own threshold is `N < 5`.

---

## 5. `<prefix>splitped.txt`

### 5.1 When it is written

Only by `--ibdseg` (including `--related --ibdseg`), never by `--related` alone, `--ibs`,
`--bysample`, `--bySNP`, `--build`, `--cluster` or `--unrelated`. 50 of the corpus's 65
`ibdseg` cases have one — all except `trio`, `pair`, `singleton`, which downgrade to
`--kinship`. It is byte-identical across `--degree` and `--seglength`.

It is written **iff** at least one family survives §5.3's drop rule. A probe fileset of six
one-member families produces neither the file nor the console line
`<prefix>splitped.txt is generated for certain pedigree plot applications.` — the line and
the file always travel together.

### 5.2 Format

Nine **space**-separated fields per line, `\n`, trailing newline, **no header**:

```
OldFID  OldIID  NewFID  NewIID  Father  Mother  Sex  Pheno  Dummy
```

* `OldFID` is the `.fam` FID; `OldIID` == `NewIID` always (only the FID is ever rewritten).
* `Father`/`Mother` are the (possibly rewritten) parent IIDs, `0` when absent.
* `Sex` is the `.fam` sex verbatim.
* `Pheno` is the `.fam` phenotype with **`-9` mapped to `0`** (other values pass through:
  `nuclear`'s `1`/`2` and `bigish`'s `1`/`2` survive).
* `Dummy` is `1` for a person KING materialised, `0` for a genotyped `.fam` record.

### 5.3 Which families appear

A family is dropped iff it has exactly one member **and** that member declares no parent.
Everything else is written, including a one-member family whose member names an absent
parent (probe `q1`, family `ONE`: one genotyped sample plus two dummies).

### 5.4 The three ways a dummy is created

1. **Named parent absent from the whole fileset** — materialised under its declared name,
   `0 0` parents, sex from the slot (father ⇒ 1, mother ⇒ 2), `Dummy = 1`.
2. **Named parent that lives in a different family** — a *copy* is imported into the
   referencing family as a founder with `0 0` parents, the **source record's sex**, and
   `Dummy = 1`. (`multifam`: `FAM3`'s `C_F` declares `A_F`/`A_M`, who live in `FAM1`; both
   appear again under `FAM3` with `Dummy = 1`.)
3. **Exactly one parent slot filled** — the other parent is invented as `KING<n>` with the
   complementary sex and `Dummy = 1`, and the child's empty slot is rewritten to that name.
   `n` is a **global counter across families, incremented in `.fam` family order** (probe
   `q1`: `.fam` order MIX, ONE, TWO, THREE, FOUR gives ONE→`KING1`, TWO→`KING2`,`KING3`,
   THREE→`KING4`,`KING5`, while the *rows* come out in sorted-FID order).

### 5.5 Splitting, family order and row order

* Families are emitted in **sorted FID order** under the `docs/BEHAVIOR.md` §Q6 comparator,
  *not* `.fam` order. (Every corpus fileset happens to be already sorted; probe `q1`, whose
  `.fam` order is MIX, ONE, TWO, THREE, FOUR, emits FOUR, MIX, ONE, THREE, TWO.)
* Within a family, build the person list = dummies + `.fam` members and sort it by
  **(generation depth, ID comparator)**, where depth(founder) = 0 and
  depth(x) = 1 + max(depth of present parents).
* Split into connected components over parent–child edges.
  * **One component** → `NewFID = OldFID`, rows in the depth/ID order above.
  * **Two or more** → the k-th component (seeded by the first still-unassigned person in
    that order) gets `NewFID = <OldFID>_S<k>`, and its rows come out in **breadth-first
    traversal order from the seed**, neighbours visited in the depth/ID order.

  The BFS is only used on split families, and it is visible: probe family `BB` = two sibs
  `N_A`,`N_B` with absent parents `NP`,`NM` plus an unrelated founder `N_C` emits
  `NM, N_A, N_B, NP` for `BB_S1` (BFS from `NM`) while the same family without `N_C`
  (no split) emits `MM, MP, M_A, M_B` (plain depth/ID). Any implementation that uses one
  order for both cases will be wrong on one of them.

### 5.6 Verification

A ~120-line reconstruction of the above (scratchpad `splitped.py`) reproduces
**byte-for-byte**:

* all ten corpus `kingsplitped.txt` files (`admixed`, `bigish`, `dups`, `missing`,
  `monomorphic`, `multifam`, `nuclear`, `sexchr`, `threegen`, `unrelated`);
* four adversarial probe pedigrees run through the reference binary: mixed
  split/non-split families, one-member families with absent parents, cross-family parent
  references with a contradicting sex, a family of seven founders whose IDs pin the ID
  comparator (`QA, QM, QP, QZ, Q_, Q_A, Q0` — letters, then `_`, then digits), and families
  whose parents sort after their children.

---

## 6. Reproducing the probes

Everything in §2, §3.3 and §5 came from three generators, kept in `tests/parity/probes/`:

* `mkpairs.py OUT targets.txt` — writes a PLINK fileset: chromosome 1, 250 Mb,
  5 000 markers at 50 kb spacing, allele frequencies U(0.10, 0.50), one two-person family
  per target line `f1 f2 [offset]`. The pair shares a contiguous IBD2 block of `f2·M`
  markers and a contiguous IBD1 block of `f1·M` markers starting at `offset`; A1 is
  re-oriented per SNP to the observed minor allele so the
  `Too many first alleles as the major allele` gate never fires. Realized `(π1, π2)` is read
  back from the reference's `.seg`, so the prescribed values only need to be roughly right.
* `probe.py` — a thin `run(name, targets)` wrapper that generates, runs
  `king -b … --ibdseg --prefix name`, and returns the `(π1, π2, π, label)` rows for the
  self-pairs.
* `splitped.py FAM [GOLDEN]` — the §5.6 reconstruction / differ.

The `Error` probe (`perr`) overwrites the generated `.fam` with a hand-written pedigree
before running `--related --degree 4`; the `splitped` probes (`q1`–`q5`) do the same before
running `--ibdseg`.

---

## 7. What this changes for the implementation

| item | status |
| --- | --- |
| `allsegs.txt` for `--ibdseg` | already correct; nothing to do (consider de-duplicating the two implementations) |
| `InfType` | add clause B (`π > 0.32 && π2 > 0.15`); everything else in `docs/research/02-ibdseg.md` §7 stands |
| 16/14-column `.kin`/`.kin0` | fully specified here: columns, formats, `HomIBS0`, inclusion (`Kinship >= 2^-(d+1.5)` **or** `PropIBD > 2^-(d+0.5)`), order, `Error` |
| `splitped.txt` | fully specified and verified; independent of the segment caller, so it can land now |
| `--related` between-family control flow | **still open** — see §3.5 |
| pedigree class of a declared 4th-degree pair | **still open** — see §3.3 |

None of it unblocks the IBD1 run-acceptance rule (`docs/PARITY.md` §11.1): the `IBD1Seg`/
`IBD2Seg`/`PropIBD` values that all of the above consume are still mis-called. But every rule that
consumes them is now pinned, so once the caller is right the downstream is not a second
research problem.
