# PLINK 1 binary fileset (.bed/.bim/.fam) — implementation spec, as KING 2.3.2 consumes it

Recon note for the clean-room MIT reimplementation of KING's relatedness inference.
Compiled 2026-08-13.

**Reference binary:** `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(Mach-O 64-bit arm64, KING 2.3.2, `(c) 2010-2023 Wei-Min Chen`).
**PLINK binaries used for empirical verification:**
`runtimes/macos-arm64/plink1.9/plink` (v1.9.0-b.7.11, 19 Aug 2025) and
`runtimes/macos-arm64/plink2/plink2`.

**Provenance tags used throughout:**
- `[VERIFIED]` — reproduced empirically in this session (hexdump / KING run / cross-check).
- `[STRINGS]`  — read from the binary's embedded string constants (facts about output, not source).
- `[SPEC]`     — PLINK 1 format documentation / long-standing public format definition.

**Clean-room note:** nothing here is transcribed from KING's C++ source. All KING-specific
claims come from running the binary and from `strings -a` on it. The `.bed` format itself is
public PLINK documentation and is independently re-derived below from hexdumps.

**Scratch workspace with every fixture cited here:**
`…/74b7491e-2c75-4297-8170-f18e23fe8596/scratchpad/plinkio/`

---

## 0. Summary of what a reimplementation must get right

| # | Rule | Consequence if wrong |
|---|---|---|
| 1 | `.bed` bytes 0–1 are `0x6c 0x1b`; byte 2 is the mode | Silent garbage / wrong file accepted |
| 2 | Mode `0x01` = SNP-major (variant-major). KING **rejects** `0x00` | — |
| 3 | 2-bit codes: `00`=hom A1, `01`=**missing**, `10`=het, `11`=hom A2 | `01` is missing, NOT "one copy" — the single most common bug |
| 4 | **Low** bits of a byte hold the **first** sample | Whole matrix transposed within each byte |
| 5 | Each variant row is padded to a whole byte; padding bits are `0`, which decodes as `00` = hom A1 | Phantom homozygotes for up to 3 non-existent samples per variant |
| 6 | Read `.fam` first (gives `n`), `.bim` second (gives `m`), then `.bed` | Cannot compute `bytes_per_variant` |
| 7 | KING's autosome set = chr `1..22` **plus chr 25 (XY)** | ~PAR SNPs silently dropped or double-counted |
| 8 | All counts restricted to the pairwise-complete SNP set | Subtly wrong kinship for differential missingness |

---

## 1. `.bed` — the genotype matrix

### 1.1 Overall layout `[SPEC]` `[VERIFIED]`

```
offset  size                       content
------  -------------------------  --------------------------------------------
0       1                          0x6c        magic byte 0
1       1                          0x1b        magic byte 1
2       1                          0x01        mode: 1 = SNP-major, 0 = individual-major
3       m * bytes_per_variant      genotype payload, variant-major
```

with

```
bytes_per_variant = ceil(n_samples / 4) = (n_samples + 3) / 4     [integer division]
expected_file_len = 3 + n_variants * bytes_per_variant
```

`n_samples` comes from the `.fam` line count, `n_variants` from the `.bim` line count. The
`.bed` itself carries **no** counts — it is not self-describing, and the header is exactly
3 bytes with no length fields.

**Verified.** For the 6-sample × 4-variant fixture `t1`:

```
n=6 m=4 bytes_per_snp=2 expected_size=3+4*2=11 actual=11
```

### 1.2 The 2-bit genotype code `[SPEC]` `[VERIFIED]`

| bits | value | meaning | IBS/counting role |
|---|---|---|---|
| `00` | 0 | homozygous **A1/A1** (A1 = the allele in `.bim` column 5) | `hom_a1` |
| `01` | 1 | **missing** | excluded from every count |
| `10` | 2 | heterozygous A1/A2 | `het` |
| `11` | 3 | homozygous **A2/A2** (A2 = `.bim` column 6) | `hom_a2` |

Note the deliberately non-monotone assignment: the "dosage" order is `00 → 10 → 11` and `01`
is the odd one out. A naive `code as dosage` read is wrong.

### 1.3 Bit packing order within a byte `[SPEC]` `[VERIFIED]`

Four samples per byte. **Sample `i` occupies bits `2*(i mod 4)` and `2*(i mod 4)+1`, i.e. the
lowest bit pair holds the lowest-numbered sample.**

```
 bit:   7   6 | 5   4 | 3   2 | 1   0
        sample 3 | sample 2 | sample 1 | sample 0      (0-based, within the byte)
        ^ highest addresses                ^ FIRST sample in .fam order
```

Extraction:

```
byte_index = 3 + variant * bytes_per_variant + (sample >> 2)
code       = (bed[byte_index] >> (2 * (sample & 3))) & 3
```

### 1.4 Empirical proof of the bit order — hexdump vs decoded matrix `[VERIFIED]`

The fixture was **designed backwards from the bit-order hypothesis**: genotypes were chosen so
that, if and only if low bits hold the first sample, variant 1 must serialize to the byte
`0xE4` and variant 2 to `0x1B`. PLINK produced exactly those bytes.

Input `t1.ped` / `t1.map` → `plink --file t1 --make-bed --out t1`.

Full `.bed`, byte for byte:

```
$ xxd t1.bed
00000000: 6c1b 01e4 0f1b 0fff 0f55 05              l........U.
          ^^^^ ^^                                  magic 0x6c 0x1b, mode 0x01
                 ^^^^ ^^                           variant 0 (rsE4)      : e4 0f
                         ^^^^ ^^                   variant 1 (rs1B)      : 1b 0f
                                 ^^^^ ^^           variant 2 (rsMONO)    : ff 0f
                                         ^^^^ ^^   variant 3 (rsALLMISS) : 55 05
```

Decoded side by side (sample order = `.fam` order `I1..I6`):

```
SNP       A1 A2  raw bytes           I1        I2        I3        I4        I5        I6
rsE4      G  A   e4 0f           hom-A1   MISSING       het    hom-A2    hom-A2    hom-A2
rs1B      T  C   1b 0f           hom-A2       het   MISSING    hom-A1    hom-A2    hom-A2
rsMONO    0  A   ff 0f           hom-A2    hom-A2    hom-A2    hom-A2    hom-A2    hom-A2
rsALLMISS 0  0   55 05          MISSING   MISSING   MISSING   MISSING   MISSING   MISSING
```

Worked derivation for variant 0, byte 0 = `0xE4` = `0b11 10 01 00`:

```
bits 1-0 = 00 -> I1 hom A1 (GG, A1=G)   <- LOWEST bits = FIRST sample
bits 3-2 = 01 -> I2 missing
bits 5-4 = 10 -> I3 het    (AG)
bits 7-6 = 11 -> I4 hom A2 (AA, A2=A)
```
byte 1 = `0x0F` = `0b00 00 11 11`:
```
bits 1-0 = 11 -> I5 hom A2
bits 3-2 = 11 -> I6 hom A2
bits 7-4 = 0000 -> PADDING (samples 7,8 do not exist) — zero-filled
```

**Padding hazard, confirmed here:** the padding is `0`, and `00` is a *valid* code meaning
hom-A1. A reader that popcounts a whole word without masking will invent up to 3 phantom
hom-A1 samples per variant. `0x0F` (not `0x55`) is the direct evidence. Always mask the tail.

Mirror check on variant 1, byte `0x1B` = `0b00 01 10 11` → `I1`=hom A2, `I2`=het, `I3`=missing,
`I4`=hom A1 — the exact reverse pattern, as designed. Bit order is therefore unambiguous.

`rsMONO` (all samples `A A`) → `0xFF 0x0F`, i.e. every sample `11` = hom A2, and PLINK writes
`A1 = 0` in the `.bim` (no second allele observed). `rsALLMISS` → `0x55 0x05`, every sample
`01`, and **both** `.bim` allele columns are `0`.

### 1.5 Independent cross-check against PLINK's own decoder `[VERIFIED]`

To rule out a self-consistent misreading, the decoded matrix was compared against
`plink --bfile t1 --recode tab` (PLINK's own round-trip back to text), reconstructing the
expected 2-bit code from the `.ped` allele pair and the `.bim` A1/A2:

```
cross-check vs plink --recode: ALL MATCH
```

### 1.6 Mode byte and magic — KING's actual behavior `[VERIFIED]`

| Input | KING exit | KING message |
|---|---|---|
| mode byte `0x01` | 0 | normal load |
| mode byte `0x00` (individual-major) | 1 | `Currently only SNP-major mode can be analyzed.` |
| byte 0 corrupted (`0x00 0x1b 0x01`) | 1 | `Please use either PLINK or KING binary format as input.` |
| byte 1 corrupted (`0x6c 0x00 0x01`) | 1 | `Please use either PLINK or KING binary format as input.` |
| zero-length `.bed` | 1 | `Please use either PLINK or KING binary format as input.` |

All fatal messages are printed as:

```
\nFATAL ERROR - \n<message>\n
```

**Individual-major (`mode = 0x00`) is not supported by KING and need not be supported by us.**
It is a PLINK 1.07-era layout in which the matrix is transposed (one row per *sample*,
`ceil(m/4)` bytes each). Recommended behavior for parity: detect it and return a typed error
rather than attempting to transpose — matching KING means refusing the file. Modern PLINK
(1.9 / 2.0) never writes it.

### 1.7 Truncation / size mismatch `[VERIFIED]`

KING does **not** stat the file or validate its length up front. It streams
`bytes_per_variant` per variant and fails when a read comes up short:

```
FATAL ERROR -
Not enough genotypes at the <k>th marker
```

`<k>` is **0-based** (verified by bisection):

| `.bed` size | complete variants present | message |
|---|---|---|
| 3 bytes (header only) | 0 | `Not enough genotypes at the 0th marker` |
| 5 bytes | 1 | `Not enough genotypes at the 1th marker` |
| 7 bytes | 2 | `Not enough genotypes at the 2th marker` |

(`"0th"`/`"1th"`/`"2th"` — the suffix is a literal `th`, not English-correct ordinals.)

---

## 2. `.bim` — the variant map

### 2.1 Columns `[SPEC]` `[VERIFIED]`

Six columns, **one line per variant, in the same order as the `.bed` variant rows**. No header
line, no comment convention.

| # | Field | Type | Notes |
|---|---|---|---|
| 1 | Chromosome code | integer (or string with `--allow-extra-chr`) | see §2.3 |
| 2 | Variant ID | string | not required to be unique by KING |
| 3 | Position | double | centimorgans; `0` is the universal "unknown" |
| 4 | Base-pair coordinate | integer | 1-based |
| 5 | Allele 1 (**A1**) | string | PLINK 1.9 text-conversion convention: the **minor** allele; `0` if not observed |
| 6 | Allele 2 (**A2**) | string | the **major** allele |

**Delimiter:** PLINK 1.9 writes **TAB**-separated with a trailing `\n`. Verified byte-exactly:

```
$ xxd t1.bim
00000000: 3109 7273 4534 0930 0931 3030 3009 4709  1.rsE4.0.1000.G.
00000010: 410a ...                                  A.
```
i.e. `1<TAB>rsE4<TAB>0<TAB>1000<TAB>G<TAB>A<LF>`.

**KING's parser is whitespace-tolerant:** a space-delimited `.bim` loads identically
(verified). Treat the file as whitespace-separated tokens, not strictly TSV.

**KING requires all six columns.** A 5-column `.bim` is **not** an error — KING silently loads
**0 SNPs** (`PLINK maps loaded: 0 SNPs`) and continues. A reimplementation should reject a
short line loudly; this KING behavior is a silent-wrong-answer trap, not a feature to copy.

### 2.2 A1/A2 orientation `[VERIFIED]` `[STRINGS]`

- PLINK **1.9**, converting from text (`.ped`/`.map`), sets **A1 = minor allele** by allele
  count. Verified: in a 5000-variant fileset, the fraction of variants with `freq(A1) > 0.5`
  was **0/5000**.
- PLINK **2**, and PLINK 1.9 with `--keep-allele-order`, do **not** reorder; A1 tracks
  REF/ALT provenance instead and may well be the major allele.
- KING carries a warning about this `[STRINGS]`:

  > `Too many first alleles as the major allele (~%.1lf%%). Please use plink1.9 --make-bed to regenerate the genotype data again.`

  Emission could not be triggered in this session: a fileset with **98.74%** of variants having
  `freq(A1) > 0.5` produced no warning under `--kinship`, `--related`, `--ibdseg`, `--ibs`,
  `--pca`, `--bySNP`, `--roh`, `--autoQC`, `--duplicate`, `--unrelated`, `--bysample`, or
  `--makeGRM`. Treat it as reachable only from a path not exercised here.

**Why this does not affect us.** Every quantity in KING's robust estimator — `N_Aa^i`,
`N_Aa^j`, `N_HetHet`, `N_AA,aa` (IBS0), `M_ij` — is **invariant under swapping A1 with A2**:
heterozygosity is symmetric, and "opposite homozygotes" is symmetric. Swapping A1/A2 permutes
codes `00 ↔ 11`, which leaves all five counts unchanged. **Do not normalize allele order and
do not read allele letters at all for relatedness.** (Frequency-based estimators and the
`--risk`/GRM paths are a different matter; they are outside this document.)

### 2.3 Chromosome codes `[VERIFIED]`

PLINK normalizes human chromosome names to integers on `--make-bed`. Verified round-trip
(input `.map` used the string names, output `.bim` used the integers):

| Input name | `.bim` code |
|---|---|
| `1` … `22` | `1` … `22` |
| `X` | `23` |
| `Y` | `24` |
| `XY` (pseudo-autosomal region) | `25` |
| `MT` | `26` |
| unplaced / unknown | `0` |

PLINK also **sorts** the output `.bim` by (chromosome, bp); the observed order was
`0, 1, 22, 23, 24, 25, 26`. Do not rely on this — see §2.5.

### 2.4 KING's chromosome partition, and `--sexchr` `[VERIFIED]`

KING partitions loaded variants and reports the split on stderr/stdout. With a fixture of
40 SNPs each on chr 1 and 22, 40 on 23, and 20 each on 24, 25, 26, 0:

```
Genotype data consist of 100 autosome SNPs (including 20 XY SNPs), 40 X-chromosome SNPs,
  20 Y-chromosome SNPs, 20 mitochondrial SNPs
  20 other SNPs are removed.
PLINK maps loaded: 180 SNPs
```

**Therefore, with the default `--sexchr 23`:**

| Class | Chromosome codes |
|---|---|
| **Autosome** (what relatedness uses) | `1..22` **and `25` (XY / PAR)** |
| X | `23` |
| Y | `24` |
| MT | `26` |
| **Dropped** ("other SNPs are removed") | `0` and anything else |

> **The XY rule is the easy one to miss:** chr 25 is folded into the autosome pool and is
> counted as `100 autosome SNPs (including 20 XY SNPs)`. It is *not* treated as a sex
> chromosome. `40 (chr1) + 40 (chr22) + 20 (chr25) = 100`. ✓

**`--sexchr N` semantics, fully determined by sweeping `N`** (default `N = 23`):

```
autosomes = chromosomes 1 .. N-1,  plus chromosome N+2
X         = N
Y         = N+1
XY        = N+2      (folded into autosomes)
MT        = N+3
everything else (incl. 0) = removed
```

Verification table (same 200-variant fixture throughout):

| `--sexchr` | autosome | X | Y | MT | removed | extra banner |
|---|---|---|---|---|---|---|
| 23 (default) | 100 (incl. 20 XY) | 40 | 20 | 20 | 20 | — |
| 22 | 60 (incl. 20 XY) | 40 | 40 | 20 | 40 | `Non-human samples are analyzed, with 22 pairs of chromosomes` |
| 26 | 160 | 20 | 0 | 0 | 20 | `… with 26 pairs of chromosomes` |
| 30 | 180 | 0 | 0 | 0 | 20 | `… with 30 pairs of chromosomes` |

The `Non-human samples are analyzed, with %d pairs of chromosomes` banner prints whenever
`--sexchr != 23`. Note the report line elides zero-count classes (`--sexchr 26` prints only
`160 autosome SNPs, 20 X-chromosome SNPs`).

### 2.5 Sort order requirements `[VERIFIED]`

KING's tolerance depends on the analysis:

| Analysis | Unsorted `.bim` |
|---|---|
| `--kinship`, `--ibs` (order-independent counting) | **Accepted.** 800 variants with chromosomes interleaved `1,2,1,2,…` loaded and ran normally. Descending bp within a chromosome also accepted. |
| `--ibdseg` (position-dependent segment detection) | **Refused**, with a non-fatal early exit: `Chromosomes unsorted: u1 on chr 2, u2 on chr 1.` followed by `Note chromosomal positions can be sorted conveniently using other tools such as PLINK.` |

**Implication for us:** the KING-robust kinship/IBS work we are reimplementing is
order-independent, so we need not require or impose a sort. If IBD-segment work is added
later, a `(chr, bp)` ascending precondition becomes necessary.

---

## 3. `.fam` — the sample table

### 3.1 Columns `[SPEC]` `[VERIFIED]`

Six columns, **one line per sample, in the same order as the `.bed` sample bit positions.**
No header line.

| # | Field | Notes |
|---|---|---|
| 1 | Family ID (**FID**) | see §4.5 for the `0` convention |
| 2 | Within-family ID (**IID**) | `(FID, IID)` is the identity key |
| 3 | Paternal IID (**PAT**) | `0` = unknown/not in dataset |
| 4 | Maternal IID (**MAT**) | `0` = unknown/not in dataset |
| 5 | **SEX** | `1` = male, `2` = female, `0` (or any other value) = unknown |
| 6 | **PHENO** | `1` = control, `2` = case, `-9` / `0` / non-numeric = missing |

**Delimiter:** PLINK 1.9 writes **SPACE**-separated with a trailing `\n`. Verified:

```
$ xxd t1.fam
00000000: 4631 2049 3120 3020 3020 3120 310a ...   F1 I1 0 0 1 1.
```
i.e. `F1<SP>I1<SP>0<SP>0<SP>1<SP>1<LF>`.

Note the asymmetry with `.bim`: **`.bim` is TAB-delimited, `.fam` is SPACE-delimited.** Both
KING and any sane reader should treat both as whitespace-separated. Verified: a TAB-delimited
`.fam` loads identically in KING.

PLINK wrote `0` for the sample declared with an unrecognized sex and `-9` for missing
phenotype (`6 people (3 males, 2 females, 1 ambiguous) loaded from .fam`).

**KING requires all six columns.** A 5-column `.fam` is **not** an error — KING silently loads
**0 samples** (`PLINK pedigrees loaded: 0 samples`). Again: reject loudly in our
implementation rather than mimic this.

### 3.2 Sample count is authoritative

`n_samples = number of lines in .fam`. This single number determines `bytes_per_variant` for
the whole `.bed`. There is no cross-check available inside the `.bed` — see §4.3.

---

## 4. Edge cases

### 4.1 Individual-major `.bed` — §1.6

Mode `0x00`. KING: `FATAL ERROR - Currently only SNP-major mode can be analyzed.` Recommended:
same refusal.

### 4.2 Missing / zero-length files `[VERIFIED]`

| Situation | KING behavior |
|---|---|
| `.bed` given without extension (`-b prefix`) | `FATAL ERROR - Genotype file prefix cannot be opened` — **KING requires the literal `.bed` filename**, it does not append the extension |
| `.bim` absent | `FATAL ERROR - Map file <prefix>.bim cannot be opened` |
| `.fam` absent | `FATAL ERROR - Pedigree file <prefix>.fam cannot be opened` |
| `.bed` zero-length | `FATAL ERROR - Please use either PLINK or KING binary format as input.` (fails the magic check) |
| `.bed` truncated mid-matrix | `FATAL ERROR - Not enough genotypes at the <k>th marker` (`k` 0-based) |

KING's load order is strictly **`.fam` → `.bim` → `.bed`**, with progress lines:

```
Loading genotype data in PLINK binary format...
Read in PLINK fam file <prefix>.fam...
  PLINK pedigrees loaded: <n> samples
Read in PLINK bim file <prefix>.bim...
  Genotype data consist of <a> autosome SNPs[, ...]
  PLINK maps loaded: <m> SNPs
Read in PLINK bed file <prefix>.bed...
  PLINK binary genotypes loaded.
  KING format genotype data successfully converted.
Autosome genotypes stored in <w> words for each of <n> individuals.
```

`--fam <path>` and `--bim <path>` override the derived sidecar paths (verified working).

### 4.3 Mismatched counts — **silent corruption in KING** `[VERIFIED]`

This is the most dangerous class, because KING has no length validation:

| Situation | KING behavior |
|---|---|
| `.fam` has **fewer** samples than the `.bed` was written for | **exit 0, no warning.** KING recomputes `bytes_per_variant` from the short `.fam`; trailing samples are dropped, and if `ceil(n'/4) != ceil(n/4)` every variant row after the first is **misaligned** and silently decodes to garbage. |
| `.bim` has **fewer** variants than the `.bed` contains | **exit 0, no warning.** Trailing variants ignored. |
| `.bim` has **more** variants than the `.bed` contains | `FATAL ERROR - Not enough genotypes at the <m_bed>th marker` |

**Recommendation for our implementation (deliberate divergence from KING):** validate up front

```
expected = 3 + n_variants * ceil(n_samples / 4)
if actual_len != expected  ->  hard error naming all three counts
```

An exact equality check catches every case above, costs one `stat`, and converts a silent
wrong answer into a diagnosable one. This is a *safety* divergence, not a numeric one: on any
well-formed fileset the results are identical.

### 4.4 Duplicate IDs `[VERIFIED]`

The identity key is the **pair** `(FID, IID)`, not IID alone.

| Situation | KING behavior |
|---|---|
| Same IID under **different** FIDs (30 samples all named `SAME`, FIDs `F1..F30`) | **Accepted**, loads 30 samples |
| Same `(FID, IID)` twice (all rows forced to `F1 I1`) | `Family F1: Person I1 is duplicated` then `FATAL ERROR - Please correct problems with pedigree structure` |

KING also carries `[STRINGS]` an ID-uniqueness note used by the pedigree-clustering path:
`Individual IDs are not unique (e.g., %s), and family IDs will be used as well.` and
`Individual IDs are unique across all families.`

### 4.5 The `"0"` family-ID convention `[VERIFIED]` `[STRINGS]`

FID `0` is the conventional "no known family" marker, and PLINK 2's `--dummy` emits it by
default (`0 per0 0 0 2 1`). KING treats it specially:

```
All individuals with family ID 0 are considered as relatives.
```

i.e. FID `0` does **not** mean "each sample is its own family" — every FID-`0` sample lands in
one pool treated as within-family for the relatedness pass. Relevant because a `.gq` bundle
whose `.fam` uses `0` throughout will be analyzed as a single family.

`PAT`/`MAT` values of `0` mean "parent unknown / not genotyped" and must not be looked up as a
sample named `"0"`.

### 4.6 Padding bits — §1.4

Zero-filled, decode as `00` (hom A1). **Must be masked off.** For `n_samples % 4 != 0` the
final byte of every variant row carries `4 - (n_samples % 4)` phantom samples.

### 4.7 Monomorphic and all-missing variants `[VERIFIED]`

- A variant monomorphic in the fileset gets `.bim` A1 = `0` (e.g. `1 rsMONO 0 3000 0 A`); all
  samples encode as `11`.
- An all-missing variant gets **both** allele columns `0` (`1 rsALLMISS 0 4000 0 0`); all
  samples encode as `01`.
- **KING keeps monomorphic variants in its IBS/kinship counts.** Verified: a 4000-variant
  fileset containing 506 monomorphic variants reproduced KING's `--ibs` counts exactly when
  our reference counter used all 4000 (§6.4). Do **not** filter monomorphic variants.

### 4.8 Sample-level SNP-count screen — observed, not fully explained `[VERIFIED]`

KING excludes samples from the kinship analysis with:

```
The following %d samples are excluded from the kinship analysis (M<%d):
```

The printed threshold was always `512`. The exclusion boundary was located precisely and is
**a function of the total autosomal variant count only**:

- excluded at `m = 544`, not excluded at `m = 545` — a sharp boundary;
- **independent of sample count** `n` (identical at `n = 4, 8, 12, 40`);
- **independent of MAF / monomorphic content** — a fixture of 640 variants that were *all*
  monomorphic passed, while 544 fully polymorphic variants failed. So `M` is not "polymorphic
  variant count" (an initially plausible reading that a partial sweep appeared to support and
  a controlled fixture then falsified).

The mapping from `m` to the compared `M` is therefore **not** `M = m` (544 ≥ 512 yet fails) and
not the polymorphic count. Flagging as **open**: it belongs to the estimator/QC document, not
to PLINK I/O, and it does not affect file parsing. It matters only for tiny filesets — every
realistic `.gq` bundle is far past this boundary. Do not encode a guessed rule.

---

## 5. Reference reader — recommended Rust

```rust
pub const BED_MAGIC: [u8; 2] = [0x6c, 0x1b];
pub const MODE_SNP_MAJOR: u8 = 0x01;
pub const MODE_INDIVIDUAL_MAJOR: u8 = 0x00;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Geno { HomA1, Missing, Het, HomA2 }

#[inline]
pub fn decode(code: u8) -> Geno {
    match code & 3 {
        0 => Geno::HomA1,
        1 => Geno::Missing,
        2 => Geno::Het,
        _ => Geno::HomA2,
    }
}

#[inline]
pub fn bytes_per_variant(n_samples: usize) -> usize { (n_samples + 3) / 4 }

/// Byte offset of `sample`'s code within a variant row, and the shift to apply.
#[inline]
pub fn locate(sample: usize) -> (usize, u32) { (sample >> 2, 2 * (sample as u32 & 3)) }
```

Validation order (mirrors KING, plus the §4.3 length check):

1. `.fam` → `n_samples`; error on any line with < 6 whitespace fields.
2. `.bim` → `n_variants`, per-variant `(chrom, id, cm, bp, a1, a2)`; error on < 6 fields.
3. `.bed`: require `len >= 3`, `buf[0..2] == BED_MAGIC`.
4. `buf[2]`: `0x01` ok; `0x00` → `Err(IndividualMajorUnsupported)`; else `Err(BadMode)`.
5. `len == 3 + n_variants * bytes_per_variant(n_samples)` → else `Err(LengthMismatch{..})`.
6. Classify chromosomes per §2.4; keep the autosome set (`1..=22` and `25` at default
   `sexchr = 23`).

`mmap` the `.bed` — it is read-only, read-once, and often multi-GB; a bundle's genome-wide
fileset should never be copied into a `Vec<u8>`.

---

## 6. The performance core: bit-plane representation and popcount expressions

### 6.1 Why transpose

KING's inner loop is pairwise over samples and sequential over variants. The `.bed` layout is
variant-major, so a pair `(i, j)` would touch two bytes per variant with shifts and masks —
`O(m)` branchy work per pair, `O(n² m)` overall. Transposing once into per-sample bitvectors
turns the inner loop into whole-word `AND`/`popcount`, processing **64 variants per
instruction**.

**KING does exactly this, and uses 64-bit words** `[VERIFIED]`. The load banner
`Autosome genotypes stored in %d words for each of %d individuals` reports the per-sample word
count; sweeping the variant count pins the word width:

| autosomal SNPs | words reported |
|---|---|
| 16, 17, 32, 33, **64** | 1 |
| **65**, 128 | 2 |
| 100 | 2 |

The 64 → 1, 65 → 2 step is conclusive: `words = ceil(m / 64)`, i.e. **64-bit words**.

### 6.2 Recommended layout

Per sample, store `W = ceil(m_autosome / 64)` words in **each of four planes**. Bit `b` of
word `w` corresponds to autosomal variant index `64*w + b`.

| Plane | Set when | Storage |
|---|---|---|
| `hom_a1[i]` | code `00` | `W` u64 |
| `het[i]`    | code `10` | `W` u64 |
| `hom_a2[i]` | code `11` | `W` u64 |
| `nonmiss[i]`| code `!= 01`, i.e. `hom_a1 \| het \| hom_a2` | `W` u64 |

`nonmiss` is derivable but is precomputed once at load: it is used in four of the five
expressions below, so materializing it removes two `OR`s from every word of every pair.

Cost: 4 bits/variant/sample = 4 MB per sample per 8M variants. For a `.gq` bundle (a handful
of genomes) this is trivial; even at biobank scale it is `n * m / 2` bytes.

**Tail masking (essential, see §4.6 and §1.4).** If `m_autosome % 64 != 0`, the high bits of
the last word must be zero in every plane. Build the planes by setting bits only for real
variants (never by copying `.bed` bytes wholesale) and the invariant holds by construction;
assert it once after load:

```rust
debug_assert_eq!(nonmiss[i][W-1] & !tail_mask, 0);
// tail_mask = if m % 64 == 0 { !0u64 } else { (1u64 << (m % 64)) - 1 }
```

With the tail zeroed in `nonmiss`, every expression below is automatically tail-safe, because
each one is `AND`ed with at least one `nonmiss`/`het`/`hom_*` plane.

### 6.3 The five counts — exact popcount expressions

All counts are taken over the **pairwise-complete** variant set `M_ij` (both samples
non-missing) — this is the parity-critical rule from `01-paper-estimators.md` §1.

```rust
let (mut m_ij, mut n_aa_i, mut n_aa_j, mut n_hethet, mut n_ibs0) = (0u64, 0u64, 0u64, 0u64, 0u64);
let mut n_homhom = 0u64; // only if KING's NHomHom column is needed

for w in 0..W {
    let (nm_i, nm_j)   = (nonmiss[i][w], nonmiss[j][w]);
    let (he_i, he_j)   = (het[i][w],     het[j][w]);
    let (a1_i, a1_j)   = (hom_a1[i][w],  hom_a1[j][w]);
    let (a2_i, a2_j)   = (hom_a2[i][w],  hom_a2[j][w]);

    // N_shared_nonmissing  ==  M_ij
    m_ij     += (nm_i & nm_j).count_ones() as u64;

    // N_Aa^i : i heterozygous, restricted to sites where j is also non-missing.
    //          (het_i already implies nm_i, so no `& nm_i` is needed.)
    n_aa_i   += (he_i & nm_j).count_ones() as u64;

    // N_Aa^j : symmetric
    n_aa_j   += (he_j & nm_i).count_ones() as u64;

    // N_hethet == N_Aa,Aa : both heterozygous (both non-missing implied)
    n_hethet += (he_i & he_j).count_ones() as u64;

    // N_IBS0 == N_AA,aa : opposite homozygotes (both non-missing implied)
    n_ibs0   += ((a1_i & a2_j) | (a2_i & a1_j)).count_ones() as u64;

    // optional: both homozygous (any combination)
    n_homhom += (((a1_i | a2_i) & (a1_j | a2_j))).count_ones() as u64;
}
```

Derived, no extra passes:

```rust
let n_ibs1 = n_aa_i + n_aa_j - 2 * n_hethet;   // exactly one of the pair is het
let n_ibs2 = m_ij - n_ibs0 - n_ibs1;
```

**Cost:** 5 `popcnt` + 7 bitwise ops per 64 variants per pair. On arm64 `count_ones()` lowers
to `cnt`+`addv` over a NEON register; on x86-64 with `popcnt` it is a single instruction. The
five accumulators are independent, so they pipeline.

**Note on the `& nm_j` in `n_aa_i`:** this is *the* subtle correctness point. Dropping it
yields each sample's marginal heterozygote count instead of the pairwise-restricted one, which
is wrong for any pair with differential missingness. It is the exact failure mode called out
in `01-paper-estimators.md` §1.

**Alternative 3-plane layout** (if memory is the binding constraint): store `nonmiss`, `hom`
(= `hom_a1 | hom_a2`), and `a2` (= `hom_a2`). Then `het_i = nm_i & !hom_i` and
`IBS0 = hom_i & hom_j & (a2_i ^ a2_j)` — 25% less memory, ~2 extra ops per word. Not worth it
at bundle scale; prefer the 4-plane version.

### 6.4 Empirical validation against KING `[VERIFIED]`

The expressions above were implemented against a fixture with real missingness
(`plink2 --dummy 8 4000 0.08`, 8 samples × 4000 variants, ~8% missing, 506 of them
monomorphic) and compared to `king -b ibs1.bed --ibs`, which reports these counts directly.

KING's `.ibs` header (tab-separated) is:

```
FID	ID1	ID2	Z0	Phi	N_SNP	N_IBS0	N_IBS1	N_IBS2	NHetHet	NHomHom	N_Het1	N_Het2	IBS	Dist	HetConc	Het2|1	Het1|2	HomConc	Kinship
```

Mapping to our symbols: `N_SNP` = `M_ij`, `N_Het1` = `N_Aa^i`, `N_Het2` = `N_Aa^j`,
`NHetHet` = `N_HetHet`, `N_IBS0` = `N_AA,aa`.

Result over all C(8,2) = 28 pairs, comparing all eight integer columns
(`N_SNP, N_IBS0, N_IBS1, N_IBS2, NHetHet, NHomHom, N_Het1, N_Het2`):

```
pairs compared: 28   mismatches: 0
sample row: ('per0','per1') -> M_ij=3378 IBS0=265 IBS1=1398 IBS2=1715
                               HetHet=423 HomHom=1557 N_Het1=1113 N_Het2=1131
```

This simultaneously validates (a) the `.bed` decoding of §1, (b) the bit-plane expressions of
§6.3, (c) the pairwise-complete restriction on `N_Aa^i`/`N_Aa^j`, and (d) that monomorphic
variants are retained (§4.7). The derived identities `N_IBS1 = N_Aa^i + N_Aa^j - 2·N_HetHet`
and `N_IBS2 = M_ij - N_IBS0 - N_IBS1` reproduced KING's own columns exactly, so they are
confirmed as KING's definitions and not merely ours.

---

## 7. Fixtures produced (all in `scratchpad/plinkio/`)

| Fixture | Purpose |
|---|---|
| `t1.{bed,bim,fam,ped,map}` | 6×4 designed bit-order proof (`0xE4`/`0x1B` signature), monomorphic + all-missing variants |
| `t1_im.bed` | mode byte forced to `0x00` (individual-major rejection) |
| `bad.bed`, `m2.bed`, `t0.bed` | corrupted magic byte 0 / byte 1 / zero-length |
| `t3.bed`, `t5.bed`, `trunc.bed` | truncation at marker 0 / 1 / 2 |
| `f5.*`, `b3.*`, `b5.*` | short `.fam`, short `.bim`, long `.bim` |
| `f5c.*`, `b5c.*` | 5-column `.fam` / 5-column `.bim` |
| `c1.{bed,bim,fam}` | chromosome-code round-trip (`X/Y/XY/MT/0` → `23/24/25/26/0`) |
| `m1.*` | 200 variants across chr 1,22,23,24,25,26,0 — `--sexchr` sweep |
| `ibs1.*` + `ibs1k.ibs` | 8×4000 with 8% missingness — the §6.4 popcount validation |
| `uns.*`, `desc.*` | unsorted `.bim` (chromosome-interleaved / descending bp) |
| `s544/s545/n{4,8,12,40}m{544,545}` | §4.8 threshold bisection |

---

## 8. Open items

1. **§4.8** — the `M<512` sample screen: boundary empirically at `m = 545` autosomal variants,
   independent of `n`, MAF, and monomorphic fraction, but the `m → M` mapping is unexplained.
   Belongs to the estimator/QC doc. Does not affect parsing.
2. **§2.2** — the `Too many first alleles as the major allele` warning could not be triggered
   through twelve analysis options. Its emission path is unknown. Irrelevant to relatedness
   (all five counts are allele-order invariant), but would matter if frequency-based
   estimators or the GRM/risk paths are ever implemented.
