# KING 2.3.2 — Command-Line Surface & Console Output (recon 06)

**Reference binary:** `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
Mach-O 64-bit executable **arm64**, 1,815,336 bytes. Version banner: `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`.

**Method:** black-box execution of the binary (≈1,200 invocations) plus extraction of
`__TEXT,__cstring` literals with `strings`/`otool`. **No source code was read or
transcribed.** Every rule below was derived from observed input→output behavior; the
string dump was used only to enumerate literal constants (facts about output format).

**Validation:** a clean-room model of the parameters-in-effect block was written from
these rules and diffed against the binary over **600 randomized flag combinations**
(random option subsets, random command-line ordering, long/short values, all numeric
magnitudes): **600/600 exact byte matches.**
Model: `…/scratchpad/probe/cli06/model.py`
Captured transcripts: `…/scratchpad/probe/cli/cli06-*.stdout|.stderr|.exit`

---

## 1. Top-level facts

| Property | Value |
|---|---|
| All console output goes to | **stdout only** — stderr is empty (0 bytes) in every observed case, including FATAL ERROR |
| Exit code, argument/input error | **1** |
| Exit code, successful analysis | **0** |
| Exit code, unknown option + valid input | **0** (unknown options do **not** abort the run) |
| Buffering | plain `printf`; progress uses `\r` |
| No `--help`, `--version`, `-h`, `-v` | all rejected as undefined/ignored |
| Locale | C numeric formatting (`.` decimal point) |

Skeleton of every invocation:

```
<banner>
<blank>
The following parameters are in effect:
<Binary File line>
<blank>
Additional Options
<group lines…>
<blank>                                  ← always emitted (options-block trailer)
[ WARNING block ]                        ← only if there were parse problems
[ "KING starts at <ctime>"  … analysis ] ← only if a genotype file was supplied
[ FATAL ERROR block ]                    ← on fatal error
```

---

## 2. Verbatim console output — no arguments

`king` (exit **1**, stderr empty, stdout 1611 bytes). `·` = space, `»` = TAB, `¶` = LF.

```
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen¶
¶
The following parameters are in effect:¶
···················Binary File :·················(-bname)¶
¶
Additional Options¶
·········Close Relative Inference : --related, --duplicate¶
···Pairwise Relatedness Inference : --kinship, --ibdseg, --ibs, --makeGRM¶
··············Inference Parameter : --degree, --noscreen [-1717986816],¶
····································--seglength, --minConc [0.80]¶
·········Relationship Application : --unrelated, --cluster, --build¶
························QC Report : --bysample, --bySNP, --roh, --autoQC¶
·····················QC Parameter : --callrateN, --callrateM¶
·············Population Structure : --pca, --mds¶
··············Structure Parameter : --projection, --pcs¶
··········Quantitative Trait GWAS : --lmm¶
················Binary Trait GWAS : --tdt, --gdt¶
················Association Model : --trait [], --covariate [], --maxP¶
·····Association Method Parameter : --invnorm¶
···············Genetic Risk Score : --risk, --model [], --prevalence, --noflip¶
··············Computing Parameter : --cpus¶
···················Optional Input : --fam [], --bim [], --phefile [],¶
····································--covfile [], --prunedsnp [],¶
····································--sexchr [23]¶
···························Output : --rplot, --pngplot, --plink¶
·················Output Parameter : --prefix [king], --rpath []¶
¶
¶
FATAL ERROR - ¶
Genotype files are required. e.g.,¶
··king -b ex.bed --related¶
¶
Please check the reference paper Manichaikul et al. 2010 Bioinformatics,¶
»»»»»Chen et al. 2024,¶
··········or the KING website at kingrelatedness.com¶
¶
```

Byte-exact notes:
- Line 1 banner is 39 chars; line 3 header is 39 chars (coincidence).
- `FATAL ERROR - ` has a **trailing space** before the LF.
- The `Chen et al. 2024,` line is indented with **exactly 5 TAB characters (0x09)**, not spaces.
- The last footer line is indented with **10 spaces**.
- Output ends with `\n\n`.

`--related` (exit 1) is byte-identical except line 7 becomes
`·········Close Relative Inference : --related [ON], --duplicate`.

---

## 3. Verbatim console output — the other probes

### 3a. `king -b nonexistent.bed --related` (exit **1**, 1491 bytes)

Options block as above (with `--related [ON]` and the Binary File line filled in), then:

```
···················Binary File : nonexistent.bed (-bname)¶
…¶
·················Output Parameter : --prefix [king], --rpath []¶
¶
KING starts at Thu Aug 13 17:55:44 2026¶
¶
FATAL ERROR - ¶
Genotype file nonexistent.bed cannot be opened¶
¶
```

Note: **no** "Please check the reference paper…" epilogue — that epilogue is part of the
*text* of the "Genotype files are required" message, not of the FATAL frame.

### 3b. `king --bogusflag` (exit **1**, 1717 bytes)

Options block (unchanged from the no-args form), then:

```
…--prefix [king], --rpath []¶
¶
<BEL>WARNING - ¶
Problems encountered parsing command line:¶
¶
Command line parameter --bogusflag is undefined¶
¶
¶
FATAL ERROR - ¶
Genotype files are required. e.g.,¶
…
```

**`<BEL>` is a literal 0x07 byte immediately preceding `WARNING`.** It is present on every
WARNING block and absent from FATAL ERROR blocks. Exact bytes:
`…--rpath []\n` `\n` `\n\x07WARNING - \nProblems encountered parsing command line:\n\n` …

### 3c. Exact frame grammar (derived, verified against all captures)

```
options_block_trailer := "\n"
warning_block         := "\n\x07WARNING - \nProblems encountered parsing command line:\n\n"
                         + problem_line*                       (each ends "\n")
                         + "\n"
starts_line           := "KING starts at " + ctime()           (ctime already ends "\n")
fatal_block           := "\nFATAL ERROR - \n" + message + "\n\n"
```

Cross-checks: no-args → `\n` + `\nFATAL…`; missing-bed → `\n` + `KING starts at…\n` +
`\nFATAL…`; bogusflag → `\n` + `\n\x07WARNING…undefined\n` + `\n` + `\nFATAL…`;
bogusflag + valid bed → `\n` + `\n\x07WARNING…undefined\n` + `\n` + `KING starts at…`.

### 3d. Problem-line forms (in command-line order, one per problem)

| Situation | Line |
|---|---|
| Unknown long option | `Command line parameter --bogusflag is undefined` |
| Ambiguous long-option prefix | `Command line parameter --r is ambiguous` |
| Unknown short option | `Command line parameter -q (#3) ignored` |
| Bare positional token | `Command line parameter positional (#2) ignored` |
| Non-numeric value token after a numeric option | `Command line parameter hello (#2) ignored` |

`#N` is the **1-based index into `argv[1..]`** of the offending token.
Real capture (`king --bogus1 positional -q --bogus2 --r`):

```
Command line parameter --bogus1 is undefined¶
Command line parameter positional (#2) ignored¶
Command line parameter -q (#3) ignored¶
Command line parameter --bogus2 is undefined¶
Command line parameter --r is ambiguous¶
```

---

## 4. The complete option table

This build defines **exactly 46 long options and 1 short option** — the ones printed in the
help block. Confirmed by sweeping **5,402 candidate names** harvested from every printable
string in the binary: every accepted token is a unique case-insensitive **prefix** of one of
these 46; nothing else is accepted.

Notably **absent** (they appear in the binary's strings as message text or as
"Options in effect:" echo formats, but are *not* accepted on the command line in 2.3.2):
`--ibdall --ibdmds --ibdmap --ibdH2 --ibdMI --ibdgdt --popibd --popdist --poproh --porel
--HEreg --herit --homog --homomap --mthomo --mtscore --grm --grm-lmm --makePC --pcgdt
--popgdt --casecontrol --ancestry --aucmap --npl --lessmem --exact --distant --merlin
--paternity --phase --linear --poodt --novclmm --strat --nperm --search --faster --slower
--errorrate --minMAF --mincons --prunesnp --position --sysbit --bysnp --model` (last two
are aliases/echo forms). Treat these as **not part of the 2.3.2 CLI**.

### 4a. Table (order = order of appearance in the help block = print order)

Column *Type* legend:
`SW` switch · `INT` integer · `SINT` "smart" integer · `DBL` double · `STR` string

| # | Section header (right-aligned label) | Option | Type | Default | Rendered when unset |
|---|---|---|---|---|---|
| 1 | `Close Relative Inference` | `--related` | SW | off | `--related` |
| 2 | | `--duplicate` | SW | off | `--duplicate` |
| 3 | `Pairwise Relatedness Inference` | `--kinship` | SW | off | `--kinship` |
| 4 | | `--ibdseg` | SW | off | `--ibdseg` |
| 5 | | `--ibs` | SW | off | `--ibs` |
| 6 | | `--makeGRM` | SW | off | `--makeGRM` |
| 7 | `Inference Parameter` | `--degree` | SINT | 0 | `--degree` |
| 8 | | `--noscreen` | INT | **uninitialised** (see §7) | `--noscreen [-1717986816]` |
| 9 | | `--seglength` | DBL | unset (NaN) | `--seglength` |
| 10 | | `--minConc` | DBL | **0.8** | `--minConc [0.80]` |
| 11 | `Relationship Application` | `--unrelated` | SW | off | `--unrelated` |
| 12 | | `--cluster` | SW | off | `--cluster` |
| 13 | | `--build` | SW | off | `--build` |
| 14 | `QC Report` | `--bysample` | SW | off | `--bysample` |
| 15 | | `--bySNP` | SW | off | `--bySNP` |
| 16 | | `--roh` | SW | off | `--roh` |
| 17 | | `--autoQC` | SW | off | `--autoQC` |
| 18 | `QC Parameter` | `--callrateN` | DBL | unset (NaN) | `--callrateN` |
| 19 | | `--callrateM` | DBL | unset (NaN) | `--callrateM` |
| 20 | `Population Structure` | `--pca` | SW | off | `--pca` |
| 21 | | `--mds` | SW | off | `--mds` |
| 22 | `Structure Parameter` | `--projection` | SINT | 0 | `--projection` |
| 23 | | `--pcs` | SINT | 0 | `--pcs` |
| 24 | `Quantitative Trait GWAS` | `--lmm` | SW | off | `--lmm` |
| 25 | `Binary Trait GWAS` | `--tdt` | SW | off | `--tdt` |
| 26 | | `--gdt` | SW | off | `--gdt` |
| 27 | `Association Model` | `--trait` | STR | `""` | `--trait []` |
| 28 | | `--covariate` | STR | `""` | `--covariate []` |
| 29 | | `--maxP` | DBL | unset (NaN) | `--maxP` |
| 30 | `Association Method Parameter` | `--invnorm` | SW | off | `--invnorm` |
| 31 | `Genetic Risk Score` | `--risk` | SW | off | `--risk` |
| 32 | | `--model` | STR | `""` | `--model []` |
| 33 | | `--prevalence` | DBL | unset (NaN) | `--prevalence` |
| 34 | | `--noflip` | SW | off | `--noflip` |
| 35 | `Computing Parameter` | `--cpus` | SINT | 0 | `--cpus` |
| 36 | `Optional Input` | `--fam` | STR | `""` | `--fam []` |
| 37 | | `--bim` | STR | `""` | `--bim []` |
| 38 | | `--phefile` | STR | `""` | `--phefile []` |
| 39 | | `--covfile` | STR | `""` | `--covfile []` |
| 40 | | `--prunedsnp` | STR | `""` | `--prunedsnp []` |
| 41 | | `--sexchr` | INT | **23** | `--sexchr [23]` |
| 42 | `Output` | `--rplot` | SW | off | `--rplot` |
| 43 | | `--pngplot` | SW | off | `--pngplot` |
| 44 | | `--plink` | SW | off | `--plink` |
| 45 | `Output Parameter` | `--prefix` | STR | **`king`** | `--prefix [king]` |
| 46 | | `--rpath` | STR | `""` | `--rpath []` |

### 4b. Short options

Only **`-b`** exists (matched **case-insensitively**, so `-B` works too).
Both forms accepted: `-b file.bed` (separate token) and `-bfile.bed` (attached).
`-b` as the last token with nothing following is silently accepted and leaves the value empty.
Every other single-dash token produces `Command line parameter -X (#N) ignored`.
The label shown in the help block for this parameter is `Binary File`, hint text `(-bname)`.

---

## 5. Long-option matching rules

1. **Case-insensitive.** `--RELATED`, `--ReLaTeD`, `--related` are equivalent.
2. **Unique-prefix matching.** The supplied token must be a prefix of exactly one option
   name. `--re`, `--rel`, `--k`, `--kin`, `--ibd`, `--SEX`, `--Make` all resolve.
3. **Ambiguous prefix → rejected** with `is ambiguous`, not silently resolved. Observed
   ambiguous 2-char prefixes: `by` (bysample/bySNP), `ca` (callrateN/M), `co`
   (covariate/covfile), `ib` (ibdseg/ibs), `ma` (makeGRM/maxP), `no` (noscreen/noflip),
   `pc` (pca/pcs), `pr` (prefix/prevalence/projection/prunedsnp), `rp` (rplot/rpath),
   `se` (seglength/sexchr). Also `--` alone is *ambiguous* (empty prefix matches everything).
4. **Longer-than-name → undefined.** `--relatedXX`, `--pca2`, `--related5` are rejected.
5. Unrecognised → `is undefined`. Neither `undefined` nor `ambiguous` aborts the run.

### 5a. Value-token consumption (per type)

| Type | Bare (no value / value not consumed) | Token accepted iff | Parse |
|---|---|---|---|
| `SW` | sets ON | never consumes a token | — |
| `SINT` | sets **1** | token matches `^[+-]?[0-9]+$` entirely | 32-bit `int`, wraps (`2147483648`→`-2147483648`); `007`→7 (decimal, not octal) |
| `INT` | sets **0** | same as SINT | same |
| `DBL` | leaves value **unchanged** | first char ∈ `[0-9] . + -` | `strtod`-like: `3x`→3.0, `1e3`→1000.0, `0x10`→16.0, `1,5`→1.0, `.5`→0.5, `-`→0.0. Trailing/leading space, `hello`, `INF`, `NaN`, `e5` are **rejected** |
| `STR` | leaves value unchanged | **always consumes** the next token, even another option | verbatim |

Consequences worth reproducing:
- `--cpus foo` → `--cpus [1]` **and** `Command line parameter foo (#2) ignored`.
- `--minConc foo` → value stays `0.80` **and** `foo (#2) ignored`.
- `--trait --related` → `--trait [--related]`, and `--related` is **not** turned on, **no warning**.
- `--degree` followed by a valid option (e.g. `--degree --related`) does not consume it:
  `--degree [1]`, `--related [ON]`, no warning.
- A value-taking option as the **last** argv token consumes nothing and warns nothing.

---

## 6. Exact layout algorithm of the parameters-in-effect block

### 6a. The "Binary File" line (top block)

```
printf("%30s : %15s (-b%s)\n", "Binary File", <bed value>, "name")
```
- Label **right-aligned in a 30-column field** (`:` lands at 1-based column 32).
- Value **right-aligned in a 15-column field**, overflowing to the right if longer.
- Followed by one space and the literal hint `(-bname)`.
- Unset value → 15 spaces, so the line is `…Binary File : ` + 15 spaces + ` (-bname)` = **57 chars**.

Verified: values of length 1/3/10/15 all yield line length 57; length 16 → 58; length 30 → 72.

### 6b. The "Additional Options" groups

Header line is the literal `Additional Options` (18 chars, no colon).

For each group, in table order:

```
LABEL_W = 33      # group label right-aligned in 33 columns  → ':' at 1-based column 35
MAXCOL  = 78      # hard wrap budget, EXCLUDING the trailing comma

line = f"{label:>33} :"                 # length 35, note: NO trailing space
for k, item in enumerate(items):        # item = the rendered "--name" / "--name [value]"
    if len(line) + 1 + len(item) > 78:
        emit(line)                      # flush current line as-is
        line = " " * 35                 # continuation stub, same width as "label :"
    line += " " + item                  # separator is always a single leading space
    if k != last:
        line += ","                     # comma may push the line to 79 chars
emit(line)
```

Key details, each confirmed by experiment:

- **Continuation indent is 36 spaces** (35-space stub + the separator space), i.e. the
  continued item starts at 1-based column 37 — exactly under the first item of line 1.
- **The comma is not counted in the fit test.** A line may therefore reach **79** characters
  when the last item that fits is followed by a comma. Observed real 79-char line:
  `                Association Model : --trait [XXXXXXXXXXXXXXXX], --covariate [],`
- **An over-long single item is not broken.** If the first item does not fit, the line is
  flushed immediately — producing a bare `{label:>33} :` line of exactly **35 characters
  with no trailing space** — and the item then overflows past column 78 on the next line.
  Example (`--trait` with a 40-char value):
  ```
  [ 35] |                Association Model :|
  [ 87] |                                    --trait [XXXX…XXXX],|
  [ 58] |                                    --covariate [], --maxP|
  ```
- The block **re-wraps dynamically** as values change. E.g. `--noscreen 0` removes the
  `[-1717986816]` and the Inference Parameter group re-flows from
  `--degree, --noscreen [-1717986816],` / `--seglength, --minConc [0.80]` to
  `--degree, --noscreen, --seglength,` / `--minConc [0.80]`.

Worked check against the shipped default block:

| Group | Arithmetic | Result |
|---|---|---|
| Genetic Risk Score | 35+1+6=42, +1+10=53(+`,`54)… final 78 | 1 line, len 78 (fits exactly) |
| Genetic Risk Score with `--risk [ON]` | …73(+`,`74), then 74+1+8=83 > 78 | wraps; line1 len 74 |
| Inference Parameter | 35+1+8=44(45), +1+24=70(71), 71+1+11=83 > 78 | wraps; lines 71 / 65 |
| Optional Input | …68(69), 69+1+12=82 > 78 → wrap; then 64(65), 65+1+13=79 > 78 → wrap | lines 69 / 65 / 49 |

The label field width 33 equals `max(len(group label)) + 3` (longest is
`Pairwise Relatedness Inference`, 30 chars); the top block's 30 does **not** follow that
rule (`Binary File` is 11), so treat both 30 and 33 as **hard-coded constants**.

### 6c. Value rendering

| Type | Rendered | Condition |
|---|---|---|
| `SW` | `--name [ON]` / `--name` | on / off |
| `INT`, `SINT` | `--name [%d]` / `--name` | **shown iff value ≠ 0** |
| `DBL` | `--name [<fmt>]` / `--name` | **shown iff the value is not NaN** — i.e. always once set, *including an explicit 0* → `[0.00]` |
| `STR` | `--name [%s]` | **always** shown, empty string → `[]` |

Unset doubles (`seglength`, `callrateN`, `callrateM`, `maxP`, `prevalence`) are NaN and
therefore hidden; `minConc` is initialised to 0.8 and always shown.

### 6d. Double number format (exact)

```
value == 0.0 or value >= 0.01   →   "%.2f"
otherwise (0 < value < 0.01, or value < 0)  →  "%.1e"
```
No `fabs()` — **every negative value uses scientific notation.** Measured:

| input | rendered | | input | rendered |
|---|---|---|---|---|
| `0` | `0.00` | | `0.009` | `9.0e-03` |
| `0.01` | `0.01` | | `0.001` | `1.0e-03` |
| `0.1` | `0.10` | | `0.0001` | `1.0e-04` |
| `0.5` | `0.50` | | `1e-5` | `1.0e-05` |
| `0.999` | `1.00` | | `1e-10` | `1.0e-10` |
| `1` | `1.00` | | `-0.5` | `-5.0e-01` |
| `9.995` | `9.99` | | `-1` | `-1.0e+00` |
| `12345.678` | `12345.68` | | `-100` | `-1.0e+02` |
| `1e20` | `100000000000000000000.00` | | `1e300` | full 301-digit `%.2f` expansion |

Integers are plain `%d` (32-bit, signed, wrapping).

---

## 7. ⚠ The `--noscreen` display bug (required for byte parity)

`--noscreen` shows a garbage default `[-1717986816]`. This is **not** a constant — it is a
memory-aliasing artifact, fully characterised by black-box probing:

> `noscreen` is an `int` whose storage **overlaps** the `double minConc`, with
> `&minConc == &noscreen + 1`. Both are written at full width, so each clobbers the other.
> `noscreen` is never initialised (its own byte is observed as `0x00`); `minConc` is
> initialised to `0.8` at startup. The status printer reads 4 bytes at `&noscreen`.

Default case: `0.8` = `0x3FE999999999999A` little-endian → bytes `9A 99 99 99 99 99 E9 3F`.
Reading `[0x00, 0x9A, 0x99, 0x99]` as little-endian `int32` = `0x99999A00` = **-1717986816**. ✔

Exact emulation (validated on 600 randomized runs, including argument ordering):

```python
buf = bytearray(16)
struct.pack_into('<d', buf, 1, 0.8)          # startup init of minConc
for name, v in options_in_command_line_order:
    if name == 'noscreen': struct.pack_into('<i', buf, 0, v)
    elif name == 'minConc': struct.pack_into('<d', buf, 1, v)
displayed_noscreen = struct.unpack_from('<i', buf, 0)[0]
displayed_minConc  = struct.unpack_from('<d', buf, 1)[0]
```

Confirmations: `--minConc 2.5` → `0` (hidden); `--minConc 0.0001` → `474164480`;
`--noscreen 5 --minConc 0.7` → `1717986821`; `--noscreen 999999 --minConc 1.0` → `63`;
`--minConc 123456.789 --noscreen 200` → `200`; `--noscreen 1 --minConc 0.9` → `-858993407`.

The reverse clobber (a `--noscreen` write corrupting the low 3 bytes of `minConc`) is
invisible at `%.2f`, e.g. `--minConc 0.7 --noscreen 5` still prints `0.70`.

**Recommendation for the reimplementation:** expose this as a bug-compatibility switch.
The default no-args banner must print `--noscreen [-1717986816]` for byte parity, but the
value is a stack/BSS artifact and may differ on other platforms/builds — do not treat it as
a documented constant.

---

## 8. Runtime console output (successful run)

`king -b t.bed --related --prefix demo` (exit **0**), after the options block:

```
¶
KING starts at Thu Aug 13 17:55:32 2026¶
Loading genotype data in PLINK binary format...¶
Read in PLINK fam file t.fam...¶
··PLINK pedigrees loaded: 6 samples¶
Read in PLINK bim file t.bim...¶
··Genotype data consist of 2000 autosome SNPs¶
··PLINK maps loaded: 2000 SNPs¶
Read in PLINK bed file t.bed...¶
0%<CR>6%<CR>13%<CR>19%<CR>25%<CR>31%<CR>38%<CR>44%<CR>50%<CR>56%<CR>63%<CR>69%<CR>75%<CR>81%<CR>88%<CR>··PLINK binary genotypes loaded.¶
94%<CR>··KING format genotype data successfully converted.¶
¶
--related is replaced with --kinship for a small sample size.¶
Autosome genotypes stored in 32 words for each of 6 individuals.¶
¶
Options in effect:¶
»--kinship¶
»--prefix demo¶
¶
Within-family kinship data saved in file demo.kin¶
Relationship inference across families starts at Thu Aug 13 17:55:32 2026¶
8 CPU cores are used.¶
·········································ends at Thu Aug 13 17:55:32 2026¶
Between-family kinship data saved in file demo.kin0¶
Note --kinship --degree <n> can filter & speed up the kinship computing.¶
KING ends at Thu Aug 13 17:55:32 2026¶
¶
```

Notes:
- `KING starts at %s` / `KING ends at %s` where `%s` is `ctime()` output, which **already
  ends in `\n`** — hence the extra blank line before a following `FATAL ERROR`.
  Format is `Thu Aug 13 17:55:32 2026` (C `ctime`, day-of-month space-padded).
- Progress is `%d%%` followed by `\r`; percentages here were 0,6,13,…,88 then 94.
  The `94%` appears *after* the `PLINK binary genotypes loaded.` line (a second counter).
- **`Options in effect:`** is a second, distinct block: literal header, then one
  **TAB-indented** line per effective option, then a blank line. It echoes the *resolved*
  options (`--related` was rewritten to `--kinship`), not the raw command line.
- The `ends at` line is padded with **41 leading spaces** so `ends at` aligns under
  `starts at` of the preceding `… starts at …` line.
- With a larger cohort the run also emits a `Relationship summary (total relatives: N by
  pedigree, M by inference)` table with TAB-separated columns
  `··Source»MZ»PO»FS»2nd»3rd»OTHER`, a `··====…` rule of 59 `=`, then `··Pedigree»…` and
  `··Inference»…` rows.

### 8a. "Options in effect:" echo formats found in `__cstring`

These printf templates exist in the binary and drive that block (many belong to options not
exposed on the 2.3.2 CLI — see §4):

`--cpus %d` · `--prefix %s` · `--degree %d` · `--minConc %g` · `--pcs %d` · `--projection %d` ·
`--trait %s` · `--covariate %s` · `--model %s` · `--prevalence %G` · `--grm %s` · `--strat %s` ·
`--prunesnp %s` · `--faster %d` · `--slower %d` · `--errorrate %lf` · `--minMAF %lf` ·
`--mincons %d` · `--nperm %d` · `--search %d-%d` · `--sysbit 64` · `--degree 5`
plus bare switch names for every switch option.

Note `--minConc` is echoed with `%g` here (→ `0.8`), **not** the `%.2f` used in the banner,
and `--prevalence` with `%G`.

---

## 9. Other exact literals harvested (relevant to CLI/console parity)

```
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen
The following parameters are in effect:
Additional Options
Genotype files are required. e.g.,
  king -b ex.bed --related
Please check the reference paper Manichaikul et al. 2010 Bioinformatics,
\t\t\t\t\tChen et al. 2024,
          or the KING website at kingrelatedness.com
KING starts at %s
KING ends at %s
Genotype file %s cannot be opened
Please use PLINK binary format as input.
Please use either PLINK or KING binary format as input.
Sex chromosome %d out of range.
Non-human samples are analyzed, with %d pairs of chromosomes
Minimum segment length is set as %d bp
KING supports minimum segment length from 1 to 10 Mb at the moment.
Default seglength of 3Mb is used.
minConc value is out of range and not specified.
Please use --model <file> to specify a risk model.
Please specify one of the following %d options:
 --%s
The following analyses will run separately: 
Trait %s cannot be found.
Covariate %s cannot be found.
Covariate %s is duplicated.
Binary trait %s is now considered as a quantitative trait.
No covariates are included in the analysis.
--related is replaced with --kinship for a small sample size.
--related is skipped for a rather small sample size.
--related is skipped.
--roh is skipped.
--ibdseg is skipped.
--kinship analysis carried out instead for such a small sample size.
  For additional relative pair plots please use --degree, --degree 2, or --degree 3.
Please do not run --related together with --autoQC
Please do not run --roh together with --autoQC
Please do not run --ibdseg together with --autoQC
R plot for --%s is not available.        (instantiated per analysis)
```

---

## 10. Reimplementation checklist (byte-parity)

1. Print everything to **stdout**; keep stderr empty. Exit 1 on argument/input failure, 0 otherwise.
2. Banner exactly `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen` + `\n\n`.
3. Top block: `%30s : %15s (-b%s)` with label `Binary File`, hint word `name`.
4. Groups: label right-aligned in 33, `" :"`, then items joined by `" "` with `","`
   appended after each non-final item; wrap when `len(line)+1+len(item) > 78`;
   continuation stub is 35 spaces.
5. Value rendering per §6c/§6d — note *negative doubles use `%.1e`* and *explicit-zero
   doubles still print `[0.00]`* while zero ints print nothing.
6. Emit exactly one blank line after the options block.
7. WARNING block is prefixed with a **0x07 BEL**; FATAL ERROR is not. `FATAL ERROR - ` has a
   trailing space. The reference footer uses **5 TABs**.
8. Long options: case-insensitive **unique-prefix** match; ambiguity is an error, not a
   resolution; longer-than-name is undefined; `-b`/`-B` is the only short option.
9. Reproduce the `--noscreen`/`--minConc` aliasing (§7) behind a compatibility flag.
10. `KING starts at`/`KING ends at` use `ctime()` output verbatim (trailing `\n` included).
