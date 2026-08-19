# KING 2.3.2 — Clean-Room Reimplementation Specification

**Target:** an MIT-licensed Rust reimplementation of KING 2.3.2's relatedness inference,
byte-compatible with the reference binary's output files on the supported flag set.

**Reference binary (oracle, read-only, never redistributed):**
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
Mach-O 64-bit arm64, 1 815 336 bytes, banner `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`.

This document is self-contained. A Rust engineer should not need to read any of the nine
recon reports (`01-…09-…`) to implement it. Where a claim is uncertain it says so and
§8 gives the experiment that settles it.

---

## 0. Preamble

### 0.1 Clean-room posture (binding — read before writing any code)

KING's C++ source **must never be fetched, opened, or read** by any human or agent on this
project. It is available at `kingrelatedness.com/KINGcode.tar.gz` and mirrored in
`statgen/topmed_variant_calling/king/`. Both are on a hard blacklist.

Rationale: the author's only published grant is the sentence *"Feel free to use KING for your
research, but please do not redistribute AND make profits."* — not an OSI licence, with no
reproduction or derivative-works grant. Bioconda simultaneously labels the same tarball
`GPL-3.0-or-later`, and the tarball statically bundles Abecasis's libStatGen, which *is*
GPL-3-or-later. Under the website reading we have no licence to copy; under the GPL reading any
derivation forces GPL-3 onto a crate we ship inside a commercial-capable desktop app. Either
way there is no MIT path through the source, and the algorithm is short and published, so the
clean-room route costs almost nothing.

**Permitted inputs (whitelist):**
1. Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen W-M. "Robust relationship
   inference in genome-wide association studies." *Bioinformatics* 2010;26(22):2867–2873
   (PMC3025716). Mathematical methods are uncopyrightable.
2. kingrelatedness.com manual / tutorial / download pages (paraphrase; do not copy prose).
3. `strings`, `nm`, `otool -L`, `otool -tV` on the shipped binary, and **running** it as a
   black box and diffing its output bytes. Format strings and column headers are facts about
   an interface, not protected expression.
4. MIT-licensed reimplementations, for cross-checking the scalar only: Hail `hl.king`, cuKING,
   somalier. Docs (not source) of GPL tools: plink2, SNPRelate.

**Prohibited:** KING source in any form; the R scripts embedded in the binary (they are
authored code even though `strings` reveals them); plink2/SNPRelate/akt source; prompting any
LLM to "recall" or "show" KING's source. Illumina **akt** is PolyForm Strict 1.0.0 — struck
from the reference set entirely.

**Process controls:** the module carries a provenance header naming this whitelist; a
`THIRD_PARTY.md` entry records that KING is neither vendored, linked, nor redistributed;
fixtures + expected-output text (program output, not source) are the committed oracle.

**Separate action item:** `runtimes/king/king` is currently committed and shipped. Delete
`runtimes/king/` once this port reaches parity; keep the binary only at the external read-only
path above, as a local oracle.

### 0.2 Evidence tags used throughout

| Tag | Meaning |
|---|---|
| **[V]** | Verified by running the reference binary and diffing bytes, in the recon phase or in this session. Highest confidence. |
| **[V-NEW]** | Verified *while writing this spec*, with the check reproduced in §A.3. Supersedes earlier recon guesses where they conflict. |
| **[S]** | Read from the binary's `__TEXT,__cstring` literals, pinned to a named writer function via `otool -tV` literal-pool annotations. Facts about the output format. |
| **[P]** | Manichaikul et al. 2010. |
| **[D]** | Public documentation (kingrelatedness.com). Weakest — the site documents ~2.2.7 and its whitespace is a rendering artifact. |
| **[?]** | Unresolved. Cross-referenced to a numbered item in §8. |

**Conflict rule:** binary observation **[V]/[V-NEW]** beats binary strings **[S]** beats the
paper **[P]** beats the website **[D]**. Every conflict actually encountered is called out
inline.

### 0.3 Notation

* `\t` = one 0x09 byte. `\n` = one 0x0A byte. Shown explicitly everywhere; never rely on
  visual alignment. **KING never pads fields** — one separator byte between fields, no
  alignment, no trailing separator (with one bug exception, §6.11).
* `%.4lf` etc. are C printf conversions. Rust's `format!("{:.4}", x)` and C's `%.4lf` both
  round the *exact binary double* correctly and agree bit-for-bit on all finite values,
  including ties. Do **not** pre-round intermediates; format once from the full-precision f64.
* `%G` appears only for the `Error` column and `MI_Removal`, whose only observed values are
  `0`, `0.5`, `1`, rendered exactly as those three strings. Implement as:
  ```rust
  fn fmt_g(x: f64) -> String {           // C "%G", 6 significant digits, trailing zeros dropped
      if x == 0.0 { return "0".into() }
      let s = format!("{:.5E}", x);      // then normalise; for the three legal values this is
      // 0 -> "0", 0.5 -> "0.5", 1 -> "1"
      shortest_g(x, 6)
  }
  ```
  A three-way match on `{0.0, 0.5, 1.0}` is acceptable for v1; assert on any other value.
* `<p>` denotes the value of `--prefix` (default `king`). Output paths are **literal string
  concatenation**, not stem+separator: `--prefix ZZ_` gives `ZZ_.kin` and `ZZ_allsegs.txt`.
  **[V]**
* Sample indices `i`, `j` are 0-based positions in `.fam` order unless stated otherwise.

---

## 1. SCOPE

### 1.1 Tier 1 — relatedness core (v1 byte-parity target)

These five analyses are the product requirement and the v1 parity target.

| Flag | Analysis | Outputs | Notes |
|---|---|---|---|
| `--related` | Integrated relationship inference (kinship + IBD segments) | `<p>.kin` (16 col), `<p>.kin0` (14 col), `<p>allsegs.txt` | **Split delivery**, see §1.1.1 |
| `--kinship` | KING-robust kinship for all pairs | `<p>.kin` (10 col), `<p>.kin0` (8 col) | Fully implementable now; no dependency on the segment engine |
| `--duplicate` | Duplicate / MZ detection by heterozygote concordance | `<p>.con` | |
| `--ibs` | Full IBS/concordance statistics per pair | `<p>.ibs`, `<p>.ibs0`, `<p>allsegs.txt` | Two trailing columns depend on the segment engine (§6.9) |
| `--unrelated` | Greedy maximal unrelated subset | `<p>unrelated.txt`, `<p>unrelated_toberemoved.txt`, `<p>allsegs.txt` | Selection algorithm is **[?]** — §8 item 15 |

Supporting flags that are **in Tier 1** because these five analyses consume them:

`-b` / `-B` (single fileset only), `--prefix`, `--degree`, `--minConc`, `--cpus`, `--fam`,
`--bim`, `--sexchr`, `--noscreen`.

**Tier 1 also includes**, as hard requirements:
* the complete console surface of §2 (banner, parameters-in-effect block, warning/fatal
  frames, exit codes) — it is how the parity harness diffs runs;
* `.bed`/`.bim`/`.fam` parsing and validation of §3;
* the bit-plane popcount engine of §4.4, which every Tier-1 and Tier-2 analysis shares.

#### 1.1.1 `--related` is delivered in two phases

`--related`'s 16-column `.kin` / 14-column `.kin0` carry `IBD1Seg`, `IBD2Seg`, `PropIBD`,
`InfType`, which come from the IBD-segment engine (Tier 2). Therefore:

* **Phase 1a (Tier 1, no Tier-2 dependency).** Implement `--related`'s *degraded* paths, which
  are byte-identical to `--kinship`: the `N < 10` sample downgrade (`--related is replaced with
  --kinship for a small sample size.` **[V]**) and the no-usable-segments path
  (`Relationship inference will be based on kinship estimation only.` **[S]**). This makes
  `--related` correct-and-parity-complete for small and segment-less inputs.
* **Phase 1b (needs Tier 2's segment engine).** The full 16/14-column form. Schedule it with
  `--ibdseg`; do not attempt it before the segment engine passes §8 item 16.

### 1.2 Tier 2 — second wave

| Flag | Analysis | Outputs | Blocking unknowns |
|---|---|---|---|
| `--ibdseg` | Pairwise IBD-segment inference | `<p>.seg`, `<p>allsegs.txt`, `<p>splitped.txt` | Segment *calling* rule (§8 item 16); usable-segment cut rule (§8 item 17) |
| `--build` | Pedigree reconstruction from SNP data | `<p>updateids.txt`, `<p>updateparents.txt`, `<p>build.log`, `<p>splitped.txt` | The data-derived FS/PO IBS0 cutoff (§8 item 18) |
| `--bysample` | Sample-level QC | `<p>bySample.txt` (6 dynamic header variants) | Mendelian-error accounting |
| `--bySNP` | SNP-level QC | `<p>bySNP.txt` (3 dynamic header variants) | Mendelian-error accounting; a non-deterministic KING fatal at tiny N (§8 item 22) |
| `--autoQC` | Call-rate + sex QC pipeline | `<p>_autoQC_Summary.txt`, `<p>_autoQC_snptoberemoved.txt`, `<p>_autoQC_sampletoberemoved.txt`, `<p>_autoQC_updatesex.txt` | **None — implemented.** Layouts and rules: VERIFIED_FORMULAS.md § `--autoQC`; byte-identical on all 13 datasets, PARITY.md §3 |

Tier-2 supporting flags: `--seglength`, `--callrateN`, `--callrateM`.

Tier 2 also includes `<p>allsegs.txt`, which several Tier-1 analyses emit as a side effect
(§6.10) — it is specified here because Tier 1 needs it, but its *content* depends on the
Tier-2 usable-segment rule.

### 1.3 Deliberately excluded product surface

Accepted on the command line so the banner and parser remain compatible, but not implemented
by this minimal relatedness/QC package. [SCOPE.md](SCOPE.md) is the binding product-scope
statement. Clear unsupported-option diagnostics are desirable, but implementing these
analysis families is not a parity requirement for the supported core.

* **Analyses:** `--makeGRM`, `--roh`, `--pca`, `--mds`, `--lmm`, `--tdt`, `--gdt`,
  `--risk`, `--invnorm`, `--plink`.
* **Parameters of Tier-3 analyses:** `--projection`, `--pcs`, `--trait`, `--covariate`,
  `--maxP`, `--model`, `--prevalence`, `--noflip`, `--phefile`, `--covfile`, `--prunedsnp`.
* **Plotting:** `--rplot`, `--pngplot`, `--rpath`. KING shells out to `R CMD BATCH` with
  scripts embedded as string constants. Those scripts are authored code — **do not
  transcribe**. GeneQuire renders its own MDX/SVG figures; the Rust port has **zero** R
  dependency, which is a product win.
* ~~**X-chromosome analysis**~~ — **brought in scope.** `<p>X.kin`, `<p>X.kin0` and
  `<p>X.seg` are all implemented and byte-identical on every capture that produces them.
  Each of the three X passes has its own gate, and they are not the same gate:
  `--kinship`'s needs 512 or more X markers, no `--degree` and more than one family;
  `--related`'s and `--ibdseg`'s need only that the X map yields a **usable segment**
  (`crate::analysis::xseg`), and `--ibdseg`'s needs `--degree` as well. (`--sexchr` was
  always Tier 1 because it changes which variants count as *autosomal* — §3.4.)
* **Multi-dataset input** (`-b a.bed,b.bed`, `--fam a.fam,b.fam`) with its auto-flip / ambiguous-SNP
  merge rules. Reject a comma in `-b` with a clear error.
* **Retired / never-exposed flags:** `--homog` and `--mtscore` (retired in 2.3.0); and the ~40
  names that exist only as strings in the binary and are *not* accepted by the 2.3.2 CLI
  (`--ibdall`, `--HEreg`, `--exact`, `--distant`, `--porel`, `--paternity`, `--lessmem`, …).
  Do not implement, do not accept.
* **KING bug-compatibility** beyond what §2 requires: the `X.kin0` OpenMP data race and
  the mis-indexed exclusion list are **not** reproduced. The `X.seg` header/row mismatch
  **is** — it moved in scope with X.seg itself, and a file whose eleven-name header sits
  over nine-value rows is what the reference writes, so parity requires writing it too
  (`crate::analysis::xseg`). The two console bugs we also reproduce are the `--noscreen`
  garbage default (§2.7) and the `3nd-degree` typo (§2.8), both behind the
  `king-bugcompat` feature flag.

### 1.4 Explicit non-goals for v1

1. Bit-identical *performance* characteristics (KING's two-stage screening); we may compute all
   pairs directly provided the numbers match (§5.5).
2. `.segments.gz` / `.rohseg.gz`. **The reference binary never writes them** — this build has
   no zlib in the segment writer (`Engine::AllIBDSegments()` contains only
   `--ibdall cannot run without ZLIB`) **[S][V]**. The manual documents `.segments.gz`; it is
   unreachable in 2.3.2. Do not target it.
3. Reproducing KING's temp files `<p>$$$.kin0` / `<p>$$$.ibs0` (visible mid-run, deleted after).
4. Parity against the website's published example outputs. IBD-segment numerics changed
   materially at 2.1.2, 2.1.3, 2.2.1, 2.2.5, 2.2.6 and 2.2.7 — parity is meaningful **only**
   against the 2.3.2 binary in hand.
5. Enforcing `--cpus` as a strict Rayon thread cap. The option remains part of the
   compatibility parser and console surface; deterministic results, not performance parity,
   are the requirement.

### 1.5 Library-first shape

The CLI is a parity harness, not the product. The crate exposes a typed API
(`Bundle -> RelatednessReport`) and the `king`-compatible CLI is a thin `bin/` on top. Every
number in every output file must be reachable from the API without going through text.

---

## 2. CLI

All console output goes to **stdout**; stderr is empty in every observed case, including fatal
errors **[V]**. Exit **1** on argument/input error, **0** on success. An unknown option does
**not** abort a run **[V]**.

### 2.1 The 46 long options + 1 short option

This build defines exactly 46 long options and `-b`. Confirmed by sweeping 5 402 candidate
names harvested from the binary's strings: every accepted token is a unique case-insensitive
prefix of one of these 46, and nothing else is accepted **[V]**.

Types: `SW` switch · `INT` integer · `SINT` "smart" integer · `DBL` double · `STR` string.

| # | Group label | Option | Type | Default | Tier |
|---|---|---|---|---|---|
| 1 | Close Relative Inference | `--related` | SW | off | 1 |
| 2 | | `--duplicate` | SW | off | 1 |
| 3 | Pairwise Relatedness Inference | `--kinship` | SW | off | 1 |
| 4 | | `--ibdseg` | SW | off | 2 |
| 5 | | `--ibs` | SW | off | 1 |
| 6 | | `--makeGRM` | SW | off | 3 |
| 7 | Inference Parameter | `--degree` | SINT | 0 (unset) | 1 |
| 8 | | `--noscreen` | INT | uninitialised — see §2.7 | 1 |
| 9 | | `--seglength` | DBL | unset (NaN) → 3 Mb | 2 |
| 10 | | `--minConc` | DBL | **0.8** | 1 |
| 11 | Relationship Application | `--unrelated` | SW | off | 1 |
| 12 | | `--cluster` | SW | off | 3 |
| 13 | | `--build` | SW | off | 2 |
| 14 | QC Report | `--bysample` | SW | off | 2 |
| 15 | | `--bySNP` | SW | off | 2 |
| 16 | | `--roh` | SW | off | 3 |
| 17 | | `--autoQC` | SW | off | 2 |
| 18 | QC Parameter | `--callrateN` | DBL | unset (NaN) → 0.95 | 2 |
| 19 | | `--callrateM` | DBL | unset (NaN) → 0.95 | 2 |
| 20 | Population Structure | `--pca` | SW | off | 3 |
| 21 | | `--mds` | SW | off | 3 |
| 22 | Structure Parameter | `--projection` | SINT | 0 | 3 |
| 23 | | `--pcs` | SINT | 0 (docs contradict: 10 vs 20) | 3 |
| 24 | Quantitative Trait GWAS | `--lmm` | SW | off | 3 |
| 25 | Binary Trait GWAS | `--tdt` | SW | off | 3 |
| 26 | | `--gdt` | SW | off | 3 |
| 27 | Association Model | `--trait` | STR | `""` | 3 |
| 28 | | `--covariate` | STR | `""` | 3 |
| 29 | | `--maxP` | DBL | unset (NaN) | 3 |
| 30 | Association Method Parameter | `--invnorm` | SW | off | 3 |
| 31 | Genetic Risk Score | `--risk` | SW | off | 3 |
| 32 | | `--model` | STR | `""` | 3 |
| 33 | | `--prevalence` | DBL | unset (NaN) | 3 |
| 34 | | `--noflip` | SW | off | 3 |
| 35 | Computing Parameter | `--cpus` | SINT | 0 → half the logical cores | 1 |
| 36 | Optional Input | `--fam` | STR | `""` | 1 |
| 37 | | `--bim` | STR | `""` | 1 |
| 38 | | `--phefile` | STR | `""` | 3 |
| 39 | | `--covfile` | STR | `""` | 3 |
| 40 | | `--prunedsnp` | STR | `""` | 3 |
| 41 | | `--sexchr` | INT | **23** | 1 |
| 42 | Output | `--rplot` | SW | off | 3 |
| 43 | | `--pngplot` | SW | off | 3 |
| 44 | | `--plink` | SW | off | 3 |
| 45 | Output Parameter | `--prefix` | STR | **`king`** | 1 |
| 46 | | `--rpath` | STR | `""` | 3 |

`-b` / `-B` (case-insensitive) is the only short option. Both `-b file.bed` and `-bfile.bed`
work. `-b` as the last token is silently accepted, leaving the value empty. Every other
single-dash token yields `Command line parameter -X (#N) ignored` **[V]**.

There is **no** `--help`, `--version`, `-h` or `-v`; all four are rejected as undefined **[V]**.

### 2.2 Matching rules **[V]**

1. **Case-insensitive.** `--RELATED` ≡ `--related`.
2. **Unique-prefix matching.** `--re`→`--related`, `--k`→`--kinship`, `--ibd`→`--ibdseg`.
3. **Ambiguous prefix is an error, not a resolution:**
   `Command line parameter --r is ambiguous`. Known ambiguous 2-char prefixes: `by`, `ca`,
   `co`, `ib`, `ma`, `no`, `pc`, `pr`, `rp`, `se`. Bare `--` is ambiguous (empty prefix).
4. **Longer than the option name is undefined:** `--relatedXX`, `--pca2` → `is undefined`.
5. Neither `undefined` nor `ambiguous` aborts the run.

### 2.3 Value consumption **[V]**

| Type | Bare (no value consumed) | Consumes the next token iff | Parse |
|---|---|---|---|
| `SW` | sets ON | never | — |
| `SINT` | sets **1** | token fully matches `^[+-]?[0-9]+$` | 32-bit `int`, wraps; `007`→7 (decimal) |
| `INT` | sets **0** | same | same |
| `DBL` | leaves value **unchanged** | first char ∈ `[0-9.+-]` | `strtod`-like: `3x`→3.0, `1e3`→1000.0, `0x10`→16.0, `1,5`→1.0, `.5`→0.5, `-`→0.0. `hello`, `INF`, `NaN`, leading space → rejected |
| `STR` | leaves value unchanged | **always**, even if the next token is an option | verbatim |

Consequences to reproduce: `--cpus foo` → `--cpus [1]` **and** `foo (#2) ignored`;
`--trait --related` → `--trait [--related]` with `--related` **not** turned on and **no
warning**; a value-taking option as the last argv token consumes nothing and warns nothing.

### 2.4 Console frame grammar **[V]**

```
banner                := "KING 2.3.2 - (c) 2010-2023 Wei-Min Chen\n\n"
params_header         := "The following parameters are in effect:\n"
binary_file_line      := sprintf("%30s : %15s (-b%s)\n", "Binary File", <-b value>, "name")
blank                 := "\n"
groups_header         := "Additional Options\n"
<group lines, §2.5>
options_block_trailer := "\n"
warning_block         := "\n\x07WARNING - \nProblems encountered parsing command line:\n\n"
                         + problem_line* + "\n"
starts_line           := "KING starts at " + ctime()          // ctime() already ends in "\n"
fatal_block           := "\nFATAL ERROR - \n" + message + "\n\n"
ends_line             := "KING ends at " + ctime()
```

Byte-exact details, each confirmed by experiment:
* `FATAL ERROR - ` has a **trailing space** before its `\n`.
* The WARNING block is preceded by a literal **0x07 BEL** byte. FATAL ERROR is not.
* The no-genotype-file message's footer indents `Chen et al. 2024,` with **exactly 5 TAB
  characters**, and the last line with **10 spaces**:
  ```
  Genotype files are required. e.g.,\n  king -b ex.bed --related\n\n
  Please check the reference paper Manichaikul et al. 2010 Bioinformatics,\n
  \t\t\t\t\tChen et al. 2024,\n
            or the KING website at kingrelatedness.com\n\n
  ```
  This epilogue belongs to *that message's text*, not to the fatal frame — a missing-`.bed`
  fatal does not print it.
* Problem lines, one per problem, in command-line order:
  | Situation | Line |
  |---|---|
  | Unknown long option | `Command line parameter --bogusflag is undefined` |
  | Ambiguous prefix | `Command line parameter --r is ambiguous` |
  | Unknown short option | `Command line parameter -q (#3) ignored` |
  | Bare positional token | `Command line parameter positional (#2) ignored` |
  | Non-numeric value after a numeric option | `Command line parameter hello (#2) ignored` |

  `#N` is the **1-based index into `argv[1..]`**.

### 2.5 Parameters-in-effect rendering algorithm **[V]**

Validated against the binary on 600 randomised option sets and orderings: **600/600 exact byte
matches.**

```python
LABEL_W_TOP = 30      # "Binary File" label, right-aligned; ':' lands at 1-based column 32
VALUE_W_TOP = 15      # -b value, right-aligned, overflows right
LABEL_W     = 33      # group label, right-aligned; ':' lands at 1-based column 35
MAXCOL      = 78      # wrap budget, EXCLUDING the trailing comma

def render_group(label, items):            # items already rendered, e.g. "--minConc [0.80]"
    line = f"{label:>33} :"                # length 35, NO trailing space
    out  = []
    for k, item in enumerate(items):
        if len(line) + 1 + len(item) > MAXCOL:
            out.append(line)               # flush as-is
            line = " " * 35                # continuation stub
        line += " " + item                 # separator is always one leading space
        if k != len(items) - 1:
            line += ","                    # the comma is NOT counted in the fit test
    out.append(line)
    return out
```

* Continuation items start at 1-based column **37**, exactly under the first item of line 1.
* Because the comma is not counted, a line legitimately reaches **79** characters.
* An over-long single item is never broken: the bare `{label:>33} :` line (exactly 35 chars,
  no trailing space) is flushed and the item overflows past column 78 on the next line.
* The block re-wraps dynamically as values change.

**Value rendering:**

| Type | Rendered | Shown when |
|---|---|---|
| `SW` | `--name [ON]` / `--name` | on / off |
| `INT`, `SINT` | `--name [%d]` / `--name` | value **≠ 0** |
| `DBL` | `--name [<fmt>]` / `--name` | value is **not NaN** — an explicit `0` prints `[0.00]` |
| `STR` | `--name [%s]` | **always**; empty → `[]` |

**Double format (exact):** `(v == 0.0 || v >= 0.01) → "%.2f"`, otherwise `"%.1e"`.
There is no `fabs()`, so **every negative value uses scientific notation** (`-1` →
`-1.0e+00`). Integers are plain `%d`, 32-bit signed, wrapping.

Group order and membership are the table of §2.1, columns 2–3, in that order.

### 2.6 Default block, byte-exact

`king` with no arguments (exit 1, stdout 1611 bytes). `·` = space:

```
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen

The following parameters are in effect:
···················Binary File :·················(-bname)

Additional Options
·········Close Relative Inference : --related, --duplicate
···Pairwise Relatedness Inference : --kinship, --ibdseg, --ibs, --makeGRM
··············Inference Parameter : --degree, --noscreen [-1717986816],
····································--seglength, --minConc [0.80]
·········Relationship Application : --unrelated, --cluster, --build
························QC Report : --bysample, --bySNP, --roh, --autoQC
·····················QC Parameter : --callrateN, --callrateM
·············Population Structure : --pca, --mds
··············Structure Parameter : --projection, --pcs
··········Quantitative Trait GWAS : --lmm
················Binary Trait GWAS : --tdt, --gdt
················Association Model : --trait [], --covariate [], --maxP
·····Association Method Parameter : --invnorm
···············Genetic Risk Score : --risk, --model [], --prevalence, --noflip
··············Computing Parameter : --cpus
···················Optional Input : --fam [], --bim [], --phefile [],
····································--covfile [], --prunedsnp [],
····································--sexchr [23]
···························Output : --rplot, --pngplot, --plink
·················Output Parameter : --prefix [king], --rpath []


FATAL ERROR - 
Genotype files are required. e.g.,
··king -b ex.bed --related

Please check the reference paper Manichaikul et al. 2010 Bioinformatics,
[5 TABs]Chen et al. 2024,
··········or the KING website at kingrelatedness.com

```

An unset `-b` renders as 15 spaces, making the Binary File line exactly 57 chars.

### 2.7 The `--noscreen` display bug — required for parity **[V]**

`--noscreen` shows a garbage default `[-1717986816]`. It is not a constant: `noscreen` is an
`int` whose storage **overlaps** the `double minConc`, with `&minConc == &noscreen + 1`. Each
write clobbers the other. `noscreen`'s own byte is 0x00; `minConc` is initialised to 0.8 at
startup, whose little-endian bytes are `9A 99 99 99 99 99 E9 3F`; reading
`[0x00, 0x9A, 0x99, 0x99]` as `i32` gives `0x99999A00` = −1 717 986 816. ✔

Exact emulation (validated over 600 randomised runs including argument ordering):

```rust
let mut buf = [0u8; 16];
buf[1..9].copy_from_slice(&0.8f64.to_le_bytes());          // startup init of minConc
for (name, v) in options_in_command_line_order {
    match name {
        "noscreen" => buf[0..4].copy_from_slice(&(v as i32).to_le_bytes()),
        "minConc"  => buf[1..9].copy_from_slice(&(v as f64).to_le_bytes()),
        _ => {}
    }
}
let shown_noscreen = i32::from_le_bytes(buf[0..4].try_into().unwrap());
let shown_minconc  = f64::from_le_bytes(buf[1..9].try_into().unwrap());
```

Spot checks: `--minConc 2.5` → `0` (hidden); `--minConc 0.0001` → `474164480`;
`--noscreen 5 --minConc 0.7` → `1717986821`; `--noscreen 1 --minConc 0.9` → `-858993407`.
The reverse clobber is invisible at `%.2f`.

Put this behind the `king-bugcompat` feature. Without the feature, print `--noscreen` bare.
Note the *functional* effect of `--noscreen` on us is nil (§5.5); only the display matters.

### 2.8 Runtime console, successful run **[V]**

```
\n
KING starts at <ctime>                       // ctime() already ends in "\n"
Loading genotype data in PLINK binary format...\n
Read in PLINK fam file <path>.fam...\n
··PLINK pedigrees loaded: <n> samples\n
Read in PLINK bim file <path>.bim...\n
··Genotype data consist of <a> autosome SNPs[, <x> X-chromosome SNPs][, ...]\n
··PLINK maps loaded: <m> SNPs\n
Read in PLINK bed file <path>.bed...\n
0%\r6%\r13%\r…88%\r··PLINK binary genotypes loaded.\n
94%\r··KING format genotype data successfully converted.\n
Autosome genotypes stored in <W> words for each of <n> individuals.\n
\n
Options in effect:\n
\t--kinship\n
\t--prefix demo\n
\n
<per-analysis body, §6>
KING ends at <ctime>
```

* `<W> = ceil(m_autosome / 64)` **[V]** (48 400 autosomal SNPs → 757 words).
* The `Options in effect:` block echoes the **resolved** options one per TAB-indented line
  (`--related` rewritten to `--kinship` on small samples). Its ordering with many
  simultaneous options is **[?]** — §8 item 20.
* `<ctime>` is C `ctime()` verbatim: `Thu Aug 13 17:55:32 2026`, day-of-month space-padded,
  and it **already ends in `\n`**.
* The "…ends at" status lines are padded with **41 leading spaces** so `ends at` aligns under
  `starts at` of the preceding line.
* Degrees ≥ 3 print with the typo `3nd-degree`, `4nd-degree`, `5nd-degree`; only `1st` and
  `2nd` are correct **[V]**. Reproduce under `king-bugcompat`.
* Relationship-summary tables are TAB-separated:
  ```
  Relationship summary (total relatives: %d by pedigree, %d by inference)\n
  ··Source\tMZ\tPO\tFS\t2nd\t3rd\tOTHER\n
  ··===========================================================\n      (59 '=')
  ··Pedigree\t%d\t%d\t%d\t%d\t%d\t%d\n
  ··Inference\t%d\t%d\t%d\t%d\t%d\t%d\n\n
  ```
  The cross-family variant uses 8 spaces instead of `Source`, a 57-`=` rule, and `4th`
  instead of `OTHER`, and prints only the `Inference` row **[V]**.

---

## 3. INPUT — PLINK 1 binary fileset

Load order is strictly **`.fam` → `.bim` → `.bed`** **[V]**. `--fam <path>` and `--bim <path>`
override the derived sidecar paths. KING requires the literal `.bed` filename (`-b prefix`
without the extension fails).

### 3.1 `.fam` — samples

Six whitespace-separated columns, one line per sample, **in the same order as the `.bed`
sample bit positions**. No header. PLINK 1.9 writes it SPACE-delimited; parse as
whitespace-separated (a TAB-delimited `.fam` loads identically in KING **[V]**).

| # | Field | Notes |
|---|---|---|
| 1 | FID | see §3.5 for the `0` convention |
| 2 | IID | `(FID, IID)` is the identity key |
| 3 | PAT | `0` = unknown / not genotyped; never look it up as a sample named `"0"` |
| 4 | MAT | same |
| 5 | SEX | `1` male, `2` female, anything else unknown |
| 6 | PHENO | `1` control, `2` case, `-9`/`0`/non-numeric missing |

`n_samples` = line count. This single number determines `bytes_per_variant` for the whole
`.bed`; the `.bed` carries no counts.

**Deliberate divergence:** a 5-column `.fam` makes KING silently load **0 samples** **[V]**.
We reject it loudly (`FamShortLine{line, fields}`).

### 3.2 `.bim` — variants

Six whitespace-separated columns, one line per variant, in `.bed` row order. No header.
PLINK 1.9 writes TAB-delimited; parse as whitespace-separated.

| # | Field | Type |
|---|---|---|
| 1 | chromosome code | integer (§3.4) |
| 2 | variant ID | string, not required unique |
| 3 | genetic position | double, cM, `0` = unknown |
| 4 | base-pair coordinate | integer, 1-based |
| 5 | **A1** | PLINK 1.9 text-conversion convention: the **minor** allele; `0` if unobserved |
| 6 | **A2** | the major allele |

**Deliberate divergence:** a 5-column `.bim` makes KING silently load **0 SNPs** **[V]**.
We reject it loudly.

**Allele orientation matters for exactly one output column.** Four of the five per-pair
counts (`M_ij`, `N_Aa^i`, `N_Aa^j`, `N_HetHet`, `N_IBS0`) are invariant under swapping A1↔A2,
because that permutes codes `00 ↔ 11` and both heterozygosity and "opposite homozygotes" are
symmetric. So `Kinship`, `HetHet`, `IBS0`, `HetConc`, `IBS`, `Dist`, `HomConc`, `Concord` are
allele-order invariant. **`HomIBS0` is not** (§4.6) — its denominator counts A1-homozygotes.
Do not normalise allele order; do not read allele letters at all except for `HomIBS0`.

### 3.3 `.bed` — the genotype matrix

```
offset  size                     content
------  -----------------------  --------------------------------------
0       1                        0x6c            magic
1       1                        0x1b            magic
2       1                        0x01            mode: 1 = SNP-major
3       m * bytes_per_variant    payload, variant-major
```

```
bytes_per_variant = (n_samples + 3) / 4          // integer division
expected_len      = 3 + n_variants * bytes_per_variant
```

**2-bit code table** (the single most common source of bugs):

| bits | meaning | plane |
|---|---|---|
| `00` | homozygous **A1/A1** | `hom_a1` |
| `01` | **missing** | excluded everywhere |
| `10` | heterozygous | `het` |
| `11` | homozygous **A2/A2** | `hom_a2` |

Note `01` is *missing*, not "one copy"; the dosage order is `00 → 10 → 11`.

**Bit order (proved, not assumed).** The lowest bit pair of a byte holds the *lowest-numbered*
sample:

```
 bit:    7 6 | 5 4 | 3 2 | 1 0
       sample3 sample2 sample1 sample0     (0-based within the byte)

byte_index = 3 + variant * bytes_per_variant + (sample >> 2)
code       = (bed[byte_index] >> (2 * (sample & 3))) & 3
```

Proof: a fileset was designed so that low-bits-first predicts bytes `0xE4` then `0x1B`; PLINK
emitted exactly `6c1b 01 e4 0f 1b 0f ff 0f 55 05`, and the decode was independently
cross-checked against `plink --recode` — all match **[V]**.

**Padding hazard.** The high bits of a variant row's last byte are zero-filled, and `00` is a
*valid* code meaning hom-A1. `0x0F` (not `0x55`) in the tail byte is the direct evidence. An
unmasked popcount invents up to 3 phantom hom-A1 samples per variant. Build the bit planes by
setting bits only for real samples/variants and the invariant holds by construction.

### 3.4 Chromosome partition and `--sexchr` **[V]**

PLINK normalises names on `--make-bed`: `X`→23, `Y`→24, `XY`(PAR)→25, `MT`→26,
unplaced→0.

**With the default `--sexchr 23`, KING's autosome set is `1..=22` PLUS `25` (XY/PAR).**
This is the easy rule to miss. KING reports it as
`Genotype data consist of 100 autosome SNPs (including 20 XY SNPs), 40 X-chromosome SNPs, …`
and `20 other SNPs are removed.`

General rule, fully determined by sweeping `N`:

```
autosomes = chromosomes 1 .. N-1,  plus chromosome N+2
X  = N        Y  = N+1        XY = N+2 (folded into autosomes)        MT = N+3
everything else, including chromosome 0, is REMOVED
```

`--sexchr != 23` additionally prints
`Non-human samples are analyzed, with %d pairs of chromosomes`. Zero-count classes are elided
from the SNP-count line.

**Sort order.** `--kinship`/`--duplicate` counting is order-independent. `--ibs` continues
its counting pass but disables its segment columns, while every other segment consumer
falls back or stops segment work. The diagnostic is
`Chromosomes unsorted: %s on chr %d, %s on chr %d.` / `Positions unsorted: %s at %d, %s at %d.`
So the relatedness counts need no sort; segment work requires ascending `(chr, bp)` **[V]**.

### 3.5 Identity, duplicates, and the FID `0` convention **[V]**

* The identity key is the **pair** `(FID, IID)`. The same IID under different FIDs is fine.
* An exact `(FID, IID)` duplicate is fatal:
  `Family F1: Person I1 is duplicated` then
  `FATAL ERROR - \nPlease correct problems with pedigree structure`.
  Comparison is **case-insensitive** — `a` and `A` in one family collide.
* **FID `0` does not mean "each sample is its own family."** KING pools them:
  `All individuals with family ID 0 are considered as relatives.` A `.gq` bundle whose `.fam`
  uses `0` throughout is analysed as a **single family**, so everything lands in `.kin` and
  `.kin0` is empty. Surface this in the GeneQuire import path.

### 3.6 Error messages (reproduce verbatim) **[V]**

| Condition | Message |
|---|---|
| Bad magic / zero-length `.bed` | `Please use either PLINK or KING binary format as input.` |
| Mode byte `0x00` (individual-major) | `Currently only SNP-major mode can be analyzed.` |
| `.bed` truncated | `Not enough genotypes at the <k>th marker` — `k` is **0-based**, suffix is a literal `th` (`0th`, `1th`, `2th`) |
| `.bed` unopenable | `Genotype file <path> cannot be opened` |
| `.bim` absent | `Map file <path>.bim cannot be opened` |
| `.fam` absent | `Pedigree file <path>.fam cannot be opened` |
| No autosomes | `No autosome SNPs are available. Please check your map file.` |

All fatals print as `\nFATAL ERROR - \n<message>\n\n` and exit 1.

### 3.7 Length validation — a deliberate safety divergence

KING never validates file length. A **short `.fam`** or **short `.bim`** produces
**exit 0 with no warning** and misaligned/dropped data — silent corruption **[V]**. A long
`.bim` is fatal.

We check up front:

```
expected = 3 + n_variants * bytes_per_variant(n_samples)
if actual_len != expected -> Err(LengthMismatch { n_samples, n_variants, expected, actual })
```

One `stat`, catches every mismatch class, converts a silent wrong answer into a diagnosable
one. On any well-formed fileset the results are identical, so this costs no parity.

`mmap` the `.bed`: read-only, read-once, often multi-GB.

---

## 4. ALGORITHMS

### 4.1 Symbols

Fix a pair `(i, j)` and the set of variants that pass the filter of §5.

| Symbol | Definition |
|---|---|
| `M_ij` | # variants with a **non-missing genotype in both** `i` and `j`. This is the `N_SNP` column. |
| `N_i` (`N_Aa^i`) | # variants where `i` is heterozygous, **restricted to the `M_ij` set** |
| `N_j` (`N_Aa^j`) | symmetric |
| `N_HetHet` | # variants where **both** are heterozygous |
| `N_IBS0` (`N_AA,aa`) | # variants where the two are **opposite homozygotes** |
| `N_HomHom` | # variants where both are homozygous (any combination); `N_IBS0 ⊆ N_HomHom` |
| `N_A1any` | # variants in the `M_ij` set where **at least one** of the pair is hom-A1 |

**The parity-critical rule, stated three ways because it is the most-botched detail:**
`N_i` and `N_j` are **recomputed per pair** over the pairwise-complete set. They are *not*
each sample's global heterozygote count. The paper states it verbatim ("excluding those SNPs
with missing genotypes in either individual of the pair"); Hail, SNPRelate, plink2 and akt all
do it; getting it wrong yields subtly wrong kinship for every pair with differential
missingness and passes every equal-call-rate test.

### 4.2 IBS bookkeeping (exact, not approximate) **[P][V]**

| geno i | geno j | IBS | `(x_i−x_j)²` |
|---|---|---|---|
| AA/AA, aa/aa, Aa/Aa | | 2 | 0 |
| any one het, other hom | | 1 | 1 |
| **AA/aa or aa/AA** | | **0** | **4** |

```
N_IBS1 = N_i + N_j - 2*N_HetHet          // exactly one of the pair is het
N_IBS2 = M_ij - N_IBS0 - N_IBS1
D_ij   = Σ_m (x_i - x_j)²  =  N_IBS1 + 4*N_IBS0  =  N_i + N_j - 2*N_HetHet + 4*N_IBS0
```

`D_ij` is the master quantity: **every KING-robust estimator is `1/2 − D_ij/(2·Denom)`.**
The identity `N_IBS1 = N_i + N_j − 2·N_HetHet` and `N_IBS2 = M_ij − N_IBS0 − N_IBS1` were
confirmed against KING's own `.ibs` integer columns with **0 mismatches** over all pairs of
two independent fixtures **[V][V-NEW]** — they are KING's definitions, not merely ours.

### 4.3 The two estimators

**Within-family — paper Eq (9). Written to `.kin`. [P][V-NEW]**

```
              N_HetHet - 2*N_IBS0                              D_ij
phi_within =  -------------------      ==      1/2  -  ----------------------
                  N_i + N_j                              2 * (N_i + N_j)
```

**Between-family — paper Eq (11). Written to `.kin0`. [P][V-NEW]**

```
                        4*N_IBS0 + (N_i - N_HetHet) + (N_j - N_HetHet)
phi_between = 1/2  -   -----------------------------------------------
                                4 * min(N_i, N_j)

            == 1/2  -  D_ij / (4 * min(N_i, N_j))

            == 1/2  +  (2*N_HetHet - 4*N_IBS0 - N_i - N_j) / (4 * min(N_i, N_j))
```

**Unified:** `phi = 1/2 − D_ij / (2·Denom)`, with `Denom = N_i + N_j` (within) or
`Denom = 2·min(N_i, N_j)` (between). They coincide exactly when `N_i == N_j`.

Both forms were confirmed against the reference binary this session: every `Kinship` value in
`.kin`, `.kin0`, `.ibs` and `.ibs0` on a 10-sample × 48 400-SNP fixture reproduced to the
printed `%.4f`, **0 mismatches** **[V-NEW]**.

**The `+1/2 − (N_i+N_j)/(4·min)` correction in Eq (11) is the single most commonly dropped
piece of this formula.** It vanishes only when `N_i == N_j`, so a naive
`(N_HetHet − 2·N_IBS0)/(2·min)` passes every MZ-twin and equal-heterozygosity test and is then
silently wrong on real pairs — golden vector B in §A.1 returns **0.5000** instead of **0.1250**,
i.e. it reports an unrelated pair as a duplicate.

**Which estimator for which pair (unanimous across the paper, SNPRelate, and the binary):**

| Pair | Estimator | File |
|---|---|---|
| same FID | Eq (9), **sum** denominator | `.kin`, `.ibs` |
| different FID | Eq (11), **min** denominator | `.kin0`, `.ibs0` |

plink2 and Hail implement Eq (11) *only* ("Pedigree information is currently ignored" /
"only implements the 'between-family' estimator"). If we validate `--kinship` against plink2
alone, every within-family pair will appear to disagree and we will fix the wrong thing.

**Rules that follow:**
* **Accumulate integers; divide once.** Both estimators are an exact-integer numerator over an
  exact-integer denominator with a single f64 division. Results are then bit-reproducible
  regardless of SIMD or threading. Never accumulate `D_ij` as a float.
* **Never clamp.** Negative kinship is meaningful — an extreme negative value signals the pair
  is drawn from two distinct populations (Eq 10). KING's manual says so explicitly. The 0.5
  upper bound is analytic, not enforced.
* Degenerate denominators: §4.9.

### 4.4 Bit-plane popcount engine

KING packs autosomal genotypes into **64-bit words**, proven by the 64→1 / 65→2 step in
`Autosome genotypes stored in %d words for each of %d individuals` **[V]**. We adopt the same
strategy independently — the paper itself describes it ("When each genotype is stored in two
bits, … can be computed using only bit operations"), and popcount genotype counting is
standard practice.

Per sample, `W = ceil(m_autosome / 64)` words in **each of four planes**. Bit `b` of word `w`
is autosomal variant `64*w + b` in `.bim` order restricted to the autosome set.

| Plane | Set when |
|---|---|
| `hom_a1[i]` | code `00` |
| `het[i]` | code `10` |
| `hom_a2[i]` | code `11` |
| `nonmiss[i]` | code `!= 01`, i.e. `hom_a1 \| het \| hom_a2` (materialised once; used in 4 of 6 expressions) |

Cost: 4 bits/variant/sample. Tail bits above `m_autosome % 64` must be zero in every plane;
build by setting only real variants, then assert once:

```rust
let tail_mask = if m % 64 == 0 { !0u64 } else { (1u64 << (m % 64)) - 1 };
debug_assert_eq!(nonmiss[i][W-1] & !tail_mask, 0);
```

The six counts:

```rust
for w in 0..W {
    let (nm_i, nm_j) = (nonmiss[i][w], nonmiss[j][w]);
    let (he_i, he_j) = (het[i][w],     het[j][w]);
    let (a1_i, a1_j) = (hom_a1[i][w],  hom_a1[j][w]);
    let (a2_i, a2_j) = (hom_a2[i][w],  hom_a2[j][w]);

    m_ij     += (nm_i & nm_j).count_ones() as u64;                       // M_ij
    n_i      += (he_i & nm_j).count_ones() as u64;                       // N_Aa^i  <- the `& nm_j` is load-bearing
    n_j      += (he_j & nm_i).count_ones() as u64;                       // N_Aa^j
    n_hethet += (he_i & he_j).count_ones() as u64;                       // N_HetHet
    n_ibs0   += ((a1_i & a2_j) | (a2_i & a1_j)).count_ones() as u64;     // N_IBS0
    n_homhom += (((a1_i|a2_i) & (a1_j|a2_j)) & nm_i & nm_j).count_ones() as u64;
    n_a1any  += (((a1_i | a1_j) & nm_i & nm_j)).count_ones() as u64;     // N_A1any, for HomIBS0
}
```

`het_i` already implies `nonmiss_i`, so `& nm_i` is redundant on the `n_i` line; the `& nm_j`
is **not** — dropping it yields the marginal het count and is exactly the failure mode of §4.1.

These expressions were validated against KING's own `--ibs` integer columns: all 28 pairs × 8
integer columns on an 8-sample × 4 000-variant fixture with 8 % missingness, **0 mismatches**
**[V]**; and re-validated this session on a 10-sample × 48 400-variant fixture **[V-NEW]**.

Threading: pairs are independent. Use `rayon` over the upper-triangular pair list. Because
every accumulator is an exact integer, thread count cannot change any output byte (unlike
KING, whose X-chromosome writer has a data race).

### 4.5 Reported rates on `.kin` / `.kin0` **[V-NEW]**

```
N_SNP  = M_ij                            %d
HetHet = N_HetHet / M_ij                 %.4lf      <- PROPORTION, not a count
IBS0   = N_IBS0   / M_ij                 %.4lf      <- PROPORTION, not a count
```

> **Resolved conflict.** The website's `--kinship` example prints `HetHet` at **3** decimals
> and its `--related` example at **4**, and recon flagged this as two different writers. It is
> a **32-bit vs 64-bit** difference: `ComputeShortRobustKinship()` (32-bit path) uses `%.3lf`
> **[S]**; the 64-bit path — the only one that runs on arm64/x86-64 — uses `%.4lf` in **both**
> writers. Confirmed by hexdumping a live `--kinship` run: `FAM1\tFA1\tKID1\t48400\t0.000\t0.2500\t0.2007\t…`
> **[V-NEW]**. **Always emit `%.4lf`.**

### 4.6 Derived columns — all definitions confirmed empirically **[V-NEW]**

Every formula below was fitted and then verified against the reference binary's own output:
254 field comparisons across `.kin`, `.kin0`, `.ibs`, `.ibs0` and 315 across `.con`,
**0 mismatches**. Reproduction in §A.3.

```
HetConc  = N_HetHet / (N_i + N_j - N_HetHet)          // symmetric union / Jaccard form
Het2|1   = N_HetHet / N_i                              // "het in 2 given het in 1"
Het1|2   = N_HetHet / N_j
HomConc  = (N_HomHom - N_IBS0) / N_HomHom
HomIBS0  = N_IBS0 / N_A1any                            // see the warning below
IBS      = (N_IBS1 + 2*N_IBS2) / M_ij                  // mean IBS allele sharing, in [0,2]
Dist     = (N_IBS1 + 4*N_IBS0) / M_ij  ==  D_ij / M_ij // mean squared dosage difference
Concord  = N_IBS2 / M_ij                               // .con only
```

> **`Dist` is NOT `2 − IBS`.** That is the obvious guess and it is wrong: on a pair with
> `M=48400, N_IBS0=3927, N_IBS1=22139, N_IBS2=22334`, KING prints `IBS 1.3803` and
> `Dist 0.7820`, whereas `2 − IBS = 0.6197`. `Dist` is the mean of `(x_i − x_j)²`, which
> weights IBS0 by 4 rather than 2. It coincides with `2 − IBS` only when `N_IBS0 = 0`, which
> is why a PO-only fixture will not catch the error. **[V-NEW]**

> **`HomIBS0` is the one allele-order-dependent number in the whole relatedness surface.**
> Its denominator `N_A1any` counts pairwise-called variants where **at least one** individual
> is homozygous for **A1** (the minor allele under the PLINK 1.9 convention) — i.e. the
> "informative" sites of KING's QC vocabulary ("at least one carries the minor homozygote").
> Fitted denominators matched exactly: 9144/9144, 9203/9203, and 7636 / 8743 within the 4-dp
> rounding bracket. Candidate definitions that are **wrong**: `N_IBS0/N_HomHom` (that is
> `1 − HomConc`), `N_IBS0/M_ij` (that is `IBS0`). This is also why KING warns
> `Too many first alleles as the major allele (~%.1lf%%). Please use plink1.9 --make-bed …`.
> **[V-NEW]**, §8 item 12.

### 4.7 Pedigree-derived `Z0` and `Phi`

`Z0` and `Phi` are **expectations from the `.fam` pedigree**, not estimates **[D][V]**:

| pedigree relationship | `Z0` (`%.3lf`) | `Phi` (`%.4lf`) |
|---|---|---|
| parent–offspring | `0.000` | `0.2500` |
| full siblings | `0.250` | `0.2500` |
| grandparent / avuncular / half-sib | `0.500` | `0.1250` |
| unrelated (including spouses) | `1.000` | `0.0000` |

General computation for non-inbred pedigrees (founders with unknown parents are mutually
unrelated):

```
phi(a,a) = 1/2 * (1 + phi(father_a, mother_a))
phi(a,b) = 1/2 * (phi(father_a, b) + phi(mother_a, b))       // a not an ancestor of b
k2(a,b)  = phi(Fa,Fb)*phi(Ma,Mb) + phi(Fa,Mb)*phi(Ma,Fb)     // both parents known for both
k1       = 4*phi(a,b) - 2*k2
Z0 = k0  = 1 - k1 - k2
```

Memoise `phi` over the topologically sorted pedigree. Behaviour on inbred loops and on
FID-`0`-pooled samples is **[?]** — §8 item 13.

### 4.8 The `Error` flag **[V-NEW, fitted]**

`Error` (`.kin` only; never in `.kin0`) flags disagreement between the estimated kinship and
the pedigree expectation. Documented values: **1 = error, 0.5 = warning**, `0` otherwise
**[D]**. The predicate is not documented. Fitted from 247 rows spanning `Phi ∈ {0, 0.125,
0.25}`, every row consistent:

```
if Phi > 0:
    r = Kinship / Phi
    Error = 0     if  2^-0.5 < r <= 2^0.5      // within a factor of sqrt(2)
    Error = 0.5   if  2^-1   < r <= 2^1        // within a factor of 2
    Error = 1     otherwise
else:  // Phi == 0, pedigree says unrelated
    Error = 0     if  Kinship <= 2^-5.5        // 0.0220970869
    Error = 0.5   if  Kinship <= 2^-4.5        // 0.0441941738
    Error = 1     otherwise
```

Observed brackets that pin it: at `Phi=0.25`, Error 1 up to `0.1235` and Error 0.5 from
`0.1285` (boundary `Phi/2 = 0.125`); Error 0 from `0.1816` to `0.3463`, Error 0.5 at `0.1638`
and `0.3680` (boundaries `0.17678 = 2^-2.5` and `0.35355 = 2^-1.5`); at `Phi=0.125`, a
`Kinship 0.258` row is Error 1 (ratio 2.064 > 2 ✓). At `Phi=0`, Error 0 max `0.0186`,
Error 0.5 spans `0.0248…0.0361`, Error 1 from `0.0482`.

Note that for `Phi` an exact power of two the inner band **coincides** with the standard degree
class, but the outer band is **not** the adjacent classes — it is `[Phi/2, 2·Phi]`, i.e. only
the upper half of the next class down. A class-distance rule reproduces `Phi=0.25` rows
incorrectly. Boundary inclusivity is unverified — §8 item 4.

### 4.9 Degenerate denominators **[V-NEW]**

Fitted against purpose-built fixtures (an all-homozygous sample; a pair with disjoint call sets):

| Case | KING's behaviour |
|---|---|
| `M_ij == 0` (no shared called variants) | The pair is **omitted entirely** from `.kin` / `.kin0` / `.ibs` |
| within-family, `N_i + N_j == 0` | The row is **omitted from `.kin`** (verified: an all-hom pair in one family produced no row, while an all-hom sample paired with a normal sample — sum > 0, min = 0 — *did* produce a row with `Kinship -0.6647`) |
| between-family, `min(N_i, N_j) == 0` | The row **is written**, with `Kinship` printed as `0.0000` (not `-inf`, not `nan`) |

For reference, plink2 emits `-inf` here and SNPRelate emits `NaN`; **KING emits `0.0000`**, so
neither third-party convention is the parity target. Implement as an explicit guard:
`if denom == 0 { 0.0 } else { … }` for the between case, and a skip for the within case.
Whether the within-family skip is keyed on `N_i+N_j == 0` or on something else is §8 item 5.

### 4.10 `--duplicate`

Compute the same per-pair counts, then:

```
Concord = N_IBS2 / M_ij                                     %.5lf
HomConc = (N_HomHom - N_IBS0) / N_HomHom                    %.5lf
HetConc = N_HetHet / (N_i + N_j - N_HetHet)                 %.5lf
```

**Selection:** a pair is written iff `HetConc > minConc` — strict `>`, default `0.80` **[V]**.
Both within- and cross-family pairs go into the single `.con` file, in `.fam` serial order.
All 45 pairs of a 10-sample fixture at `--minConc 0.01` reproduced with **0 mismatches**
across `N`, `N_IBS0/1/2`, `Concord`, `HomConc`, `HetConc` **[V-NEW]**.

Console: `%lli pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s`
where `%d` is `round(minConc * 100)` rendered without decimals (`80`, `1`). If nothing passes:
`No duplicates are found with heterozygote concordance rate > %d%%.` and the `.con` file is
still written **with its header line only** **[V]**.

### 4.11 `--unrelated`

Extracts a maximal subset containing no pair at or above the `--degree` threshold, writing the
kept list to `<p>unrelated.txt` and the complement to `<p>unrelated_toberemoved.txt`. The two
files partition the sample set (5 + 3 = 8; 116 + 48 = 164) **[V]**. Row order is neither `.fam`
order nor sorted order — it is the greedy visit order, and it is reproducible run to run
**[V]**. The algorithm is described in Manichaikul et al. 2012 (PLoS Genet 8:e1002640), which
is a permitted input. **[?]** — §8 item 15.

Gate: `N < 10` prints
`This function is currently disabled for tiny dataset with sample size < 10.` and both files
are still written, without family clustering **[V]**.

### 4.12 IBD segments (Tier 2 — sketch + the honest gap)

Definitions **[D]**, all denominators over the **autosomal** usable length:

```
D        = Σ Length over all rows of <p>allsegs.txt        (autosomes only)
IBD1Seg  = (total Mb called IBD1) / D
IBD2Seg  = (total Mb called IBD2) / D
PropIBD  = IBD2Seg + IBD1Seg / 2                            ( = 2 * phi )
```

`PropIBD` is computed in **f64 from full-precision** π1/π2 and formatted once at 4 dp — never
from the printed values **[V]**. Confirmed in the binary by an `fmov d8,#0.5` + `fmadd`
**[S]**.

Fixed filters, announced verbatim on every run **[V]**:
* a pair is emitted iff the conditioned, floor-dependent merged IBD1 or IBD2 call set has
  at least one segment **≥ 10 Mb** (not tunable). A held-out IBD1 pair and an independent
  IBD2 canvas distinguish this from filtering the unmerged calls
  (`tests/parity/probes/segment_residuals.py`);
* segments shorter than `--seglength` (default **3 Mb**, valid 1–10 Mb, units Mb stored as bp;
  out-of-range values silently revert to 3 Mb — `--seglength 0` and `11` gave md5-identical
  output) are discarded.

Usable-segment construction for `allsegs.txt` **[V, partial]**: cut each chromosome wherever
consecutive markers are **more than 1 000 000 bp apart** (exactly 1 000 000 does *not* cut —
verified to the bp), then drop pieces that are too small. The drop rule brackets to
"span > 10 Mb **and** ≥ 5 complete 64-marker words" but one counter-case (11 Mb / 352 markers)
remains unexplained — §8 item 17. This rule scales **every** IBD proportion, so it is a
first-order parity item.

**The IBD1/IBD2 calling rule itself is undocumented and unpublished.** The manual states the
manuscript is "yet to be published"; the binary cites "Chen et al. 2024", which does not exist
in any index (not on kingrelatedness.com, not in Wei-Min Chen's publication list, not in PubMed
or Crossref). Behavioural parity is the only spec. Fit it against the reference binary on
synthetic pairs — never from KING source. Observations to fit against: within a usable segment,
IBD1 stretches behave like runs with no IBS0 (an IBD1 region cannot produce IBS0), IBD2 like
runs with neither IBS0 nor het/hom mismatch, with error tolerance and word-granular (64-marker)
scanning; a truly IBD1 block of *L* Mb is recovered as ~*L*+1–2 Mb because boundaries extend to
the next IBS0. **Treat that as the acceptance test, not as the algorithm.** §8 item 16.

---

## 5. SNP FILTERING — exactly which variants enter each computation

### 5.1 The rule

For **all** Tier-1 relatedness computations (`--related`, `--kinship`, `--duplicate`, `--ibs`,
`--unrelated`) the variant set is:

```
all variants whose .bim chromosome code is in the AUTOSOME set of §3.4
   ( = 1..=22 and 25, at the default --sexchr 23 )
```

and nothing else. In particular:

* **No MAF filter.** **No call-rate filter.** **No LD pruning.** **[V-NEW]**
* **Monomorphic variants are kept and counted.** Verified: a 4 000-variant fileset containing
  506 monomorphic variants reproduced KING's `--ibs` counts exactly only when all 4 000 were
  used **[V]**. Do not filter them.
* **All-missing variants are kept** in the variant list; they simply never contribute, since
  every count is masked by `nonmiss`.
* Missingness is handled **pairwise**, per §4.1 — never by dropping a variant globally.
* Chromosome `0` and anything outside the partition are dropped ("other SNPs are removed").

> **Resolved risk.** Recon flagged the binary string
> `"%d autosome SNPs with MAF>%.3lf and call rate>%d%% are used."` as evidence that KING
> pre-filters variants by MAF and call rate by default, and ranked it the largest threat to
> `N_SNP` parity. Cross-referencing the string against its owning function shows it belongs to
> **`Engine::MakeGRM0_LMM(IntArray&, Matrix&)`** — the GRM/LMM path (`--makeGRM`, Tier 3)
> **[V-NEW]**. It is not on any relatedness path. Corroborated behaviourally: `--kinship` on
> a fixture with 48 400 autosomal variants reports `N_SNP = 48400` for every fully-called pair,
> i.e. **no variant was filtered** **[V-NEW]**. Confidence: high. Residual verification in
> §8 item 6.

The website's own guidance agrees: *"Please do not prune or filter any 'good' SNPs that pass QC
prior to any KING inference … LD pruning is not recommended in KING."* **[D]**

### 5.2 Sample-level exclusion: the `M ≤ 512` screen **[V]**

Before the **cross-family** stage, KING drops samples with too few called autosomal variants:

```
The following %d samples are excluded from the kinship analysis (M<512):\n
\t(%s %s)\t(%s %s)…\n
```

Measured with a per-sample call-count ladder: samples with **M ≥ 513** appear in `.kin0`;
samples with **M ≤ 512** do not. The real predicate is therefore `M < 513`, despite the
printed `M<512`. Reproduced this session on a fixture where one sample had 10 called variants
**[V-NEW]**.

**The printed name list is wrong.** The *count* is correct; the names are simply the first
`count` samples in `.fam` serial order, including samples demonstrably present in `.kin0`.
This is a KING display bug — emit the correct names (do not bug-match), or gate the correct
list behind `king-bugcompat`.

Excluded samples still participate in the **within-family** stage (`.kin`) **[V]**.

Unresolved wrinkle: a separate sweep found the boundary tracking the *total* autosomal variant
count (`m = 544` excluded, `m = 545` not), independent of `n`, MAF and monomorphic fraction —
so the `m → M` mapping is not simply `M = m`. It matters only for tiny filesets; every
realistic `.gq` bundle is far past it. **Do not encode a guessed rule** — §8 item 7.

### 5.3 Pairs that are skipped entirely

* `M_ij == 0` → the pair is omitted from every output file **[V]**.
* Within-family, `N_i + N_j == 0` → omitted from `.kin` **[V-NEW]**.
* A sample with no non-missing genotypes never appears in any pair **[V]**.
* A family with a single genotyped member contributes no `.kin` rows
  (`Each family consists of one individual.`) **[V]**.
* If the entire fileset is one family, KING prints `There is only one family.`, writes
  `<p>.kin` as a **0-byte file** (not even a header) and writes no `<p>.kin0` **[V]**.
  Decide deliberately whether to bug-match — §8 item 10.

### 5.4 What `--degree` filters (it is *only* a write filter for Tier 1)

| Analysis | Effect of `--degree d` | Threshold |
|---|---|---|
| `--kinship` | filters `.kin0` rows; `.kin` untouched. **Without `--degree`, all cross-family pairs are written** | `Kinship >= 2^-(d+1.5)` |
| `--related` | filters `.kin0` rows; **default is `d = 1`** so `.kin0` is *always* filtered | `Kinship >= 2^-(d+1.5)` |
| `--ibdseg` (Tier 2) | filters `.seg` rows; without `--degree`, all surviving pairs are written including `UN` | `PropIBD > 2^-(d+0.5)` |
| `--unrelated`, `--build`, `--cluster` | depth of the relatedness search | as above |

Thresholds, printed to stdout at `%.5f` **[V]**:

| `d` | 1 | 2 | 3 | 4 | 5 |
|---|---|---|---|---|---|
| `2^-(d+1.5)` kinship | 0.17678 | 0.08839 | 0.04419 | 0.02210 | 0.01105 |
| `2^-(d+0.5)` PropIBD | 0.35355 | 0.17678 | 0.08839 | 0.04419 | 0.02210 |

`--degree 0` is treated as `1` **[V]**. Row counts on a 164-sample fixture: 7 / 12 / 59 / 965
/ 2455 for d = 1…5. Note that the stdout count of "identified relatives" is **smaller** than
the number of rows written (13 identified vs 59 rows at degree 3), because rows above the raw
kinship threshold that classify as `UN` are still written.

Inclusivity (`>=` vs `>`) is taken from the console message
`Between-family relatives (kinship >= %.5lf) saved in file %s` **[S]** but is not
byte-verified — §8 item 3.

### 5.5 Two-stage screening is a performance optimisation, not a filter

KING screens on a SNP subset then confirms on all SNPs
(`A subset of informative SNPs will be used to screen close relatives.`,
`Stage 1 (with %d SNPs) screening ends at %s`, `Stage 2 (with all SNPs) inference ends at %s`).
`--noscreen` disables it. **We compute all pairs directly and never screen.** This must not
change any reported number; verify by diffing a default run against `--noscreen` on the same
data — §8 item 8. Accept `--noscreen` on the command line and make it a no-op.

### 5.6 A hard invariant worth testing

The manual guarantees that `--kinship --projection N` produces estimates *"identical (with no
numerical differences)"* to an unsplit `--kinship` run **[D]**. That is only possible if the
estimator is exactly pairwise, using no dataset-wide quantity — no allele frequencies, no
global heterozygosity. Encode it as a property test: **splitting the sample set must not change
a single digit of any pair's `Kinship`, `HetHet`, `IBS0` or `N_SNP`.**

---

## 6. OUTPUT

Universal rules **[V]**:
* LF line endings only. Every file ends with exactly one `\n`; no trailing blank line, no BOM,
  no trailing field separator — except `<p>X.seg`, whose rows genuinely end in a tab
  (§1.3, `crate::analysis::xseg`).
* **TAB** separators in `.kin`, `.kin0`, `.con`, `.ibs`, `.ibs0`, `allsegs.txt`,
  `unrelated*.txt`, `.seg`. **SPACE** separators in `bySample.txt`, `bySNP.txt`,
  `splitped.txt`, `pc.txt`. Never mix.
* No field padding or alignment anywhere. Negative values print with a leading `-` and no
  padding (`-0.0059`).
* Filenames are `<p>` + suffix with **no separator**. The dot belongs to the suffix for
  `.kin/.kin0/.con/.ibs/.ibs0/.seg`; there is **no** dot before `allsegs.txt`,
  `unrelated.txt`, `unrelated_toberemoved.txt`, `bySample.txt`, `bySNP.txt`.
  `--prefix sub/run1.` → `sub/run1..kin`; the directory must already exist.

### 6.1 Row ordering — two different orders, both mandatory **[V][V-NEW]**

| Order | Files | Rule |
|---|---|---|
| **Sorted** | `.kin`, `.ibs` (and `X.kin`) | families in sorted order; members sorted within each family by the comparator below; then `i < j` over the sorted member list |
| **Serial** | `.kin0`, `.con`, `.ibs0`, `.seg`, `bySample.txt` (and `X.kin0`) | `.fam` file order; `i < j` over the global sample index |
| **Mirrored** | `X.seg` | exactly the rows of `<p>.seg`, in exactly its order (which is the 16-sample tiling, not the serial order) |

These genuinely differ. Confirmed this session on a fixture whose `.fam` order is
`FA1, MO1, KID1, KID2`: `.kin` emits `FA1–KID1, FA1–KID2, FA1–MO1, KID1–KID2, KID1–MO1,
KID2–MO1` (sorted), while `.kin0` and `.con` emit `.fam` order (`FA1–MO1` first) **[V-NEW]**.

**The comparator** (reverse-engineered from a 40-value and an independent 36-value probe that
produced the same total order, so one comparator serves both FIDs and IIDs) **[V]**:

1. **Case-insensitive, folding to UPPERCASE.** (`zA` < `z_` because `'A'`=65 < `'_'`=95;
   `zzz…` < `z_` because `'z'`→`'Z'`=90 < 95.) Two ids differing only in case are *equal* —
   KING rejects such a fileset as a duplicate person.
2. **Digits sort after every non-digit.** `ABD` < `A0`; `zA` < `z0`; `_a` < `0`. Otherwise
   ordinary ASCII on the uppercased byte: `m!`(0x21) < `m#` < `mZ` < `m}` < `m~` < `m0`.
3. **Runs of digits compare as a unit: longer run first, then digit by digit.**
   `x1y < x1z < x2y < x01y`; `b1 < b2 < b9 < b10`; `0 < 0a < 1a < 2 < 7 < 9 < 00 < 09 < 10 < 007`.
   Leading zeros are **not** normalised.
4. **End of string sorts before any character:** `a < ab < abc`.

Whether the comparison is byte-wise or locale-aware, and the FID tie-break, are §8 item 9.

### 6.2 `<p>.kin` from `--kinship` — 10 columns

**Header (literal):**
```
FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError\n
```
**Row:**
```
%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n
```
| # | Column | Source |
|---|---|---|
| 1 | `FID` | family id |
| 2–3 | `ID1`, `ID2` | pair members, in sorted order (§6.1) |
| 4 | `N_SNP` | `M_ij` |
| 5 | `Z0` | pedigree-expected Pr(IBD=0), `%.3lf` |
| 6 | `Phi` | pedigree-expected kinship, `%.4lf` |
| 7 | `HetHet` | `N_HetHet / M_ij` |
| 8 | `IBS0` | `N_IBS0 / M_ij` |
| 9 | `Kinship` | **Eq (9)**, within-family |
| 10 | `Error` | `%G`; `0`, `0.5` or `1` (§4.8) |

Contains **all** within-family pairs (no `--degree` filter). Live sample:
```
FID→ID1→ID2→N_SNP→Z0→Phi→HetHet→IBS0→Kinship→Error
FAM1→FA1→KID1→48400→0.000→0.2500→0.2007→0.0000→0.2526→0
FAM1→FA1→MO1→48400→1.000→0.0000→0.1679→0.0808→0.0079→0
```
(`→` = TAB.)

Console: `Within-family kinship data saved in file %s\n`.

### 6.3 `<p>.kin0` from `--kinship` — 8 columns

**Header:**
```
FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship\n
```
**Row:**
```
%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\n
```
`Kinship` is **Eq (11)**, between-family. **No `Error` column.** Without `--degree` this file
contains every cross-family pair (13 120 rows = C(164,2) − 246 on a 164-sample fixture)
**[V]**.

Console: `Between-family kinship data saved in file %s\n`, or with a filter
`Between-family kinship data (up to degree %d, %lli pairs in total) saved in file %s\n`, then
`Note --kinship --degree <n> can filter & speed up the kinship computing.`

### 6.4 `<p>.kin` from `--related` — 16 columns

**Header** — assembled from three pieces; the middle piece appears only when IBD segments are
available, and the tail is an **inline 8-byte immediate** invisible to `strings` **[S]**:
```
FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship   (60 bytes)
  + \tIBD1Seg\tIBD2Seg\tPropIBD\tInfType                                 (32 bytes, conditional)
  + \tError\n
```
**Row:**
```
%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf
  then  \t%.4lf\t%.4lf\t%.4lf\t%s\t%G\n      (segment path)
  or    \t%G\n                                (no-segment path — 12 columns)
```
Live sample:
```
FID→ID1→ID2→N_SNP→Z0→Phi→HetHet→IBS0→HetConc→HomIBS0→Kinship→IBD1Seg→IBD2Seg→PropIBD→InfType→Error
FAM1→KID1→KID2→48400→0.250→0.2500→0.2385→0.0202→0.4314→0.1283→0.2502→0.5364→0.2598→0.5280→FS→0
```
`Kinship` is still **Eq (9)** (verified **[V-NEW]**). `HetConc` and `HomIBS0` per §4.6.

### 6.5 `<p>.kin0` from `--related` — 14 columns

**Header** — parallel construction, but the tail immediate is a bare `\n`: **there is no
`Error` column** **[S][V]**.
```
FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship   (59 bytes)
  + \tIBD1Seg\tIBD2Seg\tPropIBD\tInfType                               (32 bytes, conditional)
  + \n
```
**Row:**
```
%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf
  then  \t%.4lf\t%.4lf\t%.4lf\t<InfType>\n
  or    \n                                    (no-segment path — 10 columns)
```
Filtered by `--degree` (default 1). **If no cross-family pair clears the threshold the file is
not created at all**, and stdout says `No close relatives are inferred.` **[V]** — note this
differs from `--kinship`, which always writes the file, and from `.con`, which writes a
header-only file.

Console: `\nBetween-family relatives (kinship >= %.5lf) saved in file %s\n`.

### 6.6 The `--related` → `--kinship` downgrade **[V]**

With **N < 10 samples** KING prints
`\n--related is replaced with --kinship for a small sample size.` and from then on behaves
exactly as `--kinship`: 10-column `.kin`, 8-column `.kin0`, `--degree` ignored. Measured
boundary: N = 8, 9 downgrade; N = 10 does not.

Related messages in the same family **[S]**:
`\n--related is skipped for a rather small sample size.`,
`\n--kinship analysis carried out instead for such a small sample size.`,
`  Relationship inference will be based on kinship estimation only.` — their exact triggers
are §8 item 11.

### 6.7 `<p>.con` from `--duplicate`

**Header:**
```
FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tConcord\tHomConc\tHetConc\n
```
**Row:**
```
%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\n
```
`N` = `M_ij`. These three concordance columns are the **only** `%.5lf` fields KING emits.
Both within- and cross-family pairs share this one file, in `.fam` serial order. Selection and
console text per §4.10.

### 6.8 `<p>.ibs` from `--ibs` — within-family pairs

**Header** (note the literal pipe characters):
```
FID\tID1\tID2\tZ0\tPhi\tN_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\tIBS\tDist\tHetConc\tHet2|1\tHet1|2\tHomConc\tKinship\tMaxIBD2\tPr_IBD2\n
```
**Row:**
```
%s\t%s\t%s\t%.3lf\t%.4lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf
  + \t%.3lf\t%.4lf\n          (MaxIBD2, Pr_IBD2 — only when informative segments exist)
```
Column mapping: `N_SNP` = `M_ij`, `N_Het1` = `N_i`, `N_Het2` = `N_j`, `NHetHet` = `N_HetHet`,
`NHomHom` = `N_HomHom`, `N_IBS0` = `N_IBS0`. `Kinship` is **Eq (9)**. Derived columns per §4.6.
Ordering: sorted (§6.1).

### 6.9 `<p>.ibs0` from `--ibs` — cross-family pairs

**Header:** same as `.ibs` but with `FID1\tID1\tFID2\tID2` in place of `FID\tID1\tID2`, and
**without `Z0` and `Phi`**:
```
FID1\tID1\tFID2\tID2\tN_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\tIBS\tDist\tHetConc\tHet2|1\tHet1|2\tHomConc\tKinship\tMaxIBD2\tPr_IBD2\n
```
**Row:** the same conversions minus `Z0`/`Phi`, then the trailing pair. `Kinship` is
**Eq (11)**. Ordering: serial.

Two quirks **[V]**:
* `.ibs0` writes the **literal string `-9`** for both `MaxIBD2` and `Pr_IBD2` on pairs without
  an IBD-segment analysis — bare `-9`, not `-9.000`/`-9.0000`. `.ibs` writes `0.000`/`0.0000`
  in the same situation.
* When there are **no informative segments at all** (stdout `No informative IBD segments.`),
  `MaxIBD2` and `Pr_IBD2` are **removed from the header and every row**, giving 20 / 19
  columns instead of 22 / 21.
* `MaxIBD2` is a base-pair magnitude — values like `48941787.000` occur. It is emphatically
  not `%g`. Its definition is §8 item 14.

### 6.10 `<p>allsegs.txt`

Written by `--ibdseg`, `--related`, `--ibs`, `--unrelated`, `--bysample`, `--bySNP` — i.e.
whenever the segment pre-pass runs — but **only if at least one chromosomal segment is long
enough**; otherwise stdout says `No informative IBD segments.` and the file is absent **[V]**.

**Header:**
```
Segment\tChr\tStartMB\tStopMB\tLength\tN_SNP\tStartSNP\tStopSNP\n
```
**Row:**
```
%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%d\t%s\t%s\n
```
`Segment` is a **1-based running index across the whole genome**, not per chromosome.
`Length == StopMB − StartMB` exactly. Positions are bp/1e6. X segments, when present, are
appended to the same file with `Chr` = 23. Content rule: §4.12.

Console: `Total length of %d chromosomal segments usable for IBD segment analysis is %.1lf Mb.\n`
then `  Information of these chromosomal segments can be found in file %s\n\n`.

### 6.11 `<p>.seg` from `--ibdseg` (Tier 2)

**Header:**
```
FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n
```
**Row:**
```
%s\t%s\t%s\t%s\t%.4lf\t%.4lf\t%.4lf\t<InfType>\n
```
One flat file containing **both** within- and between-family pairs (unlike `--related`, which
splits them), in `.fam` serial order. Numbers are always 4 decimals including `0.0000` and
`1.0000`; no `NA` sentinel.

### 6.12 `<p>unrelated.txt` / `<p>unrelated_toberemoved.txt`

**No header line.** Two TAB-separated columns:
```
%s\t%s\n           // FID, IID
```
Console: `A list of %d unrelated individuals saved in file %s` /
`An alternative list of %d to-be-removed individuals saved in file %s`.

### 6.13 Tier-2 output shapes (for planning only)

* `<p>bySample.txt` — SPACE-delimited, **6 dynamic header variants** keyed on X-SNP presence
  and PO/trio counts. Base: `FID IID FA MO SEX N_SNP Missing Heterozygosity`; optional blocks
  `N_xSNP xHeterozygosity` (iff X SNPs exist), `N_pair N_MIp Err_MIp` (iff ≥1 PO pair),
  `N_trio N_MIt Err_MIt` (iff ≥1 complete trio), `MI_Removal` (iff ≥1 PO pair). `MI_Removal`
  is `%G`. Serial order, one row per sample **including samples with no genotypes**.
* `<p>bySNP.txt` — SPACE-delimited, **3 header variants**. Base:
  `SNP Chr Pos Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate`, plus a PO block and a trio
  block. `Chr` prints the **symbol** `X` for chromosome 23. Covers **all** chromosomes,
  unlike the relatedness modes. Rows in `.bim` order.
* `<p>updateids.txt` / `<p>updateparents.txt` must match PLINK 1.9's `--update-ids`
  (`oldFID oldIID newFID newIID`) and `--update-parents` (`FID IID newPatID newMatID`)
  exactly — that is the documented downstream use **[D]**.
* `<p>splitped.txt` — SPACE-separated, 9 columns, no header.

### 6.14 Number-format summary (exhaustive for Tiers 1–2)

| Conversion | Columns |
|---|---|
| `%d` | `N_SNP`, `N`, `N_IBS0/1/2`, `NHetHet`, `NHomHom`, `N_Het1`, `N_Het2`, `Segment`, `Chr`, `Pos`, `SEX`, and every QC count |
| `%.3lf` | `Z0`; `StartMB`, `StopMB`, `Length` in `allsegs.txt`; `MaxIBD2` in `.ibs`; `Het` in X files |
| `%.4lf` | `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `Kinship`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `IBS`, `Dist`, `Het2\|1`, `Het1\|2`, `HomConc`, `Pr_IBD2`, `Missing`, `Heterozygosity`, `Freq_A`, `CallRate`, all `Err_*` |
| `%.5lf` | `Concord`, `HomConc`, `HetConc` — **`.con` only** |
| `%.1lf` | `MaxROH` (Tier 3); Mb totals in status messages |
| `%G` | `Error` (`.kin`), `MI_Removal` (`bySample.txt`) — values `0`, `0.5`, `1` |
| `%s` | all IDs, `InfType`, SNP names |
| literal `-9` | `MaxIBD2`, `Pr_IBD2` in `.ibs0` for unanalysed pairs |

---

## 7. RELATIONSHIP INFERENCE

### 7.1 The `InfType` vocabulary

Exactly seven literals, stored contiguously in the binary in this order **[S]**:
```
PO   FS   2nd   3rd   4th   UN   Dup/MZ
```

> **The manual is wrong.** It documents `Dup/MZTwin`; the binary writes **`Dup/MZ`** **[V]**.
> Match the binary.

`InfType` appears in `.seg`, and in `.kin`/`.kin0` **only on the `--related` segment path**.
`--kinship` has no `InfType` column at all.

### 7.2 The decision tree (segment-based) — first match wins

`π1 = IBD1Seg`, `π2 = IBD2Seg`, `π = PropIBD`:

```
if   π2 > 0.7                                            -> "Dup/MZ"
elif (π1+π2) > 0.96  or  ((π1+π2) > 0.9 and π2 < 0.08)   -> "PO"
elif π > 0.35355339059327373  and  π2 >= 0.08            -> "FS"     // 2^-1.5
elif π > 0.17677669529663687                             -> "2nd"    // 2^-2.5
elif π > 0.08838834764831845                             -> "3rd"    // 2^-3.5
elif π > 0.04419417382415922                             -> "4th"    // 2^-4.5
else                                                     -> "UN"
```

Provenance: KING's own emitted R script states the cut-points verbatim (`IBD2Seg>0.7`;
`IBD1Seg+IBD2Seg>0.96`; `>0.9 & IBD2Seg<0.08`; `PropIBD>0.35355 & IBD2Seg>=0.08`; `>0.17678`;
`>0.08839`; `>0.04419`) **[S]**, and the C++ writer's labels were verified to match it exactly
on **225 synthetic pairs with prescribed (π1, π2)**, with every transition bracketed to
±0.005 **[V]**.

Details that bite:
* **`2nd` has no upper bound.** A pair with π = 0.4481, π1+π2 = 0.8962, π2 = 0 labels `2nd`,
  not FS/PO. An entire arm of the π2 = 0 axis from π ≈ 0.18 to ≈ 0.45 is `2nd`. This is real,
  reproducible behaviour.
* **`FS` requires `π2 >= 0.08`** (`>=`, not `>`); `PO`'s second clause requires `π2 < 0.08`.
  Together with the π1+π2 clauses the PO/FS split is a two-piece boundary: a pair with
  π2 ≥ 0.08 is still `PO` if π1+π2 > 0.96.
* Use `>` for every PropIBD cut. Use the **exact powers of two**, not the rounded literals
  KING's R script prints.

### 7.3 The kinship degree bands (no segments)

Used for the relationship-summary table, the `Error` flag (§4.8) and `--degree` filtering.
Boundaries are geometric midpoints `2^-(d+3/2)` between adjacent degrees:

| Class | Interval on `Kinship` | Exact boundary | Published as |
|---|---|---|---|
| Dup/MZ | `> 0.3535533906` | `2^-1.5` | 0.354 |
| 1st degree (PO or FS) | `(0.1767766953, 0.3535533906]` | `2^-2.5` | 0.177 |
| 2nd degree | `(0.0883883476, 0.1767766953]` | `2^-3.5` | 0.0884 |
| 3rd degree | `(0.0441941738, 0.0883883476]` | `2^-4.5` | 0.0442 |
| 4th degree | `(0.0220970869, 0.0441941738]` | `2^-5.5` | 0.0221 |
| Unrelated | `<= 0.0220970869` | — | — |

**Use the exact powers of two, never the rounded manual values.** Boundary inclusivity is
§8 item 4.

### 7.4 PO vs FS without segments — scope correction

Two binary messages mention a data-derived IBS0 cutoff **[S]**:
```
1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf
Cutoff value for IBS0 between FS and PO is set at %.4f
```
Recon ranked this "the single highest parity risk in the whole reimplementation". Resolving the
strings against their owning functions changes that assessment **[V-NEW]**:

| String | Owning function | Reached by |
|---|---|---|
| `Cutoff value for IBS0 between FS and PO is set at %.4f` | `Engine::internalKING(int)` | the pedigree-reconstruction engine → `--build`, `--cluster` (Tier 2/3) |
| `1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf` | `Engine::MakeTrioForTDT()`, `Engine::MakeFamilyForMI()` | `--tdt` (Tier 3) and the Mendelian-error QC of `--bysample`/`--bySNP` (Tier 2) |

Neither is on the `--kinship`, `--related` or `--ibdseg` path, and neither message appeared in
any observed Tier-1 run. **For Tier 1 there is no PO/FS discrimination at all** (`--kinship`
emits no `InfType`), and for `--related`/`--ibdseg` the split comes from the segment rule of
§7.2. The cutoff formula therefore blocks only `--build` and the QC MI paths — §8 item 18.

Biologically it rests on the fact that a true PO pair shares ≥1 allele IBD at every autosomal
locus, so `N_IBS0 = 0` up to genotyping error, whereas FS has π0 = 0.25 and a clearly non-zero
IBS0 rate **[P]**.

### 7.5 The `π0` axis (KING-homo only — not implemented)

The paper's Table 1 adds a second inference axis using `π̂0 = N_IBS0 / Σ_m 2·p̂²(1−p̂)²`, with
`π̂2 = 4φ̂ + π̂0 − 1` and `π̂1 = 2 − 2π̂0 − 4φ̂`. These require allele frequencies and belong to
the **KING-homo** branch, which was retired in KING 2.3.0 (`--homog` is not in the 2.3.2 CLI).
Recorded for completeness; **do not implement**.

---

## 8. OPEN QUESTIONS

Each item states what is unknown, why it matters, and the exact experiment against the
reference binary that settles it. Ordered by parity risk within each tier.

**Reference binary for every experiment:**
`KING="/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"`

---

### Tier-1 blockers

**1. Does any Tier-1 analysis apply a variant-level MAF or call-rate filter?**
Resolved to "no" in §5.1 by pinning `"%d autosome SNPs with MAF>%.3lf and call rate>%d%% are
used."` to `Engine::MakeGRM0_LMM`, plus the observation that `N_SNP` equals the full autosomal
count on a fully-called fixture. Residual risk: a filter that happens to be a no-op on our
fixtures.
*Experiment:* build a fileset with 2 000 clean variants plus 200 deliberately awful ones
(MAF 0.001; 40 % missing; monomorphic; all-missing). Run `king -b f.bed --kinship --prefix a`
and compare `N_SNP` for a fully-called pair against the exact autosomal count. Repeat with the
200 removed. If `N_SNP` differs by exactly the number of awful variants, a filter exists and
its thresholds must be bisected; if it tracks the full count, close this item.

**2. Is `Kinship` in `.kin0` really the between-family estimator for *every* cross-FID pair,
including pairs KING might treat as "same family" through the FID-`0` pooling rule?**
`All individuals with family ID 0 are considered as relatives.` **[S]** suggests FID-0 samples
are pooled into one family, which would route them to Eq (9).
*Experiment:* build two filesets identical except that one uses distinct FIDs and the other
uses FID `0` for everyone. Run `--kinship` on both. If the FID-0 run puts all pairs in `.kin`
with Eq-(9) values, the pooling rule is confirmed and must be implemented before the estimator
is selected.

**3. Is the `--degree` filter on `.kin0` `>=` or `>`?**
The console says `kinship >= %.5lf`; the code may differ. Affects rows exactly on the boundary.
*Experiment:* construct a pair whose `Kinship` lands exactly on `2^-2.5` to more digits than
`%.4f` shows (tune by adding/removing variants with a binary search on `N_HetHet` and
`N_IBS0`, since the estimator is an exact rational). Run `--kinship --degree 1` and check
whether the row appears. Repeat at `2^-3.5` with `--degree 2`.

**4. Boundary inclusivity everywhere else:** the `Error`-flag bands of §4.8, the degree bands
of §7.3, and the `π2 >= 0.08` / `π > 2^-1.5` cuts of §7.2.
*Experiment:* the estimator is a ratio of integers, so exact boundary values are constructible.
For each boundary, build two pairs whose exact rational `Kinship` (or `PropIBD`) is one ULP
below and one ULP above, run `--kinship` / `--ibdseg`, and read the emitted class or `Error`.
Do this once for all 12 boundaries in one fixture with 12 pairs.

**5. What exactly is the within-family skip predicate?**
Observed: a within-family pair with `N_i + N_j == 0` produced no `.kin` row, while `N_i = 0`
with `N_j > 0` did produce one (§4.9). Candidate predicates that agree with both observations:
`N_i + N_j == 0`, `M_ij == 0 || N_i + N_j == 0`, or `N_HetHet + N_IBS0 == 0`.
*Experiment:* build three within-family pairs: (a) `N_i+N_j == 0` with `M_ij > 0`;
(b) `N_HetHet == 0` and `N_IBS0 == 0` but `N_i + N_j > 0`; (c) `N_HetHet == 0` with
`N_IBS0 > 0`. Run `--kinship` and record which rows appear.

**6. Confirm the between-family zero-denominator output byte-for-byte across sign cases.**
`Kinship 0.0000` was observed when `min(N_i,N_j) == 0` with a positive numerator.
*Experiment:* also construct `min == 0` with `N_IBS0 == 0` (numerator 0) and with a huge
`N_IBS0`. If all three print `0.0000`, it is an unconditional guard; if one prints something
else, the guard is on the ratio, not the denominator.

**7. The `M ≤ 512` sample screen: what exactly is `M`?**
Per-sample called-autosomal-variant count reproduces the observed exclusions (a 10-called
sample was excluded), but an independent sweep found the boundary tracking the *total*
autosomal variant count (`m = 544` excluded, `m = 545` not), independent of `n`, MAF, and
monomorphic fraction. The two readings are not reconciled.
*Experiment:* hold `m = 4096` fixed and vary a single sample's call count across
`511, 512, 513, 514`; then hold every sample fully called and vary `m` across
`540…548`. The first sweep pins the per-sample predicate, the second pins whatever additional
global gate produces the 545 boundary. Also check whether the *within-family* stage applies the
screen (§5.2 says it does not — confirm with a sample that is excluded yet has a family
partner).

**8. Is default two-stage screening lossless?**
We never screen; if KING's screening drops a pair we would report, the files differ.
*Experiment:* `king -b big.bed --related --degree 4 --prefix s1` vs
`king -b big.bed --related --degree 4 --noscreen --prefix s2`; `diff s1.kin0 s2.kin0`. Repeat
on a fixture engineered to contain pairs just above the degree-4 threshold. Any difference
means screening is lossy and must be modelled.

**9. The sort comparator: byte-wise or locale-aware, and what is the FID tie-break?**
§6.1's rules reproduce two independent 40-value/36-value probes exactly, but non-ASCII and
FID-collision behaviour is untested.
*Experiment:* one fixture with 60 crafted IDs including UTF-8 bytes, embedded spaces (if the
parser allows), very long digit runs (`000000000001` vs `1`), and two families whose sorted
member lists interleave. Compare `.kin` row order against our comparator on all
C(60,2) orderings implied.

**10. Bug-match or fix: the 0-byte `.kin` on single-family input.**
KING prints `There is only one family.` and writes a **0-byte** `.kin` (not even a header).
*Experiment:* none needed — this is a product decision. Recommendation: emit the header line
and gate the 0-byte behaviour behind `king-bugcompat`, since a headerless empty file breaks
every downstream reader.

**11. Exact triggers for the four "small sample size" messages.**
`--related is replaced with --kinship` fires below 10 samples **[V]**; the triggers for
`--related is skipped for a rather small sample size.`,
`--kinship analysis carried out instead for such a small sample size.` and
`Relationship inference will be based on kinship estimation only.` are unknown.
*Experiment:* sweep N = 2…20 with `--related`, and separately sweep the usable-segment total
(short chromosomes) at fixed N, capturing stdout each time. The first sweep separates the
sample-count gates; the second isolates the segment-availability gate.

**12. Confirm `HomIBS0`'s allele-order dependence and its exact denominator.**
Fitted as `N_IBS0 / |{m : (g_i = homA1 or g_j = homA1), both called}|`, matching on 4/4 fitted
pairs and 10/10 verified rows.
*Experiment:* take a fixture, run `--related`, then swap the `.bim` A1/A2 columns for 100 % of
variants (which flips every `00 ↔ 11` in meaning without touching the `.bed`) and re-run. Every
column except `HomIBS0` must be byte-identical; `HomIBS0` must change to
`N_IBS0 / |{homA2 union}|`. If `HomIBS0` is *also* unchanged, the denominator is
allele-symmetric and the fit is coincidental — re-fit.

**13. `Z0`/`Phi` for pedigrees beyond the simple cases.**
Verified only for PO / FS / 2nd / unrelated. Inbred loops, half-sib chains, three-generation
pedigrees, samples whose parents are not in the `.fam`, and FID-`0` pools are untested.
*Experiment:* one fixture containing a first-cousin marriage (inbred offspring), a
double-first-cousin pair, a half-sib pair, a great-grandparent pair, and a pair whose
connecting ancestor is ungenotyped. Run `--kinship` and diff `Z0`/`Phi` against §4.7's
recursive computation.

**14. `MaxIBD2` and `Pr_IBD2` in `.ibs`/`.ibs0`.**
`MaxIBD2` is a base-pair magnitude (`48941787.000`, `%.3lf`); `Pr_IBD2` is `%.4lf`; both are
dropped entirely when no informative segments exist, and `.ibs0` prints bare `-9` for
unanalysed pairs while `.ibs` prints `0.000`/`0.0000`.
*Experiment:* on a fixture with known IBD structure, cross-reference `MaxIBD2` against the
longest IBD2 run in the corresponding `.seg`/segment analysis, and `Pr_IBD2` against candidate
probability expressions. Needs the Tier-2 segment engine first; low priority because these two
columns are absent whenever segments are unavailable.

**15. `--unrelated`'s greedy algorithm and its output order.**
Reproducible run to run but neither `.fam` order nor sorted order.
*Experiment:* read Manichaikul et al. 2012 (permitted input) for the algorithm, then validate
by running `--unrelated --degree d` for d = 1…4 on a fixture with a hand-constructible optimal
answer (a star pedigree, a chain, and a clique). Diff both files including row order.

---

### Tier-2 blockers

**16. The IBD1/IBD2 calling rule.** Genuinely unspecified — the manual says the manuscript is
unpublished and the cited "Chen et al. 2024" does not exist in any index. Error tolerance,
boundary placement, and the word-granularity of the scan are all unknown.
*Experiment:* a fitting campaign, not a single run. Generate synthetic pairs sharing a single
IBD block of prescribed length at prescribed positions, with prescribed genotyping-error rates
(0 %, 0.1 %, 1 %), and read back `IBD1Seg`/`IBD2Seg`. Sweep block length across
1–50 Mb, marker density across 1–50 markers/Mb, and error rate; fit the boundary-extension
behaviour (observed: an *L* Mb IBD1 block is recovered as ~*L*+1–2 Mb, boundaries extending to
the next IBS0) and the error tolerance. Never read KING source. Acceptance criterion:
`IBD1Seg`/`IBD2Seg` byte-identical at `%.4f` on 200 held-out synthetic pairs.

**17. The usable-segment drop rule for `allsegs.txt`.** Cutting at gaps > 1 000 000 bp is
verified to the bp. The drop rule for small pieces brackets to "span > 10 Mb **and** ≥ 5
complete 64-marker words", but an 11 Mb / 352-marker / 5-word piece was still dropped.
*Experiment:* a 2-D sweep on one synthetic chromosome — span ∈ {10.0, 10.2, 10.5, 11, 12, 15,
20, 40} Mb × marker count ∈ {256, 320, 352, 353, 360, 368, 370, 384, 400, 512} — reading
`allsegs.txt` each time, plus a repeat with the piece offset so that word alignment shifts. This
scales every IBD proportion, so it must be closed before any `.seg` parity claim.

**18. The data-derived FS/PO IBS0 cutoff formula** (needed only for `--build` and the QC MI
paths, §7.4). KING prints it at `%.4f`.
*Experiment:* run `king -b f.bed --build` on several filesets with deliberately different MAF
spectra (uniform 0.05–0.5; all MAF ≈ 0.5; a rare-variant-heavy set) and identical pedigrees,
and capture the printed cutoff. First hypothesis to test: it is a fixed fraction of the mean
IBS0 rate among inferred-unrelated pairs (FS expectation is 0.25× the unrelated rate, so a
natural cutoff sits well below that). Fit the constant, then validate on a held-out spectrum.

**19. `X.seg` emission gate — ANSWERED.** Two conditions: a non-zero `--degree`, and an X
map that yields at least one usable segment. Swept `--degree` −1…10 and the X marker count
over 319…1000 on built filesets; there is no 512-marker threshold (that one belongs to
`--kinship`'s X pass alone) and no family-count condition. Rules and evidence:
`crate::analysis::xseg`.

---

### Lower-risk / hygiene

**20. Ordering of the `Options in effect:` echo block** when many options are set. Observed
one-per-TAB-indented-line; the binary also contains multi-option-per-line literals.
*Experiment:* run `--related --degree 3 --cpus 4 --prefix p --noscreen --minConc 0.9` and
capture the block verbatim; repeat with a different argv order to check whether the echo
follows command-line order or a fixed order.

**21. Does `-0.0000` ever print?** A `Kinship` of exactly `-0.0` would print with a leading
minus in both C and Rust, but KING may compute `+0.0`.
*Experiment:* construct a pair with `N_HetHet == 2·N_IBS0` (exact zero numerator) in both a
within- and a between-family position and inspect the emitted bytes.

**22. The non-deterministic spurious fatal on tiny `--bySNP` runs.**
`Too many first alleles as the major allele (~10.8%)` fired on 9/40 runs of a byte-identical
8-sample fixture whose true A1-major fraction is 0.00 %, with the percentage varying run to run;
`--cpus 1` does not fix it (7/40); never seen at 164 samples. This is uninitialised memory in
KING.
The stable part is now resolved: the intended window is the first 4,096 retained autosomal
markers and the boundary is 410 A1-major markers. open-king does **not** reproduce the short-map
tail read; it skips the gate until a complete window exists. Operational rule for reference
research remains: **never build regression fixtures at ~8 samples for `--bySNP`**; use
≥ 100 samples, and retry on failure. Successful runs are byte-deterministic.

**23. `--cpus` default and thread-invariance.** Default is half the logical cores; autosomal
`.kin0` is byte-identical between `--cpus 1` and `--cpus 8`, but the X writer races.
*Experiment:* run every Tier-1 analysis at `--cpus 1`, `2`, `8` on the same fixture and md5 the
outputs. Our implementation must be thread-invariant by construction (integer accumulators), so
this is a guard test, not a design input.

**24. Header emission when a file has zero data rows.** Observed: `.con` → header only;
`--related`'s `.kin0` → file absent; `--kinship`'s `.kin0` → header present; single-family
`.kin` → 0 bytes.
*Experiment:* one fixture per case (nothing above `--minConc`; nothing above `--degree`; a
single family; a fully cross-family cohort) and record file existence, size, and first bytes.
Encode the four outcomes as an explicit table in the writer.

---

## Appendix A — verification assets

### A.1 Golden vectors for unit tests

Genotypes as dosages `0/1/2`, `.` = missing, one SNP per character. Counts over pairwise-called
sites. Assert on the **exact rationals** in unit tests and on the `%.4f` strings only in format
tests.

| # | g1 | g2 | M | hethet | het1hom2 | het2hom1 | ibs0 | homhom | N_i | N_j | D_ij | between Eq(11) | within Eq(9) | naive (no correction) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| **A** dup/MZ | `1111002020` | `1111002020` | 10 | 4 | 0 | 0 | 0 | 6 | 4 | 4 | 0 | `1/2` → **0.5000** | `1/2` → **0.5000** | 0.5000 |
| **B** unequal het, no IBS0 | `1111102020` | `1102002020` | 10 | 2 | 3 | 0 | 0 | 5 | 5 | 2 | 3 | `1/8` → **0.1250** | `2/7` → **0.2857** | **0.5000** ✗ |
| **C** IBS0 + missing, equal het | `021120.120` | `2011201.00` | 8 | 2 | 0 | 0 | 3 | 6 | 2 | 2 | 12 | `−1` → **−1.0000** | `−1` → **−1.0000** | −1.0000 |
| **D′** IBS0 + missing + unequal het | `02.2000110` | `11.001211.` | 8 | 2 | 0 | 3 | 2 | 3 | 2 | 5 | 11 | `−7/8` → **−0.8750** | `−2/7` → **−0.2857** | **−0.5000** ✗ |
| **E** mild unequal het + IBS0 | `1210112010` | `1200110210` | 10 | 4 | 1 | 0 | 2 | 5 | 5 | 4 | 9 | `−1/16` → **−0.0625** | `0` → **0.0000** | 0.0000 ✗ |
| **F** zero min denominator | `0022002020` | `1212112121` | 10 | 0 | 0 | 6 | 1 | 4 | 0 | 6 | 10 | min = 0 → **`0.0000`** (§4.9) | `−1/3` → **−0.3333** | undefined |
| **G** no overlap | `01....2...` | `....12..01` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | **row omitted** | **row omitted** | undefined |

How to use them:
* **B and D′ are the discriminating tests.** They are the only rows where correct-between,
  correct-within and naive-without-correction give three different plausible answers. B catches
  a dropped correction term spectacularly: the naive form reports an unrelated pair as a
  duplicate.
* **A is worthless alone** — every wrong variant returns 0.5.
* **C is a trap**: IBS0 *and* missingness but equal het counts, so between and within coincide.
  A suite of only A and C cannot distinguish min from sum.
* **E** is the gentlest realistic case where between ≠ within (−0.0625 vs 0.0000).
* **F, G** are the degenerate fixtures for §4.9.

Additional invariants to assert as property tests over random genotype vectors:
```
Σ (x_i − x_j)²                 == N_IBS1 + 4·N_IBS0
min(N_i, N_j)                  == N_HetHet + min(het1hom2, het2hom1)
M_ij                           == N_HetHet + het1hom2 + het2hom1 + N_HomHom
within count-form              == within squared-difference form   (exactly)
within                         == between   whenever N_i == N_j
all counts except N_A1any      invariant under x → 2 − x           (A1/A2 swap)
simulated MZ → 0.5;  PO → 0.25 with IBS0 == 0;  unrelated → ~0
```
All nine published/observed algebraic forms of the two estimators were cross-checked over
20 000 random pairs with **0 mismatches** — the algebra is settled; do not re-litigate it.

### A.2 Parity ladder

1. **Unit:** golden vectors §A.1, exact rationals.
2. **Estimator:** our `Kinship` vs **Hail `hl.king`** (MIT, safe to read) on the same autosomal
   fileset — between-family only. Full double precision.
3. **Counts:** our integer accumulators vs `plink2 --make-king-table counts` on the *same*
   post-filter variant list. Expect exact agreement on `NSNP`, `HETHET`, `IBS0`, `HET1_HOM2`,
   `HET2_HOM1`. Expect the header, column names, row order and number formatting to differ —
   plink2 emits lower-triangular order "rather than KING's upper-triangular order" and its own
   source concedes "the header line still doesn't perfectly match KING due to e.g.
   capitalization". **plink2 is a value oracle, never a format oracle.**
4. **Within-family:** vs `snpgdsIBDKING(..., family.id=)` — the only readable implementation of
   Eq (9). (SNPRelate's default `family.id=NULL` uses Eq (11) everywhere; passing `family.id`
   is what activates the within/between split.)
5. **Byte parity:** vs the reference binary on the fixtures below. This is the only step that
   validates column order, header text and number formatting.

**Staging rule:** match `N_SNP` first, then the count-level quantities, only then `Kinship`.
A `Kinship` mismatch with a matching `N_SNP` is a formula bug; a mismatch with a differing
`N_SNP` is a filtering bug, and the two need completely different fixes.

### A.3 Reproducing this session's new verifications

All fixtures are tiny and synthetic; none contain real genome data.

```bash
KING="/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
SP=<scratchpad>/probe/run        # 10 samples x 49,600 SNPs (48,400 autosomal), ex.{bed,bim,fam}

# (a) estimator + every derived column, .kin/.kin0/.ibs/.ibs0 — 254 fields, 0 mismatches
#     decode ex.bed with the §3.3 rules, compute the §4.4 counts, apply §4.3/§4.6,
#     format with "%.4f", diff against $SP/{r_.kin, rd_.kin0, i_.ibs, k_.kin, k_.kin0}

# (b) .con definitions — 315 fields, 0 mismatches
"$KING" -b ex.bed --duplicate --minConc 0.01 --prefix cc_
#     Concord = N_IBS2/M, HomConc = (N_HomHom-N_IBS0)/N_HomHom,
#     HetConc = N_HetHet/(N_i+N_j-N_HetHet), all at %.5f

# (c) degenerate denominators
#     fixture 1: 12 samples x 5000 SNPs, FAM4 = {HOM1, HOM2} both all-homozygous,
#                FAM5 = {X1, X2} with disjoint call sets
"$KING" -b dg.bed --kinship --prefix dg_
#     -> .kin has no FAM4 row and no FAM5 row; .kin0 rows involving HOM1/HOM2 print Kinship 0.0000
#     fixture 2: FAM4 = {HOM1 (all hom), NORM1 (normal)}
"$KING" -b dg2.bed --kinship --prefix d2_
#     -> .kin row present: FAM4 HOM1 NORM1 ... Kinship -0.6647   (N_i=0, N_j>0, sum>0)
#     -> stdout: "The following 1 samples are excluded from the kinship analysis (M<512):"

# (d) string-to-function attribution (no source read; otool annotates literal pools)
grep -B12 "MAF>" king.formatstrings.byfunc.txt        # -> Engine::MakeGRM0_LMM
grep -B12 "Cutoff value for IBS0" king.formatstrings.byfunc.txt   # -> Engine::internalKING
grep -B12 "treated as parent-offspring" king.formatstrings.byfunc.txt
                                                      # -> MakeTrioForTDT, MakeFamilyForMI
```

### A.4 Recommended fixture set to commit

| Fixture | Shape | Exercises |
|---|---|---|
| `tiny` | 8 samples × 2 000 SNPs, 2 chromosomes, 1 % missing, real meiosis | `--kinship`, `--duplicate`, `--ibs`; the `--related` downgrade |
| `big` | 164 samples × 2 000 SNPs, 40 nuclear families + a dup/GC/GGC/half-sib family | the 16/14-column `--related` forms; `--degree` sweeps; the M ≤ 512 screen |
| `bitorder` | 6 × 4, designed so low-bits-first predicts bytes `0xE4`/`0x1B` | `.bed` decoding, padding, monomorphic, all-missing |
| `degenerate` | 12 × 5 000 with an all-hom pair, an all-hom/normal pair, a disjoint-call pair, a 10-called sample | §4.9, §5.2, §5.3 |
| `boundary` | 12 pairs constructed to sit one ULP either side of each of the 12 thresholds | §8 items 3, 4 |
| `chrom` | 200 variants across chr 1, 22, 23, 24, 25, 26, 0 | §3.4, `--sexchr` sweep |
| `ordering` | crafted FID/IID corpora (40 + 36 values) | §6.1 comparator |

Commit each fixture **plus the reference binary's expected output text**. Expected-output text
is factual data produced by running a program we are licensed to use for research; it is not
KING's source. Keep every fixture tiny and synthetic — the project's data rule forbids
committing real genomes, and a `.gq` bundle is family-level PHI.

---

## Appendix B — implementation checklist

1. `.bed`/`.bim`/`.fam` reader with the exact-length check (§3.7) and the chromosome partition
   of §3.4. `mmap` the `.bed`.
2. Four bit planes per sample, tail-masked, `W = ceil(m_autosome/64)` (§4.4).
3. Per-pair counts by popcount; **integers only**, one division at the end.
4. `phi_within` = Eq (9) → `.kin`; `phi_between` = Eq (11) → `.kin0`. Never clamp.
5. Derived columns exactly as §4.6 — in particular `Dist ≠ 2 − IBS` and `HomIBS0`'s A1-based
   denominator.
6. Pedigree `Z0`/`Phi` (§4.7) and the `Error` flag (§4.8).
7. Two ordering functions (§6.1) and the KING string comparator.
8. Writers per §6, with `%.4lf` everywhere the 64-bit path uses it, `%.5lf` only in `.con`,
   `%G` only for `Error`.
9. Console surface per §2, including the 78-column wrap algorithm and the BEL byte.
10. `king-bugcompat` feature gating exactly two behaviours: the `--noscreen` display value and
    the `3nd-degree` typo.
11. Property tests of §A.1 before any parity test; parity ladder of §A.2 after.
12. Provenance header + `THIRD_PARTY.md` entry (§0.1); delete `runtimes/king/` at parity.
