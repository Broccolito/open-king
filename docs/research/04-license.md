# KING 2.3.2 — Licensing, Source Distribution, Third-Party Libraries, and Clean-Room Process

**Recon date:** 2026-08-13
**Target:** KING 2.3.2 (released 2023-09-08), Wei-Min Chen, University of Virginia
**Reference binary inspected:** `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king` (Mach-O 64-bit arm64, 1,815,336 bytes)
**Scope of this document:** legal posture only. Output formats are covered by the sibling recon docs.

> **Compliance note for this document.** No KING source code was downloaded, opened, or
> read while producing this file. `KINGcode.tar.gz` was deliberately *not* fetched. The
> `king/` subtree in `statgen/topmed_variant_calling` was identified but **not** opened
> beyond its repository-level metadata (license badge + directory presence). Everything
> below about KING's internals comes from (a) the public website, (b) the published paper,
> (c) `strings`/`otool`/`nm` on the compiled binary — i.e. facts about observable behavior
> and link-time dependencies, not source text.

---

## 1. TL;DR verdict

| Question | Answer |
| --- | --- |
| Is KING's C++ source publicly distributed? | **Yes** — `https://www.kingrelatedness.com/KINGcode.tar.gz` (and versioned `.../executables/KING2.3.2code.tar.gz`), plus a partial mirror inside `statgen/topmed_variant_calling`. |
| Under a recognized OSI license? | **No clear one.** The website states a bespoke, non-OSI research-use restriction. Third parties (Bioconda) label it GPL-3.0-or-later. The two are mutually contradictory. |
| Are derivative works / redistribution allowed? | **Ambiguous at best, prohibited at worst.** Website text forbids redistribution; a GPL-3 reading would *allow* derivatives but force GPL-3 copyleft on them. |
| Can we relicense any of it MIT? | **No, under either reading.** |
| Should we read the source for our reimplementation? | **NO. Do not fetch it, do not open it, do not have any agent open it.** |
| Is a clean-room MIT reimplementation legitimate? | **Yes** — the *method* is published (Manichaikul et al. 2010) and mathematical methods are not copyrightable; the *output format* is factual interface data. Only the C++ expression is protected. |

---

## 2. Is the source publicly distributed, and under what license?

### 2.1 Distribution — yes, three channels

1. **Primary, from the author's website** (`https://www.kingrelatedness.com/Download.shtml`).
   Verbatim from the page (HTML tags stripped, line 140-145 of the served document):

   > Source Code
   > Here is the source code in C++: KINGcode.tar.gz.
   > If your platform is not included here, you are welcome to download the source code and compile it in command line:
   > `wget https://www.kingrelatedness.com/KINGcode.tar.gz`
   > `tar -xzvf KINGcode.tar.gz`
   > `c++ -lm -lz -O2 -fopenmp -o king *.cpp`

   The versioned URL used by packagers is `https://www.kingrelatedness.com/executables/KING2.3.2code.tar.gz`
   (this is the `source.url` in the Bioconda recipe). The Download-History page (`history.shtml`)
   carries per-version tarballs, e.g. `KING2.3.1code.tar.gz`.

   Note the build line: `*.cpp` in a single flat directory — i.e. the tarball is a flat blob of
   KING's own translation units **plus** bundled third-party library translation units, all
   statically compiled into one binary. There is no separate `libsrc/` build step and no
   dynamic linkage to a system copy of the bundled library.

2. **Precompiled binaries** — `Linux-king.tar.gz` (GNU/Linux x86-64), `Windows-king.zip`.
   There is **no** official macOS binary on the download page; the page instead says:

   > For Mac users, if the executable is not working (with dyld errors), you may run "brew install gcc" first in order to run (or build) KING on your Mac.

   Consistent with that, our reference arm64 binary links `/opt/homebrew/opt/libomp/lib/libomp.dylib`,
   i.e. it was **locally compiled from `KINGcode.tar.gz` on a Homebrew Mac**, not downloaded
   from the vendor. (Relevant to our redistribution question — see §6.4.)

3. **Partial third-party mirror on GitHub.** `https://github.com/statgen/topmed_variant_calling`
   (University of Michigan Center for Statistical Genetics) contains a top-level `king/`
   directory holding KING C++ sources; its README's install steps include
   `g++ -O3 -c *.cpp; g++ -O3 -o king *.o -lz; cd ..`. That repository declares
   **Apache-2.0** at the repo level, which is almost certainly the CSG pipeline wrapper's
   license and **not** an authoritative relicensing of KING by its copyright holder
   (Wei-Min Chen is not a party to that repo). Treat this mirror as *tainted, do-not-read*
   exactly like the vendor tarball. There is no official `kingrelatedness`/Wei-Min-Chen
   GitHub organization or repository; searches surface only wrappers
   (`CDNMBioinformatics/KING_Wrappers`, `chenlab-uva/AncestryInference_KING`) and
   reimplementations (`populationgenomics/cuKING`).

### 2.2 The license — the author's own words

The **only** license statement the author publishes is a single sentence on the Download page
(verbatim, line 125 of the served HTML):

> **"Feel free to use KING for your research, but please do not redistribute AND make profits."**

Program banner, from the binary's string table (offset context in `king_strings.txt` line 1277):

> `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`

That is the entire published grant. Analysis:

- **It is not an OSI-approved license** and not any recognized public license family.
- **It is not a license text at all**, strictly speaking — it is a politely phrased request
  ("please do not…"), with no definitions, no grant clause, no warranty disclaimer, no
  termination clause, no patent grant.
- **Grant scope:** "use KING for your research." That is a *use* permission for research.
  It does not grant the copyright-relevant rights we would need: reproduction of source,
  preparation of derivative works, or distribution.
- **Restriction scope:** the conjunction is genuinely ambiguous. "do not redistribute AND
  make profits" can be read (a) conjunctively — only *commercial* redistribution is
  forbidden — or (b) as forbidding both redistribution and profiting. Courts construe
  ambiguous licenses against the drafter, but a prudent engineering org assumes the
  restrictive reading, especially for a bundled desktop product.
- **No explicit derivative-works permission exists under any reading.** Silence is not a grant.
- Neither `kingrelatedness.com/` (home), `manual.shtml`, nor `Download.shtml` contains a
  LICENSE/COPYING link, a copyright footer beyond the banner, or third-party attribution.

### 2.3 The contradicting third-party claim: GPL-3.0-or-later

Bioconda's recipe for `king` (`recipes/king/meta.yaml`, source `KING2.3.2code.tar.gz`) declares:

- `license: GPL-3.0-or-later`
- `license_family: GPL3`

That is a strong signal — Bioconda maintainers set that field by reading headers inside the
tarball, and they would not invent GPL-3. Corroborating circumstantial evidence from the
binary (see §3.1) is that KING statically bundles **Goncalo Abecasis's `libsrc` / libStatGen**
code, which *is* GPL-3-or-later. If GPL-3 headers ride along on those bundled files and KING
compiles them into the same executable (`c++ ... *.cpp`), then the distributed KING binary is,
on its face, a work whose distribution is governed by GPL-3.

**This creates a direct legal conflict, not a resolution:**

- If KING's own files carry GPL-3 headers → derivatives are permitted **but must be GPL-3**.
  We could not ship an MIT reimplementation derived from it. The website's
  "do not redistribute" sentence would additionally be a **GPL-3 §10 violation**
  (imposing a further restriction), making the license status internally inconsistent
  and litigable.
- If only the bundled Abecasis files are GPL-3 and KING's own files are all-rights-reserved →
  KING as distributed is arguably an infringing combination, and *KING's own* code remains
  proprietary. Still no MIT path.

**Either way the answer is the same for us: no readable-and-relicensable path exists.**
The ambiguity is itself a reason to stay out — a "we read it but only for X" defense is
much weaker when the license status is contested.

### 2.4 What we may and may not rely on

| Source | Read it? | Why |
| --- | --- | --- |
| Manichaikul et al. 2010, *Bioinformatics* 26(22):2867-2873 (PMC3025716) | ✅ **Yes** | Published method. Mathematical methods and algorithms are uncopyrightable (17 U.S.C. §102(b)); the paper's *expression* may not be copied verbatim, but the estimator equations may be implemented freely. |
| Chen et al. 2021 and Chen et al. 2024 (cited in the binary's own help text for `--ibdseg` / relatedness) | ✅ **Yes** | Same reasoning. |
| kingrelatedness.com manual/tutorial pages | ✅ **Yes** | Public documentation of interface and output semantics. Do not copy prose verbatim into our docs. |
| `strings`/`otool`/`nm`/behavioral runs of the shipped binary | ✅ **Yes** | Format strings, column labels, and error text are *facts about the interface*. Short functional labels ("HetHet", "PropIBD", "InfType", `%.4f` widths) are below the threshold of copyrightability (words/short phrases, `37 C.F.R. §202.1(a)`); the interface is the uncopyrightable "method of operation" under the *Lotus v. Borland* line, and *Google v. Oracle* (2021) makes reimplementing a declaring interface for interoperability strongly favored fair use. |
| `KINGcode.tar.gz` / `KING2.3.*code.tar.gz` | ❌ **NO** | See §5. |
| `statgen/topmed_variant_calling/king/**` | ❌ **NO** | Same code, mirrored. The repo's Apache-2.0 badge does not launder it. |
| Any AI/LLM output that may have been trained on and could regurgitate KING source | ⚠️ **Guard** | See §5.3. |

---

## 3. Third-party libraries KING uses

### 3.1 Statically bundled source (compiled into the binary from the tarball)

**Goncalo Abecasis's `libsrc` / libStatGen (GPL-3.0-or-later).** The compiled binary's string
table contains unmistakable libStatGen fingerprints (line numbers refer to
`scratchpad/king_strings.txt`, produced by `strings -n 5 king`):

```
   45  MathMatrix.h
   47  MathVector.h
 1414  MathFloatVector.h
 1501  OptimizerConstraints.cpp
 1114  IntArray::InnerProduct - vectors have different dimensions
 1115  IntArrays              - Left[%d] * Right[%d]
 3497  StringArray: Null String Access
 3609  11StringArray                       (Itanium-ABI mangled type name)
  519  FATAL ERROR -
 1405  Cholesky.Decompose: Matrix %s is not square
 1427  LU.Decompose: Matrix %s is not square
 1428  LU.Decompose: Matrix %s is singular
 1429  LU.Decompose: Matrix %s has zero pivot
 1433  Matrix.Add - Attempted to add incompatible matrices
 1435  Matrix.AddMultiple - Attempted to add incompatible matrices
 1437  Matrix.Multiply - Attempted to multiply incompatible matrices
 1474  No convergence in 30 SVD.Decomp iterations
 1410  FloatVector::AddMultiple - vectors are incompatible
 1477  Vector::AddMultiple - vectors are incompatible
```

`MathMatrix`, `MathVector`, `StringArray`, `IntArray`, `Pedigree*`, `Cholesky`/`LU`/`SVD`
decompositions and the `FATAL ERROR -` idiom are the University of Michigan CSG library used
by MERLIN / PEDSTATS / libStatGen. libStatGen's `StringArray.h` carries the standard
GPL-3-or-later header ("either version 3 of the License, or (at your option) any later
version"). This is almost certainly *why* Bioconda tags the package GPL-3.0-or-later.

Corroborating: the binary also carries MERLIN-format pedigree I/O text
(`Start writing reconstructed pedigrees in MERLIN format...`, `--merlin`,
`Error reading map file header...MARKER_ID, MARKER_NAME and BASE_PAIR_POSITION`) —
KING's pedigree/datafile layer *is* the Abecasis pedigree library.

**Implication for us:** the linear-algebra and pedigree-container layer of KING is GPL-3
third-party code. We are writing pure Rust with our own containers, so we neither need nor
want any of it. It does mean the "just read the source, it's basically open" temptation is
even worse than it looks: much of the source is *someone else's* GPL-3 code.

### 3.2 Dynamically linked (verified on the reference macOS arm64 binary)

```
$ otool -L king
king:
	/usr/lib/libz.1.dylib          (compatibility version 1.0.0, current version 1.2.12)
	/opt/homebrew/opt/libomp/lib/libomp.dylib (compatibility 5.0.0, current 5.0.0)
	/usr/lib/libSystem.B.dylib     (compatibility version 1.0.0, current version 1356.0.0)
	/usr/lib/libc++.1.dylib        (compatibility version 1.0.0, current version 2100.43.0)
```

| Library | Role | License | Notes |
| --- | --- | --- | --- |
| **zlib** (`-lz`, `libz.1.dylib`) | gzip I/O for IBD segment / GRM outputs | zlib license (permissive) | Compile-time optional. Binary contains `--ibdall cannot run without ZLIB` and `--ibdGRM cannot run without ZLIB`, and emits `.rohseg.gz`. Rust equivalent: `flate2` (MIT/Apache-2.0). |
| **OpenMP** (`-fopenmp`, `libomp.dylib`) | parallel kinship loops | LLVM `libomp`: Apache-2.0 WITH LLVM-exception | Confirmed by undefined symbols `_omp_get_max_threads`, `_omp_get_thread_num`; strings `--cpus %d`, `%d CPU cores are used to compute the pairwise kinship coefficients...`, `OMP loop starts.` Rust equivalent: `rayon` (MIT/Apache-2.0). |
| **libc++ / libSystem** | C++ runtime + libc | Apache-2.0 WITH LLVM-exception / APSL | Not applicable to a Rust port. |
| **LAPACK** | *optional, NOT linked here* | BSD-3-Clause | String at line 975: `  Please re-compile KING with LAPACK library.` — LAPACK is a compile-time-optional accelerator (eigen-decomposition for `--pca`/`--mds`, variance components). Our reference binary was built **without** it. `nm -u` shows no `dsyev`/`dgemm`/BLAS symbols. Relevant only if we ever port `--pca`; for relatedness we do not need it. |

`nm -u` otherwise shows only C stdlib (`_malloc`, `_memcpy`, `_fopen`, `_fread`, `_exp`,
`_log`, `_pow`, `_longjmp`, …). **No Boost, no Eigen, no Intel MKL, no htslib, no GSL.**

### 3.3 Runtime (out-of-process) dependencies — R, not linked

KING shells out to R for plotting (`--rplot` / `--pngplot`; strings show `R CMD BATCH %s` and
`--rpath`). The R scripts are **emitted from string constants embedded in the binary**, headed
with comments like:

```
## %s for KING --related, by Wei-Min Chen and Zhennan Zhu
## %s for KING Ancestry plot, by Zhennan Zhu and Wei-Min Chen
## %s for KING --ibdseg, by Wei-Min Chen
```

R packages referenced by those emitted scripts: **ggplot2** (GPL-2), **kinship2** (GPL-2+),
**igraph** (GPL-2+), **e1071** (GPL-2+), **doParallel** (GPL-2). Graceful-degradation strings
exist, e.g. `--roh is done but R plot is not available for missing R library ggplot2.` and
`Please rerun R code %s (or rerun KING) after ggplot2 is installed.`

**Two consequences.** (1) The embedded R scripts are substantial creative expression by named
authors — treat them as **source code, do not transcribe**, even though `strings` surfaces
them. Our scope is the numeric relatedness output, not KING's plots; GeneQuire renders its own
MDX/SVG figures. (2) An MIT Rust port has zero R dependency, which is a genuine product win
(no `runtimes/r` requirement for the relatedness path).

---

## 4. Open-source reimplementations of KING kinship (legitimate cross-checks)

Ranked by how safely we may read them **and** how useful they are as numeric oracles.

| Tool | License | Readable for an MIT port? | Notes / output differences |
| --- | --- | --- | --- |
| **Hail** `hl.king()` (`hail/methods/relatedness/king.py`) | **MIT** (Copyright 2015-2023, Hail Authors) | ✅ **Yes — fully. Best legal position.** MIT-in → MIT-out with attribution. | Returns a **square kinship MatrixTable/BlockMatrix**, not a pair table. Implements **KING-robust between-family only**. Diploid calls only. No `.kin0` text emission, no `InfType` classification, no `PropIBD`/IBD-segment logic. Ideal for validating the *kinship scalar* only. |
| **cuKING** (`populationgenomics/cuKING`) | **MIT** | ✅ **Yes** | CUDA reimplementation; README states results are "identical to Hail's `hl.king` implementation." Emits **Parquet** (`i`, `j`, `kin`, plus `ibs0`/`ibs1`/`ibs2`), converted to Hail tables by a helper script. Useful as a second MIT-licensed cross-check and for the IBS0/1/2 counts. |
| **somalier** (`brentp/somalier`) | **MIT** | ✅ **Yes** | Sketch-based relatedness for QC; "similar kinship estimates to KING in much less time." Uses its own selected-sites sketches, so numbers are *close but not bit-identical* to KING on the same BED. Good sanity oracle, not a parity oracle. |
| **PLINK 2.0** `--make-king` / `--make-king-table` / `--king-cutoff` | **GPL-3.0** | ⚠️ **Read docs freely; DO NOT read plink2 C source** and copy into MIT. | Best *numeric* oracle (documented as identical to KING when no chrX). Documented divergences below. |
| **SNPRelate** `snpgdsIBDKING()` | **GPL-3** | ⚠️ Docs yes, source no (copyleft). | Offers both `"KING-robust"` (returns `IBS0` proportion + `kinship`) and `"KING-homo"` (returns `k0`, `k1`). R-object output, no KING file formats. |
| **akt** (`Illumina/akt`) `akt kin` | **PolyForm Strict 1.0.0** | ❌ **NO — do not read.** | PolyForm Strict is *noncommercial + no-derivatives-flavored*; strictly worse than GPL for us. Repo archived 2026-04-20. htslib-based, VCF/BCF input, kinship in column 6 of a bespoke table. Skip entirely. |
| **NgsRelate** | **GNU (GPL)** | ⚠️ Docs only. | Genotype-likelihood MLE relatedness for low-depth NGS — a *different estimator*, not KING-robust. Not a parity oracle. |
| **KIMGENS / exKING-robust** | academic | Paper readable | Extends KING-robust to haploid-diploid systems. Method reference only; out of scope. |

### 4.1 PLINK 2.0 — the documented divergences from KING (important for our parity work)

Because plink2 is our most convenient numeric oracle, its *known* deltas from KING matter:

- **Chromosome scope.** plink2: *"Only autosomes are included in this computation."*
  KING uses autosomes **and** X/XY SNPs in some paths → `NSNP` (and, marginally, kinship) can
  differ. Where no X data is present, plink2 and KING agree exactly.
- **Pedigree handling.** plink2: *"Pedigree information is currently ignored; the
  between-family estimator is used for all pairs."* KING has both within-family
  (`--kinship` on `.kin`) and between-family (`.kin0`) estimators. Our port must implement
  both if we want KING parity; plink2 can only validate the between-family arm.
- **Multiallelics.** plink2: *"For multiallelic variants, REF allele counts are used."*
- **Scaling agrees.** Both scale so a duplicate pair is **0.5**, not 1 (1st-degree ≈ 0.25,
  2nd-degree ≈ 0.125, 3rd ≈ 0.0625).
- **`.kin0` column set differs from KING's.** plink2's `--make-king-table` columns (per its
  format spec) are:
  `FID1, ID1, SID1, FID2, ID2, SID2, NSNP, HETHET, IBS0, HET1_HOM2, HET2_HOM1, IBS, KINSHIP`
  — with `SID*` source-ID columns KING does not have, counts *or* proportions depending on
  modifiers, and an added Hamming-distance `IBS` column. plink2's own docs say it uses
  *"KING's original .kin0 text table format with minor changes,"* with *"row order … more
  friendly to incremental addition of samples."* By contrast, KING 2.3.2's own `.kin0`
  column tokens visible in the binary's string table are
  `FID1 / ID1 / FID2 / ID2 / N_SNP / HetHet / IBS0 / HomIBS0 / Kinship` and, for `--ibdseg`
  paths, `PropIBD` / `InfType`; plus an `N_IBS0` token and an `IID1` variant in other
  reports. **Therefore: use plink2 for numeric cross-checking only — never for
  column/format parity.** Exact KING formats come from the format-recon doc and from
  running the reference binary.
- **`--king-table-filter`** exists in plink2 to restrict output to high kinship; KING's
  analogous gate is `--degree`.

### 4.2 Suggested validation ladder

1. **Unit-level:** hand-computed KING-robust on tiny synthetic genotype matrices, derived
   straight from the paper's estimator equations.
2. **Estimator parity:** our Rust output vs **Hail `hl.king`** (MIT, safe to read) on the
   same autosomal BED — should match to floating-point tolerance.
3. **Numeric parity at scale:** our output vs **plink2 `--make-king-table`** restricted to
   autosomes — should match KINSHIP/IBS0/NSNP to tolerance.
4. **Byte parity:** our output vs the **shipped KING binary** on the real
   `bundled_data/` fixtures — this is the only thing that validates exact column order,
   header text, field widths, and `%.4f`/`%.3f` formatting.

Step 4 uses the binary as a black box. That is the clean-room-legal way to get byte parity:
compare bytes, never read source.

---

## 5. RECOMMENDATION — the safe process

**Confirmed: your prior is correct. Do NOT read KING's source. Implement from paper +
website + observed binary behavior only.**

### 5.1 Why, stated as a decision rather than a hedge

- Under the **website reading**, we have no reproduction/derivative-works grant at all;
  reading is at best unlicensed access and any resemblance in our code is exposure.
- Under the **GPL-3 reading**, reading is permitted but *any* derivation forces GPL-3 on
  GeneQuire's relatedness module — fatal for an MIT crate shipped inside a commercial-capable
  desktop app.
- The license is **internally contradictory**, so we cannot even get a clean opinion on which
  reading applies without counsel. An engineering team does not need to resolve that: the
  clean-room route costs us little and removes the question entirely.
- Third-party GPL-3 code (libStatGen) is statically bundled in the same tarball, so "reading
  KING's source" also means reading *Michigan's* GPL-3 source.
- The algorithm is **published, short, and closed-form**. KING-robust is a handful of counting
  statistics and a ratio. There is nothing in the source we actually need. The cost of the
  clean-room constraint here is near zero, which makes accepting any legal risk irrational.

### 5.2 Permitted inputs (whitelist — nothing else)

1. Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) *Robust relationship
   inference in genome-wide association studies.* Bioinformatics 26(22):2867-2873. PMC3025716.
2. Chen et al. 2021 and Chen et al. 2024 (the papers KING's own help text cites for
   `--ibdseg`/relatedness).
3. `kingrelatedness.com` manual, tutorial, and download pages (paraphrase; do not copy prose).
4. `strings`, `otool -L`, `nm`, `otool -tv` on the shipped binary; and **running** the binary
   on synthetic and fixture data and diffing its output files.
5. MIT-licensed reimplementations: **Hail `hl.king`**, **cuKING**, **somalier**.
6. Public docs (not source) of GPL tools: plink2, SNPRelate.

### 5.3 Prohibited inputs (blacklist — enforce on every agent, human and LLM)

1. `KINGcode.tar.gz`, `KING2.*code.tar.gz`, any `king*.cpp` / `king*.h`.
2. `statgen/topmed_variant_calling` `king/` subtree, and any other mirror of the same files.
3. plink2's C/C++ source, SNPRelate's C/R source, akt's source (any part).
4. Verbatim transcription of the R scripts embedded in KING's binary (they are authored code,
   even though `strings` reveals them). We do not need them.
5. **Prompting any LLM to "recall," "reproduce," or "show" KING's source.** Ask for
   implementations *from the paper's equations*. If a model volunteers something that looks
   like transcribed C++, discard it and restate the task from the paper.

### 5.4 Process controls to put in place

- **Written provenance header** in `crates/gq-core/src/relatedness/` (or wherever the module
  lands): a short comment stating this is a clean-room implementation from Manichaikul et al.
  2010 + observed binary I/O, that KING's source was never consulted, and listing the
  permitted inputs. Cheap, and it is the artifact that matters if the question is ever raised.
- **A `NOTICE`/`THIRD_PARTY.md` entry** recording: KING is not vendored, not linked, and not
  redistributed by GeneQuire; only its published output format is interoperated with.
- **Test fixtures, not source, are the oracle.** Commit small synthetic genotype fixtures plus
  the *expected output text* we generated by running the reference binary. Expected-output
  text is factual data produced by running a program we are licensed to use for research; it
  is not KING's source. Keep fixtures tiny and synthetic (the CLAUDE.md data rule already
  forbids committing real genomes).
- **Cite KING in our UI/docs** — "relatedness estimates computed with a KING-robust estimator
  (Manichaikul et al. 2010)" — which is scientifically correct, courteous to the authors, and
  makes clear we are not passing off their software.
- **Do not use the phrase "KING-compatible" as a product claim** in a way that implies
  endorsement; "implements the KING-robust estimator" is accurate and safe.

### 5.5 Independent, separate issue: stop shipping the KING binary

This is out of scope for the reimplementation but it is the *reason* the reimplementation
exists, so record it:

`runtimes/king/king` is currently committed to this repository. Under the vendor's own
sentence — *"please do not redistribute AND make profits"* — redistributing that binary in a
GitHub repo and in a shipped `.dmg` is at minimum against the author's stated wishes and, on
the restrictive reading, outside the granted permission. Under the GPL-3 reading it would
additionally trigger source-offer obligations (GPL-3 §6), and note the committed binary is a
**locally-built** arm64 artifact (Homebrew `libomp` path), i.e. our own build of GPL-3-flavored
sources, which is the clearest §6 trigger of all. **Once the Rust reimplementation reaches
parity, delete `runtimes/king/` from the shipped app and from git history-going-forward
(new commits; a rewrite is optional).** Keep the binary only outside the repo, as a local
parity oracle, at the read-only path
`/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`.

The same audit should be run over the rest of `runtimes/` (`plink` is GPL-3; `bcftools` is
MIT/GPL-dual; the JRE, python, node, r each carry their own terms) — but that is a separate
ticket, tracked in `docs/RUNTIMES.md`.

---

## 6. Raw evidence appendix

### 6.1 Verbatim license sentence, Download.shtml (tags stripped)

```
Feel free to use KING for your research, but please do not redistribute AND make profits.
```

Surrounding context (same page, in order): current release 2.3.2 released September 8, 2023;
pointer to Download History; "register at the KING User Forum"; "you may contact the developer
Dr. Wei-Min Chen directly"; then the license sentence; then Precompiled Binaries; then
Source Code; then the example dataset (`ex.tar.gz`, 1.35 MB, 332 HapMap samples — 165 CEU +
167 YRI — at 18,290 SNPs); then REFERENCE (Manichaikul et al. 2010).

### 6.2 Binary identity strings

```
KING 2.3.2 - (c) 2010-2023 Wei-Min Chen
Please check the reference paper Manichaikul et al. 2010 Bioinformatics,
Chen et al. 2024,
          or the KING website at kingrelatedness.com
Please check the reference paper Manichaikul et al. 2010 Bioinformatics,
Chen et al. 2021,
          or the KING website at kingrelatedness.com
```

### 6.3 Commands used for this recon (all read-only, on the local binary)

```bash
strings -n 5 king > king_strings.txt      # 3154 lines
otool -L king                              # dynamic deps
nm -u king                                 # undefined symbols
file king                                  # Mach-O 64-bit executable arm64
curl -s -L https://www.kingrelatedness.com/Download.shtml   # for verbatim license text
```

`KINGcode.tar.gz` was **not** downloaded.

### 6.4 Sources

- KING Download page (license sentence, source tarball, build line) — https://www.kingrelatedness.com/Download.shtml
- KING home — https://www.kingrelatedness.com/
- KING manual (citation, no license statement) — https://www.kingrelatedness.com/manual.shtml
- KING download history (versioned tarballs) — https://www.kingrelatedness.com/history.shtml
- Manichaikul et al. 2010, Bioinformatics 26(22):2867-2873 — https://pmc.ncbi.nlm.nih.gov/articles/PMC3025716/
- Bioconda recipe `king` (declares GPL-3.0-or-later; source `KING2.3.2code.tar.gz`) — https://bioconda.github.io/recipes/king/README.html and https://github.com/bioconda/bioconda-recipes/blob/master/recipes/king/meta.yaml
- statgen/topmed_variant_calling (contains a `king/` source subtree; repo badge Apache-2.0) — https://github.com/statgen/topmed_variant_calling
- libStatGen (GPL-3-or-later; `StringArray.h` header) — https://csg.sph.umich.edu/mktrost/doxygen/current/StringArray_8h_source.html and https://genome.sph.umich.edu/wiki/C++_Library:_libStatGen
- Hail `hl.king` source + LICENSE (MIT) — https://hail.is/docs/0.2/_modules/hail/methods/relatedness/king.html and https://github.com/hail-is/hail/blob/main/LICENSE
- Hail relatedness docs — https://hail.is/docs/0.2/methods/relatedness.html
- cuKING (MIT; "identical to Hail's hl.king") — https://github.com/populationgenomics/cuKING
- PLINK 2.0 distance/KING docs — https://www.cog-genomics.org/plink/2.0/distance
- PLINK 2.0 `.kin0` format spec — https://www.cog-genomics.org/plink/2.0/formats#kin0
- plink2-users, "Comparison between plink --make-king-table and king --kinship" (chrX/NSNP divergence) — https://groups.google.com/g/plink2-users/c/JiCBudDTwjY
- SNPRelate (GPL-3) `snpgdsIBDKING` — https://bioconductor.org/packages/release/bioc/html/SNPRelate.html and https://rdrr.io/bioc/SNPRelate/man/snpgdsIBDKING.html
- Illumina akt (PolyForm Strict 1.0.0; archived) — https://github.com/Illumina/akt
- somalier (MIT) — https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7362544/
- NgsRelate (GPL) — https://academic.oup.com/bioinformatics/article/31/24/4009/198242
