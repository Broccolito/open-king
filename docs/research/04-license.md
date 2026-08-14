# KING 2.3.2 — licensing, source availability, third-party libraries, and clean-room process

Recon note for the clean-room MIT reimplementation of KING's relatedness inference.
Compiled 2026-08-13. Not legal advice; escalate to counsel before shipping if any doubt remains.

**Reference binary examined:** `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(Mach-O 64-bit executable arm64, 1,815,336 bytes, mtime 2026-06-17 18:30).

---

## 0. Verdict in one paragraph

KING's C++ source **is publicly downloadable** (`KINGcode.tar.gz` from the author's site) but it is
**not open source**. There is no LICENSE/COPYING file anywhere in the distribution, the only
published terms are a single sentence on the download page, and **63 of the 172 source files carry
an inherited MERLIN header that explicitly forbids redistribution and forbids distributing modified
versions**. No grant of derivative-work rights exists. Therefore: **do not read, download, or open
KING's source code as part of this reimplementation.** Implement from (1) the 2010 Bioinformatics
paper, (2) the kingrelatedness.com manual/tutorial pages, (3) observed behavior and embedded string
constants of the local binary, and (4) genuinely open reimplementations — of which **Hail (MIT) is
the only one whose code we may actually copy from**.

---

## 1. Is KING's source publicly distributed?

Yes — as a plain tarball, with no version control, no issue tracker, and no license file.

| Fact | Value |
|---|---|
| Current release | **KING 2.3.2**, released **September 8, 2023** |
| Source archive (current) | `https://www.kingrelatedness.com/KINGcode.tar.gz` |
| Source archive (versioned, from history page) | `KING2.3.2code.tar.gz`, `KING2.3.1code.tar.gz`, … |
| Precompiled binaries | `Linux-king.tar.gz` (GNU/Linux x86-64), `Windows-king.zip`. **No macOS binary is published** — mac users are told to `brew install gcc` and build from source. |
| Documented build line | `c++ -lm -lz -O2 -fopenmp -o king *.cpp` |
| Author / copyright holder | **Wei-Min Chen** (Dept. of Genome Sciences / Public Health Sciences, University of Virginia) |
| Binary copyright string | `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen` |
| Page last updated | "Last updated: September 8, 2023 by Wei-Min Chen" |
| Download history | 30 releases listed, 1.0.0 (Oct 5, 2010) → 2.3.2 (Sept 8, 2023) |

There is **no** GitHub repository operated by the author for the KING core. The `chenlab-uva`
GitHub org (Chen's lab) publishes only satellite projects — `AncestryInference_KING` (R),
`InteractivePlots` / `InteractiveRelatednessPlots` (R), `king_testing` (Python),
`VirginiaKingWebserverDevelopment` (Perl) — none of which is the KING C++ core, and none of which
shows a license on the org overview page.

---

## 2. The license terms, verbatim

### 2.1 The only published terms — kingrelatedness.com/Download.shtml

The entire license, quoted byte-for-byte from the raw HTML (line 125 of the fetched page):

> `Feel free to use KING for your research, but please do not redistribute AND make profits.`

That is it. There is no other terms-of-use, copyright, or license statement on the Download page,
the homepage (`index.shtml`), the manual (`manual.shtml`), or the download-history page
(`history.shtml`). The Download page's only other legal-adjacent content is a citation request:

> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM (2010) Robust relationship
> inference in genome-wide association studies. *Bioinformatics* 26(22):2867-2873

**Correction to an existing project document.** `THIRD-PARTY-LICENSES` line 197-201 in this repo
states that the download page "carries no licence, copyright or terms statement at all." That is no
longer accurate (or was missed): the sentence above **is** on the page today, and it is a terms
statement — a permissive-sounding but restrictive one. The finding's *conclusion* is unaffected
(arguably strengthened: the one published sentence says "do not redistribute"), but the wording of
that entry should be updated to quote it.

### 2.2 How to read that sentence

- **"Feel free to use KING for your research"** — a use grant, scoped to *research*. It does not
  mention commercial use, clinical use, embedding in a product, or use by a company.
- **"but please do not redistribute AND make profits"** — grammatically ambiguous. Two readings:
  (a) conjunctive: *don't do both* (redistribute-for-profit is banned; non-profit redistribution is
  tolerated); (b) each is independently discouraged. Ambiguity in a permission grant resolves
  **against** the party relying on it — we would be relying on it, so we must assume the strict
  reading. Either way it is a *request* ("please"), not a formal grant, and it conveys **no**
  right to prepare or distribute derivative works.
- No warranty disclaimer, no patent grant, no sublicense right, no attribution formula, no
  choice of law. This is not a license in the OSI sense; it is a courtesy note.

### 2.3 The per-file headers — the real blocker

The source files carry per-file copyright headers. Two families exist.

**(a) MERLIN / libsrc files (Goncalo Abecasis).** Fetched verbatim from the public mirror
`statgen/topmed_variant_calling`, file `king/PedigreeGlobals.cpp`:

```
//////////////////////////////////////////////////////////////////////
// libsrc/PedigreeGlobals.cpp
// (c) 2000-2007 Goncalo Abecasis
//
// This file is distributed as part of the MERLIN source code package
// and may not be redistributed in any form, without prior written
// permission from the author. Permission is granted for you to
// modify this file for your own personal use, but modified versions
// must retain this copyright notice and must not be distributed.
//
// Permission is granted for you to use this file to compile MERLIN.
//
// All computer programs have bugs. Use this file at your own risk.
//
// Tuesday December 18, 2007
//
```

This repo's own prior audit (`THIRD-PARTY-LICENSES` lines 196-221) found **63 of the 172 files** in
`KINGcode.tar.gz` carry this header.

**(b) KING's own files.** The same header pattern with KING substituted, reported verbatim by web
search over the same mirror:

> "This file is distributed as part of the KING source code package and may not be redistributed in
> any form, without prior written permission from the author. Permission is granted for you to
> modify this file for your own personal use, but modified versions must retain this copyright
> notice and must not be distributed. Permission is granted for you to use this file to compile
> KING."

(Reported by search index over `github.com/statgen/topmed_variant_calling/blob/master/king/*`;
not independently re-fetched, deliberately — see §5. Its structure is identical to the MERLIN
header confirmed above, which is consistent with it being the same boilerplate adapted.)

### 2.4 Consequences — what is and is not permitted

| Action | Permitted? | Basis |
|---|---|---|
| Run KING locally for research | Yes | "Feel free to use KING for your research" |
| Run KING locally to generate reference output for our tests | Yes (research/verification use, local, not distributed) | same |
| Read the source | Legally yes (it's published) — **but see §5, we choose not to** | no confidentiality obligation attaches |
| Modify source for personal use | Yes, explicitly | per-file header |
| **Distribute modified versions** | **NO — explicitly forbidden** | per-file header: "must not be distributed" |
| **Redistribute the source** | **NO — explicitly forbidden** | per-file header: "may not be redistributed in any form, without prior written permission" |
| **Ship a compiled KING binary inside our app** | **NO** | building and shipping a binary from non-redistributable source is redistribution in object form |
| Relicense any part of it as MIT | **NO** | we are not the copyright holder; no sublicense right granted |

This is why `runtimes/` in this repo contains no `king/` directory: KING was removed 2026-08-12
(`docs/RUNTIMES.md` §"Intentionally not bundled"; `THIRD-PARTY-LICENSES` line 196). That decision
stands and this recon confirms it.

**Note the double encumbrance:** even if Wei-Min Chen granted us written permission tomorrow, ~37%
of the archive is Abecasis's MERLIN code, and Chen cannot grant rights in it. Permission would have
to come from *both* authors. That is one more reason the reimplementation route is the right one.

### 2.5 What is *not* encumbered

- **The algorithm itself.** The KING-robust and KING-homo estimators are published in
  Manichaikul et al. 2010, *Bioinformatics* 26(22):2867-2873. Mathematical methods are not
  copyrightable; an independent implementation of a published estimator is lawful.
- **Output formats.** Column names, field widths, header lines, and number formatting are facts
  about an interface, not expressive code. Extracting them from the binary's string table
  (`strings king`) and from the published manual is legitimate and is the correct source for
  byte-level parity.
- **Command-line flag names.** Documented on the manual pages.
- **Patents.** No patent covering the KING estimator surfaced in this search. Absence of evidence
  is not a clearance opinion, but the paper predates any obvious filing and none is asserted on
  the site or in the binary.

---

## 3. GitHub mirrors of KING source

| Location | What it is | License shown |
|---|---|---|
| `github.com/statgen/topmed_variant_calling` → `king/` | **A full copy of the KING C++ source (170 files)** vendored into the University of Michigan TOPMed variant-calling pipeline. Contains `KingCore.cpp/.h`, `Kinship.cpp/.h`, `KinshipX.cpp/.h`, `ShortKingCore.cpp`, `relationship.cpp`, `shortrelationship.cpp`, `ibdsegment.cpp`, `integrated.cpp`, `bigdataKING.cpp`, `SubsetKING.cpp`, plus the MERLIN libsrc set. | No repo-level license covering `king/`; the per-file headers described in §2.3 govern. **Vendoring by a third party does not create a license.** |
| `github.com/CDNMBioinformatics/KING_Wrappers` | Snakemake wrappers that *invoke* the KING binary. No KING source. | Wrapper's own; irrelevant. |
| `github.com/chenlab-uva/*` | Author's lab satellites (R plots, ancestry projection, webserver, test harness). Not the core. | Mostly unstated. |

**Do not treat any mirror as a license grant.** A third party republishing code that says "may not
be redistributed in any form" is evidence of a possible license violation upstream, not evidence
that we may copy it.

---

## 4. Third-party libraries KING uses

### 4.1 From the published build line

`c++ -lm -lz -O2 -fopenmp -o king *.cpp` →

- **libm** — C standard math library (system).
- **zlib** (`-lz`) — zlib license (permissive, BSD-like). Used for compressed I/O. Binary strings
  confirm it is *conditionally required*: `--ibdGRM cannot run without ZLIB`,
  `--ibdall cannot run without ZLIB`.
- **OpenMP** (`-fopenmp`) — the runtime, not the spec. On this arm64 build that is LLVM
  `libomp` (Apache-2.0 with LLVM Exceptions); with GCC it would be `libgomp` (GPLv3 + Runtime
  Library Exception).

### 4.2 Verified against the actual binary

`otool -L` on the local arm64 build:

```
/usr/lib/libz.1.dylib                        (current version 1.2.12)
/opt/homebrew/opt/libomp/lib/libomp.dylib    (compatibility version 5.0.0)
/usr/lib/libSystem.B.dylib
/usr/lib/libc++.1.dylib
```

Undefined-symbol audit (`nm -u`, 97 symbols total) — the complete external surface is:

- **libc/libm only**: `fopen fclose fread fwrite fprintf fscanf fgets fgetc fputs fputc fseek
  rewind feof fflush remove tmpfile printf sprintf snprintf vprintf vfprintf vsnprintf sscanf
  puts putchar getchar malloc free memcpy memmove memset bzero memset_pattern16 strcpy stpcpy
  strcat strchr strcmp strlen strtod atof atoi setjmp longjmp exit atexit system time clock
  localtime asctime exp exp2 log log10 pow sin atan __maskrune __tolower __toupper
  __DefaultRuneLocale __stdinp __stdoutp __assert_rtn __stack_chk_* __chkstk_darwin`
- **C++ runtime**: `__Unwind_Resume`, `___gxx_personality_v0`, `operator new/delete` (mangled).
- **OpenMP**: `__kmpc_barrier __kmpc_fork_call __kmpc_global_thread_num __kmpc_push_num_threads
  __kmpc_for_static_init_4 __kmpc_for_static_init_4u __kmpc_for_static_init_8
  __kmpc_for_static_fini __kmpc_reduce_nowait __kmpc_end_reduce_nowait omp_get_max_threads
  omp_get_thread_num`

**No BLAS, no LAPACK, no Eigen, no Boost, no GSL, no htslib, no Intel MKL** are linked in this
build. LAPACK is an *optional compile-time* dependency: the binary contains the string
`Please re-compile KING with LAPACK library.` — some analyses (large eigen-decompositions for
`--pca`/`--mds`, likely the `Largest %d eigenvalues:` path) degrade or refuse without it. Linear
algebra otherwise comes from KING's own vendored numeric files (`MathSVD`, `MathCholesky`,
`MathLu`, `MathMatrix`, `MathNormal`, `MathGenMin`, `MathGold`, `MathDeriv`, `MathSobol`,
`MathMiser`, `MathVegas`).

`_system` is in the undefined set, consistent with KING shelling out to **R** for `--rplot` and to
a PNG path for `--pngplot`; the binary emits `## %s for KING --…, by Wei-Min Chen and Zhennan Zhu`
banners at the top of the R scripts it writes. R is a runtime tool invocation, not a linked library.

### 4.3 Vendored source-level third-party code (from mirror file listing — names only, no code read)

The `king/` directory listing (170 files) shows a large **Abecasis MERLIN `libsrc` subset**
vendored wholesale. Identified by filename:

- Containers/strings/CLI: `StringBasics`, `StringArray`, `StringHash`, `StringMap`, `IntArray`,
  `LongArray`, `LongHash`, `LongInt`, `LongLongCounter`, `BasicHash`, `Hash`, `Sort`, `MerlinSort`,
  `QuickIndex`, `Parameters`, `Error`, `Constant.h`, `MemoryAllocators`, `MemoryInfo`,
  `InputFile`, `Input`, `FortranFormat`, `WindowsHelper`, `MiniDeflate`.
- Math: `MathVector`, `MathFloatVector`, `MathMatrix`, `MathSVD`, `MathCholesky`, `MathLu`,
  `MathStats`, `MathConstant.h`, `MathNormal`, `MathDeriv`, `MathGenMin`, `MathGold`,
  `MathSobol`, `MathMiser`, `MathVegas`, `MathAssoc`, `Random`, `OLS`, `OptimizerConstraints`.
- Pedigree/genetics: `Pedigree`, `PedigreeGlobals`, `PedigreeFamily`, `PedigreePerson`,
  `PedigreeDescription`, `PedigreeLoader`, `PedigreeAlleles.h`, `PedigreeAlleleFreq`,
  `PedigreeTrim`, `PedigreeTwin`, `Genetics`, `Matings`, `MapFunction`, `GenotypeLists`,
  `GenotypeCompressor`, `PeelerNodes`, `Intervals`, `Davies` (Davies 1980 AS 155 quadratic-form
  tail probability, used by SKAT-type tests), `BrentC`.

KING's own files are the rest: `Main.cpp`, `KingCore`, `ShortKingCore`, `Kinship`, `KinshipX`,
`relationship`, `shortrelationship`, `integrated`, `ibdsegment`, `ibdmapping`, `rohmapping`,
`bigdataKING`, `bigdataOffline`, `SubsetKING`, `ReadPLINK`, `structure`, `ancestry`, `admixture`,
`qc`, `autoQC`, `buildped`, `phase`, `poly`, `rplot`, `assoc`, `LMM`, `LMMSCORE`, `GRAMMAR`,
`FASTASSOC`, `SKAT`, `famSKAT`, `TDT.h`, `permuteTDT`, `rareTDT`, `ROADTRIPS`, `RSCORE`, `VC*`,
`RiskPrediction`, `TraitTransformations`, `diseaseGEE`, `IBD`.

The binary's own string table corroborates the MERLIN CLI heritage without any source access —
these are `Parameters.cpp`/`Error.cpp` hallmarks:

```
FATAL ERROR - 
Problems encountered parsing command line:
The following parameters are in effect:
Command line parameter %s (#%d) ignored
Command line parameter -%c%s: the option '%c' has no meaning
Command line parameter --%s is ambiguous
Command line parameter --%s is undefined
```

**Implication for us:** our Rust reimplementation needs none of this. We need a PLINK `.bed/.bim/.fam`
reader, popcount-based genotype arithmetic, and the output formatter. No SVD, no LAPACK, no pedigree
peeling for `--kinship`/`--related` kinship estimation.

---

## 5. RECOMMENDATION — the safe clean-room process

**Do not read KING's source code. Do not download `KINGcode.tar.gz`. Do not open any file under
`statgen/topmed_variant_calling/king/`.** This matches your prior. Reasons, in order of weight:

1. **No derivative-work grant exists.** Every route from "we read it" to "we ship it" requires a
   license we do not have. Even accidental structural similarity becomes hard to defend.
2. **We intend to relicense as MIT.** Publishing MIT code that a plaintiff can trace to
   all-rights-reserved source is the specific failure mode a clean room exists to prevent.
3. **Two copyright holders.** Chen and Abecasis. Contamination from either is fatal.
4. **The cost of avoiding it is near zero.** The estimators are two closed-form expressions in a
   published paper; the output format is fully recoverable from the binary's string table and the
   manual; and Hail's MIT implementation is available as a legitimate reference.

### Rules of engagement for the implementing agents

**ALLOWED (use freely, and log what you used):**

- Manichaikul et al. 2010, *Bioinformatics* 26(22):2867-2873, and its supplement.
- Every page under `https://www.kingrelatedness.com/` (manual, tutorial, index, history).
- `strings`, `nm`, `otool` on the local binary; running the binary on synthetic genotypes and
  diffing its output files byte-for-byte. **This is the primary source for format parity.**
- Hail (MIT), plink2 (GPLv3), SNPRelate (GPL-3), GENESIS (GPL) — read for *understanding*; copy
  code only from Hail, and only with attribution (see §6 for what each license permits).
- Third-party descriptions of KING's output format (e.g. GENESIS `kingToMatrix` docs, biostars
  threads, plink2 docs comparing to KING).

**FORBIDDEN:**

- `KINGcode.tar.gz`, any file inside it, any mirror of it, any GitHub blob under a `king/` source
  tree, any code-search snippet of KING source, any Stack Overflow/biostars post that pastes KING
  source lines.
- Decompilation/disassembly of the binary into pseudo-C for the purpose of transcribing logic.
  (Reading the *string table* is fine and is not decompilation; reconstructing control flow to
  copy it is not.)
- Asking an LLM to "recall" or reproduce KING source.

**Contamination note already in play:** the licensing audit recorded in this repo's
`THIRD-PARTY-LICENSES` (2026-08-12) states file counts and per-file headers from inside
`KINGcode.tar.gz` — i.e. **someone on this project has already opened the archive**, for licensing
review. That is defensible (license review is a recognized, non-implementing purpose) but it means
the clean-room record should be explicit: *the license reviewer and the implementer are separate
roles, and the implementer works only from the allowed-inputs list above.* Nothing from that audit
beyond license headers and file counts should reach the implementation.

### Record-keeping to produce alongside the code

1. A short `PROVENANCE.md` (or a section in the implementation's module doc) listing exactly which
   sources were consulted, with dates — paper, manual URLs, `strings` output, Hail file paths.
2. An explicit statement in that doc: *"KING's source code was not consulted. Output-format parity
   was derived from the KING 2.3.2 binary's embedded string constants and observed output files."*
3. MIT headers on our files; a `NOTICE`-style line crediting the *method* to Manichaikul et al.
   2010 (academic citation, not a license obligation) and any code genuinely derived from Hail
   under Hail's MIT terms (copyright line: `Copyright (c) 2015-2023, Hail Authors`).
4. Do **not** name the crate/module in a way that implies endorsement (`king-rs` invites confusion;
   something like `gq-relatedness` with "KING-robust estimator (Manichaikul et al. 2010)" in the
   docs is cleaner). KING is not a registered trademark we found, but the courtesy costs nothing.

---

## 6. Legitimately readable reimplementations (cross-check oracles)

| Tool | License | Copyable into MIT? | What it implements | How its output differs from KING's |
|---|---|---|---|---|
| **Hail** `hl.king()` | **MIT** (`Copyright (c) 2015-2023, Hail Authors`) | **YES**, with attribution | **Between-family (KING-robust) estimator only.** Documented formula: φ̂ᵢⱼ = ½ + (2N^{Aa,Aa} − 4N^{AA,aa} − N^{Aa}ᵢ − N^{Aa}ⱼ) / (4·min(N^{Aa}ᵢ, N^{Aa}ⱼ)) | Returns a **MatrixTable with one entry field `phi`** — no file format, no IBS0, no NSNP, no HetHet, no relationship inference. Zero format parity; pure numeric oracle. Explicitly does *not* implement the within-family estimator. |
| **PLINK 2.0** `--make-king`, `--make-king-table` | **GPLv3** (`plink-ng/2.0/COPYING`; `pgenlib` also LGPLv3 via `COPYING.LESSER`) | **NO** — GPL is incompatible with relicensing under MIT | KING-robust between-family estimator, heavily optimized | Files: `.king` + `.king.id` (triangular/square text matrix), `.king.bin` (float32/64), `.kin0` (table). `.kin0` columns: `#FID1 ID1 SID1 FID2 ID2 SID2 NSNP HETHET IBS0 HET1_HOM2 HET2_HOM1 IBS KINSHIP` — header line starts with `#`, tab-delimited. **Differences from KING:** (1) **autosomes only** — KING also uses X/XY in some analyses, so `NSNP` and hence estimates differ slightly when X is present; identical when no X is present. (2) **Pedigree is ignored** — the between-family estimator is used for *all* pairs, whereas KING switches to the within-family estimator for pairs in the same FID. (3) HETHET/IBS0/HET1_HOM2/HET2_HOM1 may be emitted as *proportions* or counts depending on modifiers. (4) Multiallelic variants use REF allele counts. |
| **SNPRelate** `snpgdsIBDKING()` | **GPL-3** | **NO** | Both `"KING-robust"` and `"KING-homo"` | Returns R **matrices** (`$IBS0`, `$kinship`, plus `sample.id`, `snp.id`, `afreq`); `KING-homo` returns `k0`/`k1` instead. No KING-compatible text file. Useful as a second numeric oracle, especially for the *homo* estimator which Hail doesn't implement. |
| **akt** (`akt kin -M 1`) | **PolyForm Strict License 1.0.0** — **not open source**: no distribution, **no derivative works**, non-commercial | **NO — and treat as read-with-care** | KING-style kinship on VCF/BCF | Writes a whitespace table to stdout, not KING's `.kin0`. Given the license forbids "making changes or new works based on the software," it is the *least* safe of the four to read; I recommend skipping it entirely — Hail covers the same ground under MIT. |
| **GENESIS** `kingToMatrix()` | GPL (Bioconductor) | NO | Nothing — it *parses* KING's own `.kin`/`.kin0` output | Valuable for a different reason: it is a **third-party documented description of KING's actual output columns**, i.e. a legitimate secondary source for format parity to cross-check against `strings`. |
| **sgkit** | Apache-2.0 | Yes (permissive, compatible) | `genetic_relatedness` / `pc_relate` | KING-robust support not confirmed in this search; treat as a maybe. Apache-2.0 is MIT-compatible, so if it does have it, it's a second copyable reference. |
| PLINK **1.9** | GPLv3 | NO | **Does not implement KING** (`--genome` is the method-of-moments PI_HAT estimator) | Not a cross-check for KING kinship. |

### Recommended verification strategy

1. **Byte-level format parity**: run the local `king` binary on small synthetic PLINK filesets
   (permitted: research use, local, output not redistributed), capture `.kin`, `.kin0`, `.seg`,
   `.segments.gz`, `.log`, then diff our output against it byte-for-byte. This binary is the only
   authority for field widths, `%.4f` vs `%.3f`, header spelling, and row ordering.
2. **Independent numeric oracle**: `plink2 --make-king-table` on the same fileset. Because plink2
   is autosomes-only and always between-family, restrict the comparison fileset to autosomes and
   to samples with distinct FIDs to make the two agree exactly; any residual disagreement in the
   third decimal is a bug in our estimator, not a format issue.
3. **Formula sanity check**: Hail's documented between-family formula (§6 table) is the cleanest
   written statement of the estimator outside the paper; it is MIT and quotable.
4. Keep the synthetic fixtures tiny and committed; never commit KING's outputs from real bundles
   (family-level PHI, per this repo's CLAUDE.md).

---

## 7. Sources consulted (all fetched 2026-08-13)

- https://www.kingrelatedness.com/Download.shtml — license sentence, build line, binaries, citation
- https://www.kingrelatedness.com/index.shtml , /manual.shtml , /history.shtml — no license text
- https://github.com/statgen/topmed_variant_calling — `king/` directory listing (170 filenames, via
  GitHub contents API — metadata only) and the `PedigreeGlobals.cpp` MERLIN license header
- https://github.com/chenlab-uva — author's lab repos (satellites only)
- https://raw.githubusercontent.com/hail-is/hail/main/LICENSE — MIT
- https://hail.is/docs/0.2/methods/relatedness.html — `king()` semantics and formula
- https://raw.githubusercontent.com/chrchang/plink-ng/master/2.0/COPYING — GPLv3
- https://www.cog-genomics.org/plink/2.0/distance and /formats — `--make-king*`, `.kin0` columns
- https://raw.githubusercontent.com/zhengxwen/SNPRelate/master/DESCRIPTION — `License: GPL-3`
- https://raw.githubusercontent.com/Illumina/akt/master/LICENSE — PolyForm Strict 1.0.0
- Local binary: `otool -L`, `nm -u`, `strings -a`
- This repo: `THIRD-PARTY-LICENSES` lines 196-221; `docs/RUNTIMES.md` §"Intentionally not bundled"
