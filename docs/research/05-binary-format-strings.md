# KING 2.3.2 — Output Format Strings mined from the reference binary

**Recon target:** `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
Mach-O 64-bit executable **arm64**, 1 815 336 bytes, banner `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`.

**Legal posture.** Everything below is derived from (a) NUL-terminated literal constants in the
binary's `__TEXT,__cstring` section, (b) the exported C++ symbol table (which names the function
that *references* each literal), and (c) observed byte-for-byte output from running the stock
binary on synthetic PLINK data that I generated. No KING source code was consulted or
transcribed. Format strings and column headers are **facts about the file format**, not
implementation. Where I quote a disassembled instruction it is only to establish *which literal
goes into which file*, never to describe an algorithm.

---

## 0. Artefacts produced by this recon

| Path | Contents |
|---|---|
| `…/scratchpad/probe/king.strings.txt` | `strings -n 3` of the binary (3 645 lines) — the required raw dump |
| `…/scratchpad/probe/king.strings.all.txt` | `strings -a -n 3 -t d`, all sections (23 427 lines) |
| `…/scratchpad/probe/king.cstring.esc.txt` | **`__cstring` in link order**, 2 817 entries, `idx \t vmaddr \t escaped-string` — the primary evidence file |
| `…/scratchpad/probe/king.disasm.txt` | `otool -tV` full disassembly (346 343 lines); otool annotates every `adrp/add` with `literal pool for: "…"` |
| `…/scratchpad/probe/king.func_strings.demangled.txt` | **string list per C++ function, in code order** (1 508 functions) — the map from analysis → literals |
| `…/scratchpad/probe/king.formatstrings.byfunc.txt` | every literal containing a printf conversion, grouped by owning function (1 889 pairs) |
| `…/scratchpad/probe/king.xrefs.txt` | 7 885 `code-addr → cstring-addr → string` cross-references |
| `…/scratchpad/probe/run/` | live KING runs on a synthetic 10-sample / 49 600-SNP PLINK set + all output files |

Method note: `__cstring` is laid out in translation-unit order, so neighbours are meaningful;
but I did **not** rely on adjacency — every mapping in §2–§4 is pinned to a named function via
otool's literal-pool annotations, and the headline ones are additionally confirmed by running
the binary.

---

## 1. Build facts that change what this binary can emit

* `--ibdall cannot run without ZLIB` and `--ibdGRM cannot run without ZLIB` are the *entire*
  body of `Engine::AllIBDSegments()`. **This build has no zlib**, so the gzipped segment
  outputs (`.segments.gz`, `.rohseg.gz`) are never written; the only `.gz` literal in the
  binary belongs to input handling. `.rohseg.gz` appears only *inside the embedded R plotting
  script* — i.e. KING's own ROH plot script expects a file this build cannot produce.
* `  Please re-compile KING with LAPACK library.` exists → LAPACK is optional for `--ibdmds`.
* 32-bit ("Short") and 64-bit ("Long…64Bit") code paths coexist and **emit different column
  sets for the same file name**. On arm64 the 64-bit path always runs. §4 lists both so you
  don't accidentally reimplement the 32-bit variant.

---

## 2. MASTER TABLE — relatedness analyses (the ones we must reimplement)

Filenames are `<prefix>` + suffix, **no separator**; default prefix is `king` (`--prefix`).
`→` separates the header pieces KING concatenates; `⏎` is a literal `0x0A`.
Every file is opened `"wb"` (64-bit paths) or `"wt"` (32-bit paths) — identical on POSIX.

| Option | File | Owning function | Header line (exact) | Per-row printf |
|---|---|---|---|---|
| `--related` | `<p>.kin` | `Engine::IntegratedRelationshipInference()` | `FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship` → `\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType` *(iff IBD segs available)* → `\tError⏎` | `%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf` then `\t%.4lf\t%.4lf\t%.4lf\t%s\t%G\n` *(seg path)* **or** `\t%G\n` *(no-seg path)* |
| `--related` | `<p>X.kin` | same | `FID\tID1\tID2\tSex1\tSex2\tPhiX\tIBD1Seg\tIBD2Seg\tPropIBD⏎` | `%s\t%s\t%s\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n` |
| `--related` | `<p>.kin0` | same | `FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship` → `\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType` → `⏎` | `%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf` + `\t%.4lf\t%.4lf\t%.4lf\t` + *InfType* + `\n`; no-seg path uses `%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n` |
| `--related` | `<p>X.kin0` | same | `FID1\tID1\tFID2\tID2\tSex1\tSex2\tIBD1Seg\tIBD2Seg\tPropIBD⏎` | `%s\t%s\t%s\t%s\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\n` |
| `--kinship` (64-bit) | `<p>.kin` | `ComputeLongRobustKinship64Bit()` / `…WithFilter()` | `FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tKinship\tError⏎` | `%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n` |
| `--kinship` (64-bit) | `<p>.kin0` | same | `FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tKinship⏎` | `%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\n` |
| `--kinship` X (64-bit) | `<p>X.kin` | `ComputeLongRobustXKinship64Bit()` | `FID\tID1\tID2\tSex\tN_SNP\tPhiX\tHet\tIBS0\tKinshipX⏎` | `%s\t%s\t%s\t%s\t%d\t%.4lf\t%.3lf\t%.4lf\t%.4lf\n` |
| `--kinship` X (64-bit) | `<p>X.kin0` | same | `FID1\tID1\tFID2\tID2\tSex\tN_SNP\tHet\tIBS0\tKinshipX⏎` | `%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\n` |
| `--ibdseg` | `<p>.seg` | `ComputeIBDSegmentMain64BitFast()` | `FID1\tID1\tFID2\tID2\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType⏎` | `%s\t%s\t%s\t%s\t%.4lf\t%.4lf\t%.4lf\t` + *InfType* + `\n` |
| `--ibdseg` | `<p>X.seg` | `ComputeIBDSegment64Bit()` | `FID1\tID1\tFID2\tID2\tSex1\tSex2\tMaxIBD1\tMaxIBD2\tIBD1Seg\tIBD2Seg\tPropIBD⏎` ⚠ | `%s\t%s\t%s\t%s\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t` + `\n` — **9 fields under an 11-column header, and every row ends `\t\n`.** See §6.1 |
| `--ibdseg` | `<p>allsegs.txt` | `Engine::PreSegment(bool)` | `%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n` applied to `Segment,Chr,StartMB,StopMB,Length,N_SNP,StartSNP,StopSNP` | `%d\t%d\t%.3lf\t%.3lf\t%.3lf\t%d\t%s\t%s\n` |
| `--ibs` | `<p>.ibs` | `ComputeExtendedIBS64Bit()` | `FID\tID1\tID2\tZ0\tPhi\tN_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\t` → `IBS\tDist\tHetConc\tHet2\|1\tHet1\|2\tHomConc\tKinship` → `\tMaxIBD2\tPr_IBD2` → `⏎` | `%s\t%s\t%s\t%.3lf\t%.4lf\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf` + `\t%.3lf\t%.4lf` + `\n` |
| `--ibs` | `<p>.ibs0` | same | `FID1\tID1\tFID2\tID2\tN_SNP\tN_IBS0\tN_IBS1\tN_IBS2\tNHetHet\tNHomHom\tN_Het1\tN_Het2\t` → `IBS\tDist\tHetConc\tHet2\|1\tHet1\|2\tHomConc\tKinship` → `\tMaxIBD2\tPr_IBD2` → `⏎` | `%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf` + `\t-9\t-9` + `\n` |
| `--duplicate` | `<p>.con` | `ComputeBigDataDuplicate64Bit()`, `ComputeNoScreenDuplicate64Bit()` | `FID1\tID1\tFID2\tID2\tN\tN_IBS0\tN_IBS1\tN_IBS2\tConcord\tHomConc\tHetConc⏎` | `%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%.5lf\t%.5lf\t%.5lf\n` |
| `--roh` | `<p>.roh` | `Engine::ROH()` | `FID\tID\tFA\tMO\tSEX` → `\t<traitname>` → `\tMaxROH\tF_ROH` → `\tF_ROH_X`? → `\tNGENO_Y\tSexChr`? → `\tInfSex` → `\tSexErr`(only with Y) → `⏎` | `%s\t%s\t%s\t%s\t%d` + `\t%s`(affection) + `\t%.1lf\t%.4lf` + `\t%.4lf`(F_ROH_X) + `\t%d`(NGENO_Y) + one of `\tXY\tMALE` `\tXX\tFEMALE` `\tXO\tFEMALE` `\tXXY\tMALE` `\tNA\tMALE` `\tNA\tFEMALE` + `\t%s` |
| `--unrelated` | `<p>unrelated.txt` | `Engine::ExtractUnrelated()` | *(no header)* | `%s\t%s\n` |
| `--unrelated` | `<p>unrelated_toberemoved.txt` | same | *(no header)* | `%s\t%s\n` |
| `--cluster` | `<p>cluster.kin` | `Engine::ClusterFamily(int,int)` | `FID\tID1\tID2\tSex1\tSex2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType⏎` | `%s%d\t%s\t%s\t%d\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%s\n` (FID is `"KING"` + new-family number) |
| `--cluster` / `--build` | `<p>updateids.txt` | `ClusterFamily` | *(none)* | `%s\t%s` + `\t%s\t%s\n` |
| `--build` | `<p>updateparents.txt` | pedigree builder | *(none)* | `%s\t%s\t%s\t%s\t%d\t%d\t%d\n`, `%s\t%s\t%s\t%s\n`, `%s\t%s\t%s\t%s\t%d\t0\t1\n` |
| `--build` | `<p>build.log` | — | free text, `Details of pedigree reconstruction are available in log file %s\n` | — |
| any (auto) | `<p>splitped.txt` | `Engine::…SplitPed` | *(none)* | `%s %s %s_S%d %s %s %s %d %d %d\n` (split) / `%s %s %s %s %s %s %d %d %d\n` (unsplit) — **space-separated, 9 columns** |
| `--exact` | `<p>.kin` | `ComputeBigDataSecondDegree()` | `FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\tExact⏎` | `%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%s\t%s\n` |
| `--porel` | `<p>.por` | `ComputeBigDataPO()` | `FID_OFF\tID_OFF\tFID_PAR\tID_PAR\tFID_REL\tID_REL\tN_OR\tH_P\|O\tH_P\|R\tH_P\|OR⏎` | `%s\t%s\t%s\t%s\t%s\t%s\t%d\t%.3lf\t%.3lf\t%.3lf\n` |
| `--distant` | `<p>.dis` | `ComputeBigDataDistant()` | `FID1\tID1\tFID2\tID2\tKinship\tN_H\tN_HH\tFID3\tID3\tN_HHH\tH_3\|1\tH_3\|2\tH_3\|12⏎` | `%s\t%s\t%s\t%s\t%.3lf\t%d\t%d\t%s\t%s\t%d\t%.3lf\t%.3lf\t%.3lf\n` |

### 2.1 32-bit ("Short") variants of the same filenames — different columns

| Function (32-bit path) | File | Header | Row |
|---|---|---|---|
| `ComputeShortRobustKinship()` | `.kin` | `FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomConc\tKinship\tError⏎` | `%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n` |
| `ComputeShortRobustKinship()` | `.kin0` | `FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomConc\tKinship⏎` | — |
| `ComputeShortRobustXKinship()` | `X.kin` / `X.kin0` | same as 64-bit X headers | `%s\t%s\t%s\t%s\t%d\t%.4lf\t%.3lf\t%.4lf\t%.4lf\n` |
| `ComputeShortFastHomoKinship()` (`--homo`) | `.kin` | `FID\tID1\tID2\tN_SNP\tZ0\tPhi\tIBD0\tIBD1\tIBD2\tKinship\tError⏎` | `%s\t%s\t%s\t%d\t%.3lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%G\n` |
| `ComputeShortFastHomoKinship()` | `.kin0` | `FID1\tID1\tFID2\tID2\tN_SNP\tIBD0\tIBD1\tIBD2\tKinship⏎` | `%s\t%s\t%s\t%s\t%d\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n` |
| `ComputeShortExtendedIBS()` | `.ibs` / `.ibs0` | same as 64-bit IBS headers | same |

---

## 3. Header assembly is *piecewise* — the exact byte offsets

KING does not `fprintf` headers; it `memcpy`s literal chunks into a fixed buffer at known
offsets and then appends a short inline immediate. This is where the trailing columns come
from, and it is the only place the `Error` / `\n` bytes are defined.

**`.kin` (`--related`)** — buffer offsets confirmed by disassembly at `0x10007e0b4`–`0x10007e11c`:

```
offset  0 : "FID\tID1\tID2\tN_SNP\tZ0\tPhi\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship"   (60 bytes)
offset 60 : "\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType"                                       (32 bytes, conditional)
offset 60 or 92 : 8-byte immediate 0x000A_726F_7272_4509  ==  "\tError\n\0"
```

**`.kin0` (`--related`)** — `0x10007f48c`–`0x10007f510`:

```
offset  0 : "FID1\tID1\tFID2\tID2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\tKinship"      (59 bytes)
offset 59 : "\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType"                                       (32 bytes, conditional)
offset 59 or 91 : 2-byte immediate 0x000A  ==  "\n\0"      ← NO Error column on .kin0
```

**`.roh`** — `0x1000758f0`–`0x100075994`. Two of the header tokens exist *only* as inline
immediates and are absent from `strings` output:

```
0x0078_6553_666E_4909  ==  "\tInfSex\0"
0x0072_7245_7865_5309  ==  "\tSexErr\0"
```

so the ROH header is
`FID\tID\tFA\tMO\tSEX` `[\t<affection-trait-name>]` `\tMaxROH\tF_ROH` `[\tF_ROH_X]`
`[\tNGENO_Y\tSexChr]` `\tInfSex` `[\tSexErr]` `\n`.
`\tNGENO_Y\tSexChr` + `\tSexErr` appear only when the Y-chromosome SNP count `> 100`.
The affection-trait token is `DISEASEKING` when the trait is KING's synthetic one, else
`DISEASE`, else the user trait name (written via `\t%s`).

---

## 4. Verified byte-level output (stock binary, synthetic 10-sample / 49 600-SNP PLINK set)

Generated in `…/scratchpad/probe/run/`. `^I` = TAB, `$` = end of line.

### `--related` → `r_.kin`
```
FID^IID1^IID2^IN_SNP^IZ0^IPhi^IHetHet^IIBS0^IHetConc^IHomIBS0^IKinship^IIBD1Seg^IIBD2Seg^IPropIBD^IInfType^IError$
FAM1^IFA1^IKID1^I48400^I0.000^I0.2500^I0.2007^I0.0000^I0.3380^I0.0000^I0.2526^I1.0000^I0.0000^I0.5000^IPO^I0$
FAM1^IKID1^IKID2^I48400^I0.250^I0.2500^I0.2385^I0.0202^I0.4314^I0.1283^I0.2502^I0.5364^I0.2598^I0.5280^IFS^I0$
FAM1^IFA1^IMO1^I48400^I1.000^I0.0000^I0.1679^I0.0808^I0.2678^I0.4276^I0.0079^I0.0000^I0.0000^I0.0000^IUN^I0$
FAM4^IDUP1^IHS1^I48400^I1.000^I0.0000^I0.1794^I0.0455^I0.2934^I0.2516^I0.1119^I0.4490^I0.0000^I0.2245^I2nd^I1$
```
`Error` is `%G` and takes the values `0`, `0.5`, `1` (KING's own R plot script tests
`data.all$Error==0.5|data.all$Error==0`). `Kinship` can be negative and prints as e.g. `-0.0034`.

### `--related --degree 4` → `rd_.kin0`, `rd_X.kin0`
```
FID1^IID1^IFID2^IID2^IN_SNP^IHetHet^IIBS0^IHetConc^IHomIBS0^IKinship^IIBD1Seg^IIBD2Seg^IPropIBD^IInfType$
FAM1^IKID1^IFAM4^IDUP1^I48400^I0.3974^I0.0000^I1.0000^I0.0000^I0.5000^I0.0000^I1.0000^I1.0000^IDup/MZ$
```
```
FID1^IID1^IFID2^IID2^ISex1^ISex2^IIBD1Seg^IIBD2Seg^IPropIBD$
FAM1^IKID1^IFAM4^IDUP1^I1^I1^I0.0000^I1.0000^I1.0000$
```

### `--kinship` → `k_.kin`, `k_.kin0`, `k_X.kin`, `k_X.kin0`
```
FID^IID1^IID2^IN_SNP^IZ0^IPhi^IHetHet^IIBS0^IKinship^IError$
FAM1^IFA1^IKID1^I48400^I0.000^I0.2500^I0.2007^I0.0000^I0.2526^I0$

FID1^IID1^IFID2^IID2^IN_SNP^IHetHet^IIBS0^IKinship$
FAM1^IFA1^IFAM2^IMO2^I48400^I0.1665^I0.0834^I-0.0042$

FID^IID1^IID2^ISex^IN_SNP^IPhiX^IHet^IIBS0^IKinshipX$
FAM1^IFA1^IKID1^IMM^I1200^I0.0000^I0.390^I0.3725^I-0.2051$

FID1^IID1^IFID2^IID2^ISex^IN_SNP^IHet^IIBS0^IKinshipX$
FAM1^IFA1^IFAM2^IFA2^IMM^I1200^I0.390^I0.4175^I-0.3205$
```
Note `Sex` in the X files is the **pair-type string** `MM` / `FM` / `MF` / `FF`, not an integer
(`%s`), while `Sex1`/`Sex2` in `X.kin`, `X.kin0` and `X.seg` are integers (`%d`).
`Het` prints `%.3lf` (3 dp) whereas `IBS0` and `KinshipX` print `%.4lf`.

### `--ibdseg` → `s_.seg`, `sd_X.seg`, `s_allsegs.txt`
```
FID1^IID1^IFID2^IID2^IIBD1Seg^IIBD2Seg^IPropIBD^IInfType$
FAM1^IFA1^IFAM1^IKID1^I1.0000^I0.0000^I0.5000^IPO$
FAM1^IKID1^IFAM4^IDUP1^I0.0000^I1.0000^I1.0000^IDup/MZ$
FAM1^IFA1^IFAM2^IMO2^I0.0048^I0.0000^I0.0024^IUN$
```
```
FID1^IID1^IFID2^IID2^ISex1^ISex2^IMaxIBD1^IMaxIBD2^IIBD1Seg^IIBD2Seg^IPropIBD$
FAM1^IFA1^IFAM1^IKID2^I1^I2^I0.6155^I0.0000^I0.3077^I$     ← trailing TAB, 9 fields, 11-col header
```
```
Segment^IChr^IStartMB^IStopMB^ILength^IN_SNP^IStartSNP^IStopSNP$
1^I1^I0.134^I248.921^I248.786^I2200^Irs1_0^Irs1_2199$
23^I23^I0.160^I154.934^I154.774^I1200^Irs23_0^Irs23_1199$
```
`.seg` contains **both within-family and between-family pairs in one flat file** (`FID1`
and `FID2` are equal for within-family pairs). Chromosome 23 = X gets `Segment` number 23 in
`allsegs.txt` when X is present.

### `--ibs` → `i_.ibs`, `i_.ibs0`
```
FID^IID1^IID2^IZ0^IPhi^IN_SNP^IN_IBS0^IN_IBS1^IN_IBS2^INHetHet^INHomHom^IN_Het1^IN_Het2^IIBS^IDist^IHetConc^IHet2|1^IHet1|2^IHomConc^IKinship^IMaxIBD2^IPr_IBD2$
FAM1^IFA1^IKID1^I0.000^I0.2500^I48400^I0^I19025^I29375^I9714^I19661^I19217^I19236^I1.6069^I0.3931^I0.3380^I0.5055^I0.5050^I1.0000^I0.2526^I0.000^I0.0000$

FID1^IID1^IFID2^IID2^IN_SNP^IN_IBS0^IN_IBS1^IN_IBS2^INHetHet^INHomHom^IN_Het1^IN_Het2^IIBS^IDist^IHetConc^IHet2|1^IHet1|2^IHomConc^IKinship^IMaxIBD2^IPr_IBD2$
FAM1^IFA1^IFAM2^IFA2^I48400^I3927^I22139^I22334^I8019^I18242^I19217^I18960^I1.3803^I0.7820^I0.2659^I0.4173^I0.4229^I0.7847^I0.0010^I-9^I-9$
```
**`.ibs0` always writes the literal `\t-9\t-9` for `MaxIBD2`/`Pr_IBD2`** (string `\t-9\t-9`
at `0x10016cce9`), whereas `.ibs` writes `\t%.3lf\t%.4lf`.

### `--duplicate` → `d_.con`
```
FID1^IID1^IFID2^IID2^IN^IN_IBS0^IN_IBS1^IN_IBS2^IConcord^IHomConc^IHetConc$
FAM1^IKID1^IFAM4^IDUP1^I48400^I0^I0^I48400^I1.00000^I1.00000^I1.00000$
```
The three concordance columns are the only place KING uses **`%.5lf`**.

### `--roh` → `h_.roh` (dataset had no Y SNPs)
```
FID^IID^IFA^IMO^ISEX^IDISEASE^IMaxROH^IF_ROH^IF_ROH_X^IInfSex$
FAM1^IFA1^I0^I0^I1^I1^I0.0^I0.0000^I1.0000^IMALE$
```
`MaxROH` is the only `%.1lf` field in the relatedness outputs.

### `--unrelated` → `u_unrelated.txt`, `u_unrelated_toberemoved.txt`
```
FAM1^IFA1$      ← no header line at all
```

### `--ibdseg`/`--build` side effect → `s_splitped.txt` (SPACE-separated, 9 cols, no header)
```
FAM1 FA1 FAM1 FA1 0 0 1 1 0$
FAM4 DUP1 FAM4_S1 DUP1 0 0 1 1 0$
```

---

## 5. Row ordering — verified experimentally (critical for byte parity)

I re-ran with `FA1` renamed `ZFA1` to separate "`.fam` order" from "alphabetical order":

* **Within-family files (`.kin`, `X.kin`, `.ibs`)**: pairs are `i<j` over the family's members
  **sorted by IID string**, not `.fam` order. With `ZFA1` the output became
  `KID1-KID2, KID1-MO1, KID1-ZFA1, KID2-MO1, KID2-ZFA1, MO1-ZFA1`.
* **Between-family / flat files (`.kin0`, `X.kin0`, `.seg`, `.ibs0`, `.con`)**: pairs are `i<j`
  over the **global sample index, which is the `.fam` file order** (with `ZFA1` first because
  it is line 1 of `.fam`).
* Families with a single genotyped individual are skipped for `.kin`
  (`Each family consists of one individual.`).

---

## 6. KING quirks that a byte-exact clone must decide about

### 6.1 `X.seg` header/row mismatch + trailing tab (CONFIRMED by running the binary)
`Engine::ComputeIBDSegment64Bit()` writes the 11-column header
`FID1 ID1 FID2 ID2 Sex1 Sex2 MaxIBD1 MaxIBD2 IBD1Seg IBD2Seg PropIBD` (68 bytes, incl. `\n`)
but reuses the 9-conversion row format `"%s\t%s\t%s\t%s\t%d\t%d\t%.4lf\t%.4lf\t%.4lf\t"`
(literal at `0x10016487a`, shared with `ComputeIBDSegmentMain64BitFast`). The three doubles
are actually IBD1Seg, IBD2Seg, PropIBD; the header labels them MaxIBD1, MaxIBD2, IBD1Seg.
Because the format ends in `\t` and the code then appends `\n`, **every data row of `X.seg`
ends with `TAB LF`**. Reproduced: `sd_X.seg`.

### 6.2 `X.kin0` from `--kinship` is corrupted with >1 CPU (CONFIRMED)
`ComputeLongRobustXKinship64Bit()` writes `X.kin0` from multiple OpenMP threads without
serialising. With `--cpus 8` the file had **19 interleaved/garbled lines**; with `--cpus 1`
it had the correct **36 lines**. Example of the corruption:
```
FAM1^IFA1^IFAM2^IKID3^IMM^I^I120^I^I0.4020.2025^I-0.0041$
FAM2^IFA2^IFAM4^IDUP1^IMM^I1200^I0.396^I0.4200^I-0.3111FAM1^IFA1^IFAM3^IUNR1^IMF^I…
```
`.kin0` (autosome) was byte-identical between 1 and 8 CPUs. Our clone should emit the
`--cpus 1` (correct) form; do **not** try to reproduce the race.

### 6.3 Temp files with `$$$` in the name
`<prefix>$$$.kin0` and `<prefix>$$$.ibs0` are per-block scratch files (`"wb"`→`"ab"`→`"rb"`,
or `"wt"`/`"at"`/`"rt"` on the 32-bit path) that KING concatenates into the final file and
then deletes. They are visible mid-run.

### 6.4 Progress output uses carriage returns
Genotype loading prints `%d%%\r` repeatedly, which in a redirected log appears as
`0%6%13%19%25%31%38%44%50%56%63%69%75%81%88%…94%`.

### 6.5 `X.seg` is only emitted for `--ibdseg` under a size/degree condition
`--ibdseg` alone produced no `X.seg` on my data; `--ibdseg --degree 4` did. Do not assume it
is always written.

---

## 7. Shared vocabulary (exact literals)

**Relationship labels** (`__cstring` `0x1001604d9`–`0x1001604f5`, emitted in this storage order):
`PO`, `FS`, `2nd`, `3rd`, `4th`, `UN`, `Dup/MZ`.
The `%s` InfType field in `.kin`, `.kin0`, `.seg`, `cluster.kin` takes exactly these seven values.
Additional labels used by `--exact` / pedigree reconstruction (not InfType):
`HS`, `NA/NU`, `AN/UN`, `GcGp`, `GpGc`, `AV`, `GG`, `--`, `MZ`, `DZ`, `Pat`, `Mat`.

**Degree suffixes**: `st` and `nd` are separate literals, used as `%d%s-degree`
(`  Final Stage (with %d SNPs): %lli pairs of relatives (up to %d%s-degree) are confirmed\n`).

**Sex-pair strings (X analyses)**: `MM`, `FM`, `MF`, `FF`.
**ROH sex tokens**: `\tXY\tMALE`, `\tXX\tFEMALE`, `\tXO\tFEMALE`, `\tXXY\tMALE`,
`\tNA\tMALE`, `\tNA\tFEMALE`; plus bare `MALE`, `FEMALE`, `NA`.

**Sentinels**: `\t-9\t-9` (`.ibs0` MaxIBD2/Pr_IBD2), `-9` for missing affection,
`0` for missing parent id, `-99.999`.

**Relationship-summary block** (`Engine::printRelationship(int*,int*)`), printed to stdout:
```
Relationship summary (total relatives: %d by pedigree, %d by inference)
\n  Source\tMZ\tPO\tFS\t2nd\t3rd\tOTHER
  ===========================================================
  Pedigree\t%d\t%d\t%d\t%d\t%d\t%d\n
  Inference\t%d\t%d\t%d\t%d\t%d\t%d\n\n
```
There is a second, unused-looking header variant `\n        \tMZ\tPO\tFS\t2nd\t3rd\t4th`
with a 57-`=` rule (the `Source` variant uses 59 `=`).

---

## 8. Number-formatting rules (exhaustive for the relatedness outputs)

| Conversion | Where used |
|---|---|
| `%.3lf` | `Z0` (`.kin`, `.ibs`); `Het` in X kinship files; `PosMb`/`StartMB`/`StopMB`/`Length` in `allsegs.txt`; `H_*` in `.por`/`.dis` |
| `%.4lf` | `Phi`, `HetHet`, `IBS0`, `HetConc`, `HomIBS0`, `HomConc`, `Kinship`, `KinshipX`, `PhiX`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `MaxIBD1`, `MaxIBD2`, `F_ROH`, `F_ROH_X`, `IBS`, `Dist`, `Het2|1`, `Het1|2`, `Pr_IBD2` |
| `%.5lf` | `Concord`, `HomConc`, `HetConc` in `.con`; `Between-family relatives (kinship >= %.5lf)` message |
| `%.1lf` | `MaxROH`; Mb totals in status messages |
| `%G` | the `Error` column of `.kin` (values `0`, `0.5`, `1`); `%.2G`/`%.3G` for p-values elsewhere |
| `%d` | `N_SNP`, `N`, all IBS counts, `Sex1`/`Sex2`, `Chr`, `Segment` |
| `%lli` | pair counts in status messages (`%lli pairs of relatives …`) |
| `%s` | all IDs, `InfType`, `Sex` pair-type, SNP names |
| `%.4f` (float, no `l`) | only in the two threshold messages: `Cutoff value for IBS0 between FS and PO is set at %.4f\n` and `Cutoff value between full siblings and parent-offspring is set at %.4f\n` |

Separator: **TAB everywhere** in `.kin*/.seg/.ibs*/.con/.roh/allsegs.txt/unrelated.txt`.
**SPACE** is used in `bySample.txt`, `bySNP.txt`, `splitped.txt`, `pc.txt`, `_popref.txt`,
`_Dist.txt`, `.fam`, `af.ped`. Never mix.

---

## 9. Console message catalogue — relatedness paths, in emission order

### `--related` (`Engine::IntegratedRelationshipInference`)
```
No genotype data
Autosome genotypes stored in %d words for each of %d individuals.\n      (two literals concatenated)
\nOptions in effect:
\t--related        \t--noscreen       \t--degree %d\n     \t--sysbit 64
\t--cpus %d\n      \t--lessmem        \t--rplot           \t--prefix %s\n
  Relationship inference will be based on kinship estimation only.
Please consider KING options --kinship or --ibdseg.
Family %s skipped.\n
This KING version cannot handle a single family with >65536 individuals.
Each family consists of one individual.
Within-family kinship data saved in file %s\n
All individuals with family ID 0 are considered as relatives.\n
There is only one family.
Within-family X-chr IBD-sharing inference saved in file %s\n
A subset of informative SNPs will be used to screen close relatives.
Sorting autosomes...
Relationship inference across families starts at %s
%d CPU cores are used...\n
  Stages 1&2 (with %d SNPs): %lli pairs of relatives are detected (with kinship > %.4lf)\n
                               Screening ends at %s
                                         ends at %s
No close relatives are inferred.\n
  Final Stage (with %d SNPs): %lli pairs of relatives (up to %d%s-degree) are confirmed\n
                               Inference ends at %s
  %lli pairs of relatives (up to %d%s-degree) are identified\n
\nBetween-family relatives (kinship >= %.5lf) saved in file %s\n
  X-Chr IBD-sharing inference saved in file %s\n
No cryptic relatedness (up to the %d-degree) is found.\n
\nNote only duplicates and 1st-degree relatives are included in the inference.
  Specifying '--degree 2' if a higher degree relationship inference is needed.\n
This version of KING supports up to %d samples.\n
Please contact the KING authors to allow an even larger sample size.
```

### `--ibdseg` (`Engine::ComputeIBDSegment64Bit` + `PreSegment`)
```
64-bit system is required.
Chromosomes unsorted: %s on chr %d, %s on chr %d.
Positions unsorted: %s at %d, %s at %d.
No informative IBD segments.
\nToo many first alleles as the major allele (~%.1lf%%). Please use plink1.9 --make-bed to regenerate the genotype data again.\n
Total length of %d chromosomal segments usable for IBD segment analysis is %.1lf Mb.\n
  In addition to autosomes, %d segments of length %.1lf Mb on X-chr can be further used.\n
  Information of these chromosomal segments can be found in file %s\n\n
Segments too short.
IBD segment analysis starts at %s
  Note chromosomal positions can be sorted conveniently using other tools such as PLINK.
%d CPU cores are used for %s inference...\n            (%s = "autosome" | "X-chr")
                       ends at %s
\nNote with relationship inference as the primary goal, the following filters are applied:
  Sample pairs without any long IBD segments (>10Mb) are excluded.
  Short IBD segments (<3Mb) are not reported/utilized.
Summary statistics of IBD segments for individual pairs saved in file %s\n
Additional summary statistics of X-Chr IBD segments saved in file %s\n
Sample FID=%s,IID=%s is removed with only %d non-missing SNPs.\n
Only pairs between the first %d and the last %d individuals are inferred.\n
Only pairs between %d unaffected and %d affected individuals are inferred.\n
Either affection status needs to be assigned for --projection analysis,
  or count of the first part (N) can be specified as --projection N.
No valid samples for inference.
```

### `--kinship`
```
The following %d samples are excluded from the kinship analysis (M<%d):\n   then  \t(%s %s)
Within-family kinship data saved in file %s\n
Between-family kinship data saved in file %s\n
Between-family kinship data (up to degree %d, %lli pairs in total) saved in file %s\n
Note --kinship --degree <n> can filter & speed up the kinship computing.
\nX-chromosome analysis...
No sufficient X-chromosome SNPs (%d) available for the robust analysis.\n
X-chromosome genotypes stored in %d 64-bit words for each of %d individuals.\n
```

### `--duplicate`
```
Computing pairwise genotype concordance starts at %s
  %d CPU cores are used...\n
        Stage 1 (with %d SNPs) screening ends at %s
        Stage 2 (with all SNPs) inference ends at %s
                                            ends at %s
Duplicate inference ends at %s
%lli pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s\n\n
%d pairs of duplicates with heterozygote concordance rate > %d%% are saved in file %s\n\n
No duplicates are found with heterozygote concordance rate > %d%%.\n\n
No duplicates are found with heterozygote concordance rate > 80%%.\n\n
  %d additional pairs from screening stage cannot be not confirmed in the final stage\n\n   [sic]
  %d additional pairs from screening stage not confirmed in the final stage\n\n
```

### `--roh`
```
Run of homozygosity analysis starts at %s
Run of homozygosity analysis   ends at %s          (three spaces, aligns with "starts")
  Short ROHs (<5Mb) are discarded.
\nThe following %d persons have excessive run of homozygosity (F>0.0884):\n
  Person (%s, %s) has inbreeding coefficient = %.4lf\n
Run of homozygosity summary saved in file %s\n
```

### Thresholds mentioned in messages
```
1st-degree relatives are treated as parent-offspring if IBS0 < %.4lf\n
Cutoff value for IBS0 between FS and PO is set at %.4f\n
Cutoff value between full siblings and parent-offspring is set at %.4f\n
Minimum segment length is set as %d bp\n.                       [sic: period after \n]
KING supports minimum segment length from 1 to 10 Mb at the moment.
Default seglength of 3Mb is used.
minConc value is out of range and not specified.
```

### Classification cut-points, as they appear verbatim in KING's own embedded R scripts
(these are string constants, i.e. KING's published definition of its own output):
```
d0    <- data$IBD2Seg>0.7
d1.PO <- (!d0) & data$IBD1Seg+data$IBD2Seg>0.96 | (data$IBD1Seg+data$IBD2Seg>0.9 & data$IBD2Seg<0.08)
d1.FS <- (!d0) & (!d1.PO) & data$PropIBD>0.35355 & data$IBD2Seg>=0.08
d2    <- data$PropIBD>0.17678 & data$IBD1Seg+data$IBD2Seg<=0.9 & (!d1.FS)
d3    <- data$PropIBD>0.08839 & data$PropIBD<=0.17678
d4    <- data$PropIBD>0.04419 & data$PropIBD<=0.08839
dU    <- data$PropIBD<=0.04419
```
and for kinship-space plots: `abline(h = 0.35355 / 0.17678 / 0.0884 / 0.04419 / 0.02210)`,
`abline(a = 0.3535534, b = -0.5)`, `abline(a = 0.08838835, b = -0.5)`, `abline(a = 0.04419, b = -0.5)`,
`y.cut <- 2^-2.5`, and ROH `F>0.0884` / `roh$F_ROH > 2^-4.5`.
The same constants appear as immediates in the `.kin0` writer (0.8, 0.96, 0.9, 0.08, 0.35355…),
and `PropIBD = IBD2Seg + 0.5·IBD1Seg` is confirmed by an `fmov d8, #0.5` + `fmadd` in the
`X.seg` writer.

---

## 10. Error / warning message catalogue (I/O and QC, relatedness-relevant subset)

```
\nFATAL NUMERIC ERROR -            \nFATAL ERROR -            \nWARNING -
Cannot open %s to write.          Cannot open %s to write            Cannot open %s to read
File %s cannot open to write
Genotype file %s cannot be opened          Please use PLINK binary format as input.
Please use either PLINK or KING binary format as input.
The binary file cannot open\nPlease check the reference paper Manichaikul et al. 2010 Bioinformatics,\n\t\t\t\t\tChen et al. 2021,\n          or the KING website at kingrelatedness.com
Genotype files are required. e.g.,\n  king -b ex.bed --related\n\nPlease check the reference paper Manichaikul et al. 2010 Bioinformatics,\n\t\t\t\t\tChen et al. 2024,\n          or the KING website at kingrelatedness.com
Map file %s cannot be opened               Pedigree file %s cannot be opened
The PLINK format in %s cannot be recognized. The first byte is %x
The PLINK format in %s cannot be recognized. The second byte is %x
Currently only SNP-major mode can be analyzed.
Not enough memory: %d MB is needed         Not enough memory
No autosome SNPs                           No autosome SNPs are available. Please check your map file.
No autosome genotypes available for KING inferences.
Not enough genotypes at the %dth marker\n
Not enough genotypes at the %dth marker on sex chr %d\n
Not enough genotypes at %dth marker\n      Not enough genotypes at %dth SNP %s\n
No genotype data.    No genotype data     No genotype variation.    No trait variation.
Degree of relatedness not defined.
This function is currently disabled for tiny dataset with sample size < 10.
Cluster can only run on 64-bit machine.
Paternity inference can only run on a 64-bit flatform.        [sic "flatform"]
Only available for 64-bit system.          64-bit system is required.
Cannot run --exact analysis without IBD segments.
Cannot run --exact analysis without 64-bit system.
The current KING version cannot process %d samples.
This KING version cannot handle a single family with >65536 individuals.
Sex chromosome %d out of range.\n
Non-human samples are analyzed, with %d pairs of chromosomes\n
  %d samples have identical IDs between %s.fam and %s\n
  REF_ is added to all IDs in %s.fam.\n     QRY_ is added to all IDs in %s.fam
KING cannot handle %d samples with identical IDs between %s.fam and %s\n
          %d SNPs are removed for appearing more than once\n
No sufficient X-chromosome SNPs (%d) available for the robust analysis.\n
Please do not run --related together with --autoQC          … together with --homog
Please do not run --ibdseg together with --autoQC / --homog / --pca
Please do not run --roh together with --autoQC / --homog / --pca
--related is skipped.     --ibdseg is skipped.     --roh is skipped.
\n--related is replaced with --kinship for a small sample size.
\n--related is skipped for a rather small sample size.
\n--kinship analysis carried out instead for such a small sample size.
R plot for --kinship / --ibs / --unrelated / --tdt / --bySNP / --bysample / --homog / --gdt /
  --herit / --lmm / --makegrm / --pcgdt / --popgdt is not available.
```

---

## 11. Complete output-suffix inventory (all analyses)

Concatenated onto `<prefix>` with no separator unless the literal itself starts with `_`.

**Relatedness / QC:** `.kin`, `.kin0`, `X.kin`, `X.kin0`, `$$$.kin0`, `.seg`, `X.seg`,
`allsegs.txt`, `.ibs`, `.ibs0`, `$$$.ibs0`, `.con`, `.roh`, `.por`, `.dis`,
`unrelated.txt`, `unrelated_toberemoved.txt`, `cluster.kin`, `cluster`, `updateids.txt`,
`updateparents.txt`, `build.log`, `build.king`, `build`, `splitped.txt`, `_relatives.txt`,
`bySample.txt`, `bySNP.txt`, `_autoQC_Summary.txt`, `_autoQC_snptoberemoved.txt`,
`_autoQC_sampletoberemoved.txt`, `_autoQC_updatesex.txt`, `het.txt`, `af.dat`, `af.ped`,
`wipe`, `wipe.bgeno`, `.ld`, `_relative_removed.txt`.

**Mapping / association / structure:** `.her`, `.ih2`, `.mi`, `.ibdmap`, `.aucmap`, `.npl`,
`.popibd`, `.ibdgdt`, `.dst`, `.anc`, `.homomap`, `.homomapMH`, `.mthomo`, `.poproh`,
`.rohdiff`, `_eigenvalue.txt`, `_eigenvalue_grm-lmm.txt`, `_pc.txt`, `_lmmpc.txt`, `pc.txt`,
`_popref.txt`, `_grm.txt`, `_grm_grm-lmm.txt`, `_Dist.txt`, `tdt.txt`, `TDT.info`,
`TDTuninfo.txt`, `_gdt.txt`, `_gdt_ped.txt`, `_gdt_ibdseg.txt`, `_gdt_kinship.txt`,
`_gdt_lmm.txt`, `_lmm_…​.txt`, `_lmmking_disease.txt`, `_lmm_disease.txt`, `_novclmm.txt`,
`_linear.txt`, `_poodt.txt`, `grs.txt`, `hap.dat`, `hap.ped`, `flip.txt`, `_pat.{fam,bim,bed}`,
`_mat.{fam,bim,bed}`, `$$$.bed`, `.fam`, `.bim`, `.bed`, `.cov`, `.phe`, `.dat`, `.ped`,
`.map`, `.king`, `_InferredAncestry.txt`.

**R plot scripts / images:** `_relplot.R/.ps/.pdf`, `_ibd1vsibd2.R/.ps/.pdf`,
`_duplicateplot.R/.ps/.pdf`, `_clusterplot.R/.ps/.pdf`, `_pedplot.R/.ps/.pdf`,
`_buildplot.R/.ps/.pdf`, `_MIerrorplot.R/.ps/.pdf`, `_rohplot.R/.ps/.pdf`,
`_ancestryplot.R/.ps/.png/.pdf`, `_uniqfamplot`, `_uniqfam0plot`, `_herplot.R`,
`_poprohplot.R`, `_mthomoplot.R`, `_popdistplot.R`, `_ibdmapplot.R`, `_nplplot.R`,
`_aucmapplot.R`, `_%splot.ps`.
Driver strings: `%s CMD BATCH %s`, `R CMD BATCH %s`, `ps2pdf -sPAPERSIZE=letter %s_….ps`,
`  Please rerun R code %s (or rerun KING) after ggplot2 is installed.\n`,
`--related is done but R code %s failed.\n\n`.

---

## 12. Non-relatedness headers (for completeness — one line each)

| Analysis | File | Header |
|---|---|---|
| `--bysample` | `bySample.txt` | `FID IID FA MO SEX N_SNP Missing Heterozygosity` [` N_xSNP xHeterozygosity`][` N_ySNP N_yHetero`][` N_mtSNP N_mtHetero`][` N_pair N_MIp Err_MIp`][` N_trio N_MIt Err_MIt`][` MI_Removal`] |
| `--bySNP` | `bySNP.txt` | `SNP` + ` Chr` + ` Pos` + ` Label_A Label_a Freq_A N N_AA N_Aa N_aa CallRate` [` N_MZ N_HetMZ N_errMZ Err_InMZ Err_InHetMZ`][` N_PO N_HomPO N_errPO Err_InPO Err_InHomPO`][` N_trio N_HetOff N_errTrio Err_InTrio Err_InHetTrio`] |
| LD | `.ld` | `CHR_A\tBP_A\tSNP_A\tAF_A\tCHR_B\tBP_B\tSNP_B\tAF_B\tN\tHom\tD\tD_Prime\tR2\n` |
| `--makeGRM` | `_grm.txt` | `FID1\tIID1\tFID2\tIID2\tGRM_Raw\tGRM\n`, row `%s\t%s\t%s\t%s\t%.4lf\t%.4lf\n` |
| `--tdt` | `tdt.txt` | `SNP` `\tChr` `\tPos` `\tA1\tA2\tT\tNT\tOR\tChisq\tP\n`, row `%s\t%d\t%d\t%s\t%s\t%d\t%d\t%.3lf\t%.3lf\t%.3G\n` |
| `--gdt` | `_gdt.txt` | `SNP\tChr\tPos\tAllele1\tAllele2\tFrqUnaf\tN\tBeta_GT\tSE_GT\tZ\tPvalue\n` |
| `--lmm` | `_lmm_….txt` | `SNP\tChr\tPos\tAllele1\tAllele2\tFreq\tN\tBeta\tSE\tH2\tZ\tPvalue\n` |
| `--risk` | `grs.txt` | `%15s %15s %7s %7s %7s %7s %7s %9s` over `FID IID InfSNP InfVar GRS Zscore Percent ScaledGRS` |
| `--HEreg` | `.her` | `Chr\tPos\tFlankSNP1\tFlankSNP2\tN_IBD0\tN_IBD1\tN_IBD2\tSqDiff0\tSqDiff1\tSqdiff2\th2\tAlpha\tNegBeta\tSE\tT\tPvalue\tLOD\n` *(note `Sqdiff2` lower-case d — sic)* |
| `--ibdmap` | `.ibdmap` | `Chr\tPosMb\tFlank1\tFlank2\tPI_URP\tPI_DRP\tPI_ARP\tPI_Diff\tP\n` |
| `--aucmap` | `.aucmap` | `Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_URP\tPI_DRP\tPI_ARP\tNMISS\tN_URP\tN_DRP\tN_ARP\tSuccess\tAUC\n` |
| `--npl` | `.npl` | `Chr\tPos\tFlankSNP1\tFlankSNP2\tPI_USP\tPI_DSP\tPI_ASP\tPI_AP\tLOD_Raw\tLOD_ASP\tLODwDSP\n` (or the 8-col ASP-only variant) |
| `--ibdH2` | `.ih2` | `Chr\tPos\tFlankSNP1\tFlankSNP2\tTrait\tN_Obs\tN_Pred\tR2\tH2\n` |
| `--ibdMI` | `.mi` | `SNP\tChr\tPos\tN_IBD\tN_Inf\tN_MI\tR_InfMI\tRate_MI\n` |
| `--popdist` | `.dst` | `POP1\tPOP2\tN\tPropIBD\tDistIBD\n` |
| `--ancestry` | `.anc` | `FID\tIID\tNMISS\tNMIX\tN1\tN2\tAnc_P1\tAnc_P2\tAdmix\tAncestry\n` |
| `--homomap` | `.homomap` | `Chr\tPos\tFlankSNP1\tFlankSNP2\tN_Con\tHomoCon\tN_Cas\tHomoCas\tOR\tSE\tLOD\n` |
| `--paternity` | `_relatives.txt` | `FID\tIID\tFA\tMO\tSEX\tMaxXIBD\tInfType\tRelatives\n`, row `%s\t%s\t%s\t%s\t%d\t%.4lf\t%s\t` then comma list via `%s,` |
| `--mds` projection | `_Dist.txt` | `FID IID HetProj HetRef MinDist Kinship Closest RefID`, row `%s %s %.4lf %.4lf %.4lf %.4lf %d %s` |
| `--pca` | `_pc.txt` | `FID IID FA MO SEX AFF` + ` PC%d`…, values ` %.4lf`; eigenvalues ` %.5lf\n` |

---

## 13. Reproduction recipe

```bash
KING="/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
strings -n 3 "$KING" > king.strings.txt
otool -tV "$KING" > king.disasm.txt          # otool annotates: literal pool for: "…"
nm -U "$KING" | c++filt                      # full C++ symbol table — every writer is named
```
The synthetic PLINK generator used for the live runs is inline in the session log; the key
constraint is that **A1 must be the minor allele** or KING aborts with
`Too many first alleles as the major allele (~%.1lf%%)…`.

---

## 14. Open items for other recon streams

1. `.kin`'s within-family row order is alphabetical-by-IID (verified), but I did not confirm
   whether the comparison is byte-wise `strcmp` or locale-aware, nor the FID tie-break rule.
2. The exact predicate that gates `X.seg` emission under `--ibdseg` (present with `--degree 4`,
   absent without) is unresolved.
3. `Error` column semantics (`0` / `0.5` / `1`) — needs a Mendelian-error fixture to pin down
   when `0.5` is emitted.
4. Whether we reproduce §6.1 (`X.seg` trailing tab + wrong header) and §6.2 (`X.kin0` race) or
   emit the corrected forms is a product decision, not a recon one.
