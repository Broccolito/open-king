# open-king

`open-king` reads a PLINK1 fileset and reports how every pair of samples in it is related:
kinship coefficients, duplicate and monozygotic pairs, identity-by-descent segments, a
relationship label for each pair, and per-sample and per-marker quality control. Results are
written as plain text files.

**open-king** is a clean-room, MIT-licensed Rust implementation of the relatedness and QC
core of [KING 2.3.2](https://www.kingrelatedness.com/). It takes the same command line and
the same input files, and it writes byte-identical output. All 480 captured reference
invocations reproduce exactly, across 876 output files.

**Documentation: <https://broccolito.github.io/open-king/>**

Every option with a runnable example, the output-file reference, installation for each
platform, the parity evidence, and the benchmarks.

## Why a rewrite

KING was published in 2010 and the program has accumulated since. The estimators are sound
and widely relied on, and this project does not try to improve them: it reproduces them to
the byte. What a modern implementation changes is everything around the arithmetic.

### It can be redistributed

open-king is [MIT](LICENSE) licensed. It can be vendored into a pipeline, built into a
container image, packaged by a distribution, shipped inside a commercial product, forked and
modified. KING's terms are "Feel free to use KING for your research, but please do not
redistribute AND make profits", which leaves each of those uses to be negotiated. For a
workflow that has to be handed to somebody else and run without a manual download step, the
licence is often the deciding constraint rather than any technical one.

### It runs where you put it

The binary is 1.1 MB and links one library, the operating system's own:

```
$ otool -L open-king
open-king:
	/usr/lib/libSystem.B.dylib
```

KING's published macOS binary links `libstdc++`, `libgomp` and `libgcc_s` from a Homebrew
GCC 9 installed under `/usr/local`, so on an Apple silicon machine it exits before it starts:

```
Library not loaded: /usr/local/opt/gcc/lib/gcc/9/libgomp.1.dylib
```

open-king ships prebuilt, notarized binaries for macOS arm64, macOS x86_64, Linux x86_64 and Windows
x86_64, each with a published SHA-256, and needs no compiler, no runtime and no `brew
install` on any of them.

### The same input gives the same bytes

Run KING twelve times on one fileset and hash the `X.kin0` it writes each time, and you get
more than one answer: several threads append to a single unlocked file handle, so records
interleave differently from run to run. The same experiment on open-king gives one hash
twelve times.

| | distinct `X.kin0` checksums over 12 identical runs |
| --- | --- |
| open-king | 1 |
| KING 2.3.2 | 2 |

Reproducibility is a property a pipeline can depend on. Related to it, safe Rust does not
read uninitialized memory, which KING does in three places this project measured: its own
startup banner, the tail of a marker array that is an exact multiple of 64, and the value
carried by `--noscreen`, whose effect consequently differs between two builds of the same
KING source. Those are documented in [`docs/PARITY.md`](docs/PARITY.md) §5.2, §5.3 and §5.13.

### Where it is faster

The benchmark measures 104 cells, eight analyses across thirteen filesets. On the ten
reference datasets, open-king finishes first in 72 of those 80 corpus cells, at a median of
3.7 times. On a trio it is 0.004 s against 0.015 s. Peak memory is 3.8 MB against 8.6 MB.
Small runs are dominated by process startup, and open-king starts faster and holds less.

At cohort scale the picture divides. The remaining 24 cells are three synthetic cohorts of
200, 400 and 800 samples at 100,000 markers. Counting analyses stay level: on 800 samples,
`--kinship` costs 1.43 CPU seconds against KING's 1.53, and `--ibs` 1.77 against 1.95.
Analyses that build IBD segments do not stay level, and the gap grows with sample count. At
those same 800 samples, in CPU seconds, open-king against KING: `--ibdseg` 5.13 against 0.74,
`--related` 11.12 against 0.38, `--build` 20.34 against 0.46, `--cluster` 19.76 against 0.36,
`--unrelated` 19.92 against 0.36. That is 7 to 55 times more CPU on the segment work. It is
measured, tracked as [issue #12](https://github.com/Broccolito/open-king/issues/12), and
reported in full on the
[benchmarks page](https://broccolito.github.io/open-king/benchmarks.html) rather than
averaged away.

The line falls between the two. On the nine corpus filesets of 40 samples or fewer, open-king
wins all 72 cells. The eight it loses are the whole of `bigish`, at 200 samples, and from 200
samples up KING has the lower wall clock in all 32 cells measured, `--duplicate` included at
1.29 CPU seconds against 0.25. `--kinship` and `--ibs` are the two that stay level on CPU
work at every size tested.

### It is a library as well as a program

`open-king-core` exposes the engine as typed Rust data, so a program that needs kinship does not
have to spawn a process and parse a report. See [`docs/API.md`](docs/API.md).

### It is tested in the open

355 unit and integration tests, plus a differential suite that replays all 480 captured
reference invocations and compares every output file byte for byte. The suite is committed,
needs no reference binary, runs in about four seconds, and runs in CI on Linux, macOS and
Windows on every push. It fails on any difference in either direction, so an unrecorded
improvement is a failure too.

## Provenance and license

open-king is **not affiliated with, endorsed by, or derived from the source code of** the
original KING program by Wei-Min Chen and colleagues. It was written from the published
algorithm descriptions, the publicly documented file formats, and black-box observation of
the reference binary. **No KING source code was read or copied.** The rule is stated and
justified in [`docs/MAINTAINING.md`](docs/MAINTAINING.md) §1.

open-king's own code is **MIT** licensed ([LICENSE](LICENSE)). That covers this code only and
makes no claim about the original KING, which its authors license separately.

## Install

Prebuilt binaries for macOS (arm64 and x86_64), Linux x86_64 and Windows x86_64 are attached
to each [release](https://github.com/Broccolito/open-king/releases). Download, check the
SHA-256 and extract:

```bash
curl -LO https://github.com/Broccolito/open-king/releases/download/v0.1.0/open-king-macos-arm64.zip
curl -LO https://github.com/Broccolito/open-king/releases/download/v0.1.0/SHA256SUMS.txt
shasum -a 256 -c SHA256SUMS.txt --ignore-missing
unzip open-king-macos-arm64.zip
./open-king
```

Each archive holds one file, the `open-king` executable. There is no directory to descend
into and nothing to tidy up afterwards.

`--ignore-missing` is needed because `SHA256SUMS.txt` covers all four archives and you
downloaded one.

**The macOS builds are signed and notarized by Apple**, so they run straight out of the
archive. There is no Gatekeeper warning and no `xattr -d com.apple.quarantine` step.

To build instead, a Rust toolchain (1.75 or newer) is the only requirement. Python 3,
standard library only, runs the test corpus and the parity suite.

```bash
git clone https://github.com/Broccolito/open-king
cd open-king
cargo build --release
```

The binary lands at **`target/release/open-king`**. Put it on your `PATH`, or call it by path.

`open-king` prints KING 2.3.2's own banner, deliberately, so that scripts reading stdout keep
working. That means the banner does not identify which of the two programs is running: check
the path you invoked, or the binary's checksum. The
[install page](https://broccolito.github.io/open-king/install.html) covers platform details
and setup for an automated agent.

## 60-second quickstart

You need a PLINK1 fileset. If you do not have one to hand, the repository can synthesise 13
of them, with real pedigrees and real recombination, using no external tools:

```bash
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
```

Then find the relatives, out to second degree:

```bash
target/release/open-king -b /tmp/kingdocs/bigish.bed --related --degree 2 --prefix demo
```

```
Options in effect:
	--related
	--degree 2
	--prefix demo

Total length of 22 chromosomal segments usable for IBD segment analysis is 2498.9 Mb.
  Information of these chromosomal segments can be found in file demoallsegs.txt

Within-family kinship data saved in file demo.kin

Relationship summary (total relatives: 436 by pedigree, 435 by inference)
  Source	MZ	PO	FS	2nd	3rd	OTHER
  ===========================================================
  Pedigree	0	226	111	81	18	137
  Inference	0	226	111	79	19	138
...
  Final Stage (with 50000 SNPs): 26 pairs of relatives (up to 2nd-degree) are confirmed

Between-family relatives (kinship >= 0.08839) saved in file demo.kin0
```

Three files: `demo.kin` (573 within-family pairs), `demo.kin0` (26 cross-family pairs that
cleared the 2nd-degree cut), `demoallsegs.txt` (the usable-genome denominator).

```bash
head -3 demo.kin0
```

```
FID1	ID1	FID2	ID2	N_SNP	HetHet	IBS0	HetConc	HomIBS0	Kinship	IBD1Seg	IBD2Seg	PropIBD	InfType
BF01	B01_F	BF02	B02_F	50000	0.2274	0.0132	0.4855	0.1175	0.2885	0.4575	0.3676	0.5964	FS
BF01	B01_F	BF02	B02_C1	50000	0.1560	0.0299	0.2873	0.2261	0.1368	0.5684	0.0000	0.2842	2nd
```

`B01_F` and `B02_F` are declared in different families and share half their genome with a
quarter of it doubled: full siblings the pedigree never mentioned. That is the finding
`--related` exists to produce.

Two things to know before your first real run. **`--prefix` is concatenated, not joined**
(`--prefix demo` gives `demo.kin` and `demoallsegs.txt`). And **`--degree` is what bounds
`--kinship`**: without it, `.kin0` holds every between-family pair, which on this
200-sample fileset is 19,327 rows against 24 at `--degree 2`. `--related` runs the other
way and defaults to first degree, so leaving `--degree` off narrows its `.kin0` rather than
widening it.

Sample identity is the `(FID, IID)` pair under ASCII case-folding. For example, `A_F` and
`a_f` in the same family collide and the fileset is rejected, matching KING. Make those
keys case-insensitively unique before running an analysis.

## The analyses

One or more of these must be given, or nothing is computed. Full reference, including all 46
options and the parser behaviour that is not obvious from the flag name:
[the commands page](https://broccolito.github.io/open-king/usage.html) or
[`docs/CLI.md`](docs/CLI.md).

| flag | what it does | writes |
| --- | --- | --- |
| [`--kinship`](docs/CLI.md#--kinship) | KING-robust kinship for every pair; no segments, no labels | `.kin`, `.kin0`, `X.kin`, `X.kin0` |
| [`--related`](docs/CLI.md#--related) | kinship **plus** IBD segments and a relationship label per pair | `.kin`, `.kin0`, `X.kin`, `X.kin0`, `allsegs.txt` |
| [`--duplicate`](docs/CLI.md#--duplicate) | duplicate and MZ pairs by heterozygote concordance | `.con` |
| [`--ibdseg`](docs/CLI.md#--ibdseg) | pairwise IBD-segment inference on its own | `.seg`, `X.seg`, `allsegs.txt`, `splitped.txt` |
| [`--ibs`](docs/CLI.md#--ibs) | full IBS and concordance counts for every pair | `.ibs`, `.ibs0`, `allsegs.txt` |
| [`--unrelated`](docs/CLI.md#--unrelated) | a maximal mutually unrelated subset, and its complement | `unrelated.txt`, `unrelated_toberemoved.txt` |
| [`--cluster`](docs/CLI.md#--cluster) | merge families connected by inferred relatedness | `cluster.kin`, `updateids.txt` |
| [`--build`](docs/CLI.md#--build) | reconstruct pedigrees from the genotypes | `updateparents.txt`, `updateids.txt`, `build.log` |
| [`--bysample`](docs/CLI.md#--bysample) | per-sample QC: call rate, heterozygosity, Mendelian errors | `bySample.txt` |
| [`--bySNP`](docs/CLI.md#--bysnp) | per-marker QC: frequency, genotype counts, error rates | `bySNP.txt` |
| [`--autoQC`](docs/CLI.md#--autoqc) | the packaged call-rate and sex QC pipeline | four `_autoQC_*` reports |

The common modifiers are `--degree <d>`, `--prefix`, `--seglength <Mb>`, `--minConc <x>` and
`--cpus <n>`. `--cpus` changes no printed digit in any output file; it is retained for
command-line compatibility and console reporting rather than as a guaranteed worker cap.

**Deliberately excluded:** `--pca`, `--mds`, `--roh`, `--makeGRM`, `--plink`, `--lmm`,
`--tdt`, `--gdt`, `--risk`, R plotting and comma-separated multi-fileset merging are not part
of this minimal relatedness package. Their parser spellings are retained for banner
compatibility, and requesting one is a fatal error raised before any input file is opened, so
a workflow can never appear to succeed without producing what it asked for. See
[the product-scope contract](docs/SCOPE.md) and use a dedicated tool for those analyses.

## Rust library API

`open-king-core` exposes the relatedness engine as owned Rust data. `Bundle::from_plink` loads and
validates a PLINK1 fileset; `Bundle::relatedness` returns every unordered pair's exact
counts, estimators, pedigree expectations, relationship class, and optional IBD segment
metrics, without invoking the CLI or parsing a report:

```rust
use open_king_core::{Bundle, BundleError, RelatednessOptions};

fn main() -> Result<(), BundleError> {
    let bundle = Bundle::from_plink("cohort.bed")?;
    let report = bundle.relatedness(&RelatednessOptions::default());

    for pair in &report.pairs {
        let first = &report.samples[pair.first];
        let second = &report.samples[pair.second];
        println!(
            "{}/{} {}/{} {:.6} {}",
            first.fid,
            first.iid,
            second.fid,
            second.iid,
            pair.statistics.kinship,
            pair.relationship.label(),
        );
    }
    Ok(())
}
```

See [`docs/API.md`](docs/API.md) for field semantics, validation, calibration, and the cost
of all-pairs segment scanning.

## Parity

**All 480 captured reference invocations reproduce byte for byte**, covering every output
file, plus stdout, stderr and exit status. The run compares 876 output files; eight
documented files whose reference bytes are unstable on the host are excluded symmetrically.
[`docs/PARITY.md`](docs/PARITY.md) is the authoritative statement and quantifies the
differences the 480-case corpus cannot exercise.

Byte-identical wherever the corpus produces them: `--kinship` including the X pass,
`--duplicate`, `--ibs`, `--unrelated`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster`,
`--ibdseg`, `--related`, and `--build` at every captured shape and reporting floor, across
all 23 output files the project writes. `--build` reproduces the primary `bigish`
case in full, including stdout and all four reconstruction files, covering FS0/FS1/FS2, PO.S
orientation, AV.FS/AV.HS/HS.UN2, and the unstable sibling order KING inherits from
libStatGen.

**Differences the corpus cannot see** are quantified rather than assumed. An independent
24-fileset segment battery is 68/72 whole-run exact with 0 extra and 0 missing rows across
6,713 reference rows. Its four value differences occur only on exact-40,000-marker inputs
and trace to KING's uninitialised multiple-of-64 tail read; the 39,999 and 40,001 controls
are exact, and safe Rust does not reproduce undefined behaviour. Five held-out items remain
open and are tracked under
[issue #11](https://github.com/Broccolito/open-king/issues/11): a segment acceptance gate,
the data-derived sparse PO/FS cutoff, `HomIBS0` tie rendering, `MI_Removal`, and unusual
pedigree reconstruction shapes.

Every number above is measured against one reference build: KING 2.3.2, Mach-O arm64, macOS.
KING's segment algorithm is unpublished and its release notes record repeated changes across
2.1.x and 2.2.x, so "byte-identical" means to 2.3.2.

### Run the parity suite yourself

It needs no reference binary, because the goldens are committed, and it takes about two
seconds:

```bash
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
python3 tests/parity/run_parity.py --impl target/release/open-king -q
```

```
[parity] 480 case(s), impl=target/release/open-king, jobs=8
parity: 480 PASS, 0 FAIL, 480 total (876 output file(s) byte-compared, 8 diff-excluded)
```

To diff the two binaries on your own data instead, recipe 12 of
[`docs/COOKBOOK.md`](docs/COOKBOOK.md) has the procedure and the console normalizer it needs.

## Repository documentation

The published site at <https://broccolito.github.io/open-king/> is the place to start. The
repository carries the longer material behind it.

**Using it**

| | |
| --- | --- |
| [`docs/CLI.md`](docs/CLI.md) | every option, what it affects, and how the parser really behaves |
| [`docs/OUTPUTS.md`](docs/OUTPUTS.md) | every output file: columns, formats, row order, and when it is absent |
| [`docs/COOKBOOK.md`](docs/COOKBOOK.md) | twelve task-oriented recipes, from finding duplicates to diffing against KING |
| [`docs/INTERPRETING.md`](docs/INTERPRETING.md) | what the numbers mean, where they mislead, and what they cannot tell you |
| [`docs/SCOPE.md`](docs/SCOPE.md) | the supported core and the deliberately excluded analyses |

**Working on it**

| | |
| --- | --- |
| [`docs/PARITY.md`](docs/PARITY.md) | the measured parity claim, per file and per row |
| [`docs/SPEC.md`](docs/SPEC.md) | the implementation specification |
| [`docs/BEHAVIOR.md`](docs/BEHAVIOR.md) | the black-box experiments that fixed each rule |
| [`docs/VERIFIED_FORMULAS.md`](docs/VERIFIED_FORMULAS.md) | every estimator, checked numerically against the reference |
| [`docs/MAINTAINING.md`](docs/MAINTAINING.md) | the clean-room rule, the corpus, re-capturing goldens |

[`docs/README.md`](docs/README.md) indexes all of it, including the 25 research notes and the
continuation brief.

## Citation

**Cite the method, not this implementation.** The estimators are Manichaikul *et al.*:

> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. Robust relationship
> inference in genome-wide association studies. *Bioinformatics*. 2010;26(22):2867-2873.
> doi:[10.1093/bioinformatics/btq559](https://doi.org/10.1093/bioinformatics/btq559)

**Credit for KING itself belongs to Wei-Min Chen and colleagues** (University of Virginia),
whose program at <https://www.kingrelatedness.com/> this project reimplements and measures
itself against. The IBD-segment algorithm in particular is theirs and is unpublished;
open-king recovered its behaviour by observation, not by reading their work.

If a paper needs to record which implementation produced a number,
[`CITATION.cff`](CITATION.cff) carries machine-readable metadata for open-king alongside both
references above.
