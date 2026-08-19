# open-king

`king` estimates how every pair of samples in a PLINK1 fileset is related — kinship
coefficients, duplicate and MZ pairs, IBD segments, and a relationship label per pair — and
writes the answer as plain text files. **open-king** is a clean-room, MIT-licensed Rust
reimplementation of the core relatedness and QC workflows in
[KING 2.3.2](https://www.kingrelatedness.com/). On its supported surface it aims for the
same command line, the same PLINK1 input and byte-identical output files. It deliberately
does not reproduce KING's population-structure, association, risk, ROH or R-plotting
toolbox; [`docs/SCOPE.md`](docs/SCOPE.md) defines that product boundary.

## Provenance and license

open-king is **not affiliated with, endorsed by, or derived from the source code of** the
original KING program by Wei-Min Chen and colleagues. It was written from the published
algorithm descriptions, the publicly documented file formats, and black-box observation of
the reference binary. **No KING source code was read or copied.** The rule is stated and
justified in [`docs/MAINTAINING.md`](docs/MAINTAINING.md) §1.

open-king's own code is **MIT** licensed ([LICENSE](LICENSE)). That covers this code only and
makes no claim about the original KING, which its authors license separately.

## Install

A Rust toolchain (1.75+) is the only requirement. Python 3, standard library only, runs the
test corpus and the parity suite.

```bash
git clone https://github.com/Broccolito/open-king
cd open-king
cargo build --release
```

The binary lands at **`target/release/king`**. Put it on your `PATH`, or call it by path.

## 60-second quickstart

You need a PLINK1 fileset. If you do not have one to hand, the repository can synthesise 13
of them — real pedigrees, real recombination, no external tools:

```bash
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
```

Then find the relatives, out to second degree:

```bash
target/release/king -b /tmp/kingdocs/bigish.bed --related --degree 2 --prefix demo
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
quarter of it doubled — full siblings the pedigree never mentioned. That is the finding
`--related` exists to produce.

Two things to know before your first real run: **`--prefix` is concatenated, not joined**
(`--prefix demo` gives `demo.kin` *and* `demoallsegs.txt`), and on a cohort of any size
**always pass `--degree`** — unfiltered, `.kin0` is every pair in the dataset.

Sample identity is the `(FID, IID)` pair under ASCII case-folding. For example, `A_F` and
`a_f` in the same family collide and the fileset is rejected, matching KING. Make those
keys case-insensitively unique before running an analysis.

## The analyses

One or more of these must be given, or nothing is computed. Full reference, including all 46
options and the parser's several surprises: [`docs/CLI.md`](docs/CLI.md).

| flag | what it does | writes |
| --- | --- | --- |
| [`--kinship`](docs/CLI.md#--kinship) | KING-robust kinship for every pair; no segments, no labels | `.kin`, `.kin0`, `X.kin`, `X.kin0` |
| [`--related`](docs/CLI.md#--related) | kinship **plus** IBD segments and a relationship label per pair | `.kin`, `.kin0`, `X.kin`, `allsegs.txt` |
| [`--duplicate`](docs/CLI.md#--duplicate) | duplicate and MZ pairs by heterozygote concordance | `.con` |
| [`--ibdseg`](docs/CLI.md#--ibdseg) | pairwise IBD-segment inference on its own | `.seg`, `X.seg`, `allsegs.txt`, `splitped.txt` |
| [`--ibs`](docs/CLI.md#--ibs) | full IBS and concordance counts for every pair | `.ibs`, `.ibs0`, `allsegs.txt` |
| [`--unrelated`](docs/CLI.md#--unrelated) | a maximal mutually unrelated subset, and its complement | `unrelated.txt`, `unrelated_toberemoved.txt` |
| [`--cluster`](docs/CLI.md#--cluster) | merge families connected by inferred relatedness | `cluster.kin`, `updateids.txt` |
| [`--build`](docs/CLI.md#--build) | reconstruct pedigrees from the genotypes | `updateparents.txt`, `updateids.txt`, `build.log` |
| [`--bysample`](docs/CLI.md#--bysample) | per-sample QC: call rate, heterozygosity, Mendelian errors | `bySample.txt` |
| [`--bySNP`](docs/CLI.md#--bysnp) | per-marker QC: frequency, genotype counts, error rates | `bySNP.txt` |
| [`--autoQC`](docs/CLI.md#--autoqc) | the packaged call-rate and sex QC pipeline | four `_autoQC_*` reports |

The common modifiers are `--degree <d>` (report only relatives that close), `--prefix`,
`--seglength <Mb>`, `--minConc <x>` and `--cpus <n>`. `--cpus` changes no printed digit in
any output file; it is retained for command-line compatibility and console reporting, not
as a guaranteed Rayon thread cap.

**Deliberately excluded:** `--pca`, `--mds`, `--roh`, `--makeGRM`, `--plink`, `--lmm`,
`--tdt`, `--gdt`, `--risk`, R plotting and comma-separated multi-fileset merging are not
part of this minimal relatedness package. Their parser spellings are retained for banner
compatibility but do not run an analysis. See [the product-scope contract](docs/SCOPE.md)
and use a dedicated tool for those workflows. Until unsupported-option diagnostics land,
assert on expected output files rather than exit status alone.

## Parity, honestly

**477 of the 480 captured reference invocations are byte-identical** — every output file,
plus stdout, stderr and exit status. [`docs/PARITY.md`](docs/PARITY.md) is the authoritative
statement: the full analysis × dataset matrix, per-file and per-row scorecards, and a labelled
limitations section. Everything here is a summary of it.

Byte-identical everywhere the corpus produces them: `--kinship` (including the X pass),
`--duplicate`, `--ibs`, `--unrelated`, `--bysample`, `--bySNP`, `--autoQC`, `--cluster`,
`--ibdseg` and `--related` at all three captured `--seglength` floors — 30 of the 31 output
files this project writes.

**The three failing cases, and their blast radius:**

| cases | what differs | does it affect an output file? |
| ---: | --- | --- |
| 2 | one stdout line — `--related`'s two-stage screening count on the 200-sample dataset (`36` vs `50`) | **No.** `.kin`, `.kin0` and `.seg` are byte-identical; the rows come from the exhaustive re-estimate below that line |
| 1 | `<prefix>build.log`'s `INFERENCE` half | **Yes, that one file.** Its header and `RULE` lines are byte-identical and the file is a strict subsequence of the reference's. `--build`'s `updateids.txt` and `updateparents.txt` are byte-identical |

**Differences the corpus cannot see** cost no case but a user can still hit them. The
`--related` path now detects a marker panel too sparse for the segment caller and switches to
the reference's short, kinship-only output. The same fallback now drives `--unrelated`,
`--cluster` and `--build`: every file from those three analyses is byte-identical on the
held-out sparse fixture. The shared two-stage screen still admits two extra between-family
`--related` candidates (17 rather than 15), so its console summary and `.kin0` row set are
not yet exact. The A1-major input gate, unsorted-map rejection and case-insensitive sample-ID
validation are implemented and have focused differential probes. These are measured in
[`docs/PARITY.md`](docs/PARITY.md) §5.10–§5.12 and §4.6; case-only sample-ID collisions are
now rejected like KING. The independent 24-fileset segment battery is **68/72** whole-run
exact with **0 extra and 0 missing rows** across 6,713 reference rows. Its four value
differences occur only on exact-40,000-marker inputs and are KING's uninitialised
multiple-of-64 tail read; 39,999- and 40,001-marker controls are exact, so safe Rust does
not emulate that undefined behavior.

Every number above is measured against **one** reference build: KING 2.3.2, Mach-O arm64,
macOS. KING's segment algorithm is unpublished and its release notes record repeated changes
across 2.1.x–2.2.x, so "byte-identical" means *to 2.3.2*.

### Run the parity suite yourself

It needs no reference binary — the goldens are committed — and takes about four seconds:

```bash
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs
python3 tests/parity/run_parity.py --impl target/release/king -q
```

```
[parity] 480 case(s), impl=/Users/wgu/Desktop/open-king/target/release/king, jobs=8
FAIL  apps/bigish__build                          stdout!=; kingbuild.log!=(num)
FAIL  core/bigish__related_degree2                stdout!=
FAIL  ibdseg/bigish__related_degree2_ibdseg       stdout!=

parity: 477 PASS, 3 FAIL, 480 total (2.1s wall, 876 output file(s) byte-compared, 8 diff-excluded)
```

(The three `FAIL` rows can arrive in any order — cases run in parallel — and the wall-clock
figure varies.) To diff the two binaries on your own data instead, recipe 12 of
[`docs/COOKBOOK.md`](docs/COOKBOOK.md) has the procedure and the console normalizer it needs.

`cargo test --workspace` is 330 tests. CI replays all 480 captures against a committed
baseline on every push and fails on any difference **in either direction**, so an unrecorded
improvement is a failure too.

## Documentation

**Using it**

| | |
| --- | --- |
| [`docs/CLI.md`](docs/CLI.md) | every option, what it affects, and how the parser really behaves |
| [`docs/OUTPUTS.md`](docs/OUTPUTS.md) | every output file: columns, formats, row order, and when it is absent |
| [`docs/COOKBOOK.md`](docs/COOKBOOK.md) | twelve task-oriented recipes, from "find duplicates" to "diff against KING" |
| [`docs/INTERPRETING.md`](docs/INTERPRETING.md) | what the numbers mean, where they mislead, and what they cannot tell you |
| [`docs/SCOPE.md`](docs/SCOPE.md) | the supported minimal core and deliberately excluded legacy/non-core analyses |

**Working on it**

| | |
| --- | --- |
| [`docs/PARITY.md`](docs/PARITY.md) | the measured parity claim, per file and per row |
| [`docs/SPEC.md`](docs/SPEC.md) | the implementation specification |
| [`docs/BEHAVIOR.md`](docs/BEHAVIOR.md) | the black-box experiments that fixed each rule |
| [`docs/VERIFIED_FORMULAS.md`](docs/VERIFIED_FORMULAS.md) | every estimator, checked numerically against the reference |
| [`docs/MAINTAINING.md`](docs/MAINTAINING.md) | the clean-room rule, the corpus, re-capturing goldens |

[`docs/README.md`](docs/README.md) indexes all of it, including the 25 research notes and the handoff brief.

## Citation

**Cite the method, not this reimplementation.** The estimators are Manichaikul *et al.*:

> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. Robust relationship
> inference in genome-wide association studies. *Bioinformatics*. 2010;26(22):2867–2873.
> doi:[10.1093/bioinformatics/btq559](https://doi.org/10.1093/bioinformatics/btq559)

**Credit for KING itself belongs to Wei-Min Chen and colleagues** (University of Virginia),
whose program — <https://www.kingrelatedness.com/> — this project reimplements and measures
itself against. The IBD-segment algorithm in particular is theirs and is unpublished;
open-king recovered its behaviour by observation, not by reading their work.

If a paper needs to say which implementation produced a number, [`CITATION.cff`](CITATION.cff)
carries machine-readable metadata for open-king alongside both references above.
