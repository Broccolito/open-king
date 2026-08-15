# open-king

A clean-room, MIT-licensed reimplementation of **KING** (Kinship-based INference for
Genome-wide association studies) — the relatedness-inference program described in
[Manichaikul *et al.* 2010, *Bioinformatics*](https://pmc.ncbi.nlm.nih.gov/articles/PMC3025716/).

KING answers one question about a genotype dataset: **for every pair of samples in it, how
are those two people related?** Given a PLINK1 fileset (`.bed`/`.bim`/`.fam`) it estimates
kinship coefficients from genome-wide SNPs, finds duplicate samples, calls the IBD segments
two people share, labels each pair `PO`, `FS`, `2nd`, `3rd`, `4th`, `Dup/MZ` or `UN`, flags
pairs whose inferred relationship contradicts the pedigree, reconstructs families, and picks
a maximal unrelated subset. It is a standard first step in GWAS quality control, and the
original implementation is a widely used binary with no source license permitting reuse.

open-king is that program rewritten in Rust from the published algorithms, the public
documentation and black-box observation of the reference binary — **no KING source code was
read or copied** — so that a permissively licensed implementation of these formats and
estimators can be embedded in other software.

The goal is **drop-in parity**: the same command line, the same input files, and
byte-identical output files as the reference `king` 2.3.2 binary. Not "statistically
equivalent" — the same bytes, including column widths, row order and the reference's own
rounding quirks. That is a falsifiable target, and the whole project is organised around
measuring it: 480 captured reference invocations across 13 synthetic datasets are replayed
on every change and diffed byte for byte.

```bash
cargo build --release
target/release/king -b study.bed --related --prefix study
```

Requires a Rust toolchain and nothing else. Python 3 (standard library only) runs the parity
suite.

## Status

**477 of the 480 captured reference invocations reproduce byte-identically (99.4 %)**,
including all 220 flag-plumbing and error probes. Run the suite yourself:

```bash
cargo build --release
python3 tests/parity/run_parity.py --impl target/release/king
```

`docs/PARITY.md` is the authoritative claim — the full analysis × dataset matrix, the
measured size of every remaining gap, and a labelled limitations section. Everything below
is a summary of it and says nothing it does not support.

The one-paragraph version: the relatedness estimators, the QC reports, duplicate detection,
auto-QC, unrelated-set selection, clustering, `--ibs`, the whole X-chromosome surface
(`X.kin`, `X.kin0`, `X.seg`), the whole command-line surface **and the IBD-segment engine at
every captured reporting floor** are byte-identical everywhere. **Thirty of the thirty-one
output files this project writes are byte-identical in every case that produces them**; only
`<prefix>build.log` differs anywhere, and only in 1 of the 8 cases that write it. On the
primary `--ibdseg` capture all **982**
rows are byte-exact on all four printed fields — `IBD1Seg`, `IBD2Seg`, `PropIBD` and
`InfType` — with **0 spurious and 0 missing rows on every output file in the corpus**.

The 3 cases that are not byte-identical fall under two named causes, and neither is a
segment estimate:

| cases | cause |
| ---: | --- |
| 2 | one stdout line: `--related`'s two-stage screening count on `bigish`. Every output file in both cases, `.kin`, `.kin0` and `.seg` alike, is byte-identical |
| 1 | `--build`'s `<prefix>build.log` — its `RULE` half lands, its `INFERENCE` half does not |

**Row-level, over the same 982 pairs, at all three captured floors:**

| `--seglength` | all four columns | `IBD1Seg` | `IBD2Seg` | extra | missing | MAE |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 3 Mb (default) | **982 / 982** | **982** | **982** | 0 | 0 | 0.000000 |
| 5 Mb | **982 / 982** | **982** | **982** | 0 | 0 | 0.000000 |
| 10 Mb | **982 / 982** | **982** | **982** | 0 | 0 | 0.000000 |

(`python3 tests/parity/fit/scorecard.py`.)

**Both of those numbers are saturated, and that is a fact about the corpus as much as about
the caller.** So this release also measures the segment engine on filesets it has never
seen — 24 of them, on 8 seeds used nowhere else in the repository, at all three floors,
byte-diffed against the reference:

```
72 runs, 6 713 reference rows   ->   66 byte-identical
rows: 0 extra, 2 missing, 4 value-differing
```

**66 of 72 runs and 6 707 of 6 713 rows.** The 6 wrong rows are one IBD2 deficit of
0.0182 on a full-sib pair (identical on two independent seeds — one missed piece in one
place) and one distant pair we drop that the reference reports. Both are written up with
their next experiment in `docs/PARITY.md` §4.6, and neither is visible to the 480 captures.
`python3 docs/research/fixtures/oosseg.py --ref /path/to/king` reproduces it.

**One caveat that applies to every number on this page.** All of it is measured against a
single reference build — `KING 2.3.2`, Mach-O arm64, on macOS. KING's own release notes
record repeated changes to the IBD-segment algorithm across 2.1.x–2.2.x, and that algorithm
has never been published, so "byte-identical" here means *to 2.3.2*. No cross-build or
cross-platform differential has been run. `docs/PARITY.md` opens with the full statement.

Two further differences sit outside the 480 captures entirely, so they cost no case but a
user could still hit them: `--ibdseg` does not apply the reference's 100 Mb usable-total
floor, and `splitped.txt` is written unconditionally. Both are measured and written up in
`docs/PARITY.md` §5.10.

### Two numbers, and why they are different

The headline is a **whole-file** count: a case turns `PASS` only when every row of every file
it writes is byte-exact. For most of this project's life the residual was spread thinly
across nearly every dataset, so it was routine for a large row-level gain to move the
headline by nothing — and, at the very end, for a change that moved no estimate at all to
move it by 28:

| change | row-level effect | headline |
| --- | --- | --- |
| the `.seg` IBD2 caller (`docs/research/17-seg-caller.md`) | `IBD2Seg` 822 → 896 of 982 exact; mean `PropIBD` error ÷3.7, worst row ÷24 | **+0** |
| its bridge and gate, re-bisected (§14 of the same doc) | none at all — a binary with the change reverted scores the same `.seg` scorecard to the digit; on constructed canvases it goes 5 723 → **6 000 of 6 000** | **+0** |
| the `IBD1Seg` overlap rule (`18-ibd1-caller.md`) | `IBD1Seg` 826 → **982 of 982** exact | **+5** |
| the IBD2 segment fringe (`19-ibd2seg-residual.md`) | `IBD2Seg` 896 → **982 of 982** exact | 408 → **436** |
| `.seg`'s two writer rules (`20-seg-writer.md`) | no estimate changed at all; byte-exact rows 806 → **982 of 982** | 436 → **464** |
| `<prefix>X.seg` implemented (`crate::analysis::xseg`) | a new file, 28 rows, byte-exact | 464 → **466** |
| the `--seglength` run merge (`20-seglength-floor.md`) | at the 10 Mb floor `IBD1Seg` 844 → **960 of 982** and byte-exact rows 832 → **943**; mean `PropIBD` error ÷3.2, worst row 0.0916 → 0.0111. Nothing at 3 Mb, where the rule cannot fire | 466 → **472** |
| the push and the IBD2 merge, re-measured (`21-push-merge.md`) | `--seglength 5` becomes as exact as the default floor: both estimate columns 982 of 982 and byte-exact rows 947 → **982**. At 10 Mb `IBD2Seg` 945 → **972**, byte-exact rows 943 → **970** | 472 → **475** |
| the gate window's own length bound and the IBD1 merge's budget word set (`23-gap-bound.md`) | `--seglength 10` becomes exact too: `IBD1Seg` 970 → **982**, `IBD2Seg` 972 → **982**, byte-exact rows 970 → **982**, printed-column MAE 0.000046 → **0** | 475 → **477** |
| `<prefix>build.log`'s header, `Duplicate … is removed.` and `RULE` lines (`crate::analysis::build`) | 0 → **243 of 806 bytes**, 6 of 18 lines, every one byte-identical; byte-identical on **53 of 59** held-out pedigree shapes | **+0** — a case is all-or-nothing |
| the merge queue, the clustering gate and the duplicate rule (`crate::analysis::unrelated`) | no corpus row at all: **19 of 19**, **19 of 19** and **27 of 27** held-out shapes, against 7, 18 and 21 for the rules they replaced | **+0** — invisible to the corpus |

The `20-seg-writer.md` row is the point about graders: `PropIBD` computed from the printed
columns instead of the totals, and rows listed in 16-sample blocks instead of by index.
Neither touches a segment, an estimate or a reported pair — and between them they were worth
28 cases, because the numbers underneath had finally stopped being wrong. The three
`--seglength` rows are the opposite case: real changes to the caller. The run merge was worth
6 cases and 158 exact rows (47 at the 5 Mb floor, 111 at 10); the push and IBD2-merge
corrections 3 cases and 62 more (35 at 5 Mb, 27 at 10), which took `--seglength 5` to exact;
and the window bound with the budget word set took `--seglength 10` there too.

The last three rows are the two ends of that spectrum: `23-gap-bound.md` moved both
scoreboards, while `build.log`'s `RULE` lines and the three clustering corrections moved
neither despite being right — the corpus simply has no fileset that can tell the difference,
which is why they were graded on 59, 19 and 27 held-out shapes instead.

A further track landed **nothing on purpose, and that is the result** —
`docs/research/22-screen.md` *proves* the screening count is not the kinship over any subset
of markers (the algebra is exact and covers every subset and every weighting), then closes
the two mechanisms that survived the proof, and declines to fit the one constant that would
have reproduced `bigish` and nothing else. `docs/PARITY.md` §5.7 carries the proof and a list
of what a future maintainer should *not* attempt.

So read both. `docs/PARITY.md` §4.4 is the row-level scoreboard, §3 the file-level one, §4.6
the out-of-sample one — the only one that still discriminates — and §5.0 says which grader to
use for what, what is left, and which experiment to run next.

## Scope (v1)

`byte-identical` below means exactly that — every file, every column, plus stdout, stderr
and exit status — on every dataset and flag combination the corpus captures. Counts are
cases, not files.

| Flag | Output files | Status |
| --- | --- | --- |
| `--kinship` | `.kin` (10 col), `.kin0` (8 col), `X.kin`, `X.kin0` | **byte-identical** (13/13 datasets, 220/220 param cases) |
| `--duplicate` | `.con` | **byte-identical** (13/13) |
| `--bysample` | `bySample.txt`, `allsegs.txt` | **byte-identical** (13/13) |
| `--bySNP` | `bySNP.txt`, `allsegs.txt` | **byte-identical** (13/13) |
| `--autoQC` | `_autoQC_Summary.txt`, `_autoQC_snptoberemoved.txt`, `_autoQC_sampletoberemoved.txt`, `_autoQC_updatesex.txt` | **byte-identical** (13/13) |
| `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt`, `allsegs.txt` | **byte-identical** (26/26) |
| `--ibs` | `.ibs`, `.ibs0`, `allsegs.txt` | **byte-identical** (13/13) — every column, `MaxIBD2` and `Pr_IBD2` included, on all 21 561 rows |
| `--cluster` | `allsegs.txt`, `updateids.txt`, `cluster.kin` | **byte-identical** (13/13) |
| `--build` | `updateids.txt`, `updateparents.txt`, `build.log`, `allsegs.txt` | 12/13 — `updateids.txt` and `updateparents.txt` are byte-identical in all 8 cases that write them; on `bigish` only `build.log` differs, and only in its `INFERENCE` half — its header and `RULE` lines are byte-identical (6 of 18 lines, 243 of 806 bytes, a strict subsequence of the reference's file), while the rest needs exact segment *placement* and one still-unidentified ordering (`docs/PARITY.md` §6.2) |
| `--related` | `.kin` (16 col), `.kin0` (14 col), `X.kin`, `allsegs.txt` | 64/65 — every file byte-identical in every case, all 4 805 16-column rows exact on all sixteen. The one failure is a single stdout line, the two-stage screening count on `bigish` (`docs/PARITY.md` §5.7) |
| `--ibdseg` | `.seg`, `allsegs.txt`, `splitped.txt`, `X.seg` | 64/65 — **52/52 alone**, 12/13 with `--related`. `.seg`, `allsegs.txt`, `splitped.txt` and `X.seg` are byte-identical in **every** case, on all 4 172 `.seg` rows, at 3, 5 and 10 Mb alike; the one loss is `bigish --related --degree 2 --ibdseg` on the stdout line above, with its `.seg` byte-identical. Off-corpus: `docs/PARITY.md` §4.6, and `splitped.txt` is written unconditionally where the reference does not always do so (§5.10) |

`--related` is **not** a synonym for `--kinship`: it emits six extra columns
(`HetConc`, `HomIBS0`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`), four of which come from
the IBD-segment engine, so full `--related` parity depends on `--ibdseg`. Below 10 samples
the reference itself downgrades `--related` to the `--kinship` path and emits the
10-column form; `--ibdseg` does the same below 5 samples.

The X chromosome **is** in scope, and each of its three passes has its own gate.
`--kinship` writes `<prefix>X.kin` and `<prefix>X.kin0` — with their own three sex-specific
estimators — when the map holds 512 or more X markers, no `--degree` is given and there is
more than one family; byte-identical in all 17 and all 5 diffable cases. `--related`'s
`<prefix>X.kin` and `--ibdseg`'s `<prefix>X.seg` are gated instead on the X map yielding a
**usable segment** (there is no marker-count threshold on those two), and `X.seg`
additionally on `--degree` being non-zero. `X.seg` carries the reference's own malformed
header — eleven names over nine-value rows — deliberately.

Note that `--prefix` is a plain **concatenation**, not a stem plus separator:
`--prefix ZZ_` yields `ZZ_.kin` and `ZZ_allsegs.txt`. The reference also opens
`<prefix>$TMP$.ped` for writing while it loads the `.fam`, so an unwritable prefix is a
fatal error there rather than at output time.

Out of scope for v1: `--pca`, `--mds`, `--roh`, `--lmm`, `--tdt`, `--gdt`, `--risk`,
`--makeGRM`, `--plink`, the R plotting flags (`--rplot`, `--pngplot`, `--rpath`), and
multi-dataset input. These are still *accepted* on the command line so the banner stays
byte-exact, then rejected at dispatch rather than silently ignored.

`.segments.gz` is **not** a target: the reference 2.3.2 build ships without zlib in its
segment writer, so it never produces that file despite the manual documenting it.

## Three things the reference itself gets wrong

All three are measured. `docs/PARITY.md` §4.3, §5.1 and §5.2 have the evidence — and the
third one stopped being an obstacle only when it was taken literally rather than treated as
noise.

* **`<prefix>X.kin0` is written by racing threads.** Six identical reference runs produced
  six different files, one truncated to 187 of 662 bytes, with records torn mid-number.
  `--cpus 1` makes it deterministic. The harness diffs only the `--cpus 1` captures — 5 of
  the 13 — and all 5 are byte-identical.
* **The four small datasets do not grade the segment caller.** They each report the 14
  within-family pairs of one six-person nuclear family over 5 000–10 000 markers, against
  `bigish`'s 50 000. On `monomorphic` the reference reports two full siblings as
  `IBD1Seg 0.9800 / IBD2Seg 0.0000`, labels them `PO`, and its own `--kinship` puts the same
  pair at 0.3384 where those segment numbers imply 0.2450. open-king reproduces that row
  exactly, which says nothing about either implementation recovering the underlying IBD.
* **The reference disagrees with itself about `PropIBD`.** In a single
  `--related --degree 2 --ibdseg` run on `bigish`, 147 pairs appear in both `king.kin` and
  `king.seg` — all 147 with identical `IBD1Seg` and `IBD2Seg`, and **43** with a different
  `PropIBD` in the two files (e.g. 0.5048 against 0.5049). Corpus-wide the two writers
  disagree on 54 of the 201 pairs the reference puts in both. This turned out to be the last
  thing standing between the project and an exact `.seg`: `.kin` computes `PropIBD` from the
  full-precision totals and `.seg` computes it from the four-decimal columns it is about to
  print (`i2*1e-4 + i1*5e-5`, exact on all 4 172 captured `.seg` rows including all 1 313
  that land on an exact decimal half). open-king implements both, one per writer, and
  reproduces each file. `docs/research/20-seg-writer.md`.

## Building

```bash
cargo build --release
```

The binary is emitted as `target/release/king`. It builds from a clean checkout in
**about 9.5 seconds** with no external toolchain: `Cargo.lock` has 15 packages, three of which
are this workspace. `cargo test --workspace` is 325 tests; CI additionally replays all 480
captured invocations against `tests/parity/BASELINE.txt` on every push, and fails on any
difference **in either direction** — an unrecorded improvement is a failure too, so the
committed baseline can never drift from what the tree actually does. Contributing,
regenerating the corpus, re-capturing goldens and the fixture technique the segment work
depends on: `docs/MAINTAINING.md`.

## Relationship to the original KING

This project is **not** affiliated with, endorsed by, or derived from the source code
of the original KING program by Wei-Min Chen. It is an independent implementation
written from:

* the published algorithm descriptions in the peer-reviewed literature,
* the publicly documented input/output file formats, and
* black-box observation of the reference binary's output.

No KING source code was read or copied. The original KING remains the work of its
authors under its own license terms; if you use relatedness inference in published
research, cite the original paper.

`docs/MAINTAINING.md` §1 states the clean-room rule, why it is absolute rather than
best-effort, and exactly which three sources of information are permitted. Every rule in the
segment engine names the black-box experiment that established it.

## Citation

**Cite the method, not this reimplementation.** The estimators are Manichaikul *et al.*:

> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. Robust relationship
> inference in genome-wide association studies. *Bioinformatics*. 2010;26(22):2867–2873.
> doi:10.1093/bioinformatics/btq559

**Credit for KING itself belongs to Wei-Min Chen and colleagues** (University of Virginia),
whose program — <https://www.kingrelatedness.com/> — this project reimplements and measures
itself against. The IBD-segment algorithm in particular is theirs and is unpublished;
open-king recovered its behaviour by observation, not by reading their work.

If a paper needs to say which implementation produced a number, `CITATION.cff` in this
repository carries machine-readable metadata for open-king alongside both references above.

## License

MIT — see [LICENSE](LICENSE). This license covers open-king's own code only; it makes no
claim about the original KING, which is separately licensed by its authors.
