# open-king

A clean-room, MIT-licensed reimplementation of **KING** (Kinship-based INference for
Genome-wide association studies) — the relatedness-inference program described in
[Manichaikul *et al.* 2010, *Bioinformatics*](https://pmc.ncbi.nlm.nih.gov/articles/PMC3025716/).

The goal is **drop-in parity**: the same command line, the same input files, and
byte-identical output files as the reference `king` 2.3.2 binary.

```bash
king -b study.bed --related --prefix study
```

## Status

**408 of the 480 captured reference invocations reproduce byte-identically (85.0 %)**,
including all 220 flag-plumbing and error probes. Run the suite yourself:

```bash
cargo build --release
python3 tests/parity/run_parity.py --impl target/release/king
```

`docs/PARITY.md` is the authoritative claim — the full analysis × dataset matrix, the
measured size of every remaining gap, and a labelled limitations section. Everything below
is a summary of it and says nothing it does not support.

The one-paragraph version: the relatedness estimators, the QC reports, duplicate detection,
auto-QC, unrelated-set selection, `--ibs` and the whole command-line surface are
byte-identical everywhere. What is left is the `.seg` IBD-segment **caller**, which finds
exactly the right set of segments — **0 spurious and 0 missing rows on every output file in
the corpus** — but places a called IBD2 segment's endpoints a few markers from the reference.
So the segment columns are close without being equal: on the primary `--ibdseg` capture
**896 of 982** rows carry both estimates exactly, every one of the 982 `InfType` labels is
right, and the mean absolute `PropIBD` error is **0.00006** with a worst row of **0.0042**.
**`IBD1Seg` is exact on all 982 rows**, and on every `.kin`, `.kin0`, `X.kin` and
`cluster.kin` row in the corpus; whether a pair has any IBD2 at all is decided correctly
everywhere. The residual is purely the *length* of the 86 rows whose IBD2 calls are a few
markers out.

Of the 72 cases that are not byte-identical, 71 are that one cause; the remaining one is
`--build`'s pedigree reconstruction, which is unimplemented — and which is itself blocked on
the same segment residual, so it is not a second independent problem.

### Two numbers, and why they are different

The headline is a **whole-file** count: a case turns `PASS` only when every row of every file
it writes is byte-exact. The residual is spread thinly across nearly every dataset, so it is
routine for a large row-level gain to move the headline by little or nothing. That is not a
caveat invented after the fact — it is how the last two corrections landed:

| change | row-level effect | headline |
| --- | --- | --- |
| the `.seg` IBD2 caller (`docs/research/17-seg-caller.md`) | `IBD2Seg` 822 → 896 of 982 exact; mean `PropIBD` error ÷3.7, worst row ÷24 | **+0** |
| its bridge and gate, re-bisected (§14 of the same doc) | none at all — a binary with the change reverted scores the same 408/480 and the same `.seg` scorecard to the digit; on constructed canvases it goes 5 723 → **6 000 of 6 000** | **+0** |
| the `IBD1Seg` overlap rule (`docs/research/18-ibd1-caller.md`) | `IBD1Seg` 826 → **982 of 982** exact; mean error 0.00037 → 0.00007 | **+5** |

So read both. `docs/PARITY.md` §4.4 is the row-level scoreboard, §3 the file-level one, and
§5.0 says which grader to use for what.

`--ibs` left the residual entirely when its IBD2 caller was solved outright: both of its IBD2
columns, `MaxIBD2` and `Pr_IBD2`, are exact on all 21 561 rows of every dataset. It is a
*different* caller from `.seg`'s (`docs/research/16-segment-extension.md`, `docs/PARITY.md`
§5.8), so the `.seg` caller needed its own campaign and its own instrument — a canvas that
inverts a printed segment column back to the number of calls and the number of words
(`docs/research/17-seg-caller.md`, and `18-ibd1-caller.md` for the same canvas built the other
way up, so the painted region is IBD2-free). `docs/PARITY.md` §4.4 is where both are graded,
and §5.0 says what is left, in what shape, and which experiment to run next.

Two known differences sit outside the 480 captures entirely, so they cost no case but a user
could still hit them: `--ibdseg` does not apply the reference's 100 Mb usable-total floor,
and `splitped.txt` is written unconditionally. Both are measured, localised and written up in
`docs/PARITY.md` §5.10.

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
| `--cluster` | `allsegs.txt`, `updateids.txt`, `cluster.kin` | 12/13 — on `bigish`, `cluster.kin`'s segment columns differ on 16 of 165 rows; the pair set, the ordering and every other column are exact |
| `--build` | `updateids.txt`, `updateparents.txt`, `build.log`, `allsegs.txt` | 12/13 — on `bigish` the pedigree-reconstruction rules are unimplemented, so `build.log` and `updateparents.txt` come out empty where the reference writes 18 and 33 lines; `updateids.txt` is already byte-identical |
| `--related` | `.kin` (16 col), `.kin0` (14 col), `X.kin`, `allsegs.txt` | 40/65 — the only columns that differ anywhere are `IBD2Seg` (463 of 4 805 rows) and `PropIBD` (456). `IBD1Seg`, `InfType`, `Error`, `HetConc`, `HomIBS0` and the ten `--kinship` columns are exact on every row |
| `--ibdseg` | `.seg`, `allsegs.txt`, `splitped.txt`, `X.seg` | 20/65 (16/52 alone, 4/13 with `--related`) — `allsegs.txt` byte-identical everywhere; `.seg` has the right rows (0 extra, 0 missing) with 1 086 of 4 169 differing numerically; `X.seg` is not written; `splitped.txt` is byte-identical in all 50 corpus cases but is written unconditionally, which the reference does not always do (`docs/PARITY.md` §5.10) |

`--related` is **not** a synonym for `--kinship`: it emits six extra columns
(`HetConc`, `HomIBS0`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`), four of which come from
the IBD-segment engine, so full `--related` parity depends on `--ibdseg`. Below 10 samples
the reference itself downgrades `--related` to the `--kinship` path and emits the
10-column form; `--ibdseg` does the same below 5 samples.

X-chromosome kinship **is** implemented: with 512 or more X markers, no `--degree`, and
more than one family, `--kinship` writes `<prefix>X.kin` and `<prefix>X.kin0` with their
own three sex-specific estimators. `<prefix>X.seg` — which `--ibdseg` writes only when
`--degree` is non-zero and the fileset has usable X segments — is **not** implemented.

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

All three are measured, and they are why a handful of captures cannot be graded normally.
`docs/PARITY.md` §4.3, §5.1 and §5.2 have the evidence.

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
  `PropIBD` in the two files (e.g. 0.5048 against 0.5049). open-king computes it once, from
  the unrounded estimates, which matches `.kin` almost everywhere and costs 149 `.seg` rows a
  ±1 in the fourth decimal. Fifteen arithmetic reformulations were scored and all fifteen tie
  exactly; rounding off the printed columns instead is measurably worse, in both files and in
  both rounding modes. This is not a formula that can be fixed.

## Building

```bash
cargo build --release
```

The binary is emitted as `target/release/king`. It builds offline from a clean checkout in
under ten seconds with no external toolchain: `Cargo.lock` has 15 packages, three of which
are this workspace. Contributing, regenerating the corpus, re-capturing goldens and the
fixture technique the segment work depends on: `docs/MAINTAINING.md`.

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

## Citation

> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. Robust relationship
> inference in genome-wide association studies. *Bioinformatics*. 2010;26(22):2867-2873.

## License

MIT — see [LICENSE](LICENSE).
