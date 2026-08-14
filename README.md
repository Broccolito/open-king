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

**365 of the 480 captured reference invocations reproduce byte-identically**, including
all 220 flag-plumbing and error probes. Run the suite yourself:

```bash
python3 tests/parity/run_parity.py --impl target/release/king
```

`docs/PARITY.md` is the authoritative matrix — every case group with its status, and §11
measures each remaining gap. The one-paragraph version: the relatedness estimators,
the QC reports and the whole command-line surface are done; the IBD-segment *caller* is
not, and it is what `--related`, `--ibdseg` and two columns of `--ibs` are waiting on.

## Scope (v1)

`byte-identical` below means exactly that — every file, every column, plus stdout, stderr
and exit status — on every dataset and flag combination the corpus captures.

| Flag | Output files | Status |
| --- | --- | --- |
| `--kinship` | `.kin` (10 col), `.kin0` (8 col), `X.kin`, `X.kin0` | **byte-identical** (13/13 datasets, 220/220 param cases) |
| `--duplicate` | `.con` | **byte-identical** (13/13) |
| `--bysample` | `bySample.txt`, `allsegs.txt` | **byte-identical** (13/13) |
| `--bySNP` | `bySNP.txt`, `allsegs.txt` | **byte-identical** (13/13) |
| `--cluster` | `allsegs.txt`, `updateids.txt`, `cluster.kin` | byte-identical 12/13 — merged-cluster tail missing |
| `--build` | `updateids.txt`, `updateparents.txt`, `build.log`, `allsegs.txt` | byte-identical 11/13 |
| `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt`, `allsegs.txt` | byte-identical 12/13 — 3 of 84 rows differ on the one merged dataset |
| `--ibs` | `.ibs`, `.ibs0`, `allsegs.txt` | every column byte-identical **except** `MaxIBD2`/`Pr_IBD2` (57 and 137 of 21 561 rows) |
| `--ibdseg` | `.seg`, `allsegs.txt`, `splitped.txt` | `allsegs.txt` byte-identical; `.seg` **not** (558 of 982 rows exact, 188 extra); `splitped.txt` unwritten |
| `--related` | `.kin` (16 col), `.kin0` (14 col), `allsegs.txt` | byte-identical only on the five datasets the reference downgrades to `--kinship`; the 16-column pass is unimplemented |
| `--autoQC` | `_autoQC_Summary.txt`, `_autoQC_snptoberemoved.txt`, `_autoQC_sampletoberemoved.txt` | unimplemented |

`--related` is **not** a synonym for `--kinship`: it emits six extra columns
(`HetConc`, `HomIBS0`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`) that come from the
IBD-segment engine, so full `--related` parity depends on `--ibdseg`. Below 10 samples
the reference itself downgrades `--related` to the `--kinship` path and emits the
10-column form; `--ibdseg` does the same below 5 samples.

X-chromosome kinship **is** implemented: with 512 or more X markers, no `--degree`, and
more than one family, `--kinship` writes `<prefix>X.kin` and `<prefix>X.kin0` with their
own three sex-specific estimators. `<prefix>X.seg` (the `--ibdseg` X table) is not.

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

## Building

```bash
cargo build --release
```

The binary is emitted as `target/release/king`.

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
