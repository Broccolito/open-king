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

Under active development. See [`docs/SPEC.md`](docs/SPEC.md) for the implementation
specification and [`docs/PARITY.md`](docs/PARITY.md) for the current parity matrix
against the reference binary.

## Scope (v1)

| Flag | Output files | Status |
| --- | --- | --- |
| `--kinship` | `.kin` (10 col), `.kin0` (8 col) | in progress |
| `--duplicate` | `.con` | in progress |
| `--ibs` | `.ibs`, `.ibs0`, `allsegs.txt` | in progress |
| `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt` | in progress |
| `--related` | `.kin` (16 col), `.kin0` (14 col), `allsegs.txt` | in progress |
| `--ibdseg` | `.seg`, `allsegs.txt`, `splitped.txt` | in progress |
| `--build` | `updateids.txt`, `updateparents.txt`, `build.log`, `splitped.txt` | in progress |
| `--bysample` | `bySample.txt` | in progress |
| `--bySNP` | `bySNP.txt` | in progress |

`--related` is **not** a synonym for `--kinship`: it emits six extra columns
(`HetConc`, `HomIBS0`, `IBD1Seg`, `IBD2Seg`, `PropIBD`, `InfType`) that come from the
IBD-segment engine, so full `--related` parity depends on `--ibdseg`. Below 10 samples
the reference itself downgrades `--related` to the `--kinship` path and emits the
10-column form.

Note that `--prefix` is a plain **concatenation**, not a stem plus separator:
`--prefix ZZ_` yields `ZZ_.kin` and `ZZ_allsegs.txt`.

Out of scope for v1: `--pca`, `--mds`, `--roh`, `--cluster`, `--lmm`, `--tdt`, `--gdt`,
`--risk`, `--makeGRM`, `--plink`, the R plotting flags (`--rplot`, `--pngplot`,
`--rpath`), X-chromosome analysis (`X.kin`, `X.kin0`, `X.seg`), and multi-dataset input.
These are still *accepted* on the command line so the banner stays byte-exact, then
rejected at dispatch rather than silently ignored.

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
