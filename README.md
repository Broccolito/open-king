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
| `--related` | `.kin`, `.kin0` | in progress |
| `--kinship` | `.kin`, `.kin0` | in progress |
| `--duplicate` | `.con` | in progress |
| `--ibs` | `.ibs`, `.ibs0` | in progress |
| `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt` | in progress |
| `--ibdseg` | `.seg`, `.segments.gz` | in progress |
| `--build` | rebuilt `.fam` / update files | in progress |
| `--bysample` | `bySample.txt` | in progress |
| `--bySNP` | `bySNP.txt` | in progress |

Out of scope for v1: `--pca`, `--mds`, `--roh`, `--cluster`, `--lmm`, `--tdt`,
`--gdt`, `--risk`, `--makeGRM`.

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
