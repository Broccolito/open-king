# open-king v0.1.0

First release of **open-king**, a clean-room, MIT-licensed reimplementation of
**KING 2.3.2** for relatedness inference from PLINK filesets.

The goal is drop-in replacement: the same command line, the same input files, and
byte-identical output.

```bash
king -b study.bed --related --prefix study
```

## Parity

**477 of the 480 captured reference invocations reproduce byte-identically.** "Byte
-identical" means the whole invocation — every output file, every column, plus stdout,
stderr and exit status.

The differential harness is in the repository, so you can check this yourself rather than
taking our word for it:

```bash
cargo build --release
python3 tests/parity/run_parity.py --impl ./target/release/king
```

Running the harness against the *reference* binary scores 480/480, which is what
establishes that a failure means a real difference rather than harness noise.

### Byte-identical

| Analysis | Outputs |
| --- | --- |
| `--kinship` | `.kin`, `.kin0`, `X.kin`, `X.kin0` |
| `--related` | `.kin` (16 col), `.kin0` (14 col), `allsegs.txt` |
| `--ibdseg` | `.seg`, `X.seg`, `allsegs.txt`, `splitped.txt` |
| `--duplicate` | `.con` |
| `--ibs` | `.ibs`, `.ibs0` |
| `--unrelated` | `unrelated.txt`, `unrelated_toberemoved.txt` |
| `--cluster` | `cluster.kin`, `updateids.txt` |
| `--bysample` / `--bySNP` | `bySample.txt`, `bySNP.txt` |
| `--autoQC` | the four `_autoQC_*` files |

Segment estimates are exact on all 982 corpus rows at every supported `--seglength`
floor (3, 5 and 10 Mb). The command-line surface was additionally validated by a
differential fuzz of 10,000 random command lines against the reference.

### Known gaps

Three cases, all on the largest test dataset, and none of them a wrong estimate:

1. **The two-stage screen's stdout line.** `--related --degree 2` reports a different
   count of screened pairs. **No output file is affected** — `.kin0`'s rows come from the
   exhaustive re-estimate that follows the screen and are byte-correct at every degree,
   and every pair the reference's screen drops falls below the reporting threshold
   anyway. Four rounds of measurement have ruled out the obvious mechanisms; see
   `docs/research/22-screen.md`.
2. **`build.log`.** Every line emitted is byte-identical, but the file is a subsequence:
   some lines the reference prints are still missing.
3. **`--related --degree 2` stdout** on one dataset, for the same reason as (1).

`docs/PARITY.md` is the authoritative matrix and measures every gap.

## Documentation

| Document | For |
| --- | --- |
| `docs/CLI.md` | every command-line option |
| `docs/OUTPUTS.md` | every output file and column |
| `docs/COOKBOOK.md` | task-oriented recipes |
| `docs/INTERPRETING.md` | how to read the numbers, and the pitfalls |
| `docs/PARITY.md` | the parity matrix and every measured gap |
| `docs/MAINTAINING.md` | how to continue the work |

## Install

Download the archive for your platform below, or build from source:

```bash
cargo build --release      # -> target/release/king
```

Binaries are provided for macOS (arm64 and x86_64), Linux x86_64, and Windows x86_64.
Verify with the published `SHA256SUMS.txt`.

## Provenance

This project is **not** affiliated with, endorsed by, or derived from the source code of
the original KING by Wei-Min Chen. It was written from the published algorithm
descriptions, the publicly documented file formats, and black-box observation of the
reference binary's behaviour. No KING source code was read or copied.

If you use relatedness inference in published research, cite the original work:

> Manichaikul A, Mychaleckyj JC, Rich SS, Daly K, Sale M, Chen WM. Robust relationship
> inference in genome-wide association studies. *Bioinformatics*. 2010;26(22):2867-2873.

## Caveat

Parity is measured against **one build** — KING 2.3.2, macOS arm64. KING's segment
numerics changed across earlier 2.1.x and 2.2.x releases, so agreement with a different
build is not implied.

## License

MIT.
