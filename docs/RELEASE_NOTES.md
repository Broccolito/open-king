# open-king v0.1.0

`open-king` reads a PLINK1 fileset and reports how every pair of samples in it is related:
kinship coefficients, duplicate and monozygotic pairs, IBD segments, a relationship label
for each pair, and per-sample and per-marker quality control.

It is a clean-room, MIT-licensed Rust implementation of the relatedness core of
[KING 2.3.2](https://www.kingrelatedness.com/), taking the same command line and writing
byte-identical output. All 480 captured reference invocations reproduce exactly, across 876
output files.

## Install

Each archive holds a single file, the `open-king` executable. Unzip it and run it.

| Platform | Asset |
| --- | --- |
| macOS, Apple silicon | `open-king-macos-arm64.zip` |
| macOS, Intel | `open-king-macos-x86_64.zip` |
| Linux x86_64 | `open-king-linux-x86_64.zip` |
| Windows x86_64 | `open-king-windows-x86_64.zip` |

The macOS builds are signed with a Developer ID and notarized by Apple, so they run straight
out of the archive with no Gatekeeper prompt and no quarantine step.

```bash
unzip open-king-macos-arm64.zip
./open-king -b study.bed --related --degree 2 --prefix study
```

`SHA256SUMS.txt` carries a checksum for every asset.

## Documentation

<https://broccolito.github.io/open-king/>

Every command and option with a runnable example, the output-file reference, the parity
evidence, and the benchmarks.
