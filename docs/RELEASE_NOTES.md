# open-king v0.1.1

`open-king` reads a PLINK1 fileset and reports how every pair of samples in it is related:
kinship coefficients, duplicate and monozygotic pairs, IBD segments, a relationship label
for each pair, and per-sample and per-marker quality control.

It is a clean-room, MIT-licensed Rust implementation of the relatedness core of
[KING 2.3.2](https://www.kingrelatedness.com/), taking the same command line and writing
byte-identical output. All 480 captured reference invocations reproduce exactly, across 876
output files.

## What changed in 0.1.1

**`--ibdseg` now matches KING 2.3.2 on dense panels.** Two segment-calling rules were right
for small filesets and wrong for large ones. Both were found the same way — by running the
reference and this binary over a real 663,197-marker autosomal panel of 157 samples and
comparing the `.seg` files byte for byte — and both had to be fixed together.

* **The informativeness gate is a table, not a constant.** KING picks how many informative
  markers a candidate segment must carry once per run, from the fileset's total marker
  count: 10 below 400,000 markers, 20 from there to 2,000,000, and 100 above, with the
  minimum candidate length stepping alongside it. Every fixture and every parity dataset in
  this project sits in the first row, so the 10 measured there had been generalised into a
  universal rule. Both boundaries are bisected to the marker.
* **The IBD1 merge cap of two unusable words was a first-row reading too.** It was measured
  on synthetic filesets whose markers sit 20,000 bp apart; on a panel whose markers sit
  about 4,300 bp apart the reference joins runs across four. The cap now applies only to
  filesets under 400,000 markers.

Correcting the merge cap on its own makes that panel *worse* — an eighth related pair
appears that the reference does not report. Corrected together, `--ibdseg` is byte-identical
to KING 2.3.2 on it: all seven reported pairs, all four printed columns, plus
`allsegs.txt` and `splitped.txt`.

Nothing else moves. `--kinship` was already exact on the same panel and stays exact
(11,476 `.kin` rows, 770 `.kin0` rows). The 480-case captured corpus is unchanged at
480 PASS / 0 FAIL over 876 byte-compared files — every dataset in it is far below 400,000
markers, which is exactly why it could not see any of this.

Reported as [#13](https://github.com/Broccolito/open-king/issues/13) (the merge cap) and
[#14](https://github.com/Broccolito/open-king/issues/14) (the marker-count dependence);
the measurement is
[`docs/PARITY.md`](https://github.com/Broccolito/open-king/blob/main/docs/PARITY.md) §5.14.

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
