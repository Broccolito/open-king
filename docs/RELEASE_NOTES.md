# open-king v0.1.2

`open-king` reads a PLINK1 fileset and reports how every pair of samples in it is related:
kinship coefficients, duplicate and monozygotic pairs, IBD segments, a relationship label
for each pair, and per-sample and per-marker quality control.

It is a clean-room, MIT-licensed Rust implementation of the relatedness core of
[KING 2.3.2](https://www.kingrelatedness.com/), taking the same command line and writing
byte-identical output. All 480 captured reference invocations reproduce exactly, across 876
output files.

## What changed in 0.1.2

**The binary is unchanged.** No estimator, no output format and no command-line behaviour
moves in this release: 480 PASS / 0 FAIL over 876 byte-compared files, `baseline: MATCH`,
and the `.seg` row scorecard reads 982 / 982 / 982 with MAE 0.000000 at the 3, 5 and 10 Mb
floors — the same numbers 0.1.1 published. What changed is that the macOS archives are now
actually what this page has claimed they are.

**v0.1.1's macOS archives shipped unsigned, and have been re-signed in place.** CI holds no
Apple credential by design, so the macOS pair is signed and notarized by hand afterwards
(`docs/MAINTAINING.md` §9). That step ran for 0.1.1 and failed, silently enough that the
release went out anyway while `README.md`, this page and the documentation site all said the
builds were notarized.

The cause was the guard rather than the signing. `scripts/notarize-macos-release.sh`
asserted the hardened runtime with

```bash
codesign -dv --verbose=2 "$d/open-king" 2>&1 | grep -qE 'flags=0x10000\(runtime\)' || exit 1
```

`grep -q` exits at its first match — line 4 of the 14 `codesign` writes — so `codesign` is
still writing when the read end closes, takes `SIGPIPE`, and exits 141; `set -o pipefail`
three lines above then makes 141 the status of the pipeline. The signature was correct every
time and the script rejected it every time, reporting a signing failure for what was a
plumbing failure. The assertion now captures `codesign`'s output once and matches against
that, and checks the secure timestamp and `TeamIdentifier` while it has it, so a bad
signature is caught before submission rather than twenty minutes into Apple's queue.

**If you downloaded v0.1.1's macOS archives before 2026-08-22, replace them.** The assets
and `SHA256SUMS.txt` on the v0.1.1 release have been regenerated, so the checksums you
recorded then no longer match. The executable inside is byte-for-byte the one CI built from
the tagged commit — only the signature is added. Both submissions were `Accepted` by Apple,
`scripts/verify-release-assets.sh v0.1.1` passes all 17 checks, and
`codesign --test-requirement="=notarized"` is satisfied against the published archive.

**Two documentation-only repairs in `open-king-core`.** `seg_prop_ibd`'s doc comment linked
twice to `Segments::prop_ibd`, a type that does not exist — the method is on `PairSegments`
— which broke `cargo doc` under `-D warnings`. With that fixed the same build still failed
on 29 `private_intra_doc_links` errors, so that lint is now allowed crate-wide: linking a
rule's doc comment to the private constant that pins it is the convention this crate is
documented in, and the alternative was deleting 29 working cross-references.
`RUSTDOCFLAGS="-D warnings" cargo doc -p open-king-core --no-deps` exits 0.

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
