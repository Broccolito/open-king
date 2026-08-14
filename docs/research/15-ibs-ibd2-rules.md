# The `--ibs` IBD2 columns, measured on constructed filesets

`MaxIBD2` and `Pr_IBD2` are the only place the reference prints an **individual** IBD2
segment length, so they grade an IBD2 rule one segment at a time instead of through an
aggregate. Everything below was established by building filesets in which the pair's
per-word genotype pattern is exact by construction, running the reference binary, and
reading the two columns back — never by reading KING's source.

Reproduce with `python3 docs/research/fixtures/ibs_rules.py` (needs the reference binary
at the path `fixtures/fixlab.py` names). The rig:

* **chr1** — the *carrier*: 30–60 complete words in which the pair is IBD1, which keeps
  the pair's kinship above the `--ibs` IBD2 gate (`2^-3.5`) and the usable total above the
  100 Mb floor that makes the two columns appear at all.
* **chr2** — the *canvas*: complete words whose content is written marker by marker.
  A "wall" word is 64 het-vs-hom mismatches — dirty under any rule, and worth nothing to
  any count. Blocks between walls are painted from `HetHet`, `A2A2/A2A2`, `A1A1/A1A1`,
  opposite-homozygote and missing markers in stated proportions.
* Marker spacing is uniform where only a length is needed, and **per-word** (word *w* at
  `40 000 + 137 w` bp) where the reported length has to identify *which* words were
  called: with distinct word widths the length inverts to exactly one interval `[u, e]`.

The pair is `S00`/`S01` in separate families, so its row is in `<prefix>.ibs0`.

---

## 1. The ruler, and the interval that is measured

A block of complete words `u..v` painted IBD2 inside a canvas the pair is unrelated over
reports

```
MaxIBD2 = pos[64·e + 63] − pos[64·u]        e = w1 if v + 2 ≥ w1 else v + 1
```

where `w0..w1` are the usable segment's complete words. Sliding a three-word block along a
ten-word segment and reading the length back gives spans of **3, 4, 5, 4, 4, 4** words as
the number of trailing dirty words goes 0, 1, 2, 3, 4, 5 — the call swallows *all* the
trailing words when there are at most two of them, and exactly one otherwise. The
per-word-spacing rig confirms the interval directly rather than through its length.

There is **no marker-level refinement** at either end: putting a single mismatch (or a
single opposite homozygote) at offset 0, 1, 5, 20 or 63 of the run's first or last word
does not move the reported length by one marker interval.

## 2. Which words break a run

Two clean blocks separated by a gap of *g* words, on an IBD1 canvas so that no word except
the gap carries anything unusual:

| gap content | `g = 1` | `g = 2` |
| --- | --- | --- |
| 64 het-vs-hom mismatches | merged | split |
| 5 mismatches, 59 HetHet | merged | split |
| 4 mismatches, 60 HetHet | merged | merged |
| **64 opposite homozygotes** | merged | **merged** |
| 64 missing calls | merged | merged |
| 64 `A1A1/A1A1` | merged | merged |
| 4 mismatches + 60 opposite homozygotes | merged | merged |
| 5 mismatches + 59 opposite homozygotes | merged | split |

So a word breaks an `--ibs` IBD2 run **iff it carries five or more het-vs-hom
mismatches**. Opposite homozygotes are irrelevant to this pass at any density, and so are
missing calls. One dirty word between two clean ones is absorbed; two in a row are not.

## 3. `--ibs` and `--ibdseg` do not report the same calls

Take a `W`-word all-HetHet block bounded by words in which *every* marker is an opposite
homozygote, on a chromosome that is otherwise all such words, and run both analyses on the
same fileset:

| W | `--ibs` `MaxIBD2` | interval it inverts to | `--ibdseg` `IBD2Seg` |
| ---: | ---: | --- | ---: |
| 4 | 25 953 327 | the whole usable segment `[0, 9]` | 0.0531 |
| 6 | 31 257 181 | the whole usable segment `[0, 11]` | 0.0779 |

`--ibs` runs straight through the IBS0 words; `IBD2Seg` is worth about the block alone.
The two passes agree about the *ruler* — `dups`' MZ pair prints `IBD2Seg 1.0000` and
`Pr_IBD2 0.8984`, and 0.8984 is the word-aligned total over the same `D` — but they do not
agree about the calls, so `king_core::ibdseg::Scan::ibd2_words` (this pass) and
`Scan::ibd2` (`.seg`) are separate functions over the same word masks until someone
explains the IBS0 asymmetry.

## 4. `Pr_IBD2` — the 10 Mb rule gates the pair, not the call

Sweeping the marker spacing under a fixed three-word call moves its length across 10 Mb:

| call length | `MaxIBD2` | `Pr_IBD2` |
| ---: | --- | --- |
| 9 932 000 | printed | `0.0000` |
| 9 999 423 | printed | `0.0000` |
| 10 000 187 | printed | `0.0498` |
| 10 027 500 | printed | `0.0498` |

but a fixture with **two** calls — 9.55 Mb and 19.15 Mb — reports a numerator of
28 699 910, which is both of them. So the threshold is not a filter on each call: a pair
whose longest call is under 10 Mb prints `Pr_IBD2 0.0000`, and a pair with one call at or
over it counts **every** call, short ones included. `MaxIBD2` is never gated.

The lengths summed are the word-aligned ones. A usable segment's trailing fringe (the
markers past its last complete word) does not enter the sum and does not help a short call
across the 10 Mb line: with a 40-marker fringe the same 9.55 Mb call still prints
`Pr_IBD2 0.0000`.

## 5. The acceptance count: 95 HetHet markers over the measured interval

A block of complete words holding exactly `k` HetHet markers, every other marker
`A2A2/A2A2` in the whole cohort, walled on both sides:

* **`k ≥ 95` is reported, `k ≤ 94` is not.** Eight bisections — block widths 2, 3, 4, 5
  words, block positions, 50 kb and 150 kb spacing, 6 and 12 samples, carrier lengths 10
  and 20 words — all land on exactly 95.
* **Only HetHet counts.** The same fixture with 200 `A1A1/A1A1` markers (which are
  `inf2`-informative under the `.seg` pass's rule) and no HetHet reports nothing.
* **The window is the measured interval `[u, e]`, not the run.** Loading the terminating
  word with 59 HetHet markers drops the block's own requirement from 95 to 36, and with 54
  it drops to 41; loading the word *before* the run leaves it at 95.
* **The count is waived when the run reaches the segment's last two words.** A block of
  384 markers with no HetHet at all is reported when it ends on `w1` or `w1 − 1` and
  refused one word earlier; sliding it along a 10-word and a 16-word canvas shows the
  exemption follows `w1`, and it is the same whether the canvas is the first chromosome or
  the last, so it is the *usable segment's* end and not the genotype array's.
* **The interval must span three words.** A two-word interval — only reachable against
  `w1` — is refused with 128 HetHet markers in it, while the same block one word earlier
  (three measured words) is reported.

## 6. What was still open: sustained low-grade mismatch — **SOLVED, see `16-…`**

> **Resolved by `docs/research/16-segment-extension.md`, which is the document to read
> instead of this section.** The boundary tabulated below is not an acceptance boundary at
> all: it is where the reported interval gets **cut**. The rule is a quantised confirmation
> scan — a run is confirmed in chunks of five het-vs-hom mismatches, each needing ≥ 95
> HetHet over ≥ 3 words — so `19`, `32` and `63` below are `5×19 = 95`, `3×32 = 96` and
> `2×63 = 126`, counted *within a chunk*, from wherever that chunk starts. That is also why
> the order matters. The rule is committed in `Scan::ibd2_words` and both `--ibs` IBD2
> columns are now exact on the whole corpus; the last paragraph of this section (what the
> residual costs) no longer describes the engine.


The rules above are exact for words that are either clean or hard-dirty. They are **not**
complete for words carrying one to four mismatches *and* few HetHet markers, which is what
real data looks like. Uniform blocks of eight words, `h` HetHet and `m` mismatches per
word, are reported only above a boundary that no rule tried here reproduces:

| mismatches per word | smallest HetHet per word reported |
| ---: | ---: |
| 0 | 12 (this is just the 95 total) |
| 1 | 19 |
| 2 | 32 |
| 3 | 63 |
| 4 | not reachable (≤ 60 tested) |

The boundary is not linear in `(h, m)`: `m = 1, h = 16` is refused while `m = 2, h = 32`
is reported, and the two have the same ratio. It is also **order-dependent**, which rules
out any per-word predicate: eight clean words followed by eight words at `m = 4` are
reported as one interval covering all sixteen, while the same sixteen words in the
opposite order report only the clean eight — a run that has started absorbs words that
could not have started it. Placement inside a word does not matter (three mismatches
adjacent and three spread 21 markers apart behave identically), so whatever the rule is,
it is a function of per-word counts and their order.

This is what the residual costs: on the corpus the implemented rule prints `Pr_IBD2`
slightly high on 59 of 673 `.ibs` rows (mean 0.024) and calls three `monomorphic` pairs
the reference refuses outright — all of them runs of 10–17 words carrying one mismatch per
word at 10–20 HetHet per word, which is exactly the region above.
