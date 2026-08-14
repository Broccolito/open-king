# The IBD-segment calling rule, measured on controlled fixtures

**Status:** black-box measurement. Every statement below is a reading taken off the
reference binary `/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king`
(KING 2.3.2) run on filesets whose IBD state was **constructed**, not simulated through a
pedigree. KING's C++ source was not downloaded, opened or read. Where a rule is not
pinned, it says so.

Harness: `docs/research/fixtures/fixlab.py` (fileset writer) and
`docs/research/fixtures/rig2.py` (the two-chromosome rig). Both are standard-library
Python and re-run in seconds.

---

## 1. Method

### 1.1 Construction

Six samples, each in its own family, no parents. Samples `S00`/`S01` are the test pair;
`S02`–`S05` are unrelated padding, present only because `--ibdseg` silently downgrades to
`--kinship` below five samples.

Genotypes are built from explicit haplotypes, one marker at a time. Every sample draws two
haplotype alleles from a per-marker frequency; the pair's IBD state is imposed directly —
`h[S01][0] = h[S00][0]` over an IBD1 region, both haplotypes copied over an IBD2 region,
nothing copied elsewhere. IBD status is therefore exact by construction and local to the
marker, with no recombination model in the way. Genotype = count of A1 alleles; A1/A2 are
re-oriented per marker after drawing so A1 is the observed minor allele.

Three override layers sit on top, all applied before re-orientation:

| layer | effect |
| --- | --- |
| `force_ibs0` | make the pair opposite homozygotes at a marker (an IBS0), polarity alternating by marker index so allele frequencies stay balanced |
| `pat` | set the pair's two genotypes at a marker exactly |
| `pat_all` | set all six samples' genotypes at a marker exactly |

### 1.2 The rig

Two chromosomes, both a multiple of 64 markers so that chromosome 2 starts on a boundary
of the global 64-marker word grid and local word `w` is global word `n1/64 + w`.

* **chr1 — the carrier.** Entirely IBD1. It is always called whole (verified: `IBD1Seg`
  equals its length over the denominator to the last digit, over many seeds), so it
  contributes a known constant and guarantees the pair clears the pair-level filter.
* **chr2 — the canvas.** Background forced to IBS0 at *every* marker, with test blocks
  carved out of it. A solid background makes the flanking IBS0 positions deterministic,
  which is what turns the endpoint rules from noise into arithmetic.

Readback: the only per-pair output this build writes is `<prefix>.seg`, four decimals.
Segment length in base pairs is recovered as `IBD1Seg × denominator` and rounded to the
nearest marker interval, where the denominator is the sum of `Length` over
`<prefix>allsegs.txt`. Marker spacing is chosen so one ulp of the printed value is well
under one marker: at spacing 100 kb with a 153.4 Mb denominator, 1 ulp = 15.3 kb = 0.15
marker.

### 1.3 Two calibration facts found while building the rig

* **`--ibdseg` requires a total usable autosomal length of at least 100,000,000 bp.**
  Below it the run prints `Segments too short.` and writes no `.seg` at all. Bisected to
  the marker: 99.975 Mb → `Segments too short.`, 100.000 Mb → full output. This gate is
  the same 100 Mb that switches on `MaxIBD2`/`Pr_IBD2` in `.ibs`. Any fixture smaller than
  this measures nothing.
* **The `Too many first alleles as the major allele` fatal is unseeded** (as
  `docs/PARITY.md` §11.1 records). It fires on perhaps one run in three for skewed
  fixtures and produces a *missing output file*, which reads exactly like "the pair was
  filtered out". Two early conclusions in this investigation were wrong because of it.
  The harness now retries and raises rather than returning an empty result.

---

## 2. The IBD1 rule

### 2.1 The word test — no error tolerance at all

The scan is quantised to the 64-marker words of the **global** retained-autosome array. A
word is *clean* iff it contains no IBS0 (opposite homozygote) whatsoever.

Fixture: a 12-word IBD1 block on the solid background, one forced IBS0 swept across it,
spacing 100 kb (so that no piece can ever be lost to the 3 Mb reporting floor).

```
1 forced IBS0 in block word bw   ->  reported total (clean = 831 mk)
bw = 0                831 - 1 - 8(a-1)        left end trimmed
bw = 1 .. 10          830                     run splits, pieces abut
bw = 11 (last word)   704 + 8(a-1)            right end trimmed
```

Every interior position costs exactly one marker interval — the signature of a split whose
two pieces abut — and never zero. **One opposite homozygote disqualifies its entire word.**
There is no tolerance to find: sweeping the count `a` of IBS0 within one word (1, 2, 3, 4,
8, 16, 32, 64) and their spread (consecutive vs every 8th bit) changes only *where* the
endpoints land, never *whether* the word breaks.

> An earlier reading of this same sweep at 25 kb spacing appeared to show a tolerance of
> two IBS0 per word. It was the `--seglength` floor: at 25 kb a word is 1.6 Mb, so the
> short side of a split fell under 3 Mb and vanished, and the arithmetic coincided. At
> 100 kb the same sweep is flat. Do not re-derive the tolerance; it does not exist.

### 2.2 The endpoint rule — asymmetric, and exact

For a run of clean words `w0..w1` inside a usable segment:

```
lo = 1 + (index of the LAST IBS0 in word w0-1)      or  segment start if w0 is the first word
hi =     (index of the LAST IBS0 in word w1+1)      or  segment end   if w1 is the last word
```

Consequence: a run of `W` clean words bounded by IBS0-bearing words reports exactly
`64(W+1) - 1` marker intervals — it grows one whole word to the right, and not at all to
the left. Verified for `W` = 2,3,4,5,6,8 at 25 kb and 100 kb spacing, and independently by
sliding the block's start marker by marker:

```
block = [a, 512), background solid  ->  reported span = [a, 575] for every a
a = 255 -> 320 mk   a = 320 -> 255 mk   a = 383 -> 192 mk   a = 448 -> 127 mk
```

The predicted span was correct at all 26 values of `a` tested, including every crossing of
a word boundary.

### 2.3 There is no merge rule

Two IBD1 blocks separated by `g` markers of IBS0, `g` = 0,1,2,4,8,16,32,63,64,65,96,128,
129,192,256,320: the reported total equals *sum of the two pieces computed independently
by §2.2* at every single `g`, and never the merged span. `g = 0` is the only case that
reads as one segment, because there is no gap. Nothing is bridged, at any distance.

### 2.4 Minimum run length

Two or more clean words are **always** accepted: a 2-word block placed at each of 13
positions × 4 seeds gave 191 markers in all 52 runs.

Exactly one clean word is **data-dependent**. Sweeping a 1-word block across 14 positions
gives an irregular kept/dropped pattern that is reproducible, is a function of the global
word index for a given dataset, and is *not* explained by physical length — a 1-word run
spanning 12.7 Mb is dropped while another spanning the same 12.7 Mb is kept. It is also
not the `--seglength` floor (the effect survives `--seglength 1`) and not a
`het`/`het-het`/`hom-concordant` count of the word (words with identical count triples land
on both sides). **Open.** Treating the floor as "≥ 2 clean words" is right in every case
tested except this one, which is why `MIN_RUN1 = 2` fits the corpus as well as it does.

---

## 3. The IBD2 rule

* **IBD2 segments are not extended.** A run of `W` clean IBD2 words reports exactly
  `64W - 1` marker intervals — clipped to its own words, with none of the one-word
  right-hand growth that IBD1 gets. Measured at `W` = 1,2,3,4,6,8: 63, 127, 191, 255, 383,
  511.
* **One word is enough for IBD2**, unlike IBD1: a single IBD2 word reports its 63
  intervals.
* **Heterozygote disagreements confined to one word never break IBD2.** Forcing
  `k` het-vs-hom disagreements into one interior word of an 8-word IBD2 block leaves the
  call untouched for every `k` from 1 to 64 — including a word in which *all* 64 markers
  disagree. This is consistent with, and independently confirms, the two-word contingency
  table already in `king_core::ibdseg::Scan::het_break` (`a ≥ 2` on the left of a boundary
  *and* `b ≥ 1` on the right): with `b = 0` there is no break however large `a` is.
* A forced IBS0 inside an IBD2 block splits the IBD2 run on word boundaries and leaves an
  IBD1 remainder, as expected. The remainder's exact size is off by ±1 marker from the
  §2.2 model in some configurations; see §6.

---

## 4. `--seglength`, and the pair filter

Both thresholds are **inclusive**, which the console text ("Short IBD segments (<3Mb) are
not reported", "pairs without any long IBD segments (>10Mb) are excluded") does not say.

Measured at 10 kb spacing, where a span of exactly 300 intervals is exactly 3,000,000 bp:

| segment span | `--seglength 3` |
| ---: | --- |
| 2,990,000 bp | dropped |
| **3,000,000 bp** | **kept** |
| 3,010,000 bp | kept |

and with `--seglength 10`: 9,990,000 dropped, 10,000,000 kept. So the test is
`length_bp >= seglength_bp`.

The pair-level filter behaves the same way. With the carrier removed so the pair's only
segment is the one under test, and `--seglength 1` so the reporting floor cannot interfere:

| pair's longest segment | row in `.seg` |
| ---: | --- |
| 9,990,000 bp | absent |
| **10,000,000 bp** | **present** |
| 10,010,000 bp | present |

So the pair filter is `longest_reported_segment_bp >= 10,000,000`, not `>`.

**The length compared is the reported (extended) length of §2.2, not the run's own span.**
A 2-clean-word run at 25 kb spacing has an own-span of 3.175 Mb and a reported span of
4.775 Mb; it survives `--seglength 4` and dies at `--seglength 5`. The same crossing was
checked for runs of 3, 4, 5 and 6 words: the threshold always tracks `64(W+1)-1` intervals,
never `64W-1`.

---

## 5. The rule that is actually missing from our implementation

`docs/PARITY.md` §11.1 describes the open problem as "the reference accepts far fewer runs
of IBS0-free words than any run-length or physical-length filter explains". It is not a
run-length rule. **The reference additionally requires the genotypes under a run to be
informative, and informativeness is a property of the whole sample, not of the pair.**

### 5.1 The over-call, reproduced in isolation

chr2 constructed IBD0 by construction — the pair shares nothing — with the per-marker
allele frequency swept. chr1 is the IBD1 carrier so the pair is always reported.

| maf | IBS0-free words on chr2 | our model calls | the reference calls |
| ---: | ---: | ---: | ---: |
| 0.50 | 0 / 16 | 0 | **0** |
| 0.30 | 0 / 16 | 0 | **0** |
| 0.20 | 1 / 16 | 0 | **0** |
| 0.15 | 2 / 16 | 17.4 Mb in 1 segment | **0** |
| 0.10 | 7 / 16 | 37.8 Mb in 2 segments | **0** |
| 0.05 | 12 / 16 | **101.9 Mb in 5 segments** | **0** |

"our model" is §2.1–§2.4 evaluated on the same genotypes: clean-word runs of ≥ 2, endpoints
by §2.2, 3 Mb floor. This is the whole of the over-calling gap, in one table, on a fileset
where the truth is known.

### 5.2 It is not a positive-evidence test

The obvious hypothesis — that KING demands evidence of sharing (het-het excess, a local
kinship estimate) and not merely the absence of IBS0 — is **refuted**. Take a genuinely
IBD0 six-word region and repair every IBS0 in it by turning one member of the pair into a
heterozygote, creating no sharing at all. The reference calls it exactly as it calls a true
IBD1 block, 447 markers, at maf 0.5, 0.3 and 0.2 alike. At high frequency, absence of IBS0
*is* sufficient.

### 5.3 It is per-SNP informativeness, and it comes from the other samples

A pair that is IBD1 by construction across a whole chromosome is not called at all when the
markers are near-monomorphic. The decisive experiment holds the pair's genotypes byte-for-
byte fixed and changes only the **other four samples**:

| pair's genotypes | other 4 samples | reference |
| --- | --- | --- |
| drawn at maf 0.15 | maf 0.50 | called (958 mk) |
| drawn at maf 0.15 | maf 0.15 | called (958 mk) |
| drawn at maf 0.12 | maf 0.50 | **called (958 mk)** |
| drawn at maf 0.12 | maf 0.12 | **dropped (511 mk)** |
| drawn at maf 0.10 | maf 0.50 | **called (958 mk)** |
| drawn at maf 0.10 | maf 0.10 | **dropped (511 mk)** |

Identical pair genotypes, identical map, opposite outcomes. The discriminator is per-SNP
allele frequency estimated from the sample.

`allsegs.txt` is byte-identical between the two — same `N_SNP`, same `Length` — so
uninformative markers are **not** removed from the marker array or from the word grid. The
gate is inside the per-pair scan.

### 5.4 How much informativeness

Blocks in which `n` markers per word are polymorphic (maf 0.5) and the remaining `64 - n`
are monomorphic across all six samples:

```
6-word block:  n = 7  dropped   n = 8  dropped   n = 9  CALLED   n = 10 CALLED
```

A sharp threshold at **9 informative markers per 64-marker word**. Making the filler
markers barely polymorphic instead of monomorphic — exactly one minor allele copy among
the 12 in the sample, maf 0.083 — moves nothing: the threshold is still 9 strong markers,
so those barely-polymorphic markers count as uninformative.

The threshold is **not** a per-word constant, though. Repeating the sweep at other block
widths:

| block width (clean words) | per-word threshold | markers in block |
| ---: | ---: | ---: |
| 3 | 17 | 51 |
| 4 | 14 | 56 |
| 6 | 9 | 54 |
| 10 | 7 | 70 |
| 14 | 3 | 42 |

Neither the per-word count nor the block total is constant, so the aggregation is some
window or density statistic that these five points do not determine. A single weak word
embedded among strong ones is tolerated down to 4 informative markers and only bites at 0,
which is further evidence that the test spans more than one word.

**This is the next thing to nail.** It is the whole of the 188 spurious rows and most of
the 357 rows whose `IBD1Seg` is wrong, and §5.1 is a ready-made acceptance test for any
candidate rule.

---

## 6. What is still open

1. **The exact form of the informativeness statistic** (§5.4) — window shape and
   aggregation. Existence, direction, and order of magnitude are settled; the formula is
   not.
2. **The one-clean-word case** (§2.4) — reproducible, deterministic, not explained by
   length, `--seglength`, or any per-word genotype count tried.
3. **A ±1 marker discrepancy at split points.** At 25 kb spacing, 4 of 12 interior split
   positions report one marker interval more than the abutting-pieces model of §2.2
   predicts (printed `0.7124` where the model says `0.7122`). Real, not rounding —
   1 marker is 1.95 ulp there. Under 1 % of any segment length, but it will block
   byte-parity on some rows.
4. **IBD2 endpoint placement after an interior IBS0**, off by ±1 marker from the model in
   some configurations (§3).

## 7. Actionable now

Independently of the open items, three of these are corrections to
`crates/king-core/src/ibdseg.rs` that are pinned to the byte:

* `PairSegments::reported()` uses `longest_bp > LONG_SEGMENT_BP`; the measurement says
  **`>=`** (§4).
* `keep_long()` drops calls shorter than `min_bp`; the measurement says keep iff
  **`length >= min_bp`** (§4).
* `MIN_RUN2` for IBD2 is 2; the measurement says **1** — a single clean word is a valid
  IBD2 segment — and IBD2 calls must **not** get the one-word right-hand extension that
  `Scan::right_end` applies (§3).
