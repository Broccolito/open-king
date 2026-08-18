# Continuation brief

Written for whoever picks this project up next. It says where things stand, what the
three failing cases actually are, and what is worth doing. Read
[`docs/MAINTAINING.md`](MAINTAINING.md) next — it has the working rules and the rigs.

## State as of this handoff

| | |
| --- | --- |
| Parity | **477 / 480** captured reference invocations byte-identical |
| Harness self-check | **480 / 480** (reference vs its own captures — proves no false positives) |
| Rust tests | 330 passing |
| clippy `-D warnings`, `cargo fmt --check` | clean |
| CI | green on ubuntu, macOS, Windows, and the parity job |
| Repo | `Broccolito/open-king` (private), MIT |

Reproduce in two commands:

```bash
cargo build --release
python3 tests/parity/run_parity.py --impl ./target/release/king
```

Add `--impl <path-to-reference-king>` instead to run the self-check; it must score 480/480.
If it does not, the harness is broken, not the program.

## The three failing cases, exactly

All three are on `bigish` (200 samples × 50,000 SNPs), the largest dataset. **None is a
wrong relatedness estimate.** Every kinship, IBS, concordance and segment number this
program computes is byte-identical to the reference, on every dataset, at every
`--seglength` floor.

### 1 and 2 — the two-stage screen's stdout line (2 cases)

```
core/bigish__related_degree2            stdout!=
ibdseg/bigish__related_degree2_ibdseg   stdout!=
```

One line differs, and it is the same line in both cases:

```
-  Stages 1&2 (with 32768 SNPs): 36 pairs of relatives are detected (with kinship > 0.0625)
+  Stages 1&2 (with 32768 SNPs): 50 pairs of relatives are detected (with kinship > 0.0625)
```

**Why it fails.** Before its exhaustive pass, KING runs a cheap screen over
`min(m, 32768)` markers to decide which pairs are worth computing properly. We report 50
pairs where the reference reports 36 — our screen is less selective. The *next* line,
`Final Stage (with 50000 SNPs): 26 pairs ... are confirmed`, already matches, so the
exhaustive pass is correct; only the screen's own count differs.

**Blast radius: one console line. No output file is affected.** `.kin0`'s rows come from
the exhaustive re-estimate that follows the screen and are byte-correct at every
`--degree`, and all 14 pairs the reference's screen drops sit below the 0.08839 reporting
threshold, so they would not have been written anyway.

**What happened to it.** Four rounds of measurement (`docs/research/22-screen.md`). It is
not a lack of effort; the mechanism is genuinely elusive, and the rounds are valuable
mostly for what they *exclude*:

* **Not kinship over any subset or weighting of markers.** This is an algebraic proof, not
  a failed search: at a marker of frequency `p` the `p²q²` terms cancel, so numerator and
  denominator are both proportional to `Σpq` over whatever index set you sum on — making
  every subset and every non-negative weighting unbiased for the same φ.
* **Not a merge or compression into 32,768 slots.** Appending markers one at a time shows
  no step at one marker over budget, though a block merge must step; and the same multiset
  appended, interleaved or shuffled prints an identical count, though every idempotent
  combining operation is lossless for at least one of those arrangements.
* **Not two stages intersected.** Per-pair labels are perfectly sharp — zero inversions.
* **Not rank-block or rank-stride grouping** under or/and/xor/saturating-sum.
* **Not a bound-based early exit** over the sorted map (round 4).

The decisive positive fact, which any surviving hypothesis must explain: **the screen
reads markers it does not keep.** Holding `m = 50000` as 32,768 markers at MAF 0.45 plus
17,232 at MAF `x`, the top-32,768 index set is the same markers with bit-identical
genotypes for all `x` through 0.30 — yet the printed count runs 46, 46, 45, 43, 39, 37.

What it *does* obey: `k_screen = 0.5 + R(k − 0.5)`, with `R` exactly 1 when `m ≤ 32768`
and, above that, depending on the **MAF spectrum** rather than on `(m, n)` — swinging
0.998 to 1.085. **Do not fit `R`.** It would reproduce `bigish` and nothing else, which is
the one thing this project forbids.

Round 4 also found a sharper target than the 2% deflation everyone was chasing: a **knee**,
a factor of three in kinship, reproducing out of sample at three map sizes and two
degrees, one run per point to measure. `docs/research/22-screen.md` §20 lists the next
scans.

### 3 — `build.log` is a strict subsequence (1 case)

```
apps/bigish__build   stdout!=; kingbuild.log!=(num)
```

**Why it fails.** `<prefix>build.log` narrates pedigree reconstruction. Every line we emit
is **byte-identical** to the reference's, and the file is a strict subsequence of it — we
are missing *triggers*, never formatting. On `bigish` we emit 6 of 18 lines. The same
lines are correspondingly missing from stdout.

**Blast radius: that one file.** `--build`'s actual products, `updateids.txt` and
`updateparents.txt`, are byte-identical, so the reconstruction itself is right.

**What happened to it.** Solved along the way: the log's line templates; that the same
bytes go to file and stdout; cluster **numbering** (a staged merge queue worked by
relationship type — Dup/MZ, then PO, then FS, then weaker — with file order breaking ties
inside a type and `OriginalFamID` in absorption order; 19/19 on fresh shapes); the merge
gate (`.kin0`'s own disjunction following `--degree`, not `kinship > 2^-2.5`); and
duplicate-removal (keep the copy with more declared first-degree relatives, ties to the
later id; 27/27).

Still open:

* **The remaining `INFERENCE` triggers.** Which conditions make the reference emit each
  template. Work from constructed pedigrees, not from `bigish`.
* **The `uncle|aunt` vs `grandparent|HS|nephew` cut** on `Join3/Join2`, bracketed to
  (0.846, 0.902) — needs tightening.
* **`FS1` member order** is a hash-table iteration order over a family-scoped, id-keyed
  container. It *is* a function of the id strings but *not* a per-id ranking (13 subsets of
  one 8-id pool contradict each other 91 times). Reproducing it means identifying the hash.
  Worth 3 of 59 shapes, and `bigish` has no `FS1` line — **deprioritise it.**

One retraction to respect: the claim that our sib-pair over-call fully accounts for the
`Join3/Join2` residual was measured **in sample** (39 of 39 triples) and is **11 of 34**
out of sample. Closing the over-call is necessary, not demonstrably sufficient.

## Working rules that earned their place

These are not style preferences; each was bought with wasted work.

1. **Never fit to the corpus.** A prior session reached `IBD2Seg` 982/982 by fitting, and
   it had to be thrown away and re-derived. Validate on fresh canvases with unused seeds.
2. **Land a rule only if exact rows *and* mean error both improve** (or MAE equal). This
   caught a wrong `.seg` port that looked like progress: better tail, worse mean — the
   signature of a partly-right rule.
3. **A negative result is a real result.** Several rounds deliberately shipped nothing.
   The impossibility proofs above are worth more than a fitted constant would have been.
4. **The reference wins over intuition.** Verify by running it, not by reasoning.
5. **Never read KING's source.** Its licence is not MIT-compatible; that is why this repo
   can be MIT at all.
6. **Synthetic filesets must code A1 as the minor allele**, or KING aborts with
   `Too many first alleles as the major allele` — which silently turns a bisection into
   "no bracket" and has cost hours more than once.

## The instrument you will not invent yourself

The **canvas read-back**: choose marker spacing so one ulp of a printed 4-decimal column is
~0.11 of a marker gap. The printed value then reads back the *exact* number of marker
intervals called, turning a lossy aggregate into an exact instrument. Sweeping
`--seglength` as a continuous parameter makes the jumps of `IBD2Seg(L)` the individual call
lengths. Rigs: `docs/research/fixtures/{segcanvas,ibd1canvas,fringecanvas,mergelab,push1,
screenfold}.py`, each with a JSON cache of reference answers.

The **ladder fileset** (`screenfold.py`): 48 pairs climbing through a cutoff, so one run
reads the effective threshold — 17× cheaper than bisection.

## Traps a newcomer will hit

* **Genotype files are never committed.** A fresh checkout has the golden text but no
  corpus; regenerate with `tests/parity/generate_corpus.py`. CI does this, and both
  scorecard tests assert they measured all 10 datasets so a missing corpus fails loudly
  instead of passing having measured nothing.
* **`.gitattributes` sets `* -text`.** Without it, Windows checks the byte-exact fixtures
  out as CRLF and every console test fails for a reason unrelated to the program.
* **Parity is measured against one build**: KING 2.3.2, macOS arm64. Segment numerics
  changed across 2.1.x and 2.2.x, so agreement with another build is not implied.
* **Differences the corpus cannot see** are listed in `docs/PARITY.md` §5.10–§5.12 and
  §4.6 — the sparse-panel segment fallback is the one most likely to bite a real user.

## If you want to keep going, in value order

1. **`build.log` triggers** — bounded, well-understood, closes 1 case.
2. **The screen's knee** — closes 2 cases; hardest thing in the project, and four rounds
   of exclusions mean the next attempt starts far ahead of the first.
3. **The corpus-invisible differences** in `docs/PARITY.md` §5.10–§5.12 — these affect real
   users even though they cost no case, which arguably makes them more valuable than
   either of the above.

## Open work

Tracked as GitHub issues:

| # | Title | Parity cases |
| --- | --- | --- |
| [1](https://github.com/Broccolito/open-king/issues/1) | `build.log` missing INFERENCE triggers | 1 |
| [2](https://github.com/Broccolito/open-king/issues/2) | Two-stage screen pair count | 2 |
| [3](https://github.com/Broccolito/open-king/issues/3) | Differences the corpus cannot see | 0 — but they affect real users |
| [4](https://github.com/Broccolito/open-king/issues/4) | Unimplemented analyses exit 0 silently | 0 |

Issue 3 is the one worth doing first: it costs no parity case, which is exactly why the
suite cannot protect a user from it.
