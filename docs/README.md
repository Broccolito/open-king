# open-king documentation

Everything about this project, grouped by what you are trying to do. The project itself is
introduced in the [top-level README](../README.md).

## If you want to run it

Read in this order. Each page assumes you have PLINK1 filesets and know what a kinship
coefficient is; none of them assumes anything about this codebase.

| page | answers |
| --- | --- |
| [`CLI.md`](CLI.md) | Which flag do I need, and what does it do? Every option the parser accepts, which analyses it affects, the parser quirks that will surprise you, and the exit statuses |
| [`OUTPUTS.md`](OUTPUTS.md) | I have the files — what is in them? Every output file: column meanings, numeric formats, the three different row orders, and the rules for when a file is empty, truncated or absent |
| [`INTERPRETING.md`](INTERPRETING.md) | What does 0.177 mean, and when is it lying to me? The two estimators, the cutoffs, PO vs FS, the `Error` column, and seven demonstrated ways the numbers mislead |
| [`COOKBOOK.md`](COOKBOOK.md) | Just show me the command. Twelve task-oriented recipes: find duplicates, screen a cohort, pick an unrelated set, QC a fileset, migrate from KING |
| [`SCOPE.md`](SCOPE.md) | What is this package intentionally responsible for? The minimal relatedness/QC core, deliberately excluded analyses, and computing scope |

Every command in all four pages was executed against `target/release/open-king` and its output
pasted from that run.

**Start here for common tasks:** find duplicate samples →
[COOKBOOK recipe 1](COOKBOOK.md); screen a GWAS cohort for cryptic relatedness →
[recipe 2](COOKBOOK.md); a number looks wrong →
[INTERPRETING §7](INTERPRETING.md#7-pitfalls); a file you expected
is missing → [OUTPUTS: empty, header-only, truncated, absent](OUTPUTS.md#empty-header-only-truncated-absent).

## If you want to know how faithful it is

| page | contents |
| --- | --- |
| [`PARITY.md`](PARITY.md) | **The authoritative claim.** 480 of 480 captured reference invocations byte-identical; the analysis × dataset matrix, per-file and per-row scorecards, and labelled held-out limitations (§4.6, §5.10–§5.12) |
| [`RELEASE_NOTES.md`](RELEASE_NOTES.md) | what changed in this release |
| [`CONTINUATION.md`](CONTINUATION.md) | the handoff brief: completed milestones, deliberate exclusions, and the remaining held-out work |
| [`SCOPE.md`](SCOPE.md) | the product boundary that separates deliberate exclusions from supported-core parity defects |

## If you want to work on it

| page | contents |
| --- | --- |
| [`SPEC.md`](SPEC.md) | the implementation specification: CLI surface, I/O, every analysis |
| [`BEHAVIOR.md`](BEHAVIOR.md) | the black-box experiments behind the rules — SNP inclusion, `--degree` semantics, the ID sort comparator, output-file existence, column-set variation |
| [`VERIFIED_FORMULAS.md`](VERIFIED_FORMULAS.md) | every estimator and every printed column, checked numerically against the reference |
| [`MAINTAINING.md`](MAINTAINING.md) | the clean-room rule, the test corpus, re-capturing goldens, and the fixture technique the segment work depends on |
| [`research/`](research/) | 25 investigation notes, roughly chronological: `01-paper-estimators.md` through `23-gap-bound.md`. `research/SPEC.md` is the original clean-room spec; `research/fixtures/` holds the probe scripts each note names |
| [`reference-captures/`](reference-captures/) | recorded reference-binary transcripts used by the parity harness |

**The clean-room rule** ([`MAINTAINING.md`](MAINTAINING.md) §1) is absolute: KING's source is
never read. Every rule in this project traces to a published description, a documented file
format, or a named black-box experiment.

## Reproducing anything in these docs

```bash
cargo build --release                                           # -> target/release/open-king
python3 tests/parity/generate_corpus.py --outdir /tmp/kingdocs  # 13 filesets, ~6 s
```

The corpus is 13 synthetic PLINK1 filesets with known-by-construction pedigrees, generated
from one seed with the Python standard library. `docs/CLI.md`
[§10](CLI.md#10-the-derived-filesets-used-above) and `docs/INTERPRETING.md`
[appendix](INTERPRETING.md#appendix--reshapepy) publish the two scripts that derive the
deliberately awkward filesets a few examples need.

- [CONTINUATION.md](CONTINUATION.md) — handoff brief: state, the three open cases, working rules, next steps.
