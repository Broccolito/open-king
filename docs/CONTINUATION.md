# Continuation brief

Current handoff for the minimal open-king relatedness/QC product. Read
[`SCOPE.md`](SCOPE.md) for the binding product boundary and [`PARITY.md`](PARITY.md) for
the measurements behind every claim.

## Current state

| Gate | Current evidence |
| --- | --- |
| Captured differential | **480 PASS / 0 FAIL / 480**, 876 output files byte-compared, 8 documented diff-exclusions |
| Regression baseline | `MATCH (480 cases)` |
| Rust verification | all workspace tests pass; clippy `-D warnings`, formatting, release build, and `open-king-core` docs pass |
| Library surface | typed `Bundle -> RelatednessReport` API in `open-king-core`; see [`API.md`](API.md) |
| Dense-panel `--ibdseg` | byte-identical to KING 2.3.2 on a 663,197-marker x 157-sample panel (v0.1.1; [`PARITY.md`](PARITY.md) §5.14) |
| Live issues | only [#11](https://github.com/Broccolito/open-king/issues/11), the five held-out supported-core residuals |

Reproduce the primary gates:

```bash
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo build --release
python3 tests/parity/run_parity.py --impl ./target/release/open-king --baseline
```

## Completed development-plan milestones

The implementation sequence is intentionally preserved as small commits:

| Commit | Result |
| --- | --- |
| `8723ef0` | binding minimal-product scope and exclusions |
| `c6435b9` | case-insensitive sample identity |
| `a93c33d`, `54e38ba` | sparse kinship-only fallback across related/ibdseg/unrelated/cluster/build |
| `b69d250` | chromosome/position map-order validation |
| `ba38dbd` | A1-minor orientation gate |
| `42d6c9d` | closed 100 Mb usable-segment floor |
| `8a4c12c` | conditional `splitped.txt` generation |
| `9d075ee` | held-out IBD1/IBD2 long-pair filtering and zero row-set residual |
| `d26ba24` | primary `build.log` reconstruction parity |
| `1362a54` | exact dense and sparse relatedness screening stages |
| `e89a66d` | fatal pre-I/O diagnostics for excluded product requests |
| `0c582a9` | typed owned relatedness API |
| `9863c75` | permissive nonstandard `.fam` `SEX` parsing |

Issues #1, #2, #4, #5, #6, #7, #8, #9, and #10 are closed with differential evidence.

## Deliberate product exclusions

Population structure, ROH, GRM, association/risk analyses, PLINK orchestration, R plotting,
comma-separated multi-fileset merging, and a strict `--cpus` worker cap are not part of the
minimal product. Their recognized CLI spellings fail clearly before input is opened. These
are product boundaries, not remaining parity defects.

## Settled, and not carried forward

**Exact-multiple-of-64 segment tail.** Four value differences among 6,713 rows in the
24-fileset battery, all at exactly 40,000 markers; the 39,999 and 40,001 controls are exact
and the row sets are exact. This is a deliberate safety divergence: KING reads uninitialized
memory there, which safe Rust must not emulate. It is retained on purpose and asserted by
`tests/parity/probes/segment_residuals.py::check_safe_tail_divergence`, so it is not open
work and is not part of issue #11.

## Remaining supported-core work

The 480-case corpus is saturated, so every remaining item needs a held-out discriminator.
All five are tracked by issue #11.

1. **Segment acceptance gate:** one constructed brother/sister pair is emitted by open-king
   and omitted by KING. It does not occur in the captured corpus. Since v0.1.1 the `>=10`
   rule it contradicts is known to be the *first row* of a table the reference reads off
   the total marker count (`PARITY.md` §5.14), but that does not explain this pair: its
   fileset carries 5,000 markers, which is that same first row, and open-king's answer on
   it is unchanged. The discriminator this still needs must live inside the first tier.
2. **`HomIBS0` exact ties:** zero golden rows differ; two of 1,189 random-pedigree rows and
   four of nine hand-placed exact ties differ in the last printed digit. Integer counts and
   the algebraic ratio are correct; no tested floating evaluation reproduces all reference
   perturbations.
3. **`MI_Removal`:** the greedy pair-error cover matches six of seven focused probes, while
   every golden row is zero. The seventh reference flag at a 0.18% rate refutes the current
   monotone threshold rule.
4. **Sparse PO/FS cutoff:** the segment-unavailable fallback currently uses `0.0050`, while
   reference probes show a deterministic data-derived value in roughly `[0.0035, 0.0060]`.
   The application rule is known, but the value derivation is not; a held-out fileset that
   changes the reference cutoff must pin the rule before replacing the constant.
5. **Rare pedigree reconstruction shapes:** the cached `build.log` replay is 277/347
   byte-identical; another 52 have the same distinct lines with repetition/count residue.
   Remaining semantic work includes cross-family named-parent materialization and
   `<FID>-><IID>` renaming plus rare inference-loop trigger/repetition shapes. Two cached
   reference logs are truncated evidence and should not be fitted.

Do not weaken the headline by mixing these held-out results into the 480-case denominator,
and do not hide them because the headline is green. A fix should add a probe that rejects a
plausible wrong rule, pass the held-out comparison, preserve the 480-case baseline, and then
update issue #11 and `PARITY.md`.

## Working rules

* Never read KING source; this is a clean-room MIT implementation.
* Never fit a constant only to the golden corpus. Use fresh canvases or constructed
  filesets and keep the counterexample in-tree.
* Treat a negative or unsafe-reference result as evidence when the reproducer pins it.
* Preserve exact output bytes, row sets, console, stderr, exit status, and file presence.
* `tests/parity/run_parity.py --baseline` must match, not merely improve.
