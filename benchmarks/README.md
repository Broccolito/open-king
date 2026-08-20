# open-king benchmarks

`run_benchmark.py` measures wall-clock time, peak resident set size and CPU time
for a matrix of (dataset, analysis) pairs, for one or two KING-compatible
binaries. When a second binary is given it also byte-compares every analysis
output file the two produced for the same invocation.

Python 3 standard library only. POSIX only: it uses `os.posix_spawn`, `os.wait4`
and `resource`.

## Running it

```bash
# open-king alone, results into benchmarks/results/
python3 benchmarks/run_benchmark.py --work-dir /tmp/king-bench

# open-king against a reference KING build, with output diffs
python3 benchmarks/run_benchmark.py \
    --work-dir /tmp/king-bench \
    --binary-a target/release/king      --label-a open-king \
    --binary-b /path/to/king-2.3.2      --label-b king-2.3.2

# just build the input filesets and print their sizes
python3 benchmarks/run_benchmark.py --work-dir /tmp/king-bench --gen-only
```

`--work-dir` must sit outside the repository. Generated filesets and every run's
output files go there, so the repository stays clean; the harness warns if the
work dir is inside the tree. The only things written into the repository are
`results/results.json` and `results/results.md`.

Useful flags:

| flag | effect |
| --- | --- |
| `--reps N` | timed repetitions per cell (default 3) |
| `--warmup N` | untimed runs before timing, to warm the page cache (default 1) |
| `--datasets ...` | restrict to named datasets |
| `--analyses ...` | restrict to named analyses |
| `--large NAME:SAMPLES:MARKERS` | override the synthetic size ladder |
| `--matched-cpus N` | thread count used by the `matched` calibration mode |
| `--no-calibrate` | skip the thread-count sweep |
| `--regen-large` | rebuild the synthetic filesets even if cached |

## What it measures

Each run is spawned with `os.posix_spawn` and reaped with `os.wait4`, so the
`rusage` returned belongs to that one child and nothing else. This matters:
`resource.getrusage(RUSAGE_CHILDREN)` reports `ru_maxrss` as a running maximum
over every child the process has ever reaped, so differencing it across a run
gives the wrong answer as soon as one earlier child was larger.

Per cell the harness reports:

- **Wall time**: median across the timed repetitions, with min and max.
- **Peak RSS**: `ru_maxrss` of the child, converted to MB.
- **CPU time**: `ru_utime + ru_stime`.
- **cpu/wall**: CPU seconds divided by wall seconds. On an idle host this reads
  as the average number of cores kept busy. On a loaded host it is pushed below
  1.0 by descheduling and says nothing about the binary.

### The `ru_maxrss` unit trap

`ru_maxrss` is in **bytes on macOS** and **kilobytes on Linux**. The harness
detects the platform, converts accordingly, and records which unit it used in
`results.json` under `machine.ru_maxrss_unit`. Check that field before comparing
memory numbers taken on different hosts. Getting this wrong is a factor of 1024.

### Load average

The harness samples the load average before and after every cell and records the
range in the results. If the peak exceeds half the logical core count,
`results.md` prints a contention warning and points the reader at the CPU-seconds
column instead of wall time. CPU time counts work done rather than time elapsed,
so it degrades far more gracefully under contention. Wall times are an upper
bound on what an idle host would show. They are only comparable across cells to
the extent the load held steady, so `results.json` records the load average
before and after every cell: check the spread there before comparing two rows
from a run that was made on a busy machine.

Check `uptime` before trusting wall-clock numbers from any run.

## Input filesets

**Golden corpus.** The differential-test filesets in `tests/parity/golden/`:
`trio`, `nuclear`, `multifam`, `sexchr`, `missing`, `monomorphic`, `admixed`,
`dups`, `unrelated`, `bigish`. These are small, from 1 to 200 samples. They are
not committed; the harness regenerates any that are missing by shelling out to
`tests/parity/generate_corpus.py`. They are included because they are the parity
fixtures, not because they produce interesting timings.

**Synthetic ladder.** Three larger filesets built by the harness into the work
directory:

| dataset | samples | markers |
| --- | ---: | ---: |
| `synth_s` | 200 | 100,000 |
| `synth_m` | 400 | 100,000 |
| `synth_l` | 800 | 100,000 |

### Why these sizes

Relationship inference is O(n<sup>2</sup> &middot; m): quadratic in samples,
linear in markers. Sample count is therefore the expensive axis and marker count
the cheap one. The ladder holds markers fixed at 100,000 and doubles samples, so
each rung should cost about four times the one below it. That makes the
quadratic term visible directly in the results table, and it buys meaningful
runtimes without the blowup that more samples would cause.

The sizes were chosen against a wall-clock budget: the whole suite, across both
binaries, every dataset and every repetition, should fit in roughly 40 minutes so
that it stays rerunnable. Two larger configurations were measured and rejected
for being too slow to rerun:

| configuration | slowest single run |
| --- | --- |
| 2,000 samples x 50,000 markers | `--cluster` at 79 s |
| 1,000 samples x 100,000 markers | `--unrelated` at 76 s |

At 800 x 100,000 the slowest analysis lands in the tens of seconds, which is
enough that process startup is negligible while still leaving the suite
rerunnable. Repetitions default to 3 for the same reason. Median with min and max
is reported either way, so the spread is visible.

If you have time to spend, `--large synth_xl:2000:100000` adds a rung.

### What the synthetic filesets contain

Every sample gets its own FID and no declared parents. The relatedness is real
but undeclared, which is the cohort shape KING's inference targets. Underneath,
the genotypes carry a full pedigree: parent-offspring pairs, full sibs, half
sibs, grandparents, avuncular pairs, first cousins, and duplicate samples with a
low genotyping error rate. IBD segments come from a recombination process over
the `.bim` genetic map, so `--ibdseg` has real segments to find.

Map construction, the allele frequency model and `.bim` writing are reused from
`tests/parity/generate_corpus.py`. Only the genotype inner loop is replaced. That
generator draws allele by allele in Python, which costs about 6 s for 200 samples
x 50,000 markers and does not scale to this ladder. The harness instead
simulates the whole cohort per SNP with big-integer bit vectors: bit *f* of a
vector is family *f*, so inheritance is a whole-cohort bitwise select and a
per-SNP allele draw costs 16 `getrandbits` calls regardless of sample count. The
result is roughly 300 times faster; all three rungs generate in about 18 s.

The simulator was checked against theory by comparing inferred kinship for every
designed relationship class:

| relation | expected phi | observed mean |
| --- | ---: | ---: |
| duplicate | 0.5000 | 0.4972 |
| parent-offspring | 0.2500 | 0.2487 |
| full sib | 0.2500 | 0.2636 |
| half sib | 0.1250 | 0.1186 |
| grandparent | 0.1250 | 0.1230 |
| avuncular | 0.1250 | 0.1316 |
| half avuncular | 0.0625 | 0.0559 |
| first cousin | 0.0625 | 0.0647 |

Parent-offspring shows near-zero variance and the others show real segregation
variance, which is the expected signature of correct meiosis.

Filesets are cached: a rung is rebuilt only if its size or seed changed, or
`--regen-large` is passed. The seed defaults to 20260819, so runs are
reproducible.

## Adding the reference binary

The harness is built to compare two binaries. Pass the reference with
`--binary-b` and label it with `--label-b`:

```bash
python3 benchmarks/run_benchmark.py \
    --work-dir /tmp/king-bench \
    --binary-b /path/to/king-2.3.2 --label-b king-2.3.2
```

**Do not copy the reference binary into the repository.** It is referenced by
path and executed in place. Record where the build came from in the run notes:
`results.json` captures the path and the first line of its banner under
`binaries[]`, which is enough to identify a build after the fact, but not enough
to say where it was obtained.

Nothing else changes. The comparison columns in `results.md` are always present
and simply blank when only one binary is measured, so adding the reference does
not change the table shape.

### How the published results obtained their reference

Recorded here because "KING 2.3.2" is not one binary, and the results are only
readable if you know which one was used.

KING publishes a precompiled macOS binary at `Mac-king.tar.gz`. It is x86_64,
and on an arm64 host it does not launch: it asks for `libgomp.1.dylib` from an
Intel Homebrew GCC 9. Running it under translation would have meant timing an
emulated program against a native one, which is not a comparison worth
publishing. So the reference was compiled from the published source:

```bash
curl -LO https://www.kingrelatedness.com/KINGcode.tar.gz   # sha256 b6c636ac99d2...
tar -xzf KINGcode.tar.gz
g++-16 -lm -lz -O2 -fopenmp -o king-ref *.cpp              # Homebrew GCC 16, arm64
```

The resulting binary self-reports `KING 2.3.2 - (c) 2010-2023 Wei-Min Chen`.

**The source was compiled and then deleted, and was never read.** open-king is a
clean-room implementation and `docs/MAINTAINING.md` §1 forbids reading KING's
source. Compiling a binary to measure against is not reading it, but the
distinction only holds if the source does not stay on disk waiting to be opened.
If you repeat this, delete the source when the build finishes.

**This rebuild is not the binary the parity goldens were captured from.** The two
agree on every output file in the corpus and disagree on the banner and on one
option whose value was never initialised. That is measured in `docs/PARITY.md`
§5.13, and it matters here: a rebuilt reference can answer `--noscreen` questions
differently from the capture binary.

**Do not commit the reference binary.** KING's terms are "Feel free to use KING
for your research, but please do not redistribute AND make profits." Referencing
it by path keeps this repository clear of that question entirely.

### Output comparison

For each (dataset, analysis) cell both binaries run the same invocation into
separate output directories, and every file is compared by sha256. Each file is
reported as `identical`, `different`, `only_in_a` or `only_in_b`.

Only analysis output files are compared. Each run's stdout and stderr are
captured to files prefixed with `_` and excluded, because the KING banner carries
a wall-clock timestamp and would differ on every run even when the results are
identical. If an output file itself embeds a timestamp it will show as
`different`; check the diff before reading that as a numerical disagreement.

### Thread-count calibration

Running two binaries at their own default thread counts compares two different
things at once. The harness therefore also runs a small sweep on one dataset
(`synth_m` by default) across three settings:

- `default`: no `--cpus` flag, each binary picks for itself.
- `matched`: both pinned to `--cpus N` (`--matched-cpus`, default 4), so the
  comparison is not confounded by that choice.
- `single`: `--cpus 1`, which gives the serial cost.

Dividing single-thread wall time by default wall time gives the parallel speedup
actually realised, which is a more reliable parallelism signal than `cpu/wall`
on a busy host.

The sweep is restricted to one dataset and a few analyses on purpose: a
single-threaded pass over the whole matrix would dominate the budget. Widen it
with `--calibrate-dataset`, `--calibrate-analyses` and `--calibrate-modes`, or
turn it off with `--no-calibrate`.

## Analyses

`--kinship`, `--ibs`, `--related`, `--ibdseg`, `--unrelated`, `--duplicate`,
`--build`, `--cluster`. All eight run on all datasets. Cells that exit non-zero
are listed in their own table in `results.md` rather than being dropped.

Some cells produce no output files legitimately: `--duplicate` writes nothing on
a dataset with no duplicate pairs, for instance. The "Output files produced"
table records what each cell actually wrote.

## Results

`results/results.json` is the full record: machine context, binary banners, every
individual repetition, per-cell load averages, and the per-file comparison.
`results/results.md` is the readable table. Both are overwritten on each run.

Machine context captured: CPU model, physical and logical core counts, RAM, OS
version, architecture, Python version, open-king git commit and whether the tree
was dirty, the detected `ru_maxrss` unit, and load average at start and end.
