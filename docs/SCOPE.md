# Product scope

open-king is the smallest useful, permissively licensed implementation of KING's
relatedness workflow. It is not intended to reproduce every analysis that accumulated in
the original KING executable.

## Supported core

The supported product surface is:

* PLINK1 `.bed` / `.bim` / `.fam` input;
* pairwise kinship, relatedness classification, duplicate detection, IBS statistics and
  autosomal/X-chromosome IBD segments;
* unrelated-sample selection, family clustering and pedigree reconstruction;
* per-sample, per-marker and automatic QC reports; and
* the output files and command-line behavior needed to substitute open-king for those
  workflows.

The top-level [README](../README.md) gives the task-oriented command table. [CLI.md](CLI.md)
documents every accepted option, and [PARITY.md](PARITY.md) is the authoritative measured
comparison with KING 2.3.2.

## Deliberately excluded

The following are legacy or non-core analysis families and are not planned for this
minimal package:

* population structure: `--pca`, `--mds`;
* runs of homozygosity: `--roh`;
* GRM, association and risk analysis: `--makeGRM`, `--lmm`, `--tdt`, `--gdt`, `--risk`;
* PLINK-export orchestration: `--plink`;
* R-based plotting: `--rplot`, `--pngplot`, `--rpath`; and
* KING's comma-separated multi-fileset merge mode.

Their spellings remain in the parser and banner to preserve KING 2.3.2's recognizable
command-line surface. The executable rejects an excluded analysis, its associated parameter,
or comma-separated input with exit status 1 before opening any input file. This prevents a
requested workflow from appearing to succeed without producing its promised output. Use a
dedicated population-structure, association, ROH, plotting or fileset-merging tool before or
after open-king as appropriate.

These exclusions are product boundaries, not defects in the supported relatedness core.
They must not be counted as relatedness-parity failures. A future expansion requires an
explicit scope decision, clean-room behavioral evidence and its own tests; merely accepting
an option is never evidence that the analysis exists.

## Computing scope

`--cpus <n>` is accepted, echoed and used in compatibility console text. It does not
currently promise a strict Rayon worker-thread cap. Computed files are deterministic across
thread counts, and exact performance parity with KING is not a product goal.

## What still counts as a defect

A difference remains in scope when it changes a supported command's accepted input,
relationship result, report contents, output-file set or documented console contract.
Examples include a segment acceptance-gate counterexample, rare `HomIBS0` ties,
the data-derived sparse PO/FS cutoff, `MI_Removal`, and unusual reconstruction pedigrees.
Those held-out differences are tracked in [PARITY.md](PARITY.md) and GitHub issue #11 and
remain work even though all 480 captured cases pass and the excluded analyses above do not
count as product defects.
