# Typed Rust API

The `open-king-core` crate exposes the minimal relatedness product without console output or
text-file parsing. The owned flow is:

```text
PLINK1 .bed/.bim/.fam -> Bundle -> RelatednessReport
```

```rust
use open_king_core::{Bundle, BundleError, RelatednessOptions};

fn main() -> Result<(), BundleError> {
    let bundle = Bundle::from_plink("cohort.bed")?;
    let report = bundle.relatedness(&RelatednessOptions::default());
    println!("{} pairs", report.pairs.len());
    Ok(())
}
```

`Bundle::from_plink` loads the autosomal analysis matrix. `Bundle::new` adopts an existing
`open_king_io::Fileset`, and `Bundle::fileset` keeps the validated low-level data reachable.

## Report contents

`RelatednessReport` owns the sample table and one `RelatedPair` for every unordered pair,
in ascending sample-index order. `report.pair(a, b)` performs symmetric lookup. Each pair
contains:

* sample indices and within-/between-family estimator scope;
* exact `PairCounts`, including pairwise-missing-aware SNP, heterozygote, IBS0/1/2 and
  homozygote counts;
* `PairStatistics`: kinship, mean IBS, distance, concordance, HetHet, IBS0 and HomIBS0;
* pedigree-expected kinship and Z0;
* a typed `Relationship`; and
* optional raw and normalized IBD1/IBD2 segment statistics, their denominator, long-call
  reporting gate and segment inference label.

These are the unrounded values. Formatting and row/file selection remain presentation
policy in the compatibility CLI. The lower-level `counts`, `kinship`, `infer`, and `ibdseg`
modules remain public when a caller needs a single kernel rather than an all-pairs report.

## Options

`RelatednessOptions::default()` uses the measured KING 2.3.2 relationship cutoffs, derives
the parent-offspring IBS0 threshold from declared full siblings, and scans segments at the
3 Mb default floor.

Set `segment_length_bp` to `None` when only counts and kinship are needed. This avoids the
segment caller while retaining the full pair table. `PoThreshold::Fixed` supplies an
explicit strict IBS0 boundary; `PoThreshold::FromPedigree` uses half the mean among declared
full-sibling anchors, and treats exact IBS0 zero as parent-offspring when none exist.

## Validation and failures

`BundleError` distinguishes I/O failures from malformed in-memory data. Construction checks
sample/variant dimensions, bit-plane lengths, unused tail bits, retained map indices,
recognized chromosomes and map order. Once a `Bundle` exists, `relatedness` is infallible.

The report contains `n(n-1)/2` pairs. Counts and optional segment calls are therefore
quadratic in sample count; disabling segment scanning lowers the constant but does not make
an all-pairs report subquadratic.
