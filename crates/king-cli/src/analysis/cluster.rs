//! `--cluster`: family clustering, and nothing else.
//!
//! The whole pass is the shared clustering prologue
//! ([`unrelated::clustering_prologue`]) with a short tail. On a fileset where no two
//! families are joined that tail is empty, so the body is:
//!
//! ```text
//! Family clustering starts at <t>
//! Autosome genotypes stored in <w> words for each of <n> individuals.
//! Sorting autosomes...
//! Total length of <k> chromosomal segments usable for IBD segment analysis is <d> Mb.
//!   Information of these chromosomal segments can be found in file <p>allsegs.txt
//!
//! <n> CPU cores are used to compute the pairwise kinship coefficients...
//! Clustering up to 1st-degree relatives in families...
//! Individual IDs are unique across all families.
//! No families were found to be connected.
//! ```
//!
//! and `<p>allsegs.txt` is the only file written. Under ten samples the prologue
//! collapses to `This function is currently disabled for tiny dataset with sample size <
//! 10.` and even that file is skipped.
//!
//! When families *are* joined the reference goes on to write `<p>updateids.txt` and
//! `<p>cluster.kin` and to print
//!
//! ```text
//! Update-ID information is saved in file <p>updateids.txt
//!
//! Pair-wise relatedness in newly clustered families saved in <p>cluster.kin.
//! KING cluster analysis ends at <t>
//! ```
//!
//! **That tail is unimplemented.** Merging needs a hundred samples before the reference
//! will even look at a cross-family pair, so `bigish` is the only capture that reaches
//! it; see `docs/PARITY.md` §11.

use std::io::Write;

use crate::analysis::unrelated;
use crate::cli::Options;
use crate::load::Loaded;

/// Run the pass. The caller has already printed `Options in effect:`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let _ = unrelated::clustering_prologue(opts, loaded, out);
}
