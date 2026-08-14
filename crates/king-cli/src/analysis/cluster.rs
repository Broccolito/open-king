//! `--cluster`: family clustering, and nothing else.
//!
//! The whole pass is the shared clustering prologue
//! ([`unrelated::clustering_prologue`]) with a tail that only runs when two families were
//! actually joined. On a fileset where none were the tail is empty, so the body is:
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
//! When families *are* joined the prologue's `connected` table is followed by
//!
//! ```text
//! Update-ID information is saved in file <p>updateids.txt
//!
//! Pair-wise relatedness in newly clustered families saved in <p>cluster.kin.
//! KING cluster analysis ends at <t>
//! ```
//!
//! with `<p>updateids.txt` and `<p>cluster.kin` written. Note that `--unrelated` shares
//! everything above this tail and none of it: it prints a second blank line where the
//! `Update-ID` line stands here, and writes neither file.
//!
//! # `<prefix>cluster.kin`
//!
//! Every pair inside a **merged** cluster — cross-family pairs included, which is the
//! whole point — under the new `KING<k>` family ID. Fifteen tab-separated columns: the
//! `.kin` ones without `Z0`/`Phi`, plus `Sex1`/`Sex2` and the four segment columns
//! `--related` also carries. `Kinship` is the **within-family** estimator, because after
//! the merge the pair is inside one family; that is what makes the cross-family rows of
//! `bigish`'s `KING1` print small negatives rather than the between-family estimate.
//!
//! `IBD1Seg`, `IBD2Seg` and `PropIBD` come from the IBD-segment engine, which is not yet
//! exact (`docs/PARITY.md` §11.1); everything else in the file is. That is the only thing
//! standing between `apps/bigish__cluster` and a pass.

use std::fmt::Write as _;
use std::io::Write;

use king_core::{counts, kinship as est, PairCounts, Scope};

use crate::analysis::unrelated::Clustering;
use crate::analysis::{f, ibdseg, out_path, unrelated};
use crate::cli::Options;
use crate::console;
use crate::load::Loaded;

/// Header of `<prefix>cluster.kin`.
const HEADER: &str = "FID\tID1\tID2\tSex1\tSex2\tN_SNP\tHetHet\tIBS0\tHetConc\tHomIBS0\t\
                      Kinship\tIBD1Seg\tIBD2Seg\tPropIBD\tInfType\n";

/// Run the pass. The caller has already printed `Options in effect:`.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let clustering = unrelated::clustering_prologue(opts, loaded, out);
    if !clustering.any_merged {
        return;
    }

    let ids_path = out_path(opts, "updateids.txt");
    let _ = std::fs::write(
        &ids_path,
        clustering.updateids_text(&loaded.fileset.samples),
    );
    let _ =
        out.write_all(format!("Update-ID information is saved in file {ids_path}\n\n").as_bytes());

    let kin_path = out_path(opts, "cluster.kin");
    let _ = std::fs::write(&kin_path, cluster_kin(opts, loaded, &clustering));
    let _ = out.write_all(
        format!("Pair-wise relatedness in newly clustered families saved in {kin_path}.\n")
            .as_bytes(),
    );
    let _ = out.write_all(
        format!(
            "KING cluster analysis ends at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
}

/// Render `<prefix>cluster.kin`: every within-cluster pair of every merged cluster.
fn cluster_kin(opts: &Options, loaded: &Loaded, clustering: &Clustering) -> String {
    let samples = &loaded.fileset.samples;
    let segments = ibdseg::Segments::new(opts, loaded);
    let mut text = String::from(HEADER);
    for (key, members) in clustering.merged() {
        for (n, &a) in members.iter().enumerate() {
            for &b in &members[n + 1..] {
                let c = counts::pair_counts(&loaded.fileset.genotypes, a, b);
                // Zero when the map has no usable segment, which is also when the
                // reference has nothing to report.
                let (pi1, pi2, prop) = segments
                    .as_ref()
                    .map_or((0.0, 0.0, 0.0), |s| s.of(loaded, a, b));
                let _ = writeln!(
                    text,
                    "{key}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    samples[a].iid,
                    samples[b].iid,
                    samples[a].sex,
                    samples[b].sex,
                    c.n_snp,
                    f(est::het_het_prop(&c), 4),
                    f(est::ibs0_prop(&c), 4),
                    f(est::het_concordance(&c), 4),
                    f(hom_ibs0(&c), 4),
                    f(est::kinship(&c, Scope::WithinFamily), 4),
                    f(pi1, 4),
                    f(pi2, 4),
                    f(prop, 4),
                    king_core::ibdseg::inf_type(pi1, pi2, prop),
                );
            }
        }
    }
    text
}

/// `HomIBS0` — `N_IBS0` over the number of variants at which **either** sample is
/// homozygous for A1.
///
/// The same statistic `--related` prints; it belongs in `king_core::kinship` beside the
/// other estimators once both writers are settled.
fn hom_ibs0(c: &PairCounts) -> f64 {
    f64::from(c.ibs0) / f64::from(c.hom_a1_union)
}
