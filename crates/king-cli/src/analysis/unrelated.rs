//! `--unrelated` — extract a maximal set of mutually unrelated individuals.
//!
//! Two files, both **tab** separated and both without a header line:
//! `<prefix>unrelated.txt` holds the kept `FID IID` pairs and
//! `<prefix>unrelated_toberemoved.txt` the complement. Together they partition the
//! sample set exactly.
//!
//! # What "related" means here
//!
//! A pair is an edge of the graph when **either** the declared pedigree **or** the
//! genotypes make it closer than unrelated, i.e. kinship above `2^-5.5`:
//!
//! * a mosaic fixture stepping a pair's kinship through 0.0202, 0.0213, 0.0225 and 0.0249
//!   is cut at the last two and kept at the first two — the boundary is the 4th-degree
//!   band edge, not the 1st-degree one the console line mentions;
//! * a declared parent → child → grandchild chain of unrelated genotypes still loses the
//!   two descendants, so the pedigree alone makes edges. Extending it to seven
//!   generations keeps exactly the two ends, which are `2^-7` apart: the same threshold
//!   applies to pedigree kinship.
//!
//! `--degree` does **not** move it. Running `--degree 1 … 4` over a fixture with pairs at
//! 1st, 2nd and 4th degree gives four byte-identical answers, which is why the golden
//! corpus's `*__unrelated` and `*__unrelated_degree2` captures are identical files.
//!
//! # Which pairs are looked at
//!
//! Only pairs **inside one family**, where "family" is the cluster the reference builds
//! first. Cross-family relatives are invisible to the selection: `dups`' exact duplicate
//! pair spans two FIDs and both copies survive, while its within-family MZ pair loses one.
//!
//! Clustering merges families joined by an inferred 1st-degree pair, but **only when the
//! dataset has at least 100 samples** — a 99-sample fileset carrying a 1st-degree
//! cross-family pair reports `No families were found to be connected.` and keeps both
//! members; adding one unrelated sample to make 100 merges the families and drops one.
//! That single threshold explains the whole corpus: `bigish` (200) merges three family
//! pairs, `admixed` (40), `unrelated` (30) and `multifam` (20) merge nothing despite
//! carrying real cross-family parent–offspring and full-sibling pairs.
//!
//! # Under ten samples
//!
//! The reference prints
//! `This function is currently disabled for tiny dataset with sample size < 10.` and
//! skips the load-time preamble, the segment pre-pass and the clustering entirely — but
//! still writes both files, selecting on the pedigree alone.
//!
//! # The visit order
//!
//! The kept list is written in the order the greedy picked, and that order is neither
//! `.fam` order nor sorted order. What is established:
//!
//! * it depends only on each member's **rank within its family under the ID comparator**
//!   — renaming a family's members so their ranks permute permutes the output identically,
//!   while moving their `.bed` rows (genotypes and IDs together) changes nothing;
//! * members are visited in ascending order of their **number of relatives inside the
//!   family**, which is what makes a three-generation family emit its least-connected
//!   founder first;
//! * members with equal counts are visited in a scrambled order that is a pure function
//!   of the family size when *every* member ties — measured for sizes 2…70 against
//!   families of mutually unrelated founders, and reproduced identically by two families
//!   of the same size in one dataset and by families of the same size in datasets of
//!   different totals.
//!
//! [`TIE_ORDER`] is that measurement. It is an empirical table, not a derivation: the
//! generating rule is a sort of the relative-count array whose exact algorithm remains
//! unresolved, and probe fixtures with two distinct counts show the scramble is not
//! confined to each tied run. See the note on [`visit_order`].

use std::io::Write;

use king_core::{counts, kinship, Scope};
use king_io::Sample;

use crate::analysis::{band, king_id_cmp, out_path};
use crate::cli::{Opt, Options};
use crate::console;
use crate::load::Loaded;

/// Below this sample count the reference disables clustering and selects on the pedigree.
const TINY_DATASET: usize = 10;

/// At or above this sample count the reference screens cross-family pairs and merges
/// families; below it, no pair spanning two FIDs is ever looked at.
const CLUSTERING_MIN_SAMPLES: usize = 100;

/// Kinship above which a pair is an edge: the 4th-degree band edge, `2^-5.5`.
const EDGE: f64 = band::FOURTH;

/// Kinship above which two families are merged into one cluster: the 1st-degree band
/// edge, `2^-2.5`, matching the console's `Clustering up to 1st-degree relatives`.
const CLUSTER_EDGE: f64 = band::FIRST;

/// IBS0 proportion below which a 1st-degree pair is called parent–offspring rather than
/// full siblings, for the screening summary alone.
///
/// The reference derives this from the data and never writes it to a file; the value it
/// prints on ordinary data is 0.0055 (`docs/BEHAVIOR.md` §Q2). Nothing in the corpus
/// discriminates: every 1st-degree cross-family pair the screening reports is a
/// full-sibling pair with an IBS0 proportion two orders of magnitude above this.
const PO_IBS0_CUTOFF: f64 = 0.0055;

// ---------------------------------------------------------------------------
// The tie order
// ---------------------------------------------------------------------------

/// The order the reference visits the members of an all-tied family of size *n*, as
/// 1-based ranks under the ID comparator.
///
/// Measured directly: a family of *n* mutually unrelated founders keeps all *n*, so the
/// emitted file *is* the visit order. Sizes 2…70 were swept against the reference with
/// every run verified to remove nobody. Two independent checks say the table is a
/// property of *n* alone — a dataset holding two 10-member families emits the same
/// permutation twice, and padding the dataset from 10 to 30 samples leaves a 10-member
/// family's permutation unchanged.
///
/// Indexed by `n - 2`; sizes above `TIE_ORDER.len() + 1` fall back to rank order.
const TIE_ORDER: &[&[u8]] = &[
    &[1, 2],
    &[1, 3, 2],
    &[1, 2, 4, 3],
    &[1, 2, 5, 4, 3],
    &[2, 3, 1, 6, 5, 4],
    &[2, 3, 1, 7, 6, 4, 5],
    &[3, 2, 4, 1, 8, 7, 5, 6],
    &[3, 2, 9, 4, 5, 1, 7, 8, 6],
    &[4, 3, 10, 9, 2, 5, 1, 8, 6, 7],
    &[4, 3, 11, 5, 2, 6, 1, 9, 10, 7, 8],
    &[4, 3, 5, 12, 11, 2, 6, 1, 10, 7, 8, 9],
    &[4, 3, 5, 13, 6, 2, 7, 11, 10, 1, 12, 8, 9],
    &[5, 4, 6, 14, 13, 3, 2, 7, 12, 11, 1, 8, 9, 10],
    &[5, 4, 6, 15, 7, 3, 2, 8, 13, 12, 1, 14, 9, 11, 10],
    &[5, 6, 4, 7, 16, 15, 3, 2, 8, 14, 13, 1, 9, 10, 12, 11],
    &[5, 6, 4, 7, 17, 8, 3, 2, 9, 14, 15, 13, 1, 16, 10, 12, 11],
    &[
        6, 7, 18, 5, 4, 8, 2, 17, 3, 9, 15, 16, 14, 1, 10, 11, 13, 12,
    ],
    &[
        6, 7, 19, 5, 4, 8, 2, 9, 3, 10, 16, 17, 18, 15, 14, 1, 12, 11, 13,
    ],
    &[
        6, 7, 20, 19, 8, 5, 9, 2, 4, 3, 10, 17, 18, 11, 16, 15, 1, 13, 12, 14,
    ],
    &[
        6, 7, 21, 10, 8, 5, 9, 2, 4, 3, 11, 17, 18, 20, 12, 19, 16, 1, 13, 15, 14,
    ],
    &[
        7, 8, 22, 6, 9, 5, 10, 2, 21, 4, 3, 11, 18, 19, 12, 13, 20, 17, 1, 14, 16, 15,
    ],
    &[
        7, 8, 23, 6, 9, 5, 10, 2, 11, 4, 3, 12, 19, 20, 22, 18, 21, 17, 1, 14, 13, 16, 15,
    ],
    &[
        8, 9, 7, 24, 23, 10, 6, 11, 2, 5, 4, 3, 12, 20, 21, 13, 19, 22, 18, 1, 15, 14, 17, 16,
    ],
    &[
        8, 9, 7, 25, 12, 10, 6, 11, 2, 5, 4, 3, 13, 21, 22, 20, 24, 14, 23, 19, 1, 15, 18, 17, 16,
    ],
    &[
        9, 10, 8, 26, 7, 11, 6, 2, 3, 12, 25, 5, 4, 13, 22, 23, 21, 14, 15, 24, 20, 1, 16, 19, 18,
        17,
    ],
    &[
        9, 10, 8, 27, 7, 11, 6, 2, 3, 12, 13, 5, 4, 14, 23, 24, 22, 26, 21, 25, 20, 16, 17, 1, 15,
        19, 18,
    ],
    &[
        9, 10, 8, 28, 27, 11, 12, 7, 2, 3, 13, 6, 5, 4, 14, 24, 25, 23, 15, 22, 26, 21, 17, 18, 1,
        16, 20, 19,
    ],
    &[
        9, 10, 8, 29, 14, 11, 12, 7, 2, 3, 13, 6, 5, 4, 15, 24, 25, 23, 28, 16, 26, 27, 22, 17, 18,
        1, 21, 20, 19,
    ],
    &[
        10, 11, 9, 30, 8, 12, 13, 7, 2, 3, 14, 29, 6, 4, 5, 15, 25, 26, 24, 16, 17, 27, 28, 23, 18,
        19, 1, 22, 21, 20,
    ],
];

/// The order in which a family's members are offered to the greedy.
///
/// Ascending relative count, and within one count the [`TIE_ORDER`] permutation of that
/// run. The second half is the honest weak point: the table was measured on families
/// where *every* member ties, and probe fixtures carrying two distinct counts show the
/// reference's scramble depends on the whole array rather than on each run in isolation
/// — a family of four isolated members followed by six paired ones emits the isolated
/// four as ranks 4, 3, 1, 2, which no per-run table produces. Applying the table per run
/// nevertheless reproduces every `*__unrelated` capture in the golden corpus except the
/// three merged clusters of `bigish`, where a two-member tie at ranks 4 and 10 of an
/// eleven-member cluster comes out in the opposite order.
fn visit_order(degree: &[usize]) -> Vec<usize> {
    let mut counts: Vec<usize> = degree.to_vec();
    counts.sort_unstable();
    counts.dedup();
    let mut out = Vec::with_capacity(degree.len());
    for d in counts {
        let run: Vec<usize> = (0..degree.len()).filter(|&i| degree[i] == d).collect();
        match TIE_ORDER.get(run.len().wrapping_sub(2)) {
            Some(perm) => out.extend(perm.iter().map(|&p| run[usize::from(p) - 1])),
            None => out.extend(run),
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Relatedness
// ---------------------------------------------------------------------------

/// Pedigree kinship over the whole sample set, resolving parents **within a family**.
///
/// A parent named by a row of another family is not that person: `multifam`'s `C_F` names
/// `A_F` and `A_M`, who live in `FAM1`, and the reference treats `C_F` as a founder of
/// `FAM3` whose parents merely happen to share their names. It materialises them as new,
/// ungenotyped members of `FAM3` — which is what makes it report the IDs as not unique.
fn pedigree_kinship(samples: &[Sample]) -> Vec<Vec<f64>> {
    let n = samples.len();
    let locate = |fid: &str, iid: &str| {
        samples
            .iter()
            .position(|s| s.fid == fid && s.iid == iid)
            .filter(|_| iid != "0")
    };
    let parents: Vec<(Option<usize>, Option<usize>)> = samples
        .iter()
        .map(|s| (locate(&s.fid, &s.pat), locate(&s.fid, &s.mat)))
        .collect();

    // Parents before children: a `.fam` may list them either way round.
    let mut order = Vec::with_capacity(n);
    let mut state = vec![0u8; n];
    let mut stack: Vec<usize> = Vec::new();
    for start in 0..n {
        if state[start] != 0 {
            continue;
        }
        stack.push(start);
        while let Some(&i) = stack.last() {
            if state[i] == 0 {
                state[i] = 1;
                for p in [parents[i].0, parents[i].1].into_iter().flatten() {
                    if state[p] == 0 {
                        stack.push(p);
                    }
                }
            } else {
                stack.pop();
                if state[i] == 1 {
                    state[i] = 2;
                    order.push(i);
                }
            }
        }
    }

    // `phi[i][j] = (phi[father(i)][j] + phi[mother(i)][j]) / 2`, and `phi[i][i]` carries
    // the inbreeding term. Walking `order` guarantees both parents are already filled in.
    let mut phi = vec![vec![0.0f64; n]; n];
    for &i in &order {
        let (fa, mo) = parents[i];
        let row: Vec<f64> = (0..n)
            .map(|j| {
                let a = fa.map_or(0.0, |p| phi[p][j]);
                let b = mo.map_or(0.0, |p| phi[p][j]);
                0.5 * (a + b)
            })
            .collect();
        for (j, &v) in row.iter().enumerate() {
            if j != i {
                phi[i][j] = v;
                phi[j][i] = v;
            }
        }
        phi[i][i] = 0.5 * (1.0 + fa.zip(mo).map_or(0.0, |(f, m)| phi[f][m]));
    }
    phi
}

/// The graph the greedy runs on: pedigree kinship, plus estimated kinship once the
/// genotypes are available.
struct Graph {
    n: usize,
    edge: Vec<bool>,
}

impl Graph {
    fn build(loaded: &Loaded, use_genotypes: bool) -> Self {
        let samples = &loaded.fileset.samples;
        let n = samples.len();
        let phi = pedigree_kinship(samples);
        let mut edge = vec![false; n * n];
        for i in 0..n {
            for j in i + 1..n {
                let related =
                    phi[i][j] > EDGE || (use_genotypes && Self::estimate(loaded, i, j) > EDGE);
                edge[i * n + j] = related;
                edge[j * n + i] = related;
            }
        }
        Graph { n, edge }
    }

    /// The estimator the reference would print for this pair, chosen by `FID` exactly as
    /// `.kin` and `.kin0` are.
    fn estimate(loaded: &Loaded, i: usize, j: usize) -> f64 {
        let samples = &loaded.fileset.samples;
        let scope = if samples[i].fid == samples[j].fid {
            Scope::WithinFamily
        } else {
            Scope::BetweenFamily
        };
        kinship::kinship(&counts::pair_counts(&loaded.fileset.genotypes, i, j), scope)
    }

    fn related(&self, i: usize, j: usize) -> bool {
        self.edge[i * self.n + j]
    }
}

// ---------------------------------------------------------------------------
// Clustering
// ---------------------------------------------------------------------------

/// One family as the selection sees it: the original FIDs it was built from and its
/// members, sorted by the ID comparator.
struct Cluster {
    /// The name the reference sorts families by: the FID for an unmerged family, `KING<k>`
    /// for a merged one.
    key: String,
    original: Vec<String>,
    members: Vec<usize>,
}

/// Group samples into clusters, merging families joined by an inferred 1st-degree pair.
///
/// Merging is gated on the sample count; below [`CLUSTERING_MIN_SAMPLES`] the reference
/// never looks at a cross-family pair, so every family stands alone.
fn clusters(loaded: &Loaded) -> Vec<Cluster> {
    let samples = &loaded.fileset.samples;
    let n = samples.len();
    let mut fids: Vec<String> = Vec::new();
    for s in samples {
        if !fids.iter().any(|f| f == &s.fid) {
            fids.push(s.fid.clone());
        }
    }
    let mut parent: Vec<usize> = (0..fids.len()).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }
    if n >= CLUSTERING_MIN_SAMPLES {
        let of = |fid: &str| fids.iter().position(|f| f == fid).unwrap_or(0);
        for i in 0..n {
            for j in i + 1..n {
                if samples[i].fid == samples[j].fid {
                    continue;
                }
                if Graph::estimate(loaded, i, j) > CLUSTER_EDGE {
                    let (a, b) = (
                        find(&mut parent, of(&samples[i].fid)),
                        find(&mut parent, of(&samples[j].fid)),
                    );
                    if a != b {
                        parent[a] = b;
                    }
                }
            }
        }
    }

    // Group in first-appearance order, so a merged cluster is numbered by the position of
    // its earliest family — `bigish` names BF01+BF02 `KING1` and BF25+BF26 `KING3`.
    let mut groups: Vec<(usize, Vec<String>)> = Vec::new();
    for (k, fid) in fids.iter().enumerate() {
        let root = find(&mut parent, k);
        match groups.iter_mut().find(|(r, _)| *r == root) {
            Some((_, list)) => list.push(fid.clone()),
            None => groups.push((root, vec![fid.clone()])),
        }
    }
    let mut merged = 0usize;
    let mut out: Vec<Cluster> = groups
        .into_iter()
        .map(|(_, original)| {
            let key = if original.len() > 1 {
                merged += 1;
                format!("KING{merged}")
            } else {
                original[0].clone()
            };
            let mut members: Vec<usize> = (0..n)
                .filter(|&i| original.iter().any(|f| *f == samples[i].fid))
                .collect();
            members.sort_by(|&a, &b| {
                king_id_cmp(samples[a].iid.as_bytes(), samples[b].iid.as_bytes())
            });
            Cluster {
                key,
                original,
                members,
            }
        })
        .collect();
    out.sort_by(|a, b| king_id_cmp(a.key.as_bytes(), b.key.as_bytes()));
    out
}

// ---------------------------------------------------------------------------
// The pass
// ---------------------------------------------------------------------------

/// The family-clustering prologue, shared verbatim by `--unrelated`, `--cluster` and
/// `--build`.
///
/// Everything from `Family clustering starts at` down to and including either
/// `No families were found to be connected.` or the `The following families are found to
/// be connected` table — byte for byte the same text in all three passes, which is why it
/// lives here rather than being written out three times.
///
/// Two shapes, decided by [`TINY_DATASET`]:
///
/// ```text
/// Family clustering starts at <t>              Family clustering starts at <t>
/// This function is currently disabled ...      Autosome genotypes stored in <w> words ...
///                                              Sorting autosomes...
///                                              <the allsegs.txt pre-pass>
///                                              <n> CPU cores are used to compute ...
///                                              Clustering up to <d> relatives ...
///                                              <the ID-uniqueness line>
///                                              No families were found to be connected.
/// ```
///
/// The block never ends with a blank line of its own: `--cluster` follows it straight
/// with `KING ends at`, and `--unrelated` supplies its own separator. Callers that need
/// the clusters get them back rather than recomputing.
pub fn clustering_prologue(opts: &Options, loaded: &Loaded, out: &mut dyn Write) -> Clustering {
    let samples = &loaded.fileset.samples;
    let tiny = samples.len() < TINY_DATASET;

    let _ = out.write_all(
        format!(
            "Family clustering starts at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );

    if tiny {
        let _ = out.write_all(
            b"This function is currently disabled for tiny dataset with sample size < 10.\n",
        );
    } else {
        // The preamble line without the blank line under it: `Sorting autosomes...`
        // follows immediately.
        let _ = out.write_all(
            console::autosome_words(loaded.words(), loaded.fileset.samples.len()).as_bytes(),
        );
        let _ = out.write_all(console::SORTING_AUTOSOMES.as_bytes());
        let _ = out.write_all(super::ibdseg::segment_prepass(opts, loaded).as_bytes());
        let _ = out.write_all(
            format!(
                "{} CPU cores are used to compute the pairwise kinship coefficients...\n",
                super::cpu_count(opts)
            )
            .as_bytes(),
        );
        let _ = out.write_all(
            format!(
                "Clustering up to {} relatives in families...\n",
                degree_label(opts.int(Opt::Degree))
            )
            .as_bytes(),
        );
        let _ = out.write_all(id_uniqueness(samples).as_bytes());
    }

    let graph = Graph::build(loaded, !tiny);
    let groups = clusters(loaded);
    let mut any_merged = false;

    if !tiny {
        let merged: Vec<&Cluster> = groups.iter().filter(|c| c.original.len() > 1).collect();
        any_merged = !merged.is_empty();
        if merged.is_empty() {
            let _ = out.write_all(b"No families were found to be connected.\n");
        } else {
            let _ = out.write_all(screening_summary(loaded, opts.int(Opt::Degree)).as_bytes());
            let _ = out.write_all(connected_families(&merged).as_bytes());
        }
    }
    Clustering {
        tiny,
        any_merged,
        graph,
        groups,
    }
}

/// What [`clustering_prologue`] worked out, for the pass that called it.
pub struct Clustering {
    /// Whether the dataset was too small for clustering to run at all.
    pub tiny: bool,
    /// Whether any two families were joined.
    pub any_merged: bool,
    graph: Graph,
    groups: Vec<Cluster>,
}

/// Run the `--unrelated` pass: console body plus the two lists.
pub fn run(opts: &Options, loaded: &Loaded, out: &mut dyn Write) {
    let samples = &loaded.fileset.samples;
    let Clustering { graph, groups, .. } = clustering_prologue(opts, loaded, out);
    // The prologue stops on its last line of text; the list block below opens with a
    // blank one, which is the second blank line after a `connected` table.
    let _ = out.write_all(b"\n");

    let (kept, removed) = select(&groups, &graph);
    let keep_path = out_path(opts, "unrelated.txt");
    let drop_path = out_path(opts, "unrelated_toberemoved.txt");
    write_list(&keep_path, samples, &kept);
    write_list(&drop_path, samples, &removed);

    let _ = out.write_all(
        format!(
            "A list of {} unrelated individuals saved in file {keep_path}\n",
            kept.len()
        )
        .as_bytes(),
    );
    let _ = out.write_all(
        format!(
            "An alternative list of {} to-be-removed individuals saved in file {drop_path}\n\n",
            removed.len()
        )
        .as_bytes(),
    );
    let _ = out.write_all(
        format!(
            "Extracting a subset of unrelated individuals ends at {}\n",
            console::ctime(console::now_local())
        )
        .as_bytes(),
    );
}

/// The greedy, cluster by cluster: kept in visit order, removed in ID order.
fn select(groups: &[Cluster], graph: &Graph) -> (Vec<usize>, Vec<usize>) {
    let mut kept = Vec::new();
    let mut removed = Vec::new();
    for c in groups {
        let m = &c.members;
        let degree: Vec<usize> = (0..m.len())
            .map(|i| {
                (0..m.len())
                    .filter(|&j| j != i && graph.related(m[i], m[j]))
                    .count()
            })
            .collect();
        let mut taken: Vec<usize> = Vec::new();
        for i in visit_order(&degree) {
            if taken.iter().all(|&t| !graph.related(m[i], m[t])) {
                taken.push(i);
                kept.push(m[i]);
            }
        }
        removed.extend((0..m.len()).filter(|i| !taken.contains(i)).map(|i| m[i]));
    }
    (kept, removed)
}

/// `FID\tIID\n` per line, no header. An empty list still writes an empty file.
fn write_list(path: &str, samples: &[Sample], rows: &[usize]) {
    let mut text = String::new();
    for &i in rows {
        text.push_str(&samples[i].fid);
        text.push('\t');
        text.push_str(&samples[i].iid);
        text.push('\n');
    }
    let _ = std::fs::write(path, text);
}

/// The line reporting whether IIDs identify a person on their own.
///
/// The check is over the pedigree the reference *builds*, not over the `.fam` rows: a row
/// naming a parent who is absent from its own family creates that parent there, so
/// `multifam`'s `FAM3 C_F A_F A_M` puts a second `A_F` beside `FAM1`'s and the reference
/// reports `Individual IDs are not unique (e.g., A_F)` on an otherwise duplicate-free
/// file. The example named is the first duplicate met walking the rows.
fn id_uniqueness(samples: &[Sample]) -> String {
    // Every person the pedigree names, in the order the reference meets them: each row,
    // then the parents that row creates inside its own family.
    let mut people: Vec<(&str, &str)> = Vec::new();
    for s in samples {
        people.push((&s.fid, &s.iid));
        for parent in [&s.pat, &s.mat] {
            if parent != "0" {
                people.push((&s.fid, parent));
            }
        }
    }
    let mut seen: Vec<(&str, &str)> = Vec::new();
    for (fid, iid) in people {
        if let Some((_, dup)) = seen.iter().find(|(f, i)| *i == iid && *f != fid) {
            return not_unique(dup);
        }
        if !seen.contains(&(fid, iid)) {
            seen.push((fid, iid));
        }
    }
    "Individual IDs are unique across all families.\n".to_string()
}

/// `--degree`'s echo in the clustering line: `1st`, `2nd`, `3rd`, `4th`, …
///
/// `--degree 0` is the reference's "unset", and the line then reads `1st-degree` — the
/// only thing `--degree` changes here, since the selection itself ignores it.
fn degree_label(degree: i32) -> String {
    let d = degree.max(1);
    let suffix = match (d % 10, d % 100) {
        (1, 11) | (2, 12) | (3, 13) => "th",
        (1, _) => "st",
        (2, _) => "nd",
        (3, _) => "rd",
        _ => "th",
    };
    format!("{d}{suffix}-degree")
}

fn not_unique(iid: &str) -> String {
    format!("  Individual IDs are not unique (e.g., {iid}), and family IDs will be used as well.\n")
}

/// The summary of what the cross-family screening found, printed only when it found
/// something.
///
/// A narrower table than `.kin`'s: one `Inference` row, no `Pedigree` row and no `OTHER`
/// column, because the screening only ever reports 1st-degree-or-closer pairs — `bigish`
/// shows three full-sibling pairs and zeros everywhere else even though it carries plenty
/// of more distant cross-family relatives.
fn screening_summary(loaded: &Loaded, degree: i32) -> String {
    let samples = &loaded.fileset.samples;
    let n = samples.len();
    // `--degree` widens what the screening *reports* without widening what it merges:
    // `bigish --unrelated --degree 2` adds 23 second-degree pairs to the table and still
    // connects exactly the three family pairs a bare `--unrelated` connects.
    let reported = band::MZ / f64::from(1u32 << degree.max(1));
    let mut cells = [0u64; 6];
    for i in 0..n {
        for j in i + 1..n {
            if samples[i].fid == samples[j].fid {
                continue;
            }
            let counts = counts::pair_counts(&loaded.fileset.genotypes, i, j);
            let phi = kinship::kinship(&counts, Scope::BetweenFamily);
            if phi < reported {
                continue;
            }
            let cell = if phi >= band::MZ {
                0
            } else if phi >= band::FIRST {
                usize::from(kinship::ibs0_prop(&counts) >= PO_IBS0_CUTOFF) + 1
            } else if phi >= band::SECOND {
                3
            } else if phi >= band::THIRD {
                4
            } else {
                5
            };
            cells[cell] += 1;
        }
    }
    let total: u64 = cells.iter().sum();
    let [mz, po, fs, second, third, fourth] = cells;
    format!(
        "\nRelationship summary (total relatives: 0 by pedigree, {total} by inference)\n\
         \x20       \tMZ\tPO\tFS\t2nd\t3rd\t4th\n\
         \x20 =========================================================\n\
         \x20 Inference\t{mz}\t{po}\t{fs}\t{second}\t{third}\t{fourth}\n\n"
    )
}

/// The table naming each merged cluster and the families it absorbed.
fn connected_families(merged: &[&Cluster]) -> String {
    let mut s = String::from("The following families are found to be connected\n");
    s.push_str("  NewFamID  OriginalFamID                                     \n");
    for c in merged {
        s.push_str(&format!("  {:<10}{}\n", c.key, c.original.join(",")));
    }
    // One blank line closes the table. `--cluster` puts its `Update-ID information …`
    // line straight after it and `--unrelated` adds a second blank of its own, which is
    // why the two passes show one and two blank lines here respectively.
    s.push('\n');
    s
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(fid: &str, iid: &str, pat: &str, mat: &str) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: pat.to_string(),
            mat: mat.to_string(),
            sex: 1,
            pheno: "-9".to_string(),
        }
    }

    #[test]
    fn tie_order_is_a_permutation_for_every_measured_size() {
        for (k, perm) in TIE_ORDER.iter().enumerate() {
            let n = k + 2;
            assert_eq!(perm.len(), n, "size {n}");
            let mut seen: Vec<u8> = perm.to_vec();
            seen.sort_unstable();
            assert_eq!(seen, (1..=n as u8).collect::<Vec<_>>(), "size {n}");
        }
    }

    /// The `unrelated` capture: ten mutually unrelated members of one FID, all kept.
    #[test]
    fn all_tied_family_uses_the_measured_order() {
        let order = visit_order(&[0; 10]);
        assert_eq!(order, [3, 2, 9, 8, 1, 4, 0, 7, 5, 6]);
    }

    /// `threegen`: the least-connected founder first, then the two spouses that married
    /// in, then the grandmother, then the grandfather.
    #[test]
    fn visit_order_is_ascending_relative_count() {
        // Counts for TG_C1..C4, TG_GF, TG_GM1, TG_GM2, TG_P1..P3, TG_S1, TG_S2 in ID order.
        let degree = [6, 6, 6, 6, 7, 6, 1, 8, 8, 4, 2, 2];
        let order = visit_order(&degree);
        assert_eq!(order[0], 6, "TG_GM2, the only member with one relative");
        assert_eq!(&order[1..3], [10, 11], "TG_S1 then TG_S2");
        assert_eq!(order[3], 9, "TG_P3, the only member with four");
    }

    #[test]
    fn a_family_with_one_member_needs_no_permutation() {
        assert_eq!(visit_order(&[0]), [0]);
        assert!(visit_order(&[]).is_empty());
    }

    #[test]
    fn pedigree_kinship_follows_the_declared_chain() {
        let samples = [
            sample("F", "A", "0", "0"),
            sample("F", "B", "A", "0"),
            sample("F", "C", "B", "0"),
            sample("F", "D", "C", "0"),
        ];
        let phi = pedigree_kinship(&samples);
        assert!((phi[0][1] - 0.25).abs() < 1e-12);
        assert!((phi[0][2] - 0.125).abs() < 1e-12);
        assert!((phi[0][3] - 0.0625).abs() < 1e-12);
        // All four are above the 4th-degree edge, so a chain of four keeps only its ends.
        assert!(phi[0][3] > EDGE);
    }

    #[test]
    fn a_parent_named_in_another_family_is_a_stranger() {
        // `multifam`'s shape: FAM3's C_F names FAM1's A_F and A_M.
        let samples = [
            sample("FAM1", "A_F", "0", "0"),
            sample("FAM1", "A_M", "0", "0"),
            sample("FAM3", "C_F", "A_F", "A_M"),
        ];
        let phi = pedigree_kinship(&samples);
        assert_eq!(phi[0][2], 0.0);
        assert_eq!(
            id_uniqueness(&samples),
            "  Individual IDs are not unique (e.g., A_F), and family IDs will be used as well.\n"
        );
    }

    #[test]
    fn unique_ids_report_as_unique() {
        let samples = [
            sample("F1", "A", "0", "0"),
            sample("F1", "B", "0", "0"),
            sample("F1", "C", "A", "B"),
            sample("F2", "D", "0", "0"),
        ];
        assert_eq!(
            id_uniqueness(&samples),
            "Individual IDs are unique across all families.\n"
        );
    }
}
