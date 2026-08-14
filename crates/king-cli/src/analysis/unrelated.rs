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
//! * members with equal counts come out in an order that is *not* a tie-break on rank: it
//!   is the residue of an unstable sort, and it depends on the whole count array rather
//!   than on each tied run in isolation.
//!
//! [`sort_by_count_descending`] is that sort, identified rather than tabulated. See its
//! documentation for the algorithm and for the experiments that fixed every parameter.
//!
//! # How to interrogate the order directly
//!
//! Two facts make the reference's visit order fully observable, which is worth knowing
//! before anyone re-derives any of the above by hand:
//!
//! * a family of mutually unrelated members keeps everybody, so the emitted list **is**
//!   the visit order — build one family per size and read the permutation straight off;
//! * a run of members that all share one count can be totally ordered even when the
//!   greedy drops most of them: realise the run as a perfect matching (each member
//!   related to exactly one other, so every member's count is 1 whatever the pairing is),
//!   and the count array — hence the visit order — is the same for every pairing. Pair
//!   *a* with *b* and the survivor of that pair names which of the two is visited first.
//!
//! Both work on filesets whose relatedness comes from the `.fam` alone, so the genotypes
//! can be independent random draws: a group of *m* declared full sibs is a clique of *m*,
//! and a lone founder is an isolated vertex, which is enough to realise any count array
//! that is a disjoint union of cliques.

use std::io::Write;

use king_core::{counts, kinship, Scope};
use king_io::Sample;

use crate::analysis::{band, king_id_cmp, out_path, with_phantom_parents};
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

// ---------------------------------------------------------------------------
// The visit order
// ---------------------------------------------------------------------------

/// Subfile length above which the reference's sort takes its median-of-three pivot.
///
/// The condition is `r - l > 7` on the inclusive range, i.e. median-of-three from nine
/// elements up and the plain last-element partition at eight and below. Both halves are
/// forced: a family of *n* mutually unrelated members keeps everybody, so its emitted
/// list *is* the visit order, and replaying sizes 2…70 against the reference singles this
/// cut-point out — `> 6` already disagrees at *n* = 8, `> 8` at *n* = 9, and every other
/// value in 3…11 disagrees somewhere in 2…30.
const MEDIAN_OF_THREE_CUTOFF: isize = 7;

/// The order in which a family's members are offered to the greedy.
///
/// The reference sorts the members by relative count **descending** and then walks the
/// sorted array **backwards**, which is why the visit order is by ascending count and why
/// its tie residue is the residue of an unstable sort rather than a tie-break on rank.
/// [`sort_by_count_descending`] is that sort.
fn visit_order(degree: &[usize]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..degree.len()).collect();
    sort_by_count_descending(&mut order, degree);
    order.reverse();
    order
}

/// `a` precedes `b` under the reference's ordering: a *larger* count sorts first.
fn precedes(degree: &[usize], a: usize, b: usize) -> bool {
    degree[a] > degree[b]
}

/// Sedgewick's median-of-three quicksort, which is the sort the reference uses.
///
/// ```text
/// quicksort(l, r):
///     if r <= l: return
///     if r - l > 7:                       # median of three, into the r-1 slot
///         exch(a[(l+r)/2], a[r-1])
///         compexch(a[l], a[r-1]); compexch(a[l], a[r]); compexch(a[r-1], a[r])
///         (pl, pr) = (l+1, r-1)           # a[l] and a[r] are now sentinels
///     else:
///         (pl, pr) = (l, r)               # plain last-element pivot
///     v = a[pr]; i = pl-1; j = pr
///     loop:
///         while a[++i] < v: ;
///         while v < a[--j]: if j == pl: break
///         if i >= j: break
///         exch(a[i], a[j])
///     exch(a[i], a[pr])
///     quicksort(l, i-1); quicksort(i+1, r)
/// ```
///
/// # How it was established
///
/// Nothing here is a table. Equal keys make every `<` false, so the algorithm collapses to
/// one fixed permutation per size, and a family of mutually unrelated members exposes that
/// permutation directly. Working backwards from the measured permutations for *n* = 2…30:
///
/// * assuming a quicksort and inverting the recursion recovers, for every *n* from 9 to
///   20, a partition step that is **exactly** this one — twelve independent *n*-element
///   agreements, and in each case it is also the fewest-swaps partition consistent with
///   the data;
/// * the same inversion fails at *n* ≤ 8, which is what says there is a cut-off; searching
///   17 280 parameterised two-pointer partitions for one matching all seven measured
///   permutations at *n* = 2…8 returns the plain last-element partition above, uniquely;
/// * the two branches together reproduce all 29 measured permutations, and then
///   **extrapolate**: *n* = 31, 35, 40, 55 and 70 were predicted before being measured and
///   all five matched. (The table this replaced stopped at 31 and fell back to rank order.)
///
/// # What the real-data probes fix
///
/// Equal keys cannot see the median-of-three or the `j == pl` guard, so those were checked
/// against filesets whose counts differ. Two families of *a* and *b* children joined by a
/// full-sib link between the fathers give the count array
/// `[a+b+2]*a, a+b+1, a, [a+b+2]*b, a+b+1, b` — `bigish`'s merged clusters are the
/// `a = 3, b = 4` case. All 49 shapes for *a*, *b* ∈ 1…7 were run against the reference:
/// this function reproduces every one, and dropping the median-of-three breaks 30 of them.
/// A second, structurally unrelated sweep — 45 random families of up to 26 members built
/// as random disjoint cliques, so the count multiset and its rank layout are arbitrary —
/// also matches on all 45.
///
/// The `j == pl` guard never fires in any of those 94 probes; it is Sedgewick's own bound
/// on the descending scan and is kept because without it the scan can run off the array on
/// input the corpus does not contain.
fn sort_by_count_descending(order: &mut [usize], degree: &[usize]) {
    fn rec(a: &mut [usize], degree: &[usize], l: isize, r: isize) {
        if r <= l {
            return;
        }
        let (pl, pr) = if r - l > MEDIAN_OF_THREE_CUTOFF {
            let (lo, hi) = (l as usize, r as usize);
            a.swap((lo + hi) / 2, hi - 1);
            if precedes(degree, a[hi - 1], a[lo]) {
                a.swap(lo, hi - 1);
            }
            if precedes(degree, a[hi], a[lo]) {
                a.swap(lo, hi);
            }
            if precedes(degree, a[hi], a[hi - 1]) {
                a.swap(hi - 1, hi);
            }
            (l + 1, r - 1)
        } else {
            (l, r)
        };
        let v = a[pr as usize];
        let (mut i, mut j) = (pl - 1, pr);
        loop {
            i += 1;
            while precedes(degree, a[i as usize], v) {
                i += 1;
            }
            j -= 1;
            while precedes(degree, v, a[j as usize]) {
                if j == pl {
                    break;
                }
                j -= 1;
            }
            if i >= j {
                break;
            }
            a.swap(i as usize, j as usize);
        }
        a.swap(i as usize, pr as usize);
        rec(a, degree, l, i - 1);
        rec(a, degree, i + 1, r);
    }

    let last = order.len() as isize - 1;
    rec(order, degree, 0, last);
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
///
/// Those materialised parents are **relatives**, not just names, which is why the matrix
/// is built over [`with_phantom_parents`] and then cut back to the genotyped rows: two
/// rows of one family naming the same absent parents are declared full sibs, and the
/// reference drops one of them. No corpus fileset has such a pair — the rule was measured
/// on a two-row fixture run against the reference, where it keeps one of the two.
fn pedigree_kinship(samples: &[Sample]) -> Vec<Vec<f64>> {
    let genotyped = samples.len();
    let samples = &with_phantom_parents(samples)[..];
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
    phi.truncate(genotyped);
    for row in &mut phi {
        row.truncate(genotyped);
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
            let _ =
                out.write_all(screening_summary(opts, loaded, opts.int(Opt::Degree)).as_bytes());
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

impl Clustering {
    /// The clusters that absorbed more than one family, in the order the console table
    /// lists them: `(KING<k>, its members in ID order)`.
    ///
    /// The unmerged families are not here. `--cluster` and `--build` only ever rename,
    /// re-analyse and reconstruct the *newly clustered* ones.
    pub fn merged(&self) -> impl Iterator<Item = (&str, &[usize])> {
        self.groups
            .iter()
            .filter(|c| c.original.len() > 1)
            .map(|c| (c.key.as_str(), c.members.as_slice()))
    }

    /// `<prefix>updateids.txt` — `OLDFID OLDIID NEWFID NEWIID`, tab separated, no header.
    ///
    /// Only the merged clusters get a row; the IID never changes, so the fourth column
    /// repeats the second and the file is really the FID rename map. Rows follow the
    /// cluster order and, inside a cluster, the ID comparator — which for `bigish` is
    /// also original-family order, the two being indistinguishable there.
    pub fn updateids_text(&self, samples: &[Sample]) -> String {
        let mut s = String::new();
        for (key, members) in self.merged() {
            for &i in members {
                let (fid, iid) = (&samples[i].fid, &samples[i].iid);
                s.push_str(&format!("{fid}\t{iid}\t{key}\t{iid}\n"));
            }
        }
        s
    }
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
/// column. It is the **same screening `--related` runs**, so the two rules that pass are
/// taken from there rather than re-derived:
///
/// * **which pairs are reported** is a disjunction — kinship at or above `2^-(d+1.5)`
///   *or* `PropIBD` above `2^-(d+0.5)`. `bigish --unrelated --degree 2` is exactly the
///   case that separates the two halves: the reference counts 23 second-degree pairs
///   where the kinship cut alone finds 21, the two extra ones sitting at 0.0882 and
///   0.0870, just under `2^-3.5`, and reaching the band on their segment sharing;
/// * **which column a pair lands in** is its segment `InfType`, not its kinship band.
///   Those same two pairs go in `2nd`, which their kinship (below `2^-3.5`) would not do,
///   and the `4th` column stays zero throughout because `4th` and `UN` are `OTHER` here.
///
/// `bigish` is the only corpus fileset whose families merge, so it is the only one that
/// prints this table at all; both of its captures — bare and `--degree 2` — are
/// byte-identical under these rules.
fn screening_summary(opts: &Options, loaded: &Loaded, degree: i32) -> String {
    let samples = &loaded.fileset.samples;
    let n = samples.len();
    let d = f64::from(degree.max(1));
    let kin_cut = 2f64.powf(-(d + 1.5));
    let prop_cut = 2f64.powf(-(d + 0.5));
    let engine = super::related::Engine::autosomes(
        &loaded.fileset.variants,
        &loaded.fileset.kept,
        i64::from(opts.int(Opt::Sexchr)),
        super::ibdseg::seglength_bp(opts),
    );
    let genotypes = &loaded.fileset.genotypes;

    let mut cells = [0u64; 6];
    for i in 0..n {
        for j in i + 1..n {
            if samples[i].fid == samples[j].fid {
                continue;
            }
            let counts = counts::pair_counts(genotypes, i, j);
            if counts.n_snp == 0 {
                continue;
            }
            let phi = kinship::kinship(&counts, Scope::BetweenFamily);
            let ibd = engine.pair(genotypes, i, j);
            if !(phi >= kin_cut || ibd.prop_ibd > prop_cut) {
                continue;
            }
            let cell = match ibd.inf_type(kinship::het_concordance(&counts)) {
                "Dup/MZ" => 0,
                "PO" => 1,
                "FS" => 2,
                "2nd" => 3,
                "3rd" => 4,
                // `4th` and `UN` are the table's `OTHER`: neither printed nor totalled.
                _ => continue,
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
    fn the_visit_order_is_always_a_permutation() {
        for n in 0..80usize {
            for degree in [vec![0; n], (0..n).map(|i| i % 5).collect()] {
                let mut seen = visit_order(&degree);
                assert_eq!(seen.len(), n, "size {n}");
                seen.sort_unstable();
                assert_eq!(seen, (0..n).collect::<Vec<_>>(), "size {n}");
            }
        }
    }

    /// The order really is by ascending count; only its tie residue is subtle.
    #[test]
    fn the_visit_order_is_sorted_by_count() {
        let degree: Vec<usize> = (0..40).map(|i| (i * 7) % 6).collect();
        let counts: Vec<usize> = visit_order(&degree).iter().map(|&i| degree[i]).collect();
        assert!(counts.windows(2).all(|w| w[0] <= w[1]), "{counts:?}");
    }

    /// The `unrelated` capture: ten mutually unrelated members of one FID, all kept, so
    /// the emitted file *is* the visit order. Measured against the reference.
    #[test]
    fn all_tied_family_matches_the_reference() {
        // 1-based ranks 4 3 10 9 2 5 1 8 6 7, as the reference emits them.
        assert_eq!(visit_order(&[0; 10]), [3, 2, 9, 8, 1, 4, 0, 7, 5, 6]);
        // Two more sizes from the same sweep, one either side of the median-of-three
        // cut-off at nine.
        assert_eq!(visit_order(&[0; 8]), [2, 1, 3, 0, 7, 6, 4, 5]);
        assert_eq!(visit_order(&[0; 9]), [2, 1, 8, 3, 4, 0, 6, 7, 5]);
    }

    /// `bigish`'s merged clusters: eleven members, and the count-8 tie at ranks 4 and 10
    /// comes out 10 first — the case a per-run tie-break table gets backwards.
    #[test]
    fn a_merged_cluster_visits_the_later_of_the_tied_pair_first() {
        let degree = [9, 9, 9, 8, 3, 9, 9, 9, 9, 8, 4];
        let order = visit_order(&degree);
        assert_eq!(order[0], 4, "rank 5, the only member with three relatives");
        assert_eq!(order[1], 10, "rank 11, the only member with four");
        assert_eq!(order[2], 9, "rank 10 before rank 4 in the count-8 tie");
        assert_eq!(order[3], 3);
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
        assert_eq!(phi.len(), 3, "the matrix is cut back to the genotyped rows");
        assert_eq!(phi[0][2], 0.0);
        assert_eq!(
            id_uniqueness(&samples),
            "  Individual IDs are not unique (e.g., A_F), and family IDs will be used as well.\n"
        );
    }

    /// A parent the `.fam` never lists is still a person: two rows of one family naming
    /// the same absent parents are full sibs, and the reference keeps only one of them.
    #[test]
    fn rows_sharing_an_absent_parent_are_relatives() {
        let samples = [
            sample("F", "M01", "PP", "QQ"),
            sample("F", "M02", "PP", "QQ"),
            sample("F", "M03", "0", "0"),
        ];
        let phi = pedigree_kinship(&samples);
        assert!((phi[0][1] - 0.25).abs() < 1e-12, "declared full sibs");
        assert!(phi[0][1] > EDGE);
        assert_eq!(phi[0][2], 0.0);
    }

    /// `bigish`'s `KING1`: two families renamed, the IID carried through unchanged.
    #[test]
    fn updateids_names_only_the_merged_clusters() {
        let samples = [
            sample("BF01", "B01_F", "0", "0"),
            sample("BF02", "B02_F", "0", "0"),
            sample("BF03", "B03_F", "0", "0"),
        ];
        let clustering = Clustering {
            tiny: false,
            any_merged: true,
            graph: Graph {
                n: 3,
                edge: vec![false; 9],
            },
            groups: vec![
                Cluster {
                    key: "KING1".to_string(),
                    original: vec!["BF01".to_string(), "BF02".to_string()],
                    members: vec![0, 1],
                },
                Cluster {
                    key: "BF03".to_string(),
                    original: vec!["BF03".to_string()],
                    members: vec![2],
                },
            ],
        };
        assert_eq!(
            clustering.updateids_text(&samples),
            "BF01\tB01_F\tKING1\tB01_F\nBF02\tB02_F\tKING1\tB02_F\n"
        );
        assert_eq!(clustering.merged().count(), 1);
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
