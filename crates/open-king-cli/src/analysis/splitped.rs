//! `<prefix>splitped.txt` — the pedigree-splitting artefact `--ibdseg` leaves behind.
//!
//! # What it is
//!
//! The reference announces it on an `--ibdseg` run only when [`is_generated`] is true
//! (`<prefix>splitped.txt is generated for certain pedigree plot applications.`) and
//! writes it before any segment work happens. It is a *pedigree* file, not a segment
//! file: nothing in it depends on genotypes, on `--degree`, or on `--seglength`, and the
//! captured corpus confirms that — the file is byte-identical across all five `--ibdseg`
//! captures of every dataset.
//!
//! # Format (verified)
//!
//! Nine **space**-separated fields, no header, `\n` terminated:
//!
//! ```text
//! OldFID OldIID NewFID NewIID Father Mother Sex Pheno Dummy
//! ```
//!
//! `Pheno` is the `.fam` phenotype with `-9` rewritten to `0`. `Dummy` is `1` for a
//! person the reference invented and `0` for a genotyped one.
//!
//! # The rules
//!
//! Established black-box, by running the reference on the ten corpus datasets and on four
//! adversarial probe pedigrees (`tests/parity/probes/splitped.py` is the reference
//! implementation this module was ported from; it reproduces all fourteen byte for byte).
//!
//! 1. **A family of one parentless founder is dropped entirely.** Any other family is
//!    emitted, including a singleton that names a parent.
//! 2. **Absent parents are materialised.** A `PAT`/`MAT` naming somebody the family does
//!    not contain becomes a dummy founder row inside that family — taking its `SEX` from
//!    that ID's own `.fam` row if one exists anywhere in the file (a cross-family
//!    import), otherwise from the slot it fills.
//! 3. **A half-specified parentage is completed with an invented mate.** A person with
//!    exactly one of `PAT`/`MAT` set gets a `KING<n>` founder for the empty slot, and the
//!    person's own row is rewritten to point at it. `<n>` is a **global** counter running
//!    in `.fam` family order — not in output order, which is why it has to be assigned in
//!    a first pass.
//! 4. **Families are written in [`king_id_cmp`] order of the FID**, not `.fam` order.
//! 5. **Members are ordered by (generation depth, [`king_id_cmp`] of the IID)**, dummies
//!    included, where depth is 0 for a founder and `1 + max(parent depth)` otherwise.
//! 6. **A family that turns out to be several disconnected pedigrees is split**: each
//!    connected component becomes `<FID>_S<k>`, `k` counting from 1 in the order the
//!    components' first members appear in rule 5's order, and each component's members
//!    are then listed in **breadth-first** order from that first member. A family with a
//!    single component keeps its FID *and* rule 5's order — one order for both cases is
//!    wrong, which is the subtlety this module exists to record.

use std::collections::{HashMap, HashSet, VecDeque};

use open_king_io::Sample;

use crate::analysis::king_id_cmp;

/// Whether KING writes and announces `splitped.txt` for this pedigree.
///
/// A lone parentless founder contributes nothing. The file exists iff at least one family
/// has two members or a singleton names a parent; this is the same rule [`text`] applies
/// while deciding which families to retain.
pub fn is_generated(samples: &[Sample]) -> bool {
    let mut family_sizes: HashMap<&str, usize> = HashMap::new();
    for sample in samples {
        *family_sizes.entry(sample.fid.as_str()).or_default() += 1;
    }
    samples.iter().any(|sample| {
        family_sizes.get(sample.fid.as_str()).copied().unwrap_or(0) >= 2
            || sample.pat != "0"
            || sample.mat != "0"
    })
}

/// One row of the file, before ordering.
#[derive(Clone)]
struct Person {
    iid: String,
    fa: String,
    mo: String,
    sex: String,
    pheno: String,
    /// Rule 2/3: invented by the reference rather than read from the `.fam`.
    dummy: bool,
}

impl Person {
    fn from_sample(s: &Sample) -> Person {
        Person {
            iid: s.iid.clone(),
            fa: s.pat.clone(),
            mo: s.mat.clone(),
            sex: s.sex.to_string(),
            // `-9` — PLINK's "missing phenotype" — is written as `0`.
            pheno: if s.pheno == "-9" {
                "0".to_string()
            } else {
                s.pheno.clone()
            },
            dummy: false,
        }
    }

    fn founder(iid: String, sex: String) -> Person {
        Person {
            iid,
            fa: "0".to_string(),
            mo: "0".to_string(),
            sex,
            pheno: "0".to_string(),
            dummy: true,
        }
    }

    fn is_parentless(&self) -> bool {
        self.fa == "0" && self.mo == "0"
    }
}

/// Render the whole file. Empty when every family is a lone parentless founder.
pub fn text(samples: &[Sample]) -> String {
    // Families in `.fam` first-appearance order — the order rule 3's counter runs in.
    let mut order: Vec<&str> = Vec::new();
    let mut members: HashMap<&str, Vec<Person>> = HashMap::new();
    for s in samples {
        let e = members.entry(s.fid.as_str()).or_insert_with(|| {
            order.push(s.fid.as_str());
            Vec::new()
        });
        e.push(Person::from_sample(s));
    }
    // Rule 2's sex lookup is global and last-writer-wins, matching a plain `.fam` scan.
    let mut sex_of: HashMap<&str, String> = HashMap::new();
    for s in samples {
        sex_of.insert(s.iid.as_str(), s.sex.to_string());
    }

    let mut kept: Vec<&str> = Vec::new();
    let mut counter = 0usize;
    for &fid in &order {
        let ms = members.get_mut(fid).expect("family present");
        // Rule 1.
        if ms.len() < 2 && ms.iter().all(Person::is_parentless) {
            continue;
        }
        let mut have: HashSet<String> = ms.iter().map(|m| m.iid.clone()).collect();
        let mut invented: Vec<Person> = Vec::new();
        for m in ms.iter_mut() {
            // Rule 2, father slot then mother slot.
            for (named, slot_sex) in [(m.fa.clone(), "1"), (m.mo.clone(), "2")] {
                if named == "0" || have.contains(&named) {
                    continue;
                }
                let sex = sex_of
                    .get(named.as_str())
                    .cloned()
                    .unwrap_or_else(|| slot_sex.to_string());
                have.insert(named.clone());
                invented.push(Person::founder(named, sex));
            }
            // Rule 3 — tested on the *original* parent fields, so a person missing both
            // parents is left alone and one missing exactly one gains a mate.
            let (fa_missing, mo_missing) = (m.fa == "0", m.mo == "0");
            if fa_missing != mo_missing {
                counter += 1;
                let name = format!("KING{counter}");
                let sex = if fa_missing { "1" } else { "2" }.to_string();
                if fa_missing {
                    m.fa = name.clone();
                } else {
                    m.mo = name.clone();
                }
                have.insert(name.clone());
                invented.push(Person::founder(name, sex));
            }
        }
        // Dummies lead; rule 5 reorders everything a moment later anyway, but the
        // starting order decides ties the comparator calls equal.
        invented.extend(ms.iter().cloned());
        *ms = invented;
        kept.push(fid);
    }

    // Rule 4.
    kept.sort_by(|a, b| king_id_cmp(a.as_bytes(), b.as_bytes()));

    let mut out = String::new();
    for fid in kept {
        let people = order_family(members.get(fid).expect("family present"));
        let comps = components(&people);
        let split = comps.len() > 1;
        for (ci, comp) in comps.iter().enumerate() {
            let new_fid = if split {
                format!("{fid}_S{}", ci + 1)
            } else {
                fid.to_string()
            };
            for &p in comp {
                let p = &people[p];
                out.push_str(&format!(
                    "{} {} {} {} {} {} {} {} {}\n",
                    fid,
                    p.iid,
                    new_fid,
                    p.iid,
                    p.fa,
                    p.mo,
                    p.sex,
                    p.pheno,
                    u8::from(p.dummy),
                ));
            }
        }
    }
    out
}

/// Rule 5: `(generation depth, id)`.
fn order_family(people: &[Person]) -> Vec<Person> {
    let index: HashMap<&str, usize> = people
        .iter()
        .enumerate()
        .map(|(i, p)| (p.iid.as_str(), i))
        .collect();
    let mut depth = vec![usize::MAX; people.len()];
    for i in 0..people.len() {
        generation(people, &index, &mut depth, i);
    }
    let mut ord: Vec<usize> = (0..people.len()).collect();
    ord.sort_by(|&a, &b| {
        depth[a]
            .cmp(&depth[b])
            .then_with(|| king_id_cmp(people[a].iid.as_bytes(), people[b].iid.as_bytes()))
    });
    ord.into_iter().map(|i| people[i].clone()).collect()
}

/// Depth of person `i`: 0 for a founder, else one past its deepest parent.
///
/// Iterative so a pedigree loop cannot blow the stack; a person reached while its own
/// depth is still being computed contributes nothing, which makes a cycle terminate at
/// the depth of whatever else it hangs from.
fn generation(
    people: &[Person],
    index: &HashMap<&str, usize>,
    depth: &mut [usize],
    start: usize,
) -> usize {
    const IN_PROGRESS: usize = usize::MAX - 1;
    let mut stack = vec![start];
    while let Some(&i) = stack.last() {
        if depth[i] != usize::MAX && depth[i] != IN_PROGRESS {
            stack.pop();
            continue;
        }
        let parents: Vec<usize> = [&people[i].fa, &people[i].mo]
            .into_iter()
            .filter(|p| *p != "0")
            .filter_map(|p| index.get(p.as_str()).copied())
            .collect();
        let pending: Vec<usize> = parents
            .iter()
            .copied()
            .filter(|&p| depth[p] == usize::MAX)
            .collect();
        if pending.is_empty() {
            depth[i] = parents
                .iter()
                .filter(|&&p| depth[p] != IN_PROGRESS)
                .map(|&p| depth[p] + 1)
                .max()
                .unwrap_or(0);
            stack.pop();
        } else {
            depth[i] = IN_PROGRESS;
            stack.extend(pending);
        }
    }
    depth[start]
}

/// Rule 6: connected components over parent–child edges, each in breadth-first order.
///
/// Seeds and adjacency lists both follow the rule-5 order, so the traversal is
/// deterministic. A family with one component is returned in rule-5 order untouched —
/// breadth-first from its first member is a *different* order and is not what the
/// reference writes.
fn components(people: &[Person]) -> Vec<Vec<usize>> {
    let index: HashMap<&str, usize> = people
        .iter()
        .enumerate()
        .map(|(i, p)| (p.iid.as_str(), i))
        .collect();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); people.len()];
    for (i, p) in people.iter().enumerate() {
        for parent in [&p.fa, &p.mo] {
            if parent == "0" {
                continue;
            }
            if let Some(&j) = index.get(parent.as_str()) {
                adj[j].push(i);
                adj[i].push(j);
            }
        }
    }
    for a in &mut adj {
        a.sort_unstable();
    }

    let mut seen = vec![false; people.len()];
    let mut comps: Vec<Vec<usize>> = Vec::new();
    for i in 0..people.len() {
        if seen[i] {
            continue;
        }
        seen[i] = true;
        let mut queue = VecDeque::from([i]);
        let mut comp = Vec::new();
        while let Some(cur) = queue.pop_front() {
            comp.push(cur);
            for &nb in &adj[cur] {
                if !seen[nb] {
                    seen[nb] = true;
                    queue.push_back(nb);
                }
            }
        }
        comps.push(comp);
    }
    if comps.len() == 1 {
        return vec![(0..people.len()).collect()];
    }
    comps
}

#[cfg(test)]
mod tests {
    use super::*;

    fn sample(fid: &str, iid: &str, pat: &str, mat: &str, sex: u8) -> Sample {
        Sample {
            fid: fid.to_string(),
            iid: iid.to_string(),
            pat: pat.to_string(),
            mat: mat.to_string(),
            sex,
            pheno: "-9".to_string(),
        }
    }

    #[test]
    fn a_lone_parentless_founder_is_dropped() {
        assert_eq!(text(&[sample("F", "A", "0", "0", 1)]), "");
    }

    #[test]
    fn generation_starts_at_family_size_two_or_a_named_parent() {
        let singleton_families = [
            sample("F1", "A", "0", "0", 1),
            sample("F2", "B", "0", "0", 2),
            sample("F3", "C", "0", "0", 1),
        ];
        assert!(!is_generated(&singleton_families));

        let mut pair = singleton_families.clone();
        pair[1].fid = "F1".to_string();
        assert!(is_generated(&pair));

        let named_parent = [sample("F1", "A", "DAD", "0", 1)];
        assert!(is_generated(&named_parent));
    }

    #[test]
    fn a_trio_keeps_its_family_and_orders_parents_first() {
        let fam = [
            sample("F", "KID", "DAD", "MOM", 1),
            sample("F", "DAD", "0", "0", 1),
            sample("F", "MOM", "0", "0", 2),
        ];
        assert_eq!(
            text(&fam),
            "F DAD F DAD 0 0 1 0 0\nF MOM F MOM 0 0 2 0 0\nF KID F KID DAD MOM 1 0 0\n"
        );
    }

    #[test]
    fn a_half_specified_parentage_invents_a_mate() {
        let fam = [
            sample("F", "KID", "DAD", "0", 1),
            sample("F", "DAD", "0", "0", 1),
        ];
        // KING1 is the invented mother; it is a founder, so it sorts with DAD.
        assert_eq!(
            text(&fam),
            "F DAD F DAD 0 0 1 0 0\nF KING1 F KING1 0 0 2 0 1\nF KID F KID DAD KING1 1 0 0\n"
        );
    }

    #[test]
    fn an_absent_parent_is_materialised_with_its_own_sex() {
        // GRANDPA is genotyped in another family, so his sex comes from there.
        let fam = [
            sample("G", "GRANDPA", "0", "0", 1),
            sample("G", "GRANDMA", "0", "0", 2),
            sample("F", "KID", "GRANDPA", "GRANDMA", 2),
            sample("F", "SIB", "GRANDPA", "GRANDMA", 1),
        ];
        let out = text(&fam);
        assert!(out.contains("F GRANDPA F GRANDPA 0 0 1 0 1\n"), "{out}");
        assert!(out.contains("F GRANDMA F GRANDMA 0 0 2 0 1\n"), "{out}");
    }

    #[test]
    fn two_disconnected_pedigrees_in_one_family_are_split() {
        let fam = [
            sample("F", "A1", "0", "0", 1),
            sample("F", "A2", "0", "0", 2),
            sample("F", "A3", "A1", "A2", 1),
            sample("F", "B1", "0", "0", 1),
            sample("F", "B2", "0", "0", 2),
            sample("F", "B3", "B1", "B2", 1),
        ];
        let out = text(&fam);
        assert!(out.contains("F A1 F_S1 A1 "), "{out}");
        assert!(out.contains("F B1 F_S2 B1 "), "{out}");
    }
}
