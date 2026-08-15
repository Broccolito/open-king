#!/usr/bin/env python3
"""`Duplicate <a> (of <b>) is removed.` — which copy goes, and when the line is raised.

Two questions, both opened by one shape (`clusternum.py`'s `mixed_cluster`) and neither
answerable from it:

1. **Is the line rule-half or inference-half?**  It matters because this binary writes the
   rule half only.  The companion question for `Reconstruct parent-offspring pair` came
   back *inference-half* — 42 of 42 clusters that print it also print an `INFERENCE` line,
   and a `PO` merge between two families with no sibship anywhere prints nothing at all —
   so the same test is run here.
2. **Which of the two copies is removed?**  `mixed_cluster` removes the later one under the
   ID comparator; the singleton shapes below remove the *earlier*.  The candidate rule is
   that the reference keeps the copy with more **declared 1st-degree relatives present in
   the fileset** (parents plus full sibs the `.fam` names), ties going to the later id.

    python3 dupkeep.py            # the scorecard over every shape and seed

Nothing here reads KING's source.  Filesets land in `work/dupkeep/` (gitignored).
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
WORK = os.path.join(_HERE, "work", "dupkeep")
KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
)
DUP = re.compile(r"Duplicate (\S+) \(of (\S+)\) is removed\.")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


BS = _load("bs", os.path.join(_HERE, "build_shapes.py"))
G = BS.G


def _fs_singletons(ped, tags, ph):
    for k, t in enumerate(tags):
        ped.add(t, "F" + t, father=ph[0], mother=ph[1], sex=1 + (k % 2))


# --- shapes -----------------------------------------------------------------
# Each returns (pedigree, expected-tie-state) and puts exactly one duplicate pair in a
# cluster that also carries an FS edge, so the cluster raises rule lines either way.

def solo(clone_of="Q1", dup_id="Q4"):
    """Four singleton families: three undeclared full sibs and a duplicate of one.

    Neither copy declares a parent or a sib, so the connectivity statistic ties.
    """
    def build():
        ped = G.Ped()
        ph = G.add_couple(ped, "PHs", "PHs", emit=False)
        _fs_singletons(ped, ["Q1", "Q2", "Q3"], ph)
        ped.add(dup_id, "F" + dup_id, clone_of=clone_of, sex=1)
        return ped
    return build


def parented(first_has_parents):
    """One copy declares a couple and a sib, the other declares nothing.

    `first_has_parents` decides whether the connected copy sorts first or last, which is
    what separates "keep the connected one" from "keep the later one".
    """
    def build():
        ped = G.Ped()
        zf, zm, kids = G.add_nuclear(ped, "FZ", "Z", 2)
        lone = "A1" if not first_has_parents else "ZZ9"
        ped.add(lone, "F" + lone, clone_of=kids[0], sex=1)
        sf = ped.add("SHZF", "SHZF", sex=1, clone_of=zf, emit=False)
        sm = ped.add("SHZM", "SHZM", sex=2, clone_of=zm, emit=False)
        ped.add("A2", "FA2", father=sf, mother=sm, sex=2)
        return ped
    return build


def both_parented(sib_on_first):
    """Both copies declare a couple; only one of them also has a declared sib.

    `sib_on_first` puts the sib on the copy that sorts first, so "keep the connected one"
    and "keep the later one" again disagree.
    """
    def build():
        ped = G.Ped()
        n_first, n_second = (2, 1) if sib_on_first else (1, 2)
        af, am, akids = G.add_nuclear(ped, "FA", "A", n_first)
        bf = ped.add("B_F", "FB", sex=1)
        bm = ped.add("B_M", "FB", sex=2)
        ped.add("B_C1", "FB", father=bf, mother=bm, clone_of=akids[0], sex=1)
        for k in range(n_second - 1):
            ped.add("B_C%d" % (k + 2), "FB", father=bf, mother=bm, sex=2)
        sf = ped.add("SHAF", "SHAF", sex=1, clone_of=af, emit=False)
        sm = ped.add("SHAM", "SHAM", sex=2, clone_of=am, emit=False)
        ped.add("W1", "FW", father=sf, mother=sm, sex=2)
        return ped
    return build


def dup_only():
    """A duplicate and nothing else: the control that says the line needs company."""
    def build():
        ped = G.Ped()
        _f, _m, kids = G.add_nuclear(ped, "FA", "A", 2)
        G.add_nuclear(ped, "FB", "B", 2)
        ped.add("B_D", "FB", clone_of=kids[0], sex=1)
        return ped
    return build


def with_children(kids_on_first):
    """One copy declares no parents but *is* a declared parent of two children.

    Children are the third kind of declared 1st-degree relative and the shapes above
    exercise none, so this is the one that says whether they count.
    """
    def build():
        ped = G.Ped()
        ph = G.add_couple(ped, "PHk", "PHk", emit=False)
        # the two copies: one a lone singleton, one the father of a nuclear family
        dad_fid, lone = ("FD", "ZL1") if kids_on_first else ("FZ", "AL1")
        dad = ped.add(dad_fid + "_F", dad_fid, father=ph[0], mother=ph[1], sex=1)
        mo = ped.add(dad_fid + "_M", dad_fid, sex=2)
        for k in range(2):
            ped.add("%s_C%d" % (dad_fid, k + 1), dad_fid, father=dad, mother=mo,
                    sex=1 + (k % 2))
        ped.add(lone, "F" + lone, clone_of=dad, sex=1)
        ped.add("OS1", "FOS", father=ph[0], mother=ph[1], sex=2)   # the FS edge
        return ped
    return build


SHAPES = {
    "kids_first": with_children(True),
    "kids_last": with_children(False),
    "solo_q1": solo("Q1", "Q4"),
    "solo_q2": solo("Q2", "Q4"),
    "solo_q3": solo("Q3", "Q9"),
    "lone_first": parented(False),      # unparented copy sorts first
    "lone_last": parented(True),        # unparented copy sorts last
    "sib_first": both_parented(True),
    "sib_last": both_parented(False),
    "dup_only": dup_only(),
}


def declared_degree1(fam, iid):
    """How many declared 1st-degree relatives of `iid` the `.fam` actually carries."""
    rows = [l.split() for l in open(fam)]
    me = next(r for r in rows if r[1] == iid)
    ids = {r[1] for r in rows}
    n = sum(1 for p in (me[2], me[3]) if p != "0" and p in ids)
    if me[2] != "0":
        n += sum(1 for r in rows if r[1] != iid and (r[2], r[3]) == (me[2], me[3]))
    n += sum(1 for r in rows if iid in (r[2], r[3]))
    return n


def run(name, seed):
    tag = "%s_%d" % (name, seed)
    d = os.path.join(WORK, tag)
    bed = os.path.join(d, tag + ".bed")
    if not os.path.exists(os.path.join(d, "kingbuild.log")):
        ped = SHAPES[name]()
        BS.pad(ped, max(0, 105 - ped.n_emitted()))
        spec = G.Spec(tag, ped, G.AUTOSOMES, 40000, notes="duplicate keep/remove probe")
        os.makedirs(d, exist_ok=True)
        G.simulate(spec, seed, d)
        subprocess.run([KING, "-b", bed, "--build", "--cpus", "1"], cwd=d,
                       check=True, capture_output=True)
    return d, bed, open(os.path.join(d, "kingbuild.log")).read()


def main(seeds=(4242, 11, 909)):
    tally = {"connectivity_then_later": [0, 0], "later_id": [0, 0], "earlier_id": [0, 0]}
    half = [0, 0]
    for name in SHAPES:
        for seed in seeds:
            d, bed, log = run(name, seed)
            m = DUP.search(log)
            if not m:
                print("%-12s seed=%-5d no Duplicate line (log %d bytes, %s)"
                      % (name, seed, len(log),
                         "rules only" if log.strip() else "empty"))
                continue
            removed, kept = m.group(1), m.group(2)
            fam = bed[:-4] + ".fam"
            dr, dk = declared_degree1(fam, removed), declared_degree1(fam, kept)
            preds = {
                "connectivity_then_later": kept if dk != dr else max(removed, kept),
                "later_id": max(removed, kept),
                "earlier_id": min(removed, kept),
            }
            # "connectivity" predicts the KEPT id; compare directly.
            for h, who in preds.items():
                ok = who == kept
                tally[h][0] += ok
                tally[h][1] += not ok
            half[0] += "INFERENCE" in log
            half[1] += "INFERENCE" not in log
            print("%-12s seed=%-5d removed=%-8s kept=%-8s  deg1 %d vs %d   inference=%s"
                  % (name, seed, removed, kept, dr, dk, "yes" if "INFERENCE" in log else "NO"))
    print("\nrule                       predicts the kept copy   right / wrong")
    for h, (a, b) in tally.items():
        print("  %-38s %3d / %3d" % (h, a, b))
    print("\nDuplicate line seen with an INFERENCE line %d time(s), without one %d time(s)"
          % tuple(half))


if __name__ == "__main__":
    main() if len(sys.argv) == 1 else sys.exit(__doc__)
