#!/usr/bin/env python3
"""How the reference numbers merged clusters `KING1`, `KING2`, … — measured, not guessed.

`build_shapes.py` found, on three of its twenty shapes, that merged clusters are *not*
numbered in family order, and recorded a hypothesis: that they are numbered by the
**relationship type** of the pair that joined them, `Dup/MZ` first, then `PO`, then `FS`,
ties broken by family order.  Three shapes is not a rule.  This file builds nineteen that
*discriminate* it from the hypotheses it is confounded with on those three — the largest
kinship of a joining pair, the cluster's size, and plain family order — and then pins the
clauses none of them could see.

**What it found.**  The type hypothesis survives (19 of 19, against 11 for kinship and 7
each for size and family order), and the mechanism behind it is sharper than a sort: the
reference works through the joining pairs **by relationship type**, and a cluster takes its
number the first time the queue *creates* it.  The giveaway is the `OriginalFamID` list,
which is in **absorption** order — a cluster whose `Dup/MZ` edge is `QBB–QBC` and whose
`FS` edge is `QBA–QBB` prints `QBB,QBC,QBA`, not the file order `QBA,QBB,QBC`.
`seeds` kills the kinship hypothesis outright: over eight fresh seeds of one two-type
shape, **4** have the `FS` pair scoring a higher kinship than the `PO` pair and the `PO`
cluster is `KING1` in all eight.

**And a second bug, from `gate`.**  Our merge gate was `kinship > 2^-2.5`.  The reference
uses the disjunction `--related` uses to decide on a `.kin0` row at `--degree d`,
`kinship >= 2^-(d+1.5) || PropIBD > 2^-(d+0.5)` — 19 of 19 against the kinship rule's 18.
The shape that separates them is `threequarter`, whose 3/4-sib pair sits at `kinship
0.1749`, *under* the cut, and merges on `PropIBD 0.3646`.

    python3 clusternum.py            # the full scorecard
    python3 clusternum.py seeds      # the kinship-vs-type discriminator, many seeds
    python3 clusternum.py gate       # what predicate actually merges two families
    python3 clusternum.py dump NAME  # one shape's reference numbering

Nothing here reads KING's source.  Filesets land in `work/clusternum/` (gitignored).
"""

from __future__ import annotations

import importlib.util
import itertools
import os
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
WORK = os.path.join(_HERE, "work", "clusternum")

KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
)
IMPL = os.environ.get("OPEN_KING", os.path.join(_ROOT, "target", "release", "open-king"))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G = _load("gc", os.path.join(_ROOT, "tests", "parity", "generate_corpus.py"))
BS = _load("bs", os.path.join(_HERE, "build_shapes.py"))


def shadow(ped, who, tag):
    return ped.add(tag, tag, sex=ped.get(who).sex, clone_of=who, emit=False)


def fs_pair(ped, tag, fid_a, fid_b, nk=2):
    """Two families whose fathers are undeclared full sibs."""
    ph = G.add_couple(ped, "PH" + tag, "PH" + tag, emit=False)
    G.add_nuclear(ped, fid_a, fid_a, nk, father_parents=ph)
    G.add_nuclear(ped, fid_b, fid_b, nk, father_parents=ph)


def po_pair(ped, tag, fid_a, fid_b, nk=2):
    """Two families joined by one undeclared parent-offspring link."""
    fa, _ma, _k = G.add_nuclear(ped, fid_a, fid_a, nk)
    sf = shadow(ped, fa, "SH" + tag + "F")
    pm = ped.add("SH" + tag + "M", "SH" + tag + "M", sex=2, emit=False)
    ped.add(fid_b + "_F", fid_b, father=sf, mother=pm, sex=1)
    mb = ped.add(fid_b + "_M", fid_b, sex=2)
    for k in range(nk):
        ped.add("%s_C%d" % (fid_b, k + 1), fid_b, father=fid_b + "_F", mother=mb,
                sex=1 + (k % 2))


def dup_pair(ped, _tag, fid_a, fid_b, nk=2):
    """Two families sharing one individual entered twice under different ids."""
    _fa, _ma, kids = G.add_nuclear(ped, fid_a, fid_a, nk)
    G.add_nuclear(ped, fid_b, fid_b, nk)
    ped.add(fid_b + "_D", fid_b, clone_of=kids[0], sex=1)


MAKERS = {"FS": fs_pair, "PO": po_pair, "DUP": dup_pair}


def chain(types):
    """One cluster per entry of `types`, laid out in the given family order.

    Family ids are `Q<i>A` / `Q<i>B`, so file order and FID sort order agree; the
    `sorted_fids` shapes below break that on purpose.
    """
    def build():
        ped = G.Ped()
        for i, t in enumerate(types):
            MAKERS[t](ped, "%d" % i, "Q%dA" % i, "Q%dB" % i)
        return ped
    return build


def mixed_cluster():
    """One cluster carrying BOTH an FS edge and a Dup edge, and one carrying only FS.

    Family order puts the FS-only cluster first.  If a cluster's number follows the
    strongest edge inside it, the mixed one is KING1; if it follows the weakest, or the
    first edge found in family order, it is KING2.
    """
    ped = G.Ped()
    fs_pair(ped, "x", "QAA", "QAB")            # cluster 1: FS only
    ph = G.add_couple(ped, "PHm", "PHm", emit=False)
    G.add_nuclear(ped, "QBA", "QBA", 2, father_parents=ph)
    _f, _m, kids = G.add_nuclear(ped, "QBB", "QBB", 2, father_parents=ph)
    G.add_nuclear(ped, "QBC", "QBC", 2)
    ped.add("QBC_D", "QBC", clone_of=kids[0], sex=1)   # Dup edge QBB–QBC
    return ped


def mixed_po_fs():
    """A cluster with both a PO edge and an FS edge, against a pure-FS cluster first."""
    ped = G.Ped()
    fs_pair(ped, "y", "RAA", "RAB")
    ph = G.add_couple(ped, "PHn", "PHn", emit=False)
    fa, _ma, _k = G.add_nuclear(ped, "RBA", "RBA", 2, father_parents=ph)
    G.add_nuclear(ped, "RBB", "RBB", 2, father_parents=ph)
    sf = shadow(ped, fa, "SHnF")
    pm = ped.add("SHnM", "SHnM", sex=2, emit=False)
    ped.add("RBC_F", "RBC", father=sf, mother=pm, sex=1)
    mo = ped.add("RBC_M", "RBC", sex=2)
    for k in range(2):
        ped.add("RBC_C%d" % (k + 1), "RBC", father="RBC_F", mother=mo, sex=1 + (k % 2))
    return ped


def sorted_fids():
    """Three FS clusters whose FILE order and FID SORT order disagree.

    File order Z*, M*, A*; sort order A*, M*, Z*.  Family order alone cannot say which
    of the two the tie-break inside a type uses.
    """
    ped = G.Ped()
    for tag in ("Z", "M", "A"):
        fs_pair(ped, tag, tag + "1", tag + "2")
    return ped


def sorted_fids_po():
    """The same disagreement, with PO edges, so the tie-break is measured inside PO too."""
    ped = G.Ped()
    for tag in ("Z", "M", "A"):
        po_pair(ped, tag, tag + "1", tag + "2")
    return ped


def threequarter():
    """A cluster joined by a **3/4 sib** pair (φ=0.1875) — above the merge cut but neither
    `PO` nor `FS` — against an `FS` cluster first and a `PO` cluster after it.

    The question is where a type outside {Dup/MZ, PO, FS} lands in the priority list.
    """
    ped = G.Ped()
    fs_pair(ped, "q", "TAA", "TAB")
    gf, gm = G.add_couple(ped, "PHt", "PHt", emit=False)          # the shared grandparents
    dad = ped.add("TB_DAD", "PHt", sex=1, emit=False)             # the shared father
    m1 = ped.add("TB_M1", "PHt", sex=2, father=gf, mother=gm, emit=False)
    m2 = ped.add("TB_M2", "PHt", sex=2, father=gf, mother=gm, emit=False)
    for fid, mo in (("TBA", m1), ("TBB", m2)):
        kid = ped.add(fid + "_K", fid, father=dad, mother=mo, sex=1)
        sp = ped.add(fid + "_S", fid, sex=2)
        for k in range(2):
            ped.add("%s_C%d" % (fid, k + 1), fid, father=kid, mother=sp, sex=1 + (k % 2))
    po_pair(ped, "q2", "TCA", "TCB")
    return ped


SHAPES = {
    # all one type: the tie-break inside a type
    "fs_fs_fs": chain(["FS", "FS", "FS"]),
    "fs_dup_fs_dup": chain(["FS", "DUP", "FS", "DUP"]),
    "po_fs_po_fs": chain(["PO", "FS", "PO", "FS"]),
    "threequarter": threequarter,
    "po_po_po": chain(["PO", "PO", "PO"]),
    "dup_dup_dup": chain(["DUP", "DUP", "DUP"]),
    # two types, in the order that makes family order and type priority disagree
    "po_then_dup": chain(["PO", "DUP"]),
    "fs_then_po": chain(["FS", "PO"]),
    "fs_then_dup": chain(["FS", "DUP"]),
    "dup_then_fs": chain(["DUP", "FS"]),          # the two agree: a control
    "po_fs": chain(["PO", "FS"]),                 # the two agree: a control
    # three types, every permutation that discriminates
    "dup_fs_po": chain(["DUP", "FS", "PO"]),
    "fs_po_dup2": chain(["FS", "PO", "DUP"]),
    "po_dup_fs": chain(["PO", "DUP", "FS"]),
    "fs_dup_po": chain(["FS", "DUP", "PO"]),
    # a cluster with two different edge types inside it
    "mixed_cluster": mixed_cluster,
    "mixed_po_fs": mixed_po_fs,
    # file order vs FID sort order
    "sorted_fids": sorted_fids,
    "sorted_fids_po": sorted_fids_po,
}


def make(name, seed, nsnp=40000, total=110):
    d = os.path.join(WORK, "%s_%d" % (name, seed))
    stem = "%s_%d" % (name, seed)
    bed = os.path.join(d, stem + ".bed")
    if not os.path.exists(os.path.join(d, "kingupdateids.txt")) and not os.path.exists(
            os.path.join(d, ".done")):
        ped = SHAPES[name]()
        BS.pad(ped, max(0, total - ped.n_emitted()))
        spec = G.Spec(stem, ped, G.AUTOSOMES, nsnp, notes="cluster numbering probe")
        os.makedirs(d, exist_ok=True)
        G.simulate(spec, seed, d)
        subprocess.run([KING, "-b", bed, "--build", "--cpus", "1"], cwd=d,
                       check=True, capture_output=True)
        open(os.path.join(d, ".done"), "w").close()
    return d, bed


def numbering(d):
    """`{KING<k>: [original FIDs, in file order]}` read off the reference's updateids."""
    p = os.path.join(d, "kingupdateids.txt")
    if not os.path.exists(p):
        return {}
    out = {}
    for line in open(p):
        fid, _iid, new, _new_iid = line.rstrip("\n").split("\t")
        out.setdefault(new, [])
        if fid not in out[new]:
            out[new].append(fid)
    return out


def our_numbering(d, bed):
    o = os.path.join(d, "impl")
    os.makedirs(o, exist_ok=True)
    subprocess.run([IMPL, "-b", bed, "--build", "--cpus", "1"], cwd=o,
                   check=True, capture_output=True)
    return numbering(o)


def edges(d, bed):
    """Every cross-family pair the reference itself calls 1st-degree or duplicate.

    Read from `--related`'s `.kin0`, which is the same screen clustering uses, so the
    joining edge of each cluster is observed rather than assumed from the pedigree.
    """
    p = os.path.join(d, "king.kin0")
    if not os.path.exists(p):
        subprocess.run([KING, "-b", bed, "--related", "--degree", "1", "--cpus", "1"],
                       cwd=d, check=True, capture_output=True)
    rows = []
    if not os.path.exists(p):
        return rows
    with open(p) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            rows.append((f["FID1"], f["FID2"], f.get("InfType", "?"),
                         float(f["Kinship"])))
    return rows


ORDER = {"Dup/MZ": 0, "PO": 1, "FS": 2}


def predict(fids, groups, ev):
    """The numbering each hypothesis predicts, as `{hypothesis: [cluster, …]}`.

    `groups` is the list of merged clusters as lists of FIDs, in file order; `ev` maps a
    frozenset of two FIDs to `(InfType, kinship)`.
    """
    pos = {f: i for i, f in enumerate(fids)}

    def cluster_edges(g):
        return [ev[frozenset(p)] for p in itertools.combinations(g, 2)
                if frozenset(p) in ev]

    keyed = []
    for g in groups:
        es = cluster_edges(g)
        best_type = min((ORDER.get(t, 9) for t, _k in es), default=9)
        worst_type = max((ORDER.get(t, 9) for t, _k in es), default=9)
        best_kin = max((k for _t, k in es), default=0.0)
        keyed.append((g, best_type, worst_type, best_kin, min(pos[f] for f in g)))

    out = {}
    out["family"] = [g for g, *_ in sorted(keyed, key=lambda r: r[4])]
    out["type"] = [g for g, *_ in sorted(keyed, key=lambda r: (r[1], r[4]))]
    out["type_weakest"] = [g for g, *_ in sorted(keyed, key=lambda r: (r[2], r[4]))]
    out["kinship"] = [g for g, *_ in sorted(keyed, key=lambda r: (-r[3], r[4]))]
    out["size"] = [g for g, *_ in sorted(keyed, key=lambda r: (-len(g), r[4]))]
    return out


def score(names=None, seed=4242, verbose=True):
    names = names or list(SHAPES)
    tally = {h: [0, 0] for h in ("family", "type", "type_weakest", "kinship", "size")}
    ours = [0, 0]
    for name in names:
        d, bed = make(name, seed)
        ref = numbering(d)
        if not ref:
            print("%-16s (no merge)" % name)
            continue
        fids = []
        for line in open(bed[:-4] + ".fam"):
            f = line.split()[0]
            if f not in fids:
                fids.append(f)
        ev = {}
        for a, b, t, k in edges(d, bed):
            key = frozenset((a, b))
            if key not in ev or ORDER.get(t, 9) < ORDER.get(ev[key][0], 9):
                ev[key] = (t, k)
        groups = [ref[k] for k in sorted(ref, key=lambda s: int(s[4:]))]
        got = groups                              # reference order, by construction
        pred = predict(fids, groups, ev)
        mine = our_numbering(d, bed)
        mine_groups = [mine[k] for k in sorted(mine, key=lambda s: int(s[4:]))]
        ours[0] += mine_groups == got
        ours[1] += mine_groups != got
        if verbose:
            desc = []
            for k in sorted(ref, key=lambda s: int(s[4:])):
                es = sorted({ev[frozenset(p)][0]
                             for p in itertools.combinations(ref[k], 2)
                             if frozenset(p) in ev})
                desc.append("%s=%s[%s]" % (k, "+".join(ref[k]), ",".join(es) or "?"))
            print("%-16s %s" % (name, "  ".join(desc)))
        for h, order in pred.items():
            ok = order == got
            tally[h][0] += ok
            tally[h][1] += not ok
            if verbose and not ok:
                print("      %-13s predicts %s" % (
                    h, " ".join("+".join(g) for g in order)))
    print("\nhypothesis        agrees / disagrees")
    for h, (a, b) in tally.items():
        print("  %-14s %3d / %3d" % (h, a, b))
    print("  %-14s %3d / %3d  (our binary)" % ("open-king", ours[0], ours[1]))


def seeds(name="fs_then_po", how_many=8, base=7100):
    """The kinship-vs-type discriminator: many seeds of one two-type shape.

    Prints, per seed, the joining kinship of each cluster and which hypothesis the
    reference's numbering agrees with.  If numbering followed the kinship value the
    order would flip on the seeds where the FS pair outscores the PO pair.
    """
    flips = 0
    for s in range(base, base + how_many):
        d, bed = make(name, s)
        ref = numbering(d)
        if not ref:
            print("seed %-6d (no merge)" % s)
            continue
        ev = {}
        for a, b, t, k in edges(d, bed):
            key = frozenset((a, b))
            if key not in ev or ORDER.get(t, 9) < ORDER.get(ev[key][0], 9):
                ev[key] = (t, k)
        cells = []
        for k in sorted(ref, key=lambda s2: int(s2[4:])):
            es = [ev[frozenset(p)] for p in itertools.combinations(ref[k], 2)
                  if frozenset(p) in ev]
            t = min(es, key=lambda e: ORDER.get(e[0], 9))[0] if es else "?"
            kin = max((x[1] for x in es), default=0.0)
            cells.append("%s=%s(%s φ=%.4f)" % (k, "+".join(ref[k]), t, kin))
        # would sorting by kinship give the same order?
        kins = []
        for k in sorted(ref, key=lambda s2: int(s2[4:])):
            es = [ev[frozenset(p)] for p in itertools.combinations(ref[k], 2)
                  if frozenset(p) in ev]
            kins.append(max((x[1] for x in es), default=0.0))
        by_kin = kins == sorted(kins, reverse=True)
        flips += not by_kin
        print("seed %-6d %s   kinship-descending? %s" % (s, "  ".join(cells),
                                                         "yes" if by_kin else "NO"))
    print("\n%d of %d seeds contradict the kinship ordering" % (flips, how_many))


def gate(names=None, seed=4242):
    """What predicate actually merges two families?

    Our binary merges on `kinship > 2^-2.5` alone.  The reference's `.kin0` admits a pair
    at `--degree d` on a **disjunction** — `kinship >= 2^-(d+1.5)` OR
    `PropIBD > 2^-(d+0.5)` — and this checks whether the merge follows that same rule, by
    testing whether the reference's clusters are exactly the connected components of the
    cross-family pairs its own `.kin0` lists.
    """
    names = names or list(SHAPES)
    kin_cut, prop_cut = 2 ** -2.5, 2 ** -1.5
    agree = {"kin_only": [0, 0], "disjunction": [0, 0], "any_kin0_row": [0, 0]}
    for name in names:
        d, bed = make(name, seed)
        ref = numbering(d)
        merged = {frozenset(g) for g in ref.values()}
        rows = []
        p = os.path.join(d, "king.kin0")
        if not os.path.exists(p):
            edges(d, bed)
        with open(p) as fh:
            hdr = next(fh).rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                rows.append((f["FID1"], f["FID2"], float(f["Kinship"]),
                             float(f["PropIBD"])))

        def components(pred):
            par = {}

            def find(x):
                par.setdefault(x, x)
                while par[x] != x:
                    par[x] = par[par[x]]
                    x = par[x]
                return x
            touched = set()
            for a, b, k, pr in rows:
                if a == b or not pred(k, pr):
                    continue
                touched |= {a, b}
                ra, rb = find(a), find(b)
                if ra != rb:
                    par[ra] = rb
            out = {}
            for f in touched:
                out.setdefault(find(f), set()).add(f)
            return {frozenset(v) for v in out.values()}

        cases = {
            "kin_only": lambda k, pr: k > kin_cut,
            "disjunction": lambda k, pr: k >= kin_cut or pr > prop_cut,
            "any_kin0_row": lambda k, pr: True,
        }
        line = []
        for h, pred in cases.items():
            ok = components(pred) == merged
            agree[h][0] += ok
            agree[h][1] += not ok
            line.append("%s=%s" % (h, "ok" if ok else "NO"))
        print("%-16s %d cluster(s)  %s" % (name, len(merged), "  ".join(line)))
    print("\npredicate       agrees / disagrees")
    for h, (a, b) in agree.items():
        print("  %-13s %3d / %3d" % (h, a, b))


def dump(name, seed=4242):
    d, bed = make(name, seed)
    print("--- fam (merged families only)")
    fids = {f for g in numbering(d).values() for f in g}
    for line in open(bed[:-4] + ".fam"):
        if line.split()[0] in fids:
            print("   ", line.rstrip())
    print("--- reference numbering")
    for k, g in sorted(numbering(d).items(), key=lambda kv: int(kv[0][4:])):
        print("   ", k, "=", "+".join(g))
    print("--- cross-family 1st-degree edges")
    for a, b, t, k in edges(d, bed):
        if a != b:
            print("    %-6s %-6s %-8s %.4f" % (a, b, t, k))
    print("--- reference kingbuild.log")
    p = os.path.join(d, "kingbuild.log")
    print(open(p).read() if os.path.exists(p) else "(none)")


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "score"
    if mode == "score":
        score(sys.argv[2:] or None)
    elif mode == "seeds":
        seeds(*(sys.argv[2:3] or ["fs_then_po"]))
    elif mode == "gate":
        gate(sys.argv[2:] or None)
    elif mode == "dump":
        dump(sys.argv[2])
    else:
        sys.exit(__doc__)
