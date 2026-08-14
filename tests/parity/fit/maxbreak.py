"""Turn the 154 localised MaxIBD2 intervals into a labelled break / no-break dataset.

Every called interval [a..b] inside a usable segment implies, for whatever rule the
reference uses, that some word boundaries are *breaks* and every boundary strictly inside
is *not*.  Collecting both classes over the corpus turns rule-guessing into a
separability question that can be answered instead of argued.

Two start conventions are collected side by side:

    C1  segment = [p+1 .. q]   break at (a-1, a)   and (b, b+1)
    C2  segment = [p+2 .. q]   break at (a-2, a-1) and (b, b+1)

    python3 maxbreak.py
"""

from collections import Counter

import numpy as np

import kingdata as kd
import maxfit as M
import maxdump as D

WORD = 64


def features(ds, i, j):
    """Per-word feature table for a pair: ibs0, ibs1, hethet, called."""
    p0i, p1i = ds.p0[i], ds.p1[i]
    p0j, p1j = ds.p0[j], ds.p1[j]
    het_i = ~p0i & p1i
    het_j = ~p0j & p1j
    pc = np.bitwise_count
    return {
        "ibs0": pc(p0i & p0j & (p1i ^ p1j)).astype(np.int32),
        "ibs1": pc((het_i & p0j) | (p0i & het_j)).astype(np.int32),
        "hethet": pc(het_i & het_j).astype(np.int32),
        "called": pc((p0i | p1i) & (p0j | p1j)).astype(np.int32),
        "het1": pc(het_i).astype(np.int32),
        "het2": pc(het_j).astype(np.int32),
    }


def collect():
    """[(conv, label, left_feats, right_feats)] over every localised interval."""
    rows = []
    for ds, i, j, t in M.all_targets():
        cs = D.candidates(ds, i, j, t)
        if len(cs) != 1:
            continue
        a, b, w0, w1 = cs[0]
        f = features(ds, i, j)

        def at(w):
            return {k: int(v[w]) for k, v in f.items()}

        for conv, s in (("C1", a - 1), ("C2", a - 2)):
            if s >= w0:
                rows.append((conv, 1, at(s), at(s + 1), (ds.name, i, j, "start")))
        if b + 1 <= w1:
            for conv in ("C1", "C2"):
                rows.append((conv, 1, at(b), at(b + 1), (ds.name, i, j, "end")))
        for w in range(a, b):
            for conv in ("C1", "C2"):
                rows.append((conv, 0, at(w), at(w + 1), (ds.name, i, j, f"in{w}")))
        if a - 1 >= w0:                     # (a-1, a) is interior under C2 only
            rows.append(("C2", 0, at(a - 1), at(a), (ds.name, i, j, "a-1")))
    return rows


def report(rows, conv):
    r = [x for x in rows if x[0] == conv]
    brk = [x for x in r if x[1] == 1]
    non = [x for x in r if x[1] == 0]
    print(f"\n===== convention {conv}: {len(brk)} breaks, {len(non)} non-breaks")
    # IBS0 first: does any IBS0 force a break?
    print("  breaks    with ibs0>0 on either side:",
          sum(x[2]["ibs0"] or x[3]["ibs0"] for x in brk), f"/ {len(brk)}")
    print("  nonbreaks with ibs0>0 on either side:",
          sum(x[2]["ibs0"] or x[3]["ibs0"] for x in non), f"/ {len(non)}")
    # restrict to IBS0-free boundaries and look at the ibs1 pair
    b2 = [x for x in brk if not x[2]["ibs0"] and not x[3]["ibs0"]]
    n2 = [x for x in non if not x[2]["ibs0"] and not x[3]["ibs0"]]
    print(f"  IBS0-free boundaries: {len(b2)} breaks, {len(n2)} non-breaks")
    for name in ("ibs1", "hethet"):
        L = [x[2][name] for x in b2]
        R = [x[3][name] for x in b2]
        Ln = [x[2][name] for x in n2]
        Rn = [x[3][name] for x in n2]
        print(f"  {name:7s} break  left  min {min(L):3d} max {max(L):3d} | "
              f"right min {min(R):3d} max {max(R):3d}")
        print(f"  {name:7s} nobrk  left  min {min(Ln):3d} max {max(Ln):3d} | "
              f"right min {min(Rn):3d} max {max(Rn):3d}")
    # search the (A, B) contingency rule: break iff l>=A and r>=B (either orientation)
    best = []
    for A in range(0, 40):
        for B in range(0, 40):
            def hit(x):
                return ((x[2]["ibs1"] >= A and x[3]["ibs1"] >= B)
                        or (x[3]["ibs1"] >= A and x[2]["ibs1"] >= B))
            tp = sum(hit(x) for x in b2)
            fp = sum(hit(x) for x in n2)
            best.append((tp - fp, tp, fp, A, B))
    best.sort(reverse=True)
    print("  best symmetric contingency (score, tp, fp, A, B):")
    for row in best[:5]:
        print("     ", row, f"of {len(b2)} breaks / {len(n2)} non-breaks")
    # and the directed version
    best = []
    for A in range(0, 40):
        for B in range(0, 40):
            def hit(x):
                return x[2]["ibs1"] >= A and x[3]["ibs1"] >= B
            tp = sum(hit(x) for x in b2)
            fp = sum(hit(x) for x in n2)
            best.append((tp - fp, tp, fp, A, B))
    best.sort(reverse=True)
    print("  best directed contingency (score, tp, fp, A=left, B=right):")
    for row in best[:5]:
        print("     ", row)
    # single-word rule: break iff the right word is bad
    best = []
    for B in range(0, 40):
        tp = sum(x[3]["ibs1"] >= B for x in b2)
        fp = sum(x[3]["ibs1"] >= B for x in n2)
        best.append((tp - fp, tp, fp, B))
    best.sort(reverse=True)
    print("  best single-word (right) rule (score, tp, fp, B):", best[:3])
    best = []
    for B in range(0, 40):
        tp = sum(x[2]["ibs1"] >= B for x in b2)
        fp = sum(x[2]["ibs1"] >= B for x in n2)
        best.append((tp - fp, tp, fp, B))
    best.sort(reverse=True)
    print("  best single-word (left)  rule (score, tp, fp, B):", best[:3])


def main():
    rows = collect()
    for conv in ("C1", "C2"):
        report(rows, conv)


if __name__ == "__main__":
    main()
