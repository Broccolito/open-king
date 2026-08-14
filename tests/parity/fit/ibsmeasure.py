"""Score `--ibs`'s MaxIBD2 / Pr_IBD2 columns under the shared IBD2 caller.

The caller is `king_core::ibdseg::Scan::ibd2` transcribed; the only thing varied here is
how its calls are *measured* (`.seg`'s own-ends ruler vs `--ibs`'s word grid) and whether
the informativeness gate of `docs/research/13-informativeness-gate.md` is applied, and
over which words.

    python3 ibsmeasure.py [gate_scope ...]     # none | core | reported
"""

import os
import sys

import numpy as np

import kingdata as kd

WORD = 64
PC = np.bitwise_count
DIRTY = 5           # IBD2_HET_DIRTY
GATE = 10           # MIN_INFORMATIVE

IBSDIR = os.path.join(kd.GOLDEN, "core")


def targets(ds):
    """{(i, j): (MaxIBD2, Pr_IBD2)} as printed, by sample index."""
    idx = {iid: k for k, (fid, iid) in enumerate(ds.fam)}
    out = {}
    for suffix in (".ibs", ".ibs0"):
        p = os.path.join(IBSDIR, f"{ds.name}__ibs", "king" + suffix)
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            head = next(fh).rstrip("\n").split("\t")
            if "MaxIBD2" not in head:
                continue
            cm, cp = head.index("MaxIBD2"), head.index("Pr_IBD2")
            i1, i2 = head.index("ID1"), head.index("ID2")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                a, b = idx[f[i1]], idx[f[i2]]
                out[(min(a, b), max(a, b))] = (f[cm], f[cp])
    return out


def call(ds, i, j, gate_scope="none"):
    """(max_bp, total_bp) measured on the word grid, under one gate variant."""
    ibs0, ibs1, _, _ = ds.masks(i, j)
    n0, n1 = PC(ibs0), PC(ibs1)
    m2 = ds.p1[i] & ds.p1[j]              # inf2: both carry A1
    cum = np.concatenate(([0], np.cumsum(PC(m2).astype(np.int64))))
    pos = ds.pos
    best = total = 0
    for _, lo, hi in ds.segs:
        w0, w1 = -(-lo // WORD), (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        clean = [(n0[w] == 0) and (n1[w] < DIRTY) for w in range(w0, w1 + 1)]
        ok = list(clean)
        for k in range(1, len(clean) - 1):
            if not clean[k] and clean[k - 1] and clean[k + 1] and n0[w0 + k] == 0:
                ok[k] = True
        prev = None
        k = 0
        while k < len(ok):
            if not ok[k]:
                k += 1
                continue
            k0 = k
            while k < len(ok) and ok[k]:
                k += 1
            u, v = w0 + k0, w0 + k - 1
            e = w1 if v + 2 >= w1 else v + 1
            if gate_scope != "none":
                last = v if gate_scope == "core" else e
                if int(cum[last + 1] - cum[u]) < GATE:
                    continue
            a = u if prev is None else max(u, prev + 1)
            if a > e:
                continue
            prev = e
            ln = int(pos[WORD * e + WORD - 1] - pos[WORD * a])
            best = max(best, ln)
            total += ln
    return best, total


def main():
    scopes = sys.argv[1:] or ["none", "core", "reported"]
    for scope in scopes:
        n = mok = pok = 0
        per = {}
        for name in kd.DATASETS:
            ds = kd.load(name)
            for (i, j), (tm, tp) in sorted(targets(ds).items()):
                if float(tm) <= 0:
                    continue
                b, t = call(ds, i, j, scope)
                n += 1
                gm = f"{b:.3f}"
                gp = "%.4f" % (t / ds.denom)
                mok += gm == tm
                pok += gp == tp
                r = per.setdefault(name, [0, 0, 0])
                r[0] += 1
                r[1] += gm == tm
                r[2] += gp == tp
        print(f"gate={scope:9s}  n={n}  MaxIBD2 {mok}  Pr_IBD2 {pok}")
        for name, r in sorted(per.items()):
            print(f"    {name:12s} n={r[0]:4d} max={r[1]:4d} pr={r[2]:4d}")


if __name__ == "__main__":
    main()
