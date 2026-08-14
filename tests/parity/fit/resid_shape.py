"""Which *shape* do the reference-only IBD2 calls have?

For every reference IBD2 call under 10 Mb that our caller does not reproduce, invert the
length over the pair's own usable segments and tag each candidate interval by its position
in the 64-marker word grid — and, separately, by its position relative to *our* calls.
A shape that recurs across pairs and datasets is the clause to interrogate.
"""

import json
import os
from collections import Counter

import numpy as np

import engine as E
import kingdata as kd
import resid19 as R
import resid_calls as RC

WORD = 64


def tag(ds, i, j, a, b, seg, ours):
    """A structural label for the candidate interval [a, b]."""
    sc = E.SegScan(ds, i, j, seg, E.BASE)
    ua, ub = a // WORD, b // WORD
    z = [int(sc.n0[k]) != 0 for k in range(ds.nwords)]
    mis = [int(sc.n1[k]) for k in range(ds.nwords)]
    t = []
    t.append("a=%s" % ("seg.lo" if a == seg[1] else
                       ("64u" if a % WORD == 0 else
                        ("64u+63" if a % WORD == 63 else "mid"))))
    t.append("b=%s" % ("seg.hi" if b == seg[2] else
                       ("64u+63" if b % WORD == WORD - 1 else
                        ("64u" if b % WORD == 0 else "mid"))))
    t.append("words=%d" % (ub - ua + 1))
    inner = range(ua + 1, ub)
    nz = sum(1 for k in inner if z[k])
    nd = sum(1 for k in inner if not z[k] and mis[k] >= 2)
    t.append("innerIBS0=%d" % nz)
    t.append("innerDirty=%d" % nd)
    hit = [(x, y) for x, y in ours if x == a or y == b]
    if any(x == a and y == b for x, y in ours):
        t.append("ours=exact")
    elif hit:
        t.append("ours=shares-end")
    return " ".join(t)


def main():
    shapes = Counter()
    lines = []
    for name, n1, n2 in RC.PAIRS:
        ds = kd.load(name)
        idx = {f[1]: k for k, f in enumerate(ds.fam)}
        i, j = sorted((idx[n1], idx[n2]))
        pos = ds.pos
        ours = R.calls_of(ds, i, j)
        oursp = [(a, b) for _s, a, b in ours]
        want = json.load(open(os.path.join(
            R.SEGLEN, "%s.IBD2Seg.json" % name))).get("%d,%d" % (i, j), [])
        gshort = sorted(int(pos[b] - pos[a]) for _s, a, b in ours
                        if int(pos[b] - pos[a]) < 10_000_000)
        _m, _gl, wleft = R._match(gshort, sorted(want), RC.TOL)
        for t in wleft:
            cands = RC.invert(ds, t)
            for a, b, seg in cands:
                lab = tag(ds, i, j, a, b, seg, oursp)
                shapes[lab] += 1
                lines.append((name, n1, n2, t, a, b, lab, len(cands)))
    print("=== candidate shapes for reference-only IBD2 calls (%d candidates)"
          % sum(shapes.values()))
    for k, v in shapes.most_common():
        print("  %4d  %s" % (v, k))
    print()
    print("=== the wall-spanning shape, pair by pair")
    for name, n1, n2, t, a, b, lab, nc in lines:
        if "a=64u+63" in lab and "b=64u" in lab and "innerIBS0=1" in lab:
            print("  %-10s %-8s %-8s %8.3f Mb  [%d,%d]  (1 of %d candidates)  %s"
                  % (name, n1, n2, t / 1e6, a, b, nc, lab))
    print()
    print("=== per target: how many candidates, and does exactly one carry the shape?")
    seen = {}
    for name, n1, n2, t, a, b, lab, nc in lines:
        key = (name, n1, n2, t)
        seen.setdefault(key, []).append(lab)
    hit = uniq = 0
    for key, labs in sorted(seen.items()):
        wall = [x for x in labs
                if "a=64u+63" in x and "b=64u" in x and "innerIBS0=1" in x]
        hit += bool(wall)
        uniq += len(wall) == 1
        print("  %-10s %-8s %-8s %8.3f Mb  %2d candidates, %d wall-spanning"
              % (key[0], key[1], key[2], key[3] / 1e6, len(labs), len(wall)))
    print("  -> %d of %d targets have a wall-spanning candidate (%d uniquely)"
          % (hit, len(seen), uniq))


if __name__ == "__main__":
    main()
