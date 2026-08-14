"""Aligned context around the 154 localised MaxIBD2 intervals.

Prints, for the start and end boundary and for the interior, the per-word IBS1 counts
relative to the boundary, so the run predicate can be read as a joint condition on
neighbouring words rather than guessed one word at a time.

    python3 maxctx.py [thresh]
"""

import sys
from collections import Counter

import numpy as np

import kingdata as kd
import maxfit as M
import maxdump as D

WORD = 64


def main():
    T = int(sys.argv[1]) if len(sys.argv) > 1 else 5      # "dirty" means ibs1 >= T
    tg = M.all_targets()
    starts, ends = [], []
    interior_dirty = []
    inter_runs = Counter()
    for ds, i, j, t in tg:
        cs = D.candidates(ds, i, j, t)
        if len(cs) != 1:
            continue
        a, b, w0, w1 = cs[0]
        _, n0, n1 = M.counts(ds, i, j)

        def g(w):
            return int(n1[w]) if 0 <= w < len(n1) else None

        if a > w0:
            starts.append((g(a - 3), g(a - 2), g(a - 1), g(a), g(a + 1)))
        if b < w1:
            ends.append((g(b - 1), g(b), g(b + 1), g(b + 2), g(b + 3)))
        # interior structure: longest run of dirty words strictly inside [a..b]
        d = (n1[a:b + 1] >= T)
        best = cur = 0
        for x in d.tolist():
            cur = cur + 1 if x else 0
            best = max(best, cur)
        inter_runs[best] += 1
        for k in range(a, b + 1):
            if int(n1[k]) >= T:
                interior_dirty.append((g(k - 1), g(k), g(k + 1), k == b, k == a))

    def show(name, rows, labels):
        print(f"\n{name} ({len(rows)} boundaries), IBS1 by offset")
        for c, lab in enumerate(labels):
            vals = [r[c] for r in rows if r[c] is not None]
            if not vals:
                continue
            print(f"   {lab:>5}: min {min(vals):3d}  max {max(vals):3d}  "
                  f"<{T}: {sum(v < T for v in vals):3d}  >={T}: {sum(v >= T for v in vals):3d}")

    show("START", starts, ["a-3", "a-2", "a-1", "a", "a+1"])
    show("END", ends, ["b-1", "b", "b+1", "b+2", "b+3"])
    print(f"\nlongest run of words with IBS1>={T} strictly inside a called interval:",
          sorted(inter_runs.items()))
    print(f"\ninterior words with IBS1>={T}: {len(interior_dirty)}")
    ctx = Counter()
    for lft, own, rgt, is_b, is_a in interior_dirty:
        ctx[(lft is not None and lft >= T, rgt is not None and rgt >= T, is_b)] += 1
    print("  (left dirty, right dirty, is the last word):")
    for k, v in sorted(ctx.items()):
        print("   ", k, v)


if __name__ == "__main__":
    main()
