"""Per-pair dump of the MaxIBD2 inversion, with usable-segment context.

    python3 maxdump.py [t2] [ds...] [-v]

For every pair carrying a MaxIBD2 this prints the matching word-aligned interval and the
per-word (IBS0, IBS1) counts around it *inside its usable segment*, so the run predicate
and the two endpoint rules can be read straight off.  `|` marks a usable-segment edge.
"""

import sys
from collections import Counter

import numpy as np

import kingdata as kd
import maxfit as M

WORD = 64


def seg_of(ds, lo_m):
    for _, lo, hi in ds.segs:
        if lo <= lo_m <= hi:
            return lo, hi
    return None


def candidates(ds, i, j, t):
    """Every word-aligned, IBS0-free interval whose span equals the target."""
    _, n0, n1 = M.counts(ds, i, j)
    out = []
    for lo_m, hi_m, seglo, seghi in M.matching_intervals(ds, t):
        if lo_m % WORD or hi_m % WORD != WORD - 1:
            continue
        a, b = lo_m // WORD, (hi_m + 1) // WORD - 1
        if int(n0[a:b + 1].max()) != 0:
            continue
        w0 = -(-seglo // WORD)
        w1 = (seghi + 1) // WORD - 1
        out.append((a, b, w0, w1))
    return out


def main():
    argv = [a for a in sys.argv[1:] if a != "-v"]
    verbose = "-v" in sys.argv
    t2 = int(argv[0]) if argv else 0
    names = argv[1:] or None
    tg = M.all_targets(names)
    nuniq = Counter()
    kinds = Counter()
    starts = []
    ends = []
    for ds, i, j, t in tg:
        cs = candidates(ds, i, j, t)
        nuniq[len(cs)] += 1
        if len(cs) != 1:
            continue
        a, b, w0, w1 = cs[0]
        _, n0, n1 = M.counts(ds, i, j)
        # what the two endpoints demand of the predicate
        sedge = a == w0
        eedge = b == w1
        if not sedge:
            starts.append((int(n1[a]), int(n0[a - 1]), int(n1[a - 1])))
        if not eedge:
            ends.append((int(n1[b]), int(n0[b + 1]), int(n1[b + 1])))
        kinds[(sedge, eedge)] += 1
        if verbose:
            lo = max(w0 - 1, a - 3)
            hi = min(w1 + 1, b + 3)
            cells = []
            for w in range(lo, hi + 1):
                s = ""
                if w == w0:
                    s += "|"
                if w == a:
                    s += "["
                s += f"{int(n0[w])}/{int(n1[w])}"
                if w == b:
                    s += "]"
                if w == w1:
                    s += "|"
                cells.append(s)
            print(f"{ds.name:12s} {i:3d},{j:3d} w{a}..{b} seg w{w0}..{w1}  " + " ".join(cells))
    print("\ncandidates per pair:", sorted(nuniq.items()))
    print("(start at seg edge, end at seg edge):", dict(kinds))

    def q(name, rows):
        print(f"\n{name}: {len(rows)} interior endpoints")
        own = Counter(r[0] for r in rows)
        print("  IBS1 of the boundary word itself:", sorted(own.items()))
        nb0 = Counter(r[1] for r in rows)
        print("  IBS0 of the word just outside:   ", sorted(nb0.items()))
        nb1 = Counter(r[2] for r in rows)
        print("  IBS1 of the word just outside:   ", sorted(nb1.items())[:12], "...")
        print("  min IBS1 just outside:", min(r[2] for r in rows),
              "  max IBS1 inside:", max(r[0] for r in rows))

    q("START", starts)
    q("END", ends)


if __name__ == "__main__":
    main()
