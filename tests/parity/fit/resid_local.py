"""Run `seglocal.py` over every wrong `.seg` IBD2 row and profile the guilty segments.

The output is the residual restated per *usable segment* instead of per pair: how many
segments each wrong row actually disagrees on, how big the disagreement is in marker
intervals, and what the segment looks like — its word count, its predicate string, how
many calls each side makes, and whether the disagreement sits at a segment fringe, at a
bridge, at a push, or at a plain interior run.

    python3 resid_local.py                 # every wrong row in the probeable datasets
    python3 resid_local.py --all           # every row, wrong or not (the control)
"""

import sys
from collections import Counter

import engine as E
import kingdata as kd
import resid19 as R
import seglocal as SL

WORD = 64
# `bigish` is 40 samples x 50 000 markers; localising all 763 of its rows is a different
# scale of job. The nine smaller datasets carry 20 of the 86 wrong rows.
SMALL = [n for n in kd.DATASETS if n != "bigish"]


def seg_shape(sc, ds, i, j):
    """A compact description of one usable segment for this pair."""
    n = sc.n
    w0 = sc.w0
    z = [int(sc.n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(sc.n1[w0 + k]) for k in range(n)]
    s = "".join("W" if z[k] else ("C" if mis[k] == 0 else ("x" if mis[k] == 1 else "y"))
                for k in range(n))
    return s


def classify(ds, i, j, seg, ourcalls):
    """Which clauses this segment's calls exercise."""
    sc = E.SegScan(ds, i, j, seg, E.BASE)
    tags = []
    if not ourcalls:
        return ["no-call"]
    for a, b in ourcalls:
        if a == sc.lo:
            tags.append("fringeL")
        if b == sc.hi:
            tags.append("fringeR")
    if len(ourcalls) > 1:
        tags.append("push")
    s = seg_shape(sc, ds, i, j)
    if "W" in s:
        tags.append("hasW")
    if "y" in s:
        tags.append("hasY")
    if all(c == "C" for c in s):
        tags.append("all-clean")
    return sorted(set(tags))


def main(only_wrong=True):
    rows = []
    for name in SMALL:
        ds = kd.load(name)
        ref = ds.ref
        for i, j in ds.pairs():
            if (i, j) not in ref:
                continue
            r = R.row(ds, i, j, E.SEGLEN, ref)
            if only_wrong and r["ok2"]:
                continue
            rows.append((name, i, j, r))
    print("=== localising %d rows" % len(rows))
    tally = Counter()
    guilty = []
    for name, i, j, r in rows:
        ds = kd.load(name)
        bad, _ = SL.per_segment(name, ds.fam[i][1], ds.fam[j][1], quiet=True)
        tally["rows"] += 1
        tally["rows_with_a_bad_segment"] += bool(bad)
        tally["bad_segments"] += len(bad)
        for k, seg, c, tot, refv, delta in bad:
            sc = E.SegScan(ds, i, j, seg, E.BASE)
            shape = seg_shape(sc, ds, i, j)
            cls = classify(ds, i, j, seg, c)
            gap = R._median_gap(ds)
            guilty.append((name, ds.fam[i][1], ds.fam[j][1], k, seg, shape, c, tot,
                           refv, delta, cls, sc))
            for t in cls:
                tally["cls:" + t] += 1
            tally["markers:%+d" % round(delta / gap)] += 1
    print("  " + "  ".join("%s=%d" % kv for kv in sorted(tally.items())
                           if not kv[0].startswith(("cls:", "markers:"))))
    print("  clause exposure of the guilty segments:")
    for k, v in sorted(tally.items()):
        if k.startswith("cls:"):
            print("     %-14s %d" % (k[4:], v))
    print("  size of the disagreement, in median marker gaps:")
    print("     " + "  ".join("%s:%d" % (k.split(":")[1], v)
                              for k, v in sorted(tally.items(),
                                                 key=lambda kv: _key(kv[0]))
                              if k.startswith("markers:")))
    print()
    print("=== the guilty segments")
    for (name, a, b, k, seg, shape, c, tot, refv, delta, cls, sc) in guilty:
        gap = R._median_gap(kd.load(name))
        print("  %-10s %-8s %-8s seg%2d chr%-2d w%d..%d n=%d  d=%+9d bp (%+6.1f markers)"
              % (name, a, b, k, seg[0], sc.w0, sc.w1, sc.n, delta, delta / gap))
        print("      shape %s" % shape)
        print("      lo=%d(w%d+%d) hi=%d(w%d+%d)  our calls %s"
              % (sc.lo, sc.lo // WORD, sc.lo % WORD, sc.hi, sc.hi // WORD, sc.hi % WORD,
                 ["[%d,%d] w%d+%d..w%d+%d" % (x, y, x // WORD, x % WORD,
                                              y // WORD, y % WORD) for x, y in c]))
        print("      tags %s" % ",".join(cls))


def _key(s):
    if not s.startswith("markers:"):
        return (0, s)
    return (1, int(s.split(":")[1]))


if __name__ == "__main__":
    main(only_wrong="--all" not in sys.argv)
