"""Row-level `.seg` diff between a candidate rule and the captured reference.

    python3 segdiff.py [dataset ...] [-p KEY=VALUE]...

Prints, for every row that is not exact, the reference's four columns against the
candidate's, the per-segment calls behind them, and the sign of each error — so a residual
can be attributed to a specific called segment instead of to a dataset.
"""

import sys
from dataclasses import replace

import kingdata as kd
import engine as E

WORD = 64


def pair_calls(ds, i, j, p=E.BASE, min_bp=E.SEGLEN):
    """Every call for one pair as (kind, seg_index, lo_word, hi_word, bp)."""
    out = []
    for n, seg in enumerate(ds.segs):
        sc = E.SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        for kind, calls in (("2", sc.ibd2(ds.pos, min_bp)), ("1", sc.ibd1(ds.pos, min_bp))):
            for lo, hi in calls:
                out.append((kind, n, lo, hi, int(ds.pos[hi] - ds.pos[lo])))
    return out


def main():
    argv = sys.argv[1:]
    kw = {}
    while "-p" in argv:
        k = argv.index("-p")
        key, _, val = argv[k + 1].partition("=")
        kw[key] = eval(val)  # noqa: S307
        del argv[k:k + 2]
    p = replace(E.BASE, **kw) if kw else E.BASE
    names = [a for a in argv if not a.startswith("-")] or kd.DATASETS
    verbose = "-v" in argv
    for name in names:
        ds = kd.load(name)
        d = ds.denom
        bad = 0
        for i, j in ds.pairs():
            if (i, j) not in ds.ref:
                continue
            i1, i2, lg, _ = E.call_pair(ds, i, j, p)
            a1, a2, ap, at = ds.ref[(i, j)]
            g1, g2 = i1 / d, i2 / d
            gp = g2 + g1 / 2
            if (kd.fmt4(g1), kd.fmt4(g2), kd.fmt4(gp),
                    kd.inf_type(g1, g2, gp)) == (a1, a2, ap, at):
                continue
            bad += 1
            print("%-12s %-9s/%-9s  ref %.4f %.4f %.4f %-6s | got %.4f %.4f %.4f %-6s"
                  " | d1 %+.4f d2 %+.4f"
                  % (name, ds.fam[i][1], ds.fam[j][1], a1, a2, ap, at,
                     g1, g2, gp, kd.inf_type(g1, g2, gp), g1 - a1, g2 - a2))
            if verbose:
                for kind, n, lo, hi, bp in pair_calls(ds, i, j, p):
                    print("      IBD%s seg%-3d w%d..w%d  markers %d..%d  %.3f Mb"
                          % (kind, n, lo // WORD, hi // WORD, lo, hi, bp / 1e6))
        print("%-12s %d inexact of %d" % (name, bad, len(ds.ref)))


if __name__ == "__main__":
    main()
