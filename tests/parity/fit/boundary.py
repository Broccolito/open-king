"""Fit the IBD2 word predicate and boundary convention against MaxIBD2.

Each `--ibs` MaxIBD2 value pins one segment exactly (see maxibd2.py).  For a candidate
word predicate `good(k)` and a candidate boundary convention, this script re-derives the
segment from the observed start word and checks whether it lands on the observed end word,
and whether the observed start is where the predicate says a run begins.  Sweeping the
candidate space over ~150 pinned segments is what selects the rule.

    python3 boundary.py <dir-with-ibs_*.ibs>
"""

import sys
from collections import Counter

import numpy as np

import kingdata as kd
import maxibd2 as M

PC = np.bitwise_count


def cases(dirname, datasets=None):
    """(ds, i, j, u, v, w0, w1) for every unambiguously located MaxIBD2 segment."""
    out = []
    for name in datasets or kd.DATASETS:
        ds = kd.load(name)
        segw = [(-(-lo // 64), (hi + 1) // 64 - 1) for _, lo, hi in ds.segs]
        m = M.collect(dirname, ds)
        for (i, j), t in sorted(m.items()):
            n0, n1 = M.profile(ds, i, j)
            hits = M.endpoints_for(ds, i, j, t)
            good = [(u, v) for u, v in hits
                    if n0[u:v + 1].sum() == 0 and (n1[u:v + 1] <= 2).mean() > 0.5]
            if len(good) != 1:
                continue
            u, v = good[0]
            sw = [k for k in segw if k[0] <= u <= k[1]]
            if not sw:
                continue
            out.append((ds, i, j, u, v, sw[0][0], sw[0][1]))
    return out


def predicate(n0, n1, kind, t):
    """Candidate `good word` masks."""
    if kind == "count":
        return (n0 == 0) & (n1 <= t)
    if kind == "count_only1":        # ignore IBS0, threshold IBS1 alone
        return n1 <= t
    if kind == "pair":               # break when this word and the next are both busy
        busy = n1 > t
        g = ~busy
        g[:-1] |= busy[:-1] & ~busy[1:]
        return (n0 == 0) & g
    raise ValueError(kind)


def check(dirname, datasets=None):
    cs = cases(dirname, datasets)
    print(f"{len(cs)} pinned segments")
    print(f"{'kind':<12}{'t':>3}{'end ok':>8}{'start ok':>10}{'both':>7}")
    best = None
    for kind in ("count", "count_only1", "pair"):
        for t in range(0, 9):
            end_ok = start_ok = both = 0
            for ds, i, j, u, v, w0, w1 in cs:
                n0, n1 = M.profile(ds, i, j)
                g = predicate(n0, n1, kind, t)
                if not g[u]:
                    continue
                k = u
                while k + 1 <= w1 and g[k + 1]:
                    k += 1
                v_pred = min(k + 1, w1)
                s_ok = (u == w0) or (not g[u - 1])
                e_ok = v_pred == v
                end_ok += e_ok
                start_ok += s_ok
                both += e_ok and s_ok
            print(f"{kind:<12}{t:3d}{end_ok:8d}{start_ok:10d}{both:7d}")
            if best is None or both > best[0]:
                best = (both, kind, t)
    print("best:", best)
    return cs


def residuals(dirname, kind, t, datasets=None):
    """Where the winning candidate still misses."""
    cs = cases(dirname, datasets)
    bad = Counter()
    for ds, i, j, u, v, w0, w1 in cs:
        n0, n1 = M.profile(ds, i, j)
        g = predicate(n0, n1, kind, t)
        k = u
        while k + 1 <= w1 and g[k + 1]:
            k += 1
        v_pred = min(k + 1, w1)
        if v_pred != v or not (u == w0 or not g[u - 1]):
            bad[(v_pred - v, u == w0, int(n1[v]), int(n1[min(v + 1, len(n1) - 1)]))] += 1
            if sum(bad.values()) <= 12:
                lo = " ".join(f"{w}:{int(n0[w])},{int(n1[w])}" for w in range(max(u - 2, 0), u + 2))
                hi = " ".join(f"{w}:{int(n0[w])},{int(n1[w])}"
                              for w in range(v - 2, min(v + 3, len(n1))))
                print(f"  {ds.name} {ds.fam[i][1]}x{ds.fam[j][1]} seg {w0}..{w1} "
                      f"obs {u}..{v} pred ..{v_pred}\n     L {lo}\n     R {hi}")
    print("residual signature counts:", bad.most_common(10))


if __name__ == "__main__":
    d = sys.argv[1]
    if len(sys.argv) > 2 and sys.argv[2] == "residuals":
        residuals(d, sys.argv[3], int(sys.argv[4]))
    else:
        check(d)
