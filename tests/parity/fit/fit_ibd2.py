"""Fit the IBD2 rule directly against `--ibs`'s MaxIBD2, one exact number per pair.

`MaxIBD2` is the length in bp of the longest IBD2 segment the reference called for a pair.
It is a far sharper target than `.seg`'s totals: a candidate rule either reproduces the
number exactly or it does not, and every pair with kinship >= 2^-3.5 supplies one.

    python3 fit_ibd2.py <dir-with-ibs_*.ibs files>
"""

import sys

import numpy as np

import kingdata as kd
import maxibd2 as M
import rules2 as R

WORD = 64


def longest_ibd2(ds, i, j, t2, end_mode, start_mode, bridge, minw):
    """Longest IBD2 segment in bp under one candidate rule."""
    ibs0, n0, n1 = R.counts(ds, i, j)
    pos = ds.pos
    best = 0
    for _, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        ok = (n0[w0:w1 + 1] == 0) & (n1[w0:w1 + 1] <= t2)
        if bridge:
            ok = R._bridge(ok)
        for a, b in R._runs(ok):
            if b - a + 1 < minw:
                continue
            u, v = w0 + a, w0 + b
            # left
            if u == w0:
                lo_m = lo
            else:
                m = int(ibs0[u - 1])
                lo_m = WORD * (u - 1) + (63 - R._clz(m)) + 1 if m else WORD * u
                if start_mode == "back":
                    lo_m = WORD * (u - 1) if not m else lo_m
            # right
            if v >= w1:
                hi_m = hi
            else:
                m = int(ibs0[v + 1])
                if m:
                    hi_m = WORD * (v + 1) + (63 - R._clz(m))
                else:
                    hi_m = WORD * (v + 1) + 63 if end_mode == "next" else WORD * v + 63
                hi_m = min(hi_m, hi)
            best = max(best, int(pos[hi_m] - pos[lo_m]))
    return best


def main():
    dirname = sys.argv[1]
    targets = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        for (i, j), t in M.collect(dirname, ds).items():
            targets.append((ds, i, j, t))
    print(f"{len(targets)} pairs carry a MaxIBD2")
    print(f"{'rule':<44}{'exact':>7}{'mean |rel err|':>16}")
    for end_mode in ("next", "same"):
        for start_mode in ("run", "back"):
            for t2 in (0, 2, 3, 4, 5, 6):
                for bridge in (0, 1):
                    ok = 0
                    err = 0.0
                    for ds, i, j, t in targets:
                        g = longest_ibd2(ds, i, j, t2, end_mode, start_mode, bridge, 1)
                        ok += g == t
                        err += abs(g - t) / t
                    tag = f"end={end_mode} start={start_mode} t2={t2} bridge={bridge}"
                    print(f"{tag:<44}{ok:7d}{err / len(targets):16.4f}")


if __name__ == "__main__":
    main()
