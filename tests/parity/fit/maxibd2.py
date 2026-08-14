"""Solve the IBD2 boundary rule from `--ibs`'s MaxIBD2 column.

`--ibs` prints `MaxIBD2`, the **length in base pairs of one single segment** — the longest
IBD2 segment the same engine called — for every pair with kinship >= 2^-3.5 on a map with
>= 100 Mb usable.  That is far more information than `.seg`'s two aggregates: one exact
segment length per pair, ~1000 of them across the corpus.

This module inverts it: given the length, find every word-aligned marker interval of a
usable segment whose span equals it exactly, then report the per-word IBS0/IBS1 profile at
and around the two endpoints.  The offsets that repeat across a thousand pairs are the
boundary rule.

Usage:
    python3 maxibd2.py <dir-with-ibs_*.ibs files>
"""

import os
import sys
from collections import Counter

import numpy as np

import kingdata as kd

PC = np.bitwise_count


def read_maxibd2(path):
    """{(iid1, iid2): maxibd2_bp} from a .ibs / .ibs0 file."""
    out = {}
    with open(path) as fh:
        head = next(fh).rstrip("\n").split("\t")
        c = head.index("MaxIBD2")
        id1 = head.index("ID1")
        id2 = head.index("ID2")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[(f[id1], f[id2])] = int(float(f[c]))
    return out


def collect(dirname, ds):
    m = {}
    for suffix in (".ibs", ".ibs0"):
        p = os.path.join(dirname, f"ibs_{ds.name}{suffix}")
        if os.path.exists(p):
            m.update(read_maxibd2(p))
    idx = {iid: k for k, (fid, iid) in enumerate(ds.fam)}
    out = {}
    for (a, b), v in m.items():
        if v > 0:
            i, j = idx[a], idx[b]
            out[(min(i, j), max(i, j))] = v
    return out


def word_runs(ok):
    d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()))


def endpoints_for(ds, i, j, target):
    """Every word-aligned interval inside a usable segment whose span == target."""
    pos = ds.pos.astype(np.int64)
    hits = []
    for _, lo, hi in ds.segs:
        w0 = -(-lo // 64)
        w1 = (hi + 1) // 64 - 1
        if w1 < w0:
            continue
        starts = np.arange(w0, w1 + 1)
        sp = pos[64 * starts]
        ep = pos[64 * starts + 63]
        # for each start word u find v with pos[64v+63] - pos[64u] == target
        k = np.searchsorted(ep, sp + target)
        for a, v in zip(starts.tolist(), k.tolist()):
            if v < len(starts) and ep[v] - pos[64 * a] == target:
                hits.append((a, int(starts[v])))
    return hits


def profile(ds, i, j):
    ibs0, ibs1, _, _ = ds.masks(i, j)
    return (PC(ibs0), PC(ibs1))


def main():
    dirname = sys.argv[1]
    tol = int(sys.argv[2]) if len(sys.argv) > 2 else 2
    datasets = sys.argv[3:] or kd.DATASETS
    left = Counter()
    right = Counter()
    amb = 0
    tot = 0
    for name in datasets:
        ds = kd.load(name)
        m = collect(dirname, ds)
        for (i, j), target in sorted(m.items()):
            n0, n1 = profile(ds, i, j)
            hits = endpoints_for(ds, i, j, target)
            # the intended hit is the one whose interior is IBS0-free and IBS1-poor
            good = []
            for (u, v) in hits:
                interior = slice(u, v + 1)
                if n0[interior].sum() == 0 and (n1[interior] <= tol).mean() > 0.5:
                    good.append((u, v))
            tot += 1
            if len(good) != 1:
                amb += 1
                continue
            u, v = good[0]
            # where does the strict ibs1-free/ibs0-free run inside this interval start/end?
            ok = (n0 == 0) & (n1 <= tol)
            runs = [r for r in word_runs(ok) if r[0] <= v and r[1] >= u]
            if not runs:
                amb += 1
                continue
            r0 = min(r[0] for r in runs)
            r1 = max(r[1] for r in runs)
            left[u - r0] += 1
            right[v - r1] += 1
    print(f"pairs with MaxIBD2 {tot}, unambiguous {tot - amb}")
    print("start-word offset (segment start - run start):", sorted(left.items()))
    print("end-word offset   (segment end   - run end):  ", sorted(right.items()))


if __name__ == "__main__":
    main()
