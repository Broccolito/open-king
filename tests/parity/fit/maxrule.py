"""The fitted IBD2 caller, and its residuals against MaxIBD2.

Rule (all three parts read off the inversion in `maxdump.py` / `maxbreak.py`):

    dirty(w) = IBS0(w) > 0 or IBS1(w) >= DIRTY
    a run is a maximal stretch of non-dirty words, absorbing a *single* dirty word that
    has a non-dirty word on both sides (two in a row always break it)
    segment = [64u, 64(v+EXT)+63], clipped to the usable segment's complete-word grid

    python3 maxrule.py [dirty] [bridge] [ext]      # score + list every residual
"""

import sys
from collections import Counter

import numpy as np

import kingdata as kd
import maxfit as M
import maxdump as D

WORD = 64
DIRTY, BRIDGE, EXT = 5, 1, 1


def segments(ds, i, j, dirty_t=DIRTY, bridge=BRIDGE, ext=EXT, min_run=1, min_bp=0):
    _, n0, n1 = M.counts(ds, i, j)
    pos = ds.pos
    clean = (n0 == 0) & (n1 < dirty_t)
    out = []
    for _, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        ok = clean[w0:w1 + 1].copy()
        if bridge:
            d = np.diff(np.concatenate(([False], ~ok, [False])).astype(np.int8))
            for a, b in zip(np.flatnonzero(d == 1).tolist(),
                            (np.flatnonzero(d == -1) - 1).tolist()):
                if b - a + 1 <= bridge and a > 0 and b + 1 < ok.size:
                    ok[a:b + 1] = True
        d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
        for a, b in zip(np.flatnonzero(d == 1).tolist(),
                        (np.flatnonzero(d == -1) - 1).tolist()):
            if b - a + 1 < min_run:
                continue
            u, v = w0 + a, min(w0 + b + ext, w1)
            lo_m, hi_m = WORD * u, WORD * v + 63
            ln = int(pos[hi_m] - pos[lo_m])
            if ln >= min_bp:
                out.append((u, v, ln))
    return out


def longest(ds, i, j, **kw):
    return max((s[2] for s in segments(ds, i, j, **kw)), default=0)


def main():
    a = sys.argv[1:]
    kw = dict(dirty_t=int(a[0]) if a else DIRTY,
              bridge=int(a[1]) if len(a) > 1 else BRIDGE,
              ext=int(a[2]) if len(a) > 2 else EXT)
    tg = M.all_targets()
    ok = 0
    per = Counter()
    bad = []
    for ds, i, j, t in tg:
        g = longest(ds, i, j, **kw)
        if g == t:
            ok += 1
            per[ds.name] += 1
        else:
            bad.append((ds, i, j, t, g))
    print(f"MaxIBD2 exact: {ok} / {len(tg)}   {kw}")
    print("  by dataset:", dict(per))
    print(f"\n{len(bad)} residuals:")
    for ds, i, j, t, g in bad:
        cs = D.candidates(ds, i, j, t)
        _, n0, n1 = M.counts(ds, i, j)
        loc = ""
        if len(cs) == 1:
            aa, bb, w0, w1 = cs[0]
            mine = segments(ds, i, j, **kw)
            best = max(mine, key=lambda s: s[2], default=None)
            same = [s for s in mine if s[0] == aa]
            loc = (f"ref w{aa}..{bb} (seg w{w0}..{w1})  "
                   f"ours-at-same-start={same}  our-max=w{best[0]}..{best[1]}")
            ctx = " ".join(
                ("[" if w == aa else "") + f"{int(n1[w])}"
                + ("*" if int(n0[w]) else "") + ("]" if w == bb else "")
                for w in range(max(w0, aa - 4), min(w1, bb + 4) + 1))
            loc += "\n        " + ctx
        else:
            loc = f"NOT LOCALISED ({len(cs)} candidates)"
        print(f"  {ds.name:12s} {i:3d},{j:3d} target={t:>10d} got={g:>10d} "
              f"d={g - t:>+9d}\n        {loc}")


if __name__ == "__main__":
    main()
