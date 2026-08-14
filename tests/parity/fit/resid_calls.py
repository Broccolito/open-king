"""Localise a per-call disagreement to its markers.

`resid19.py` reports, per pair, which reference IBD2 calls under 10 Mb we fail to
reproduce.  This inverts each such length back to the marker interval(s) that could have
produced it — `seglen_invert.py`'s trick — and prints them beside our own calls over the
same usable segment, with the per-word predicate spelled out.

    python3 resid_calls.py multifam A_C2 C_F
    python3 resid_calls.py                       # every disagreeing pair
"""

import json
import os
import sys

import numpy as np

import engine as E
import kingdata as kd
import resid19 as R

WORD = 64
TOL = 2


def invert(ds, target, tol=TOL):
    out = []
    for chrom, lo, hi in ds.segs:
        p = ds.pos[lo:hi + 1]
        k = np.searchsorted(p, p + target - tol)
        for a in range(len(p)):
            for b in range(int(k[a]), min(int(k[a]) + 2, len(p))):
                if abs(int(p[b]) - int(p[a]) - target) <= tol:
                    out.append((lo + a, lo + b, (chrom, lo, hi)))
    return out


def wordline(sc, ds, i, j, a, b):
    """The per-word predicate string over words [a, b] of the usable segment."""
    w0 = sc.w0
    s = []
    for k in range(max(0, a), min(sc.n, b + 1)):
        z = int(sc.n0[w0 + k])
        m = int(sc.n1[w0 + k])
        inf = int(sc.cum2s[w0 + k + 1] - sc.cum2s[w0 + k])
        tag = "W" if z else ("C" if m == 0 else ("x" if m == 1 else "y"))
        s.append("%s%d[%d/%d/%d]" % (tag, w0 + k, z, m, inf))
    return " ".join(s)


def show_pair(name, n1, n2):
    ds = kd.load(name)
    idx = {f[1]: k for k, f in enumerate(ds.fam)}
    i, j = sorted((idx[n1], idx[n2]))
    pos = ds.pos
    print("=== %s %s/%s" % (name, n1, n2))
    ours = R.calls_of(ds, i, j)
    path = os.path.join(R.SEGLEN, "%s.IBD2Seg.json" % name)
    want = json.load(open(path)).get("%d,%d" % (i, j), [])
    gshort = sorted(int(pos[b] - pos[a]) for _s, a, b in ours if
                    int(pos[b] - pos[a]) < 10_000_000)
    matched, gleft, wleft = R._match(gshort, sorted(want), TOL)
    print("  reference short calls: %s" % [round(v / 1e6, 3) for v in sorted(want)])
    print("  our short calls:       %s" % [round(v / 1e6, 3) for v in gshort])
    print("  matched %d  ours-only %s  ref-only %s"
          % (matched, [round(v / 1e6, 3) for v in gleft],
             [round(v / 1e6, 3) for v in wleft]))
    segs = {}
    for t in wleft:
        for a, b, seg in invert(ds, t):
            segs.setdefault(seg, []).append(("REF", t, a, b))
    for _s, a, b in ours:
        ln = int(pos[b] - pos[a])
        if ln in gleft or any(abs(ln - v) <= TOL for v in gleft):
            segs.setdefault(_s, []).append(("OURS*", ln, a, b))
    for seg, items in sorted(segs.items()):
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        print("  -- usable segment chr%d markers %d..%d words %d..%d"
              % (seg[0], seg[1], seg[2], sc.w0, sc.w1))
        allours = [(a, b) for s, a, b in ours if s == seg]
        print("     our calls: %s"
              % ["[%d,%d] w%d..%d %.3fMb" % (a, b, a // WORD, b // WORD,
                                             (pos[b] - pos[a]) / 1e6)
                 for a, b in allours])
        for tag, t, a, b in sorted(items, key=lambda v: v[2]):
            print("     %-5s %8.3f Mb  [%d,%d]  w%d+%d .. w%d+%d"
                  % (tag, t / 1e6, a, b, a // WORD, a % WORD, b // WORD, b % WORD))
            print("        %s" % wordline(sc, ds, i, j, a // WORD - 2 - sc.w0,
                                          b // WORD + 2 - sc.w0))


PAIRS = [("nuclear", "N_C1", "N_C4"), ("threegen", "TG_C3", "TG_C4"),
         ("multifam", "C_M", "D_M"), ("multifam", "D_C1", "D_C2"),
         ("multifam", "D_C2", "D_C3"), ("multifam", "A_C1", "A_C3"),
         ("multifam", "A_C2", "C_F"), ("multifam", "A_C3", "C_F"),
         ("multifam", "B_C1", "B_C2"), ("dups", "MZ_1", "MZ_2"),
         ("missing", "M_C1", "M_C3"), ("missing", "M_C2", "M_C3"),
         ("admixed", "X_C1", "X_C2"), ("admixed", "Z_C1", "Z_C2")]


if __name__ == "__main__":
    if len(sys.argv) > 3:
        show_pair(sys.argv[1], sys.argv[2], sys.argv[3])
    else:
        for t in PAIRS:
            show_pair(*t)
            print()
