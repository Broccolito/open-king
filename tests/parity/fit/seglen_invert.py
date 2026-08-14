"""Turn each `--seglength`-recovered segment length back into its marker endpoints.

`seglen_probe.py` gives the exact base-pair length of every called segment under 10 Mb.
Each length is inverted here over all marker pairs of the pair's own usable segments, and
the surviving candidates are printed with their offsets in the 64-marker word grid and
their relation to the usable segment's own ends — which is what says whether a `.seg`
IBD2 segment is word-aligned, fringe-extended, or neither.

    python3 seglen_invert.py [dataset ...] [--col IBD1Seg]
"""

import json
import os
import sys

import numpy as np

import kingdata as kd
import engine as E

WORD = 64
OUT = os.path.join(kd.ROOT, "tests", "parity", "work", "seglen")


TOL = 1
"""Base pairs of slack allowed when matching a probed length to a marker span.

`--seglength X` is converted to base pairs through `X * 1e6` in double, and a six-decimal
`X` is never exactly representable, so the comparison flips one base pair either side of
the integer the banner prints. Calibrating the probe against open-king's own binary — whose
filter is exactly `pos[hi] - pos[lo] >= min_bp` — recovers its segment lengths precisely,
while the reference's come back one base pair short; every length here is therefore only
determined to ±1.
"""


def intervals(ds, target, tol=TOL):
    """(a, b, usable segment) with `pos[b] - pos[a] ~= target`, inside one segment."""
    out = []
    for chrom, lo, hi in ds.segs:
        p = ds.pos[lo:hi + 1]
        k = np.searchsorted(p, p + target - tol)
        for a in range(len(p)):
            for b in range(int(k[a]), min(int(k[a]) + 2, len(p))):
                if abs(int(p[b]) - int(p[a]) - target) <= tol:
                    out.append((lo + a, lo + b, (chrom, lo, hi)))
    return out


def describe(ds, a, b, seg, ibs0):
    _chrom, lo, hi = seg
    w0, w1 = -(-lo // WORD), (hi + 1) // WORD - 1
    tags = []
    tags.append("a=seg.lo" if a == lo else ("a=64w%d" % (a // WORD) if a % WORD == 0
                                            else "a=w%d+%d" % (a // WORD, a % WORD)))
    tags.append("b=seg.hi" if b == hi else ("b=64w%d+63" % (b // WORD)
                                            if b % WORD == WORD - 1
                                            else "b=w%d+%d" % (b // WORD, b % WORD)))
    n = int(sum(int(x).bit_count() for x in ibs0[a // WORD:b // WORD + 1]))
    tags.append("words %d..%d of %d..%d" % (a // WORD, b // WORD, w0, w1))
    tags.append("IBS0 in those words=%d" % n)
    return "  ".join(tags)


def main():
    argv = sys.argv[1:]
    col = "IBD2Seg"
    if "--col" in argv:
        k = argv.index("--col")
        col = argv[k + 1]
        del argv[k:k + 2]
    for name in argv:
        path = os.path.join(OUT, "%s.%s.json" % (name, col))
        if not os.path.exists(path):
            continue
        ref = json.load(open(path))
        ds = kd.load(name)
        print("=== %s" % name)
        for key, lens in sorted(ref.items()):
            if not lens:
                continue
            i, j = map(int, key.split(","))
            ibs0, _n0, _n1, _c1, _c2 = E.masks(ds, i, j)
            print("  %s/%s" % (ds.fam[i][1], ds.fam[j][1]))
            for t in sorted(lens):
                cands = intervals(ds, t)
                print("    %10d bp -> %d interval(s)" % (t, len(cands)))
                for a, b, seg in cands:
                    print("        [%d,%d] %s" % (a, b, describe(ds, a, b, seg, ibs0)))


if __name__ == "__main__":
    main()
