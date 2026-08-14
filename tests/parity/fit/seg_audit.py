"""Per-**segment** audit: every pair, every usable segment, ours against the reference.

The `.seg` row is a sum over 5-22 usable segments, so a row can be right by cancellation.
This grades the summands instead, using `seglocal.py` to make the reference print one
segment's IBD2 total at a time — about 4 000 reference invocations over the nine smaller
datasets, and the strongest statement the corpus can make about the caller's geometry.

    python3 seg_audit.py            # every dataset but `bigish`
    python3 seg_audit.py dups       # one dataset
"""

import sys
from collections import Counter

import engine as E
import kingdata as kd
import resid19 as R
import seglocal as SL

SMALL = [n for n in kd.DATASETS if n != "bigish"]


def audit(names):
    tot = Counter()
    worst = []
    for name in names:
        ds = kd.load(name)
        ulp = ds.denom / 10000.0
        for (i, j) in sorted(ds.ref):
            for k, seg in enumerate(ds.segs):
                c, ours = SL.ours(ds, i, j, seg)
                val, dg = SL.run(ds, i, j, [(seg[1], seg[2])], None,
                                 SL.pick_carrier(ds, k))
                if dg != SL.base_digest(ds):
                    tot["allsegs moved"] += 1
                    continue
                if val is None:
                    tot["pair not reported"] += 1
                    continue
                d = ours - float(val[1]) * ds.denom
                tot["segments"] += 1
                tot["exact(<0.5ulp)"] += abs(d) <= 0.5 * ulp
                tot["within 0.25ulp"] += abs(d) <= 0.25 * ulp
                if abs(d) > 0.5 * ulp:
                    worst.append((abs(d) / ulp, name, ds.fam[i][1], ds.fam[j][1], k,
                                  seg, c, ours, val[1], d))
        print("  %-12s done (%d segments graded)" % (name, tot["segments"]))
    print("=== per-segment audit")
    for k in ("segments", "exact(<0.5ulp)", "within 0.25ulp", "pair not reported",
              "allsegs moved"):
        print("  %-20s %d" % (k, tot[k]))
    worst.sort(reverse=True)
    for w in worst[:40]:
        ds = kd.load(w[1])
        print("  %+6.2f ulp  %-10s %-8s %-8s seg%2d chr%-2d  ref %s  ours %d  d=%+.0f"
              % (w[9] / (ds.denom / 10000.0), w[1], w[2], w[3], w[4], w[5][0], w[8],
                 w[7], w[9]))


if __name__ == "__main__":
    audit(sys.argv[1:] or SMALL)
