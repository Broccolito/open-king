"""**SUPERSEDED — read this first.** The premise below is wrong, and the way it is wrong is
worth keeping. It assumes `.seg`'s `PropIBD` is a third *rounding* of the same two base-pair
totals, so a disagreement there must mean our totals are slightly off. They are not:
`docs/research/20-seg-writer.md` shows `.seg` computes `PropIBD` from the **printed**
columns (`i2*1e-4 + i1*5e-5`), so it carries no information about the totals at all and the
"hexagon" this module derives is empty of signal. Every correction it reports was an artefact
of comparing two different formulas. Kept because `19-ibd2seg-residual.md` quotes its
numbers and because it is a clean example of a well-built instrument pointed at the wrong
question.

Read the reference's `.seg` totals to **finer than the printed ulp**, using PropIBD.

`prof19.py` says the whole residual at the default floor is `PropIBD`: both estimate
columns are exact on all 982 rows and 176 rows still differ, always by 0.5-1.2 ulp of
the printed `%.4lf`.  That is only possible because `PropIBD = IBD2Seg + IBD1Seg/2` is a
*third* rounding of the same two base-pair totals, and it reads about one extra bit of
them than either column does alone.

This module turns that into an instrument.  Each printed column is an interval:

    ibd1 in [ (r1 - h) D, (r1 + h) D )      h = 0.00005
    ibd2 in [ (r2 - h) D, (r2 + h) D )
    ibd2 + ibd1/2 in [ (rp - h) D, (rp + h) D )

so the reference's true `(ibd1, ibd2)` lies in a **hexagon** — the intersection of a
square with a diagonal band a third as wide.  `feasible()` returns it, `deficit()`
returns the smallest correction that would move our own totals into it, and `profile()`
reports how those corrections are distributed.

The point is that the corrections are *directed*: a row can only be fixed by adding
base pairs, or only by removing them, and the amount is bounded on both sides.  That is
a much stronger statement about a wrong rule than "the row differs".

    python3 prop19.py            # the decomposition, at 3 / 5 / 10 Mb
    python3 prop19.py rows       # ...plus every row whose correction is forced
    python3 prop19.py gaps       # corrections expressed in marker gaps
"""

import sys
from collections import Counter

import kingdata as kd
import seg19 as S19

H = 0.00005          # half an ulp of the printed %.4lf
FLOORS = [(3_000_000, "__ibdseg"), (5_000_000, "__ibdseg_seglength5"),
          (10_000_000, "__ibdseg_seglength10")]
RULE = S19.R19()


def feasible(r1, r2, rp, d):
    """The reference's `(ibd1, ibd2)` region, in base pairs, as three intervals.

    Returns `(lo1, hi1, lo2, hi2, lop, hip)` — the two axis-aligned bounds and the
    bound on `ibd2 + ibd1/2`.  All are half-open on the right.
    """
    return ((r1 - H) * d, (r1 + H) * d,
            (r2 - H) * d, (r2 + H) * d,
            (rp - H) * d, (rp + H) * d)


def deficit(a, b, box):
    """The signed correction our `(ibd1=a, ibd2=b)` needs, on the PropIBD axis only.

    Positive = we are short by that many base pairs of `ibd2 + ibd1/2`; negative = long.
    Zero when the row already satisfies the PropIBD constraint.
    """
    lo1, hi1, lo2, hi2, lop, hip = box
    p = b + a / 2.0
    if p < lop:
        return lop - p
    if p >= hip:
        return hip - p          # negative, and reaching it exactly is not enough
    return 0.0


def rows_of(min_bp, suffix, rule=RULE):
    out = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        ref = ds._read_seg(suffix)
        gap = ds.median_gap if hasattr(ds, "median_gap") else None
        for i, j in ds.pairs():
            if (i, j) not in ref:
                continue
            a, b, _lg = S19.call_pair(ds, i, j, rule, min_bp)
            r1, r2, rp, rt = ref[(i, j)]
            box = feasible(r1, r2, rp, d)
            out.append(dict(ds=ds.name, i=i, j=j, a=a, b=b, d=d, gap=gap,
                            r1=r1, r2=r2, rp=rp, rt=rt, box=box,
                            dp=deficit(a, b, box),
                            ok1=kd.fmt4(a / d) == r1, ok2=kd.fmt4(b / d) == r2,
                            okp=kd.fmt4(b / d + a / 2.0 / d) == rp,
                            zero2=(b == 0 and r2 == 0.0)))
    return out


def profile(min_bp, suffix, verbose=False):
    rs = rows_of(min_bp, suffix)
    bad = [r for r in rs if not r["okp"]]
    print("=== --seglength %d Mb: %d rows, %d wrong on PropIBD"
          % (min_bp // 1_000_000, len(rs), len(bad)))
    print("  of those, columns already exact: IBD1 %d   IBD2 %d   both %d"
          % (sum(r["ok1"] for r in bad), sum(r["ok2"] for r in bad),
             sum(r["ok1"] and r["ok2"] for r in bad)))

    # The decisive split: rows where BOTH we and the reference report no IBD2 at all.
    # There `PropIBD = IBD1Seg/2`, so the correction is unambiguously the IBD1 pass's.
    z = [r for r in bad if r["zero2"]]
    zall = [r for r in rs if r["zero2"]]
    print("  rows with IBD2 == 0 on both sides: %d of %d wrong  "
          "(their PropIBD is IBD1Seg/2 alone, so the fault is the IBD1 pass)"
          % (len(z), len(zall)))
    nz = [r for r in bad if not r["zero2"]]
    print("  rows carrying IBD2:                %d of %d wrong"
          % (len(nz), len(rs) - len(zall)))

    for tag, sub in (("IBD2 == 0 rows", z), ("IBD2 > 0 rows", nz), ("all wrong", bad)):
        if not sub:
            continue
        short = [r for r in sub if r["dp"] > 0]
        long_ = [r for r in sub if r["dp"] < 0]
        print("  %-16s short %3d  long %3d   |correction| bp: min %8.0f  "
              "median %8.0f  max %9.0f"
              % (tag, len(short), len(long_),
                 min(abs(r["dp"]) for r in sub),
                 sorted(abs(r["dp"]) for r in sub)[len(sub) // 2],
                 max(abs(r["dp"]) for r in sub)))

    # In marker gaps: is the correction a whole number of markers?
    g = Counter()
    for r in bad:
        ds = kd.load(r["ds"])
        mg = float(kd.np.median(kd.np.diff(ds.pos)))
        g[round(abs(r["dp"]) / mg * 2) / 2.0] += 1      # in half-gaps
    print("  |correction| in median marker gaps (of ibd2+ibd1/2): "
          + "  ".join("%.1f:%d" % kv for kv in sorted(g.items())[:14]))

    print("  by dataset: " + "  ".join(
        "%s %d/%d" % (k, sum(1 for r in bad if r["ds"] == k),
                      sum(1 for r in rs if r["ds"] == k)) for k in kd.DATASETS))
    if verbose:
        for r in sorted(bad, key=lambda r: -abs(r["dp"]))[:40]:
            ds = kd.load(r["ds"])
            print("    %-12s %-9s %-9s  ibd1 %10d  ibd2 %10d  need %+9.0f bp  "
                  "(%s)" % (r["ds"], ds.fam[r["i"]][1], ds.fam[r["j"]][1],
                            r["a"], r["b"], r["dp"],
                            "IBD1 only" if r["zero2"] else "mixed"))
    return rs, bad


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    for bp, sfx in FLOORS:
        profile(bp, sfx, verbose=(mode == "rows"))
        print()
