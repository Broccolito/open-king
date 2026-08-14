"""Profile the residual left by `seg19.py`'s fringe rule, field by field.

`seg19.py` makes both printed columns exact on all 982 rows at the default floor, yet only
806 rows are byte-identical: the disagreement has moved into `PropIBD`, which is
`IBD2Seg + IBD1Seg/2` and therefore reads about one extra digit of the underlying base-pair
totals than either column does on its own.  This script grades every field separately and
reports where the surviving error lives.

    python3 prof19.py            # the field-by-field profile at 3 / 5 / 10 Mb
    python3 prof19.py rows       # ...plus every disagreeing row
"""

import sys
from collections import Counter

import engine as E
import kingdata as kd
import seg19 as S19

FLOORS = [(3_000_000, "__ibdseg"), (5_000_000, "__ibdseg_seglength5"),
          (10_000_000, "__ibdseg_seglength10")]

RULE = S19.R19()


def rows_of(min_bp, suffix, rule=RULE):
    out = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            if (i, j) not in ref:
                continue
            a, b, lg = S19.call_pair(ds, i, j, rule, min_bp)
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = a / d, b / d
            gp = g2 + g1 / 2
            out.append(dict(
                ds=ds.name, i=i, j=j, bp1=a, bp2=b, denom=d,
                g1=g1, g2=g2, gp=gp, r1=a1, r2=a2, rp=ap, rt=at,
                ok1=kd.fmt4(g1) == a1, ok2=kd.fmt4(g2) == a2,
                okp=kd.fmt4(gp) == ap, okt=kd.inf_type(g1, g2, gp) == at,
                d1=a - a1 * d, d2=b - a2 * d, dp=gp - ap,
            ))
    return out


def profile(min_bp, suffix):
    rs = rows_of(min_bp, suffix)
    n = len(rs)
    bad = [r for r in rs if not (r["ok1"] and r["ok2"] and r["okp"] and r["okt"])]
    print("=== --seglength %d Mb: %d rows, %d not byte-identical"
          % (min_bp // 1_000_000, n, len(bad)))
    print("  wrong IBD1Seg %d   wrong IBD2Seg %d   wrong PropIBD %d   wrong InfType %d"
          % (sum(not r["ok1"] for r in rs), sum(not r["ok2"] for r in rs),
             sum(not r["okp"] for r in rs), sum(not r["okt"] for r in rs)))
    only_p = [r for r in bad if r["ok1"] and r["ok2"]]
    print("  rows wrong ONLY on PropIBD/InfType (both columns exact): %d" % len(only_p))
    ulp = 1e-4
    sgn = Counter("+" if r["dp"] > 0 else "-" for r in bad if not r["okp"])
    print("  PropIBD sign: over %d  under %d" % (sgn["+"], sgn["-"]))
    h = Counter(round(abs(r["dp"]) / ulp, 1) for r in bad if not r["okp"])
    print("  |PropIBD delta| in ulps: " + "  ".join("%s:%d" % kv
                                                    for kv in sorted(h.items())[:12]))
    print("  by dataset: " + "  ".join(
        "%s %d/%d" % (k, sum(1 for r in bad if r["ds"] == k),
                      sum(1 for r in rs if r["ds"] == k)) for k in kd.DATASETS))
    print("  by InfType: " + "  ".join(
        "%s %d/%d" % (k, sum(1 for r in bad if r["rt"] == k),
                      sum(1 for r in rs if r["rt"] == k))
        for k in sorted({r["rt"] for r in rs})))
    # Is the reference's own PropIBD consistent with its own printed columns?
    cons = sum(1 for r in rs if kd.fmt4(r["r2"] + r["r1"] / 2) == r["rp"])
    print("  reference rows where fmt4(r2 + r1/2) == printed PropIBD: %d / %d"
          % (cons, n))
    # How close is our PropIBD to the half-ulp boundary on the rows that pass?
    near = sum(1 for r in rs if r["okp"]
               and abs(r["gp"] * 10000 - round(r["gp"] * 10000)) > 0.45)
    print("  passing rows within 0.05 ulp of a rounding boundary: %d" % near)
    return rs, bad


if __name__ == "__main__":
    verbose = len(sys.argv) > 1 and sys.argv[1] == "rows"
    for bp, sfx in FLOORS:
        rs, bad = profile(bp, sfx)
        if verbose:
            for r in sorted(bad, key=lambda r: -abs(r["dp"])):
                ds = kd.load(r["ds"])
                print("    %-12s %-9s %-9s  ibd1 %s%.4f/%.4f  ibd2 %s%.4f/%.4f  "
                      "prop %s%.6f/%.4f  type %s%s/%s"
                      % (r["ds"], ds.fam[r["i"]][1], ds.fam[r["j"]][1],
                         " " if r["ok1"] else "*", r["g1"], r["r1"],
                         " " if r["ok2"] else "*", r["g2"], r["r2"],
                         " " if r["okp"] else "*", r["gp"], r["rp"],
                         " " if r["okt"] else "*",
                         kd.inf_type(r["g1"], r["g2"], r["gp"]), r["rt"]))
        print()
