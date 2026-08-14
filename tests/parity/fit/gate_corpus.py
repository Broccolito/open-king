"""Score the informativeness gate against all 982 captured `.seg` rows.

    python3 gate_corpus.py            # scorecard for the fitted gate
    python3 gate_corpus.py sweep      # sweep the constant C and the counting scope
"""

import sys

import numpy as np

import kingdata as kd
import rules2 as R2
import rules3 as R3

BASE = dict(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")
POISON = {"nuclear", "missing", "monomorphic"}


def score(p, datasets=kd.DATASETS, gate=True):
    tot = dict(rows=0, exact=0, ibd1=0, ibd2=0, extra=0, missing=0, mae=0.0, worst=0.0)
    per = {}
    for name in datasets:
        ds = kd.load(name)
        s = dict(rows=0, exact=0, ibd1=0, ibd2=0, extra=0, missing=0, mae=0.0, worst=0.0)
        for (i, j) in ds.pairs():
            if gate:
                a, b, longest = R3.call_pair(ds, i, j, p)
            else:
                a, b, longest = R2.call_pair(ds, i, j, p)
            got = longest >= p.long_bp
            ref = ds.ref.get((i, j))
            if not got and ref is None:
                continue
            if got and ref is None:
                s["extra"] += 1
                continue
            if not got and ref is not None:
                s["missing"] += 1
                s["rows"] += 1
                s["mae"] += abs(ref[2])
                s["worst"] = max(s["worst"], abs(ref[2]))
                continue
            s["rows"] += 1
            i1 = kd.fmt4(a / ds.denom)
            i2 = kd.fmt4(b / ds.denom)
            pr = kd.fmt4(b / ds.denom + a / ds.denom / 2)
            s["ibd1"] += i1 == ref[0]
            s["ibd2"] += i2 == ref[1]
            s["exact"] += (i1 == ref[0] and i2 == ref[1] and pr == ref[2])
            d = abs((b / ds.denom + a / ds.denom / 2) - ref[2])
            s["mae"] += d
            s["worst"] = max(s["worst"], d)
        s["mae"] /= max(s["rows"], 1)
        per[name] = s
        for k in tot:
            if k == "mae":
                tot[k] += s[k] * s["rows"]
            elif k == "worst":
                tot[k] = max(tot[k], s[k])
            else:
                tot[k] += s[k]
    tot["mae"] /= max(tot["rows"], 1)
    return per, tot


def show(per, tot, title):
    print(f"\n=== {title} ===")
    print(f"{'dataset':<14}{'rows':>6}{'exact':>7}{'IBD1':>7}{'IBD2':>7}"
          f"{'extra':>7}{'miss':>6}{'MAE':>10}{'worst':>9}")
    for name, s in per.items():
        mark = "  (poisoned)" if name in POISON else ""
        print(f"{name:<14}{s['rows']:>6}{s['exact']:>7}{s['ibd1']:>7}{s['ibd2']:>7}"
              f"{s['extra']:>7}{s['missing']:>6}{s['mae']:>10.5f}{s['worst']:>9.4f}{mark}")
    print(f"{'ALL':<14}{tot['rows']:>6}{tot['exact']:>7}{tot['ibd1']:>7}{tot['ibd2']:>7}"
          f"{tot['extra']:>7}{tot['missing']:>6}{tot['mae']:>10.5f}{tot['worst']:>9.4f}")
    clean = [k for k in per if k not in POISON]
    c = {k: sum(per[n][k] for n in clean) for k in ("rows", "exact", "extra", "missing")}
    print(f"{'CLEAN':<14}{c['rows']:>6}{c['exact']:>7}{'':>7}{'':>7}"
          f"{c['extra']:>7}{c['missing']:>6}   (bigish/admixed/multifam/threegen/dups"
          f"/sexchr/unrelated)")


def main():
    per, tot = score(R3.P3(**BASE, gate=0), gate=False)
    show(per, tot, "baseline: no gate, min run 2 words (the committed rule)")
    per, tot = score(R3.P3(**BASE, gate=10, gate_scope="core"))
    show(per, tot, "gate C=10 over the run's own words, min run still 2 words")
    per, tot = score(R3.P3(**{**BASE, "min1": 1}, gate=10, gate_scope="core"))
    show(per, tot, "gate C=10, min run 1 word  <- the fitted rule")


def sweep():
    print("constant C x counting scope, over all 982 rows")
    print(f"{'C':>4} {'scope':>6} {'exact':>7} {'extra':>7} {'missing':>8} {'MAE':>10}")
    for scope in ("core", "ext"):
        for C in (0, 4, 6, 8, 9, 10, 11, 12, 14, 16, 20):
            _, tot = score(R3.P3(**BASE, gate=C, gate_scope=scope))
            print(f"{C:>4} {scope:>6} {tot['exact']:>7} {tot['extra']:>7} "
                  f"{tot['missing']:>8} {tot['mae']:>10.5f}")


if __name__ == "__main__":
    (sweep if len(sys.argv) > 1 and sys.argv[1] == "sweep" else main)()
