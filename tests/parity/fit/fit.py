"""Score a candidate rule against every captured reference .seg row.

    python3 tests/parity/fit/fit.py baseline
    python3 tests/parity/fit/fit.py sweep <knob> <values...>

Scoring is on the *printed* values: a row counts as exact only when IBD1Seg and IBD2Seg
both agree at four decimals, which is what byte parity needs.
"""

import sys
import numpy as np

import kingdata as kd
import rules as R


def score(p=R.Params(), datasets=None, verbose=False, collect=None):
    datasets = datasets or kd.DATASETS
    tot = dict(exact=0, ref_rows=0, got_rows=0, missing=0, extra=0, both=0,
               inftype_ok=0, abs_err=0.0, worst=0.0, exact1=0, exact2=0)
    for name in datasets:
        ds = kd.load(name)
        for (i, j) in ds.pairs():
            ibd1, ibd2, longest = R.call_pair(ds, i, j, p)
            rep = R.reported(ibd1, ibd2, longest, ds.denom, p)
            ref = ds.ref.get((i, j))
            if ref is not None:
                tot["ref_rows"] += 1
            if rep:
                tot["got_rows"] += 1
            if ref is None and rep:
                tot["extra"] += 1
                if collect is not None:
                    collect.append((name, i, j, "extra", ibd1, ibd2, longest, None))
                continue
            if ref is not None and not rep:
                tot["missing"] += 1
                if collect is not None:
                    collect.append((name, i, j, "missing", ibd1, ibd2, longest, ref))
                continue
            if ref is None:
                continue
            tot["both"] += 1
            pi1 = ibd1 / ds.denom
            pi2 = ibd2 / ds.denom
            prop = pi2 + pi1 / 2.0
            g1, g2 = kd.fmt4(pi1), kd.fmt4(pi2)
            tot["exact1"] += g1 == ref[0]
            tot["exact2"] += g2 == ref[1]
            ok = g1 == ref[0] and g2 == ref[1]
            tot["exact"] += ok
            tot["inftype_ok"] += kd.inf_type(pi1, pi2, prop) == ref[3]
            e = abs(prop - ref[2])
            tot["abs_err"] += e
            tot["worst"] = max(tot["worst"], e)
            if collect is not None and not ok:
                collect.append((name, i, j, "value", ibd1, ibd2, longest, ref))
    tot["mae"] = tot["abs_err"] / max(1, tot["both"])
    if verbose:
        report(tot)
    return tot


def report(t):
    print(f"  rows ref={t['ref_rows']} got={t['got_rows']} "
          f"missing={t['missing']} extra={t['extra']}")
    print(f"  exact(both cols)={t['exact']}/{t['ref_rows']}  "
          f"IBD1 ok={t['exact1']}  IBD2 ok={t['exact2']}  "
          f"InfType ok={t['inftype_ok']}")
    print(f"  MAE(PropIBD)={t['mae']:.5f}  worst={t['worst']:.4f}")


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "baseline"
    if what == "baseline":
        print("baseline (mirrors crates/king-core/src/ibdseg.rs):")
        score(verbose=True)
    elif what == "sweep":
        knob = sys.argv[2]
        for v in sys.argv[3:]:
            try:
                val = int(v)
            except ValueError:
                try:
                    val = float(v)
                except ValueError:
                    val = v
            p = R.Params(**{knob: val})
            t = score(p)
            print(f"{knob}={v!r:>12}  exact={t['exact']:4d}  extra={t['extra']:4d}  "
                  f"missing={t['missing']:4d}  mae={t['mae']:.5f}")
    else:
        raise SystemExit("usage: fit.py [baseline|sweep knob v1 v2 ...]")


if __name__ == "__main__":
    main()
