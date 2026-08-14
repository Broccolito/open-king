"""Final scorecard for the fitted rule: per dataset, per flag variant, per residual class.

    python3 report.py
"""

import kingdata as kd
import rules2 as R
import search

BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")


def per_dataset(p):
    print(f"{'dataset':<12}{'rows':>5}{'exact':>7}{'IBD1':>6}{'IBD2':>6}"
          f"{'extra':>7}{'miss':>6}{'inf':>5}{'mae':>9}{'worst':>8}")
    for name in kd.DATASETS:
        t = search.score(p, datasets=[name])
        print(f"{name:<12}{t['ref']:5d}{t['exact']:7d}{t['exact1']:6d}{t['exact2']:6d}"
              f"{t['extra']:7d}{t['missing']:6d}{t['inf']:5d}{t['mae']:9.5f}{t['worst']:8.4f}")
    t = search.score(p)
    print(f"{'ALL':<12}{t['ref']:5d}{t['exact']:7d}{t['exact1']:6d}{t['exact2']:6d}"
          f"{t['extra']:7d}{t['missing']:6d}{t['inf']:5d}{t['mae']:9.5f}{t['worst']:8.4f}")


def seglength_variants(p):
    """The same rule against the --seglength 5 / 10 captures."""
    for suffix, bp in (("__ibdseg_seglength5", 5_000_000),
                       ("__ibdseg_seglength10", 10_000_000)):
        q = R.P(**{**p.__dict__, "seglength_bp": bp})
        ref_rows = exact = extra = miss = 0
        err = 0.0
        for name in kd.DATASETS:
            ds = kd.load(name)
            try:
                ref = ds._read_seg(suffix)
            except FileNotFoundError:
                continue
            for (i, j) in ds.pairs():
                a, b, longest = R.call_pair(ds, i, j, q)
                rep = longest > q.long_bp
                r = ref.get((i, j))
                ref_rows += r is not None
                if r is None:
                    extra += rep
                    continue
                if not rep:
                    miss += 1
                    continue
                pi1, pi2 = a / ds.denom, b / ds.denom
                exact += kd.fmt4(pi1) == r[0] and kd.fmt4(pi2) == r[1]
                err += abs(pi2 + pi1 / 2 - r[2])
        print(f"--seglength {bp // 10**6:<3} rows={ref_rows} exact={exact} "
              f"extra={extra} missing={miss} mae={err / max(1, ref_rows):.5f}")


if __name__ == "__main__":
    print("rule:", BEST)
    per_dataset(BEST)
    print()
    seglength_variants(BEST)
