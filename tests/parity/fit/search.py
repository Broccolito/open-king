"""Score / search the word-aligned caller (`rules2.P`) against every captured .seg row.

    python3 search.py score                  # the current default knob-set
    python3 search.py sweep t2 0 2 4 6       # one knob
    python3 search.py grid                   # the full cross-product
"""

import itertools
import sys
from dataclasses import replace

import kingdata as kd
import rules2 as R


def score(p=R.P(), datasets=None, collect=None):
    datasets = datasets or kd.DATASETS
    t = dict(ref=0, got=0, exact=0, exact1=0, exact2=0, missing=0, extra=0,
             both=0, inf=0, err=0.0, worst=0.0)
    for name in datasets:
        ds = kd.load(name)
        for (i, j) in ds.pairs():
            a, b, longest = R.call_pair(ds, i, j, p)
            rep = longest > p.long_bp
            ref = ds.ref.get((i, j))
            t["ref"] += ref is not None
            t["got"] += rep
            if ref is None:
                t["extra"] += rep
                if collect is not None and rep:
                    collect.append((name, i, j, "extra", a, b, longest, None))
                continue
            if not rep:
                t["missing"] += 1
                if collect is not None:
                    collect.append((name, i, j, "missing", a, b, longest, ref))
                continue
            t["both"] += 1
            pi1, pi2 = a / ds.denom, b / ds.denom
            prop = pi2 + pi1 / 2
            g1, g2 = kd.fmt4(pi1), kd.fmt4(pi2)
            t["exact1"] += g1 == ref[0]
            t["exact2"] += g2 == ref[1]
            ok = g1 == ref[0] and g2 == ref[1]
            t["exact"] += ok
            t["inf"] += kd.inf_type(pi1, pi2, prop) == ref[3]
            e = abs(prop - ref[2])
            t["err"] += e
            t["worst"] = max(t["worst"], e)
            if collect is not None and not ok:
                collect.append((name, i, j, "value", a, b, longest, ref))
    t["mae"] = t["err"] / max(1, t["both"])
    return t


def line(tag, t):
    return (f"{tag:<34} exact={t['exact']:4d}/{t['ref']:<4d} "
            f"IBD1={t['exact1']:4d} IBD2={t['exact2']:4d} "
            f"extra={t['extra']:4d} miss={t['missing']:4d} "
            f"inf={t['inf']:4d} mae={t['mae']:.5f} worst={t['worst']:.4f}")


def main():
    what = sys.argv[1] if len(sys.argv) > 1 else "score"
    if what == "score":
        print(line("default", score()))
    elif what == "sweep":
        knob = sys.argv[2]
        for v in sys.argv[3:]:
            p = replace(R.P(), **{knob: int(v)})
            print(line(f"{knob}={v}", score(p)))
    elif what == "grid":
        best = []
        space = dict(t1=[0, 1], t2=[0, 2, 4], min1=[1, 2, 3], min2=[1, 2],
                     bridge1=[0, 1], bridge2=[0, 1],
                     edge=["edge", "fringe"], end2=["next", "same"])
        keys = list(space)
        for combo in itertools.product(*(space[k] for k in keys)):
            p = R.P(**dict(zip(keys, combo)))
            t = score(p)
            best.append((t["exact"], -t["extra"], dict(zip(keys, combo)), t))
        best.sort(reverse=True, key=lambda x: (x[0], x[1]))
        for e, ne, cfg, t in best[:15]:
            print(line(str(cfg), t))
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
