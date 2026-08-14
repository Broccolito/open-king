"""Profile the wrong `.seg` rows at `--seglength` 5 and 10 Mb — `20-seglength-floor.md`.

`seg19.py` scores 982/982 on both estimate columns at the default 3 Mb floor and 910/946
(5 Mb) / 844/937 (10 Mb) above it.  This says *which* rows, on *which* column, with what
sign, and how big the miss is in units of the floor — so the shape of the missing clause
is read off the residual instead of guessed.

    python3 where20.py            # per-dataset summary at 5 and 10 Mb
    python3 where20.py rows       # every wrong row
"""

import sys
from collections import Counter

import engine as E
import kingdata as kd
import seg19 as S19

RULE = S19.R19()


def rows(min_bp, suffix, p=RULE, caller=None):
    """Every graded row: (ds, i, j, d1, d2, want1, want2, got1, got2)."""
    out = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            a, b, lg = (caller or S19.call_pair)(ds, i, j, p, min_bp)
            if (i, j) not in ref:
                if lg >= E.LONG:
                    out.append((name, i, j, "EXTRA", None, None, None, None, None))
                continue
            if lg < E.LONG:
                out.append((name, i, j, "MISS", None, None, None, None, None))
                continue
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = a / d, b / d
            out.append((name, i, j, kd.fmt4(g1) - a1, kd.fmt4(g2) - a2,
                        a1, a2, g1, g2))
    return out


def summary(min_bp, suffix, p=RULE):
    print("=== --seglength %g Mb ===" % (min_bp / 1e6))
    rr = rows(min_bp, suffix, p)
    per = {}
    for name, i, j, d1, d2, w1, w2, g1, g2 in rr:
        s = per.setdefault(name, Counter())
        s["rows"] += 1
        if isinstance(d1, str):
            s[d1] += 1
            continue
        if abs(d1) > 5e-5:
            s["bad1"] += 1
            s["bad1+" if d1 > 0 else "bad1-"] += 1
        if abs(d2) > 5e-5:
            s["bad2"] += 1
            s["bad2+" if d2 > 0 else "bad2-"] += 1
    tot = Counter()
    for name in kd.DATASETS:
        s = per.get(name, Counter())
        tot.update(s)
        if s["bad1"] or s["bad2"]:
            print("  %-12s rows %4d  ibd1 wrong %3d (+%d/-%d)  ibd2 wrong %3d (+%d/-%d)"
                  % (name, s["rows"], s["bad1"], s["bad1+"], s["bad1-"],
                     s["bad2"], s["bad2+"], s["bad2-"]))
    print("  %-12s rows %4d  ibd1 wrong %3d (+%d/-%d)  ibd2 wrong %3d (+%d/-%d)"
          % ("TOTAL", tot["rows"], tot["bad1"], tot["bad1+"], tot["bad1-"],
             tot["bad2"], tot["bad2+"], tot["bad2-"]))
    return rr


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    for bp, sfx in S19.FLOORS:
        rr = summary(bp, sfx)
        if mode == "rows":
            for name, i, j, d1, d2, w1, w2, g1, g2 in rr:
                if isinstance(d1, str) or abs(d1) > 5e-5 or abs(d2) > 5e-5:
                    if isinstance(d1, str):
                        print("    %-11s %3d %3d  %s" % (name, i, j, d1))
                    else:
                        print("    %-11s %3d %3d  ibd1 %+.4f (want %.4f got %.4f)"
                              "   ibd2 %+.4f (want %.4f got %.4f)"
                              % (name, i, j, d1, w1, g1, d2, w2, g2))
