"""Final scorecard for the informativeness gate, including the held-out --seglength runs.

    python3 gate_final.py           # scorecard at the default --seglength 3
    python3 gate_final.py sweeplen  # the same rule against the seglength 5 / 10 captures
    python3 gate_final.py knobs     # re-search the remaining knobs with the gate in place
"""

import sys

import kingdata as kd
import rules3 as R3
from gate_corpus import score, show, BASE, POISON

FIT = dict(BASE)
FIT["min1"] = 1                      # the gate replaces the "2 clean words" rule
BEST = R3.P3(**FIT, gate=10, gate_scope="core")


def _score_suffix(p, suffix, seglength_bp):
    """Score against a capture taken with a different --seglength."""
    tot = dict(rows=0, exact=0, extra=0, missing=0)
    for name in kd.DATASETS:
        ds = kd.load(name)
        ref = ds._read_seg(suffix)
        q = R3.P3(**{**FIT, "seglength_bp": seglength_bp}, gate=p.gate,
                  gate_scope=p.gate_scope)
        for (i, j) in ds.pairs():
            a, b, longest = R3.call_pair(ds, i, j, q)
            got = longest >= q.long_bp
            r = ref.get((i, j))
            if not got and r is None:
                continue
            if got and r is None:
                tot["extra"] += 1
                continue
            tot["rows"] += 1
            if not got:
                tot["missing"] += 1
                continue
            i1, i2 = kd.fmt4(a / ds.denom), kd.fmt4(b / ds.denom)
            pr = kd.fmt4(b / ds.denom + a / ds.denom / 2)
            tot["exact"] += (i1 == r[0] and i2 == r[1] and pr == r[2])
    return tot


def main():
    per, tot = score(R3.P3(**BASE, gate=0), gate=False)
    show(per, tot, "committed engine (no gate, min run 2 words)")
    per, tot = score(BEST)
    show(per, tot, "gate C=10 over the run's own words, min run 1 word")


def sweeplen():
    print("held-out captures: the same rule at other --seglength values")
    print(f"{'capture':<28}{'rows':>6}{'exact':>7}{'extra':>7}{'miss':>6}")
    for suffix, bp in (("__ibdseg", 3_000_000), ("__ibdseg_seglength5", 5_000_000),
                       ("__ibdseg_seglength10", 10_000_000)):
        t = _score_suffix(BEST, suffix, bp)
        print(f"{suffix:<28}{t['rows']:>6}{t['exact']:>7}{t['extra']:>7}{t['missing']:>6}")


def knobs():
    print("remaining knobs, with the gate fixed at C=10 / core / min1=1")
    print(f"{'knob':<34}{'exact':>7}{'IBD1':>7}{'IBD2':>7}{'extra':>7}{'miss':>6}{'MAE':>10}")
    cands = [("as fitted", {}),
             ("min2=2", dict(min2=2)),
             ("t2=1", dict(t2=1)),
             ("t2=4", dict(t2=4)),
             ("end2=same (IBD2 clipped)", dict(end2="same")),
             ("edge=edge", dict(edge="edge")),
             ("edge=clamp", dict(edge="clamp")),
             ("bridge1=1", dict(bridge1=1))]
    for label, kw in cands:
        p = R3.P3(**{**FIT, **kw}, gate=10, gate_scope="core")
        _, t = score(p)
        print(f"{label:<34}{t['exact']:>7}{t['ibd1']:>7}{t['ibd2']:>7}"
              f"{t['extra']:>7}{t['missing']:>6}{t['mae']:>10.5f}")


if __name__ == "__main__":
    globals()[sys.argv[1] if len(sys.argv) > 1 else "main"]()
