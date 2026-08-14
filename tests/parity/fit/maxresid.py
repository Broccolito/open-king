"""Residual analysis for the MaxIBD2 fit: what the misses have in common.

    python3 maxresid.py [t2]

For every pair the candidate rule gets wrong, this prints the interval the reference must
have used (from the inversion) next to the one the rule produced, plus the per-word
(IBS0/IBS1) profile around both, so the failure mode can be named rather than guessed.
"""

import sys
from collections import Counter

import numpy as np

import kingdata as kd
import maxfit as M
import maxdump as D

WORD = 64
RULE = dict(t2=4, end="next", start="aligned", edge="word", min_run=1)


def profile(ds, i, j, lo_w, hi_w, mark):
    _, n0, n1 = M.counts(ds, i, j)
    cells = []
    for w in range(max(0, lo_w), min(len(n0), hi_w + 1)):
        s = "".join(mark.get(w, ""))
        cells.append(f"{s}{int(n0[w])}/{int(n1[w])}")
    return " ".join(cells)


def main():
    t2 = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    rule = dict(RULE, t2=t2)
    tg = M.all_targets()
    kinds = Counter()
    nbad = 0
    for ds, i, j, t in tg:
        segs = M.ibd2_segments(ds, i, j, **rule)
        got = max((s[2] for s in segs), default=0)
        if got == t:
            continue
        nbad += 1
        cs = D.candidates(ds, i, j, t)
        best = max(segs, key=lambda s: s[2], default=None)
        # which of our segments is closest to the target?
        near = sorted(segs, key=lambda s: abs(s[2] - t))[:1]
        if len(cs) != 1:
            kinds["target not localised"] += 1
            print(f"\n{ds.name} {i},{j} target={t} got={got}  NOT LOCALISED "
                  f"({len(cs)} candidates)")
            if best:
                a, b = best[0] // WORD, (best[1] + 1) // WORD - 1
                print("   ours w%d..%d  %s" % (a, b, profile(ds, i, j, a - 2, b + 2,
                                                             {a: "[", b: "]"})))
            continue
        a, b, w0, w1 = cs[0]
        ours = None
        for lo_m, hi_m, ln in segs:
            if lo_m == WORD * a:
                ours = (lo_m // WORD, (hi_m + 1) // WORD - 1, ln)
        tag = []
        if ours is None:
            tag.append("start differs / run not found")
        else:
            if ours[1] != b:
                tag.append(f"end off by {ours[1] - b} words")
            if ours[2] > t and best and best[2] > ours[2]:
                tag.append("another segment is longer")
        if best and best[0] // WORD != a:
            tag.append("max is a different segment")
        k = "; ".join(tag) or "?"
        kinds[k] += 1
        mark = {WORD and a: "[", b: "]"}
        print(f"\n{ds.name} {i},{j} target={t} got={got}  [{k}]")
        print(f"   ref  w{a}..{b}   seg w{w0}..{w1}   {profile(ds, i, j, a - 2, b + 2, {a: '[', b: ']'})}")
        if ours:
            print(f"   ours w{ours[0]}..{ours[1]} len={ours[2]}")
        if best:
            print(f"   our max w{best[0] // WORD}..{(best[1] + 1) // WORD - 1} len={best[2]}")
    print(f"\n{nbad} misses of {len(tg)}")
    for k, v in kinds.most_common():
        print(f"  {v:4d}  {k}")


if __name__ == "__main__":
    main()
