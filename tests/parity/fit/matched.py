"""Search for ANY feature that separates the over-called pairs from reported ones.

`docs/research/11-segment-rule-fit.md` §7 leaves the 188 extra rows unexplained: matched
on core size, no summary statistic told the two groups apart.  This script rebuilds the
labelled set and brute-forces a decision stump over a much wider feature list, printing
the best achievable purity so the failure (or success) is quantified rather than asserted.

    python3 matched.py [dataset]
"""

import sys

import numpy as np

import kingdata as kd
import rules2 as R

BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")
PC = np.bitwise_count


def seg_of(ds, lo, hi):
    for k, (_, a, b) in enumerate(ds.segs):
        if a <= lo <= b:
            return k, a, b
    return -1, lo, hi


def feats(ds, i, j, lo, hi):
    ibs0, ibs1, hethet, called = ds.masks(i, j)
    k, slo, shi = seg_of(ds, lo, hi)
    w0, w1 = -(-slo // 64), (shi + 1) // 64 - 1
    u, v = (lo + 63) // 64, (hi + 1) // 64 - 1     # complete words inside the call
    pos = ds.pos
    core_lo, core_hi = 64 * u, min(64 * v + 63, shi)
    f = {
        "len": int(pos[hi] - pos[lo]),
        "core_len": int(pos[core_hi] - pos[core_lo]) if core_hi > core_lo else 0,
        "nsnp": hi - lo + 1,
        "words": v - u + 1,
        "left_ext": int(pos[core_lo] - pos[lo]),
        "right_ext": int(pos[hi] - pos[core_hi]),
        "at_seg_start": int(u <= w0),
        "at_seg_end": int(v >= w1),
        "seg_words": w1 - w0 + 1,
        "chr": ds.chr[lo],
        "u_mod2": u % 2,
        "u_mod4": u % 4,
        "u": u,
        "seg_idx": k,
        "dist_from_seg_start": u - w0,
        "dist_to_seg_end": w1 - v,
        "ibs0_prev": int(PC(ibs0[u - 1])) if u > 0 else -1,
        "ibs0_next": int(PC(ibs0[v + 1])) if v + 1 < ds.nwords else -1,
        "ibs1_core": int(PC(ibs1[u:v + 1]).sum()),
        "hethet_core": int(PC(hethet[u:v + 1]).sum()),
        "called_core": int(PC(called[u:v + 1]).sum()),
        "miss_core": 64 * (v - u + 1) - int(PC(called[u:v + 1]).sum()),
        "ibs1_prev": int(PC(ibs1[u - 1])) if u > 0 else -1,
        "ibs1_next": int(PC(ibs1[v + 1])) if v + 1 < ds.nwords else -1,
        "spacing": int(pos[hi] - pos[lo]) // max(1, hi - lo),
    }
    return f


def main(name="bigish"):
    ds = kd.load(name)
    rows, labels = [], []
    for (i, j) in ds.pairs():
        a, b, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        f = feats(ds, i, j, best[1], best[2])
        f["kind"] = best[0]
        f["nseg"] = len(detail)
        rows.append(f)
        labels.append(1 if (i, j) in ds.ref else 0)
    labels = np.array(labels)
    keep = np.array([r["words"] == 2 for r in rows])
    print(f"{name}: reported={labels.sum()} extra={(labels == 0).sum()}  "
          f"(2-word core: {labels[keep].sum()} reported / {(labels[keep] == 0).sum()} extra)")

    lab = labels[keep]
    sub = [r for r, k in zip(rows, keep) if k]
    print("\nquartiles over the 2-word-core group")
    print(f"{'feature':<22} {'group':<8} {'min':>10} {'q1':>10} {'med':>10} "
          f"{'q3':>10} {'max':>10}")
    best_split = []
    for key in sub[0]:
        v = np.array([r[key] for r in sub], dtype=float)
        for g, name_g in ((1, "reported"), (0, "extra")):
            x = v[lab == g]
            if x.size:
                qs = np.percentile(x, [0, 25, 50, 75, 100])
                print(f"{key:<22} {name_g:<8} " + " ".join(f"{q:10.5g}" for q in qs))
        # best single threshold
        order = np.unique(v)
        acc = 0.0
        thr = None
        for t in order:
            for sign in (1, -1):
                pred = (v > t) if sign > 0 else (v <= t)
                a = (pred == (lab == 1)).mean()
                if a > acc:
                    acc, thr = a, (t, sign)
        best_split.append((acc, key, thr))
    print("\nbest single-threshold accuracy (0.583 = predict-majority baseline)")
    for acc, key, thr in sorted(best_split, reverse=True)[:10]:
        print(f"  {key:<22} acc={acc:.4f}  {thr}")


if __name__ == "__main__":
    main(*sys.argv[1:])
