"""Dense word-grid sweep: one pair, 64 alignments, the reference's verdict for each.

`gridshift.py` showed the verdict on a borderline two-word run moves when the global
64-marker grid moves.  Genotypes are constant along this sweep, so everything the rule can
depend on is a function of the alignment alone — which makes 64 rows per pair a far
sharper instrument than 437 pairs with one alignment each.

For every shift `m` (the first `m` markers of the fileset are deleted) this prints the
run's word span, our predicted boundaries and length, and whether the reference emitted a
row at all.

    python3 sweep.py [n_reported] [n_extra] [max_shift]
"""

import sys
import tempfile

import numpy as np

import kingdata as kd
import gridshift as G
import probe_seg as PS
import rules2 as R

BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")
PC = np.bitwise_count


def predict(ds, i, j, anchor, drop):
    """Our call for the run containing marker `anchor`, on the grid shifted by `drop`.

    Returns (u, v, lo, hi, length_bp) in *shifted* word indices but original marker
    indices, or None when the run has no complete clean word left.
    """
    ibs0, _, _, _ = ds.masks(i, j)

    def is0(t):
        return 0 <= t < ds.pos.size and bool((int(ibs0[t // 64]) >> (t % 64)) & 1)

    # the clean marker stretch around the anchor
    left = anchor
    while left > 0 and not is0(left - 1):
        left -= 1
    right = anchor
    while right + 1 < ds.pos.size and not is0(right + 1):
        right += 1
    # complete words of the shifted grid inside it
    u = -(-(left - drop) // 64)
    v = (right - drop + 1) // 64 - 1
    if v < u:
        return None
    lo_m, hi_m = 64 * u + drop, 64 * v + 63 + drop
    # extend: back to one past the previous word's last IBS0, out to the next word's last
    prev = [t for t in range(64 * (u - 1) + drop, lo_m) if is0(t)]
    nxt = [t for t in range(hi_m + 1, min(64 * (v + 2) + drop, ds.pos.size)) if is0(t)]
    lo = prev[-1] + 1 if prev else lo_m
    hi = nxt[-1] if nxt else hi_m
    return u, v, lo, hi, int(ds.pos[hi] - ds.pos[lo])


def main(n_rep=2, n_extra=2, max_shift=64):
    n_rep, n_extra, max_shift = int(n_rep), int(n_extra), int(max_shift)
    ds = kd.load("bigish")
    groups = {"reported": [], "extra": []}
    for (i, j) in ds.pairs():
        _, _, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        u, v = (best[1] + 63) // 64, (best[2] + 1) // 64 - 1
        if len(detail) != 1:
            continue
        ibs0, _, _, _ = ds.masks(i, j)
        n0 = PC(ibs0)
        uu = u
        while uu - 1 >= 0 and n0[uu - 1] == 0:
            uu -= 1
        vv = (best[2] + 1) // 64 - 1
        while n0[vv] != 0:
            vv -= 1
        while vv + 1 < ds.nwords and n0[vv + 1] == 0:
            vv += 1
        if vv - uu + 1 != 2:
            continue
        groups["reported" if (i, j) in ds.ref else "extra"].append(
            (i, j, 64 * uu, int(ds.chr[best[1]])))

    probes = ([("reported", *t) for t in groups["reported"][:n_rep]]
              + [("extra", *t) for t in groups["extra"][:n_extra]])
    used = {x for t in probes for x in t[1:3]}
    pad = [s for s in range(len(ds.fam)) if s not in used][:28]

    for label, i, j, anchor, chrom in probes:
        name = f"{ds.fam[i][1]}/{ds.fam[j][1]}"
        print(f"\n=== {label} {name} anchor={anchor} chr={chrom}")
        print(f"{'m':>3} {'words':>5} {'u':>6} {'lo':>7} {'hi':>7} {'ours Mb':>8} "
              f"{'ref':>6} {'ref Mb':>8}")
        keep = sorted(set([i, j] + pad))
        for m in range(max_shift):
            p = predict(ds, i, j, anchor, m)
            with tempfile.TemporaryDirectory() as td:
                G.write_shifted(ds, keep, m, chrom, 2, td)
                rows, denom = PS.run_king(td)
            key = (ds.fam[i][1], ds.fam[j][1])
            r = rows.get(key) or rows.get((key[1], key[0]))
            refmb = f"{r[0] * denom / 2 / 1e6:8.3f}" if r else "       -"
            if p is None:
                print(f"{m:>3} {'-':>5} {'-':>6} {'-':>7} {'-':>7} {'-':>8} "
                      f"{'YES' if r else 'no':>6} {refmb}")
            else:
                u, v, lo, hi, ln = p
                print(f"{m:>3} {v - u + 1:>5} {u:>6} {lo:>7} {hi:>7} {ln / 1e6:8.3f} "
                      f"{'YES' if r else 'no':>6} {refmb}")


if __name__ == "__main__":
    main(*sys.argv[1:])
