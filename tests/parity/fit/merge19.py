"""Does the reference report two *touching* IBD1 calls as one segment?

`prop19.py` localises the residual: at the default floor both estimate columns are exact
on all 982 rows and 176 rows still miss `PropIBD`, 115 of them on rows where neither side
reports any IBD2 at all — so those 115 are the IBD1 pass's, and 101 of the 115 need us to
report **more** base pairs, by a median of about one marker gap.

One marker gap is exactly what a *merge* is worth.  Two IBD1 calls separated by a single
bad word come out adjacent — the first ends on that word's last opposite homozygote and
the second starts one marker later — so reporting them as one segment rather than two adds
`pos[lo2] - pos[hi1]`, one marker interval, and nothing else.  The corpus has 287 such
touching pairs.

This scores that hypothesis, and the two neighbouring ones, against `PropIBD`.  It is a
*screen*, not evidence: whatever wins here has to be confirmed on the canvas, where the
call count is read back directly, before it goes anywhere near the engine.

    python3 merge19.py
"""

import sys

import engine as E
import kingdata as kd
import seg19 as S19

FLOORS = [(3_000_000, "__ibdseg"), (5_000_000, "__ibdseg_seglength5"),
          (10_000_000, "__ibdseg_seglength10")]
RULE = S19.R19()


def merged_len(pos, calls, mode):
    """Total base pairs of a list of IBD1 calls under one reporting convention."""
    if mode == "split":                      # the committed rule: each call on its own
        return sum(int(pos[b] - pos[a]) for a, b in calls)
    out = []
    for a, b in calls:
        if out and ((mode == "touch" and a == out[-1][1] + 1)
                    or (mode == "any" and a <= out[-1][1] + 1)):
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return sum(int(pos[y] - pos[x]) for x, y in out)


def call_pair(ds, i, j, min_bp, mode):
    pos = ds.pos
    ibd1 = ibd2 = longest = 0
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = S19.ibd2_19(sc, ds, i, j, RULE, pos, min_bp)
        c1 = sc.ibd1(pos, min_bp)
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2 += ln
            longest = max(longest, ln)
        keep = []
        for lo, hi in c1:
            longest = max(longest, int(pos[hi] - pos[lo]))
            for x, y in E._pieces((lo, hi), c2):
                if int(pos[y] - pos[x]) >= min_bp:
                    keep.append((x, y))
        ibd1 += merged_len(pos, keep, mode)
    return ibd1, ibd2, longest


def score(mode, min_bp, suffix):
    rows = exact = i1 = i2 = ip = 0
    err = 0.0
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            if (i, j) not in ref:
                continue
            a, b, _lg = call_pair(ds, i, j, min_bp, mode)
            r1, r2, rp, rt = ref[(i, j)]
            g1, g2 = a / d, b / d
            gp = g2 + g1 / 2
            rows += 1
            ok1, ok2 = kd.fmt4(g1) == r1, kd.fmt4(g2) == r2
            okp = kd.fmt4(gp) == rp
            i1 += ok1
            i2 += ok2
            ip += okp
            exact += ok1 and ok2 and okp and kd.inf_type(g1, g2, gp) == rt
            err += abs(gp - rp)
    return dict(rows=rows, exact=exact, ibd1=i1, ibd2=i2, prop=ip, mae=err / rows)


if __name__ == "__main__":
    modes = sys.argv[1:] or ["split", "touch", "any"]
    for bp, sfx in FLOORS:
        print("--seglength %d Mb" % (bp // 1_000_000))
        for m in modes:
            s = score(m, bp, sfx)
            print("  %-8s exact %4d  ibd1 %4d  ibd2 %4d  prop %4d  MAE %.6f"
                  % (m, s["exact"], s["ibd1"], s["ibd2"], s["prop"], s["mae"]))
