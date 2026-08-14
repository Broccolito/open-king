"""Where does the reference put an IBD2 segment's right end? A labelled table.

For every localised `MaxIBD2` target this pairs the reference's own interval `[a, b]`
(from `invert.py`) with the run `[u, v]` the current word predicate produces, and prints
`b - v` — 1 when the call reaches one word past its last clean word, 0 when it does not.
Every candidate separating feature is printed alongside, so the question "what makes the
reference decline the extension" is answered by reading a column rather than by guessing.

    python3 endfit.py            # the contingency tables
    python3 endfit.py rows       # one tab-separated line per target
"""

import sys
from collections import Counter, defaultdict

import kingdata as kd
import engine as E
import invert

WORD = 64


def runs_of(ds, i, j, p=E.BASE):
    """Every IBD2 word run the predicate produces, as (u, v, w0, w1)."""
    _, n0, n1, _c1, _c2 = E.masks(ds, i, j)
    out = []
    for _chrom, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        sl = slice(w0, w1 + 1)
        clean = (n0[sl] == 0) & (n1[sl] < p.ibd2_dirty_ibs1)
        ok = clean.copy()
        if p.bridge:
            nn = n0[sl]
            for k in range(1, len(clean) - 1):
                if not clean[k] and clean[k - 1] and clean[k + 1] and nn[k] == 0:
                    ok[k] = True
        for a, b in E._runs(ok):
            out.append((w0 + a, w0 + b, w0, w1))
    return out


def rows(p=E.BASE):
    """One record per localised target: the reference interval and its run context."""
    out = []
    for name, i, j, t in E.max_targets():
        ds = kd.load(name)
        cs = invert.candidates(ds, i, j, t)
        if len(cs) != 1:
            continue
        a, b, w0, w1 = cs[0]
        rs = runs_of(ds, i, j, p)
        # the run that owns the reference's start word
        own = [r for r in rs if r[0] <= a <= r[1]]
        if not own:
            out.append(dict(ds=name, i=i, j=j, a=a, b=b, w0=w0, w1=w1, u=None, v=None))
            continue
        u, v, _, _ = own[0]
        _, n0, n1, _c1, c2 = E.masks(ds, i, j)

        def nx(w, arr):
            return int(arr[w]) if 0 <= w < len(arr) else -1

        nxt = [r for r in rs if r[0] > v]
        out.append(dict(
            ds=name, i=i, j=j, a=a, b=b, w0=w0, w1=w1, u=u, v=v,
            db=b - v, da=a - u,
            runlen=v - u + 1,
            atedge_end=int(b == w1), atedge_start=int(a == w0),
            v_is_w1=int(v == w1), v1_is_w1=int(v + 1 == w1), v2_is_w1=int(v + 2 == w1),
            ibs1_v=nx(v, n1), ibs1_v1=nx(v + 1, n1), ibs1_v2=nx(v + 2, n1),
            ibs0_v1=nx(v + 1, n0), ibs0_v2=nx(v + 2, n0),
            inf_v1=int(c2[v + 2] - c2[v + 1]) if v + 2 <= len(n1) else -1,
            gap_next=(nxt[0][0] - v) if nxt else -1,
            nextlen=(nxt[0][1] - nxt[0][0] + 1) if nxt else -1,
            bridged=int(any(n1[w] >= p.ibd2_dirty_ibs1 for w in range(u, v + 1))),
        ))
    return out


def table(recs, key, label):
    c = defaultdict(Counter)
    for r in recs:
        c[r[key]][r["db"]] += 1
    print("\n%-14s   " % label + "  ".join("db=%s" % d for d in (0, 1, 2)))
    for k in sorted(c, key=lambda x: (x is None, x)):
        row = c[k]
        print("  %-12s " % k + "  ".join("%4d " % row.get(d, 0) for d in (0, 1, 2)))


def main():
    recs = [r for r in rows() if r["u"] is not None]
    print("localised targets with a matching run: %d" % len(recs))
    print("db = b - v distribution:", Counter(r["db"] for r in recs))
    print("da = a - u distribution:", Counter(r["da"] for r in recs))
    if "rows" in sys.argv:
        keys = ["ds", "i", "j", "a", "b", "u", "v", "db", "da", "runlen", "ibs1_v1",
                "ibs1_v2", "ibs0_v1", "v_is_w1", "v1_is_w1", "v2_is_w1", "inf_v1",
                "gap_next", "nextlen", "bridged"]
        print("\t".join(keys))
        for r in sorted(recs, key=lambda r: (-r["db"], r["ds"])):
            print("\t".join(str(r[k]) for k in keys))
        return
    for k, lab in [("v1_is_w1", "v+1 == w1"), ("v2_is_w1", "v+2 == w1"),
                   ("v_is_w1", "v == w1"), ("bridged", "run bridged"),
                   ("gap_next", "gap to next run"), ("nextlen", "next run len")]:
        table(recs, k, lab)
    for k, lab in [("ibs1_v1", "IBS1(v+1)"), ("ibs1_v2", "IBS1(v+2)")]:
        buckets = []
        for r in recs:
            x = r[k]
            buckets.append(dict(r, **{k: "%2d-%2d" % (x // 5 * 5, x // 5 * 5 + 4)}))
        table(buckets, k, lab)


if __name__ == "__main__":
    main()
