"""What single-marker change would make a `.seg` row's *third* column right?

After `19-…` §3-§5 the `IBD1Seg` and `IBD2Seg` columns are exact on all 982 corpus rows,
but `PropIBD = IBD2Seg + IBD1Seg/2` — printed from the unrounded totals — still differs on
some, which means base-pair errors too small for either column to show.  A base pair is
not a free parameter: every total is a sum of `pos[b] - pos[a]` over marker indices, so an
error is a whole number of marker gaps.

This enumerates every one-marker move of every call end, keeps the moves that make **all
four** columns match, and tags them, so the fix is read off a pattern rather than fitted.
"""

import sys
from collections import Counter

import engine as E
import kingdata as kd
import seg19 as S

WORD = 64


def totals(ds, c1, c2, min_bp=E.SEGLEN):
    pos = ds.pos
    ibd2 = sum(int(pos[b] - pos[a]) for a, b in c2)
    ibd1 = 0
    for a, b in c1:
        ibd1 += sum(v for v in (int(pos[y] - pos[x])
                                for x, y in E._pieces((a, b), c2)) if v >= min_bp)
    return ibd1, ibd2


def row_ok(ds, ibd1, ibd2, ref):
    d = ds.denom
    a1, a2, ap, at = ref
    g1, g2 = ibd1 / d, ibd2 / d
    gp = g2 + g1 / 2
    return (kd.fmt4(g1) == a1 and kd.fmt4(g2) == a2 and kd.fmt4(gp) == ap
            and kd.inf_type(g1, g2, gp) == at)


def calls(ds, i, j, min_bp=E.SEGLEN):
    c1, c2, owner = [], [], []
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        a2 = S.ibd2_19(sc, ds, i, j, S.R19(), ds.pos, min_bp)
        a1 = sc.ibd1(ds.pos, min_bp)
        for t, (x, y) in enumerate(a1):
            c1.append((x, y))
            owner.append(("ibd1", sc, t, len(a1)))
        for t, (x, y) in enumerate(a2):
            c2.append((x, y))
            owner.append(("ibd2", sc, t, len(a2)))
    return c1, c2, owner


def tag(kind, sc, idx, ntot, end, step, x, y):
    """A structural label for one candidate move."""
    t = [kind, end, "%+d" % step]
    t.append("first" if idx == 0 else ("last" if idx == ntot - 1 else "mid"))
    v = x if end == "left" else y
    if v == sc.lo:
        t.append("at seg.lo")
    elif v == sc.hi:
        t.append("at seg.hi")
    elif v % WORD == 0:
        t.append("at 64u")
    elif v % WORD == WORD - 1:
        t.append("at 64u+63")
    else:
        t.append("interior")
    if v < WORD * sc.w0 or v > WORD * (sc.w1 + 1) - 1:
        t.append("in fringe")
    return " ".join(t)


def analyse(names=None, verbose=False):
    hits = Counter()
    rows = fixed = 0
    for name in (names or [n for n in kd.DATASETS if n != "bigish"]):
        ds = kd.load(name)
        for (i, j), ref in sorted(ds.ref.items()):
            c1, c2, owner = calls(ds, i, j)
            base = totals(ds, c1, c2)
            if row_ok(ds, base[0], base[1], ref):
                continue
            rows += 1
            found = []
            for k in range(len(c1) + len(c2)):
                kind, sc, idx, ntot = owner[k]
                src = c1 if kind == "ibd1" else c2
                pos_in = k if kind == "ibd1" else k - len(c1)
                # `owner` is built in the same order the two lists are appended, so
                # recover the index inside its own list.
                pos_in = ([q for q, o in enumerate(owner) if o[0] == kind]
                          .index(k))
                x, y = src[pos_in]
                for end in ("left", "right"):
                    for step in (-1, 1):
                        nx, ny = (x + step, y) if end == "left" else (x, y + step)
                        if nx > ny or nx < sc.lo or ny > sc.hi:
                            continue
                        n1 = list(c1)
                        n2 = list(c2)
                        (n1 if kind == "ibd1" else n2)[pos_in] = (nx, ny)
                        t1, t2 = totals(ds, n1, n2)
                        if row_ok(ds, t1, t2, ref):
                            found.append(tag(kind, sc, idx, ntot, end, step, x, y))
            if found:
                fixed += 1
                for f in set(found):
                    hits[f] += 1
            if verbose:
                print("  %-11s %-9s %-9s  %d single-marker fixes: %s"
                      % (name, ds.fam[i][1], ds.fam[j][1], len(found),
                         sorted(set(found))[:4]))
    print("=== rows whose PropIBD is wrong: %d; fixable by ONE marker: %d" % (rows, fixed))
    for k, v in hits.most_common(30):
        print("  %4d  %s" % (v, k))


if __name__ == "__main__":
    analyse(sys.argv[1:] or None, verbose="-v" in sys.argv)
