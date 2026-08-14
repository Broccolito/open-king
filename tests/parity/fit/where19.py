"""**SUPERSEDED — read this first.** The premise below is wrong, and the way it is wrong is
worth keeping. It assumes `.seg`'s `PropIBD` is a third *rounding* of the same two base-pair
totals, so a disagreement there must mean our totals are slightly off. They are not:
`docs/research/20-seg-writer.md` shows `.seg` computes `PropIBD` from the **printed**
columns (`i2*1e-4 + i1*5e-5`), so it carries no information about the totals at all and the
"hexagon" this module derives is empty of signal. Every correction it reports was an artefact
of comparing two different formulas. Kept because `19-ibd2seg-residual.md` quotes its
numbers and because it is a clean example of a well-built instrument pointed at the wrong
question.

Which IBD1 endpoint is short?  One-marker perturbations, row by row.

`prop19.py` says 115 of the 176 wrong rows at the default floor carry no IBD2 on either
side, so their `PropIBD` is `IBD1Seg/2` and the fault is the IBD1 pass alone; 101 of them
need us to report *more*, by about one marker gap.  `merge19.py` ruled out the obvious
source (reporting touching calls as one), and `18-…` §5 had already measured that dead.

So the missing base pairs are at an **endpoint of one call**.  This enumerates the
single-marker perturbations of our own call set — move one call's left end out or in by
one marker, or its right end — and reports which ones would land `PropIBD` inside the
band the three printed columns allow.  If one *kind* of endpoint explains most rows, that
names the clause; if the working perturbations are scattered, the residual is not a single
one-marker rule and this says so.

Endpoints are classified by what stopped them, so the answer is a clause and not a row
number:

    seg.lo / seg.hi      the usable segment's own first/last marker (no fringe breaker)
    fringe               a fringe stop: the nearest IBS0 in the partial word
    grid                 the word-grid edge with no fringe at all (segment word-aligned)
    ibs0                 the interior rule: the flanking word's last opposite homozygote

    python3 where19.py            # the summary
    python3 where19.py rows       # ...plus each row's working perturbations
"""

import sys
from collections import Counter

import engine as E
import kingdata as kd
import prop19 as P
import seg19 as S19

WORD = E.WORD
RULE = S19.R19()


def classify(ds, sc, seg, m, side):
    """What stopped an IBD1 call at marker `m` on the given side."""
    lo, hi, w0, w1 = sc.lo, sc.hi, sc.w0, sc.w1
    if side == "lo":
        if m == lo:
            return "seg.lo"
        if m < WORD * w0:
            return "fringe.lo"
        if m == WORD * w0:
            return "grid.lo"
        return "ibs0.lo"
    if m == hi:
        return "seg.hi"
    if m > WORD * (w1 + 1) - 1:
        return "fringe.hi"
    if m == WORD * (w1 + 1) - 1:
        return "grid.hi"
    return "ibs0.hi"


def row_calls(ds, i, j, min_bp):
    """Our IBD1 and IBD2 calls, with the scan that produced each, per usable segment."""
    out = []
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = S19.ibd2_19(sc, ds, i, j, RULE, ds.pos, min_bp)
        c1 = sc.ibd1(ds.pos, min_bp)
        out.append((seg, sc, c1, c2))
    return out


def analyse(min_bp=3_000_000, suffix="__ibdseg", verbose=False):
    hit = Counter()
    tally = Counter()
    rows = fixable = 0
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        pos = ds.pos
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            if (i, j) not in ref:
                continue
            per = row_calls(ds, i, j, min_bp)
            ibd1 = ibd2 = 0
            for _seg, _sc, c1, c2 in per:
                ibd2 += sum(int(pos[b] - pos[a]) for a, b in c2)
                for lo, hi in c1:
                    ibd1 += sum(v for v in (int(pos[y] - pos[x])
                                            for x, y in E._pieces((lo, hi), c2))
                                if v >= min_bp)
            r1, r2, rp, _rt = ref[(i, j)]
            if b_ok(ibd1, ibd2, d, r1, r2, rp):
                continue
            # only the rows this instrument can speak about
            if not (ibd2 == 0 and r2 == 0.0):
                continue
            rows += 1
            box = P.feasible(r1, r2, rp, d)
            need = P.deficit(ibd1, ibd2, box)      # bp of ibd2 + ibd1/2
            works = []
            for _seg, sc, c1, _c2 in per:
                for idx, (lo, hi) in enumerate(c1):
                    for side, m2, tag in (("lo", lo - 1, "extend lo"),
                                          ("lo", lo + 1, "shrink lo"),
                                          ("hi", hi + 1, "extend hi"),
                                          ("hi", hi - 1, "shrink hi")):
                        if not (sc.lo <= m2 <= sc.hi):
                            continue
                        delta = ((pos[lo] - pos[m2]) if side == "lo"
                                 else (pos[m2] - pos[hi]))
                        newp = ibd2 + (ibd1 + int(delta)) / 2.0
                        if box[4] <= newp < box[5]:
                            k = classify(ds, sc, _seg, lo if side == "lo" else hi, side)
                            works.append((tag, k))
            if works:
                fixable += 1
                for tag, k in set(works):
                    hit[(tag, k)] += 1
                tally[tuple(sorted({t for t, _ in works}))] += 1
            else:
                tally[("none",)] += 1
            if verbose:
                ds_ = ds
                print("  %-12s %-9s %-9s need %+8.0f bp   %s"
                      % (name, ds_.fam[i][1], ds_.fam[j][1], need,
                         ", ".join(sorted({"%s@%s" % (t, k) for t, k in works}))
                         or "NO single-marker fix"))
    print("=== IBD2-free rows wrong on PropIBD at %d Mb: %d"
          % (min_bp // 1_000_000, rows))
    print("  fixable by moving ONE endpoint by ONE marker: %d" % fixable)
    print("  which perturbation, and what stopped that endpoint:")
    for (tag, k), n in sorted(hit.items(), key=lambda kv: -kv[1]):
        print("    %-10s at %-10s %4d" % (tag, k, n))
    print("  rows by the set of perturbations that work:")
    for k, n in sorted(tally.items(), key=lambda kv: -kv[1])[:10]:
        print("    %-40s %4d" % (",".join(k), n))


def b_ok(a, b, d, r1, r2, rp):
    return (kd.fmt4(a / d) == r1 and kd.fmt4(b / d) == r2
            and kd.fmt4(b / d + a / 2.0 / d) == rp)


if __name__ == "__main__":
    analyse(verbose=len(sys.argv) > 1 and sys.argv[1] == "rows")
