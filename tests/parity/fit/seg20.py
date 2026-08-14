"""The `--seglength` run merge — `docs/research/20-seglength-floor.md`.

`seg19.py` closed both estimate columns at the default 3 Mb floor (982 / 982) and left
`IBD1Seg` at 910 / 844 and `IBD2Seg` at 946 / 937 at `--seglength 5` and `10`.  The whole
residual is one clause: **above the floor the reference joins two runs across a short
interruption**.  `18-ibd1-caller.md` §9 measured two of its five conditions and
deliberately did not implement it, because the two-condition version makes those floors
much worse; the three missing conditions are what this module adds.

Every constant is bisected on `docs/research/fixtures/mergelab.py` canvases against the
KING 2.3.2 reference and validated out of sample on unused seeds — never fitted to the
corpus.  Both passes merge, with the same geometry:

    the runs of one pass, **after the gate has refused what it refuses**, are walked in
    order; two adjacent survivors are joined iff

      (a) at most 2 *unusable* words lie between them                     §3
          (a gate-refused run lies inside the interruption and is stepped over: it can
          never end a merged segment, but it does not stop one either)    §6
      (b) pos[first marker of the later run] - pos[last marker of the earlier run]
            <  --seglength,  strictly                                     §2
      (c) cost * (bad - 2) <= X, summed over those unusable words         §4, §5, §7

              IBD1   bad = opposite homozygotes                cost = 4
                     X   = A1A1/A1A1 markers  (`inf1 & ~ibs1`)
                           ...unless the het-vs-A1A1 markers (`inf1 & ibs1`) alone
                           reach 10, in which case X is those instead
              IBD2   bad = opposite homozygotes + het-vs-hom mismatches   cost = 3
                     X   = `inf2` = HetHet + A1A1/A1A1  (`share & ~ibs1`)

The merged run then takes the gate, the endpoints and the `--seglength` floor exactly as
an unmerged one does.  A merged call may **not** satisfy the ">10 Mb" pair-reporting
filter (§8).

    python3 seg20.py            # the corpus scorecard at 3 / 5 / 10 Mb, before and after
    python3 seg20.py grid       # every knob of the merge swept

**Read the `exact` column with care.** Like `seg19.py`, which this module extends, it
grades `PropIBD` with the retired **`.kin`** rule (`IBD2Seg + IBD1Seg/2`, unrounded), so
that the "before" and "after" rows are comparable with that write-up's. The committed
engine prints `.seg`'s own rule instead (`engine.seg_prop_ibd`), under which the same
runs score 947 and 943 exact rather than 795 and 793. The two `IBD1Seg`/`IBD2Seg` columns
are unaffected and are the ones this write-up is about. For the committed engine's own
scorecard, measured from the binary against the goldens, use `scorecard.py`.
"""

import sys
from dataclasses import dataclass, replace

import numpy as np

import engine as E
import kingdata as kd
import seg19 as S19

WORD = E.WORD
PC = np.bitwise_count

MERGE_GATE = 10          # the het-informative count at which `X` switches from U to V
MERGE_FREE = 2           # bad markers a merge gets for nothing
MERGE_MULT1 = 4          # informative markers one further bad marker costs, IBD1 pass
MERGE_MULT2 = 3          # ...and IBD2
MERGE_WORDS = 2          # unusable words a merge may bridge


@dataclass(frozen=True)
class R20:
    """`19-…`'s geometry plus the merge.  Every field is a knob a section bisects."""

    base: S19.R19 = S19.R19()
    merge1: bool = True          # merge IBD1 runs
    merge2: bool = True          # merge IBD2 runs
    words: int = MERGE_WORDS     # how many unusable words a merge may bridge
    free: int = MERGE_FREE
    gate: int = MERGE_GATE
    mult1: int = MERGE_MULT1
    mult2: int = MERGE_MULT2
    stat: str = "switch"         # "switch" | "u" | "v" | "sum" — what `X` is
    span: str = "unusable"       # which words between two runs the budget is summed over
    gap: str = "runs"            # "runs" | "calls" — what the floor is compared to
    # Whether a merged call may satisfy the ">10 Mb" pair-reporting filter.  It may not:
    # the reference reports the **same pair set** at 3, 5 and 10 Mb on all ten corpus
    # datasets, and the merge is floor-dependent, so whatever the filter reads is not it.
    filter_merged: bool = False


# ---------------------------------------------------------------------------
# the three per-word counts the merge test reads
# ---------------------------------------------------------------------------

_C = {}


def counts(ds, i, j):
    """`(Z, U, V, U2, V2)` per-word popcount vectors for one pair.

    `U`/`V` split `inf1` (both carry A1, at least one homozygous) by whether the marker is
    also a het-vs-hom mismatch: `U` is A1A1/A1A1, `V` is het-vs-A1A1.  `U2`/`V2` do the
    same to the IBD2 pass's `share = p1i & p1j`, so `U2` is HetHet + A1A1/A1A1 — the very
    `inf2` the `.seg` gate counts — and `V2` is again het-vs-A1A1.
    """
    key = (ds.name, i, j)
    v = _C.get(key)
    if v is None:
        p0i, p1i = ds.p0[i], ds.p1[i]
        p0j, p1j = ds.p0[j], ds.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        ibs0 = p0i & p0j & (p1i ^ p1j)
        ibs1 = (het_i & p0j) | (p0i & het_j)
        share = p1i & p1j
        inf1 = share & (p0i | p0j)
        v = (PC(ibs0).astype(np.int32), PC(inf1 & ~ibs1).astype(np.int32),
             PC(inf1 & ibs1).astype(np.int32), PC(share & ~ibs1).astype(np.int32),
             PC(ibs1).astype(np.int32))
        _C[key] = v
    return v


def stat_x(U, V, p):
    """`X`, the informative count the bad-marker budget is measured against."""
    if p.stat == "u":
        return U
    if p.stat == "v":
        return V
    if p.stat == "sum":
        return U + V
    return V if V >= p.gate else U


def mergeable(bad, U, V, mult, p):
    """`mult * (bad - FREE) <= X` — two bad markers are free, each further one costs
    `mult` informative markers.  `bad` is the opposite homozygotes on the IBD1 pass and
    the opposite homozygotes plus het-vs-hom mismatches on the IBD2 one."""
    return mult * (bad - p.free) <= stat_x(U, V, p)


# ---------------------------------------------------------------------------
# the merge, over a run list
# ---------------------------------------------------------------------------

def join(runs, usable, pos, min_bp, uc, vc, mult, p, badc):
    """Runs (global word index pairs) with the merge of `20-…` applied, left to right.

    `runs` are the runs that **survived the gate** — the gate is asked first, and a run it
    refuses lies inside the interruption rather than ending it.  Only *unusable* words
    count toward the two-word cap, and only they are summed for the budget test.
    """
    out = []
    for a, b in runs:
        if out:
            pa, pb = out[-1]
            mid = [k for k in range(pb + 1, a) if not usable[k]]
            if 1 <= len(mid) <= p.words:
                lo = WORD * (pb + 1) - 1          # last marker of the earlier run
                hi = WORD * a                     # first marker of the later run
                if int(pos[hi] - pos[lo]) < min_bp:
                    if p.span == "all":
                        mid = list(range(pb + 1, a))
                    bad = int(sum(badc[k] for k in mid))
                    U = int(sum(uc[k] for k in mid))
                    V = int(sum(vc[k] for k in mid))
                    if mergeable(bad, U, V, mult, p):
                        out[-1] = (pa, b)
                        continue
        out.append((a, b))
    return out


# ---------------------------------------------------------------------------
# the two callers
# ---------------------------------------------------------------------------

def ibd1_20(sc, ds, i, j, pos, min_bp, p):
    """`Scan::ibd1` with the merge.  Everything else is `18-…` §7, unchanged."""
    if sc.n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    ok = sc.n0[w0:w1 + 1] == 0
    runs = [(w0 + a, w0 + b) for a, b in E._runs(ok)]
    # The gate first — `20-…` §6: a run under `inf1 >= 10` is refused outright and can
    # neither end a merge nor stop one.
    runs = [r for r in runs if sc.informative(sc.cum1, *r)]
    if p.merge1:
        zc, uc, vc, _u2, _mis = counts(ds, i, j)
        usable = {w0 + k: bool(v) for k, v in enumerate(ok)}
        runs = join(runs, usable, pos, min_bp, uc, vc, p.mult1, p, zc)
    out = []
    for u, v in runs:
        hi = sc.right_end(v)
        lo = sc.left_end(u)
        out = sc._emit(out, lo, hi, pos, min_bp)
    return out


def ibd2_20(sc, ds, i, j, p, pos, min_bp):
    """`seg19.ibd2_19` with the merge applied to its run list.

    The bridge of `17-…` §14 runs first and is untouched: it absorbs a *mismatch-only*
    unusable word at any floor, where this clause joins across words the bridge refused
    (an opposite homozygote makes a word un-bridgeable) and only under the floor.
    """
    n = sc.n
    if n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    ibs0, ibs1, n0, n1, cum = S19.masks(ds, i, j)
    z = [int(n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < p.base.dirty for k in range(n)]
    head, tail = S19.fringe_masks(sc, ibs1, ibs0, p.base)

    def ge_of(b):
        return b + 1 if (b + 1 < n and not z[b + 1] and mis[b + 1]) else b

    def gate_ok(g, b):
        return int(cum[w0 + ge_of(b) + 1] - cum[w0 + g]) >= p.base.gate

    ok = list(usable)
    gs0 = None
    for k in range(n):
        if usable[k]:
            if gs0 is None and mis[k] == 0:
                gs0 = k
            continue
        bridged = False
        if (gs0 is not None and k > 0 and not z[k] and k + 1 < n
                and usable[k + 1] and mis[k + 1] == 0):
            b2 = k + 1
            while b2 + 1 < n and usable[b2 + 1]:
                b2 += 1
            bridged = gate_ok(gs0, k - 1) and gate_ok(k + 1, b2)
        if bridged:
            ok[k] = True
        else:
            gs0 = None

    def gs_of(a, b):
        return next((t for t in range(a, b + 1) if mis[t] == 0), None)

    runs = []
    for a, b in E._runs(np.array(ok)):
        g = gs_of(a, b)
        if g is not None and gate_ok(g, b):
            runs.append((w0 + a, w0 + b))
    if p.merge2:
        zc, _u, _v, u2c, misc = counts(ds, i, j)
        # The IBD2 pass counts a het-vs-hom mismatch as a bad marker alongside an
        # opposite homozygote, and its informative count is `inf2 = share & ~ibs1`.
        usable = {w0 + k: bool(v) for k, v in enumerate(ok)}
        runs = join(runs, usable, pos, min_bp, u2c, np.zeros_like(u2c), p.mult2, p,
                    zc + misc)

    head_stop = (WORD * (w0 - 1) + S19._last(head) + p.base.fringe_off) if head else sc.lo
    tail_stop = (WORD * (w1 + 1) + S19._first(tail) - p.base.fringe_off) if tail else sc.hi
    head_stop = max(head_stop, sc.lo)
    tail_stop = min(tail_stop, sc.hi)

    out, emitted = [], 0
    for u, v in runs:
        a, b = u - w0, v - w0
        left = WORD * u
        if a > 0 and not z[a - 1] and int(ibs1[u - 1]):
            left = max(0, WORD * (u - 1) + S19._last(int(ibs1[u - 1])) - p.base.reach)
            if a < 2 or z[a - 2]:
                left = max(left, WORD * (u - 1))
        if left <= WORD * w0:
            left = head_stop
        right = WORD * v + WORD - 1
        if b + 1 < n and not z[b + 1] and int(ibs1[v + 1]):
            right = WORD * (v + 1) + S19._first(int(ibs1[v + 1])) + p.base.reach
            if b + 2 >= n or z[b + 2]:
                right = min(right, WORD * (v + 2) - 1)
        if right >= WORD * (w1 + 1) - 1:
            right = tail_stop
        gs = gs_of(a, b)
        if gs is None or not gate_ok(gs, b):
            continue
        if emitted:
            left = max(left, WORD * (w0 + gs + 1))
        emitted += 1
        left, right = max(left, sc.lo), min(right, sc.hi)
        if out:
            left = max(left, out[-1][1] + p.base.clip)
        if left <= right and pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


# ---------------------------------------------------------------------------
# aggregation and scoring
# ---------------------------------------------------------------------------

def call_pair(ds, i, j, p, min_bp=E.SEGLEN):
    pos = ds.pos
    ibd1 = ibd2 = longest = 0
    q = replace(p, merge1=False, merge2=False) if not p.filter_merged else p
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = ibd2_20(sc, ds, i, j, p, pos, min_bp)
        c1 = ibd1_20(sc, ds, i, j, pos, min_bp, p)
        # The pair-reporting filter reads the *unmerged* calls (see `R20.filter_merged`).
        f2 = c2 if q is p else ibd2_20(sc, ds, i, j, q, pos, min_bp)
        f1 = c1 if q is p else ibd1_20(sc, ds, i, j, pos, min_bp, q)
        for lo, hi in c2:
            ibd2 += int(pos[hi] - pos[lo])
        for lo, hi in f2 + f1:
            longest = max(longest, int(pos[hi] - pos[lo]))
        for lo, hi in c1:
            ibd1 += sum(v for v in (int(pos[y] - pos[x])
                                    for x, y in E._pieces((lo, hi), c2))
                        if v >= min_bp)
    return ibd1, ibd2, longest


def score(p, min_bp=E.SEGLEN, suffix="__ibdseg", datasets=None):
    rows = exact = i1 = i2 = extra = missing = 0
    err = worst = 0.0
    for name in (datasets or kd.DATASETS):
        ds = kd.load(name)
        d = ds.denom
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            a, b, lg = call_pair(ds, i, j, p, min_bp)
            got, want = lg >= E.LONG, (i, j) in ref
            if not want:
                extra += got
                continue
            if not got:
                missing += 1
                continue
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = a / d, b / d
            gp = g2 + g1 / 2
            rows += 1
            ok1, ok2 = kd.fmt4(g1) == a1, kd.fmt4(g2) == a2
            i1 += ok1
            i2 += ok2
            exact += (ok1 and ok2 and kd.fmt4(gp) == ap
                      and kd.inf_type(g1, g2, gp) == at)
            err += abs(gp - ap)
            worst = max(worst, abs(gp - ap))
    return dict(rows=rows, exact=exact, ibd1=i1, ibd2=i2, extra=extra, missing=missing,
                mae=err / max(rows, 1), worst=worst)


# ---------------------------------------------------------------------------
# the canvas mirror
# ---------------------------------------------------------------------------

def predict(info, pos, min_bp, w0=0, w1=None, lo=None, hi=None, p=R20()):
    """`.seg` IBD1 calls as marker intervals — the canvas mirror of `ibd1_20`.

    `info[k]` describes global word `k` as `(ibs0_count, U, V, last_ibs0_bit, inf1)`:

        ibs0_count     opposite homozygotes in the word
        U              A1A1/A1A1 markers            (`inf1 & ~ibs1`)
        V              het-vs-A1A1 markers          (`inf1 &  ibs1`)
        last_ibs0_bit  bit index of the last opposite homozygote, or None
        inf1           U + V, the gate's own count

    `pos` is the marker position vector and `min_bp` is `--seglength`; with `min_bp = 0`
    no merge can fire and this reduces to `18-…` §7 exactly.  `mergelab.ibd1()` is the
    same function over word *compositions* and is what the reference is graded against.
    """
    n = len(info)
    w1 = n - 1 if w1 is None else w1
    lo = WORD * w0 if lo is None else lo
    hi = WORD * (w1 + 1) - 1 if hi is None else hi
    usable = [info[k][0] == 0 for k in range(n)]

    runs, k = [], w0
    while k <= w1:
        if not usable[k]:
            k += 1
            continue
        a = k
        while k <= w1 and usable[k]:
            k += 1
        runs.append((a, k - 1))
    # The gate first — a run it refuses can neither end a merge nor stop one.
    runs = [r for r in runs
            if sum(info[t][4] for t in range(r[0], r[1] + 1)) >= p.gate]

    joined = []
    for a, b in runs:
        if joined and p.merge1:
            pa, pb = joined[-1]
            mid = [t for t in range(pb + 1, a) if not usable[t]]
            if (1 <= len(mid) <= p.words
                    and pos[WORD * a] - pos[WORD * (pb + 1) - 1] < min_bp
                    and mergeable(sum(info[t][0] for t in mid),
                                  sum(info[t][1] for t in mid),
                                  sum(info[t][2] for t in mid), p.mult1, p)):
                joined[-1] = (pa, b)
                continue
        joined.append((a, b))

    out = []
    for a, b in joined:
        left = lo if a == w0 else (WORD * (a - 1) + info[a - 1][3] + 1
                                   if info[a - 1][3] is not None else WORD * a)
        right = hi if b == w1 else (WORD * (b + 1) + info[b + 1][3]
                                    if info[b + 1][3] is not None
                                    else WORD * (b + 2) - 1)
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1] + 1)
        if left <= right and pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


FLOORS = S19.FLOORS


def show(tag, s):
    print("%-34s exact %4d  ibd1 %4d  ibd2 %4d  extra %3d  miss %3d  MAE %.6f  "
          "worst %.4f"
          % (tag, s["exact"], s["ibd1"], s["ibd2"], s["extra"], s["missing"],
             s["mae"], s["worst"]))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    OFF = R20(merge1=False, merge2=False)
    for bp, sfx in FLOORS:
        print("--seglength %d Mb" % (bp // 1_000_000))
        show("  19 (no merge)", score(OFF, bp, sfx))
        show("  20 (merge, both passes)", score(R20(), bp, sfx))
        if mode == "grid":
            show("  merge on IBD1 only", score(replace(R20(), merge2=False), bp, sfx))
            show("  merge on IBD2 only", score(replace(R20(), merge1=False), bp, sfx))
            for w in (1, 3, 99):
                show("  words<=%d" % w, score(replace(R20(), words=w), bp, sfx))
            for st in ("u", "v", "sum"):
                show("  stat=%s" % st, score(replace(R20(), stat=st), bp, sfx))
            for m in (2, 3, 5, 6, 999):
                show("  mult1=%d" % m, score(replace(R20(), mult1=m), bp, sfx))
            for m in (2, 4, 5, 999):
                show("  mult2=%d" % m, score(replace(R20(), mult2=m), bp, sfx))
            for f in (0, 1, 3):
                show("  free=%d" % f, score(replace(R20(), free=f), bp, sfx))
