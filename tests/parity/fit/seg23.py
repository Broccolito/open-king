"""The gate-window length bound of `docs/research/23-gap-bound.md`.

`21-push-merge.md` §8.1 left `--seglength 10` wrong on 12 rows and named the merge's gap
as the suspect.  It is not the gap.  Muting all but one chromosome of a real corpus pair
(`docs/research/fixtures/chrprobe.py`) makes the reference print that chromosome's own
called length, and sweeping `--seglength` against it shows the fault as a single IBD2
call that the reference stops reporting at a floor **far below its own length**:

    multifam 2/4, chr13, call 11.2066 Mb     kept at 6 290 751 bp, dropped at 6 290 752
    bigish 66/69, chr6,  call 10.1511 Mb     kept at 6 297 345 bp, dropped at 6 297 346

Both flips are `2 * (the span of the call's gate window) + 1`, so the clause is

    pos[WORD * ge_of(e) + 63] - pos[WORD * gs]  >=  seglength / 2      (integer division)

— the same window the informativeness gate already reads (`13-…`, `19-…`) and the same
`seglength / 2` with the same integer division the push uses (`21-…` §2.5).  It is a
second condition on the gate, not on the reported call: the call itself is 11.2 Mb long
and its 3.1 Mb window is what fails.

    python3 seg23.py            # the corpus scorecard at 3 / 5 / 10 Mb, 21 vs 23
    python3 seg23.py grid       # where the bound is asked, and on which pass
"""

import sys
from dataclasses import dataclass, replace

import numpy as np

import engine as E
import kingdata as kd
import seg19 as S19
import seg20 as S20
import seg21 as S21

WORD = E.WORD

#: The fraction of `--seglength` a gate window must span. Bisected to the base pair on
#: two independent corpus calls; the division is the reference's own integer one.
WINDOW_FRACTION = 2


@dataclass(frozen=True)
class R23:
    """`21-…`'s rule plus the window bound. Every field is a knob a section drops."""

    base: S21.R21 = S21.R21()
    window: bool = True        # the bound, on the IBD2 pass
    pre_merge: bool = False    # asked with the gate rather than only at emit
    window1: bool = True       # the bound, on the IBD1 pass (strict, `> L/2`)
    span_all: bool = True      # IBD1 merge budget over every word between the two runs


def ibd1_23(sc, ds, i, j, pos, min_bp, p):
    """`seg20.ibd1_20` with the IBD1 pass's own window bound.

    The window is the run's own complete words — the very span `Scan::informative`
    counts over — and the comparison is **strict**, one unit of `min_bp / 2` tighter
    than the IBD2 pass's (`window1.py` §4, bisected at three spacings).
    """
    b = replace(p.base.base, span="all" if p.span_all else "unusable")
    if sc.n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    ok = sc.n0[w0:w1 + 1] == 0
    runs = [(w0 + a, w0 + b2) for a, b2 in E._runs(ok)]
    runs = [r for r in runs if sc.informative(sc.cum1, *r)]
    if b.merge1:
        zc, uc, vc, _u2, _mis = S20.counts(ds, i, j)
        usable = {w0 + k: bool(v) for k, v in enumerate(ok)}
        runs = S20.join(runs, usable, pos, min_bp, uc, vc, b.mult1, b, zc)
    out = []
    for u, v in runs:
        if p.window1 and not (int(pos[WORD * v + WORD - 1] - pos[WORD * u])
                              > min_bp // WINDOW_FRACTION):
            continue
        out = sc._emit(out, sc.left_end(u), sc.right_end(v), pos, min_bp)
    return out


def ibd2_23(sc, ds, i, j, p, pos, min_bp):
    """`seg21.ibd2_21` with the gate window required to span `min_bp / 2`."""
    n = sc.n
    if n == 0:
        return []
    q = p.base
    w0, w1 = sc.w0, sc.w1
    b = q.base
    ibs0, ibs1, n0, n1, cum = S19.masks(ds, i, j)
    z = [int(n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < b.base.dirty for k in range(n)]
    head, tail = S19.fringe_masks(sc, ibs1, ibs0, b.base)

    def ge_of(e):
        return e + 1 if (e + 1 < n and not z[e + 1] and mis[e + 1]) else e

    def gate_ok(g, e):
        return int(cum[w0 + ge_of(e) + 1] - cum[w0 + g]) >= b.base.gate

    def wide_ok(g, e):
        """The window bound: `pos[last marker of ge_of(e)] - pos[first of gs] >= L/2`."""
        if not p.window:
            return True
        lo = WORD * (w0 + g)
        hi = WORD * (w0 + ge_of(e)) + WORD - 1
        return int(pos[hi] - pos[lo]) >= min_bp // WINDOW_FRACTION

    # the `17-…` §14 bridge, unchanged
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

    def gs_of(a, e):
        return next((t for t in range(a, e + 1) if mis[t] == 0), None)

    runs = []
    for a, e in E._runs(np.array(ok)):
        g = gs_of(a, e)
        if g is None or not gate_ok(g, e):
            continue
        if p.pre_merge and not wide_ok(g, e):
            continue
        runs.append((a, e))

    if q.base.merge2:
        badc, hc, uc = S21.counts2(ds, i, j)
        joined = []
        for a, e in runs:
            if joined:
                pa, pe = joined[-1]
                qq = ge_of(pe) if q.reach else pe
                g2 = gs_of(a, e) if q.reach else a
                mid = [t for t in range(qq + 1, a) if not ok[t]
                       and not (q.reach and t > 0 and ok[t - 1] and not z[t])]
                bad = int(sum(badc[w0 + t] for t in mid))
                H = int(sum(hc[w0 + t] for t in mid))
                U = int(sum(uc[w0 + t] for t in mid))
                X = (H if H >= S21.MERGE_GATE else U) if q.hethet else H + U
                gap = int(pos[WORD * (w0 + g2)] - pos[WORD * (w0 + qq + 1) - 1])
                cap = 10 ** 9 if q.no_cap else q.base.words
                if (any(not ok[t] for t in range(pe + 1, a)) and len(mid) <= cap
                        and gap < min_bp
                        and S21.MERGE_MULT2 * max(0, bad - S21.MERGE_FREE) <= X):
                    joined[-1] = (pa, e)
                    continue
            joined.append((a, e))
        runs = joined

    head_stop = (WORD * (w0 - 1) + S19._last(head) + b.base.fringe_off) if head else sc.lo
    tail_stop = (WORD * (w1 + 1) + S19._first(tail) - b.base.fringe_off) if tail else sc.hi
    head_stop = max(head_stop, sc.lo)
    tail_stop = min(tail_stop, sc.hi)

    out, armed = [], False
    for a, e in runs:
        u, v = w0 + a, w0 + e
        left = WORD * u
        if a > 0 and not z[a - 1] and int(ibs1[u - 1]):
            left = max(0, WORD * (u - 1) + S19._last(int(ibs1[u - 1])) - b.base.reach)
            if a < 2 or z[a - 2]:
                left = max(left, WORD * (u - 1))
        if left <= WORD * w0:
            left = head_stop
        right = WORD * v + WORD - 1
        if e + 1 < n and not z[e + 1] and int(ibs1[v + 1]):
            right = WORD * (v + 1) + S19._first(int(ibs1[v + 1])) + b.base.reach
            if e + 2 >= n or z[e + 2]:
                right = min(right, WORD * (v + 2) - 1)
        if right >= WORD * (w1 + 1) - 1:
            right = tail_stop
        gs = gs_of(a, e)
        if gs is None or not gate_ok(gs, e) or not wide_ok(gs, e):
            continue
        left, right = max(left, sc.lo), min(right, sc.hi)
        if armed:
            left = max(left, WORD * (w0 + gs + 1))
        if out:
            left = max(left, out[-1][1] + b.base.clip)
        if left > right:
            continue
        gsm = min(max(WORD * (w0 + gs), sc.lo), right)
        armed = armed or int(pos[right] - pos[gsm]) >= min_bp // S21.PUSH_FRACTION
        if int(pos[right] - pos[left]) >= min_bp:
            out.append((left, right))
    return out


def call_pair(ds, i, j, p, min_bp=E.SEGLEN):
    pos = ds.pos
    ibd1 = ibd2 = longest = 0
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = ibd2_23(sc, ds, i, j, p, pos, min_bp)
        c1 = ibd1_23(sc, ds, i, j, pos, min_bp, p)
        for lo, hi in c2:
            ibd2 += int(pos[hi] - pos[lo])
        for lo, hi in c2 + c1:
            longest = max(longest, int(pos[hi] - pos[lo]))
        for lo, hi in c1:
            ibd1 += sum(w for w in (int(pos[y] - pos[x])
                                    for x, y in E._pieces((lo, hi), c2))
                        if w >= min_bp)
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


def predict(st, pos, min_bp, p=R23()):
    """The IBD2 calls of `23-…` over per-word composition stats — the canvas mirror.

    `st[k]` is a dict with `z`, `mis`, `mf`, `ml`, `U2` (= `inf2`), `U1` (= A1A1/A1A1);
    `mergelab.wstat` produces exactly that.  Identical to `seg21.predict` apart from the
    window bound.
    """
    n = len(st)
    lo, hi = 0, WORD * n - 1
    z = [st[k]["z"] > 0 for k in range(n)]
    mis = [st[k]["mis"] for k in range(n)]
    i2 = [st[k]["U2"] for k in range(n)]
    usable = [(not z[k]) and mis[k] < 2 for k in range(n)]

    def ge_of(e):
        return e + 1 if (e + 1 < n and not z[e + 1] and mis[e + 1]) else e

    def gate_ok(g, e):
        return sum(i2[t] for t in range(g, ge_of(e) + 1)) >= S21.MERGE_GATE

    def wide_ok(g, e):
        if not p.window:
            return True
        return (pos[min(WORD * ge_of(e) + WORD - 1, hi)] - pos[WORD * g]
                >= min_bp // WINDOW_FRACTION)

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

    runs, k = [], 0
    while k < n:
        if not ok[k]:
            k += 1
            continue
        a = k
        while k < n and ok[k]:
            k += 1
        runs.append((a, k - 1))

    def gs_of(a, e):
        return next((t for t in range(a, e + 1) if mis[t] == 0), None)

    kept = []
    for a, e in runs:
        g = gs_of(a, e)
        if g is None or not gate_ok(g, e):
            continue
        if p.pre_merge and not wide_ok(g, e):
            continue
        kept.append((a, e))

    joined = []
    for a, e in kept:
        if joined:
            pa, pe = joined[-1]
            qq = ge_of(pe)
            g2 = gs_of(a, e)
            mid = [t for t in range(qq + 1, a) if not ok[t]
                   and not (t > 0 and ok[t - 1] and not z[t])]
            bad = sum(st[t]["z"] + st[t]["mis"] for t in mid)
            H = sum(st[t]["U2"] - st[t]["U1"] for t in mid)
            U = sum(st[t]["U1"] for t in mid)
            X = H if H >= S21.MERGE_GATE else U
            if (any(not ok[t] for t in range(pe + 1, a))
                    and pos[WORD * g2] - pos[WORD * (qq + 1) - 1] < min_bp
                    and S21.MERGE_MULT2 * max(0, bad - S21.MERGE_FREE) <= X):
                joined[-1] = (pa, e)
                continue
        joined.append((a, e))

    out, armed = [], False
    for a, e in joined:
        left = WORD * a
        if a > 0 and not z[a - 1] and st[a - 1]["ml"] is not None:
            left = max(0, WORD * (a - 1) + st[a - 1]["ml"] - 63)
            if a < 2 or z[a - 2]:
                left = max(left, WORD * (a - 1))
        right = WORD * e + WORD - 1
        if e + 1 < n and not z[e + 1] and st[e + 1]["mf"] is not None:
            right = WORD * (e + 1) + st[e + 1]["mf"] + 63
            if e + 2 >= n or z[e + 2]:
                right = min(right, WORD * (e + 2) - 1)
        gs = gs_of(a, e)
        if gs is None or not gate_ok(gs, e) or not wide_ok(gs, e):
            continue
        if armed:
            left = max(left, WORD * (gs + 1))
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1])
        if left > right:
            continue
        armed = armed or (pos[right] - pos[min(WORD * gs, right)]
                          >= min_bp // S21.PUSH_FRACTION)
        if pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


FLOORS = S19.FLOORS


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    for bp, sfx in FLOORS:
        print("--seglength %d Mb" % (bp // 1_000_000))
        S20.show("  21 (committed)", S21.score(S21.R21(), bp, sfx))
        S20.show("  23 (window bound)", score(R23(), bp, sfx))
        if mode == "grid":
            for f in ("window", "window1", "span_all"):
                S20.show("  23 without %s" % f, score(replace(R23(), **{f: False}), bp, sfx))
            S20.show("  23 window pre-merge", score(R23(pre_merge=True), bp, sfx))
