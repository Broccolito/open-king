"""The push/merge corrections of `docs/research/21-push-merge.md`.

`seg20.py` left five parity cases open, all at `--seglength` 5 and 10.  Its residual was
one-sided at 5 Mb (every wrong `IBD2Seg` too low, every wrong `IBD1Seg` too high), which
said the reference merges IBD2 where that caller does not.  This module is the campaign
that measured why.  Four clauses change, all on the **IBD2** pass; the IBD1 pass of
`20-…` is imported unchanged.

    1. THE ONE-WORD PUSH IS CONDITIONAL.  `17-…` §6 has every call after the first in a
       usable segment starting one word later.  It is not every call: a call arms the
       push only when it is at least **half** the floor long, measured from its own
       **gate-start word** rather than from its left end.

           armed |= pos[right] - pos[WORD * gs]  >=  seglength / 2

       Once armed it stays armed for the rest of the segment.  At the default floor this
       is almost always true, which is why `17-…` §6 saw the unconditional version.

    2. THE IBD2 MERGE HAS NO WORD CAP.  `20-…` §3 bisected "at most two unusable words"
       on the IBD1 pass and the IBD2 pass inherited it.  The IBD2 pass has no cap at all:
       fifteen unusable words merge when the gap and the budget allow.

    3. THE INTERRUPTION IS BETWEEN THE TWO GATE WINDOWS, NOT THE TWO RUNS.  A run's gate
       window runs from its gate-start word `gs` through `ge_of(b)`, the one word its
       right end reaches into.  The gap is measured from the end of the earlier window to
       the start of the later one, and the word the earlier run reaches into is not part
       of the interruption.  The same exclusion applies after any usable word, so a
       gate-refused run's own reach word is skipped too.

    4. `X` IS THE HetHet COUNT, WITH `20-…` §5's SWITCH.  `20-…` §7 read `X` as `inf2`
       (HetHet + A1A1/A1A1); its fixtures could not separate the two.  It is

           X = HetHet   if HetHet >= 10   else   A1A1/A1A1

       — the IBD1 pass's own clause, bisected at 9/10 on this pass as well.

Every constant is measured on `docs/research/fixtures/push1.py` and `mergelab.py`
canvases against the KING 2.3.2 reference and validated out of sample on unused seeds.

    python3 seg21.py            # the corpus scorecard at 3 / 5 / 10 Mb, 20 vs 21
    python3 seg21.py grid       # each clause dropped in turn
"""

import sys
from dataclasses import dataclass, replace

import numpy as np

import engine as E
import kingdata as kd
import seg19 as S19
import seg20 as S20

WORD = E.WORD
PC = np.bitwise_count

MERGE_FREE = 2           # bad markers a merge gets for nothing
MERGE_MULT2 = 3          # informative markers one further bad marker costs, IBD2 pass
MERGE_GATE = 10          # the count at which `X` switches from HetHet to A1A1/A1A1
PUSH_FRACTION = 2        # a call arms the push at `seglength / PUSH_FRACTION`


@dataclass(frozen=True)
class R21:
    """`20-…`'s rule plus the four corrections.  Every field is a knob a section drops."""

    base: S20.R20 = S20.R20()
    push_half: bool = True       # clause 1 — the push is armed at half the floor
    no_cap: bool = True          # clause 2 — the IBD2 merge has no word cap
    reach: bool = True           # clause 3 — gate windows, not runs
    hethet: bool = True          # clause 4 — `X` is HetHet with the switch at 10


# ---------------------------------------------------------------------------
# the per-word counts the IBD2 merge reads
# ---------------------------------------------------------------------------

_C2 = {}


def counts2(ds, i, j):
    """`(bad, H, U)` per-word popcount vectors for the IBD2 merge.

    `bad` is opposite homozygotes plus het-vs-hom mismatches (`20-…` §7), `H` is the
    HetHet count and `U` the A1A1/A1A1 count — the two halves of `inf2`.
    """
    key = (ds.name, i, j)
    v = _C2.get(key)
    if v is None:
        p0i, p1i = ds.p0[i], ds.p1[i]
        p0j, p1j = ds.p0[j], ds.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        ibs0 = p0i & p0j & (p1i ^ p1j)
        ibs1 = (het_i & p0j) | (p0i & het_j)
        share = p1i & p1j
        v = ((PC(ibs0) + PC(ibs1)).astype(np.int32),
             PC(share & ~p0i & ~p0j).astype(np.int32),
             PC(share & p0i & p0j).astype(np.int32))
        _C2[key] = v
    return v


# ---------------------------------------------------------------------------
# the IBD2 caller
# ---------------------------------------------------------------------------

def ibd2_21(sc, ds, i, j, p, pos, min_bp):
    """`seg20.ibd2_20` with the four corrections of `21-…`."""
    n = sc.n
    if n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    b = p.base
    ibs0, ibs1, n0, n1, cum = S19.masks(ds, i, j)
    z = [int(n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < b.base.dirty for k in range(n)]
    head, tail = S19.fringe_masks(sc, ibs1, ibs0, b.base)

    def ge_of(e):
        return e + 1 if (e + 1 < n and not z[e + 1] and mis[e + 1]) else e

    def gate_ok(g, e):
        return int(cum[w0 + ge_of(e) + 1] - cum[w0 + g]) >= b.base.gate

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
        if g is not None and gate_ok(g, e):
            runs.append((a, e))

    # --- the merge -------------------------------------------------------
    if p.base.merge2:
        badc, hc, uc = counts2(ds, i, j)
        joined = []
        for a, e in runs:
            if joined:
                pa, pe = joined[-1]
                # clause 3: the earlier run's gate window ends at `ge_of`, the later
                # one's begins at `gs`; the words each window covers are not part of
                # the interruption.
                q = ge_of(pe) if p.reach else pe
                g2 = gs_of(a, e) if p.reach else a
                mid = [t for t in range(q + 1, a) if not ok[t]
                       and not (p.reach and t > 0 and ok[t - 1] and not z[t])]
                bad = int(sum(badc[w0 + t] for t in mid))
                H = int(sum(hc[w0 + t] for t in mid))
                U = int(sum(uc[w0 + t] for t in mid))
                X = (H if H >= MERGE_GATE else U) if p.hethet else H + U
                gap = int(pos[WORD * (w0 + g2)] - pos[WORD * (w0 + q + 1) - 1])
                cap = 10 ** 9 if p.no_cap else p.base.words
                if (any(not ok[t] for t in range(pe + 1, a)) and len(mid) <= cap
                        and gap < min_bp
                        and MERGE_MULT2 * max(0, bad - MERGE_FREE) <= X):
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
        if gs is None or not gate_ok(gs, e):
            continue
        if armed:
            left = max(left, WORD * (w0 + gs + 1))
        left, right = max(left, sc.lo), min(right, sc.hi)
        if out:
            left = max(left, out[-1][1] + b.base.clip)
        if left > right:
            continue
        # clause 1: the push is armed by a call at least half the floor long, measured
        # from its gate-start word.  Unconditional in `20-…`.
        if p.push_half:
            gsm = min(max(WORD * (w0 + gs), sc.lo), right)
            armed = armed or int(pos[right] - pos[gsm]) >= min_bp // PUSH_FRACTION
        else:
            armed = True
        if int(pos[right] - pos[left]) >= min_bp:
            out.append((left, right))
    return out


# ---------------------------------------------------------------------------
# aggregation and scoring — `seg20.call_pair` with the new IBD2 caller
# ---------------------------------------------------------------------------

def call_pair(ds, i, j, p, min_bp=E.SEGLEN):
    pos = ds.pos
    ibd1 = ibd2 = longest = 0
    b = p.base
    q = replace(b, merge1=False, merge2=False) if not b.filter_merged else b
    qp = replace(p, base=q)
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = ibd2_21(sc, ds, i, j, p, pos, min_bp)
        c1 = S20.ibd1_20(sc, ds, i, j, pos, min_bp, b)
        f2 = c2 if q is b else ibd2_21(sc, ds, i, j, qp, pos, min_bp)
        f1 = c1 if q is b else S20.ibd1_20(sc, ds, i, j, pos, min_bp, q)
        for lo, hi in c2:
            ibd2 += int(pos[hi] - pos[lo])
        for lo, hi in f2 + f1:
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


# ---------------------------------------------------------------------------
# the canvas mirror
# ---------------------------------------------------------------------------

def predict(st, pos, min_bp, p=R21()):
    """The IBD2 calls of `21-…` over per-word composition stats — the canvas mirror.

    `st[k]` is a dict with `z`, `mis`, `mf`, `ml`, `U2` (= `inf2`), `U1` (= A1A1/A1A1);
    `mergelab.wstat` produces exactly that.  `lab21.ibd2` in the scratch lab and this
    function are the same rule; this is the one the write-up quotes.
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
        return sum(i2[t] for t in range(g, ge_of(e) + 1)) >= MERGE_GATE

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
        if g is not None and gate_ok(g, e):
            kept.append((a, e))

    joined = []
    for a, e in kept:
        if joined:
            pa, pe = joined[-1]
            q = ge_of(pe)
            g2 = gs_of(a, e)
            mid = [t for t in range(q + 1, a) if not ok[t]
                   and not (t > 0 and ok[t - 1] and not z[t])]
            bad = sum(st[t]["z"] + st[t]["mis"] for t in mid)
            H = sum(st[t]["U2"] - st[t]["U1"] for t in mid)
            U = sum(st[t]["U1"] for t in mid)
            X = H if H >= MERGE_GATE else U
            if (any(not ok[t] for t in range(pe + 1, a))
                    and pos[WORD * g2] - pos[WORD * (q + 1) - 1] < min_bp
                    and MERGE_MULT2 * max(0, bad - MERGE_FREE) <= X):
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
        if gs is None or not gate_ok(gs, e):
            continue
        if armed:
            left = max(left, WORD * (gs + 1))
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1])
        if left > right:
            continue
        armed = armed or pos[right] - pos[min(WORD * gs, right)] >= min_bp // PUSH_FRACTION
        if pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


FLOORS = S19.FLOORS


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    for bp, sfx in FLOORS:
        print("--seglength %d Mb" % (bp // 1_000_000))
        S20.show("  20 (committed)", S20.score(S20.R20(), bp, sfx))
        S20.show("  21 (push + merge)", score(R21(), bp, sfx))
        if mode == "grid":
            for f in ("push_half", "no_cap", "reach", "hethet"):
                S20.show("  21 without %s" % f,
                         score(replace(R21(), **{f: False}), bp, sfx))
