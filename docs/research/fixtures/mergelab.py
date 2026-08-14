#!/usr/bin/env python3
"""Scratch lab for the `--seglength` run merge — `20-seglength-floor.md`.

Reuses `ibd1canvas`'s alphabet and `segcanvas`'s Canvas, with its own answer cache so the
other two campaigns' caches stay exactly as they were.
"""
from __future__ import annotations
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segcanvas as S
import ibd1canvas as I

ROOT = os.path.dirname(os.path.abspath(__file__))
S.CACHE = os.path.join(ROOT, "mergelab_measured.json")
S.JOBS = int(os.environ.get("MERGELAB_JOBS", "10"))

WORD = S.WORD
K, K0, WALL, B, B0, Z = I.K, I.K0, S.WALL, I.B, I.B0, I.Z
CLEAN = S.CLEAN


def cv(name, block, nw2=36, spacing=40_000, **kw):
    return S.Canvas(name, block, nw2=nw2, spacing=spacing, **kw)


def sweep(items, Ls, what=1):
    """items = [(label, canvas)]; Ls = seglengths in Mb. Returns {(label,L): markers}."""
    reqs = [(c, ("--seglength", "%.6f" % L)) for _, c in items for L in Ls]
    res = S.many(reqs)
    out = {}
    k = 0
    for lab, c in items:
        for L in Ls:
            out[(lab, L)] = (S.mk(c, res[k], 1), S.mk(c, res[k], 2))
            k += 1
    return out


# ---------------------------------------------------------------------------
# the model — `20-seglength-floor.md`
# ---------------------------------------------------------------------------

GATE = 10        # informativeness gate, both passes (`13-…` §5)
FREE = 2         # bad markers a merge gets for nothing
MULT1 = 4        # informative markers one further bad marker costs — IBD1
MULT2 = 3        # ...IBD2
WORDS = 2        # unusable words a merge may bridge


def wstat(spec):
    """Every per-word count the two callers and the merge read."""
    ks = S.expand(spec)
    z = [i for i, k in enumerate(ks) if k in ("ibs0", "ibs0b")]
    mis = [i for i, k in enumerate(ks) if k in ("ibs1", "ibs1b", "ibs1c")]
    U1 = sum(1 for k in ks if k == "hom1")
    V1 = sum(1 for k in ks if k in ("ibs1b", "ibs1c"))
    U2 = sum(1 for k in ks if k in ("hethet", "hom1"))
    return dict(z=len(z), zf=z[0] if z else None, zl=z[-1] if z else None,
                mis=len(mis), mf=mis[0] if mis else None, ml=mis[-1] if mis else None,
                U1=U1, V1=V1, inf1=U1 + V1, U2=U2, inf2=U2)


def mergeable(bad, U, V, mult):
    """`mult * (bad - FREE) <= X`, with `X = V if V >= GATE else U` (IBD1) or `U` (IBD2).

    On the IBD2 pass `V` is passed as 0, so `X` is always `U2 = inf2`; its het-vs-A1A1
    markers are counted on the *left*, inside `bad`, instead.
    """
    X = V if V >= GATE else U
    return mult * (bad - FREE) <= X


def _join(runs, st, pos, min_bp, mult, pass2, usable=None, span="unusable"):
    """Merge adjacent *surviving* runs.  Only the **unusable** words between them count
    toward the two-word cap: a run refused by the gate lies inside the interruption and
    is stepped over, so it can never be the endpoint of a merged segment but does not
    stop one either (`20-…` §6)."""
    out = []
    for a, b in runs:
        if out:
            pa, pb = out[-1]
            mid = list(range(pb + 1, a))
            badw = [k for k in mid if not usable[k]]
            if 1 <= len(badw) <= WORDS and \
               pos[WORD * a] - pos[WORD * (pb + 1) - 1] < min_bp:
                g = [st[k] for k in (mid if span == "all" else badw)]
                bad = sum(w["z"] for w in g) + (sum(w["mis"] for w in g) if pass2 else 0)
                U = sum(w["U2" if pass2 else "U1"] for w in g)
                V = 0 if pass2 else sum(w["V1"] for w in g)
                if mergeable(bad, U, V, mult):
                    out[-1] = (pa, b)
                    continue
        out.append((a, b))
    return out


def ibd1(st, pos, min_bp, lo=None, hi=None, span="unusable"):
    """`18-…` §7 with the merge — IBD1 calls as marker intervals."""
    n = len(st)
    lo = 0 if lo is None else lo
    hi = WORD * n - 1 if hi is None else hi
    ok = [st[k]["z"] == 0 for k in range(n)]
    runs, k = [], 0
    while k < n:
        if not ok[k]:
            k += 1
            continue
        a = k
        while k < n and ok[k]:
            k += 1
        runs.append((a, k - 1))
    # The gate runs FIRST: a run that carries under `GATE` informative markers is
    # refused outright and cannot take part in a merge (`20-…` §5).
    runs = [r for r in runs if sum(st[t]["inf1"] for t in range(r[0], r[1] + 1)) >= GATE]
    runs = _join(runs, st, pos, min_bp, MULT1, False, usable=ok, span=span)
    out = []
    for a, b in runs:
        left = lo if a == 0 else (WORD * (a - 1) + st[a - 1]["zl"] + 1
                                 if st[a - 1]["zl"] is not None else WORD * a)
        right = hi if b == n - 1 else (WORD * (b + 1) + st[b + 1]["zl"]
                                       if st[b + 1]["zl"] is not None
                                       else WORD * (b + 2) - 1)
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1] + 1)
        if left <= right and pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


def model1(cv, seglen_bp, span="unusable"):
    """`IBD1Seg`'s chr2 numerator in marker intervals, for an IBD2-free canvas."""
    _, p2 = cv.positions()
    st = [wstat(x) for x in cv.words]
    return sum(p2[b] - p2[a] for a, b in ibd1(st, p2, seglen_bp, span=span)) / cv.s


# ---------------------------------------------------------------------------
# the out-of-sample battery
# ---------------------------------------------------------------------------

S.PAIR.setdefault("ibs1c", [2, 1])       # het vs A1A1, the other way round


def rword(rng):
    """One chr2 word for the IBD1 battery: IBD2-dead (>= 2 mismatches), else free."""
    ks = ["zero"] * WORD
    slots = list(range(WORD))
    rng.shuffle(slots)
    n = 0

    def put(kind, cnt):
        nonlocal n
        for _ in range(cnt):
            if n < WORD:
                ks[slots[n]] = kind
                n += 1

    if rng.random() < 0.45:
        put("ibs0", rng.choice([1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 12, 30, 64]))
    if rng.random() < 0.15:
        put("ibs0b", rng.choice([1, 2, 3]))
    put("hom1", rng.choice([0, 2, 4, 6, 8, 10, 12, 16, 20, 24, 30]))
    put(rng.choice(["ibs1b", "ibs1c"]), rng.choice([0, 0, 2, 4, 8, 9, 10, 11, 12, 16, 20]))
    put("ibs1", rng.choice([2, 2, 3, 6, 12, 20, 34]))
    put("hethet", rng.choice([0, 0, 4, 12, 30]))
    put("miss", rng.choice([0, 0, 0, 2, 6]))
    return ks


def battery(seed, count, seglen, width=12, nw2=70, spacing=20_000, verbose=False,
            span="unusable"):
    """`model1()` against the reference on random IBD2-free canvases at one floor."""
    import random
    rng = random.Random(seed)
    cvs = [cv("mb%d_%d_%d" % (seed, int(seglen * 10), t),
              [rword(rng) for _ in range(width)], nw2=nw2, spacing=spacing)
           for t in range(count)]
    res = S.many([(c, ("--seglength", "%.6f" % seglen)) for c in cvs])
    ok, bad, dirty = 0, [], 0
    for c, r in zip(cvs, res):
        if float(r["row"]["IBD2Seg"]) != 0.0:
            dirty += 1
            continue
        got, want = S.mk(c, r, 1), model1(c, int(round(seglen * 1e6)), span=span)
        if abs(got - want) <= 0.3:
            ok += 1
        else:
            bad.append((c.name, round(got, 1), round(want, 1)))
    if verbose:
        for b in bad:
            print("      miss", b)
    return ok, len(cvs) - dirty, dirty, bad


def ibd2(st, pos, min_bp, lo=None, hi=None, span="unusable"):
    """`19-…`'s IBD2 caller with the merge — the mirror of `seg20.ibd2_20`."""
    n = len(st)
    lo = 0 if lo is None else lo
    hi = WORD * n - 1 if hi is None else hi
    z = [st[k]["z"] > 0 for k in range(n)]
    m = [st[k]["mis"] for k in range(n)]
    i2 = [st[k]["inf2"] for k in range(n)]
    usable = [(not z[k]) and m[k] < 2 for k in range(n)]

    def ge_of(b):
        return b + 1 if (b + 1 < n and not z[b + 1] and m[b + 1]) else b

    def gate_ok(g, b):
        return sum(i2[t] for t in range(g, ge_of(b) + 1)) >= GATE

    ok = list(usable)                       # the `17-…` §14 bridge, unchanged
    gs0 = None
    for k in range(n):
        if usable[k]:
            if gs0 is None and m[k] == 0:
                gs0 = k
            continue
        bridged = False
        if (gs0 is not None and k > 0 and not z[k] and k + 1 < n
                and usable[k + 1] and m[k + 1] == 0):
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

    def gs_of(a, b):
        return next((t for t in range(a, b + 1) if m[t] == 0), None)

    kept = []
    for a, b in runs:
        g = gs_of(a, b)
        if g is not None and gate_ok(g, b):
            kept.append((a, b))
    kept = _join(kept, st, pos, min_bp, MULT2, True, usable=ok, span=span)

    out, emitted = [], 0
    for a, b in kept:
        left = WORD * a
        if a > 0 and not z[a - 1] and st[a - 1]["ml"] is not None:
            left = max(0, WORD * (a - 1) + st[a - 1]["ml"] - 63)
            if a < 2 or z[a - 2]:
                left = max(left, WORD * (a - 1))
        right = WORD * b + WORD - 1
        if b + 1 < n and not z[b + 1] and st[b + 1]["mf"] is not None:
            right = WORD * (b + 1) + st[b + 1]["mf"] + 63
            if b + 2 >= n or z[b + 2]:
                right = min(right, WORD * (b + 2) - 1)
        gs = gs_of(a, b)
        if gs is None or not gate_ok(gs, b):
            continue
        if emitted:
            left = max(left, WORD * (gs + 1))
        emitted += 1
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1])
        if left <= right and pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


def model2(cv, seglen_bp, span="unusable"):
    _, p2 = cv.positions()
    st = [wstat(x) for x in cv.words]
    return sum(p2[b] - p2[a] for a, b in ibd2(st, p2, seglen_bp, span=span)) / cv.s


def rword2(rng):
    """One chr2 word for the IBD2 battery — the `17-…` axes, IBS0 kept sparse."""
    ks = ["zero"] * WORD
    slots = list(range(WORD))
    rng.shuffle(slots)
    n = 0

    def put(kind, cnt):
        nonlocal n
        for _ in range(cnt):
            if n < WORD:
                ks[slots[n]] = kind
                n += 1

    if rng.random() < 0.40:
        put("ibs0", rng.choice([1, 1, 2, 3, 4, 6, 10, 30, 64]))
    put("ibs1", rng.choice([0, 0, 1, 1, 2, 3, 5, 10]))
    if rng.random() < 0.3:
        put(rng.choice(["ibs1b", "ibs1c"]), rng.choice([1, 2, 4]))
    put("hethet", rng.choice([0, 4, 10, 20, 30, 40, 55]))
    put("hom1", rng.choice([0, 0, 4, 10, 20]))
    put("miss", rng.choice([0, 0, 0, 2, 6]))
    return ks


def battery2(seed, count, seglen, width=12, nw2=70, spacing=20_000, span="unusable"):
    """`model2()` against the reference's `IBD2Seg` on random canvases at one floor."""
    import random
    rng = random.Random(seed)
    cvs = [cv("m2b%d_%d_%d" % (seed, int(seglen * 10), t),
              [rword2(rng) for _ in range(width)], nw2=nw2, spacing=spacing)
           for t in range(count)]
    res = S.many([(c, ("--seglength", "%.6f" % seglen)) for c in cvs])
    ok, bad = 0, []
    for c, r in zip(cvs, res):
        got, want = S.mk(c, r, 2), model2(c, int(round(seglen * 1e6)), span=span)
        if abs(got - want) <= 0.3:
            ok += 1
        else:
            bad.append((c.name, round(got, 1), round(want, 1)))
    return ok, len(cvs), bad
