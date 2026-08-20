"""Score the full `.seg` corpus with the committed IBD1 rule and a candidate IBD2 rule.

IBD1 is mirrored from `crates/open-king-core/src/ibdseg.rs` unchanged (IBS0-free words, run of
>= 2, boundaries refined by the flanking word's last IBS0, usable-segment fringes).  IBD2
is the rule fitted against `--ibs`'s MaxIBD2 in `maxfit.py` / `maxrule.py`:

    dirty(w) = IBS0(w) > 0 or IBS1(w) >= 5
    a single dirty word flanked by clean words does not break a run; two in a row do
    segment = [64u, 64(v+1)+63] on the usable segment's complete-word grid, and the last
    word boundary of a usable segment never breaks, so a run reaching it takes w1

    python3 score_new.py [old|new] [dirty] [tailwindow]
"""

import sys
from collections import Counter

import numpy as np

import kingdata as kd

WORD = 64
PC = np.bitwise_count
POISONED = {"nuclear", "missing", "monomorphic"}

DIRTY_T = 5
TAILW = 2
SEGLEN = 3_000_000
LONG = 10_000_000

_C = {}


def counts(ds, i, j):
    key = (ds.name, i, j)
    v = _C.get(key)
    if v is None:
        ibs0, ibs1, _, _ = ds.masks(i, j)
        v = (ibs0, PC(ibs0).astype(np.int32), PC(ibs1).astype(np.int32))
        _C[key] = v
    return v


def _runs(ok):
    d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()))


def _hi_bit(m):
    return m.bit_length() - 1


def ibd1_calls(ds, i, j, lo, hi, seglen=SEGLEN):
    """The committed engine's IBD1 pass over one usable segment."""
    ibs0, n0, _ = counts(ds, i, j)
    pos = ds.pos
    w0 = -(-lo // WORD)
    w1 = (hi + 1) // WORD - 1
    if w1 < w0:
        return []
    head = 0
    if lo != WORD * w0 and w0 > 0:
        head = int(ibs0[w0 - 1]) & ~((1 << (lo - WORD * (w0 - 1))) - 1)
    tail = 0
    if hi != WORD * (w1 + 1) - 1 and w1 + 1 < len(ibs0):
        keep = hi - WORD * (w1 + 1) + 1
        tail = int(ibs0[w1 + 1]) & ((1 << keep) - 1 if keep < 64 else ~0)
    out = []
    for a, b in _runs(n0[w0:w1 + 1] == 0):
        if b - a + 1 < 2:
            continue
        u, v = w0 + a, w0 + b
        if v < w1:
            m = int(ibs0[v + 1])
            hi_m = min(WORD * (v + 1) + (_hi_bit(m) if m else 63), hi)
        elif tail:
            hi_m = WORD * (w1 + 1) + (tail & -tail).bit_length() - 2
        else:
            hi_m = hi
        if u > w0:
            m = int(ibs0[u - 1])
            lo_m = WORD * (u - 1) + _hi_bit(m) + 1 if m else WORD * u
        elif head:
            lo_m = WORD * (w0 - 1) + _hi_bit(head) + 1
        else:
            lo_m = lo
        if out:
            lo_m = max(lo_m, out[-1][1] + 1)
        if lo_m <= hi_m and pos[hi_m] - pos[lo_m] >= seglen:
            out.append((lo_m, hi_m))
    return out


EDGE = "fringe"


def ibd2_words(ds, i, j, lo, hi, dirty_t=DIRTY_T, tailw=TAILW):
    """The word runs an IBD2 segment covers: (u, e) inclusive, on the global grid.

    Fitted against `--ibs`'s MaxIBD2 (145/158 exact).  `--ibs` reports the *word-aligned*
    span of these same runs, which is why a duplicate pair's MaxIBD2 is the usable
    segment's aligned span while its `.seg` row still reads 1.0000.
    """
    _, n0, n1 = counts(ds, i, j)
    w0 = -(-lo // WORD)
    w1 = (hi + 1) // WORD - 1
    if w1 < w0:
        return []
    ok = ((n0[w0:w1 + 1] == 0) & (n1[w0:w1 + 1] < dirty_t))
    d = np.diff(np.concatenate(([False], ~ok, [False])).astype(np.int8))
    for a, b in zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()):
        if a == b and a > 0 and b + 1 < ok.size:
            ok[a] = True                      # a lone dirty word does not break a run
    out = []
    for a, b in _runs(ok):
        u, v = w0 + a, w0 + b
        out.append((u, w1 if v >= w1 - tailw else v + 1))
    return out


def ibd2_calls_new(ds, i, j, lo, hi, seglen=SEGLEN, dirty_t=DIRTY_T, tailw=TAILW):
    ibs0, _, _ = counts(ds, i, j)
    pos = ds.pos
    w0 = -(-lo // WORD)
    w1 = (hi + 1) // WORD - 1
    out = []
    for u, e in ibd2_words(ds, i, j, lo, hi, dirty_t, tailw):
        lo_m = WORD * u
        hi_m = WORD * e + 63
        if EDGE == "fringe":
            if u == w0:
                lo_m = lo
            if e == w1:
                hi_m = hi
        elif EDGE == "refine":
            if u == w0:
                lo_m = lo
            elif int(ibs0[u - 1]):
                lo_m = WORD * (u - 1) + _hi_bit(int(ibs0[u - 1])) + 1
            if e == w1:
                hi_m = hi
            elif int(ibs0[e]):
                hi_m = WORD * e + _hi_bit(int(ibs0[e]))
        if out:
            lo_m = max(lo_m, out[-1][1] + 1)
        if lo_m <= hi_m and pos[hi_m] - pos[lo_m] >= seglen:
            out.append((lo_m, hi_m))
    return out


def ibd2_calls_old(ds, i, j, lo, hi, seglen=SEGLEN):
    """The committed engine: IBD1's geometry with a zero-tolerance IBS1 word test."""
    ibs0, n0, n1 = counts(ds, i, j)
    pos = ds.pos
    w0 = -(-lo // WORD)
    w1 = (hi + 1) // WORD - 1
    if w1 < w0:
        return []
    tail = 0
    if hi != WORD * (w1 + 1) - 1 and w1 + 1 < len(ibs0):
        keep = hi - WORD * (w1 + 1) + 1
        tail = int(ibs0[w1 + 1]) & ((1 << keep) - 1 if keep < 64 else ~0)
    head = 0
    if lo != WORD * w0 and w0 > 0:
        head = int(ibs0[w0 - 1]) & ~((1 << (lo - WORD * (w0 - 1))) - 1)
    out = []
    for a, b in _runs((n0[w0:w1 + 1] == 0) & (n1[w0:w1 + 1] == 0)):
        u, v = w0 + a, w0 + b
        if v < w1:
            m = int(ibs0[v + 1])
            hi_m = min(WORD * (v + 1) + (_hi_bit(m) if m else 63), hi)
        elif tail:
            hi_m = WORD * (w1 + 1) + (tail & -tail).bit_length() - 2
        else:
            hi_m = hi
        if u > w0:
            m = int(ibs0[u - 1])
            lo_m = WORD * (u - 1) + _hi_bit(m) + 1 if m else WORD * u
        elif head:
            lo_m = WORD * (w0 - 1) + _hi_bit(head) + 1
        else:
            lo_m = lo
        if out:
            lo_m = max(lo_m, out[-1][1] + 1)
        if lo_m <= hi_m and pos[hi_m] - pos[lo_m] >= seglen:
            out.append((lo_m, hi_m))
    return out


def pair(ds, i, j, ibd2fn):
    pos = ds.pos
    ibd1_bp = ibd2_bp = longest = 0
    for _, lo, hi in ds.segs:
        c2 = ibd2fn(ds, i, j, lo, hi)
        c1 = ibd1_calls(ds, i, j, lo, hi)
        for a, b in c2:
            ln = int(pos[b] - pos[a])
            ibd2_bp += ln
            longest = max(longest, ln)
        for a, b in c1:
            ln = int(pos[b] - pos[a])
            longest = max(longest, ln)
            ov = 0
            for x, y in c2:
                p, q = max(a, x), min(b, y)
                if p < q:
                    ov += int(pos[q] - pos[p])
            ibd1_bp += ln - ov
    return ibd1_bp, ibd2_bp, longest


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "new"
    global EDGE
    EDGE = sys.argv[4] if len(sys.argv) > 4 else EDGE
    fn = ibd2_calls_new if which == "new" else ibd2_calls_old
    if len(sys.argv) > 2:
        global DIRTY_T
        DIRTY_T = int(sys.argv[2])
    if len(sys.argv) > 3:
        global TAILW
        TAILW = int(sys.argv[3])
    print(f"IBD2 rule = {which}  dirty>={DIRTY_T} tailwindow={TAILW}\n")
    print(f"{'dataset':<13}{'rows':>6}{'exact':>7}{'IBD1':>6}{'IBD2':>6}"
          f"{'extra':>7}{'miss':>6}{'MAE':>10}{'worst':>9}")
    tot = np.zeros(6, dtype=np.int64)
    allerr = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        rows = ex = o1 = o2 = extra = miss = 0
        errs = []
        for (i, j) in ds.pairs():
            b1, b2, lg = pair(ds, i, j, fn)
            got = lg >= LONG
            ref = ds.ref.get((i, j))
            if ref is None:
                extra += got
                continue
            rows += 1
            if not got:
                miss += 1
                continue
            p1 = kd.fmt4(b1 / ds.denom)
            p2 = kd.fmt4(b2 / ds.denom)
            prop = b2 / ds.denom + b1 / ds.denom / 2
            o1 += p1 == ref[0]
            o2 += p2 == ref[1]
            ex += (p1 == ref[0]) and (p2 == ref[1])
            errs.append(abs(prop - ref[2]))
        allerr += errs
        mae = np.mean(errs) if errs else 0.0
        worst = max(errs) if errs else 0.0
        flag = "  (poisoned)" if name in POISONED else ""
        print(f"{name:<13}{rows:6d}{ex:7d}{o1:6d}{o2:6d}{extra:7d}{miss:6d}"
              f"{mae:10.5f}{worst:9.4f}{flag}")
        tot += np.array([rows, ex, o1, o2, extra, miss])
    print(f"{'ALL':<13}{tot[0]:6d}{tot[1]:7d}{tot[2]:6d}{tot[3]:6d}{tot[4]:7d}"
          f"{tot[5]:6d}{np.mean(allerr):10.5f}{max(allerr):9.4f}")


if __name__ == "__main__":
    main()
