"""Invert `--ibs`'s MaxIBD2 column and score candidate IBD2 segment rules against it.

`--ibs` prints, for every pair whose kinship clears 2^-3.5 on a map with >= 100 Mb usable,
`MaxIBD2 <bp>.000` — the length in base pairs of **one single segment**, the longest IBD2
segment the same engine called.  158 corpus pairs carry one.  That is the sharpest
instrument available: a candidate rule reproduces the number exactly or it does not.

Two modes:

    python3 maxfit.py invert [ds...]   # enumerate every marker interval matching the
                                       # target and classify its two endpoints
    python3 maxfit.py score  [ds...]   # score rule variants, exact-match counts

Everything is read off the corpus + the reference's own output.  No KING source involved.
"""

import os
import sys
from collections import Counter, defaultdict

import numpy as np

import kingdata as kd

WORD = 64
PC = np.bitwise_count

IBSDIR = os.path.join(kd.ROOT, "tests", "parity", "work", "ibs")


# --------------------------------------------------------------------------- targets

def read_maxibd2(path):
    out = {}
    with open(path) as fh:
        head = next(fh).rstrip("\n").split("\t")
        if "MaxIBD2" not in head:
            return out
        c, i1, i2 = head.index("MaxIBD2"), head.index("ID1"), head.index("ID2")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            v = int(float(f[c]))
            if v > 0:
                out[(f[i1], f[i2])] = v
    return out


def targets(ds):
    """{(i, j): maxibd2_bp} by sample index, for one dataset."""
    m = {}
    for suffix in (".ibs", ".ibs0"):
        p = os.path.join(IBSDIR, f"ibs_{ds.name}{suffix}")
        if os.path.exists(p):
            m.update(read_maxibd2(p))
    idx = {iid: k for k, (fid, iid) in enumerate(ds.fam)}
    return {(min(idx[a], idx[b]), max(idx[a], idx[b])): v for (a, b), v in m.items()}


def all_targets(names=None):
    out = []
    for name in (names or kd.DATASETS):
        ds = kd.load(name)
        for (i, j), t in sorted(targets(ds).items()):
            out.append((ds, i, j, t))
    return out


# ----------------------------------------------------------------------- per-pair data

_C = {}


def counts(ds, i, j):
    """(ibs0 mask, ibs0 count, ibs1 count) per word of the global grid."""
    key = (ds.name, i, j)
    v = _C.get(key)
    if v is None:
        ibs0, ibs1, _, _ = ds.masks(i, j)
        v = (ibs0, PC(ibs0).astype(np.int16), PC(ibs1).astype(np.int16))
        _C[key] = v
    return v


def last_ibs0(m):
    """Highest set bit index of a word mask, or None."""
    return (m.bit_length() - 1) if m else None


def first_ibs0(m):
    return ((m & -m).bit_length() - 1) if m else None


# -------------------------------------------------------------------------- inversion

def matching_intervals(ds, target):
    """Every (lo, hi) marker interval inside a usable autosomal segment with span target."""
    pos = ds.pos
    hits = []
    for _, lo, hi in ds.segs:
        p = pos[lo:hi + 1]
        k = np.searchsorted(p, p + target)
        ok = np.flatnonzero((k < p.size) & (p[np.minimum(k, p.size - 1)] - p == target))
        for a in ok.tolist():
            hits.append((lo + a, lo + int(k[a]), lo, hi))
    return hits


def classify_start(ds, i, j, lo_m, seglo):
    """How could the reference have produced start marker lo_m?  -> list of tags."""
    ibs0, _, _ = counts(ds, i, j)
    tags = []
    if lo_m == seglo:
        tags.append(("segstart", lo_m // WORD + (0 if lo_m % WORD == 0 else 1)))
    if lo_m % WORD == 0:
        u = lo_m // WORD
        tags.append(("aligned", u))
    else:
        u = lo_m // WORD + 1
        w = u - 1
        m = int(ibs0[w]) if w < len(ibs0) else 0
        b = lo_m - WORD * w - 1
        if m and (m >> b) & 1:
            if last_ibs0(m) == b:
                tags.append(("lastIBS0", u))
            if first_ibs0(m) == b:
                tags.append(("firstIBS0", u))
            tags.append(("someIBS0", u))
    return tags


def classify_end(ds, i, j, hi_m, seghi):
    ibs0, _, _ = counts(ds, i, j)
    tags = []
    if hi_m == seghi:
        tags.append(("segend", (hi_m + 1) // WORD - 1))
    if hi_m % WORD == WORD - 1:
        # either v = hi/64 (run ends on its own last word) or v = hi/64 - 1 (+1 word)
        tags.append(("aligned_same", hi_m // WORD))
        tags.append(("aligned_next", hi_m // WORD - 1))
    else:
        w = hi_m // WORD
        m = int(ibs0[w]) if w < len(ibs0) else 0
        b = hi_m - WORD * w
        if m and (m >> b) & 1:
            v = w - 1
            if last_ibs0(m) == b:
                tags.append(("lastIBS0", v))
            if first_ibs0(m) == b:
                tags.append(("firstIBS0", v))
            tags.append(("someIBS0", v))
    return tags


def run_profile(ds, i, j, u, v):
    """IBS0/IBS1 counts of the implied core [u..v] and its two flanking words."""
    _, n0, n1 = counts(ds, i, j)
    core0 = int(n0[u:v + 1].max()) if v >= u else -1
    core1 = int(n1[u:v + 1].max()) if v >= u else -1
    lft = (int(n0[u - 1]), int(n1[u - 1])) if u > 0 else None
    rgt = (int(n0[v + 1]), int(n1[v + 1])) if v + 1 < len(n0) else None
    return core0, core1, lft, rgt


def invert(names=None):
    tg = all_targets(names)
    print(f"{len(tg)} pairs carry a MaxIBD2\n")
    startc, endc = Counter(), Counter()
    corec1 = Counter()
    nloc = 0
    unloc = []
    for ds, i, j, t in tg:
        hits = matching_intervals(ds, t)
        _, n0, n1 = counts(ds, i, j)
        cands = []
        for lo_m, hi_m, seglo, seghi in hits:
            st = classify_start(ds, i, j, lo_m, seglo)
            en = classify_end(ds, i, j, hi_m, seghi)
            for sname, u in st:
                for ename, v in en:
                    if v < u:
                        continue
                    if int(n0[u:v + 1].max()) != 0:
                        continue          # interior must be IBS0-free
                    cands.append((sname, u, ename, v, lo_m, hi_m))
        if not cands:
            unloc.append((ds.name, i, j, t))
            continue
        nloc += 1
        us = {(u, v) for _, u, _, v, _, _ in cands}
        if len(us) == 1:
            u, v = next(iter(us))
            corec1[int(n1[u:v + 1].max())] += 1
        for sname, u, ename, v, _, _ in cands:
            startc[sname] += 1
            endc[ename] += 1
    print("localised (>=1 IBS0-free candidate):", nloc, " unlocalised:", len(unloc))
    print("start tags:", dict(startc))
    print("end tags:  ", dict(endc))
    print("max IBS1 in core, when the core is unique:", sorted(corec1.items()))
    if unloc:
        print("\nunlocalised pairs:")
        for row in unloc:
            print("   ", row)


# ---------------------------------------------------------------------------- callers

def ibd2_segments(ds, i, j, *, t2=4, min_run=1, end="next", start="aligned",
                  edge="word", min_bp=0, bridge=False):
    """Every IBD2 segment (lo_m, hi_m, bp) under one candidate rule.

    `edge="word"` keeps the segment strictly on the complete-word grid of the usable
    segment; `edge="fringe"` lets a run touching an end of the usable segment reach that
    segment's own first/last marker, the way the IBD1 rule does.
    """
    ibs0, n0, n1 = counts(ds, i, j)
    pos = ds.pos
    out = []
    for _, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        ok = (n0[w0:w1 + 1] == 0) & (n1[w0:w1 + 1] <= t2)
        if bridge and ok.size >= 3:
            ok = ok.copy()
            ok[1:-1] |= ~ok[1:-1] & ok[:-2] & ok[2:]
        d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
        runs = list(zip(np.flatnonzero(d == 1).tolist(),
                        (np.flatnonzero(d == -1) - 1).tolist()))
        keep = []
        for a, b in runs:
            if b - a + 1 < min_run:
                continue
            u, v = w0 + a, w0 + b
            # ---- start
            if u == w0:
                lo_m = lo if edge == "fringe" else WORD * w0
            elif start == "refine" and int(ibs0[u - 1]):
                lo_m = WORD * (u - 1) + last_ibs0(int(ibs0[u - 1])) + 1
            else:
                lo_m = WORD * u
            # ---- end
            if v >= w1:
                hi_m = hi if edge == "fringe" else WORD * w1 + 63
            elif end == "same":
                hi_m = WORD * v + 63
            else:
                m = int(ibs0[v + 1])
                hi_m = WORD * (v + 1) + (last_ibs0(m) if (m and start == "refine") else 63)
                if v + 1 == w1 and edge == "fringe":
                    hi_m = max(hi_m, hi)
                hi_m = min(hi_m, WORD * w1 + 63 if edge != "fringe" else hi)
            if keep and lo_m <= keep[-1][1]:
                lo_m = keep[-1][1] + 1
            if lo_m > hi_m:
                continue
            ln = int(pos[hi_m] - pos[lo_m])
            if ln < min_bp:
                continue
            keep.append((lo_m, hi_m, ln))
        out.extend(keep)
    return out


def ibd2_pairbreak(ds, i, j, *, t=5, hard="dirty", endoff=0, skip=2, min_run=1,
                   min_bp=0):
    """IBD2 segments under the two-consecutive-dirty-words break rule.

    A word is *dirty* when its IBS1 (het-vs-hom) count reaches `t`; `hard` says what an
    IBS0 in the word does — "dirty" makes it dirty like any other, "break" makes it break
    on its own.  A break at word boundary (w, w+1) closes the running segment at
    `w + endoff` and opens the next at `w + skip`.
    """
    ibs0, n0, n1 = counts(ds, i, j)
    pos = ds.pos
    out = []
    for _, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        dirty = (n1 >= t) | (n0 > 0)
        cuts = []
        s = w0
        w = w0
        while w < w1:
            solo = hard == "break" and (n0[w] > 0 or n0[w + 1] > 0)
            if solo or (dirty[w] and dirty[w + 1]):
                cuts.append((s, w + endoff))
                s = w + skip
                w += skip
            else:
                w += 1
        cuts.append((s, w1))
        for a, b in cuts:
            if b < a or b - a + 1 < min_run:
                continue
            lo_m, hi_m = WORD * a, WORD * b + 63
            ln = int(pos[hi_m] - pos[lo_m])
            if ln < min_bp:
                continue
            out.append((lo_m, hi_m, ln))
    return out


def longest(ds, i, j, **kw):
    segs = ibd2_segments(ds, i, j, **kw)
    return max((s[2] for s in segs), default=0)


def longest_pb(ds, i, j, **kw):
    segs = ibd2_pairbreak(ds, i, j, **kw)
    return max((s[2] for s in segs), default=0)


def score_pb(names=None):
    tg = all_targets(names)
    print(f"{len(tg)} pairs carry a MaxIBD2")
    print(f"{'rule':<58}{'exact':>7}{'mean|rel|':>11}")
    best = None
    for hard in ("dirty", "break"):
        for t in range(2, 12):
            for endoff in (0, -1):
                for skip in (1, 2):
                    ok = 0
                    err = 0.0
                    for ds, i, j, tt in tg:
                        g = longest_pb(ds, i, j, t=t, hard=hard, endoff=endoff,
                                       skip=skip)
                        ok += (g == tt)
                        err += abs(g - tt) / tt
                    tag = f"hard={hard} t={t} endoff={endoff} skip={skip}"
                    print(f"{tag:<58}{ok:7d}{err / len(tg):11.4f}")
                    if best is None or ok > best[0]:
                        best = (ok, tag)
    print("\nbest:", best)


def score(names=None):
    tg = all_targets(names)
    print(f"{len(tg)} pairs carry a MaxIBD2")
    print(f"{'rule':<64}{'exact':>7}{'mean|rel|':>11}")
    best = None
    for end in ("next", "same"):
        for start in ("aligned", "refine"):
            for edge in ("word", "fringe"):
                for t2 in range(0, 9):
                    ok = 0
                    err = 0.0
                    for ds, i, j, t in tg:
                        g = longest(ds, i, j, t2=t2, end=end, start=start, edge=edge)
                        ok += (g == t)
                        err += abs(g - t) / t
                    tag = f"end={end} start={start} edge={edge} t2={t2}"
                    print(f"{tag:<64}{ok:7d}{err / len(tg):11.4f}")
                    if best is None or ok > best[0]:
                        best = (ok, tag)
    print("\nbest:", best)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "invert"
    ns = sys.argv[2:] or None
    {"invert": invert, "score": score, "pb": score_pb}[mode](ns)
