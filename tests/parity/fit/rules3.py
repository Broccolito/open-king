"""The fitted segment caller of `rules2.py` plus the informativeness gate.

The gate is measured in `docs/research/fixtures/gate{2,3,4,5,6}.py` and written up in
`docs/research/13-informativeness-gate.md`.  In one line: a run of clean words is only
callable when its **own complete words** carry at least `C` markers at which the pair
shares the A1 allele, and (for IBD1) at least one of them is homozygous for it.

    inf1(marker) = carriesA1_i & carriesA1_j & (hom_i | hom_j)      # IBD1
    inf2(marker) = carriesA1_i & carriesA1_j                        # IBD2
    a run [u..v] is callable iff  popcount(inf & words u..v) >= C,  C = 10

`carriesA1` and `hom` are the two PLINK bit planes `kingdata.read_bed_planes` builds, so
`inf1 = p1_i & p1_j & (p0_i | p0_j)` and `inf2 = p1_i & p1_j` are single word ops.
"""

from dataclasses import dataclass

import numpy as np

import rules2 as R2

WORD = 64
PC = np.bitwise_count


@dataclass(frozen=True)
class P3(R2.P):
    gate: int = 10          # minimum informative markers under an IBD1 run
    gate2: int = -1         # same for IBD2 runs; -1 means "use `gate`"
    gate2_mask: str = "share"   # share = both carry A1 | hom = the IBD1 mask
    gate_scope: str = "core"    # core | ext : count over the run's own words, or over
                                # every word the reported segment touches


_CACHE = {}


def _partial(mask, lo, hi):
    """Popcount of `mask` over marker indices [lo, hi), which may straddle words."""
    n = 0
    a, b = lo, hi
    while a < b:
        w = a // WORD
        end = min(b, WORD * (w + 1))
        bits = int(mask[w]) >> (a - WORD * w)
        bits &= (1 << (end - a)) - 1
        n += bin(bits).count("1")
        a = end
    return n


def inf_masks(ds, i, j):
    """(IBD1 informative mask, IBD2 informative mask) per word, memoised per pair."""
    key = (ds.name, i, j)
    v = _CACHE.get(key)
    if v is None:
        p0i, p1i = ds.p0[i], ds.p1[i]
        p0j, p1j = ds.p0[j], ds.p1[j]
        share = p1i & p1j                 # both carry the A1 allele
        v = (share & (p0i | p0j), share)
        _CACHE[key] = v
    return v


def call_pair(ds, i, j, p=P3(), want=False):
    """`rules2.call_pair` with the informativeness gate applied to each run."""
    ibs0, n0, n1 = R2.counts(ds, i, j)
    m1, m2 = inf_masks(ds, i, j)
    c1, c2 = PC(m1).astype(np.int64), PC(m2).astype(np.int64)
    k1 = np.concatenate(([0], np.cumsum(c1)))
    k2 = np.concatenate(([0], np.cumsum(c2 if p.gate2_mask == "share" else c1)))
    g1 = p.gate
    g2 = p.gate if p.gate2 < 0 else p.gate2
    pos = ds.pos
    ibd1_bp = ibd2_bp = 0
    longest = 0
    detail = []
    for _, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        a0 = n0[w0:w1 + 1]
        a1 = n1[w0:w1 + 1]
        ok1 = a0 <= p.t1
        ok2 = ok1 & (a1 <= p.t2)
        if p.bridge1:
            ok1 = R2._bridge(ok1)
        if p.bridge2:
            ok2 = R2._bridge(ok2)
        span = {}
        for tag, ok, minrun, cum, gate in ((2, ok2, p.min2, k2, g2),
                                           (1, ok1, p.min1, k1, g1)):
            keep = []
            for a, b in R2._runs(ok):
                if b - a + 1 < minrun:
                    continue
                u, v = w0 + a, w0 + b
                lo_m, hi_m = R2._bounds(ibs0, u, v, w0, w1, lo, hi, p.edge,
                                        p.end2 if tag == 2 else "next")
                if gate:
                    mk = m1 if tag == 1 or p.gate2_mask != "share" else m2
                    if p.gate_scope == "core":
                        got = int(cum[v + 1] - cum[u])
                    elif p.gate_scope == "fringe":
                        # the run's own words, plus the usable segment's own markers in
                        # the incomplete word beyond each end the run actually touches
                        got = int(cum[v + 1] - cum[u])
                        if u == w0 and lo < WORD * w0:
                            got += _partial(mk, lo, WORD * w0)
                        if v == w1 and hi > WORD * (w1 + 1) - 1:
                            got += _partial(mk, WORD * (w1 + 1), hi + 1)
                    else:
                        got = int(cum[min(hi_m // WORD + 1, cum.size - 1)]
                                  - cum[lo_m // WORD])
                    if got < gate:
                        continue
                if keep and lo_m <= keep[-1][1]:
                    lo_m = keep[-1][1] + 1
                    if lo_m > hi_m:
                        continue
                ln = int(pos[hi_m] - pos[lo_m])
                if ln < p.seglength_bp:
                    continue
                keep.append((lo_m, hi_m, ln))
                longest = max(longest, ln)
            span[tag] = keep
        for lo_m, hi_m, ln in span[2]:
            ibd2_bp += ln
            if want:
                detail.append((2, lo_m, hi_m, ln))
        for lo_m, hi_m, ln in span[1]:
            ov = 0
            for a, b, _ln in span[2]:
                x, y = max(lo_m, a), min(hi_m, b)
                if x < y:
                    ov += int(pos[y] - pos[x])
            ibd1_bp += ln - ov
            if want:
                detail.append((1, lo_m, hi_m, ln - ov))
    if want:
        return ibd1_bp, ibd2_bp, longest, detail
    return ibd1_bp, ibd2_bp, longest
