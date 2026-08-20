"""The fitted segment caller, parameterised for search.

Rewritten from `rules.py` around what the `--ibs` MaxIBD2 column proves (`maxibd2.py`):
segment boundaries live on the **64-marker word grid**, refined only by the IBS0 markers
of the flanking word.  Every MaxIBD2 in the corpus that can be located resolves to an
interval whose ends are `64u` and `64v+63`, and the ones that do not are exactly the ones
whose flanking word carries an IBS0.

One boundary rule serves both segment types (`_bounds`):

    lo = segment start           if the run starts at the first complete word
       = last IBS0 of word u-1 + 1   if that word has one
       = 64u                         otherwise
    hi = segment end             if the run ends at the last complete word
       = last IBS0 of word v+1       if that word has one
       = 64(v+1)+63                  otherwise

so a segment always reaches into the word that ended it, stopping on that word's **last**
IBS0 — and when the word that ended it holds no IBS0 at all (the IBD2 case, where runs end
on het-mismatch, not IBS0) it swallows that word whole.  That single asymmetry is what
makes a parent-offspring pair read exactly `1.0000` while an IBD2 run reads one word longer
than its own words.

Knobs (searched in `search.py`):

* `t1`, `t2` — per-word tolerance: a word joins an IBD1 run if its IBS0 count is `<= t1`,
  an IBD2 run if additionally its IBS1 (het-vs-hom) count is `<= t2`;
* `min1`, `min2` — minimum run length in complete words;
* `bridge1`, `bridge2` — absorb a single bad word flanked by good ones;
* `seglength_bp`, `long_bp` — the reporting floor and the fixed pair filter.
"""

from dataclasses import dataclass

import numpy as np

WORD = 64
PC = np.bitwise_count


@dataclass(frozen=True)
class P:
    t1: int = 0
    min1: int = 2
    bridge1: int = 0
    t2: int = 4
    min2: int = 1
    bridge2: int = 1
    ibd2_mode: str = "count"
    edge: str = "edge"      # edge | clamp | fringe
    end2: str = "next"      # next | same: where an IBD2 run stops when the word that
                            # ended it holds no IBS0 at all
    seglength_bp: int = 3_000_000
    long_bp: int = 10_000_000


def _clz(m):
    return 64 - m.bit_length()


def _runs(ok):
    d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()))


def _bridge(ok):
    """Absorb single bad words that have a good word on both sides."""
    out = ok.copy()
    if ok.size >= 3:
        out[1:-1] |= ~ok[1:-1] & ok[:-2] & ok[2:]
    return out


_CACHE = {}


def counts(ds, i, j):
    """Per-word (IBS0 mask, IBS0 count, IBS1 count), memoised per pair."""
    key = (ds.name, i, j)
    v = _CACHE.get(key)
    if v is None:
        ibs0, ibs1, _, _ = ds.masks(i, j)
        v = (ibs0, PC(ibs0).astype(np.uint8), PC(ibs1).astype(np.uint8))
        _CACHE[key] = v
    return v


def _bounds(ibs0, u, v, w0, w1, lo, hi, edge="edge", end="next"):
    """Marker span of a word run, per the rule in the module docstring.

    `edge` decides what happens where the run touches the end of the usable segment:
    "edge" takes the segment's own end; "fringe"/"fringe_last" additionally pull the end in
    to the segment's own markers in the incomplete flanking word, to the first / last IBS0
    marker there; "clamp" applies the ordinary flanking-word formula to the word outside
    the segment and clamps.
    """
    if u == w0:
        lo_m = lo
        if edge in ("fringe", "fringe_last", "clamp") and u > 0:
            m = int(ibs0[u - 1])
            if edge != "clamp" and lo > WORD * (u - 1):
                m &= ~((1 << (lo - WORD * (u - 1))) - 1)
            if m:
                lo_m = max(lo, WORD * (u - 1) + (63 - _clz(m)) + 1)
    else:
        m = int(ibs0[u - 1])
        lo_m = WORD * (u - 1) + (63 - _clz(m)) + 1 if m else WORD * u
    if v >= w1:
        hi_m = hi
        if edge in ("fringe", "fringe_last", "clamp") and v + 1 < len(ibs0):
            m = int(ibs0[v + 1])
            keep = hi - WORD * (v + 1) + 1
            if edge != "clamp" and keep < 64:
                m &= (1 << keep) - 1 if keep > 0 else 0
            if m:
                hi_m = min(hi, WORD * (v + 1) + ((63 - _clz(m)) if edge == "fringe_last"
                                                 else ((m & -m).bit_length() - 1) - 1))
    else:
        m = int(ibs0[v + 1])
        if m:
            hi_m = WORD * (v + 1) + (63 - _clz(m))
        else:
            hi_m = WORD * (v + 1) + 63 if end == "next" else WORD * v + 63
        hi_m = min(hi_m, hi)
    return lo_m, hi_m


def call_pair(ds, i, j, p=P(), want=False):
    """(ibd1_bp, ibd2_bp, longest_bp[, detail]) for one pair under knob-set `p`."""
    ibs0, n0, n1 = counts(ds, i, j)
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
        if p.ibd2_mode == "count":
            ok2 = ok1 & (a1 <= p.t2)
        elif p.ibd2_mode == "hetbreak":
            # the two-word contingency rule pinned in crates/open-king-core/src/ibdseg.rs
            brk = np.zeros(a1.size, dtype=bool)
            if a1.size > 1:
                brk[:-1] = (a1[:-1] >= 2) & (a1[1:] > 0)
            ok2 = ok1 & ~brk
            ok2[1:] &= ~brk[:-1]
        elif p.ibd2_mode == "pairsum":
            s2 = a1.astype(np.int32).copy()
            s2[:-1] += a1[1:]
            ok2 = ok1 & (s2 <= p.t2)
        else:
            raise ValueError(p.ibd2_mode)
        if p.bridge1:
            ok1 = _bridge(ok1)
        if p.bridge2:
            ok2 = _bridge(ok2)
        span = {}
        for tag, ok, minrun in ((2, ok2, p.min2), (1, ok1, p.min1)):
            keep = []
            for a, b in _runs(ok):
                if b - a + 1 < minrun:
                    continue
                lo_m, hi_m = _bounds(ibs0, w0 + a, w0 + b, w0, w1, lo, hi, p.edge,
                                     p.end2 if tag == 2 else "next")
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
