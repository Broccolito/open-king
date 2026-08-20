"""Candidate IBD-segment callers, parameterised so the rule space can be searched.

`call_pair` is one knob-set applied to one pair; `Params` is the knob-set.  The default
`Params()` reproduces the Rust engine in `crates/open-king-core/src/ibdseg.rs` exactly (the
`baseline_matches_rust` check in fit.py pins that), so every experiment is a delta from
a known point.
"""

from dataclasses import dataclass
import numpy as np

WORD = 64
MB = 1_000_000


@dataclass(frozen=True)
class Params:
    # --- IBD1 run formation -------------------------------------------------
    min_run1: int = 2            # complete words a zero-IBS0 run must span
    tol1: int = 0                # IBS0 markers tolerated inside an IBD1 word
    # --- IBD2 run formation -------------------------------------------------
    min_run2: int = 2
    ibd2_mode: str = "hetbreak"  # hetbreak | ibs1free | ibs1rate
    ibs1_max: int = 0            # for ibd2_mode == "ibs1rate": IBS1 tolerated per word
    # --- boundaries ---------------------------------------------------------
    boundary: str = "asym"       # asym (Rust) | clip | mid
    # --- reporting ----------------------------------------------------------
    seglength_bp: int = 3 * MB   # --seglength
    long_bp: int = 10 * MB       # the fixed pair filter
    long_on: str = "any"         # any | ibd1 | total  (what must exceed long_bp)
    min_total_bp: int = 0        # extra pair-level floor on IBD1+IBD2 base pairs
    min_prop: float = 0.0        # extra pair-level floor on PropIBD


POPCNT = np.bitwise_count if hasattr(np, "bitwise_count") else None


def _popcount(x):
    if POPCNT is not None:
        return POPCNT(x)
    return np.array([bin(int(v)).count("1") for v in np.atleast_1d(x)])


def _runs(ok):
    """Maximal runs of True in a bool array, as (start, stop_inclusive) pairs."""
    if ok.size == 0:
        return []
    d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
    starts = np.flatnonzero(d == 1)
    stops = np.flatnonzero(d == -1) - 1
    return list(zip(starts.tolist(), stops.tolist()))


class SegScan:
    """One pair over one usable segment: the per-word masks plus the two fringes."""

    def __init__(self, ds, ibs0, ibs1, seg):
        _, lo, hi = seg
        self.lo, self.hi = lo, hi
        self.w0 = -(-lo // WORD)               # first complete word
        self.w1 = (hi + 1) // WORD - 1         # last complete word
        self.n = max(0, self.w1 - self.w0 + 1)
        self.ibs0 = ibs0[self.w0:self.w1 + 1] if self.n else np.zeros(0, dtype=np.uint64)
        self.ibs1 = ibs1[self.w0:self.w1 + 1] if self.n else np.zeros(0, dtype=np.uint64)
        # Fringe markers: in a word the segment only partly owns.
        self.head = np.uint64(0)
        self.tail = np.uint64(0)
        if self.n:
            if lo != WORD * self.w0:
                keep = lo - WORD * (self.w0 - 1)
                self.head = ibs0[self.w0 - 1] & ~np.uint64((1 << keep) - 1)
            if hi != WORD * (self.w1 + 1) - 1:
                keep = hi - WORD * (self.w1 + 1) + 1
                self.tail = ibs0[self.w1 + 1] & np.uint64((1 << keep) - 1)

    def marker(self, k, b):
        return WORD * (self.w0 + k) + b

    def right_end(self, k1, mode):
        if mode == "clip":
            return min(self.marker(k1, 63), self.hi)
        if k1 + 1 < self.n:
            m = int(self.ibs0[k1 + 1])
            if m == 0:
                return self.marker(k1, 63)
            return self.marker(k1 + 1, 63 - _clz(m))
        if int(self.tail):
            return WORD * (self.w1 + 1) + _ctz(int(self.tail)) - 1
        return self.hi

    def left_end(self, k0, mode):
        if mode == "clip":
            return max(self.marker(k0, 0), self.lo)
        if k0 > 0:
            m = int(self.ibs0[k0 - 1])
            if m == 0:
                return self.marker(k0, 0)
            return self.marker(k0 - 1, 63 - _clz(m)) + 1
        if int(self.head):
            return WORD * (self.w0 - 1) + (63 - _clz(int(self.head))) + 1
        return self.lo


def _clz(m):
    return 64 - m.bit_length()


def _ctz(m):
    return (m & -m).bit_length() - 1


def _segments_from_runs(scan, ok, min_run, mode):
    out = []
    for k0, k1 in _runs(ok):
        if k1 - k0 + 1 < min_run:
            continue
        hi = scan.right_end(k1, mode)
        lo = scan.left_end(k0, mode)
        if out:
            lo = max(lo, out[-1][1] + 1)
        if lo <= hi:
            out.append((lo, hi))
    return out


def _ibd2_ok(scan, p):
    z = scan.ibs0 == 0
    if p.ibd2_mode == "ibs1free":
        return z & (scan.ibs1 == 0)
    if p.ibd2_mode == "ibs1rate":
        return z & (_popcount(scan.ibs1) <= p.ibs1_max)
    # hetbreak: the measured two-word contingency table
    n = scan.n
    if n == 0:
        return z
    brk = np.zeros(n, dtype=bool)
    if n > 1:
        brk[:-1] = (_popcount(scan.ibs1[:-1]) >= 2) & (scan.ibs1[1:] != 0)
    ok = z & ~brk
    ok[1:] &= ~brk[:-1]
    return ok


def call_pair(ds, i, j, p=Params(), want_segments=False):
    """Return (ibd1_bp, ibd2_bp, longest_bp[, segment detail]) for one pair."""
    ibs0, ibs1, _, _ = ds.masks(i, j)
    pos = ds.pos
    ibd1_bp = ibd2_bp = 0
    longest = 0
    detail = []
    for seg in ds.segs:
        scan = SegScan(ds, ibs0, ibs1, seg)
        if scan.n == 0:
            continue
        if p.tol1 == 0:
            ok1 = scan.ibs0 == 0
        else:
            ok1 = _popcount(scan.ibs0) <= p.tol1
        c1 = _segments_from_runs(scan, ok1, p.min_run1, p.boundary)
        c2 = _segments_from_runs(scan, _ibd2_ok(scan, p), p.min_run2, p.boundary)
        c1 = [c for c in c1 if pos[c[1]] - pos[c[0]] >= p.seglength_bp]
        c2 = [c for c in c2 if pos[c[1]] - pos[c[0]] >= p.seglength_bp]
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2_bp += ln
            longest = max(longest, ln)
            if want_segments:
                detail.append((seg[0], 2, lo, hi, ln))
        for lo, hi in c1:
            ln = int(pos[hi] - pos[lo])
            longest = max(longest, ln)
            ov = 0
            for a, b in c2:
                x, y = max(lo, a), min(hi, b)
                if x < y:
                    ov += int(pos[y] - pos[x])
            ibd1_bp += ln - ov
            if want_segments:
                detail.append((seg[0], 1, lo, hi, ln - ov))
    if want_segments:
        return ibd1_bp, ibd2_bp, longest, detail
    return ibd1_bp, ibd2_bp, longest


def reported(ibd1_bp, ibd2_bp, longest, denom, p=Params()):
    if p.long_on == "any":
        ok = longest > p.long_bp
    elif p.long_on == "total":
        ok = (ibd1_bp + ibd2_bp) > p.long_bp
    else:
        ok = longest > p.long_bp
    if not ok:
        return False
    if p.min_total_bp and (ibd1_bp + ibd2_bp) < p.min_total_bp:
        return False
    if p.min_prop:
        prop = (ibd2_bp + ibd1_bp / 2.0) / denom
        if prop < p.min_prop:
            return False
    return True
