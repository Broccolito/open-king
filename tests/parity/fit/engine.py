"""A faithful Python mirror of `crates/king-core/src/ibdseg.rs`, made variable.

Why this exists: the remaining parity gap is *boundary geometry* on a minority of called
segments, and iterating on it in Rust costs a rebuild plus a full corpus replay per idea.
This module reimplements the committed engine exactly — same word rule, same bridging,
same informativeness gate, same asymmetric refinement — with every disputed decision
exposed as a field of `Params`, so a candidate rule can be scored against the whole
corpus in a second.

**It is a mirror, not a second source of truth.** `check_mirror.py` asserts that with
`Params()` (the committed defaults) it reproduces the Rust binary's own `.seg` columns on
all 982 corpus rows and the Rust binary's `MaxIBD2` on all 158 non-zero rows. Any
divergence there is a bug in this file, never a discovery.

Two rulers are implemented over one caller, exactly as `analysis/segments.rs` describes:

* `.seg` measures a call from its refined `lo` to its refined `hi`;
* `--ibs`'s `MaxIBD2` measures the same call **word-aligned**, `pos[64e+63] - pos[64u]`,
  recovering `(u, e)` from the refined endpoints by integer division.
"""

from dataclasses import dataclass

import numpy as np

import kingdata as kd

WORD = 64
PC = np.bitwise_count

SEGLEN = 3_000_000
LONG = 10_000_000


@dataclass(frozen=True)
class Params:
    """Every knob the geometry investigation wants to turn.

    Defaults are the committed engine. Anything not `None`/default here changes the
    rule, so a scorecard line always names what it changed.
    """

    # --- word predicates -------------------------------------------------
    ibd2_dirty_ibs1: int = 5        # IBS1 count that makes a word too dirty for IBD2
    bridge: bool = True             # a lone dirty word between clean ones is absorbed
    gate: int = 10                  # MIN_INFORMATIVE
    min_run1: int = 1
    min_run2: int = 1

    # --- IBD2 geometry ---------------------------------------------------
    ibd2_tail: int = 2              # `v + tail >= w1` snaps the end to the segment end
    ibd2_ext: int = 1               # words the run reaches past its last clean word
    ibd2_ext_last: bool = True      # ...also when the extending word would be `w1`
    ibd2_start_refine: bool = False  # refine the start by the previous word's last IBS0
    ibd2_geom: str = "word"         # "word" = the committed geometry; "ibd1" = borrow
    #                                 Scan.left_end / Scan.right_end wholesale
    ibs_pad: int = 0                # words `--ibs`'s ruler adds past the call's own end

    # --- ordering --------------------------------------------------------
    clip_before_len: bool = True    # clip against the previous call before the length test
    ibd1_clip_ibd2: bool = False    # IBD1 calls are clipped off IBD2 territory, not
    #                                 just subtracted from the total

    # --- marker-level boundary refinement --------------------------------
    # Where inside the flanking word a call stops. `last`/`first` name which IBS0 of that
    # word is used and the integer is the offset added to it, so `("last", 0)` is "end on
    # the flanking word's last IBS0" — the committed rule.
    ibd1_right: tuple = ("last", 0)
    ibd1_left: tuple = ("last", 1)
    ibd2_right: tuple = ("last", 0)
    # The fringe rules, for a run that reaches the usable segment's own first/last word.
    fringe_right: tuple = ("first", -1)
    fringe_left: tuple = ("last", 1)

    def label(self):
        d = Params()
        bits = [f"{k}={getattr(self, k)!r}" for k in self.__dataclass_fields__
                if getattr(self, k) != getattr(d, k)]
        return "baseline" if not bits else " ".join(bits)


BASE = Params()


# ---------------------------------------------------------------------------
# per-pair word masks
# ---------------------------------------------------------------------------

_MASKS = {}


def masks(ds, i, j):
    """(ibs0, ibs1, inf1, inf2) word masks plus the popcount vectors.

    `nhh` is the per-word HetHet popcount, which only `--ibs`'s caller reads
    (`Scan::ibd2_words`); the `.seg` caller never looks at it.
    """
    key = (ds.name, i, j)
    v = _MASKS.get(key)
    if v is None:
        p0i, p1i = ds.p0[i], ds.p1[i]
        p0j, p1j = ds.p0[j], ds.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        ibs0 = p0i & p0j & (p1i ^ p1j)
        ibs1 = (het_i & p0j) | (p0i & het_j)
        share = p1i & p1j
        inf1 = share & (p0i | p0j)
        v = (ibs0, PC(ibs0).astype(np.int32), PC(ibs1).astype(np.int32),
             np.concatenate(([0], np.cumsum(PC(inf1).astype(np.int64)))),
             np.concatenate(([0], np.cumsum(PC(share).astype(np.int64)))),
             PC(het_i & het_j).astype(np.int32))
        _MASKS[key] = v
    return v


def _last_bit(m):
    """Index of the highest set bit of a 64-bit mask (mask must be non-zero)."""
    return int(m).bit_length() - 1


def _first_bit(m):
    return (int(m) & -int(m)).bit_length() - 1


# ---------------------------------------------------------------------------
# the caller, one usable segment at a time
# ---------------------------------------------------------------------------

def _runs(ok):
    """Maximal runs of True in a bool array, as (start, stop) inclusive index pairs."""
    d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()))


class SegScan:
    """One pair over one usable segment — the Rust `Scan`, with `Params` applied."""

    def __init__(self, ds, i, j, seg, p=BASE):
        self.ds, self.p = ds, p
        chrom, lo, hi = seg
        self.lo, self.hi = lo, hi
        self.w0 = -(-lo // WORD)
        self.w1 = (hi + 1) // WORD - 1
        self.n = max(0, self.w1 - self.w0 + 1)
        ibs0, n0, n1, k1, k2, nhh = masks(ds, i, j)
        self.ibs0 = ibs0
        self.n0 = n0
        self.n1 = n1
        self.nhh = nhh
        self.cum1, self.cum2 = k1, k2
        # fringe IBS0 masks: the segment's own markers in the words it does not own
        self.head = 0
        self.tail = 0
        if self.n > 0:
            if lo != WORD * self.w0:
                keep = lo - WORD * (self.w0 - 1)
                self.head = int(ibs0[self.w0 - 1]) & ~((1 << keep) - 1)
            if hi != WORD * (self.w1 + 1) - 1:
                keep = hi - WORD * (self.w1 + 1) + 1
                self.tail = int(ibs0[self.w1 + 1]) & ((1 << keep) - 1)

    # --- gates ---------------------------------------------------------
    def informative(self, cum, u, v):
        return int(cum[v + 1] - cum[u]) >= self.p.gate

    # --- IBD1 ----------------------------------------------------------
    @staticmethod
    def _pick(mask, rule):
        """Bit of `mask` named by `rule = (which, offset)`; `mask` must be non-zero."""
        which, off = rule
        return (_last_bit(mask) if which == "last" else _first_bit(mask)) + off

    def right_end(self, v):
        """Right end of an IBD1 run whose last good word is `v` (global word index)."""
        if v + 1 <= self.w1:
            m = int(self.ibs0[v + 1])
            if m == 0:
                return min(WORD * (v + 1) + 63, self.hi)
            return WORD * (v + 1) + self._pick(m, self.p.ibd1_right)
        if self.tail:
            return WORD * (self.w1 + 1) + self._pick(self.tail, self.p.fringe_right)
        return self.hi

    def left_end(self, u):
        if u - 1 >= self.w0:
            m = int(self.ibs0[u - 1])
            if m == 0:
                return WORD * u
            return WORD * (u - 1) + self._pick(m, self.p.ibd1_left)
        if self.head:
            return WORD * (self.w0 - 1) + self._pick(self.head, self.p.fringe_left)
        return self.lo

    def ibd1(self, pos, min_bp):
        p = self.p
        if self.n == 0:
            return []
        ok = self.n0[self.w0:self.w1 + 1] == 0
        out = []
        for a, b in _runs(ok):
            if b - a + 1 < p.min_run1:
                continue
            u, v = self.w0 + a, self.w0 + b
            if not self.informative(self.cum1, u, v):
                continue
            hi = self.right_end(v)
            lo = self.left_end(u)
            out = self._emit(out, lo, hi, pos, min_bp)
        return out

    # --- IBD2 ----------------------------------------------------------
    def ibd2(self, pos, min_bp):
        p = self.p
        if self.n == 0:
            return []
        sl = slice(self.w0, self.w1 + 1)
        clean = (self.n0[sl] == 0) & (self.n1[sl] < p.ibd2_dirty_ibs1)
        ok = clean.copy()
        if p.bridge:
            n0 = self.n0[sl]
            for k in range(1, self.n - 1):
                if not clean[k] and clean[k - 1] and clean[k + 1] and n0[k] == 0:
                    ok[k] = True
        out = []
        for a, b in _runs(ok):
            if b - a + 1 < p.min_run2:
                continue
            u, v = self.w0 + a, self.w0 + b
            if not self.informative(self.cum2, u, v):
                continue
            if p.ibd2_geom == "ibd1":
                lo, hi = self.left_end(u), self.right_end(v)
                out = self._emit(out, lo, hi, pos, min_bp)
                continue
            e = self._ibd2_end_word(v)
            lo = WORD * u if u != self.w0 else self.lo
            if p.ibd2_start_refine and u != self.w0:
                m = int(self.ibs0[u - 1])
                if m:
                    lo = WORD * (u - 1) + _last_bit(m) + 1
            if e == self.w1:
                hi = self.hi
            else:
                m = int(self.ibs0[e])
                hi = WORD * e + (63 if m == 0 else self._pick(m, self.p.ibd2_right))
            out = self._emit(out, lo, hi, pos, min_bp)
        return out

    def _ibd2_end_word(self, v):
        p = self.p
        if v + p.ibd2_tail >= self.w1:
            return self.w1
        e = v + p.ibd2_ext
        if not p.ibd2_ext_last and e >= self.w1:
            e = self.w1 - 1 if self.w1 - 1 >= v else v
        return min(e, self.w1)

    # --- shared emit ---------------------------------------------------
    def _emit(self, out, lo, hi, pos, min_bp):
        p = self.p
        if p.clip_before_len:
            if out:
                lo = max(lo, out[-1][1] + 1)
            if lo <= hi and pos[hi] - pos[lo] >= min_bp:
                out.append((lo, hi))
        else:
            if lo <= hi and pos[hi] - pos[lo] >= min_bp:
                if out:
                    lo = max(lo, out[-1][1] + 1)
                if lo <= hi:
                    out.append((lo, hi))
        return out


# ---------------------------------------------------------------------------
# pair aggregation — the two rulers
# ---------------------------------------------------------------------------

def call_pair(ds, i, j, p=BASE, min_bp=SEGLEN):
    """Returns (ibd1_bp, ibd2_bp, longest_bp, max_ibd2_wordaligned)."""
    pos = ds.pos
    ibd1_bp = ibd2_bp = longest = 0
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        c2 = sc.ibd2(pos, min_bp)
        c1 = sc.ibd1(pos, min_bp)
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2_bp += ln
            longest = max(longest, ln)
        for lo, hi in c1:
            ln = int(pos[hi] - pos[lo])
            longest = max(longest, ln)
            ibd1_bp += ln - _overlap((lo, hi), c2, pos)
    return ibd1_bp, ibd2_bp, longest, max_ibd2(ds, i, j, p)


def _overlap(c, others, pos):
    tot = 0
    for lo, hi in others:
        a, b = max(c[0], lo), min(c[1], hi)
        if a < b:
            tot += int(pos[b] - pos[a])
    return tot


IBS_IBD2_DIRTY = 5      # het-vs-hom mismatches that make a word break an --ibs run
IBS_IBD2_HETHET = 95    # HetHet markers a call's measured words must hold
IBS_IBD2_MIN_WORDS = 3  # words a call's measured interval must span


def ibd2_words(sc):
    """`--ibs`'s own IBD2 caller — the mirror of `Scan::ibd2_words`.

    NOT the `.seg` caller: opposite homozygotes and missing calls are irrelevant here at
    any density, only het-vs-hom mismatches break a run, and the call runs straight
    through IBS0 words that `Scan::ibd2` would split on.  Returns word intervals
    `(lo, hi)` inclusive, in global word indices.
    """
    n = sc.n
    if n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    clean = [int(sc.n1[w0 + k]) < IBS_IBD2_DIRTY for k in range(n)]
    # A lone dirty word between two clean ones is absorbed; two in a row are not. Read
    # from `clean`, never from the running copy, so dirty words cannot chain in.
    ok = list(clean)
    for k in range(1, max(0, n - 1)):
        if not clean[k] and clean[k - 1] and clean[k + 1]:
            ok[k] = True

    out = []
    k = 0
    while k < n:
        if not ok[k]:
            k += 1
            continue
        k0 = k
        while k < n and ok[k]:
            k += 1
        u, v = w0 + k0, w0 + k - 1
        hi = w1 if v + 2 >= w1 else v + 1
        lo = u
        if out:
            lo = max(lo, out[-1][1] + 1)
        if lo > hi or hi + 1 - lo < IBS_IBD2_MIN_WORDS:
            continue
        # The segment's own tail is exempt from the HetHet count.
        if v + 1 < w1:
            if int(sc.nhh[lo:hi + 1].sum()) < IBS_IBD2_HETHET:
                continue
        out.append((lo, hi))
    return out


def max_ibd2(ds, i, j, p=BASE):
    """`--ibs`'s `MaxIBD2`: the longest `ibd2_words` call, measured word-aligned."""
    pos = ds.pos
    best = 0
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        for u, e in ibd2_words(sc):
            best = max(best, int(pos[WORD * e + 63] - pos[WORD * u]))
    return best


def pr_ibd2(ds, i, j, p=BASE):
    """`--ibs`'s `Pr_IBD2`: the word-aligned `ibd2_words` total over `D`.

    A second aggregate over the same calls as `MaxIBD2`, which grades only the longest
    member.  The 10 Mb rule gates the **pair**, not the call: if no single call reaches
    `LONG`, `Pr_IBD2` is 0 however much shorter material was called.
    """
    pos = ds.pos
    tot = 0
    best = 0
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        for u, e in ibd2_words(sc):
            ln = int(pos[WORD * e + 63] - pos[WORD * u])
            tot += ln
            best = max(best, ln)
    if best < LONG:
        return 0.0
    return tot / ds.denom


_PRT = None


def pr_targets():
    """Reference `Pr_IBD2` for every pair whose `MaxIBD2` is non-zero."""
    global _PRT
    if _PRT is not None:
        return _PRT
    import os
    out = []
    base = os.path.join(kd.ROOT, "tests", "parity", "work", "ibs")
    for name in kd.DATASETS:
        ds = kd.load(name)
        idx = {(f, i): k for k, (f, i) in enumerate(ds.fam)}
        for ext, wide in ((".ibs", False), (".ibs0", True)):
            path = os.path.join(base, "ibs_%s%s" % (name, ext))
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                head = fh.readline().split()
                if "Pr_IBD2" not in head:
                    continue
                c, cm = head.index("Pr_IBD2"), head.index("MaxIBD2")
                for line in fh:
                    f = line.split()
                    if f[cm] == "-9" or float(f[cm]) <= 0:
                        continue
                    if wide:
                        i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
                    else:
                        i, j = idx[(f[0], f[1])], idx[(f[0], f[2])]
                    out.append((name, min(i, j), max(i, j), f[c]))
    _PRT = out
    return out


def score_pr(p=BASE):
    ok = 0
    err = 0.0
    tg = pr_targets()
    for name, i, j, want in tg:
        ds = kd.load(name)
        g = pr_ibd2(ds, i, j, p)
        if "%.4f" % g == want:
            ok += 1
        err += g - float(want)
    return ok, len(tg), err / len(tg)


def max_ibd2_words(ds, i, j, p=BASE):
    """Same as `max_ibd2` but returning `(u, e, bp)` of the winning call."""
    pos = ds.pos
    best = (None, None, 0)
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        prev = None
        for lo, hi in sc.ibd2(pos, 0):
            e = min(hi // WORD + p.ibs_pad, sc.w1)
            u = max(lo // WORD, sc.w0)
            if prev is not None:
                u = max(u, prev + 1)
            if u > e:
                continue
            prev = e
            ln = int(pos[WORD * e + 63] - pos[WORD * u])
            if ln > best[2]:
                best = (u, e, ln)
    return best


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

def score_seg(p=BASE, datasets=None, suffix="__ibdseg", min_bp=SEGLEN, per_ds=False):
    """`.seg` scorecard: exact rows on all four printed columns, plus the pair set."""
    rows = exact = both = ibd1ok = ibd2ok = extra = missing = 0
    err = 0.0
    worst = 0.0
    by = {}
    for name in (datasets or kd.DATASETS):
        ds = kd.load(name)
        ref = ds._read_seg(suffix) if suffix != "__ibdseg" else ds.ref
        d = ds.denom
        e = b = r = 0
        for i, j in ds.pairs():
            i1, i2, lg, _ = call_pair(ds, i, j, p, min_bp)
            got = lg >= LONG
            want = (i, j) in ref
            if not want:
                if got:
                    extra += 1
                continue
            if not got:
                missing += 1
                continue
            rows += 1
            r += 1
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = i1 / d, i2 / d
            gp = g2 + g1 / 2
            ok1 = kd.fmt4(g1) == a1
            ok2 = kd.fmt4(g2) == a2
            ibd1ok += ok1
            ibd2ok += ok2
            if ok1 and ok2:
                both += 1
                b += 1
            if ok1 and ok2 and kd.fmt4(gp) == ap and kd.inf_type(g1, g2, gp) == at:
                exact += 1
                e += 1
            err += abs(gp - ap)
            worst = max(worst, abs(gp - ap))
        by[name] = (e, b, r)
    out = dict(rows=rows, exact=exact, both=both, ibd1=ibd1ok, ibd2=ibd2ok,
               extra=extra, missing=missing,
               mae=err / rows if rows else 0.0, worst=worst)
    if per_ds:
        out["by"] = by
    return out


def max_targets(nonzero=True):
    """Reference `MaxIBD2` values: [(dataset, i, j, bp)] over `.ibs` **and** `.ibs0`.

    `.ibs` is the within-family table (`FID ID1 ID2`) and `.ibs0` the between-family one
    (`FID1 ID1 FID2 ID2`); both carry the column and together they cover every pair the
    reference grades.
    """
    import os
    out = []
    base = os.path.join(kd.ROOT, "tests", "parity", "work", "ibs")
    for name in kd.DATASETS:
        ds = kd.load(name)
        idx = {(f, i): k for k, (f, i) in enumerate(ds.fam)}
        for ext, wide in ((".ibs", False), (".ibs0", True)):
            path = os.path.join(base, "ibs_%s%s" % (name, ext))
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                head = fh.readline().split()
                if "MaxIBD2" not in head:
                    continue
                c = head.index("MaxIBD2")
                for line in fh:
                    f = line.split()
                    v = float(f[c])
                    if nonzero and v <= 0:
                        continue
                    if wide:
                        i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
                    else:
                        i, j = idx[(f[0], f[1])], idx[(f[0], f[2])]
                    out.append((name, min(i, j), max(i, j), int(round(v))))
    return out


def score_max(p=BASE, targets=None):
    tg = targets if targets is not None else max_targets()
    ok = 0
    bad = []
    for name, i, j, t in tg:
        ds = kd.load(name)
        g = max_ibd2(ds, i, j, p)
        if g == t:
            ok += 1
        else:
            bad.append((name, i, j, t, g))
    return ok, len(tg), bad


def main():
    import sys
    tg = max_targets()
    print("MaxIBD2 targets:", len(tg))
    ok, n, bad = score_max(BASE, tg)
    print("MaxIBD2 exact: %d/%d" % (ok, n))
    s = score_seg(BASE, per_ds=True)
    print(".seg: exact %(exact)d  both %(both)d  ibd1 %(ibd1)d  ibd2 %(ibd2)d  "
          "of %(rows)d   extra %(extra)d missing %(missing)d  MAE %(mae).5f "
          "worst %(worst).4f" % s)
    for k, v in s["by"].items():
        print("   %-12s exact %3d  both %3d  of %3d" % (k, v[0], v[1], v[2]))
    if "-v" in sys.argv:
        for row in bad:
            print("   MISS %-12s %3d,%-3d want %10d got %10d  d=%+d"
                  % (row[0], row[1], row[2], row[3], row[4], row[4] - row[3]))


if __name__ == "__main__":
    main()
