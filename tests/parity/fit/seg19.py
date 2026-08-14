"""The `.seg` IBD2 caller with the **fringe** clause corrected — `19-ibd2seg-residual.md`.

`17-seg-caller.md` §5 stated the fringe as *"a call that reaches the usable segment's own
first or last complete word runs on to the segment's first or last marker"*, measured from
the corpus alone because the canvas's chr2 is word-aligned and cannot see it.  It is right
about which markers are *reachable* and wrong about where the call stops inside them: the
partial word beyond the segment's word grid is read **marker by marker, and a het-vs-hom
mismatch there ends the call**, exactly as an opposite homozygote ends an IBD1 call in the
same place.

    left  = (last  mismatch in markers [lo, 64*w0 - 1]) + 1   , else lo
    right = (first mismatch in markers [64*(w1+1), hi])  - 1   , else hi

Everything else is `17-…` §7 with §14's bridge, unchanged.

    python3 seg19.py              # the corpus scorecard at 3 / 5 / 10 Mb
    python3 seg19.py grid         # every knob of the fringe swept
    python3 seg19.py rust         # the numbers the Rust port must reproduce
"""

import sys
from dataclasses import dataclass, replace

import numpy as np

import engine as E
import kingdata as kd

WORD = E.WORD
PC = np.bitwise_count

#: `engine.py`'s committed geometry with the `20-…` run merge switched **off** — the
#: IBD1 pass exactly as it stood when this write-up was measured.
PRE20 = replace(E.BASE, merge=False)


@dataclass(frozen=True)
class R19:
    """`17-…` §7 + §14, with the fringe of `19-…` §3."""

    dirty: int = 2           # het-vs-hom mismatches that make a word unusable
    gate: int = 10           # informative markers a run needs
    reach: int = 63          # markers a call reaches past the flanking word's mismatch
    push: str = "always"     # "never" | "always"
    clip: int = 0            # calls may touch (0) or must not (1)
    # --- the fringe, the clause this write-up measures -------------------
    #   "extend"  — `17-…` §5: run to the segment's own first/last marker (retired)
    #   "mis"     — stop one marker short of the nearest het-vs-hom mismatch there
    #   "mis0"    — ...or opposite homozygote, whichever is nearer
    #   "none"    — do not leave the word grid at all
    fringe: str = "mis"
    fringe_off: int = 1      # markers between the call's end and that mismatch
    # What stops the §5 reach when the *second* word out is the partial word beyond the
    # segment's word grid.  "grid" is `17-…` §5 verbatim — the reach is capped at the
    # grid's own edge; "ibs0" caps only on a real opposite homozygote.  Under
    # `reach_fringe="snap"` the two are **provably identical** (`19-…` §5.2) and the
    # corpus confirms it row for row, so §5 is left exactly as it was.
    edge_cap: str = "grid"
    # What happens once an end leaves the word grid and lands in the partial word.
    #   "snap" — the marker scan takes over: the end *is* the fringe stop, whether that
    #            is nearer than the reach or further (`19-…` §5). Subsumes the §3 fringe.
    #   "stop" — the fringe stop only ever shortens the reach.
    #   "free" — the reach is not limited by the fringe at all.
    reach_fringe: str = "snap"


RETIRED = R19(fringe="extend")


_M = {}


def masks(ds, i, j):
    """`(ibs0, ibs1, n0, n1, cum_inf2)` for one pair — `engine.masks`, by another name."""
    key = (ds.name, i, j)
    v = _M.get(key)
    if v is None:
        ibs0, n0, n1, _c1, _c2, _nhh, ibs1, k2s = E.masks(ds, i, j)
        v = (ibs0, ibs1, n0, n1, k2s)
        _M[key] = v
    return v


def _last(mask):
    return int(mask).bit_length() - 1


def _first(mask):
    return (int(mask) & -int(mask)).bit_length() - 1


def fringe_masks(sc, ibs1, ibs0, p):
    """The two partial words' breaking markers, restricted to the segment's own.

    Returns `(head, tail)` 64-bit masks over word `w0-1` and word `w1+1`.  A set bit is a
    marker that stops an IBD2 call: a het-vs-hom mismatch always, and an opposite
    homozygote too under `fringe="mis0"`.
    """
    lo, hi, w0, w1 = sc.lo, sc.hi, sc.w0, sc.w1
    head = tail = 0
    if p.fringe in ("extend", "none"):
        return 0, 0
    if lo != WORD * w0:
        keep = lo - WORD * (w0 - 1)
        m = int(ibs1[w0 - 1])
        if p.fringe == "mis0":
            m |= int(ibs0[w0 - 1])
        head = m & ~((1 << keep) - 1)
    if hi != WORD * (w1 + 1) - 1:
        keep = hi - WORD * (w1 + 1) + 1
        m = int(ibs1[w1 + 1])
        if p.fringe == "mis0":
            m |= int(ibs0[w1 + 1])
        tail = m & ((1 << keep) - 1)
    return head, tail


def ibd2_19(sc, ds, i, j, p, pos, min_bp):
    """The `.seg` IBD2 caller — `17-…` §7 + §14, with the fringe of `19-…` §3."""
    n = sc.n
    if n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    ibs0, ibs1, n0, n1, cum = masks(ds, i, j)
    z = [int(n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < p.dirty for k in range(n)]
    head, tail = fringe_masks(sc, ibs1, ibs0, p)

    def ge_of(b):
        return b + 1 if (b + 1 < n and not z[b + 1] and mis[b + 1]) else b

    def gate_ok(g, b):
        return int(cum[w0 + ge_of(b) + 1] - cum[w0 + g]) >= p.gate

    # --- the §14 bridge: the ordinary gate, asked twice
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

    # The two marker-level stops beyond the word grid (`19-…` §3, §4).
    head_stop = (WORD * (w0 - 1) + _last(head) + p.fringe_off) if head else sc.lo
    tail_stop = (WORD * (w1 + 1) + _first(tail) - p.fringe_off) if tail else sc.hi
    head_stop = max(head_stop, sc.lo)
    tail_stop = min(tail_stop, sc.hi)

    out, emitted = [], 0
    for a, b in E._runs(ok):
        u, v = w0 + a, w0 + b
        left = WORD * u
        if a > 0 and not z[a - 1] and int(ibs1[u - 1]):
            left = max(0, WORD * (u - 1) + _last(int(ibs1[u - 1])) - p.reach)
            capped = (a < 2 or z[a - 2]) if p.edge_cap == "grid" else (a >= 2
                                                                      and z[a - 2])
            if capped:
                left = max(left, WORD * (u - 1))
            elif a < 2 and p.reach_fringe == "stop":
                left = max(left, head_stop)
        if p.reach_fringe == "snap" and left <= WORD * w0:
            left = head_stop
        right = WORD * v + WORD - 1
        if b + 1 < n and not z[b + 1] and int(ibs1[v + 1]):
            right = WORD * (v + 1) + _first(int(ibs1[v + 1])) + p.reach
            capped = ((b + 2 >= n or z[b + 2]) if p.edge_cap == "grid"
                      else (b + 2 < n and z[b + 2]))
            if capped:
                right = min(right, WORD * (v + 2) - 1)
            elif b + 2 >= n and p.reach_fringe == "stop":
                right = min(right, tail_stop)
        if p.reach_fringe == "snap" and right >= WORD * (w1 + 1) - 1:
            right = tail_stop
        gs = next((t for t in range(a, b + 1) if mis[t] == 0), None)
        if gs is None or not gate_ok(gs, b):
            continue
        if emitted and p.push == "always":
            left = max(left, WORD * (w0 + gs + 1))
        emitted += 1
        # --- the fringe (`19-…` §3). The partial word beyond the segment's grid is read
        # marker by marker: the call runs into it only as far as the nearest het-vs-hom
        # mismatch, and to the segment's own end when there is none.
        if p.reach_fringe != "snap" and p.fringe != "none":
            if u == w0:
                left = min(left, head_stop)
            if v == w1:
                right = max(right, tail_stop)
        left, right = max(left, sc.lo), min(right, sc.hi)
        if out:
            left = max(left, out[-1][1] + p.clip)
        if left <= right and pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


# ---------------------------------------------------------------------------
# aggregation and scoring — the `.seg` row, both columns
# ---------------------------------------------------------------------------

def call_pair(ds, i, j, p, min_bp=E.SEGLEN):
    """`(IBD1 bp, IBD2 bp, longest bp)` — `18-…` §6 assembles IBD1 from the pieces."""
    pos = ds.pos
    ibd1 = ibd2 = longest = 0
    for seg in ds.segs:
        # This module is the `19-…` era engine: its IBD1 pass is `engine.py`'s, but
        # pinned **before** the run merge of `20-…`, which is what makes the 5 and 10 Mb
        # rows of its scorecard the "before" column that write-up is measured against.
        # Without the pin it would silently track whatever `E.BASE` becomes.
        sc = E.SegScan(ds, i, j, seg, PRE20)
        if sc.n == 0:
            continue
        c2 = ibd2_19(sc, ds, i, j, p, pos, min_bp)
        c1 = sc.ibd1(pos, min_bp)
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2 += ln
            longest = max(longest, ln)
        for lo, hi in c1:
            longest = max(longest, int(pos[hi] - pos[lo]))
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


def show(tag, s):
    print("%-34s exact %4d  ibd1 %4d  ibd2 %4d  extra %3d  miss %3d  MAE %.6f  "
          "worst %.4f"
          % (tag, s["exact"], s["ibd1"], s["ibd2"], s["extra"], s["missing"],
             s["mae"], s["worst"]))


FLOORS = [(3_000_000, "__ibdseg"), (5_000_000, "__ibdseg_seglength5"),
          (10_000_000, "__ibdseg_seglength10")]


# ---------------------------------------------------------------------------
# the canvas mirror
# ---------------------------------------------------------------------------

def predict(info, w0=0, w1=None, lo=None, hi=None, head=0, tail=0, p=R19()):
    """`.seg` IBD2 calls as marker intervals, over `segcanvas.wordinfo` tuples.

    `info[k] = (ibs0?, mismatches, first mismatch bit, last mismatch bit, inf2)`.
    `head`/`tail` are the **het-vs-hom mismatch** bit masks of the partial words `w0-1`
    and `w1+1`, restricted to the segment's own markers.  Opposite homozygotes there do
    *not* belong in these masks: `19-…` §2 measures that an IBS0 in a fringe does not stop
    an IBD2 call, though one anywhere in a whole word inside the grid disqualifies it.

    This is the exact mirror of `ibd2_19` above and of `Scan::ibd2` in the Rust engine,
    written over word *descriptions* rather than genotypes so the canvas rigs can grade it
    without a fileset.  `fringecanvas.py` §6 runs it against the reference on random
    canvases whose segment does not start on a word boundary.
    """
    n = len(info)
    w1 = n - 1 if w1 is None else w1
    lo = WORD * w0 if lo is None else lo
    hi = WORD * (w1 + 1) - 1 if hi is None else hi
    z = [info[k][0] for k in range(n)]
    m = [info[k][1] for k in range(n)]
    fb = [info[k][2] for k in range(n)]
    lb = [info[k][3] for k in range(n)]
    i2 = [info[k][4] for k in range(n)]
    usable = [(not z[k]) and m[k] < p.dirty for k in range(n)]

    def ge_of(b):
        return b + 1 if (b + 1 <= w1 and not z[b + 1] and m[b + 1]) else b

    def gate_ok(g, b):
        return sum(i2[t] for t in range(g, ge_of(b) + 1)) >= p.gate

    ok = list(usable)
    gs0 = None
    for k in range(w0, w1 + 1):
        if usable[k]:
            if gs0 is None and m[k] == 0:
                gs0 = k
            continue
        bridged = False
        if (gs0 is not None and k > w0 and not z[k] and k + 1 <= w1
                and usable[k + 1] and m[k + 1] == 0):
            b2 = k + 1
            while b2 + 1 <= w1 and usable[b2 + 1]:
                b2 += 1
            bridged = gate_ok(gs0, k - 1) and gate_ok(k + 1, b2)
        if bridged:
            ok[k] = True
        else:
            gs0 = None

    runs, k = [], w0
    while k <= w1:
        if not ok[k]:
            k += 1
            continue
        a = k
        while k <= w1 and ok[k]:
            k += 1
        runs.append((a, k - 1))

    # The two marker-level stops beyond the word grid (`19-…` §1, §3).
    head_stop = max(WORD * (w0 - 1) + _last(head) + p.fringe_off, lo) if head else lo
    tail_stop = min(WORD * (w1 + 1) + _first(tail) - p.fringe_off, hi) if tail else hi

    out, emitted = [], 0
    for a, b in runs:
        left = WORD * a
        if a - 1 >= w0 and not z[a - 1] and lb[a - 1] is not None:
            left = WORD * (a - 1) + lb[a - 1] - p.reach
            if a - 2 < w0 or z[a - 2]:
                left = max(left, WORD * (a - 1))
        # Once an end reaches the grid's own edge the word scan is over and the marker
        # scan takes over, whether that moves the end out or pulls it back in (`19-…` §4).
        if left <= WORD * w0:
            left = head_stop
        right = WORD * b + WORD - 1
        if b + 1 <= w1 and not z[b + 1] and fb[b + 1] is not None:
            right = WORD * (b + 1) + fb[b + 1] + p.reach
            if b + 2 > w1 or z[b + 2]:
                right = min(right, WORD * (b + 2) - 1)
        if right >= WORD * (w1 + 1) - 1:
            right = tail_stop
        gs = next((t for t in range(a, b + 1) if m[t] == 0), None)
        if gs is None or not gate_ok(gs, b):
            continue
        if emitted and p.push == "always":
            left = max(left, WORD * (gs + 1))
        emitted += 1
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1] + p.clip)
        if left <= right:
            out.append((left, right))
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    for bp, sfx in FLOORS:
        print("--seglength %d Mb" % (bp // 1_000_000))
        # NB: this is *not* the committed `17-…`/`18-…` rule — it is this file's geometry
        # with the fringe knob turned off, which also moves where the retired "extend"
        # clause is applied relative to the push. `seg18.py` scores the real baseline
        # (747 / 982 / 896 at 3 Mb); see `19-…` §8.
        show("  fringe=extend (not the 18 rule)", score(RETIRED, bp, sfx))
        show("  19 (fringe=mis)", score(R19(), bp, sfx))
        if mode == "grid":
            for f in ("none", "mis0"):
                show("  fringe=%s" % f, score(replace(R19(), fringe=f), bp, sfx))
            for o in (0, 2):
                show("  fringe_off=%d" % o, score(replace(R19(), fringe_off=o), bp, sfx))
