"""Score candidate `.seg` IBD2 callers from `docs/research/17-seg-caller.md`.

Reuses `engine.py` for everything already verified (bit planes, usable segments, the
IBD1 caller, `PropIBD`/`InfType`) and replaces only `Scan::ibd2`.

    python3 seg17.py            # the scorecard
    python3 seg17.py grid       # the knob-by-knob grid
"""

import sys
from dataclasses import dataclass, replace

import numpy as np

import engine as E
import kingdata as kd

WORD = E.WORD
LONG = E.LONG
SEGLEN = E.SEGLEN
PC = np.bitwise_count


@dataclass(frozen=True)
class R17:
    """The rule of `docs/research/17-seg-caller.md` §6, with every knob exposed."""
    dirty: int = 2          # het-vs-hom mismatches that make a word unusable
    bridge: str = "clean"   # "clean" | an int mismatch cap | "none"
    gate: int = 10          # `inf2` markers, counted from the run's first clean word
    gate_from_clean: bool = True
    gate_right: bool = True  # ...through the words the right extension reaches
    reach: int = 63         # markers the call reaches past the flanking word's mismatch
    push: str = "always"      # "never" | "always" | "ibs0" (2nd+ call after an IBS0 word)
    clip: int = 0           # calls may touch (0) or must not (1)
    tail_snap: int = 0      # v + tail_snap >= w1 snaps the end to the segment end


_M1 = {}


def ibs1_mask(ds, i, j):
    key = (ds.name, i, j)
    v = _M1.get(key)
    if v is None:
        p0i, p1i, p0j, p1j = ds.p0[i], ds.p1[i], ds.p0[j], ds.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        v = (het_i & p0j) | (p0i & het_j)
        _M1[key] = v
    return v


def _runs(ok):
    out, a = [], None
    for k, v in enumerate(ok):
        if v and a is None:
            a = k
        elif not v and a is not None:
            out.append((a, k - 1))
            a = None
    if a is not None:
        out.append((a, len(ok) - 1))
    return out


def ibd2_17(sc, ds, i, j, p, pos, min_bp):
    """The `.seg` IBD2 caller of `17-seg-caller.md`, as marker intervals."""
    n = sc.n
    if n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    m1 = ibs1_mask(ds, i, j)
    z = [int(sc.n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(sc.n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < p.dirty for k in range(n)]
    ok = list(usable)
    for k in range(n):
        if usable[k] or z[k]:
            continue
        if not (usable[k - 1] if k else False):
            continue
        if not (usable[k + 1] if k + 1 < n else False):
            continue
        if p.bridge == "clean":
            # absorbed only if the run picks up cleanly: the next word carries no
            # mismatch at all, and the usable words after it carry `gate` between them.
            if mis[k + 1] != 0:
                continue
            t, acc = k + 1, 0
            while t < n and usable[t]:
                acc += int(sc.cum2[w0 + t + 1] - sc.cum2[w0 + t])
                t += 1
            if acc < p.gate:
                continue
        elif p.bridge == "none" or mis[k] > int(p.bridge):
            continue
        ok[k] = True

    def firstbit(k):
        v = int(m1[w0 + k])
        return (v & -v).bit_length() - 1 if v else None

    def lastbit(k):
        v = int(m1[w0 + k])
        return v.bit_length() - 1 if v else None

    out, emitted = [], 0
    for a, b in _runs(ok):
        u, v = w0 + a, w0 + b
        left = WORD * u
        after_ibs0 = a > 0 and z[a - 1]
        if a > 0 and not z[a - 1] and lastbit(a - 1) is not None:
            left = WORD * (u - 1) + lastbit(a - 1) - p.reach
            if a < 2 or z[a - 2]:
                left = max(left, WORD * (u - 1))
        right = WORD * v + WORD - 1
        if b + 1 < n and not z[b + 1] and firstbit(b + 1) is not None:
            right = WORD * (v + 1) + firstbit(b + 1) + p.reach
            if b + 2 >= n or z[b + 2]:
                right = min(right, WORD * (v + 2) - 1)
        if p.tail_snap and v + p.tail_snap >= w1:
            right = WORD * w1 + WORD - 1
        # --- the gate: `inf2` from the run's first mismatch-free word, through the
        # words the right end reaches into.
        gs = a
        if p.gate_from_clean:
            gs = next((t for t in range(a, b + 1) if mis[t] == 0), None)
            if gs is None:
                continue
        ge = b
        if p.gate_right and right > WORD * v + WORD - 1:
            ge = min(w1 - w0, right // WORD - w0)
        if int(sc.cum2[w0 + ge + 1] - sc.cum2[w0 + gs]) < p.gate:
            continue
        if emitted and (p.push == "always" or (p.push == "ibs0" and after_ibs0)):
            left = max(left, WORD * (w0 + gs + 1))
        emitted += 1
        # A call that reaches the usable segment's own first/last complete word runs on
        # to the segment's first/last marker — the fringe the word grid does not cover.
        # (`dups`' duplicate pair reads `IBD2Seg 1.0000`, not the word-aligned 0.8984.)
        if u == w0:
            left = min(left, sc.lo)
        if v == w1:
            right = max(right, sc.hi)
        left, right = max(left, sc.lo), min(right, sc.hi)
        if out:
            left = max(left, out[-1][1] + p.clip)
        if left <= right and pos[right] - pos[left] >= min_bp:
            out.append((left, right))
    return out


def call_pair(ds, i, j, p, min_bp=SEGLEN):
    pos = ds.pos
    ibd1_bp = ibd2_bp = longest = 0
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = ibd2_17(sc, ds, i, j, p, pos, min_bp)
        c1 = sc.ibd1(pos, min_bp)
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2_bp += ln
            longest = max(longest, ln)
        for lo, hi in c1:
            ln = int(pos[hi] - pos[lo])
            longest = max(longest, ln)
            ibd1_bp += ln - E._overlap((lo, hi), c2, pos)
    return ibd1_bp, ibd2_bp, longest


def score(p, datasets=None, min_bp=SEGLEN, per_ds=False):
    rows = exact = both = i1ok = i2ok = extra = missing = 0
    err = bias = worst = 0.0
    by = {}
    for name in (datasets or kd.DATASETS):
        ds = kd.load(name)
        ref, d = ds.ref, ds.denom
        e = r = 0
        for i, j in ds.pairs():
            a, b, lg = call_pair(ds, i, j, p, min_bp)
            got, want = lg >= LONG, (i, j) in ref
            if not want:
                extra += got
                continue
            if not got:
                missing += 1
                continue
            rows += 1
            r += 1
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = a / d, b / d
            gp = g2 + g1 / 2
            ok1, ok2 = kd.fmt4(g1) == a1, kd.fmt4(g2) == a2
            i1ok += ok1
            i2ok += ok2
            both += ok1 and ok2
            allok = (ok1 and ok2 and kd.fmt4(gp) == ap
                     and kd.inf_type(g1, g2, gp) == at)
            exact += allok
            e += allok
            err += abs(gp - ap)
            bias += gp - ap
            worst = max(worst, abs(gp - ap))
        by[name] = (e, r)
    out = dict(rows=rows, exact=exact, both=both, ibd1=i1ok, ibd2=i2ok, extra=extra,
               missing=missing, mae=err / max(rows, 1), bias=bias / max(rows, 1),
               worst=worst)
    if per_ds:
        out["by"] = by
    return out


def show(tag, s):
    print("%-30s exact %4d  both %4d  ibd2 %4d  extra %3d  miss %3d  MAE %.5f  "
          "bias %+.5f  worst %.4f"
          % (tag, s["exact"], s["both"], s["ibd2"], s["extra"], s["missing"],
             s["mae"], s.get("bias", 0), s["worst"]))


if __name__ == "__main__":
    show("retired geometry", E.score_seg(E.RETIRED))
    show("committed (engine.py)", E.score_seg())
    base = R17()
    show("17 fitted", score(base))
    if len(sys.argv) > 1 and sys.argv[1] == "grid":
        sys.stdout.reconfigure(line_buffering=True)
        for t in (1, 2, 3, 5):
            show("  dirty=%d" % t, score(replace(base, dirty=t)))
        for b in ("clean", "none", 0, 1, 64):
            show("  bridge=%s" % b, score(replace(base, bridge=b)))
        for r in (0, 32, 62, 63, 64):
            show("  reach=%d" % r, score(replace(base, reach=r)))
        for pu in ("never", "always", "ibs0"):
            show("  push=%s" % pu, score(replace(base, push=pu)))
        for c in (0, 1):
            show("  clip=%d" % c, score(replace(base, clip=c)))
        for g in (False, True):
            show("  gate_from_clean=%s" % g, score(replace(base, gate_from_clean=g)))
            show("  gate_right=%s" % g, score(replace(base, gate_right=g)))
        for ts in (0, 1, 2, 3):
            show("  tail_snap=%d" % ts, score(replace(base, tail_snap=ts)))
