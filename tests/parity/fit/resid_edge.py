"""Dump the endpoint arithmetic of every still-disagreeing usable segment.

After `19-…` §3 the residual is entirely **under-call**, a handful of markers per segment.
This prints, for each guilty segment, the exact reference total, our calls, and every
quantity the endpoint rule reads at the two ends — the flanking words' mismatch bits, the
second word out, and the partial words beyond the grid — so the offending clause can be
read off rather than guessed.
"""

import sys

import numpy as np

import engine as E
import kingdata as kd
import resid19 as R
import seg19 as S19
import seglocal as SL

WORD = 64


def bitpos(mask, base, lo=None, hi=None):
    m = int(mask)
    out = [base + b for b in range(64) if (m >> b) & 1]
    if lo is not None:
        out = [v for v in out if v >= lo]
    if hi is not None:
        out = [v for v in out if v <= hi]
    return out


def dump(name, n1, n2, k=None):
    ds = kd.load(name)
    idx = {f[1]: t for t, f in enumerate(ds.fam)}
    i, j = sorted((idx[n1], idx[n2]))
    pos, d = ds.pos, ds.denom
    bad, _ = SL.per_segment(name, n1, n2, only=k, quiet=True)
    ibs0, ibs1, n0, n1w, _cum = S19.masks(ds, i, j)
    for t, seg, calls, tot, refv, delta in bad:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        w0, w1, n = sc.w0, sc.w1, sc.n
        ref = float(refv) * d
        print("=== %s %s/%s seg%d chr%d  w%d..%d (n=%d)  lo=%d(w%d+%d) hi=%d(w%d+%d)"
              % (name, ds.fam[i][1], ds.fam[j][1], t, seg[0], w0, w1, n,
                 sc.lo, sc.lo // WORD, sc.lo % WORD, sc.hi, sc.hi // WORD, sc.hi % WORD))
        print("    ref %s = %.0f bp   ours %d bp   d = %+.0f bp (%+.2f markers)"
              % (refv, ref, tot, delta, delta / R._median_gap(ds)))
        shape = "".join("W" if int(n0[w0 + q]) else
                        ("C" if int(n1w[w0 + q]) == 0 else
                         ("x" if int(n1w[w0 + q]) == 1 else "y")) for q in range(n))
        print("    shape %s" % shape)
        print("    head part w%d markers>=%d: mismatches %s  IBS0 %s"
              % (w0 - 1, sc.lo, bitpos(ibs1[w0 - 1], WORD * (w0 - 1), lo=sc.lo),
                 bitpos(ibs0[w0 - 1], WORD * (w0 - 1), lo=sc.lo)))
        print("    tail part w%d markers<=%d: mismatches %s  IBS0 %s"
              % (w1 + 1, sc.hi, bitpos(ibs1[w1 + 1], WORD * (w1 + 1), hi=sc.hi),
                 bitpos(ibs0[w1 + 1], WORD * (w1 + 1), hi=sc.hi)))
        for a, b in calls:
            print("    call [%d,%d]  w%d+%d .. w%d+%d   %.4f Mb"
                  % (a, b, a // WORD, a % WORD, b // WORD, b % WORD,
                     (pos[b] - pos[a]) / 1e6))
        # the two ends of the run that owns each call, spelled out
        ok = _runs_of(sc, ds, i, j)
        for (ra, rb) in ok:
            u, v = w0 + ra, w0 + rb
            print("    run w%d..w%d:" % (u, v))
            if ra > 0:
                print("       left  flank w%d mis@%s ibs0=%d   second-out w%d %s"
                      % (u - 1, bitpos(ibs1[u - 1], 0), int(n0[u - 1]), u - 2,
                         "outside grid" if ra < 2 else
                         ("IBS0" if int(n0[u - 2]) else "clean")))
            if rb < n - 1:
                print("       right flank w%d mis@%s ibs0=%d   second-out w%d %s"
                      % (v + 1, bitpos(ibs1[v + 1], 0), int(n0[v + 1]), v + 2,
                         "outside grid" if rb + 2 >= n else
                         ("IBS0" if int(n0[v + 2]) else "clean")))
        # what a few candidate right/left ends would measure
        _candidates(ds, i, j, sc, calls, tot, ref, ibs1, ibs0, n0, n1w)
        print()


def _runs_of(sc, ds, i, j, p=S19.R19()):
    n = sc.n
    w0 = sc.w0
    _ibs0, _ibs1, n0, n1, cum = S19.masks(ds, i, j)
    z = [int(n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < p.dirty for k in range(n)]

    def ge_of(b):
        return b + 1 if (b + 1 < n and not z[b + 1] and mis[b + 1]) else b

    def gate_ok(g, b):
        return int(cum[w0 + ge_of(b) + 1] - cum[w0 + g]) >= p.gate

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
    return E._runs(ok)


def _candidates(ds, i, j, sc, calls, tot, ref, ibs1, ibs0, n0, n1w):
    """How many markers each end would have to move to close the gap."""
    pos = ds.pos
    need = ref - tot
    if not calls:
        return
    print("    to close %+.0f bp the ends would move:" % need)
    for a, b in calls:
        for step in range(1, 8):
            if a - step >= sc.lo:
                dv = int(pos[a] - pos[a - step])
                if abs(dv - need) < 0.4 * ds.denom / 10000.0:
                    print("       left of [%d,%d] back %d markers to %d (w%d+%d) = %+d bp"
                          % (a, b, step, a - step, (a - step) // WORD,
                             (a - step) % WORD, dv))
            if b + step <= sc.hi:
                dv = int(pos[b + step] - pos[b])
                if abs(dv - need) < 0.4 * ds.denom / 10000.0:
                    print("       right of [%d,%d] on %d markers to %d (w%d+%d) = %+d bp"
                          % (a, b, step, b + step, (b + step) // WORD,
                             (b + step) % WORD, dv))


CASES = [("nuclear", "N_C1", "N_C3"), ("nuclear", "N_C1", "N_C4"),
         ("nuclear", "N_C2", "N_C4"), ("nuclear", "N_C3", "N_C4"),
         ("multifam", "A_F", "B_F"), ("multifam", "A_C1", "A_C2"),
         ("multifam", "A_C1", "C_F"), ("multifam", "A_C2", "C_F"),
         ("multifam", "A_C3", "C_F"), ("multifam", "B_C1", "B_C2"),
         ("multifam", "B_C2", "B_C3"), ("multifam", "C_C1", "C_C3"),
         ("missing", "M_C2", "M_C4"), ("sexchr", "S_SON2", "S_DAU1"),
         ("admixed", "Z_C1", "Z_C2")]


if __name__ == "__main__":
    if len(sys.argv) > 3:
        dump(sys.argv[1], sys.argv[2], sys.argv[3],
             int(sys.argv[4]) if len(sys.argv) > 4 else None)
    else:
        for c in CASES:
            dump(*c)
