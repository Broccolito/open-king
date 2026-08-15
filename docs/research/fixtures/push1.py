#!/usr/bin/env python3
"""Fixtures behind `docs/research/21-push-merge.md` — the push, and the IBD2 merge.

`20-seglength-floor.md` §11 left the IBD2 merge open at 56-58 of 60 on random canvases
and named the one-word push of `17-seg-caller.md` §6 as the first suspect.  This rig
separates them.  Everything here is built out of two words:

    CLEAN  64 HetHet          — usable, `inf2` 64, no mismatch
    WALL   64 opposite homs   — unusable, unbridgeable, and never absorbed by anything

so a canvas is a sentence in runs and walls, and the printed `IBD2Seg` decodes to the
number of marker intervals called (`segcanvas.py`'s ruler).  A word-aligned call over
`k` words measures `64k - 1`; a call pushed one word measures `64k - 65`.

Sections
--------
`push()`   §2 — `--seglength` swept past a call's own length, which separates "the push
                is armed by every call" from "by a call at least half the floor long".
`cap()`    §4 — how many unusable words an IBD2 merge may bridge (the answer is: any
                number), against the same sweep on the IBD1 pass (the answer is two).
`window()` §3 — whether the interruption is measured between the two runs or between
                their gate windows.
`stat()`   §5 — whether `X` is `inf2` or the HetHet count.

Nothing here reads KING's source; every number is a reading off the reference binary.

    python3 push1.py             # every section
    python3 push1.py 2 5         # only §2 and §5
"""
from __future__ import annotations

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segcanvas as S  # noqa: E402
import mergelab as M  # noqa: E402  (points S.CACHE at mergelab_measured.json)

WORD = S.WORD
WALL, CLEAN = S.WALL, S.CLEAN

#: One chr2 word, by marker kind counts.  `z` opposite homozygotes, `m` het-vs-hom
#: mismatches, `h` HetHet, `u` A1A1/A1A1, the rest A2A2/A2A2 (inert).
def word(z=0, m=0, h=0, u=0):
    return (["ibs0"] * z + ["ibs1"] * m + ["hethet"] * h + ["hom1"] * u
            + ["zero"] * (64 - z - m - h - u))


def canvas(block, s=20_000, nw2=70, pad=(3, 3), tag=""):
    """A chr2 painted with `block`, walled and padded out to `nw2` words."""
    import hashlib
    import json
    key = hashlib.sha1(json.dumps([S.expand(x) for x in block]).encode()).hexdigest()[:10]
    return M.cv("push1%s_%s_%d_%d" % (tag, key, s, nw2), block, nw2=nw2, spacing=s,
                pad=pad)


def read(c, seglen):
    """`IBD2Seg` as marker intervals at one floor."""
    r = S.many([(c, ("--seglength", "%.6f" % seglen))])[0]
    return round(float(r["row"]["IBD2Seg"]) * c.denom / c.s)


def steps(c, lo, hi, step):
    """Every `--seglength` at which the printed `IBD2Seg` changes, with its new value.

    The step function's jumps are the individual call lengths, so this reads a canvas's
    whole behaviour in one pass.  **Keep `hi` at or under 10**: above roughly 10 Mb the
    reference stops behaving like a floor at all (§7).
    """
    xs = [round(lo + step * i, 4) for i in range(int(round((hi - lo) / step)) + 1)]
    res = S.many([(c, ("--seglength", "%.6f" % x)) for x in xs])
    out, prev = [], None
    for x, r in zip(xs, res):
        raw = r["row"]["IBD2Seg"] if r["row"] else None
        if raw != prev:
            out.append((x, round(float(raw) * c.denom / c.s)))
        prev = raw
    return out


def _show(label, rows):
    print("  %-46s %s" % (label, " ".join("%.2f:%d" % t for t in rows)))


# ---------------------------------------------------------------------------
# §2 — the push
# ---------------------------------------------------------------------------

def push():
    """Runs of `k1`, `k2`, ... clean words, one wall between each, swept over the floor.

    `c1 = 64*k1 - 1` markers; the next call measures `64*k2 - 1` unpushed and 64 markers
    less pushed.  Raising the floor past `c1` kills the first call, and what the second
    one does next says what arms the push.
    """
    print("§2 the push: runs separated by walls, --seglength swept 1.0 .. 9.0")
    print("   (a call's own length in markers is 64k-1; pushed, 64k-65)")
    for ks in [(1, 4), (2, 5), (4, 2, 6), (1, 2, 6), (2, 2, 6), (6, 2, 6), (6, 1, 6)]:
        block = []
        for i, k in enumerate(ks):
            if i:
                block.append(WALL)
            block += [CLEAN] * k
        _show("runs %s  unpushed %s" % (list(ks), [64 * k - 1 for k in ks]),
              steps(canvas(block, nw2=80), 1.0, 9.0, 0.02))
    print("   the un-push threshold is 2 x (call length measured from its gate-start")
    print("   word) in every row; §2.4 of the write-up bisects it to the base pair.")


# ---------------------------------------------------------------------------
# §4 — the width of the interruption
# ---------------------------------------------------------------------------

Z1 = word(z=1, h=63)              # one opposite homozygote: unusable, unbridgeable
RUN1 = {"hom1": 20, "ibs1": 34, "zero": 10}   # IBD1-live, IBD2-dead
INT1 = {"ibs0": 1, "hom1": 20, "ibs1": 43}    # an IBD1 interruption the budget allows


def cap():
    """`j` unusable words between two four-word runs, at a floor above the gap."""
    print("§4 how many unusable words a merge may bridge")
    print("   IBD2: 4 CLEAN, j x (1 IBS0 + 63 HetHet), 4 CLEAN;  L = 9")
    for j in range(1, 8):
        c = canvas([CLEAN] * 4 + [Z1] * j + [CLEAN] * 4)
        m = read(c, 9.0)
        print("     j=%d  gap %5.2f Mb  merged = %3d  ->  %3d  %s"
              % (j, (64 * j + 1) * 0.02, (8 + j) * 64 - 1, m,
                 "MERGED" if m > 100 else "split"))
    print("   IBD1: the same shape on the IBD1 pass, read off IBD1Seg")
    for j in range(1, 8):
        c = M.cv("push1_i1_%d" % j, [RUN1] * 4 + [INT1] * j + [RUN1] * 4,
                 nw2=70, spacing=20_000)
        r = S.many([(c, ("--seglength", "9.000000"))])[0]
        m = S.mk(c, r, 1)
        print("     j=%d  gap %5.2f Mb  ->  %5.0f  %s"
              % (j, (64 * j + 1) * 0.02, m, "MERGED" if m > 100 else "split"))


# ---------------------------------------------------------------------------
# §3 — runs, or gate windows
# ---------------------------------------------------------------------------

MISW = word(m=4, h=60)            # mismatch-only: the word a run's right end reaches into
U1 = word(m=1, h=63)              # usable, but not mismatch-free: not a gate-start word


def window():
    """Where the interruption begins and ends."""
    print("§3 the interruption runs between the two gate windows")
    print("   4 CLEAN, [MISW, Z1], 4 CLEAN — the earlier run reaches into MISW")
    print("     gap from the run's last word 2.58 Mb, from MISW 1.30 Mb:")
    _show("", steps(canvas([CLEAN] * 4 + [MISW, Z1] + [CLEAN] * 4), 1.2, 2.7, 0.01))
    print("   4 CLEAN, [Z1, MISW], 4 CLEAN — nothing to reach into (Z1 carries an IBS0):")
    _show("", steps(canvas([CLEAN] * 4 + [Z1, MISW] + [CLEAN] * 4), 1.2, 2.7, 0.01))
    print("   and the later run's window opens at its gate-start word:")
    for lead, want in [([], "1.30"), ([U1], "2.58"), ([U1, U1], "3.86")]:
        c = canvas([CLEAN] * 4 + [Z1] + lead + [CLEAN] * (4 - len(lead)))
        _show("run2 opens with %d non-mismatch-free words (gap %s Mb)"
              % (len(lead), want), steps(c, 1.0, 4.0, 0.01))


# ---------------------------------------------------------------------------
# §5 — what `X` is
# ---------------------------------------------------------------------------

def stat():
    """A one-word interruption of 10 opposite homozygotes, HetHet against A1A1/A1A1.

    `3 * (10 - 2) = 24`, so `X = inf2` merges whenever `h + u >= 24` and `X = HetHet`
    only when `h >= 24` — except that at `h < 10` the switch of `20-…` §5 hands `X`
    back to `u`, which is what makes the merge region non-convex.
    """
    print("§5 `X` is the HetHet count, with the switch at 10")
    print("   one word: 10 IBS0, h HetHet, u A1A1/A1A1.  3*(10-2) = 24")
    print("   %4s %4s %4s   %-8s %-8s %-8s" % ("h", "u", "h+u", "ref", "HetHet", "inf2"))
    for h, u in [(0, 54), (5, 49), (9, 45), (10, 44), (12, 12), (20, 34), (23, 31),
                 (24, 30), (24, 0), (23, 0), (30, 20)]:
        c = canvas([CLEAN] * 4 + [word(z=10, h=h, u=u)] + [CLEAN] * 4)
        m = read(c, 9.0)
        x = h if h >= 10 else u
        print("   %4d %4d %4d   %-8s %-8s %-8s"
              % (h, u, h + u, "MERGED" if m > 400 else "split",
                 "M" if 24 <= x else ".", "M" if 24 <= h + u else "."))


SECTIONS = {"2": push, "3": window, "4": cap, "5": stat}

if __name__ == "__main__":
    want = sys.argv[1:] or list(SECTIONS)
    for k in want:
        SECTIONS[k]()
        print()
