#!/usr/bin/env python3
"""Fixtures behind `docs/research/23-gap-bound.md` - the gate window's length bound, and
what the IBD1 merge's budget is summed over.

`21-push-merge.md` section 8.1 left `--seglength 10` wrong on twelve corpus rows and named
the IBD2 merge's *gap* as the suspect.  It is not the gap.  `chrprobe.py` reads the
reference's answer one chromosome at a time on the corpus's own data and localises two
faults, both floor-dependent and neither in the merge:

    1. an IBD2 call the reference stops reporting at a floor **far below its own length**
       - 11.2066 Mb, kept at `--seglength 6.290751` and gone at 6.290752;
    2. an IBD1 merge the reference makes across an interruption that contains a
       gate-refused run, and this caller does not.

This rig measures both on purpose-built canvases, away from the corpus.

Sections
--------
`window()`    2 - sweep `--seglength` against a canvas whose IBD2 call is three words long
                  but whose **gate window** is one word.  The call dies at twice the
                  window, not at its own length.
`ends()`      3 - which two words the window runs between: `gs` on the left (not the run's
                  first word) and `ge_of(e)` on the right (the reach word).
`ibd1w()`     4 - the same bound on the IBD1 pass, over the run's own words, with a
                  **strict** comparison - one unit of `seglength / 2` tighter.
`budget()`    5 - an IBD1 interruption holding a gate-refused run: whether that run's own
                  markers enter the merge budget (`20-...` section 11 item 4).  They do.
`when()`      6 - whether the bound is asked with the gate or at emit.  At emit: a run it
                  refuses still merges with its neighbour.
`held_out()`  7 - random IBD2 canvases on unused seeds, in the two spacings where the
                  bound can decide a call.
`held_out1()` 8 - the same for the IBD1 pass.

Nothing here reads KING's source; every number is a reading off the reference binary.

    python3 window1.py             # every section
    python3 window1.py 2 5         # only sections 2 and 5
"""
from __future__ import annotations

import hashlib
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segcanvas as S  # noqa: E402
import mergelab as M  # noqa: E402  (points S.CACHE at mergelab_measured.json)

WORD = S.WORD
WALL, CLEAN = S.WALL, S.CLEAN


def word(z=0, m=0, h=0, u=0, v=0):
    """`z` opposite homs, `m` het-vs-A2A2, `h` HetHet, `u` A1A1/A1A1, `v` het-vs-A1A1.

    The kinds are laid down in that order, so `m` mismatches with no `z` before them sit
    at bits 0..m-1 - which is what fixes `mf`/`ml` and therefore the endpoint reach.
    """
    return (["ibs0"] * z + ["ibs1"] * m + ["hethet"] * h + ["hom1"] * u
            + ["ibs1b"] * v + ["zero"] * (64 - z - m - h - u - v))


def canvas(block, s=50_000, nw2=45, pad=(3, 3), tag=""):
    key = hashlib.sha1(json.dumps([S.expand(x) for x in block]).encode()).hexdigest()[:10]
    return M.cv("win1%s_%s_%d_%d" % (tag, key, s, nw2), block, nw2=nw2, spacing=s,
                pad=pad)


def read(c, seglen, what=2):
    r = S.many([(c, ("--seglength", "%.6f" % (seglen / 1e6)))])[0]
    v = S.mk(c, r, what)
    return None if v is None else round(v)


def steps(c, lo, hi, step, what=2):
    """Every `--seglength` (in bp) at which the printed column changes."""
    n = int(round((hi - lo) / step))
    xs = [lo + step * k for k in range(n + 1)]
    res = S.many([(c, ("--seglength", "%.6f" % (x / 1e6))) for x in xs])
    out, prev = [], object()
    for x, r in zip(xs, res):
        v = S.mk(c, r, what)
        k = None if v is None else round(v)
        if k != prev:
            out.append((x, k))
        prev = k
    return out


def flip(c, lo, hi, what=2):
    """Bisect the single `--seglength` at which the printed column changes."""
    a, b = read(c, lo, what), read(c, hi, what)
    assert a != b, (a, b)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if read(c, mid, what) == a:
            lo = mid
        else:
            hi = mid
    return lo, hi


# ---------------------------------------------------------------------------
# the words the IBD2 sections are built from
# ---------------------------------------------------------------------------

M0 = word(m=2, h=62)        # unusable (2 mismatches), mismatches at bits 0..1: ml = 1
MISW = word(m=4, h=60)      # unusable, mismatches at bits 0..3: mf = 0
U1 = word(m=1, h=63)        # usable, but not mismatch-free: never a gate-start word


# ---------------------------------------------------------------------------
# §2 — the bound
# ---------------------------------------------------------------------------

def window():
    """One CLEAN word, walled on the right, with two unusable words on its left.

    The left end reaches back past the whole of `M0` into the word before it, so the
    reported call is 189 marker intervals — but the run, and therefore the gate window,
    is one word: 63 intervals.  Under `21-…` the call survives every floor under its own
    length; under the bound it dies at twice the window.
    """
    print("§2 an IBD2 call three times its own gate window")
    c = canvas([M0, M0, CLEAN])
    print("   block [M0, M0, CLEAN]; call = 189 markers = 9.450 Mb, window = 63 = 3.150 Mb")
    print("   %-22s %s" % ("--seglength sweep",
                           " ".join("%.3f:%s" % (x / 1e6, v)
                                    for x, v in steps(c, 1_000_000, 10_000_000, 100_000))))
    lo, hi = flip(c, 1_000_000, 10_000_000)
    print("   bisected: kept at %d bp, dropped at %d bp" % (lo, hi))
    print("   2 x window + 1 = %d   (the test is `window >= seglength / 2`, integer)"
          % (2 * 63 * c.s + 1))


# ---------------------------------------------------------------------------
# §3 — which two words the window runs between
# ---------------------------------------------------------------------------

def ends():
    """`gs` on the left and `ge_of(e)` on the right, each separated from the alternative."""
    print("§3 the window's two ends")
    print("   right end: [M0, M0, CLEAN, MISW] — the run reaches into MISW, so the")
    print("   window is two words (127 = 6.350 Mb) and nothing under 10 Mb can kill it")
    c = canvas([M0, M0, CLEAN, MISW])
    print("   %-22s %s" % ("", " ".join("%.3f:%s" % (x / 1e6, v) for x, v in
                                        steps(c, 1_000_000, 10_000_000, 250_000))))
    print("   left end: [M0, M0, U1, CLEAN] — the run is two words but its gate-start")
    print("   word is the second, so the window is one word again (63 = 3.150 Mb)")
    c = canvas([M0, M0, U1, CLEAN])
    print("   %-22s %s" % ("", " ".join("%.3f:%s" % (x / 1e6, v) for x, v in
                                        steps(c, 1_000_000, 10_000_000, 100_000))))
    lo, hi = flip(c, 1_000_000, 10_000_000)
    print("   bisected: kept at %d bp, dropped at %d bp  (2 x 63 x %d + 1 = %d)"
          % (lo, hi, c.s, 2 * 63 * c.s + 1))


# ---------------------------------------------------------------------------
# §4 — the IBD1 pass has no such bound
# ---------------------------------------------------------------------------

#: IBD1-live (no opposite homozygote), IBD2-dead (44 mismatches), `inf1` = 20.
RUN1 = {"hom1": 20, "ibs1": 44}
#: one opposite homozygote at bit 0 — the left flank an IBD1 run reaches back over
LFLANK = S.at({"ibs1": 64}, ibs0=0)
#: one opposite homozygote at bit 63 — the right flank
RFLANK = S.at({"ibs1": 64}, ibs0=63)


def ibd1w():
    """A one-word IBD1 run whose reported call is 189 markers, swept over the floor.

    The IBD1 pass has the bound too, over the run's own complete words, but its
    comparison is one unit of `seglength / 2` tighter: the IBD2 pass keeps a window of
    exactly `min_bp // 2` and this one does not.  Bisected at four spacings, the IBD2
    flip is always `(2w + 1, 2w + 2)` and the IBD1 flip `(2w - 1, 2w)`.
    """
    print("\u00a74 the same bound on the IBD1 pass, one unit tighter")
    for s2, nw2 in ((30_000, 50), (45_000, 50), (50_000, 45), (70_000, 25)):
        c1 = canvas([LFLANK, RUN1, RFLANK], s=s2, nw2=nw2, tag="_i1")
        c2 = canvas([M0, M0, CLEAN], s=s2, nw2=nw2)
        w = 63 * s2
        print("   spacing %5d  window %8d   IBD2 flip %s (want %s)   IBD1 flip %s"
              " (want %s)"
              % (s2, w, flip(c2, 1_000_000, 10_000_000), (2 * w + 1, 2 * w + 2),
                 flip(c1, 1_000_000, 10_000_000, what=1), (2 * w - 1, 2 * w)))
    c = canvas([LFLANK, RUN1, RUN1, RFLANK], tag="_i1")
    print("   two-word run (window 127 = 6.350 Mb) survives every floor under 10 Mb:")
    print("     %s" % [(x / 1e6, v) for x, v in
                       steps(c, 1_000_000, 10_000_000, 1_000_000, what=1)])


# ---------------------------------------------------------------------------
# §5 — what the IBD1 merge's budget is summed over
# ---------------------------------------------------------------------------

#: interruption word: 2 opposite homozygotes, 4 het-vs-A1A1, no A1A1/A1A1
INTW = word(z=2, m=58, v=4)


def refused(vr):
    """A word with no opposite homozygote and `vr` het-vs-A1A1 markers.

    Between two interruption words it forms a run of its own, which the gate refuses
    while `vr < 10` — `20-…` §6's stepped-over run.
    """
    return word(m=64 - vr, v=vr)


def budget():
    """Sweep the gate-refused run's own informative load with the budget on the boundary.

    Four-word runs either side; the interruption is `[INTW, refused(vr), INTW]`, so
    `bad = 4` and the two unusable words carry `V = 8` and `U = 0`.  `4 * (4 - 2) = 8`,
    and `X = V if V >= 10 else U`, so summed over the unusable words alone `X = 0` and no
    `vr` can merge.  Summed over every word between the runs, `V = 8 + vr` and the merge
    turns on at `vr = 2` and stays on until the gate stops refusing the run at `vr = 10`.
    """
    print("§5 does a stepped-over gate-refused run's own markers enter the budget?")
    print("   4 runs, [INTW, refused(vr), INTW], 4 runs; gap 3.860 Mb, L = 5")
    print("   bad = 4, unusable words carry U = 0 and V = 8; 4*(4-2) = 8")
    print("   %3s  %8s  %-12s   %-14s %-14s"
          % ("vr", "markers", "(words,calls)", "span=unusable", "span=all"))
    for vr in range(0, 12):
        c = M.cv("win1_bud_%d" % vr,
                 [RUN1] * 4 + [INTW, refused(vr), INTW] + [RUN1] * 4,
                 nw2=70, spacing=20_000)
        r = S.many([(c, ("--seglength", "5.000000"))])[0]
        m = S.mk(c, r, 1)
        _, p2 = c.positions()
        st = [M.wstat(x) for x in c.words]
        pred = {}
        for span in ("unusable", "all"):
            calls = M.ibd1(st, p2, 5_000_000, span=span)
            pred[span] = sum(p2[b] - p2[a] for a, b in calls) / c.s
        print("   %3d  %8.0f  %-12s   %-14.0f %-14.0f"
              % (vr, m, S.wc(m), pred["unusable"], pred["all"]))
    print("   one call = merged, two = split.  The flip is at vr = 2, which is where")
    print("   the refused run's own het-vs-A1A1 markers take V from 8 to 10 and hand")
    print("   `X` the switch.  At vr >= 10 the middle run passes the gate itself, so")
    print("   it is an endpoint and both readings merge for a different reason.")


# ---------------------------------------------------------------------------
# §6 — where the bound is asked
# ---------------------------------------------------------------------------

Z1 = word(z=1, h=63)        # unusable, unbridgeable, and cheap enough to merge across


def when():
    """A run that fails the bound, next to one that passes, with a mergeable gap.

    `[CLEAN, Z1, CLEAN x4]`: the first run's window is one word (63 = 3.150 Mb) and the
    second's is four (255 = 12.750 Mb).  At `--seglength 8` the bound refuses the first.
    If it is asked with the gate — before the merge — the refused run cannot be an
    endpoint and the answer is the second run alone, 255 markers.  If it is asked only
    when the call is emitted, the two merge first and the merged window passes: 383.
    """
    print("§6 the bound is asked at emit, not with the gate")
    c = canvas([S.CLEAN, Z1] + [S.CLEAN] * 4)
    print("   [CLEAN, Z1, CLEAN x4]  merged = 383 markers, second run alone = 255")
    for L in (1_000_000, 5_000_000, 6_000_000, 7_000_000, 8_000_000, 9_000_000):
        print("     L=%5.1f Mb -> %s markers" % (L / 1e6, read(c, L)))
    print("   (pre-merge predicts 255 from L = 6.3 Mb up; emit-only predicts 383)")


# ---------------------------------------------------------------------------
# §7 — held out
# ---------------------------------------------------------------------------

def rword3(rng):
    """One random chr2 word, drawn over the axes the window bound lives on.

    `mergelab.rword2` was drawn for the merge and puts its features at shuffled bits with
    a fixed budget; what this needs instead is a healthy supply of **mismatch-only** words
    with the mismatch bits placed anywhere, because `ml` and `mf` are what decide how far
    a call reaches past its own run — and therefore how far the reported length can run
    ahead of the gate window.
    """
    ks = ["zero"] * WORD
    t = rng.random()
    if t < 0.25:
        kind, cnt = "ibs0", rng.choice([1, 1, 2, 4, 10, 64])
    elif t < 0.55:
        kind, cnt = "ibs1", rng.choice([2, 2, 3, 5, 9])
    else:
        kind, cnt = "ibs1", rng.choice([0, 0, 1])
    for b in rng.sample(range(WORD), cnt):
        ks[b] = kind
    for b in rng.sample(range(WORD), rng.choice([0, 8, 16, 30, 50, 64])):
        if ks[b] == "zero":
            ks[b] = rng.choice(["hethet", "hethet", "hethet", "hom1"])
    return ks


def rblock4(rng, width):
    """A whole random canvas, drawn as a sentence of run/interruption *units*.

    `rword3` draws each word independently, and the configuration the bound decides — a
    short run whose call reaches two words back — then turns up on about one canvas in a
    hundred.  Drawing units instead ("a couple of mismatch-only words, then a one- or
    two-word run") puts it on most of them, while every word's content, every mismatch
    bit and the order of the units stay random.
    """
    ks = []
    while len(ks) < width:
        t = rng.random()
        if t < 0.20:
            ks.append(rword3(rng))
            continue
        if t < 0.45:
            ks.append(_mk(rng, "ibs0", rng.choice([1, 2, 8, 64])))
            continue
        if t < 0.75:
            for _ in range(rng.choice([1, 2, 2])):
                ks.append(_mk(rng, "ibs1", rng.choice([2, 3, 5, 9])))
        for _ in range(rng.choice([1, 1, 2, 3])):
            ks.append(_mk(rng, "ibs1", rng.choice([0, 0, 1])))
    return ks[:width]


def _mk(rng, kind, cnt):
    """One word: `cnt` markers of `kind`, the rest random filler.

    The bits are usually random, but two draws in five put them in one contiguous block
    at a random offset. That is a property of the *sampler*, not of the rule: how far a
    call reaches past its own run is set by the flanking word's first and last mismatch
    bit, and scattering `cnt` bits over 64 puts the last one near bit 63 almost always,
    so a purely uniform draw almost never builds the long reach the bound decides.
    """
    ks = ["zero"] * WORD
    if rng.random() < 0.4 and cnt:
        off = rng.randrange(WORD - cnt + 1)
        bits = range(off, off + cnt)
    else:
        bits = rng.sample(range(WORD), cnt)
    for b in bits:
        ks[b] = kind
    for b in rng.sample(range(WORD), rng.choice([16, 30, 50, 64])):
        if ks[b] == "zero":
            ks[b] = rng.choice(["hethet", "hethet", "hethet", "hom1"])
    return ks


def battery(seed=20260901, count=60, seglen_bp=10_000_000, width=24, nw2=30,
            spacing=50_000, rules=(), draw=None):
    """Random IBD2-native canvases graded against one or more `predict`-style models.

    **The spacing matters.** At `mergelab`'s 20 kb a word spans 1.26 Mb, an endpoint can
    reach at most one word either side, and a call long enough to clear a 10 Mb floor is
    eight words wide — so its window clears `L / 2` for free and the bound never bites.
    At 50 kb, the corpus's own spacing, one word is 3.15 Mb and a one-word run with both
    reaches spent measures 9.5 Mb: exactly the regime the bound decides.
    """
    import random
    rng = random.Random(seed)
    cvs = [M.cv("w1d%d_%d_%d_%d" % (seed, seglen_bp // 1000, spacing, t),
                (draw(rng, width) if draw else rblock4(rng, width)),
                nw2=nw2, spacing=spacing)
           for t in range(count)]
    res = S.many([(c, ("--seglength", "%.6f" % (seglen_bp / 1e6))) for c in cvs])
    score = [0] * len(rules)
    live = 0
    for c, r in zip(cvs, res):
        _, p2 = c.positions()
        st = [M.wstat(x) for x in c.words]
        got = S.mk(c, r, 2)
        wants = [sum(p2[b] - p2[a] for a, b in f(st, p2, seglen_bp)) / c.s
                 for f in rules]
        live += len(set(round(w, 1) for w in wants)) > 1
        for k, w in enumerate(wants):
            score[k] += abs(got - w) <= 0.3
    return score, len(cvs), live


def held_out():
    """`23` against `21` on three unused seeds, in the two regimes where the bound bites.

    A one-word window spans `63s` and can carry at most two more words of reach, so the
    bound decides a call only when `126s < seglength <= 189s`.  At the corpus's 50 kb
    that window straddles the 10 Mb floor, which is exactly why the corpus sees the
    clause at 10 Mb and at no smaller floor; here it is reproduced twice, at 30 kb
    against `--seglength 5` and at 60 kb against `--seglength 10`.
    """
    root = os.path.dirname(os.path.dirname(os.path.dirname(
        os.path.dirname(os.path.abspath(__file__)))))
    sys.path.insert(0, os.path.join(root, "tests", "parity", "fit"))
    import seg21 as S21  # noqa: E402
    import seg23 as S23  # noqa: E402
    print("§7 held out: random IBD2 canvases on unused seeds")
    tot = [0, 0]
    n = live = 0
    for seed in (20260901, 55512345, 987654321):
        for L, s, nw2, width in ((5_000_000, 30_000, 50, 42),
                                 (10_000_000, 60_000, 26, 18)):
            sc, t, lv = battery(seed, 60, L, width=width, nw2=nw2, spacing=s,
                                rules=(S21.predict, S23.predict))
            tot = [a + b for a, b in zip(tot, sc)]
            n += t
            live += lv
            print("   seed %-10d L=%2d Mb  %2d kb    21: %2d/%d    23: %2d/%d"
                  "   (bound live on %2d)"
                  % (seed, L // 10 ** 6, s // 1000, sc[0], t, sc[1], t, lv))
    print("   TOTAL                          21: %d/%d   23: %d/%d   (live on %d)"
          % (tot[0], n, tot[1], n, live))


def held_out1():
    """The IBD1 pass, held out: `20-…`'s battery in the regime where the bound bites."""
    print("§8 held out: random IBD1 canvases on unused seeds")
    tot = [0, 0, 0]
    n = 0
    for seed in (20260901, 55512345, 987654321):
        for L, s, nw2, width in ((5.0, 30_000, 50, 42), (10.0, 60_000, 26, 18)):
            got = []
            for span, win in (("unusable", None), ("all", None), ("all", 2)):
                ok, t, _d, _b = M.battery(seed, 60, L, width=width, nw2=nw2, spacing=s,
                                          span=span, window=win)
                got.append(ok)
            tot = [a + b for a, b in zip(tot, got)]
            n += t
            print("   seed %-10d L=%2d Mb  %2d kb    20: %2d/%d   +span: %2d/%d"
                  "   +window: %2d/%d" % (seed, int(L), s // 1000,
                                          got[0], t, got[1], t, got[2], t))
    print("   TOTAL                          20: %d/%d   +span: %d/%d   +window: %d/%d"
          % (tot[0], n, tot[1], n, tot[2], n))


SECTIONS = {"2": window, "3": ends, "4": ibd1w, "5": budget, "6": when, "7": held_out,
            "8": held_out1}

if __name__ == "__main__":
    want = sys.argv[1:] or list(SECTIONS)
    for k in want:
        SECTIONS[k]()
        print()
