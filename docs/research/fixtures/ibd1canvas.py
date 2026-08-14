#!/usr/bin/env python3
"""The IBD1-native canvas — fixtures behind `docs/research/18-ibd1-caller.md`.

`docs/research/17-seg-caller.md` built a `.seg`-native canvas whose read-back column is
`IBD2Seg`; every word meant to be *called* there had to be IBD2-clean, which meant no
opposite homozygote and at most one het-vs-hom mismatch.  This rig is the mirror image:
the read-back column is **`IBD1Seg`**, so the painted region must be **IBD2-free** while
IBS0 placement stays under control.  A word carrying thirty-four het-vs-hom mismatches is
exactly that — perfectly usable to the IBD1 pass, refused outright by the IBD2 one — so
IBD2 is silenced by paint rather than by walls, and every fixture below is checked to
report `IBD2Seg 0.0000` before its `IBD1Seg` is read.

The instrument is `segcanvas.Canvas` unchanged: chromosome 1 carries the same 5-word
`inf1` carrier (it is what earns the pair a `.seg` row), chromosome 2 is painted one
complete word at a time and walled with all-IBS0 words, and its uniform spacing puts `D`
just over the reference's 100 Mb floor so one ulp of the printed column is about a ninth
of a marker gap.  `segcanvas.mk(cv, res, 1)` subtracts the carrier's fixed contribution
and divides by the spacing, so the column reads back **the number of marker intervals
called on chromosome 2**.

Decoding is the same trick: IBD1 calls inside one usable segment are emitted adjacent
(each new call starts one marker past the previous call's end), so `c` calls covering `w`
whole words measure `64w - c` markers and `c = (-M) mod 64` recovers both counts.

    python3 ibd1canvas.py           # every section
    python3 ibd1canvas.py 2 5       # only sections 2 and 5
    IBD1CANVAS_JOBS=12 python3 ibd1canvas.py

Nothing here reads KING's source.
"""

from __future__ import annotations

import itertools
import os
import random
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import segcanvas as S  # noqa: E402

WORD = S.WORD
ROOT = os.path.dirname(os.path.abspath(__file__))

# A cache of this campaign's reference answers, kept apart from `segcanvas`'s own.
S.CACHE = os.path.join(ROOT, "ibd1canvas_measured.json")
S.JOBS = int(os.environ.get("IBD1CANVAS_JOBS", os.environ.get("SEGCANVAS_JOBS", "8")))

# ---------------------------------------------------------------------------
# the alphabet
# ---------------------------------------------------------------------------
#
# Every letter that is meant to be *inside* an IBD1 call carries at least two het-vs-hom
# mismatches, which is what keeps the IBD2 pass off the canvas (`17-…` §3).  `B` and `b`
# are the deliberate exceptions: they are the only words the IBD2 pass can use, and they
# exist for §6, where the two passes have to meet.

WALL = S.WALL                                        # 64 opposite homozygotes
K = {"ibs1": 34, "hom1": 12, "zero": 18}             # IBD1-usable, IBD2-dead, inf1 = 12
K0 = {"ibs1": 34, "zero": 30}                        # ...and worth no `inf1` at all
B = {"hom1": 12, "hethet": 52}                       # usable to BOTH passes, inf1 = 12
B0 = {"hethet": 64}                                  # usable to IBD2 only (inf1 = 0)


def kw(**counts):
    """A word composition from explicit marker counts, padded out with A2A2/A2A2."""
    n = sum(counts.values())
    assert n <= WORD, counts
    d = dict(counts)
    d["zero"] = d.get("zero", 0) + WORD - n
    return d


def Z(bits, base=None):
    """`base` (default `K`) with opposite homozygotes forced at the named bit positions."""
    return S.at(base or K, ibs0=bits)


def Y(bits, base=None):
    """`base` with het-vs-hom mismatches forced at the named bits — an IBD2 boundary."""
    return S.at(base or B, ibs1=bits)


ALPHA = {
    "K": K,          # IBD1 only
    "k": K0,         # IBD1 only, no `inf1`
    "W": WALL,       # breaks both passes
    "B": B,          # both passes
    "b": B0,         # IBD2 only
}


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------

GATE1 = 10        # `inf1` markers a run must carry over its own complete words (§4)
MIN_RUN1 = 1      # complete words a run must span (§4)


def wordinfo1(spec):
    """`(has_ibs0, first_ibs0_bit, last_ibs0_bit, inf1)` for one word."""
    ks = S.expand(spec)
    z = [i for i, k in enumerate(ks) if k in ("ibs0", "ibs0b")]
    inf1 = sum(1 for k in ks if k in ("hom1", "ibs1b"))
    return (bool(z), z[0] if z else None, z[-1] if z else None, inf1)


def predict1(info, w0=0, w1=None, lo=None, hi=None, gate=GATE1, min_run=MIN_RUN1,
             left_rule="last+1", right_rule="last", push="never", bridge=False,
             clip=1):
    """`.seg` IBD1 calls as marker intervals — the rule of `18-ibd1-caller.md` §7.

    `info[k]` is `wordinfo1` for global word `k`; the usable segment is words `w0..w1`
    covering markers `lo..hi`.  Every keyword is a knob a section below bisects.
    """
    n = len(info)
    w1 = n - 1 if w1 is None else w1
    lo = WORD * w0 if lo is None else lo
    hi = WORD * (w1 + 1) - 1 if hi is None else hi
    z = [info[k][0] for k in range(n)]
    fb = [info[k][1] for k in range(n)]
    lb = [info[k][2] for k in range(n)]
    i1 = [info[k][3] for k in range(n)]

    ok = [not z[k] for k in range(n)]
    if bridge:                                    # §5 says this never happens
        for k in range(w0 + 1, w1):
            if not ok[k] and ok[k - 1] and ok[k + 1]:
                ok[k] = True

    runs, k = [], w0
    while k <= w1:
        if not ok[k]:
            k += 1
            continue
        a = k
        while k <= w1 and ok[k]:
            k += 1
        runs.append((a, k - 1))

    out, emitted = [], 0
    for a, b in runs:
        if b - a + 1 < min_run:
            continue
        if sum(i1[t] for t in range(a, b + 1)) < gate:
            continue
        # left: one marker past the *last* opposite homozygote of the word before the run
        if a - 1 >= w0:
            bit = fb[a - 1] if left_rule.startswith("first") else lb[a - 1]
            left = WORD * a if bit is None else WORD * (a - 1) + bit + 1
        else:
            left = lo
        # right: out to the *last* opposite homozygote of the word after the run
        if b + 1 <= w1:
            bit = fb[b + 1] if right_rule.startswith("first") else lb[b + 1]
            right = WORD * (b + 1) + WORD - 1 if bit is None else WORD * (b + 1) + bit
        else:
            right = hi
        if push == "always" and emitted:
            left = max(left, WORD * (a + 1))
        emitted += 1
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1] + clip)
        if left <= right:
            out.append((left, right))
    return out


def pieces(c, c2, cut="exclusive"):
    """One IBD1 call with the IBD2 calls cut out of it, as marker intervals (§6).

    `"exclusive"` — the measured convention — removes the IBD2 call's own end markers
    too, so `[lo, hi]` cut by `[a, b]` leaves `[lo, a-1]` and `[b+1, hi]`.
    `"inclusive"` is the naive "length minus overlap" the engine used to compute.
    """
    d = 1 if cut == "exclusive" else 0
    lo, hi = c
    out, cur = [], lo
    for a, b in sorted(c2):
        if b < lo or a > hi:
            continue
        if a - d > cur:
            out.append((cur, a - d))
        cur = max(cur, b + d)
    if cur + (0 if d else 1) <= hi:
        out.append((cur, hi))
    return out


def subtract(c1, c2, pos, seglen_bp, frag="drop", cut="exclusive"):
    """`IBD1Seg`'s numerator: the IBD1 calls with the IBD2 calls taken out (§6).

    `frag` says what happens to the pieces: `"drop"` applies the `--seglength` floor to
    each piece (the measured rule), `"keep"` counts them all, `"whole"` applies the floor
    to their total.
    """
    tot = 0
    for c in c1:
        lens = [int(pos[b] - pos[a]) for a, b in pieces(c, c2, cut)]
        if frag == "drop":
            tot += sum(v for v in lens if v >= seglen_bp)
        elif frag == "whole":
            tot += sum(lens) if sum(lens) >= seglen_bp else 0
        else:
            tot += sum(lens)
    return tot


def model(cv, seglen_bp=1_000_000, frag="drop", cut="exclusive", **kw1):
    """`(IBD1 markers, IBD2 markers)` the model predicts for a canvas's chromosome 2.

    The `--seglength` floor is applied to the IBD2 calls **before** they are subtracted:
    a dropped IBD2 call is not taken out of `IBD1Seg` at all (§6.3).
    """
    _, p2 = cv.positions()
    i1 = [wordinfo1(x) for x in cv.words]
    i2 = [S.wordinfo(x) for x in cv.words]
    c2 = [c for c in S.predict(i2) if p2[c[1]] - p2[c[0]] >= seglen_bp]
    c1 = [c for c in predict1(i1, **kw1) if p2[c[1]] - p2[c[0]] >= seglen_bp]
    return (subtract(c1, c2, p2, seglen_bp, frag, cut) / cv.s,
            sum(p2[b] - p2[a] for a, b in c2) / cv.s)


# ---------------------------------------------------------------------------
# running a family of canvases
# ---------------------------------------------------------------------------

def run(items, extra=("--seglength", "1"), what=1):
    """`[(label, canvas)]` -> `[(label, canvas, result)]`, cached and parallel."""
    res = S.many([(cv, extra) for _, cv in items])
    return [(lab, cv, r) for (lab, cv), r in zip(items, res)]


def line(lab, cv, r, expect=None, width=26):
    """One printed row: the read-back, the decode, and the IBD2 isolation check."""
    m = S.mk(cv, r, 1)
    d = S.wc(m, tol=0.3)
    s = "%8.2f mk %-22s ibd2 %s" % (
        m, "" if d is None else "= %d words / %d calls" % d, r["row"]["IBD2Seg"])
    if expect is not None:
        s += "   model %8.2f %s" % (expect, "ok" if abs(m - expect) <= 0.3 else "MISS")
    return "  %-*s %s" % (width, lab, s)


def seq_canvas(s, **kwargs):
    kwargs.setdefault("nw2", max(16, len(s) + 8))
    return S.Canvas("q1_" + s, [ALPHA[c] for c in s], **kwargs)


# ---------------------------------------------------------------------------
# sections
# ---------------------------------------------------------------------------

def hdr(t):
    print("\n== %s" % t)


def section0():
    hdr("0. the rig — chromosome 2 reads back IBD1, and IBD2 stays at zero")
    items = [("all wall", S.Canvas("i1_empty", []))]
    for W in (1, 2, 3, 4, 8):
        items.append(("K x %d" % W, S.Canvas("i1_kw%d" % W, [K] * W)))
    print("    a block of W IBD1-callable words, walled by all-IBS0 words:")
    for lab, cv, r in run(items):
        exp = 0.0 if lab == "all wall" else WORD * (int(lab.split()[-1]) + 1) - 1
        print(line(lab, cv, r, exp))
    print("    -> 64(W+1) - 1: the call swallows the whole trailing wall word, because")
    print("       the wall's last IBS0 sits at bit 63 (§2).")


def section1():
    hdr("1. the word predicate — one opposite homozygote, no tolerance, no bridging")
    zs = (0, 1, 2, 3, 5, 64)
    items = []
    for j in (1, 2, 3):
        for z in zs:
            body = [K] * 8
            for t in range(3, 3 + j):
                body[t] = K if z == 0 else (WALL if z == 64 else Z(list(range(z))))
            items.append((("z", j, z), S.Canvas("i1z_%d_%d" % (j, z), body)))
    other = [("ibs1=64", kw(ibs1=64)), ("ibs1b=64", kw(ibs1b=64)),
             ("miss=64", kw(miss=64)), ("missA=64", kw(missA=64)),
             ("hethet=62 + 2 mis", kw(hethet=62, ibs1=2)),
             ("hom1=62 + 2 mis", kw(hom1=62, ibs1=2)),
             ("zero=62 + 2 mis", kw(ibs1=2)),
             ("ibs1=63 + ibs0=1", S.at(kw(ibs1=63), ibs0=[63]))]
    for lab, spec in other:
        for j in (1, 2):
            body = [K] * 8
            for t in range(3, 3 + j):
                body[t] = spec
            items.append((("o", j, lab), S.Canvas("i1o_%d_%s" % (j, S.slug(lab)), body)))
    got = {k: (S.wc(S.mk(cv, r, 1), tol=0.3), r["row"]["IBD2Seg"])
           for k, cv, r in run(items)}
    print("    8 IBD1-callable words with `j` consecutive words replaced by one carrying")
    print("    `z` opposite homozygotes, read back as words called / calls:")
    print("     j\\z  " + "".join("%-9d" % z for z in zs))
    for j in (1, 2, 3):
        print("     %-4d " % j + "".join("%-9s" % ("%d/%d" % got[("z", j, z)][0])
                                         for z in zs))
    print("    -> one IBS0 splits the block at every j; a *lone* bad word is never")
    print("       absorbed (the j = 1 row is 2 calls, not 1).")
    print()
    print("    the same block with words 3.. replaced by other content:")
    for lab, _ in other:
        print("    %-20s j=1 %-9s j=2 %-9s   ibd2 %s"
              % (lab, "%d/%d" % got[("o", 1, lab)][0], "%d/%d" % got[("o", 2, lab)][0],
                 got[("o", 1, lab)][1]))
    print("    -> nothing but an opposite homozygote breaks an IBD1 word.")


BITS = [("{0}", [0]), ("{1}", [1]), ("{5}", [5]), ("{31}", [31]), ("{32}", [32]),
        ("{62}", [62]), ("{63}", [63]), ("{0,1}", [0, 1]), ("{62,63}", [62, 63]),
        ("{0,63}", [0, 63]), ("{0..19}", list(range(20))),
        ("{44..63}", list(range(44, 64))), ("all 64", list(range(64)))]


def section2():
    hdr("2. the endpoints — how far a call reaches into the words that bound it")
    print("    block = 6 IBD1-callable words with one boundary word beside them, whose")
    print("    opposite homozygotes sit at the named bits.  6 words walled = 447 mk;")
    print("    the run's own word-aligned span is 383 mk.")
    left = [(lab, S.Canvas("i2L_" + S.slug(lab), [Z(b)] + [K] * 6, nw2=20))
            for lab, b in BITS]
    right = [(lab, S.Canvas("i2R_" + S.slug(lab), [K] * 6 + [Z(b)], nw2=20))
             for lab, b in BITS]
    rl, rr = run(left), run(right)
    print("    %-10s %-26s %s" % ("IBS0 bits", "left of the block", "right of it"))
    for (lab, b), (_, ca, ra), (_, cb, rb) in zip(BITS, rl, rr):
        ma, mb = S.mk(ca, ra, 1), S.mk(cb, rb, 1)
        print("    %-10s %8.1f mk (+%3.0f)      %8.1f mk (+%3.0f)"
              % (lab, ma, ma - 447, mb, mb - 383))
    print("    left  = 64(a-1) + (last IBS0 bit) + 1     -> extension 63 - lastbit")
    print("    right = 64(b+1) + (last IBS0 bit)         -> extension  1 + lastbit")
    print("    Both ends read the *last* opposite homozygote of the immediately")
    print("    flanking word.  There is no 63-marker reach and no whole-word blocking:")
    print("    this is not the IBD2 geometry of `17-…` §5.")
    scope = [("R: K6 Z{5}", [K] * 6 + [Z([5])]),
             ("R: K6 Z{5} Z{10}", [K] * 6 + [Z([5]), Z([10])]),
             ("R: K6 Z{5} WALL", [K] * 6 + [Z([5]), WALL]),
             ("R: K6 WALL Z{5}", [K] * 6 + [WALL, Z([5])]),
             ("L: Z{10} Z{5} K6", [Z([10]), Z([5])] + [K] * 6),
             ("L: WALL Z{5} K6", [WALL, Z([5])] + [K] * 6),
             ("L: Z{5} WALL K6", [Z([5]), WALL] + [K] * 6)]
    print()
    print("    ...and it reads one word and no further:")
    for lab, cv, r in run([(lab, S.Canvas("i2s_" + S.slug(lab), body, nw2=20))
                           for lab, body in scope]):
        print(line(lab, cv, r, width=18))


def section3():
    hdr("3. the push — there isn't one")
    items = [("K2 W K2", [K, K, WALL, K, K]),
             ("K2 W K2 W K2", [K, K, WALL, K, K, WALL, K, K]),
             ("K2 W2 K2 W2 K2", [K, K, WALL, WALL, K, K, WALL, WALL, K, K]),
             ("K W K4", [K, WALL] + [K] * 4),
             ("k W K4", [K0, WALL] + [K] * 4),
             ("K4 W K", [K] * 4 + [WALL, K])]
    for lab, cv, r in run([(lab, S.Canvas("i3_" + S.slug(lab), body, nw2=20))
                           for lab, body in items]):
        info = [wordinfo1(x) for x in cv.words]
        pn = S.predict_bp(cv, predict1(info, push="never")) / cv.s
        pa = S.predict_bp(cv, predict1(info, push="always")) / cv.s
        print("    %-16s ref %7.1f mk    push=never %6.0f   push=always %6.0f"
              % (lab, S.mk(cv, r, 1), pn, pa))
    print("    -> every call starts exactly where the previous IBS0 left off. The")
    print("       one-word push that `17-…` §6 measured on the IBD2 pass does not")
    print("       exist on this one, at any number of calls.")


def section4():
    hdr("4. the gate — `inf1` >= 10 over the run's own complete words")
    items = []
    for k in (0, 8, 9, 10, 11, 20):
        items.append(("1 word, %d A1A1/A1A1" % k, [kw(hom1=k, ibs1=34)]))
    for k in (9, 10, 11):
        items.append(("1 word, %d het-vs-A1A1" % k, [kw(ibs1b=k, ibs1=34)]))
    for k in (4, 5):
        items.append(("2 words, %d each" % k, [kw(hom1=k, ibs1=34)] * 2))
    for a, b in ((4, 5), (5, 4), (5, 5)):
        items.append(("1 word, %d A1A1 + %d het-vs-A1A1" % (a, b),
                      [kw(hom1=a, ibs1b=b, ibs1=34)]))
    items.append(("1 word, 62 HetHet", [kw(hethet=62, ibs1=2)]))
    rich_r = S.at(kw(hom1=40, ibs1=20), ibs0=[63])
    rich_l = S.at(kw(hom1=40, ibs1=20), ibs0=[0])
    items.append(("9 + rich right flank", [kw(hom1=9, ibs1=34), rich_r]))
    items.append(("9 + rich left flank", [rich_l, kw(hom1=9, ibs1=34)]))
    items.append(("10 + rich right flank", [kw(hom1=10, ibs1=34), rich_r]))
    for lab, cv, r in run([(lab, S.Canvas("i4_" + S.slug(lab), body, nw2=20))
                           for lab, body in items]):
        print(line(lab, cv, r, width=24))
    print("    -> 9 refused / 10 accepted on A1A1/A1A1 and on het-vs-A1A1 alike, split")
    print("       across two words, and mixed; HetHet is worth nothing. The markers the")
    print("       call reaches into do NOT count: 9 + a flank carrying 40 is refused.")


def section5():
    hdr("5. bridging — a lone bad word is never absorbed on the IBD1 pass")
    items = []
    for lab, mid in (("1 IBS0 at bit 0", Z([0])), ("1 IBS0 at bit 63", Z([63])),
                     ("1 IBS0 at bit 31", Z([31])), ("all 64 IBS0", WALL)):
        items.append((lab, [K] * 3 + [mid] + [K] * 3))
    items.append(("no bad word", [K] * 7))
    for lab, cv, r in run([(lab, S.Canvas("i5_" + S.slug(lab), body, nw2=20))
                           for lab, body in items]):
        info = [wordinfo1(x) for x in cv.words]
        pn = S.predict_bp(cv, predict1(info)) / cv.s
        pb = S.predict_bp(cv, predict1(info, bridge=True)) / cv.s
        print("    %-18s ref %7.1f mk   no bridge %6.0f   bridge %6.0f"
              % (lab, S.mk(cv, r, 1), pn, pb))


FAM6 = [("K B3 K", [K, B, B, B, K]),
        ("K B6 K", [K] + [B] * 6 + [K]),
        ("K4 B4 K", [K] * 4 + [B] * 4 + [K]),
        ("K B4 K B4 K", [K] + [B] * 4 + [K] + [B] * 4 + [K]),
        ("k b4 k", [K0] + [B0] * 4 + [K0]),
        ("K b4 K", [K] + [B0] * 4 + [K])]


def section6():
    hdr("6. the IBD1/IBD2 overlap — what `IBD1Seg` subtracts, and in what order")
    print("    `B` is the only word both passes can use (no IBS0, no mismatch, inf1=12);")
    print("    `K` bounds it, so the IBD2 call spills one word out of the B block while")
    print("    the IBD1 call runs to the walls.  `b` is `B` without `inf1`, so the IBD2")
    print("    pass takes it and the IBD1 gate does not.")
    items = [(lab, S.Canvas("i6_" + S.slug(lab), body, nw2=20)) for lab, body in FAM6]
    print("\n    6.1 the cut is at marker granularity and EXCLUDES the IBD2 end markers")
    for lab, cv, r in run(items):
        e = model(cv, 1_000_000, cut="exclusive")[0]
        i = model(cv, 1_000_000, cut="inclusive")[0]
        print("      %-12s ref ibd1 %7.1f mk   exclusive %6.1f   inclusive %6.1f"
              % (lab, S.mk(cv, r, 1), e, i))
    print("      (`k b4 k` reads 0: the IBD1 gate refuses the run, and the IBD2 call over")
    print("       it is not subtracted from anywhere — the subtraction is per IBD1 call.)")

    print("\n    6.2 every piece faces the `--seglength` floor on its own")
    cv = items[0][1]
    frag_bp = 63 * cv.s
    ls = [(frag_bp - 1000) / 1e6, (frag_bp - 1) / 1e6, frag_bp / 1e6,
          (frag_bp + 1) / 1e6, (frag_bp + 1000) / 1e6]
    for L, r in zip(ls, S.many([(cv, ("--seglength", "%.6f" % L)) for L in ls])):
        kept = S.mk(cv, r, 1) > 1
        print("      K B3 K, piece = %d bp, --seglength %.6f Mb -> piece %s"
              % (frag_bp, L, "counted" if kept else "dropped"))
    print("      -> the floor is `>=`, bisected to the base pair, and it is applied to")
    print("         the piece, not to the 26.8 Mb IBD1 call it came from.")

    print("\n    6.3 a dropped IBD2 call is not subtracted (floor first, then subtract)")
    cv2 = S.Canvas("i6b_K_B_K", [K, B, K], nw2=40)
    ls2 = [1, 2.2, 2.3, 6.6, 6.7, 8, 9, 10]
    res = S.many([(cv2, ("--seglength", "%.6f" % L)) for L in ls2])
    print("      [K B K] at %d bp spacing: IBD1 call 8.921 Mb, IBD2 call 6.686 Mb,"
          % cv2.s)
    print("      the piece left over 2.205 Mb.")
    for L, r in zip(ls2, res):
        m1, m2 = S.mk(cv2, r, 1), S.mk(cv2, r, 2)
        d1, d2 = model(cv2, int(L * 1e6))
        print("      L=%-5s ref ibd1 %7.2f mk  ibd2 %7.2f mk   model %7.1f / %7.1f %s"
              % (L, m1, m2, d1, d2,
                 "ok" if abs(m1 - d1) <= 0.3 and abs(m2 - d2) <= 0.3 else "MISS"))
    print("      At L = 6.7..8 the IBD2 call is under the floor, vanishes from `IBD2Seg`,")
    print("      and the whole IBD1 call comes back — so the floor runs before the")
    print("      subtraction. At L = 9..10 the piece (now the whole call) is under the")
    print("      floor too and `IBD1Seg` goes to zero.")

    print("\n    6.4 the whole model against the reference, over the family and floors")
    miss = 0
    for L in (1, 2, 3, 4, 5, 6, 8, 10):
        for lab, cv, r in run(items, extra=("--seglength", "%d" % L)):
            m1, m2 = S.mk(cv, r, 1), S.mk(cv, r, 2)
            d1, d2 = model(cv, L * 1_000_000)
            if abs(m1 - d1) > 0.3 or abs(m2 - d2) > 0.3:
                miss += 1
                print("      MISS L=%d %-12s ref %7.1f/%7.1f model %7.1f/%7.1f"
                      % (L, lab, m1, m2, d1, d2))
    print("      %d of %d canvas x floor combinations reproduced"
          % (8 * len(items) - miss, 8 * len(items)))

    print("\n    6.5 the reference's own range check on `--seglength` (a cross-check)")
    for L in ("0.990000", "0.991000", "10.001000", "10.010000"):
        r = S.probe(cv, ("--seglength", L))
        print("      --seglength %-10s IBD1Seg %s   (%s)"
              % (L, r["row"]["IBD1Seg"],
                 "in range" if 0.99 < float(L) < 10.01 else "reverts to the 3 Mb default"))


def section7():
    hdr("7. the exhaustive word-sequence battery")
    letters = "KkW"
    seqs = ["".join(t) for n in (1, 2, 3, 4) for t in itertools.product(letters, repeat=n)]
    res = run([(s, seq_canvas(s)) for s in seqs])
    bad = []
    for s, cv, r in res:
        got = S.mk(cv, r, 1)
        want = S.predict_bp(cv, predict1([wordinfo1(x) for x in cv.words])) / cv.s
        if abs(got - want) > 0.3:
            bad.append((s, round(got, 1), round(want, 1)))
    print("    K = IBD1-callable (inf1 12), k = the same with inf1 0, W = all-IBS0")
    print("    predict1() reproduces %d of %d" % (len(seqs) - len(bad), len(seqs)))
    for b in bad:
        print("      miss", b)


def random_word(rng, mixed=False):
    """One chromosome-2 word, drawn to exercise every axis the IBD1 rule reads.

    `mixed` lets a word carry fewer than two het-vs-hom mismatches, which is what makes
    it usable to the IBD2 pass as well — the battery that grades both columns.
    """
    ks = ["zero"] * WORD
    slots = list(range(WORD))
    rng.shuffle(slots)
    n = 0
    if rng.random() < 0.35:
        for _ in range(rng.choice([1, 1, 2, 3, 8, 64])):
            if n < WORD:
                ks[slots[n]] = "ibs0"
                n += 1
    for _ in range(rng.choice([0, 2, 4, 6, 9, 10, 11, 14, 30])):
        if n < WORD:
            ks[slots[n]] = rng.choice(["hom1", "ibs1b"])
            n += 1
    mis = rng.choice([0, 1, 1, 2, 3, 8] if mixed else [2, 3, 8, 20, 34])
    for _ in range(mis):
        if n < WORD:
            ks[slots[n]] = "ibs1"
            n += 1
    for _ in range(rng.choice([0, 4, 12, 30, 64] if mixed else [0, 4, 12, 30])):
        if n < WORD:
            ks[slots[n]] = rng.choice(["hethet", "miss", "zero"])
            n += 1
    return ks


def random_canvases(seed, count, width=10, nw2=16, mixed=False):
    rng = random.Random(seed)
    return [S.Canvas("i1%s%d_%d" % ("mix" if mixed else "rnd", seed, t),
                     [random_word(rng, mixed) for _ in range(width)], nw2=nw2)
            for t in range(count)]


def battery(seed, count, width=10, mixed=False, seglen=1, verbose=False):
    """`predict1()` (and, when `mixed`, the whole `model()`) against the reference."""
    cvs = random_canvases(seed, count, width, mixed=mixed)
    res = S.many([(c, ("--seglength", "%d" % seglen)) for c in cvs])
    agree, bad, dirty = 0, [], 0
    for cv, r in zip(cvs, res):
        if float(r["row"]["IBD2Seg"]) != 0.0:
            dirty += 1
        g1, g2 = S.mk(cv, r, 1), S.mk(cv, r, 2)
        w1, w2 = model(cv, seglen * 1_000_000)
        if abs(g1 - w1) <= 0.3 and (not mixed or abs(g2 - w2) <= 0.3):
            agree += 1
        else:
            bad.append((cv.name, (round(g1, 1), round(g2, 1)), (round(w1, 1),
                                                                round(w2, 1))))
    if verbose:
        for b in bad:
            print("      miss", b)
    return agree, len(cvs), dirty, bad


def section8():
    hdr("8. random canvases — the model against the reference, out of sample")
    print("    IBD2-free canvases (only `IBD1Seg` is graded):")
    tot = ok = 0
    for seed, n in ((201, 80), (5150, 80), (99991, 80)):
        a, t, dirty, _ = battery(seed, n)
        tot += t
        ok += a
        print("      seed %-6d %3d / %3d agree (%.0f %%)   ibd2-contaminated %d"
              % (seed, a, t, 100 * a / t, dirty))
    print("      overall %d / %d (%.0f %%)" % (ok, tot, 100 * ok / tot))
    print("    mixed canvases — both columns graded, so the §6 subtraction is on trial:")
    tot2 = ok2 = 0
    for seed, n, L in ((3300, 60, 1), (61803, 60, 1), (3300, 60, 4), (61803, 60, 8)):
        a, t, _, bad = battery(seed, n, mixed=True, seglen=L)
        tot2 += t
        ok2 += a
        print("      seed %-6d --seglength %d   %3d / %3d agree (%.0f %%)"
              % (seed, L, a, t, 100 * a / t))
    print("      overall %d / %d (%.0f %%)" % (ok2, tot2, 100 * ok2 / tot2))


def section9():
    hdr("9. the `--seglength`-dependent run merge — measured, NOT modelled")
    print("    §5 says a lone bad word is never absorbed. That is true at the default")
    print("    3 Mb and it is not the whole story: raise `--seglength` past the gap the")
    print("    bad word leaves and the two runs merge into one.")
    print()
    print("    9.1 the gap threshold, bisected to the base pair at two spacings")
    for nw2 in (16, 30):
        cv = S.Canvas("i10_K3_Z0_K4_%d" % nw2, [K] * 3 + [Z([0])] + [K] * 4, nw2=nw2)
        # run 1 is words 0..2 (last marker 191), run 2 is words 4..7 (first marker 256)
        gap = 65 * cv.s
        ls = [(gap - 1) / 1e6, gap / 1e6, (gap + 1) / 1e6]
        for L, r in zip(ls, S.many([(cv, ("--seglength", "%.6f" % L)) for L in ls])):
            v = S.mk(cv, r, 1)
            print("      s=%d  gap = %d bp   --seglength %10.6f Mb -> %s"
                  % (cv.s, gap, L, "MERGED" if v > 574.5 else "split"))
    print("      -> merged iff `pos[first marker of the later run]`")
    print("         `- pos[last marker of the earlier run] < --seglength`, strictly.")
    print()
    print("    9.2 how much IBS0 the interruption may carry, at --seglength 8")
    tr = [(z, S.Canvas("i11_z%d" % z, [K] * 3 + [Z(list(range(z)))] + [K] * 4, nw2=16))
          for z in (4, 5, 6, 7, 8, 16, 64)]
    for (z, cv), r in zip(tr, S.many([(cv, ("--seglength", "8")) for _, cv in tr])):
        v = S.mk(cv, r, 1)
        print("      %2d opposite homozygotes -> %s" % (z, "MERGED" if v > 574.5 else "split"))
    print("      -> 5 merges, 6 does not, wherever the bits sit; and it is the total over")
    print("         the interrupting words, which may be more than one (§9.3).")
    print()
    print("    9.3 more than one interrupting word, and the negative")
    alt = [("two words, 1 IBS0 each, L=8", [K] * 3 + [Z([0]), Z([0])] + [K] * 3, 30, "8"),
           ("two words, 1 IBS0 each, L=6", [K] * 3 + [Z([0]), Z([0])] + [K] * 3, 30, "6"),
           ("three words, 1 IBS0 each, L=8",
            [K] * 3 + [Z([0])] * 3 + [K] * 2, 30, "8")]
    its = [(lab, S.Canvas("i12_" + S.slug(lab.split(",")[0] + lab[-2:]), b, nw2=nw), L)
           for lab, b, nw, L in alt]
    for (lab, cv, L), r in zip(its, S.many([(cv, ("--seglength", L))
                                            for _, cv, L in its])):
        m = sum(b - a for a, b in predict1([wordinfo1(x) for x in cv.words]))
        print("      %-32s ref %6.1f   unmerged %4d" % (lab, S.mk(cv, r, 1), m))
    print("      The clause is real and bisected, and it is deliberately NOT in the")
    print("      engine: the obvious generalisation (merge two runs when the words")
    print("      between them carry <= 5 IBS0 and the gap is under `--seglength`) makes")
    print("      the corpus much worse at 5 and 10 Mb — `IBD1Seg` 795 of 982 against 909")
    print("      at 5 Mb — because unrelated pairs are full of one-IBS0 interruptions.")
    print("      At the default 3 Mb it never fires (the gap is 65 marker intervals, and")
    print("      at real spacings that is over 3 Mb), which is exactly why `IBD1Seg` is")
    print("      exact on all 982 corpus rows there and is not at 5 and 10 Mb.")


SECTIONS = {"0": section0, "1": section1, "2": section2, "3": section3, "4": section4,
            "5": section5, "6": section6, "7": section7, "8": section8, "9": section9}


def main(argv):
    for k in (argv[1:] or sorted(SECTIONS)):
        SECTIONS[k]()
    S.save_cache()


if __name__ == "__main__":
    main(sys.argv)
