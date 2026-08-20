#!/usr/bin/env python3
"""The **fringe** canvas — a usable segment that does not start on a word boundary.

`segcanvas.py` and `ibd1canvas.py` both lay chromosome 2 out so that its first marker is
also the first marker of a 64-marker word of the global grid.  Every usable segment they
build is therefore word-aligned at both ends, and the clause that says what happens in
the **partial word** beyond a segment's word grid is invisible to them.  Both campaigns
say so in as many words (`17-…` §5: *"the canvas cannot see this; the corpus can"*), and
both then took that clause from the corpus — which is a fit, not a measurement.

This rig makes it visible.  The trick is to spend markers, not words:

    chr1   64*nw1 - f markers   the carrier, so the pair still earns a `.seg` row
    chr2   f + 64*nw2 + t       f head-fringe markers, nw2 complete words, t tail

Because chr1 is **short by exactly f markers**, chromosome 2 opens `f` markers before a
word boundary, so its usable segment is

    lo = 64*nw1 - f     first_word = nw1        head fringe = chr2 markers [0, f)
    hi = 64*(nw1+nw2) + t - 1                   tail fringe = the last t markers
                        last_word  = nw1+nw2-1

and the `nw2` painted words stay exactly word-aligned to the global grid, so every rule
`17-…` and `18-…` measured still applies to them unchanged.  The two fringes are painted
marker by marker and are the only thing that varies.

The head-fringe word is *shared*: its low bits are chromosome 1's last markers and its
high bits are chromosome 2's first.  Each segment can only own its own side, and the
carrier's mismatches all sit at bits below the fringe, so they can never pull a stop past
`seg.lo` — masked or not, the reading is the same.  That is a limitation of the rig, not
a result.

The ruler is `segcanvas`'s: chr2's spacing is uniform and chosen so one ulp of the printed
`%.4lf` is about an eighth of a marker gap, so the column reads back **the number of chr2
marker intervals called**, exactly.  `mk2` does that for `IBD2Seg`; `IBD1Seg` also carries
the carrier chromosome, so it is read *differentially* against a canvas with the same
chr1 (`mk1`, and §5's baseline) rather than by subtracting an assumed carrier length.

    python3 fringecanvas.py            # every section
    python3 fringecanvas.py 0 1        # only sections 0 and 1
    FRINGECANVAS_JOBS=12 python3 fringecanvas.py

Nothing here reads KING's source.
"""

from __future__ import annotations

import hashlib
import json
import os
import shutil
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixlab as F  # noqa: E402
import segcanvas as S  # noqa: E402

WORD = S.WORD
ROOT = os.path.dirname(os.path.abspath(__file__))

# This campaign's own cache, kept apart from the other two rigs'.
S.CACHE = os.path.join(ROOT, "fringecanvas_measured.json")
S.WORK = os.path.join(ROOT, "work", "fringecanvas")
S.JOBS = int(os.environ.get("FRINGECANVAS_JOBS", os.environ.get("SEGCANVAS_JOBS", "8")))

PAIR = S.PAIR
WALL = S.WALL          # 64 opposite homozygotes
CLEAN = S.CLEAN        # 64 HetHet: callable by both passes
CARRIER = S.CARRIER    # IBD1-clean, IBD2-dead

#: The IBD1 canvas letter (`18-…`): 34 het-vs-hom mismatches makes a word useless to the
#: IBD2 pass (`17-…` §3) and perfect for the IBD1 one, so an IBD1 fringe can be read with
#: IBD2 held at zero.  12 `hom1` keeps `inf1` over the gate.
K = {"ibs1": 34, "hom1": 12, "zero": 18}

NW2 = 16               # complete words on chr2, so the 100 Mb floor is cleared cheaply


class FringeCanvas(S.Canvas):
    """A canvas whose chr2 segment opens `f` markers and closes `t` markers off-grid."""

    def __init__(self, name, words, head=(), tail=(), nw1=6, sp1=33_000,
                 spacing=None, nsample=6, seed=5, carrier=None):
        self.name = name
        self.words = list(words)
        self.head = list(head)
        self.tail = list(tail)
        self.nw1 = nw1
        self.sp1 = sp1
        self.nsample = nsample
        self.seed = seed
        self.carrier = carrier or CARRIER
        self.gaps = None
        self.pad = (0, 0)
        f, t = len(self.head), len(self.tail)
        assert 0 <= f < WORD and 0 <= t < WORD, (f, t)
        # chr1 loses f markers, so it owns nw1-1 complete words once f > 0.
        self.n1 = WORD * nw1 - f
        self.n2 = f + WORD * self.nw2 + t
        assert self.n1 // WORD >= 5, "chr1 must keep 5 complete words"
        assert self.nw2 >= 5, "chr2 must keep 5 complete words"
        need = S.MIN_TOTAL_BP - (self.n1 - 1) * sp1
        auto = -(-need // (self.n2 - 1))
        self.s = spacing if spacing else int(-(-auto // 1000) * 1000)
        assert WORD * self.s <= 10_000_000, "a chr2 word spans over 10 Mb"
        assert (WORD * self.nw2 - 1) * self.s > 10_000_000, "chr2 segment too short"
        assert (WORD * (self.n1 // WORD) - 1) * sp1 > 10_000_000, "chr1 segment too short"

    # ---- geometry ----------------------------------------------------
    @property
    def nw2(self):
        return len(self.words)

    @property
    def f(self):
        return len(self.head)

    @property
    def t(self):
        return len(self.tail)

    def positions(self):
        p1 = [1_000_000 + self.sp1 * i for i in range(self.n1)]
        p2 = [1_000_000 + self.s * i for i in range(self.n2)]
        return p1, p2

    @property
    def denom(self):
        return (self.n1 - 1) * self.sp1 + (self.n2 - 1) * self.s

    def kinds2(self):
        """Every chr2 marker's kind, in chr2-local order."""
        out = list(self.head)
        for spec in self.words:
            out += S.expand(spec)
        return out + list(self.tail)

    # ---- build -------------------------------------------------------
    def build(self, wd):
        F.SPACING = 50_000
        fix = F.Fixture(self.name, [(1, self.n1), (2, self.n2)],
                        nsample=self.nsample, seed=self.seed, maf=0.5)
        pad = [0] * (self.nsample - 2)
        ck = S.expand(self.carrier)
        lo1, _ = fix.chrom_span(0)
        for i in range(self.n1):
            fix.pat_all[lo1 + i] = PAIR[ck[i % WORD]] + pad
            fix.noflip.add(lo1 + i)
        lo2, _ = fix.chrom_span(1)
        for i, kind in enumerate(self.kinds2()):
            g = PAIR[kind]
            fix.pat_all[lo2 + i] = (g[:self.nsample] if len(g) >= self.nsample
                                    else g + pad)
            fix.noflip.add(lo2 + i)
        prefix = fix.build(wd)
        p1, p2 = self.positions()
        out = []
        with open(prefix + ".bim") as fh:
            for n, line in enumerate(fh):
                v = line.rstrip("\n").split("\t")
                p = p1[n] if n < len(p1) else p2[n - len(p1)]
                v[2], v[3] = f"{p / 1e6:.6f}", str(p)
                out.append("\t".join(v))
        with open(prefix + ".bim", "w") as fh:
            fh.write("\n".join(out) + "\n")
        return prefix

    def key(self, extra):
        spec = json.dumps(["fringe", [S.expand(x) for x in self.words],
                           list(self.head), list(self.tail), self.nw1, self.sp1,
                           self.s, self.nsample, self.seed, S.expand(self.carrier),
                           list(extra)], sort_keys=True, default=str)
        return hashlib.sha1(spec.encode()).hexdigest()


# ---------------------------------------------------------------------------
# rigs: the block sits against the head (so a call opens on the grid's first word)
# or against the tail
# ---------------------------------------------------------------------------

def head_rig(name, block, head, nw2=NW2, **kw):
    """`block` at chr2 words 0.., walled off behind; `head` is the partial word."""
    words = list(block) + [WALL] * (nw2 - len(block))
    return FringeCanvas(name, words, head=head, **kw)


def tail_rig(name, block, tail, nw2=NW2, **kw):
    """`block` against chr2's last word, walled off in front."""
    words = [WALL] * (nw2 - len(block)) + list(block)
    return FringeCanvas(name, words, tail=tail, **kw)


# ---------------------------------------------------------------------------
# reading a result back, in chr2 marker intervals
# ---------------------------------------------------------------------------

def mk2(cv, res):
    """`IBD2Seg` as the number of chr2 marker intervals called."""
    row = res["row"]
    if row is None:
        return None
    return float(row["IBD2Seg"]) * cv.denom / cv.s


def raw1(cv, res):
    """`IBD1Seg` in chr2 marker intervals, carrier included (read differentially)."""
    row = res["row"]
    if row is None:
        return None
    return float(row["IBD1Seg"]) * cv.denom / cv.s


def show(tag, val, cands):
    """One measured value against named candidate predictions."""
    if val is None:
        return "  %-30s refused" % tag
    best = [k for k, v in cands.items() if v is not None and abs(val - v) <= 0.35]
    return ("  %-30s %8.3f   %s   [%s]"
            % (tag, val, "  ".join("%s=%s" % (k, v) for k, v in cands.items()),
               ",".join(best) if best else "NONE"))


# ---------------------------------------------------------------------------
# sections
# ---------------------------------------------------------------------------

def hdr(t):
    print("\n" + "=" * 78 + "\n" + t + "\n" + "=" * 78)


def run(items, extra=()):
    res = S.many([(cv, extra) for _lab, cv in items])
    return [(lab, cv, r) for (lab, cv), r in zip(items, res)]


NB = 8                 # callable words in the block
BASE = WORD * NB - 1   # markers a word-aligned call over NB words measures


def section0():
    """Controls: does a call leave the word grid at all, and does the ruler read true?"""
    hdr("0. controls — the instrument itself")
    cv = head_rig("f0_plain", [CLEAN] * NB, [])
    r = S.probe(cv)
    print("  no fringe, %d clean words then walls" % NB)
    print("    nseg %d (want 2)   D %d   Dref %d   one ulp = %.3f markers"
          % (r["nseg"], cv.denom, r["Dref"], (cv.denom / 1e4) / cv.s))
    print(show("IBD2Seg", mk2(cv, r), {"word-aligned": BASE}))

    print("\n  a clean head fringe of f markers — is it reached?")
    items = [("head=%2d clean" % f,
              head_rig("f0_h%d" % f, [CLEAN] * NB, ["zero"] * f)) for f in (1, 8, 32, 63)]
    for lab, cv, r in run(items):
        print(show(lab, mk2(cv, r), {"extend": BASE + cv.f, "none": BASE}))

    print("\n  a clean tail fringe of t markers")
    items = [("tail=%2d clean" % t,
              tail_rig("f0_t%d" % t, [CLEAN] * NB, ["zero"] * t)) for t in (1, 8, 32, 63)]
    for lab, cv, r in run(items):
        print(show(lab, mk2(cv, r), {"extend": BASE + cv.t, "none": BASE}))


def section1():
    """The head fringe: one mismatch at a known marker, swept."""
    hdr("1. the head fringe — one het-vs-hom mismatch, position swept")
    f = 32
    items = []
    for q in range(0, f, 2):
        head = ["zero"] * f
        head[q] = "ibs1"
        items.append((q, head_rig("f1_m%d" % q, [CLEAN] * NB, head)))
    print("  f = %d;  extend = %d (ignores the mismatch), none = %d" % (f, BASE + f, BASE))
    for q, cv, r in run(items):
        print(show("mismatch at head[%2d]" % q, mk2(cv, r),
                   {"off=1": BASE + f - q - 1, "off=0": BASE + f - q,
                    "extend": BASE + f, "none": BASE}))


def section2():
    """Is it the *last* mismatch, and does an opposite homozygote break it too?"""
    hdr("2. the head fringe — two breakers, and which kinds break")
    f = 32
    specs = [
        ("mismatches at 5 and 20", {5: "ibs1", 20: "ibs1"}),
        ("ibs0 at 20 only", {20: "ibs0"}),
        ("ibs0 at 20, mismatch at 5", {5: "ibs1", 20: "ibs0"}),
        ("mismatch at 20, ibs0 at 5", {5: "ibs0", 20: "ibs1"}),
        ("het-vs-A1A1 at 20", {20: "ibs1b"}),
        ("both missing at 20", {20: "miss"}),
        ("A1A1/A1A1 at 20", {20: "hom1"}),
    ]
    items = []
    for lab, at in specs:
        head = ["zero"] * f
        for i, k in at.items():
            head[i] = k
        items.append((lab, head_rig("f2_" + S.slug(lab), [CLEAN] * NB, head)))
    for lab, cv, r in run(items):
        print(show(lab, mk2(cv, r),
                   {"stop after 20": BASE + f - 21, "stop after 5": BASE + f - 6,
                    "extend": BASE + f, "none": BASE}))


def section3():
    """The tail fringe: the mirror."""
    hdr("3. the tail fringe — one mismatch, position swept")
    t = 32
    items = []
    for q in range(0, t, 2):
        tail = ["zero"] * t
        tail[q] = "ibs1"
        items.append((q, tail_rig("f3_m%d" % q, [CLEAN] * NB, tail)))
    print("  t = %d;  extend = %d, none = %d" % (t, BASE + t, BASE))
    for q, cv, r in run(items):
        print(show("mismatch at tail[%2d]" % q, mk2(cv, r),
                   {"off=1": BASE + q, "off=0": BASE + q + 1,
                    "extend": BASE + t, "none": BASE}))

    print("\n  two breakers, and which kinds break")
    specs = [("mismatches at 8 and 24", {8: "ibs1", 24: "ibs1"}),
             ("ibs0 at 8 only", {8: "ibs0"}),
             ("ibs0 at 24, mismatch at 8", {8: "ibs1", 24: "ibs0"}),
             ("het-vs-A1A1 at 8", {8: "ibs1b"})]
    items = []
    for lab, at in specs:
        tail = ["zero"] * t
        for i, k in at.items():
            tail[i] = k
        items.append((lab, tail_rig("f3b_" + S.slug(lab), [CLEAN] * NB, tail)))
    for lab, cv, r in run(items):
        print(show(lab, mk2(cv, r),
                   {"stop before 8": BASE + 8, "stop before 24": BASE + 24,
                    "extend": BASE + t, "none": BASE}))


def section4():
    """Both fringes at once, and whether the §5 reach *snaps* out to the fringe."""
    hdr("4. both fringes, and a run that opens at the second complete word")
    f = t = 24
    cv = FringeCanvas("f4_both", [CLEAN] * NB + [WALL] * (NW2 - NB),
                      head=["zero"] * f, tail=[])
    # A call can only touch both fringes when the block spans the whole segment.
    both = FringeCanvas("f4_span", [CLEAN] * NW2, head=["zero"] * f, tail=["zero"] * t)
    r = S.probe(both)
    span = WORD * NW2 - 1
    print(show("clean head AND tail, full span", mk2(both, r),
               {"both": span + f + t, "head only": span + f, "tail only": span + t,
                "neither": span}))
    both2 = FringeCanvas("f4_span_m", [CLEAN] * NW2,
                         head=["ibs1" if i == 10 else "zero" for i in range(f)],
                         tail=["ibs1" if i == 9 else "zero" for i in range(t)])
    r = S.probe(both2)
    print(show("mismatch head[10] and tail[9]", mk2(both2, r),
               {"both stop": span + (f - 11) + 9, "extend both": span + f + t,
                "neither": span}))

    print("\n  a run that cannot start at block word 0 (2 mismatches there, no IBS0).")
    print("  The §5 reach carries its left end back over word 0 to the grid's own edge,")
    print("  where 'snap' hands it to the fringe and 'stop'/§5-as-written does not.")
    dirty = S.w(m=2, h=62)
    right = f + WORD * NB - 1       # the run ends at block word NB-1's last marker
    for lab, head in (("clean head", ["zero"] * f),
                      ("mismatch at head[10]",
                       ["ibs1" if i == 10 else "zero" for i in range(f)])):
        cv = head_rig("f4_reach_" + S.slug(lab), [dirty] + [CLEAN] * (NB - 1), head)
        r = S.probe(cv)
        stop = 0 if "clean" in lab else 11
        print(show("run from word 1, " + lab, mk2(cv, r),
                   {"snap to fringe": right - stop, "stop at grid": right - f,
                    "no reach (word 1 start)": right - f - WORD}))


def section5():
    """The IBD1 fringe — the same questions, with IBD2 silenced by paint."""
    hdr("5. the IBD1 fringe (IBD2 held at zero by 34 mismatches per word)")
    f = 32
    specs = [("clean head", {}),
             ("ibs0 at head[ 0]", {0: "ibs0"}),
             ("ibs0 at head[ 8]", {8: "ibs0"}),
             ("ibs0 at head[20]", {20: "ibs0"}),
             ("ibs0 at head[31]", {31: "ibs0"}),
             ("ibs0 at 5 and 20", {5: "ibs0", 20: "ibs0"}),
             ("mismatch at head[20]", {20: "ibs1"})]
    items = []
    for lab, at in specs:
        head = ["zero"] * f
        for i, k in at.items():
            head[i] = k
        items.append((lab, head_rig("f5_" + S.slug(lab), [K] * NB, head)))
    got = run(items)
    # Read differentially against the walls-only canvas: same chr1, no chr2 IBD1 at all.
    zero = head_rig("f5_zero", [WALL] * NB, ["zero"] * f)
    rz = S.probe(zero)
    b0 = raw1(zero, rz)
    # An IBD1 run reaches into the word that ended it, out to that word's *last* IBS0
    # (`18-…`), and here that word is a full wall — so the call swallows it whole and the
    # right end is one word past the block.
    top = f + WORD * (NB + 1) - 1
    print("  carrier baseline (chr2 all walls): IBD1Seg reads %8.3f markers" % b0)
    for (lab, at), (_l, cv, r) in zip(specs, got):
        i2 = float(r["row"]["IBD2Seg"]) if r["row"] else -1
        v = raw1(cv, r)
        # Candidates built from this fringe's own breakers, so the row is self-checking.
        cands = {"extend": top, "none": top - f}
        for q in sorted(at):
            cands["stop after %d" % q] = top - (q + 1)
        print(show(lab + "  (IBD2 %.4f)" % i2, None if v is None else v - b0, cands))

    print("\n  the tail, same rig (chr1 word-aligned here, so the baseline is clean)")
    t = 32
    specs = [("clean tail", {}),
             ("ibs0 at tail[ 0]", {0: "ibs0"}),
             ("ibs0 at tail[12]", {12: "ibs0"}),
             ("ibs0 at tail[24]", {24: "ibs0"}),
             ("ibs0 at tail[31]", {31: "ibs0"}),
             ("ibs0 at 12 and 24", {12: "ibs0", 24: "ibs0"}),
             ("mismatch at tail[12]", {12: "ibs1"})]
    items = []
    for lab, at in specs:
        tail = ["zero"] * t
        for i, k in at.items():
            tail[i] = k
        items.append((lab, tail_rig("f5t_" + S.slug(lab), [K] * NB, tail)))
    got = run(items)
    zero = tail_rig("f5t_zero", [WALL] * NB, ["zero"] * t)
    rz = S.probe(zero)
    b0 = raw1(zero, rz)
    print("  carrier baseline (chr2 all walls): IBD1Seg reads %8.3f markers" % b0)
    for (lab, at), (_l, cv, r) in zip(specs, got):
        i2 = float(r["row"]["IBD2Seg"]) if r["row"] else -1
        v = raw1(cv, r)
        cands = {"extend": BASE + t, "none": BASE}
        for q in sorted(at):
            cands["stop before %d" % q] = BASE + q
        print(show(lab + "  (IBD2 %.4f)" % i2, None if v is None else v - b0, cands))


# ---------------------------------------------------------------------------
# §6 — grading the rule out of sample
# ---------------------------------------------------------------------------

def _fit_dir():
    repo = os.path.dirname(os.path.dirname(os.path.dirname(ROOT)))
    return os.path.join(repo, "tests", "parity", "fit")


sys.path.insert(0, _fit_dir())


def wordinfo(spec):
    """`(ibs0?, mismatches, first mismatch bit, last mismatch bit, inf2)` of one word."""
    ks = S.expand(spec)
    mis = [b for b, k in enumerate(ks) if k in ("ibs1", "ibs1b")]
    return (any(k in ("ibs0", "ibs0b") for k in ks), len(mis),
            mis[0] if mis else None, mis[-1] if mis else None,
            sum(1 for k in ks if k in ("hethet", "hom1", "ibs1b")))


def fringe_mask(kinds, side, f_or_t):
    """The mismatch bit mask of a fringe, in the *global word's* bit numbering."""
    m = 0
    for i, k in enumerate(kinds):
        if k in ("ibs1", "ibs1b"):
            b = (WORD - f_or_t + i) if side == "head" else i
            m |= 1 << b
    return m


def predicted_markers(cv):
    """`seg19.predict` applied to a FringeCanvas, as chr2 marker intervals called."""
    import seg19 as S19
    w0, w1 = cv.nw1, cv.nw1 + cv.nw2 - 1
    lo = WORD * cv.nw1 - cv.f
    hi = WORD * (cv.nw1 + cv.nw2) + cv.t - 1
    # `predict` indexes `info` by absolute word number, so the words before the segment
    # are present but never read.
    info = [(False, 0, None, None, 0)] * w0 + [wordinfo(spec) for spec in cv.words]
    calls = S19.predict(info, w0=w0, w1=w1, lo=lo, hi=hi,
                        head=fringe_mask(cv.head, "head", cv.f),
                        tail=fringe_mask(cv.tail, "tail", cv.t))
    return sum(b - a for a, b in calls)


def random_word(rng):
    """A word from the alphabet the two other rigs use, with IBS0 kept rare."""
    r = rng.random()
    if r < 0.18:
        return WALL if rng.random() < 0.4 else S.w(m=rng.randint(0, 3),
                                                   h=rng.randint(0, 20), z=1)
    m = rng.choice([0, 0, 1, 1, 2, 3, 5])
    h = rng.randint(0, WORD - m)
    rest = WORD - m - h
    hom = rng.randint(0, min(rest, 12))
    return {"ibs1": m, "hethet": h, "hom1": hom, "zero": rest - hom}


def random_fringe(rng, n):
    out = []
    for _ in range(n):
        r = rng.random()
        out.append("ibs1" if r < 0.12 else "ibs0" if r < 0.18
                   else "hethet" if r < 0.5 else "zero")
    return out


def battery(seed, count, width=10, verbose=False):
    """Random fringe canvases: the reference against `seg19.predict`."""
    import random
    rng = random.Random(seed)
    items = []
    for c in range(count):
        f = rng.randint(0, 40)
        t = rng.randint(0, 40)
        words = [random_word(rng) for _ in range(width)]
        words += [WALL] * (NW2 - width)
        cv = FringeCanvas("f6_%d_%d" % (seed, c), words,
                          head=random_fringe(rng, f), tail=random_fringe(rng, t))
        items.append(cv)
    res = S.many([(cv, ()) for cv in items])
    ok = bad = 0
    for cv, r in zip(items, res):
        got = mk2(cv, r)
        want = predicted_markers(cv)
        if got is None:
            got = 0.0
        if abs(got - want) <= 0.35:
            ok += 1
        else:
            bad += 1
            if verbose:
                print("    MISS %-14s reference %8.3f  predict %6d  f=%d t=%d"
                      % (cv.name, got, want, cv.f, cv.t))
    print("  seed %-5d %4d canvases:  %4d agree, %3d differ" % (seed, count, ok, bad))
    return ok, bad


def section6():
    """Out of sample: fresh random fringe canvases, and an exhaustive short sweep."""
    hdr("6. out-of-sample batteries — the reference against `seg19.predict`")
    print("  random canvases with random head/tail fringes (unused seeds):")
    tot = Counter()
    for seed in (4919, 6271):
        ok, bad = battery(seed, 60, verbose=True)
        tot["ok"] += ok
        tot["bad"] += bad
        S.save_cache()

    print("\n  exhaustive: every 3-word block over {clean, quiet, 1-mis, 2-mis}")
    print("  crossed with 6 fringe shapes on each side")
    alpha = {"C": CLEAN, "q": {"zero": WORD}, "y": S.w(m=1, h=63),
             "d": S.w(m=2, h=62)}
    shapes = {"-": [], "z8": ["zero"] * 8, "m8": ["ibs1"] + ["zero"] * 7,
              "z20": ["zero"] * 20, "m20": ["zero"] * 10 + ["ibs1"] + ["zero"] * 9,
              "x20": ["ibs0"] + ["zero"] * 19}
    items = []
    for a in alpha:
        for b in alpha:
            for c in alpha:
                for hs, hv in shapes.items():
                    body = [alpha[a], alpha[b], alpha[c]] + [WALL] * (NW2 - 3)
                    items.append(FringeCanvas("f6e_%s%s%s_%s" % (a, b, c, hs), body,
                                              head=hv, tail=[]))
    res = S.many([(cv, ()) for cv in items])
    ok = bad = 0
    for cv, r in zip(items, res):
        got = mk2(cv, r)
        want = predicted_markers(cv)
        if got is None:
            got = 0.0
        if abs(got - want) <= 0.35:
            ok += 1
        else:
            bad += 1
            print("    MISS %-18s reference %8.3f  predict %6d"
                  % (cv.name, got, want))
    tot["ok"] += ok
    tot["bad"] += bad
    print("  %d compositions x fringes: %d agree, %d differ" % (len(items), ok, bad))
    print("\n  TOTAL %d canvases: %d agree, %d differ" % (tot["ok"] + tot["bad"],
                                                          tot["ok"], tot["bad"]))
    S.save_cache()


def section7():
    """The Rust port against the same rig: our binary must read what the reference reads."""
    hdr("7. the port — `target/release/open-king` on the fringe canvases")
    repo = os.path.dirname(os.path.dirname(os.path.dirname(ROOT)))
    ours = os.environ.get("OPENKING", os.path.join(repo, "target", "release", "open-king"))
    if not os.path.exists(ours):
        print("  %s not built — skipping" % ours)
        return

    # Everything §0-§4 and §6's exhaustive sweep put in front of the reference.
    alpha = {"C": CLEAN, "q": {"zero": WORD}, "y": S.w(m=1, h=63), "d": S.w(m=2, h=62)}
    shapes = {"-": [], "z8": ["zero"] * 8, "m8": ["ibs1"] + ["zero"] * 7,
              "z20": ["zero"] * 20, "m20": ["zero"] * 10 + ["ibs1"] + ["zero"] * 9,
              "x20": ["ibs0"] + ["zero"] * 19}
    items = []
    for q in range(0, 32, 2):
        items.append(head_rig("f1_m%d" % q, [CLEAN] * NB,
                              ["ibs1" if i == q else "zero" for i in range(32)]))
        items.append(tail_rig("f3_m%d" % q, [CLEAN] * NB,
                              ["ibs1" if i == q else "zero" for i in range(32)]))
    for a in alpha:
        for b in alpha:
            for c in alpha:
                for hs, hv in shapes.items():
                    items.append(FringeCanvas(
                        "f6e_%s%s%s_%s" % (a, b, c, hs),
                        [alpha[a], alpha[b], alpha[c]] + [WALL] * (NW2 - 3),
                        head=hv, tail=[]))

    ref = [mk2(cv, r) for cv, r in zip(items, S.many([(cv, ()) for cv in items]))]
    S.save_cache()

    # Our binary's answers are deliberately **not** cached: a stale cache here would
    # report agreement after an engine change that broke the rule.  Every run of §7 puts
    # the current build in front of every canvas again.
    keep = (F.KING, S.CACHE, S._cache, S.USE_CACHE)
    F.KING, S.USE_CACHE = ours, False
    S.CACHE = os.path.join(S.WORK, "ours-scratch.json")
    S._cache = None
    try:
        got = [mk2(cv, r) for cv, r in zip(items, S.many([(cv, ()) for cv in items]))]
    finally:
        F.KING, S.CACHE, S._cache, S.USE_CACHE = keep
        shutil.rmtree(os.path.join(S.WORK, "ours-scratch.json"), ignore_errors=True)

    ok = bad = 0
    for cv, a, b in zip(items, ref, got):
        a = 0.0 if a is None else a
        b = 0.0 if b is None else b
        if abs(a - b) <= 0.35:
            ok += 1
        else:
            bad += 1
            print("    MISS %-18s reference %8.3f  ours %8.3f" % (cv.name, a, b))
    print("  %d canvases: %d agree, %d differ" % (len(items), ok, bad))


SECTIONS = {0: section0, 1: section1, 2: section2, 3: section3, 4: section4,
            5: section5, 6: section6, 7: section7}


def main(argv):
    want = [int(a) for a in argv if a.isdigit()] or sorted(SECTIONS)
    for k in want:
        SECTIONS[k]()
        S.save_cache()


if __name__ == "__main__":
    main(sys.argv[1:])
