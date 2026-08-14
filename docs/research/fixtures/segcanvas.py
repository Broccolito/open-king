#!/usr/bin/env python3
"""The `.seg`-native canvas — fixtures behind `docs/research/17-seg-caller.md`.

`docs/research/16-segment-extension.md` §10 named the blocker: the `--ibs` rig varies
het-vs-hom mismatches, and the `.seg` IBD2 pass refuses any word that carries an
**opposite homozygote**, so every `--ibs` fixture reports `IBD2Seg 0.0000` under
`--ibdseg` and the `.seg` constants cannot be measured through it.  This rig is built the
other way round: IBS0 is the paint, and every block meant to be *called* carries none.

The instrument
--------------

* **chr1** — the carrier.  `nw1` complete words that are IBD1-clean (no IBS0) and
  IBD2-dirty (20 het-vs-hom mismatches per word), so the pair always owns one IBD1
  segment over 10 Mb — which is what earns it a `.seg` row — and contributes nothing to
  `IBD2Seg`.
* **chr2** — the canvas.  `nw2` complete words painted marker by marker from an explicit
  composition, walled at both ends by all-IBS0 words, and padded out to a fixed word
  count so `D` is the same across a family.
* **the ruler** — chr2's marker spacing `s` is uniform and chosen so `D` lands just over
  the reference's 100 Mb floor (§0).  One ulp of the printed `%.4lf` is then `D/10000`,
  about a fifth of one marker gap, so `IBD2Seg * D / s` reads back the **number of marker
  intervals called** to better than a tenth of a marker.  A word-aligned call over `n`
  words contributes `64n - 1`, so a total `M` decodes to `(words, calls)` by
  `calls = (-M) mod 64` — the count of calls *and* the count of words, exactly.
* **`--seglength`** separates the calls when there is more than one: `IBD2Seg(L)` is a
  step function whose jumps are individual segment lengths.

Nothing here reads KING's source.

    python3 segcanvas.py            # every section
    python3 segcanvas.py 1 3        # only sections 1 and 3
    SEGCANVAS_JOBS=12 python3 segcanvas.py
    KING=../../../target/release/king python3 segcanvas.py 9   # our binary, same rig
"""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import random
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixlab as F  # noqa: E402

WORD = 64
ROOT = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(ROOT, "work", "segcanvas")
CACHE = os.path.join(ROOT, "segcanvas_measured.json")
JOBS = int(os.environ.get("SEGCANVAS_JOBS", "8"))
USE_CACHE = os.environ.get("SEGCANVAS_NOCACHE", "") == ""

# The reference refuses to analyse a fileset whose usable total is under this. Bisected
# to the base pair in §0: 100000000 is analysed, 99999999 prints "Segments too short."
MIN_TOTAL_BP = 100_000_000

# per-marker patterns, as (sample 0, sample 1) A1 dosages; 3 is a missing call
PAIR = {
    "hethet": [1, 1],   # both het        — HetHet, and `inf2`-informative
    "hom1": [2, 2],     # both A1A1       — `inf1`/`inf2`-informative, not HetHet
    "zero": [0, 0],     # both A2A2       — clean and worth nothing
    "ibs1": [1, 0],     # het vs A2A2     — het-vs-hom mismatch, uninformative
    "ibs1b": [1, 2],    # het vs A1A1     — het-vs-hom mismatch, `inf2`-informative
    "ibs0": [2, 0],     # opposite homs   — THE paint
    "ibs0b": [0, 2],    # opposite homs, other polarity
    "miss": [3, 3],     # both missing
    "missA": [3, 1],    # one missing, one het
    # Polymorphic controls: the pair reads exactly as `zero`/`hethet` above, but the
    # padding samples carry A1, so the marker is not monomorphic in the cohort. Used to
    # rule out any special handling of monomorphic markers.
    "zeroP": [0, 0, 1, 1, 0, 0],
    "hetP": [1, 1, 1, 1, 0, 0],
}

WALL = {"ibs0": 64}                  # all opposite homozygotes: dirty for IBD1 and IBD2
CLEAN = {"hethet": 64}               # the canonical callable IBD2 word
HOM = {"hom1": 64}                   # clean, `inf2`-informative, no HetHet
# The carrier: IBD1-clean (no IBS0), IBD2-dirty (34 mismatches per word), and just
# `inf1`-rich enough to clear the gate (4 x 5 words = 20 markers).  A1A1/A1A1 markers are
# kept scarce on purpose: KING's "Too many first alleles as the major allele" QC samples
# individuals with an unseeded RNG, and a marker where the test pair is A1A1/A1A1 reads as
# A1-major whenever the pair is what it sampled.
CARRIER = {"ibs1": 34, "hom1": 4, "zero": 26}


def w(m=0, h=0, z=0, **kw):
    """A word carrying `m` het-vs-hom mismatches, `h` HetHet and `z` IBS0; rest A2A2."""
    d = {"ibs1": m, "hethet": h, "ibs0": z}
    d.update(kw)
    return {k: v for k, v in d.items() if v}


def at(spec, **places):
    """`spec` as an explicit 64-list with `kind=bit` overrides, e.g. `at(CLEAN, ibs0=7)`."""
    ks = expand(spec)
    for kind, bits in places.items():
        for b in (bits if isinstance(bits, (list, tuple)) else [bits]):
            ks[b] = kind
    return ks


def expand(spec):
    """A composition dict (or explicit 64-list) as a list of 64 marker kinds."""
    if not isinstance(spec, dict):
        ks = list(spec)
        assert len(ks) == WORD, len(ks)
        return ks
    ks = []
    for kind, cnt in spec.items():
        ks.extend([kind] * cnt)
    assert len(ks) <= WORD, spec
    return (ks + ["zero"] * WORD)[:WORD]


def counts(spec):
    """`(mis, het, ibs0, inf2, inf1)` of one word composition."""
    ks = expand(spec)
    return (sum(1 for k in ks if k in ("ibs1", "ibs1b")),
            sum(1 for k in ks if k == "hethet"),
            sum(1 for k in ks if k in ("ibs0", "ibs0b")),
            sum(1 for k in ks if k in ("hethet", "hom1", "ibs1b")),
            sum(1 for k in ks if k in ("hom1", "ibs1b")))


class Canvas:
    """One fixture: a carrier chromosome plus a painted chr2.

    `block` is the interesting part; `pad[0]` walls precede it and the rest of chr2 is
    walls, out to `nw2` words total.  Block word 0 is chr2 word `pad[0]`.
    """

    def __init__(self, name, block, nw1=5, sp1=33_000, spacing=None,
                 pad=(3, 3), nw2=16, nsample=6, seed=5, carrier=None, gaps=None):
        body = [WALL] * pad[0] + list(block) + [WALL] * pad[1]
        assert len(body) <= nw2, (name, len(body), nw2)
        body += [WALL] * (nw2 - len(body))
        self.name = name
        self.pad = pad
        self.words = body
        self.nw1 = nw1
        self.sp1 = sp1
        self.nsample = nsample
        self.seed = seed
        self.carrier = carrier or CARRIER
        # The smallest uniform chr2 spacing that clears the 100 Mb floor, rounded up to
        # a round number so the ruler is legible; capped so no word spans over 10 Mb.
        need = MIN_TOTAL_BP - (WORD * nw1 - 1) * sp1
        auto = -(-need // (WORD * self.nw2 - 1))
        self.s = spacing if spacing else int(-(-auto // 1000) * 1000)
        # An optional per-word gap vector: the ruler stops being uniform, so a marker
        # count no longer decodes on its own, but two candidate intervals of the same
        # width become distinguishable by length. `decode` is the reader for it.
        self.gaps = list(gaps) if gaps else None
        assert WORD * max(self.gaps or [self.s]) <= 10_000_000, "a word spans over 10 Mb"

    # ---- geometry ----------------------------------------------------
    @property
    def nw2(self):
        return len(self.words)

    def positions(self):
        p1 = [1_000_000 + self.sp1 * i for i in range(WORD * self.nw1)]
        if self.gaps is None:
            p2 = [1_000_000 + self.s * i for i in range(WORD * self.nw2)]
        else:
            p2, y = [], 1_000_000
            for k in range(self.nw2):
                for _ in range(WORD):
                    p2.append(y)
                    y += self.gaps[k]
        return p1, p2

    @property
    def denom(self):
        _, p2 = self.positions()
        return (WORD * self.nw1 - 1) * self.sp1 + p2[-1] - p2[0]

    def wspan(self, u, e):
        """bp of the word-aligned interval `[u, e]` in *block* word coordinates."""
        _, p2 = self.positions()
        o = WORD * self.pad[0]
        return p2[WORD * e + WORD - 1 + o] - p2[WORD * u + o]

    def wbp(self, calls):
        """bp a list of block word intervals would measure."""
        return sum(self.wspan(u, e) for u, e in calls)

    # ---- build -------------------------------------------------------
    def build(self, wd):
        F.SPACING = 50_000
        fix = F.Fixture(self.name, [(1, WORD * self.nw1), (2, WORD * self.nw2)],
                        nsample=self.nsample, seed=self.seed, maf=0.5)
        pad = [0] * (self.nsample - 2)
        lo1, _ = fix.chrom_span(0)
        ck = expand(self.carrier)
        for k in range(self.nw1):
            for idx, kind in enumerate(ck):
                fix.pat_all[lo1 + WORD * k + idx] = PAIR[kind] + pad
                fix.noflip.add(lo1 + WORD * k + idx)
        lo2, _ = fix.chrom_span(1)
        for k, spec in enumerate(self.words):
            for idx, kind in enumerate(expand(spec)):
                g = PAIR[kind]
                fix.pat_all[lo2 + WORD * k + idx] = (
                    g[:self.nsample] if len(g) >= self.nsample else g + pad)
                fix.noflip.add(lo2 + WORD * k + idx)
        prefix = fix.build(wd)
        p1, p2 = self.positions()
        out = []
        with open(prefix + ".bim") as f:
            for n, line in enumerate(f):
                v = line.rstrip("\n").split("\t")
                p = p1[n] if n < len(p1) else p2[n - len(p1)]
                v[2], v[3] = f"{p / 1e6:.6f}", str(p)
                out.append("\t".join(v))
        with open(prefix + ".bim", "w") as f:
            f.write("\n".join(out) + "\n")
        return prefix

    def key(self, extra):
        spec = json.dumps([[expand(x) for x in self.words], self.nw1, self.sp1, self.s,
                           self.gaps,
                           self.pad, self.nsample, self.seed, expand(self.carrier),
                           list(extra)], sort_keys=True, default=str)
        return hashlib.sha1(spec.encode()).hexdigest()

    # ---- the model's view of this canvas -----------------------------
    def wordstats(self):
        """Per-word `(mis, het, ibs0, inf2, inf1)` over all of chr2."""
        return [counts(x) for x in self.words]


# ---------------------------------------------------------------------------
# running the reference
# ---------------------------------------------------------------------------

_cache = None
_dirty = False


def cache():
    global _cache
    if _cache is None:
        try:
            with open(CACHE) as f:
                _cache = json.load(f)
        except (OSError, ValueError):
            _cache = {}
    return _cache


def save_cache():
    global _dirty
    if _dirty and _cache is not None:
        tmp = CACHE + ".tmp"
        with open(tmp, "w") as f:
            json.dump(_cache, f, indent=0, sort_keys=True)
        os.replace(tmp, CACHE)
        _dirty = False


def slug(s):
    return "".join(c if c.isalnum() else "_" for c in str(s))[:100]


def _invoke(cv, extra, keep=False):
    wd = os.path.join(WORK, slug(cv.name) + "_" + cv.key(extra)[:8])
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd, exist_ok=True)
    prefix = cv.build(wd)
    cmd = [F.KING, "-b", prefix + ".bed", "--ibdseg", *extra,
           "--prefix", os.path.join(wd, "k")]
    out, rows, segs = "", {}, []
    for attempt in range(24):
        if attempt:
            time.sleep(1.05)        # the QC gate's RNG is seeded from the clock

        r = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
        out = r.stdout
        segs = F.parse_allsegs(os.path.join(wd, "kallsegs.txt"))
        rows = F.parse_seg(os.path.join(wd, "k.seg"))
        # The carrier guarantees a >10 Mb IBD1 segment, so a missing row means the run
        # itself failed (a fatal on the unseeded A1-major gate, or a truncated write).
        if "FATAL ERROR" not in out and (("S00", "S01") in rows
                                         or "Segments too short." in out):
            break
    if "FATAL ERROR" in out or (("S00", "S01") not in rows
                                and "Segments too short." not in out):
        raise RuntimeError("KING produced no row for %s" % cv.name)
    res = dict(row=rows.get(("S00", "S01")), nseg=len(segs),
               Dref=int(round(sum(float(s["Length"]) for s in segs
                                  if int(s["Chr"]) < 23) * 1e6)))
    if keep:
        res["wd"] = wd
        res["stdout"] = out
    else:
        shutil.rmtree(wd, ignore_errors=True)
    return res


def probe(cv, extra=(), keep=False):
    extra = tuple(extra)
    k = cv.key(extra)
    if USE_CACHE and not keep and k in cache():
        return cache()[k]
    res = _invoke(cv, extra, keep)
    if not keep:
        global _dirty
        cache()[k] = res
        _dirty = True
    return res


def many(items):
    """`[(canvas, extra), ...]` -> results in order; cached, and run in parallel."""
    items = [(cv, tuple(ex)) for cv, ex in items]
    todo = [(cv, ex) for cv, ex in items
            if not (USE_CACHE and cv.key(ex) in cache())]
    seen, uniq = set(), []
    for cv, ex in todo:
        if cv.key(ex) not in seen:
            seen.add(cv.key(ex))
            uniq.append((cv, ex))
    if uniq:
        with ThreadPoolExecutor(max_workers=JOBS) as pool:
            got = list(pool.map(lambda t: _invoke(t[0], t[1]), uniq))
        global _dirty
        for (cv, ex), res in zip(uniq, got):
            cache()[cv.key(ex)] = res
            _dirty = True
        save_cache()
    return [cache()[cv.key(ex)] for cv, ex in items]


# ---------------------------------------------------------------------------
# reading a result back
# ---------------------------------------------------------------------------

def mk(cv, res, what=2):
    """Marker intervals called on chr2, as a float. `what` is 1 for IBD1, 2 for IBD2."""
    row = res["row"]
    if row is None:
        return None
    v = float(row["IBD1Seg" if what == 1 else "IBD2Seg"]) * cv.denom
    if what == 1:                       # IBD1Seg always carries the carrier chromosome
        v -= (WORD * cv.nw1 - 1) * cv.sp1
    return v / cv.s


def wc(m, tol=0.25):
    """`(words, calls)` behind a marker total `m`, or None if it is not word-aligned.

    A word-aligned call over `n` words measures `64n - 1` markers, so a total of `c`
    calls over `n` words总 measures `64n - c`.
    """
    if m is None:
        return None
    r = int(round(m))
    if abs(m - r) > tol:
        return None
    calls = (-r) % WORD
    words = (r + calls) // WORD
    return (words, calls)


def pick(cv, res, cands, what=2):
    """Which of `cands` (each a list of block word intervals) the printed total is."""
    row = res["row"]
    if row is None:
        return [c for c in cands if not c]
    want = float(row["IBD1Seg" if what == 1 else "IBD2Seg"]) * cv.denom
    if what == 1:
        want -= (WORD * cv.nw1 - 1) * cv.sp1
    tol = cv.denom / 2e4 + 1
    return [c for c in cands if abs(cv.wbp(c) - want) <= tol]


def fmt(cv, res, what=2):
    m = mk(cv, res, what)
    if m is None:
        return "refused"
    d = wc(m)
    return "%8.2f mk %s" % (m, "" if d is None else "= %d words / %d calls" % d)


def seglen_lengths(cv, extra=(), lo=1.0, hi=10.0, what=2, steps=None):
    """Individual call lengths, in markers, by sweeping `--seglength` over a grid.

    A jump in `IBD2Seg` between `L` and `L'` means a call of length in `[L, L')` was
    dropped; with the grid fine enough (one marker gap) each jump names one call.
    """
    s = cv.s
    if steps is None:
        steps = []
        L = int(lo * 1e6)
        while L <= hi * 1e6:
            steps.append(L)
            L += s
    res = many([(cv, tuple(extra) + ("--seglength", "%.6f" % (L / 1e6)))
                for L in steps])
    vals = [mk(cv, r, what) for r in res]
    out = []
    for (L, a), (L2, b) in zip(zip(steps, vals), zip(steps[1:], vals[1:])):
        if a is None or b is None:
            continue
        if abs(a - b) > 0.5:
            out.append((round(a - b), L2 / s))
    return out, list(zip(steps, vals))


# ---------------------------------------------------------------------------
# sections
# ---------------------------------------------------------------------------

def hdr(t):
    print("\n== %s" % t)


def section0():
    hdr("0. the rig — the 100 Mb floor, the denominator, the marker ruler")
    print("    The reference refuses a fileset whose usable total is under 100 Mb")
    print("    ('Segments too short.'), bisected to the base pair by solving")
    print("    D = 319*sp1 + 895*spacing for exact totals:")
    for tgt, sp1, sp in ((99_999_999, 32_141, 100_276),
                         (100_000_000, 32_040, 100_312),
                         (100_000_001, 32_834, 100_029)):
        cv = Canvas("floor_%d" % tgt, [CLEAN] * 8, nw1=5, sp1=sp1, spacing=sp,
                    pad=(3, 3), nw2=14)
        r = probe(cv)
        assert cv.denom == tgt, (cv.denom, tgt)
        print("      D = %d  ->  %s" % (tgt, "row" if r["row"] else "no row"))
    empty = Canvas("rig_empty", [])
    hh = Canvas("rig_hh8", [CLEAN] * 8)
    r = many([(empty, ()), (hh, ())])
    print("    geometry: chr1 %d words @ %d bp, chr2 %d words @ %d bp, D = %d"
          % (empty.nw1, empty.sp1, empty.nw2, empty.s, empty.denom))
    print("    one ulp of IBD2Seg = %.0f bp = %.3f marker gaps"
          % (empty.denom / 1e4, empty.denom / 1e4 / empty.s))
    print("    all-wall canvas : IBD1 %s  IBD2 %s   (carrier only, chr2 silent)"
          % (r[0]["row"]["IBD1Seg"], r[0]["row"]["IBD2Seg"]))
    print("    8 HetHet words  : IBD1 %s  IBD2 %s -> %s"
          % (r[1]["row"]["IBD1Seg"], r[1]["row"]["IBD2Seg"], fmt(hh, r[1])))
    print("    ...and its IBD1 column reads %s: a pure-HetHet block carries no `inf1`"
          % fmt(hh, r[1], 1))
    print("       at all, so the IBD1 pass never calls it. Use CARRIER words (section 3)")
    print("       to read IBD1 back.")


def section1():
    hdr("1. a pure-HetHet block of W words, walled by all-IBS0 words")
    ws = list(range(1, 11))
    cvs = [Canvas("hhW%d" % W, [CLEAN] * W) for W in ws]
    res = many([(c, ("--seglength", "1")) for c in cvs])
    print("    W    IBD2Seg   markers called   decoded            expected 64W-1")
    for W, cv, r in zip(ws, cvs, res):
        row = r["row"]
        print("    %-4d %-9s %s   %d"
              % (W, "-" if row is None else row["IBD2Seg"], fmt(cv, r), 64 * W - 1))


def _edge(name, spec, side, body=6, filler=None, **kw):
    """A block whose `side` end is the test word, the rest `filler` (default HetHet)."""
    b = [filler or CLEAN] * body
    return Canvas(name, ([spec] + b) if side == "L" else (b + [spec]), **kw)


def section2():
    hdr("2. the IBD2 word predicate")
    print("    8 pure-HetHet words with `j` consecutive words replaced by")
    print("    (m mismatches, 64-m HetHet), read back as words called / calls:")
    ms = (0, 1, 2, 3, 4, 5, 8, 16, 64)
    items = []
    for j in (1, 2, 3):
        for m in ms:
            body = [CLEAN] * 8
            for t in range(3, 3 + j):
                body[t] = w(m=m, h=64 - m)
            items.append(((j, m), Canvas("jm_%d_%d" % (j, m), body)))
    res = many([(cv, ("--seglength", "1")) for _, cv in items])
    d = {k: wc(mk(cv, r)) for (k, cv), r in zip(items, res)}
    print("     j\\m  " + "".join("%-8d" % m for m in ms))
    for j in (1, 2, 3):
        print("     %-4d " % j
              + "".join("%-8s" % ("%d/%d" % d[(j, m)]) for m in ms))
    print("    -> a word is unusable at DIRTY_MIS=%d mismatches. A *lone* unusable word"
          % DIRTY_MIS)
    print("       is absorbed whatever it carries (the j = 1 row); two in a row never")
    print("       are, and both are then dropped from the call.")
    print()
    print("    the same 8-word block with words 3 and 4 replaced by other content:")
    trials = []
    for z in (1, 2, 64):
        trials.append(("ibs0=%d" % z, w(h=64 - z, z=z)))
    for q in (1, 2, 3):
        trials.append(("ibs1b=%d" % q, {"ibs1b": q, "hethet": 64 - q}))
    for q in (32, 64):
        trials.append(("miss=%d" % q, {"miss": q, "hethet": 64 - q}))
    trials.append(("hom1=64", HOM))
    trials.append(("zero=64", {"zero": 64}))
    trials.append(("ibs1=1 + ibs0=1", w(m=1, h=62, z=1)))
    cvs = []
    for lab, spec in trials:
        body = [CLEAN] * 8
        body[3] = body[4] = spec
        cvs.append((lab, Canvas("pr2_" + slug(lab), body)))
    res = many([(cv, ("--seglength", "1")) for _, cv in cvs])
    for (lab, cv), r in zip(cvs, res):
        print("    %-18s %s" % (lab, fmt(cv, r)))


def section3():
    hdr("3. the IBD1 word predicate and endpoints")
    print("    block = [X] + 6 carrier words (IBD1-clean, IBD2-dirty), walled;")
    print("    IBD1Seg reads the chr2 call directly because IBD2 is empty.")
    trials = [("clean", CARRIER),
              ("ibs0=1", {"ibs0": 1, "ibs1": 34, "hom1": 4}),
              ("ibs0=2", {"ibs0": 2, "ibs1": 34, "hom1": 4}),
              ("ibs0=64", WALL),
              ("ibs1=34", {"ibs1": 34, "hom1": 4}),
              ("ibs1=64", {"ibs1": 64}),
              ("miss=64", {"miss": 64}),
              ("zero=64", {"zero": 64})]
    cvs = [(lab, _edge("p1L_" + slug(lab), spec, "L", body=6, filler=CARRIER),
            _edge("p1R_" + slug(lab), spec, "R", body=6, filler=CARRIER))
           for lab, spec in trials]
    res = many([(c, ("--seglength", "1")) for _, a, b in cvs for c in (a, b)])
    print("    %-14s %-30s %s" % ("X", "left edge (IBD1)", "right edge (IBD1)"))
    for n, (lab, a, b) in enumerate(cvs):
        print("    %-14s %-30s %s"
              % (lab, fmt(a, res[2 * n], 1), fmt(b, res[2 * n + 1], 1)))


def section4():
    hdr("4. the endpoints — how far a call reaches into the words that bound it")
    print("    block = 6 HetHet words with one boundary word beside them; the boundary")
    print("    word's mismatches sit at the named bits.  6 words alone = 383 markers.")
    def M(bits):
        return at(CLEAN, ibs1=bits)
    rows = []
    for name, bits in (("bits 0,1", [0, 1]), ("bits 30,31", [30, 31]),
                       ("bits 62,63", [62, 63]), ("bits 0..19", list(range(20))),
                       ("bits 44..63", list(range(44, 64))), ("all 64", list(range(64)))):
        rows.append((name, Canvas("e4L_" + slug(name), [M(list(range(64))), M(bits)]
                                  + [CLEAN] * 6, nw2=24),
                     Canvas("e4R_" + slug(name), [CLEAN] * 6
                            + [M(bits), M(list(range(64)))], nw2=24)))
    res = many([(c, ("--seglength", "1")) for _, a, b in rows for c in (a, b)])
    print("    %-14s %-28s %s" % ("flanking word", "left of the block", "right of it"))
    for n, (name, a, b) in enumerate(rows):
        ma, mb = mk(a, res[2 * n]), mk(b, res[2 * n + 1])
        print("    %-14s %-28s %s"
              % (name, "%7.1f mk (+%.0f)" % (ma, ma - 383),
                 "%7.1f mk (+%.0f)" % (mb, mb - 383)))
    print("    left  = 127 - (last mismatch bit);  right = 64 + (first mismatch bit)")
    print("    -> the call reaches REACH=%d markers past the nearest mismatch." % REACH)
    blk = []
    for name, body in (("b+2 is a wall", [CLEAN] * 6 + [M([62, 63]), WALL]),
                       ("b+2 has 1 IBS0", [CLEAN] * 6 + [M([62, 63]), at(CLEAN, ibs0=10)]),
                       ("b+1 has 1 IBS0", [CLEAN] * 6 + [at(CLEAN, ibs1=0, ibs0=40),
                                                         M([0, 1])])):
        blk.append((name, Canvas("e4c_" + slug(name), body, nw2=24)))
    res = many([(c, ("--seglength", "1")) for _, c in blk])
    for (name, cv), r in zip(blk, res):
        v = mk(cv, r)
        print("    %-16s %7.1f mk (+%.0f)  <- an IBS0 word blocks the reach whole-word"
              % (name, v, v - 383))


def section5():
    hdr("5. the exhaustive word-sequence battery")
    import itertools
    letters = "Czxy"
    seqs = ["".join(t) for n in (1, 2, 3, 4)
            for t in itertools.product(letters, repeat=n)]
    res = many([(seq_canvas(s), ("--seglength", "1")) for s in seqs])
    bad = []
    for s, r in zip(seqs, res):
        cv = seq_canvas(s)
        if abs(mk(cv, r) - predict_mk(cv)) > 0.3:
            bad.append((s, round(mk(cv, r), 1), round(predict_mk(cv), 1)))
    print("    C = 64 HetHet, z = 64 A2A2, x = 1 mismatch + 63 HetHet, y = 2 + 62")
    print("    predict() reproduces %d of %d" % (len(seqs) - len(bad), len(seqs)))
    old = sum(1 for s, r in zip(seqs, res)
              if abs(mk(seq_canvas(s), r)
                     - predict_bp(seq_canvas(s),
                                  predict_old([wordinfo(w) for w in seq_canvas(s).words]))
                     / seq_canvas(s).s) <= 0.3)
    print("    the committed geometry reproduces %d of %d" % (old, len(seqs)))
    for b in bad:
        print("      miss", b)


def section6():
    hdr("6. the informativeness gate")
    items = []
    for k in (0, 8, 9, 10, 11, 20):
        items.append(("1 word, %d HetHet" % k, Canvas("g6h%d" % k, [w(h=k)])))
    for k in (9, 10, 11, 64):
        items.append(("1 word, %d A1A1/A1A1" % k,
                      Canvas("g6m%d" % k, [{"hom1": k, "zero": 64 - k}])))
    for k in (4, 5):
        items.append(("2 words, %d HetHet each" % k, Canvas("g6t%d" % k, [w(h=k)] * 2)))
    res = many([(c, ("--seglength", "1")) for _, c in items])
    for (lab, cv), r in zip(items, res):
        print("    %-24s %s" % (lab, fmt(cv, r)))
    print("    -> `inf2` (HetHet + A1A1/A1A1) >= %d, and the two are interchangeable."
          % GATE)


def section7():
    hdr("7. the push — every call after the first starts one word late")
    C, Y = CLEAN, w(m=2, h=62)
    items = [("C2 y2 C2", Canvas("p7a", [C, C, Y, Y, C, C])),
             ("C2 y3 C2", Canvas("p7b", [C, C, Y, Y, Y, C, C])),
             ("C2 W2 C2", Canvas("p7c", [C, C, WALL, WALL, C, C])),
             ("C W C6", Canvas("p7d", [C, WALL] + [C] * 6, pad=(0, 3))),
             ("Z0 W C6", Canvas("p7e", [{"zero": 64}, WALL] + [C] * 6, pad=(0, 3))),
             ("C W C6 @6Mb", Canvas("p7f", [C, WALL] + [C] * 6, pad=(0, 3)))]
    extra = [("--seglength", "1")] * 5 + [("--seglength", "6")]
    res = many(list(zip([c for _, c in items], extra)))
    for (lab, cv), r in zip(items, res):
        info = [wordinfo(x) for x in cv.words]
        pn = predict_bp(cv, predict(info, push="never")) / cv.s
        pa = predict_bp(cv, predict(info, push="always")) / cv.s
        print("    %-14s ref %7.1f mk   push=never %6.0f   push=always %6.0f"
              % (lab, mk(cv, r), pn, pa))
    print("    'Z0 W C6' shows the push follows an *emitted* call, not the break. The")
    print("    6 Mb row is the same canvas at --seglength 6, which drops the 1-word first")
    print("    call (5.54 Mb): the reference still reports the second call one word short")
    print("    (319 = 5 words), so the clip runs before the length filter. `predict` does")
    print("    not model --seglength, which is why its column reads 382 there.")


def section8():
    hdr("8. random canvases — the model against the reference")
    for seed, n in ((101, 80), (777, 80), (8081, 80)):
        a, t, bad = battery(seed, n)
        print("    seed %-5d %3d / %3d agree (%.0f %%)" % (seed, a, t, 100 * a / t))
    new, old, tot = compare()
    print("    over all three: 17-rule %d/%d (%.0f %%), committed geometry %d/%d (%.0f %%)"
          % (new, tot, 100 * new / tot, old, tot, 100 * old / tot))


ALPHA = {
    # letter: (composition, (mis, hethet, ibs0, inf2))
    "C": {"hethet": 64},                      # clean, HetHet-rich
    "z": {"zero": 64},                        # clean, worth nothing
    "x": {"ibs1": 1, "hethet": 63},           # one mismatch, informative
    "p": {"ibs1": 1, "zero": 63},             # one mismatch, uninformative
    "y": {"ibs1": 2, "hethet": 62},           # two mismatches
    "W": {"ibs0": 1, "hethet": 63},           # one opposite homozygote
}


def seq_canvas(s, **kw):
    return Canvas("sq_" + s, [ALPHA[c] for c in s], **kw)


def seq_stats(s):
    return [counts(ALPHA[c]) for c in s]


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------

DIRTY_MIS = 2      # het-vs-hom mismatches that make a word unusable for a `.seg` IBD2 run
BRIDGE_MIS = "clean"   # how a lone unusable word is absorbed — see `predict`
GATE = 10          # `inf2` markers, counted from the run's first mismatch-free word
REACH = 63         # markers the call reaches past the flanking word's nearest mismatch


def wordinfo(spec):
    """`(has_ibs0, mis, first_mis_bit, last_mis_bit, inf2)` for one word."""
    ks = expand(spec)
    mis = [i for i, k in enumerate(ks) if k in ("ibs1", "ibs1b")]
    z = any(k in ("ibs0", "ibs0b") for k in ks)
    inf2 = sum(1 for k in ks if k in ("hethet", "hom1", "ibs1b"))
    return (z, len(mis), mis[0] if mis else None, mis[-1] if mis else None, inf2)


def predict(info, w0=0, w1=None, lo=None, hi=None, dirty=DIRTY_MIS,
            bridge=BRIDGE_MIS, gate=GATE, reach=REACH, push="always", clip=0):
    """`.seg` IBD2 calls as marker intervals — the rule of `17-seg-caller.md` §6.

    `info[k]` is `wordinfo` for global word `k`; the usable segment is words `w0..w1`
    covering markers `lo..hi`.
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
    usable = [(not z[k]) and m[k] < dirty for k in range(n)]
    ok = list(usable)
    for k in range(w0, w1 + 1):
        if usable[k] or z[k]:
            continue
        if not (usable[k - 1] if k > w0 else False):
            continue
        if not (usable[k + 1] if k < w1 else False):
            continue
        # A lone unusable word is absorbed only if the run picks up again cleanly: the
        # very next word must carry no mismatch at all, and the usable words after it
        # must carry `gate` informative markers between them.
        if bridge == "clean":
            if m[k + 1] != 0:
                continue
            t, acc = k + 1, 0
            while t <= w1 and usable[t]:
                acc += i2[t]
                t += 1
            if acc < gate:
                continue
        elif m[k] > bridge:
            continue
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
        left = WORD * a
        after_ibs0 = a - 1 >= w0 and z[a - 1]
        if a - 1 >= w0 and not z[a - 1] and lb[a - 1] is not None:
            left = WORD * (a - 1) + lb[a - 1] - reach
            if a - 2 < w0 or z[a - 2]:
                left = max(left, WORD * (a - 1))
        right = WORD * b + WORD - 1
        if b + 1 <= w1 and not z[b + 1] and fb[b + 1] is not None:
            right = WORD * (b + 1) + fb[b + 1] + reach
            if b + 2 > w1 or z[b + 2]:
                right = min(right, WORD * (b + 2) - 1)
        gs = next((t for t in range(a, b + 1) if m[t] == 0), None)
        if gs is None:
            continue
        ge = min(w1, right // WORD) if right > WORD * b + WORD - 1 else b
        if sum(i2[t] for t in range(gs, ge + 1)) < gate:
            continue
        if emitted and (push == "always" or (push == "ibs0" and after_ibs0)):
            left = max(left, WORD * (gs + 1))
        emitted += 1
        # A call touching the usable segment's own first/last complete word runs on to
        # the segment's first/last marker — the fringe the word grid does not cover.
        if a == w0:
            left = min(left, lo)
        if b == w1:
            right = max(right, hi)
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1] + clip)
        if left <= right:
            out.append((left, right))
    return out
def predict_bp(cv, calls=None):
    """The `IBD2Seg` numerator `predict()` implies for a canvas, in base pairs."""
    info = [wordinfo(x) for x in cv.words]
    _, p2 = cv.positions()
    calls = predict(info) if calls is None else calls
    return sum(p2[b] - p2[a] for a, b in calls)


def predict_mk(cv):
    """...in marker intervals (uniform ruler only)."""
    return predict_bp(cv) / cv.s


# ---------------------------------------------------------------------------
# random canvases — the out-of-sample battery
# ---------------------------------------------------------------------------

def random_word(rng):
    """One chr2 word, drawn to exercise every axis the rule reads."""
    ks = ["zero"] * WORD
    slots = list(range(WORD))
    rng.shuffle(slots)
    n = 0
    if rng.random() < 0.30:                      # opposite homozygotes
        for _ in range(rng.choice([1, 1, 2, 5, 64])):
            if n < WORD:
                ks[slots[n]] = "ibs0"
                n += 1
    for _ in range(rng.choice([0, 0, 0, 1, 1, 2, 3, 5, 12, 40, 64])):
        if n < WORD:
            ks[slots[n]] = "ibs1"
            n += 1
    fill = rng.choice(["hethet", "hom1", "zero", "hethet"])
    for _ in range(rng.choice([0, 2, 8, 12, 30, 64])):
        if n < WORD:
            ks[slots[n]] = fill
            n += 1
    return ks


def random_canvases(seed, count, width=10, nw2=16):
    rng = random.Random(seed)
    out = []
    for t in range(count):
        block = [random_word(rng) for _ in range(width)]
        out.append(Canvas("rnd%d_%d" % (seed, t), block, nw2=nw2))
    return out


def battery(seed, count, width=10, verbose=False):
    cvs = random_canvases(seed, count, width)
    res = many([(c, ("--seglength", "1")) for c in cvs])
    agree, bad = 0, []
    for cv, r in zip(cvs, res):
        got, pred = mk(cv, r), predict_mk(cv)
        if abs(got - pred) <= 0.3:
            agree += 1
        else:
            bad.append((cv.name, round(got, 1), round(pred, 1)))
    if verbose:
        for b in bad:
            print("      miss", b)
    return agree, len(cvs), bad


def predict_old(info, w0=0, w1=None, lo=None, hi=None):
    """The **committed** `Scan::ibd2` geometry, for the before/after comparison."""
    n = len(info)
    w1 = n - 1 if w1 is None else w1
    lo = WORD * w0 if lo is None else lo
    hi = WORD * (w1 + 1) - 1 if hi is None else hi
    z = [info[k][0] for k in range(n)]
    m = [info[k][1] for k in range(n)]
    i2 = [info[k][4] for k in range(n)]
    clean = [(not z[k]) and m[k] < 5 for k in range(n)]
    ok = list(clean)
    for k in range(w0 + 1, w1):
        if not clean[k] and clean[k - 1] and clean[k + 1] and not z[k]:
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
    out = []
    for a, b in runs:
        if sum(i2[t] for t in range(a, b + 1)) < GATE:
            continue
        e = w1 if b + 2 >= w1 else b + 1
        left = lo if a == w0 else WORD * a
        right = hi if e == w1 else WORD * e + WORD - 1
        if out:
            left = max(left, out[-1][1] + 1)
        if left <= right:
            out.append((left, right))
    return out


# ---------------------------------------------------------------------------
# §14 — the corrected caller
#
# `predict` above is §7 as it was first fitted, kept byte-for-byte so the numbers §10
# quotes still reproduce.  `predict19` is the same geometry with the one clause §11.1
# flagged replaced by the bisected one: the bridge is the ordinary gate, asked twice.
# ---------------------------------------------------------------------------

def wordinfo2(spec, ibs1b_inf2=False):
    """`wordinfo` with the corrected `inf2` — a het-vs-A1A1 marker is not informative.

    §14.3: 20 words carrying one het-vs-A1A1 marker each are refused by the gate, where
    10 HetHet or 10 A1A1/A1A1 markers pass it.  Pass `ibs1b_inf2=True` for the old
    `p1 & p1` statistic, which §14.5's batteries reject.
    """
    ks = expand(spec)
    mis = [i for i, k in enumerate(ks) if k in ("ibs1", "ibs1b")]
    z = any(k in ("ibs0", "ibs0b") for k in ks)
    kinds = ("hethet", "hom1", "ibs1b") if ibs1b_inf2 else ("hethet", "hom1")
    inf2 = sum(1 for k in ks if k in kinds)
    return (z, len(mis), mis[0] if mis else None, mis[-1] if mis else None, inf2)


def predict19(info, w0=0, w1=None, lo=None, hi=None, dirty=DIRTY_MIS, gate=GATE,
              reach=REACH, bridge="both", push="always", clip=0):
    """The corrected `.seg` IBD2 caller — `17-seg-caller.md` §14.

    Identical to `predict` outside the bridge.  `bridge` selects how a lone unusable
    word is absorbed: `"both"` (the measured rule — the run so far must already pass the
    gate *and* the continuation must pass it from the very next word), `"left"` or
    `"cont"` for either half alone, `"none"` for no absorption at all.
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
    usable = [(not z[k]) and m[k] < dirty for k in range(n)]

    def ge_of(b):
        """The last word a run ending at `b` reaches into — the §5 reach, whole-word."""
        return b + 1 if (b + 1 <= w1 and not z[b + 1] and m[b + 1]) else b

    def gate_ok(g, b):
        return sum(i2[t] for t in range(g, ge_of(b) + 1)) >= gate

    ok = list(usable)
    a = gs = None
    for k in range(w0, w1 + 1):
        if usable[k]:
            if a is None:
                a = k
            if gs is None and m[k] == 0:
                gs = k
            continue
        # An IBS0 word, a run that never found a gate-start, and two unusable words in a
        # row all end the run outright.
        if z[k] or a is None or gs is None or k == w1 or not usable[k + 1]:
            a = gs = None
            continue
        left_ok = bridge in ("both", "left") and gate_ok(gs, k - 1)
        if bridge == "cont":
            left_ok = True
        if not left_ok or m[k + 1] != 0:
            a = gs = None
            continue
        if bridge in ("both", "cont"):
            b2 = k + 1
            while b2 + 1 <= w1 and usable[b2 + 1]:
                b2 += 1
            if not gate_ok(k + 1, b2):
                a = gs = None
                continue
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
        left = WORD * a
        if a - 1 >= w0 and not z[a - 1] and lb[a - 1] is not None:
            left = WORD * (a - 1) + lb[a - 1] - reach
            if a - 2 < w0 or z[a - 2]:
                left = max(left, WORD * (a - 1))
        right = WORD * b + WORD - 1
        if b + 1 <= w1 and not z[b + 1] and fb[b + 1] is not None:
            right = WORD * (b + 1) + fb[b + 1] + reach
            if b + 2 > w1 or z[b + 2]:
                right = min(right, WORD * (b + 2) - 1)
        gs = next((t for t in range(a, b + 1) if m[t] == 0), None)
        if gs is None or not gate_ok(gs, b):
            continue
        if emitted and push == "always":
            left = max(left, WORD * (gs + 1))
        emitted += 1
        if a == w0:
            left = min(left, lo)
        if b == w1:
            right = max(right, hi)
        left, right = max(left, lo), min(right, hi)
        if out:
            left = max(left, out[-1][1] + clip)
        if left <= right:
            out.append((left, right))
    return out


def predict19_mk(cv, ibs1b_inf2=False, **kw):
    """`predict19`'s answer for a canvas, in marker intervals."""
    info = [wordinfo2(x, ibs1b_inf2) for x in cv.words]
    _, p2 = cv.positions()
    return sum(p2[b] - p2[a] for a, b in predict19(info, **kw)) / cv.s


def compare(seeds=(101, 777, 8081), count=80, width=10):
    """`predict` and `predict_old` against the reference on random canvases."""
    tot = new = old = 0
    for seed in seeds:
        cvs = random_canvases(seed, count, width)
        res = many([(c, ("--seglength", "1")) for c in cvs])
        for cv, r in zip(cvs, res):
            got = mk(cv, r)
            info = [wordinfo(x) for x in cv.words]
            tot += 1
            new += abs(got - predict_bp(cv, predict(info)) / cv.s) <= 0.3
            old += abs(got - predict_bp(cv, predict_old(info)) / cv.s) <= 0.3
    return new, old, tot


# ---------------------------------------------------------------------------
# §14 — the campaign that bisected the bridge
# ---------------------------------------------------------------------------

#: The alphabet §14 works in.  `d` and `q` are the two letters §5's alphabet was missing:
#: an unusable word that carries *no* informative markers, and a clean word that carries
#: exactly the gate.  Between them they separate every candidate bridging clause.
ALPHA9 = dict(ALPHA)
ALPHA9.update({
    "d": {"ibs1": 2, "zero": 62},          # unusable and worth nothing
    "q": {"hethet": 10, "zero": 54},       # clean, exactly clears the gate
})


def seq_canvas9(s, **kw):
    return Canvas("sq2_" + s, [ALPHA9[c] for c in s], **kw)


def _b9(items, tag, one=255.0, other=253.0, labels=("BRIDGED", "refused")):
    res = many([(cv, ("--seglength", "1")) for _, cv in items])
    print("    %s" % tag)
    for (lab, cv), r in zip(items, res):
        v = mk(cv, r)
        if one is None:
            which = ""
        else:
            which = labels[0] if abs(v - one) < .5 else (
                labels[1] if abs(v - other) < .5 else "?")
        print("      %-26s %7.1f mk   %-8s  model %6.1f"
              % (lab, v, which, predict19_mk(cv)))


def section9():
    hdr("9. §14 — the bridge, bisected: the gate asked twice")
    print("    A lone unusable word is absorbed iff BOTH halves pass the ordinary gate:")
    print("    the run so far (from its gate-start, through this word), and the")
    print("    continuation (from the very next word, which must be mismatch-free).")

    print("\n  9.1 the left half — [Q(j), d, C, C]; bridged 255, refused 253")
    items = [("hethet=%d" % j,
              Canvas("s9L_h%d" % j, [{"hethet": j, "zero": 64 - j}, ALPHA9["d"],
                                     ALPHA["C"], ALPHA["C"]]))
             for j in (0, 8, 9, 10, 11, 20)]
    _b9(items, "j informative markers before the dirty word:")
    items = [("hom1=%d" % j, Canvas("s9L_m%d" % j,
                                    [{"hom1": j, "zero": 64 - j}, ALPHA9["d"],
                                     ALPHA["C"], ALPHA["C"]])) for j in (9, 10)]
    items += [("hethet=%d+hom1=%d" % (a, b),
               Canvas("s9L_x%d_%d" % (a, b),
                      [{"hethet": a, "hom1": b, "zero": 64 - a - b}, ALPHA9["d"],
                       ALPHA["C"], ALPHA["C"]])) for a, b in ((5, 4), (5, 5), (4, 5))]
    _b9(items, "...the same 9/10 bisection, and the two kinds interchangeable:")

    print("\n  9.2 the right half — [C, y, Q(j)]; bridged 191, refused 127")
    items = [("hethet=%d" % j, Canvas("s9R_h%d" % j,
                                      [ALPHA["C"], ALPHA["y"], {"hethet": j}]))
             for j in (0, 9, 10, 11)]
    _b9(items, "j informative markers after the dirty word:", 191.0, 126.9)

    print("\n  9.3 the window the right half counts over — [C, y, z*r, Q10] and stops")
    items = [("z*%d then q" % r, Canvas("s9W_%d" % r, [ALPHA["C"], ALPHA["y"]]
                                        + [ALPHA["z"]] * r + [ALPHA9["q"]]))
             for r in (0, 1, 2, 3, 4)]
    _b9(items, "distance: the continuation is a whole usable run, however long:", None)
    stops = [("then z, then q", [ALPHA["z"], ALPHA9["q"]]),
             ("then 1 IBS0, then q", [{"ibs0": 1, "zero": 63}, ALPHA9["q"]]),
             ("then 2 mis, then q", [ALPHA9["d"], ALPHA9["q"]]),
             ("then 1 mis, then q", [ALPHA["p"], ALPHA9["q"]]),
             ("then a 2-mis word worth 62", [{"ibs1": 2, "hethet": 62}]),
             ("then a 2-mis word worth 9", [{"ibs1": 2, "hethet": 9, "zero": 53}]),
             ("then an IBS0 word worth 63", [{"ibs0": 1, "hethet": 63}])]
    items = [(lab, Canvas("s9S_" + slug(lab), [ALPHA["C"], ALPHA["y"], ALPHA["z"]] + tail))
             for lab, tail in stops]
    _b9(items, "what ends it: the first unusable word is counted, an IBS0 word is not:",
        None)

    print("\n  9.4 the informativeness correction — one informative marker per word,")
    print("      k words: 10 HetHet or 10 A1A1/A1A1 pass the gate, 20 het-vs-A1A1 do not")
    kinds = [("hethet", {"hethet": 1, "zero": 63}), ("hom1", {"hom1": 1, "zero": 63}),
             ("ibs1b (het/A1A1)", {"ibs1b": 1, "zero": 63}),
             ("ibs1  (het/A2A2)", {"ibs1": 1, "zero": 63})]
    ks = (9, 10, 20)
    cvs = [(lab, k, Canvas("s9I_%s_%d" % (slug(lab), k),
                           [{"zero": 64}] + [spec] * k, nw2=30))
           for lab, spec in kinds for k in ks]
    res = many([(cv, ("--seglength", "1")) for _, _, cv in cvs])
    print("      %-20s %s" % ("one marker per word", "".join("%10s" % ("k=%d" % k)
                                                             for k in ks)))
    row, cur = [], None
    for (lab, k, cv), r in zip(cvs, res):
        if lab != cur:
            if row:
                print("      %-20s %s" % (cur, "".join(row)))
            cur, row = lab, []
        v = mk(cv, r)
        row.append("%10s" % ("called" if abs(v - (64 * (k + 1) - 1)) < .5 else "refused"))
    print("      %-20s %s" % (cur, "".join(row)))

    print("\n  9.5 two unusable words in a row are never absorbed, at any composition")
    def dword(bits, het=0):
        ks = ["zero"] * WORD
        for b in bits:
            ks[b] = "ibs1"
        free = [i for i in range(WORD) if ks[i] == "zero"]
        for i in free[:het]:
            ks[i] = "hethet"
        return ks
    ds = [dword([0, 1], 62), dword([62, 63], 62), dword([0, 1], 0), dword([31, 32], 62),
          dword([0, 1, 2], 61), dword(list(range(8)), 56), dword(list(range(64))),
          dword([0, 1], 10)]
    pairs = [(u, v) for u in range(len(ds)) for v in range(len(ds))]
    cvs = [Canvas("s9T_%d_%d" % (u, v),
                  [ALPHA["C"], ALPHA["C"], ds[u], ds[v], ALPHA["C"], ALPHA["C"]])
           for u, v in pairs]
    res = many([(c, ("--seglength", "1")) for c in cvs])
    merged = sum(1 for c, r in zip(cvs, res) if abs(mk(c, r) - 383.0) < .5)
    okc = sum(1 for c, r in zip(cvs, res) if abs(mk(c, r) - predict19_mk(c)) <= .3)
    print("      %d compositions of [C C d1 d2 C C]: %d merged into one call, "
          "model exact on %d" % (len(cvs), merged, okc))

    print("\n  9.6 the batteries")
    fams = [("length<=4 over {C,z,x,y}",
             ["".join(t) for n in (1, 2, 3, 4)
              for t in itertools.product("Czxy", repeat=n)], seq_canvas),
            ("length-5 over {C,z,x,y}",
             ["".join(t) for t in itertools.product("Czxy", repeat=5)], seq_canvas),
            ("length-4 over {C,z,x,p,y,d,W,q}",
             ["".join(t) for t in itertools.product("CzxpydWq", repeat=4)],
             seq_canvas9)]
    for tag, seqs, mkcv in fams:
        cvs = [mkcv(s) for s in seqs]
        res = many([(c, ("--seglength", "1")) for c in cvs])
        got = [mk(c, r) for c, r in zip(cvs, res)]
        new = sum(1 for c, g in zip(cvs, got) if abs(g - predict19_mk(c)) <= .3)
        nob = sum(1 for c, g in zip(cvs, got)
                  if abs(g - predict19_mk(c, bridge="none")) <= .3)
        b1b = sum(1 for c, g in zip(cvs, got)
                  if abs(g - predict19_mk(c, ibs1b_inf2=True)) <= .3)
        old = sum(1 for c, g in zip(cvs, got)
                  if abs(g - predict_bp(c, predict([wordinfo(x) for x in c.words]))
                         / c.s) <= .3)
        print("    %-34s §14 %4d/%4d   §7 %4d   no bridge %4d   inf2 with ibs1b %4d"
              % (tag, new, len(cvs), old, nob, b1b))
    for seed in (101, 777, 8081):
        cvs = random_canvases(seed, 80, 10)
        res = many([(c, ("--seglength", "1")) for c in cvs])
        got = [mk(c, r) for c, r in zip(cvs, res)]
        new = sum(1 for c, g in zip(cvs, got) if abs(g - predict19_mk(c)) <= .3)
        old = sum(1 for c, g in zip(cvs, got)
                  if abs(g - predict_bp(c, predict([wordinfo(x) for x in c.words]))
                         / c.s) <= .3)
        print("    %-34s §14 %4d/  80   §7 %4d" % ("random seed %d" % seed, new, old))
    for seed in (31337, 424242, 90210):
        cvs = rich_canvases(seed, 100)
        res = many([(c, ("--seglength", "1")) for c in cvs])
        got = [mk(c, r) for c, r in zip(cvs, res)]
        new = sum(1 for c, g in zip(cvs, got) if abs(g - predict19_mk(c)) <= .3)
        old = sum(1 for c, g in zip(cvs, got)
                  if abs(g - predict_bp(c, predict([wordinfo(x) for x in c.words]))
                         / c.s) <= .3)
        b1b = sum(1 for c, g in zip(cvs, got)
                  if abs(g - predict19_mk(c, ibs1b_inf2=True)) <= .3)
        print("    %-34s §14 %4d/ 100   §7 %4d   inf2 with ibs1b %4d"
              % ("rich random seed %d" % seed, new, old, b1b))


    print("\n  9.7 §14.8 — `--seglength` is not a post-filter")
    y = ALPHA["y"]
    cvL = Canvas("s9len", [ALPHA["C"], ALPHA["C"], y, y, ALPHA["C"], ALPHA["C"]],
                 nw1=5, sp1=100_000, nw2=16)     # a 31.9 Mb carrier: the row survives
    steps = [1.0, 4.0, 4.25, 4.5, 8.0, 10.0, 10.25, 15.0, 20.0]
    res = many([(cvL, ("--seglength", "%.6f" % L)) for L in steps])
    print("      [C C y y C C], one word = %.3f Mb" % (WORD * cvL.s / 1e6))
    print("      --seglength Mb %s" % "".join("%8.2f" % L for L in steps))
    print("      IBD2Seg, mk    %s"
          % "".join("%8.0f" % mk(cvL, r) for r in res))
    print("      -> it RISES past one word and falls back past 10 Mb; three re-runs on")
    print("         fresh working directories give the same numbers, so this is not the")
    print("         unseeded QC RNG. A minimum length that only dropped short calls")
    print("         could not do this.")


def rich_word(rng):
    """A random word drawing on *every* marker class the caller can read."""
    ks = ["zero"] * WORD
    slots = list(range(WORD))
    rng.shuffle(slots)
    n = 0

    def put(kind, cnt):
        nonlocal n
        for _ in range(cnt):
            if n < WORD:
                ks[slots[n]] = kind
                n += 1
    if rng.random() < 0.28:
        put("ibs0", rng.choice([1, 1, 2, 3, 9, 64]))
    put("ibs1", rng.choice([0, 0, 0, 1, 1, 1, 2, 2, 3, 6, 20, 64]))
    put("ibs1b", rng.choice([0, 0, 0, 0, 1, 1, 2, 5, 30]))
    put(rng.choice(["hethet", "hom1", "hethet", "zero"]),
        rng.choice([0, 1, 4, 9, 10, 11, 25, 64]))
    put("miss", rng.choice([0, 0, 0, 3, 20]))
    return ks


def rich_canvases(seed, count, width=12, nw2=20):
    rng = random.Random(seed)
    return [Canvas("r2_%d_%d" % (seed, t), [rich_word(rng) for _ in range(width)],
                   nw2=nw2) for t in range(count)]


SECTIONS = {"0": section0, "1": section1, "2": section2, "3": section3,
            "4": section4, "5": section5, "6": section6, "7": section7,
            "8": section8, "9": section9}


def main(argv):
    for k in (argv[1:] or sorted(SECTIONS)):
        SECTIONS[k]()
    save_cache()


if __name__ == "__main__":
    main(sys.argv)


# ---------------------------------------------------------------------------
# the alphabet used by the exhaustive sequence battery (§5)
# ---------------------------------------------------------------------------

