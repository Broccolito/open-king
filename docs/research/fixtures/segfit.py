#!/usr/bin/env python3
"""Fixtures behind `docs/research/16-segment-extension.md` — how an IBD2 run is
started, extended and ended.

The instrument: chr1 is an IBD1 carrier that keeps the pair above `--ibs`'s gates;
chr2 is a canvas painted one *complete word* at a time from an explicit composition
`(mismatch, hethet, ...)`.  Every word of chr2 gets its own marker spacing
(`40000 + 137 w` bp), so the reported `MaxIBD2` length inverts to exactly one word
interval `[u, e]` — the called segment's endpoints, not a total.

Nothing here reads KING's source.

    python3 segfit.py            # all sections
    python3 segfit.py 1 3        # only sections 1 and 3
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixlab as F  # noqa: E402

WORD = 64
ROOT = os.path.dirname(os.path.abspath(__file__))
WORK = os.path.join(ROOT, "work", "segfit")
JOBS = int(os.environ.get("SEGFIT_JOBS", "8"))

# per-marker patterns, as (sample 0, sample 1) A1 dosages; 3 is a missing call
PAIR = {
    "hethet": [1, 1],   # both het          — the IBD2-informative marker
    "zero": [0, 0],     # both A2A2         — clean, uninformative
    "hom1": [2, 2],     # both A1A1         — clean, uninformative to this pass
    "ibs1": [1, 0],     # het vs hom        — the mismatch that breaks a run at five
    "ibs0": [2, 0],     # opposite homs     — irrelevant to this pass
    "miss": [3, 3],
}

WALL = {"ibs1": 64}          # hard-dirty word: breaks any run, worth nothing
CLEAN = {"hethet": 64}       # all-HetHet word


def w(m=0, h=0, **kw):
    """A word carrying `m` mismatches and `h` HetHet; the rest is A2A2."""
    d = {"ibs1": m, "hethet": h}
    d.update(kw)
    return d


def word_positions(nwords, base=40_000, step=137):
    pos, x = [], 1_000_000
    for k in range(nwords):
        for _ in range(WORD):
            pos.append(x)
            x += base + step * k
    return pos


def marker_positions(nmark, base=40_000, step=17):
    """A ruler with a distinct gap after every marker, so a *marker* interval
    inverts uniquely from its length (not just a word interval)."""
    pos, x = [], 1_000_000
    for i in range(nmark):
        pos.append(x)
        x += base + step * i
    return pos


class Probe:
    """One chr2 canvas: `words` is a list of composition dicts, one per complete word."""

    def __init__(self, name, words, nw1=60, nsample=6, seed=3, pad=(2, 4),
                 permarker=False):
        self.name = name
        self.pad = pad
        self.words = [WALL] * pad[0] + list(words) + [WALL] * pad[1]
        self.nw1 = nw1
        self.nsample = nsample
        self.seed = seed
        self.permarker = permarker

    @property
    def nw2(self):
        return len(self.words)

    def build(self, wd):
        F.SPACING = 50_000
        fix = F.Fixture(self.name, [(1, WORD * self.nw1), (2, WORD * self.nw2)],
                        nsample=self.nsample, seed=self.seed, maf=0.5)
        fix.set_state(0, 0, WORD * self.nw1, F.IBD1)
        lo, _ = fix.chrom_span(1)
        pad = [0] * (self.nsample - 2)
        for k, spec in enumerate(self.words):
            if isinstance(spec, dict):
                kinds = []
                for kind, cnt in spec.items():
                    kinds.extend([kind] * cnt)
                kinds = (kinds + ["zero"] * WORD)[:WORD]
            else:                       # explicit 64-marker pattern
                kinds = list(spec)
                assert len(kinds) == WORD, len(kinds)
            for idx, kind in enumerate(kinds):
                fix.pat_all[lo + WORD * k + idx] = PAIR[kind] + pad
                fix.noflip.add(lo + WORD * k + idx)
        prefix = fix.build(wd)
        if self.permarker:
            p1 = [1_000_000 + 50_000 * i for i in range(WORD * self.nw1)]
            p2 = marker_positions(WORD * self.nw2)
        else:
            p1, p2 = word_positions(self.nw1), word_positions(self.nw2)
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

    def decode(self, maxbp):
        """The word interval [u, e] of chr2 whose word-aligned span is `maxbp`."""
        p = word_positions(self.nw2)
        hits = [(u, e) for u in range(self.nw2) for e in range(u, self.nw2)
                if abs(p[WORD * e + WORD - 1] - p[WORD * u] - maxbp) < 1]
        return hits

    def decode_markers(self, maxbp):
        """Every *marker* interval [i, j] of chr2 whose span is `maxbp`, in block
        marker coordinates (0 = first marker of the first non-wall word)."""
        p = marker_positions(WORD * self.nw2)
        off = WORD * self.pad[0]
        want, hits = maxbp, []
        n = len(p)
        j = 0
        for i in range(n):
            while j < n and p[j] - p[i] < want:
                j += 1
            if j < n and abs(p[j] - p[i] - want) < 1:
                hits.append((i - off, j - off))
        return hits


def slug(s):
    return "".join(c if c.isalnum() else "_" for c in str(s))


def run(pr, extra=(), keep=False):
    """(MaxIBD2, Pr_IBD2, stdout, workdir) for the test pair, or None."""
    wd = os.path.join(WORK, slug(pr.name))
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd, exist_ok=True)
    prefix = pr.build(wd)
    cmd = [F.KING, "-b", prefix + ".bed", "--ibs", *extra,
           "--prefix", os.path.join(wd, "k")]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
    out = r.stdout
    got = None
    path = os.path.join(wd, "k.ibs0")
    if "FATAL ERROR" not in out and os.path.exists(path):
        with open(path) as f:
            hdr = f.readline().rstrip("\n").split("\t")
            for line in f:
                d = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if {d.get("ID1"), d.get("ID2")} == {"S00", "S01"} and "MaxIBD2" in d:
                    got = (float(d["MaxIBD2"]), float(d["Pr_IBD2"]))
    if not keep:
        shutil.rmtree(wd, ignore_errors=True)
    return got, out, wd


def call(pr, extra=()):
    """The called word interval, in *block* coordinates (0 = first non-wall word).

    Returns None when nothing is reported, or (u, e) relative to the block start.
    """
    got, _, _ = run(pr, extra)
    if got is None or got[0] <= 0:
        return None
    hits = pr.decode(got[0])
    if len(hits) != 1:
        return ("ambiguous", got[0], hits)
    u, e = hits[0]
    return (u - pr.pad[0], e - pr.pad[0])


def dump(pr, extra=()):
    """Everything the reference says about one canvas, with the decode spelled out."""
    got, out, wd = run(pr, extra, keep=True)
    segs = F.parse_allsegs(os.path.join(wd, "kallsegs.txt"))
    D = sum(float(s["Length"]) for s in segs if int(s["Chr"]) < 23) * 1e6
    shutil.rmtree(wd, ignore_errors=True)
    if got is None:
        return dict(max=None, pr=None, D=D, hits=[], total=0.0)
    p = word_positions(pr.nw2)
    hits = [(u - pr.pad[0], e - pr.pad[0])
            for u in range(pr.nw2) for e in range(u, pr.nw2)
            if abs(p[WORD * e + WORD - 1] - p[WORD * u] - got[0]) < 1]
    return dict(max=got[0], pr=got[1], D=D, hits=hits, total=got[1] * D,
                words=[(u, e, (p[WORD * (e + pr.pad[0]) + WORD - 1]
                               - p[WORD * (u + pr.pad[0])])) for u, e in hits])


def many(probes, extra=()):
    with ThreadPoolExecutor(max_workers=JOBS) as ex:
        return list(ex.map(lambda p: call(p, extra), probes))


# ----------------------------------------------------------------------------
# the model
# ----------------------------------------------------------------------------

DIRTY = 5        # het-vs-hom mismatches that make a word break an IBD2 run
CHUNK_MIS = 5    # mismatches that close a chunk
CHUNK_HET = 95   # HetHet a chunk must carry to be confirmed
MIN_WORDS = 3    # words the measured interval must span
MIN_CHUNK = 3    # words a chunk must span to be confirmable
EXT_MIS = 1      # mismatches the measured interval may reach past the confirmed end


def predict(words, w0=0, w1=None):
    """The rule of `docs/research/16-segment-extension.md`, over per-word `(m, h)`.

    `words[k] = (mismatches, hethet)` for every complete word of one usable segment.
    Returns the reported word intervals `[u, e]`, inclusive.
    """
    n = len(words)
    w1 = n - 1 if w1 is None else w1
    clean = [words[k][0] < DIRTY for k in range(n)]
    ok = list(clean)
    for k in range(1, max(0, n - 1)):
        if not clean[k] and clean[k - 1] and clean[k + 1]:
            ok[k] = True

    def extend(conf, b):
        cum, e = 0, conf
        for k in range(conf + 1, b + 1):
            cum += words[k][0]
            if cum > EXT_MIS:
                break
            e = k
        return e

    raw = []
    k0 = w0
    while k0 <= w1:
        if not ok[k0]:
            k0 += 1
            continue
        a = k0
        while k0 <= w1 and ok[k0]:
            k0 += 1
        b = k0 - 1
        scan_last = min(b + 1, w1)
        exempt = b + 1 >= w1
        u, mis, het, conf, last_mis, cstart = a, 0, 0, None, None, a
        k = a
        while k <= scan_last:
            m, h = words[k]
            before = mis
            mis += m
            het += h
            if mis >= CHUNK_MIS:
                good = het >= CHUNK_HET and k - cstart + 1 >= MIN_CHUNK
                if good or (exempt and k >= scan_last):
                    conf, mis, het, last_mis, cstart = k, 0, 0, None, k + 1
                else:
                    if conf is not None:
                        raw.append((u, extend(conf, b)))
                    if before == CHUNK_MIS - 1 and m == 1 and last_mis is not None:
                        u = last_mis + 1
                    else:
                        u = k + 1
                    mis, het, conf, last_mis, cstart = 0, 0, None, None, u
                    k = u
                    continue
            if m:
                last_mis = k
            k += 1
        if exempt:
            raw.append((u, w1))
        elif conf is not None:
            raw.append((u, w1 if b + 2 >= w1 else extend(conf, b)))

    out = []
    for u, e in raw:
        if out:
            u = max(u, out[-1][1] + 1)
        if u > e or e + 1 - u < MIN_WORDS:
            continue
        out.append((u, e))
    return out


def predict_call(pr):
    """`predict` on a `Probe`, reported the way `call` reports it: the longest interval
    in block coordinates, or None."""
    words = []
    for spec in pr.words:
        if isinstance(spec, dict):
            m = spec.get("ibs1", 0)
            h = spec.get("hethet", 0)
            m, h = min(m, WORD), min(h, max(0, WORD - m))
        else:
            m = sum(1 for x in spec if x == "ibs1")
            h = sum(1 for x in spec if x == "hethet")
        words.append((m, h))
    calls = predict(words)
    if not calls:
        return None
    p = word_positions(pr.nw2)
    best = max(calls, key=lambda c: p[WORD * c[1] + WORD - 1] - p[WORD * c[0]])
    return (best[0] - pr.pad[0], best[1] - pr.pad[0])


# ----------------------------------------------------------------------------
# sections
# ----------------------------------------------------------------------------

def section0():
    print("\n== 0. the grid is not disturbed by monomorphic markers")
    # 8 words, half of every word monomorphic A2A2, the rest HetHet
    pr = Probe("grid_mono", [w(h=32)] * 8)
    got, out, wd = run(pr, keep=True)
    for line in out.splitlines():
        if "words" in line or "SNPs" in line.lower() or "Autosome" in line:
            print("   ", line.strip())
    with open(os.path.join(wd, "kallsegs.txt")) as f:
        for line in f:
            print("   allsegs:", line.rstrip())
    shutil.rmtree(wd, ignore_errors=True)
    print("    MaxIBD2", got, "->", pr.decode(got[0]) if got else None,
          " (block occupies words %d..%d)" % (pr.pad[0], pr.pad[0] + 7))


def section1():
    print("\n== 1. uniform blocks: the reported interval, not just the verdict")
    print(f"{'m':>3} {'h':>3} {'W':>3}  called interval (block coords)")
    for m, h, W in [(0, 12, 8), (0, 11, 8), (0, 13, 8), (1, 19, 8), (1, 18, 8),
                    (1, 20, 8), (1, 24, 8), (1, 32, 8), (1, 64 - 1, 8),
                    (2, 32, 8), (2, 31, 8), (2, 40, 8), (2, 62, 8),
                    (3, 61, 8), (3, 48, 8), (4, 60, 8), (5, 59, 8)]:
        pr = Probe(f"u_{m}_{h}_{W}", [w(m=m, h=h)] * W)
        print(f"{m:>3} {h:>3} {W:>3}  {call(pr)}")


def section2():
    print("\n== 2. calibration: K clean words, and the block-length dependence")
    probes = [Probe(f"cal_{K}", [CLEAN] * K) for K in range(1, 11)]
    for K, got in zip(range(1, 11), many(probes)):
        print(f"  {K:>2} clean words -> {got}")
    print("  -- same (m,h), different block widths --")
    for m, h in ((1, 19), (2, 32), (1, 24), (2, 40), (0, 12)):
        row = []
        for W in (6, 8, 12, 16, 24):
            row.append((W, call(Probe(f"bw_{m}_{h}_{W}", [w(m=m, h=h)] * W))))
        print(f"  m={m} h={h}: " + "  ".join(f"W={x[0]}:{x[1]}" for x in row))


def sweep(m, hs, W=16, tag="s"):
    probes = [Probe(f"{tag}_{m}_{h}_{W}", [w(m=m, h=h)] * W) for h in hs]
    return list(zip(hs, many(probes)))


def section3():
    print("\n== 3. run length against per-word composition (W=16 blocks)")
    for m in (0, 1, 2, 3, 4):
        hs = [h for h in range(4, 64 - m + 1, 1)]
        res = sweep(m, hs)
        # compress: report the interval for each h, collapsing equal runs
        out, prev = [], object()
        for h, v in res:
            if v != prev:
                out.append((h, v))
                prev = v
        print(f"  m={m}: " + "  ".join(f"h>={h}:{v}" for h, v in out))


def section4():
    print("\n== 4. the trim: 20 clean words + j words at (m, h), then a wall")
    P, js = 20, list(range(0, 21))
    print("       j: " + " ".join(f"{j:>2}" for j in js))
    for m, h in [(1, 0), (1, 18), (1, 19), (1, 23), (1, 24), (1, 31), (1, 32),
                 (1, 47), (1, 48), (1, 63), (2, 31), (2, 32), (3, 61), (4, 60)]:
        res = many([Probe(f"S4_{m}_{h}_{j}", [CLEAN] * P + [w(m=m, h=h)] * j)
                    for j in js])
        row = []
        for j, v in zip(js, res):
            e = (v[1] if isinstance(v, tuple) else None) if v else None
            row.append(" -" if e is None else f"{min(e, P + j - 1) - P + 1:>2}")
        print(f"m={m} h={h:>2}: " + " ".join(row))
    print("  (the figure is how many of the j words the reported interval covers;"
          " 5h, 4h, 3h, 2h >= 95 are the h = 19, 24, 32, 48 thresholds)")


def section5():
    print("\n== 5. the start: j words at (m, h) then 20 clean words")
    js = list(range(0, 13))
    print("       j: " + " ".join(f"{j:>2}" for j in js))
    for m, h in [(1, 0), (1, 19), (2, 31), (2, 32), (3, 61), (4, 60), (5, 59)]:
        res = many([Probe(f"S5_{m}_{h}_{j}", [w(m=m, h=h)] * j + [CLEAN] * 20)
                    for j in js])
        print(f"m={m} h={h:>2}: "
              + " ".join((" -" if not v else f"{v[0]:>2}") for v in res))


def section6():
    print("\n== 6. the model against the reference, on random word sequences")
    import random
    random.seed(5)
    probes = []
    for t in range(60):
        ws = []
        for _ in range(random.randint(4, 22)):
            m = random.choice([0, 0, 0, 1, 1, 2, 3, 4, 5, 8, 20, 64])
            ws.append(w(m=m, h=random.randint(0, 64 - min(m, 64))))
        probes.append(Probe(f"S6_{t}", ws))
    got = many(probes)
    ok = bad = 0
    for pr, g in zip(probes, got):
        if isinstance(g, tuple) and len(g) == 3:
            continue
        if g == predict_call(pr):
            ok += 1
        else:
            bad += 1
    print(f"  agree {ok}   disagree {bad}")


SECTIONS = {0: section0, 1: section1, 2: section2, 3: section3,
            4: section4, 5: section5, 6: section6}

if __name__ == "__main__":
    want = [int(a) for a in sys.argv[1:]] or sorted(SECTIONS)
    for n in want:
        SECTIONS[n]()
