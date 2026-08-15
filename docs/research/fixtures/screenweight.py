#!/usr/bin/env python3
"""What the `--related` screening statistic weights, and what it cannot be.

Second instrument for `docs/PARITY.md` §5.7 --- the one stdout line `--related` still
gets wrong:

    Stages 1&2 (with 32768 SNPs): <d> pairs of relatives are detected (with kinship > <t>)

`screencanvas.py` established *that* the stage is a per-pair threshold whose effective
cut-point sits above the printed one.  This file asks the next question --- **what
functional of the pair is being thresholded** --- and answers it well enough to retire
two whole families of candidate rules, including the `mean(z_i z_j)/2` lead that §5.7
used to record as suggestive.

Nothing here reads KING's source and nothing is fitted to `bigish`'s own answer: every
number below is a bisection against the reference binary on a constructed canvas.

    python3 screenweight.py replicate            # R, both degrees, several m and n
    python3 screenweight.py bands                # boundary confined to one MAF band
    python3 screenweight.py curve --m 50000      # per-marker weight vs. MAF
    python3 screenweight.py curve --m 32768      #   ... and its control
    python3 screenweight.py superpose            # is the rule additive over markers?
    python3 screenweight.py subset               # offline: can any subset fit? (no runs)

# 0. The canvas

`screencanvas.py`'s: `bigish`'s 167 unrelated fillers, which on their own make the
reference print `No close relatives are inferred.`, plus one pair (`BSNG001`,`BSNG002`)
whose kinship is set by the marker set `C` on which B is overwritten with A's genotypes.
Bisecting `|C|` on a nested family locates the acceptance boundary to one marker, which
pins the boundary kinship to about 1.5e-5.  `R` is the deflation of `screencanvas.py`
§2, `k_screen = 0.5 + R*(k - 0.5)`, read off the boundary as `(cut - 0.5)/(k_b - 0.5)`.

# 1. Below 32768 markers the screen is the exact robust kinship

12 bisections, m = 20000 / 30000 / 32768, 4 seeds each:

    R = 1.00010 +- 0.00005,  boundary kinship 0.06252 against a printed cut of 0.0625

So there the statistic is the robust kinship over the whole map and the test is a strict
`>`.  This is the constraint every candidate rule has to survive, and most do not: it
leaves **no room** for a different estimator, a different denominator or a frequency
standardisation, because any of those would still be in force at m <= 32768.

# 2. Above it the deflation is affine in `k - 0.5`, not multiplicative

At m = 50000, n = 167, over independent random clone families:

    degree 2 (cut 0.0625):  R = 1.01779 +- 0.00062   (5 seeds, boundary k 0.0692..0.0708)
    degree 1 (cut 0.1250):  R = 1.01806 +- 0.00053   (4 seeds, boundary k 0.1312..0.1321)

The degree-1 boundary is a **prediction test**, not a second fit.  Carrying degree 2's
`R` across predicts a degree-1 boundary of 0.13161; a multiplicative deflation of the
same size predicts 0.13979.  Measured: 0.1316.  The affine reading is right and the
multiplicative one is dead by ~20 sd --- which also kills every rule of the shape
"numerator over a selected marker set, denominator over the whole map", since those are
multiplicative by construction.

# 3. The rule is additive over markers

Superposition, with a parameter-free prediction.  Pure-band boundaries are 3978 markers
for MAF `[0.35,0.50)` and 8130 for MAF `[0.147,0.249)`; mixing a fraction `x` of the
first should then need `(1-x)` of the second:

    x = 0.50 of band-4's boundary -> band-2 needed 3906, linear prediction 4065  (0.961)
    x = 0.25 of band-4's boundary -> band-2 needed 6001, linear prediction 6098  (0.984)

So the statistic is a weighted marker count, and "what does it weight" is a well posed
question.  (The residual 2-4 % is the draw-to-draw variation in how many of the cloned
markers A is actually heterozygous at, which is what a clone marker contributes.)

# 4. The weight is 2p(1-p) below 32768 markers, and is not above it

Differential probe: hold a base clone set of top-MAF markers just below the boundary,
then bisect how many extra clone markers from a narrow in-sample MAF window are needed
to cross.  The reciprocal of that count is the marker's weight.

**Control, m = 32768** (subsetting inactive).  Ten bands, MAF 0.08 to 0.35:

    w / 2p(1-p) = 1.09 1.12 0.97 1.07 1.07 1.04 0.95 0.91 1.07 1.00      (mean 1.03)

Flat, as it must be: a clone marker moves the robust kinship only when A is heterozygous
there, and `4pq(p^2+q^2) + 8p^2q^2 = 4pq` exactly, so the weight of a marker in *any*
consistent kinship estimator is proportional to `2pq`.

**The measurement, m = 50000** (same markers, same pair, longer map):

    MAF          0.05 0.08 0.11 0.14 | 0.16  0.18  0.20  0.22  0.25  0.30  0.35
    markers needed  -    -    -    - | 2350  1814  1444  1171  1040   765   661
    w / 2p(1-p)     0    0    0    0 | 0.48  0.56  0.65  0.75  0.77  0.94  1.00

The first four bands cannot cross **at any size** --- 2133, 3278, 3772 and 3739 markers
all cloned and still rejected --- so markers below MAF ~0.15 carry no weight at all, and
between there and 0.35 the weight climbs smoothly to `2pq`.

Two families die here:

* **Equal-weight / frequency-standardised estimators**, `mean(z_i z_j)/2` included.  They
  predict the boundary at a fixed *marker count* whatever the band.  Measured counts are
  3978, 4552 and 8130 --- a factor of two.  Het-weighting predicts a fixed het *mass*,
  and the top two bands give 1910 and 1898, agreeing to 0.6 %.  This is the out-of-sample
  refutation of §5.7's recorded lead, and it does not depend on which frequencies the
  lead is computed with.
* **A step-shaped marker subset.**  The probe bands are +-0.015 wide, so a hard MAF
  cut-off would smear over 0.03 of MAF; this ramp spans 0.20.  It is also far broader
  than the ramp any *top-32768-by-X* rule produces: selecting on the in-sample
  heterozygote count over all 167 samples saturates by MAF 0.22 (predicted 0.10 0.45 0.75
  0.92 0.99 1.00 1.00 against the measured row above).  Estimating the selector from a
  16-32 sample subset broadens it about far enough, which is recorded here as the one
  surviving hint and nothing more.

# 5. And it is not the robust kinship over *any* marker subset

The decisive negative, and the reason the placeholder in `related.rs` cannot simply be
swapped for a better subset.  Eight measured brackets --- 5 random clone families, 3
confined to a MAF band --- each say `k_S(size-1) <= 0.0625 < k_S(size)` for the true
screening set `S`.  Scanning `S` over 3 selectors (in-sample MAF, in-sample heterozygote
count, heterozygote count over 24 samples) x 11 sizes from 500 to 50000:

    none of the 33 satisfies all eight.

The reason is structural, and `subset` prints it: a clone family drawn uniformly over the
map meets any subset in its own proportion, so `k_S` tracks `k_all` and the five random
boundaries read **0.069 to 0.072 for every subset tried**, never the cutoff.  A subset can
bend the *band* boundaries around --- that is §4's ramp --- but it cannot supply §2's
uniform deflation.  Whatever the screen computes, it is not this estimator on a subset of
these markers, and no amount of searching for the right subset will find it.

# 6. What is left standing

A statistic that (a) coincides with the robust kinship on the whole map when m <= 32768,
(b) is additive over markers, (c) weights them by `2pq` times a factor that ramps from 0
at MAF ~0.15 to 1 at MAF ~0.35, and (d) is uniformly deflated by an affine `R ~ 1.018`
that varies pair to pair and falls with the sample count.  (a) and (d) together are the
hard part: whatever produces the deflation has to switch off entirely below 32768
markers.  One shape that does all four is the robust kinship over a subset with a
*pair-dependent* denominator --- `R = a_S/a_all` with `a_X` the pair's heterozygosity
asymmetry `(h_i + h_j)/(2*min(h_i,h_j))` over `X` --- which is exactly affine, is exactly
1 when `S` is the whole map, and is pair-specific.  It is written down here as the shape
the constraints imply, **not** as a rule: for `S` = top-32768 by in-sample MAF this pair's
`a_S/a_all` is 0.998, not 1.018, so the subset that would carry it has not been found.
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from screencanvas import CORE_FAMS, kinship, load, run, write  # noqa: E402

WORK = os.path.join(HERE, "work", "screenweight")
CACHE = os.path.join(HERE, "screenweight_measured.json")

# Boundaries measured against the reference by this file; `subset` replays them offline.
MEASURED = {
    "random_d2": [(101, 6659), (102, 6696), (103, 6652), (104, 6765), (105, 6707)],
    "random_d1": [(101, 12797), (102, 13069), (103, 12647), (104, 13088)],
    "bands_d2": [((0.1467, 0.2485), 8130), ((0.2485, 0.3503), 4552),
                 ((0.3503, 0.5001), 3978)],
}


class Canvas:
    """`screencanvas.Canvas` with the in-sample per-marker statistics it needs here."""

    def __init__(self, n_fillers=None, markers=None, a="BSNG001", b="BSNG002"):
        fam, bim, code = load()
        fid = [l.split()[0] for l in fam]
        iid = [l.split()[1] for l in fam]
        keep = [i for i in range(len(fam)) if fid[i] not in CORE_FAMS]
        if n_fillers is not None:
            ia, ib = iid.index(a), iid.index(b)
            keep = sorted([i for i in keep if i not in (ia, ib)][: n_fillers - 2] + [ia, ib])
        self.markers = np.arange(code.shape[0]) if markers is None else np.asarray(markers)
        self.fam = [fam[i] for i in keep]
        self.bim = [bim[i] for i in self.markers]
        self.code = code[np.ix_(self.markers, keep)]
        self.m, self.n = self.code.shape
        self.iA = [k for k, i in enumerate(keep) if iid[i] == a][0]
        self.iB = [k for k, i in enumerate(keep) if iid[i] == b][0]
        self.runs = 0

    @staticmethod
    def freq(code):
        dos = np.select([code == 0, code == 2, code == 3], [0, 1, 2], default=0)
        called = (code != 1).sum(axis=1)
        return dos.sum(axis=1) / (2 * np.maximum(called, 1))

    def maf(self):
        p = self.freq(self.code)
        return np.minimum(p, 1 - p)

    def build(self, clone):
        cc = self.code.copy()
        cc[clone, self.iB] = cc[clone, self.iA]
        return cc

    def stats(self, clone):
        cc = self.build(clone)
        p = self.freq(cc)
        h = 2 * p * (1 - p)
        return kinship(cc, self.iA, self.iB), len(clone), h[clone].sum(), h.sum()

    def detect(self, clone, degree=2, binary=None):
        bed = write(WORK, self.fam, self.bim, self.build(clone))
        self.runs += 1
        return run(bed, degree, binary)[1] > 0

    def bisect(self, pool, degree=2, lo=None, hi=None, binary=None):
        """Smallest accepted prefix of `pool`, or (None, why) if it does not bracket."""
        pool = np.asarray(pool)
        lo = 1 if lo is None else lo
        hi = len(pool) if hi is None else hi
        if not self.detect(np.sort(pool[:hi]), degree, binary):
            return None, "upper bracket rejects"
        if self.detect(np.sort(pool[:lo]), degree, binary):
            return None, "lower bracket accepts"
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.detect(np.sort(pool[:mid]), degree, binary):
                hi = mid
            else:
                lo = mid
        return hi, None


def cmd_replicate(a):
    """§1 and §2: `R` over independent clone families."""
    cv = Canvas(markers=np.arange(a.m) if a.m else None,
                n_fillers=a.n if a.n else None)
    cut = 2.0 ** -(a.degree + 2)
    print(f"# m={cv.m} n={cv.n} degree={a.degree} cutoff={cut:.4f}")
    print("seed  boundary |C|      k_true         R")
    rs = []
    for seed in range(a.seed, a.seed + a.reps):
        pool = np.random.default_rng(seed).permutation(cv.m)
        s, err = cv.bisect(pool, a.degree, int(0.03 * cv.m), int(0.55 * cv.m), a.impl)
        if s is None:
            print(f"{seed:4d}  -- {err}")
            continue
        k = cv.stats(np.sort(pool[:s]))[0]
        rs.append((cut - 0.5) / (k - 0.5))
        print(f"{seed:4d} {s:12d} {k:12.5f} {rs[-1]:9.5f}")
    if rs:
        print(f"mean R {np.mean(rs):.5f}  sd {np.std(rs):.5f}  n={len(rs)}")
        if a.degree == 2:
            r = np.mean(rs)
            print(f"  -> degree-1 boundary predicted affine {0.5 + (0.125-0.5)/r:.5f}, "
                  f"multiplicative {0.125/(0.0625/ (0.5+(0.0625-0.5)/r)):.5f}")
    print(f"[{cv.runs} reference runs]")


def cmd_bands(a):
    """§4: the boundary with the clone set confined to one in-sample MAF quartile."""
    cv = Canvas(markers=np.arange(a.m) if a.m else None)
    maf = cv.maf()
    qs = np.percentile(maf, [0, 25, 50, 75, 100])
    print(f"# m={cv.m} n={cv.n} degree={a.degree}")
    print("band                 avail  boundary |C|    k_true      H(C)   mean het")
    for lo, hi in zip(qs[:-1], qs[1:]):
        idx = np.where((maf >= lo) & (maf < hi + (1e-9 if hi == qs[-1] else 0)))[0]
        pool = np.random.default_rng(a.seed).permutation(idx)
        s, err = cv.bisect(pool, a.degree, binary=a.impl)
        if s is None:
            print(f"MAF[{lo:.4f},{hi:.4f}) {len(idx):6d}   -- {err}")
            continue
        k, c, hc, _ = cv.stats(np.sort(pool[:s]))
        print(f"MAF[{lo:.4f},{hi:.4f}) {len(idx):6d} {c:11d} {k:11.5f} {hc:9.1f} {hc/c:9.4f}")
    print("\n(equal-weight estimators predict one |C| for every band; het-weighted ones\n"
          " predict one H(C).)")
    print(f"[{cv.runs} reference runs]")


def cmd_curve(a):
    """§4: per-marker weight against in-sample MAF, by differential probe."""
    cv = Canvas(markers=np.arange(a.m) if a.m else None)
    maf, rng = cv.maf(), np.random.default_rng(a.seed)
    het = 2 * cv.freq(cv.code) * (1 - cv.freq(cv.code))
    top = rng.permutation(np.where(maf >= a.base)[0])
    s, err = cv.bisect(top, a.degree, binary=a.impl)
    if s is None:
        return print(f"base pool (MAF>={a.base}, {len(top)} markers): {err}")
    back = int(round(a.deficit * 2 * het.sum() / het[top[:s]].mean()))
    base = np.sort(top[: max(s - back, 1)])
    used = set(base.tolist())
    print(f"# m={cv.m} n={cv.n} degree={a.degree}; base pool MAF>={a.base} "
          f"({len(top)} markers, boundary {s}), base {len(base)} markers, "
          f"k={cv.stats(base)[0]:.5f}")
    print(f"{'MAF':>6} {'avail':>6} {'needed':>8} {'2p(1-p)':>9} {'w/2pq':>8}")
    rows = []
    for v in (0.05, 0.08, 0.11, 0.14, 0.16, 0.18, 0.20, 0.22, 0.25, 0.30, 0.35):
        idx = np.where((maf >= v - 0.015) & (maf < v + 0.015) & (maf < a.base))[0]
        idx = np.array([i for i in idx if i not in used])
        if len(idx) < 50:
            continue
        pool = rng.permutation(idx)
        probe = lambda t: np.sort(np.concatenate([base, pool[:t]]))  # noqa: E731
        if not cv.detect(probe(len(pool)), a.degree, a.impl):
            print(f"{v:6.2f} {len(idx):6d} {'none':>8} {2*v*(1-v):9.4f}     0.00"
                  "   (whole band cloned and still rejected)")
            rows.append((v, None))
            continue
        lo, hi = 0, len(pool)
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if cv.detect(probe(mid), a.degree, a.impl):
                hi = mid
            else:
                lo = mid
        rows.append((v, hi))
        print(f"{v:6.2f} {len(idx):6d} {hi:8d} {2*v*(1-v):9.4f}", flush=True)
    ok = [(v, h) for v, h in rows if h]
    if ok:
        v0, h0 = ok[-1]
        print(f"\n{'MAF':>6} {'w/2pq, 1.0 = het-weighted':>28}")
        for v, h in ok:
            print(f"{v:6.2f} {(h0/h)*(2*v0*(1-v0))/(2*v*(1-v)):28.3f}")
    print(f"[{cv.runs} reference runs]")


def cmd_superpose(a):
    """§3: is the acceptance rule additive over markers?"""
    cv = Canvas()
    maf = cv.maf()
    b2 = np.random.default_rng(101).permutation(
        np.where((maf >= 0.1467) & (maf < 0.2485))[0])
    b4 = np.random.default_rng(101).permutation(
        np.where((maf >= 0.3503) & (maf < 0.5001))[0])
    B2, B4 = 8130, 3978
    print(f"# pure-band boundaries: band2={B2} band4={B4}")
    for frac in (0.5, 0.25):
        x, pred = int(round(frac * B4)), int(round((1 - frac) * B2))
        probe = lambda y: np.sort(np.concatenate([b4[:x], b2[:y]]))  # noqa: E731
        lo, hi = 0, min(len(b2), 2 * pred)
        if not cv.detect(probe(hi), a.degree, a.impl):
            print("  upper bracket rejects")
            continue
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if cv.detect(probe(mid), a.degree, a.impl):
                hi = mid
            else:
                lo = mid
        print(f"  band4 x={x:5d}  ->  band2 y measured {hi:5d}, linear prediction "
              f"{pred:5d}, ratio {hi/pred:.3f}")
    print(f"[{cv.runs} reference runs]")


def cmd_subset(a):
    """§5: can the robust kinship over ANY marker subset fit every bracket? No runs."""
    cv = Canvas()
    maf, code = cv.maf(), cv.code
    brackets = []
    for seed, s in MEASURED["random_d2"]:
        brackets.append((f"rand{seed}", np.random.default_rng(seed).permutation(cv.m), s))
    for (lo, hi), s in MEASURED["bands_d2"]:
        idx = np.where((maf >= lo) & (maf < hi))[0]
        brackets.append((f"band{lo:.2f}", np.random.default_rng(101).permutation(idx), s))
    built = [(cv.build(np.sort(p[: s - 1])), cv.build(np.sort(p[:s])))
             for _, p, s in brackets]
    sels = {"in-sample MAF": maf,
            "obs het, all 167": (code == 2).sum(1).astype(float),
            "obs het, 24 samples": (code[:, :24] == 2).sum(1).astype(float)}
    print(f"{'selector':>21} {'|S|':>6}  " + " ".join(f"{n:>8}" for n, _, _ in brackets))
    hits = []
    for name, st in sels.items():
        order = np.argsort(-st, kind="stable")
        for K in (500, 1000, 2000, 4000, 8000, 16000, 24000, 32768, 40000, 50000):
            sel = np.sort(order[:K])
            v = [kinship(h, cv.iA, cv.iB, sel) for _, h in built]
            ok = all(kinship(l, cv.iA, cv.iB, sel) <= 0.0625 < x
                     for (l, _), x in zip(built, v))
            if K in (500, 8000, 32768, 50000) or ok:
                print(f"{name:>21} {K:6d}  " + " ".join(f"{x:8.4f}" for x in v)
                      + ("  ALL OK" if ok else ""))
            if ok:
                hits.append((name, K))
    print(f"\ncutoff is 0.0625; subsets satisfying all {len(brackets)} brackets: {hits}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["replicate", "bands", "curve", "superpose", "subset"])
    ap.add_argument("--degree", type=int, default=2)
    ap.add_argument("--seed", type=int, default=101)
    ap.add_argument("--reps", type=int, default=5)
    ap.add_argument("--m", type=int, default=0, help="truncate the map to this many markers")
    ap.add_argument("--n", type=int, default=0, help="use this many samples")
    ap.add_argument("--base", type=float, default=0.36, help="curve: base pool MAF floor")
    ap.add_argument("--deficit", type=float, default=0.010,
                    help="curve: how far below the boundary the base sits, in kinship")
    ap.add_argument("--impl", help="binary to drive instead of the reference")
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    {"replicate": cmd_replicate, "bands": cmd_bands, "curve": cmd_curve,
     "superpose": cmd_superpose, "subset": cmd_subset}[a.cmd](a)


if __name__ == "__main__":
    main()
