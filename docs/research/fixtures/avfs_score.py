#!/usr/bin/env python3
"""Scorer for `--build`'s `INFERENCE AV.FS` statistic — the rig `avfs.py` was missing.

`avfs.py` builds the pedigree shapes that make the reference log `AV.FS` lines.  This
file is the other half: it drives the reference over those shapes, reads the
`Join3/Join2` value out of `kingbuild.log`, recomputes the same statistic from **our**
committed segment engine, and prints the residual next to the amount of it the segment
caller's known over-call can account for.

Nothing here reads KING's source.  Everything is black-box: build a fileset, run the
reference, read its log.

    python3 avfs_score.py                 # the default shape sweep + corpus `bigish`
    python3 avfs_score.py --shapes 4:4 6:6 --seeds 3
    python3 avfs_score.py --pairs         # only the named-pair shape map (§3)

# 1. The statistic

For an ordered triple `(R; N1, N2)` write `IBD(x, y)` for the union of that pair's
called IBD1 and IBD2 segments, as a set of base pairs on the usable-segment map:

    Join2 = | IBD(R,N1) ∩ IBD(R,N2) |
    Join3 = | IBD(R,N1) ∩ IBD(R,N2) ∩ IBD(N1,N2) |

printed at `%.3lf`.  Measuring in SNP counts instead of base pairs, and word-aligning the
intervals instead of using their refined endpoints, were both tried and are worse — see
§2 of the output.

# 2. Why the residual is not this file's fault

Our value is **always** high, by a mean of about +0.004.  The cause is not the formula:

* `Join2` depends only on the two `R`-to-nephew pairs, which are avuncular, so the
  reference reports `IBD2Seg 0.0000` for them — and the reported union `IBD1Seg+IBD2Seg`
  is exact on **all 823** corpus rows whose reference `IBD2Seg` is zero.  Both inputs to
  the denominator are right.
* `Join3` additionally intersects with the *sib* pair, whose reference `IBD2Seg` is not
  zero — and there the union is exact on only **3 of 159** corpus rows, always because
  ours is too big.

So the excess `ΔS` in our sib-pair union can inflate `Join3` by at most `ΔS`, and the
ratio by at most `ΔS / Join2`.  The scorecard prints that bound: every triple measured
so far falls inside `[0, ΔS / Join2]`, i.e. the sib-pair over-call accounts for the whole
residual with nothing left over.  `--build` closes when `docs/PARITY.md` §4.1 closes, and
not before.

# 3. What names the two nephews

The reference prints one `AV.FS` line per `(R, sibship)` and names two of the sibship's
members.  Which two, and in which order, is **a function of the pedigree shape alone**:
`--pairs` re-measures the map.  Holding the shape fixed and changing the genotype seed
never moves it (9 seeds on 4:4, 3 seeds on seven other shapes), and neither does changing
every child's sex (5 sex patterns).  Sliding the whole pedigree down the `.fam` by
prepending singletons moves the named *individuals* but not their positions inside the
sibship, so it is positional, not an absolute sample index.  The positions themselves are
not a fixed permutation — a four-child second family names positions (1,2), (3,4) or
(3,2) depending on how many children the first family has — and the generating rule is
still unidentified.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
_FIT = os.path.join(_ROOT, "tests", "parity", "fit")
sys.path.insert(0, _FIT)

import engine as E  # noqa: E402
import kingdata as kd  # noqa: E402

# The reference. `$KING` repoints it, exactly as `fixlab.py` does.
KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
)
IMPL = os.environ.get("OPEN_KING", os.path.join(_ROOT, "target", "release", "king"))
WORK = os.path.join(_HERE, "work", "avfs")

AV = re.compile(
    r"INFERENCE AV\.FS: (\S+) is [\w, ]+ of (\S+) and (\S+), Join3/Join2=([\d.]+)"
)

# The five `AV.FS` lines the corpus's own `bigish` capture carries.
BIGISH = [
    ("B02_F", "B01_C2", "B01_C3", 0.778),
    ("B01_F", "B02_C3", "B02_C4", 0.801),
    ("B14_F", "B13_C2", "B13_C1", 0.779),
    ("B13_F", "B14_C1", "B14_C2", 0.827),
    ("B25_F", "B26_C3", "B26_C1", 0.803),
]


# ---------------------------------------------------------------------------
# interval algebra over marker indices
# ---------------------------------------------------------------------------

def merge(ivs):
    out = []
    for a, b in sorted(ivs):
        if out and a <= out[-1][1]:
            out[-1] = (out[-1][0], max(out[-1][1], b))
        else:
            out.append((a, b))
    return out


def inter(a, b):
    out = []
    for lo1, hi1 in a:
        for lo2, hi2 in b:
            lo, hi = max(lo1, lo2), min(hi1, hi2)
            if lo < hi:
                out.append((lo, hi))
    return merge(out)


# ---------------------------------------------------------------------------
# datasets
# ---------------------------------------------------------------------------

class Fileset:
    """Any PLINK fileset, in the shape `engine.SegScan` wants.

    The usable-segment list comes from our own binary's `allsegs.txt`, which is
    byte-identical to the reference's in all 163 captures that produce one.
    """

    def __init__(self, bed):
        import numpy as np

        base = bed[:-4]
        self.name = os.path.basename(base)
        chrom, snpname, pos = kd.read_bim(base + ".bim")
        self.fam = kd.read_fam(base + ".fam")
        keep = np.array([c is not None and (1 <= c <= 22 or c == 25) for c in chrom])
        self.pos = pos[keep]
        self.chr = np.array([c for c, k in zip(chrom, keep) if k], dtype=np.int64)
        self.index = {s: i for i, s in enumerate(
            [s for s, k in zip(snpname, keep) if k])}
        self.p0, self.p1 = kd.read_bed_planes(base + ".bed", len(self.fam), keep)
        self.nwords = self.p0.shape[1]
        d = os.path.join(os.path.dirname(os.path.abspath(bed)), "ours")
        os.makedirs(d, exist_ok=True)
        if not os.path.exists(os.path.join(d, "king.seg")):
            subprocess.run([IMPL, "-b", bed, "--ibdseg", "--cpus", "1"], cwd=d,
                           check=True, capture_output=True)
        self.segs = []
        with open(os.path.join(d, "kingallsegs.txt")) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if int(f[1]) != 23:
                    self.segs.append((int(f[1]), self.index[f[6]], self.index[f[7]]))
        self.denom = int(sum(self.pos[hi] - self.pos[lo] for _, lo, hi in self.segs))
        self.ours = read_seg(os.path.join(d, "king.seg"))


def read_seg(path):
    out = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[frozenset((f[1], f[3]))] = (float(f[4]), float(f[5]))
    return out


class Scorer:
    def __init__(self, ds, refseg):
        self.ds = ds
        self.refseg = refseg
        self.idx = {iid: k for k, (_fid, iid) in enumerate(ds.fam)}
        self._calls = {}

    def calls(self, a, b):
        key = frozenset((a, b))
        if key not in self._calls:
            i, j = sorted((self.idx[a], self.idx[b]))
            got = []
            for seg in self.ds.segs:
                sc = E.SegScan(self.ds, i, j, seg, E.BASE)
                if sc.n == 0:
                    continue
                got += list(sc.ibd2(self.ds.pos, E.SEGLEN))
                got += list(sc.ibd1(self.ds.pos, E.SEGLEN))
            self._calls[key] = merge(got)
        return self._calls[key]

    def bp(self, ivs):
        return sum(int(self.ds.pos[b] - self.ds.pos[a]) for a, b in ivs)

    def ratio(self, R, n1, n2):
        j2 = inter(self.calls(R, n1), self.calls(R, n2))
        j3 = inter(j2, self.calls(n1, n2))
        return (self.bp(j3) / self.bp(j2) if self.bp(j2) else 0.0), self.bp(j2)

    def union_error(self, a, b):
        """Ours minus the reference's reported `IBD1Seg + IBD2Seg`, as a fraction of D."""
        key = frozenset((a, b))
        g = self.refseg[key]
        o = self.ds.ours[key]
        return (o[0] + o[1]) - (g[0] + g[1])


# ---------------------------------------------------------------------------
# driving the reference
# ---------------------------------------------------------------------------

def shape_fixture(spec, seed):
    """`spec` is "na:nb" (two families) or "m*k" (m families of k). Returns its dir."""
    tag = spec.replace(":", "x").replace("*", "m")
    name = "av%s_%d" % (tag, seed)
    d = os.path.join(WORK, name)
    bed = os.path.join(d, name + ".bed")
    if not os.path.exists(bed):
        mod = _load("avfs")
        if "*" in spec:
            nfam, kids = (int(x) for x in spec.split("*"))
            mod.multi(name, nfam, kids, max(4, 100 - nfam * (kids + 2)), d, seed)
        else:
            na, nb = (int(x) for x in spec.split(":"))
            mod.two(name, na, nb, max(4, 100 - na - nb - 4), 50000, d, seed)
    if not os.path.exists(os.path.join(d, "kingbuild.log")):
        subprocess.run([KING, "-b", bed, "--build", "--cpus", "1"], cwd=d,
                       check=True, capture_output=True)
    rd = os.path.join(d, "ref")
    os.makedirs(rd, exist_ok=True)
    if not os.path.exists(os.path.join(rd, "king.seg")):
        subprocess.run([KING, "-b", bed, "--ibdseg", "--cpus", "1"], cwd=rd,
                       check=True, capture_output=True)
    return d, bed


def _load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(_HERE, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def lines_of(d):
    seen, out = set(), []
    with open(os.path.join(d, "kingbuild.log")) as fh:
        for R, n1, n2, val in AV.findall(fh.read()):
            if (R, n1, n2) not in seen:
                seen.add((R, n1, n2))
                out.append((R, n1, n2, float(val)))
    return out


# ---------------------------------------------------------------------------
# sections
# ---------------------------------------------------------------------------

def score(shapes, seeds):
    print("Reference `Join3/Join2` against ours, and the residual the sib pair's")
    print("union over-call can account for.  `dU` columns are ours minus the")
    print("reference's reported IBD1Seg+IBD2Seg, as a fraction of D.\n")
    print("%-20s %6s %7s %8s %9s %9s %9s %8s %s"
          % ("case", "ref", "ours", "err", "dU sib", "dU R-N1", "dU R-N2", "J2/D", "bound"))
    rows = []
    ds = kd.load("bigish")
    ds.ours = read_seg(_our_seg("bigish"))
    sc = Scorer(ds, {frozenset((ds.fam[i][1], ds.fam[j][1])): (v[0], v[1])
                     for (i, j), v in ds.ref.items()})
    for R, n1, n2, val in BIGISH:
        rows.append(_row(sc, "bigish", R, n1, n2, val))
    for spec in shapes:
        for k in range(seeds):
            seed = 9000 + 97 * k + sum(ord(c) for c in spec)
            d, bed = shape_fixture(spec, seed)
            got = lines_of(d)
            if not got:
                continue
            ds = Fileset(bed)
            sc = Scorer(ds, read_seg(os.path.join(d, "ref", "king.seg")))
            for R, n1, n2, val in got:
                rows.append(_row(sc, "%s/%d" % (spec, seed), R, n1, n2, val))
    inside = sum(1 for r in rows if 0.0 <= r[0] <= r[1] + 1e-9)
    exact = sum(1 for r in rows if r[2])
    print("\n%d triples: %d inside [0, dU_sib/(J2/D)], %d round to the printed 3 dp"
          % (len(rows), inside, exact))
    print("mean residual %+.4f" % (sum(r[0] for r in rows) / max(1, len(rows))))


def _our_seg(name):
    """Our own `king.seg` for a corpus dataset, run on demand into `work/avfs`."""
    d = os.path.join(WORK, "corpus", name)
    os.makedirs(d, exist_ok=True)
    out = os.path.join(d, "king.seg")
    if not os.path.exists(out):
        bed = os.path.join(_ROOT, "tests", "parity", "work", "data", name + ".bed")
        subprocess.run([IMPL, "-b", bed, "--ibdseg", "--cpus", "1"], cwd=d,
                       check=True, capture_output=True)
    return out


def _row(sc, tag, R, n1, n2, val):
    ours, j2 = sc.ratio(R, n1, n2)
    du_s = sc.union_error(n1, n2)
    du_1 = sc.union_error(R, n1)
    du_2 = sc.union_error(R, n2)
    f2 = j2 / sc.ds.denom
    err = ours - val
    bound = du_s / f2 if f2 else 0.0
    print("%-20s %6.3f %7.4f %+8.4f %+9.4f %+9.4f %+9.4f %8.4f  <=%.4f %s"
          % (tag, val, ours, err, du_s, du_1, du_2, f2, bound,
             "OK" if 0.0 <= err <= bound + 1e-9 else "OUTSIDE"))
    return err, bound, "%.3f" % ours == "%.3f" % val


def pairs(shapes, seeds):
    """§3: which two nephews get named, as a function of shape alone."""
    print("%-14s %-5s %-8s %-26s %s" % ("shape", "seed", "sibship", "members", "named"))
    for spec in shapes:
        for k in range(seeds):
            seed = 9000 + 97 * k + sum(ord(c) for c in spec)
            d, bed = shape_fixture(spec, seed)
            rows = [ln.split() for ln in open(bed[:-4] + ".fam")]
            parents = {r[1]: (r[2], r[3]) for r in rows}
            seen = set()
            for R, n1, n2, _v in lines_of(d):
                key = parents[n1]
                if key in seen:
                    continue
                seen.add(key)
                sibs = [r[1] for r in rows if (r[2], r[3]) == key and r[2] != "0"]
                pos = (sibs.index(n1) + 1, sibs.index(n2) + 1) if n1 in sibs else None
                print("%-14s %-5d %-8s %-26s (%s,%s) pos %s"
                      % (spec, seed, key[0], ",".join(sibs), n1, n2, pos))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--shapes", nargs="*", default=["2:3", "3:3", "3:4", "4:4", "5:5"],
                    help='"na:nb" for two families, "m*k" for m families of k children')
    ap.add_argument("--seeds", type=int, default=2)
    ap.add_argument("--pairs", action="store_true", help="only the named-pair shape map")
    a = ap.parse_args()
    os.makedirs(WORK, exist_ok=True)
    if a.pairs:
        pairs(a.shapes, a.seeds)
    else:
        score(a.shapes, a.seeds)


if __name__ == "__main__":
    main()
