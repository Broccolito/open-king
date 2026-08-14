#!/usr/bin/env python3
"""Canvas rig for `--related`'s two-stage screening count (`docs/PARITY.md` §5.7).

The line under investigation is the one stdout line `--related` still gets wrong:

    Stages 1&2 (with 32768 SNPs): <d> pairs of relatives are detected (with kinship > <t>)

`bigish --degree 2` prints **36**; open-king prints **50**.  This file is the instrument
that turns that single aggregate number into per-pair, per-marker measurements.

Nothing here reads KING's source.  Everything is black box: build a fileset, run the
reference, read the count off its stdout.

    python3 screencanvas.py --facts        # re-measure every headline number below
    python3 screencanvas.py --onepair      # which of bigish's candidate pairs pass
    python3 screencanvas.py --boundary     # effective threshold, in true kinship
    python3 screencanvas.py --sweep-n      # threshold vs. sample count
    python3 screencanvas.py --sweep-m      # threshold vs. marker count
    python3 screencanvas.py --informative  # in-sample marker informativeness gate

`KING` repoints the rig at another binary (ours, to compare).

# 1. The two canvases

**Single-pair probe.**  Take `bigish`, drop the six families that carry every cross-family
relative, and keep the remaining 167 samples as a fixed background — on its own the
reference reports `No close relatives are inferred.`, so the count is a clean zero.  Add
back exactly one candidate pair and the printed count is 0 or 1: a direct read of whether
the screen accepts *that pair*.  Summed over `bigish`'s candidates — 47 by the whole-map
estimate, 50 by the map's first 32768 markers, and the three that differ are rejected
either way — this reproduces the reference's **36** exactly, which is what proves the
stage is a **per-pair** rule and not a budget, a cap or a ranking.

**Clone canvas.**  Sample `B'` equals filler `B` outside a chosen marker set `C` and is an
exact clone of filler `A` inside it, so the pair's kinship is tunable continuously by
`|C|` while nothing else about the fileset moves.  Bisecting `|C|` locates the screen's
acceptance boundary to one marker; the exact robust kinship of the constructed pair is
computed here, over any marker subset, so the boundary is read in kinship units.

# 2. What the rig measured (all of it out of sample, on constructed filesets)

* The stage is **per-pair**: 36 = the sum of the single-pair runs.  Sample order does not
  change it (6 permutations); which fillers, only through how many.
* It does **not** use the map's first 32768 markers, which is what this repo implements.
  A clone window placed at markers [40000, 50000) is accepted at the same kinship as one
  at [0, 10000) (`--boundary --windows`), and every stride/offset subset behaves the same.
* The decision is a threshold on the pair's kinship, but the effective threshold sits
  **above** the printed cutoff: 0.0700 against 0.0625 at n=167, m=50000 — lossy in exactly
  the direction that makes 36 < 50.  Reading it as a scaling of the deviation from 0.5,
  `k_screen = 0.5 + R*(k - 0.5)`, gives the same `R` at both degrees (1.0186 at cutoff
  0.1250, 1.0176 at 0.0625); an additive offset disagrees by 18 % between the two and a
  plain multiplicative one by a factor of two, so this is the reading the rig quotes.
* `R` is **1.0000** when `m <= 32768` (measured 0.99995 at 32768, 0.99993 at 33280) and
  grows with the map: 1.0106 at 36864, 1.0128 at 40000, 1.0204 at 45000, 1.0176 at 50000.
* `R` **falls with the sample size** — 1.034 at n=100 to 1.018 at n=167, smoothly, with no
  block structure at 16/32 — and varies pair to pair (1.018…1.026 over six filler pairs)
  but was never below 1.
* The screen's estimate for a pair therefore depends on the **other samples' genotypes**.
  Markers made uninformative in the sample (fillers driven to one homozygote) are
  invisible to it: a pair related *only* on such markers is not detected at kinship 0.154,
  while the same pair related on the informative markers is detected at 0.062.  With the
  whole background replaced by HWE-consistent random genotypes the threshold moves past
  0.25.  No marker subset can produce that — a uniformly spread clone set meets any subset
  in its own proportion — so the estimator itself reads sample-level allele frequencies.

# 3. The lead this leaves, and why it is not landed

A frequency-standardised estimate over a MAF-selected subset — the GRM inner product
`mean_m (x_i - 2p)(x_j - 2p) / (2p(1-p)) / 2` over the 32768 highest-MAF markers — gives
**36** on `bigish` at degree 2, with the whole fileset's 200-sample frequencies, and
agrees with 48 of the 50 single-pair labels.  It is **not** the rule: recomputed with each
single-pair run's own 169-sample frequencies, which is what the reference actually saw in
those runs, it predicts 44 and agrees with 42 of 50, and it gives 16 at degree 1 where the
reference gives 18.  Landing it would be fitting to `bigish`.  It is recorded as a lead,
not a rule.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
DATA = os.path.join(ROOT, "tests", "parity", "work", "data")
BIGISH = os.path.join(DATA, "bigish")
REF = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")

RE_STAGES = re.compile(r"Stages 1&2 \(with (\d+) SNPs\): (\d+) pairs of relatives")

# The six families that hold every cross-family relative in `bigish`.
CORE_FAMS = {"BF01", "BF02", "BF13", "BF14", "BF25", "BF26"}


# --------------------------------------------------------------------------
# PLINK I/O — genotype codes are 0 hom A1, 1 missing, 2 het, 3 hom A2
# --------------------------------------------------------------------------

def load(base=BIGISH):
    fam = open(base + ".fam").read().splitlines()
    bim = open(base + ".bim").read().splitlines()
    n, m = len(fam), len(bim)
    raw = np.fromfile(base + ".bed", dtype=np.uint8)
    assert tuple(raw[:3]) == (0x6C, 0x1B, 0x01), "not a PLINK1 SNP-major .bed"
    body = raw[3:].reshape(m, (n + 3) // 4)
    bits = np.unpackbits(body, axis=1, bitorder="little").reshape(m, -1, 2)
    code = (bits[:, :, 0] | (bits[:, :, 1] << 1))[:, :n]
    return fam, bim, code


def write(out, fam, bim, code):
    m, n = code.shape
    os.makedirs(os.path.dirname(out), exist_ok=True)
    open(out + ".fam", "w").write("\n".join(fam) + "\n")
    open(out + ".bim", "w").write("\n".join(bim) + "\n")
    bits = np.empty((m, n, 2), dtype=np.uint8)
    bits[:, :, 0] = code & 1
    bits[:, :, 1] = (code >> 1) & 1
    flat = bits.reshape(m, n * 2)
    pad = (-flat.shape[1]) % 8
    if pad:
        flat = np.concatenate([flat, np.zeros((m, pad), np.uint8)], axis=1)
    with open(out + ".bed", "wb") as fh:
        fh.write(bytes((0x6C, 0x1B, 0x01)))
        fh.write(np.packbits(flat, axis=1, bitorder="little").tobytes())
    return out + ".bed"


def run(bed, degree=2, binary=None):
    """Run `--related --degree d` and return (screen_snps, detected, stdout)."""
    d = tempfile.mkdtemp(prefix="screencanvas")
    try:
        p = subprocess.run([binary or REF, "-b", bed, "--related", "--degree", str(degree)],
                           cwd=d, capture_output=True, text=True)
        mo = RE_STAGES.search(p.stdout)
        return (int(mo.group(1)) if mo else None,
                int(mo.group(2)) if mo else 0,
                p.stdout)
    finally:
        shutil.rmtree(d, ignore_errors=True)


# --------------------------------------------------------------------------
# Kinship, exactly as `--related`'s between-family stage computes it
# --------------------------------------------------------------------------

def kinship(code, i, j, sub=None):
    a = code[:, i] if sub is None else code[sub, i]
    b = code[:, j] if sub is None else code[sub, j]
    het_i = int(((a == 2) & (b != 1)).sum())
    het_j = int(((b == 2) & (a != 1)).sum())
    hh = int(((a == 2) & (b == 2)).sum())
    ibs0 = int((((a == 0) & (b == 3)) | ((a == 3) & (b == 0))).sum())
    if min(het_i, het_j) == 0:
        return 0.0
    return 0.5 + (2 * hh - 4 * ibs0 - het_i - het_j) / (4 * min(het_i, het_j))


# --------------------------------------------------------------------------
# The clone canvas
# --------------------------------------------------------------------------

class Canvas:
    """Fillers from `bigish` plus one pair whose kinship the cloned marker set sets."""

    def __init__(self, n_fillers=None, markers=None, a="BSNG001", b="BSNG002",
                 work=None):
        fam, bim, code = load()
        fid = [line.split()[0] for line in fam]
        iid = [line.split()[1] for line in fam]
        keep = [i for i in range(len(fam)) if fid[i] not in CORE_FAMS]
        if n_fillers is not None:
            ia, ib = iid.index(a), iid.index(b)
            keep = sorted([i for i in keep if i not in (ia, ib)][:n_fillers - 2] + [ia, ib])
        self.markers = np.arange(code.shape[0]) if markers is None else np.asarray(markers)
        self.fam = [fam[i] for i in keep]
        self.bim = [bim[i] for i in self.markers]
        self.code = code[np.ix_(self.markers, keep)]
        self.m, self.n = self.code.shape
        self.iA = [k for k, i in enumerate(keep) if iid[i] == a][0]
        self.iB = [k for k, i in enumerate(keep) if iid[i] == b][0]
        self.work = work or os.path.join(HERE, "work", "screencanvas")

    def build(self, clone):
        cc = self.code.copy()
        cc[clone, self.iB] = cc[clone, self.iA]
        return cc

    def kinship(self, clone, sub=None):
        return kinship(self.build(clone), self.iA, self.iB, sub)

    def detect(self, clone, degree=2, binary=None, code=None):
        cc = self.build(clone) if code is None else code
        bed = write(self.work, self.fam, self.bim, cc)
        snps, det, out = run(bed, degree, binary)
        return det > 0, out

    def boundary(self, degree=2, seed=1, lo=None, hi=None, binary=None):
        """Bisect the nested clone family; returns (s*, kinship just below/above)."""
        perm = np.random.default_rng(seed).permutation(self.m)
        lo = int(0.03 * self.m) if lo is None else lo
        hi = int(0.45 * self.m) if hi is None else hi
        assert self.detect(np.sort(perm[:hi]), degree, binary)[0], "upper bracket fails"
        assert not self.detect(np.sort(perm[:lo]), degree, binary)[0], "lower bracket fails"
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if self.detect(np.sort(perm[:mid]), degree, binary)[0]:
                hi = mid
            else:
                lo = mid
        return hi, self.kinship(np.sort(perm[:lo])), self.kinship(np.sort(perm[:hi]))


def deflation(k_boundary, degree=2):
    """`R` in `k_screen = 0.5 + R*(k - 0.5)` from a measured boundary kinship."""
    cut = 2.0 ** -(degree + 2)
    return (cut - 0.5) / (k_boundary - 0.5)


# --------------------------------------------------------------------------
# The single-pair probe
# --------------------------------------------------------------------------

def candidate_pairs():
    """`bigish`'s between-family pairs over the degree-2 screen cutoff, best first."""
    fam, bim, code = load()
    fid = [line.split()[0] for line in fam]
    iid = [line.split()[1] for line in fam]
    core = [i for i in range(len(fam)) if fid[i] in CORE_FAMS]
    out = []
    for x, i in enumerate(core):
        for j in core[x + 1:]:
            if fid[i] == fid[j]:
                continue
            k = kinship(code, i, j)
            if k > 0.0625:
                out.append((k, i, j, iid[i], iid[j]))
    return sorted(out, reverse=True), fam, bim, code, fid


def onepair(binary=None):
    pairs, fam, bim, code, fid = candidate_pairs()
    fill = [i for i in range(len(fam)) if fid[i] not in CORE_FAMS]
    work = os.path.join(HERE, "work", "screenone")
    accepted = 0
    print("rank  pair                      kinship  screen")
    for rank, (k, i, j, ni, nj) in enumerate(pairs, 1):
        sel = sorted(fill + [i, j])
        bed = write(work, [fam[s] for s in sel], bim, code[:, sel])
        _, det, _ = run(bed, 2, binary)
        accepted += det > 0
        print("%4d  %-9s %-9s %11.5f  %s" % (rank, ni, nj, k, "accept" if det else "reject"))
    print("\n%d of %d candidate pairs accepted (reference prints 36 for the whole fileset)"
          % (accepted, len(pairs)))
    return accepted


# --------------------------------------------------------------------------
# Sweeps
# --------------------------------------------------------------------------

def sweep_windows(binary=None):
    cv = Canvas()
    print("clone window of 10000 markers slid across the map (prefix model predicts")
    print("detection only while the window overlaps [0, 32768)):")
    for w0 in range(0, cv.m - 10000 + 1, 5000):
        clone = np.arange(w0, w0 + 10000)
        ok, _ = cv.detect(clone, 2, binary)
        print("  W=[%5d,%5d)  k(all)=%.4f  k(first 32768)=%.4f  %s"
              % (w0, w0 + 10000, cv.kinship(clone),
                 cv.kinship(clone, np.arange(32768)), "accept" if ok else "reject"))


def sweep_boundary(binary=None):
    cv = Canvas()
    for degree in (2, 1):
        s, klo, khi = cv.boundary(degree=degree, binary=binary)
        k = (klo + khi) / 2
        print("degree %d (cutoff %.4f): boundary at kinship %.5f  ->  R=%.5f"
              % (degree, 2.0 ** -(degree + 2), k, deflation(k, degree)))


def sweep_n(binary=None):
    for n in (100, 110, 120, 130, 140, 150, 167):
        cv = Canvas(n_fillers=n)
        s, klo, khi = cv.boundary(degree=2, lo=5000, hi=9000, binary=binary)
        k = (klo + khi) / 2
        print("n=%3d  boundary kinship %.5f  R=%.5f" % (n, k, deflation(k)))


def sweep_m(binary=None):
    for m in (32768, 33280, 36864, 40000, 45000, 50000):
        cv = Canvas(markers=np.arange(m))
        s, klo, khi = cv.boundary(degree=2, lo=int(0.04 * m), hi=int(0.40 * m),
                                  binary=binary)
        k = (klo + khi) / 2
        print("m=%5d  boundary kinship %.5f  R=%.5f" % (m, k, deflation(k)))


def sweep_informative(binary=None):
    """Markers the *sample* cannot see are invisible to the screen."""
    cv = Canvas()
    rng = np.random.default_rng(31)
    inf = np.zeros(cv.m, bool)
    inf[rng.choice(cv.m, 32768, replace=False)] = True
    others = [k for k in range(cv.n) if k not in (cv.iA, cv.iB)]
    base = cv.code.copy()
    cv.code = base.copy()
    # Drive the fillers to one homozygote outside `inf`; the pair keeps its genotypes.
    # Use hom A2 (code 3): forcing hom A1 trips the reference's major-allele check.
    cv.code[np.ix_(np.where(~inf)[0], others)] = 3
    for tag, pool in (("informative", np.where(inf)[0]), ("uninformative", np.where(~inf)[0])):
        perm = np.random.default_rng(5).permutation(pool)
        for s in (4000, 6000, 9000, 12000):
            clone = np.sort(perm[:s])
            ok, _ = cv.detect(clone, 2, binary)
            print("  clones in %-13s |C|=%5d  k(all)=%.4f  k(informative)=%.4f  %s"
                  % (tag, s, cv.kinship(clone), cv.kinship(clone, np.where(inf)[0]),
                     "accept" if ok else "reject"))


def facts(binary=None):
    print("== window slide (is the screen the map's first 32768 markers?) ==")
    sweep_windows(binary)
    print("\n== acceptance boundary, both degrees ==")
    sweep_boundary(binary)
    print("\n== boundary vs. sample count ==")
    sweep_n(binary)
    print("\n== boundary vs. marker count ==")
    sweep_m(binary)
    print("\n== in-sample informativeness gate ==")
    sweep_informative(binary)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--facts", action="store_true")
    ap.add_argument("--onepair", action="store_true")
    ap.add_argument("--boundary", action="store_true")
    ap.add_argument("--windows", action="store_true")
    ap.add_argument("--sweep-n", action="store_true")
    ap.add_argument("--sweep-m", action="store_true")
    ap.add_argument("--informative", action="store_true")
    ap.add_argument("--impl", help="binary to drive instead of the reference")
    args = ap.parse_args(argv)
    if not os.path.exists(BIGISH + ".bed"):
        sys.exit("corpus missing: run tests/parity/run_parity.py once to generate it")
    b = args.impl
    if args.onepair:
        onepair(b)
    elif args.boundary:
        sweep_boundary(b)
    elif args.windows:
        sweep_windows(b)
    elif args.sweep_n:
        sweep_n(b)
    elif args.sweep_m:
        sweep_m(b)
    elif args.informative:
        sweep_informative(b)
    else:
        facts(b)


if __name__ == "__main__":
    main()
