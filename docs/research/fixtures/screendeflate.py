#!/usr/bin/env python3
"""Third instrument for `--related`'s screening count (`docs/research/22-screen.md`).

The line under investigation is the one stdout line `--related` still gets wrong:

    Stages 1&2 (with 32768 SNPs): <d> pairs of relatives are detected (with kinship > <t>)

`screencanvas.py` established that the stage is a per-pair threshold sitting above the
printed cutoff; `screenweight.py` asked what the statistic weights.  This file answers a
different and sharper question --- **is the screen the kinship over a marker subset at
all** --- and the answer is a clean no, established four independent ways.  It also
measures the law the screen *does* obey to 0.1 %.

Nothing here reads KING's source.  Every number is a bisection or a count read off the
reference binary's stdout, on filesets constructed here.

    python3 screendeflate.py labels      # per-pair accept/reject on `bigish`, both degrees
    python3 screendeflate.py boundary    # the affine law: R at both cutoffs
    python3 screendeflate.py seeds       # is the deflation systematic or realisation noise?
    python3 screendeflate.py order       # does the marker file order matter?
    python3 screendeflate.py subset      # offline: is any marker subset consistent? (no runs)
    python3 screendeflate.py replicate   # r copies of one map: kinship fixed, screen moves
    python3 screendeflate.py append      # appending uninformative markers costs nothing
    python3 screendeflate.py spectrum    # R against the MAF spectrum, synthetic filesets
    python3 screendeflate.py tiegroup    # R against how a tied group is split
    python3 screendeflate.py blocks      # is the loss spread evenly over the map?
    python3 screendeflate.py facts       # all of the above

`KING` repoints the rig at another binary.

# The instrument

**Dilution bisection.**  Take one `bigish` between-family pair in a 169-sample fileset
(the 167 fillers that on their own print `No close relatives are inferred.`, plus the
pair).  Replace one member's genotypes at a growing random marker set with genotypes
drawn from the fileset's own allele frequencies --- a synthetic *unrelated* replacement,
not another sample's, so no new relative pair appears.  The pair's kinship falls
continuously; bisecting the replaced-marker count locates the screen's acceptance
boundary to one marker, i.e. to about 1e-5 in kinship.  The boundary kinship `k*` is
reported as

    R = (cut - 0.5) / (k* - 0.5)          `k_screen = 0.5 + R*(k - 0.5)`

For synthetic filesets the same bisection runs on a cloned pair instead.

# What it measures

1. **The law.**  `0.5 - k_screen = R*(0.5 - k)`, one `R` per fileset.  On `bigish`
   (m=50000, n=169) 36 bisections at each cutoff give `R = 1.02257 +- 0.00065` at 0.0625
   and `1.02079 +- 0.00062` at 0.1250 --- equal to 0.2 %.  On a flat-MAF synthetic map
   where the deflation is four times larger the two agree just as well (1.0798 / 1.0838),
   which is a 25x lever arm on the alternative: a multiplicative `k_screen = c*k` would
   need `c` = 0.659 and 0.812 at the two cutoffs.
2. **`R` is exactly 1 when m <= 32768** --- 0.99999 +- 0.00001 on `bigish`, and 1.00000
   on three synthetic MAF spectra.  Above it `R` grows with `m` and depends strongly on
   the spectrum.
3. **The deflation is systematic, not noise.**  Realisation-to-realisation jitter of the
   boundary is 0.0018 in kinship (8 dilution seeds x 4 pairs); the deflation is 0.0089.
   The per-pair labels agree: every `bigish` pair above kinship 0.0731 is accepted and
   every one below 0.0718 rejected, a sharp threshold, not a 50/50 zone.
4. **It is not the kinship over any marker subset or weighting.**  Algebraically the KING
   estimator is exactly scale-free in `p(1-p)` --- `E[N_l] = 4p_l q_l (1-2phi)` against
   `E[het_l] = 2p_l q_l`, and the `p^2 q^2` terms cancel --- so *every* subset and every
   weighting is unbiased.  Measured: replicating a map r times leaves every kinship
   bit-identical and still moves the count (41 -> 36 at r=3), while `kinship` over
   top-K-by-MAF subsets of the same map stays at 41.
5. **Selecting on an in-sample statistic does not rescue it.**  Simulated ascertainment
   for `top-32768 by in-sample MAF` gives R = 0.995 +- 0.002 (no bias) and for
   `top-32768 by in-sample heterozygote count` R = 0.920 (bias of the *wrong sign*).
6. **The marker file order matters, but only through ties.**  The printed count is 36/18
   under nine permutations of `bigish`'s marker order; the boundary bisection, forty
   times finer, moves by 0.0004 in kinship --- about 5 % of the deflation.
7. **What the deflation tracks.**  It appears exactly when the map holds more
   equally-informative markers than the 32768 budget.  Two-point MAF maps at m=65536:
   with 32768 informative markers (the budget lands between the two groups) R = 1.0081;
   with 40000 R = 1.0596; with 50000 R = 1.0776.  Appending 17232 MAF-0.02 markers to a
   32768-marker map leaves the count at exactly its m=32768 value.
8. **The loss is spread evenly over the map.**  On a flat-MAF map a contiguous clone block
   grown from the head, from the middle, from marker 32768 or backwards from the tail hits
   the boundary at the same kinship to ~2 %.  Whatever the budget costs, it costs every
   marker alike --- there is no favoured region and in particular no favoured prefix.

# A trap worth knowing

KING refuses filesets whose A1 is the major allele --- `FATAL ERROR - Too many first
alleles as the major allele (~100.0%)`.  A synthetic generator must follow PLINK and
code A1 as the *minor* allele, or every run dies and every bisection silently reports
"no bracket".  It also rejects a 50/50 mix (~49.2 %).
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
WORK = os.path.join(HERE, "work", "screendeflate")
REF = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")

RE_STAGES = re.compile(r"Stages 1&2 \(with (\d+) SNPs\): (\d+) pairs of relatives")

# The six families that hold every cross-family relative in `bigish`.
CORE_FAMS = {"BF01", "BF02", "BF13", "BF14", "BF25", "BF26"}

# Rough chromosome lengths in Mb, for synthetic maps that look like a real genome.
CHRLEN = [249, 243, 198, 191, 181, 171, 159, 146, 141, 136, 135,
          134, 115, 107, 102, 90, 81, 78, 59, 63, 48, 51]


# --------------------------------------------------------------------------
# PLINK I/O --- genotype codes are 0 hom A1, 1 missing, 2 het, 3 hom A2
# --------------------------------------------------------------------------

def load(base=BIGISH):
    fam = open(base + ".fam").read().splitlines()
    bim = open(base + ".bim").read().splitlines()
    n, m = len(fam), len(bim)
    raw = np.fromfile(base + ".bed", dtype=np.uint8)
    assert tuple(raw[:3]) == (0x6C, 0x1B, 0x01), "not a PLINK1 SNP-major .bed"
    body = raw[3:].reshape(m, (n + 3) // 4)
    bits = np.unpackbits(body, axis=1, bitorder="little").reshape(m, -1, 2)
    return fam, bim, (bits[:, :, 0] | (bits[:, :, 1] << 1))[:, :n]


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
    """`--related --degree d`; returns (screen SNP count, detected pairs, stdout)."""
    d = tempfile.mkdtemp(prefix="screendeflate")
    try:
        p = subprocess.run([binary or REF, "-b", bed, "--related", "--degree", str(degree)],
                           cwd=d, capture_output=True, text=True)
        mo = RE_STAGES.search(p.stdout)
        return (int(mo.group(1)) if mo else None,
                int(mo.group(2)) if mo else 0, p.stdout)
    finally:
        shutil.rmtree(d, ignore_errors=True)


def kinship(code, i, j, sub=None):
    """KING's reported kinship, verified against `king.kin0` on every map used here."""
    a = code[:, i] if sub is None else code[sub, i]
    b = code[:, j] if sub is None else code[sub, j]
    het_i = int(((a == 2) & (b != 1)).sum())
    het_j = int(((b == 2) & (a != 1)).sum())
    hh = int(((a == 2) & (b == 2)).sum())
    ibs0 = int((((a == 0) & (b == 3)) | ((a == 3) & (b == 0))).sum())
    if min(het_i, het_j) == 0:
        return 0.0
    return 0.5 + (2 * hh - 4 * ibs0 - het_i - het_j) / (4 * min(het_i, het_j))


def mafs(code):
    d = np.where(code == 1, 0, np.where(code == 0, 0, np.where(code == 2, 1, 2)))
    cnt = (code != 1).sum(1)
    p = np.where(code == 1, 0, d).sum(1) / np.maximum(2 * cnt, 1)
    return np.minimum(p, 1 - p)


def deflation(k_boundary, degree=2):
    """`R` in `k_screen = 0.5 + R*(k - 0.5)` from a measured boundary kinship."""
    return (2.0 ** -(degree + 2) - 0.5) / (k_boundary - 0.5)


# --------------------------------------------------------------------------
# Synthetic filesets
# --------------------------------------------------------------------------

def synth_bim(m):
    tot = sum(CHRLEN)
    per = [max(1, round(m * c / tot)) for c in CHRLEN]
    while sum(per) > m:
        per[per.index(max(per))] -= 1
    while sum(per) < m:
        per[per.index(max(per))] += 1
    out = []
    for c, (ln, k) in enumerate(zip(CHRLEN, per), 1):
        step = ln * 1000000 // (k + 1)
        for t in range(1, k + 1):
            out.append("%d rs%d_%d %.6f %d A G" % (c, c, t * step, t * step / 1e6, t * step))
    return out


def synth(n, m, maf, seed=7):
    """`n` unrelated samples over `m` markers. `maf` is a scalar or a length-m array.

    A1 is coded as the **minor** allele, as PLINK writes it; coding it the other way
    round trips KING's major-allele check and every run dies (see the docstring).
    """
    rng = np.random.default_rng(seed)
    p = np.full(m, maf, float) if np.isscalar(maf) else np.asarray(maf, float)[:m]
    g = ((rng.random((m, n)) < p[:, None]).astype(np.int8)
         + (rng.random((m, n)) < p[:, None]).astype(np.int8))
    code = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
    fam = ["F%03d I%03d 0 0 1 0" % (k, k) for k in range(n)]
    return fam, synth_bim(m), code


# --------------------------------------------------------------------------
# `bigish` canvases
# --------------------------------------------------------------------------

def candidates(min_k=0.05):
    """`bigish`'s between-family pairs over `min_k`, best first, plus the fillers."""
    fam, bim, code = load()
    fid = [line.split()[0] for line in fam]
    iid = [line.split()[1] for line in fam]
    core = [i for i in range(len(fam)) if fid[i] in CORE_FAMS]
    fill = [i for i in range(len(fam)) if fid[i] not in CORE_FAMS]
    out = []
    for x, i in enumerate(core):
        for j in core[x + 1:]:
            if fid[i] != fid[j]:
                k = kinship(code, i, j)
                if k > min_k:
                    out.append((k, i, j, iid[i], iid[j]))
    return sorted(out, reverse=True), fill, fam, bim, code


def dilution_boundary(i, j, fill, fam, bim, code, markers=None, degree=2, seed=1,
                      order=None, binary=None):
    """Bisect a dilution of `j` toward an unrelated synthetic sample; boundary kinship."""
    m = code.shape[0] if markers is None else markers
    sel = sorted([f for f in fill if f not in (i, j)] + [i, j])
    ii, jj = sel.index(i), sel.index(j)
    base, bb = code[:m][:, sel], bim[:m]
    perm = np.random.default_rng(seed).permutation(m)
    d = np.where(base == 1, 0, np.where(base == 0, 0, np.where(base == 2, 1, 2)))
    pf = np.where(base == 1, 0, d).sum(1) / np.maximum(2 * (base != 1).sum(1), 1)
    rng = np.random.default_rng(seed * 7919 + 13)
    g = ((rng.random(m) < pf).astype(np.int8) + (rng.random(m) < pf).astype(np.int8))
    synthetic = np.where(g == 0, 0, np.where(g == 1, 2, 3)).astype(base.dtype)

    def build(t):
        c = base.copy()
        c[np.sort(perm[:t]), jj] = synthetic[np.sort(perm[:t])]
        return c

    def ok(t):
        c = build(t)
        bed = write(os.path.join(WORK, "dil"), [fam[s] for s in sel], bb,
                    c if order is None else c[order])
        return run(bed, degree, binary)[1] > 0

    lo, hi = 0, m
    if not ok(lo) or ok(hi - 1):
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            lo = mid
        else:
            hi = mid
    return (kinship(build(lo), ii, jj) + kinship(build(hi), ii, jj)) / 2


def clone_boundary(fam, bim, code, degree=2, seed=1, binary=None):
    """Bisect a clone set on a synthetic fileset; boundary kinship of the last two."""
    m, n = code.shape
    iA, iB = n - 2, n - 1
    perm = np.random.default_rng(seed).permutation(m)

    def build(t):
        c = code.copy()
        c[np.sort(perm[:t]), iB] = c[np.sort(perm[:t]), iA]
        return c

    def ok(t):
        return run(write(os.path.join(WORK, "cln"), fam, bim, build(t)), degree, binary)[1] > 0

    lo, hi = 0, m
    if ok(lo) or not ok(hi - 1):
        return None
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if ok(mid):
            hi = mid
        else:
            lo = mid
    return (kinship(build(lo), iA, iB) + kinship(build(hi), iA, iB)) / 2


# --------------------------------------------------------------------------
# The measurements
# --------------------------------------------------------------------------

def m_labels(binary=None):
    """Per-pair accept/reject, both degrees. The degree-2 column sums to the printed 36."""
    pairs, fill, fam, bim, code = candidates()
    print("rank  pair                    kinship  deg2  deg1")
    a2 = a1 = 0
    for rank, (k, i, j, ni, nj) in enumerate(pairs, 1):
        sel = sorted(fill + [i, j])
        bed = write(os.path.join(WORK, "one"), [fam[s] for s in sel], bim, code[:, sel])
        d2 = run(bed, 2, binary)[1] > 0
        d1 = run(bed, 1, binary)[1] > 0
        a2 += d2
        a1 += d1
        print("%4d  %-9s %-9s %9.5f  %-5s %s" % (rank, ni, nj, k,
              "ACCEPT" if d2 else ".", "ACCEPT" if d1 else "."))
    print("\ndegree 2: %d accepted (reference prints 36 for the whole fileset)" % a2)
    print("degree 1: %d accepted (reference prints 18; the stage is per-pair only to +-1)" % a1)


def m_boundary(binary=None, pairs_n=6, seeds=6):
    """The affine law: `R` at both cutoffs, on `bigish` and on a flat-MAF synthetic map."""
    cand, fill, fam, bim, code = candidates()
    for degree in (2, 1):
        rs = [deflation(b, degree)
              for k, i, j, _, _ in cand[:pairs_n]
              for s in range(1, seeds + 1)
              for b in [dilution_boundary(i, j, fill, fam, bim, code, degree=degree,
                                          seed=s, binary=binary)] if b]
        rs = np.array(rs)
        print("bigish m=50000 n=169 cut %.4f:  R = %.5f +- %.5f  (%d bisections, k* = %.5f)"
              % (2.0 ** -(degree + 2), rs.mean(), rs.std() / np.sqrt(len(rs)), len(rs),
                 0.5 - (0.5 - 2.0 ** -(degree + 2)) / rs.mean()))
    fam2, bim2, code2 = synth(140, 50000, 0.25)
    for degree in (2, 1):
        rs = [deflation(b, degree) for s in (1, 2, 3)
              for b in [clone_boundary(fam2, bim2, code2, degree, s, binary)] if b]
        cut = 2.0 ** -(degree + 2)
        k = 0.5 - (0.5 - cut) / np.mean(rs)
        print("flat-MAF m=50000 n=140 cut %.4f:  R = %.5f +- %.5f   k* = %.5f   cut/k* = %.4f"
              % (cut, np.mean(rs), np.std(rs) / np.sqrt(len(rs)), k, cut / k))
    print("  the two `cut/k*` differ by 23 %; the two `R` agree to 0.4 % -> affine, not"
          " multiplicative")


def m_seeds(binary=None):
    """Systematic or noise? Spread of the boundary over dilution realisations."""
    cand, fill, fam, bim, code = candidates()
    for k, i, j, ni, nj in cand[:4]:
        bs = [b for s in range(1, 9)
              for b in [dilution_boundary(i, j, fill, fam, bim, code, seed=s, binary=binary)]
              if b]
        bs = np.array(bs)
        print("%-9s %-9s  boundary %s  mean %.5f sd %.5f"
              % (ni, nj, " ".join("%.5f" % v for v in bs), bs.mean(), bs.std()))
    print("  deflation 0.0089 against a realisation sd of 0.0018 -> systematic")


def m_order(binary=None):
    """The printed count under marker permutations, then the finer boundary probe."""
    fam, bim, code = load()
    m = code.shape[0]
    mf = mafs(code)
    for tag, order in ([("identity", np.arange(m))]
                       + [("random %d" % s, np.random.default_rng(s).permutation(m))
                          for s in range(1, 7)]
                       + [("MAF desc", np.argsort(-mf, kind="stable")),
                          ("MAF asc", np.argsort(mf, kind="stable"))]):
        bed = write(os.path.join(WORK, "perm"), fam, bim, code[order])
        print("  %-12s deg2 = %2d   deg1 = %2d"
              % (tag, run(bed, 2, binary)[1], run(bed, 1, binary)[1]))
    print("  ... count invariant; now the boundary, forty times finer:")
    cand, fill, fam, bim, code = candidates()
    k, i, j, _, _ = cand[0]
    for tag, order in ([("identity", np.arange(m))]
                       + [("random %d" % s, np.random.default_rng(s).permutation(m))
                          for s in (1, 2, 3)]):
        b = dilution_boundary(i, j, fill, fam, bim, code, seed=1, order=order, binary=binary)
        print("  %-12s boundary %.5f   R = %.5f" % (tag, b, deflation(b)))
    print("  spread 0.0004 in kinship, 5 % of the deflation -> ties in the ranking, only")


def m_subset(binary=None):
    """Offline. Every subset and every weighting of markers is an unbiased estimator."""
    fam, bim, code = load()
    fid = [line.split()[0] for line in fam]
    n = len(fam)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if fid[i] != fid[j]]
    print("count of between-family pairs over the two cutoffs, by marker subset:")
    for m0, Ks in ((50000, (50000, 32768, 25000, 16384)),
                   (16384, (16384, 10922, 8192, 5461))):
        c = code[:m0]
        order = np.argsort(-mafs(c), kind="stable")
        for K in Ks:
            sub = np.sort(order[:K])
            ks = np.array([kinship(c, i, j, sub) for i, j in pairs])
            print("  m=%5d  top-%-6d by MAF:  >0.0625: %2d   >0.1250: %2d"
                  % (m0, K, (ks > 0.0625).sum(), (ks > 0.125).sum()))
    print("  the counts do not move: E[N_l] = 4 p q (1-2phi) and E[het_l] = 2 p q are both")
    print("  proportional to p(1-p), so the ratio --- and the estimate --- is subset-free.")
    print("\nascertainment: does selecting on an in-sample statistic bias it? (simulated)")
    rng = np.random.default_rng(1)
    for name, key in (("in-sample MAF", "maf"), ("in-sample het count", "het")):
        rs = []
        for seed in range(1, 6):
            _, _, c = synth(140, 50000, 0.25, seed)
            iA, iB = 138, 139
            d = rng.permutation(50000)[:int(2 * 0.095 * 50000)]
            c = c.copy()
            c[d, iB] = c[d, iA]
            sel = mafs(c) if key == "maf" else (c == 2).sum(1).astype(float)
            sub = np.sort(np.argsort(-sel, kind="stable")[:32768])
            ka, ks = kinship(c, iA, iB), kinship(c, iA, iB, sub)
            rs.append((0.5 - ks) / (0.5 - ka))
        print("  top-32768 by %-20s R = %.4f +- %.4f"
              % (name, np.mean(rs), np.std(rs) / np.sqrt(len(rs))))
    print("  neither is the 1.02 the reference shows; the het-count one leans the wrong way")


def m_replicate(binary=None):
    """r exact copies of one map: every kinship is bit-identical, the screen still moves."""
    fam, bim, code = load()
    fid = [line.split()[0] for line in fam]
    n = len(fam)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if fid[i] != fid[j]]
    for m0 in (8192, 16384):
        base, bb = code[:m0], bim[:m0]
        ks = np.sort(np.array([kinship(base, i, j) for i, j in pairs]))[::-1]
        print("base m0=%d: %d pairs over 0.0625, %d over 0.1250 (unchanged by replication)"
              % (m0, (ks > 0.0625).sum(), (ks > 0.125).sum()))
        for r in [r for r in range(1, 10) if r * m0 >= 32768][:5]:
            b2 = []
            for q in range(r):
                for line in bb:
                    f = line.split()
                    f[1], f[3] = "c%d_%s" % (q, f[1]), str(int(f[3]) + q)
                    f[2] = "%.6f" % (int(f[3]) / 1e6)
                    b2.append(" ".join(f))
            bed = write(os.path.join(WORK, "rep"), fam, b2, np.concatenate([base] * r, 0))
            d2, d1 = run(bed, 2, binary)[1], run(bed, 1, binary)[1]
            t2 = ks[d2 - 1] if 0 < d2 <= len(ks) else float("nan")
            print("  x%d  m=%6d  deg2 = %2d (implies a threshold of %.5f, R = %.4f)"
                  % (r, r * m0, d2, t2, deflation(t2)))


def m_append(binary=None):
    """Appending uninformative markers to a 32768-marker map costs nothing at all."""
    fam, bim, code = load()
    base = 32768
    bed = write(os.path.join(WORK, "app"), fam, bim[:base], code[:base])
    print("  base 32768 alone            deg2 = %2d  deg1 = %2d"
          % (run(bed, 2, binary)[1], run(bed, 1, binary)[1]))
    k = 50000 - base
    xbim = ["22 rsX_%d %.6f %d A G" % (t, (200000000 + t * 137) / 1e6, 200000000 + t * 137)
            for t in range(k)]
    for maf, seed in ((0.02, 1), (0.05, 2), (0.10, 3), (0.25, 4), (0.50, 5)):
        rng = np.random.default_rng(seed)
        g = ((rng.random((k, code.shape[1])) < maf).astype(np.int8)
             + (rng.random((k, code.shape[1])) < maf).astype(np.int8))
        extra = np.where(g == 2, 0, np.where(g == 1, 2, 3)).astype(np.uint8)
        bed = write(os.path.join(WORK, "app"), fam, bim[:base] + xbim,
                    np.concatenate([code[:base], extra], 0))
        print("  + %d markers at MAF %.2f  deg2 = %2d  deg1 = %2d"
              % (k, maf, run(bed, 2, binary)[1], run(bed, 1, binary)[1]))
    bed = write(os.path.join(WORK, "app"), fam, bim, code)
    print("  + the real bigish tail       deg2 = %2d  deg1 = %2d"
          % (run(bed, 2, binary)[1], run(bed, 1, binary)[1]))


def m_spectrum(binary=None):
    """`R` against the MAF spectrum. Flat spectra deflate hard, skewed ones barely."""
    rng = np.random.default_rng(3)
    specs = (("flat 0.25", 0.25), ("flat 0.45", 0.45),
             ("uniform(0.05,0.5)", rng.uniform(0.05, 0.5, 200000)),
             ("beta(0.6,2.2)", np.clip(rng.beta(0.6, 2.2, 200000), 0.01, 0.5)))
    for name, sp in specs:
        for m in (32768, 40960, 50000, 65536):
            fam, bim, code = synth(140, m, sp)
            rs = [deflation(b) for s in (1, 2, 3)
                  for b in [clone_boundary(fam, bim, code, 2, s, binary)] if b]
            print("  %-18s m=%6d  R = %s" % (name, m,
                  "%.5f +- %.5f" % (np.mean(rs), np.std(rs) / np.sqrt(len(rs)))
                  if rs else "no bracket"))


def m_tiegroup(binary=None):
    """Two-point MAF maps: the deflation switches on when the budget splits a tied group."""
    m, n = 65536, 140
    print("m=65536, MAF is 0.45 on K markers and 0.12 on the rest; the screen keeps 32768")
    for K in (20000, 28000, 32768, 33536, 40000, 50000, 65536):
        rng = np.random.default_rng(11)
        p = np.where(np.arange(m) < K, 0.45, 0.12)[rng.permutation(m)]
        fam, bim, code = synth(n, m, p)
        rs = [deflation(b) for s in (1, 2, 3)
              for b in [clone_boundary(fam, bim, code, 2, s, binary)] if b]
        print("  informative markers K = %6d  R = %s" % (K,
              "%.5f +- %.5f" % (np.mean(rs), np.std(rs) / np.sqrt(len(rs)))
              if rs else "no bracket"))


def m_blocks(binary=None):
    """Is the loss spread evenly? Grow a contiguous clone block from four places."""
    m, n = 50000, 140
    fam, bim, code = synth(n, m, 0.25)
    iA, iB = n - 2, n - 1
    print("flat MAF 0.25, m=50000, n=140; contiguous clone block grown from:")
    for tag, start, step in (("marker 0", 0, +1), ("marker 20000", 20000, +1),
                             ("marker 32768", 32768, +1), ("the tail, backwards", m, -1)):
        def build(t):
            c = code.copy()
            d = (np.arange(start, start + t) if step > 0 else np.arange(start - t, start)) % m
            c[d, iB] = c[d, iA]
            return c

        def ok(t):
            return run(write(os.path.join(WORK, "blk"), fam, bim, build(t)), 2, binary)[1] > 0

        lo, hi = 0, 20000
        if ok(lo) or not ok(hi):
            print("  %-20s no bracket" % tag)
            continue
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if ok(mid):
                hi = mid
            else:
                lo = mid
        b = (kinship(build(lo), iA, iB) + kinship(build(hi), iA, iB)) / 2
        print("  %-20s boundary kinship %.5f   R = %.4f" % (tag, b, deflation(b)))
    print("  the same to ~2 %: no favoured region, and in particular no favoured prefix")


MEASUREMENTS = {"labels": m_labels, "boundary": m_boundary, "seeds": m_seeds,
                "order": m_order, "subset": m_subset, "replicate": m_replicate,
                "append": m_append, "spectrum": m_spectrum, "tiegroup": m_tiegroup,
                "blocks": m_blocks}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("what", nargs="?", default="facts", choices=list(MEASUREMENTS) + ["facts"])
    ap.add_argument("--impl", help="binary to drive instead of the reference")
    args = ap.parse_args(argv)
    if not os.path.exists(BIGISH + ".bed"):
        sys.exit("corpus missing: run tests/parity/run_parity.py once to generate it")
    for name in (MEASUREMENTS if args.what == "facts" else [args.what]):
        print("== %s ==" % name)
        MEASUREMENTS[name](args.impl)
        print()


if __name__ == "__main__":
    main()
