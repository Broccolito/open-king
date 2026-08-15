#!/usr/bin/env python3
"""Fourth instrument for `--related`'s screening count (`docs/research/22-screen.md`).

`screendeflate.py` closed the "which 32 768 markers?" question and left two candidate
mechanisms standing: **merging** markers into 32 768 slots, and **two stages intersected**.
This file kills the first outright, quantitatively as well as structurally, and replaces
the second with two much harder facts:

  * the screen carries a **second necessary condition** that has nothing to do with the
    budget --- it refuses pairs at m = 32 768, where the doc proves the screen is the exact
    whole-map kinship, when the pair's sharing avoids informative markers;
  * above the budget the screen's statistic **depends on the markers it discarded**, which
    no subset of any kind can do, and it does so *deterministically* --- the per-pair
    labels are a sharp threshold with zero inversions.

    python3 screenfold.py step        # one marker over budget: a ramp, not a step
    python3 screenfold.py arrange     # same multiset, three arrangements, same count
    python3 screenfold.py merge       # merge models scored against an 11-point curve
    python3 screenfold.py gate        # the second condition, at m = 32768
    python3 screenfold.py region      # the whole accept region in a two-stratum plane
    python3 screenfold.py strata      # MAF-stratified dilution of a real bigish pair
    python3 screenfold.py separation  # R against how resolvable the ranking is
    python3 screenfold.py sharp       # per-pair labels: sharp, or scattered?
    python3 screenfold.py keys        # in-sample ascertainment: which ranking key?
    python3 screenfold.py facts       # all of the above

`KING` repoints the rig at another binary.

# The instrument

**The ladder.** `screendeflate.py`'s bisection costs ~17 reference runs per number. A
ladder fileset costs one: build `npair` pairs whose kinships climb through the cutoff in
even steps, plus unrelated fillers, and read the printed count. The count *is* the
effective threshold, to the ladder's step (~0.0015 in kinship, ~0.004 in `R`). Everything
below that scans a parameter uses it.

The trap of `screendeflate.py` still applies: A1 must be coded as the **minor** allele.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import screendeflate as S                                        # noqa: E402

WORK = os.path.join(HERE, "work", "screenfold")
BUDGET = 32768
CUT = 0.0625


# --------------------------------------------------------------------------
# The ladder
# --------------------------------------------------------------------------

def ladder(spec, m, n=200, npair=48, lo=0.115, hi=0.34, seed=7, cseed=3):
    """`npair` pairs with clone fractions spread over [lo, hi], plus fillers."""
    fam, bim, code = S.synth(n, m, spec, seed)
    rng = np.random.default_rng(cseed)
    pairs = []
    for t, f in enumerate(np.linspace(lo, hi, npair)):
        iA, iB = 2 * t, 2 * t + 1
        d = np.sort(rng.permutation(m)[:int(round(f * m))])
        code[d, iB] = code[d, iA]
        pairs.append((iA, iB))
    return fam, bim, code, pairs


def kinships(code, pairs, sub=None):
    return np.array([S.kinship(code, i, j, sub) for i, j in pairs])


def read_threshold(fam, bim, code, pairs, degree=2, binary=None, tag="lad"):
    """One reference run -> (printed count, true count, implied threshold, R)."""
    bed = S.write(os.path.join(WORK, tag), fam, bim, code)
    det = S.run(bed, degree, binary)[1]
    ks = np.sort(kinships(code, pairs))[::-1]
    cut = 2.0 ** -(degree + 2)
    thr = ks[det - 1] if 0 < det <= len(ks) else float("nan")
    return det, int((ks > cut).sum()), thr, (cut - 0.5) / (thr - 0.5)


def twopoint(nhi, mafhi, nlo, maflo, shuffle=1234):
    p = np.concatenate([np.full(nhi, mafhi), np.full(nlo, maflo)])
    return p[np.random.default_rng(shuffle).permutation(len(p))]


def between_pairs(fam):
    fid = [line.split()[0] for line in fam]
    n = len(fam)
    return [(i, j) for i in range(n) for j in range(i + 1, n) if fid[i] != fid[j]]


def count_over(code, pairs, cut, sub=None):
    return int((kinships(code, pairs, sub) > cut).sum())


# --------------------------------------------------------------------------
# 1. one marker over budget
# --------------------------------------------------------------------------

def m_step(binary=None):
    """A block merge (blockSize = ceil(m/32768)) predicts a step at m = 32769."""
    fam, bim, code = S.load()
    pairs = between_pairs(fam)
    print("bigish's first 32768 markers, then the real tail one marker at a time:")
    for t in (0, 1, 2, 4, 16, 64, 256, 1024, 2048, 4096, 8192, 12288, 16384, 17232):
        m = BUDGET + t
        bed = S.write(os.path.join(WORK, "step"), fam, bim[:m], code[:m])
        d2 = S.run(bed, 2, binary)[1]
        true = count_over(code[:m], pairs, CUT)
        print("  m=%6d (+%6d)  printed %2d   true %2d   %s"
              % (m, t, d2, true, "" if d2 == true else "<- loses %d" % (true - d2)),
              flush=True)
    print("  no step at the first overflowing marker: exact through +256, then a ramp.")


# --------------------------------------------------------------------------
# 2. the same multiset, three arrangements
# --------------------------------------------------------------------------

def m_arrange(binary=None):
    """Merging needs a grouping; every grouping is lossless for one of these."""
    fam, bim, code = S.load()
    pairs = between_pairs(fam)
    base, dup = BUDGET, 8192
    rng = np.random.default_rng(5)
    adj = []
    for j in range(base):
        adj.append(j)
        if j < dup:
            adj.append(j)
    arrangements = {
        "self-aligned": np.concatenate([np.arange(base), np.arange(dup)]),
        "shuffled": np.concatenate([np.arange(base), rng.permutation(dup)]),
        "adjacent": np.array(adj),
        "offset block": np.concatenate([np.arange(base), 16384 + np.arange(dup)]),
    }
    print("bigish's first 32768 markers + 8192 duplicates, arranged three ways (plus one")
    print("control with a different multiset).")
    print("  self-aligned: slot j = j mod 32768 would merge a marker with its own copy;")
    print("  adjacent:     consecutive-pair merging would; both are lossless for any")
    print("                idempotent op (or / and / max / take-the-more-informative).")
    for name, idx in arrangements.items():
        c = code[idx]
        bb = ["22 x%d %.6f %d A G" % (k, (2e8 + k * 13) / 1e6, 2e8 + k * 13)
              for k in range(len(idx))]
        bed = S.write(os.path.join(WORK, "arr"), fam, bb, c)
        print("  %-14s m=%d  printed deg2 %2d deg1 %2d   true %2d / %2d"
              % (name, len(idx), S.run(bed, 2, binary)[1], S.run(bed, 1, binary)[1],
                 count_over(c, pairs, CUT), count_over(c, pairs, 0.125)), flush=True)
    print("  identical printed counts, and all of them lossy: no index- or rank-adjacent")
    print("  grouping survives.")


# --------------------------------------------------------------------------
# 3. merge models, scored
# --------------------------------------------------------------------------

def _planes(code):
    return (code == 2) | (code == 0), (code == 0) | (code == 1)


def _decode(b1, b2):
    out = np.full(b1.shape, 3, np.uint8)
    out[b1 & ~b2] = 2
    out[b1 & b2] = 0
    out[~b1 & b2] = 1
    return out


def fold(code, key, grouping="stride", op="or"):
    """Fold m markers into 32768 slots. Sparse encoding: hom-major = 00, so a slot that
    merges junk keeps its partner unchanged --- the only encoding in which merging can be
    free, which §4.1 requires."""
    m, n = code.shape
    if m <= BUDGET:
        return code
    order = np.argsort(-key, kind="stable")
    g = (m + BUDGET - 1) // BUDGET
    blocks = ([order[i * BUDGET:(i + 1) * BUDGET] for i in range(g)] if grouping == "stride"
              else [order[i::g][:BUDGET] for i in range(g)])
    if op == "sum":
        dose = np.where(code == 0, 2, np.where(code == 2, 1, 0)).astype(np.int16)
        acc = np.zeros((BUDGET, n), np.int16)
        for blk in blocks:
            a = np.zeros((BUDGET, n), np.int16)
            a[:len(blk)] = dose[blk]
            acc += a
        acc = np.minimum(acc, 2)
        return np.where(acc == 2, 0, np.where(acc == 1, 2, 3)).astype(np.uint8)
    b1, b2 = _planes(code)
    a1 = a2 = None
    for blk in blocks:
        x1 = np.zeros((BUDGET, n), bool); x1[:len(blk)] = b1[blk]
        x2 = np.zeros((BUDGET, n), bool); x2[:len(blk)] = b2[blk]
        if a1 is None:
            a1, a2 = x1, x2
            continue
        if op == "or":
            a1, a2 = a1 | x1, a2 | x2
        elif op == "and":
            pad = np.zeros(BUDGET, bool); pad[:len(blk)] = True
            a1 = np.where(pad[:, None], a1 & x1, a1)
            a2 = np.where(pad[:, None], a2 & x2, a2)
        elif op == "xor":
            a1, a2 = a1 ^ x1, a2 ^ x2
    return _decode(a1, a2)


MERGES = [("stride/or", "stride", "or"), ("stride/and", "stride", "and"),
          ("stride/xor", "stride", "xor"), ("stride/sum", "stride", "sum"),
          ("block/or", "block", "or")]


def m_merge(binary=None):
    """Score every merge model against the reference on the separation scan."""
    print("m=50000 = 32768@MAF 0.45 + 17232@MAF x, n=200, ladder of 48 pairs")
    print("%5s %8s %6s | %s" % ("x", "printed", "true",
                                " ".join("%-10s" % t for t, _, _ in MERGES)))
    for x in (0.05, 0.15, 0.25, 0.35, 0.45):
        fam, bim, code, pairs = ladder(twopoint(BUDGET, 0.45, 17232, x), 50000)
        det = read_threshold(fam, bim, code, pairs, binary=binary, tag="mrg")[0]
        het = (code == 2).sum(1)
        key = (2 * (code == 0).sum(1) + het).astype(float)
        cells = [count_over(fold(code, key, gr, op), pairs, CUT) for _, gr, op in MERGES]
        print("%5.2f %8d %6d | %s" % (x, det, count_over(code, pairs, CUT),
                                      " ".join("%-10d" % c for c in cells)), flush=True)
    print("  or/and inflate (they accept every pair), xor/sum destroy --- and both are")
    print("  already wrong at x = 0.05, where the reference is exact.")


# --------------------------------------------------------------------------
# 4. the second condition
# --------------------------------------------------------------------------

def two_stratum(nA, a, nB, b, n=140, gseed=7):
    m = nA + nB
    p = np.concatenate([np.full(nA, a), np.full(nB, b)])
    lab = np.concatenate([np.zeros(nA, int), np.ones(nB, int)])
    o = np.random.default_rng(1234).permutation(m)
    fam, bim, code = S.synth(n, m, p[o], gseed)
    return fam, bim, code, lab[o]


def m_gate(binary=None):
    """m = 32768, where the screen is the exact whole-map kinship for uniform pairs."""
    print("m=32768 (budget not binding). Pair cloned at every marker of one stratum:")
    for nA, a, nB, b in ((16384, 0.20, 16384, 0.45), (16384, 0.10, 16384, 0.45),
                         (16384, 0.25, 16384, 0.45), (24576, 0.15, 8192, 0.45)):
        fam, bim, code, lab = two_stratum(nA, a, nB, b)
        n = code.shape[1]
        iA, iB = n - 2, n - 1
        c = code.copy()
        d = np.flatnonzero(lab == 0)
        c[d, iB] = c[d, iA]
        bed = S.write(os.path.join(WORK, "gate"), fam, bim, c)
        det = S.run(bed, 2, binary)[1]
        print("  %5d@%.2f cloned + %5d@%.2f untouched:  kinship %.5f   screen: %s"
              % (nA, a, nB, b, S.kinship(c, iA, iB),
                 "%d pair" % det if det else "NO CLOSE RELATIVES"), flush=True)
    print("  KING's own --kinship agrees with each of those numbers; the screen refuses")
    print("  them anyway. No budget is involved: m is exactly 32768.")


def m_region(binary=None):
    """The whole accept region in (fA, fB), so the boundary's shape is visible."""
    fam, bim, code, lab = two_stratum(16384, 0.20, 16384, 0.45)
    n = code.shape[1]
    iA, iB = n - 2, n - 1
    pA = np.random.default_rng(101).permutation(np.flatnonzero(lab == 0))
    pB = np.random.default_rng(102).permutation(np.flatnonzero(lab == 1))
    fs = np.linspace(0, 1, 11)
    print("16384 @ MAF 0.20 (fA cloned) + 16384 @ MAF 0.45 (fB cloned), m=32768")
    print("fA\\fB " + " ".join("%5.2f" % f for f in fs))
    for fA in fs:
        row, ks = [], []
        for fB in fs:
            c = code.copy()
            d = np.sort(pA[:int(round(fA * len(pA)))]); c[d, iB] = c[d, iA]
            d = np.sort(pB[:int(round(fB * len(pB)))]); c[d, iB] = c[d, iA]
            det = S.run(S.write(os.path.join(WORK, "reg"), fam, bim, c), 2, binary)[1]
            row.append("  Y  " if det else "  .  ")
            ks.append(S.kinship(c, iA, iB))
        print("%5.2f %s  | k %s" % (fA, " ".join(row), " ".join("%.3f" % v for v in ks)),
              flush=True)
    print("  the boundary is flat in fA at fB ~ 0.045: a pair at kinship 0.20 is refused")
    print("  for want of sharing among the informative markers.")


def m_strata(binary=None):
    """The same thing on real data: dilute a bigish pair inside one MAF stratum."""
    cand, fill, fam, bim, code = S.candidates()
    for M in (BUDGET, 50000):
        c, b = code[:M], bim[:M]
        order = np.argsort(-S.mafs(c), kind="stable")
        pools = (("everywhere", np.arange(M)),
                 ("top 3/4 by MAF", np.sort(order[:3 * M // 4])),
                 ("top half by MAF", np.sort(order[:M // 2])))
        print("bigish's first %d markers:" % M)
        for k, i, j, ni, nj in cand[:2]:
            for name, pool in pools:
                st, t, kb = _strat_boundary(i, j, fill, fam, b, c, pool, binary=binary)
                print("   %s %s  dilute %-16s -> %s"
                      % (ni, nj, name,
                         "boundary kinship %.5f   R = %.4f" % (kb, S.deflation(kb))
                         if st == "ok" else
                         "%s (kinship %s at full dilution)"
                         % (st, "%.5f" % kb if kb is not None else "-")), flush=True)
    print("  at m=32768 both brackets land on 0.06250 exactly; at m=50000 the deflation")
    print("  more than quadruples when the dilution is aimed at the informative half.")


def _strat_boundary(i, j, fill, fam, bim, code, pool, degree=2, seed=1, binary=None):
    m = code.shape[0]
    sel = sorted([f for f in fill if f not in (i, j)] + [i, j])
    ii, jj = sel.index(i), sel.index(j)
    base = code[:, sel]
    perm = np.random.default_rng(seed).permutation(np.asarray(pool))
    d = np.where(base == 1, 0, np.where(base == 0, 0, np.where(base == 2, 1, 2)))
    pf = np.where(base == 1, 0, d).sum(1) / np.maximum(2 * (base != 1).sum(1), 1)
    rng = np.random.default_rng(seed * 7919 + 13)
    g = ((rng.random(m) < pf).astype(np.int8) + (rng.random(m) < pf).astype(np.int8))
    syn = np.where(g == 0, 0, np.where(g == 1, 2, 3)).astype(base.dtype)

    def build(t):
        c = base.copy()
        d0 = np.sort(perm[:t])
        c[d0, jj] = syn[d0]
        return c

    def ok(t):
        bed = S.write(os.path.join(WORK, "str"), [fam[s] for s in sel], bim, build(t))
        return S.run(bed, degree, binary)[1] > 0

    lo, hi = 0, len(perm)
    if not ok(lo):
        return ("start rejected", None, None)
    if ok(hi - 1):
        return ("never rejected", None, S.kinship(build(hi - 1), ii, jj))
    while hi - lo > 1:
        mid = (lo + hi) // 2
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return ("ok", lo, (S.kinship(build(lo), ii, jj) + S.kinship(build(hi), ii, jj)) / 2)


# --------------------------------------------------------------------------
# 5. resolvability
# --------------------------------------------------------------------------

def keyvals(code, name):
    het = (code == 2).sum(1).astype(np.int64)
    hommin = (code == 0).sum(1).astype(np.int64)
    return {"allele count": (2 * hommin + het).astype(float),
            "hom-minor": hommin.astype(float),
            "het count": het.astype(float)}[name]


def model_count(code, pairs, key, cut=CUT):
    m = code.shape[0]
    sub = (np.arange(m) if m <= BUDGET else
           np.sort(np.argsort(-keyvals(code, key), kind="stable")[:BUDGET]))
    return count_over(code, pairs, cut, sub)


KEYS = ("allele count", "hom-minor", "het count")


def m_separation(binary=None):
    """R against how resolvable the ranking at rank 32768 is, at fixed m and spectrum."""
    print("m=50000 = 32768@MAF 0.45 + 17232@MAF x, n=200.")
    print("%5s %8s %6s %8s %7s %8s | %s"
          % ("x", "printed", "true", "R", "swaps", "kept ==", " ".join("%-13s" % k
                                                                       for k in KEYS)))
    ref = None
    for x in (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.34, 0.38, 0.41, 0.43, 0.45):
        spec = twopoint(BUDGET, 0.45, 17232, x)
        fam, bim, code, pairs = ladder(spec, 50000)
        det, true, thr, R = read_threshold(fam, bim, code, pairs, binary=binary, tag="sep")
        # is the top-32768 by allele count exactly the MAF-0.45 group, and are those
        # markers' genotypes the same map they were at x = 0.05?
        hi = np.flatnonzero(spec > 0.44)
        key = (2 * (code == 0).sum(1) + (code == 2).sum(1)).astype(float)
        kept = np.sort(np.argsort(-key, kind="stable")[:BUDGET])
        swaps = BUDGET - len(np.intersect1d(kept, hi))
        same = "-" if ref is None else str(np.array_equal(code[hi], ref))
        if ref is None:
            ref = code[hi].copy()
        mods = [model_count(code, pairs, k) for k in KEYS]
        print("%5.2f %8d %6d %8.4f %7d %8s | %s"
              % (x, det, true, R, swaps, same, " ".join("%-13d" % v for v in mods)),
              flush=True)
    print("  `swaps` is how many of the kept 32768 come from the low group; `kept ==` is")
    print("  whether those 32768 markers' genotypes are bit-identical to the x = 0.05 map.")
    print("  Through x = 0.30 the kept set is the same index set, holding the same bits,")
    print("  giving every pair the same kinship -- and the printed count still falls from")
    print("  46 to 37. The screen is reading markers it did not keep.")


def m_sharp(binary=None):
    """Per-pair labels where the ranking resolves: systematic shift, or scatter?"""
    x = 0.25
    fam, bim, code, pairs = ladder(twopoint(BUDGET, 0.45, 17232, x), 50000)
    n = code.shape[1]
    key = (2 * (code == 0).sum(1) + (code == 2).sum(1)).astype(float)
    kept = np.sort(np.argsort(-key, kind="stable")[:BUDGET])
    fillers = [i for i in range(n) if i >= 2 * len(pairs)]
    rows = []
    for iA, iB in pairs:
        sel = sorted(fillers + [iA, iB])
        ii, jj = sel.index(iA), sel.index(iB)
        sub = code[:, sel]
        bed = S.write(os.path.join(WORK, "shp"), [fam[s] for s in sel], bim, sub)
        rows.append((S.kinship(sub, ii, jj), S.kinship(sub, ii, jj, kept),
                     S.run(bed, 2, binary)[1] > 0))
    rows.sort(reverse=True)
    print("48 ladder pairs on 32768@0.45 + 17232@0.25, one pair per run, sorted by kinship:")
    print("  " + " ".join("Y" if d else "." for _, _, d in rows))
    acc = [r for r in rows if r[2]]
    rej = [r for r in rows if not r[2]]
    print("  lowest accepted  kinship %.5f (over the kept 32768: %.5f)" % acc[-1][:2])
    print("  highest rejected kinship %.5f (over the kept 32768: %.5f)" % rej[0][:2])
    print("  inversions: %d in the whole-map kinship, %d in the kept-subset kinship"
          % (sum(1 for a in acc for r in rej if r[0] > a[0]),
             sum(1 for a in acc for r in rej if r[1] > a[1])))
    print("  a sharp threshold, displaced: not an intersection of noisy estimates.")


# --------------------------------------------------------------------------
# 6. ascertainment
# --------------------------------------------------------------------------

def m_keys(binary=None):
    """Which ranking key? Only keys that lean against the pair's own heterozygosity can
    deflate at all --- but none of them is the reference."""
    fam, bim, code = S.load()
    pairs = between_pairs(fam)
    print("bigish, m=50000, n=200; the reference prints deg2 = 36, deg1 = 18")
    het = (code == 2).sum(1).astype(np.int64)
    hommin = (code == 0).sum(1).astype(np.int64)
    cands = {"whole map": None,
             "allele count": (2 * hommin + het).astype(float),
             "het count": het.astype(float),
             "hom-minor count": hommin.astype(float),
             "carriers": (hommin + het).astype(float)}
    for name, kv in cands.items():
        sub = (np.arange(code.shape[0]) if kv is None
               else np.sort(np.argsort(-kv, kind="stable")[:BUDGET]))
        print("  top-32768 by %-18s deg2 = %2d  deg1 = %2d"
              % (name, count_over(code, pairs, CUT, sub),
                 count_over(code, pairs, 0.125, sub)), flush=True)
    print()
    print("and against the sample count, on a flat map where ascertainment is the only")
    print("channel a subset has (every marker is exchangeable):")
    for n in (110, 200, 400, 700):
        fam, bim, code, pairs = ladder(0.25, 50000, n=n)
        det, true, thr, R = read_threshold(fam, bim, code, pairs, binary=binary, tag="key")
        print("  n=%4d  printed %2d  true %2d  R = %.4f   (R-1)*n = %5.2f"
              % (n, det, true, R, (R - 1) * n), flush=True)
    print("  (R-1) falls with n but nothing like as fast as 1/n, so the pair's 2-in-n")
    print("  share of the ranking is not the whole effect either.")


MEASUREMENTS = {"step": m_step, "arrange": m_arrange, "merge": m_merge, "gate": m_gate,
                "region": m_region, "strata": m_strata, "separation": m_separation,
                "sharp": m_sharp, "keys": m_keys}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("what", nargs="?", default="facts",
                    choices=list(MEASUREMENTS) + ["facts"])
    ap.add_argument("--impl", help="binary to drive instead of the reference")
    args = ap.parse_args(argv)
    if not os.path.exists(os.path.join(S.DATA, "bigish.bed")):
        sys.exit("corpus missing: run tests/parity/run_parity.py once to generate it")
    os.makedirs(WORK, exist_ok=True)
    for name in (MEASUREMENTS if args.what == "facts" else [args.what]):
        print("== %s ==" % name)
        MEASUREMENTS[name](args.impl)
        print()


if __name__ == "__main__":
    main()
