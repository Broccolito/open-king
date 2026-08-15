#!/usr/bin/env python3
"""Fifth instrument for `--related`'s screening count (`docs/research/22-screen.md`).

Round 2 (`screenfold.py`) left one hypothesis standing: a **deterministic bound-based
early exit** over the informativeness-sorted map.  That shape predicts everything Round 2
measured --- it reads markers it does not keep, it is deterministic, it is exact below the
budget, and a pair sharing only at the bottom of the ranking looks unrelated in every
prefix.  This file tests it, and refutes it.

The refutation is a single measurement repeated many ways: **rank is not what matters, an
absolute MAF band is.**  A block of unshared markers placed *above* the pair's sharing in
the ranking costs the pair nothing at all when the block sits at MAF 0.40 --- at any block
size up to 20 480, i.e. five eighths of the map --- and is lethal at MAF 0.45.  A prefix
scan cannot tell those two blocks apart: both are unshared, both are ranked above the
sharing, and they differ by 3 % in informativeness.

    python3 screenexit.py tail     # the screen reads the pair's genotypes at discarded markers
    python3 screenexit.py below    # R < 1 is real: sharing high in the spectrum is a bargain
    python3 screenexit.py band     # the threshold against the unshared block's MAF
    python3 screenexit.py size     # ... and against its size, at MAF 0.45 and 0.40
    python3 screenexit.py rank     # rank position, held apart from MAF: nothing
    python3 screenexit.py coding   # control: recoding A1/A2 changes nothing
    python3 screenexit.py oos      # out-of-sample: fresh seeds, sizes, spectra, both degrees
    python3 screenexit.py facts    # all of the above

`KING` repoints the rig at another binary; `--impl` drives a different one for one run.

# The instruments

**The kept/discarded canvas.**  Build `m` markers as `nA` at MAF 0.45 plus the rest at MAF
`x < 0.42`, so the informativeness ranking separates the two groups with zero swaps, and
clone each ladder pair at fraction `fA` inside the high group and `fB` inside the low one,
drawing the two clone sets from *separate* generators keyed on the rung.  Every arm of a
scan over `fB` then holds the high group's genotypes **bit-identical**, so any movement in
the printed count is the screen reading the markers it discarded.

**The block canvas.**  `u` markers at MAF `a` that the pair never shares, above a window of
`w` markers at a low MAF that it clones at fraction `f`, with junk at MAF 0.03 filling out
to `m`.  A 48-rung ladder over `f` reads the accept threshold in whole-map kinship in one
run.  Quoting the threshold *in kinship* is what makes the scan fair: it already charges
the pair for however much information the block adds to the map.

**Read the threshold, not the label.**  The accept set is the top `det` pairs by realised
kinship, not by ladder index --- near `f = 1` successive rungs differ by less in kinship
than realisation noise.  `implied()` therefore reports the `det`-th largest kinship, as
`screenfold.py` does.  Single-pair labels away from the boundary agree with it.

The `screendeflate.py` trap still applies: A1 must be coded as the **minor** allele.
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import screendeflate as S                                        # noqa: E402

WORK = os.path.join(HERE, "work", "screenexit")
BUDGET = 32768
CUT = 0.0625


# --------------------------------------------------------------------------
# canvases
# --------------------------------------------------------------------------

def twogroup(nhi, mafhi, nlo, maflo, shuffle=1234):
    p = np.concatenate([np.full(nhi, mafhi), np.full(nlo, maflo)])
    o = np.random.default_rng(shuffle).permutation(len(p))
    return p[o], (p[o] > (mafhi + maflo) / 2)


def kept_ladder(spec, hi, fB, lo, hi_f, n=200, npair=48, seed=7):
    """Ladder in fA inside the high group; fB fixed inside the low one.  The high
    group's bits do not depend on fB."""
    m = len(spec)
    fam, bim, code = S.synth(n, m, spec, seed)
    A, B = np.flatnonzero(hi), np.flatnonzero(~hi)
    pairs, fs = [], []
    for t, f in enumerate(np.linspace(lo, hi_f, npair)):
        iA, iB = 2 * t, 2 * t + 1
        d = np.sort(np.random.default_rng(10_000 + t + 97 * seed).permutation(A)
                    [:int(round(f * len(A)))])
        code[d, iB] = code[d, iA]
        if fB > 0:
            d2 = np.sort(np.random.default_rng(20_000 + t + 97 * seed).permutation(B)
                         [:int(round(fB * len(B)))])
            code[d2, iB] = code[d2, iA]
        pairs.append((iA, iB))
        fs.append(f)
    return fam, bim, code, pairs, np.array(fs), A, B


def block_ladder(u, a, w, bmaf, m=32768, junk=0.03, n=300, npair=48, seed=7, shuffle=1234):
    """u unshared markers @ MAF a, above a w-marker window @ bmaf cloned at the rung's f."""
    rest = m - w - u
    p = np.concatenate([np.full(u, a), np.full(w, bmaf), np.full(rest, junk)])
    lab = np.concatenate([np.zeros(u, int), np.ones(w, int), np.full(rest, 2)])
    o = np.random.default_rng(shuffle).permutation(m)
    p, lab = p[o], lab[o]
    fam, bim, code = S.synth(n, m, p, seed)
    pool = np.flatnonzero(lab == 1)
    for t, f in enumerate(np.linspace(0.0, 1.0, npair)):
        d = np.sort(np.random.default_rng(10_000 + t + 97 * seed).permutation(pool)
                    [:int(round(f * len(pool)))])
        code[d, 2 * t + 1] = code[d, 2 * t]
    return fam, bim, code, [(2 * t, 2 * t + 1) for t in range(npair)], lab


def implied(fam, bim, code, pairs, degree=2, binary=None, tag="x"):
    """(det, threshold in whole-map kinship, largest kinship on the ladder)."""
    bed = S.write(os.path.join(WORK, tag), fam, bim, code)
    det = S.run(bed, degree, binary)[1]
    ks = np.sort(np.array([S.kinship(code, i, j) for i, j in pairs]))[::-1]
    if det == 0:
        return det, float("nan"), ks[0]
    if det >= len(ks):
        return det, float("-inf"), ks[0]
    return det, ks[det - 1], ks[0]


def block_threshold(u, a, w=8192, bmaf=0.15, m=32768, n=300, seeds=(7, 11, 23),
                    degree=2, binary=None, tag="blk"):
    out, top = [], []
    for sd in seeds:
        fam, bim, code, pairs, _ = block_ladder(u, a, w, bmaf, m, n=n, seed=sd)
        det, thr, kmax = implied(fam, bim, code, pairs, degree, binary, tag)
        out.append(thr)
        top.append(kmax)
    return np.array(out), float(np.mean(top))


def fmt(v):
    return "  n/a" if np.isnan(v) else ("free" if v == float("-inf") else "%.4f" % v)


# --------------------------------------------------------------------------
# 1. the screen reads markers it discards
# --------------------------------------------------------------------------

def m_tail(binary=None):
    """Round 2 §10 held the kept bits fixed and changed the discarded markers' *frequencies*.
    This holds them fixed and changes the discarded markers' *genotypes* for the pair."""
    m, nA, x = 50000, BUDGET, 0.25
    spec, hi = twogroup(nA, 0.45, m - nA, x)
    print("m=%d = %d @ MAF 0.45 (kept by any informativeness rank) + %d @ MAF %.2f."
          % (m, nA, m - nA, x))
    print("48-rung ladder in the clone fraction fA *inside the kept group*; fB is the clone")
    print("fraction inside the discarded group and is the only thing that changes.")
    print()
    print("  %6s %8s %10s %10s %10s %12s" %
          ("fB", "printed", "thr(whole)", "thr(kept)", "R", "kept bits =="))
    ref = None
    for fB in (0.0, 0.05, 0.10, 0.20, 0.40, 0.70, 1.00):
        fam, bim, code, pairs, fs, A, B = kept_ladder(spec, hi, fB, 0.0, 0.42)
        same = "-" if ref is None else str(np.array_equal(code[A], ref))
        if ref is None:
            ref = code[A].copy()
        det, thr, _ = implied(fam, bim, code, pairs, 2, binary, "tail")
        kk = np.sort(np.array([S.kinship(code, i, j, A) for i, j in pairs]))[::-1]
        print("  %6.2f %8d %10s %10s %10s %12s"
              % (fB, det, fmt(thr), fmt(kk[det - 1] if det else float("nan")),
                 fmt((CUT - 0.5) / (thr - 0.5)) if det else "  n/a", same), flush=True)
    print()
    print("  The kept group's genotypes are bit-identical down the column, so every pair's")
    print("  kinship over them is too -- and the screen's demand on the kept group falls by")
    print("  a factor of two as the pair shares more of what the budget throws away.  The")
    print("  screen reads the pair's genotypes at markers it does not keep.")


# --------------------------------------------------------------------------
# 2. R below 1
# --------------------------------------------------------------------------

def m_below(binary=None):
    """Round 1 saw one beta spectrum at R = 0.998 and called it slack.  It is not."""
    print("Sharing confined to the most informative markers, m = 50000 > budget:")
    print("  %8s %8s %8s %10s %10s" % ("nA", "MAF A", "MAF B", "thr(whole)", "R"))
    for nA, a, b in ((BUDGET, 0.45, 0.25), (BUDGET, 0.45, 0.15), (24576, 0.45, 0.20)):
        spec, hi = twogroup(nA, a, 50000 - nA, b)
        fam, bim, code, pairs, fs, A, B = kept_ladder(spec, hi, 0.0, 0.0, 0.42)
        det, thr, _ = implied(fam, bim, code, pairs, 2, binary, "below")
        print("  %8d %8.2f %8.2f %10s %10s"
              % (nA, a, b, fmt(thr), fmt((CUT - 0.5) / (thr - 0.5))), flush=True)
    print()
    print("  R < 1: the screen accepts pairs whose whole-map kinship is *below* the printed")
    print("  cutoff when their sharing sits high in the spectrum.  A conservative bound can")
    print("  only ever reject pairs a full computation would accept, so this alone retires")
    print("  the 'bound-based early exit' reading of the deflation.")


# --------------------------------------------------------------------------
# 3. the MAF band
# --------------------------------------------------------------------------

def m_band(binary=None):
    """The threshold against the MAF of an unshared block of fixed size and rank."""
    print("m=32768: u=4096 unshared markers @ MAF a, above 8192 @ MAF 0.15 that the pair")
    print("clones at the rung's f, junk @ MAF 0.03 filling the rest.  The block outranks")
    print("the window at every a below, and its informativeness spans only 11 %.")
    for n, degree, w, bmaf in ((150, 2, 8192, 0.15), (300, 2, 8192, 0.15),
                               (600, 2, 8192, 0.15), (300, 2, 16384, 0.20),
                               (300, 1, 16384, 0.20)):
        print()
        print("  n = %d, degree %d (cutoff %.4f)%s:"
              % (n, degree, 2.0 ** -(degree + 2),
                 "" if degree == 2 else ", window widened to %d @ MAF %.2f so the ladder"
                 " can clear the higher cutoff" % (w, bmaf)))
        thr, kmax = block_threshold(0, 0.45, w=w, bmaf=bmaf, n=n, degree=degree,
                                    binary=binary, tag="band")
        print("    no block   %11s  thr = %s   (ladder reaches %.3f)"
              % ("", " ".join(fmt(t) for t in thr), kmax))
        for a in (0.34, 0.38, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.48):
            thr, kmax = block_threshold(4096, a, w=w, bmaf=bmaf, n=n, degree=degree,
                                        binary=binary, tag="band")
            print("    a=%.2f  2pq=%.4f  thr = %s   (ladder reaches %.3f)"
                  % (a, 2 * a * (1 - a), " ".join(fmt(t) for t in thr), kmax), flush=True)
    print()
    print("  Flat at the cutoff to a = 0.42, then a knee: by a = 0.45 the same 4096 unshared")
    print("  markers have more than doubled what the pair must score.  The knee does not")
    print("  move between n = 150 and n = 600, so it is not a binomially smeared cut on the")
    print("  in-sample frequency sitting somewhere else.")


def m_size(binary=None):
    """The same threshold against the block's size, on each side of the knee."""
    print("m=32768, window 8192 @ MAF 0.15, junk @ 0.03; block size u swept:")
    print("  %7s | %-26s | %-26s" % ("u", "a = 0.40 (below the knee)", "a = 0.45 (above)"))
    for u in (0, 1024, 2048, 4096, 8192, 12288, 16384, 20480):
        cells = []
        for a in (0.40, 0.45):
            thr, kmax = block_threshold(u, a, binary=binary, tag="size")
            cells.append("%-26s" % (" ".join(fmt(t) for t in thr)))
        print("  %7d | %s | %s" % (u, cells[0], cells[1]), flush=True)
    print()
    print("  At MAF 0.40 the block is free at every size -- 20480 unshared markers, five")
    print("  eighths of the map, all of them ranked above the pair's sharing, cost nothing.")
    print("  At MAF 0.45 the threshold is already 1.5x the cutoff at 2048 markers and out of")
    print("  reach by 8192.  No scan over a ranking can distinguish those two blocks.")


def m_rank(binary=None):
    """Rank position on its own, with MAF held below the knee."""
    print("Where the sharing sits in the ranking, with every block held at MAF 0.40:")
    print("  %7s %10s %s" % ("u", "win rank", "thr, 3 seeds"))
    for u in (0, 4096, 8192, 16384, 20480):
        thr, kmax = block_threshold(u, 0.40, binary=binary, tag="rank")
        print("  %7d %10d %s   (ladder reaches %.3f)"
              % (u, u, " ".join(fmt(t) for t in thr), kmax), flush=True)
    print()
    print("and the same sweep with the window's own MAF raised into the band, so that its")
    print("markers interleave with the block instead of sitting below it:")
    print("  %7s %10s %s" % ("u", "win MAF", "thr, 3 seeds"))
    for wmaf in (0.15, 0.30, 0.38, 0.42, 0.45):
        thr, kmax = block_threshold(12288, 0.45, bmaf=wmaf, binary=binary, tag="rank2")
        print("  %7d %10.2f %s   (ladder reaches %.3f)"
              % (12288, wmaf, " ".join(fmt(t) for t in thr), kmax), flush=True)
    print()
    print("  The pair is rescued exactly when its own sharing enters the band, not when it")
    print("  climbs the ranking: MAF 0.38 is ranked above two thirds of the map and is still")
    print("  refused, MAF 0.42 and 0.45 are inside the band and pass at the cutoff.")


def m_coding(binary=None):
    """Control: swapping the two homozygote codes is invisible to every pairwise count."""
    print("Swapping hom-A1 and hom-A2 at a marker leaves het, HetHet and IBS0 -- and so")
    print("every kinship -- bit-identical.  If the screen keyed on the unfolded A1")
    print("frequency rather than on MAF, that recoding would move it.")
    fam, bim, code, pairs, lab = block_ladder(8192, 0.45, 8192, 0.15)
    blk = np.flatnonzero(lab == 0)
    d = np.where(code == 0, 2, np.where(code == 2, 1, 0))
    a1 = d.sum(1) / np.maximum(2 * (code != 1).sum(1), 1)
    over = np.intersect1d(np.flatnonzero(a1 > 0.5), blk)
    for name, idx in (("as generated", np.array([], int)),
                      ("A1-major markers flipped", over),
                      ("whole block flipped", blk)):
        c = code.copy()
        c[idx] = np.where(c[idx] == 0, 3, np.where(c[idx] == 3, 0, c[idx]))
        det, thr, kmax = implied(fam, bim, c, pairs, 2, binary, "cod")
        print("  %-26s %4d markers recoded   printed %2d   thr %s"
              % (name, len(idx), det, fmt(thr)), flush=True)
    print("  Identical, as it must be.  The band is a band in MAF, not in allele coding.")


# --------------------------------------------------------------------------
# 4. out of sample
# --------------------------------------------------------------------------

def m_oos(binary=None):
    """Fresh seeds, fresh sizes, fresh spectra, both degrees, above and below the budget."""
    print("Every number below uses seeds, map sizes and MAFs that appear nowhere above.")
    print()
    print("the band, at m = 45056 (above the budget, block 8192, window 16384) and")
    print("m = 20480 (well below it, block 2048, window 8192); window @ MAF 0.22,")
    print("junk @ MAF 0.06, n = 260, seeds 101/202/303:")
    print("  %8s %8s | %s" % ("m", "a", "thr, 3 seeds"))
    for m, u, w in ((45056, 8192, 16384), (20480, 2048, 8192)):
        for a in (0.36, 0.40, 0.43, 0.45, 0.48):
            thr, kmax = block_threshold(u, a, w=w, bmaf=0.22, m=m, n=260,
                                        seeds=(101, 202, 303), binary=binary, tag="oos1")
            print("  %8d %8.2f | %s   (ladder reaches %.3f)"
                  % (m, a, " ".join(fmt(t) for t in thr), kmax), flush=True)
    print()
    print("degree 1 (cutoff 0.1250), m = 40960, block 6144, window 20480 @ MAF 0.22,")
    print("junk @ MAF 0.06, n = 220:")
    print("  %8s | %s" % ("a", "thr, 3 seeds"))
    for a in (0.0, 0.30, 0.34, 0.38, 0.45):
        thr, kmax = block_threshold(0 if a == 0 else 6144, a if a else 0.45, w=20480,
                                    bmaf=0.22, m=40960, n=220, seeds=(101, 202, 303),
                                    degree=1, binary=binary, tag="oos2")
        print("  %8s | %s   (ladder reaches %.3f)"
              % ("no block" if a == 0 else "%.2f" % a,
                 " ".join(fmt(t) for t in thr), kmax), flush=True)
    print()
    print("the discarded-marker read, on a fresh canvas: m = 45056 = 30000 @ MAF 0.45 +")
    print("15056 @ MAF 0.20, n = 240, seed 404:")
    spec, hi = twogroup(30000, 0.45, 15056, 0.20)
    print("  %6s %8s %10s %12s" % ("fB", "printed", "thr(whole)", "kept bits =="))
    ref = None
    for fB in (0.0, 0.10, 0.30, 0.60, 1.00):
        fam, bim, code, pairs, fs, A, B = kept_ladder(spec, hi, fB, 0.0, 0.45,
                                                      n=240, seed=404)
        same = "-" if ref is None else str(np.array_equal(code[A], ref))
        if ref is None:
            ref = code[A].copy()
        det, thr, _ = implied(fam, bim, code, pairs, 2, binary, "oos3")
        print("  %6.2f %8d %10s %12s" % (fB, det, fmt(thr), same), flush=True)


def m_bigish(binary=None):
    """Is the band the deflation?  Delete it from `bigish` and see."""
    fam, bim, code = S.load()
    fid = [line.split()[0] for line in fam]
    n = len(fam)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n) if fid[i] != fid[j]]
    mf = S.mafs(code)
    print("bigish, m = 50000, n = 200.  Drop every marker with in-sample MAF >= c:")
    print("  %6s %8s %9s %7s %9s %7s %9s"
          % ("c", "m left", "printed2", "true2", "printed1", "true1", "R(deg2)"))
    for c in (1.01, 0.45, 0.44, 0.42, 0.40, 0.38, 0.35, 0.30):
        keep = np.flatnonzero(mf < c)
        sub, bb = code[keep], [bim[k] for k in keep]
        bed = S.write(os.path.join(WORK, "bandcut"), fam, bb, sub)
        d2, d1 = S.run(bed, 2, binary)[1], S.run(bed, 1, binary)[1]
        ks = np.sort(np.array([S.kinship(sub, i, j) for i, j in pairs]))[::-1]
        thr = ks[d2 - 1] if 0 < d2 <= len(ks) else float("nan")
        print("  %6.2f %8d %9d %7d %9d %7d %9s"
              % (c, len(keep), d2, int((ks > CUT).sum()), d1, int((ks > 0.125).sum()),
                 fmt((CUT - 0.5) / (thr - 0.5))), flush=True)
    print()
    print("  The deflation tracks m against the 32768 budget and nothing else: it is")
    print("  undisturbed by deleting the whole band (c = 0.45 -> 0.38, R still 1.015-1.022)")
    print("  and it vanishes only when m itself falls under the budget.  Compare §4.3's")
    print("  prefix truncation at the same m -- 1.0087 at 36864, 1.0128 at 40960 -- which")
    print("  the band-deleted maps reproduce.  **The band and the deflation are two")
    print("  different phenomena**, and `bigish`'s uniformly-related pairs share inside the")
    print("  band in proportion, so the band never binds there.")


MEASUREMENTS = {"tail": m_tail, "below": m_below, "band": m_band, "size": m_size,
                "rank": m_rank, "coding": m_coding, "bigish": m_bigish, "oos": m_oos}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("what", nargs="?", default="facts",
                    choices=list(MEASUREMENTS) + ["facts"])
    ap.add_argument("--impl", help="binary to drive instead of the reference")
    args = ap.parse_args(argv)
    os.makedirs(WORK, exist_ok=True)
    for name in (MEASUREMENTS if args.what == "facts" else [args.what]):
        print("== %s ==" % name)
        MEASUREMENTS[name](args.impl)
        print()


if __name__ == "__main__":
    main()
