#!/usr/bin/env python3
"""Read the reference's answer one chromosome at a time, on the corpus's own data.

The instrument behind `docs/research/23-gap-bound.md`.  Every earlier campaign had to
guess which segment of which pair a wrong `.seg` row came from, because `--ibdseg` prints
only two genome-wide proportions and one bad call of 11 Mb in 690 Mb is four decimal
places away from visible.  This makes it visible: run the reference twice on the same
fileset with all but one chromosome **muted** for the probe pair, and the difference in
the printed `IBD1Seg` / `IBD2Seg` is that chromosome's own called length in base pairs.

Two things make it exact rather than approximate.

* **Mute, do not subset.**  Dropping the other chromosomes' `.bim` rows changes the
  answer: KING packs the whole retained marker list into 64-marker words, so deleting
  markers re-phases every later usable segment and moves its calls.  Muting keeps every
  `.bim` row and only rewrites the *probe pair's* genotypes — one sample to `A1A1`, the
  other to `A2A2`, so every muted marker is an opposite homozygote.  A chromosome muted
  that way carries no usable word on either pass and contributes exactly zero, while the
  usable-segment map, the denominator and the word phase stay byte-identical.  (Muting
  nothing reproduces the golden `.seg` row exactly; that is the rig's own control.)
* **Recompute the denominator.**  The printed total is `%.1lf Mb`, rounded to 100 kb.
  `Pallsegs.txt` names each usable segment's first and last SNP, so the exact denominator
  comes back out of the `.bim`.  What is left is the 4-dp rounding of `IBD?Seg` itself,
  about `0.00005 * denom` — 35 kb on `multifam`, well under one marker gap.

Two rig hazards, both bisected and both silent if ignored:

* `--seglength` is **in Mb** and KING clamps it to 1..10 Mb, falling back to its 3 Mb
  default outside that range with no diagnostic.  Passing base pairs by mistake reads as
  "the 3 Mb answer at every floor", which looks exactly like a floor-independent rule.
  `run()` takes base pairs and asserts the range.
* `--ibdseg` refuses a fileset whose usable total is under 100 Mb ("Segments too short"),
  and the `>10 Mb` pair filter drops a pair whose kept chromosomes hold no long segment.
  Muting solves the first; `sweep()` takes a **carrier** chromosome for the second.

    python3 chrprobe.py multifam 2 4 5 10        # per-chromosome table at two floors
    python3 chrprobe.py multifam 2 4 --sweep 13 1 --lo 1e6 --hi 1e7
"""
from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))))
sys.path.insert(0, os.path.join(ROOT, "tests", "parity", "fit"))

import engine as E  # noqa: E402
import kingdata as kd  # noqa: E402
import seg20 as S20  # noqa: E402
import seg21 as S21  # noqa: E402

KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
JOBS = int(os.environ.get("CHRPROBE_JOBS", "8"))


# ---------------------------------------------------------------------------
# the fileset
# ---------------------------------------------------------------------------

_CODES: dict[str, np.ndarray] = {}
_BIMC: dict[str, list] = {}


def codes(ds):
    """The dataset's genotypes as an `(n_variant, n_sample)` matrix of PLINK codes."""
    v = _CODES.get(ds.name)
    if v is not None:
        return v
    src = os.path.join(kd.DATA, ds.name)
    n = len(ds.fam)
    bpv = (n + 3) // 4
    raw = np.fromfile(src + ".bed", dtype=np.uint8)
    body = raw[3:].reshape(-1, bpv)
    c = np.empty((body.shape[0], bpv * 4), dtype=np.uint8)
    for k in range(4):
        c[:, k::4] = (body >> (2 * k)) & 3
    _CODES[ds.name] = v = c[:, :n]
    return v


def bimchr(ds):
    """KING's chromosome code for every `.bim` row, in file order."""
    v = _BIMC.get(ds.name)
    if v is None:
        v = []
        with open(os.path.join(kd.DATA, ds.name + ".bim")) as fh:
            for line in fh:
                v.append(kd.king_chrom_code(line.split()[0]))
        _BIMC[ds.name] = v
    return v


def chrom_mask(ds, chroms):
    return np.array([c in chroms for c in bimchr(ds)])


def write(ds, mat, out):
    """`<out>.{bed,bim,fam}` from a code matrix; `.bim` and `.fam` are copied verbatim."""
    src = os.path.join(kd.DATA, ds.name)
    for ext in (".bim", ".fam"):
        with open(out + ext, "wb") as fh:
            fh.write(open(src + ext, "rb").read())
    nvar, m = mat.shape
    obpv = (m + 3) // 4
    packed = np.zeros((nvar, obpv), dtype=np.uint8)
    for k in range(m):
        packed[:, k >> 2] |= mat[:, k] << (2 * (k & 3))
    with open(out + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        packed.tofile(fh)


def run(prefix, workdir, seglen_bp=None):
    """`(rows keyed by (ID1, ID2), exact denominator in bp, stdout)`.

    `seglen_bp` is `--seglength` in **base pairs**; the flag is in Mb and KING silently
    falls back to 3 Mb outside 1..10 Mb, so the range is asserted here.
    """
    cmd = [KING, "-b", prefix + ".bed", "--ibdseg", "--prefix", "P"]
    if seglen_bp is not None:
        assert 1e6 <= seglen_bp <= 1e7, "--seglength outside KING's 1..10 Mb range"
        cmd += ["--seglength", "%.6f" % (seglen_bp / 1e6)]
    r = subprocess.run(cmd, cwd=workdir, capture_output=True, text=True, check=False)
    pos = {}
    with open(prefix + ".bim") as fh:
        for line in fh:
            f = line.split()
            pos[f[1]] = int(f[3])
    denom = 0
    allsegs = os.path.join(workdir, "Pallsegs.txt")
    if os.path.exists(allsegs):
        with open(allsegs) as fh:
            next(fh)
            for line in fh:
                f = line.split()
                denom += pos[f[7]] - pos[f[6]]
    rows = {}
    seg = os.path.join(workdir, "P.seg")
    if os.path.exists(seg):
        with open(seg) as fh:
            hdr = next(fh).split()
            for line in fh:
                f = line.split()
                rows[(f[1], f[3])] = dict(zip(hdr, f))
    return rows, denom, r.stdout


def measure(dsname, i, j, keep_mask, seglen_bp):
    """`(IBD1 bp, IBD2 bp, denom)` the reference reports for the pair on the kept markers."""
    ds = kd.load(dsname)
    mat = codes(ds).copy()
    off = ~np.asarray(keep_mask)
    mat[off, i] = 0        # A1A1
    mat[off, j] = 3        # A2A2
    with tempfile.TemporaryDirectory() as d:
        pre = os.path.join(d, "sub")
        write(ds, mat, pre)
        rows, denom, _ = run(pre, d, seglen_bp)
        key = (ds.fam[i][1], ds.fam[j][1])
        row = rows.get(key) or rows.get((key[1], key[0]))
        if row is None or not denom:
            return None
        return (float(row["IBD1Seg"]) * denom, float(row["IBD2Seg"]) * denom, denom)


# ---------------------------------------------------------------------------
# what we say
# ---------------------------------------------------------------------------

def ours(dsname, i, j, seglen_bp, p=S21.R21(), caller=None):
    """This caller's per-chromosome `(IBD1, IBD2)` in bp, for the same pair."""
    ds = kd.load(dsname)
    pos = ds.pos
    b = p.base
    out = {}
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = (caller or S21.ibd2_21)(sc, ds, i, j, p, pos, seglen_bp)
        c1 = S20.ibd1_20(sc, ds, i, j, pos, seglen_bp, b)
        a2 = sum(int(pos[hi] - pos[lo]) for lo, hi in c2)
        a1 = 0
        for lo, hi in c1:
            a1 += sum(w for w in (int(pos[y] - pos[x])
                                  for x, y in E._pieces((lo, hi), c2)) if w >= seglen_bp)
        x, y = out.get(seg[0], (0, 0))
        out[seg[0]] = (x + a1, y + a2)
    return out


def carrier_for(dsname, i, j, seglen_bp):
    """A chromosome whose own calls are long enough to keep the pair past the filter."""
    us = ours(dsname, i, j, seglen_bp)
    return max(us, key=lambda c: us[c][0] + us[c][1])


def perchrom(dsname, i, j, seglen_bp, carrier=None, chroms=None, jobs=JOBS):
    """Every chromosome's own `(IBD1, IBD2)` from the reference, by difference."""
    ds = kd.load(dsname)
    carrier = carrier or carrier_for(dsname, i, j, seglen_bp)
    allch = sorted({c for c in bimchr(ds) if c is not None and 1 <= c <= 22})
    plan = [c for c in (chroms or allch) if c != carrier]
    reqs = [{carrier}] + [{carrier, c} for c in plan]
    with ThreadPoolExecutor(jobs) as ex:
        res = list(ex.map(lambda s: measure(dsname, i, j, chrom_mask(ds, s), seglen_bp),
                          reqs))
    base = res[0]
    out = {carrier: (base[0], base[1]) if base else None}
    for c, v in zip(plan, res[1:]):
        out[c] = None if (v is None or base is None) else (v[0] - base[0], v[1] - base[1])
    return out


def report(dsname, i, j, floors_bp, tol=80_000):
    """Print, per floor, every chromosome where the reference and this caller differ."""
    print("%s %d/%d" % (dsname, i, j))
    car = carrier_for(dsname, i, j, max(floors_bp))
    for L in floors_bp:
        ref, us = perchrom(dsname, i, j, L, car), ours(dsname, i, j, L)
        print("--- floor %.1f Mb (carrier chr%d) ---   chr   ref_ibd1   our_ibd1 |"
              "   ref_ibd2   our_ibd2" % (L / 1e6, car))
        hit = False
        for c in sorted(ref):
            r, o = ref[c], us.get(c, (0, 0))
            if r is None:
                print("   %3d   ?" % c)
                continue
            f1, f2 = abs(r[0] - o[0]) > tol, abs(r[1] - o[1]) > tol
            if not f1 and not f2:
                continue
            hit = True
            print("   %3d  %9.4f  %9.4f %s|  %9.4f  %9.4f %s"
                  % (c, r[0] / 1e6, o[0] / 1e6, "<<<" if f1 else "   ",
                     r[1] / 1e6, o[1] / 1e6, "<<<" if f2 else ""))
        if not hit:
            print("     (every chromosome agrees)")


def sweep(dsname, i, j, chrom, carrier, lo, hi, step=None, jobs=JOBS):
    """`--seglength` swept against one chromosome; the jumps are its call lengths."""
    ds = kd.load(dsname)
    mask = chrom_mask(ds, {chrom, carrier})
    base = chrom_mask(ds, {carrier})
    xs = ([lo + step * k for k in range(int(round((hi - lo) / step)) + 1)]
          if step else [lo, hi])
    with ThreadPoolExecutor(jobs) as ex:
        got = list(ex.map(lambda L: (measure(dsname, i, j, mask, L),
                                     measure(dsname, i, j, base, L)), xs))
    out, prev = [], object()
    for L, (v, b) in zip(xs, got):
        d = None if v is None else (v[0] - (b[0] if b else 0), v[1] - (b[1] if b else 0))
        k = None if d is None else (round(d[0] / 1e4), round(d[1] / 1e4))
        if k != prev:
            out.append((L, d))
        prev = k
    return out


def flip(dsname, i, j, chrom, carrier, lo, hi, what=1):
    """Bisect the single `--seglength` at which that chromosome's answer changes."""
    def val(L):
        v = measure(dsname, i, j, chrom_mask(kd.load(dsname), {chrom, carrier}), L)
        return None if v is None else round(v[what] / 1e4)
    a = val(lo)
    assert a != val(hi), "no transition in [%d, %d]" % (lo, hi)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if val(mid) == a:
            lo = mid
        else:
            hi = mid
    return lo, hi


if __name__ == "__main__":
    a = sys.argv[1:]
    dsname, i, j = a[0], int(a[1]), int(a[2])
    if "--sweep" in a:
        k = a.index("--sweep")
        chrom, carrier = int(a[k + 1]), int(a[k + 2])
        lo = int(float(a[a.index("--lo") + 1])) if "--lo" in a else 1_000_000
        hi = int(float(a[a.index("--hi") + 1])) if "--hi" in a else 10_000_000
        for L, v in sweep(dsname, i, j, chrom, carrier, lo, hi, (hi - lo) // 90 or 1):
            print("  L=%9d  IBD1 %10.4f  IBD2 %10.4f"
                  % (L, v[0] / 1e6, v[1] / 1e6) if v else "  L=%9d  ABSENT" % L)
    else:
        report(dsname, i, j, [int(float(x) * 1e6) for x in a[3:]] or [5_000_000,
                                                                     10_000_000])
