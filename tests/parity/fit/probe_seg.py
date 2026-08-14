"""Make the reference print the length of a segment it would normally hide.

`subset.py` shows the over-called pairs are excluded for a reason carried entirely by
their own genotypes: a 30-sample subset reproduces the full run's verdict exactly.  So the
reference is either calling a *shorter* segment there than we do, or not calling it at all
— and a pair whose longest segment is under 10 Mb prints nothing, which is why the gap has
been invisible.

The instrument removes the 10 Mb floor without touching a single genotype:

* keep every marker, so the 64-marker word grid — which is indexed by position in the
  retained-autosome array, not by base pair — is byte-identical to the full run;
* **stretch** the target chromosome's base-pair positions by `K`, so a segment of length
  `L` measures `K·L` and clears the fixed 10 Mb pair filter;
* **compress** every other chromosome to 1 kb spacing so it falls out of the usable-segment
  list entirely, leaving the denominator equal to the target chromosome alone. That buys
  back the print resolution the stretch costs: one ulp of `%.4lf` is `D/10000`.

`K = 2` is the largest safe stretch for this corpus: a usable segment is cut wherever one
64-marker word spans over 10 Mb, and 64 x 2 x 50 kb = 6.4 Mb still clears it.

    python3 probe_seg.py [n_probes] [K]
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

import kingdata as kd
import rules2 as R
import subset as S

KING = S.KING
BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")


def write_bim(ds, target_chr, k, out):
    """Target chromosome stretched by `k`; every other chromosome squashed to 1 kb."""
    src = os.path.join(kd.DATA, ds.name + ".bim")
    lines = open(src).read().splitlines()
    n = 0
    with open(out, "w") as fh:
        for line in lines:
            f = line.split()
            chrom = kd.king_chrom_code(f[0])
            if chrom == target_chr:
                f[3] = str(int(f[3]) * k)
            else:
                n += 1
                f[3] = str(n * 1000)
            fh.write("\t".join(f) + "\n")


def run_king(workdir):
    subprocess.run([KING, "-b", "S.bed", "--ibdseg", "--prefix", "P", "--seglength", "1"],
                   cwd=workdir, capture_output=True, check=False)
    seg = os.path.join(workdir, "P.seg")
    allsegs = os.path.join(workdir, "Pallsegs.txt")
    denom = 0
    if os.path.exists(allsegs):
        for line in open(allsegs).readlines()[1:]:
            f = line.rstrip("\n").split("\t")
            if int(f[1]) < 23:
                denom += round(float(f[4]) * 1e6)
    rows = {}
    if os.path.exists(seg):
        for line in open(seg).readlines()[1:]:
            f = line.rstrip("\n").split("\t")
            rows[(f[1], f[3])] = (float(f[4]), float(f[5]), float(f[6]))
    return rows, denom


def main(n_probes=10, k=2):
    n_probes, k = int(n_probes), int(k)
    ds = kd.load("bigish")
    groups = {"reported": [], "extra": []}
    for (i, j) in ds.pairs():
        a, b, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        u, v = (best[1] + 63) // 64, (best[2] + 1) // 64 - 1
        if v - u + 1 != 2 or len(detail) != 1:
            continue
        groups["reported" if (i, j) in ds.ref else "extra"].append(
            (i, j, best[1], best[2], u, v))

    used = set()
    for g in groups.values():
        for t in g[:n_probes]:
            used |= {t[0], t[1]}
    pad = [s for s in range(len(ds.fam)) if s not in used][:28]

    print(f"K={k}.  lengths in Mb on the ORIGINAL scale (ref value divided by K).")
    print(f"{'group':<9} {'pair':<24} {'core':>7} {'ext':>7} {'ref':>8} {'refIBD2':>8}"
          f" {'ref-core':>9}")
    for label, group in (("reported", groups["reported"]), ("extra", groups["extra"])):
        for (i, j, lo, hi, u, v) in group[:n_probes]:
            chrom = int(ds.chr[lo])
            keep = sorted(set([i, j] + pad))
            core = int(ds.pos[min(64 * v + 63, hi)] - ds.pos[64 * u]) / 1e6
            ext = int(ds.pos[hi] - ds.pos[lo]) / 1e6
            with tempfile.TemporaryDirectory() as td:
                S.write_subset(ds, keep, os.path.join(td, "S"))
                write_bim(ds, chrom, k, os.path.join(td, "S.bim"))
                rows, denom = run_king(td)
            key = (ds.fam[i][1], ds.fam[j][1])
            r = rows.get(key) or rows.get((key[1], key[0]))
            if r is None:
                print(f"{label:<9} {'/'.join(key):<24} {core:7.3f} {ext:7.3f} "
                      f"{'absent':>8}")
                continue
            l1 = r[0] * denom / k / 1e6
            l2 = r[1] * denom / k / 1e6
            print(f"{label:<9} {'/'.join(key):<24} {core:7.3f} {ext:7.3f} {l1:8.3f} "
                  f"{l2:8.3f} {l1 - core:9.3f}")
        print()


if __name__ == "__main__":
    main(*sys.argv[1:])
