"""Ask the reference directly whether an over-called pair's segment really is > 10 Mb.

The 188 extra rows are indistinguishable from the reported ones on every summary of the
called segment (`matched.py`).  Two explanations survive: either the reference measures
those segments shorter than we do, or it excludes the pair for a reason that is not about
that segment at all.

This separates them.  For each probe pair we write a fileset containing that pair plus a
fixed padding cohort — same `.bim`, so the usable segments and the denominator are
byte-identical to the full run — and re-run the reference.  If a pair the full run omitted
comes back with a row here, the omission was never about its segment length.

    python3 subset.py [n_probes] [n_pad]
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

import kingdata as kd
import rules2 as R

KING = "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")


def write_subset(ds, keep, out):
    """Write <out>.{bed,bim,fam} holding only samples `keep`, in that order."""
    src = os.path.join(kd.DATA, ds.name)
    fam = open(src + ".fam").read().splitlines()
    with open(out + ".fam", "w") as fh:
        for s in keep:
            fh.write(fam[s] + "\n")
    with open(out + ".bim", "wb") as fh:
        fh.write(open(src + ".bim", "rb").read())
    n = len(fam)
    bpv = (n + 3) // 4
    raw = np.fromfile(src + ".bed", dtype=np.uint8)
    body = raw[3:].reshape(-1, bpv)
    nvar = body.shape[0]
    codes = np.empty((nvar, bpv * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (body >> (2 * k)) & 3
    sub = codes[:, keep]
    m = len(keep)
    obpv = (m + 3) // 4
    packed = np.zeros((nvar, obpv), dtype=np.uint8)
    for k in range(m):
        packed[:, k >> 2] |= sub[:, k] << (2 * (k & 3))
    with open(out + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        packed.tofile(fh)


def run_king(prefix, workdir):
    subprocess.run([KING, "-b", prefix + ".bed", "--ibdseg", "--prefix", "P"],
                   cwd=workdir, capture_output=True, check=False)
    path = os.path.join(workdir, "P.seg")
    rows = {}
    if os.path.exists(path):
        for line in open(path).readlines()[1:]:
            f = line.rstrip("\n").split("\t")
            rows[(f[1], f[3])] = (float(f[4]), float(f[5]), float(f[6]))
    return rows


def main(n_probes=12, n_pad=28):
    n_probes, n_pad = int(n_probes), int(n_pad)
    ds = kd.load("bigish")
    reported, extra = [], []
    for (i, j) in ds.pairs():
        a, b, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        u = (best[1] + 63) // 64
        v = (best[2] + 1) // 64 - 1
        if v - u + 1 != 2 or len(detail) != 1:
            continue
        (reported if (i, j) in ds.ref else extra).append((i, j, longest))
    reported.sort(key=lambda t: t[2])
    extra.sort(key=lambda t: t[2])
    print(f"candidates: {len(reported)} reported, {len(extra)} extra "
          f"(single 2-word-core segment)")

    # A padding cohort that is in no probe pair, taken from the far end of the .fam.
    used = {i for i, j, _ in reported[:n_probes] + extra[:n_probes]}
    used |= {j for i, j, _ in reported[:n_probes] + extra[:n_probes]}
    pad = [s for s in range(len(ds.fam)) if s not in used][:n_pad]

    print(f"{'group':<9} {'pair':<22} {'ours Mb':>9} {'ref row':>9} {'ref Mb':>9}")
    for label, group in (("reported", reported), ("extra", extra)):
        hit = 0
        for (i, j, longest) in group[:n_probes]:
            keep = sorted(set([i, j] + pad))
            with tempfile.TemporaryDirectory() as td:
                write_subset(ds, keep, os.path.join(td, "S"))
                rows = run_king("S", td)
            k = (ds.fam[i][1], ds.fam[j][1])
            k2 = (ds.fam[j][1], ds.fam[i][1])
            r = rows.get(k) or rows.get(k2)
            hit += r is not None
            ref_mb = f"{r[0] * ds.denom / 1e6:9.3f}" if r else "        -"
            print(f"{label:<9} {k[0]}/{k[1]:<12} {longest / 1e6:9.3f} "
                  f"{'YES' if r else 'no':>9} {ref_mb}")
        print(f"  -> {label}: {hit}/{min(n_probes, len(group))} reported in subset\n")


if __name__ == "__main__":
    main(*sys.argv[1:])
