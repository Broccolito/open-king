"""Does the reference's verdict on a borderline run move with the 64-marker word grid?

`probe_seg.py` establishes that the pairs we over-call are ones the reference does not
call *at all* — not ones it measures shorter — and `matched.py` establishes that no local
summary of the run tells the two groups apart.  What is left to test is the grid itself:
every one of these runs is exactly two complete words wide, but the underlying IBS0-free
stretch of markers is ~165 long, so shifting the global word grid changes how many
complete words fit inside it (2 or 3).

Deleting the first `m` markers of the fileset shifts the grid by `m` and changes nothing
else about the pair.  If the verdict tracks the complete-word count, a word-count rule is
alive and our count is what is wrong; if the verdict never moves, the rule is not a word
count at all.

    python3 gridshift.py [n_pairs]
"""

import os
import subprocess
import sys
import tempfile

import numpy as np

import kingdata as kd
import probe_seg as PS
import rules2 as R
import subset as S

BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")


def write_shifted(ds, keep, drop, target_chr, k, outdir):
    """Subset to `keep` samples, delete the first `drop` markers, stretch `target_chr`."""
    src = os.path.join(kd.DATA, ds.name)
    fam = open(src + ".fam").read().splitlines()
    with open(os.path.join(outdir, "S.fam"), "w") as fh:
        for s in keep:
            fh.write(fam[s] + "\n")
    lines = open(src + ".bim").read().splitlines()[drop:]
    n = 0
    with open(os.path.join(outdir, "S.bim"), "w") as fh:
        for line in lines:
            f = line.split()
            if kd.king_chrom_code(f[0]) == target_chr:
                f[3] = str(int(f[3]) * k)
            else:
                n += 1
                f[3] = str(n * 1000)
            fh.write("\t".join(f) + "\n")
    nfam = len(fam)
    bpv = (nfam + 3) // 4
    raw = np.fromfile(src + ".bed", dtype=np.uint8)
    body = raw[3:].reshape(-1, bpv)[drop:]
    codes = np.empty((body.shape[0], bpv * 4), dtype=np.uint8)
    for b in range(4):
        codes[:, b::4] = (body >> (2 * b)) & 3
    sub = codes[:, keep]
    m = len(keep)
    obpv = (m + 3) // 4
    packed = np.zeros((sub.shape[0], obpv), dtype=np.uint8)
    for c in range(m):
        packed[:, c >> 2] |= sub[:, c] << (2 * (c & 3))
    with open(os.path.join(outdir, "S.bed"), "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        packed.tofile(fh)


def complete_words(ds, i, j, lo, hi, drop):
    """Complete clean words inside the IBS0-free stretch once the grid shifts by `drop`."""
    ibs0, _, _, _ = ds.masks(i, j)

    def is0(t):
        return bool((int(ibs0[t // 64]) >> (t % 64)) & 1)

    left = lo
    while left > 0 and not is0(left - 1):
        left -= 1
    right = hi
    while right + 1 < ds.pos.size and not is0(right + 1):
        right += 1
    a, b = left - drop, right - drop           # indices in the shifted array
    return max(0, (b + 1) // 64 - -(-a // 64))


def main(n_pairs=3):
    n_pairs = int(n_pairs)
    ds = kd.load("bigish")
    groups = {"reported": [], "extra": []}
    for (i, j) in ds.pairs():
        _, _, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        u, v = (best[1] + 63) // 64, (best[2] + 1) // 64 - 1
        if v - u + 1 != 2 or len(detail) != 1:
            continue
        groups["reported" if (i, j) in ds.ref else "extra"].append((i, j, best[1], best[2]))

    probes = groups["reported"][:n_pairs] + groups["extra"][:n_pairs]
    labels = ["reported"] * n_pairs + ["extra"] * n_pairs
    used = {x for t in probes for x in t[:2]}
    pad = [s for s in range(len(ds.fam)) if s not in used][:28]

    shifts = [0, 8, 16, 24, 32, 40, 48, 56]
    print(f"{'group':<9} {'pair':<22} " + " ".join(f"m={m:<5}" for m in shifts))
    for label, (i, j, lo, hi) in zip(labels, probes):
        chrom = int(ds.chr[lo])
        keep = sorted(set([i, j] + pad))
        cells = []
        for m in shifts:
            with tempfile.TemporaryDirectory() as td:
                write_shifted(ds, keep, m, chrom, 2, td)
                rows, denom = PS.run_king(td)
            key = (ds.fam[i][1], ds.fam[j][1])
            r = rows.get(key) or rows.get((key[1], key[0]))
            w = complete_words(ds, i, j, lo, hi, m)
            mb = f"{r[0] * denom / 2 / 1e6:.2f}" if r else "-"
            cells.append(f"{w}w/{mb}")
        print(f"{label:<9} {'/'.join((ds.fam[i][1], ds.fam[j][1])):<22} "
              + " ".join(f"{c:<7}" for c in cells))


if __name__ == "__main__":
    main(*sys.argv[1:])
