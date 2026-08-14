"""Batch word-grid sweep -> a labelled table the acceptance rule has to explain.

For each probe pair and each shift `m` of the global 64-marker grid, this records the
geometry of the pair's borderline IBS0-free stretch *in that grid* together with the
reference's verdict.  Genotypes never change along a sweep, so any feature that varies is
a pure function of the alignment — which is exactly the kind of feature the rule has to be
built from.

Output: CSV on stdout, one row per (pair, shift).

    python3 sweepdata.py OUT.csv [n_reported] [n_extra] [shifts]
"""

import sys
import tempfile

import numpy as np

import kingdata as kd
import gridshift as G
import probe_seg as PS
import rules2 as R

BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")
PC = np.bitwise_count


def geometry(ds, i, j, anchor, m):
    """Everything about the clean stretch around `anchor` on the grid shifted by `m`."""
    ibs0, _, _, _ = ds.masks(i, j)

    def is0(t):
        return 0 <= t < ds.pos.size and bool((int(ibs0[t // 64]) >> (t % 64)) & 1)

    left = anchor
    while left > 0 and not is0(left - 1):
        left -= 1
    right = anchor
    while right + 1 < ds.pos.size and not is0(right + 1):
        right += 1
    u = -(-(left - m) // 64)                 # first complete word (shifted index)
    v = (right - m + 1) // 64 - 1            # last complete word
    g = {"clean_n": right - left + 1, "W": v - u + 1, "u": u,
         "clean_bp": int(ds.pos[right] - ds.pos[left])}
    if v < u:
        return g
    g["head"] = 64 * u + m - left            # clean markers in word u-1
    g["tail"] = right - (64 * (v + 1) + m) + 1   # clean markers in word v+1
    prev = [t for t in range(64 * (u - 1) + m, 64 * u + m) if is0(t)]
    nxt = [t for t in range(64 * (v + 1) + m, 64 * (v + 2) + m) if is0(t)]
    g["n0_prev"] = len(prev)
    g["n0_next"] = len(nxt)
    g["lastbit_prev"] = (prev[-1] - m) % 64 if prev else -1
    g["firstbit_next"] = (nxt[0] - m) % 64 if nxt else -1
    g["lastbit_next"] = (nxt[-1] - m) % 64 if nxt else -1
    lo = prev[-1] + 1 if prev else 64 * u + m
    hi = nxt[-1] if nxt else 64 * (v + 1) + 63 + m
    g["lo"], g["hi"] = lo, hi
    g["len"] = int(ds.pos[hi] - ds.pos[lo])
    g["core_len"] = int(ds.pos[min(64 * (v + 1) + 63 + m, ds.pos.size - 1)]
                        - ds.pos[64 * u + m])
    return g


def main(out="sweep.csv", n_rep=3, n_extra=3, shifts=64):
    n_rep, n_extra, shifts = int(n_rep), int(n_extra), int(shifts)
    ds = kd.load("bigish")
    groups = {"reported": [], "extra": []}
    for (i, j) in ds.pairs():
        _, _, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp or len(detail) != 1:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        ibs0, _, _, _ = ds.masks(i, j)
        n0 = PC(ibs0)
        u = -(-best[1] // 64)
        while u - 1 >= 0 and n0[u - 1] == 0:
            u -= 1
        v = (best[2] + 1) // 64 - 1
        while n0[v] != 0:
            v -= 1
        while v + 1 < ds.nwords and n0[v + 1] == 0:
            v += 1
        if v - u + 1 != 2:
            continue
        groups["reported" if (i, j) in ds.ref else "extra"].append(
            (i, j, 64 * u, int(ds.chr[best[1]])))

    probes = ([("reported", *t) for t in groups["reported"][:n_rep]]
              + [("extra", *t) for t in groups["extra"][:n_extra]])
    used = {x for t in probes for x in t[1:3]}
    pad = [s for s in range(len(ds.fam)) if s not in used][:28]

    cols = ["pair", "label", "m", "clean_n", "clean_bp", "W", "u", "head", "tail",
            "n0_prev", "n0_next", "lastbit_prev", "firstbit_next", "lastbit_next",
            "lo", "hi", "len", "core_len", "called", "reflen"]
    fh = open(out, "w")
    fh.write(",".join(cols) + "\n")
    for label, i, j, anchor, chrom in probes:
        name = f"{ds.fam[i][1]}_{ds.fam[j][1]}"
        keep = sorted(set([i, j] + pad))
        for m in range(shifts):
            g = geometry(ds, i, j, anchor, m)
            with tempfile.TemporaryDirectory() as td:
                G.write_shifted(ds, keep, m, chrom, 2, td)
                rows, denom = PS.run_king(td)
            key = (ds.fam[i][1], ds.fam[j][1])
            r = rows.get(key) or rows.get((key[1], key[0]))
            g.update(pair=name, label=label, m=m, called=int(r is not None),
                     reflen=int(r[0] * denom / 2) if r else -1)
            fh.write(",".join(str(g.get(c, -1)) for c in cols) + "\n")
            print(".", end="", flush=True)
    fh.close()
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main(*sys.argv[1:])
