"""Does the gate explain the word-grid sweep?  (The sharpest test available.)

`docs/research/13-segment-acceptance.md` §5 shows the reference's verdict on a borderline
run **moves with the alignment of the 64-marker grid** while the pair's genotypes, the
refined endpoints `lo`/`hi` and the segment's length in bp all stay fixed.  No count of
words, markers or base pairs can explain that.

The informativeness gate can: it counts informative markers over the run's *own complete
words*, and shifting the grid moves markers across those word boundaries even when the
reported endpoints do not move.

Deleting the first `m` markers of the fileset shifts the global grid by `m` and changes
nothing else, so each pair yields 64 labelled observations with genotypes held constant.

    python3 gate_sweep.py [n_reported] [n_extra] [shifts]
"""

import sys
import tempfile

import numpy as np

import kingdata as kd
import gridshift as G
import probe_seg as PS
import rules2 as R2
import rules3 as R3
import sweepdata as SW

BEST = R2.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")
PC = np.bitwise_count
GATE = 10


def informative_markers(ds, i, j):
    """Marker-level bool: both carry A1 and at least one of them is homozygous."""
    m1, _ = R3.inf_masks(ds, i, j)
    bits = np.unpackbits(m1.view(np.uint8), bitorder="little").astype(bool)
    return bits[:ds.pos.size]


def borderline(ds, n_rep, n_extra):
    groups = {"reported": [], "extra": []}
    for (i, j) in ds.pairs():
        _, _, longest, detail = R2.call_pair(ds, i, j, BEST, want=True)
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
        g = groups["reported" if (i, j) in ds.ref else "extra"]
        if len(g) < max(n_rep, n_extra):
            g.append((i, j, 64 * u, int(ds.chr[best[1]])))
    return groups["reported"][:n_rep], groups["extra"][:n_extra]


def main(n_rep=4, n_extra=4, shifts=64):
    n_rep, n_extra, shifts = int(n_rep), int(n_extra), int(shifts)
    ds = kd.load("bigish")
    rep, ext = borderline(ds, n_rep, n_extra)
    probes = [("reported", *t) for t in rep] + [("extra", *t) for t in ext]
    used = {x for t in probes for x in t[1:3]}
    pad = [s for s in range(len(ds.fam)) if s not in used][:28]

    n_ok = n_tot = 0
    confusion = {(0, 0): 0, (0, 1): 0, (1, 0): 0, (1, 1): 0}
    for label, i, j, anchor, chrom in probes:
        inf = informative_markers(ds, i, j)
        cum = np.concatenate(([0], np.cumsum(inf)))
        name = f"{ds.fam[i][1]}/{ds.fam[j][1]}"
        keep = sorted(set([i, j] + pad))
        line = []
        for m in range(shifts):
            g = SW.geometry(ds, i, j, anchor, m)
            with tempfile.TemporaryDirectory() as td:
                G.write_shifted(ds, keep, m, chrom, 2, td)
                rows, _denom = PS.run_king(td)
            key = (ds.fam[i][1], ds.fam[j][1])
            called = int((rows.get(key) or rows.get((key[1], key[0]))) is not None)
            if g["W"] < 2:
                pred, cnt = 0, 0
            else:
                a = 64 * g["u"] + m
                b = 64 * (g["u"] + g["W"]) + m          # exclusive
                cnt = int(cum[min(b, cum.size - 1)] - cum[max(a, 0)])
                pred = int(cnt >= GATE)
            confusion[(called, pred)] += 1
            n_ok += called == pred
            n_tot += 1
            line.append(("Y" if called else ".") + ("Y" if pred else "."))
        print(f"{label:<9} {name:<20} " + " ".join(line))
        print(f"{'':<9} {'':<20} " + " ".join(
            f"{SW.geometry(ds, i, j, anchor, m).get('W', 0)}{'':1}" for m in range(shifts)))
    print(f"\nrefVerdict/prediction agreement: {n_ok}/{n_tot} = {n_ok / n_tot:.4f}")
    print(f"confusion (reference, predicted): {confusion}")


if __name__ == "__main__":
    main(*sys.argv[1:])
