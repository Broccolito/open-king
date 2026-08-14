"""Does marker informativeness separate the 188 extra pairs from the reported ones?

`docs/research/10-segment-rule-fixtures.md` §5 shows on constructed filesets that the
reference refuses to call an IBS0-free run when the markers under it are uninformative,
and that the discriminator is the **sample** allele frequency, not the pair's genotypes.
`docs/research/11-segment-rule-fit.md` §7 shows that no *pair* feature separates the 182
`bigish` pairs we over-call from the 255 comparable ones the reference does report.

This script tests the frequency hypothesis on the real corpus: for every pair whose
longest called segment is the deciding one, it measures candidate informativeness
statistics over that segment and prints their distribution split by whether the reference
reported the pair.
"""

import sys

import numpy as np

import kingdata as kd
import rules2 as R

BEST = R.P(t1=0, min1=2, bridge1=0, t2=0, min2=1, bridge2=0, edge="fringe", end2="next")


def allele_freq(ds):
    """Per-marker frequency of the A1 allele, over all samples, non-missing only."""
    n = len(ds.fam)
    nmark = ds.pos.size
    bits = np.zeros((n, ds.nwords * 64), dtype=np.uint8)
    for k in range(n):
        b0 = np.unpackbits(ds.p0[k].view(np.uint8), bitorder="little")
        b1 = np.unpackbits(ds.p1[k].view(np.uint8), bitorder="little")
        # plane encoding: (1,1)=hom A1A1 -> 2 copies, (0,1)=het -> 1, (1,0)=hom A2A2 -> 0,
        # (0,0) = missing
        called = (b0 | b1).astype(bool)
        copies = (b0 & b1) * 2 + ((~b0.astype(bool)) & b1.astype(bool)) * 1
        bits[k] = copies + 4 * (~called)
    copies = np.where(bits[:, :nmark] < 4, bits[:, :nmark], 0).astype(np.float64)
    called = (bits[:, :nmark] < 4).astype(np.float64)
    tot = called.sum(axis=0) * 2
    p = np.divide(copies.sum(axis=0), np.maximum(tot, 1e-9))
    return p


def features(ds, p, lo, hi):
    q = 1 - p
    sl = slice(lo, hi + 1)
    pp, qq = p[sl], q[sl]
    maf = np.minimum(pp, qq)
    return {
        "n": hi - lo + 1,
        "E_ibs0": float(np.sum(2 * pp * pp * qq * qq)),
        "sum2pq": float(np.sum(2 * pp * qq)),
        "n_maf05": int(np.sum(maf >= 0.05)),
        "n_maf10": int(np.sum(maf >= 0.10)),
        "n_maf20": int(np.sum(maf >= 0.20)),
        "mean_maf": float(np.mean(maf)),
    }


def main(name="bigish"):
    ds = kd.load(name)
    p = allele_freq(ds)
    groups = {"reported": [], "extra": []}
    for (i, j) in ds.pairs():
        a, b, longest, detail = R.call_pair(ds, i, j, BEST, want=True)
        if longest < BEST.long_bp:
            continue
        best = max(detail, key=lambda d: int(ds.pos[d[2]] - ds.pos[d[1]]))
        f = features(ds, p, best[1], best[2])
        f["len"] = int(ds.pos[best[2]] - ds.pos[best[1]])
        f["words"] = (best[2] // 64) - (best[1] + 63) // 64 + 1
        key = "reported" if (i, j) in ds.ref else "extra"
        groups[key].append(f)

    keys = ["len", "n", "words", "E_ibs0", "sum2pq", "n_maf05", "n_maf10", "n_maf20",
            "mean_maf"]
    print(f"{name}: reported={len(groups['reported'])} extra={len(groups['extra'])}")
    print(f"{'feature':<10} {'group':<9} {'min':>10} {'q1':>10} {'med':>10} "
          f"{'q3':>10} {'max':>10}")
    for k in keys:
        for g in ("reported", "extra"):
            v = np.array([f[k] for f in groups[g]], dtype=float)
            if v.size == 0:
                continue
            qs = np.percentile(v, [0, 25, 50, 75, 100])
            print(f"{k:<10} {g:<9} " + " ".join(f"{x:10.4g}" for x in qs))
        print()


if __name__ == "__main__":
    main(*sys.argv[1:])
