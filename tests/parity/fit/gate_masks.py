"""Does the corpus independently confirm the *shape* of the informative mask?

The fixtures say the IBD1 mask is "both carry A1 AND at least one of them is homozygous".
Several plausible alternatives are scored here over all 982 rows, each with its own best
constant, so the comparison is not rigged by holding C at the fitted value.

For the extra/missing metric only the *pair-inclusion* decision matters, and that decision
is monotone in C: a pair is reported iff some run that clears the 10 Mb filter carries at
least C informative markers.  So one pass per mask records

    M(pair) = max over runs longer than 10 Mb of (informative markers under that run)

and every constant is then scored by comparing `M >= C` against the reference.  That also
prints the *margin*: how close the corpus comes to the fitted threshold from each side.
"""

import numpy as np

import kingdata as kd
import rules2 as R2
import rules3 as R3
from gate_corpus import BASE

PC = np.bitwise_count
WORD = 64
FIT = dict(BASE)
FIT["min1"] = 1
P = R3.P3(**FIT, gate=0)

MASKS = {
    "both carry A1 & >=1 hom  (fitted)":
        lambda p0i, p1i, p0j, p1j: p1i & p1j & (p0i | p0j),
    "both carry A1":
        lambda p0i, p1i, p0j, p1j: p1i & p1j,
    "both A1A1":
        lambda p0i, p1i, p0j, p1j: p0i & p1i & p0j & p1j,
    "HetHet":
        lambda p0i, p1i, p0j, p1j: (~p0i & p1i) & (~p0j & p1j),
    ">=1 hom, either allele":
        lambda p0i, p1i, p0j, p1j: (p0i | p0j) & (p0i | p1i) & (p0j | p1j),
    "both non-missing (N_SNP)":
        lambda p0i, p1i, p0j, p1j: (p0i | p1i) & (p0j | p1j),
    "both carry A1 & both hom":
        lambda p0i, p1i, p0j, p1j: p1i & p1j & p0i & p0j,
    ">=1 carries A1":
        lambda p0i, p1i, p0j, p1j: (p1i | p1j) & (p0i | p1i) & (p0j | p1j),
}


def maxinf(ds, i, j, fn):
    """max informative markers over the runs whose reported segment clears 10 Mb."""
    ibs0, n0, n1 = R2.counts(ds, i, j)
    m = fn(ds.p0[i], ds.p1[i], ds.p0[j], ds.p1[j])
    cum = np.concatenate(([0], np.cumsum(PC(m).astype(np.int64))))
    best = -1
    for _, lo, hi in ds.segs:
        w0, w1 = -(-lo // WORD), (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        for ok, tag in ((n0[w0:w1 + 1] <= P.t1, 1),
                        ((n0[w0:w1 + 1] <= P.t1) & (n1[w0:w1 + 1] <= P.t2), 2)):
            for a, b in R2._runs(ok):
                u, v = w0 + a, w0 + b
                lo_m, hi_m = R2._bounds(ibs0, u, v, w0, w1, lo, hi, P.edge,
                                        "next" if tag == 1 else P.end2)
                if int(ds.pos[hi_m] - ds.pos[lo_m]) >= P.long_bp:
                    best = max(best, int(cum[v + 1] - cum[u]))
    return best


def main():
    head = "margin (max dropped / min called)"
    print(f"{'mask':<36}{'bestC':>6}{'extra':>7}{'miss':>6}   {head:>34}")
    for label, fn in MASKS.items():
        obs = []
        for name in kd.DATASETS:
            ds = kd.load(name)
            for (i, j) in ds.pairs():
                M = maxinf(ds, i, j, fn)
                if M < 0:
                    continue
                obs.append((M, (i, j) in ds.ref))
        arr = np.array([[m, int(r)] for m, r in obs])
        best = None
        for C in range(0, int(arr[:, 0].max()) + 2):
            pred = arr[:, 0] >= C
            extra = int(((pred == 1) & (arr[:, 1] == 0)).sum())
            miss = int(((pred == 0) & (arr[:, 1] == 1)).sum())
            key = (extra + miss, C)
            if best is None or key < best[0]:
                best = (key, C, extra, miss)
        _, C, extra, miss = best
        hi_drop = arr[arr[:, 1] == 0][:, 0]
        lo_call = arr[arr[:, 1] == 1][:, 0]
        margin = (f"{hi_drop.max() if hi_drop.size else '-'} / "
                  f"{lo_call.min() if lo_call.size else '-'}")
        print(f"{label:<36}{C:>6}{extra:>7}{miss:>6}   {margin:>28}")


if __name__ == "__main__":
    main()
