#!/usr/bin/env python3
"""Bracket the IBD2 clause of the `--ibdseg --degree` reporting filter.

The corpus says the filter is *not* `PropIBD > 2^-(d+0.5)` alone: the reference keeps
some pairs that fail it.  Every such pair has a non-zero `IBD2Seg` and every pair the
reference drops has `IBD2Seg` exactly 0.0000, but the corpus (and a 53 000-case
`--seglength` x `--degree` sweep over it) never produces an `IBD2Seg` strictly between
0 and 0.1089, so it cannot separate "any IBD2 sharing" from a threshold at 0.08.

This fixture can.  One pair, one IBD2 block of `L` markers on chromosome 1 and nothing
else shared, so `IBD2Seg` is a dial:

    IBD2Seg ~ L * spacing / D          PropIBD ~ IBD2Seg (the pair is IBD0 elsewhere)

with `D` ~ 300 Mb, so a 15 Mb block — comfortably over the fixed 10 Mb "long segment"
filter that decides whether the pair is reported at all — is only IBD2Seg 0.05.  Sweeping
`L` walks IBD2Seg across 0.08 while PropIBD stays far below every degree cut-point.

    python3 gate8.py            # the sweep
"""

import sys

import fixlab as F

CHROMS = [(c, 3000) for c in range(1, 11)]  # 10 x 30 Mb at 10 kb spacing


def probe_block(nmark, degree):
    """(IBD1Seg, IBD2Seg, PropIBD, reported?) for a pair sharing one IBD2 block."""
    fix = F.Fixture(f"g8_{nmark}", CHROMS, nsample=6, maf=0.3, seed=11)
    fix.set_state(0, 500, 500 + nmark, F.IBD2)
    args = [] if degree is None else ["--degree", str(degree)]
    row, _segs, _denom, _out, _wd = F.probe(fix, args, tag=f"_d{degree}")
    if row is None:
        return None
    return (
        float(row["IBD1Seg"]),
        float(row["IBD2Seg"]),
        float(row["PropIBD"]),
    )


def main():
    print(f"{'markers':>8} {'IBD1Seg':>8} {'IBD2Seg':>8} {'PropIBD':>8}  base  d1  d2")
    for nmark in (1200, 1600, 2000, 2200, 2400, 2500, 2600, 2800, 3000):
        base = probe_block(nmark, None)
        if base is None:
            print(f"{nmark:8d}  not reported at all")
            continue
        d1 = probe_block(nmark, 1)
        d2 = probe_block(nmark, 2)
        print(
            f"{nmark:8d} {base[0]:8.4f} {base[1]:8.4f} {base[2]:8.4f}"
            f"   yes  {'yes' if d1 else ' no'} {'yes' if d2 else ' no'}"
        )


if __name__ == "__main__":
    sys.exit(main())
