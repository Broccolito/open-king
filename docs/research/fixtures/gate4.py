"""Pin the constant, the counting window and the per-genotype weight.

`gate3.py` gives a pass/fail ladder whose boundary is a *total* over the run, not a rate:

    W       2   3   4   6   8  10  14
    k_min   5   4   3   2   2   1   1
    k*W    10  12  12  12  16  10  14      <- smallest passing total is 10 everywhere

This script settles (1) that the constant is exactly 10, (2) whether the count is taken
over the run's own words or over the reported (extended) interval, (3) what each pair
genotype is worth, and (4) whether missing calls matter.

    python3 gate4.py exact
    python3 gate4.py window
    python3 gate4.py weight
    python3 gate4.py miss
"""

import sys

import gatelab as G

R = G.GateRig(spacing=100_000, n1=640, n2=1600, nsample=6, seed=1)

INFO = [2, 1, 0, 0, 0, 0]      # A1A1 / het   -> both carry A1
FILL = [1, 0, 1, 1, 1, 1]      # het  / A2A2  -> only one carries A1
IBS0 = [2, 0, 1, 1, 0, 0]      # A1A1 / A2A2  -> opposite homozygotes


def full(W):
    return 64 * (W + 1) - 1


def place(W, n, positions=None):
    """n informative markers at explicit block-local positions (default: evenly spread)."""
    if positions is None:
        positions = [(i * 64 * W) // max(n, 1) for i in range(n)]
    s = set(positions)
    return lambda t: list(INFO) if t in s else list(FILL)


def run(name, W, kinds, start_word=1, pre=None, post=None):
    """pre/post: optional overrides for the flanking background words."""
    f = R.new(name)
    R.block(f, start_word, W, kinds)
    lo, _ = f.chrom_span(1)
    for wordoff, fn in ((start_word - 1, pre), (start_word + W, post)):
        if fn is None:
            continue
        for b in range(64):
            m = lo + wordoff * 64 + b
            f.force_ibs0.discard(m)
            f.pat_all[m] = fn(b)
    r = R.read(f)
    return r["chr2_mk"] if r else 0


def exact():
    """Total informative markers vs the verdict, at several widths and layouts."""
    print("informative markers placed by explicit total; filler = het/A2A2")
    print(f"{'W':>3} {'layout':>10} " + " ".join(f"{n:>4}" for n in range(6, 15)))
    for W in (2, 3, 4, 8, 14):
        for layout in ("spread", "packed", "split"):
            cells = []
            for n in range(6, 15):
                if layout == "spread":
                    pos = [(i * 64 * W) // n for i in range(n)]
                elif layout == "packed":
                    pos = list(range(n))
                else:                       # half at each end of the run
                    pos = list(range(n // 2)) + \
                        list(range(64 * W - (n - n // 2), 64 * W))
                got = run(f"g4_e_{W}_{layout}_{n}", W, place(W, n, pos))
                cells.append("Y" if got else ".")
            print(f"{W:>3} {layout:>10} " + " ".join(f"{c:>4}" for c in cells))


def window():
    """Do informative markers in the *flanking* (IBS0-bearing) words count?

    The left flanking word is IBS0 at bits 0..63-j and informative at the last j bits, so
    the refined start `lo` = 64u - j puts those j markers inside the reported segment.
    The right flanking word is informative at its first j bits and IBS0 afterwards, so its
    last IBS0 is bit 63 and the whole word — including those j markers — is inside.
    """
    W = 2
    print(f"W={W}: core carries `core` informative markers, each flanking word `j`")
    print(f"{'core':>5} {'j':>3} {'total in core':>14} {'total in [lo,hi]':>17} {'chr2':>6}")
    for core in (6, 8, 9, 10):
        for j in (0, 1, 2, 4):
            def pre(b, j=j):
                return list(INFO) if b >= 64 - j else list(IBS0)

            def post(b, j=j):
                return list(INFO) if b < j else list(IBS0)
            got = run(f"g4_w_{core}_{j}", W, place(W, core), pre=pre, post=post)
            print(f"{core:>5} {j:>3} {core:>14} {core + 2 * j:>17} {got:>6}")
        print()


def weight():
    """What is each pair genotype worth?  Threshold at W=4 is 10 markers of weight 1."""
    cands = {
        "A1A1/het   (2,1)": [2, 1, 0, 0, 0, 0],
        "het/A1A1   (1,2)": [1, 2, 0, 0, 0, 0],
        "het/het    (1,1)": [1, 1, 0, 0, 0, 0],
        "A1A1/A1A1  (2,2)": [2, 2, 1, 0, 0, 0],
        "het/A2A2   (1,0)": [1, 0, 1, 1, 1, 1],
        "A2A2/A2A2  (0,0)": [0, 0, 2, 2, 1, 1],
        "A1A1/miss  (2,-)": [2, 3, 0, 0, 0, 0],
        "het/miss   (1,-)": [1, 3, 1, 1, 0, 0],
        "miss/miss  (-,-)": [3, 3, 1, 1, 1, 1],
    }
    W = 4
    print(f"W={W}: n candidate markers among het/A2A2 filler.  A weight-1 marker "
          f"first passes at n=10.")
    print(f"{'candidate':>20} " + " ".join(f"{n:>4}" for n in range(3, 13)))
    for label, vec in cands.items():
        cells = []
        for n in range(3, 13):
            pos = [(i * 64 * W) // n for i in range(n)]
            s = set(pos)
            got = run(f"g4_wt_{abs(hash(label)) % 100000}_{n}", W,
                      lambda t, v=vec, s=s: list(v) if t in s else list(FILL))
            cells.append("Y" if got else ".")
        print(f"{label:>20} " + " ".join(f"{c:>4}" for c in cells))


def miss():
    """Does a missing call inside the run break the clean-word test or just not count?"""
    W = 4
    print(f"W={W}: 12 informative markers, then `m` markers of the run set to miss/miss")
    for m in (0, 1, 8, 32, 64, 128):
        pos = set((i * 64 * W) // 12 for i in range(12))

        def kinds(t, m=m, pos=pos):
            if t in pos:
                return list(INFO)
            if t >= 64 * W - m:
                return [3, 3, 1, 1, 1, 1]
            return list(FILL)
        got = run(f"g4_m_{m}", W, kinds)
        print(f"  missing={m:>4}  chr2={got:>5}   (full={full(W)})")


if __name__ == "__main__":
    globals()[sys.argv[1] if len(sys.argv) > 1 else "exact"]()
