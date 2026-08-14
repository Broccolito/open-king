"""How many "informative" markers does a clean word need, and over what window?

`gate2.py` establishes the per-marker predicate: a marker only supports a call when
**both** members of the pair carry at least one copy of the A1 allele.  Everything here
counts those markers deterministically -- no allele frequencies are drawn, so a threshold
measured in `k` is a threshold in the statistic itself.

    INFO   = A1A1 / het        (both carry A1; IBS1, so the run stays IBD1)
    FILL   = het / A2A2        (polymorphic, pair carries a het -- but only one member
                                carries A1, so it must not count)
    MONO   = A2A2 / A2A2 with the whole cohort A2A2 (monomorphic control)

    python3 gate3.py perword [W ...]
    python3 gate3.py layout
    python3 gate3.py window
"""

import sys

import gatelab as G

R = G.GateRig(spacing=100_000, n1=640, n2=1600, nsample=6, seed=1)

# dosages of A1: [pair_a, pair_b, o, o, o, o].  Chosen so fixlab's A1-minor
# re-orientation never fires (A1 count <= half of 12).
INFO = [2, 1, 0, 0, 0, 0]      # A1A1 / het      -> both carry A1     (3/12)
FILL = [1, 0, 1, 1, 1, 1]      # het  / A2A2     -> only one carries  (5/12)
MONO = [0, 0, 0, 0, 0, 0]      # A2A2 / A2A2, monomorphic             (0/12)
HH = [1, 1, 0, 0, 0, 0]        # het / het       -> both carry A1     (2/12)


def full(W):
    return 64 * (W + 1) - 1


def run(name, W, kinds, start_word=1):
    f = R.new(name)
    R.block(f, start_word, W, kinds)
    r = R.read(f)
    return r["chr2_mk"] if r else 0


def spread(k, fill=FILL):
    """k informative markers spaced evenly inside each 64-marker word."""
    def kinds(t):
        b = t % 64
        if k <= 0:
            return list(fill)
        return list(INFO) if (b * k) // 64 != ((b - 1) * k) // 64 or b == 0 else list(fill)
    return kinds


def perword(widths=(2, 3, 4, 6, 8, 10, 14)):
    print("k informative markers per word, evenly spread; filler = het/A2A2")
    print(f"{'W':>3} " + " ".join(f"k={k:<2}" for k in range(0, 17)))
    for W in widths:
        cells = []
        for k in range(0, 17):
            got = run(f"g3_pw_{W}_{k}", W, spread(k))
            cells.append("Y" if got else ".")
        print(f"{W:>3} " + " ".join(f"{c:<4}" for c in cells) + f"   full={full(W)}")


def layout():
    """Does placement inside the word matter, or only the count?"""
    W = 4
    print(f"W={W}, full={full(W)}; k informative markers laid out three ways")
    print(f"{'k':>3} {'spread':>8} {'packed-lo':>10} {'packed-hi':>10}")
    for k in range(0, 13):
        a = run(f"g3_l_s_{k}", W, spread(k))

        def packed_lo(t, k=k):
            return list(INFO) if (t % 64) < k else list(FILL)

        def packed_hi(t, k=k):
            return list(INFO) if (t % 64) >= 64 - k else list(FILL)
        b = run(f"g3_l_a_{k}", W, packed_lo)
        c = run(f"g3_l_b_{k}", W, packed_hi)
        print(f"{k:>3} {a:>8} {b:>10} {c:>10}")


def window():
    """Give one word k informative markers and its neighbours none (or vice versa)."""
    W = 6
    print(f"W={W}, full={full(W)}: informative markers confined to some words")
    print(f"{'pattern':>28} {'k':>3} {'chr2':>6}")
    patterns = {
        "all 6 words": [1] * 6,
        "words 0-1 only": [1, 1, 0, 0, 0, 0],
        "words 2-3 only": [0, 0, 1, 1, 0, 0],
        "word 2 only": [0, 0, 1, 0, 0, 0],
        "words 0,2,4": [1, 0, 1, 0, 1, 0],
        "words 0-2": [1, 1, 1, 0, 0, 0],
        "words 0-4": [1, 1, 1, 1, 1, 0],
    }
    for k in (4, 8, 12, 16, 24, 32):
        for label, mask in patterns.items():
            def kinds(t, k=k, mask=mask):
                w, b = t // 64, t % 64
                if not mask[w]:
                    return list(FILL)
                return list(INFO) if (b * k) // 64 != ((b - 1) * k) // 64 or b == 0 \
                    else list(FILL)
            got = run(f"g3_w_{k}_{abs(hash(label)) % 10000}", W, kinds)
            print(f"{label:>28} {k:>3} {got:>6}")
        print()


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "perword"
    if cmd == "perword":
        perword(tuple(int(x) for x in sys.argv[2:]) or (2, 3, 4, 6, 8, 10, 14))
    else:
        globals()[cmd]()
