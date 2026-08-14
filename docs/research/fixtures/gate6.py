"""The IBD2 side of the gate, and the classification of hom-concordant runs.

    python3 gate6.py fine
    python3 gate6.py mix
    python3 gate6.py homrun
"""

import sys

import gatelab as G

R = G.GateRig(spacing=100_000, n1=640, n2=1600, nsample=6, seed=1)

HET2 = [1, 1, 0, 0, 0, 0]      # het / het        -> HetHet
A1A1 = [2, 2, 1, 0, 0, 0]      # A1A1 / A1A1      -> one member hom for A1
A1HE = [2, 1, 0, 0, 0, 0]      # A1A1 / het       -> one member hom for A1, but IBS1
HOM0 = [0, 0, 0, 0, 0, 0]      # A2A2 / A2A2, monomorphic; neither IBS0 nor IBS1
FILL = [1, 0, 1, 1, 1, 1]      # het / A2A2, an IBS1 that forbids IBD2


def run(name, W, kinds, start_word=1):
    f = R.new(name)
    R.block(f, start_word, W, kinds)
    r = R.read(f)
    if not r:
        return 0, 0
    return r["ibd1_mk"] - (R.n1 - 1) if r["ibd2_mk"] < R.n1 // 2 else r["ibd1_mk"], \
        r["ibd2_mk"]


def mixplace(W, spec, fill):
    """spec: list of (vector, count).  Markers are interleaved deterministically."""
    slots = {}
    total = sum(c for _, c in spec)
    idx = 0
    for vec, c in spec:
        for i in range(c):
            slots[((idx + i * len(spec)) * 64 * W) // max(total, 1) % (64 * W)] = vec
        idx += 1
    # re-lay them on a clean stride to avoid collisions
    slots = {}
    t = 0
    for vec, c in spec:
        for _ in range(c):
            while t in slots:
                t += 1
            slots[t] = vec
            t += max(1, (64 * W) // max(total, 1))
    return lambda x: list(slots.get(x, fill))


def fine():
    print("IBD2-eligible run (filler A2A2/A2A2 so IBS1 == 0); evidence markers swept.")
    for label, vec in (("het/het", HET2), ("A1A1/A1A1", A1A1)):
        print(f"\n  evidence = {label}")
        print(f"  {'W':>3} " + " ".join(f"k={k:<3}" for k in range(6, 14)))
        for W in (1, 2, 4):
            cells = []
            for k in range(6, 14):
                a, b = run(f"g6_f_{label[0]}{W}_{k}", W, mixplace(W, [(vec, k)], HOM0))
                cells.append(f"{a}/{b}")
            print(f"  {W:>3} " + " ".join(f"{c:<5}" for c in cells)
                  + f"   ibd1 full={64 * (W + 1) - 1}, ibd2 full={64 * W - 1}")


def mix():
    """Do HetHet and A1A1 markers pool into one count, or are they separate?"""
    W = 2
    print(f"W={W}, IBD2-eligible filler: a HetHet markers + b A1A1/A1A1 markers")
    print(f"{'a(HetHet)':>10} {'b(A1A1)':>9} {'a+b':>5} {'ibd1':>6} {'ibd2':>6}")
    for a, b in ((0, 9), (0, 10), (9, 0), (10, 0), (5, 5), (5, 4), (4, 5), (6, 4),
                 (9, 1), (1, 9), (12, 0), (0, 12)):
        i1, i2 = run(f"g6_m_{a}_{b}", W, mixplace(W, [(HET2, a), (A1A1, b)], HOM0))
        print(f"{a:>10} {b:>9} {a + b:>5} {i1:>6} {i2:>6}")


def homrun():
    """A run with no heterozygote at all: IBD1 or IBD2?"""
    W = 4
    for label, fill, evid in (("A1A1/A1A1 everywhere", A1A1, None),
                              ("A2A2/A2A2 everywhere", HOM0, None),
                              ("A1A1/het everywhere (IBS1)", A1HE, None)):
        i1, i2 = run(f"g6_h_{abs(hash(label)) % 99999}", W, lambda t, f=fill: list(f))
        print(f"{label:>30}  ibd1={i1:>5} ibd2={i2:>5}   "
              f"(ibd1 full={64 * (W + 1) - 1}, ibd2 full={64 * W - 1})")


if __name__ == "__main__":
    globals()[sys.argv[1] if len(sys.argv) > 1 else "fine"]()
