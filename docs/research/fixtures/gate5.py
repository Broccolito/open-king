"""Is the gate keyed on the .bim's A1 column, or on the observed minor allele?

`gate4.py` shows a marker only counts when at least one member of the pair is homozygous
for the **A1** allele -- but every fixture so far re-orients A1 to the observed minor
allele, so "A1A1" and "hom for the minor allele" have never been separated.  Here the
block's markers are written with the re-orientation switched off, so the two come apart:

    A1major   pair = A2A2 / A2A2, cohort = A1A1  -> pair is hom for the MINOR allele,
                                                    but A1A1 is false
    A1minor   pair = A1A1 / het,  cohort = A2A2  -> the gate4 control

Also measures the IBD2 side, which `gate2.py` shows is *not* subject to this gate.

    python3 gate5.py allele
    python3 gate5.py ibd2
"""

import sys

import gatelab as G

R = G.GateRig(spacing=100_000, n1=640, n2=1600, nsample=6, seed=1)

FILL = [1, 0, 1, 1, 1, 1]      # het / A2A2, weight 0
HOM0 = [0, 0, 0, 0, 0, 0]      # A2A2 / A2A2 everywhere: monomorphic, IBS1-free


def run(name, W, kinds, noflip_all=False, start_word=1):
    f = R.new(name)
    R.block(f, start_word, W, kinds)
    if noflip_all:
        lo, _ = f.chrom_span(1)
        f.noflip = set(range(lo + start_word * 64, lo + (start_word + W) * 64))
    r = R.read(f)
    return r["chr2_mk"] if r else 0


def place(W, n, vec, fill=FILL):
    pos = set((i * 64 * W) // n for i in range(n)) if n else set()
    return lambda t: list(vec) if t in pos else list(fill)


def allele():
    W = 4
    cands = {
        # written verbatim (no re-orientation): A1 is the MAJOR allele here, and the
        # pair is homozygous for the minor one
        ("A1major: pair A2A2/A2A2, cohort A1A1", True): [0, 0, 2, 2, 2, 2],
        ("A1major: pair A2A2/het,  cohort A1A1", True): [0, 1, 2, 2, 2, 2],
        # controls, same fixture path, A1 left as the minor allele
        ("A1minor: pair A1A1/het,  cohort A2A2", True): [2, 1, 0, 0, 0, 0],
        ("A1minor: pair A1A1/A1A1, cohort A2A2", True): [2, 2, 0, 0, 0, 0],
    }
    print(f"W={W}: n candidate markers among het/A2A2 filler; re-orientation OFF.")
    print("A weight-1 marker first passes at n=10.")
    print(f"{'candidate':>40} " + " ".join(f"{n:>4}" for n in range(8, 13)))
    for (label, nf), vec in cands.items():
        cells = []
        for n in range(8, 13):
            got = run(f"g5_a_{abs(hash(label)) % 100000}_{n}", W, place(W, n, vec),
                      noflip_all=nf)
            cells.append("Y" if got else ".")
        print(f"{label:>40} " + " ".join(f"{c:>4}" for c in cells))


def ibd2():
    """An IBD2-eligible run (no het-vs-hom anywhere): how much evidence does it need?"""
    print("IBD2 run: k het/het markers, filler A2A2/A2A2 (so IBS1 == 0 everywhere).")
    print(f"{'W':>3} " + " ".join(f"k={k:<2}" for k in
                                  (0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32)))
    for W in (1, 2, 4):
        cells = []
        for k in (0, 1, 2, 3, 4, 6, 8, 10, 12, 16, 24, 32):
            got = run(f"g5_2_{W}_{k}", W, place(W, k, [1, 1, 0, 0, 0, 0], fill=HOM0))
            cells.append(str(got))
        print(f"{W:>3} " + " ".join(f"{c:<4}" for c in cells)
              + f"   ibd2 full={64 * W - 1}, ibd1 full={64 * (W + 1) - 1}")


if __name__ == "__main__":
    globals()[sys.argv[1] if len(sys.argv) > 1 else "allele"]()
