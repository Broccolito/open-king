"""Sanity + first discrimination sweep for the acceptance gate.

    python3 gate1.py sanity
    python3 gate1.py sweep [W ...]
"""

import sys

import gatelab as G

R = GR = G.GateRig(spacing=100_000, n1=640, n2=640, nsample=6, seed=1)

# marker vocabulary: (pair_a, pair_b, others...) as A1 dosages
PAT = {
    # pair het/het, marker polymorphic
    "HH":   [1, 1, 1, 1, 1, 1],
    # pair het vs hom (IBS1), marker polymorphic
    "H0":   [1, 0, 1, 1, 1, 1],
    # pair hom-concordant on the major allele, marker polymorphic in the others
    "CC0":  [0, 0, 2, 2, 1, 1],
    # pair hom-concordant on the minor allele, marker polymorphic
    "CC2":  [2, 2, 1, 1, 0, 0],
    # everybody homozygous for the same allele: a monomorphic marker
    "MONO": [0, 0, 0, 0, 0, 0],
    # pair het/het but nobody else carries the allele: high pair evidence, tiny MAF
    "HHmono": [1, 1, 0, 0, 0, 0],
}


def expected(W):
    """Marker intervals a fully-called W-word block reports on the solid background."""
    return 64 * (W + 1) - 1


def run(name, W, kinds, start_word=1, tag=""):
    """kinds: callable t -> key of PAT, for t in [0, 64W)."""
    f = R.new(name)
    R.block(f, start_word, W, lambda t: PAT[kinds(t)])
    r = R.read(f, tag=tag)
    return (r["chr2_mk"] if r else None), r


def sanity():
    for W in (2, 4):
        for k in ("HH", "H0", "CC0", "CC2", "MONO", "HHmono"):
            got, _ = run(f"g1_{k}_{W}", W, lambda t, k=k: k)
            print(f"W={W} all-{k:<7} chr2={got!s:>5}  (full block = {expected(W)})")
        print()


def sweep(widths=(2, 3, 4, 6, 8)):
    """n strong markers per word (spread evenly), rest monomorphic."""
    print("n = HetHet markers per word, evenly spread; fillers monomorphic")
    print(f"{'W':>3} {'n':>3} {'total':>5} {'chr2_mk':>8} {'full':>5}")
    for W in widths:
        for n in range(0, 33):
            def kinds(t, n=n):
                return "HH" if (t % 64) * n // 64 != ((t % 64) - 1) * n // 64 or (
                    t % 64) == 0 and n > 0 else "MONO"
            got, _ = run(f"g1_sw_{W}_{n}", W, kinds)
            print(f"{W:>3} {n:>3} {n * W:>5} {got!s:>8} {expected(W):>5}")
            if got and got >= expected(W) - 1:
                break


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "sanity"
    if cmd == "sanity":
        sanity()
    else:
        sweep(tuple(int(x) for x in sys.argv[2:]) or (2, 3, 4, 6, 8))
