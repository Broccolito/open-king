"""Which pair genotypes let a clean word be called?  Full vocabulary enumeration.

Every marker inside the block is set explicitly for every sample, and the fileset that
reaches the reference is read back off disk, so the PLINK codes under test are the codes
KING saw -- fixlab's A1-minor re-orientation cannot silently change the question.

    python3 gate2.py vocab
    python3 gate2.py freq
"""

import os
import sys

import gatelab as G

R = G.GateRig(spacing=100_000, n1=640, n2=640, nsample=6, seed=1)

CODE_NAME = {0: "A1A1", 1: "het", 2: "A2A2", 3: "miss"}


def read_back(wd, name, nsample, marker_indices):
    """(codes per sample) for the given global marker indices, straight out of the .bed."""
    bed = os.path.join(wd, name + ".bed")
    raw = open(bed, "rb").read()
    bpv = (nsample + 3) // 4
    body = raw[3:]
    out = {}
    for m in marker_indices:
        row = body[m * bpv:(m + 1) * bpv]
        codes = []
        for s in range(nsample):
            codes.append((row[s >> 2] >> (2 * (s & 3))) & 3)
        # PLINK: 00 = A1A1, 01 = missing, 10 = het, 11 = A2A2  -> dosage of A1
        out[m] = [{0: 0, 1: 3, 2: 1, 3: 2}[c] for c in codes]
    return out


def probe(name, W, vec, start_word=1):
    """Fill a W-word block with one repeated genotype vector; return (chr2_mk, codes)."""
    f = R.new(name)
    R.block(f, start_word, W, lambda t: list(vec))
    r = R.read(f, tag="")
    lo, _ = f.chrom_span(1)
    back = read_back(r["wd"] if r else _last_wd(name), name,
                     R.nsample, [lo + start_word * 64])
    return (r["chr2_mk"] if r else None), back[lo + start_word * 64]


def _last_wd(name):
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "work", name)


def vocab():
    """Sweep the pair's genotypes x the padding cohort's genotypes."""
    others_sets = {
        "o=het":   [1, 1, 1, 1],
        "o=hom0":  [0, 0, 0, 0],
        "o=hom2":  [2, 2, 2, 2],
        "o=mix":   [2, 2, 0, 0],
    }
    pairs = [(0, 0), (1, 1), (2, 2), (1, 0), (0, 1), (2, 1), (1, 2)]
    W = 4
    print(f"W={W}; a full block reports {64 * (W + 1) - 1} intervals, "
          f"IBD2 reports {64 * W - 1}")
    print(f"{'pair':>8} {'others':>8} {'written pair codes':>22} {'maf':>6} {'chr2':>6}")
    for (a, b) in pairs:
        for oname, o in others_sets.items():
            vec = [a, b] + o
            nm = f"g2_{a}{b}_{oname.replace('=', '')}"
            got, codes = probe(nm, W, vec)
            n1 = sum(codes)
            maf = min(n1, 12 - n1) / 12
            print(f"{str((a, b)):>8} {oname:>8} "
                  f"{CODE_NAME[0] if False else str(codes[:2]):>22} {maf:>6.3f} "
                  f"{str(got):>6}")
        print()


def freq():
    """Hold the pair at het/het and slide the padding cohort's allele frequency."""
    W = 4
    print("pair = het/het at every marker in the block; padding cohort varied")
    print(f"{'others':>18} {'A1 count/12':>12} {'chr2':>6}")
    for o in ([0, 0, 0, 0], [1, 0, 0, 0], [1, 1, 0, 0], [1, 1, 1, 0], [1, 1, 1, 1],
              [2, 1, 1, 1], [2, 2, 1, 1]):
        vec = [1, 1] + o
        nm = "g2_f_" + "".join(map(str, o))
        got, codes = probe(nm, W, vec)
        print(f"{str(o):>18} {sum(codes):>12} {str(got):>6}")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "vocab"
    globals()[cmd]()
