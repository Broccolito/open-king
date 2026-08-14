#!/usr/bin/env python3
"""Build a PLINK fileset of independent pairs with prescribed IBD1/IBD2 fractions.

Usage: mkpairs.py OUT targets.txt   where targets.txt has "f1 f2" per line.
Each target becomes one 2-person family (FID Pk, IIDs Pk_A / Pk_B).
"""
import random

import sys

L_BP = 250_000_000
SPACING = 50_000
NSNP = L_BP // SPACING  # 5000
SEED = 20260813


def main():
    out = sys.argv[1]
    targets = []
    for line in open(sys.argv[2]):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        p = line.split()
        targets.append((float(p[0]), float(p[1]),
                        int(p[2]) if len(p) > 2 else 0))
    rnd = random.Random(SEED)

    freq = [rnd.uniform(0.10, 0.50) for _ in range(NSNP)]
    npair = len(targets)
    nsamp = 2 * npair
    # geno[s][m] in 0,1,2 = count of allele "1"
    geno = [[0] * NSNP for _ in range(nsamp)]

    for k, (f1, f2, off) in enumerate(targets):
        n2 = int(round(f2 * NSNP))
        n1 = int(round(f1 * NSNP))
        a = 2 * k
        b = 2 * k + 1
        for m0 in range(NSNP):
            m = m0
            p = freq[m]
            mm = (m0 - off) % NSNP
            if mm < n2:           # IBD2: identical genotypes
                g = draw2(rnd, p)
                geno[a][m] = g
                geno[b][m] = g
            elif mm < n2 + n1:    # IBD1: one shared haplotype
                hs = 1 if rnd.random() < p else 0
                ha = 1 if rnd.random() < p else 0
                hb = 1 if rnd.random() < p else 0
                geno[a][m] = hs + ha
                geno[b][m] = hs + hb
            else:                 # IBD0
                geno[a][m] = draw2(rnd, p)
                geno[b][m] = draw2(rnd, p)

    # orient A1 as the observed minor allele (what plink --make-bed does)
    for m in range(NSNP):
        cnt1 = sum(geno[s][m] for s in range(nsamp))
        if cnt1 > nsamp:  # allele "1" is the major allele -> swap
            for s in range(nsamp):
                geno[s][m] = 2 - geno[s][m]

    with open(out + ".fam", "w") as f:
        for k in range(npair):
            f.write(f"P{k:03d} P{k:03d}_A 0 0 1 -9\n")
            f.write(f"P{k:03d} P{k:03d}_B 0 0 2 -9\n")
    with open(out + ".bim", "w") as f:
        for m in range(NSNP):
            bp = (m + 1) * SPACING
            f.write(f"1\trs{m}\t{bp / 1e6:.6f}\t{bp}\tA\tG\n")
    # .bed SNP-major.  bit code: 0=hom A1(2 copies of allele1), 2=het, 3=hom A2
    code = {2: 0b00, 1: 0b10, 0: 0b11}
    with open(out + ".bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        nb = (nsamp + 3) // 4
        for m in range(NSNP):
            buf = bytearray(nb)
            for s in range(nsamp):
                buf[s >> 2] |= code[geno[s][m]] << (2 * (s & 3))
            f.write(bytes(buf))



def draw2(rnd, p):
    return (1 if rnd.random() < p else 0) + (1 if rnd.random() < p else 0)


main()
