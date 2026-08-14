#!/usr/bin/env python3
"""Controlled IBD fixture lab for open-king.

Builds PLINK1 binary filesets in which the IBD state of ONE designated pair is
exact by construction (explicit shared haplotypes, no pedigree transmission),
runs the KING 2.3.2 reference binary on them, and parses back .seg/allsegs.txt.

Sample 0 = "A", sample 1 = "B" are the test pair. Samples 2..n-1 are unrelated
padding, present only because --ibdseg downgrades to --kinship below 5 samples.
"""

from __future__ import annotations

import os
import random
import shutil
import subprocess
import struct
import sys

# The binary every fixture in this directory drives. Defaults to the KING 2.3.2
# reference; set `$KING` to point a fixture family at *our* build instead, which is how
# a rule fitted against the reference is then checked against the Rust port.
KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
)
BED_MAGIC = bytes((0x6C, 0x1B, 0x01))
# index by A1 dosage 0,1,2 ; 3 == missing
BED_CODE = (0b11, 0b10, 0b00, 0b01)

SPACING = 10_000  # bp between consecutive markers, constant everywhere

IBD0, IBD1, IBD2 = 0, 1, 2


class Fixture:
    """chroms: list of (chrom_code, n_markers). state: per-global-marker IBD code."""

    def __init__(self, name, chroms, nsample=6, maf=0.5, seed=7):
        self.name = name
        self.chroms = chroms
        self.nsample = nsample
        self.maf = maf
        self.seed = seed
        self.nm = sum(n for _, n in chroms)
        self.state = [IBD0] * self.nm
        self.force_ibs0 = set()   # global marker idx: force A=hom A1A1, B=hom A2A2
        self.force_ibs1 = set()   # global marker idx: force A=het, B=hom  (het/hom)
        self.missing = set()      # global marker idx: both of the pair missing
        self.maf_of = None        # optional per-marker maf override list
        self.pat = {}             # global marker idx -> (g0, g1) exact pair genotypes
        self.pat_all = {}         # global marker idx -> full list of genotypes, all samples
        self.noflip = set()       # markers exempt from the A1-minor re-orientation, so a
                                  # fixture can put the MAJOR allele in the A1 column

    # ---- span helpers -------------------------------------------------
    def chrom_span(self, k):
        """(lo, hi) global marker indices of the k-th chromosome, hi exclusive."""
        lo = sum(n for _, n in self.chroms[:k])
        return lo, lo + self.chroms[k][1]

    def set_state(self, k, a, b, code):
        """Set state over local markers [a, b) of chromosome k."""
        lo, _ = self.chrom_span(k)
        for m in range(lo + a, lo + b):
            self.state[m] = code

    # ---- generation ---------------------------------------------------
    def build(self, outdir):
        os.makedirs(outdir, exist_ok=True)
        rng = random.Random(self.seed)
        nm, ns = self.nm, self.nsample
        maf = self.maf_of or [self.maf] * nm

        # positions / chromosome codes
        chrom = []
        pos = []
        for code, n in self.chroms:
            for i in range(n):
                chrom.append(code)
                pos.append((i + 1) * SPACING)

        # haplotypes: hap[s][0], hap[s][1] as lists of A1-allele indicators
        # generated marker by marker so the shared-haplotype rule is local.
        geno = [[0] * nm for _ in range(ns)]
        for m in range(nm):
            p = maf[m]
            h = [[rng.random() < p, rng.random() < p] for _ in range(ns)]
            st = self.state[m]
            if st >= IBD1:
                h[1][0] = h[0][0]
            if st == IBD2:
                h[1][1] = h[0][1]
            for s in range(ns):
                geno[s][m] = int(h[s][0]) + int(h[s][1])

        # forced overrides on the test pair
        # alternate the polarity so the forced background does not skew allele
        # frequencies (KING's A1-major fatal samples markers with an unseeded RNG)
        for m in self.force_ibs0:
            geno[0][m] = 2 if (m & 1) == 0 else 0
            geno[1][m] = 0 if (m & 1) == 0 else 2
        for m in self.force_ibs1:
            geno[0][m] = 1
            geno[1][m] = 0
        for m in self.missing:
            geno[0][m] = 3
            geno[1][m] = 3
        for m, (a, b) in self.pat.items():
            geno[0][m] = a
            geno[1][m] = b
        for m, gl in self.pat_all.items():
            for s_, g_ in enumerate(gl):
                geno[s_][m] = g_

        # re-orient so A1 is the observed minor allele (KING's --make-bed gate)
        flipped = [False] * nm
        for m in range(nm):
            a1 = sum(geno[s][m] for s in range(ns) if geno[s][m] != 3)
            called = sum(1 for s in range(ns) if geno[s][m] != 3)
            if a1 * 2 > called * 2:  # a1 count > a2 count
                pass
            if a1 > (2 * called - a1) and m not in self.noflip:
                flipped[m] = True
                for s in range(ns):
                    if geno[s][m] != 3:
                        geno[s][m] = 2 - geno[s][m]

        prefix = os.path.join(outdir, self.name)
        with open(prefix + ".fam", "w") as f:
            for s in range(ns):
                iid = "S%02d" % s
                f.write("F%02d\t%s\t0\t0\t1\t-9\n" % (s, iid))
        with open(prefix + ".bim", "w") as f:
            for m in range(nm):
                a1, a2 = ("A", "G") if not flipped[m] else ("G", "A")
                f.write("%d\trs%d\t%.6f\t%d\t%s\t%s\n"
                        % (chrom[m], m, pos[m] / 1e6, pos[m], a1, a2))
        with open(prefix + ".bed", "wb") as f:
            f.write(BED_MAGIC)
            per = (ns + 3) // 4
            for m in range(nm):
                buf = bytearray(per)
                for s in range(ns):
                    buf[s >> 2] |= BED_CODE[geno[s][m]] << (2 * (s & 3))
                f.write(bytes(buf))
        self.chrom_arr, self.pos_arr, self.geno = chrom, pos, geno
        return prefix


def run_king(prefix, args, workdir):
    cmd = [KING, "-b", prefix + ".bed"] + args + ["--prefix", os.path.join(workdir, "k")]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=workdir)
    return r


def parse_seg(path):
    if not os.path.exists(path):
        return {}
    rows = {}
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            v = line.rstrip("\n").split("\t")
            rows[(v[1], v[3])] = dict(zip(hdr, v))
    return rows


def parse_allsegs(path):
    if not os.path.exists(path):
        return []
    out = []
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            out.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    return out


def probe(fix, args=(), tag=""):
    """Build, run, and return (segrow_for_pair, allsegs, denom_bp, stdout)."""
    root = os.path.dirname(os.path.abspath(__file__))
    wd = os.path.join(root, "work", fix.name + tag)
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    prefix = fix.build(wd)
    r = run_king(prefix, list(args) + ["--ibdseg"], wd)
    for _ in range(6):
        if "FATAL ERROR" not in r.stdout:
            break
        r = run_king(prefix, list(args) + ["--ibdseg"], wd)   # unseeded A1-major gate
    if "FATAL ERROR" in r.stdout:
        raise RuntimeError("KING fatal on %s: %s" % (prefix,
                           [l for l in r.stdout.splitlines() if "FATAL" in l]))
    segs = parse_allsegs(os.path.join(wd, "kallsegs.txt"))
    denom = sum(float(s["Length"]) for s in segs if int(s["Chr"]) < 23) * 1e6
    rows = parse_seg(os.path.join(wd, "k.seg"))
    return rows.get(("S00", "S01")), segs, denom, r.stdout, wd
