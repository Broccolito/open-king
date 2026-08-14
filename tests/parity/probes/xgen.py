#!/usr/bin/env python3
"""Build X-bearing PLINK filesets with a controlled pedigree.

Support code for `xseg_probe.py`; `generate_corpus.py` builds the committed corpus and
this builds throwaway filesets to hold rules out against. Haplotype simulation with
recombination and correct X transmission:

  * a mother transmits one recombined X to every child;
  * a father transmits his single X unchanged to daughters, and nothing to sons;
  * a male is hemizygous on X and is written **homozygous**, which is both what PLINK
    does and the representational fact the reference's X segment caller relies on.

`make_map` lays out an autosomal map plus an X arm of a chosen marker count and span;
`build(out, people, rows, seed)` writes `out.{bed,bim,fam}`. `people` is a list of
`(fid, iid, pat, mat, sex)` with parents before children.

Nothing here is committed as a fixture: every caller regenerates from a seed.
"""
import random
import sys

MISS = None


class Spec:
    def __init__(self, name, seed=1):
        self.name = name
        self.rnd = random.Random(seed)
        self.people = []          # (fid, iid, pat, mat, sex)
        self.index = {}


def make_map(auto_chroms, auto_n, auto_span, x_n, x_span, x_chr=23,
             extra=()):
    """Return list of (chrom, snp_id, bp)."""
    rows = []
    for c in auto_chroms:
        for m in range(auto_n):
            bp = 1_000_000 + m * (auto_span // max(auto_n - 1, 1))
            rows.append((str(c), f"rs{c}_{bp}", bp))
    for m in range(x_n):
        bp = 1_000_000 + m * (x_span // max(x_n - 1, 1))
        rows.append((str(x_chr), f"rsX_{bp}", bp))
    for (c, n, span) in extra:
        for m in range(n):
            bp = 1_000_000 + m * (span // max(n - 1, 1))
            rows.append((str(c), f"rs{c}_{bp}", bp))
    return rows


class Sim:
    """Haplotype simulator over one map."""

    def __init__(self, rows, seed):
        self.rows = rows
        self.rnd = random.Random(seed)
        self.freq = [self.rnd.uniform(0.15, 0.5) for _ in rows]
        # per-chrom marker index lists
        self.by_chr = {}
        for i, (c, _, _) in enumerate(rows):
            self.by_chr.setdefault(c, []).append(i)

    def founder_hap(self):
        return [1 if self.rnd.random() < p else 0 for p in self.freq]

    def meiosis(self, h1, h2, chrom=None):
        """Recombine two haplotypes; restrict to `chrom` if given."""
        out = list(h1)
        idxs = self.by_chr[chrom] if chrom else range(len(self.rows))
        cur = self.rnd.random() < 0.5
        # one crossover per ~50Mb
        last_chr = None
        for i in idxs:
            c, _, bp = self.rows[i]
            if c != last_chr:
                cur = self.rnd.random() < 0.5
                last_chr = c
                next_x = bp + self.rnd.expovariate(1 / 50_000_000)
            if bp > next_x:
                cur = not cur
                next_x = bp + self.rnd.expovariate(1 / 50_000_000)
            out[i] = (h2 if cur else h1)[i]
        return out


def build(out, people, rows, seed=7, x_chr="23", missing=()):
    """people: list of (fid, iid, pat, mat, sex).  Parents must precede children
    unless they are '0' (founder)."""
    sim = Sim(rows, seed)
    n = len(rows)
    xidx = set(i for i, (c, _, _) in enumerate(rows) if c == x_chr)

    hap = {}   # iid -> (h_pat, h_mat) autosomal+X; for males h_pat is None on X
    sexof = {}
    order = []
    for (fid, iid, pat, mat, sex) in people:
        sexof[iid] = sex
        order.append(iid)

    for (fid, iid, pat, mat, sex) in people:
        if pat in hap and mat in hap:
            fp, fm = hap[pat]
            mp, mm = hap[mat]
            from_dad = sim.meiosis(fp, fm)
            from_mom = sim.meiosis(mp, mm)
            # X: mother always transmits a recombined X
            xmom = sim.meiosis(mp, mm, x_chr)
            for i in xidx:
                from_mom[i] = xmom[i]
            if sex == 1:
                # no paternal X
                for i in xidx:
                    from_dad[i] = None
            else:
                # father's single X, unchanged (fp is his hemizygous hap)
                for i in xidx:
                    from_dad[i] = fp[i]
            hap[iid] = (from_dad, from_mom)
        else:
            h1 = sim.founder_hap()
            h2 = sim.founder_hap()
            if sex == 1:
                for i in xidx:
                    h2[i] = None      # hemizygous: only h1 carries the X
            hap[iid] = (h1, h2)

    # genotype: count of allele 1
    geno = {}
    for iid in order:
        h1, h2 = hap[iid]
        g = []
        for i in range(n):
            a, b = h1[i], h2[i]
            if a is None and b is None:
                g.append(MISS)
            elif a is None:
                g.append(2 * b)
            elif b is None:
                g.append(2 * a)
            else:
                g.append(a + b)
        geno[iid] = g

    for (iid, frac) in missing:
        rnd = random.Random(hash(iid) & 0xffff)
        for i in range(n):
            if rnd.random() < frac:
                geno[iid][i] = MISS

    # orient A1 as minor allele
    ns = len(order)
    for i in range(n):
        cnt = sum(geno[s][i] for s in order if geno[s][i] is not None)
        tot = sum(2 for s in order if geno[s][i] is not None)
        if tot and cnt * 2 > tot:
            for s in order:
                if geno[s][i] is not None:
                    geno[s][i] = 2 - geno[s][i]

    with open(out + ".fam", "w") as f:
        for (fid, iid, pat, mat, sex) in people:
            f.write(f"{fid} {iid} {pat} {mat} {sex} -9\n")
    with open(out + ".bim", "w") as f:
        for (c, sid, bp) in rows:
            f.write(f"{c}\t{sid}\t{bp / 1e6:.6f}\t{bp}\tA\tG\n")
    code = {2: 0b00, 1: 0b10, 0: 0b11, None: 0b01}
    with open(out + ".bed", "wb") as f:
        f.write(bytes([0x6C, 0x1B, 0x01]))
        nb = (ns + 3) // 4
        for i in range(n):
            buf = bytearray(nb)
            for s, iid in enumerate(order):
                buf[s // 4] |= code[geno[iid][i]] << (2 * (s % 4))
            f.write(bytes(buf))
