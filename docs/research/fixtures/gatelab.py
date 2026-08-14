"""Exact-genotype rig for the acceptance ("informativeness") gate.

`rig2.py` draws the block's genotypes from allele frequencies, which confounds four
things at once: the per-variant allele frequency across the sample, the pair's
heterozygosity, the pair's HetHet count and the pair's het-vs-hom count.  This rig sets
**every genotype of every sample explicitly**, so those four can be moved one at a time.

Layout (identical to `rig2.py` so results are comparable):

    chr1  = carrier, n1 markers, entirely IBD1  -> the pair always clears the 10 Mb
            pair filter, so the chr2 block's verdict is readable as a length increment
    chr2  = canvas, background forced to opposite homozygotes (IBS0) at every marker;
            a block of `W` complete words is carved out of it starting on a word boundary

Both chromosomes are multiples of 64 markers, so chr2's local word `w` is global word
`n1/64 + w`.

Read-out: `IBD1Seg * denominator` in markers, minus the carrier's `n1 - 1`, is the number
of marker intervals the reference called on chr2.  A run of `W` clean words on a solid
IBS0 background reports exactly `64(W+1) - 1` intervals (10-segment-rule-fixtures §2.2),
so the increment says not just *whether* the block was called but *which words* were.
"""

from __future__ import annotations

import fixlab as L

# genotype vectors are (pair_a, pair_b, other0, other1, ...) counts of the A1 allele.
# fixlab re-orients A1 to the observed minor allele, which maps g -> 2-g for every
# sample at once; that leaves every pair-level IBS category invariant.


class GateRig:
    def __init__(self, spacing=100_000, n1=640, n2=640, nsample=6, seed=1):
        assert n1 % 64 == 0 and n2 % 64 == 0
        self.sp, self.n1, self.n2 = spacing, n1, n2
        self.nsample, self.seed = nsample, seed
        self.carrier = (n1 - 1) * spacing
        self.denom = self.carrier + (n2 - 1) * spacing
        assert self.denom >= 100_000_000, "usable genome must be >= 100 Mb"

    # ---- genotype vocabulary -----------------------------------------
    def vec(self, a, b, others):
        """One marker's genotypes: pair (a, b) plus `others` recycled over the padding."""
        n = self.nsample - 2
        o = [others[k % len(others)] for k in range(n)]
        return [a, b] + o

    def new(self, name, seed=None):
        L.SPACING = self.sp
        f = L.Fixture(name, [(1, self.n1), (2, self.n2)], nsample=self.nsample,
                      maf=0.5, seed=self.seed if seed is None else seed)
        f.set_state(0, 0, self.n1, L.IBD1)              # carrier
        lo, hi = f.chrom_span(1)
        f.force_ibs0 = set(range(lo, hi))               # solid background
        f._rig = self
        return f

    def block(self, f, start_word, W, marker):
        """Carve W complete words out of the background and set every genotype in them.

        `marker(local_index_within_block) -> list of genotypes (len == nsample)`.
        """
        lo, _ = f.chrom_span(1)
        a = start_word * 64
        for t in range(W * 64):
            m = lo + a + t
            f.force_ibs0.discard(m)
            f.pat_all[m] = marker(t)

    def read(self, f, args=(), tag=""):
        L.SPACING = self.sp
        row, segs, denom, out, wd = L.probe(f, args, tag)
        assert abs(denom - self.denom) < 1000, (denom, self.denom)
        if row is None:
            return None
        b1 = float(row["IBD1Seg"]) * self.denom
        b2 = float(row["IBD2Seg"]) * self.denom
        mk1 = int(round(b1 / self.sp))
        mk2 = int(round(b2 / self.sp))
        return dict(row=row, ibd1_mk=mk1, ibd2_mk=mk2,
                    chr2_mk=mk1 + mk2 - (self.n1 - 1), out=out, wd=wd)
