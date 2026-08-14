"""Parameterised two-chromosome rig.

chr1 = carrier, fully IBD1, always called whole (verified).
chr2 = test canvas, background forced to opposite homozygotes (IBS0) everywhere
       unless a block is carved out of it.

Both chromosomes are multiples of 64 markers, so chr2 starts on a word boundary
of the global grid and local word w == global word (n1/64 + w).
"""
import fixlab as L


class Rig:
    def __init__(self, spacing, n1=640, n2=640, nsample=6, maf=0.5, seed=1):
        self.sp = spacing
        self.n1, self.n2 = n1, n2
        self.nsample, self.maf, self.seed = nsample, maf, seed
        assert n1 % 64 == 0 and n2 % 64 == 0
        self.carrier = (n1 - 1) * spacing
        self.denom = self.carrier + (n2 - 1) * spacing
        assert self.denom >= 100_000_000, "usable genome must be >= 100 Mb"

    def new(self, name, solid=True, seed=None):
        L.SPACING = self.sp
        f = L.Fixture(name, [(1, self.n1), (2, self.n2)],
                      nsample=self.nsample, maf=self.maf,
                      seed=self.seed if seed is None else seed)
        f.set_state(0, 0, self.n1, L.IBD1)          # carrier
        if solid:
            lo, hi = f.chrom_span(1)
            f.force_ibs0 = set(range(lo, hi))
        f._rig = self
        return f

    def block(self, f, a, b, code=L.IBD1):
        lo, _ = f.chrom_span(1)
        for m in range(lo + a, lo + b):
            f.force_ibs0.discard(m)
        f.set_state(1, a, b, code)

    def poke(self, f, local_idx, kind="ibs0"):
        """Force one marker of chr2 back to a disagreement."""
        lo, _ = f.chrom_span(1)
        (f.force_ibs0 if kind == "ibs0" else f.force_ibs1).add(lo + local_idx)
        if kind != "ibs0":
            f.force_ibs0.discard(lo + local_idx)

    def read(self, f, args=(), tag=""):
        L.SPACING = self.sp
        row, segs, denom, out, wd = L.probe(f, args, tag)
        assert abs(denom - self.denom) < 1000, (denom, self.denom)
        if row is None:
            return None
        b1 = float(row["IBD1Seg"]) * self.denom
        b2 = float(row["IBD2Seg"]) * self.denom
        r = dict(row=row, raw1=float(row["IBD1Seg"]), raw2=float(row["IBD2Seg"]),
                 ibd1_mk=int(round(b1 / self.sp)), ibd2_mk=int(round(b2 / self.sp)),
                 out=out, wd=wd)
        # carrier is IBD1 unless the pair reads as IBD2 there too
        r["test1_mk"] = r["ibd1_mk"] - (self.n1 - 1 if r["ibd2_mk"] < self.n1 // 2 else 0)
        r["test2_mk"] = r["ibd2_mk"]
        return r
