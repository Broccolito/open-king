#!/usr/bin/env python3
"""Regenerate the probes behind `--related`'s `Error` rule and its `Dup/MZ` gate.

Both rules are invisible in the golden corpus and were fitted against the reference
binary on filesets built here.  Each probe is a set of independent families; the
*target* pair of every family is one of `mkpairs.py`'s prescribed-IBD pairs, and the
rest of the family exists only to give that pair a declared relationship.

    python3 pederr.py error   OUT   # declared PO/FS/2nd/3rd/4th/UN x prescribed IBD
    python3 pederr.py dupgate OUT   # an IBD2 ladder across the HetConc = 0.8 boundary
    python3 pederr.py homibs0 OUT   # pairs whose HomIBS0 lands on an exact 4-dp tie

then run the reference and read the `.kin`:

    king -b OUT.bed --related --degree 4 --prefix OUT

`error` pins `Error` (see `related::error_flag`): `InfType` in {2nd, 3rd, 4th} is graded
by `kinship::error_flag(PropIBD / 2, Phi)`, every other label by exact agreement.
`dupgate` pins the `Dup/MZ` clause of `.kin`'s `InfType` (see `PairIbd::inf_type`): it
needs `HetConc > 0.8` on top of `IBD2Seg > 0.7`, which `.seg` does not.

Python 3 standard library only; writes `<OUT>.txt` and shells out to `mkpairs.py`.
"""

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


# --- family shapes -------------------------------------------------------
# Each takes the family's IIDs in .fam order -- which must stay mkpairs' own order,
# P<k>_A, P<k>_B, P<k+1>_A, ... -- and returns (iid, father, mother, sex) rows.
# The target pair is always the first two, and the shape fixes its (Z0, Phi).

def unrel(n):                                    # Z0 1.000  Phi 0.0000
    x, y = n
    return [(x, "0", "0", "1"), (y, "0", "0", "2")]


def po(n):                                       # Z0 0.000  Phi 0.2500
    x, y, m, z = n
    return [(x, "0", "0", "1"), (y, x, m, "2"),
            (m, "0", "0", "2"), (z, "0", "0", "1")]


def fs(n):                                       # Z0 0.250  Phi 0.2500
    x, y, f, m = n
    return [(x, f, m, "1"), (y, f, m, "2"),
            (f, "0", "0", "1"), (m, "0", "0", "2")]


def halfsib(n):                                  # Z0 0.500  Phi 0.1250
    x, y, p, m1, m2, z = n
    return [(x, p, m1, "1"), (y, p, m2, "2"), (p, "0", "0", "1"),
            (m1, "0", "0", "2"), (m2, "0", "0", "2"), (z, "0", "0", "1")]


def cousin(n):                                   # Z0 0.750  Phi 0.0625
    x, y, g1, g2, a, b, sa, sb = n
    return [(x, a, sa, "1"), (y, sb, b, "2"), (g1, "0", "0", "1"), (g2, "0", "0", "2"),
            (a, g1, g2, "1"), (b, g1, g2, "2"), (sa, "0", "0", "2"), (sb, "0", "0", "1")]


def c1r(n):                                      # Z0 0.875  Phi 0.0312
    x, y, g1, g2, a, b, sa, sb, d, sd = n
    return [(x, a, sa, "1"), (y, d, sd, "2"), (g1, "0", "0", "1"), (g2, "0", "0", "2"),
            (a, g1, g2, "1"), (b, g1, g2, "2"), (sa, "0", "0", "2"), (sb, "0", "0", "1"),
            (d, sb, b, "1"), (sd, "0", "0", "2")]


SHAPES = {"unrel": (1, unrel), "po": (2, po), "fs": (2, fs),
          "halfsib": (3, halfsib), "cousin": (4, cousin), "c1r": (5, c1r)}


# --- the two probes ------------------------------------------------------
# (shape, IBD1 fraction, IBD2 fraction) per family.  PropIBD ~= f2 + f1/2, so the
# fractions below place each target pair in a chosen InfType band; the realized values
# are read back from the reference's own output and need only be roughly right.

ERROR_SPEC = (
    # declared unrelated against inferred UN / 4th / 3rd / 2nd
    [("unrel", f1, 0.0) for f1 in (0.05, 0.08, 0.10, 0.125, 0.15, 0.17, 0.25, 0.35)]
    # declared 4th (first cousin once removed) right across the 4th/UN boundary
    + [("c1r", f1, 0.0) for f1 in (0.062, 0.064, 0.066, 0.070, 0.074, 0.076,
                                   0.080, 0.084, 0.086, 0.090, 0.10, 0.13,
                                   0.16, 0.25, 0.50)]
    # declared 3rd and 2nd, sweeping the graded band and past it into PO/FS
    + [("cousin", f1, f2) for f1, f2 in ((0.04, 0.0), (0.12, 0.0), (0.24, 0.0),
                                         (0.30, 0.0), (0.36, 0.0), (0.50, 0.0),
                                         (0.90, 0.0), (0.0, 0.85))]
    + [("halfsib", f1, f2) for f1, f2 in ((0.04, 0.0), (0.12, 0.0), (0.24, 0.0),
                                          (0.50, 0.0), (0.82, 0.0), (0.90, 0.0),
                                          (0.98, 0.0), (0.70, 0.10), (0.50, 0.25),
                                          (0.0, 0.85))]
    # declared PO and FS against every close label, including the ratio band where a
    # multiplicative rule alone would say 0.5 and the reference says 0 or 1
    + [(shape, f1, f2) for shape in ("po", "fs")
       for f1, f2 in ((0.98, 0.0), (0.90, 0.0), (0.70, 0.10), (0.60, 0.20),
                      (0.50, 0.25), (0.50, 0.50), (0.40, 0.55), (0.30, 0.35),
                      (0.15, 0.75), (0.02, 0.80), (0.0, 0.85))]
)

# One two-person family per point, IBD2 stepped finely through the gate.  Four offsets
# per fraction spread the block over the map so the realized HetConc varies within a
# step; the boundary lands between 0.7986 (FS) and 0.8004 (Dup/MZ).
DUPGATE_SPEC = [("unrel", 0.0, f2 / 1000.0)
                for f2 in range(760, 832, 4) for _ in range(4)]
DUPGATE_OFFSETS = [0, 500, 1000, 1500]


def build(spec, out, offsets=None):
    start, targets, fams = 0, {}, []
    for q, (shape, f1, f2) in enumerate(spec):
        npairs, make = SHAPES[shape]
        off = 0 if offsets is None else offsets[q % len(offsets)]
        targets[start] = (f1, f2, off)
        names = []
        for k in range(start, start + npairs):
            names += [f"P{k:03d}_A", f"P{k:03d}_B"]
        fams.append((f"{shape[:2].upper()}{q:02d}", make(names)))
        start += npairs

    Path(out + ".txt").write_text("\n".join(
        "%.6f %.6f %d" % targets.get(k, (0.0, 0.0, 0))
        for k in range(start)) + "\n")
    subprocess.run([sys.executable, str(HERE / "mkpairs.py"), out, out + ".txt"],
                   check=True)
    # mkpairs writes its own one-family-per-pair .fam; replace it, keeping row order.
    Path(out + ".fam").write_text("".join(
        f"{fid} {iid} {fa} {mo} {sex} -9\n"
        for fid, members in fams for iid, fa, mo, sex in members))
    print(f"{out}: {start} pairs, {2 * start} samples, {len(fams)} families")


# --- the HomIBS0 tie probe ----------------------------------------------
# `HomIBS0 = N_IBS0 / |{i hom-A1} u {j hom-A1}|`, so a pair with exactly `a` opposite
# homozygotes and exactly `b` sites where either is hom-A1 prints `a/b`.  Every (a, b)
# here makes `a/b` an exact four-decimal tie, which is the only place the reference and
# a plain double division disagree.  Genotypes are placed by hand rather than drawn:
# markers [o, o+a) are hom-A1 / hom-A2, [o+a, o+b) are hom-A1 / hom-A1, and everything
# else is hom-A2 / hom-A2 -- so A1 stays the minor allele fileset-wide and the
# `Too many first alleles as the major allele` gate never fires.  Each pair also gets a
# private 200-marker block where both members are heterozygous, without which
# `Het1 + Het2` is zero and the row is dropped.
HOMIBS0_TARGETS = [(3, 32), (31, 32), (15, 96), (51, 96), (9, 160),
                   (17, 160), (13, 160), (7, 160), (414, 960)]
HOMIBS0_M = 5000
HOM_A1, HET, HOM_A2 = 0, 2, 3


def build_homibs0(out, targets=HOMIBS0_TARGETS, m=HOMIBS0_M):
    n = 2 * len(targets)
    geno = [[HOM_A2] * m for _ in range(n)]
    off = 0
    for f, (a, b) in enumerate(targets):
        lo, hi = 2 * f, 2 * f + 1
        for k in range(off, off + a):
            geno[lo][k] = HOM_A1
        for k in range(off + a, off + b):
            geno[lo][k] = HOM_A1
            geno[hi][k] = HOM_A1
        off += b
    het_at = off + 64
    assert het_at + 200 * len(targets) <= m, "blocks overflow the map"
    for f in range(len(targets)):
        for k in range(het_at + 200 * f, het_at + 200 * (f + 1)):
            geno[2 * f][k] = geno[2 * f + 1][k] = HET

    Path(out + ".fam").write_text("".join(
        f"H{f:02d} H{f:02d}_A 0 0 1 -9\nH{f:02d} H{f:02d}_B 0 0 2 -9\n"
        for f in range(len(targets))))
    Path(out + ".bim").write_text("".join(
        f"1\trs{k}\t{(k + 1) * 0.05:.6f}\t{(k + 1) * 50000}\tA\tG\n" for k in range(m)))
    nb = (n + 3) // 4
    with open(out + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        for k in range(m):
            buf = bytearray(nb)
            for s in range(n):
                buf[s >> 2] |= geno[s][k] << (2 * (s & 3))
            fh.write(bytes(buf))
    print(f"{out}: {n} samples, {m} markers, ties " +
          ", ".join(f"{a}/{b}" for a, b in targets))


def main():
    if len(sys.argv) != 3 or sys.argv[1] not in ("error", "dupgate", "homibs0"):
        sys.exit(__doc__)
    which, out = sys.argv[1], sys.argv[2]
    if which == "error":
        build(ERROR_SPEC, out)
    elif which == "dupgate":
        build(DUPGATE_SPEC, out, DUPGATE_OFFSETS)
    else:
        build_homibs0(out)


main()
