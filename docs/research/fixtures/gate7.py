"""Re-examine `10-segment-rule-fixtures.md` §2.4 — the "one clean word is data-dependent"
open item — now that the gate is known.

That section swept a 1-word IBD1 block across positions of the solid IBS0 background and
found an irregular, reproducible kept/dropped pattern that no length or count explained,
and concluded the minimum run must be 2 words.  If the gate is the whole story then
(a) the minimum run is really 1, and (b) each verdict is predicted by the number of
markers in that one word at which the pair shares A1 with at least one of them homozygous.

    python3 gate7.py sweep [maf ...]
"""

import sys

import fixlab as L
import rig2

WORD = 64


def informative(f, lo, hi):
    """Markers in [lo, hi) where both carry A1 and at least one is homozygous.

    `f.geno` holds A1 dosages after fixlab's re-orientation, so dosage >= 1 is
    "carries A1" and dosage == 2 is "homozygous for A1".
    """
    n = 0
    for m in range(lo, hi):
        a, b = f.geno[0][m], f.geno[1][m]
        if a == 3 or b == 3:
            continue
        if a >= 1 and b >= 1 and (a == 2 or b == 2):
            n += 1
    return n


def sweep(mafs=(0.5, 0.3, 0.2)):
    for maf in mafs:
        rig = rig2.Rig(spacing=100_000, n1=640, n2=640, nsample=6, maf=maf, seed=1)
        print(f"\nmaf={maf}: one clean word swept across the canvas "
              f"(a 1-word run reports {2 * WORD - 1} intervals)")
        print(f"{'word':>5} {'informative':>12} {'predict':>8} {'chr2_mk':>8}")
        agree = tot = 0
        for w in range(1, 9):
            f = rig.new(f"g7_{str(maf).replace('.', '')}_{w}", solid=True)
            rig.block(f, w * WORD, (w + 1) * WORD, L.IBD1)
            r = rig.read(f)
            lo, _ = f.chrom_span(1)
            n = informative(f, lo + w * WORD, lo + (w + 1) * WORD)
            got = 0 if r is None else r["ibd1_mk"] + r["ibd2_mk"] - (rig.n1 - 1)
            pred = "call" if n >= 10 else "drop"
            agree += (n >= 10) == (got > 0)
            tot += 1
            print(f"{w:>5} {n:>12} {pred:>8} {got:>8}")
        print(f"  agreement: {agree}/{tot}")


if __name__ == "__main__":
    sweep(tuple(float(x) for x in sys.argv[2:]) or (0.5, 0.3, 0.2))
