"""Which degrees of freedom in the gate does the corpus actually constrain?"""

import gate_corpus as GC
import rules3 as R3

BASE = GC.BASE

VARIANTS = [
    ("IBD1 gate only (IBD2 ungated)", dict(gate=10, gate2=0)),
    ("both gated, C=10, IBD2 mask = share", dict(gate=10, gate2=-1)),
    ("both gated, IBD2 mask = IBD1 mask", dict(gate=10, gate2=-1, gate2_mask="hom")),
    ("IBD2 gate C=8", dict(gate=10, gate2=8)),
    ("IBD2 gate C=12", dict(gate=10, gate2=12)),
    ("IBD2 gate C=16", dict(gate=10, gate2=16)),
    ("IBD2 gate C=24", dict(gate=10, gate2=24)),
]

print(f"{'variant':<40}{'exact':>7}{'IBD1':>7}{'IBD2':>7}{'extra':>7}{'miss':>6}{'MAE':>10}")
for label, kw in VARIANTS:
    _, t = GC.score(R3.P3(**BASE, **kw))
    print(f"{label:<40}{t['exact']:>7}{t['ibd1']:>7}{t['ibd2']:>7}"
          f"{t['extra']:>7}{t['missing']:>6}{t['mae']:>10.5f}")
