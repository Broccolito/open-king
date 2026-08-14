"""With the gate in place, is the "at least 2 clean words" rule still needed?"""

import gate_corpus as GC
import rules3 as R3

BASE = dict(GC.BASE)

print(f"{'variant':<44}{'exact':>7}{'IBD1':>7}{'IBD2':>7}{'extra':>7}{'miss':>6}{'MAE':>10}")
for min1 in (1, 2, 3):
    for gate in (0, 10):
        kw = dict(BASE)
        kw["min1"] = min1
        _, t = GC.score(R3.P3(**kw, gate=gate))
        label = f"min1={min1}, gate={'10' if gate else 'off'}"
        print(f"{label:<44}{t['exact']:>7}{t['ibd1']:>7}{t['ibd2']:>7}"
              f"{t['extra']:>7}{t['missing']:>6}{t['mae']:>10.5f}")
