"""Characterise where a candidate rule disagrees with the reference.

    python3 tests/parity/fit/analyze.py extras     # the over-called pairs
    python3 tests/parity/fit/analyze.py values     # the wrong-value rows
    python3 tests/parity/fit/analyze.py po         # PO rows: does the denominator close?
    python3 tests/parity/fit/analyze.py perds
"""

import sys
from collections import Counter

import kingdata as kd
import rules as R
import fit


def extras(p=R.Params()):
    got = []
    fit.score(p, collect=got)
    ex = [g for g in got if g[3] == "extra"]
    print(f"{len(ex)} extra rows")
    print("  by dataset:", Counter(g[0] for g in ex).most_common())
    print(f"{'ds':<12}{'longestMb':>10}{'ibd1Mb':>9}{'ibd2Mb':>9}{'prop':>9}{'nseg':>6}")
    rows = []
    for name, i, j, _, ibd1, ibd2, longest, _ in ex:
        ds = kd.load(name)
        prop = (ibd2 + ibd1 / 2) / ds.denom
        rows.append((longest / 1e6, ibd1 / 1e6, ibd2 / 1e6, prop, name))
    rows.sort()
    for r in rows[:15] + [("...",) * 5] + rows[-10:]:
        if r[0] == "...":
            print("   ...")
            continue
        print(f"{r[4]:<12}{r[0]:10.2f}{r[1]:9.2f}{r[2]:9.2f}{r[3]:9.4f}")
    print("\n  longest-segment quantiles (Mb):",
          [round(sorted(x[0] for x in rows)[k * (len(rows) - 1) // 10], 2)
           for k in range(11)])
    print("  PropIBD quantiles:",
          [round(sorted(x[3] for x in rows)[k * (len(rows) - 1) // 10], 4)
           for k in range(11)])
    # what does the reference report at the low end?
    lo = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        for v in ds.ref.values():
            lo.append((v[2], name))
    lo.sort()
    print("\n  lowest 15 PropIBD the reference DOES report:",
          [(round(x[0], 4), x[1]) for x in lo[:15]])


def values(p=R.Params()):
    got = []
    fit.score(p, collect=got)
    bad = [g for g in got if g[3] == "value"]
    print(f"{len(bad)} rows present in both but numerically wrong")
    print("  by dataset:", Counter(g[0] for g in bad).most_common())
    print(f"{'ds':<10}{'refIBD1':>9}{'gotIBD1':>9}{'d1':>9}"
          f"{'refIBD2':>9}{'gotIBD2':>9}{'d2':>9}  ref")
    rows = []
    for name, i, j, _, ibd1, ibd2, longest, ref in bad:
        ds = kd.load(name)
        g1, g2 = ibd1 / ds.denom, ibd2 / ds.denom
        rows.append((abs(g1 - ref[0]) + abs(g2 - ref[1]), name, ref, g1, g2))
    rows.sort(reverse=True)
    for e, name, ref, g1, g2 in rows[:20]:
        print(f"{name:<10}{ref[0]:9.4f}{g1:9.4f}{g1-ref[0]:+9.4f}"
              f"{ref[1]:9.4f}{g2:9.4f}{g2-ref[1]:+9.4f}  {ref[3]}")
    d1 = [r[3] - r[2][0] for r in rows]
    d2 = [r[4] - r[2][1] for r in rows]
    print(f"\n  IBD1 delta: mean {sum(d1)/len(d1):+.4f}  "
          f"neg {sum(1 for x in d1 if x < 0)}  pos {sum(1 for x in d1 if x > 0)}")
    print(f"  IBD2 delta: mean {sum(d2)/len(d2):+.4f}  "
          f"neg {sum(1 for x in d2 if x < 0)}  pos {sum(1 for x in d2 if x > 0)}")
    # split by reference InfType
    by = {}
    for e, name, ref, g1, g2 in rows:
        by.setdefault(ref[3], []).append((g1 - ref[0], g2 - ref[1]))
    for k, v in sorted(by.items()):
        m1 = sum(x[0] for x in v) / len(v)
        m2 = sum(x[1] for x in v) / len(v)
        print(f"  {k:<8} n={len(v):4d}  dIBD1 {m1:+.4f}  dIBD2 {m2:+.4f}")


def po():
    """Solve for the denominator: which D makes the reference's PO rows exactly 1.0000?"""
    for name in kd.DATASETS:
        ds = kd.load(name)
        n_po = sum(1 for v in ds.ref.values() if v[3] == "PO")
        exact = sum(1 for v in ds.ref.values()
                    if v[3] == "PO" and v[0] == 1.0 and v[1] == 0.0 and v[2] == 0.5)
        tot_len = ds.denom
        alt = sum(int(ds.pos[hi]) - int(ds.pos[lo]) + 1 for _, lo, hi in ds.segs)
        print(f"{name:<12} segs={len(ds.segs):3d}  D(last-first)={tot_len/1e6:9.3f}Mb  "
              f"D(+1 per seg)={alt/1e6:9.3f}Mb  PO rows={n_po:3d} at 1.0000/0.0000/0.5000: {exact}")


def perds(p=R.Params()):
    print(f"{'dataset':<12}{'ref':>5}{'got':>5}{'extra':>6}{'miss':>5}"
          f"{'exact':>6}{'mae':>9}{'worst':>8}")
    for name in kd.DATASETS:
        t = fit.score(p, datasets=[name])
        print(f"{name:<12}{t['ref_rows']:5d}{t['got_rows']:5d}{t['extra']:6d}"
              f"{t['missing']:5d}{t['exact']:6d}{t['mae']:9.5f}{t['worst']:8.4f}")


if __name__ == "__main__":
    globals()[sys.argv[1] if len(sys.argv) > 1 else "perds"]()
