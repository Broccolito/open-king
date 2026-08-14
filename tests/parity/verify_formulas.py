#!/usr/bin/env python3
"""Check the documented estimator formulas against captured reference output.

Every row below was read out of output produced by the reference KING 2.3.2 binary
(see docs/reference-captures/tiny/). If this script passes, the formulas in
docs/VERIFIED_FORMULAS.md reproduce the reference to the precision it prints.

    python3 tests/parity/verify_formulas.py
"""

import sys

# Rows captured from docs/reference-captures/tiny/out_ibs/king.ibs and king.ibs0,
# cross-referenced with out_kinship/king.kin and king.kin0 for the Kinship column.
ROWS = [
    ("f1dad/f1kid1", "within", dict(
        N=1970, ibs0=0, ibs1=593, ibs2=1377, hh=302, homhom=1075, h1=602, h2=595,
        IBS=1.6990, Dist=0.3010, HetConc=0.3374, H21=0.5017, H12=0.5076,
        HomConc=1.0000, Kin=0.2523)),
    ("f1dad/f1kid2", "within", dict(
        N=1942, ibs0=0, ibs1=594, ibs2=1348, hh=269, homhom=1079, h1=595, h2=537,
        IBS=1.6941, Dist=0.3059, HetConc=0.3117, H21=0.4521, H12=0.5009,
        HomConc=1.0000, Kin=0.2376)),
    ("f1dad/f2dad", "between", dict(
        N=1960, ibs0=28, ibs1=479, ibs2=1453, hh=337, homhom=1144, h1=595, h2=558,
        IBS=1.7270, Dist=0.3015, HetConc=0.4130, H21=0.5664, H12=0.6039,
        HomConc=0.9755, Kin=0.2352)),
    ("f1dad/f2mom", "between", dict(
        N=1803, ibs0=110, ibs1=673, ibs2=1020, hh=221, homhom=909, h1=549, h2=566,
        IBS=1.5047, Dist=0.6173, HetConc=0.2472, H21=0.4026, H12=0.3905,
        HomConc=0.8790, Kin=-0.0068)),
]


def r4(x):
    return float("%.4f" % x)


def main():
    failures = 0
    for name, kind, d in ROWS:
        N, ibs0, ibs1, ibs2 = d["N"], d["ibs0"], d["ibs1"], d["ibs2"]
        hh, homhom, h1, h2 = d["hh"], d["homhom"], d["h1"], d["h2"]

        checks = [
            ("N_IBS1 identity", h1 + h2 - 2 * hh, ibs1, True),
            ("N_IBS2 identity", N - ibs0 - (h1 + h2 - 2 * hh), ibs2, True),
            ("IBS", (ibs1 + 2 * ibs2) / N, d["IBS"], False),
            ("Dist", (ibs1 + 4 * ibs0) / N, d["Dist"], False),
            ("HetConc", hh / (h1 + h2 - hh), d["HetConc"], False),
            ("Het2|1", hh / h1, d["H21"], False),
            ("Het1|2", hh / h2, d["H12"], False),
            ("HomConc", (homhom - ibs0) / homhom, d["HomConc"], False),
        ]
        if kind == "within":
            checks.append(("Kinship", (hh - 2 * ibs0) / (h1 + h2), d["Kin"], False))
        else:
            checks.append(("Kinship",
                           0.5 + (2 * hh - 4 * ibs0 - h1 - h2) / (4 * min(h1, h2)),
                           d["Kin"], False))

        print("--- %s (%s-family) ---" % (name, kind))
        for label, got, exp, exact in checks:
            ok = (round(got) == exp) if exact else (r4(got) == r4(exp))
            failures += 0 if ok else 1
            print("  %s %-16s computed=%.6f reference=%s"
                  % ("OK  " if ok else "FAIL", label, got, exp))

    if failures:
        print("\n%d MISMATCH(ES) — docs/VERIFIED_FORMULAS.md is wrong or stale" % failures)
        return 1
    print("\nAll formulas reproduce the reference output.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
