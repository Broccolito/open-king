#!/usr/bin/env python3
"""Verify the row-ordering rules for .kin0 / .ibs0 against captured golden output.

Between-family pairs are emitted by a block-tiled loop keyed on
(i // B, j // B, i, j) over .fam index order, with B = 32 for .kin0 and B = 8 for
.ibs0. Plain ascending (i, j) coincides with this whenever n <= B, which is why it
must be checked on a large fileset.

    python3 tests/parity/verify_row_order.py
"""
import os
import sys

GOLDEN = os.path.join(os.path.dirname(__file__), "golden")
CASES = [
    ("dups .ibs0", "core/dups__ibs/king.ibs0", "dups", 8),
    ("unrelated .ibs0", "core/unrelated__ibs/king.ibs0", "unrelated", 8),
    ("bigish .ibs0", "core/bigish__ibs/king.ibs0", "bigish", 8),
    ("bigish .kin0", "core/bigish__kinship/king.kin0", "bigish", 32),
]


def fam_order(name):
    with open(os.path.join(GOLDEN, name + ".fam")) as fh:
        return [(l.split()[0], l.split()[1]) for l in fh if l.strip()]


def emitted(path):
    out = []
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out.append(((f[0], f[1]), (f[2], f[3])))
    return out


def tiled(order, b):
    n = len(order)
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    pairs.sort(key=lambda t: (t[0] // b, t[1] // b, t[0], t[1]))
    return [(order[i], order[j]) for i, j in pairs if order[i][0] != order[j][0]]


def main():
    failures = 0
    for label, rel, ds, b in CASES:
        path = os.path.join(GOLDEN, rel)
        if not os.path.exists(path):
            print("SKIP %-18s (golden missing: %s)" % (label, rel))
            continue
        order, obs = fam_order(ds), emitted(path)
        ok = tiled(order, b) == obs
        plain_ok = tiled(order, 10**9) == obs
        failures += 0 if ok else 1
        note = "" if not plain_ok else "  (n <= B: plain order also matches)"
        print("%s %-18s n=%-4d rows=%-6d B=%d%s"
              % ("OK  " if ok else "FAIL", label, len(order), len(obs), b, note))
    if failures:
        print("\n%d ordering rule(s) do not match golden output" % failures)
        return 1
    print("\nRow-ordering rules reproduce the reference.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
