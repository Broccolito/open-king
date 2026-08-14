#!/usr/bin/env python3
"""probe.py NAME  <<targets on stdin>>  -> prints realized (pi1,pi2,pi,label) rows."""
import subprocess
import sys

KING = "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"


def run(name, targets, extra=()):
    with open(name + ".txt", "w") as f:
        f.write("\n".join(f"{t[0]:.6f} {t[1]:.6f} {t[2] if len(t)>2 else 0}" for t in targets))
    subprocess.run([sys.executable, "mkpairs.py", name, name + ".txt"], check=True)
    subprocess.run([KING, "-b", name + ".bed", "--ibdseg", "--prefix", name, *extra],
                   check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    rows = []
    for line in open(name + ".seg").readlines()[1:]:
        r = line.rstrip("\n").split("\t")
        if r[0] == r[2]:
            rows.append((float(r[4]), float(r[5]), float(r[6]), r[7]))
    return rows


def bracket(rows, key, lo_label, hi_label):
    lo = max((key(x) for x in rows if x[3] == lo_label), default=None)
    hi = min((key(x) for x in rows if x[3] == hi_label), default=None)
    return lo, hi
