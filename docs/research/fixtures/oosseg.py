#!/usr/bin/env python3
"""Out-of-sample `.seg` differential — fresh filesets, unused seeds, both binaries.

The capture corpus is **saturated**: `tests/parity/fit/scorecard.py` reads 982 of 982 rows
byte-exact at 3, 5 and 10 Mb, so it can no longer grade the segment caller in either
direction (`docs/PARITY.md` §4.4, "How to grade further work on it"). The canvas rigs grade
constructed word sequences. This rig grades the thing in between and the thing a *user*
actually runs: whole simulated filesets the caller has never seen, on seeds that appear
nowhere else in this repository, compared **byte for byte** against KING 2.3.2.

    python3 oosseg.py --ref /path/to/king                 # the committed draw
    python3 oosseg.py --ref /path/to/king --impl ./target/release/king
    python3 oosseg.py --ref ... --seeds 5 6 7 --shapes twofam

Every fileset is built by `tests/parity/generate_corpus.py`'s own simulator, so the marker
map, the allele frequencies and the transmission model are the corpus's — only the pedigree
shapes and the seeds are new. Nothing is written into the repository: filesets live in a
temporary directory and are deleted.

# What it reports, and why the three numbers are separate

Per run it reports **extra** rows (pairs we report and the reference does not), **missing**
rows (the reverse) and **value-differing** rows (a pair both report, on which some printed
column differs). They fail for different reasons and a single "rows differ" count hides
that: a whole-file line diff charges a single dropped row against every row after it —
`twofam31415926` reads "39 of 106 rows differ" that way and is in fact **one** missing pair
with all 105 shared rows byte-identical.

# The measured result at the time of writing (the committed draw, 24 filesets x 3 floors)

    72 runs, 66 byte-identical; 6 of 6 713 rows wrong: 0 extra, 2 missing, 4 value-differing

and both residuals have a shape:

* **4 value rows** — one full-sib pair (`FA A_C2 / A_C3`) on two filesets at 3 and 5 Mb.
  `IBD1Seg` is exact, `IBD2Seg` is **low** by 0.0181-0.0182 on all four, and `PropIBD`
  follows it. The deficit is the same to a printed ulp across two independent seeds, which
  points at one missed IBD2 piece in one place on the shared map rather than at a
  data-dependent error.
* **2 missing rows** — the same distant pair of padding singletons, dropped at 5 and 10 Mb
  but reported correctly at 3. Only the `.seg` **pair filter** (longest single segment
  > 10 Mb, `king_core::ibdseg::pair_segments`) can drop a pair, and ours reads the unmerged
  calls *at the requested floor*. The reference reports this pair at all three floors, so
  whatever it reads is floor-independent. That is the sharpest open hypothesis this rig has
  produced and it is **not** landed here: it is one pair.

# It is also how the window bound was validated out of sample

`docs/research/23-gap-bound.md`'s window bound is invisible to the corpus at 3 and 5 Mb and
worth two cases at 10. Re-running this rig against a binary with `WINDOW_FRACTION` disabled
scores **60 of 72** where the shipped binary scores **66 of 72** — six further filesets
wrong, every one of them at `--seglength 10`, and none right that the shipped binary gets
wrong. Out-of-sample and in one direction only, which is the bar of `MAINTAINING.md` §8.6.

Exit status is 0 iff every run is byte-identical.
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
import subprocess
import sys
import tempfile

_HERE = os.path.dirname(os.path.abspath(__file__))
_GC = os.path.normpath(os.path.join(_HERE, "..", "..", "..", "tests", "parity",
                                    "generate_corpus.py"))

# Seeds used nowhere else in this repository — not by `generate_corpus.py` (master
# 20260813), not by the canvas rigs, not by `avfs.py`, not by `window1.py`. Change them
# only by ADDING; a seed that has ever graded a landed change is no longer out of sample.
SEEDS = [424242424, 13572468, 8675309, 20260814, 1010101, 777000777, 31415926, 27182818]

FLOORS = ["3", "5", "10"]


def _corpus():
    spec = importlib.util.spec_from_file_location("gc", _GC)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _pad(gc, ped, n):
    """Singletons, to push the sample count past the reference's own gates."""
    for k in range(n):
        ped.add("SG%03d" % k, "SF%03d" % k, sex=1 + (k % 2))


def _build(gc, outdir, name, seed, shape):
    """One fileset. The shapes differ in which relationships the caller has to find."""
    ped = gc.Ped()
    if shape == "twofam":
        # Two nuclear families whose fathers are undeclared full sibs: avuncular,
        # cousin and full-sib pairs at once.
        phantom = gc.add_couple(ped, "PH", "PH", emit=False)
        gc.add_nuclear(ped, "FA", "A", 4, father_parents=phantom)
        gc.add_nuclear(ped, "FB", "B", 4, father_parents=phantom)
        _pad(gc, ped, 90)
        nsnp = 40000
    elif shape == "threegen":
        # Three generations in one family: PO, FS and grandparental pairs.
        gp = gc.add_couple(ped, "G", "G", emit=True)
        gc.add_nuclear(ped, "G", "P", 3, father_parents=gp)
        _pad(gc, ped, 60)
        nsnp = 30000
    elif shape == "multi":
        # Three sibships sharing one phantom couple: many 2nd-degree pairs.
        phantom = gc.add_couple(ped, "PH", "PH", emit=False)
        for f in range(3):
            gc.add_nuclear(ped, "F%d" % f, "K%d" % f, 3, father_parents=phantom)
        _pad(gc, ped, 80)
        nsnp = 50000
    else:
        raise SystemExit("unknown shape %r" % shape)
    spec = gc.Spec(name, ped, gc.AUTOSOMES, nsnp, notes="out-of-sample rig")
    os.makedirs(outdir, exist_ok=True)
    gc.simulate(spec, seed, outdir)
    return os.path.join(outdir, name + ".bed")


def _run(binary, bed, floor, wd):
    os.makedirs(wd, exist_ok=True)
    subprocess.run([binary, "-b", bed, "--cpus", "1", "--ibdseg", "--seglength", floor],
                   cwd=wd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   check=False)
    path = os.path.join(wd, "king.seg")
    return open(path, "rb").read() if os.path.exists(path) else b""


def _rows(blob):
    """Data rows keyed on the four identifier columns — the same match `measure_gaps.py`
    makes, so extra/missing/value-differing are counted the way §3 counts them."""
    if not blob:
        return {}
    out = {}
    for line in blob.decode().splitlines()[1:]:
        f = line.split()
        if len(f) >= 4:
            out[tuple(f[:4])] = f
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--ref", required=True, help="path to the KING 2.3.2 reference binary")
    ap.add_argument("--impl", default=os.path.normpath(
        os.path.join(_HERE, "..", "..", "..", "target", "release", "king")))
    ap.add_argument("--seeds", type=int, nargs="*", default=SEEDS)
    ap.add_argument("--shapes", nargs="*", default=["twofam", "threegen", "multi"])
    ap.add_argument("--floors", nargs="*", default=FLOORS)
    args = ap.parse_args()

    for p in (args.ref, args.impl):
        if not os.path.exists(p):
            raise SystemExit("not found: %s" % p)

    gc = _corpus()
    root = tempfile.mkdtemp(prefix="oosseg")
    runs = same = rows = extra = missing = valdiff = 0
    notes: list[str] = []
    try:
        for shape in args.shapes:
            for seed in args.seeds:
                name = "%s%d" % (shape, seed)
                d = os.path.join(root, name)
                bed = _build(gc, d, name, seed, shape)
                for fl in args.floors:
                    a = _run(args.ref, bed, fl, os.path.join(d, "ref" + fl))
                    b = _run(args.impl, bed, fl, os.path.join(d, "our" + fl))
                    ra, rb = _rows(a), _rows(b)
                    runs += 1
                    rows += len(ra)
                    if a == b:
                        same += 1
                        continue
                    ex = sorted(set(rb) - set(ra))
                    mi = sorted(set(ra) - set(rb))
                    vd = [k for k in sorted(set(ra) & set(rb)) if ra[k] != rb[k]]
                    extra += len(ex)
                    missing += len(mi)
                    valdiff += len(vd)
                    notes.append("%-22s L=%-3s of %3d rows: extra %d, missing %d, "
                                 "value-differing %d"
                                 % (name, fl, len(ra), len(ex), len(mi), len(vd)))
                    for k in mi:
                        notes.append("      MISSING  ref: " + " ".join(ra[k]))
                    for k in ex:
                        notes.append("      EXTRA    our: " + " ".join(rb[k]))
                    for k in vd:
                        notes.append("      ref: " + " ".join(ra[k]))
                        notes.append("      our: " + " ".join(rb[k]))
                shutil.rmtree(d, ignore_errors=True)
    finally:
        shutil.rmtree(root, ignore_errors=True)

    print("ref : %s" % args.ref)
    print("impl: %s" % args.impl)
    print("%d run(s) over %d fileset(s) x %d floor(s), %d reference row(s)"
          % (runs, runs // max(1, len(args.floors)), len(args.floors), rows))
    print("byte-identical: %d/%d" % (same, runs))
    print("rows: %d extra, %d missing, %d value-differing" % (extra, missing, valdiff))
    for n in notes:
        print("  " + n)
    return 0 if same == runs else 1


if __name__ == "__main__":
    sys.exit(main())
