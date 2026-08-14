#!/usr/bin/env python3
"""Differential probe for `<prefix>X.seg`, out of sample.

Establishes, against the reference binary, everything implemented in
`crate::analysis::xseg`:

    written  <->  --degree != 0  AND  the X map yields a usable segment
    rows     ==   <prefix>.seg's rows, in <prefix>.seg's order
    columns  ==   FID1 ID1 FID2 ID2 Sex1 Sex2 <IBD1Seg> <IBD2Seg> <PropIBD> TAB
    values   ==   the autosomal segment caller over the X array, PropIBD at full
                  precision (not `.seg`'s printed-column rule)

The corpus captures only two X.seg files — 28 rows of one 6-sample family — so the
rules cannot be validated there without fitting to them. This probe builds fresh
filesets instead: eight pedigree shapes crossed with five X maps, two seeds each,
thirteen flag combinations, 1 040 reference-vs-open-king runs.

    python3 tests/parity/probes/xseg_probe.py --impl target/release/king

Reports four numbers:

    presence                    the emission gate, both directions.  Must be N/N.
    default                     X.seg byte-identical at the default 3 Mb floor.
    default_given_autosome_ok   the same, restricted to runs whose autosomal .seg
                                is byte-identical.  Must be N/N — anything else is
                                an X.seg bug rather than an inherited one.
    nondefault                  --seglength 5 / 10, which inherit the open run-merge
                                residual (docs/PARITY.md §5.0) and are informational.

Two conditions are deliberately kept out of the maps, both documented in
docs/PARITY.md §5.11 and neither an X.seg property:

  * an X array whose length is an exact multiple of 64 — the reference reads past the
    end and adds an absolute coordinate to the pair's IBD2 total, in X.kin as well;
  * a `.fam` SEX field outside {0, 1, 2}, which the reference and `king-io` read
    differently.

Python 3 standard library only.  The genotype simulator is `xgen.py` beside this file.
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import xgen  # noqa: E402

REFERENCE = "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"

# --------------------------------------------------------------------------
# Pedigrees
# --------------------------------------------------------------------------


def nuclear(fid, p, nson, ndau):
    out = [(fid, p + "FA", "0", "0", 1), (fid, p + "MO", "0", "0", 2)]
    out += [(fid, f"{p}S{k + 1}", p + "FA", p + "MO", 1) for k in range(nson)]
    out += [(fid, f"{p}D{k + 1}", p + "FA", p + "MO", 2) for k in range(ndau)]
    return out


def threegen(fid, p):
    out = [(fid, p + "GF", "0", "0", 1), (fid, p + "GM", "0", "0", 2),
           (fid, p + "FA", p + "GF", p + "GM", 1),
           (fid, p + "MO", "0", "0", 2),
           (fid, p + "AU", p + "GF", p + "GM", 2)]
    out += [(fid, f"{p}C{k + 1}", p + "FA", p + "MO", 1 + k % 2) for k in range(4)]
    return out


HUGE = [r for k in range(6) for r in nuclear(f"H{k}", f"H{k}_", 3, 3)]  # 48 samples

PEDS = {
    # three families, one of them a pair of unrelated founders
    "twofam": (nuclear("F1", "A_", 2, 2) + nuclear("F2", "B_", 1, 2)
               + [("F3", "C_1", "0", "0", 1), ("F3", "C_2", "0", "0", 2)]),
    "onefam": nuclear("F1", "A_", 3, 3),
    "threegen": threegen("G1", "T_") + nuclear("G2", "N_", 2, 1),
    # 28 samples: crosses the 16-sample tiling boundary of <prefix>.seg's row order
    "wide": (nuclear("F1", "A_", 3, 3) + nuclear("F2", "B_", 3, 3)
             + nuclear("F3", "C_", 2, 2) + nuclear("F4", "D_", 2, 2)),
    # 48 samples: three tiles
    "huge48": HUGE,
    # sex 0 in the middle of a reported family — not excluded, unlike in X.kin
    "sex0": ([("F1", "A_FA", "0", "0", 1), ("F1", "A_MO", "0", "0", 2),
              ("F1", "A_S1", "A_FA", "A_MO", 1), ("F1", "A_S2", "A_FA", "A_MO", 0),
              ("F1", "A_D1", "A_FA", "A_MO", 2), ("F1", "A_D2", "A_FA", "A_MO", 0)]
             + nuclear("F2", "B_", 1, 1)),
    # no family of two: the within-family stages do not run
    "unrel": [(f"U{k}", f"U{k}_1", "0", "0", 1 + k % 2) for k in range(12)],
    # below the <5 sample downgrade: no .seg at all, so no X.seg either
    "tiny4": nuclear("F1", "A_", 1, 1),
}

# --------------------------------------------------------------------------
# Maps.  No X array length is a multiple of 64 (see the module docstring).
# --------------------------------------------------------------------------

AUTO = [1, 2]
MAPS = {
    "x1000": xgen.make_map(AUTO, 2000, 100_000_000, 1000, 60_000_000),
    # 333 markers: above the five-word floor, far below --kinship's 512
    "x333": xgen.make_map(AUTO, 2000, 100_000_000, 333, 25_000_000),
    "x1500": xgen.make_map(AUTO, 2000, 100_000_000, 1500, 75_000_000),
    # three usable X segments, cut by >1 Mb gaps
    "x3seg": (xgen.make_map(AUTO, 2000, 100_000_000, 0, 1)
              + [("23", f"rsXa{m}", 1_000_000 + m * 60_000) for m in range(401)]
              + [("23", f"rsXb{m}", 40_000_000 + m * 60_000) for m in range(403)]
              + [("23", f"rsXc{m}", 80_000_000 + m * 60_000) for m in range(407)]),
    # X alongside Y, XY and MT, as the `sexchr` corpus fixture has it
    "xfull": xgen.make_map(AUTO, 2000, 100_000_000, 901, 60_000_000,
                           extra=[(24, 300, 20_000_000), (25, 150, 10_000_000),
                                  (26, 50, 100_000)]),
}

DEFAULT_FLOOR = [
    ["--ibdseg"],
    ["--ibdseg", "--degree", "0"],
    ["--ibdseg", "--degree", "-1"],
    ["--ibdseg", "--degree", "1"],
    ["--ibdseg", "--degree", "2"],
    ["--ibdseg", "--degree", "3"],
    ["--ibdseg", "--degree", "4"],
    ["--ibdseg", "--degree", "2", "--prefix", "Q"],
    ["--ibdseg", "--degree", "2", "--cpus", "1"],
    ["--related", "--degree", "1", "--ibdseg"],
    ["--related", "--degree", "2", "--ibdseg"],
]
RAISED_FLOOR = [
    ["--ibdseg", "--degree", "2", "--seglength", "5"],
    ["--ibdseg", "--degree", "2", "--seglength", "10"],
]


def run(binary, workdir, args):
    """Run one invocation and return every `.seg` file it wrote."""
    for f in os.listdir(workdir):
        if not f.startswith("in."):
            os.remove(os.path.join(workdir, f))
    subprocess.run([binary, "-b", "in.bed", *args], cwd=workdir,
                   capture_output=True, text=True)
    return {f: (Path(workdir) / f).read_bytes()
            for f in sorted(os.listdir(workdir))
            if not f.startswith("in.") and f.endswith(".seg")}


def unified(ref: bytes | None, ours: bytes | None) -> list[str]:
    rl = (ref or b"").decode().splitlines()
    ol = (ours or b"").decode().splitlines()
    out = []
    for k in range(max(len(rl), len(ol))):
        a = rl[k] if k < len(rl) else "<none>"
        b = ol[k] if k < len(ol) else "<none>"
        if a != b:
            out.append(f"      ref  {a}")
            out.append(f"      ours {b}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--impl", required=True, help="the open-king binary under test")
    ap.add_argument("--ref", default=REFERENCE, help="the KING 2.3.2 reference binary")
    ap.add_argument("--seed", type=int, default=500_000)
    ap.add_argument("--reps", type=int, default=2, help="seeds per (pedigree, map)")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="print every differing row, not just the tag")
    a = ap.parse_args()
    # Both binaries run with the fileset directory as cwd, so a relative --impl would
    # not resolve there.
    a.impl = str(Path(a.impl).resolve())
    a.ref = str(Path(a.ref).resolve())

    stats = {k: [0, 0] for k in
             ("presence", "default", "default_given_autosome_ok", "nondefault")}
    fails: list[str] = []
    root = tempfile.mkdtemp(prefix="xseg_probe")
    seed = a.seed

    for pname, ped in PEDS.items():
        for mname, rows in MAPS.items():
            for rep in range(a.reps):
                seed += 13
                ref_dir = os.path.join(root, f"{pname}_{mname}_{rep}")
                our_dir = ref_dir + "_o"
                os.makedirs(ref_dir, exist_ok=True)
                os.makedirs(our_dir, exist_ok=True)
                cwd = os.getcwd()
                os.chdir(ref_dir)
                try:
                    xgen.build("in", ped, rows, seed=seed)
                finally:
                    os.chdir(cwd)
                for f in ("in.bed", "in.bim", "in.fam"):
                    shutil.copy(os.path.join(ref_dir, f), os.path.join(our_dir, f))

                for kind, flags in (("default", DEFAULT_FLOOR),
                                    ("nondefault", RAISED_FLOOR)):
                    for args in flags:
                        r = run(a.ref, ref_dir, args)
                        o = run(a.impl, our_dir, args)
                        tag = f"{pname}/{mname}#{rep} {' '.join(args)}"

                        rx = {f: v for f, v in r.items() if f.endswith("X.seg")}
                        ox = {f: v for f, v in o.items() if f.endswith("X.seg")}
                        auto_ok = all(r.get(f) == o.get(f)
                                      for f in set(r) | set(o)
                                      if not f.endswith("X.seg"))

                        stats["presence"][1] += 1
                        if set(rx) == set(ox):
                            stats["presence"][0] += 1
                        else:
                            fails.append(f"PRESENCE {tag}: "
                                         f"ref={sorted(rx)} ours={sorted(ox)}")

                        stats[kind][1] += 1
                        if rx == ox:
                            stats[kind][0] += 1
                        elif set(rx) == set(ox):
                            fails.append(f"BYTES {tag} autosome_ok={auto_ok}")
                            if a.verbose:
                                for f in sorted(rx):
                                    fails += unified(rx[f], ox.get(f))

                        if kind == "default" and rx and auto_ok:
                            stats["default_given_autosome_ok"][1] += 1
                            if rx == ox:
                                stats["default_given_autosome_ok"][0] += 1

    for k in ("presence", "default", "default_given_autosome_ok", "nondefault"):
        got, want = stats[k]
        print(f"{k:26s} {got}/{want}")
    for line in fails:
        print(line)
    print("workdir:", root)

    hard = stats["presence"][0] == stats["presence"][1] and (
        stats["default_given_autosome_ok"][0]
        == stats["default_given_autosome_ok"][1])
    return 0 if hard else 1


if __name__ == "__main__":
    sys.exit(main())
