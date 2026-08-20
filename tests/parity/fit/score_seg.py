#!/usr/bin/env python3
"""Score our `.seg` against the captured reference `.seg`, row by row.

Runs the implementation under test on every corpus dataset that has an
``ibdseg/<ds>__ibdseg*`` capture and reports, per dataset and overall:

    rows      reference rows
    exact     rows whose IBD1Seg *and* IBD2Seg print identically
    ibd1/ibd2 per-column agreement
    extra     pairs we report and the reference does not
    missing   pairs the reference reports and we do not
    MAE       mean |our PropIBD - reference PropIBD| over shared rows

Usage:  python3 score_seg.py [--impl PATH] [--variant ibdseg|ibdseg_seglength5|...]
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARITY = HERE.parent
GOLDEN = PARITY / "golden" / "ibdseg"
DATA = PARITY / "work" / "data"

DATASETS = ["nuclear", "threegen", "multifam", "dups", "missing", "monomorphic",
            "sexchr", "unrelated", "admixed", "bigish"]

EXTRA_ARGS = {
    "ibdseg": [],
    "ibdseg_seglength5": ["--seglength", "5"],
    "ibdseg_seglength10": ["--seglength", "10"],
    "ibdseg_degree2": ["--degree", "2"],
}


def read_seg(path: Path) -> dict:
    rows = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[(f[0], f[1], f[2], f[3])] = (f[4], f[5], f[6], f[7])
    return rows


def run_impl(impl: str, ds: str, args: list[str], workdir: Path) -> dict:
    cmd = [impl, "-b", str(DATA / f"{ds}.bed"), "--ibdseg", *args]
    subprocess.run(cmd, cwd=workdir, capture_output=True, check=False)
    out = workdir / "king.seg"
    return read_seg(out) if out.exists() else {}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--impl", default=str(PARITY.parent.parent / "target/release/open-king"))
    ap.add_argument("--variant", default="ibdseg", choices=sorted(EXTRA_ARGS))
    ap.add_argument("--dump", help="dataset whose per-row diff to print")
    ap.add_argument("--limit", type=int, default=40)
    a = ap.parse_args()

    tot = dict(rows=0, exact=0, ibd1=0, ibd2=0, inf=0, extra=0, missing=0, err=0.0, worst=0.0)
    print(f"{'dataset':<13} {'rows':>5} {'exact':>6} {'ibd1':>6} {'ibd2':>6} "
          f"{'inf':>5} {'extra':>6} {'miss':>5} {'MAE':>9} {'worst':>7}")
    for ds in DATASETS:
        gold_dir = GOLDEN / f"{ds}__{a.variant}"
        gold_file = gold_dir / "king.seg"
        if not gold_file.exists():
            continue
        want = read_seg(gold_file)
        with tempfile.TemporaryDirectory() as td:
            got = run_impl(a.impl, ds, EXTRA_ARGS[a.variant], Path(td))
        shared = [k for k in want if k in got]
        exact = sum(1 for k in shared if got[k][0] == want[k][0] and got[k][1] == want[k][1])
        ibd1 = sum(1 for k in shared if got[k][0] == want[k][0])
        ibd2 = sum(1 for k in shared if got[k][1] == want[k][1])
        inf = sum(1 for k in shared if got[k][3] == want[k][3])
        errs = [abs(float(got[k][2]) - float(want[k][2])) for k in shared]
        mae = sum(errs) / len(errs) if errs else 0.0
        worst = max(errs) if errs else 0.0
        extra = len(set(got) - set(want))
        missing = len(set(want) - set(got))
        print(f"{ds:<13} {len(want):>5} {exact:>6} {ibd1:>6} {ibd2:>6} {inf:>5} "
              f"{extra:>6} {missing:>5} {mae:>9.5f} {worst:>7.4f}")
        tot["rows"] += len(want)
        tot["exact"] += exact
        tot["ibd1"] += ibd1
        tot["ibd2"] += ibd2
        tot["inf"] += inf
        tot["extra"] += extra
        tot["missing"] += missing
        tot["err"] += sum(errs)
        tot["worst"] = max(tot["worst"], worst)
        if a.dump == ds:
            n = 0
            for k in sorted(set(want) | set(got)):
                w, g = want.get(k), got.get(k)
                if w == g:
                    continue
                n += 1
                if n <= a.limit:
                    print(f"    {'/'.join(k):<28} ref={w} ours={g}")
    n = max(1, tot["rows"] - tot["missing"])
    print(f"{'ALL':<13} {tot['rows']:>5} {tot['exact']:>6} {tot['ibd1']:>6} "
          f"{tot['ibd2']:>6} {tot['inf']:>5} {tot['extra']:>6} {tot['missing']:>5} "
          f"{tot['err'] / n:>9.5f} {tot['worst']:>7.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
