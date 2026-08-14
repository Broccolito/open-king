#!/usr/bin/env python3
"""Differential probe for the `--ibdseg --degree` reporting filter.

Establishes, against the reference binary, the rule implemented as
`king_core::ibdseg::reported_at_degree`:

    d == 0  ->  every pair with a >10 Mb segment is reported
    d >  0  ->  PropIBD > 2^-(d+0.5),  or, at d == 1 only, IBD2Seg >= 0.08
    d <  0  ->  PropIBD <= 2^-(|d|+0.5)   -- the comparison inverts

For every (dataset, --seglength, --degree) it runs the reference twice, once
unfiltered and once filtered, and checks the predicate against the reference's own
verdict for every pair. Reports false-keeps and false-drops; both must be zero.

The `d == 1` IBD2 clause is invisible here — a real first-degree pair has
IBD1Seg ~ 0.5, so its PropIBD clears 2^-1.5 on the first clause anyway, and no corpus
pair has IBD2Seg strictly between 0 and 0.1089. It is pinned by the constructed
fixture `docs/research/fixtures/gate8.py` instead. The two together are what the rule
rests on.

    python3 tests/parity/probes/degree_filter.py --ref /path/to/reference/king
    python3 tests/parity/probes/degree_filter.py --ref ... --seglengths 3 8 15

Python 3 standard library only. The full sweep is ~2 500 reference invocations.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

DATASETS = [
    "bigish", "multifam", "threegen", "admixed", "nuclear", "missing",
    "monomorphic", "sexchr", "dups", "unrelated", "trio", "pair", "singleton",
]

FIRST_DEGREE_IBD2 = 0.08


def reported_at_degree(degree: int, ibd2_seg: float, prop_ibd: float) -> bool:
    """The predicate under test — mirrors `king_core::ibdseg::reported_at_degree`."""
    if degree == 0:
        return True
    cut = 2.0 ** -(abs(degree) + 0.5)
    if degree < 0:
        return prop_ibd <= cut
    if degree == 1:
        return prop_ibd > cut or ibd2_seg >= FIRST_DEGREE_IBD2
    return prop_ibd > cut


def seg(ref: str, data: Path, dataset: str, args: list[str]) -> dict:
    """Run the reference and parse `king.seg` into {pair: (IBD1Seg, IBD2Seg, PropIBD)}."""
    tmp = Path(tempfile.mkdtemp(prefix="degree-"))
    try:
        subprocess.run(
            [ref, "-b", str(data / f"{dataset}.bed"), "--ibdseg", *args],
            cwd=tmp, capture_output=True, check=False,
        )
        out = tmp / "king.seg"
        if not out.is_file():
            return {}
        rows = {}
        for line in out.read_text().split("\n")[1:]:
            if not line:
                continue
            f = line.split("\t")
            rows[tuple(f[:4])] = (float(f[4]), float(f[5]), float(f[6]))
        return rows
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def main() -> int:
    here = Path(__file__).resolve().parents[1]
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="path to the KING 2.3.2 reference binary")
    ap.add_argument("--data", type=Path, default=here / "work" / "data")
    ap.add_argument("--degrees", type=int, nargs="+",
                    default=list(range(-6, 7)))
    ap.add_argument("--seglengths", nargs="+", default=["3", "8", "15"],
                    help="--seglength values in Mb; 3 is the default floor")
    args = ap.parse_args()

    if not Path(args.ref).is_file():
        print(f"no such binary: {args.ref}", file=sys.stderr)
        return 2
    if not args.data.is_dir():
        print(f"no corpus at {args.data}; run generate_corpus.py first", file=sys.stderr)
        return 2

    total = false_keep = false_drop = 0
    for sl in args.seglengths:
        base_args = ["--seglength", sl] if sl != "3" else []
        for dataset in DATASETS:
            base = seg(args.ref, args.data, dataset, base_args)
            if not base:
                continue
            for d in args.degrees:
                kept = seg(args.ref, args.data, dataset, base_args + ["--degree", str(d)])
                for pair, (_pi1, pi2, prop) in base.items():
                    total += 1
                    want = pair in kept
                    got = reported_at_degree(d, pi2, prop)
                    if got and not want:
                        false_keep += 1
                        print(f"FALSE-KEEP {dataset} seglength={sl} degree={d} "
                              f"{'/'.join(pair)} IBD2Seg={pi2} PropIBD={prop}")
                    elif want and not got:
                        false_drop += 1
                        print(f"FALSE-DROP {dataset} seglength={sl} degree={d} "
                              f"{'/'.join(pair)} IBD2Seg={pi2} PropIBD={prop}")
    print(f"\n{total} (dataset, seglength, degree, pair) cases: "
          f"false-keep {false_keep}, false-drop {false_drop}")
    return 1 if (false_keep or false_drop) else 0


if __name__ == "__main__":
    sys.exit(main())
