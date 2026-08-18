#!/usr/bin/env python3
"""Differential probe for KING's case-insensitive `(FID, IID)` identity key."""

import argparse
import shutil
import subprocess
import tempfile
from pathlib import Path


def run(binary: Path, bed: Path, fam: Path, cwd: Path):
    proc = subprocess.run(
        [str(binary), "-b", str(bed), "--fam", str(fam), "--kinship", "--cpus", "1"],
        cwd=cwd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    text = proc.stdout.decode("utf-8", "replace")
    marker = text.find("Family ")
    tail = text[marker:] if marker >= 0 else text
    return proc.returncode, tail, proc.stderr


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[3]
    golden = root / "tests" / "parity" / "golden"
    rows = [line.split() for line in (golden / "multifam.fam").read_text().splitlines()]

    with tempfile.TemporaryDirectory(prefix="open-king-idcase-") as tmp:
        work = Path(tmp)
        bed = work / "case.bed"
        shutil.copyfile(golden / "multifam.bed", bed)
        shutil.copyfile(golden / "multifam.bim", work / "case.bim")

        # Duplicate one unreferenced child as another. Replacing a parent would make the
        # reference emit additional, unrelated parental-sex diagnostics.
        source = 2
        target = 3
        variants = {
            "exact": (rows[source][0], rows[source][1]),
            "iid-case": (rows[source][0], rows[source][1].lower()),
            "fid-case": (rows[source][0].lower(), rows[source][1]),
        }
        for name, (fid, iid) in variants.items():
            changed = [row[:] for row in rows]
            changed[target][0] = fid
            changed[target][1] = iid
            fam = work / f"{name}.fam"
            fam.write_text("".join(" ".join(row) + "\n" for row in changed))

            ref = run(args.ref.resolve(), bed, fam, work)
            impl = run(args.impl.resolve(), bed, fam, work)
            if ref != impl:
                print(f"FAIL {name}\nreference={ref!r}\nimplementation={impl!r}")
                return 1
            if ref[0] != 1 or b"" != ref[2] or "is duplicated" not in ref[1]:
                print(f"FAIL {name}: unexpected reference behavior {ref!r}")
                return 1
            print(f"PASS {name}: {ref[1].splitlines()[0]}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
