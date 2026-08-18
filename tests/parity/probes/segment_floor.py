#!/usr/bin/env python3
"""Differential probe for --ibdseg's closed 100,000,000 bp floor."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


TOTALS = {"below": 99_999_999, "exact": 100_000_000, "above": 100_000_001}
TIME = re.compile(rb"[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}")
PROGRESS = re.compile(rb"(?:[0-9]+%\r)+")


def fixtures(golden: Path, root: Path):
    fam = (golden / "multifam.fam").read_bytes()
    rows = [line.split() for line in (golden / "multifam.bim").read_text().splitlines()[:2001]]
    bed = (golden / "multifam.bed").read_bytes()
    samples = sum(bool(line.strip()) for line in fam.splitlines())
    bytes_per_variant = (samples + 3) // 4
    body = bed[: 3 + len(rows) * bytes_per_variant]

    result = {}
    for name, total in TOTALS.items():
        prefix = root / name
        shaped = []
        for index, row in enumerate(rows):
            variant = row[:]
            variant[0] = "1"
            variant[1] = f"v{index}"
            variant[2] = "0"
            variant[3] = str(1_000_000 + round(index * total / 2000))
            shaped.append(variant)
        prefix.with_suffix(".fam").write_bytes(fam)
        prefix.with_suffix(".bim").write_text(
            "".join(" ".join(variant) + "\n" for variant in shaped)
        )
        prefix.with_suffix(".bed").write_bytes(body)
        result[name] = prefix
    return result


def run(binary: Path, prefix: Path, work: Path):
    proc = subprocess.run(
        [
            str(binary),
            "-b",
            str(prefix.with_suffix(".bed")),
            "--ibdseg",
            "--cpus",
            "1",
        ],
        cwd=work,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    files = {path.name: path.read_bytes() for path in work.iterdir() if path.is_file()}
    stdout = TIME.sub(b"<TIME>", PROGRESS.sub(b"", proc.stdout))
    return proc.returncode, stdout, proc.stderr, files


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()
    ref_binary = args.ref.resolve()
    impl_binary = args.impl.resolve()
    golden = Path(__file__).resolve().parents[1] / "golden"

    with tempfile.TemporaryDirectory(prefix="open-king-segment-floor-") as tmp:
        root = Path(tmp)
        for name, prefix in fixtures(golden, root).items():
            ref_work = root / f"ref-{name}"
            impl_work = root / f"impl-{name}"
            ref_work.mkdir()
            impl_work.mkdir()
            reference = run(ref_binary, prefix, ref_work)
            implementation = run(impl_binary, prefix, impl_work)
            if reference != implementation:
                print(f"FAIL {name} ({TOTALS[name]} bp)")
                print(f"reference process={reference[:3]!r}")
                print(f"implementation process={implementation[:3]!r}")
                print(f"reference files={sorted(reference[3])}")
                print(f"implementation files={sorted(implementation[3])}")
                return 1
            print(f"PASS {name:5} {TOTALS[name]} bp files={','.join(sorted(reference[3]))}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
