#!/usr/bin/env python3
"""Differential matrix for KING's pre-segment BIM ordering validation."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


ANALYSES = (
    "--related",
    "--ibs",
    "--unrelated",
    "--build",
    "--bysample",
    "--bySNP",
    "--cluster",
    "--ibdseg",
)
TIME = re.compile(rb"[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}")
PROGRESS = re.compile(rb"(?:[0-9]+%\r)+")


def groups(rows):
    result = []
    start = 0
    while start < len(rows):
        stop = start + 1
        while stop < len(rows) and rows[stop][0] == rows[start][0]:
            stop += 1
        result.append((start, stop))
        start = stop
    return result


def write_fileset(prefix: Path, rows, fam: bytes, bed: bytes) -> None:
    prefix.with_suffix(".fam").write_bytes(fam)
    prefix.with_suffix(".bim").write_text("".join(" ".join(row) + "\n" for row in rows))
    prefix.with_suffix(".bed").write_bytes(bed)


def fixtures(golden: Path, root: Path):
    fam = (golden / "multifam.fam").read_bytes()
    rows = [line.split() for line in (golden / "multifam.bim").read_text().splitlines()]
    bed = (golden / "multifam.bed").read_bytes()
    blocks = groups(rows)

    positions = [row[:] for row in rows]
    for start, stop in blocks:
        reversed_bp = [row[3] for row in positions[start:stop]][::-1]
        for offset, bp in enumerate(reversed_bp):
            positions[start + offset][3] = bp
    position_prefix = root / "positions"
    write_fileset(position_prefix, positions, fam, bed)

    samples = sum(bool(line.strip()) for line in fam.splitlines())
    bytes_per_variant = (samples + 3) // 4
    order = [index for start, stop in reversed(blocks) for index in range(start, stop)]
    chromosome_bed = bed[:3] + b"".join(
        bed[3 + index * bytes_per_variant : 3 + (index + 1) * bytes_per_variant]
        for index in order
    )
    chromosome_prefix = root / "chromosomes"
    write_fileset(chromosome_prefix, [rows[index] for index in order], fam, chromosome_bed)
    return {"positions": position_prefix, "chromosomes": chromosome_prefix}


def normalized(stdout: bytes) -> bytes:
    return TIME.sub(b"<TIME>", PROGRESS.sub(b"", stdout))


def run(binary: Path, prefix: Path, analysis: str, work: Path):
    proc = subprocess.run(
        [
            str(binary),
            "-b",
            str(prefix.with_suffix(".bed")),
            analysis,
            "--cpus",
            "1",
        ],
        cwd=work,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    files = {path.name: path.read_bytes() for path in work.iterdir() if path.is_file()}
    return proc.returncode, normalized(proc.stdout), proc.stderr, files


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()
    ref_binary = args.ref.resolve()
    impl_binary = args.impl.resolve()
    golden = Path(__file__).resolve().parents[1] / "golden"

    with tempfile.TemporaryDirectory(prefix="open-king-map-order-") as tmp:
        root = Path(tmp)
        shaped = fixtures(golden, root)
        for shape, prefix in shaped.items():
            for analysis in ANALYSES:
                ref_work = root / f"ref-{shape}-{analysis[2:]}"
                impl_work = root / f"impl-{shape}-{analysis[2:]}"
                ref_work.mkdir()
                impl_work.mkdir()
                reference = run(ref_binary, prefix, analysis, ref_work)
                implementation = run(impl_binary, prefix, analysis, impl_work)
                if reference != implementation:
                    print(f"FAIL {shape} {analysis}")
                    if reference[:3] != implementation[:3]:
                        print(f"reference process={reference[:3]!r}")
                        print(f"implementation process={implementation[:3]!r}")
                    print(f"reference files={sorted(reference[3])}")
                    print(f"implementation files={sorted(implementation[3])}")
                    for name in reference[3].keys() & implementation[3].keys():
                        if reference[3][name] != implementation[3][name]:
                            print(f"different file: {name}")
                    return 1
                print(f"PASS {shape:11} {analysis}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
