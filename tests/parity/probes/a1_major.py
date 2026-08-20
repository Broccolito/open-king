#!/usr/bin/env python3
"""Differential regression for KING's A1-major orientation gate.

This pins the stable, observable contract without reproducing KING's uninitialized
short-map tail: the first 4,096 retained autosomal markers, a strict 10% boundary,
analysis-specific application, fatal placement, percentage, exit code and file set.

Usage:
    python3 a1_major.py --ref /path/to/king --impl target/release/open-king
"""

from __future__ import annotations

import argparse
import difflib
import re
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path


N_MARKERS = 5000
FATAL = "Too many first alleles as the major allele"
TIME = re.compile(r"KING (starts|ends) at [^\n]+")


@dataclass(frozen=True)
class Result:
    code: int
    output: str
    files: tuple[tuple[str, bytes], ...]

    @property
    def rejected(self) -> bool:
        return FATAL in self.output


def write_fileset(
    prefix: Path, sample_count: int, major: set[int], grouped_family: bool = False
) -> None:
    prefix.with_suffix(".fam").write_text(
        "".join(
            f"{'F' if grouped_family else f'F{s:02d}'} I{s:02d} 0 0 {1 + s % 2} -9\n"
            for s in range(sample_count)
        )
    )
    prefix.with_suffix(".bim").write_text(
        "".join(
            f"1\trs{m}\t{(m + 1) * 0.05:.6f}\t{(m + 1) * 50000}\tA\tG\n"
            for m in range(N_MARKERS)
        )
    )
    row_bytes = (sample_count + 3) // 4
    body = b"".join(
        bytes([0x00 if m in major else 0xFF]) * row_bytes for m in range(N_MARKERS)
    )
    prefix.with_suffix(".bed").write_bytes(bytes([0x6C, 0x1B, 0x01]) + body)


def run(binary: Path, prefix: Path, analysis: str, work: Path) -> Result:
    proc = subprocess.run(
        [str(binary), "-b", str(prefix.with_suffix(".bed")), analysis, "--prefix", "out"],
        cwd=work,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        check=False,
    )
    return Result(
        proc.returncode,
        proc.stdout,
        tuple((p.name, p.read_bytes()) for p in sorted(work.iterdir())),
    )


def normalized(output: str) -> str:
    return TIME.sub(lambda m: f"KING {m.group(1)} at <TIME>", output)


def compare(
    ref: Path,
    impl: Path,
    prefix: Path,
    analysis: str,
    root: Path,
    tag: str,
) -> tuple[Result, Result]:
    ref_work = root / f"ref-{tag}"
    impl_work = root / f"impl-{tag}"
    ref_work.mkdir()
    impl_work.mkdir()
    reference = run(ref, prefix, analysis, ref_work)
    implementation = run(impl, prefix, analysis, impl_work)
    assert (implementation.code, implementation.rejected) == (
        reference.code,
        reference.rejected,
    ), f"{tag}: outcome differs\nref={reference}\nimpl={implementation}"
    if reference.rejected:
        expected = normalized(reference.output)
        actual = normalized(implementation.output)
        assert actual == expected, f"{tag}: fatal console differs\n" + "".join(
            difflib.unified_diff(
                expected.splitlines(keepends=True),
                actual.splitlines(keepends=True),
                fromfile="reference",
                tofile="implementation",
            )
        )
        assert implementation.files == reference.files, (
            f"{tag}: fatal artifacts differ: "
            f"ref={[name for name, _ in reference.files]}, "
            f"impl={[name for name, _ in implementation.files]}"
        )
    print(f"PASS {tag}: rejected={reference.rejected}")
    return reference, implementation


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()
    ref = args.ref.resolve()
    impl = args.impl.resolve()

    with tempfile.TemporaryDirectory(prefix="a1-major-") as tmp:
        root = Path(tmp)
        fixtures: dict[tuple[int, str], Path] = {}

        def fixture(
            samples: int, name: str, major: set[int], grouped_family: bool = False
        ) -> Path:
            key = (samples, name)
            if key not in fixtures:
                prefix = root / f"data-{samples}-{name}"
                write_fileset(prefix, samples, major, grouped_family)
                fixtures[key] = prefix
            return fixtures[key]

        below = fixture(20, "below", set(range(409)))
        above = fixture(20, "above", set(range(410)))
        trailing = fixture(20, "trailing", set(range(4096, 5000)))
        compare(ref, impl, below, "--related", root, "boundary-409-of-4096")
        compare(ref, impl, above, "--related", root, "boundary-410-of-4096")
        compare(ref, impl, trailing, "--related", root, "markers-after-window")

        for analysis in (
            "--kinship",
            "--related",
            "--ibs",
            "--unrelated",
            "--build",
            "--bysample",
            "--bySNP",
            "--cluster",
            "--ibdseg",
            "--duplicate",
            "--autoQC",
        ):
            compare(ref, impl, above, analysis, root, f"surface-{analysis[2:]}")

        for analysis, passing_n, rejecting_n in (
            ("--related", 9, 10),
            ("--ibdseg", 4, 5),
            ("--cluster", 9, 10),
        ):
            passing = fixture(passing_n, "above", set(range(410)))
            rejecting = fixture(rejecting_n, "above", set(range(410)))
            compare(ref, impl, passing, analysis, root, f"{analysis[2:]}-n{passing_n}")
            compare(ref, impl, rejecting, analysis, root, f"{analysis[2:]}-n{rejecting_n}")

        grouped = fixture(20, "above-grouped", set(range(410)), grouped_family=True)
        compare(ref, impl, grouped, "--ibdseg", root, "ibdseg-splitped-artifact")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
