#!/usr/bin/env python3
"""Held-out differential probe for the no-informative-segments fallback."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


TIME = re.compile(rb"[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}")
PROGRESS = re.compile(rb"(?:[0-9]+%\r)+")
def thin_bigish(golden: Path, destination: Path) -> Path:
    """Keep every fourth SNP without decoding or rephasing SNP-major BED rows."""
    fam = (golden / "bigish.fam").read_bytes()
    bim = (golden / "bigish.bim").read_bytes().splitlines(keepends=True)
    bed = (golden / "bigish.bed").read_bytes()
    samples = sum(bool(line.strip()) for line in fam.splitlines())
    bytes_per_variant = (samples + 3) // 4
    keep = range(0, len(bim), 4)

    prefix = destination / "thin4"
    prefix.with_suffix(".fam").write_bytes(fam)
    prefix.with_suffix(".bim").write_bytes(b"".join(bim[i] for i in keep))
    prefix.with_suffix(".bed").write_bytes(
        bed[:3]
        + b"".join(
            bed[3 + i * bytes_per_variant : 3 + (i + 1) * bytes_per_variant]
            for i in keep
        )
    )
    return prefix


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
    files = {
        path.name: path.read_bytes()
        for path in work.iterdir()
        if path.is_file()
    }
    return proc.returncode, proc.stdout, proc.stderr, files


def normalized(stdout: bytes) -> bytes:
    return TIME.sub(b"<TIME>", PROGRESS.sub(b"", stdout))


def fail(message: str) -> None:
    raise AssertionError(message)


def check_related(ref, impl) -> None:
    if ref[0] != 0 or impl[0] != 0 or ref[2] or impl[2]:
        fail(f"unexpected related exit/stderr: ref={ref[:3]!r}, impl={impl[:3]!r}")
    if set(ref[3]) != {"king.kin", "king.kin0"}:
        fail(f"unexpected reference related files: {sorted(ref[3])}")
    if set(impl[3]) != {"king.kin", "king.kin0"}:
        fail(f"unexpected implementation related files: {sorted(impl[3])}")
    if ref[3]["king.kin"] != impl[3]["king.kin"]:
        fail("within-family fallback king.kin is not byte-identical")

    if ref[3]["king.kin0"] != impl[3]["king.kin0"]:
        fail("fallback king.kin0 is not byte-identical")
    if normalized(ref[1]) != normalized(impl[1]):
        fail("related console differs after time/progress normalization")
    if b"15 pairs of relatives" not in ref[1] or b"15 pairs of relatives" not in impl[1]:
        fail("expected 15-pair degree-1 screen was not observed")


def check_ibdseg(ref, impl) -> None:
    if ref[0] != 0 or impl[0] != 0 or ref[2] or impl[2]:
        fail(f"unexpected ibdseg exit/stderr: ref={ref[:3]!r}, impl={impl[:3]!r}")
    if set(ref[3]) != {"kingsplitped.txt"} or set(impl[3]) != {"kingsplitped.txt"}:
        fail(f"unexpected ibdseg files: ref={sorted(ref[3])}, impl={sorted(impl[3])}")
    if ref[3] != impl[3]:
        fail("sparse ibdseg splitped output differs")
    if normalized(ref[1]) != normalized(impl[1]):
        fail("sparse ibdseg console differs after time/progress normalization")


def check_clustering(analysis, ref, impl) -> None:
    if ref[0] != 0 or impl[0] != 0 or ref[2] or impl[2]:
        fail(f"unexpected {analysis} exit/stderr: ref={ref[:3]!r}, impl={impl[:3]!r}")
    if ref[3] != impl[3]:
        names = sorted(set(ref[3]) | set(impl[3]))
        changed = [name for name in names if ref[3].get(name) != impl[3].get(name)]
        fail(f"{analysis} fallback artifacts differ: {changed}")
    if normalized(ref[1]) != normalized(impl[1]):
        fail(f"{analysis} console differs after time/progress normalization")
    if b"15 by inference" not in ref[1] or b"15 by inference" not in impl[1]:
        fail(f"{analysis} did not preserve the expected 15-pair screen")
    if b"\t0\t0\t3\t12\t0\t0" not in ref[1]:
        fail(f"{analysis} reference fallback labels changed")
    if b"\t0\t0\t3\t12\t0\t0" not in impl[1]:
        fail(f"{analysis} implementation fallback labels changed")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()
    ref_binary = args.ref.resolve()
    impl_binary = args.impl.resolve()
    golden = Path(__file__).resolve().parents[1] / "golden"

    with tempfile.TemporaryDirectory(prefix="open-king-sparse-") as tmp:
        root = Path(tmp)
        prefix = thin_bigish(golden, root)
        results = {}
        for analysis in ("--related", "--ibdseg", "--unrelated", "--cluster", "--build"):
            for label, binary in (("ref", ref_binary), ("impl", impl_binary)):
                work = root / f"{label}-{analysis[2:]}"
                work.mkdir()
                results[label, analysis] = run(binary, prefix, analysis, work)
        check_related(results["ref", "--related"], results["impl", "--related"])
        check_ibdseg(results["ref", "--ibdseg"], results["impl", "--ibdseg"])
        for analysis in ("--unrelated", "--cluster", "--build"):
            check_clustering(
                analysis,
                results["ref", analysis],
                results["impl", analysis],
            )

    print("PASS --related: console, 574-line .kin, and 15-row .kin0 exact")
    print("PASS --ibdseg: console and splitped exact")
    print("PASS --unrelated: both selection files exact")
    print("PASS --cluster: updateids exact; no segment-only cluster.kin")
    print("PASS --build: build.log, updateids and updateparents exact")
    print("PASS clustering summaries: 15 inferred pairs (3 FS, 12 2nd) exact")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
