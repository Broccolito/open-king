#!/usr/bin/env python3
"""Held-out differential probe for the no-informative-segments fallback."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


TIME = re.compile(rb"[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}")
PROGRESS = re.compile(rb"(?:[0-9]+%\r)+")
COUNT_LINE = re.compile(rb"^  (?:Stages 1&2|Final Stage).*$", re.MULTILINE)
SUMMARY_TOTAL = re.compile(
    rb"^Relationship summary \(total relatives: 0 by pedigree, [0-9]+ by inference\)$",
    re.MULTILINE,
)
SUMMARY_ROW = re.compile(rb"^  Inference\t.*$", re.MULTILINE)


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


def normalized(stdout: bytes, omit_counts: bool = False) -> bytes:
    text = TIME.sub(b"<TIME>", PROGRESS.sub(b"", stdout))
    if omit_counts:
        text = COUNT_LINE.sub(b"<SCREEN-COUNT>", text)
        text = SUMMARY_TOTAL.sub(b"<SCREEN-SUMMARY-TOTAL>", text)
        text = SUMMARY_ROW.sub(b"<SCREEN-SUMMARY-ROW>", text)
    return text


def rows(data: bytes):
    lines = data.decode().splitlines()
    return lines[0], {tuple(line.split("\t")[:4]): line for line in lines[1:]}


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

    ref_header, ref_rows = rows(ref[3]["king.kin0"])
    impl_header, impl_rows = rows(impl[3]["king.kin0"])
    if ref_header != impl_header:
        fail("fallback king.kin0 header differs")
    missing = set(ref_rows) - set(impl_rows)
    changed = {key for key in ref_rows.keys() & impl_rows.keys() if ref_rows[key] != impl_rows[key]}
    if missing or changed:
        fail(f"shared fallback rows differ: missing={sorted(missing)}, changed={sorted(changed)}")

    expected_extras = {
        ("BF01", "B01_C2", "BF02", "B02_F"),
        ("BF13", "B13_C2", "BF14", "B14_F"),
    }
    extras = set(impl_rows) - set(ref_rows)
    if extras != expected_extras:
        fail(f"screen residual changed: expected {sorted(expected_extras)}, got {sorted(extras)}")
    if normalized(ref[1], omit_counts=True) != normalized(impl[1], omit_counts=True):
        fail("related console differs outside the two localised screening-count lines")
    if b"15 pairs of relatives" not in ref[1] or b"17 pairs of relatives" not in impl[1]:
        fail("expected 15-vs-17 screen residual was not observed")


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
    if normalized(ref[1], omit_counts=True) != normalized(impl[1], omit_counts=True):
        fail(f"{analysis} console differs outside the localized screen summary")
    if b"15 by inference" not in ref[1] or b"17 by inference" not in impl[1]:
        fail(f"{analysis} did not preserve the expected 15-vs-17 screen residual")
    if b"\t0\t0\t3\t12\t0\t0" not in ref[1]:
        fail(f"{analysis} reference fallback labels changed")
    if b"\t0\t0\t3\t14\t0\t0" not in impl[1]:
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

    print("PASS --related: 574-line .kin exact; 15 shared .kin0 rows exact")
    print("KNOWN SCREEN RESIDUAL --related: 2 extra candidates (17 vs 15)")
    print("PASS --ibdseg: console and splitped exact")
    print("PASS --unrelated: both selection files exact")
    print("PASS --cluster: updateids exact; no segment-only cluster.kin")
    print("PASS --build: build.log, updateids and updateparents exact")
    print("KNOWN SCREEN RESIDUAL clustering summaries: 17 vs 15 (3 FS exact, 14 vs 12 2nd)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
