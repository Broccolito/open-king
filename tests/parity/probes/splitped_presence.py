#!/usr/bin/env python3
"""Differential family-size sweep for conditional splitped generation."""

import argparse
import re
import subprocess
import tempfile
from pathlib import Path


TIME = re.compile(rb"[A-Z][a-z]{2} [A-Z][a-z]{2} [ 0-9][0-9] [0-9:]{8} [0-9]{4}")
PROGRESS = re.compile(rb"(?:[0-9]+%\r)+")


def fixtures(golden: Path, root: Path):
    original_fam = [line.split() for line in (golden / "multifam.fam").read_text().splitlines()]
    rows = [line.split() for line in (golden / "multifam.bim").read_text().splitlines()[:2001]]
    bed = (golden / "multifam.bed").read_bytes()
    bytes_per_variant = (len(original_fam) + 3) // 4
    body = bed[: 3 + len(rows) * bytes_per_variant]
    for index, row in enumerate(rows):
        row[0] = "1"
        row[1] = f"v{index}"
        row[2] = "0"
        row[3] = str(1_000_000 + index * 50_000)

    result = {}
    for maximum in (1, 2, 3):
        prefix = root / f"max{maximum}"
        fam = []
        for index, source in enumerate(original_fam):
            fid = "GROUP" if index < maximum else f"F{index:02}"
            fam.append([fid, source[1], "0", "0", source[4], source[5]])
        prefix.with_suffix(".fam").write_text("".join(" ".join(row) + "\n" for row in fam))
        prefix.with_suffix(".bim").write_text("".join(" ".join(row) + "\n" for row in rows))
        prefix.with_suffix(".bed").write_bytes(body)
        result[maximum] = prefix
    return result


def run(binary: Path, prefix: Path, analysis: str, work: Path):
    proc = subprocess.run(
        [str(binary), "-b", str(prefix.with_suffix(".bed")), analysis, "--cpus", "1"],
        cwd=work,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    stdout = TIME.sub(b"<TIME>", PROGRESS.sub(b"", proc.stdout))
    files = {path.name: path.read_bytes() for path in work.iterdir() if path.is_file()}
    return proc.returncode, stdout, proc.stderr, files


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()
    ref_binary = args.ref.resolve()
    impl_binary = args.impl.resolve()
    golden = Path(__file__).resolve().parents[1] / "golden"

    with tempfile.TemporaryDirectory(prefix="open-king-splitped-") as tmp:
        root = Path(tmp)
        for maximum, prefix in fixtures(golden, root).items():
            ref_work = root / f"ref-{maximum}"
            impl_work = root / f"impl-{maximum}"
            ref_work.mkdir()
            impl_work.mkdir()
            reference = run(ref_binary, prefix, "--ibdseg", ref_work)
            implementation = run(impl_binary, prefix, "--ibdseg", impl_work)
            if reference != implementation:
                print(f"FAIL maximum-family-size={maximum}")
                print(f"reference process={reference[:3]!r}")
                print(f"implementation process={implementation[:3]!r}")
                print(f"reference files={sorted(reference[3])}")
                print(f"implementation files={sorted(implementation[3])}")
                return 1
            present = "kingsplitped.txt" in reference[3]
            if present != (maximum >= 2):
                print(f"FAIL maximum-family-size={maximum}: unexpected presence={present}")
                return 1

            # No other supported analysis owns this artefact.
            for label, binary in (("ref", ref_binary), ("impl", impl_binary)):
                work = root / f"{label}-related-{maximum}"
                work.mkdir()
                related = run(binary, prefix, "--related", work)
                if "kingsplitped.txt" in related[3] or b"splitped.txt is generated" in related[1]:
                    print(f"FAIL {label} --related maximum-family-size={maximum}")
                    return 1
            print(f"PASS maximum-family-size={maximum} splitped={present}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
