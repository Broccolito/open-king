#!/usr/bin/env python3
"""Pin the two formerly conflated held-out segment residuals against KING 2.3.2.

The safe rule is that the >10 Mb pair filter reads the conditioned merged calls from
both passes.  The remaining value residual is not a caller rule: KING reads beyond an
exact-multiple-of-64 marker array.  Removing or appending one marker makes it disappear.
"""

from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
FIXTURES = ROOT / "docs" / "research" / "fixtures"
sys.path.insert(0, str(FIXTURES))

import segcanvas as canvas  # noqa: E402


def fail(message: str) -> None:
    raise SystemExit("FAIL: " + message)


def load_oosseg():
    path = FIXTURES / "oosseg.py"
    spec = importlib.util.spec_from_file_location("oosseg_probe", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def run(binary: Path, bed: Path, floor: int, work: Path) -> bytes:
    work.mkdir(parents=True)
    subprocess.run(
        [str(binary), "-b", str(bed), "--ibdseg", "--cpus", "1", "--seglength", str(floor)],
        cwd=work,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        check=False,
    )
    path = work / "king.seg"
    return path.read_bytes() if path.exists() else b""


def row(blob: bytes, ids: tuple[str, str, str, str]) -> list[str] | None:
    for line in blob.decode().splitlines()[1:]:
        fields = line.split()
        if tuple(fields[:4]) == ids:
            return fields
    return None


def resize_map(source: Path, destination: Path, marker_count: int) -> Path:
    fam = source.with_suffix(".fam").read_bytes()
    bim = source.with_suffix(".bim").read_text().splitlines()
    bed = source.with_suffix(".bed").read_bytes()
    samples = sum(bool(line.strip()) for line in fam.splitlines())
    bytes_per_variant = (samples + 3) // 4
    if len(bim) != 40_000 or len(bed) != 3 + len(bim) * bytes_per_variant:
        fail("unexpected twofam fixture dimensions")

    if marker_count == 39_999:
        bim = bim[:-1]
        bed = bed[:-bytes_per_variant]
    elif marker_count == 40_001:
        fields = bim[-1].split()
        fields[1] += "_missing_pad"
        fields[3] = str(int(fields[3]) + 1)
        bim.append("\t".join(fields))
        bed += bytes([0x55]) * bytes_per_variant  # PLINK missing call for every sample
    elif marker_count != 40_000:
        fail(f"unsupported marker count {marker_count}")

    destination.with_suffix(".fam").write_bytes(fam)
    destination.with_suffix(".bim").write_text("\n".join(bim) + "\n")
    destination.with_suffix(".bed").write_bytes(bed)
    return destination.with_suffix(".bed")


def check_merged_ibd2(ref: Path, impl: Path, root: Path) -> None:
    clean = canvas.CLEAN
    interruption = {"ibs0": 1, "hethet": 63}
    fixture = canvas.Canvas(
        "merged_pair_filter",
        [clean] * 4 + [interruption] * 2 + [clean] * 4,
        nw1=1,
        sp1=30_000,
        pad=(3, 3),
        nw2=60,
        spacing=30_000,
    )
    source = root / "merged-ibd2-source"
    source.mkdir(parents=True)
    bed = Path(fixture.build(str(source)) + ".bed")
    ids = ("F00", "S00", "F01", "S01")
    for floor in (3, 5, 10):
        a = run(ref, bed, floor, root / f"merged-ref-{floor}")
        b = run(impl, bed, floor, root / f"merged-impl-{floor}")
        if a != b:
            fail(f"merged-IBD2 fixture differs at {floor} Mb")
        present = row(a, ids) is not None
        if present != (floor >= 5):
            fail(f"merged-IBD2 pair presence changed at {floor} Mb")


def check_merged_ibd1(ref: Path, impl: Path, root: Path, oos) -> None:
    name = "twofam31415926"
    bed = Path(oos._build(oos._corpus(), str(root), name, 31_415_926, "twofam"))
    ids = ("SF013", "SG013", "SF024", "SG024")
    for floor in (3, 5, 10):
        a = run(ref, bed, floor, root / f"ibd1-ref-{floor}")
        b = run(impl, bed, floor, root / f"ibd1-impl-{floor}")
        if a != b:
            fail(f"held-out merged-IBD1 file differs at {floor} Mb")
        present = row(a, ids) is not None
        if present != (floor >= 5):
            fail(f"merged-IBD1 pair presence changed at {floor} Mb")


def check_safe_tail_divergence(ref: Path, impl: Path, root: Path, oos) -> None:
    ids = ("FA", "A_C2", "FA", "A_C3")
    differences = 0
    for seed in (13_572_468, 20_260_814):
        name = f"twofam{seed}"
        source = Path(oos._build(oos._corpus(), str(root), name, seed, "twofam")).with_suffix("")
        for marker_count in (39_999, 40_000, 40_001):
            prefix = root / f"{name}_{marker_count}"
            bed = resize_map(source, prefix, marker_count)
            for floor in (3, 5):
                a = row(run(ref, bed, floor, root / f"tail-ref-{seed}-{marker_count}-{floor}"), ids)
                b = row(run(impl, bed, floor, root / f"tail-impl-{seed}-{marker_count}-{floor}"), ids)
                if a is None or b is None:
                    fail("tail-boundary target row disappeared")
                if marker_count != 40_000:
                    if a != b:
                        fail(f"safe {marker_count}-marker control differs")
                else:
                    if a == b or a[4] != b[4] or a[5] == b[5]:
                        fail("exact-64 divergence no longer has the localized IBD2 shape")
                    differences += 1
    if differences != 4:
        fail(f"expected four exact-64 value divergences, found {differences}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref", type=Path, required=True)
    parser.add_argument("--impl", type=Path, required=True)
    args = parser.parse_args()
    ref, impl = args.ref.resolve(), args.impl.resolve()
    oos = load_oosseg()
    with tempfile.TemporaryDirectory(prefix="segment-residuals-") as tmp:
        root = Path(tmp)
        check_merged_ibd2(ref, impl, root / "ibd2")
        check_merged_ibd1(ref, impl, root / "ibd1", oos)
        check_safe_tail_divergence(ref, impl, root / "tail", oos)
    print("PASS merged IBD1/IBD2 calls feed the >10 Mb pair filter")
    print("PASS 39,999/40,001-marker controls are exact")
    print("EXPECTED SAFE DIVERGENCE four IBD2 rows only at exactly 40,000 markers")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
