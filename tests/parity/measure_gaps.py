#!/usr/bin/env python3
"""Measure — not merely detect — every remaining parity gap.

``run_parity.py`` answers "does this case match?".  This script answers "by how
much does it not match?", which is what ``docs/PARITY.md`` quotes.  It replays the
same captured invocations through the same case discovery, and for every output
file that differs it reports:

* for a headered table: rows the reference emits, rows only one side emits, and
  per column the number of common rows that disagree plus — for numeric columns —
  the mean and worst absolute difference;
* for anything else: how many lines differ.

Rows are matched on the identifier columns (``FID``/``FID1``/``FID2``/``ID``/
``ID1``/``ID2``), so a pair only one side reports is never mistaken for a numeric
disagreement, and column error is measured only over pairs both sides report.

    python3 tests/parity/measure_gaps.py --impl target/release/open-king
    python3 tests/parity/measure_gaps.py --impl target/release/open-king --filter ibdseg/
    python3 tests/parity/measure_gaps.py --impl target/release/open-king --by-dataset king.seg

Python 3 standard library only.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_parity as RP  # noqa: E402  (same directory, deliberate)

KEY_COLUMNS = {"FID", "FID1", "FID2", "ID", "ID1", "ID2", "IID"}
# Headered tables whose leading column is not an identifier.
OTHER_HEADERS = {"Segment", "SNP", "Sample"}


def is_number(s: str) -> bool:
    try:
        float(s)
    except ValueError:
        return False
    return True


class Table:
    """A headered tab-separated output file, keyed on its identifier columns."""

    def __init__(self, text: str):
        lines = [ln for ln in text.split("\n") if ln != ""]
        self.header = lines[0].split("\t") if lines else []
        self.nkey = 0
        while self.nkey < len(self.header) and self.header[self.nkey] in KEY_COLUMNS:
            self.nkey += 1
        if self.nkey == 0:  # keyless table: the row's ordinal is its key
            self.rows = {(str(n),): ln.split("\t") for n, ln in enumerate(lines[1:])}
            self.order = [(str(n),) for n in range(len(lines) - 1)]
            return
        self.rows = {}
        self.order = []
        for ln in lines[1:]:
            f = ln.split("\t")
            k = tuple(f[: self.nkey])
            self.rows[k] = f
            self.order.append(k)

    @staticmethod
    def looks_tabular(text: str) -> bool:
        first = text.split("\n", 1)[0]
        if "\t" not in first:
            return False
        head = first.split("\t")[0]
        return head in KEY_COLUMNS or head in OTHER_HEADERS


def compare_tables(ref: Table, got: Table) -> dict:
    """Row and column deltas between two parses of the same output file."""
    common = [k for k in ref.order if k in got.rows]
    out = {
        "ref_rows": len(ref.order),
        "got_rows": len(got.order),
        "missing": len(ref.order) - len(common),  # reference reports it, we do not
        "extra": len(got.order) - len(common),  # we report it, reference does not
        "common": len(common),
        "cols": {},
        "rows_differing": 0,
    }
    ncol = min(len(ref.header), len(got.header))
    for c in range(ref.nkey, ncol):
        name = ref.header[c]
        bad, num, total, worst, worst_at = 0, True, 0.0, 0.0, None
        for k in common:
            a, b = ref.rows[k][c], got.rows[k][c]
            if a == b:
                continue
            bad += 1
            if is_number(a) and is_number(b):
                d = abs(float(a) - float(b))
                total += d
                if d > worst:
                    worst, worst_at = d, (k, a, b)
            else:
                num = False
        if bad:
            out["cols"][name] = {
                "differing": bad,
                "numeric": num,
                "mean_abs": total / bad if num else None,
                "worst": worst if num else None,
                "worst_at": worst_at,
            }
    for k in common:
        if ref.rows[k][:ncol] != got.rows[k][:ncol]:
            out["rows_differing"] += 1
    return out


def replay(case, impl: str, data: Path, alt: Path, timeout: float) -> dict:
    """Replay one captured invocation with our binary; return the files it wrote."""
    tmp = Path(tempfile.mkdtemp(prefix="measure-"))
    try:
        subprocess.run(
            case.argv(impl, data, alt),
            cwd=tmp,
            capture_output=True,
            timeout=timeout,
            check=False,
        )
        return {
            str(p.relative_to(tmp)): p.read_bytes()
            for p in sorted(tmp.rglob("*"))
            if p.is_file()
        }
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def new_slot():
    return {
        "cases": 0,
        "differing_cases": 0,
        "ref_rows": 0,
        "rows_differing": 0,
        "missing": 0,
        "extra": 0,
        "cols": defaultdict(lambda: {"differing": 0, "sum": 0.0, "worst": 0.0}),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--impl", required=True)
    ap.add_argument("--golden", type=Path, default=RP.DEFAULT_GOLDEN)
    ap.add_argument("--data", type=Path, default=RP.DEFAULT_DATA)
    ap.add_argument("--alt", type=Path, default=RP.DEFAULT_ALT)
    ap.add_argument("--filter", dest="filt", default="")
    ap.add_argument("--timeout", type=float, default=300.0)
    ap.add_argument("-q", "--quiet", action="store_true", help="totals only")
    ap.add_argument(
        "--by-dataset",
        metavar="FILENAME",
        help="also roll the totals for this output file up per dataset",
    )
    args = ap.parse_args()

    impl = RP.resolve_binary(args.impl)
    cases = RP.discover(args.golden, args.filt or None, False)
    if not cases:
        print("no cases matched", file=sys.stderr)
        return 2
    RP.ensure_inputs(cases, args.data, args.alt, False)

    per_file = defaultdict(new_slot)
    by_dataset = defaultdict(
        lambda: {"ref_rows": 0, "rows_differing": 0, "missing": 0, "extra": 0}
    )

    for case in cases:
        dataset = case.name.split("__", 1)[0]
        produced = replay(case, impl, args.data, args.alt, args.timeout)
        for rel, path in sorted(case.golden_files().items()):
            # Same exclusions run_parity.py applies: files the *reference* writes
            # non-deterministically are not diffable goldens.
            if RP.racy_reason(case, rel):
                continue
            ref_bytes = path.read_bytes()
            got_bytes = produced.get(rel)
            slot = per_file[rel]
            slot["cases"] += 1
            if got_bytes == ref_bytes:
                continue
            slot["differing_cases"] += 1
            ref_text = ref_bytes.decode("utf-8", "replace")
            if got_bytes is None:
                if not args.quiet:
                    print(f"{case.key}: {rel} NOT PRODUCED")
                continue
            got_text = got_bytes.decode("utf-8", "replace")
            if not Table.looks_tabular(ref_text):
                a, b = ref_text.split("\n"), got_text.split("\n")
                d = sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))
                slot["ref_rows"] += len(a)
                slot["rows_differing"] += d
                if not args.quiet:
                    print(f"{case.key}: {rel} differs, {d} of {len(a)} line(s)")
                continue
            r = compare_tables(Table(ref_text), Table(got_text))
            slot["ref_rows"] += r["ref_rows"]
            slot["rows_differing"] += r["rows_differing"]
            slot["missing"] += r["missing"]
            slot["extra"] += r["extra"]
            for cname, c in r["cols"].items():
                acc = slot["cols"][cname]
                acc["differing"] += c["differing"]
                if c["numeric"]:
                    acc["sum"] += c["mean_abs"] * c["differing"]
                    acc["worst"] = max(acc["worst"], c["worst"])
            if args.by_dataset == rel:
                d = by_dataset[dataset]
                for k in ("ref_rows", "rows_differing", "missing", "extra"):
                    d[k] += r[k]
            if not args.quiet:
                cols = ", ".join(
                    f"{k}:{v['differing']}"
                    + (
                        f"(mae {v['mean_abs']:.5f} max {v['worst']:.4f})"
                        if v["numeric"]
                        else ""
                    )
                    for k, v in r["cols"].items()
                )
                print(
                    f"{case.key}: {rel} {r['rows_differing']}/{r['common']} common row(s)"
                    f" differ, +{r['extra']} extra, -{r['missing']} missing | {cols}"
                )

    print("\n=== totals by output file ===")
    for fname in sorted(per_file):
        s = per_file[fname]
        if not s["differing_cases"]:
            print(f"{fname:28s} byte-identical in all {s['cases']} case(s)")
            continue
        print(
            f"{fname:28s} differs in {s['differing_cases']}/{s['cases']} case(s); "
            f"{s['rows_differing']} of {s['ref_rows']} reference row(s) differ, "
            f"+{s['extra']} extra, -{s['missing']} missing"
        )
        for cname in sorted(s["cols"], key=lambda c: -s["cols"][c]["differing"]):
            c = s["cols"][cname]
            tail = (
                f", mae {c['sum'] / c['differing']:.6f}, worst {c['worst']:.4f}"
                if c["worst"] or c["sum"]
                else ""
            )
            print(f"    {cname:12s} {c['differing']} row(s){tail}")

    if args.by_dataset:
        print(f"\n=== {args.by_dataset} by dataset ===")
        for ds in sorted(by_dataset):
            d = by_dataset[ds]
            print(
                f"{ds:14s} {d['rows_differing']}/{d['ref_rows']} row(s) differ, "
                f"+{d['extra']} extra, -{d['missing']} missing"
            )
    return 0


if __name__ == "__main__":
    sys.exit(main())
