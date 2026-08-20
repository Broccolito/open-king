"""The row-level `.seg` scorecard, measured from the built binary against the goldens.

`run_parity.py` grades a *case*: 480 captured reference invocations, each pass/fail on
byte equality of every output file. That is the headline, and it is deliberately brutal —
one wrong digit in one row of a 763-row file fails the whole case.

This script grades the *rows underneath* those cases, so the two numbers can be read
together without either misleading. It replays `--ibdseg` at each of the three floors the
corpus captures (`--seglength` 3 — the default — 5 and 10) over the ten datasets that
carry `.seg` goldens for all three, and reports, per floor:

    rows      reference rows the pair filter agrees on
    exact     rows where **all four** printed columns match (IBD1Seg, IBD2Seg,
              PropIBD, InfType)
    ibd1      rows where IBD1Seg matches
    ibd2      rows where IBD2Seg matches
    extra     pairs we report that the reference does not
    missing   pairs the reference reports that we do not
    MAE       mean |PropIBD - reference PropIBD| over the graded rows
    worst     the largest such difference

`docs/PARITY.md` §4.4 quotes this table. It reads the goldens directly and needs no
reference binary.

    python3 scorecard.py [path/to/open-king/open-king]        # the table
    python3 scorecard.py --per-dataset                   # split by dataset
    python3 scorecard.py --residual                      # every non-exact row
"""

import os
import subprocess
import sys
import tempfile

import kingdata as kd

#: (floor in Mb, the golden-case suffix that captured it).
FLOORS = [(3, "__ibdseg"), (5, "__ibdseg_seglength5"), (10, "__ibdseg_seglength10")]

IMPL = os.path.join(kd.ROOT, "target", "release", "open-king")


def read_seg(path):
    """`{(fid1, id1, fid2, id2): (ibd1, ibd2, prop, inftype)}` from a `.seg` file."""
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        c = {n: k for k, n in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            key = (f[c["FID1"]], f[c["ID1"]], f[c["FID2"]], f[c["ID2"]])
            out[key] = (f[c["IBD1Seg"]], f[c["IBD2Seg"]],
                        f[c["PropIBD"]], f[c["InfType"]])
    return out


def golden(name, suffix):
    return os.path.join(kd.ROOT, "tests", "parity", "golden", "ibdseg",
                        name + suffix, "king.seg")


def ours(impl, name, mb, tmp):
    """Run the binary at floor `mb` and return its `.seg` rows."""
    args = [impl, "-b", os.path.join(kd.DATA, name + ".bed"), "--ibdseg", "--cpus", "1"]
    if mb != 3:
        args += ["--seglength", str(mb)]
    subprocess.run(args, cwd=tmp, check=True, capture_output=True)
    return read_seg(os.path.join(tmp, "king.seg"))


def score(impl, mb, suffix, datasets, residual=False):
    tot = dict(rows=0, exact=0, ibd1=0, ibd2=0, extra=0, missing=0)
    err = 0.0
    worst = 0.0
    by = {}
    bad = []
    for name in datasets:
        ref = read_seg(golden(name, suffix))
        with tempfile.TemporaryDirectory() as tmp:
            got = ours(impl, name, mb, tmp)
        e = r = 0
        for key in got:
            if key not in ref:
                tot["extra"] += 1
        for key, want in ref.items():
            g = got.get(key)
            if g is None:
                tot["missing"] += 1
                continue
            tot["rows"] += 1
            r += 1
            ok1 = g[0] == want[0]
            ok2 = g[1] == want[1]
            tot["ibd1"] += ok1
            tot["ibd2"] += ok2
            if ok1 and ok2 and g[2] == want[2] and g[3] == want[3]:
                tot["exact"] += 1
                e += 1
            elif residual:
                bad.append((name, key, want, g))
            d = abs(float(g[2]) - float(want[2]))
            err += d
            worst = max(worst, d)
        by[name] = (e, r)
    tot["mae"] = err / tot["rows"] if tot["rows"] else 0.0
    tot["worst"] = worst
    tot["by"] = by
    tot["bad"] = bad
    return tot


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    impl = args[0] if args else IMPL
    per_ds = "--per-dataset" in sys.argv
    residual = "--residual" in sys.argv
    # Only the datasets that carry a golden at all three floors, so the three columns
    # grade the same rows.
    datasets = [d for d in kd.DATASETS
                if all(os.path.exists(golden(d, s)) for _, s in FLOORS)]

    print("impl: %s" % impl)
    print("datasets (%d): %s\n" % (len(datasets), " ".join(datasets)))
    print("  floor   rows  exact   ibd1   ibd2  extra  missing        MAE    worst")
    rows = []
    for mb, suffix in FLOORS:
        s = score(impl, mb, suffix, datasets, residual)
        rows.append((mb, s))
        print("%5d Mb  %5d  %5d  %5d  %5d  %5d  %7d  %9.6f  %7.4f"
              % (mb, s["rows"], s["exact"], s["ibd1"], s["ibd2"],
                 s["extra"], s["missing"], s["mae"], s["worst"]))
    if per_ds:
        print("\nexact / graded rows, by dataset")
        print("  %-13s %s" % ("dataset", "  ".join("%9d Mb" % mb for mb, _ in FLOORS)))
        for name in datasets:
            cells = []
            for _mb, s in rows:
                e, r = s["by"][name]
                cells.append("%5d /%4d" % (e, r))
            print("  %-13s %s" % (name, "  ".join(cells)))
    if residual:
        for mb, s in rows:
            if not s["bad"]:
                continue
            print("\n%d Mb — %d non-exact row(s)" % (mb, len(s["bad"])))
            for name, key, want, g in s["bad"]:
                print("  %-12s %-10s %-10s want %s got %s"
                      % (name, key[1], key[3], " ".join(want), " ".join(g)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
