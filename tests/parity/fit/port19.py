"""Is the Rust engine the same rule as `seg19.py`?  Row by row, on the corpus.

The scorecards in this directory grade a *Python* model against the reference.  That is
only useful if the shipped engine computes the same thing, so this runs
`target/release/open-king --ibdseg` on every corpus fileset and compares its four printed `.seg`
columns against `seg19.call_pair` for the same pair — not against the reference, against
the model.  Any disagreement is a porting bug and is reported as such.

    python3 port19.py             # all three floors
    python3 port19.py 3           # one floor
"""

import os
import subprocess
import sys
import tempfile

import kingdata as kd
import seg19 as S19

KING = os.environ.get("OPENKING",
                      os.path.join(kd.ROOT, "target", "release", "open-king"))
FLOORS = [(3_000_000, []), (5_000_000, ["--seglength", "5"]),
          (10_000_000, ["--seglength", "10"])]


def run_ours(name, extra):
    bed = os.path.join(kd.DATA, name + ".bed")
    with tempfile.TemporaryDirectory() as wd:
        subprocess.run([KING, "-b", bed, "--ibdseg", *extra, "--prefix",
                        os.path.join(wd, "k")],
                       capture_output=True, text=True, cwd=wd, check=False)
        path = os.path.join(wd, "k.seg")
        if not os.path.exists(path):
            return {}
        out = {}
        with open(path) as fh:
            head = fh.readline().split()
            for line in fh:
                f = line.rstrip("\n").split("\t")
                r = dict(zip(head, f))
                out[(r["ID1"], r["ID2"])] = (r["IBD1Seg"], r["IBD2Seg"],
                                             r["PropIBD"], r["InfType"])
    return out


def check(min_bp, extra):
    rows = agree = 0
    bad = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        ours = run_ours(name, extra)
        for i, j in ds.pairs():
            key = (ds.fam[i][1], ds.fam[j][1])
            a, b, lg = S19.call_pair(ds, i, j, S19.R19(), min_bp)
            want_reported = lg >= 10_000_000
            got = ours.get(key)
            if not want_reported and got is None:
                continue
            rows += 1
            if got is None or not want_reported:
                bad.append((name, key, "reported?", want_reported, got is not None))
                continue
            g1, g2 = a / d, b / d
            gp = g2 + g1 / 2
            mine = ("%.4f" % g1, "%.4f" % g2, "%.4f" % gp,
                    kd.inf_type(g1, g2, gp))
            theirs = tuple(got)
            if all(abs(float(x) - float(y)) < 5e-9 for x, y in zip(mine[:3], theirs[:3])) \
                    and mine[3] == theirs[3]:
                agree += 1
            else:
                bad.append((name, key, mine, theirs))
    print("--seglength %d Mb: %d rows, engine agrees with seg19 on %d"
          % (min_bp // 1_000_000, rows, agree))
    for b in bad[:12]:
        print("   ", b)
    return rows, agree


if __name__ == "__main__":
    want = {int(a) for a in sys.argv[1:] if a.isdigit()}
    for bp, extra in FLOORS:
        if want and bp // 1_000_000 not in want:
            continue
        check(bp, extra)
