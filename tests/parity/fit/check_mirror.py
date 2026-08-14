"""Assert that `engine.py` with default `Params` *is* the committed Rust engine.

Runs `king --ibdseg` and `king --ibs` from the built binary over every corpus dataset and
compares, row by row, the four `.seg` columns and the `MaxIBD2` column against the
mirror's own output. Any difference is a bug in `engine.py`.

**The `.seg` pass is checked at all three floors the corpus captures** — `--seglength` 3
(the default), 5 and 10 — and that is not cosmetic. The run merge of
`docs/research/20-seglength-floor.md` is floor-dependent by construction: it cannot fire
at the default floor on the corpus's marker spacings. A single-floor check therefore
cannot see it at all, and for a while did not: the merge was committed to `Scan` while
`engine.py` still had the pre-merge caller, and this script passed anyway. Anything whose
rule reads `--seglength` must be exercised away from the default or the mirror is only
asserted where the rule is dormant.

    python3 check_mirror.py [path/to/open-king/king]
"""

import os
import subprocess
import sys
import tempfile

import kingdata as kd
import engine as E

IMPL = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    kd.ROOT, "target", "release", "king")

#: The `--seglength` floors to check the `.seg` pass at, in Mb. 3 is the default and is
#: passed implicitly; 5 and 10 are where the run merge is live.
FLOORS = [3, 5, 10]


def run(ds, args, out):
    subprocess.run([IMPL, "-b", os.path.join(kd.DATA, ds + ".bed")] + args,
                   cwd=out, check=True, capture_output=True)


def check_seg(ds, name, idx, mb, tmp):
    """Compare the binary's `.seg` against the mirror at one floor. Returns (rows, bad)."""
    args = ["--ibdseg", "--cpus", "1"]
    if mb != 3:
        args += ["--seglength", str(mb)]
    run(name, args, tmp)
    min_bp = mb * 1_000_000
    n = bad = 0
    with open(os.path.join(tmp, "king.seg")) as fh:
        next(fh)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
            i, j = min(i, j), max(i, j)
            a, b, _lg, _m = E.call_pair(ds, i, j, min_bp=min_bp)
            g1, g2 = a / ds.denom, b / ds.denom
            want = (f[4], f[5], f[6], f[7])
            # `.seg` prints PropIBD from its own two printed columns; `InfType`
            # still reads the full-precision value. Both mirror `Scan`/the writer.
            got = ("%.4f" % g1, "%.4f" % g2,
                   "%.4f" % E.seg_prop_ibd(g1, g2),
                   kd.inf_type(g1, g2, g2 + g1 / 2))
            n += 1
            if want != got:
                bad += 1
                print("SEG  %2d Mb  %-12s %s/%s want %s got %s"
                      % (mb, name, f[1], f[3], want, got))
    return n, bad


def main():
    bad = 0
    for name in kd.DATASETS:
        ds = kd.load(name)
        idx = {(f, i): k for k, (f, i) in enumerate(ds.fam)}
        with tempfile.TemporaryDirectory() as tmp:
            run(name, ["--ibs", "--cpus", "1"], tmp)
            n = 0
            for mb in FLOORS:
                rows, nbad = check_seg(ds, name, idx, mb, tmp)
                n += rows
                bad += nbad
            m = 0
            for ext, wide in ((".ibs", False), (".ibs0", True)):
                p = os.path.join(tmp, "king" + ext)
                if not os.path.exists(p):
                    continue
                with open(p) as fh:
                    head = fh.readline().split()
                    if "MaxIBD2" not in head:
                        continue
                    c = head.index("MaxIBD2")
                    for line in fh:
                        f = line.split()
                        if wide:
                            i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
                        else:
                            i, j = idx[(f[0], f[1])], idx[(f[0], f[2])]
                        i, j = min(i, j), max(i, j)
                        if f[c] == "-9":
                            continue
                        got = "%.3f" % E.max_ibd2(ds, i, j)
                        m += 1
                        if got != f[c]:
                            bad += 1
                            print("MAX  %-12s %s/%s want %s got %s"
                                  % (name, f[1 if not wide else 1], f[2 if not wide else 3],
                                     f[c], got))
            print("%-12s  %4d seg rows (%s Mb), %5d MaxIBD2 values checked"
                  % (name, n, "/".join(str(x) for x in FLOORS), m))
    print("MIRROR OK" if bad == 0 else "MIRROR DIVERGES on %d values" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
