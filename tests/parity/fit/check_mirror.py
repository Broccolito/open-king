"""Assert that `engine.py` with default `Params` *is* the committed Rust engine.

Runs `king --ibdseg` and `king --ibs` from the built binary over every corpus dataset and
compares, row by row, the four `.seg` columns and the `MaxIBD2` column against the
mirror's own output. Any difference is a bug in `engine.py`.

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


def run(ds, args, out):
    subprocess.run([IMPL, "-b", os.path.join(kd.DATA, ds + ".bed")] + args,
                   cwd=out, check=True, capture_output=True)


def main():
    bad = 0
    for name in kd.DATASETS:
        ds = kd.load(name)
        idx = {(f, i): k for k, (f, i) in enumerate(ds.fam)}
        with tempfile.TemporaryDirectory() as tmp:
            run(name, ["--ibdseg", "--cpus", "1"], tmp)
            run(name, ["--ibs", "--cpus", "1"], tmp)
            seg = os.path.join(tmp, "king.seg")
            n = 0
            with open(seg) as fh:
                next(fh)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
                    i, j = min(i, j), max(i, j)
                    a, b, _lg, _m = E.call_pair(ds, i, j)
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
                        print("SEG  %-12s %s/%s want %s got %s"
                              % (name, f[1], f[3], want, got))
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
            print("%-12s  %4d seg rows, %5d MaxIBD2 values checked" % (name, n, m))
    print("MIRROR OK" if bad == 0 else "MIRROR DIVERGES on %d values" % bad)
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
