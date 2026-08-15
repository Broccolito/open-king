#!/usr/bin/env python3
"""Score the Join3/Join2 formula against every AV.FS / AV.HS line the rigs have captured."""
import os, re, subprocess, sys, collections

HERE = os.path.dirname(os.path.abspath(__file__))
SP = os.path.join(HERE, "segprobe", "target", "release", "segprobe")
AV = re.compile(r"INFERENCE AV\.FS: (\S+) is [\w, ]+? of (\S+) and (\S+), Join3/Join2=([\d.]+)")


def main(roots, variant):
    env = dict(os.environ, VARIANT=str(variant))
    exact = off = 0
    resid = []
    worst = []
    for root in roots:
        for base, _d, files in sorted(os.walk(root)):
            if "impl" in base or "kin2" in base or "kingbuild.log" not in files:
                continue
            name = os.path.basename(base)
            bed = os.path.join(base, name + ".bed")
            if not os.path.exists(bed):
                continue
            seen, trip = set(), []
            for r, n1, n2, val in AV.findall(open(os.path.join(base, "kingbuild.log")).read()):
                if (r, n1, n2) in seen:
                    continue
                seen.add((r, n1, n2))
                trip.append((r, n1, n2, float(val)))
            if not trip:
                continue
            res = subprocess.run([SP, "join", bed] + ["%s,%s,%s" % t[:3] for t in trip],
                                 capture_output=True, text=True, env=env)
            if res.returncode != 0:
                print("  !! %s: %s" % (name, res.stderr.strip()[:90]))
                continue
            for (r, n1, n2, ref), line in zip(trip, res.stdout.splitlines()):
                ours = float(line.split("ratio=")[1].split("\t")[0])
                same = abs(round(ours, 3) - ref) < 1e-9
                exact, off = exact + same, off + (not same)
                resid.append(ours - ref)
                if not same:
                    worst.append((abs(ours - ref), name, r, n1, n2, ref, ours))
    n = len(resid)
    print("variant %d: %d of %d exact at %%.3lf   mean residual %+.5f  range %+.5f … %+.5f"
          % (variant, exact, n, sum(resid) / n, min(resid), max(resid)))
    for w in sorted(worst, reverse=True)[:12]:
        print("   %-18s %-8s %-8s %-8s ref=%.3f ours=%.4f  (%+.4f)"
              % (w[1], w[2], w[3], w[4], w[5], w[6], w[6] - w[5]))


if __name__ == "__main__":
    v = int(sys.argv[1])
    main(sys.argv[2:], v)
