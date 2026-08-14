"""Recover the reference's **per-segment** IBD2 lengths by bisecting `--seglength`.

`--seglength` takes a float in Mb and is honoured to the base pair (`--seglength 6.250001`
prints `Minimum segment length is set as 6250001 bp`), and a segment is reported iff its
length is `>= ` that value. `IBD2Seg(L)` is therefore a monotone step function of `L` whose
jumps sit exactly one base pair above each called IBD2 segment's own length — so an
adaptive bisection over `L` turns an aggregate column into the **multiset of segment
lengths**, per pair, exactly.

That is the instrument this investigation was missing. `MaxIBD2` grades one segment per
pair; this grades every one of them.

    python3 seglen_probe.py nuclear missing            # write work/seglen/<ds>.json
    python3 seglen_probe.py --col IBD1Seg nuclear      # the same for the IBD1 column

Each dataset costs a few thousand reference invocations at ~20 ms each. `bigish` is far
larger than the rest; run it alone and expect an hour.
"""

import json
import os
import subprocess
import sys
import tempfile

import kingdata as kd

KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
OUT = os.path.join(kd.ROOT, "tests", "parity", "work", "seglen")


class Probe:
    def __init__(self, name, col="IBD2Seg"):
        self.ds = kd.load(name)
        self.name = name
        self.col = col
        self.idx = {(f, i): k for k, (f, i) in enumerate(self.ds.fam)}
        self.tmp = tempfile.mkdtemp(prefix="seglen_")
        self.calls = 0
        self.cache = {}

    def values(self, bp):
        """`{(i, j): printed column}` at `--seglength bp`, absent rows reading 0.0000."""
        v = self.cache.get(bp)
        if v is not None:
            return v
        self.calls += 1
        subprocess.run(
            [KING, "-b", os.path.join(kd.DATA, self.name + ".bed"), "--ibdseg",
             "--cpus", "1", "--seglength", "%.6f" % (bp / 1e6), "--prefix", "p"],
            cwd=self.tmp, check=True, capture_output=True)
        out = {}
        path = os.path.join(self.tmp, "p.seg")
        if os.path.exists(path):
            with open(path) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                c = head.index(self.col)
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    i, j = self.idx[(f[0], f[1])], self.idx[(f[2], f[3])]
                    out[(min(i, j), max(i, j))] = f[c]
        v = {k: out.get(k, "0.0000") for k in self.ds.ref}
        self.cache[bp] = v
        return v

    def jumps(self, lo, hi):
        """Every `L` at which some pair's column changes, to the base pair."""
        found = []
        stack = [(lo, self.values(lo), hi, self.values(hi))]
        while stack:
            a, va, b, vb = stack.pop()
            if va == vb:
                continue
            if b - a <= 1:
                found.append(b)
                continue
            m = (a + b) // 2
            vm = self.values(m)
            stack.append((a, va, m, vm))
            stack.append((m, vm, b, vb))
        return sorted(found)

    def lengths(self, lo=3_000_000, hi=300_000_000):
        """`{(i, j): [segment lengths]}` — a jump at `L` means a segment of `L - 1` bp."""
        js = self.jumps(lo, hi)
        out = {k: [] for k in self.ds.ref}
        prev = self.values(lo)
        for j in js:
            cur = self.values(j)
            for k in out:
                if cur[k] != prev[k]:
                    out[k].append(j - 1)
            prev = cur
        return out, js


def main():
    argv = sys.argv[1:]
    col = "IBD2Seg"
    if "--col" in argv:
        k = argv.index("--col")
        col = argv[k + 1]
        del argv[k:k + 2]
    # `--seglength` above 10 Mb is silently ignored by the reference (the
    # "Minimum segment length is set as" line disappears and the default 3 Mb is used),
    # so 10 Mb is the hard ceiling of this instrument.
    hi, lo = 10_000_000, 100_000
    if "--hi" in argv:
        k = argv.index("--hi")
        hi = int(argv[k + 1])
        del argv[k:k + 2]
    if "--lo" in argv:
        k = argv.index("--lo")
        lo = int(argv[k + 1])
        del argv[k:k + 2]
    os.makedirs(OUT, exist_ok=True)
    for name in argv:
        p = Probe(name, col)
        lens, js = p.lengths(lo=lo, hi=hi)
        n = sum(len(v) for v in lens.values())
        path = os.path.join(OUT, "%s.%s.json" % (name, col))
        with open(path, "w") as fh:
            json.dump({"%d,%d" % k: v for k, v in lens.items()}, fh)
        print("%-12s %5d reference calls, %4d jump points, %4d segment lengths -> %s"
              % (name, p.calls, len(js), n, path))


if __name__ == "__main__":
    main()
