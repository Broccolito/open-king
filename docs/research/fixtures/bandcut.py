#!/usr/bin/env python3
"""Walk `Join3/Join2` through the AV.FS verdict boundary and read the verdict off the log.

The knob is genotype surgery on the *named* sib pair: over the first `f` of the markers,
overwrite N2's calls with N1's.  Inside that stretch the pair is IBD everywhere, so it is
added to `Join3` wherever `R` already shared with N1 — the ratio walks smoothly from its
natural avuncular value up to 1 while the pedigree, the .fam and every id stay put.

    python3 bandcut.py sweep <srcdir> <R> <N1> <N2> f0 f1 step
    python3 bandcut.py bisect <srcdir> <R> <N1> <N2> lo hi   # find a verdict change
"""
import os, re, shutil, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
KING = os.environ.get("KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
SP = os.path.join(HERE, "segprobe", "target", "release", "segprobe")
OUT = os.path.join(HERE, "work", "bandcut")

AV = re.compile(r"INFERENCE AV\.FS: (\S+) is ([\w, ]+?) of (\S+) and (\S+), Join3/Join2=([\d.]+)")


def bed_dims(fam, bim):
    return sum(1 for _ in open(fam)), sum(1 for _ in open(bim))


def graft(srcdir, name, dst, i_src, i_dst, frac):
    """Copy the fileset, overwriting sample `i_dst`'s calls with `i_src`'s on the first
    `frac` of markers."""
    os.makedirs(dst, exist_ok=True)
    for ext in (".bim", ".fam"):
        shutil.copy(os.path.join(srcdir, name + ext), os.path.join(dst, name + ext))
    n, m = bed_dims(os.path.join(dst, name + ".fam"), os.path.join(dst, name + ".bim"))
    per = (n + 3) // 4
    raw = bytearray(open(os.path.join(srcdir, name + ".bed"), "rb").read())
    upto = int(m * frac)
    for v in range(upto):
        base = 3 + v * per
        b_src = (raw[base + i_src // 4] >> (2 * (i_src % 4))) & 3
        off = base + i_dst // 4
        sh = 2 * (i_dst % 4)
        raw[off] = (raw[off] & ~(3 << sh)) | (b_src << sh)
    open(os.path.join(dst, name + ".bed"), "wb").write(bytes(raw))


def sample_index(fam, iid):
    for k, line in enumerate(open(fam)):
        if line.split()[1] == iid:
            return k
    raise KeyError(iid)


def one(srcdir, name, R, N1, N2, frac):
    d = os.path.join(OUT, "%s_%s_%06d" % (name, N2, round(frac * 1e5)))
    bed = os.path.join(d, name + ".bed")
    if not os.path.exists(os.path.join(d, "kingbuild.log")):
        fam = os.path.join(srcdir, name + ".fam")
        graft(srcdir, name, d, sample_index(fam, N1), sample_index(fam, N2), frac)
        r = subprocess.run([KING, "-b", bed, "--cpus", "1", "--build"], cwd=d,
                           capture_output=True, text=True)
        if "FATAL" in r.stdout or r.returncode not in (0, 1):
            return frac, None, None, r.stdout.strip().splitlines()[-1:]
    log = open(os.path.join(d, "kingbuild.log")).read()
    verdict = None
    for m in AV.findall(log):
        if m[0] == R and {m[2], m[3]} == {N1, N2}:
            verdict = (m[1], m[4])
    res = subprocess.run([SP, "join", bed, "%s,%s,%s" % (R, N1, N2)],
                         capture_output=True, text=True, check=True)
    ours = float(res.stdout.split("ratio=")[1].split("\t")[0])
    return frac, verdict, ours, None


def sweep(srcdir, R, N1, N2, f0, f1, step):
    name = os.path.basename(srcdir.rstrip("/"))
    f = f0
    while f <= f1 + 1e-12:
        frac, verdict, ours, err = one(srcdir, name, R, N1, N2, f)
        print("f=%.5f  ours=%s  ref=%s" %
              (frac, "%.4f" % ours if ours is not None else "?",
               ("%-38s %s" % verdict) if verdict else ("SILENT" if not err else err)))
        sys.stdout.flush()
        f += step


if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "sweep":
        sweep(sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5],
              float(sys.argv[6]), float(sys.argv[7]), float(sys.argv[8]))
