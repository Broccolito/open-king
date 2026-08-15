#!/usr/bin/env python3
"""Trigger probe: for every (sibship, R) the reference *could* have inferred over,
print the 2nd-degree evidence and every Join3/Join2, beside the line it actually printed.
"""
import os, re, subprocess, sys, itertools

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
KING = os.environ.get("KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
SP = os.path.join(HERE, "segprobe", "target", "release", "segprobe")

AV = re.compile(r"INFERENCE AV\.FS: (\S+) is ([\w, ]+?) of (\S+) and (\S+), Join3/Join2=([\d.]+)")


def run(binary, bed, cwd, *flags):
    os.makedirs(cwd, exist_ok=True)
    subprocess.run([binary, "-b", bed, "--cpus", "1", *flags], cwd=cwd, check=True,
                   capture_output=True)


def kin_all(d, bed, degree=2):
    """{frozenset(pair): InfType} over .kin and .kin0 at the given degree."""
    sub = os.path.join(d, "kin%d" % degree)
    p0, p1 = os.path.join(sub, "king.kin0"), os.path.join(sub, "king.kin")
    if not os.path.exists(p0) and not os.path.exists(p1):
        run(KING, bed, sub, "--related", "--degree", str(degree))
    out = {}
    for p in (p1, p0):
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            hdr = next(fh).rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                out[frozenset((f["ID1"], f["ID2"]))] = (f["InfType"], float(f["PropIBD"]),
                                                        float(f["Kinship"]))
    return out


def fam(bed):
    rows = [l.split() for l in open(bed[:-4] + ".fam")]
    return rows


def sibships(rows):
    out = {}
    for r in rows:
        if r[2] != "0":
            out.setdefault((r[0], r[2], r[3]), []).append(r[1])
    return {k: v for k, v in out.items() if len(v) > 1}


def join(bed, triples):
    if not triples:
        return {}
    res = subprocess.run([SP, "join", bed] + ["%s,%s,%s" % t for t in triples],
                         check=True, capture_output=True, text=True)
    out = {}
    for line in res.stdout.splitlines():
        f = line.split("\t")
        t = tuple(f[0].split(","))
        out[t] = (int(f[1].split("=")[1]), int(f[2].split("=")[1]), float(f[3].split("=")[1]))
    return out


def report(d):
    name = os.path.basename(d.rstrip("/"))
    bed = os.path.join(d, name + ".bed")
    log = os.path.join(d, "kingbuild.log")
    if not os.path.exists(bed):
        return
    if not os.path.exists(log):
        run(KING, bed, d, "--build")
    text = open(log).read()
    rows = fam(bed)
    sibs = sibships(rows)
    rel = kin_all(d, bed)
    everyone = [r[1] for r in rows]
    print("=" * 78)
    print(name)
    print("-- reference log --")
    print(text.rstrip("\n") if text.strip() else "(empty)")
    printed = {(m[0], frozenset((m[2], m[3]))): (m[1], m[4]) for m in AV.findall(text)}
    print("-- candidates --")
    triples, meta = [], []
    for key, group in sorted(sibs.items()):
        for r in everyone:
            if r in group:
                continue
            second = [s for s in group if rel.get(frozenset((r, s)), ("UN",))[0] == "2nd"]
            if len(second) < 2:
                continue
            for a, b in itertools.combinations(group, 2):
                triples.append((r, a, b))
                meta.append((key, group, r, second, a, b))
    vals = join(bed, triples)
    for key, group, r, second, a, b in meta:
        j2, j3, ratio = vals[(r, a, b)]
        hit = printed.get((r, frozenset((a, b))))
        print("  sib=%-22s R=%-8s 2nd=%d/%d  pair=(%s,%s) ratio=%.4f %s"
              % (",".join(group), r, len(second), len(group), a, b, ratio,
                 ("<== PRINTED %s %s" % hit) if hit else ""))


if __name__ == "__main__":
    for d in sys.argv[1:]:
        report(d)
