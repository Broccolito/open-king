#!/usr/bin/env python3
"""Held-out battery for the `--build` log's INFERENCE triggers.

Everything here is a *fresh* shape or a fresh seed: none of it took part in finding any
of the rules it scores.

    python3 battery.py band  <nseeds>   # 2-child sibships: verdict vs Join3/Join2
    python3 battery.py hs    <nseeds>   # which 2nd-degree cousin pairs raise HS.UN2
    python3 battery.py rep   <nseeds>   # how many times one AV.FS line repeats
    python3 battery.py fs2              # can RULE FS2 be made to fire
    python3 battery.py hsrel            # HS candidate whose parent check fails
"""
import importlib.util, itertools, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, "..", "..", ".."))
KING = os.environ.get("KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
SP = os.path.join(HERE, "segprobe", "target", "release", "segprobe")
WORK = os.path.join(HERE, "work", "battery")

AV = re.compile(r"INFERENCE AV\.FS: (\S+) is ([\w, ]+?) of (\S+) and (\S+), Join3/Join2=([\d.]+)")
UN2 = re.compile(r"INFERENCE HS\.UN2: (\S+) and (\S+) are HS")
HSLINE = re.compile(r"^    HS (\S+) unrelated to (\S+)$", re.M)


def gc():
    spec = importlib.util.spec_from_file_location(
        "gc", os.path.join(ROOT, "tests", "parity", "generate_corpus.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def run(binary, bed, cwd, *flags):
    os.makedirs(cwd, exist_ok=True)
    return subprocess.run([binary, "-b", bed, "--cpus", "1", *flags], cwd=cwd,
                          capture_output=True, text=True)


def build(name, ped, seed, nsnp=50000):
    g = gc()
    d = os.path.join(WORK, name)
    bed = os.path.join(d, name + ".bed")
    if not os.path.exists(os.path.join(d, "kingbuild.log")):
        os.makedirs(d, exist_ok=True)
        g.simulate(g.Spec(name, ped, g.AUTOSOMES, nsnp, notes="battery"), seed, d)
        run(KING, bed, d, "--build")
    return d, bed


def pad(g, ped, total=104):
    """Top the fileset up past the hundred-sample clustering gate."""
    have = sum(1 for p in ped.people if p.emit)
    for k in range(max(0, total - have)):
        ped.add("SG%03d" % (k + 1), "SF%03d" % (k + 1), sex=1 + (k % 2))


def shape_two(kids_a, kids_b, total=104):
    g = gc()
    ped = g.Ped()
    ph = g.add_couple(ped, "PH", "PH", emit=False)
    g.add_nuclear(ped, "FA", "A", kids_a, father_parents=ph)
    g.add_nuclear(ped, "FB", "B", kids_b, father_parents=ph)
    pad(g, ped, total)
    return ped


def shape_multi(nfam, kids, total=104):
    g = gc()
    ped = g.Ped()
    ph = g.add_couple(ped, "PH", "PH", emit=False)
    for f in range(nfam):
        t = chr(ord("A") + f)
        g.add_nuclear(ped, "F" + t, t, kids, father_parents=ph)
    pad(g, ped, total)
    return ped


def kin(d, bed, degree=2):
    sub = os.path.join(d, "kin%d" % degree)
    p0 = os.path.join(sub, "king.kin0")
    if not os.path.exists(p0):
        run(KING, bed, sub, "--related", "--degree", str(degree))
    out = {}
    if not os.path.exists(p0):
        return out
    with open(p0) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            f = dict(zip(hdr, line.rstrip("\n").split("\t")))
            out[frozenset((f["ID1"], f["ID2"]))] = f
    return out


def ratio(bed, r, n1, n2):
    res = subprocess.run([SP, "join", bed, "%s,%s,%s" % (r, n1, n2)],
                         capture_output=True, text=True, check=True,
                         env=dict(os.environ, VARIANT="1"))
    return float(res.stdout.split("ratio=")[1].split("\t")[0])


# ---------------------------------------------------------------- band
C1 = (0.848718, 0.851164)              # bracket on the uncle cut
C2 = (0.896895, 0.903318)              # bracket on the ambiguous cut


def band(nseeds):
    ok = bad = 0
    print("%-14s %-6s %-9s %-9s %-30s %s" % ("fileset", "sib", "ours", "ref", "verdict", "test"))
    for s in range(int(nseeds)):
        seed = 700000 + 137 * s
        name = "bnd%02d" % s
        d, bed = build(name, shape_two(2, 2), seed)
        log = open(os.path.join(d, "kingbuild.log")).read()
        got = {frozenset((m[2], m[3])): (m[1], float(m[4])) for m in AV.findall(log)}
        for tag, other in (("A", "B_F"), ("B", "A_F")):
            n1, n2 = tag + "_C1", tag + "_C2"
            o = ratio(bed, other, n1, n2)
            hit = got.get(frozenset((n1, n2)))
            if hit is None:
                # silent: the value has to be inside the dead band
                good = o >= C1[0] and o <= C2[1]
                verdict, ref = "SILENT", float("nan")
            else:
                verdict, ref = hit
                good = (o <= C1[1]) if verdict in ("uncle", "aunt") else (o >= C2[0])
            ok, bad = ok + good, bad + (not good)
            print("%-14s %-6s %-9.4f %-9s %-30s %s"
                  % (name, tag, o, "%.3f" % ref if hit else "-", verdict,
                     "ok" if good else "REFUTES"))
    print("\nband rule: %d consistent, %d refuting" % (ok, bad))


# ---------------------------------------------------------------- hs
def hs(nseeds):
    used, skipped = [], []
    for s in range(int(nseeds)):
        seed = 810000 + 211 * s
        name = "hsb%02d" % s
        d, bed = build(name, shape_two(5, 5), seed)
        log = open(os.path.join(d, "kingbuild.log")).read()
        fired = {frozenset(m) for m in UN2.findall(log)}
        for pair, row in sorted(kin(d, bed).items(), key=lambda kv: sorted(kv[0])):
            a, b = sorted(pair)
            if row["InfType"] != "2nd" or "_C" not in a or "_C" not in b:
                continue
            if a.split("_")[0] == b.split("_")[0]:
                continue
            rec = (name, a, b, float(row["PropIBD"]), float(row["Kinship"]),
                   float(row["IBD1Seg"]), float(row["IBD2Seg"]))
            (used if pair in fired else skipped).append(rec)
    for tag, rows in (("HS.UN2 fired", used), ("silent", skipped)):
        print("== %s (%d) ==" % (tag, len(rows)))
        for r in sorted(rows, key=lambda r: r[3]):
            print("   %-8s %-8s %-8s PropIBD=%.4f Kinship=%.4f IBD1=%.4f IBD2=%.4f" % r)
    if used and skipped:
        lo, hi = max(r[3] for r in skipped), min(r[3] for r in used)
        print("\nPropIBD bracket on the HS gate: (%.4f, %.4f)" % (lo, hi))
        lo1, hi1 = max(r[5] for r in skipped), min(r[5] for r in used)
        print("IBD1Seg bracket:                (%.4f, %.4f)" % (lo1, hi1))
        lok, hik = max(r[4] for r in skipped), min(r[4] for r in used)
        print("Kinship  bracket:               (%.4f, %.4f)  %s"
              % (lok, hik, "consistent" if lok < hik else "REFUTED (overlaps)"))


# ---------------------------------------------------------------- rep
def rep(nseeds):
    for s in range(int(nseeds)):
        seed = 920000 + 173 * s
        for nfam in (2, 3, 4):
            name = "rep%d_%02d" % (nfam, s)
            d, bed = build(name, shape_multi(nfam, 2), seed)
            log = open(os.path.join(d, "kingbuild.log")).read()
            lines = [l for l in log.splitlines() if "INFERENCE AV.FS" in l]
            counts = {}
            for l in lines:
                counts[l] = counts.get(l, 0) + 1
            print("%-10s fams=%d  distinct=%d  repeats=%s  blanks=%d"
                  % (name, nfam, len(counts), sorted(counts.values(), reverse=True),
                     sum(1 for l in log.splitlines() if not l.strip())))


# ---------------------------------------------------------------- fs2
def fs2():
    """Two *declared* sibships joined by an inferred FS pair — the RULE FS2 shape."""
    g = gc()
    ped = g.Ped()
    ph = g.add_couple(ped, "PH", "PH", emit=False)
    # FA: declared sibs A_C1..A_C3, children of A_F/A_M; FB likewise.
    # A_F and B_F are undeclared full sibs, so the two *parent* generations merge,
    # and A_C1 is additionally simulated as a full sib of B_C1 through a second
    # phantom couple, which puts two declared sibships in one component.
    g.add_nuclear(ped, "FA", "A", 3, father_parents=ph)
    g.add_nuclear(ped, "FB", "B", 3, father_parents=ph)
    ph2 = g.add_couple(ped, "QH", "QH", emit=False)
    ped.add("A_X1", "FA", father=ph2[0], mother=ph2[1], sex=1)
    ped.add("A_X2", "FA", father=ph2[0], mother=ph2[1], sex=2)
    ped.add("B_X1", "FB", father=ph2[0], mother=ph2[1], sex=1)
    ped.add("B_X2", "FB", father=ph2[0], mother=ph2[1], sex=2)
    pad(g, ped)
    # Declare A_X1/A_X2 as children of A_F/A_M and B_X1/B_X2 of B_F/B_M: two declared
    # sibships that the inference has to combine.
    d, bed = build("fs2probe", ped, 55555)
    print(open(os.path.join(d, "kingbuild.log")).read())


# ---------------------------------------------------------------- hsrel
def hsrel():
    """Double first cousins: both the fathers and the mothers are undeclared sibs, so
    an HS candidate's `unrelated to <mother>` check has to fail."""
    g = gc()
    ped = g.Ped()
    phf = g.add_couple(ped, "PF", "PF", emit=False)
    phm = g.add_couple(ped, "PM", "PM", emit=False)
    g.add_nuclear(ped, "FA", "A", 4, father_parents=phf, mother_parents=phm)
    g.add_nuclear(ped, "FB", "B", 4, father_parents=phf, mother_parents=phm)
    pad(g, ped)
    d, bed = build("hsrel", ped, 31337)
    print(open(os.path.join(d, "kingbuild.log")).read())
    for pair, row in sorted(kin(d, bed).items(), key=lambda kv: sorted(kv[0])):
        a, b = sorted(pair)
        if row["InfType"] in ("2nd", "FS", "PO") and a.split("_")[0] != b.split("_")[0]:
            print("   %-8s %-8s %-5s PropIBD=%s" % (a, b, row["InfType"], row["PropIBD"]))


if __name__ == "__main__":
    m = sys.argv[1]
    {"band": band, "hs": hs, "rep": rep}.get(m, lambda *_: None)(*sys.argv[2:]) \
        if m in ("band", "hs", "rep") else {"fs2": fs2, "hsrel": hsrel}[m]()
