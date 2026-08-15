#!/usr/bin/env python3
"""Instruments for `--build`'s `<prefix>build.log`, the pedigree-reconstruction log.

`avfs.py` builds the shapes that make the reference log `INFERENCE AV.FS`; `avfs_score.py`
scores the `Join3/Join2` statistic; `build_shapes.py` scores `updateparents.txt`.  This
file is the fourth piece: it scores the **log itself**, and it carries the two rigs that
closed the candidate space for the one rule that is still open — *which two members of a
sibship an `AV.FS` line names*.

Nothing here reads KING's source.  Everything is black-box: build a fileset, run the
reference, read its log.  Filesets land in `work/buildlog/` (gitignored).

    python3 buildlog.py rules     # our RULE/header lines vs the reference's, all shapes
    python3 buildlog.py order     # the internal sibship order, over seeds and sizes
    python3 buildlog.py pairs     # the named pair vs every pairwise statistic

# What each section is evidence for

**`rules`** replays every shape in `build_shapes.py` plus the `avfs.py` two- and
multi-family shapes and compares the lines this binary writes — `Family KING<k>:`,
`RULE FS0`, `RULE FS1` — with the reference's, after dropping the lines we do not build.
23 of 30 match byte for byte; the seven that do not are three cluster-numbering
differences, two `<FID>-><IID>` renaming shapes that are out of scope, and two whose only
difference is the sibship member *order* that `order` below is about.

**`order`** reads that order straight off the reference's `RULE FS1: X joins in sibship
(…)` line, which prints the whole list.  Run over four seeds at each of three sibship
sizes it prints the same order every time while the sibship's own kinships move by 0.06 —
i.e. under complete genotype reseeding — which is what rules genotypes out as an input.

**`pairs`** tabulates, for every `AV.FS` line, the named pair's rank inside its sibship
under each column `--related` prints, over the sibships of three or more children (a
two-child sibship has one candidate pair and carries no information).  No column's
`argmin` or `argmax` picks the named pair more than **11 of 27** times, against a
1-in-3-or-worse chance baseline.  Between them, `order` and `pairs` are why `crates/king-cli/src/analysis/
build.rs` writes the rule half of the log and stops: the ordering is not genotype-derived,
not `.fam` order, not the sample index, not the sibship's size or position, and not any
pairwise statistic — and inventing one would be a fitted fiction.
"""

from __future__ import annotations

import importlib.util
import itertools
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
WORK = os.path.join(_HERE, "work", "buildlog")

KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
)
IMPL = os.environ.get("OPEN_KING", os.path.join(_ROOT, "target", "release", "king"))

AV = re.compile(
    r"INFERENCE AV\.FS: (\S+) is ([\w, ]+?) of (\S+) and (\S+), Join3/Join2=([\d.]+)"
)
FS1 = re.compile(r"RULE FS1: (\S+) joins in sibship \(([^)]*)\)")
COLS = ["HetHet", "IBS0", "HetConc", "HomIBS0", "Kinship", "IBD1Seg", "IBD2Seg", "PropIBD"]


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _run(binary, bed, cwd, *flags):
    os.makedirs(cwd, exist_ok=True)
    subprocess.run([binary, "-b", bed, "--cpus", "1", *flags], cwd=cwd,
                   check=True, capture_output=True)


def _kin(d, bed):
    """`--related`'s within-family rows, keyed by the unordered pair of IIDs.

    Empty when the fileset has no within-family pair at all, which is what a run of
    singletons produces: the reference writes no `.kin` in that case.
    """
    p = os.path.join(d, "king.kin")
    if not os.path.exists(p):
        _run(KING, bed, d, "--related")
    if not os.path.exists(p):
        return {}
    rows = {}
    with open(p) as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows[frozenset((f[1], f[2]))] = dict(zip(hdr, f))
    return rows


def _sibships(fam_path):
    """The declared sibships of a `.fam`, as `{(father, mother): [iid, …]}`."""
    out = {}
    for line in open(fam_path):
        f = line.split()
        if f[2] != "0":
            out.setdefault((f[2], f[3]), []).append(f[1])
    return {k: v for k, v in out.items() if len(v) > 1}


# ---------------------------------------------------------------------------
# rules: our header/RULE lines against the reference's
# ---------------------------------------------------------------------------

def _rule_lines(text):
    """The log lines this binary builds: headers and `RULE FS*`, nothing else."""
    keep = []
    for ln in text.splitlines():
        s = ln.strip()
        if not s or "INFERENCE" in ln or s.startswith(("HS ", "Reconstruct", "Duplicate")):
            continue
        keep.append(ln)
    return keep


def rules(dirs):
    ok = bad = 0
    for d in dirs:
        name = os.path.basename(d.rstrip("/"))
        bed = os.path.join(d, name + ".bed")
        if not os.path.exists(bed) or not os.path.exists(os.path.join(d, "kingbuild.log")):
            continue
        o = os.path.join(d, "impl")
        _run(IMPL, bed, o, "--build")
        ref = _rule_lines(open(os.path.join(d, "kingbuild.log")).read())
        ours = open(os.path.join(o, "kingbuild.log")).read().splitlines()
        same = ref == ours
        ok, bad = ok + same, bad + (not same)
        print("%-6s %-28s ref=%2d ours=%2d" % ("OK" if same else "DIFF", name,
                                               len(ref), len(ours)))
        if not same:
            for a, b in itertools.zip_longest(ref, ours, fillvalue=""):
                if a != b:
                    print("        ref : %r\n        ours: %r" % (a, b))
    print("\n%d match, %d differ" % (ok, bad))


# ---------------------------------------------------------------------------
# order: the sibship's internal member order, read off the `RULE FS1` line
# ---------------------------------------------------------------------------

def order(sizes=(3, 4, 5), seeds=(11, 4242, 777, 90210)):
    bs = _load("bs", os.path.join(_HERE, "build_shapes.py"))
    seen = {}
    for na in sizes:
        for seed in seeds:
            d = os.path.join(WORK, "ord%d_%d" % (na, seed))
            bed = os.path.join(d, os.path.basename(d) + ".bed")
            if not os.path.exists(os.path.join(d, "kingbuild.log")):
                ped = bs.s_fs_kids(na, 3)
                bs.pad(ped, max(0, 100 - ped.n_emitted()))
                spec = bs.G.Spec(os.path.basename(d), ped, bs.G.AUTOSOMES, 50000,
                                 notes="sibship order probe")
                os.makedirs(d, exist_ok=True)
                bs.G.simulate(spec, seed, d)
                _run(KING, bed, d, "--build")
            log = open(os.path.join(d, "kingbuild.log")).read()
            m = FS1.search(log)
            # The sibship's own kinships, as the proof that the genotypes really did
            # change: the order above must not move while these do.
            rows = _kin(d, bed)
            sibs = ["A_C%d" % (k + 1) for k in range(na)]
            ks = sorted(float(rows[frozenset(p2)]["Kinship"])
                        for p2 in itertools.combinations(sibs, 2)
                        if frozenset(p2) in rows)
            print("size=%d seed=%-6d order=%-28s kinships %.4f…%.4f"
                  % (na, seed, " ".join(m.group(2).split()) if m else "(none)",
                     ks[0], ks[-1]))
            if m:
                seen.setdefault(na, set()).add(m.group(2))
    print()
    for na, orders in sorted(seen.items()):
        print("size %d: %d distinct order(s) over %d seeds -> %s"
              % (na, len(orders), len(seeds), " | ".join(sorted(orders))))


# ---------------------------------------------------------------------------
# pairs: the named pair's rank under every printed pairwise statistic
# ---------------------------------------------------------------------------

def pairs(dirs):
    obs = []
    for d in dirs:
        name = os.path.basename(d.rstrip("/"))
        bed = os.path.join(d, name + ".bed")
        if not os.path.exists(bed) or not os.path.exists(os.path.join(d, "kingbuild.log")):
            continue
        sibs = _sibships(bed[:-4] + ".fam")
        par = {}
        for key, group in sibs.items():
            for iid in group:
                par[iid] = key
        rows, seen = _kin(d, bed), set()
        if not rows:
            continue
        for _R, _v, n1, n2, _val in AV.findall(open(os.path.join(d, "kingbuild.log")).read()):
            key = par.get(n1)
            # A two-child sibship has one candidate pair and so ranks it first and last
            # at once; it carries no information about the rule and is left out.
            if key is None or len(sibs[key]) < 3 or tuple(sibs[key]) in seen:
                continue
            seen.add(tuple(sibs[key]))
            obs.append((name, sibs[key], (n1, n2), rows))
    print("%-8s %-24s %-16s %s"
          % ("set", "sibship", "named", " ".join("%9s" % c for c in COLS)))
    hits = {c: [0, 0] for c in COLS}
    for name, group, named, rows in obs:
        cand = list(itertools.combinations(group, 2))
        cells = []
        for c in COLS:
            vals = sorted((float(rows[frozenset(p)][c]), p) for p in cand)
            rank = next(i for i, (_x, p) in enumerate(vals) if set(p) == set(named)) + 1
            cells.append("%4d/%-4d" % (rank, len(cand)))
            hits[c][0] += rank == 1
            hits[c][1] += rank == len(cand)
        print("%-8s %-24s %-16s %s"
              % (name, ",".join(x.split("_")[-1] for x in group),
                 "%s,%s" % named, " ".join(cells)))
    print("\n%d sibships" % len(obs))
    for c in COLS:
        print("  %-9s argmin picks it %2d/%2d, argmax %2d/%2d"
              % (c, hits[c][0], len(obs), hits[c][1], len(obs)))


# The 30 shapes `rules` and `pairs` score: `build_shapes.py`'s twenty merge shapes and
# ten `avfs.py` ones whose sibships are big enough to raise several `AV.FS` lines.
AVFS_SHAPES = [("two", 2, 2), ("two", 2, 3), ("two", 2, 6), ("two", 3, 3), ("two", 3, 4),
               ("two", 4, 4), ("two", 4, 6), ("two", 5, 5), ("multi", 3, 4), ("multi", 4, 3)]


def _shape_dirs():
    out = []
    bs = os.path.join(_HERE, "work", "buildshapes")
    if os.path.isdir(bs):
        out += [os.path.join(bs, d) for d in sorted(os.listdir(bs))]
    avfs = _load("avfs", os.path.join(_HERE, "avfs.py"))
    for kind, a, b in AVFS_SHAPES:
        name = "%s%d%d" % (kind[0], a, b)
        d = os.path.join(WORK, name)
        bed = os.path.join(d, name + ".bed")
        if not os.path.exists(bed):
            if kind == "two":
                avfs.two(name, a, b, 90, 50000, d, 7)
            else:
                avfs.multi(name, a, b, max(4, 100 - a * (b + 2)), d, 11)
        if not os.path.exists(os.path.join(d, "kingbuild.log")):
            _run(KING, bed, d, "--build")
        out.append(d)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    dirs = sys.argv[2:] or (_shape_dirs() if mode in ("rules", "pairs") else [])
    if mode == "rules":
        rules(dirs)
    elif mode == "order":
        order()
    elif mode == "pairs":
        pairs(dirs)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
