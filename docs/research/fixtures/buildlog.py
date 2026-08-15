#!/usr/bin/env python3
"""Instruments for `--build`'s `<prefix>build.log`, the pedigree-reconstruction log.

`avfs.py` builds the shapes that make the reference log `INFERENCE AV.FS`; `avfs_score.py`
scores the `Join3/Join2` statistic; `build_shapes.py` scores `updateparents.txt`;
`clusternum.py` pins the merge queue that numbers the clusters; `dupkeep.py` pins the
duplicate rule; `siborder.py` chases the sibship order.  This file scores the **log
itself** and carries the rigs for the three questions the log alone can answer.

Nothing here reads KING's source.  Everything is black-box: build a fileset, run the
reference, read its log.  Filesets land in `work/buildlog/` (gitignored).

    python3 buildlog.py rules     # our rule-half lines vs the reference's, 59 shapes
    python3 buildlog.py blanks    # what predicts a family block's blank-line count
    python3 buildlog.py cut       # the uncle/grandparent cut on Join3/Join2, bracketed
    python3 buildlog.py order     # the internal sibship order, over seeds and sizes
    python3 buildlog.py pairs     # the named pair vs every pairwise statistic

# What each section is evidence for

**`rules`** replays 59 held-out shapes and compares the lines this binary writes — the
`Family KING<k>:` header, `Duplicate … is removed.`, `RULE FS0` and `RULE FS1` — with the
reference's, after dropping the inference half from both sides (`_INFERENCE_HALF`).  53
match byte for byte; the six that do not are two `<FID>-><IID>` renaming shapes that are
out of scope, three that differ only in the sibship member *order* that `order` is about,
and one where the unimplemented `PO.S` branch consumes a synthetic id.

**`blanks`** scores the two candidate rules for a family block's blank-line count over
every cluster whose sibships are all pairs, which is where the named pair is forced and the
candidate `R` set can be read off a degree-2 `.kin0`.  Both score 107 of 113, on different
failure sets; what separates them by hand is `three_fs` and `ord3`.  See
`crates/king-cli/src/analysis/build.rs`.

**`cut`** re-reads every `AV.FS`/`AV.HS` line under `work/` and brackets the
`uncle|aunt` / `grandfather|…` cut: 259 values, `(0.8495, 0.9005)` once the `%.3lf`
rounding is accounted for.

**`order`** reads the sibship order off the reference's `RULE FS1: X joins in sibship (…)`
line, which prints the whole list.  Four seeds at each of three sibship sizes give the same
order while the sibship's own kinships move by 0.06, which is what rules genotypes out.
`siborder.py` takes it from there.

**`pairs`** tabulates, for every `AV.FS` line, the named pair's rank inside its sibship
under each column `--related` prints.  No column's `argmin` or `argmax` picks the named
pair more than **11 of 27** times, against a 1-in-3-or-worse chance baseline.
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

#: Lines that belong to the log's **inference** half and so are dropped from both sides.
#:
#: The blank lines go because their count is a function of the inference half.  So does
#: everything about a parent-offspring pair: `Reconstruct parent-offspring pair (X, Y)...`
#: was long assumed to be a rule line, and it is not — 42 of 42 clusters that print it also
#: print an `INFERENCE` line, and a `PO` merge whose cluster has no sibship anywhere (two
#: shapes, three seeds each) prints nothing at all.  Its `PO.S` follow-ups go with it.  The
#: `Duplicate … is removed.` line is the opposite case and is **kept**: 23 of 27 runs print
#: it with no inference line in the file (`dupkeep.py`).
_INFERENCE_HALF = ("INFERENCE", "HS ", "Reconstruct parent-offspring pair",
                   "is used to determine the parent/offspring", "RULE PO.",
                   "is created as")


def _rule_lines(text):
    """The reference log reduced to the lines this binary is expected to build.

    Headers left with nothing under them are dropped too: a cluster whose only line is an
    inference-half one does not appear in our log at all.
    """
    keep = []
    for ln in text.splitlines():
        s = ln.strip()
        if not s or any(t in ln for t in _INFERENCE_HALF):
            continue
        if s.startswith("Family KING") and s.endswith(":") and keep and \
                keep[-1].strip().startswith("Family KING") and keep[-1].strip().endswith(":"):
            keep.pop()
        keep.append(ln)
    if keep and keep[-1].strip().startswith("Family KING") and keep[-1].strip().endswith(":"):
        keep.pop()
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
    # `clusternum.py`'s nineteen shapes, which is where the merge queue, the duplicate
    # removal and the `Reconstruct parent-offspring pair` line were pinned.
    cn = _load("cn", os.path.join(_HERE, "clusternum.py"))
    for name in cn.SHAPES:
        out.append(cn.make(name, 4242)[0])
    # `dupkeep.py`'s ten shapes, which pin which copy of a duplicate survives.
    dk = _load("dk", os.path.join(_HERE, "dupkeep.py"))
    for name in dk.SHAPES:
        out.append(dk.run(name, 4242)[0])
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


def cut(root=None):
    """The `uncle|aunt` vs `grandfather|…` cut on `Join3/Join2`, bracketed.

    Scans *every* reference `kingbuild.log` under `work/` — whatever any rig has left
    behind — and reports the largest `uncle|aunt` value and the smallest ambiguous one.
    The gap between them is the bracket; nothing inside it can be ruled out.
    """
    root = root or os.path.join(_HERE, "work")
    lo, hi, n = [], [], 0
    verdicts = {}
    for base, _dirs, files in os.walk(root):
        if "impl" in base or "kingbuild.log" not in files:
            continue
        for ln in open(os.path.join(base, "kingbuild.log")):
            m = re.search(r"INFERENCE AV\.(?:FS|HS): \S+ is ([\w, ]+?) of ", ln)
            v = re.search(r"Join3/Join2=([\d.]+)", ln)
            if not m or not v:
                continue
            n += 1
            val, word = float(v.group(1)), m.group(1)
            verdicts[word] = verdicts.get(word, 0) + 1
            (lo if word in ("uncle", "aunt") else hi).append(val)
    lo.sort()
    hi.sort()
    print("%d AV lines over %s" % (n, root))
    for w, c in sorted(verdicts.items(), key=lambda kv: -kv[1]):
        print("   %-42s %4d" % (w, c))
    if lo and hi:
        # The log prints `%.3lf`, so a printed p stands for a true value in
        # [p - 0.0005, p + 0.0005).  A cut c is compatible with a printed `uncle` p iff
        # c > p - 0.0005, and with a printed ambiguous q iff c < q + 0.0005.
        half = 0.0005
        print("\nlargest  uncle|aunt          %.3f  (%d values, %.3f … %.3f)"
              % (lo[-1], len(lo), lo[0], lo[-1]))
        print("smallest grandfather|HS|…    %.3f  (%d values, %.3f … %.3f)"
              % (hi[0], len(hi), hi[0], hi[-1]))
        print("bracket on the true cut  (%.4f, %.4f)   width %.4f"
              % (lo[-1] - half, hi[0] + half, hi[0] - lo[-1] + 2 * half))
        for c in (0.85, 0.875, 0.9):
            ok = all(x - half < c for x in lo) and all(c < x + half for x in hi)
            print("   cut %.3f: %s" % (c, "survives" if ok else "REFUTED"))


def blanks(root=None):
    """How many blank lines a family block carries, and what predicts the count.

    Two rules are scored, on the clusters where the count is fully computable — every
    sibship exactly two members, so the pair an `AV.FS` line would name is forced and the
    candidate `R` set can be read straight off `.kin0`:

    * **block** (the rule this file used to carry): one blank before each sibship block
      until the family prints its first inference; if it never prints one, one per block.
    * **reject**: one blank opens the section, and one more for every candidate `R` that is
      *examined and turned down* before the first line prints.

    They differ wherever a block faces no candidate at all, or faces several — `three_fs`,
    whose first sibship faces two candidate uncles and prints three blanks where `block`
    says two, and `ord3`, whose two sibships face none and prints one where `block` says
    two.
    """
    root = root or os.path.join(_HERE, "work")
    score = {"block": [0, 0], "reject": [0, 0]}
    for base, _dirs, files in sorted(os.walk(root)):
        if "impl" in base or "kingbuild.log" not in files:
            continue
        name = os.path.basename(base)
        bed = os.path.join(base, name + ".bed")
        log = open(os.path.join(base, "kingbuild.log")).read()
        ids = os.path.join(base, "kingupdateids.txt")
        if not log.strip() or not os.path.exists(bed) or not os.path.exists(ids):
            continue
        cluster_of = {}
        for line in open(ids):
            fid, _iid, key, _n = line.rstrip("\n").split("\t")
            cluster_of[fid] = key
        rows = [l.split() for l in open(bed[:-4] + ".fam")]
        sibs = {}
        for r in rows:
            if r[2] != "0":
                sibs.setdefault((r[0], r[2], r[3]), []).append(r[1])
        sibs = {k: v for k, v in sibs.items() if len(v) > 1}
        # A dedicated degree-2 run: several rigs leave a degree-1 `king.kin0` in the
        # fileset directory, and that one carries no 2nd-degree row at all.
        second = set()
        p = os.path.join(base, "kin2", "king.kin0")
        if not os.path.exists(p):
            _run(KING, bed, os.path.join(base, "kin2"), "--related", "--degree", "2")
        if not os.path.exists(p):
            continue
        with open(p) as fh:
            hdr = next(fh).rstrip("\n").split("\t")
            for line in fh:
                f = dict(zip(hdr, line.rstrip("\n").split("\t")))
                if f.get("InfType") == "2nd":
                    second.add(frozenset((f["ID1"], f["ID2"])))
        everyone = {r[1] for r in rows}
        blocks, key, cur = [], None, None
        for ln in log.splitlines():
            if ln.startswith("Family KING"):
                if cur is not None:
                    blocks.append((key, cur))
                key, cur = ln[len("Family "):].rstrip(":"), []
            elif cur is not None:
                cur.append(ln)
        if cur is not None:
            blocks.append((key, cur))
        for key, body in blocks:
            # the cluster's declared sibships in .fam order, then the founder sibship an
            # FS0 line names
            mine = [v for k, v in sibs.items() if cluster_of.get(k[0]) == key]
            mine += [rex.split() for l in body if "RULE FS0" in l
                     for rex in re.findall(r"Sibship \((.*?)\)", l)]
            if not mine or any(len(s) != 2 for s in mine):
                continue
            printed = [l for l in body if "INFERENCE AV" in l]
            first_at, rejects, seen, nblocks = None, 0, False, 0
            for s in mine:
                nblocks += 1
                cands = sorted(x for x in everyone if x not in s
                               and frozenset((x, s[0])) in second
                               and frozenset((x, s[1])) in second)
                for r in cands:
                    words = set()
                    hit = False
                    for l in printed:
                        toks = l.replace(",", " ").split()
                        if r in toks and set(s) <= set(toks):
                            hit = True
                    if hit:
                        seen = True
                        first_at = first_at or nblocks
                        break
                    rejects += 1
                if seen:
                    break
            got = sum(1 for l in body if not l.strip())
            pred = {"block": first_at or nblocks, "reject": 1 + rejects}
            for h, v in pred.items():
                score[h][got != v] += 1
            if got != pred["reject"]:
                print("%-22s %-7s blanks=%d  block=%d  reject=%d"
                      % (name, key, got, pred["block"], pred["reject"]))
    print("\nrule     right / wrong")
    for h, (a, b) in score.items():
        print("  %-7s %3d / %3d" % (h, a, b))


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    dirs = sys.argv[2:] or (_shape_dirs() if mode in ("rules", "pairs") else [])
    if mode == "rules":
        rules(dirs)
    elif mode == "order":
        order()
    elif mode == "blanks":
        blanks()
    elif mode == "cut":
        cut()
    elif mode == "pairs":
        pairs(dirs)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
