#!/usr/bin/env python3
"""What orders a sibship's members inside `RULE FS1: X joins in sibship (…)`.

`buildlog.py order` established what the order is *not*: not genotype-derived (four fresh
seeds give the same order), not the `.fam` row order (permuting the rows moves the same
individuals), not the absolute sample index, not the sibship's size or position, and not
any pairwise statistic.  What was left standing, unmeasured, is that the order is a
function of the **id strings** — the signature of iteration over a hash-keyed container.
This file tests exactly that, and then bounds the container.

**What it found.**

* `names` — the order **is** a function of the ids.  Eight id sets in one fixed pedigree,
  same seed, same sample count, give three distinct position-orders: `A_C1..A_C3` and
  `A_D1..A_D3` give `(2,3,1)`, `K1..K3` / `1001..1003` / `B01_C1..B01_C3` the identity,
  and `zeta,alpha,mu` gives `(3,2,1)`.
* `subsets` — it is **not** a per-id ranking.  Over thirteen subsets of one eight-id pool
  the pairwise precedences contradict each other **91** times, so no `sort by f(id)`
  reproduces them.  The order moves when the *set* changes, which is what a hash table's
  capacity does, and it is why `{P1,P2,P3,P4}`, `{P1,P3,P5,P7}`, `{P2,P4,P6,P8}` and
  `{P3,P4,P5,P6}` all come out permuted `(2,1,4,3)` while `A_C1..A_C4` comes out
  `(3,4,2,1)`.
* `popn`, `setsize` — the container is **not** global.  Nine total sample counts from 102
  to 142, and four padding vocabularies at a fixed count, all leave the order
  byte-identical.
* `family` — it **is** scoped to the individual's own family.  Renaming the sibship's
  parents changes the kids' order; renaming the joiner or the other family's children does
  not.
* `perm` — four `.fam` permutations of the same three ids give one order, so the container
  is keyed by the strings, not by position.

So it is an iteration order over a family-scoped, id-keyed, capacity-sensitive container.
Reproducing it means identifying the hash, which is the next piece of work; inventing one
would be a fitted fiction.

    python3 siborder.py names     # same structure, different child ids
    python3 siborder.py popn      # same ids, different TOTAL sample count
    python3 siborder.py setsize   # same ids, other ids changed around them
    python3 siborder.py perm      # the id set fixed, its .fam order permuted
    python3 siborder.py subsets   # is the order a per-id total order?
    python3 siborder.py family    # do ids outside the sibship but inside its family matter?
    python3 siborder.py sizes     # sibship sizes 2…7 under one id vocabulary

Every run reads the order straight off the reference's own `kingbuild.log`, so nothing is
inferred from our binary.  Filesets land in `work/siborder/` (gitignored).
"""

from __future__ import annotations

import importlib.util
import os
import re
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
WORK = os.path.join(_HERE, "work", "siborder")

KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"
)

FS1 = re.compile(r"RULE FS1: (\S+) joins in sibship \(([^)]*)\)")
FS0 = re.compile(r"RULE FS0: Sibship \(([^)]*)\)")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


G = _load("gc", os.path.join(_ROOT, "tests", "parity", "generate_corpus.py"))


def kids_shape(kid_ids, pad_ids, joiner="B_X", nb=3, pad_prefix="SG"):
    """`build_shapes.s_fs_kids` with the child ids and the padding under our control.

    Family FA holds `kid_ids` as one declared sibship; family FB holds `nb` children of
    its own plus `joiner`, who is genotypically a full sib of FA's children but declares
    no parents.  Reconstruction therefore joins `joiner` to FA's sibship and prints the
    sibship's members in whatever order it holds them.
    """
    ped = G.Ped()
    fa = ped.add("A_F", "FA", sex=1)
    ma = ped.add("A_M", "FA", sex=2)
    for k, iid in enumerate(kid_ids):
        ped.add(iid, "FA", father=fa, mother=ma, sex=1 + (k % 2))
    sf = ped.add("SH_F", "SH", sex=1, clone_of=fa, emit=False)
    sm = ped.add("SH_M", "SH", sex=2, clone_of=ma, emit=False)
    fb = ped.add("B_F", "FB", sex=1)
    mb = ped.add("B_M", "FB", sex=2)
    for k in range(nb):
        ped.add("B_C%d" % (k + 1), "FB", father=fb, mother=mb, sex=1 + (k % 2))
    ped.add(joiner, "FB", father=sf, mother=sm, sex=1)
    for iid in pad_ids:
        ped.add(iid, "S" + iid, sex=1)
    return ped


def run(tag, ped, seed=4242, nsnp=30000):
    d = os.path.join(WORK, tag)
    bed = os.path.join(d, tag + ".bed")
    if not os.path.exists(os.path.join(d, "kingbuild.log")):
        spec = G.Spec(tag, ped, G.AUTOSOMES, nsnp, notes="sibship order probe")
        os.makedirs(d, exist_ok=True)
        G.simulate(spec, seed, d)
        subprocess.run([KING, "-b", bed, "--build", "--cpus", "1"], cwd=d,
                       check=True, capture_output=True)
    log = open(os.path.join(d, "kingbuild.log")).read()
    m = FS1.search(log)
    if m:
        return m.group(2).split()
    m = FS0.search(log)
    return ["FS0:"] + m.group(1).split() if m else []


def _pad(n, prefix="SG"):
    return ["%s%03d" % (prefix, k + 1) for k in range(n)]


def names():
    """Same pedigree, same seed, same sample count — only the child ids change."""
    sets = {
        "A_C":    ["A_C1", "A_C2", "A_C3"],
        "A_D":    ["A_D1", "A_D2", "A_D3"],
        "B01_C":  ["B01_C1", "B01_C2", "B01_C3"],
        "B13_C":  ["B13_C1", "B13_C2", "B13_C3"],
        "bare":   ["K1", "K2", "K3"],
        "digits": ["1001", "1002", "1003"],
        "long":   ["AAAAAAAAAAAA1", "AAAAAAAAAAAA2", "AAAAAAAAAAAA3"],
        "mixed":  ["zeta", "alpha", "mu"],
    }
    print("%-8s %-42s %s" % ("set", "ids (.fam order)", "sibship order"))
    seen = set()
    for tag, ids in sets.items():
        ped = kids_shape(ids, _pad(91))
        got = run("nm_" + tag, ped)
        # express the order as positions in the .fam order, so different id sets compare
        pos = tuple(ids.index(x) + 1 for x in got if x in ids)
        seen.add(pos)
        print("%-8s %-42s %-28s  positions %s"
              % (tag, " ".join(ids), " ".join(got), pos))
    print("\n%d distinct position-order(s) over %d id sets" % (len(seen), len(sets)))


def popn():
    """Same child ids, different TOTAL sample counts: does the container's size matter?"""
    ids = ["A_C1", "A_C2", "A_C3"]
    print("%-6s %-8s %s" % ("total", "padding", "sibship order"))
    seen = {}
    for extra in (0, 1, 2, 3, 5, 8, 13, 21, 40):
        ped = kids_shape(ids, _pad(91 + extra))
        got = run("pop_%d" % extra, ped)
        pos = tuple(ids.index(x) + 1 for x in got if x in ids)
        seen.setdefault(pos, []).append(ped.n_emitted())
        print("%-6d %-8d %-28s positions %s" % (ped.n_emitted(), 91 + extra,
                                                " ".join(got), pos))
    print("\n%d distinct order(s):" % len(seen))
    for pos, totals in seen.items():
        print("   %s at totals %s" % (pos, totals))


def setsize():
    """Same child ids, padding ids CHANGED but the count kept: does content matter?"""
    ids = ["A_C1", "A_C2", "A_C3"]
    print("%-10s %s" % ("padding", "sibship order"))
    seen = set()
    for tag, prefix in (("SG", "SG"), ("ZZ", "ZZ"), ("Q", "Q"), ("PAD", "PAD")):
        ped = kids_shape(ids, _pad(91, prefix))
        got = run("set_" + tag, ped)
        pos = tuple(ids.index(x) + 1 for x in got if x in ids)
        seen.add(pos)
        print("%-10s %-28s positions %s" % (prefix, " ".join(got), pos))
    print("\n%d distinct order(s) over 4 padding vocabularies" % len(seen))


def perm():
    """The same id set, written to the .fam in different orders."""
    base = ["A_C1", "A_C2", "A_C3"]
    orders = [base, [base[2], base[0], base[1]], [base[1], base[2], base[0]],
              list(reversed(base))]
    print("%-26s %s" % (".fam order", "sibship order"))
    seen = set()
    for k, ids in enumerate(orders):
        ped = kids_shape(ids, _pad(91))
        got = run("perm_%d" % k, ped)
        seen.add(tuple(got))
        print("%-26s %s" % (" ".join(ids), " ".join(got)))
    print("\n%d distinct order(s) over 4 .fam permutations" % len(seen))


def sizes():
    """Sibship sizes 2…7 under one id vocabulary, for the `fit` section to chew on."""
    print("%-5s %s" % ("n", "sibship order"))
    out = {}
    for n in range(2, 8):
        ids = ["A_C%d" % (k + 1) for k in range(n)]
        ped = kids_shape(ids, _pad(94 - n))
        got = run("sz_%d" % n, ped)
        out[n] = [x for x in got if x in ids]
        print("%-5d %s" % (n, " ".join(got)))
    return out


POOL = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8"]


def subsets():
    """Is the order a per-id total order, or does it move when the set changes?

    A container that sorts by a key computed from the id alone gives orders that are
    *consistent*: if `a` precedes `b` in one sibship it precedes it in every sibship
    holding both.  A hash table whose capacity follows the entry count does not — change
    the set and the buckets move.  This runs the full pool and then subsets of it and
    reports every contradiction.
    """
    obs = []
    runs = [("all", POOL)]
    runs += [("first%d" % k, POOL[:k]) for k in (2, 3, 4, 5, 6, 7)]
    runs += [("odd", POOL[0::2]), ("even", POOL[1::2]),
             ("drop1", POOL[1:]), ("drop8", POOL[:-1]),
             ("mid", POOL[2:6]), ("ends", [POOL[0], POOL[3], POOL[7]])]
    for tag, ids in runs:
        got = [x for x in run("sub_" + tag, kids_shape(list(ids), _pad(94 - len(ids))))
               if x in ids]
        obs.append((tag, list(ids), got))
        print("%-8s set %-28s order %s" % (tag, " ".join(ids), " ".join(got)))
    before = {}
    bad = 0
    for tag, _ids, got in obs:
        for i, a in enumerate(got):
            for b in got[i + 1:]:
                if (b, a) in before:
                    print("  CONTRADICTION %s before %s in %s, opposite in %s"
                          % (a, b, tag, before[(b, a)]))
                    bad += 1
                before.setdefault((a, b), tag)
    print("\n%d contradiction(s): the order is %s"
          % (bad, "NOT a per-id key" if bad else "consistent with a per-id key"))
    if not bad:
        # A consistent order means one global ranking; print it.
        rank = sorted(POOL, key=lambda x: sum(1 for y in POOL if (y, x) in before))
        print("implied ranking: %s" % " ".join(rank))


def family():
    """Do ids OUTSIDE the sibship, but inside its family, move the order?"""
    ids = ["A_C1", "A_C2", "A_C3"]
    print("%-24s %s" % ("parents / joiner", "sibship order"))
    seen = set()
    for tag, fa, ma, jn in (("A_F/A_M/B_X", "A_F", "A_M", "B_X"),
                            ("ZF/ZM/B_X", "ZF", "ZM", "B_X"),
                            ("A_F/A_M/QQQ", "A_F", "A_M", "QQQ"),
                            ("W1/W2/W3", "W1", "W2", "W3")):
        ped = G.Ped()
        f = ped.add(fa, "FA", sex=1)
        m = ped.add(ma, "FA", sex=2)
        for k, iid in enumerate(ids):
            ped.add(iid, "FA", father=f, mother=m, sex=1 + (k % 2))
        sf = ped.add("SH_F", "SH", sex=1, clone_of=f, emit=False)
        sm = ped.add("SH_M", "SH", sex=2, clone_of=m, emit=False)
        fb = ped.add("B_F", "FB", sex=1)
        mb = ped.add("B_M", "FB", sex=2)
        for k in range(3):
            ped.add("B_K%d" % (k + 1), "FB", father=fb, mother=mb, sex=1 + (k % 2))
        ped.add(jn, "FB", father=sf, mother=sm, sex=1)
        for iid in _pad(91):
            ped.add(iid, "S" + iid, sex=1)
        got = run("fam_" + tag.replace("/", "_"), ped)
        seen.add(tuple(x for x in got if x in ids))
        print("%-24s %s" % (tag, " ".join(got)))
    print("\n%d distinct order(s) over 4 surroundings" % len(seen))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    fn = {"names": names, "popn": popn, "setsize": setsize, "perm": perm,
          "sizes": sizes, "subsets": subsets, "family": family}.get(mode)
    if not fn:
        sys.exit(__doc__)
    fn()
