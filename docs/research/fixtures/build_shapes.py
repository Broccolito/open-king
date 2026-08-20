#!/usr/bin/env python3
"""Held-out `--build` merge shapes, and the scorecard for `<prefix>updateparents.txt`.

`avfs.py` builds the shapes that make the reference *log* `INFERENCE AV.FS`.  This file
builds the shapes that make it *reconstruct*, which is a different question: each builder
below gets two or more declared families merged by clustering for a **different reason** —
undeclared full sibs between founders, between mothers, or between children; a
parent-offspring link; a duplicate; a half-sib link that is too weak to merge at all;
sibships that already declare a couple; several independent merges in one run.  Running
the reference over the set is what pinned every clause of the writer in
`crates/open-king-cli/src/analysis/build.rs`, none of which is visible in `bigish` alone.

    python3 build_shapes.py                 # score our binary against the reference
    python3 build_shapes.py --dump fs_kids  # show one shape's reference outputs

Nothing here reads KING's source and nothing is committed as data: the filesets land in
`work/buildshapes/` (gitignored) and are regenerated on demand.

# The scorecard, at the time of writing

20 shapes.  Two are **out of scope**: when a `.fam` names a parent living in another
family the reference materialises a phantom for it and then, the id no longer being
unique, renames every individual to `<FID>-><IID>`; this binary implements no such
renaming, so those two fail on `updateids.txt` first and are skipped.

Of the remaining 18, **15 are byte-identical** on `updateparents.txt` *and* on the console
tail (`Update-parent information is saved…` versus `No pedigrees can be reconstructed.`).
The other three — `mixed_fs_po`, `dup_plus_fs`, `fs_po_dup` — carry the **identical**
`(IID, FATHER, MOTHER)` rows and differ only in which cluster is called `KING1`: those
three are the shapes that expose the separate cluster-numbering bug recorded in
`build.rs`'s module doc, where the reference numbers merged clusters by the relationship
type that joined them (`Dup/MZ`, then `PO`, then `FS`) rather than by family order.
"""

import importlib.util
import os
import shutil
import subprocess
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))

import importlib.util, os, subprocess, sys

ROOT = os.path.normpath(os.path.join(_HERE, "..", "..", ".."))
KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
IMPL = os.environ.get("OPEN_KING", os.path.join(ROOT, "target", "release", "open-king"))
OUT = os.path.join(_HERE, "work", "buildshapes")

def gc():
    spec = importlib.util.spec_from_file_location(
        "gc", os.path.join(ROOT, "tests", "parity", "generate_corpus.py"))
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m)
    return m

G = gc()

def pad(ped, n):
    for k in range(n):
        ped.add("SG%03d" % (k + 1), "SF%03d" % (k + 1), sex=1 + (k % 2))

# --- shapes -----------------------------------------------------------------

def s_fs_fathers(na=3, nb=3):
    """Two nuclear families whose fathers are undeclared full sibs (bigish's shape)."""
    ped = G.Ped(); ph = G.add_couple(ped, "PH", "PH", emit=False)
    G.add_nuclear(ped, "FA", "A", na, father_parents=ph)
    G.add_nuclear(ped, "FB", "B", nb, father_parents=ph)
    return ped

def s_fs_mothers(na=3, nb=3):
    """Same, but the undeclared full sibs are the two mothers."""
    ped = G.Ped(); ph = G.add_couple(ped, "PH", "PH", emit=False)
    G.add_nuclear(ped, "FA", "A", na, mother_parents=ph)
    G.add_nuclear(ped, "FB", "B", nb, mother_parents=ph)
    return ped

def shadow(ped, who, tag):
    """A non-emitted genotype-identical copy of `who`, so a child of the shadow is
    genotypically the child of `who` while the .fam declares no parents at all."""
    return ped.add(tag, "SH", sex=ped.get(who).sex, clone_of=who, emit=False)

def s_fs_kids(na=3, nb=3):
    """One child of family B is an undeclared full sib of family A's children."""
    ped = G.Ped()
    fa, ma, _ = G.add_nuclear(ped, "FA", "A", na)
    sf, sm = shadow(ped, fa, "SH_F"), shadow(ped, ma, "SH_M")
    fb, mb, _ = G.add_nuclear(ped, "FB", "B", nb)
    ped.add("B_X", "FB", father=sf, mother=sm, sex=1)
    return ped

def s_po_across(na=3, nb=3):
    """A_F is the undeclared father of B_F: a PO link across two families."""
    ped = G.Ped()
    fa, ma, _ = G.add_nuclear(ped, "FA", "A", na)
    sf = shadow(ped, fa, "SH_F")
    pm = ped.add("SH_M", "SH", sex=2, emit=False)
    ped.add("B_F", "FB", father=sf, mother=pm, sex=1)
    mb = ped.add("B_M", "FB", sex=2)
    for k in range(nb):
        ped.add("B_C%d" % (k + 1), "FB", father="B_F", mother=mb, sex=1 + (k % 2))
    return ped

def s_po_singleton(na=3):
    """One extra genotyped child of family A, declared as its own one-person family."""
    ped = G.Ped()
    fa, ma, _ = G.add_nuclear(ped, "FA", "A", na)
    sf, sm = shadow(ped, fa, "SH_F"), shadow(ped, ma, "SH_M")
    ped.add("Z_1", "FZ", father=sf, mother=sm, sex=2)
    return ped

def s_phantom_parents(na=3, nb=3):
    """B_F declares FA's couple as its parents, so PLINK materialises two ungenotyped
    phantoms inside FB. Held out for the phantom branch."""
    ped = G.Ped()
    fa, ma, _ = G.add_nuclear(ped, "FA", "A", na)
    ped.add("BB_F", "FB", father=fa, mother=ma, sex=1)
    mb = ped.add("BB_M", "FB", sex=2)
    for k in range(nb):
        ped.add("BB_C%d" % (k + 1), "FB", father="BB_F", mother=mb, sex=1 + (k % 2))
    return ped

def s_fs_singletons(n=3):
    """`n` singleton families that are really one undeclared sibship."""
    ped = G.Ped(); ph = G.add_couple(ped, "PH", "PH", emit=False)
    for k in range(n):
        ped.add("Q%d" % (k + 1), "FQ%d" % (k + 1), father=ph[0], mother=ph[1],
                sex=1 + (k % 2))
    return ped

def s_hs_across(na=3, nb=3):
    """Two families whose fathers are undeclared HALF sibs (shared phantom father)."""
    ped = G.Ped()
    gf = ped.add("PH_F", "PH", sex=1, emit=False)
    g1 = ped.add("PH_M1", "PH", sex=2, emit=False)
    g2 = ped.add("PH_M2", "PH", sex=2, emit=False)
    G.add_nuclear(ped, "FA", "A", na, father_parents=(gf, g1))
    G.add_nuclear(ped, "FB", "B", nb, father_parents=(gf, g2))
    return ped

def s_dup_across(na=3, nb=3):
    """Two families sharing one individual, entered twice under different ids."""
    ped = G.Ped()
    fa, ma, kids = G.add_nuclear(ped, "FA", "A", na)
    G.add_nuclear(ped, "FB", "B", nb)
    ped.add("B_D", "FB", clone_of=kids[0], sex=1)
    return ped

def s_three_fs(nfam=3, nkids=2):
    """`nfam` families whose fathers are one undeclared sibship."""
    ped = G.Ped(); ph = G.add_couple(ped, "PH", "PH", emit=False)
    for f in range(nfam):
        t = chr(ord("A") + f)
        G.add_nuclear(ped, "F" + t, t, nkids, father_parents=ph)
    return ped

def s_fs_declared_parents(na=3, nb=3):
    """The merged full sibs already DECLARE parents in the .fam (a genotyped couple),
    so the FS0 rule has nothing to invent."""
    ped = G.Ped()
    gf, gm = G.add_couple(ped, "FA", "A_G")
    fa = ped.add("A_F", "FA", father=gf, mother=gm, sex=1)
    ma = ped.add("A_M", "FA", sex=2)
    for k in range(na):
        ped.add("A_C%d" % (k + 1), "FA", father=fa, mother=ma, sex=1 + (k % 2))
    fb = ped.add("B_F", "FB", father=gf, mother=gm, sex=1)
    mb = ped.add("B_M", "FB", sex=2)
    for k in range(nb):
        ped.add("B_C%d" % (k + 1), "FB", father=fb, mother=mb, sex=1 + (k % 2))
    return ped

def s_two_groups(na=2, nb=2, nc=2, nd=2):
    """Two INDEPENDENT merges in one run: (A,B) and (C,D). Pins the synthetic-id
    counter's scope -- per run or per family."""
    ped = G.Ped()
    ph1 = G.add_couple(ped, "P1", "P1", emit=False)
    G.add_nuclear(ped, "FA", "A", na, father_parents=ph1)
    G.add_nuclear(ped, "FB", "B", nb, father_parents=ph1)
    ph2 = G.add_couple(ped, "P2", "P2", emit=False)
    G.add_nuclear(ped, "FC", "C", nc, father_parents=ph2)
    G.add_nuclear(ped, "FD", "D", nd, father_parents=ph2)
    return ped

def s_mixed_fs_po(na=2, nb=2, nc=3, nd=3):
    """One merged family from an FS pair (A,B) and a second from a PO link (C,D).
    Pins whether an unreconstructable family consumes synthetic ids or gets rows."""
    ped = G.Ped()
    ph = G.add_couple(ped, "P1", "P1", emit=False)
    G.add_nuclear(ped, "FA", "A", na, father_parents=ph)
    G.add_nuclear(ped, "FB", "B", nb, father_parents=ph)
    fc, mc, _ = G.add_nuclear(ped, "FC", "C", nc)
    sf = shadow(ped, fc, "SH_C")
    pm = ped.add("SH_CM", "SH", sex=2, emit=False)
    ped.add("D_F", "FD", father=sf, mother=pm, sex=1)
    md = ped.add("D_M", "FD", sex=2)
    for k in range(nd):
        ped.add("D_C%d" % (k + 1), "FD", father="D_F", mother=md, sex=1 + (k % 2))
    return ped

def s_two_sibships(na=2, nb=2):
    """One merged family containing TWO independent undeclared sibships: the two
    fathers are full sibs and so are the two mothers."""
    ped = G.Ped()
    p1 = G.add_couple(ped, "P1", "P1", emit=False)
    p2 = G.add_couple(ped, "P2", "P2", emit=False)
    G.add_nuclear(ped, "FA", "A", na, father_parents=p1, mother_parents=p2)
    G.add_nuclear(ped, "FB", "B", nb, father_parents=p1, mother_parents=p2)
    return ped

def s_fs_one_declared(na=3, nb=3):
    """The FS0 group's first member already declares a genotyped couple as parents;
    the other declares nothing. Does the group inherit that couple or get (1 2)?"""
    ped = G.Ped()
    gf, gm = G.add_couple(ped, "FA", "A_G")
    fa = ped.add("A_F", "FA", father=gf, mother=gm, sex=1)
    ma = ped.add("A_M", "FA", sex=2)
    for k in range(na):
        ped.add("A_C%d" % (k + 1), "FA", father=fa, mother=ma, sex=1 + (k % 2))
    sgf, sgm = shadow(ped, gf, "SH_GF"), shadow(ped, gm, "SH_GM")
    fb = ped.add("B_F", "FB", father=sgf, mother=sgm, sex=1)
    mb = ped.add("B_M", "FB", sex=2)
    for k in range(nb):
        ped.add("B_C%d" % (k + 1), "FB", father=fb, mother=mb, sex=1 + (k % 2))
    return ped

def s_two_sibships_rev(na=2, nb=2):
    """Two sibships again, but each family's MOTHER is written to the .fam first and is
    named so that it sorts AFTER the father. Separates .fam order from ID order as the
    thing that sequences synthetic parent ids."""
    ped = G.Ped()
    p1 = G.add_couple(ped, "P1", "P1", emit=False)
    p2 = G.add_couple(ped, "P2", "P2", emit=False)
    for tag, fid in (("A", "FA"), ("B", "FB")):
        mo = ped.add("%s_Z" % tag, fid, sex=2, father=p2[0], mother=p2[1])
        fa = ped.add("%s_F" % tag, fid, sex=1, father=p1[0], mother=p1[1])
        for k in range(na if tag == "A" else nb):
            ped.add("%s_C%d" % (tag, k + 1), fid, father=fa, mother=mo, sex=1 + (k % 2))
    return ped

def s_dup_plus_fs(na=2, nb=2, nc=3, nd=3):
    """A duplicate-only merge alongside an FS merge: does the unproductive family still
    get identity rows once the run reconstructs something?"""
    ped = G.Ped()
    ph = G.add_couple(ped, "P1", "P1", emit=False)
    G.add_nuclear(ped, "FA", "A", na, father_parents=ph)
    G.add_nuclear(ped, "FB", "B", nb, father_parents=ph)
    _fc, _mc, kids = G.add_nuclear(ped, "FC", "C", nc)
    G.add_nuclear(ped, "FD", "D", nd)
    ped.add("D_D", "FD", clone_of=kids[0], sex=1)
    return ped

def s_three_clusters():
    """Three merged clusters whose member counts (8, 12, 10) are in neither ascending
    nor appearance order. Discriminates cluster numbering by size from by position."""
    ped = G.Ped()
    for tag, kids in (("A", 2), ("C", 4), ("E", 3)):
        ph = G.add_couple(ped, "P" + tag, "P" + tag, emit=False)
        t2 = chr(ord(tag) + 1)
        G.add_nuclear(ped, "F" + tag, tag, kids, father_parents=ph)
        G.add_nuclear(ped, "F" + t2, t2, kids, father_parents=ph)
    return ped

def s_wide_vs_deep():
    """One cluster of THREE small families (9 members) and one of TWO larger ones (12).
    Separates 'most members' from 'most original families'."""
    ped = G.Ped()
    ph = G.add_couple(ped, "PW", "PW", emit=False)
    for t in ("A", "B", "C"):
        G.add_nuclear(ped, "F" + t, t, 1, father_parents=ph)
    ph2 = G.add_couple(ped, "PD", "PD", emit=False)
    for t in ("X", "Y"):
        G.add_nuclear(ped, "F" + t, t, 4, father_parents=ph2)
    return ped

def s_fs_po_dup():
    """Three clusters in appearance order FS, PO, Dup. If cluster numbering follows the
    relationship type rather than the family order, the Dup cluster is KING1."""
    ped = G.Ped()
    ph = G.add_couple(ped, "P1", "P1", emit=False)
    G.add_nuclear(ped, "FA", "A", 2, father_parents=ph)
    G.add_nuclear(ped, "FB", "B", 2, father_parents=ph)
    fc, _mc, _k = G.add_nuclear(ped, "FC", "C", 2)
    sf = shadow(ped, fc, "SH_C"); pm = ped.add("SH_CM", "SH", sex=2, emit=False)
    ped.add("D_F", "FD", father=sf, mother=pm, sex=1)
    md = ped.add("D_M", "FD", sex=2)
    for k in range(2):
        ped.add("D_C%d" % (k + 1), "FD", father="D_F", mother=md, sex=1 + (k % 2))
    _fe, _me, ke = G.add_nuclear(ped, "FE", "E", 2)
    G.add_nuclear(ped, "FG", "G", 2)
    ped.add("G_D", "FG", clone_of=ke[0], sex=1)
    return ped

SHAPES = {
    "fs_po_dup": s_fs_po_dup,
    "three_clusters": s_three_clusters, "wide_vs_deep": s_wide_vs_deep,
    "two_sibships_rev": s_two_sibships_rev, "dup_plus_fs": s_dup_plus_fs,
    "mixed_fs_po": s_mixed_fs_po, "two_sibships": s_two_sibships,
    "fs_one_declared": s_fs_one_declared,
    "fs_fathers": s_fs_fathers, "fs_mothers": s_fs_mothers, "fs_kids": s_fs_kids,
    "po_across": s_po_across, "po_singleton": s_po_singleton,
    "fs_singletons": s_fs_singletons, "hs_across": s_hs_across,
    "dup_across": s_dup_across, "three_fs": s_three_fs,
    "fs_declared_parents": s_fs_declared_parents, "two_groups": s_two_groups,
    "phantom_parents": s_phantom_parents,
}

def make(name, seed, nsnp=50000, total=100):
    d = os.path.join(OUT, "%s_%d" % (name, seed))
    bed = os.path.join(d, "%s_%d.bed" % (name, seed))
    if not os.path.exists(os.path.join(d, "kingbuild.log")):
        ped = SHAPES[name]()
        pad(ped, max(0, total - ped.n_emitted()))
        spec = G.Spec("%s_%d" % (name, seed), ped, G.AUTOSOMES, nsnp, notes="build shape")
        os.makedirs(d, exist_ok=True)
        G.simulate(spec, seed, d)
        subprocess.run([KING, "-b", bed, "--build", "--cpus", "1"], cwd=d,
                       check=True, capture_output=True)
    return d, bed


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

def console_tail(text):
    """The three lines of the reconstruction tail whose choice this rig is about."""
    return "\n".join(l for l in text.splitlines() if l.startswith(
        ("Update-ID information", "Update-parent information", "No pedigrees can be")))


def run_impl(d, bed):
    o = os.path.join(d, "impl")
    shutil.rmtree(o, ignore_errors=True)
    os.makedirs(o)
    r = subprocess.run([IMPL, "-b", bed, "--build", "--cpus", "1"], cwd=o,
                       check=True, capture_output=True, text=True)
    return o, r.stdout


def score(names):
    ok = bad = skipped = 0
    for name in names:
        d, bed = make(name, SEED)
        want = open(os.path.join(d, "kingupdateparents.txt")).read()
        if "->" in want:
            print("%-20s SKIP  reference renamed ids (unimplemented)" % name)
            skipped += 1
            continue
        o, stdout = run_impl(d, bed)
        got = open(os.path.join(o, "kingupdateparents.txt")).read()
        ref = subprocess.run([KING, "-b", bed, "--build", "--cpus", "1"], cwd=d,
                             check=True, capture_output=True, text=True).stdout
        tw, tg = console_tail(ref), console_tail(stdout)
        if want == got and tw == tg:
            print("%-20s OK" % name)
            ok += 1
            continue
        bad += 1
        rows = lambda t: {tuple(l.split("\t")[1:]) for l in t.splitlines(1)}
        same = rows(want) == rows(got)
        print("%-20s MISMATCH  parents=%s tail=%s%s"
              % (name, want == got, tw == tg,
                 "  (same rows, different KING<k> label)" if same and want != got else ""))
    print("\n%d OK, %d MISMATCH, %d skipped" % (ok, bad, skipped))


def dump(names):
    for name in names:
        d, _bed = make(name, SEED)
        print("=== %s ===" % name)
        for f in ("kingupdateids.txt", "kingupdateparents.txt", "kingbuild.log"):
            p = os.path.join(d, f)
            print("-- %s --\n%s" % (f, open(p).read().rstrip() if os.path.exists(p)
                                     else "(absent)"))


SEED = 4242

if __name__ == "__main__":
    args = sys.argv[1:]
    os.makedirs(OUT, exist_ok=True)
    if args and args[0] == "--dump":
        dump(args[1:] or sorted(SHAPES))
    else:
        score(args or sorted(SHAPES))
