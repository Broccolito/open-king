#!/usr/bin/env python3
"""Held-out pedigree shapes that raise `INFERENCE AV.FS` in the reference's `--build` log.

These are the fixtures behind the `Join3/Join2` finding recorded in
`crates/king-cli/src/analysis/build.rs` and `docs/PARITY.md` §6.2: nuclear families whose
*fathers* are undeclared full sibs, so clustering merges them and pedigree reconstruction
has to decide avuncular-vs-grandparent for each (R; N1, N2) triple. Padding singletons push
the sample count past the hundred-sample clustering gate, which is what makes the reference
log rules at all.

    python3 avfs.py two   <name> <kids_a> <kids_b> <singletons> <nsnp> <outdir> <seed>
    python3 avfs.py multi <name> <nfam> <kids> <singletons> <outdir> <seed>

`two` builds two sibships (the shape `bigish` has); `multi` builds `nfam` of them sharing
one phantom father-couple, which is how the "one synthetic parent pair per mutually-full-sib
group regardless of group size" rule was pinned.

Writes a PLINK fileset into `<outdir>`; run the reference's `--build` on it and read
`kingbuild.log`. Nothing here is committed as data — `docs/research/fixtures/work/` is
gitignored and these filesets are regenerated on demand.

**The scorer is preserved now** — losing it once cost a whole round. `segprobe/` is a
standalone crate that dumps each pair's reported IBD intervals through
`king_core::ibdseg::Scan` and intersects them; `avscore.py` drives it over every
`kingbuild.log` any rig has left behind and scores `Join3/Join2` at `%.3lf`:

    cd segprobe && cargo build --release
    python3 avscore.py 1 work/avfs work/clusternum work/battery ...

The `1` is the interval-set variant, and it is the answer: **reported** segments (IBD2
calls plus the IBD1 pieces that clear `--seglength` once the IBD2 calls are cut out), which
scores 296 of 297. Variant 0, the raw calls, is the reading that was wrong and scores 13.
"""

import importlib.util
import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
_GC = os.path.join(_HERE, "..", "..", "..", "tests", "parity", "generate_corpus.py")


def _corpus():
    spec = importlib.util.spec_from_file_location("gc", os.path.normpath(_GC))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _emit(gc, name, ped, nsnp, outdir, seed, notes):
    spec = gc.Spec(name, ped, gc.AUTOSOMES, nsnp, notes=notes)
    os.makedirs(outdir, exist_ok=True)
    gc.simulate(spec, seed, outdir)
    return os.path.join(outdir, name + ".bed")


def _pad(gc, ped, nsing):
    for k in range(nsing):
        ped.add("SG%03d" % (k + 1), "SF%03d" % (k + 1), sex=1 + (k % 2))


def two(name, na, nb, nsing, nsnp, outdir, seed):
    """Two nuclear families whose fathers are undeclared full sibs."""
    gc = _corpus()
    ped = gc.Ped()
    phantom = gc.add_couple(ped, "PH", "PH", emit=False)
    gc.add_nuclear(ped, "FA", "A", na, father_parents=phantom)
    gc.add_nuclear(ped, "FB", "B", nb, father_parents=phantom)
    _pad(gc, ped, nsing)
    return _emit(gc, name, ped, nsnp, outdir, seed, "AV.FS probe")


def multi(name, nfam, nkids, nsing, outdir, seed):
    """`nfam` nuclear families sharing one undeclared full-sib father group."""
    gc = _corpus()
    ped = gc.Ped()
    phantom = gc.add_couple(ped, "PH", "PH", emit=False)
    for f in range(nfam):
        tag = chr(ord("A") + f)
        gc.add_nuclear(ped, "F" + tag, tag, nkids, father_parents=phantom)
    _pad(gc, ped, nsing)
    return _emit(gc, name, ped, 50000, outdir, seed, "n-father sibship probe")


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("two", "multi"):
        sys.exit(__doc__)
    mode, rest = sys.argv[1], sys.argv[2:]
    if mode == "two":
        name, na, nb, nsing, nsnp, outdir, seed = rest
        print(two(name, int(na), int(nb), int(nsing), int(nsnp), outdir, int(seed)))
    else:
        name, nfam, nkids, nsing, outdir, seed = rest
        print(multi(name, int(nfam), int(nkids), int(nsing), outdir, int(seed)))
