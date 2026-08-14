#!/usr/bin/env python3
"""Build deterministic alternate .fam/.bim inputs for KING's --fam/--bim overrides.

The parity corpus (generate_corpus.py) writes <ds>.{bed,bim,fam}.  KING can be
pointed at a *different* .fam or .bim with --fam/--bim while still reading the
same .bed.  To test that path we need alternate files that are (a) byte-for-byte
reproducible and (b) provoke an observable change in KING's output.

Nothing here uses randomness: every alternate file is a pure function of the
corresponding corpus file, so re-running this script always yields identical
bytes.

Variants produced per dataset
-----------------------------
<ds>.altfam.fam   Same sample order and IIDs.  FID -> "AF<n>" (every sample in
                  its own family), parents zeroed, sex flipped (1<->2, 0 kept),
                  phenotype preserved.  Effect: all pedigree structure is
                  destroyed, so every pair must migrate from .kin to .kin0.

<ds>.altbim.bim   Same SNP order, positions and alleles.  Chromosome forced to
                  1 and SNP IDs renamed "alt<n>".  Effect: on autosome-only sets
                  this should be numerically inert (a pure relabelling); on
                  sexchr it promotes X/Y/MT markers to autosomes.

<ds>.badfam.fam   altfam minus its last line -> sample count disagrees with the
                  .bed.  Error-path probe (undercount).

<ds>.badbim.bim   corpus .bim minus its last line -> SNP count disagrees with
                  the .bed.  Error-path probe (undercount).

<ds>.bigfam.fam   altfam plus a duplicate of its last row (IID suffixed "_X")
                  -> one MORE sample than the .bed holds.  Error-path probe
                  (overcount; forces KING to read past the end of the .bed).

<ds>.bigbim.bim   corpus .bim plus one synthesised trailing SNP -> one MORE SNP
                  than the .bed holds.  Error-path probe (overcount).

Usage
-----
    python3 make_alt_inputs.py --datadir <dir with corpus files> \
                               --outdir  <dir to write alternates into>
"""

import argparse
import os
import sys

DATASETS = [
    "trio", "nuclear", "threegen", "multifam", "dups", "missing",
    "monomorphic", "sexchr", "unrelated", "admixed", "singleton",
    "pair", "bigish",
]

# .fam columns: FID IID PAT MAT SEX PHENO
FLIP_SEX = {"1": "2", "2": "1", "0": "0"}


def read_rows(path):
    with open(path, "r", newline="") as fh:
        return [ln.rstrip("\n").split() for ln in fh if ln.strip()]


def write_rows(path, rows):
    # KING's corpus files are space-delimited with a trailing newline per row.
    with open(path, "w", newline="") as fh:
        for row in rows:
            fh.write(" ".join(row) + "\n")


def make_altfam(rows):
    out = []
    for i, r in enumerate(rows, start=1):
        fid, iid, _pat, _mat, sex, pheno = r[0], r[1], r[2], r[3], r[4], r[5]
        out.append(["AF%d" % i, iid, "0", "0", FLIP_SEX.get(sex, sex), pheno])
    return out


def make_altbim(rows):
    out = []
    for i, r in enumerate(rows, start=1):
        _chr, _snp, cm, bp, a1, a2 = r[0], r[1], r[2], r[3], r[4], r[5]
        out.append(["1", "alt%d" % i, cm, bp, a1, a2])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    todo = args.only if args.only else DATASETS

    made = 0
    for ds in todo:
        fam_in = os.path.join(args.datadir, ds + ".fam")
        bim_in = os.path.join(args.datadir, ds + ".bim")
        if not (os.path.exists(fam_in) and os.path.exists(bim_in)):
            print("skip %s (corpus files absent)" % ds, file=sys.stderr)
            continue

        fam = read_rows(fam_in)
        bim = read_rows(bim_in)

        altfam = make_altfam(fam)
        altbim = make_altbim(bim)

        write_rows(os.path.join(args.outdir, ds + ".altfam.fam"), altfam)
        write_rows(os.path.join(args.outdir, ds + ".altbim.bim"), altbim)
        write_rows(os.path.join(args.outdir, ds + ".badfam.fam"), altfam[:-1])
        write_rows(os.path.join(args.outdir, ds + ".badbim.bim"), bim[:-1])

        # Overcount variants: one more record than the .bed actually holds.
        extra_fam = list(altfam[-1])
        extra_fam[0] = "AF%d" % (len(altfam) + 1)
        extra_fam[1] = extra_fam[1] + "_X"
        write_rows(os.path.join(args.outdir, ds + ".bigfam.fam"), altfam + [extra_fam])

        last = bim[-1]
        extra_bim = [last[0], "extra_snp", last[2], str(int(last[3]) + 1000),
                     last[4], last[5]]
        write_rows(os.path.join(args.outdir, ds + ".bigbim.bim"), bim + [extra_bim])

        made += 1
        print("%-12s fam=%d bim=%d" % (ds, len(fam), len(bim)))

    print("\n%d datasets -> %s" % (made, args.outdir))


if __name__ == "__main__":
    main()
