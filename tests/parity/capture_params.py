#!/usr/bin/env python3
"""Capture reference-KING golden output for the *parameter* flag group.

Flag group covered here:
    --kinship --prefix custom | --kinship --cpus N | --duplicate --minConc X
    --kinship --sexchr N      | --fam/--bim overrides pointing at alt files
plus --degree (analysis probe) and bare --kinship/--duplicate baselines so the
effect of each parameter can be isolated by diffing.

For every (dataset x flag-combination) a clean directory is created under
    <golden>/params/<dataset>__<flagslug>/
KING is run there with cwd set to that directory (KING writes relative to cwd),
and we save stdout.txt, stderr.txt, exitcode.txt, cmd.txt and every file the run
produced.

cmd.txt holds the exact argv with absolute input paths rewritten to the
placeholders {DATA} and {ALT} so it is machine-replayable.
"""

import argparse
import os
import shutil
import subprocess
import sys
import time

KING = "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king"

DATASETS = [
    "trio", "nuclear", "threegen", "multifam", "dups", "missing",
    "monomorphic", "sexchr", "unrelated", "admixed", "singleton",
    "pair", "bigish",
]

# Combinations applied to EVERY dataset.  "{ALT}/<ds>..." tokens are expanded
# per dataset; DS is replaced by the dataset name.
COMMON = [
    ("kinship",                ["--kinship"]),
    ("kinship_prefix_custom",  ["--kinship", "--prefix", "custom"]),
    ("kinship_cpus1",          ["--kinship", "--cpus", "1"]),
    ("kinship_cpus4",          ["--kinship", "--cpus", "4"]),
    ("kinship_sexchr23",       ["--kinship", "--sexchr", "23"]),
    ("kinship_degree1",        ["--kinship", "--degree", "1"]),
    ("kinship_degree2",        ["--kinship", "--degree", "2"]),
    ("kinship_degree3",        ["--kinship", "--degree", "3"]),
    ("kinship_altfam",         ["--kinship", "--fam", "{ALT}/DS.altfam.fam"]),
    ("kinship_altbim",         ["--kinship", "--bim", "{ALT}/DS.altbim.bim"]),
    ("kinship_altfam_altbim",  ["--kinship", "--fam", "{ALT}/DS.altfam.fam",
                                            "--bim", "{ALT}/DS.altbim.bim"]),
    ("duplicate",              ["--duplicate"]),
    ("duplicate_minConc0.9",   ["--duplicate", "--minConc", "0.9"]),
]

# Dataset-specific edge / error-path probes.
EXTRA = {
    "trio": [
        ("kinship_badfam",        ["--kinship", "--fam", "{ALT}/DS.badfam.fam"]),
        ("kinship_badbim",        ["--kinship", "--bim", "{ALT}/DS.badbim.bim"]),
        ("kinship_bigfam",        ["--kinship", "--fam", "{ALT}/DS.bigfam.fam"]),
        ("kinship_bigbim",        ["--kinship", "--bim", "{ALT}/DS.bigbim.bim"]),
        ("kinship_famnotfound",   ["--kinship", "--fam", "{ALT}/DS.no_such_file.fam"]),
        ("kinship_bimnotfound",   ["--kinship", "--bim", "{ALT}/DS.no_such_file.bim"]),
        ("kinship_cpus0",         ["--kinship", "--cpus", "0"]),
        ("kinship_prefix_dotted", ["--kinship", "--prefix", "cus.tom"]),
        ("kinship_prefix_subdir", ["--kinship", "--prefix", "sub/pre"]),
        ("kinship_prefix_traildot", ["--kinship", "--prefix", "custom."]),
        ("bysample",              ["--bysample"]),
        ("bysample_prefix_custom", ["--bysample", "--prefix", "custom"]),
        ("duplicate_minConc0",    ["--duplicate", "--minConc", "0"]),
        ("duplicate_minConc1",    ["--duplicate", "--minConc", "1"]),
    ],
    "multifam": [
        ("kinship_badfam",        ["--kinship", "--fam", "{ALT}/DS.badfam.fam"]),
        ("kinship_badbim",        ["--kinship", "--bim", "{ALT}/DS.badbim.bim"]),
        ("kinship_bigfam",        ["--kinship", "--fam", "{ALT}/DS.bigfam.fam"]),
        ("kinship_bigbim",        ["--kinship", "--bim", "{ALT}/DS.bigbim.bim"]),
        ("bysample",              ["--bysample"]),
        ("bysample_prefix_custom", ["--bysample", "--prefix", "custom"]),
        ("kinship_degree4",       ["--kinship", "--degree", "4"]),
        ("kinship_cpus2",         ["--kinship", "--cpus", "2"]),
    ],
    "dups": [
        ("duplicate_minConc0",    ["--duplicate", "--minConc", "0"]),
        ("duplicate_minConc0.5",  ["--duplicate", "--minConc", "0.5"]),
        ("duplicate_minConc0.99", ["--duplicate", "--minConc", "0.99"]),
        ("duplicate_minConc1",    ["--duplicate", "--minConc", "1"]),
        ("duplicate_cpus1",       ["--duplicate", "--cpus", "1"]),
        ("duplicate_cpus4",       ["--duplicate", "--cpus", "4"]),
    ],
    # NOTE: KING races when writing kingX.kin0 with >1 thread (see
    # NONDETERMINISTIC.txt).  Every X-chromosome golden must be pinned to
    # --cpus 1; the cpus1_* runs below are the authoritative X references.
    "sexchr": [
        ("kinship_sexchr1",       ["--kinship", "--sexchr", "1"]),
        ("kinship_sexchr25",      ["--kinship", "--sexchr", "25"]),
        ("kinship_sexchr24",      ["--kinship", "--sexchr", "24"]),
        ("kinship_cpus1_sexchr23", ["--kinship", "--cpus", "1", "--sexchr", "23"]),
        ("kinship_cpus1_sexchr24", ["--kinship", "--cpus", "1", "--sexchr", "24"]),
        ("kinship_cpus1_sexchr25", ["--kinship", "--cpus", "1", "--sexchr", "25"]),
        ("kinship_cpus1_sexchr2",  ["--kinship", "--cpus", "1", "--sexchr", "2"]),
        ("kinship_cpus1_prefix_custom", ["--kinship", "--cpus", "1", "--prefix", "custom"]),
        ("kinship_cpus1_degree1", ["--kinship", "--cpus", "1", "--degree", "1"]),
        ("kinship_cpus1_degree2", ["--kinship", "--cpus", "1", "--degree", "2"]),
        ("kinship_cpus1_degree3", ["--kinship", "--cpus", "1", "--degree", "3"]),
        ("kinship_cpus1_altfam",  ["--kinship", "--cpus", "1",
                                   "--fam", "{ALT}/DS.altfam.fam"]),
        ("kinship_cpus1_altbim",  ["--kinship", "--cpus", "1",
                                   "--bim", "{ALT}/DS.altbim.bim"]),
        ("duplicate_cpus1",       ["--duplicate", "--cpus", "1"]),
        ("duplicate_cpus1_minConc0.9", ["--duplicate", "--cpus", "1",
                                        "--minConc", "0.9"]),
    ],
    "bigish": [
        ("kinship_cpus8",         ["--kinship", "--cpus", "8"]),
        ("kinship_cpus16",        ["--kinship", "--cpus", "16"]),
    ],
    "unrelated": [
        ("kinship_cpus2",         ["--kinship", "--cpus", "2"]),
    ],
    # N=1 -> N=2 keeps ceil(N/4)==1, so an overcounted .fam makes KING read the
    # .bed's pad bits as a real sample.  Worth pinning down exactly.
    "singleton": [
        ("kinship_bigfam",        ["--kinship", "--fam", "{ALT}/DS.bigfam.fam"]),
        ("kinship_bigbim",        ["--kinship", "--bim", "{ALT}/DS.bigbim.bim"]),
        ("kinship_badbim",        ["--kinship", "--bim", "{ALT}/DS.badbim.bim"]),
    ],
    "pair": [
        ("kinship_bigfam",        ["--kinship", "--fam", "{ALT}/DS.bigfam.fam"]),
        ("kinship_badfam",        ["--kinship", "--fam", "{ALT}/DS.badfam.fam"]),
    ],
}

BOOKKEEPING = {"cmd.txt", "stdout.txt", "stderr.txt", "exitcode.txt"}


def expand(args, ds, datadir, altdir):
    out = []
    for a in args:
        a = a.replace("DS", ds) if "{ALT}" in a else a
        a = a.replace("{ALT}", altdir).replace("{DATA}", datadir)
        out.append(a)
    return out


def placeholderize(args, datadir, altdir):
    out = []
    for a in args:
        a = a.replace(altdir, "{ALT}").replace(datadir, "{DATA}")
        out.append(a)
    return out


def run_one(ds, slug, extra_args, datadir, altdir, outroot, timeout):
    d = os.path.join(outroot, "%s__%s" % (ds, slug))
    if os.path.isdir(d):
        shutil.rmtree(d)
    os.makedirs(d)

    bed = os.path.join(datadir, ds + ".bed")
    argv = [KING, "-b", bed] + expand(extra_args, ds, datadir, altdir)

    t0 = time.time()
    try:
        p = subprocess.run(argv, cwd=d, capture_output=True, timeout=timeout)
        out, err, code = p.stdout, p.stderr, p.returncode
        timed_out = False
    except subprocess.TimeoutExpired as e:
        out = e.stdout or b""
        err = (e.stderr or b"") + b"\n[HARNESS] TIMEOUT after %ds\n" % timeout
        code = -1
        timed_out = True
    dt = time.time() - t0

    with open(os.path.join(d, "stdout.txt"), "wb") as fh:
        fh.write(out)
    with open(os.path.join(d, "stderr.txt"), "wb") as fh:
        fh.write(err)
    with open(os.path.join(d, "exitcode.txt"), "w") as fh:
        fh.write("%d\n" % code)

    pretty = ["king", "-b", "{DATA}/%s.bed" % ds] + \
        placeholderize(expand(extra_args, ds, datadir, altdir), datadir, altdir)
    with open(os.path.join(d, "cmd.txt"), "w") as fh:
        fh.write(" ".join(pretty) + "\n")

    produced = []
    for root, _dirs, files in os.walk(d):
        for f in files:
            rel = os.path.relpath(os.path.join(root, f), d)
            if rel not in BOOKKEEPING:
                produced.append(rel)
    produced.sort()
    return {"dir": os.path.basename(d), "dataset": ds, "slug": slug,
            "exit": code, "secs": round(dt, 3), "produced": produced,
            "timeout": timed_out}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--datadir", required=True)
    ap.add_argument("--altdir", required=True)
    ap.add_argument("--outroot", required=True)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--timeout", type=int, default=180)
    args = ap.parse_args()

    datadir = os.path.abspath(args.datadir)
    altdir = os.path.abspath(args.altdir)
    outroot = os.path.abspath(args.outroot)
    os.makedirs(outroot, exist_ok=True)

    todo = args.only if args.only else DATASETS
    results = []
    for ds in todo:
        combos = list(COMMON) + EXTRA.get(ds, [])
        for slug, extra in combos:
            r = run_one(ds, slug, extra, datadir, altdir, outroot, args.timeout)
            results.append(r)
            flag = "TIMEOUT" if r["timeout"] else ""
            print("%-46s exit=%-4d %6.2fs  %2d files %s" %
                  (r["dir"], r["exit"], r["secs"], len(r["produced"]), flag))
            sys.stdout.flush()

    import json
    with open(os.path.join(outroot, "runs.json"), "w") as fh:
        json.dump(results, fh, indent=1, sort_keys=True)
    print("\n%d runs -> %s" % (len(results), outroot))


if __name__ == "__main__":
    main()
