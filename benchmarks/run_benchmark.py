#!/usr/bin/env python3
"""
open-king benchmark harness.

Measures wall-clock time, peak resident set size and CPU time for a matrix of
(dataset, analysis) pairs, for one or two KING-compatible binaries.  When two
binaries are given every output file produced by the same invocation is
byte-compared and reported per file.

Python 3 standard library only.  POSIX only (uses os.posix_spawn, os.wait4 and
resource); tested on macOS and Linux.

Measurement method
------------------
Each run is spawned with os.posix_spawn and reaped with os.wait4, so the rusage
returned belongs to that one child and nothing else.  This is more precise than
resource.getrusage(RUSAGE_CHILDREN), whose ru_maxrss is a running maximum over
every child the process has ever reaped and therefore cannot be differenced.

ru_maxrss units differ by platform: BYTES on macOS (darwin), KILOBYTES on Linux.
The detected unit is recorded in the results as "ru_maxrss_unit".

Usage
-----
    # measure open-king alone
    ./run_benchmark.py --binary-a ../target/release/open-king

    # measure open-king against the reference KING 2.3.2 build and diff outputs
    ./run_benchmark.py --binary-a ../target/release/open-king \\
                       --binary-b /path/to/king-2.3.2 \\
                       --label-a open-king --label-b king-2.3.2

    # generate the input filesets and stop
    ./run_benchmark.py --gen-only

Sizing
------
Relationship inference is O(n^2 * m) in samples and markers, so sample count is
the expensive axis.  The default ladder holds markers at 100k and doubles
samples (250 / 500 / 1000), which exposes the quadratic term while keeping the
slowest single run in the tens of seconds.  The whole suite, both binaries and
every repetition, is meant to fit inside roughly 40 minutes so that it stays
rerunnable.  Repetitions default to 3 for the same reason; median with min and
max is reported either way.

The harness never copies either binary into the repository.  Both are referenced
by path and executed in place.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import random
import resource
import shutil
import statistics
import subprocess
import sys
import time

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GOLDEN_DIR = os.path.join(REPO_ROOT, "tests", "parity", "golden")
CORPUS_GEN = os.path.join(REPO_ROOT, "tests", "parity", "generate_corpus.py")

# Analysis name -> the flags handed to the binary.
ANALYSES = {
    "kinship": ["--kinship"],
    "ibs": ["--ibs"],
    "related": ["--related"],
    "ibdseg": ["--ibdseg"],
    "unrelated": ["--unrelated"],
    "duplicate": ["--duplicate"],
    "build": ["--build"],
    "cluster": ["--cluster"],
}

DEFAULT_ANALYSES = list(ANALYSES)

# Filesets shipped by tests/parity/generate_corpus.py.
GOLDEN_DATASETS = [
    "trio", "nuclear", "multifam", "sexchr", "missing", "monomorphic",
    "admixed", "dups", "unrelated", "bigish",
]

# The synthetic size ladder, as name: (n_samples, n_markers).  Overridable with
# --large NAME:SAMPLES:MARKERS.
#
# Sizes are chosen so the whole suite stays rerunnable.  Relationship inference
# is O(n^2 * m), so sample count is the expensive axis and marker count the
# cheap one; the ladder therefore holds markers at 100k and doubles samples
# (200 / 400 / 800), which shows the quadratic term directly while keeping the
# slowest single run in the tens of seconds rather than minutes.  Two measured
# points that were rejected for being too slow to rerun: 2000 x 50k put
# --cluster at 79 s per run, and 1000 x 100k put --unrelated at 76 s.
DEFAULT_LARGE = {
    "synth_s": (200, 100000),
    "synth_m": (400, 100000),
    "synth_l": (800, 100000),
}

# Thread modes for the calibration sweep.  "default" lets each binary pick, so
# the two may disagree; "matched" pins both to the same count so a comparison is
# not confounded by that choice; "single" removes parallelism entirely.
THREAD_MODES = ("default", "matched", "single")

LARGE_SEED = 20260819


# --------------------------------------------------------------------------
# machine context
# --------------------------------------------------------------------------

def _sh(cmd):
    try:
        out = subprocess.run(cmd, capture_output=True, text=True, timeout=20)
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def loadavg():
    """1/5/15 minute load average, or None where unavailable."""
    try:
        return list(os.getloadavg())
    except (OSError, AttributeError):
        return None


def maxrss_unit():
    """ru_maxrss is bytes on macOS, kilobytes on Linux."""
    return "bytes" if sys.platform == "darwin" else "kilobytes"


def maxrss_to_mb(value):
    if maxrss_unit() == "bytes":
        return value / (1024.0 * 1024.0)
    return value / 1024.0


def machine_context():
    ctx = {
        "platform": sys.platform,
        "machine": platform.machine(),
        "python": sys.version.split()[0],
        "ru_maxrss_unit": maxrss_unit(),
        "timestamp_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    if sys.platform == "darwin":
        ctx["cpu_model"] = _sh(["sysctl", "-n", "machdep.cpu.brand_string"])
        ctx["cores_physical"] = _sh(["sysctl", "-n", "hw.physicalcpu"])
        ctx["cores_logical"] = _sh(["sysctl", "-n", "hw.logicalcpu"])
        mem = _sh(["sysctl", "-n", "hw.memsize"])
        ctx["ram_bytes"] = int(mem) if mem and mem.isdigit() else None
        ver = _sh(["sw_vers", "-productVersion"])
        build = _sh(["sw_vers", "-buildVersion"])
        ctx["os"] = "macOS %s (%s)" % (ver, build)
    else:
        ctx["cores_logical"] = str(os.cpu_count() or "")
        ctx["os"] = " ".join(platform.uname()[:3])
        try:
            with open("/proc/cpuinfo") as fh:
                for line in fh:
                    if line.startswith("model name"):
                        ctx["cpu_model"] = line.split(":", 1)[1].strip()
                        break
            with open("/proc/meminfo") as fh:
                for line in fh:
                    if line.startswith("MemTotal"):
                        ctx["ram_bytes"] = int(line.split()[1]) * 1024
                        break
        except OSError:
            pass
    if ctx.get("ram_bytes"):
        ctx["ram_gb"] = round(ctx["ram_bytes"] / (1024.0 ** 3), 1)
    ctx["loadavg_at_start"] = loadavg()
    ctx["git_commit"] = _sh(["git", "-C", REPO_ROOT, "rev-parse", "--short", "HEAD"])
    ctx["git_dirty"] = bool(_sh(["git", "-C", REPO_ROOT, "status", "--porcelain"]))
    return ctx


def binary_banner(path):
    """First non-empty line the binary prints when invoked with no input."""
    try:
        out = subprocess.run([path], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError) as exc:
        return "unavailable (%s)" % exc
    for line in (out.stdout + "\n" + out.stderr).splitlines():
        if line.strip():
            return line.strip()
    return "unknown"


# --------------------------------------------------------------------------
# input filesets
# --------------------------------------------------------------------------

def _load_corpus_module():
    if not os.path.exists(CORPUS_GEN):
        raise SystemExit("corpus generator not found at %s" % CORPUS_GEN)
    spec = importlib.util.spec_from_file_location("generate_corpus", CORPUS_GEN)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def ensure_golden(names):
    """Regenerate any golden fileset whose .bed/.bim/.fam is missing."""
    missing = [n for n in names
               if not all(os.path.exists(os.path.join(GOLDEN_DIR, n + ext))
                          for ext in (".bed", ".bim", ".fam"))]
    if not missing:
        return []
    print("regenerating golden filesets: %s" % ", ".join(missing))
    subprocess.run([sys.executable, CORPUS_GEN, "--outdir", GOLDEN_DIR,
                    "--only", *missing], check=True)
    return missing


def _spread_table():
    """byte -> 2 bytes, with input bit i placed at output bit 2*i."""
    tab = []
    for b in range(256):
        v = 0
        for i in range(8):
            if (b >> i) & 1:
                v |= 1 << (2 * i)
        tab.append(v.to_bytes(2, "little"))
    return tab


SPREAD = _spread_table()
_SPREAD_GET = SPREAD.__getitem__


def _bits_p(rng, width, p, rounds=16):
    """An integer of `width` bits, each independently 1 with probability p.

    Uses the binary expansion of p, so it costs `rounds` calls to getrandbits
    regardless of width.  That is what makes whole-cohort simulation affordable
    in pure Python: one call covers every sample at once.
    """
    if width <= 0 or p <= 0.0:
        return 0
    mask = (1 << width) - 1
    if p >= 1.0:
        return mask
    digits = []
    x = p
    for _ in range(rounds):
        x *= 2.0
        if x >= 1.0:
            digits.append(1)
            x -= 1.0
        else:
            digits.append(0)
    r = 0
    for b in reversed(digits):
        u = rng.getrandbits(width)
        r = (u | r) if b else (u & r)
    return r


class _Block:
    """A group of samples that share a role, held as two haplotype bit vectors.

    Bit f of a vector is family f, so every block of the same width is index
    aligned with every other, and inheritance is a whole-block bitwise select.
    """

    def __init__(self, name, size, offset, sex, kind, pa=None, ma=None):
        self.name = name
        self.size = size
        self.offset = offset
        self.sex = sex
        self.kind = kind
        self.pa = pa
        self.ma = ma
        self.mask = (1 << size) - 1
        self.h1 = 0
        self.h2 = 0
        self.sel_p = 0
        self.sel_m = 0


def generate_large(outdir, name, n_samples, n_markers, seed=LARGE_SEED,
                   miss_rate=0.002, dup_err=0.002, quiet=False):
    """Write a large PLINK 1 fileset with real, undeclared relatedness.

    Every sample gets its own FID and no declared parents, which is the cohort
    shape KING's relationship inference is aimed at.  The genotypes carry a real
    pedigree underneath: parent-offspring pairs, full sibs, half sibs,
    grandparents, avuncular pairs, first cousins and duplicate samples, with IBD
    segments produced by a recombination process over the .bim genetic map.

    Map construction, allele frequency model and .bim writing are reused from
    tests/parity/generate_corpus.py.  Only the genotype inner loop is replaced,
    with a bit-parallel version that simulates the whole cohort per SNP.
    """
    gc = _load_corpus_module()
    rng = random.Random(seed)

    snps = gc.build_map(gc.AUTOSOMES, n_markers, rng)
    n_snps = len(snps)
    snpinfo = gc.model_common(rng, n_snps, snps)

    fam_n = max(4, n_samples // 12)
    dup_n = max(2, fam_n // 2)
    filler = n_samples - 10 * fam_n - dup_n
    if filler < 0:
        fam_n = max(4, (n_samples - 2) // 11)
        dup_n = max(2, fam_n // 2)
        filler = max(0, n_samples - 10 * fam_n - dup_n)
    total = 10 * fam_n + dup_n + filler

    off = 0
    blocks = []

    def add(nm, size, sex, kind, pa=None, ma=None):
        nonlocal off
        blk = _Block(nm, size, off, sex, kind, pa, ma)
        off += size
        blocks.append(blk)
        return blk

    g1 = add("g1", fam_n, 1, "founder")
    g2 = add("g2", fam_n, 2, "founder")
    g3 = add("g3", fam_n, 2, "founder")
    x1 = add("x1", fam_n, 2, "founder")
    x2 = add("x2", fam_n, 1, "founder")
    a = add("a", fam_n, 1, "child", g1, g2)
    b = add("b", fam_n, 2, "child", g1, g2)
    h = add("h", fam_n, 1, "child", g1, g3)
    c = add("c", fam_n, 0, "child", a, x1)
    e = add("e", fam_n, 0, "child", x2, b)
    dup = add("dup", dup_n, 1, "dup", g1)
    unrel = add("u", filler, 0, "founder") if filler > 0 else None

    founders = [blk for blk in blocks if blk.kind == "founder"]
    children = [blk for blk in blocks if blk.kind == "child"]

    base = os.path.join(outdir, name)

    # .fam: distinct FID per sample, no declared parents.
    with open(base + ".fam", "w") as fh:
        idx = 0
        for blk in blocks:
            for k in range(blk.size):
                sex = blk.sex if blk.sex else (1 + (k % 2))
                iid = "S%06d" % idx
                fh.write("%s %s 0 0 %d -9\n" % (iid, iid, sex))
                idx += 1

    row_bytes = (total + 3) // 4
    vec_bytes = (total + 7) // 8
    row_mask = (1 << total) - 1

    founder_w = 2 * sum(blk.size for blk in founders)
    sel_w = 2 * fam_n * len(children)

    alleles = []
    out = []
    flush_every = 512
    started = time.perf_counter()

    with open(base + ".bed", "wb") as fh:
        fh.write(gc.BED_MAGIC)
        prev_chrom = None
        prev_cm = 0.0

        for j in range(n_snps):
            chrom, _bp, cm = snps[j]
            if chrom != prev_chrom:
                recomb = None
            else:
                recomb = 0.5 * (1.0 - math.exp(-2.0 * (cm - prev_cm) / 100.0))
            prev_chrom, prev_cm = chrom, cm
            q = snpinfo[j]["freqs"][0]

            # founder haplotypes: one draw covers every founder in the cohort
            fbits = _bits_p(rng, founder_w, q)
            pos = 0
            for blk in founders:
                blk.h1 = (fbits >> pos) & blk.mask
                pos += blk.size
                blk.h2 = (fbits >> pos) & blk.mask
                pos += blk.size

            # meiosis: which parental haplotype each child carries
            if recomb is None:
                sbits = rng.getrandbits(sel_w)
                pos = 0
                for blk in children:
                    blk.sel_p = (sbits >> pos) & blk.mask
                    pos += blk.size
                    blk.sel_m = (sbits >> pos) & blk.mask
                    pos += blk.size
            elif recomb > 0.0:
                sbits = _bits_p(rng, sel_w, recomb)
                pos = 0
                for blk in children:
                    blk.sel_p ^= (sbits >> pos) & blk.mask
                    pos += blk.size
                    blk.sel_m ^= (sbits >> pos) & blk.mask
                    pos += blk.size

            for blk in children:
                pa, ma = blk.pa, blk.ma
                sp, sm = blk.sel_p, blk.sel_m
                blk.h1 = (pa.h2 & sp) | (pa.h1 & (sp ^ blk.mask))
                blk.h2 = (ma.h2 & sm) | (ma.h1 & (sm ^ blk.mask))

            # duplicate samples: same haplotypes as the first dup_n of g1,
            # plus a low rate of independent genotyping error
            dup.h1 = g1.h1 & dup.mask
            dup.h2 = g1.h2 & dup.mask
            if dup_err > 0.0:
                ebits = _bits_p(rng, 2 * dup.size, dup_err)
                dup.h1 ^= ebits & dup.mask
                dup.h2 ^= (ebits >> dup.size) & dup.mask

            row_h1 = 0
            row_h2 = 0
            for blk in blocks:
                row_h1 |= blk.h1 << blk.offset
                row_h2 |= blk.h2 << blk.offset

            miss = _bits_p(rng, total, miss_rate) if miss_rate > 0.0 else 0
            keep = miss ^ row_mask

            called = total - miss.bit_count()
            a1n = (row_h1 & keep).bit_count() + (row_h2 & keep).bit_count()
            if called and a1n > called:
                # PLINK --make-bed keeps A1 as the observed minor allele.
                row_h1 ^= row_mask
                row_h2 ^= row_mask
                al = (snpinfo[j]["a2"], snpinfo[j]["a1"])
            else:
                al = (snpinfo[j]["a1"], snpinfo[j]["a2"])
            alleles.append(al)

            # 2-bit codes: bit0 set for hom-A2 or missing, bit1 set unless hom-A1
            bit0 = (((row_h1 | row_h2) ^ row_mask) & keep) | miss
            bit1 = ((row_h1 & row_h2) ^ row_mask) & keep

            s0 = int.from_bytes(
                b"".join(map(_SPREAD_GET, bit0.to_bytes(vec_bytes, "little"))),
                "little")
            s1 = int.from_bytes(
                b"".join(map(_SPREAD_GET, bit1.to_bytes(vec_bytes, "little"))),
                "little")
            out.append((s0 | (s1 << 1)).to_bytes(row_bytes, "little"))

            if len(out) >= flush_every:
                fh.write(b"".join(out))
                out.clear()
                if not quiet and j % 10000 < flush_every:
                    print("    %s: %d/%d SNPs (%.0fs)"
                          % (name, j + 1, n_snps, time.perf_counter() - started))
        if out:
            fh.write(b"".join(out))

    gc.write_bim(base + ".bim", snps, alleles)
    if not quiet:
        print("  %s: %d samples, %d markers, %.1f MB bed, %.0fs"
              % (name, total, n_snps, os.path.getsize(base + ".bed") / 1e6,
                 time.perf_counter() - started))
    return {"name": name, "n_samples": total, "n_markers": n_snps,
            "seed": seed, "fam_blocks": fam_n, "duplicate_pairs": dup_n,
            "unrelated_filler": filler}


def ensure_large(work_dir, large_spec, seed=LARGE_SEED, force=False):
    os.makedirs(work_dir, exist_ok=True)
    meta = {}
    for name, (n_samples, n_markers) in large_spec.items():
        base = os.path.join(work_dir, name)
        stamp = base + ".meta.json"
        if not force and os.path.exists(stamp):
            with open(stamp) as fh:
                cached = json.load(fh)
            if (cached.get("n_samples") == n_samples
                    and cached.get("n_markers") == n_markers
                    and cached.get("seed") == seed
                    and all(os.path.exists(base + ext)
                            for ext in (".bed", ".bim", ".fam"))):
                meta[name] = cached
                continue
        print("generating %s (%d samples x %d markers)"
              % (name, n_samples, n_markers))
        info = generate_large(work_dir, name, n_samples, n_markers, seed=seed)
        with open(stamp, "w") as fh:
            json.dump(info, fh, indent=2, sort_keys=True)
            fh.write("\n")
        meta[name] = info
    return meta


def _is_golden(name):
    return name in GOLDEN_DATASETS


def dataset_facts(name, path_base):
    with open(path_base + ".fam") as fh:
        n_samples = sum(1 for _ in fh)
    with open(path_base + ".bim") as fh:
        n_markers = sum(1 for _ in fh)
    return {"name": name,
            "n_samples": n_samples,
            "n_markers": n_markers,
            "bed_bytes": os.path.getsize(path_base + ".bed"),
            "path": path_base + ".bed"}


# --------------------------------------------------------------------------
# measurement
# --------------------------------------------------------------------------

def run_once(binary, argv_tail, out_dir, prefix="king_"):
    """Spawn one run and reap it, returning wall time, CPU time and peak RSS.

    os.wait4 gives the rusage of this one child, so peak RSS is that process's
    own maximum and nothing else's.
    """
    os.makedirs(out_dir, exist_ok=True)
    argv = [binary] + argv_tail + ["--prefix", os.path.join(out_dir, prefix)]
    stdout_path = os.path.join(out_dir, "_stdout.txt")
    stderr_path = os.path.join(out_dir, "_stderr.txt")
    wflags = os.O_WRONLY | os.O_CREAT | os.O_TRUNC
    actions = [
        (os.POSIX_SPAWN_OPEN, 1, stdout_path, wflags, 0o644),
        (os.POSIX_SPAWN_OPEN, 2, stderr_path, wflags, 0o644),
    ]
    t0 = time.perf_counter()
    pid = os.posix_spawn(binary, argv, os.environ, file_actions=actions)
    _pid, status, ru = os.wait4(pid, 0)
    wall = time.perf_counter() - t0
    cpu = ru.ru_utime + ru.ru_stime
    return {
        "wall_s": wall,
        "user_s": ru.ru_utime,
        "sys_s": ru.ru_stime,
        "cpu_s": cpu,
        "cpu_per_wall": (cpu / wall) if wall > 0 else 0.0,
        "peak_rss_raw": ru.ru_maxrss,
        "peak_rss_mb": maxrss_to_mb(ru.ru_maxrss),
        "exit_code": os.waitstatus_to_exitcode(status),
        "argv": argv,
    }


def summarise(runs):
    walls = [r["wall_s"] for r in runs]
    cpus = [r["cpu_s"] for r in runs]
    rss = [r["peak_rss_mb"] for r in runs]
    return {
        "reps": len(runs),
        "wall_median_s": statistics.median(walls),
        "wall_min_s": min(walls),
        "wall_max_s": max(walls),
        "cpu_median_s": statistics.median(cpus),
        "user_median_s": statistics.median(r["user_s"] for r in runs),
        "sys_median_s": statistics.median(r["sys_s"] for r in runs),
        "cpu_per_wall_median": statistics.median(r["cpu_per_wall"] for r in runs),
        "peak_rss_mb_median": statistics.median(rss),
        "peak_rss_mb_max": max(rss),
        "exit_codes": sorted({r["exit_code"] for r in runs}),
    }


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def output_files(out_dir, prefix="king_"):
    """Files the run produced, keyed by name with the prefix stripped."""
    found = {}
    if not os.path.isdir(out_dir):
        return found
    for entry in sorted(os.listdir(out_dir)):
        if entry.startswith("_"):
            continue
        path = os.path.join(out_dir, entry)
        if not os.path.isfile(path):
            continue
        key = entry[len(prefix):] if entry.startswith(prefix) else entry
        found[key] = path
    return found


def compare_outputs(dir_a, dir_b, prefix="king_"):
    """Byte-compare every output file the two binaries produced."""
    fa = output_files(dir_a, prefix)
    fb = output_files(dir_b, prefix)
    files = []
    identical = differing = 0
    for key in sorted(set(fa) | set(fb)):
        pa, pb = fa.get(key), fb.get(key)
        if pa is None:
            files.append({"file": key, "status": "only_in_b",
                          "size_b": os.path.getsize(pb)})
            differing += 1
            continue
        if pb is None:
            files.append({"file": key, "status": "only_in_a",
                          "size_a": os.path.getsize(pa)})
            differing += 1
            continue
        ha, hb = sha256_file(pa), sha256_file(pb)
        same = ha == hb
        files.append({"file": key,
                      "status": "identical" if same else "different",
                      "size_a": os.path.getsize(pa),
                      "size_b": os.path.getsize(pb),
                      "sha256_a": ha, "sha256_b": hb})
        if same:
            identical += 1
        else:
            differing += 1
    return {"n_files": len(files), "identical": identical,
            "different": differing, "files": files}


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def mode_flags(mode, matched_cpus):
    if mode == "matched":
        return ["--cpus", str(matched_cpus)]
    if mode == "single":
        return ["--cpus", "1"]
    return []


def bench_cell(binary, label, ds, analysis, work_dir, reps, warmup,
               mode="default", matched_cpus=4):
    tail = (["-b", ds["path"]] + ANALYSES[analysis]
            + mode_flags(mode, matched_cpus))
    out_dir = os.path.join(work_dir, "out", label, mode, ds["name"], analysis)
    shutil.rmtree(out_dir, ignore_errors=True)
    load_before = loadavg()
    for _ in range(warmup):
        run_once(binary, tail, out_dir)
    runs = []
    for _ in range(reps):
        shutil.rmtree(out_dir, ignore_errors=True)
        runs.append(run_once(binary, tail, out_dir))
    summary = summarise(runs)
    summary["output_files"] = sorted(output_files(out_dir))
    summary["out_dir"] = out_dir
    summary["mode"] = mode
    summary["cpus_flag"] = mode_flags(mode, matched_cpus)
    summary["loadavg_before"] = load_before
    summary["loadavg_after"] = loadavg()
    return summary


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Benchmark one or two KING-compatible binaries.")
    ap.add_argument("--binary-a", default=os.path.join(REPO_ROOT, "target",
                                                       "release", "open-king"),
                    help="first binary (default: the repo release build)")
    ap.add_argument("--binary-b", default=None,
                    help="optional second binary, for timing and output diffs")
    ap.add_argument("--label-a", default="open-king")
    ap.add_argument("--label-b", default="reference")
    ap.add_argument("--work-dir", required=False,
                    default=os.environ.get("OPEN_KING_BENCH_WORKDIR",
                                           "/tmp/open-king-bench"),
                    help="scratch directory for generated inputs and outputs; "
                         "keep it outside the repository")
    ap.add_argument("--results-dir",
                    default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "results"))
    ap.add_argument("--reps", type=int, default=3,
                    help="timed repetitions per cell (default 3, which keeps "
                         "the whole suite inside a ~40 minute budget)")
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--matched-cpus", type=int, default=4,
                    help="thread count both binaries are pinned to in the "
                         "'matched' calibration mode")
    ap.add_argument("--calibrate-dataset", default="synth_m")
    ap.add_argument("--calibrate-analyses", nargs="+",
                    default=["kinship", "ibdseg", "cluster"])
    ap.add_argument("--calibrate-modes", nargs="+", default=list(THREAD_MODES),
                    choices=list(THREAD_MODES))
    ap.add_argument("--no-calibrate", action="store_true",
                    help="skip the thread-count calibration sweep")
    ap.add_argument("--datasets", nargs="+", default=None,
                    help="dataset names (default: the golden corpus plus the "
                         "generated large filesets)")
    ap.add_argument("--analyses", nargs="+", default=DEFAULT_ANALYSES,
                    choices=sorted(ANALYSES))
    ap.add_argument("--large", nargs="+", default=None,
                    metavar="NAME:SAMPLES:MARKERS",
                    help="override the generated large filesets")
    ap.add_argument("--seed", type=int, default=LARGE_SEED)
    ap.add_argument("--gen-only", action="store_true",
                    help="prepare input filesets and exit")
    ap.add_argument("--regen-large", action="store_true")
    ap.add_argument("--skip-large", action="store_true")
    args = ap.parse_args(argv)

    if os.name != "posix":
        raise SystemExit("this harness needs a POSIX platform")

    large_spec = dict(DEFAULT_LARGE)
    if args.large:
        large_spec = {}
        for item in args.large:
            name, n_s, n_m = item.split(":")
            large_spec[name] = (int(n_s), int(n_m))
    if args.skip_large:
        large_spec = {}

    work_dir = os.path.abspath(args.work_dir)
    os.makedirs(work_dir, exist_ok=True)
    if os.path.commonpath([work_dir, REPO_ROOT]) == REPO_ROOT:
        print("warning: work dir is inside the repository (%s)" % work_dir)

    ensure_golden(GOLDEN_DATASETS)
    large_meta = ensure_large(work_dir, large_spec, seed=args.seed,
                              force=args.regen_large)

    names = args.datasets
    if names is None:
        names = GOLDEN_DATASETS + list(large_spec)

    datasets = []
    for name in names:
        root = GOLDEN_DIR if _is_golden(name) else work_dir
        base = os.path.join(root, name)
        if not os.path.exists(base + ".bed"):
            raise SystemExit("no fileset for %s at %s.bed" % (name, base))
        datasets.append(dataset_facts(name, base))

    if args.gen_only:
        for ds in datasets:
            print("%-10s %6d samples %8d markers %10d bed bytes"
                  % (ds["name"], ds["n_samples"], ds["n_markers"],
                     ds["bed_bytes"]))
        return 0

    binaries = [(args.label_a, os.path.abspath(args.binary_a))]
    if args.binary_b:
        binaries.append((args.label_b, os.path.abspath(args.binary_b)))
    for label, path in binaries:
        if not os.path.isfile(path) or not os.access(path, os.X_OK):
            raise SystemExit("%s is not an executable file: %s" % (label, path))

    results = {
        "machine": machine_context(),
        "binaries": [{"label": label, "path": path,
                      "banner": binary_banner(path)}
                     for label, path in binaries],
        "config": {"reps": args.reps, "warmup": args.warmup,
                   "reference_label": args.label_b,
                   "analyses": args.analyses, "work_dir": work_dir,
                   "large_spec": {k: list(v) for k, v in large_spec.items()},
                   "large_meta": large_meta, "seed": args.seed,
                   "matched_cpus": args.matched_cpus,
                   "calibrate_dataset": args.calibrate_dataset,
                   "calibrate_analyses": args.calibrate_analyses,
                   "calibrate_modes": ([] if args.no_calibrate
                                       else args.calibrate_modes)},
        "datasets": datasets,
        "cells": [],
    }

    total_cells = len(datasets) * len(args.analyses)
    done = 0
    t_start = time.perf_counter()
    for ds in datasets:
        for analysis in args.analyses:
            done += 1
            cell = {"dataset": ds["name"], "analysis": analysis, "by_binary": {}}
            for label, path in binaries:
                cell["by_binary"][label] = bench_cell(
                    path, label, ds, analysis, work_dir, args.reps, args.warmup)
            if len(binaries) == 2:
                la, lb = binaries[0][0], binaries[1][0]
                cell["comparison"] = compare_outputs(
                    cell["by_binary"][la]["out_dir"],
                    cell["by_binary"][lb]["out_dir"])
            results["cells"].append(cell)
            first = cell["by_binary"][binaries[0][0]]
            print("[%3d/%3d] %-10s %-10s median %8.3fs  rss %7.1f MB  "
                  "cpu/wall %5.2f  exit %s  (%.0fs elapsed)"
                  % (done, total_cells, ds["name"], analysis,
                     first["wall_median_s"], first["peak_rss_mb_median"],
                     first["cpu_per_wall_median"], first["exit_codes"],
                     time.perf_counter() - t_start))

    # Thread-count calibration.  Run on one dataset and a few analyses only:
    # the point is the parallelism factor, and a single-threaded sweep of the
    # whole matrix would dominate the budget.
    results["calibration"] = []
    if not args.no_calibrate and args.calibrate_modes:
        cal_ds = next((d for d in datasets
                       if d["name"] == args.calibrate_dataset), None)
        if cal_ds is None:
            print("calibration dataset %s not in the run; skipping"
                  % args.calibrate_dataset)
        else:
            cal_analyses = [a for a in args.calibrate_analyses if a in ANALYSES]
            for analysis in cal_analyses:
                for mode in args.calibrate_modes:
                    entry = {"dataset": cal_ds["name"], "analysis": analysis,
                             "mode": mode, "by_binary": {}}
                    for label, path in binaries:
                        entry["by_binary"][label] = bench_cell(
                            path, label + "-cal", cal_ds, analysis, work_dir,
                            args.reps, args.warmup, mode=mode,
                            matched_cpus=args.matched_cpus)
                    results["calibration"].append(entry)
                    f = entry["by_binary"][binaries[0][0]]
                    print("[cal ] %-10s %-10s %-8s median %8.3fs  "
                          "cpu/wall %5.2f  (%.0fs elapsed)"
                          % (cal_ds["name"], analysis, mode,
                             f["wall_median_s"], f["cpu_per_wall_median"],
                             time.perf_counter() - t_start))

    results["config"]["total_runtime_s"] = time.perf_counter() - t_start
    results["machine"]["loadavg_at_end"] = loadavg()

    os.makedirs(args.results_dir, exist_ok=True)
    json_path = os.path.join(args.results_dir, "results.json")
    with open(json_path, "w") as fh:
        json.dump(results, fh, indent=2, sort_keys=False)
        fh.write("\n")
    md_path = os.path.join(args.results_dir, "results.md")
    with open(md_path, "w") as fh:
        fh.write(render_markdown(results))
    print("wrote %s and %s" % (json_path, md_path))
    return 0


# --------------------------------------------------------------------------
# reporting
# --------------------------------------------------------------------------

def _fmt(value, spec="%.3f"):
    return "" if value is None else spec % value


def render_markdown(results):
    labels = [b["label"] for b in results["binaries"]]
    two = len(labels) == 2
    m = results["machine"]
    cfg = results["config"]
    out = []
    w = out.append

    w("# open-king benchmark results\n")
    w("Generated %s by `benchmarks/run_benchmark.py`.\n"
      % m.get("timestamp_utc", ""))

    w("\n## Machine\n")
    w("| item | value |")
    w("| --- | --- |")
    w("| CPU | %s |" % (m.get("cpu_model") or "unknown"))
    w("| Cores | %s physical / %s logical |"
      % (m.get("cores_physical", "?"), m.get("cores_logical", "?")))
    w("| RAM | %s GB |" % m.get("ram_gb", "?"))
    w("| OS | %s |" % m.get("os", "?"))
    w("| Arch | %s |" % m.get("machine", "?"))
    w("| Python | %s |" % m.get("python", "?"))
    w("| open-king commit | %s%s |"
      % (m.get("git_commit") or "?", " (dirty tree)" if m.get("git_dirty") else ""))
    w("| `ru_maxrss` unit detected | %s |" % m.get("ru_maxrss_unit"))
    for key, lbl in (("loadavg_at_start", "load average at start"),
                     ("loadavg_at_end", "load average at end")):
        la = m.get(key)
        if la:
            w("| %s | %.2f / %.2f / %.2f (1/5/15 min) |"
              % (lbl, la[0], la[1], la[2]))
    for b in results["binaries"]:
        w("| binary `%s` | `%s` |" % (b["label"], b["path"]))
        w("| banner `%s` | %s |" % (b["label"], b["banner"]))
    w("| repetitions | %d timed, %d warmup |" % (cfg["reps"], cfg["warmup"]))
    if cfg.get("total_runtime_s"):
        w("| harness runtime | %.0f s |" % cfg["total_runtime_s"])

    loads = [c["by_binary"][labels[0]].get("loadavg_before")
             for c in results["cells"]]
    loads = [l[0] for l in loads if l]
    try:
        cores = int(m.get("cores_logical") or 0)
    except (TypeError, ValueError):
        cores = 0
    if loads and cores:
        lo, hi = min(loads), max(loads)
        w("\n### Measurement conditions\n")
        w("Load average during the run ranged %.1f to %.1f on %d logical cores."
          % (lo, hi, cores))
        if hi > cores * 0.5:
            w("\n**The machine was not idle.** Other work on this host was "
              "competing for cores throughout, so wall-clock times are inflated "
              "and `cpu/wall` is pushed down by descheduling rather than by any "
              "property of the binary. Treat the CPU-seconds column as the "
              "load-resistant measure: it counts work done, not time elapsed. "
              "Wall times are comparable within this run, since every cell "
              "faced the same contention, but they are an upper bound on what "
              "an idle host would show.\n")
        else:
            w(" The host was substantially idle, so wall times are "
              "representative.\n")

    w("\n## Datasets\n")
    w("| dataset | samples | markers | .bed size |")
    w("| --- | ---: | ---: | ---: |")
    for ds in results["datasets"]:
        size = ds["bed_bytes"]
        human = ("%.1f MB" % (size / 1e6)) if size >= 1e6 else ("%.0f KB" % (size / 1e3))
        w("| %s | %d | %d | %s |"
          % (ds["name"], ds["n_samples"], ds["n_markers"], human))

    if cfg.get("large_meta"):
        w("\nLarge filesets are generated by the harness (seed %d) into the work "
          "directory, not committed. Every sample carries its own FID and no "
          "declared parents; the relatedness is real but undeclared, which is "
          "the cohort shape KING's inference targets.\n" % cfg["seed"])
        w("| dataset | family blocks | duplicate pairs | unrelated filler |")
        w("| --- | ---: | ---: | ---: |")
        for name, info in cfg["large_meta"].items():
            w("| %s | %d | %d | %d |"
              % (name, info["fam_blocks"], info["duplicate_pairs"],
                 info["unrelated_filler"]))

    a = labels[0]
    b = labels[1] if two else results["config"].get("reference_label", "reference")

    w("\n## Timings\n")
    w("Wall time is the median of %d timed repetitions, with min and max "
      "alongside. Peak RSS is the maximum resident set size of the child "
      "process. `cpu/wall` is (user + sys) CPU seconds divided by wall seconds, "
      "so it reads as the average number of cores kept busy.\n" % cfg["reps"])

    # The second binary's columns are always present. With one binary they are
    # blank, so adding the reference later does not change the table shape.
    head = ["dataset", "analysis",
            "%s wall med (s)" % a, "%s wall min/max (s)" % a,
            "%s peak RSS (MB)" % a, "%s CPU (s)" % a, "%s cpu/wall" % a,
            "%s wall med (s)" % b, "%s peak RSS (MB)" % b, "%s cpu/wall" % b,
            "%s/%s wall" % (b, a), "files", "same", "differ"]
    align = ["---", "---"] + ["---:"] * (len(head) - 2)
    w("| " + " | ".join(head) + " |")
    w("| " + " | ".join(align) + " |")

    for cell in results["cells"]:
        sa = cell["by_binary"][a]
        row = [cell["dataset"], cell["analysis"],
               "%.3f" % sa["wall_median_s"],
               "%.3f / %.3f" % (sa["wall_min_s"], sa["wall_max_s"]),
               "%.1f" % sa["peak_rss_mb_median"],
               "%.3f" % sa["cpu_median_s"],
               "%.2f" % sa["cpu_per_wall_median"]]
        sb = cell["by_binary"].get(b) if two else None
        cmpn = cell.get("comparison", {})
        if sb:
            speed = (sb["wall_median_s"] / sa["wall_median_s"]
                     if sa["wall_median_s"] > 0 else None)
            row += ["%.3f" % sb["wall_median_s"],
                    "%.1f" % sb["peak_rss_mb_median"],
                    "%.2f" % sb["cpu_per_wall_median"],
                    _fmt(speed, "%.2fx")]
        else:
            row += ["", "", "", ""]
        row += [str(cmpn.get("n_files", "")),
                str(cmpn.get("identical", "")),
                str(cmpn.get("different", ""))]
        w("| " + " | ".join(row) + " |")

    nonzero = [c for c in results["cells"]
               if c["by_binary"][a]["exit_codes"] != [0]]
    if nonzero:
        w("\n### Non-zero exits\n")
        w("| dataset | analysis | exit codes |")
        w("| --- | --- | --- |")
        for c in nonzero:
            w("| %s | %s | %s |"
              % (c["dataset"], c["analysis"],
                 c["by_binary"][a]["exit_codes"]))

    cal = results.get("calibration") or []
    if cal:
        w("\n## Thread-count calibration\n")
        w("Three thread settings on one dataset. `default` lets the binary "
          "choose, `matched` pins it to `--cpus %d` so two binaries cannot "
          "differ on that choice, and `single` is `--cpus 1`. The single-thread "
          "column is the honest serial cost; dividing it by the default wall "
          "time gives the parallel speedup actually realised.\n"
          % cfg.get("matched_cpus", 4))
        head = ["dataset", "analysis", "mode", "%s wall med (s)" % a,
                "%s cpu/wall" % a, "%s peak RSS (MB)" % a,
                "%s wall med (s)" % b, "%s cpu/wall" % b]
        w("| " + " | ".join(head) + " |")
        w("| " + " | ".join(["---", "---", "---"] + ["---:"] * 5) + " |")
        for entry in cal:
            sa = entry["by_binary"].get(a)
            sb = entry["by_binary"].get(b) if two else None
            row = [entry["dataset"], entry["analysis"], entry["mode"],
                   "%.3f" % sa["wall_median_s"],
                   "%.2f" % sa["cpu_per_wall_median"],
                   "%.1f" % sa["peak_rss_mb_median"]]
            row += (["%.3f" % sb["wall_median_s"],
                     "%.2f" % sb["cpu_per_wall_median"]] if sb else ["", ""])
            w("| " + " | ".join(row) + " |")

        by = {}
        for entry in cal:
            by.setdefault(entry["analysis"], {})[entry["mode"]] = \
                entry["by_binary"][a]["wall_median_s"]
        speedups = [(k, v["single"] / v["default"])
                    for k, v in by.items()
                    if v.get("single") and v.get("default")]
        if speedups:
            w("\nSerial-to-default speedup, %s: %s.\n"
              % (a, ", ".join("%s %.1fx" % (k, r) for k, r in speedups)))

    w("\n## Output files produced\n")
    w("| dataset | analysis | files |")
    w("| --- | --- | --- |")
    for cell in results["cells"]:
        files = cell["by_binary"][a]["output_files"]
        w("| %s | %s | %s |" % (cell["dataset"], cell["analysis"],
                                ", ".join(files) if files else "(none)"))

    w("\n## Output comparison\n")
    if not two:
        w("Only one binary was measured, so there is nothing to compare. Rerun "
          "with `--binary-b /path/to/king-2.3.2` to fill the comparison columns "
          "in the timings table and the per-file table below.\n")
        w("| dataset | analysis | file | status | sha256 %s | sha256 reference |"
          % a)
        w("| --- | --- | --- | --- | --- | --- |")
    else:
        w("Only the analysis output files are compared, byte for byte by "
          "sha256. The harness captures each run's stdout and stderr to "
          "separate files that are excluded from the comparison, because a KING "
          "banner carries a wall-clock timestamp and would differ on every run "
          "even for identical results.\n")
        w("| dataset | analysis | file | status | sha256 %s | sha256 %s |" % (a, b))
        w("| --- | --- | --- | --- | --- | --- |")
        for cell in results["cells"]:
            for f in cell.get("comparison", {}).get("files", []):
                w("| %s | %s | %s | %s | %s | %s |"
                  % (cell["dataset"], cell["analysis"], f["file"], f["status"],
                     (f.get("sha256_a") or "")[:12],
                     (f.get("sha256_b") or "")[:12]))

    return "\n".join(out) + "\n"


if __name__ == "__main__":
    sys.exit(main())
