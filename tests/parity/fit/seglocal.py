"""Localise one pair's `.seg` IBD2 total to **one usable segment** of the real corpus.

`seglen_probe.py` recovers the multiset of called segment *lengths* per pair; inverting a
length back to its markers is ambiguous (a random target already matches about one
interval per dataset, `resid_shape.py`), so a disagreement cannot be pinned to a place.
This pins it.

The trick costs one reference invocation. Rewrite the fileset so that, **for the two
samples of the pair only**, every retained marker outside the segment of interest carries
the `segcanvas` carrier pattern — 34 het-vs-A2A2 mismatches, 4 A1A1/A1A1 and 26 A2A2/A2A2
per 64-marker word. Those words are IBD1-clean and IBD2-unusable, so:

* the pair keeps long IBD1 segments everywhere and is still reported in `.seg`;
* **`IBD2Seg` is then exactly the pair's IBD2 calls inside the kept segment**, over the
  unchanged denominator `D`;
* nothing else in the cohort moves — the `.bim` is byte-identical, so `kingallsegs.txt`
  and `D` are too (asserted on every run).

The printed `%.4lf` resolves `D / 10000` — 40 kb on `dups`, 69 kb on `multifam`, about one
marker gap — so a per-segment total identifies a call set almost uniquely, and
`--seglength` on the same fileset splits it into individual calls.

    python3 seglocal.py dups MZ_1 MZ_2                # every segment, ours vs reference
    python3 seglocal.py dups MZ_1 MZ_2 --seg 8        # one segment, with a length ladder
    python3 seglocal.py --check                       # the instrument's own controls
"""

import hashlib
import json
import os
import subprocess
import sys
import tempfile

import numpy as np

import engine as E
import kingdata as kd
import resid19 as R

KING = os.environ.get(
    "KING", "/Users/wgu/Desktop/GeneQuire Project/GeneQuire/software/king/king")
CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "seglocal_measured.json")
WORD = 64

# A1 dosage -> PLINK1 two-bit code (0 = A1A1 hom-major slot in king-io's reading)
CODE = {0: 0b11, 1: 0b10, 2: 0b00, 3: 0b01}

#: The carrier word: IBD1-clean (no opposite homozygote), IBD2-unusable (34 het-vs-hom
#: mismatches) and `inf1`-rich enough that IBD1 still calls it. Bit `r` of every word.
def carrier(r):
    if r < 34:
        return (1, 0)        # het vs A2A2 — a mismatch, uninformative
    if r < 38:
        return (2, 2)        # A1A1/A1A1 — informative for IBD1 and IBD2
    return (0, 0)            # A2A2/A2A2 — clean, worthless


_MEAS = None


def measured():
    global _MEAS
    if _MEAS is None:
        _MEAS = json.load(open(CACHE)) if os.path.exists(CACHE) else {}
    return _MEAS


def save():
    with open(CACHE, "w") as fh:
        json.dump(measured(), fh, indent=0, sort_keys=True)


def write_fileset(ds, i, j, keep_segs, out, carrier_segs=()):
    """`out.{bed,bim,fam}`: `ds` with the pair neutralised outside `keep_segs`.

    Markers outside `keep_segs` become A2A2/A2A2 for the pair — no opposite homozygote,
    no mismatch, nothing informative — except inside `carrier_segs`, which get the
    mismatch-heavy carrier so the pair still owns a long IBD1 segment.  Keeping the
    carrier out of the words that flank the kept segment is what makes the measurement
    a measurement of the kept segment alone.
    """
    src = os.path.join(kd.DATA, ds.name)
    for ext in (".bim", ".fam"):
        with open(out + ext, "wb") as fh:
            fh.write(open(src + ext, "rb").read())
    n = len(ds.fam)
    bpv = (n + 3) // 4
    raw = np.fromfile(src + ".bed", dtype=np.uint8)
    body = raw[3:].reshape(-1, bpv).copy()
    nvar = body.shape[0]
    codes = np.empty((nvar, bpv * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (body >> (2 * k)) & 3
    # retained-marker index -> original bim row
    orig = np.flatnonzero(ds.keep)
    live = np.zeros(len(orig), dtype=bool)
    for lo, hi in keep_segs:
        live[lo:hi + 1] = True
    carry = np.zeros(len(orig), dtype=bool)
    for lo, hi in carrier_segs:
        carry[lo:hi + 1] = True
    carry &= ~live
    pat_i = np.array([CODE[carrier(r % WORD)[0]] for r in range(WORD)], dtype=np.uint8)
    pat_j = np.array([CODE[carrier(r % WORD)[1]] for r in range(WORD)], dtype=np.uint8)
    dead = ~live & ~carry
    codes[orig[dead], i] = CODE[0]
    codes[orig[dead], j] = CODE[0]
    idx = np.flatnonzero(carry)
    codes[orig[idx], i] = pat_i[idx % WORD]
    codes[orig[idx], j] = pat_j[idx % WORD]
    packed = np.zeros((nvar, bpv), dtype=np.uint8)
    for k in range(4):
        packed |= codes[:, k::4][:, :bpv] << (2 * k)
    with open(out + ".bed", "wb") as fh:
        fh.write(bytes([0x6C, 0x1B, 0x01]))
        packed.tofile(fh)


def pick_carrier(ds, k):
    """A usable segment to carry the pair's long IBD1, never adjacent to segment `k`."""
    best = None
    for t, s in enumerate(ds.segs):
        if abs(t - k) < 2:
            continue
        ln = int(ds.pos[s[2]] - ds.pos[s[1]])
        if best is None or ln > best[0]:
            best = (ln, (s[1], s[2]))
    return [best[1]] if best else []


def run(ds, i, j, keep_segs, seglength=None, carrier_segs=()):
    """`(IBD1Seg, IBD2Seg, D_ok)` printed by the reference for this pair."""
    key = "%s|%d,%d|%s|%s|%s" % (ds.name, i, j,
                                 ",".join("%d-%d" % s for s in keep_segs), seglength,
                                 ",".join("%d-%d" % s for s in carrier_segs))
    m = measured()
    if key in m:
        return m[key]
    with tempfile.TemporaryDirectory() as tmp:
        pre = os.path.join(tmp, "L")
        write_fileset(ds, i, j, keep_segs, pre, carrier_segs)
        args = [KING, "-b", pre + ".bed", "--ibdseg", "--cpus", "1", "--prefix", "P"]
        if seglength is not None:
            args += ["--seglength", "%.6f" % (seglength / 1e6)]
        subprocess.run(args, cwd=tmp, check=True, capture_output=True)
        segs = _allsegs_digest(os.path.join(tmp, "Pallsegs.txt"))
        path = os.path.join(tmp, "P.seg")
        val = None
        if os.path.exists(path):
            with open(path) as fh:
                head = fh.readline().rstrip("\n").split("\t")
                c1, c2 = head.index("IBD1Seg"), head.index("IBD2Seg")
                want = {ds.fam[i][1], ds.fam[j][1]}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    if {f[1], f[3]} == want:
                        val = (f[c1], f[c2])
        m[key] = [val, segs]
        save()
    return m[key]


def _allsegs_digest(path):
    if not os.path.exists(path):
        return None
    return hashlib.md5(open(path, "rb").read()).hexdigest()


_BASE_DIGEST = {}


def base_digest(ds):
    """`kingallsegs.txt` of the untouched fileset — the control the rig must preserve."""
    if ds.name not in _BASE_DIGEST:
        p = os.path.join(kd.GOLDEN, "ibdseg", ds.name + "__ibdseg", "kingallsegs.txt")
        _BASE_DIGEST[ds.name] = hashlib.md5(open(p, "rb").read()).hexdigest()
    return _BASE_DIGEST[ds.name]


def ours(ds, i, j, seg, min_bp=E.SEGLEN):
    """Our IBD2 calls inside one usable segment, and their total base pairs."""
    sc = E.SegScan(ds, i, j, seg, E.BASE)
    if sc.n == 0:
        return [], 0
    c = (sc.ibd2(ds.pos, min_bp) if R.RULE is None
         else __import__("seg19").ibd2_19(sc, ds, i, j, R.RULE, ds.pos, min_bp))
    return c, sum(int(ds.pos[b] - ds.pos[a]) for a, b in c)


def per_segment(name, n1, n2, only=None, min_bp=E.SEGLEN, quiet=False):
    """Reference IBD2 total against ours, segment by segment. Returns the disagreements."""
    ds = kd.load(name)
    idx = {f[1]: k for k, f in enumerate(ds.fam)}
    i, j = sorted((idx[n1], idx[n2]))
    pos, d = ds.pos, ds.denom
    ulp = d / 10000.0
    bad = []
    if not quiet:
        print("=== %s %s/%s   D = %d, 1 ulp = %.0f bp, median gap = %d bp"
              % (name, n1, n2, d, ulp, R._median_gap(ds)))
    for k, seg in enumerate(ds.segs):
        if only is not None and k != only:
            continue
        c, tot = ours(ds, i, j, seg, min_bp)
        val, dg = run(ds, i, j, [(seg[1], seg[2])],
                      None if min_bp == E.SEGLEN else min_bp, pick_carrier(ds, k))
        if dg != base_digest(ds):
            print("  !! segment %d: allsegs changed - instrument invalid" % k)
            continue
        if val is None:
            print("  segment %2d chr%-2d: pair not reported" % (k, seg[0]))
            continue
        ref2 = float(val[1]) * d
        delta = tot - ref2
        flag = "" if abs(delta) <= 0.5 * ulp else "   <<<"
        if flag:
            bad.append((k, seg, c, tot, val[1], delta))
        if not quiet:
            print("  seg %2d chr%-2d w%d..%d  ref %s (%11.0f bp)  ours %11d bp  "
                  "d=%+8.0f (%+5.2f ulp)  %d calls%s"
                  % (k, seg[0], -(-seg[1] // WORD), (seg[2] + 1) // WORD - 1,
                     val[1], ref2, tot, delta, delta / ulp, len(c), flag))
    return bad, (ds, i, j)


def ladder(ds, i, j, seg, lo=3_000_000, hi=10_000_000):
    """The reference's own call lengths inside one segment, by bisecting --seglength."""
    keep = [(seg[1], seg[2])]

    def v(bp):
        return run(ds, i, j, keep, bp)[0]

    found = []
    stack = [(lo, v(lo), hi, v(hi))]
    while stack:
        a, va, b, vb = stack.pop()
        if va == vb:
            continue
        if b - a <= 1:
            found.append(b - 1)
            continue
        m = (a + b) // 2
        vm = v(m)
        stack.append((a, va, m, vm))
        stack.append((m, vm, b, vb))
    return sorted(found)


def check(names=("dups", "multifam", "nuclear")):
    """The instrument's three controls.

    1. `kingallsegs.txt` is byte-identical to the untouched run, on every fileset written.
    2. Keeping every segment reproduces the pair's corpus row exactly.
    3. **Additivity** — the per-segment values sum to the whole-genome value.  This is the
       one that matters: it is the test of whether repainting the markers *outside* a
       segment disturbs the calls *inside* it, which is exactly what would invalidate the
       localisation. It is checked with the neutral fill and again with the mismatch-heavy
       carrier fill, two different contaminations of the flanking words.
    """
    for name in names:
        ds = kd.load(name)
        for (i, j) in list(ds.ref)[:3]:
            allsegs = [(s[1], s[2]) for s in ds.segs]
            val, dg = run(ds, i, j, allsegs)
            ref = ds.ref[(i, j)]
            print("%-10s %-8s/%-8s  whole-genome keep: ref %s/%s  corpus %.4f/%.4f  %s"
                  % (name, ds.fam[i][1], ds.fam[j][1], val[0], val[1], ref[0], ref[1],
                     "allsegs OK" if dg == base_digest(ds) else "ALLSEGS MOVED"))
            for tag, carr in (("neutral", False), ("carrier", True)):
                tot = 0.0
                ok = True
                for k, seg in enumerate(ds.segs):
                    cs = pick_carrier(ds, k) if carr else []
                    v, d2 = run(ds, i, j, [(seg[1], seg[2])], None, cs)
                    ok &= d2 == base_digest(ds)
                    if v is not None:
                        tot += float(v[1])
                print("     %-8s fill: sum of segments IBD2Seg = %.4f   (row %.4f)  %s"
                      % (tag, tot, ref[1], "allsegs OK" if ok else "ALLSEGS MOVED"))


if __name__ == "__main__":
    argv = sys.argv[1:]
    if "--check" in argv:
        check()
        sys.exit()
    only = None
    if "--seg" in argv:
        k = argv.index("--seg")
        only = int(argv[k + 1])
        del argv[k:k + 2]
    lad = "--ladder" in argv
    if lad:
        argv.remove("--ladder")
    bad, (ds, i, j) = per_segment(argv[0], argv[1], argv[2], only)
    if lad:
        for k, seg, c, tot, ref, delta in bad:
            print("  -- segment %d ladder" % k)
            print("     reference call lengths < 10 Mb: %s"
                  % [round(v / 1e6, 3) for v in ladder(ds, i, j, seg)])
            print("     ours: %s"
                  % [round(float(ds.pos[b] - ds.pos[a]) / 1e6, 3) for a, b in c])
