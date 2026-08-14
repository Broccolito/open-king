"""Loader + shared plumbing for the IBD-segment rule fit.

Reads the parity corpus (`tests/parity/work/data/*.bed|bim|fam`) and the captured
reference answers (`tests/parity/golden/ibdseg/<ds>__ibdseg/{king.seg,kingallsegs.txt}`)
into numpy structures the candidate segment callers in `rules.py` operate on.

Nothing here is a candidate rule: this module only encodes facts that are already
byte-verified elsewhere (the PLINK bit-plane encoding, KING's retained-autosome set, and
the usable-segment list, which is read straight out of the reference's own
`allsegs.txt` rather than recomputed).
"""

import os
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))
DATA = os.path.join(ROOT, "tests", "parity", "work", "data")
GOLDEN = os.path.join(ROOT, "tests", "parity", "golden")

WORD = 64

# Datasets whose --ibdseg capture actually contains a .seg (>= 5 samples).
DATASETS = ["nuclear", "threegen", "multifam", "dups", "missing", "monomorphic",
            "sexchr", "unrelated", "admixed", "bigish"]


def king_chrom_code(s):
    """KING's chromosome mapping: numeric 1..26, plus X/Y/XY/MT. None = variant dropped."""
    try:
        n = int(s)
        return n if 1 <= n <= 26 else None
    except ValueError:
        pass
    return {"X": 23, "x": 23, "Y": 24, "y": 24,
            "XY": 25, "xy": 25, "MT": 26, "mt": 26}.get(s)


def read_bim(path):
    chrom, name, pos = [], [], []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            chrom.append(king_chrom_code(f[0]))
            name.append(f[1])
            pos.append(int(f[3]))
    return chrom, name, np.asarray(pos, dtype=np.int64)


def read_fam(path):
    out = []
    with open(path) as fh:
        for line in fh:
            f = line.split()
            out.append((f[0], f[1]))
    return out


def read_bed_planes(path, n_samples, keep):
    """Bit planes over the retained variants, sample-major, 64 markers per uint64 word.

    Returns (plane0, plane1), each (n_samples, n_words) uint64.  plane0 = "homozygous",
    plane1 = "carries A1 when homozygous / het" — the encoding king-io/bed.rs uses:
    code 00 -> (1,1), 01 (missing) -> (0,0), 10 (het) -> (0,1), 11 -> (1,0).
    """
    bpv = (n_samples + 3) // 4
    raw = np.fromfile(path, dtype=np.uint8)
    assert raw[0] == 0x6C and raw[1] == 0x1B and raw[2] == 1, "not a SNP-major PLINK bed"
    body = raw[3:]
    n_var = body.size // bpv
    body = body.reshape(n_var, bpv)
    codes = np.empty((n_var, bpv * 4), dtype=np.uint8)
    for k in range(4):
        codes[:, k::4] = (body >> (2 * k)) & 3
    codes = codes[:, :n_samples]
    codes = codes[keep]                       # retained variants only
    p0 = (codes == 0) | (codes == 3)
    p1 = (codes == 0) | (codes == 2)
    return _pack(p0.T), _pack(p1.T)


def _pack(bits):
    """(n_samples, n_markers) bool -> (n_samples, n_words) uint64, bit k = marker 64w+k."""
    n, m = bits.shape
    pad = (-m) % WORD
    if pad:
        bits = np.hstack([bits, np.zeros((n, pad), dtype=bool)])
    packed = np.packbits(bits, axis=1, bitorder="little")
    return packed.view(np.uint64)


class Dataset:
    """One corpus fileset plus the reference's answers for it."""

    def __init__(self, name):
        self.name = name
        chrom, snpname, pos = read_bim(os.path.join(DATA, name + ".bim"))
        self.fam = read_fam(os.path.join(DATA, name + ".fam"))
        n = len(self.fam)
        keep = np.array([c is not None and (1 <= c <= 22 or c == 25) for c in chrom])
        self.keep = keep
        self.pos = pos[keep]
        self.chr = np.array([c for c, k in zip(chrom, keep) if k], dtype=np.int64)
        self.names = [s for s, k in zip(snpname, keep) if k]
        self.index = {s: i for i, s in enumerate(self.names)}
        self.p0, self.p1 = read_bed_planes(os.path.join(DATA, name + ".bed"), n, keep)
        self.nwords = self.p0.shape[1]
        self.segs = self._read_allsegs()
        self.denom = int(sum(self.pos[hi] - self.pos[lo] for _, lo, hi in self.segs))
        self.ref = self._read_seg()

    def _case(self, suffix="__ibdseg"):
        return os.path.join(GOLDEN, "ibdseg", self.name + suffix)

    def _read_allsegs(self):
        """Usable segments as (chr, lo, hi) marker-index triples, autosomes only.

        Read from the reference's own allsegs.txt: it is reproduced byte-for-byte by the
        Rust engine already, so re-deriving it here would only add a failure mode.
        """
        out = []
        with open(os.path.join(self._case(), "kingallsegs.txt")) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                chrom = int(f[1])
                if chrom == 23:
                    continue
                out.append((chrom, self.index[f[6]], self.index[f[7]]))
        return out

    def _read_seg(self, suffix="__ibdseg"):
        """Reference .seg rows: {(i, j): (ibd1, ibd2, prop, inftype)} by sample index."""
        pos = {(f, i): k for k, (f, i) in enumerate(self.fam)}
        out = {}
        with open(os.path.join(self._case(suffix), "king.seg")) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip("\n").split("\t")
                i, j = pos[(f[0], f[1])], pos[(f[2], f[3])]
                out[(min(i, j), max(i, j))] = (float(f[4]), float(f[5]), float(f[6]), f[7])
        return out

    def pairs(self):
        n = len(self.fam)
        return [(i, j) for i in range(n) for j in range(i + 1, n)]

    def masks(self, i, j):
        """(ibs0, ibs1, hethet, both_called) per-word uint64 masks for one pair."""
        p0i, p1i = self.p0[i], self.p1[i]
        p0j, p1j = self.p0[j], self.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        ibs0 = p0i & p0j & (p1i ^ p1j)
        ibs1 = (het_i & p0j) | (p0i & het_j)
        hethet = het_i & het_j
        called = (p0i | p1i) & (p0j | p1j)
        return ibs0, ibs1, hethet, called


_CACHE = {}


def load(name):
    if name not in _CACHE:
        _CACHE[name] = Dataset(name)
    return _CACHE[name]


def fmt4(x):
    """KING's %.4lf, so comparisons are made on the printed value."""
    return float("%.4f" % x)


def inf_type(pi1, pi2, prop):
    if pi2 > 0.7:
        return "Dup/MZ"
    if pi1 + pi2 > 0.96 or (pi1 + pi2 > 0.9 and pi2 < 0.08):
        return "PO"
    if prop > 0.3535533905932738 and pi2 >= 0.08:
        return "FS"
    if prop > 0.1767766952966369:
        return "2nd"
    if prop > 0.08838834764831845:
        return "3rd"
    if prop > 0.04419417382415922:
        return "4th"
    return "UN"
