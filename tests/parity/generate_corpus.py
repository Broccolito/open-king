#!/usr/bin/env python3
"""
open-king differential-test corpus generator.

Synthesises PLINK 1 binary filesets (.bed/.bim/.fam) directly -- no external
tools, no third-party Python packages, standard library only.  Everything is
derived from a single master seed, so the corpus is byte-for-byte reproducible:
the genotype files are NOT committed, this script is.

Usage
-----
    generate_corpus.py --outdir DIR [--seed 20260813] [--only NAME ...]

Genetics
--------
Founder haplotypes are drawn allele-by-allele from a per-SNP allele frequency.
Non-founders inherit one haplotype from each parent through a simulated meiosis
with recombination (Haldane map function over the .bim genetic map), so IBD
segments are real, contiguous and detectable -- not per-SNP coin flips.
Sex chromosomes follow their real transmission rules (X: sons from mother only,
fathers pass their single X to daughters; Y: father to son; MT: maternal), and
haploid genotypes are emitted as homozygous, which is PLINK's convention.

Encoding
--------
.bed is SNP-major (magic 6c 1b 01).  Allele bit 1 is the .bim A1 allele (column
5, simulated as the minor allele), bit 0 is A2 (column 6).  With g = the number
of A1 copies:

    g == 2  hom A1   -> 0b00
    g == 1  het      -> 0b10
    g == 0  hom A2   -> 0b11
    missing          -> 0b01

Two bits per sample, four samples per byte, the first sample in the low bits.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import random
import sys

GENERATOR_VERSION = 1

BED_MAGIC = bytes((0x6C, 0x1B, 0x01))

# Indexed by A1 dosage: 0,1,2 copies of A1, and 3 == missing.
BED_CODE = (0b11, 0b10, 0b00, 0b01)

MISSING = 3

# GRCh38 chromosome lengths.  24 == Y, 25 == XY (pseudoautosomal), 26 == MT.
CHROM_LEN = {
    1: 248956422, 2: 242193529, 3: 198295559, 4: 190214555, 5: 181538259,
    6: 170805979, 7: 159345973, 8: 145138636, 9: 138394717, 10: 133797422,
    11: 135086622, 12: 133275309, 13: 114364328, 14: 107043718, 15: 101991189,
    16: 90338345, 17: 83257441, 18: 80373285, 19: 58617616, 20: 64444167,
    21: 46709983, 22: 50818468, 23: 156040895, 24: 57227415, 25: 2781479,
    26: 16569,
}

AUTOSOMES = list(range(1, 23))

ALLELE_PAIRS = (("A", "G"), ("C", "T"), ("A", "C"), ("G", "T"),
                ("A", "T"), ("C", "G"))

# Nominal inter-SNP spacing; shrunk automatically on short chromosomes.
SPACING_BP = 50_000

# Genetic map: 1 cM per Mb.  .bim column 3 is written in centimorgans.
CM_PER_BP = 1.0 / 1_000_000.0


# --------------------------------------------------------------------------
# pedigree model
# --------------------------------------------------------------------------

class Person:
    __slots__ = ("idx", "fid", "iid", "sex", "pheno", "fidx", "midx",
                 "pop", "alpha", "clone_of", "error_rate", "emit", "tag")

    def __init__(self, idx, fid, iid, sex, pheno, fidx, midx,
                 pop, alpha, clone_of, error_rate, emit, tag):
        self.idx = idx
        self.fid = fid
        self.iid = iid
        self.sex = sex
        self.pheno = pheno
        self.fidx = fidx
        self.midx = midx
        self.pop = pop
        self.alpha = alpha
        self.clone_of = clone_of
        self.error_rate = error_rate
        self.emit = emit
        self.tag = tag

    @property
    def key(self):
        return "%s:%s" % (self.fid, self.iid)


class Ped:
    """Builder for a genetic pedigree.

    Individuals may be marked ``emit=False``; they take part in transmission
    and in the expected-kinship arithmetic but are not written to the .fam.
    That is how cross-family relatives get created without declaring them.
    """

    def __init__(self):
        self.people = []
        self._by_iid = {}

    def add(self, iid, fid, father=None, mother=None, sex=0, pheno=-9,
            pop=0, alpha=None, clone_of=None, error_rate=0.0, emit=True,
            tag=None):
        if iid in self._by_iid:
            raise ValueError("duplicate individual id: %s" % iid)
        idx = len(self.people)
        p = Person(
            idx=idx, fid=fid, iid=iid, sex=sex, pheno=pheno,
            fidx=self._by_iid[father].idx if father else None,
            midx=self._by_iid[mother].idx if mother else None,
            pop=pop, alpha=alpha,
            clone_of=self._by_iid[clone_of].idx if clone_of else None,
            error_rate=error_rate, emit=emit, tag=tag,
        )
        self.people.append(p)
        self._by_iid[iid] = p
        return iid

    def get(self, iid):
        return self._by_iid[iid]

    def n_emitted(self):
        return sum(1 for p in self.people if p.emit)

    def emitted(self):
        return [p for p in self.people if p.emit]


def topo_order(people):
    """Indices ordered so parents (and clone sources) precede dependants."""
    state = [0] * len(people)
    order = []

    def visit(i):
        if state[i] == 2:
            return
        if state[i] == 1:
            raise ValueError("pedigree loop at %s" % people[i].iid)
        state[i] = 1
        p = people[i]
        for dep in (p.fidx, p.midx, p.clone_of):
            if dep is not None:
                visit(dep)
        state[i] = 2
        order.append(i)

    for i in range(len(people)):
        visit(i)
    return order


def kinship_matrix(people):
    """Exact pedigree kinship over every individual, phantoms included."""
    n = len(people)
    phi = [[0.0] * n for _ in range(n)]
    for i in topo_order(people):
        p = people[i]
        if p.clone_of is not None:
            s = p.clone_of
            for j in range(n):
                if j == i:
                    continue
                v = phi[s][s] if j == s else phi[s][j]
                phi[i][j] = phi[j][i] = v
            phi[i][i] = phi[s][s]
            continue
        f, m = p.fidx, p.midx
        for j in range(n):
            if j == i:
                continue
            a = phi[f][j] if f is not None else 0.0
            b = phi[m][j] if m is not None else 0.0
            v = 0.5 * (a + b)
            if v:
                phi[i][j] = phi[j][i] = v
        if f is not None and m is not None:
            phi[i][i] = 0.5 * (1.0 + phi[f][m])
        else:
            phi[i][i] = 0.5
    return phi


def _parents(people, i):
    p = people[i]
    return [x for x in (p.fidx, p.midx) if x is not None]


def _full_sibs(people, a, b):
    if a == b:
        return False
    pa, pb = _parents(people, a), _parents(people, b)
    return len(pa) == 2 and len(pb) == 2 and set(pa) == set(pb)


def _half_sibs(people, a, b):
    if a == b:
        return False
    pa, pb = set(_parents(people, a)), set(_parents(people, b))
    return len(pa & pb) == 1


def _ancestors(people, i, depth):
    cur = {i}
    for _ in range(depth):
        nxt = set()
        for x in cur:
            nxt.update(_parents(people, x))
        cur = nxt
    return cur


def classify(people, i, j, phi):
    """Human-readable relationship label for a pair."""
    a, b = people[i], people[j]
    if a.clone_of == j or b.clone_of == i:
        src = people[a.clone_of if a.clone_of == j else b.clone_of]
        other = a if a.clone_of is not None else b
        return "MZ" if other.error_rate > 0.0 else "DUP"
    if j in (a.fidx, a.midx) or i in (b.fidx, b.midx):
        return "PO"
    if _full_sibs(people, i, j):
        return "FS"
    if _half_sibs(people, i, j):
        return "HS"
    if j in _ancestors(people, i, 2) or i in _ancestors(people, j, 2):
        return "GG"
    if j in _ancestors(people, i, 3) or i in _ancestors(people, j, 3):
        return "GGG"
    for x, y in ((i, j), (j, i)):
        for par in _parents(people, x):
            if _full_sibs(people, par, y):
                return "AV"
    for x, y in ((i, j), (j, i)):
        for par in _parents(people, x):
            if _half_sibs(people, par, y):
                return "HAV"
    for pi in _parents(people, i):
        for pj in _parents(people, j):
            if _full_sibs(people, pi, pj):
                return "FC"
    for pi in _parents(people, i):
        for pj in _parents(people, j):
            if _half_sibs(people, pi, pj):
                return "HFC"
    k = phi[i][j]
    if k <= 0.0:
        return "UN"
    deg = int(round(-math.log(k, 2.0))) - 1
    return "%ddeg" % max(deg, 1)


# --------------------------------------------------------------------------
# pedigree shorthands
# --------------------------------------------------------------------------

def add_couple(ped, fid, prefix, emit=True, pop=0):
    f = ped.add(prefix + "_F", fid, sex=1, emit=emit, pop=pop)
    m = ped.add(prefix + "_M", fid, sex=2, emit=emit, pop=pop)
    return f, m


def add_nuclear(ped, fid, prefix, n_kids, father_parents=None,
                mother_parents=None, pheno_parents=-9, pheno_kids=-9, pop=0):
    fa = ped.add(prefix + "_F", fid, sex=1, pheno=pheno_parents, pop=pop,
                 father=father_parents[0] if father_parents else None,
                 mother=father_parents[1] if father_parents else None)
    mo = ped.add(prefix + "_M", fid, sex=2, pheno=pheno_parents, pop=pop,
                 father=mother_parents[0] if mother_parents else None,
                 mother=mother_parents[1] if mother_parents else None)
    kids = []
    for k in range(n_kids):
        kids.append(ped.add("%s_C%d" % (prefix, k + 1), fid, father=fa,
                            mother=mo, sex=1 + (k % 2), pheno=pheno_kids))
    return fa, mo, kids


def add_threegen9(ped, fid, prefix):
    """Nine-person three-generation unit: PO, FS, GG, AV and FC pairs."""
    gf, gm = add_couple(ped, fid, prefix + "_G")
    p1 = ped.add(prefix + "_P1", fid, father=gf, mother=gm, sex=1)
    p2 = ped.add(prefix + "_P2", fid, father=gf, mother=gm, sex=2)
    s1 = ped.add(prefix + "_S1", fid, sex=2)
    s2 = ped.add(prefix + "_S2", fid, sex=1)
    ped.add(prefix + "_C1", fid, father=p1, mother=s1, sex=1)
    ped.add(prefix + "_C2", fid, father=p1, mother=s1, sex=2)
    ped.add(prefix + "_C3", fid, father=s2, mother=p2, sex=1)
    return 9


# --------------------------------------------------------------------------
# SNP map
# --------------------------------------------------------------------------

def allocate_counts(chroms, n_snps):
    """Largest-remainder split of n_snps across chroms, proportional to length."""
    total = float(sum(CHROM_LEN[c] for c in chroms))
    exact = {c: n_snps * CHROM_LEN[c] / total for c in chroms}
    base = {c: int(math.floor(exact[c])) for c in chroms}
    left = n_snps - sum(base.values())
    for c in sorted(chroms, key=lambda c: (-(exact[c] - base[c]), c))[:left]:
        base[c] += 1
    # Never leave a requested chromosome empty.
    for c in chroms:
        if base[c] == 0:
            donor = max(chroms, key=lambda x: base[x])
            if base[donor] > 1:
                base[donor] -= 1
                base[c] = 1
    return base


def build_map(chrom_spec, n_snps, rnd):
    """Return [(chrom, bp, cm), ...] sorted by chromosome then position."""
    if isinstance(chrom_spec, dict):
        counts = dict(chrom_spec)
        chroms = sorted(counts)
    else:
        chroms = sorted(chrom_spec)
        counts = allocate_counts(chroms, n_snps)
    snps = []
    for c in chroms:
        m = counts[c]
        if m <= 0:
            continue
        length = CHROM_LEN[c]
        margin = min(1_000_000, length // 10)
        usable = length - 2 * margin
        step = max(1, min(SPACING_BP, usable // m))
        jitter = max(1, step // 5)
        for k in range(m):
            bp = margin + k * step + rnd.randrange(jitter)
            snps.append((c, bp, bp * CM_PER_BP))
    return snps


# --------------------------------------------------------------------------
# per-SNP allele / frequency models
# --------------------------------------------------------------------------

def _alleles(rnd):
    a1, a2 = ALLELE_PAIRS[rnd.randrange(len(ALLELE_PAIRS))]
    return (a1, a2) if rnd.random() < 0.5 else (a2, a1)


def model_common(rnd, n_snps, snps):
    """MAF ~ U(0.05, 0.45) on A1, one population."""
    out = []
    for _ in range(n_snps):
        a1, a2 = _alleles(rnd)
        out.append({"freqs": (rnd.uniform(0.05, 0.45),),
                    "a1": a1, "a2": a2, "force": None})
    return out


def model_monomorphic(rnd, n_snps, snps):
    """Mix of ordinary, monomorphic, missing-allele and ultra-rare SNPs."""
    out = []
    for j in range(n_snps):
        a1, a2 = _alleles(rnd)
        cat = j % 10
        force = None
        if cat in (0, 1, 2, 3, 4):            # 50% ordinary
            p = rnd.uniform(0.05, 0.45)
        elif cat == 5:                        # monomorphic, both alleles named
            p = 0.0
        elif cat == 6:                        # monomorphic, A1 column is "0"
            p = 0.0
            a1 = "0"
        elif cat == 7:                        # MAF 0.001, usually absent
            p = 0.001
        elif cat == 8:                        # MAF 0.001, at least one carrier
            p = 0.001
            force = "singleton"
        else:                                 # exactly 0.5
            p = 0.5
        out.append({"freqs": (p,), "a1": a1, "a2": a2, "force": force})
    return out


def make_model_bn(fst, n_pops=2):
    """Balding-Nichols two-population frequencies at the given Fst."""
    def model(rnd, n_snps, snps):
        c = (1.0 - fst) / fst
        out = []
        for _ in range(n_snps):
            a1, a2 = _alleles(rnd)
            p = rnd.uniform(0.05, 0.95)
            fr = tuple(min(0.999, max(0.001, rnd.betavariate(p * c, (1.0 - p) * c)))
                       for _ in range(n_pops))
            out.append({"freqs": fr, "a1": a1, "a2": a2, "force": None})
        return out
    return model


def no_snp_missing(rnd, n_snps):
    return [0.0] * n_snps


def make_snp_missing(n_high, high_lo, high_hi, n_all_missing):
    """Sprinkle high-missingness SNPs and a few with zero call rate."""
    def plan(rnd, n_snps):
        rates = [0.0] * n_snps
        picked = set()
        while len(picked) < min(n_high, n_snps):
            picked.add(rnd.randrange(n_snps))
        for j in sorted(picked):
            rates[j] = rnd.uniform(high_lo, high_hi)
        for k in range(n_all_missing):
            rates[(k * 977 + 13) % n_snps] = 1.0
        return rates
    return plan


# --------------------------------------------------------------------------
# dataset specification
# --------------------------------------------------------------------------

class Spec:
    def __init__(self, name, ped, chrom_spec, n_snps, snp_model=model_common,
                 snp_missing=no_snp_missing, sample_missing=None, notes=""):
        self.name = name
        self.ped = ped
        self.chrom_spec = chrom_spec
        self.n_snps = n_snps
        self.snp_model = snp_model
        self.snp_missing = snp_missing
        self.sample_missing = sample_missing or {}
        self.notes = notes


def dataset_seed(master_seed, name):
    h = hashlib.sha256(("%d:%s" % (master_seed, name)).encode()).digest()
    return int.from_bytes(h[:8], "big")


# --------------------------------------------------------------------------
# simulation + writing
# --------------------------------------------------------------------------

def write_fam(path, people):
    with open(path, "w") as fh:
        for p in people:
            if p.emit:
                fh.write("%s %s %s %s %d %s\n" % (
                    p.fid, p.iid,
                    people[p.fidx].iid if p.fidx is not None and people[p.fidx].emit else "0",
                    people[p.midx].iid if p.midx is not None and people[p.midx].emit else "0",
                    p.sex, p.pheno))


def write_bim(path, snps, snpinfo):
    with open(path, "w") as fh:
        for j, (chrom, bp, cm) in enumerate(snps):
            info = snpinfo[j]
            fh.write("%d rs%d_%d %.6f %d %s %s\n" % (
                chrom, chrom, bp, cm, bp, info["a1"], info["a2"]))


def simulate(spec, seed, outdir, emit_freq=False):
    rnd = random.Random(seed)
    rr = rnd.random
    people = spec.ped.people
    n_all = len(people)
    order = topo_order(people)
    emit_idx = [p.idx for p in people if p.emit]
    n_emit = len(emit_idx)

    snps = build_map(spec.chrom_spec, spec.n_snps, rnd)
    n_snps = len(snps)
    snpinfo = spec.snp_model(rnd, n_snps, snps)
    snp_miss = spec.snp_missing(rnd, n_snps)
    samp_miss = [spec.sample_missing.get(people[i].iid, 0.0) for i in emit_idx]

    base = os.path.join(outdir, spec.name)
    write_fam(base + ".fam", people)
    write_bim(base + ".bim", snps, snpinfo)

    hap_a = [0] * n_all
    hap_b = [0] * n_all
    sel_f = [0] * n_all
    sel_m = [0] * n_all
    geno = [0] * n_all
    n_bytes = (n_emit + 3) // 4
    obs = []

    with open(base + ".bed", "wb") as fh:
        fh.write(BED_MAGIC)
        prev_chrom = None
        prev_cm = 0.0
        for j in range(n_snps):
            chrom, bp, cm = snps[j]
            new_chrom = chrom != prev_chrom
            if new_chrom:
                recomb = 0.0
            else:
                recomb = 0.5 * (1.0 - math.exp(-2.0 * (cm - prev_cm) / 100.0))
            prev_chrom, prev_cm = chrom, cm
            info = snpinfo[j]
            fr = info["freqs"]
            is_x = chrom == 23
            is_y = chrom == 24
            is_mt = chrom == 26

            for i in order:
                p = people[i]

                if p.clone_of is not None:
                    s = p.clone_of
                    hap_a[i] = hap_a[s]
                    hap_b[i] = hap_b[s]
                    g = geno[s]
                    if p.error_rate > 0.0 and g != MISSING and rr() < p.error_rate:
                        g = [x for x in (0, 1, 2) if x != g][rnd.randrange(2)]
                    geno[i] = g
                    continue

                # founder allele draw, honouring population / admixture
                if p.alpha is None:
                    q = fr[p.pop]

                    def draw(_q=q):
                        return 1 if rr() < _q else 0
                else:
                    alpha = p.alpha
                    f0, f1 = fr[0], fr[1]

                    def draw(_a=alpha, _f0=f0, _f1=f1):
                        return 1 if rr() < (_f0 if rr() < _a else _f1) else 0

                f, m = p.fidx, p.midx
                if f is not None or m is not None:
                    if new_chrom:
                        sel_f[i] = rnd.getrandbits(1)
                        sel_m[i] = rnd.getrandbits(1)
                    else:
                        if rr() < recomb:
                            sel_f[i] ^= 1
                        if rr() < recomb:
                            sel_m[i] ^= 1

                if is_y:
                    if p.sex != 1:
                        hap_a[i] = hap_b[i] = 0
                        geno[i] = MISSING
                        continue
                    a = hap_a[f] if f is not None else draw()
                    hap_a[i] = hap_b[i] = a
                    geno[i] = 2 * a
                    continue

                if is_mt:
                    a = hap_a[m] if m is not None else draw()
                    hap_a[i] = hap_b[i] = a
                    geno[i] = 2 * a
                    continue

                if is_x and p.sex == 1:
                    if m is not None:
                        a = hap_a[m] if sel_m[i] == 0 else hap_b[m]
                    else:
                        a = draw()
                    hap_a[i] = hap_b[i] = a
                    geno[i] = 2 * a
                    continue

                if is_x:
                    # female / unknown sex: father contributes his single X
                    a = hap_a[f] if f is not None else draw()
                    if m is not None:
                        b = hap_a[m] if sel_m[i] == 0 else hap_b[m]
                    else:
                        b = draw()
                else:
                    if f is not None:
                        a = hap_a[f] if sel_f[i] == 0 else hap_b[f]
                    else:
                        a = draw()
                    if m is not None:
                        b = hap_a[m] if sel_m[i] == 0 else hap_b[m]
                    else:
                        b = draw()
                hap_a[i] = a
                hap_b[i] = b
                geno[i] = a + b

            row = [geno[i] for i in emit_idx]

            if info["force"] == "singleton" and not any(g in (1, 2) for g in row):
                row[j % n_emit] = 1

            snp_rate = snp_miss[j]
            for k in range(n_emit):
                sr = samp_miss[k]
                pe = 1.0 - (1.0 - sr) * (1.0 - snp_rate)
                if pe > 0.0 and rr() < pe:
                    row[k] = MISSING

            called = 0
            a1n = 0
            for g in row:
                if g != MISSING:
                    called += 1
                    a1n += g
            obs.append((a1n, 2 * called))

            packed = bytearray(n_bytes)
            k = 0
            for byte_i in range(n_bytes):
                v = 0
                for shift in (0, 2, 4, 6):
                    if k < n_emit:
                        v |= BED_CODE[row[k]] << shift
                        k += 1
                packed[byte_i] = v
            fh.write(packed)

    if emit_freq:
        with open(base + ".expected_freq.tsv", "w") as fh:
            fh.write("SNP\tA1\tA2\tSIM_A1_FREQ\tOBS_A1_FREQ\tNCHROBS\n")
            for j in range(n_snps):
                chrom, bp, _cm = snps[j]
                a1n, nch = obs[j]
                info = snpinfo[j]
                sim = info["freqs"][0]
                fh.write("rs%d_%d\t%s\t%s\t%.6f\t%.6f\t%d\n" % (
                    chrom, bp, info["a1"], info["a2"], sim,
                    (a1n / nch) if nch else 0.0, nch))

    return snps, n_snps, n_emit


# --------------------------------------------------------------------------
# datasets
# --------------------------------------------------------------------------

def build_trio():
    ped = Ped()
    add_nuclear(ped, "TRIO", "T", 1)
    return Spec("trio", ped, [1, 2, 3], 5000,
                notes="Father, mother, one child. No missingness.")


def build_nuclear():
    ped = Ped()
    add_nuclear(ped, "NUC", "N", 4, pheno_parents=1, pheno_kids=2)
    return Spec("nuclear", ped, [1, 2, 3, 4, 5], 10000,
                notes="Father, mother, four full sibs. Case/control phenotypes.")


def build_threegen():
    ped = Ped()
    fid = "TG"
    gf = ped.add("TG_GF", fid, sex=1)
    gm1 = ped.add("TG_GM1", fid, sex=2)
    gm2 = ped.add("TG_GM2", fid, sex=2)
    p1 = ped.add("TG_P1", fid, father=gf, mother=gm1, sex=1)
    p2 = ped.add("TG_P2", fid, father=gf, mother=gm1, sex=2)
    p3 = ped.add("TG_P3", fid, father=gf, mother=gm2, sex=1)
    s1 = ped.add("TG_S1", fid, sex=2)
    s2 = ped.add("TG_S2", fid, sex=1)
    ped.add("TG_C1", fid, father=p1, mother=s1, sex=1)
    ped.add("TG_C2", fid, father=p1, mother=s1, sex=2)
    ped.add("TG_C3", fid, father=s2, mother=p2, sex=1)
    ped.add("TG_C4", fid, father=s2, mother=p2, sex=2)
    ped.add("TG_C5", fid, father=p3, mother=ped.add("TG_S3", fid, sex=2), sex=1)
    return Spec("threegen", ped, AUTOSOMES, 20000,
                notes=("Three generations. Contains PO, FS, half sibs (2nd), "
                       "grandparent-grandchild (2nd), avuncular (2nd), "
                       "half-avuncular (3rd) and first cousins (3rd)."))


def build_multifam():
    ped = Ped()
    # Phantom couples: parents of cross-family sib pairs, never written out.
    ph1f, ph1m = add_couple(ped, "PH1", "PH1", emit=False)
    ph2f, ph2m = add_couple(ped, "PH2", "PH2", emit=False)

    # FAM1 father is a full sib of FAM2 father (undeclared).
    a_f, a_m, _ = add_nuclear(ped, "FAM1", "A", 3,
                              father_parents=(ph1f, ph1m))
    add_nuclear(ped, "FAM2", "B", 3, father_parents=(ph1f, ph1m))

    # FAM3 father is a child of the FAM1 couple (undeclared): PO + FS in .kin0.
    c_f = ped.add("C_F", "FAM3", father=a_f, mother=a_m, sex=1)
    c_m = ped.add("C_M", "FAM3", sex=2, father=ph2f, mother=ph2m)
    for k in range(3):
        ped.add("C_C%d" % (k + 1), "FAM3", father=c_f, mother=c_m,
                sex=1 + (k % 2))

    # FAM4 mother is a full sib of FAM3 mother (undeclared).
    add_nuclear(ped, "FAM4", "D", 3, mother_parents=(ph2f, ph2m))

    return Spec("multifam", ped, AUTOSOMES, 15000,
                notes=("Four declared families plus undeclared cross-family "
                       "relatives (PO, FS, avuncular, first cousins) so .kin0 "
                       "carries real relatives, not only unrelated pairs."))


def build_dups():
    ped = Ped()
    for k in range(4):
        ped.add("U%d" % (k + 1), "UNR%d" % (k + 1), sex=1 + (k % 2))
    ped.add("DUP_A", "DUPA", sex=1)
    ped.add("DUP_A_COPY", "DUPB", sex=1, clone_of="DUP_A")
    ped.add("MZ_1", "MZFAM", sex=2)
    ped.add("MZ_2", "MZFAM", sex=2, clone_of="MZ_1", error_rate=0.002)
    par = ped.add("PO_P", "POFAM", sex=1)
    mom = ped.add("PO_M", "POFAM", sex=2, emit=False)
    ped.add("PO_C", "POFAM", father=par, mother=mom, sex=2)
    return Spec("dups", ped, AUTOSOMES, 10000,
                notes=("Exact duplicate pair across two FIDs, an MZ-twin-like "
                       "pair within one FID at a 0.2%% per-genotype error rate, "
                       "one parent-offspring pair and four unrelated samples."))


def build_missing():
    ped = Ped()
    add_nuclear(ped, "MIS", "M", 4)
    rates = {"M_F": 0.00, "M_M": 0.01, "M_C1": 0.05, "M_C2": 0.20,
             "M_C3": 0.50, "M_C4": 0.00}
    return Spec("missing", ped, [1, 2, 3, 4, 5], 10000,
                snp_missing=make_snp_missing(300, 0.40, 0.90, 5),
                sample_missing=rates,
                notes=("Nuclear family with per-sample missingness of 0/1/5/20/"
                       "50%%, 300 high-missingness SNPs and 5 SNPs missing in "
                       "every sample."))


def build_monomorphic():
    ped = Ped()
    add_nuclear(ped, "MONO", "P", 4)
    for k in range(6):
        ped.add("MZERO%d" % (k + 1), "MZ%d" % (k + 1), sex=1 + (k % 2))
    return Spec("monomorphic", ped, [1, 2], 5000, snp_model=model_monomorphic,
                notes=("Every tenth SNP cycles a category: ordinary, "
                       "monomorphic with both alleles named, monomorphic with "
                       "A1 written as PLINK's missing allele '0', MAF 0.001, "
                       "MAF 0.001 with a forced singleton carrier, and MAF "
                       "exactly 0.5."))


def build_sexchr():
    ped = Ped()
    fid = "SEX"
    fa = ped.add("S_F", fid, sex=1)
    mo = ped.add("S_M", fid, sex=2)
    ped.add("S_SON1", fid, father=fa, mother=mo, sex=1)
    ped.add("S_SON2", fid, father=fa, mother=mo, sex=1)
    ped.add("S_DAU1", fid, father=fa, mother=mo, sex=2)
    ped.add("S_DAU2", fid, father=fa, mother=mo, sex=2)
    ped.add("S_U0A", "SU1", sex=0)
    ped.add("S_U0B", "SU2", sex=0)
    ped.add("S_UM", "SU3", sex=1)
    ped.add("S_UF", "SU4", sex=2)
    counts = {1: 2000, 2: 2000, 23: 1500, 24: 300, 25: 150, 26: 50}
    return Spec("sexchr", ped, counts, sum(counts.values()),
                notes=("Autosomes plus X (23), Y (24), XY (25) and MT (26). "
                       "Both sexes plus sex-0 individuals. Haploid genotypes "
                       "are emitted homozygous; Y is missing for non-males."))


def build_unrelated():
    ped = Ped()
    for k in range(10):
        ped.add("P%02d" % (k + 1), "POOL", sex=1 + (k % 2))
    for k in range(20):
        ped.add("S%02d" % (k + 1), "SNG%02d" % (k + 1), sex=1 + (k % 2))
    return Spec("unrelated", ped, AUTOSOMES, 20000,
                notes=("30 mutually unrelated founders. Ten share one FID (so "
                       ".kin is exercised with unrelated pairs), twenty are in "
                       "their own FIDs (.kin0). Expected kinship 0 throughout."))


def build_admixed():
    ped = Ped()
    for k in range(14):
        ped.add("A1_%02d" % (k + 1), "AP1_%02d" % (k + 1), sex=1 + (k % 2), pop=0)
    for k in range(14):
        ped.add("A2_%02d" % (k + 1), "AP2_%02d" % (k + 1), sex=1 + (k % 2), pop=1)
    for k, alpha in enumerate((0.10, 0.25, 0.40, 0.50, 0.60, 0.90)):
        ped.add("ADM_%d" % (k + 1), "ADM%d" % (k + 1), sex=1 + (k % 2), alpha=alpha)
    xf = ped.add("X_F", "FAMX", sex=1, pop=0)
    xm = ped.add("X_M", "FAMX", sex=2, pop=0)
    ped.add("X_C1", "FAMX", father=xf, mother=xm, sex=1)
    ped.add("X_C2", "FAMX", father=xf, mother=xm, sex=2)
    yf = ped.add("Y_F", "FAMY", sex=1, pop=0)
    ym = ped.add("Y_M", "FAMY", sex=2, pop=1)
    ped.add("Y_C1", "FAMY", father=yf, mother=ym, sex=1)
    ped.add("Y_C2", "FAMY", father=yf, mother=ym, sex=2)
    zf = ped.add("Z_F", "FAMZ", sex=1, pop=1)
    zm = ped.add("Z_M", "FAMZ", sex=2, pop=1)
    ped.add("Z_C1", "FAMZ", father=zf, mother=zm, sex=1)
    ped.add("Z_C2", "FAMZ", father=zf, mother=zm, sex=2)
    return Spec("admixed", ped, AUTOSOMES, 20000,
                snp_model=make_model_bn(0.10),
                notes=("Two Balding-Nichols populations at Fst 0.10 plus six "
                       "admixed founders and three nuclear families (within "
                       "pop1, cross-population, within pop2). KING-robust and "
                       "KING-homo diverge most here."))


def build_singleton():
    ped = Ped()
    ped.add("ONLY", "SOLO", sex=1)
    return Spec("singleton", ped, AUTOSOMES, 5000,
                notes="Single individual: no pairs at all.")


def build_pair():
    ped = Ped()
    pf, pm = add_couple(ped, "PHP", "PHP", emit=False)
    ped.add("PR_A", "PFAM1", father=pf, mother=pm, sex=1)
    ped.add("PR_B", "PFAM2", father=pf, mother=pm, sex=2)
    return Spec("pair", ped, AUTOSOMES, 5000,
                notes=("Exactly two individuals, in different FIDs, that are "
                       "genetically full sibs through undeclared parents: .kin "
                       "is empty and .kin0 holds one related pair."))


def build_bigish():
    ped = Ped()
    target = 200
    fam = 0
    prev_parents = None
    while True:
        fam += 1
        n_kids = 2 + (fam % 4)
        size = 9 if fam % 3 == 0 else 2 + n_kids
        if ped.n_emitted() + size > 190:
            break
        if fam % 3 == 0:
            add_threegen9(ped, "BF%02d" % fam, "B%02d" % fam)
            prev_parents = None
        else:
            # Every fourth family shares undeclared phantom grandparents with
            # the previous one, giving cross-FID relatives for .kin0.
            fp = None
            if fam % 4 == 2 and prev_parents is not None:
                fp = prev_parents
            elif fam % 4 == 1:
                fp = add_couple(ped, "BPH%02d" % fam, "BPH%02d" % fam, emit=False)
                prev_parents = fp
            fa, mo, _ = add_nuclear(ped, "BF%02d" % fam, "B%02d" % fam, n_kids,
                                    father_parents=fp,
                                    pheno_parents=1 if fam % 2 else 2,
                                    pheno_kids=2 if fam % 2 else 1)
    k = 0
    while ped.n_emitted() < target:
        k += 1
        ped.add("BSNG%03d" % k, "BSG%03d" % k, sex=1 + (k % 2))
    return Spec("bigish", ped, AUTOSOMES, 50000,
                notes=("200 individuals: nuclear families, three-generation "
                       "units, undeclared cross-family sibships and unrelated "
                       "singletons. Performance and accumulation-order case."))


BUILDERS = {
    "trio": build_trio,
    "nuclear": build_nuclear,
    "threegen": build_threegen,
    "multifam": build_multifam,
    "dups": build_dups,
    "missing": build_missing,
    "monomorphic": build_monomorphic,
    "sexchr": build_sexchr,
    "unrelated": build_unrelated,
    "admixed": build_admixed,
    "singleton": build_singleton,
    "pair": build_pair,
    "bigish": build_bigish,
}

DATASET_ORDER = list(BUILDERS)

# Above this many samples the manifest lists only pairs with nonzero expected
# kinship; below it, every pair is listed (unrelated ones included).
ALL_PAIRS_MAX_SAMPLES = 40


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pair_records(spec):
    people = spec.ped.people
    phi = kinship_matrix(people)
    emitted = [p.idx for p in people if p.emit]
    all_pairs = len(emitted) <= ALL_PAIRS_MAX_SAMPLES
    recs = []
    for x in range(len(emitted)):
        for y in range(x + 1, len(emitted)):
            i, j = emitted[x], emitted[y]
            k = phi[i][j]
            rel = classify(people, i, j, phi)
            if not all_pairs and k <= 0.0:
                continue
            recs.append([people[i].key, people[j].key, rel, round(k, 10)])
    return recs, all_pairs


def main(argv=None):
    ap = argparse.ArgumentParser(
        description="Generate the open-king differential-test corpus.")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--seed", type=int, default=20260813)
    ap.add_argument("--only", nargs="+", metavar="NAME",
                    help="generate only these datasets (default: all)")
    ap.add_argument("--emit-freq", action="store_true",
                    help="also write <name>.expected_freq.tsv (verification aid)")
    args = ap.parse_args(argv)

    names = DATASET_ORDER
    if args.only:
        unknown = [n for n in args.only if n not in BUILDERS]
        if unknown:
            ap.error("unknown dataset(s) %s; valid: %s"
                     % (", ".join(unknown), ", ".join(DATASET_ORDER)))
        names = [n for n in DATASET_ORDER if n in set(args.only)]

    outdir = os.path.abspath(args.outdir)
    os.makedirs(outdir, exist_ok=True)
    manifest_path = os.path.join(outdir, "MANIFEST.json")

    manifest = {"generator_version": GENERATOR_VERSION,
                "master_seed": args.seed,
                "generator": "tests/parity/generate_corpus.py",
                "bed_format": "PLINK 1 SNP-major (6c 1b 01); "
                              "00=hom A1, 01=missing, 10=het, 11=hom A2; "
                              "first sample in the low bits",
                "genetic_map": "bim column 3 in centimorgans at 1 cM per Mb",
                "datasets": {}}
    if os.path.exists(manifest_path):
        try:
            with open(manifest_path) as fh:
                old = json.load(fh)
            if old.get("master_seed") == args.seed:
                manifest["datasets"] = old.get("datasets", {})
        except (ValueError, OSError):
            pass

    rows = []
    for name in names:
        spec = BUILDERS[name]()
        seed = dataset_seed(args.seed, name)
        snps, n_snps, n_emit = simulate(spec, seed, outdir,
                                        emit_freq=args.emit_freq)
        chroms = sorted({c for c, _bp, _cm in snps})
        recs, all_pairs = pair_records(spec)
        base = os.path.join(outdir, name)
        manifest["datasets"][name] = {
            "seed": seed,
            "n_samples": n_emit,
            "n_snps": n_snps,
            "chromosomes": chroms,
            "snps_per_chromosome": {str(c): sum(1 for x, _b, _m in snps if x == c)
                                    for c in chroms},
            "notes": spec.notes,
            "pairs_listed": "all" if all_pairs else "related_only",
            "pairs": recs,
            "samples": [{"fid": p.fid, "iid": p.iid, "sex": p.sex,
                         "father": (spec.ped.people[p.fidx].iid
                                    if p.fidx is not None
                                    and spec.ped.people[p.fidx].emit else "0"),
                         "mother": (spec.ped.people[p.midx].iid
                                    if p.midx is not None
                                    and spec.ped.people[p.midx].emit else "0")}
                        for p in spec.ped.emitted()],
            "sha256": {"%s.bed" % name: sha256_file(base + ".bed"),
                       "%s.bim" % name: sha256_file(base + ".bim"),
                       "%s.fam" % name: sha256_file(base + ".fam")},
        }
        n_rel = sum(1 for r in recs if r[3] > 0)
        rows.append((name, n_emit, n_snps, len(chroms), n_rel))
        print("  %-12s %4d samples  %6d SNPs  %2d chrom  %5d related pairs"
              % (name, n_emit, n_snps, len(chroms), n_rel))

    with open(manifest_path, "w") as fh:
        json.dump(manifest, fh, indent=2, sort_keys=True)
        fh.write("\n")
    print("wrote %s (%d datasets)" % (manifest_path, len(manifest["datasets"])))
    return 0


if __name__ == "__main__":
    sys.exit(main())
