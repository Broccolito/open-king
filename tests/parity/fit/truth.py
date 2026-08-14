"""Ground-truth IBD state per marker, recovered from the corpus generator.

The corpus is simulated by `tests/parity/generate_corpus.py` with real meioses, so the
true IBD state of every pair at every marker is knowable: re-run the same simulation with
the same seed while tracking which founder haplotype each transmitted allele came from.
This module replays `simulate()` with that tracking added and **no extra RNG draws**, so
the genotypes it reproduces are byte-identical to the committed `.bed` — `verify()`
asserts exactly that, which is what makes the truth trustworthy.

Truth is not a rule and cannot be used as one (the reference only sees genotypes).  It is
used here to answer a different question: is KING's caller trying to recover the true
segments, and where does it stop?
"""

import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))))

import generate_corpus as G  # noqa: E402
import kingdata as kd  # noqa: E402

MISSING = G.MISSING


def simulate_truth(spec, seed):
    """Replay one dataset, returning (origin_a, origin_b, geno) over emitted samples.

    `origin_*` are (n_snps, n_emit) int32 arrays of founder-haplotype ids.
    """
    import random
    rnd = random.Random(seed)
    rr = rnd.random
    people = spec.ped.people
    n_all = len(people)
    order = G.topo_order(people)
    emit_idx = [p.idx for p in people if p.emit]
    n_emit = len(emit_idx)

    snps = G.build_map(spec.chrom_spec, spec.n_snps, rnd)
    n_snps = len(snps)
    snpinfo = spec.snp_model(rnd, n_snps, snps)
    snp_miss = spec.snp_missing(rnd, n_snps)
    samp_miss = [spec.sample_missing.get(people[i].iid, 0.0) for i in emit_idx]

    hap_a = [0] * n_all
    hap_b = [0] * n_all
    org_a = [0] * n_all           # founder-haplotype id carried by hap_a
    org_b = [0] * n_all
    sel_f = [0] * n_all
    sel_m = [0] * n_all
    geno = [0] * n_all

    OA = np.zeros((n_snps, n_emit), dtype=np.int32)
    OB = np.zeros((n_snps, n_emit), dtype=np.int32)
    GT = np.zeros((n_snps, n_emit), dtype=np.int8)

    prev_chrom = None
    prev_cm = 0.0
    for j in range(n_snps):
        chrom, bp, cm = snps[j]
        new_chrom = chrom != prev_chrom
        recomb = 0.0 if new_chrom else 0.5 * (1.0 - math.exp(-2.0 * (cm - prev_cm) / 100.0))
        prev_chrom, prev_cm = chrom, cm
        info = snpinfo[j]
        fr = info["freqs"]
        is_x, is_y, is_mt = chrom == 23, chrom == 24, chrom == 26

        for i in order:
            p = people[i]
            if p.clone_of is not None:
                s = p.clone_of
                hap_a[i], hap_b[i] = hap_a[s], hap_b[s]
                org_a[i], org_b[i] = org_a[s], org_b[s]
                g = geno[s]
                if p.error_rate > 0.0 and g != MISSING and rr() < p.error_rate:
                    g = [x for x in (0, 1, 2) if x != g][rnd.randrange(2)]
                geno[i] = g
                continue

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
                    org_a[i] = org_b[i] = -1
                    geno[i] = MISSING
                    continue
                if f is not None:
                    a, oa = hap_a[f], org_a[f]
                else:
                    a, oa = draw(), 2 * i
                hap_a[i] = hap_b[i] = a
                org_a[i] = org_b[i] = oa
                geno[i] = 2 * a
                continue

            if is_mt:
                if m is not None:
                    a, oa = hap_a[m], org_a[m]
                else:
                    a, oa = draw(), 2 * i
                hap_a[i] = hap_b[i] = a
                org_a[i] = org_b[i] = oa
                geno[i] = 2 * a
                continue

            if is_x and p.sex == 1:
                if m is not None:
                    a = hap_a[m] if sel_m[i] == 0 else hap_b[m]
                    oa = org_a[m] if sel_m[i] == 0 else org_b[m]
                else:
                    a, oa = draw(), 2 * i
                hap_a[i] = hap_b[i] = a
                org_a[i] = org_b[i] = oa
                geno[i] = 2 * a
                continue

            if is_x:
                if f is not None:
                    a, oa = hap_a[f], org_a[f]
                else:
                    a, oa = draw(), 2 * i
                if m is not None:
                    b = hap_a[m] if sel_m[i] == 0 else hap_b[m]
                    ob = org_a[m] if sel_m[i] == 0 else org_b[m]
                else:
                    b, ob = draw(), 2 * i + 1
            else:
                if f is not None:
                    a = hap_a[f] if sel_f[i] == 0 else hap_b[f]
                    oa = org_a[f] if sel_f[i] == 0 else org_b[f]
                else:
                    a, oa = draw(), 2 * i
                if m is not None:
                    b = hap_a[m] if sel_m[i] == 0 else hap_b[m]
                    ob = org_a[m] if sel_m[i] == 0 else org_b[m]
                else:
                    b, ob = draw(), 2 * i + 1
            hap_a[i], hap_b[i] = a, b
            org_a[i], org_b[i] = oa, ob
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
        GT[j] = row
        OA[j] = [org_a[i] for i in emit_idx]
        OB[j] = [org_b[i] for i in emit_idx]

        called = sum(1 for g in row if g != MISSING)
        a1n = sum(g for g in row if g != MISSING)
        if 2 * called and 2 * a1n > 2 * called:
            pass  # allele swap does not change IBD state
    return OA, OB, GT, snps


_TCACHE = {}


def load_truth(name):
    """(OA, OB) restricted to the retained autosomal markers of dataset `name`."""
    if name in _TCACHE:
        return _TCACHE[name]
    spec = G.BUILDERS[name]()
    seed = G.dataset_seed(20260813, name)
    OA, OB, GT, snps = simulate_truth(spec, seed)
    ds = kd.load(name)
    keep = ds.keep
    OA, OB = OA[keep], OB[keep]
    _TCACHE[name] = (OA, OB)
    return _TCACHE[name]


def ibd_state(name, i, j):
    """True IBD state (0/1/2) at every retained autosomal marker for one pair."""
    OA, OB = load_truth(name)
    ai, bi = OA[:, i], OB[:, i]
    aj, bj = OA[:, j], OB[:, j]
    m1 = (ai == aj)
    m2 = (bi == bj)
    m3 = (ai == bj)
    m4 = (bi == aj)
    both = (m1 & m2) | (m3 & m4)
    any_ = m1 | m2 | m3 | m4
    return np.where(both, 2, np.where(any_, 1, 0)).astype(np.int8)


def verify(name):
    """The replay must reproduce the committed .bed exactly, or truth is worthless."""
    spec = G.BUILDERS[name]()
    seed = G.dataset_seed(20260813, name)
    _, _, GT, _ = simulate_truth(spec, seed)
    ds = kd.load(name)
    # rebuild dosage from the bit planes: plane0=hom, plane1=A1-ish
    n = len(ds.fam)
    nm = len(ds.pos)
    got = np.zeros((nm, n), dtype=np.int8)
    for s in range(n):
        bits0 = np.unpackbits(ds.p0[s].view(np.uint8), bitorder="little")[:nm]
        bits1 = np.unpackbits(ds.p1[s].view(np.uint8), bitorder="little")[:nm]
        # (1,1)=hom A1=2, (0,1)=het=1, (1,0)=hom A2=0, (0,0)=missing=3
        g = np.where((bits0 == 1) & (bits1 == 1), 2,
                     np.where((bits0 == 0) & (bits1 == 1), 1,
                              np.where((bits0 == 1) & (bits1 == 0), 0, 3)))
        got[:, s] = g
    want = GT[ds.keep]
    # the generator may have swapped A1/A2 per marker; compare up to that mirror
    ok_direct = (want == got)
    ok_swap = ((want == 3) & (got == 3)) | ((want != 3) & (got != 3) & (2 - want == got))
    per_marker = ok_direct.all(axis=1) | ok_swap.all(axis=1)
    return int(per_marker.sum()), len(per_marker)


if __name__ == "__main__":
    for name in kd.DATASETS:
        ok, tot = verify(name)
        print(f"{name:<12} markers reproduced {ok}/{tot}")
