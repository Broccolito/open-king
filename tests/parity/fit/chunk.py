"""The chunk-scan IBD2 caller for `--ibs`, and its scorecard against the corpus.

The rule was measured on constructed filesets in `docs/research/fixtures/segfit.py`;
`docs/research/16-segment-extension.md` derives it. In one sentence: an IBD2 run is
*confirmed in chunks of five het-vs-hom mismatches, each of which must carry at least 95
HetHet markers over at least three words*, and the reported segment stops at the last
confirmed chunk.

**It is now the committed rule.** It lives in `engine.py::ibd2_words` (mirroring
`Scan::ibd2_words` in `crates/king-core/src/ibdseg.rs`) and this file only re-exports it
and keeps the superseded rule alive, so the before/after table of §8.1 still reproduces.

    python3 chunk.py            # scorecard: MaxIBD2 and Pr_IBD2 against the reference
    python3 chunk.py -v         # ...plus every miss
"""

import sys

import engine as E
import kingdata as kd

WORD = 64

DIRTY = E.IBS_IBD2_DIRTY        # het-vs-hom mismatches that break an --ibs IBD2 run
CHUNK_MIS = E.IBS_IBD2_CHUNK_MIS  # mismatches that close a chunk
CHUNK_HET = E.IBS_IBD2_HETHET   # HetHet a chunk must carry to be confirmed
MIN_WORDS = E.IBS_IBD2_MIN_WORDS  # words the measured interval must span
MIN_CHUNK = E.IBS_IBD2_CHUNK_WORDS  # words a chunk must span to be confirmable
EXT_MIS = E.IBS_IBD2_EXT_MIS    # mismatches the interval may reach past the confirmation

# The committed caller, under the name the write-up uses.
ibd2_words_chunk = E.ibd2_words


def ibd2_words_prechunk(sc):
    """The rule the chunk scan replaced, kept only so the scorecard below can show what
    the change bought: accept-or-refuse the whole run, on one HetHet count over the
    measured interval. Scores MaxIBD2 148/158 and Pr_IBD2 100/158.
    """
    n, w0, w1 = sc.n, sc.w0, sc.w1
    if n == 0:
        return []
    clean = [int(sc.n1[w0 + k]) < DIRTY for k in range(n)]
    ok = list(clean)
    for k in range(1, max(0, n - 1)):
        if not clean[k] and clean[k - 1] and clean[k + 1]:
            ok[k] = True
    out = []
    k = 0
    while k < n:
        if not ok[k]:
            k += 1
            continue
        k0 = k
        while k < n and ok[k]:
            k += 1
        u, v = w0 + k0, w0 + k - 1
        hi = w1 if v + 2 >= w1 else v + 1
        lo = u
        if out:
            lo = max(lo, out[-1][1] + 1)
        if lo > hi or hi + 1 - lo < MIN_WORDS:
            continue
        if v + 1 < w1 and int(sc.nhh[lo:hi + 1].sum()) < CHUNK_HET:
            continue
        out.append((lo, hi))
    return out


# ---------------------------------------------------------------------------
# scorecard
# ---------------------------------------------------------------------------

def _calls(fn, ds, i, j, **kw):
    best = tot = 0.0
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg)
        if sc.n == 0:
            continue
        for u, e in fn(sc, **kw):
            ln = float(ds.pos[WORD * e + 63] - ds.pos[WORD * u])
            tot += ln
            best = max(best, ln)
    return best, tot


def score(fn, verbose=False, **kw):
    okmax = nmax = okpr = npr = 0
    bias = 0.0
    misses = []
    for name, i, j, want in E.max_targets():
        ds = kd.load(name)
        best, _ = _calls(fn, ds, i, j, **kw)
        nmax += 1
        if abs(best - float(want)) < 0.5:
            okmax += 1
        elif verbose:
            misses.append((name, ds.fam[i][1], ds.fam[j][1], float(want), best))
    for name, i, j, want in E.pr_targets():
        ds = kd.load(name)
        best, tot = _calls(fn, ds, i, j, **kw)
        pr = 0.0 if best < E.LONG else tot / ds.denom
        npr += 1
        if "%.4f" % pr == want:
            okpr += 1
        bias += pr - float(want)
    return okmax, nmax, okpr, npr, bias / max(1, npr), misses


if __name__ == "__main__":
    verbose = "-v" in sys.argv
    rows = [("pre-chunk  ", ibd2_words_prechunk, {}),
            ("committed  ", ibd2_words_chunk, {}),
            ("chunk after", ibd2_words_chunk, dict(restart="after")),
            ("chunk at   ", ibd2_words_chunk, dict(restart="at"))]
    print(f"{'rule':<12} {'MaxIBD2':>12} {'Pr_IBD2':>12} {'Pr bias':>10}")
    for label, fn, kw in rows:
        okm, nm, okp, npn, bias, misses = score(fn, verbose, **kw)
        print(f"{label:<12} {okm:>6}/{nm:<5} {okp:>6}/{npn:<5} {bias:>+10.4f}")
        if verbose and misses:
            for m in misses[:25]:
                print("      miss %-11s %-8s %-8s ref %12.0f  got %12.0f" % m)
