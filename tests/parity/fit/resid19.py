"""Characterise the residual of the `.seg` IBD2 caller, row by row.

`seg18.py` says `IBD2Seg` is exact on 896 of 982 corpus rows at the default floor and on
880/877 at 5/10 Mb — a systematic residual, not a floor effect.  This script does not
propose a rule: it profiles the 86 wrong rows so the next clause can be chosen rather than
guessed.

Three views, all printed by `python3 resid19.py`:

* **the row profile** — sign and size of the disagreement per row, in base pairs and in
  units of the printed ulp (`D / 10000`) and of the dataset's median marker gap, broken
  down by dataset, by `InfType`, by call count, and by whether the row's `IBD1Seg` is
  nonetheless exact;
* **the per-call diff** — for the nine datasets `seglen_probe.py` has inverted, the
  reference's exact multiset of called IBD2 segment lengths under 10 Mb against ours, so a
  disagreement is localised to an individual call rather than a total;
* **the wrong rows themselves**, one line each, with the calls we make.

    python3 resid19.py            # the profile
    python3 resid19.py rows       # ...plus every wrong row's calls
    python3 resid19.py calls      # ...plus the per-call diff, pair by pair
"""

import json
import os
import sys
from collections import Counter, defaultdict

import engine as E
import kingdata as kd
import seg19 as S19

WORD = E.WORD

#: The IBD2 caller under test. `None` = `engine.SegScan.ibd2` (the committed rule);
#: an `S19.R19` = the candidate of `19-ibd2seg-residual.md`.
RULE = S19.R19()
SEGLEN = os.path.join(kd.ROOT, "tests", "parity", "work", "seglen")

FLOORS = [(3_000_000, "__ibdseg"), (5_000_000, "__ibdseg_seglength5"),
          (10_000_000, "__ibdseg_seglength10")]


def calls_of(ds, i, j, p=E.BASE, min_bp=E.SEGLEN):
    """Our IBD2 calls for one pair, as [(seg, lo, hi)] over marker indices."""
    out = []
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        c = (sc.ibd2(ds.pos, min_bp) if RULE is None
             else S19.ibd2_19(sc, ds, i, j, RULE, ds.pos, min_bp))
        for lo, hi in c:
            out.append((seg, lo, hi))
    return out


def ibd1_calls_of(ds, i, j, p=E.BASE, min_bp=E.SEGLEN):
    out = []
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        for lo, hi in sc.ibd1(ds.pos, min_bp):
            out.append((seg, lo, hi))
    return out


def row(ds, i, j, min_bp, ref):
    """One graded row: our numbers, the reference's, and the disagreement."""
    pos, d = ds.pos, ds.denom
    c2 = calls_of(ds, i, j, min_bp=min_bp)
    c1 = ibd1_calls_of(ds, i, j, min_bp=min_bp)
    ibd2 = sum(int(pos[b] - pos[a]) for _s, a, b in c2)
    plain = [(a, b) for _s, a, b in c2]
    ibd1 = 0
    for _s, a, b in c1:
        ibd1 += sum(v for v in (int(pos[y] - pos[x])
                                for x, y in E._pieces((a, b), plain))
                    if v >= min_bp)
    a1, a2, _ap, _at = ref[(i, j)]
    return dict(ds=ds.name, i=i, j=j, calls=c2, ibd1_calls=c1,
                ibd2=ibd2, ibd1=ibd1,
                g2=ibd2 / d, g1=ibd1 / d, r2=a2, r1=a1,
                ok2=kd.fmt4(ibd2 / d) == a2, ok1=kd.fmt4(ibd1 / d) == a1,
                dbp=ibd2 - a2 * d, inftype=ref[(i, j)][3])


def profile(min_bp, suffix, verbose=False):
    rows = []
    for name in kd.DATASETS:
        ds = kd.load(name)
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            if (i, j) not in ref:
                continue
            rows.append(row(ds, i, j, min_bp, ref))
    bad = [r for r in rows if not r["ok2"]]
    print("=== --seglength %d Mb: %d rows, %d wrong IBD2Seg, %d wrong IBD1Seg"
          % (min_bp // 1_000_000, len(rows), len(bad),
             sum(1 for r in rows if not r["ok1"])))
    over = [r for r in bad if r["dbp"] > 0]
    under = [r for r in bad if r["dbp"] < 0]
    print("  over-call %d   under-call %d" % (len(over), len(under)))
    print("  wrong rows whose IBD1Seg is nonetheless exact: %d / %d"
          % (sum(1 for r in bad if r["ok1"]), len(bad)))
    print("  rows wrong on IBD1 but right on IBD2: %d"
          % sum(1 for r in rows if r["ok2"] and not r["ok1"]))

    by_ds = Counter(r["ds"] for r in bad)
    tot_ds = Counter(r["ds"] for r in rows)
    print("  by dataset: " + "  ".join("%s %d/%d" % (k, by_ds.get(k, 0), tot_ds[k])
                                       for k in kd.DATASETS))
    by_t = Counter(r["inftype"] for r in bad)
    tot_t = Counter(r["inftype"] for r in rows)
    print("  by InfType: " + "  ".join("%s %d/%d" % (k, by_t.get(k, 0), tot_t[k])
                                       for k in sorted(tot_t)))

    # size of the error, in ulps of the printed column and in median marker gaps
    print("  |delta| in ulps of the printed column (1 ulp = D/10000):")
    hist = Counter()
    for r in bad:
        d = kd.load(r["ds"]).denom
        hist[min(int(abs(r["dbp"]) / (d / 10000.0)), 20)] += 1
    print("     " + "  ".join("%s:%d" % ("20+" if k == 20 else k, v)
                              for k, v in sorted(hist.items())))
    print("  |delta| in median marker gaps of the dataset:")
    hist = Counter()
    for r in bad:
        ds = kd.load(r["ds"])
        gap = _median_gap(ds)
        hist[min(int(round(abs(r["dbp"]) / gap)), 200)] += 1
    print("     " + "  ".join("%s:%d" % (k, v) for k, v in sorted(hist.items())[:25]))

    # call-count / geometry profile
    print("  calls per pair (wrong rows):   " + _dist(len(r["calls"]) for r in bad))
    print("  calls per pair (right rows):   "
          + _dist(len(r["calls"]) for r in rows if r["ok2"] and r["r2"] > 0))
    print("  usable segments touched (wrong): "
          + _dist(len({id(s) for s, _a, _b in r["calls"]}) for r in bad))
    longest = []
    for r in bad:
        ds = kd.load(r["ds"])
        if r["calls"]:
            longest.append(max(int(ds.pos[b] - ds.pos[a]) for _s, a, b in r["calls"]))
    if longest:
        print("  longest call on a wrong row: median %.1f Mb, min %.1f, max %.1f"
              % (sorted(longest)[len(longest) // 2] / 1e6, min(longest) / 1e6,
                 max(longest) / 1e6))
    print("  reference IBD2Seg on wrong rows: " + _dist(round(r["r2"], 1) for r in bad))
    if verbose:
        print("  wrong rows:")
        for r in sorted(bad, key=lambda r: -abs(r["dbp"])):
            ds = kd.load(r["ds"])
            print("    %-12s %-10s %-10s ref %.4f got %.4f  d=%+9d bp (%+6.1f ulp)  "
                  "calls %d  ibd1 %s"
                  % (r["ds"], ds.fam[r["i"]][1], ds.fam[r["j"]][1], r["r2"], r["g2"],
                     r["dbp"], r["dbp"] / (ds.denom / 10000.0), len(r["calls"]),
                     "ok" if r["ok1"] else "WRONG"))
    return rows, bad


_GAP = {}


def _median_gap(ds):
    if ds.name not in _GAP:
        gaps = []
        for _c, lo, hi in ds.segs:
            gaps.extend((ds.pos[lo + 1:hi + 1] - ds.pos[lo:hi]).tolist())
        gaps.sort()
        _GAP[ds.name] = gaps[len(gaps) // 2]
    return _GAP[ds.name]


def _dist(it):
    c = Counter(it)
    return "  ".join("%s:%d" % (k, v) for k, v in sorted(c.items()))


# ---------------------------------------------------------------------------
# the per-call diff, against `seglen_probe.py`'s inverted segment lengths
# ---------------------------------------------------------------------------

def percall(verbose=False):
    """Reference call lengths (< 10 Mb, exact) against ours, pair by pair."""
    print("=== per-call diff against the --seglength inversion (calls < 10 Mb)")
    tot = Counter()
    detail = []
    for name in kd.DATASETS:
        path = os.path.join(SEGLEN, "%s.IBD2Seg.json" % name)
        if not os.path.exists(path):
            tot["no probe: " + name] += 1
            continue
        want = json.load(open(path))
        ds = kd.load(name)
        for key, lens in sorted(want.items()):
            i, j = map(int, key.split(","))
            got = [int(ds.pos[b] - ds.pos[a])
                   for _s, a, b in calls_of(ds, i, j)]
            gshort = sorted(v for v in got if v < 10_000_000)
            wshort = sorted(lens)
            # the probe reads one base pair short (seglen_invert.TOL)
            matched, gleft, wleft = _match(gshort, wshort, 2)
            tot["pairs"] += 1
            tot["calls_ref"] += len(wshort)
            tot["calls_ours"] += len(gshort)
            tot["matched"] += matched
            tot["ours_only"] += len(gleft)
            tot["ref_only"] += len(wleft)
            if gleft or wleft:
                tot["pairs_disagreeing"] += 1
                detail.append((name, ds.fam[i][1], ds.fam[j][1], gleft, wleft))
    for k in ("pairs", "pairs_disagreeing", "calls_ref", "calls_ours", "matched",
              "ours_only", "ref_only"):
        print("  %-18s %d" % (k, tot[k]))
    for k, v in tot.items():
        if k.startswith("no probe"):
            print("  %s" % k)
    if verbose or detail:
        for name, a, b, gleft, wleft in detail:
            print("    %-12s %-10s %-10s ours-only %s  ref-only %s"
                  % (name, a, b, [round(v / 1e6, 3) for v in gleft],
                     [round(v / 1e6, 3) for v in wleft]))


def _match(got, want, tol):
    """Greedy nearest match of two sorted length lists within `tol` base pairs."""
    g = list(got)
    w = list(want)
    n = 0
    for v in list(w):
        hit = next((x for x in g if abs(x - v) <= tol), None)
        if hit is not None:
            g.remove(hit)
            w.remove(v)
            n += 1
    return n, g, w


# ---------------------------------------------------------------------------
# structural fingerprints of the wrong rows
# ---------------------------------------------------------------------------

def fingerprint(bad):
    """What is structurally special about the calls on a wrong row.

    Every count is over the *calls*, not the rows: the question is which clause of the
    rule each wrong row's calls actually exercised.
    """
    print("=== clause exposure, wrong rows against right ones")
    tags = defaultdict(Counter)
    for tag, rows in (("wrong", bad[0]), ("right", bad[1])):
        for r in rows:
            ds = kd.load(r["ds"])
            i, j = r["i"], r["j"]
            for seg in ds.segs:
                sc = E.SegScan(ds, i, j, seg, E.BASE)
                if sc.n == 0:
                    continue
                for k, v in _seg_tags(sc, ds.pos).items():
                    tags[tag][k] += v
            tags[tag]["rows"] += 1
    keys = sorted(set(tags["wrong"]) | set(tags["right"]))
    print("  %-26s %10s %10s %10s %10s" % ("", "wrong", "/row", "right", "/row"))
    for k in keys:
        w, r = tags["wrong"][k], tags["right"][k]
        print("  %-26s %10d %10.3f %10d %10.3f"
              % (k, w, w / max(tags["wrong"]["rows"], 1),
                 r, r / max(tags["right"]["rows"], 1)))


def _seg_tags(sc, pos):
    """Per usable segment: how many times each clause of the rule fired."""
    p = sc.p
    n = sc.n
    w0, w1 = sc.w0, sc.w1
    cum = sc.cum2s
    z = [int(sc.n0[w0 + k]) != 0 for k in range(n)]
    mis = [int(sc.n1[w0 + k]) for k in range(n)]
    usable = [(not z[k]) and mis[k] < p.ibd2_dirty_ibs1 for k in range(n)]
    out = Counter()

    def ge_of(b):
        return b + 1 if (b + 1 < n and not z[b + 1] and mis[b + 1]) else b

    def gate_ok(g, b):
        return int(cum[w0 + ge_of(b) + 1] - cum[w0 + g]) >= p.gate

    ok = list(usable)
    gs0 = None
    for k in range(n):
        if usable[k]:
            if gs0 is None and mis[k] == 0:
                gs0 = k
            continue
        bridged = False
        if (gs0 is not None and k > 0 and not z[k] and k + 1 < n
                and usable[k + 1] and mis[k + 1] == 0):
            b2 = k + 1
            while b2 + 1 < n and usable[b2 + 1]:
                b2 += 1
            bridged = gate_ok(gs0, k - 1) and gate_ok(k + 1, b2)
            out["bridge_considered"] += 1
        if bridged:
            ok[k] = True
            out["bridge_taken"] += 1
        else:
            gs0 = None

    emitted = 0
    for a, b in E._runs(ok):
        u, v = w0 + a, w0 + b
        gs = next((t for t in range(a, b + 1) if mis[t] == 0), None)
        if gs is None:
            out["run_no_gatestart"] += 1
            continue
        if not gate_ok(gs, b):
            out["run_gate_refused"] += 1
            continue
        out["run_emitted"] += 1
        if emitted:
            out["push_applied"] += 1
            if gs != a:
                out["push_from_late_gs"] += 1
        if a > 0 and not z[a - 1] and int(sc.ibs1[u - 1]):
            out["reach_left"] += 1
        if b + 1 < n and not z[b + 1] and int(sc.ibs1[v + 1]):
            out["reach_right"] += 1
        if u == w0:
            out["fringe_left"] += 1
        if v == w1:
            out["fringe_right"] += 1
        if u == w0 and v == w1:
            out["fringe_both"] += 1
        if any(ok[t] and not usable[t] for t in range(a, b + 1)):
            out["run_has_bridge"] += 1
            if a == 0 or b == n - 1:
                out["bridge_at_segment_edge"] += 1
        emitted += 1
    if emitted > 1:
        out["multi_call_segment"] += 1
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    all_rows = None
    for bp, sfx in FLOORS:
        rows, bad = profile(bp, sfx, verbose=(mode == "rows"))
        if bp == 3_000_000:
            all_rows = (bad, [r for r in rows if r["ok2"] and r["r2"] > 0])
        print()
    percall(verbose=(mode == "calls"))
    print()
    fingerprint(all_rows)
