"""Score the `.seg` IBD1 rule of `docs/research/18-ibd1-caller.md` on the corpus.

`seg17.py` graded the IBD2 column; this grades the **IBD1** one. It reuses `engine.py`
for the IBD1 caller itself (whose geometry `18-…` §§1-5 re-measured and confirmed
unchanged) and `seg17.py` for the IBD2 caller, and varies only what `18-…` §6 measured:
how an IBD1 call and the IBD2 calls inside it combine into `IBD1Seg`.

    python3 seg18.py            # the scorecard, at 3 / 5 / 10 Mb
    python3 seg18.py grid       # every knob of §6 swept against the corpus
"""

import sys
from dataclasses import dataclass, replace

import engine as E
import kingdata as kd
import seg17 as S17

WORD = E.WORD
IBD2 = S17.R17()


@dataclass(frozen=True)
class R18:
    """How `IBD1Seg` is assembled from an IBD1 call and the IBD2 calls inside it."""

    # `"exclusive"` — the measured rule — drops the IBD2 call's own end markers, so an
    # IBD1 call `[lo, hi]` cut by `[a, b]` leaves `[lo, a-1]` and `[b+1, hi]`.
    # `"inclusive"` is the retired "length minus overlap".
    cut: str = "exclusive"
    # `--seglength` floor on each surviving piece (`"drop"`), on their total
    # (`"whole"`), or not at all (`"keep"` — the retired behaviour).
    frag: str = "drop"
    # subtract the IBD2 calls that survived `--seglength` (None) or all of them (0).
    sub_bp: int | None = None


RETIRED = R18(cut="inclusive", frag="keep")


def pieces(call, others, cut="exclusive"):
    """One IBD1 call with the IBD2 calls cut out of it (`18-…` §6.1)."""
    d = 1 if cut == "exclusive" else 0
    lo, hi = call
    out, cur = [], lo
    for a, b in sorted(others):
        if b < lo or a > hi:
            continue
        if a - d > cur:
            out.append((cur, a - d))
        cur = max(cur, b + d)
    if cur + (0 if d else 1) <= hi:
        out.append((cur, hi))
    return out


def call_pair(ds, i, j, p, min_bp=E.SEGLEN):
    pos = ds.pos
    ibd1 = ibd2 = longest = 0
    for seg in ds.segs:
        sc = E.SegScan(ds, i, j, seg, E.BASE)
        if sc.n == 0:
            continue
        c2 = S17.ibd2_17(sc, ds, i, j, IBD2, pos, min_bp)
        c1 = sc.ibd1(pos, min_bp)
        sub = c2 if p.sub_bp is None else S17.ibd2_17(sc, ds, i, j, IBD2, pos, p.sub_bp)
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2 += ln
            longest = max(longest, ln)
        for lo, hi in c1:
            longest = max(longest, int(pos[hi] - pos[lo]))
            lens = [int(pos[b] - pos[a]) for a, b in pieces((lo, hi), sub, p.cut)]
            if p.frag == "drop":
                ibd1 += sum(v for v in lens if v >= min_bp)
            elif p.frag == "whole":
                ibd1 += sum(lens) if sum(lens) >= min_bp else 0
            else:
                ibd1 += sum(lens)
    return ibd1, ibd2, longest


def score(p, min_bp=E.SEGLEN, suffix="__ibdseg"):
    rows = exact = i1 = i2 = extra = missing = 0
    err = worst = 0.0
    for name in kd.DATASETS:
        ds = kd.load(name)
        d = ds.denom
        ref = ds._read_seg(suffix)
        for i, j in ds.pairs():
            a, b, lg = call_pair(ds, i, j, p, min_bp)
            got, want = lg >= E.LONG, (i, j) in ref
            if not want:
                extra += got
                continue
            if not got:
                missing += 1
                continue
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = a / d, b / d
            gp = g2 + g1 / 2
            rows += 1
            ok1, ok2 = kd.fmt4(g1) == a1, kd.fmt4(g2) == a2
            i1 += ok1
            i2 += ok2
            exact += (ok1 and ok2 and kd.fmt4(gp) == ap
                      and kd.inf_type(g1, g2, gp) == at)
            err += abs(gp - ap)
            worst = max(worst, abs(gp - ap))
    return dict(rows=rows, exact=exact, ibd1=i1, ibd2=i2, extra=extra, missing=missing,
                mae=err / max(rows, 1), worst=worst)


def show(tag, s):
    print("%-34s exact %4d  ibd1 %4d  ibd2 %4d  extra %3d  miss %3d  MAE %.6f  "
          "worst %.4f"
          % (tag, s["exact"], s["ibd1"], s["ibd2"], s["extra"], s["missing"],
             s["mae"], s["worst"]))


FLOORS = [(3_000_000, "__ibdseg"), (5_000_000, "__ibdseg_seglength5"),
          (10_000_000, "__ibdseg_seglength10")]


if __name__ == "__main__":
    base = R18()
    for bp, sfx in FLOORS:
        print("--seglength %d Mb" % (bp // 1_000_000))
        show("  retired (length minus overlap)", score(RETIRED, bp, sfx))
        show("  18 measured", score(base, bp, sfx))
    if len(sys.argv) > 1 and sys.argv[1] == "grid":
        sys.stdout.reconfigure(line_buffering=True)
        print("\nknob by knob, at the default 3 Mb")
        for cut in ("exclusive", "inclusive"):
            for frag in ("drop", "whole", "keep"):
                show("  cut=%s frag=%s" % (cut, frag),
                     score(replace(base, cut=cut, frag=frag)))
        for sub in (None, 0):
            show("  sub_bp=%s" % sub, score(replace(base, sub_bp=sub)))
