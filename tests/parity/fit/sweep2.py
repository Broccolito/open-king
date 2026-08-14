"""Score a set of `engine.Params` variants against all four graders at once.

The four are deliberately different instruments over the same calls:

* `.seg` exact rows — the parity target, but an aggregate of two rulers and a filter;
* `MaxIBD2` — one exact segment length per pair, word-aligned, no filter;
* `Pr_IBD2` — the word-aligned **total**, no filter, so it grades the whole IBD2 call set
  rather than its longest member (and it is where the current engine is worst: 7/158);
* the `--seglength 5` and `10` captures, held out from every fit.

    python3 sweep2.py                       # the named variants below
    python3 sweep2.py -p ibd2_ext=0 ...     # one ad-hoc variant
"""

import sys
from dataclasses import replace

import engine as E


def line(tag, p):
    s = E.score_seg(p)
    ok, n, _ = E.score_max(p)
    pok, pn, perr = E.score_pr(p)
    print("%-34s seg %3d/%3d (both %3d) x%-3d m%-3d mae %.5f | Max %3d/%d | Pr %3d/%d "
          "bias %+.4f"
          % (tag, s["exact"], s["rows"], s["both"], s["extra"], s["missing"], s["mae"],
             ok, n, pok, pn, perr))
    return s


VARIANTS = [
    ("baseline", {}),
    ("no bridge", dict(bridge=False)),
    ("ext=0", dict(ibd2_ext=0)),
    ("tail=0", dict(ibd2_tail=0)),
    ("tail=1", dict(ibd2_tail=1)),
    ("tail=3", dict(ibd2_tail=3)),
    ("dirty>=3", dict(ibd2_dirty_ibs1=3)),
    ("dirty>=4", dict(ibd2_dirty_ibs1=4)),
    ("dirty>=6", dict(ibd2_dirty_ibs1=6)),
    ("dirty>=7", dict(ibd2_dirty_ibs1=7)),
    ("dirty>=8", dict(ibd2_dirty_ibs1=8)),
    ("dirty>=10", dict(ibd2_dirty_ibs1=10)),
    ("min_run2=2", dict(min_run2=2)),
    ("min_run2=3", dict(min_run2=3)),
    ("start refined", dict(ibd2_start_refine=True)),
    ("clip after length test", dict(clip_before_len=False)),
]


def main():
    argv = sys.argv[1:]
    kw = {}
    while "-p" in argv:
        k = argv.index("-p")
        key, _, val = argv[k + 1].partition("=")
        kw[key] = eval(val)  # noqa: S307
        del argv[k:k + 2]
    if kw:
        p = replace(E.BASE, **kw)
        line(p.label()[:34], p)
        return
    for tag, kw in VARIANTS:
        line(tag, replace(E.BASE, **kw) if kw else E.BASE)


if __name__ == "__main__":
    main()
