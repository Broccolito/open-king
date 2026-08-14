"""Invert `--ibs`'s `MaxIBD2` back to the word interval the reference must have called.

`MaxIBD2` is the sharpest instrument in the corpus: it is the length in base pairs of one
single IBD2 segment, measured word-aligned (`pos[64e+63] - pos[64u]`), so searching every
`(u, e)` inside every usable segment for that exact span recovers the reference's own
endpoints — no aggregate, no cancellation. 158 pairs carry a non-zero value.

    python3 invert.py            # localisation census + every case the engine misses
    python3 invert.py -a         # ...including the ones it already gets right
    python3 invert.py -p PARAM=V # score a variant of engine.Params first

Word cells print as `IBS0/IBS1(inf2)`; `|` marks a usable-segment edge, `[`/`]` the
interval, `<`/`>` the interval the engine produced.
"""

import sys

import kingdata as kd
import engine as E

WORD = 64


def candidates(ds, i, j, target, ext_max=2):
    """Word intervals `[u, e]` inside one usable segment whose word-aligned span is
    `target` and which an IBD2 caller could plausibly have produced.

    Plausible means: the interval opens on a word with no IBS0, and IBS0 is confined to
    its last `ext_max` words — an IBD2 run itself never contains an IBS0 word, so any
    IBS0 inside a called interval has to live in the words the call *extends* into.
    """
    _, n0, _n1, _c1, _c2 = E.masks(ds, i, j)
    pos = ds.pos
    out = []
    for _chrom, lo, hi in ds.segs:
        w0 = -(-lo // WORD)
        w1 = (hi + 1) // WORD - 1
        if w1 < w0:
            continue
        for u in range(w0, w1 + 1):
            if n0[u] != 0:
                continue
            base = int(pos[WORD * u])
            for e in range(u, w1 + 1):
                d = int(pos[WORD * e + 63]) - base
                if d > target:
                    break
                if d != target:
                    continue
                cut = max(u, e - ext_max + 1)
                if int(n0[u:cut].max(initial=0)) != 0:
                    continue
                out.append((u, e, w0, w1))
    return out


def profile(ds, i, j, a, b, marks):
    _, n0, n1, _c1, c2 = E.masks(ds, i, j)
    cells = []
    for w in range(max(0, a), min(len(n0), b + 1)):
        pre = marks.get(("pre", w), "")
        post = marks.get(("post", w), "")
        inf = int(c2[w + 1] - c2[w])
        cells.append("%s%d/%d(%d)%s" % (pre, int(n0[w]), int(n1[w]), inf, post))
    return " ".join(cells)


def report(p=E.BASE, show_all=False):
    tg = E.max_targets()
    census = {}
    misses = []
    for name, i, j, t in tg:
        ds = kd.load(name)
        got = E.max_ibd2(ds, i, j, p)
        cs = candidates(ds, i, j, t)
        census[len(cs)] = census.get(len(cs), 0) + 1
        if got != t or show_all:
            misses.append((name, i, j, t, got, cs))
    print("localisation: " + "  ".join("%d cand -> %d pairs" % (k, v)
                                       for k, v in sorted(census.items())))
    for name, i, j, t, got, cs in misses:
        ds = kd.load(name)
        u, e, ln = E.max_ibd2_words(ds, i, j, p)
        print("\n%-12s %s/%s   want %d  got %d  (%+d)"
              % (name, ds.fam[i][1], ds.fam[j][1], t, got, got - t))
        if u is not None:
            print("   engine  w%d..w%d  %d bp" % (u, e, ln))
        for (a, b, w0, w1) in cs:
            marks = {("pre", a): "[", ("post", b): "]", ("pre", w0): "|"}
            marks[("post", w1)] = marks.get(("post", w1), "") + "|"
            if u is not None:
                marks[("pre", u)] = marks.get(("pre", u), "") + "<"
                marks[("post", e)] = ">" + marks.get(("post", e), "")
            lo = max(w0 - 1, min(a, u if u is not None else a) - 3)
            hi = min(w1 + 1, max(b, e if e is not None else b) + 3)
            print("   target  w%d..w%d %s%s   %s"
                  % (a, b, " START-EDGE" if a == w0 else "",
                     " END-EDGE" if b == w1 else "",
                     profile(ds, i, j, lo, hi, marks)))
        if not cs:
            print("   target  NOT LOCALISED")


def main():
    p = E.BASE
    argv = sys.argv[1:]
    kw = {}
    while "-p" in argv:
        k = argv.index("-p")
        key, _, val = argv[k + 1].partition("=")
        kw[key] = eval(val)  # noqa: S307 - developer tool, values are literals
        del argv[k:k + 2]
    if kw:
        from dataclasses import replace
        p = replace(E.BASE, **kw)
        print("params:", p.label())
        ok, n, _ = E.score_max(p)
        print("MaxIBD2 exact %d/%d" % (ok, n))
    report(p, show_all="-a" in argv)


if __name__ == "__main__":
    main()
