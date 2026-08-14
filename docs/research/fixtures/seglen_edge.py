"""Is the `--seglength` floor inclusive? Bisected on a segment of known exact length.

`10-segment-rule-fixtures.md` §4 says the test is `length >= seglength`, measured with a
segment spanning exactly 3 000 000 bp. Probing the corpus with `--seglength` at base-pair
resolution says otherwise: open-king's own binary (whose filter is literally
`pos[hi] - pos[lo] >= min_bp`) drops a 7 901 886 bp call at `--seglength 7.901887`, while
the reference drops what looks like the same call one base pair earlier.

This settles it on a fixture whose segment length is an exact multiple of the spacing: a
one-word IBD1 block on the solid IBS0 canvas reports `64*2 - 1 = 127` marker intervals, so
at 50 kb spacing it is exactly 6 350 000 bp.

    python3 seglen_edge.py
"""

import rig2

L = rig2.L


def main():
    sp = 50_000
    rig = rig2.Rig(spacing=sp, n1=1280, n2=1280)
    f = rig.new("seglen_edge")
    rig.block(f, 64 * 3, 64 * 4, L.IBD1)          # one clean word
    base = rig.read(f)
    n = base["test1_mk"]
    exact = n * sp
    print("block reports %d marker intervals = %d bp" % (n, exact))
    for bp in (exact - 1, exact, exact + 1):
        r = rig.read(f, ("--seglength", "%.6f" % (bp / 1e6)))
        got = None if r is None else r["test1_mk"]
        print("  --seglength %d bp -> block %s"
              % (bp, "kept (%d)" % got if got == n else "dropped (%s)" % got))


if __name__ == "__main__":
    main()
