"""Re-measure the IBD2 right-hand endpoint on the solid-IBS0 canvas.

`10-segment-rule-fixtures.md` §3 states that a run of `W` clean IBD2 words reports exactly
`64W - 1` marker intervals — the run's own words and nothing more — while an IBD1 run of
`W` words reports `64(W + 1) - 1`. The committed engine ends an IBD2 call on the flanking
word's **last** IBS0, which on this canvas (every marker outside the block an opposite
homozygote) is that word's bit 63 and so reports `64(W + 1) - 1`, one whole word too many.
This re-runs the measurement rather than trusting the write-up, and prints what each
candidate endpoint rule predicts.

    python3 ibd2end.py
"""

import rig2

WIDTHS = (1, 2, 3, 4, 6, 8)


def main():
    rig = rig2.Rig(spacing=100_000, n1=640, n2=640)
    print("W   reported IBD2 markers    64W-1   64(W+1)-1")
    for w in WIDTHS:
        f = rig.new("ibd2end_w%d" % w)
        # a block of `w` complete words of the global grid, IBD2, on the solid canvas
        rig.block(f, 64 * 2, 64 * (2 + w), rig2.L.IBD2)
        r = rig.read(f)
        got = None if r is None else r["test2_mk"]
        print("%-3d %-24s %-7d %d" % (w, got, 64 * w - 1, 64 * (w + 1) - 1))


if __name__ == "__main__":
    main()
