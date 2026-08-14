"""Where exactly does an IBD2 call stop when the word that ends it carries one IBS0?

chr2 is made **entirely IBD2**, so the pair has no natural opposite homozygote anywhere on
it, and then a single IBS0 is forced at a chosen bit of a chosen word. That word is the
only thing that can break the IBD2 run, and the two pieces' total length says, to the
marker, which IBS0-relative position each endpoint took.

    python3 ibd2gap.py

Predictions printed alongside:

* `last`  — the committed rule: the left piece ends **on** the IBS0.
* `first-1` — the rule `10-segment-rule-fixtures.md` §3 implies: it stops one marker
  before the first IBS0 of that word, so on a word that is entirely IBS0 the call keeps
  none of it.
"""

import rig2

L = rig2.L


def main():
    sp = 100_000
    n2 = 640
    rig = rig2.Rig(spacing=sp, n1=640, n2=n2)
    print("poke word.bit   reported IBD2 markers   'last' pred   'first-1' pred")
    for word, bit in ((6, 20), (6, 0), (6, 63), (4, 32), (3, 5)):
        f = rig.new("ibd2gap_%d_%d" % (word, bit), solid=False)
        f.set_state(1, 0, n2, L.IBD2)
        rig.poke(f, word * 64 + bit, "ibs0")
        r = rig.read(f)
        got = None if r is None else r["test2_mk"]
        # left piece [0 .. 64*word+bit] or [0 .. 64*word+bit-1]; right piece
        # [64*(word+1) .. n2-1] — the run after the break starts on its own word.
        left_last = 64 * word + bit
        right = (n2 - 1) - 64 * (word + 1)
        print("  w%-2d.%-3d      %-23s %-13d %d"
              % (word, bit, got, left_last + right, left_last - 1 + right))


if __name__ == "__main__":
    main()
