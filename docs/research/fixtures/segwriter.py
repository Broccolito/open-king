"""How `<prefix>.seg` prints `PropIBD`, and what order it prints its rows in.

Both questions are answered from the **captured reference output alone** — no genotypes,
no reference binary, no engine. That is what makes this rig different from every other one
in this directory, and it is the whole reason the two rules were findable at all:

* `PropIBD` turned out to be a function of the two columns printed beside it, so the
  hypothesis can be tested on any `.seg` file that exists;
* the row order is a permutation of a known set, so it can be tested the same way.

Both were invisible for as long as the *numbers* were wrong. Once `19-…` closed `IBD2Seg`
the residual was 176 rows that had both estimate columns exact and still printed a
different `PropIBD`, and that is a shape no segment-caller change can produce.

    python3 segwriter.py            # every table in docs/research/20-seg-writer.md
    python3 segwriter.py prop       # only the PropIBD rule
    python3 segwriter.py order      # only the row order

Reads `tests/parity/golden/`; needs the input corpus only for `.fam` sample order, which
`tests/parity/run_parity.py` regenerates on demand.
"""

import glob
import os
import sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
GOLDEN = os.path.join(ROOT, "tests", "parity", "golden")
DATA = os.path.join(ROOT, "tests", "parity", "work", "data")


# ---------------------------------------------------------------------------
# reading captures
# ---------------------------------------------------------------------------

def read_table(path, keys):
    """`{(id columns): (IBD1Seg, IBD2Seg, PropIBD)}` for any file carrying those three."""
    out = {}
    with open(path) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        ix = {n: i for i, n in enumerate(head)}
        if not {"IBD1Seg", "IBD2Seg", "PropIBD"} <= set(ix):
            return {}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            try:
                v = (f[ix["IBD1Seg"]], f[ix["IBD2Seg"]], f[ix["PropIBD"]])
            except IndexError:
                continue
            if "" in v:
                continue        # kingX.seg's header names more columns than it fills
            out[tuple(f[ix[c]] for c in keys)] = v
    return out


W4 = ["FID1", "ID1", "FID2", "ID2"]
W3 = ["FID", "ID1", "ID2"]


def every(name, keys=W4):
    """Every capture of `name`, as a list of `(case, table)`."""
    out = []
    for d in sorted(glob.glob(os.path.join(GOLDEN, "*", "*"))):
        p = os.path.join(d, name)
        if os.path.exists(p):
            t = read_table(p, keys)
            if t:
                out.append((os.path.relpath(d, GOLDEN), t))
    return out


def units(s):
    """The integer a `%.4lf` column shows, scaled by 10 000."""
    return int(round(float(s) * 10000))


# ---------------------------------------------------------------------------
# §1  PropIBD: the reference disagrees with itself
# ---------------------------------------------------------------------------

def section_self_disagreement():
    print("=== §1  the reference's two writers disagree, on its own output")
    tot = same = diff = 0
    cases = 0
    for d in sorted(glob.glob(os.path.join(GOLDEN, "*", "*"))):
        seg = read_table(os.path.join(d, "king.seg"), W4) \
            if os.path.exists(os.path.join(d, "king.seg")) else {}
        if not seg:
            continue
        kin = {}
        for fn, keys in (("king.kin", W3), ("king.kin0", W4)):
            p = os.path.join(d, fn)
            if os.path.exists(p):
                for k, v in read_table(p, keys).items():
                    kin[(k[1], k[2]) if len(k) == 3 else (k[1], k[3])] = v
        if not kin:
            continue
        cases += 1
        for k, v in seg.items():
            w = kin.get((k[1], k[3]))
            if w is None or w[0] != v[0] or w[1] != v[1]:
                continue
            tot += 1
            if v[2] == w[2]:
                same += 1
            else:
                diff += 1
    print(f"  captures writing both a .kin/.kin0 and a .seg : {cases}")
    print(f"  pairs in both files with identical IBD1Seg and IBD2Seg : {tot}")
    print(f"    same PropIBD in the two files      : {same}")
    print(f"    DIFFERENT PropIBD in the two files : {diff}"
          f"  ({100 * diff / max(tot, 1):.1f} %)")
    print("  -> no single expression can reproduce both files; each writer needs its own.")


# ---------------------------------------------------------------------------
# §2  PropIBD is a function of the printed columns
# ---------------------------------------------------------------------------

def section_printed_columns():
    print("\n=== §2  .seg's PropIBD is a function of its own two printed columns")
    print("  For each row, is the printed PropIBD the 4-dp rounding of")
    print("  (printed IBD2Seg) + (printed IBD1Seg)/2 ?  In half-ulp integers that value")
    print("  is n = 2*i2 + i1 and the printed cell is m, so a consistent row has")
    print("  n == 2m (unambiguous) or |n - 2m| == 1 (an exact decimal tie).")
    rows = []
    for _, t in every("king.seg"):
        rows += list(t.values())
    ins = up = down = out = 0
    for s1, s2, sp in rows:
        n = 2 * units(s2) + units(s1)
        m = units(sp)
        if n % 2 == 0 and n == 2 * m:
            ins += 1
        elif n % 2 == 1 and n == 2 * m - 1:
            up += 1
        elif n % 2 == 1 and n == 2 * m + 1:
            down += 1
        else:
            out += 1
    print(f"  .seg rows                              : {len(rows)}")
    print(f"    unambiguous, and correct             : {ins}")
    print(f"    exact tie, the reference rounded UP   : {up}")
    print(f"    exact tie, the reference rounded DOWN : {down}")
    print(f"    INCONSISTENT (would refute the rule)  : {out}")
    print("  The ties go both ways, so the tie-break is arithmetic, not a convention.")

    print("\n  The same test on the other files that carry the three columns:")
    for name, keys in (("king.kin", W3), ("king.kin0", W4), ("kingX.kin", W3),
                       ("kingcluster.kin", W3)):
        rs = []
        for _, t in every(name, keys):
            rs += list(t.values())
        if not rs:
            continue
        bad = 0
        for s1, s2, sp in rs:
            n, m = 2 * units(s2) + units(s1), units(sp)
            if not ((n % 2 == 0 and n == 2 * m) or (n % 2 == 1 and abs(n - 2 * m) == 1)):
                bad += 1
        verdict = "consistent" if bad == 0 else f"REFUTED on {bad}"
        print(f"    {name:<16} {len(rs):>5} rows   {verdict}")
    print("  -> the rule is .seg's alone. open-king reproduces the others byte for byte")
    print("     with the full-precision value, which is what they actually use.")


# ---------------------------------------------------------------------------
# §3  which expression, exactly
# ---------------------------------------------------------------------------

def section_expression():
    print("\n=== §3  which expression — the ties decide, and only one survives them")
    rows = []
    for _, t in every("king.seg"):
        rows += list(t.values())
    import decimal

    def half_up(x):
        return "%.4f" % float(decimal.Decimal(x)
                              .quantize(decimal.Decimal("0.0001"),
                                        rounding=decimal.ROUND_HALF_UP))

    cand = {
        "i2*1e-4 + i1*5e-5": lambda a, b: "%.4f" % (b * 1e-4 + a * 5e-5),
        "i2*0.0001 + i1*0.00005": lambda a, b: "%.4f" % (b * 0.0001 + a * 0.00005),
        "(i1 + 2*i2) * 5e-5": lambda a, b: "%.4f" % ((a + 2 * b) * 5e-5),
        "integer round-half-up": lambda a, b: "%.4f" % (((a + 2 * b) + 1) // 2 / 10000.0),
        "i2/10000 + i1/20000": lambda a, b: "%.4f" % (b / 10000.0 + a / 20000.0),
        "(i1 + 2*i2)/20000": lambda a, b: "%.4f" % ((a + 2 * b) / 20000.0),
        "(i1/2 + i2)/10000": lambda a, b: "%.4f" % ((a / 2 + b) / 10000.0),
        "printed doubles, b + a/2": lambda a, b: "%.4f" % (b / 1e4 + (a / 1e4) / 2),
        "printed doubles, half-up": lambda a, b: half_up(b / 1e4 + (a / 1e4) / 2),
    }
    for name, fn in cand.items():
        ok = sum(1 for s1, s2, sp in rows if fn(units(s1), units(s2)) == sp)
        mark = "   <== exact" if ok == len(rows) else ""
        print(f"    {name:<26} {ok:>5} / {len(rows)}{mark}")


# ---------------------------------------------------------------------------
# §4  the row order
# ---------------------------------------------------------------------------

def sample_index(ds):
    fam = [l.split()[:2] for l in open(os.path.join(DATA, ds + ".fam")) if l.strip()]
    return {(f, i): k for k, (f, i) in enumerate(fam)}


def section_row_order():
    print("\n=== §4  the row order: pairs are listed by 16-sample block")
    files = []
    for d in sorted(glob.glob(os.path.join(GOLDEN, "*", "*"))):
        p = os.path.join(d, "king.seg")
        if os.path.exists(p):
            files.append((os.path.basename(d).split("__")[0], p))
    if not files:
        print("  no captures found")
        return
    try:
        seqs = []
        for ds, p in files:
            idx = sample_index(ds)
            with open(p) as fh:
                fh.readline()
                seqs.append([(idx[(r[0], r[1])], idx[(r[2], r[3])])
                             for r in (l.rstrip("\n").split("\t") for l in fh)])
    except FileNotFoundError:
        print("  input corpus not generated; run tests/parity/run_parity.py once first")
        return

    def works(b):
        return all(s == sorted(s, key=lambda q: (q[0] // b, q[1] // b, q[0], q[1]))
                   for s in seqs)

    good = [b for b in range(2, 80) if works(b)]
    print(f"  .seg captures graded            : {len(files)}")
    print(f"  block sizes reproducing them all: {good}")
    print("  Plain index order is block size = infinity, and it fails: multifam (20")
    print("  samples) finishes its first block before writing any pair that reaches into")
    print("  the second. threegen (12 samples) rules out everything below 12 from the")
    print("  other side, and inside that window only 16 survives bigish (200 samples).")


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "prop"):
        section_self_disagreement()
        section_printed_columns()
        section_expression()
    if which in ("all", "order"):
        section_row_order()
