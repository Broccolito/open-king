#!/usr/bin/env python3
"""Fixtures behind `docs/research/15-ibs-ibd2-rules.md` — the `--ibs` IBD2 columns.

Every experiment builds a fileset in which the test pair's per-word genotype pattern is
exact by construction, runs the reference binary with `--ibs`, and reads `MaxIBD2` /
`Pr_IBD2` back for that pair.  Nothing here reads KING's source.

    python3 ibs_rules.py            # all sections
    python3 ibs_rules.py 2 4        # only sections 2 and 4

Sections match the write-up: 1 ruler, 2 which words break a run, 3 `--ibs` vs `--ibdseg`,
4 the 10 Mb pair rule, 5 the 95-HetHet acceptance count, 6 the open boundary.
"""

import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixlab as F  # noqa: E402

WORD = 64
WORK = os.path.join(os.path.dirname(os.path.abspath(__file__)), "work", "ibs_rules")

# per-marker patterns, as (sample 0, sample 1) A1 dosages; 3 is a missing call
PAIR = {
    "hethet": [1, 1],      # both heterozygous: the only thing the acceptance count counts
    "zero": [0, 0],        # both homozygous for A2: clean, uninformative
    "hom1": [2, 2],        # both homozygous for A1: clean, and still uninformative here
    "ibs1": [1, 0],        # het against hom: the mismatch that breaks a run at five
    "ibs0": [2, 0],        # opposite homozygotes: irrelevant to this pass
    "miss": [3, 3],
}


def build(name, words, nw1=30, spacing=50_000, per_word_spacing=False, fringe=0,
          nsample=6, seed=3):
    """chr1 = `nw1` words of IBD1 carrier, chr2 = one entry of `words` per complete word.

    Each entry is a list of `(pattern, count)`; whatever is left of the word is filled
    with `hethet`.  `per_word_spacing` gives word *w* its own marker spacing so that a
    reported length identifies the called word interval uniquely.
    """
    F.SPACING = spacing
    nw2 = len(words)
    fix = F.Fixture(name, [(1, WORD * nw1), (2, WORD * nw2 + fringe)], nsample=nsample,
                    seed=seed, maf=0.5)
    fix.set_state(0, 0, WORD * nw1, F.IBD1)
    fix.set_state(1, 0, WORD * nw2 + fringe, F.IBD1)
    lo, _ = fix.chrom_span(1)
    for w, spec in enumerate(words):
        idx = 0
        for kind, cnt in spec:
            for _ in range(cnt):
                if idx >= WORD:
                    break
                m = lo + WORD * w + idx
                fix.pat_all[m] = PAIR[kind] + [0] * (nsample - 2)
                fix.noflip.add(m)
                idx += 1
        while idx < WORD:
            m = lo + WORD * w + idx
            fix.pat_all[m] = PAIR["hethet"] + [0] * (nsample - 2)
            fix.noflip.add(m)
            idx += 1
    fix._ibs = dict(nw1=nw1, nw2=nw2, per_word_spacing=per_word_spacing)
    return fix


def word_positions(nwords, base=40_000, step=137):
    pos, x = [], 1_000_000
    for w in range(nwords):
        for _ in range(WORD):
            pos.append(x)
            x += base + step * w
    return pos


def rewrite_bim(prefix, nw1, nw2, fringe=0):
    p1, p2 = word_positions(nw1), word_positions(nw2 + (fringe + WORD - 1) // WORD)
    out = []
    with open(prefix + ".bim") as f:
        for n, line in enumerate(f):
            v = line.rstrip("\n").split("\t")
            pos = p1[n] if n < len(p1) else p2[n - len(p1)]
            v[2], v[3] = f"{pos / 1e6:.6f}", str(pos)
            out.append("\t".join(v))
    with open(prefix + ".bim", "w") as f:
        f.write("\n".join(out) + "\n")


def probe(fix, extra=()):
    """(MaxIBD2, Pr_IBD2, workdir, prefix) for the test pair, or None on a FATAL."""
    wd = os.path.join(WORK, fix.name)
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd, exist_ok=True)
    prefix = fix.build(wd)
    if fix._ibs["per_word_spacing"]:
        rewrite_bim(prefix, fix._ibs["nw1"], fix._ibs["nw2"])
    cmd = [F.KING, "-b", prefix + ".bed", "--ibs", *extra,
           "--prefix", os.path.join(wd, "k")]
    r = subprocess.run(cmd, capture_output=True, text=True, cwd=wd)
    if "FATAL ERROR" in r.stdout:
        return None
    path = os.path.join(wd, "k.ibs0")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            d = dict(zip(hdr, line.rstrip("\n").split("\t")))
            if {d.get("ID1"), d.get("ID2")} == {"S00", "S01"} and "MaxIBD2" in d:
                return (float(d["MaxIBD2"]), float(d["Pr_IBD2"]), wd, prefix)
    return None


def decode(maxbp, nw2):
    """Which word interval has this length, under `per_word_spacing`?"""
    p = word_positions(nw2)
    return [(u, e) for u in range(nw2) for e in range(u, nw2)
            if abs(p[WORD * e + WORD - 1] - p[WORD * u] - maxbp) < 1]


WALL = [("ibs1", 64)]
CLEAN = []


def slug(name):
    """A fixture name that is also a directory name."""
    return "".join(c if c.isalnum() else "_" for c in name)


def section1():
    print("\n== 1. the ruler: trailing dirty words are swallowed when at most two remain")
    print(f"{'trailing dirty words':>22} {'MaxIBD2':>12} {'words spanned':>14}")
    for trail in range(6):
        words = [WALL] * 2 + [CLEAN] * 3 + [WALL] * trail
        r = probe(build(f"ruler_{trail}", words))
        span = None if not r or not r[0] else round((r[0] / 50_000 + 1) / WORD, 2)
        print(f"{trail:>22} {r and r[0]:>12} {str(span):>14}")


def section2():
    print("\n== 2. which words break a run (two-word gap unless stated)")
    gaps = [("64 mismatches", [("ibs1", 64)]),
            ("5 mismatches", [("ibs1", 5)]),
            ("4 mismatches", [("ibs1", 4)]),
            ("64 opposite homozygotes", [("ibs0", 64)]),
            ("64 missing", [("miss", 64)]),
            ("64 A1A1/A1A1", [("hom1", 64)]),
            ("4 mism + 60 opp hom", [("ibs1", 4), ("ibs0", 60)]),
            ("5 mism + 59 opp hom", [("ibs1", 5), ("ibs0", 59)])]
    for name, gap in gaps:
        row = []
        for g in (1, 2):
            words = [WALL] * 2 + [CLEAN] * 3 + [gap] * g + [CLEAN] * 3 + [WALL] * 2
            r = probe(build(f"gap_{slug(name)}_{g}", words))
            whole = (WORD * (len(words) - 2) - 1) * 50_000
            row.append("merged" if r and r[0] >= whole - 1 else "split")
        print(f"{name:>26}: g=1 {row[0]:>7}   g=2 {row[1]:>7}")


def section3():
    print("\n== 3. --ibs and --ibdseg disagree about IBS0 words")
    for W in (4, 6):
        words = [[("ibs0", 64)]] * 3 + [CLEAN] * W + [[("ibs0", 64)]] * 3
        fix = build(f"ibs0_bounded_{W}", words, nw1=60, per_word_spacing=True)
        r = probe(fix)
        if r is None:
            print(f"  W={W}: FATAL")
            continue
        maxbp, pr, wd, prefix = r
        subprocess.run([F.KING, "-b", prefix + ".bed", "--ibdseg",
                        "--prefix", os.path.join(wd, "s")],
                       capture_output=True, text=True, cwd=wd)
        seg = F.parse_seg(os.path.join(wd, "s.seg")).get(("S00", "S01"))
        print(f"  W={W}: --ibs MaxIBD2={maxbp:.0f} -> words {decode(maxbp, len(words))}"
              f"   --ibdseg IBD2Seg={seg and seg['IBD2Seg']}")


def section4():
    print("\n== 4. the 10 Mb rule gates the pair, not the call")
    for s in (45_000, 50_000, 52_353, 52_357, 55_000):
        words = [WALL] * 2 + [CLEAN] * 2 + [WALL] * 6
        r = probe(build(f"long_{s}", words, nw1=40, spacing=s))
        print(f"  spacing {s}: call {191 * s:>10} bp   MaxIBD2={r and r[0]:>12} "
              f"Pr_IBD2={r and r[1]}")
    words = ([WALL] * 2 + [CLEAN] * 2 + [WALL] * 2 + [CLEAN] * 4 + [WALL] * 2)
    r = probe(build("long_two_calls", words, nw1=30))
    print(f"  9.55 Mb call beside a 19.15 Mb one: MaxIBD2={r[0]:.0f} Pr_IBD2={r[1]} "
          f"-> numerator {r[1] * 134_300_000:.0f} (both calls)")


def hh_block(k, width):
    """`k` HetHet markers spread over `width` words, everything else A2A2 everywhere."""
    words = []
    left = k
    for _ in range(width):
        take = min(left, WORD)
        left -= take
        words.append([("hethet", take), ("zero", WORD - take)])
    return words


def section5():
    print("\n== 5. the acceptance count is 95 HetHet markers over the measured interval")
    for width in (2, 3, 4, 5):
        lo, hi = -1, 200
        while hi - lo > 1:
            mid = (lo + hi) // 2
            words = [WALL] * 2 + hh_block(mid, width) + [WALL] * 3
            r = probe(build(f"count_{width}_{mid}", words, nw1=40))
            if r and r[0] > 0:
                hi = mid
            else:
                lo = mid
        print(f"  block width {width} words: smallest k reported = {hi}")
    words = [WALL] * 2 + [[("hom1", 64)]] * 4 + [WALL] * 3
    r = probe(build("count_hom1", words, nw1=40))
    print(f"  256 A1A1/A1A1 markers and no HetHet: MaxIBD2={r and r[0]}")
    for inf in (0, 54, 59):
        lo, hi = -1, 130
        while hi - lo > 1:
            mid = (lo + hi) // 2
            words = ([WALL] * 2 + hh_block(mid, 2)
                     + [[("ibs1", 5), ("hethet", inf), ("zero", 59 - inf)]] + [WALL] * 3)
            r = probe(build(f"win_{inf}_{mid}", words, nw1=40))
            if r and r[0] > 0:
                hi = mid
            else:
                lo = mid
        print(f"  terminating word carrying {inf:>2} HetHet: block needs k = {hi}")


def section6():
    print("\n== 6. open: sustained low-grade mismatch")
    print(f"{'mismatches/word':>16} {'smallest HetHet/word reported':>30}")
    for m in (0, 1, 2, 3):
        lo, hi = -1, 63
        while hi - lo > 1:
            mid = (lo + hi) // 2
            spec = [("ibs1", m), ("hethet", mid), ("zero", max(0, WORD - m - mid))]
            words = [WALL] * 2 + [spec] * 8 + [WALL] * 4
            r = probe(build(f"open_{m}_{mid}", words, nw1=60))
            if r and r[0] > 0:
                hi = mid
            else:
                lo = mid
        print(f"{m:>16} {hi:>30}")
    for tag, block in (("8 clean then 8 at m=4",
                        [[("ibs1", 0)]] * 8 + [[("ibs1", 4)]] * 8),
                       ("8 at m=4 then 8 clean",
                        [[("ibs1", 4)]] * 8 + [[("ibs1", 0)]] * 8)):
        words = [WALL] * 2 + block + [WALL] * 4
        r = probe(build(f"order_{slug(tag)}", words, nw1=60, per_word_spacing=True))
        print(f"  {tag}: words {decode(r[0], len(words)) if r else None}")


SECTIONS = {1: section1, 2: section2, 3: section3, 4: section4, 5: section5, 6: section6}

if __name__ == "__main__":
    want = [int(a) for a in sys.argv[1:]] or sorted(SECTIONS)
    for n in want:
        SECTIONS[n]()
