#!/usr/bin/env python3
"""Grade **our binary** against the reference on the `.seg` canvas batteries.

`segcanvas.py` and `ibd1canvas.py` measure the *reference* and cache every answer in
`segcanvas_measured.json` / `ibd1canvas_measured.json`. This replays exactly those canvases
with an open-king build and compares marker-interval readings, so what is graded is the Rust
engine against KING 2.3.2 — no Python model on either side. It is the instrument
`docs/PARITY.md` §5.0 quotes, and it is how `17-seg-caller.md` §14.10 was landed: a rule the
corpus cannot see at all is separated here, 5 723 against 6 000.

    python3 gradebinary.py ../../../target/release/open-king           # the IBD2 column
    python3 gradebinary.py ../../../target/release/open-king --ibd1    # the IBD1 column

**IBD2 mode** (default) reads `IBD2Seg` on the six families of `17-seg-caller.md` §14.6: the
exhaustive word sequences of length <= 4 and 5 over `{C, z, x, y}`, length 4 over the
eight-letter alphabet, three random seeds and three "rich" random seeds — 6 000 canvases.

**IBD1 mode** reads `IBD1Seg` on the families of `18-ibd1-caller.md` §7-§8: the exhaustive
sequences of length <= 4 over `{K, k, W}` and three IBD2-free random seeds, all of which
isolate the IBD1 caller, then the *mixed* seeds, on which both columns are graded at once so
the `IBD1Seg`/`IBD2Seg` overlap rule of §6 is on trial as well. The last mixed family runs at
`--seglength 8`, which is where the measured-but-unmodelled run merge of §9 fires; it is
listed separately and is expected to miss.

Every canvas must already have a cached reference reading; families that do not are reported
and skipped rather than silently re-measured.

**It never writes either measurement cache.** Those files are the reference's own answers and
a non-reference binary must not touch them: this reads them, and keeps its own answers in
`$TMPDIR` (override with `$GRADE_CACHE`). Working directories are temporary too.

Exit status is 0 iff every canvas outside the known-open family matches.
"""

from __future__ import annotations

import itertools
import json
import os
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fixlab as F      # noqa: E402
import segcanvas as S   # noqa: E402

JOBS = int(os.environ.get("GRADE_JOBS", "8"))
#: markers of disagreement tolerated — the ruler is good to about a tenth of one
TOL = 0.3
#: the default invocation every family uses unless it names its own
L1 = ("--seglength", "1")


def families():
    """IBD2 mode: `(name, extra, columns, [canvas, ...])`, in `17-…` §14.6's order."""
    out = [
        ("length<=4 over {C,z,x,y}", L1, (2,),
         [S.seq_canvas("".join(t)) for n in (1, 2, 3, 4)
          for t in itertools.product("Czxy", repeat=n)]),
        ("length-5 over {C,z,x,y}", L1, (2,),
         [S.seq_canvas("".join(t)) for t in itertools.product("Czxy", repeat=5)]),
        ("length-4 over {C,z,x,p,y,d,W,q}", L1, (2,),
         [S.seq_canvas9("".join(t)) for t in itertools.product("CzxpydWq", repeat=4)]),
    ]
    out += [("random seed %d" % s, L1, (2,), S.random_canvases(s, 80, 10))
            for s in (101, 777, 8081)]
    out += [("rich random seed %d" % s, L1, (2,), S.rich_canvases(s, 100))
            for s in (31337, 424242, 90210)]
    return out


def families_ibd1():
    """IBD1 mode: the families of `18-ibd1-caller.md` §7 and §8.

    The IBD2-free families grade `IBD1Seg` alone; the mixed ones grade both columns, which
    is what puts the §6 overlap rule on trial.
    """
    import ibd1canvas as I  # noqa: E402  (sets S.CACHE — the caller re-reads it)

    out = [
        ("length<=4 over {K,k,W}", L1, (1,),
         [I.seq_canvas("".join(t)) for n in (1, 2, 3, 4)
          for t in itertools.product("KkW", repeat=n)]),
    ]
    out += [("IBD2-free random seed %d" % s, L1, (1,), I.random_canvases(s, 80, 10))
            for s in (201, 5150, 99991)]
    out += [("mixed seed %d, L=%d" % (s, L), ("--seglength", "%d" % L), (1, 2),
             I.random_canvases(s, 60, 10, mixed=True))
            for s, L in ((3300, 1), (61803, 1), (3300, 4))]
    out += [("mixed seed 61803, L=8  [18-… §9]", ("--seglength", "8"), (1, 2),
             I.random_canvases(61803, 60, 10, mixed=True))]
    return out


class Runner:
    """Runs `binary` over canvases, caching answers outside the repo."""

    def __init__(self, binary, cache_path):
        self.binary = os.path.abspath(binary)
        self.cache_path = cache_path
        try:
            with open(cache_path) as fh:
                self.cache = json.load(fh)
        except (OSError, ValueError):
            self.cache = {}

    def _invoke(self, cv, extra):
        wd = tempfile.mkdtemp(prefix="gradebinary_")
        try:
            prefix = cv.build(wd)
            subprocess.run([self.binary, "-b", prefix + ".bed", "--ibdseg", *extra,
                            "--prefix", os.path.join(wd, "k")],
                           cwd=wd, capture_output=True, text=True)
            rows = F.parse_seg(os.path.join(wd, "k.seg"))
            return {"row": rows.get(("S00", "S01"))}
        finally:
            shutil.rmtree(wd, ignore_errors=True)

    def __call__(self, items):
        todo, seen = [], set()
        for cv, ex in items:
            k = cv.key(ex)
            if k not in self.cache and k not in seen:
                seen.add(k)
                todo.append((cv, ex))
        if todo:
            with ThreadPoolExecutor(max_workers=JOBS) as pool:
                got = list(pool.map(lambda t: self._invoke(*t), todo))
            for (cv, ex), res in zip(todo, got):
                self.cache[cv.key(ex)] = res
            with open(self.cache_path, "w") as fh:
                json.dump(self.cache, fh, indent=0, sort_keys=True)
        return [self.cache[cv.key(ex)] for cv, ex in items]


def disagrees(cv, ref_res, our_res, columns):
    """True iff the two readings of this canvas differ on any graded column."""
    for what in columns:
        want, mine = S.mk(cv, ref_res, what), S.mk(cv, our_res, what)
        if want is None or mine is None:
            if want is not mine:
                return True
        elif abs(want - mine) > TOL:
            return True
    return False


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    ibd1 = "--ibd1" in argv[1:]
    binary = args[0] if args else os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        "..", "..", "..", "target", "release", "open-king")
    tag = os.environ.get("GRADE_TAG", "openking")
    cache = os.environ.get("GRADE_CACHE") or os.path.join(
        tempfile.gettempdir(), "openking_gradebinary_%s%s.json"
        % (tag, "_ibd1" if ibd1 else ""))

    fams = families_ibd1() if ibd1 else families()
    # `families_ibd1()` repoints `segcanvas.CACHE` at this campaign's own reference file.
    with open(S.CACHE) as fh:
        ref = json.load(fh)
    run = Runner(binary, cache)

    print("binary: %s" % os.path.abspath(binary))
    print("column: %s" % ("IBD1Seg (+ IBD2Seg on the mixed families)" if ibd1
                          else "IBD2Seg"))
    print("reference readings: %s (%d cached)" % (S.CACHE, len(ref)))
    total = bad_total = 0
    open_total = open_bad = 0
    misses = []
    for name, extra, columns, cvs in fams:
        absent = [c for c in cvs if c.key(extra) not in ref]
        if absent:
            print("  %-38s SKIPPED — %d of %d canvases have no reference reading"
                  % (name, len(absent), len(cvs)))
            continue
        bad = [cv for cv, got in zip(cvs, run([(c, extra) for c in cvs]))
               if disagrees(cv, ref[cv.key(extra)], got, columns)]
        print("  %-38s %4d/%4d exact" % (name, len(cvs) - len(bad), len(cvs)))
        if "§9" in name:            # measured, deliberately unmodelled — not a regression
            open_total += len(cvs)
            open_bad += len(bad)
            continue
        total += len(cvs)
        bad_total += len(bad)
        misses += [(name, c, extra, columns) for c in bad]
    print("  %-38s %4d/%4d exact" % ("TOTAL (closed families)", total - bad_total, total))
    if open_total:
        print("  %-38s %4d/%4d exact   (known open, `18-…` §9)"
              % ("known-open family", open_total - open_bad, open_total))
    if misses:
        print("\n  misses (canvas, reference markers, ours):")
        fmt = lambda v: "-" if v is None else "%.2f" % v  # noqa: E731
        for name, c, extra, columns in misses[:60]:
            for what in columns:
                print("    %-34s %-22s IBD%d  ref %8s  ours %8s"
                      % (name, c.name, what,
                         fmt(S.mk(c, ref[c.key(extra)], what)),
                         fmt(S.mk(c, run([(c, extra)])[0], what))))
        if len(misses) > 60:
            print("    ... and %d more" % (len(misses) - 60))
    return 1 if bad_total else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
