#!/usr/bin/env python3
"""Differential parity harness for open-king vs. the KING 2.3.2 reference binary.

Every directory under ``tests/parity/golden/<group>/<case>/`` is one captured
reference invocation:

    cmd.txt        the argv (placeholders {KING} / {DATA} / {ALT})
    exitcode.txt   reference exit status
    stdout.txt     reference stdout, verbatim (contains \\r progress tokens)
    stderr.txt     reference stderr, verbatim
    <everything else>  the output files the reference wrote into its cwd

This harness replays each cmd.txt with *our* binary in a fresh temp directory
and diffs the result against the capture: exit status, stdout, stderr, the set
of files produced, and the bytes of every file.  Known-nondeterministic stdout
lines (timestamps, thread counts, progress percentages, absolute input paths)
are normalized on BOTH sides first; everything else must match byte for byte.

Run ``--impl`` against the reference binary itself to prove the normalization is
complete: reference vs. reference must be 100% PASS.

Python 3 standard library only.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import difflib
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

# --------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------

PARITY_DIR = Path(__file__).resolve().parent
DEFAULT_GOLDEN = PARITY_DIR / "golden"
DEFAULT_DATA = PARITY_DIR / "work" / "data"
DEFAULT_ALT = PARITY_DIR / "work" / "alt"
GENERATOR = PARITY_DIR / "generate_corpus.py"
ALT_MAKER = PARITY_DIR / "make_alt_inputs.py"

# Files inside a case directory that describe the run rather than being output.
META_FILES = {"cmd.txt", "stdout.txt", "stderr.txt", "exitcode.txt",
              "MD5SUMS.txt", "README.txt"}

DATASETS = ["trio", "nuclear", "threegen", "multifam", "dups", "missing",
            "monomorphic", "sexchr", "unrelated", "admixed", "singleton",
            "pair", "bigish"]

# --------------------------------------------------------------------------
# Normalization
#
# Each rule is (name, kind, fn).  `kind` is "line" (applied to every physical
# line) or "block" (applied to the whole line list).  Rules run in list order.
# --------------------------------------------------------------------------

RE_TS = re.compile(
    r"(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun) [A-Z][a-z]{2} ?[ \d]?\d "
    r"\d{2}:\d{2}:\d{2} \d{4}"
)
RE_PROGRESS_ONLY = re.compile(r"^(?:\d{1,3}%\r?)+$")
RE_PROGRESS_HEAD = re.compile(r"^(?:\d{1,3}%\r?)+")
RE_NCPU = re.compile(r"^(\s*)\d+(\s*CPU cores are used)")
RE_NCPU_MID = re.compile(r"(with )\d+( CPU cores)")
RE_PATH = re.compile(r"/(?:[^\s()\[\],]+/)+([^/\s()\[\],]+\.(?:bed|bim|fam))")
RE_NOSCREEN = re.compile(r"(--noscreen \[)-?\d+(\])")

BANNER_START = "The following parameters are in effect:"
RE_BANNER_ENTRY = re.compile(r"^(\s*)([^:]*?)\s+:\s?(.*)$")
RE_BANNER_CONT = re.compile(r"^\s{20,}\S")


def _n_timestamp(line: str) -> str:
    """R1: wall-clock ctime() stamps -> <TS>, anywhere in the line."""
    return RE_TS.sub("<TS>", line)


def _n_ncpu(line: str) -> str:
    """R3: host/thread-count dependent 'N CPU cores are used' -> <NCPU>."""
    return RE_NCPU_MID.sub(r"\1<NCPU>\2", RE_NCPU.sub(r"\1<NCPU>\2", line))


def _n_paths(line: str) -> str:
    """R4: absolute .bed/.bim/.fam paths -> {DIR}/<basename>."""
    return RE_PATH.sub(r"{DIR}/\1", line)


def _n_noscreen(line: str) -> str:
    """R5: uninitialized --noscreen default -> <NOSCREEN>."""
    return RE_NOSCREEN.sub(r"\1<NOSCREEN>\2", line)


LINE_RULES = [
    ("R1 timestamp", _n_timestamp),
    ("R3 cpu-count", _n_ncpu),
    ("R4 input-path", _n_paths),
    ("R5 noscreen", _n_noscreen),
]


def _n_progress(lines: list[str]) -> list[str]:
    """R2: drop pure-progress lines, strip leading progress runs from the rest.

    KING writes "%d%%\\r" with no newline, so progress tokens glue onto the head
    of whatever message follows.  How many tokens appear depends on the thread
    count, so they carry no signal.  Only leading runs are stripped, which keeps
    a literal value like "rate > 100%" intact.
    """
    out = []
    for line in lines:
        if RE_PROGRESS_ONLY.match(line):
            continue
        out.append(RE_PROGRESS_HEAD.sub("", line))
    return out


def _n_banner(lines: list[str]) -> list[str]:
    """R6: unwrap the "parameters are in effect" echo block.

    KING word-wraps that block at a fixed column, so the *length of the input
    path* decides where the line breaks fall and which entries share a line.
    Continuation lines are folded back into their entry and interior whitespace
    is collapsed, making the block compare structurally instead of by column.
    """
    try:
        start = next(i for i, l in enumerate(lines) if l.strip() == BANNER_START)
    except StopIteration:
        return lines
    end = len(lines)
    for i in range(start + 1, len(lines)):
        if lines[i].startswith("KING starts at") or lines[i].startswith("KING 2."):
            end = i
            break

    folded: list[str] = []
    for line in lines[start:end]:
        if folded and RE_BANNER_CONT.match(line) and " : " not in line:
            folded[-1] = folded[-1].rstrip() + " " + line.strip()
        else:
            folded.append(line.rstrip())

    normalized = []
    for line in folded:
        m = RE_BANNER_ENTRY.match(line)
        if m and line.strip():
            label, body = m.group(2).strip(), " ".join(m.group(3).split())
            normalized.append(f"{label} : {body}".rstrip())
        else:
            normalized.append(line)
    return lines[:start] + normalized + lines[end:]


BLOCK_RULES = [
    ("R2 progress-tokens", _n_progress),
    ("R6 banner-unwrap", _n_banner),
]


def normalize_stream(raw: bytes) -> str:
    """Apply every normalization rule to a captured stdout/stderr stream."""
    text = raw.decode("utf-8", errors="replace")
    lines = text.split("\n")
    for _, fn in BLOCK_RULES:
        lines = fn(lines)
    out = []
    for line in lines:
        for _, fn in LINE_RULES:
            line = fn(line)
        out.append(line)
    return "\n".join(out)


# --------------------------------------------------------------------------
# Known-bad reference output (excluded from byte-diff)
# --------------------------------------------------------------------------

def racy_reason(case: "Case", relpath: str) -> str | None:
    """KING 2.3.2 races on the X between-family writer whenever threads > 1.

    Multiple threads append to one unlocked FILE*, so records tear mid-field and
    the file differs on every run.  Any capture of a "<prefix>X.kin0" made
    without --cpus 1 is therefore not a diffable golden.  Match on the SUFFIX:
    --prefix renames the file.
    """
    name = os.path.basename(relpath)
    if name.endswith("X.kin0") and "--cpus 1" not in case.display_cmd:
        return "X.kin0 race (KING 2.3.2, threads>1)"
    if relpath in case.diff_exclude:
        return "flagged diff_exclude in runs.json"
    return None


# --------------------------------------------------------------------------
# Case discovery
# --------------------------------------------------------------------------

class Case:
    __slots__ = ("group", "name", "dir", "tokens", "display_cmd",
                 "expected_exit", "diff_exclude")

    def __init__(self, group, name, cdir, tokens, expected_exit, diff_exclude):
        self.group = group
        self.name = name
        self.dir = cdir
        self.tokens = tokens
        self.display_cmd = " ".join(tokens)
        self.expected_exit = expected_exit
        self.diff_exclude = diff_exclude

    @property
    def key(self) -> str:
        return f"{self.group}/{self.name}"

    def golden_files(self) -> dict:
        out = {}
        for root, _dirs, files in os.walk(self.dir):
            for f in files:
                p = Path(root) / f
                rel = str(p.relative_to(self.dir))
                if rel in META_FILES:
                    continue
                out[rel] = p
        return out

    def argv(self, binary: str, data: Path, alt: Path) -> list[str]:
        argv = []
        for i, tok in enumerate(self.tokens):
            if i == 0 and tok in ("king", "{KING}"):
                argv.append(str(binary))
                continue
            tok = tok.replace("{KING}", str(binary))
            tok = tok.replace("{DATA}", str(data))
            tok = tok.replace("{ALT}", str(alt))
            argv.append(tok)
        return argv


def parse_cmd(path: Path) -> list[str]:
    """cmd.txt comes in two captured shapes: one shell line, or one token/line."""
    raw = path.read_text()
    lines = [l for l in raw.split("\n") if l.strip()]
    if len(lines) > 1:
        return [l.strip() for l in lines]
    return shlex.split(lines[0]) if lines else []


def discover(golden: Path, filt: str | None, include_analysis: bool) -> list[Case]:
    cases: list[Case] = []
    if not golden.is_dir():
        die(f"golden dir not found: {golden}")
    for group_dir in sorted(p for p in golden.iterdir() if p.is_dir()):
        group = group_dir.name
        excl_map = {}
        runs_json = group_dir / "runs.json"
        if runs_json.is_file():
            try:
                for rec in json.loads(runs_json.read_text()):
                    if rec.get("diff_exclude"):
                        excl_map[rec["dir"]] = set(rec["diff_exclude"])
            except (ValueError, KeyError):
                pass
        for cdir in sorted(p for p in group_dir.rglob("*") if p.is_dir()):
            if not (cdir / "cmd.txt").is_file():
                continue
            rel = cdir.relative_to(group_dir)
            if not include_analysis and any(part.startswith("_") for part in rel.parts):
                continue
            name = str(rel)
            key = f"{group}/{name}"
            if filt and filt not in key:
                continue
            try:
                expected_exit = int((cdir / "exitcode.txt").read_text().strip())
            except (OSError, ValueError):
                expected_exit = 0
            cases.append(Case(group, name, cdir, parse_cmd(cdir / "cmd.txt"),
                              expected_exit, excl_map.get(name, set())))
    return cases


# --------------------------------------------------------------------------
# Input corpus
# --------------------------------------------------------------------------

def ensure_inputs(cases: list[Case], data: Path, alt: Path, verbose: bool) -> None:
    missing = [ds for ds in DATASETS
               if not all((data / f"{ds}.{ext}").is_file()
                          for ext in ("bed", "bim", "fam"))]
    if missing:
        say(f"regenerating {len(missing)} dataset(s) into {data}: {' '.join(missing)}")
        data.mkdir(parents=True, exist_ok=True)
        cmd = [sys.executable, str(GENERATOR), "--outdir", str(data)]
        if len(missing) < len(DATASETS):
            cmd += ["--only"] + missing
        subprocess.run(cmd, check=True,
                       stdout=None if verbose else subprocess.DEVNULL)

    wanted = set()
    for case in cases:
        for tok in case.tokens:
            if tok.startswith("{ALT}/"):
                wanted.add(tok[len("{ALT}/"):])
    absent = [w for w in sorted(wanted)
              if not (alt / w).is_file() and "no_such_file" not in w]
    if absent:
        say(f"regenerating alternate inputs into {alt} ({len(absent)} missing)")
        alt.mkdir(parents=True, exist_ok=True)
        subprocess.run([sys.executable, str(ALT_MAKER),
                        "--datadir", str(data), "--outdir", str(alt)],
                       check=True,
                       stdout=None if verbose else subprocess.DEVNULL)


# --------------------------------------------------------------------------
# Numeric delta analysis
# --------------------------------------------------------------------------

def _split_row(line: str) -> list[str]:
    return line.split("\t") if "\t" in line else line.split()


def _num(tok: str):
    try:
        v = float(tok)
    except ValueError:
        return None
    return v


class ColStat:
    __slots__ = ("name", "n_diff", "max_abs", "max_abs_at", "max_rel",
                 "max_rel_at", "n_text", "decimals")

    def __init__(self, name):
        self.name = name
        self.n_diff = 0
        self.max_abs = 0.0
        self.max_abs_at = None
        self.max_rel = 0.0
        self.max_rel_at = None
        self.n_text = 0
        self.decimals = 0

    @property
    def ulp(self) -> float:
        """One unit in the last printed place of this column."""
        return 10.0 ** -self.decimals

    @property
    def ulps(self) -> float:
        return self.max_abs / self.ulp if self.decimals else float("inf")


def _decimals(tok: str) -> int:
    """Printed decimal places, e.g. '0.2505' -> 4.  KING uses fixed widths."""
    if "." not in tok or "e" in tok or "E" in tok:
        return 0
    return len(tok.rsplit(".", 1)[1])


def numeric_report(golden: bytes, ours: bytes, max_cols: int = 40):
    """Per-column worst absolute/relative delta, to separate rounding from bugs.

    Returns (verdict, [lines]) or (None, []) when the files are not tabular.
    """
    try:
        gt = golden.decode("utf-8")
        ot = ours.decode("utf-8")
    except UnicodeDecodeError:
        return None, []
    gl = [l for l in gt.split("\n") if l != ""]
    ol = [l for l in ot.split("\n") if l != ""]
    if not gl or not ol:
        return None, []

    head = _split_row(gl[0])
    has_header = bool(head) and all(_num(t) is None for t in head)
    names = head if has_header else [f"c{i + 1}" for i in range(len(_split_row(gl[0])))]
    gbody, obody = (gl[1:], ol[1:]) if has_header else (gl, ol)

    notes = []
    if len(gbody) != len(obody):
        notes.append(f"row count {len(gbody)} (golden) vs {len(obody)} (ours)"
                     f" - aligning by index over the common prefix")
    if has_header and gl[0] != ol[0]:
        notes.append(f"header differs: {gl[0]!r} vs {ol[0]!r}")

    stats = {}
    n_rows = min(len(gbody), len(obody))
    for r in range(n_rows):
        gr, orow = _split_row(gbody[r]), _split_row(obody[r])
        for c in range(min(len(gr), len(orow), max_cols)):
            if gr[c] == orow[c]:
                continue
            name = names[c] if c < len(names) else f"c{c + 1}"
            st = stats.setdefault(c, ColStat(name))
            st.n_diff += 1
            gv, ov = _num(gr[c]), _num(orow[c])
            if gv is None or ov is None:
                st.n_text += 1
                if st.max_abs_at is None:
                    st.max_abs_at = (r + 1, gr[c], orow[c])
                continue
            st.decimals = max(st.decimals, _decimals(gr[c]))
            ad = abs(gv - ov)
            if ad > st.max_abs:
                st.max_abs, st.max_abs_at = ad, (r + 1, gr[c], orow[c])
            if gv != 0.0:
                rd = ad / abs(gv)
                if rd > st.max_rel:
                    st.max_rel, st.max_rel_at = rd, (r + 1, gr[c], orow[c])
        if len(gr) != len(orow):
            st = stats.setdefault(-1, ColStat("<field count>"))
            st.n_diff += 1
            st.n_text += 1
            if st.max_abs_at is None:
                st.max_abs_at = (r + 1, str(len(gr)), str(len(orow)))

    if not stats and not notes:
        return None, []

    lines = list(notes)
    any_text = any(s.n_text for s in stats.values())
    worst_rel = max((s.max_rel for s in stats.values()), default=0.0)
    worst_abs = max((s.max_abs for s in stats.values()), default=0.0)
    worst_ulps = max((s.ulps for s in stats.values()), default=0.0)
    for c in sorted(stats):
        s = stats[c]
        parts = [f"col {c + 1:>2} {s.name:<12} {s.n_diff:>6} row(s) differ"]
        if s.max_abs_at:
            r, g, o = s.max_abs_at
            parts.append(f"max|d| {s.max_abs:.3e} @row {r} ({g} -> {o})")
        if s.decimals and s.max_abs:
            parts.append(f"= {s.ulps:.1f} ulp")
        if s.max_rel_at:
            parts.append(f"max rel {s.max_rel:.3e}")
        if s.n_text:
            parts.append(f"{s.n_text} non-numeric")
        lines.append("  " + " | ".join(parts))

    # Relative delta explodes near zero (0.0001 vs 0.0002 is rel 1.0 but only one
    # unit in the last printed place), so the rounding-vs-algorithm call is made
    # on ULPs of KING's fixed-width printing, with abs/rel reported alongside.
    if notes or any_text:
        verdict = "STRUCTURAL - row/field/text mismatch, not a rounding issue"
    elif worst_ulps <= 1.0 + 1e-9:
        verdict = "ROUNDING - every delta within 1 ulp of the printed precision"
    elif worst_ulps <= 2.0 + 1e-9:
        verdict = f"NEAR-ROUNDING - worst delta {worst_ulps:.1f} ulp (check tie-breaking)"
    else:
        verdict = f"ALGORITHMIC - worst delta {worst_ulps:.1f} ulp, too large for rounding"
    verdict += f" [max|d| {worst_abs:.3e}, max rel {worst_rel:.3e}]"
    lines.insert(0, verdict)
    return verdict, lines


def unified(golden: str, ours: str, label: str, max_lines: int) -> list[str]:
    """Unified diff of the first `max_lines` mismatching lines.

    When both sides have the same line count (the common case for a numeric
    regression) the lines are paired up positionally, which reads far better
    than difflib's one-giant-hunk output on a 19,000-row file.
    """
    gl, ol = golden.split("\n"), ours.split("\n")
    if len(gl) == len(ol):
        out = [f"    --- golden/{label}", f"    +++ ours/{label}"]
        shown = 0
        for i, (g, o) in enumerate(zip(gl, ol)):
            if g == o:
                continue
            if shown * 2 >= max_lines:
                out.append(f"    ... (diff truncated at {max_lines} lines)")
                break
            out.append(f"    @@ line {i + 1} @@")
            out.append(f"    -{g}")
            out.append(f"    +{o}")
            shown += 1
        return out
    diff = difflib.unified_diff(gl, ol, fromfile=f"golden/{label}",
                                tofile=f"ours/{label}", lineterm="", n=1)
    out = []
    for i, line in enumerate(diff):
        if i >= max_lines:
            out.append(f"    ... (diff truncated at {max_lines} lines)")
            break
        out.append("    " + line)
    return out


# --------------------------------------------------------------------------
# Running one case
# --------------------------------------------------------------------------

class Result:
    __slots__ = ("case", "status", "notes", "detail", "seconds", "skipped",
                 "n_compared")

    def __init__(self, case):
        self.case = case
        self.status = "PASS"
        self.notes = []
        self.detail = []
        self.seconds = 0.0
        self.skipped = []
        self.n_compared = 0

    def fail(self, note):
        self.status = "FAIL"
        self.notes.append(note)


def collect_files(root: Path) -> dict:
    out = {}
    for dirpath, _dirs, files in os.walk(root):
        for f in files:
            p = Path(dirpath) / f
            out[str(p.relative_to(root))] = p
    return out


def run_case(case: Case, binary: str, data: Path, alt: Path, timeout: float,
             max_diff_lines: int, keep: bool) -> Result:
    res = Result(case)
    tmp = Path(tempfile.mkdtemp(prefix=f"parity-{case.group}-"))
    try:
        argv = case.argv(binary, data, alt)
        t0 = time.time()
        try:
            proc = subprocess.run(argv, cwd=str(tmp), capture_output=True,
                                  timeout=timeout)
            rc, out, err = proc.returncode, proc.stdout, proc.stderr
        except subprocess.TimeoutExpired:
            res.fail(f"TIMEOUT after {timeout:g}s")
            return res
        except OSError as exc:
            res.fail(f"cannot execute: {exc}")
            return res
        res.seconds = time.time() - t0

        if rc != case.expected_exit:
            res.fail(f"exit {rc} != {case.expected_exit}")

        for stream, got in (("stdout", out), ("stderr", err)):
            gpath = case.dir / f"{stream}.txt"
            gold_raw = gpath.read_bytes() if gpath.is_file() else b""
            g, o = normalize_stream(gold_raw), normalize_stream(got)
            if g != o:
                res.fail(f"{stream}!=")
                res.detail.append(f"  --- {stream} (normalized) ---")
                res.detail += unified(g, o, f"{stream}.txt", max_diff_lines)

        golden_files = case.golden_files()
        our_files = collect_files(tmp)
        missing = sorted(set(golden_files) - set(our_files))
        extra = sorted(set(our_files) - set(golden_files))
        for m in missing:
            res.fail(f"missing:{m}")
        for e in extra:
            res.fail(f"extra:{e}")

        for rel in sorted(set(golden_files) & set(our_files)):
            reason = racy_reason(case, rel)
            gb = golden_files[rel].read_bytes()
            ob = our_files[rel].read_bytes()
            if reason:
                res.skipped.append(f"{rel} [{reason}]")
                continue
            res.n_compared += 1
            if gb == ob:
                continue
            verdict, nlines = numeric_report(gb, ob)
            tag = "num" if verdict else "bytes"
            res.fail(f"{rel}!=({tag})")
            res.detail.append(f"  --- {rel} ({len(gb)} B golden vs {len(ob)} B ours) ---")
            if nlines:
                res.detail += ["    " + l for l in nlines]
            res.detail += unified(gb.decode("utf-8", "replace"),
                                  ob.decode("utf-8", "replace"),
                                  rel, max_diff_lines)
        return res
    finally:
        if keep:
            say(f"kept work dir for {case.key}: {tmp}")
        else:
            shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# --update: re-capture goldens from the reference binary
# --------------------------------------------------------------------------

def update_case(case: Case, ref: str, data: Path, alt: Path, timeout: float) -> str:
    tmp = Path(tempfile.mkdtemp(prefix="parity-update-"))
    try:
        proc = subprocess.run(case.argv(ref, data, alt), cwd=str(tmp),
                              capture_output=True, timeout=timeout)
        produced = collect_files(tmp)
        for rel in sorted(set(case.golden_files()) - set(produced)):
            (case.dir / rel).unlink()
        for rel, path in produced.items():
            dest = case.dir / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(path, dest)
        (case.dir / "stdout.txt").write_bytes(proc.stdout)
        (case.dir / "stderr.txt").write_bytes(proc.stderr)
        (case.dir / "exitcode.txt").write_text(f"{proc.returncode}\n")
        return f"{len(produced)} file(s), exit {proc.returncode}"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# --------------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------------

def say(msg: str) -> None:
    sys.stdout.flush()
    print(f"[parity] {msg}", file=sys.stderr, flush=True)


def resolve_binary(spec: str) -> str:
    """Absolutize a binary path: cases run with cwd=<temp dir>, so a relative
    --impl like 'target/release/king' would otherwise not be found."""
    found = shutil.which(spec)
    if found:
        return str(Path(found).resolve())
    return str(Path(spec).expanduser().resolve())


def die(msg: str, code: int = 2):
    print(f"[parity] ERROR: {msg}", file=sys.stderr)
    raise SystemExit(code)


def main() -> int:
    ap = argparse.ArgumentParser(
        prog="run_parity.py",
        description="Differential parity harness: replay every captured KING "
                    "invocation with our binary and diff the results.")
    ap.add_argument("--impl", required=True,
                    help="path to the binary under test (point it at the "
                         "reference binary to self-check the harness)")
    ap.add_argument("--ref", help="path to the KING 2.3.2 reference binary; "
                                  "required by --update")
    ap.add_argument("--data", type=Path, default=DEFAULT_DATA,
                    help=f"dataset dir (default {DEFAULT_DATA})")
    ap.add_argument("--alt", type=Path, default=DEFAULT_ALT,
                    help=f"alternate --fam/--bim input dir (default {DEFAULT_ALT})")
    ap.add_argument("--golden", type=Path, default=DEFAULT_GOLDEN,
                    help=f"golden capture root (default {DEFAULT_GOLDEN})")
    ap.add_argument("--filter", dest="filt",
                    help="only run cases whose 'group/name' contains this substring")
    ap.add_argument("--update", action="store_true",
                    help="re-capture goldens from --ref instead of testing")
    ap.add_argument("--include-analysis", action="store_true",
                    help="also replay the _analysis/ side runs (their output "
                         "files were pruned at capture time; expect failures)")
    ap.add_argument("--jobs", type=int, default=min(8, os.cpu_count() or 4),
                    help="parallel case runners (default: min(8, ncpu))")
    ap.add_argument("--timeout", type=float, default=300.0,
                    help="per-case timeout in seconds (default 300)")
    ap.add_argument("--max-diff-lines", type=int, default=20,
                    help="max unified-diff lines printed per file with -v")
    ap.add_argument("--keep", action="store_true", help="keep per-case temp dirs")
    ap.add_argument("--json", type=Path, help="write machine-readable results here")
    ap.add_argument("-q", "--quiet", action="store_true",
                    help="print only the summary and failing rows")
    ap.add_argument("-v", "--verbose", action="count", default=0,
                    help="-v prints diffs for failing cases")
    args = ap.parse_args()

    # Cases run with cwd set to a temp dir, so the binary MUST be absolute.
    impl = resolve_binary(args.impl)
    if not args.update and not Path(impl).is_file():
        die(f"--impl not found: {args.impl}")

    cases = discover(args.golden, args.filt, args.include_analysis)
    if not cases:
        die("no cases matched")

    ensure_inputs(cases, args.data.resolve(), args.alt.resolve(), args.verbose > 0)

    if args.update:
        if not args.ref:
            die("--update refuses to run without --ref "
                "(it overwrites goldens; name the reference binary explicitly)")
        ref = resolve_binary(args.ref)
        if not Path(ref).is_file():
            die(f"--ref not found: {args.ref}")
        say(f"RE-CAPTURING {len(cases)} case(s) from {ref}")
        for case in cases:
            info = update_case(case, ref, args.data.resolve(), args.alt.resolve(),
                               args.timeout)
            print(f"UPDATED  {case.key:<50} {info}", flush=True)
        say(f"re-captured {len(cases)} case(s); review with 'git diff' before committing")
        return 0

    say(f"{len(cases)} case(s), impl={impl}, jobs={args.jobs}")
    results: list[Result] = []
    t0 = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        futs = {pool.submit(run_case, c, impl, args.data.resolve(),
                            args.alt.resolve(), args.timeout,
                            args.max_diff_lines, args.keep): c for c in cases}
        for fut in concurrent.futures.as_completed(futs):
            results.append(fut.result())
    results.sort(key=lambda r: r.case.key)
    elapsed = time.time() - t0

    width = max(len(r.case.key) for r in results)
    n_pass = n_fail = 0
    n_skipped_files = n_compared = 0
    for r in results:
        n_skipped_files += len(r.skipped)
        n_compared += r.n_compared
        if r.status == "PASS":
            n_pass += 1
            if not args.quiet:
                note = ""
                if r.skipped:
                    note = f"  [{len(r.skipped)} file(s) diff-excluded]"
                print(f"PASS  {r.case.key:<{width}}{note}")
        else:
            n_fail += 1
            print(f"FAIL  {r.case.key:<{width}}  {'; '.join(r.notes[:6])}"
                  + (" ..." if len(r.notes) > 6 else ""))
            if args.verbose:
                print(f"      $ {r.case.display_cmd}")
                for line in r.detail:
                    print(line)

    print()
    print(f"parity: {n_pass} PASS, {n_fail} FAIL, {len(results)} total "
          f"({elapsed:.1f}s wall, {n_compared} output file(s) byte-compared, "
          f"{n_skipped_files} diff-excluded)")

    if args.json:
        args.json.write_text(json.dumps([
            {"case": r.case.key, "status": r.status, "notes": r.notes,
             "cmd": r.case.display_cmd, "seconds": round(r.seconds, 3),
             "diff_excluded": r.skipped}
            for r in results], indent=1) + "\n")
        say(f"wrote {args.json}")

    return 1 if n_fail else 0


if __name__ == "__main__":
    sys.exit(main())
