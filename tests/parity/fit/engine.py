"""A faithful Python mirror of `crates/king-core/src/ibdseg.rs`, made variable.

Why this exists: the remaining parity gap is *boundary geometry* on a minority of called
segments, and iterating on it in Rust costs a rebuild plus a full corpus replay per idea.
This module reimplements the committed engine exactly — same word rule, same bridging,
same informativeness gate, same asymmetric refinement — with every disputed decision
exposed as a field of `Params`, so a candidate rule can be scored against the whole
corpus in a second.

**It is a mirror, not a second source of truth.** `check_mirror.py` asserts that with
`Params()` (the committed defaults) it reproduces the Rust binary's own `.seg` columns on
all 982 corpus rows **at each of the three captured floors** (`--seglength` 3, 5 and 10 —
2 946 rows in all) and the Rust binary's `MaxIBD2` on all 158 non-zero rows. Any
divergence there is a bug in this file, never a discovery.

The three floors are load-bearing, not thoroughness for its own sake: the run merge
(`Params.merge`) cannot fire at the default floor on the corpus's spacings, so a
default-only check is blind to it — and was, for the interval when `Scan` had the merge
and this file did not.

The `.seg` IBD2 caller is `SegScan.ibd2_17` — the rule of
`docs/research/17-seg-caller.md` §7 with §14's corrected bridge, gate window and `inf2`,
and `docs/research/19-ibd2seg-residual.md` §1-§4's fringe, committed to `Scan::ibd2`.
Every clause a later write-up replaced is still reachable as a knob (`bridge_rule="17"`,
`gate_end="right"`, `inf2_ibs1b=True`, `ibd2_fringe="extend"`), and the geometry all of
them replaced is `RETIRED` (`seg_rule="word"`), because the sweeps in `sweep2.py`,
`segtry.py` and `rules*.py` were scored against them and their numbers are quoted in the
research write-ups; passing those reproduces them. Four named bundles are pinned, scored
here as `exact / IBD1Seg / IBD2Seg` of 982 rows at the **default 3 Mb floor**:

    Params()   BASE       982 exact / 982 IBD1Seg / 982 IBD2Seg, MAE 0.000017  (20-…)
    PROP19                806        / 982        / 982        , MAE 0.000023  (19-…)
    FRINGE18              747        / 982        / 896        , MAE 0.000067  (18-…)
    RETIRED               705        / 822        / 822        , MAE 0.001376  (17-…)

`BASE` and `PROP19` differ only in how `PropIBD` is printed, which is why their two
estimate columns agree and their `exact` counts do not. All three retired bundles pin
`merge=False`; away from the default floor that is the whole difference between `PROP19`
and `BASE` (at 5 Mb: 910/946 against 959/947 — see `scorecard.py` for the same numbers
measured from the binary rather than from here).

Two rulers are implemented over one caller, exactly as `analysis/segments.rs` describes:

* `.seg` measures a call from its refined `lo` to its refined `hi`;
* `--ibs`'s `MaxIBD2` measures the same call **word-aligned**, `pos[64e+63] - pos[64u]`,
  recovering `(u, e)` from the refined endpoints by integer division.
"""

from dataclasses import dataclass

import numpy as np

import kingdata as kd

WORD = 64
PC = np.bitwise_count

SEGLEN = 3_000_000

#: Markers a `.seg` IBD2 call reaches past the nearest bounding mismatch
#: (`docs/research/17-seg-caller.md` §5) — `IBD2_REACH` in the Rust engine.
REACH = 63
LONG = 10_000_000

# --- the `--seglength` run merge (`docs/research/20-seglength-floor.md`) ----------
# Mirrors `MERGE_MAX_WORDS` / `MERGE_FREE` / `MERGE_COST1` / `MERGE_COST2` in
# `crates/king-core/src/ibdseg.rs`. Every one is bisected on `mergelab.py` canvases
# against the reference binary, never fitted to the corpus.
MERGE_WORDS = 2     # unusable words a merge may bridge
MERGE_FREE = 2      # bad markers a merge gets for nothing
MERGE_COST1 = 4     # informative markers one further bad marker costs, IBD1 pass
MERGE_COST2 = 3     # ...and IBD2


@dataclass(frozen=True)
class Params:
    """Every knob the geometry investigation wants to turn.

    Defaults are the committed engine. Anything not `None`/default here changes the
    rule, so a scorecard line always names what it changed.
    """

    # --- word predicates -------------------------------------------------
    seg_rule: str = "17"            # "17" = the committed `.seg` caller
    #                                 (`docs/research/17-seg-caller.md`); "word" = the
    #                                 retired word-aligned geometry, kept so the sweeps
    #                                 scored against it still run — see `RETIRED`
    ibd2_dirty_ibs1: int = 2        # IBS1 count that makes a word too dirty for IBD2
    bridge: bool = True             # a lone dirty word between clean ones is absorbed
    #                                 (`seg_rule="word"` only — the `.seg` caller's own
    #                                 bridge is `bridge_rule` below)
    gate: int = 10                  # MIN_INFORMATIVE
    # How a lone unusable word between usable ones is absorbed by the `.seg` caller.
    # "19" is the committed rule (`docs/research/17-seg-caller.md` §14): the ordinary
    # gate, asked twice — once over the run so far, from its gate-start through this
    # word, and once over the continuation, from the very next word. "17" is the fitted
    # lookahead it replaced (§7), kept so the sweeps scored against it still run.
    bridge_rule: str = "19"
    # Where the gate's window stops on the right. "next" is the committed rule — the one
    # word the reach touches (§14.2); "right" is §7's, which followed the reach into the
    # word after next when the mismatch sat late in its word.
    gate_end: str = "next"
    # `inf2`, the gate's statistic. False is the committed rule — HetHet + A1A1/A1A1
    # (§14.3); True is the retired `p1 & p1`, which also counts het-vs-A1A1.
    inf2_ibs1b: bool = False
    min_run1: int = 1
    min_run2: int = 1

    # --- IBD2 geometry ---------------------------------------------------
    ibd2_tail: int = 2              # `v + tail >= w1` snaps the end to the segment end
    ibd2_ext: int = 1               # words the run reaches past its last clean word
    ibd2_ext_last: bool = True      # ...also when the extending word would be `w1`
    ibd2_start_refine: bool = False  # refine the start by the previous word's last IBS0
    ibd2_geom: str = "word"         # "word" = the committed geometry; "ibd1" = borrow
    #                                 Scan.left_end / Scan.right_end wholesale
    ibs_pad: int = 0                # words `--ibs`'s ruler adds past the call's own end

    # --- ordering --------------------------------------------------------
    clip_before_len: bool = True    # clip against the previous call before the length test
    ibd1_clip_ibd2: bool = False    # IBD1 calls are clipped off IBD2 territory, not
    #                                 just subtracted from the total
    # How an IBD1 call and the IBD2 calls inside it make `IBD1Seg`
    # (`docs/research/18-ibd1-caller.md` §6). "pieces" is the committed rule: cut at
    # marker granularity, excluding the IBD2 call's own end markers, and apply the
    # `--seglength` floor to each piece. "overlap" is the retired `length - overlap`.
    ibd1_sub: str = "pieces"
    # Which IBD2 calls the IBD1 pass is cut by: only those that survived `--seglength`
    # ("kept", committed) or every call the IBD2 caller made ("all"). Identical at the
    # default floor, where nothing is dropped.
    ibd1_cut: str = "kept"

    # --- the `--seglength` run merge -------------------------------------
    # `docs/research/20-seglength-floor.md`, committed as `Scan::join_runs`/`merge_ok`.
    # Two runs of the same pass, **after the gate has refused what it refuses**, are
    # joined iff at most `merge_words` unusable words lie between them, the gap from the
    # earlier run's last marker to the later run's first is **strictly** under
    # `--seglength`, and `cost * (bad - merge_free) <= X` over those words. The merged
    # run then takes the gate, the endpoints and the floor exactly as an unmerged one.
    #
    # False is the engine as it stood after `19-…` and before `20-…`; the three retired
    # bundles below pin it so the numbers those write-ups quote still reproduce. It
    # cannot fire at the default floor on the corpus's spacings, which is why turning it
    # off moves nothing at 3 Mb and a great deal at 5 and 10.
    merge: bool = True
    merge_words: int = MERGE_WORDS
    merge_free: int = MERGE_FREE
    merge_cost1: int = MERGE_COST1
    merge_cost2: int = MERGE_COST2
    # The het-vs-A1A1 count at which the IBD1 budget's `X` switches from the A1A1/A1A1
    # markers to the het-vs-A1A1 ones. Bisected at `MIN_INFORMATIVE`: against A1A1/A1A1
    # loads of 16, 24, 30 and 40, nine het-vs-A1A1 markers join and ten do not.
    merge_gate: int = 10
    # --- `docs/research/21-push-merge.md` --------------------------------
    # The IBD2 pass's own merge, re-measured. `merge21=False` is the engine as it stood
    # after `20-…`: a two-word cap shared with the IBD1 pass, the interruption taken
    # between the two *runs*, and `X = inf2`. True is the committed rule: no cap, the
    # interruption taken between the two **gate windows**, and `X = HetHet` with the same
    # switch at `merge_gate` the IBD1 pass has (§3-§5).
    merge21: bool = True
    # The one-word push (`17-…` §6) is conditional: a call arms it only when it reaches
    # `min_bp // push_fraction` measured from its own gate-start word (§2). `None` is the
    # retired unconditional form — every call after the first is pushed.
    push_fraction: int = 2
    # --- `docs/research/23-gap-bound.md` ---------------------------------
    # The floor is asked twice, and the second question is about the run's **gate
    # window** rather than the reported call: a call whose window spans under
    # `min_bp // window_fraction` is dropped however long the call itself measures.
    # The IBD2 pass keeps a window of exactly that span and the IBD1 pass does not, so
    # the IBD1 comparison is the strict one (§2, §4). `None` retires the clause.
    window_fraction: int = 2
    # Which words between two runs the IBD1 merge's budget is summed over: "all" is the
    # committed rule (a gate-refused run stepped over by the cap is still paid for by
    # the budget, `23-…` §5), "unusable" is `20-…`'s reading.
    merge_span: str = "all"

    # --- marker-level boundary refinement --------------------------------
    # Where inside the flanking word a call stops. `last`/`first` name which IBS0 of that
    # word is used and the integer is the offset added to it, so `("last", 0)` is "end on
    # the flanking word's last IBS0" — the committed rule.
    ibd1_right: tuple = ("last", 0)
    ibd1_left: tuple = ("last", 1)
    ibd2_right: tuple = ("last", 0)
    # The fringe rules, for a run that reaches the usable segment's own first/last word.
    # These two are the **IBD1** pass's fringe (they read the IBS0 masks); the IBD2
    # pass's own fringe is `ibd2_fringe` below and reads the mismatch masks.
    fringe_right: tuple = ("first", -1)
    fringe_left: tuple = ("last", 1)
    # The IBD2 fringe (`docs/research/19-ibd2seg-residual.md` §1-§4). Once a computed end
    # lands on the word grid's own edge the word scan is over and a **marker** scan takes
    # over across the partial word beyond it.
    #   "mis"    — the committed rule: run out to the segment's own first/last marker, or
    #              only as far as the nearest het-vs-hom **mismatch** there, stopping one
    #              marker short of it (the last such marker on the left, the first on the
    #              right). An opposite homozygote in a fringe does not stop the call.
    #   "extend" — the retired `17-…` §5 clause: a run whose own first/last *word* is the
    #              segment's runs unconditionally to the segment's first/last marker, and
    #              it is applied after the push rather than before. This is exactly the
    #              engine of `docs/research/18-ibd1-caller.md` (747 exact / 896 `IBD2Seg`
    #              at 3 Mb), kept so those numbers still reproduce.
    #   "none"   — never leave the word grid.
    ibd2_fringe: str = "mis"
    ibd2_fringe_off: int = 1   # markers between the call's end and that mismatch

    # --- how `<prefix>.seg` prints PropIBD ---------------------------------
    # "printed"   — the committed rule (`king_core::ibdseg::seg_prop_ibd`): the two
    #               columns are rounded to the four decimals the file shows, then
    #               combined as `i2*1e-4 + i1*5e-5` and printed. Exact on all 4 172
    #               reference `.seg` rows.
    # "unrounded" — `IBD2Seg + IBD1Seg/2` at full precision, which is what the `.kin`
    #               family prints and what `.seg` was assumed to print through
    #               `docs/research/17-` and `18-`. Retired for `.seg`; keep it to
    #               reproduce those write-ups' scorecards.
    seg_prop: str = "printed"

    def label(self):
        d = Params()
        bits = [f"{k}={getattr(self, k)!r}" for k in self.__dataclass_fields__
                if getattr(self, k) != getattr(d, k)]
        return "baseline" if not bits else " ".join(bits)


BASE = Params()

#: The geometry `BASE` replaced: word-aligned ends, a five-mismatch word predicate, the
#: two-word tail snap and `length - overlap`. Scores 705 exact `.seg` rows against
#: `BASE`'s 806, MAE 0.001376 against 0.000023. Kept so the sweeps scored against it
#: still reproduce.
#: All three retired bundles pin `merge=False`: each is "the engine as it stood at
#: write-up N", and the run merge did not land until `20-…`.
RETIRED = Params(seg_rule="word", ibd2_dirty_ibs1=5, ibd1_sub="overlap",
                 seg_prop="unrounded", merge=False,
                 merge21=False, push_fraction=None,
                 window_fraction=None, merge_span="unusable")

#: The engine of `docs/research/18-ibd1-caller.md`: everything `BASE` has except the
#: IBD2 fringe of `19-…` (still `17-…` §5's unconditional "extend") and `.seg`'s own
#: `PropIBD` rule (still the `.kin` one). Scores 747 exact rows / 896 exact `IBD2Seg` /
#: MAE 0.000067 at 3 Mb. Kept so `18-…`'s numbers reproduce from this file too.
FRINGE18 = Params(ibd2_fringe="extend", seg_prop="unrounded", merge=False,
                 merge21=False, push_fraction=None,
                 window_fraction=None, merge_span="unusable")

#: `BASE`'s caller with the retired **`.kin`** `PropIBD` rule on `.seg` — the tree as it
#: stood after `19-…` and before `20-…`. 806 exact rows against `BASE`'s 982.
PROP19 = Params(seg_prop="unrounded", merge=False,
                 merge21=False, push_fraction=None,
                 window_fraction=None, merge_span="unusable")


def seg_prop_ibd(ibd1_seg, ibd2_seg):
    """`PropIBD` as `<prefix>.seg` prints it — `king_core::ibdseg::seg_prop_ibd`.

    The two columns are first rounded to the four decimals the file shows; the printed
    value is then `i2*1e-4 + i1*5e-5` on those integers. The 1 313 of 4 172 reference rows
    where that lands on an exact decimal half are resolved by which side of the half the
    double falls on, and this expression agrees with the reference on every one of them —
    `(i1+2*i2)/20000`, `i2/10000 + i1/20000` and integer round-half-up do not.
    """
    i1 = int(round(float("%.4f" % ibd1_seg) * 10000))
    i2 = int(round(float("%.4f" % ibd2_seg) * 10000))
    return i2 * 1e-4 + i1 * 5e-5


# ---------------------------------------------------------------------------
# per-pair word masks
# ---------------------------------------------------------------------------

_MASKS = {}


def masks(ds, i, j):
    """(ibs0, ibs1, inf1, inf2) word masks plus the popcount vectors.

    `nhh` is the per-word HetHet popcount, which only `--ibs`'s caller reads
    (`Scan::ibd2_words`); the `.seg` caller never looks at it.
    """
    key = (ds.name, i, j)
    v = _MASKS.get(key)
    if v is None:
        p0i, p1i = ds.p0[i], ds.p1[i]
        p0j, p1j = ds.p0[j], ds.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        ibs0 = p0i & p0j & (p1i ^ p1j)
        ibs1 = (het_i & p0j) | (p0i & het_j)
        share = p1i & p1j
        inf1 = share & (p0i | p0j)
        # The `.seg` gate's `inf2` is `share & ~ibs1` — HetHet + A1A1/A1A1, a het-vs-A1A1
        # marker being uninformative (§14.3). `share` alone is the retired statistic,
        # kept as its own cumulative so `inf2_ibs1b=True` and the older scripts that read
        # `cum2` directly still reproduce.
        v = (ibs0, PC(ibs0).astype(np.int32), PC(ibs1).astype(np.int32),
             np.concatenate(([0], np.cumsum(PC(inf1).astype(np.int64)))),
             np.concatenate(([0], np.cumsum(PC(share).astype(np.int64)))),
             PC(het_i & het_j).astype(np.int32), ibs1,
             np.concatenate(([0], np.cumsum(PC(share & ~ibs1).astype(np.int64)))))
        _MASKS[key] = v
    return v


_MERGE_COUNTS = {}


def merge_counts(ds, i, j):
    """Per-word popcounts the run merge reads — `(Z, U, V, U2, M)`.

    Kept apart from `masks()` because that function's tuple is unpacked positionally by
    `seg19.py`, `endfit.py`, `invert.py` and others; appending to it would break them.

    * `Z` — opposite homozygotes (`ibs0`), the IBD1 pass's bad markers.
    * `U` — A1A1/A1A1 (`inf1 & ~ibs1`), the IBD1 budget's `X`.
    * `V` — het-vs-A1A1 (`inf1 & ibs1`), which replaces `U` once it reaches
      `Params.merge_gate` on its own.
    * `U2` — `inf2 = share & ~ibs1` (HetHet + A1A1/A1A1), the IBD2 budget's `X` and the
      very count the `.seg` gate uses.
    * `M` — het-vs-hom mismatches (`ibs1`), which the IBD2 pass adds to `Z` to get its
      own bad-marker count.
    """
    key = (ds.name, i, j)
    v = _MERGE_COUNTS.get(key)
    if v is None:
        p0i, p1i = ds.p0[i], ds.p1[i]
        p0j, p1j = ds.p0[j], ds.p1[j]
        het_i = ~p0i & p1i
        het_j = ~p0j & p1j
        ibs0 = p0i & p0j & (p1i ^ p1j)
        ibs1 = (het_i & p0j) | (p0i & het_j)
        share = p1i & p1j
        inf1 = share & (p0i | p0j)
        v = (PC(ibs0).astype(np.int32), PC(inf1 & ~ibs1).astype(np.int32),
             PC(inf1 & ibs1).astype(np.int32), PC(share & ~ibs1).astype(np.int32),
             PC(ibs1).astype(np.int32))
        _MERGE_COUNTS[key] = v
    return v


def _last_bit(m):
    """Index of the highest set bit of a 64-bit mask (mask must be non-zero)."""
    return int(m).bit_length() - 1


def _first_bit(m):
    return (int(m) & -int(m)).bit_length() - 1


# ---------------------------------------------------------------------------
# the caller, one usable segment at a time
# ---------------------------------------------------------------------------

def _runs(ok):
    """Maximal runs of True in a bool array, as (start, stop) inclusive index pairs."""
    d = np.diff(np.concatenate(([False], ok, [False])).astype(np.int8))
    return list(zip(np.flatnonzero(d == 1).tolist(),
                    (np.flatnonzero(d == -1) - 1).tolist()))


class SegScan:
    """One pair over one usable segment — the Rust `Scan`, with `Params` applied."""

    def __init__(self, ds, i, j, seg, p=BASE):
        self.ds, self.p = ds, p
        chrom, lo, hi = seg
        self.lo, self.hi = lo, hi
        self.w0 = -(-lo // WORD)
        self.w1 = (hi + 1) // WORD - 1
        self.n = max(0, self.w1 - self.w0 + 1)
        ibs0, n0, n1, k1, k2, nhh, ibs1, k2s = masks(ds, i, j)
        self.ibs0 = ibs0
        self.ibs1 = ibs1
        self.n0 = n0
        self.n1 = n1
        self.nhh = nhh
        self.cum1, self.cum2 = k1, k2
        # the `.seg` gate's own cumulative — see `Params.inf2_ibs1b`
        self.cum2s = k2 if p.inf2_ibs1b else k2s
        # per-word counts for the run merge (`Params.merge`)
        self.mz, self.mu, self.mv, self.mu2, self.mm = merge_counts(ds, i, j)
        # Fringe masks: the segment's own markers in the two words it does not own.
        # The IBD1 pass reads the IBS0 ones, the IBD2 pass the het-vs-hom mismatch ones —
        # each pass stops at its own breaking marker (`19-…` §2, §5).
        self.head = 0
        self.tail = 0
        self.head_mis = 0
        self.tail_mis = 0
        if self.n > 0:
            if lo != WORD * self.w0:
                keep = lo - WORD * (self.w0 - 1)
                drop = ~((1 << keep) - 1)
                self.head = int(ibs0[self.w0 - 1]) & drop
                self.head_mis = int(ibs1[self.w0 - 1]) & drop
            if hi != WORD * (self.w1 + 1) - 1:
                keep = hi - WORD * (self.w1 + 1) + 1
                self.tail = int(ibs0[self.w1 + 1]) & ((1 << keep) - 1)
                self.tail_mis = int(ibs1[self.w1 + 1]) & ((1 << keep) - 1)

    # --- the IBD2 fringe stops ------------------------------------------
    def head_stop(self):
        """How far left an IBD2 call may creep past the word grid — `Scan::head_stop`."""
        if self.p.ibd2_fringe != "mis" or not self.head_mis:
            return self.lo
        last = WORD * (self.w0 - 1) + _last_bit(self.head_mis)
        return max(last + self.p.ibd2_fringe_off, self.lo)

    def tail_stop(self):
        """The mirror — `Scan::tail_stop`."""
        if self.p.ibd2_fringe != "mis" or not self.tail_mis:
            return self.hi
        first = WORD * (self.w1 + 1) + _first_bit(self.tail_mis)
        return min(first - self.p.ibd2_fringe_off, self.hi)

    # --- gates ---------------------------------------------------------
    def informative(self, cum, u, v):
        return int(cum[v + 1] - cum[u]) >= self.p.gate

    # --- the `--seglength` run merge ------------------------------------
    def merge_ok(self, mid, pass2):
        """Whether two runs may be joined across the **unusable** words `mid`.

        The mirror of `Scan::merge_ok`. `mid` holds global word indices; a gate-refused
        run's words are not in it — it is stepped over, not counted.
        """
        p = self.p
        bad = int(sum(int(self.mz[k]) for k in mid))
        if pass2:
            bad += int(sum(int(self.mm[k]) for k in mid))
            cost = p.merge_cost2
            if p.merge21:
                # `21-…` §5: `X` is the HetHet count, with the IBD1 pass's own switch.
                # `mu` is the A1A1/A1A1 half of `inf2`; the rest of `mu2` is HetHet.
                x = int(sum(int(self.mu[k]) for k in mid))
                v = int(sum(int(self.mu2[k]) - int(self.mu[k]) for k in mid))
                if v >= p.merge_gate:
                    x = v
            else:
                x = int(sum(int(self.mu2[k]) for k in mid))
        else:
            x = int(sum(int(self.mu[k]) for k in mid))
            v = int(sum(int(self.mv[k]) for k in mid))
            if v >= p.merge_gate:
                x = v
            cost = p.merge_cost1
        return cost * max(0, bad - p.merge_free) <= x

    def join_runs2(self, runs, ok, mis, pos, min_bp):
        """The IBD2 merge of `docs/research/21-push-merge.md` §3-§4 — `Scan::join_runs2`.

        No cap on the interruption's width, and what separates two runs is the space
        between their **gate windows**: the earlier ends at `ge_of(b)`, the later opens
        at its gate-start word. A word covered by any window — which is any unusable,
        IBS0-free word right after a usable one — is not part of the interruption.
        """
        n = self.n

        def ge_of(b):
            return b + 1 if (b + 1 < n and int(self.n0[self.w0 + b + 1]) == 0
                             and mis[b + 1]) else b

        out = []
        for a, b in runs:
            if out:
                pa, pb = out[-1]
                q = ge_of(pb)
                g2 = next((t for t in range(a, b + 1) if mis[t] == 0), a)
                mid = [k for k in range(q + 1, a)
                       if not ok[k] and not (k > 0 and ok[k - 1]
                                             and int(self.n0[self.w0 + k]) == 0)]
                gap = int(pos[WORD * (self.w0 + g2)]
                          - pos[WORD * (self.w0 + q + 1) - 1])
                if (any(not ok[k] for k in range(pb + 1, a)) and gap < min_bp
                        and self.merge_ok([self.w0 + k for k in mid], True)):
                    out[-1] = (pa, b)
                    continue
            out.append((a, b))
        return out

    def join_runs(self, runs, usable, pos, min_bp, pass2):
        """Join adjacent gate-passing runs across a short interruption.

        The mirror of `Scan::join_runs`. `runs` are global word-index pairs that already
        cleared the gate, in order; `usable` is indexed by scan word (0-based within the
        segment). Two runs join iff at most `merge_words` unusable words lie between
        them, the run-to-run gap is **strictly** under `--seglength`, and `merge_ok`
        passes over those words.
        """
        p = self.p
        out = []
        for u, v in runs:
            if out:
                pu, pv = out[-1]
                bad_words = [k for k in range(pv + 1, u) if not usable[k - self.w0]]
                # The cap counts only the unusable words, so a gate-refused run between
                # the two is stepped over (`20-…` §6); the budget is summed over every
                # word in the interruption, that run included (`23-…` §5).
                mid = (list(range(pv + 1, u)) if p.merge_span == "all" and not pass2
                       else bad_words)
                # `pos[first marker of the later run] - pos[last marker of the earlier]`,
                # strictly under the floor. `WORD*(pv+1) - 1` is the earlier run's last
                # marker; `WORD*u` the later run's first.
                if (bad_words and len(bad_words) <= p.merge_words
                        and int(pos[WORD * u] - pos[WORD * (pv + 1) - 1]) < min_bp
                        and self.merge_ok(mid, pass2)):
                    out[-1] = (pu, v)
                    continue
            out.append((u, v))
        return out

    # --- IBD1 ----------------------------------------------------------
    @staticmethod
    def _pick(mask, rule):
        """Bit of `mask` named by `rule = (which, offset)`; `mask` must be non-zero."""
        which, off = rule
        return (_last_bit(mask) if which == "last" else _first_bit(mask)) + off

    def right_end(self, v):
        """Right end of an IBD1 run whose last good word is `v` (global word index)."""
        if v + 1 <= self.w1:
            m = int(self.ibs0[v + 1])
            if m == 0:
                return min(WORD * (v + 1) + 63, self.hi)
            return WORD * (v + 1) + self._pick(m, self.p.ibd1_right)
        if self.tail:
            return WORD * (self.w1 + 1) + self._pick(self.tail, self.p.fringe_right)
        return self.hi

    def left_end(self, u):
        if u - 1 >= self.w0:
            m = int(self.ibs0[u - 1])
            if m == 0:
                return WORD * u
            return WORD * (u - 1) + self._pick(m, self.p.ibd1_left)
        if self.head:
            return WORD * (self.w0 - 1) + self._pick(self.head, self.p.fringe_left)
        return self.lo

    def ibd1(self, pos, min_bp):
        p = self.p
        if self.n == 0:
            return []
        ok = self.n0[self.w0:self.w1 + 1] == 0
        # The gate is asked **first** (`20-…` §6): a run under `MIN_INFORMATIVE` is
        # refused outright, and then lies inside a later interruption rather than ending
        # one — it can never be the endpoint of a merged segment, but does not stop one.
        kept = []
        for a, b in _runs(ok):
            if b - a + 1 < p.min_run1:
                continue
            u, v = self.w0 + a, self.w0 + b
            if not self.informative(self.cum1, u, v):
                continue
            kept.append((u, v))
        if p.merge:
            kept = self.join_runs(kept, ok, pos, min_bp, False)
        out = []
        for u, v in kept:
            # `23-…` §4: the run's own complete words must span more than half the floor.
            # Strict, one unit of `min_bp // 2` tighter than the IBD2 pass's test.
            if p.window_fraction is not None and not (
                    int(pos[WORD * v + WORD - 1] - pos[WORD * u])
                    > min_bp // p.window_fraction):
                continue
            hi = self.right_end(v)
            lo = self.left_end(u)
            out = self._emit(out, lo, hi, pos, min_bp)
        return out

    # --- IBD2 ----------------------------------------------------------
    def ibd2(self, pos, min_bp):
        if self.p.seg_rule == "17":
            return self.ibd2_17(pos, min_bp)
        return self.ibd2_word(pos, min_bp)

    def ibd2_17(self, pos, min_bp):
        """The committed `.seg` IBD2 caller — the mirror of `Scan::ibd2`.

        `docs/research/17-seg-caller.md` §7, with §14's corrected bridge. A word is usable
        iff it carries no opposite homozygote and at most one het-vs-hom mismatch; a lone
        unusable word is absorbed iff the gate passes on both sides of it; the gate is
        `inf2` counted from the run's first mismatch-free word through the one word the
        reach touches; the endpoints reach `REACH` markers past the nearest mismatch and
        are blocked whole-word by an IBS0; every call after the first in a usable segment
        starts one word past its gate-start word; and the segment's own fringes take a
        touching call to its first or last marker.
        """
        p = self.p
        n = self.n
        if n == 0:
            return []
        w0, w1 = self.w0, self.w1
        cum = self.cum2s
        z = [int(self.n0[w0 + k]) != 0 for k in range(n)]
        mis = [int(self.n1[w0 + k]) for k in range(n)]
        usable = [(not z[k]) and mis[k] < p.ibd2_dirty_ibs1 for k in range(n)]
        # The two marker-level stops beyond the word grid (`19-…` §3).
        head_stop, tail_stop = self.head_stop(), self.tail_stop()

        def ge_of(b):
            """The one word a run ending at `b` reaches into — whole-word."""
            return b + 1 if (b + 1 < n and not z[b + 1] and mis[b + 1]) else b

        def gate_ok(g, b):
            return int(cum[w0 + ge_of(b) + 1] - cum[w0 + g]) >= p.gate

        ok = list(usable)
        if p.bridge_rule == "19":
            # The bridge is the gate asked twice: over the run so far, from its
            # gate-start through this word, and over the continuation, from the very
            # next word — which must therefore be mismatch-free — through the words its
            # own right end reaches.
            gs0 = None
            for k in range(n):
                if usable[k]:
                    if gs0 is None and mis[k] == 0:
                        gs0 = k
                    continue
                bridged = False
                if (gs0 is not None and k > 0 and not z[k] and k + 1 < n
                        and usable[k + 1] and mis[k + 1] == 0):
                    b2 = k + 1
                    while b2 + 1 < n and usable[b2 + 1]:
                        b2 += 1
                    bridged = gate_ok(gs0, k - 1) and gate_ok(k + 1, b2)
                if bridged:
                    ok[k] = True
                else:
                    gs0 = None
        else:
            for k in range(1, max(0, n - 1)):
                if usable[k] or z[k] or not usable[k - 1] or not usable[k + 1]:
                    continue
                if mis[k + 1] != 0:
                    continue
                acc, t = 0, k + 1
                while t < n and usable[t] and acc < p.gate:
                    acc += int(cum[w0 + t + 1] - cum[w0 + t])
                    t += 1
                if acc >= p.gate:
                    ok[k] = True

        # The gate is asked **first**, so a run it refuses can sit inside a later merge's
        # interruption without ending it (`20-…` §6); the emit loop below then re-asks it
        # on the merged run, exactly as `Scan::ibd2` does. The pre-filter runs only when
        # the merge does: it uses `gate_ok(gs, b)`, which is `gate_end="next"`, so
        # applying it unconditionally would change the retired `gate_end="right"` variant
        # — whose own window is not known until `right` has been computed below.
        kept = []
        for a, b in _runs(ok):
            if b - a + 1 < p.min_run2:
                continue
            if p.merge:
                g = next((t for t in range(a, b + 1) if mis[t] == 0), None)
                if g is None or not gate_ok(g, b):
                    continue
            kept.append((w0 + a, w0 + b))
        if p.merge:
            if p.merge21:
                kept = self.join_runs2([(u - w0, v - w0) for u, v in kept],
                                       ok, mis, pos, min_bp)
                kept = [(w0 + a, w0 + b) for a, b in kept]
            else:
                kept = self.join_runs(kept, ok, pos, min_bp, True)

        out, armed = [], False
        for u, v in kept:
            a, b = u - w0, v - w0
            left = WORD * u
            if a > 0 and not z[a - 1] and int(self.ibs1[u - 1]):
                last = WORD * (u - 1) + _last_bit(int(self.ibs1[u - 1]))
                left = max(0, last - REACH)
                if a < 2 or z[a - 2]:
                    left = max(left, WORD * (u - 1))
            # Once a computed end lands on the grid's own edge the word scan is over and
            # the marker scan takes over, whether that moves the end out (the fringe
            # carries no mismatch) or pulls it back in (`19-…` §4). It is the computed
            # end that snaps, not the run: `left <= WORD*w0`, not `u == w0`.
            if p.ibd2_fringe == "mis" and left <= WORD * w0:
                left = head_stop
            right = WORD * v + WORD - 1
            if b + 1 < n and not z[b + 1] and int(self.ibs1[v + 1]):
                right = WORD * (v + 1) + _first_bit(int(self.ibs1[v + 1])) + REACH
                if b + 2 >= n or z[b + 2]:
                    right = min(right, WORD * (v + 2) - 1)
            if p.ibd2_fringe == "mis" and right >= WORD * (w1 + 1) - 1:
                right = tail_stop
            gs = next((t for t in range(a, b + 1) if mis[t] == 0), None)
            if gs is None:
                continue
            if p.gate_end == "next":
                if not gate_ok(gs, b):
                    continue
                # `23-…` §2: the same window, now for its length. Asked here, after the
                # merge, so a run the bound refuses still merges with its neighbour.
                if p.window_fraction is not None and (
                        int(pos[WORD * (w0 + ge_of(b)) + WORD - 1]
                            - pos[WORD * (w0 + gs)]) < min_bp // p.window_fraction):
                    continue
            else:
                ge = b
                if right > WORD * v + WORD - 1:
                    ge = min(n - 1, right // WORD - w0)
                if int(cum[w0 + ge + 1] - cum[w0 + gs]) < p.gate:
                    continue
            if armed:
                left = max(left, WORD * (w0 + gs + 1))
            if p.ibd2_fringe == "extend":
                # The retired `17-…` §5 clause, in the place it used to occupy.
                if u == w0:
                    left = min(left, self.lo)
                if v == w1:
                    right = max(right, self.hi)
            left, right = max(left, self.lo), min(right, self.hi)
            if out:
                left = max(left, out[-1][1])   # calls may touch, but not overlap
            if left > right:
                continue
            # `21-…` §2: the push is armed by a call reaching half the floor, measured
            # from its own gate-start word. `push_fraction=None` is the retired
            # unconditional form.
            if p.push_fraction is None:
                armed = True
            else:
                g0 = min(max(WORD * (w0 + gs), self.lo), right)
                armed = armed or pos[right] - pos[g0] >= min_bp // p.push_fraction
            if pos[right] - pos[left] >= min_bp:
                out.append((left, right))
        return out

    def ibd2_word(self, pos, min_bp):
        """The **retired** word-aligned geometry — reachable as `RETIRED`, not committed."""
        p = self.p
        if self.n == 0:
            return []
        sl = slice(self.w0, self.w1 + 1)
        clean = (self.n0[sl] == 0) & (self.n1[sl] < p.ibd2_dirty_ibs1)
        ok = clean.copy()
        if p.bridge:
            n0 = self.n0[sl]
            for k in range(1, self.n - 1):
                if not clean[k] and clean[k - 1] and clean[k + 1] and n0[k] == 0:
                    ok[k] = True
        out = []
        for a, b in _runs(ok):
            if b - a + 1 < p.min_run2:
                continue
            u, v = self.w0 + a, self.w0 + b
            if not self.informative(self.cum2, u, v):
                continue
            if p.ibd2_geom == "ibd1":
                lo, hi = self.left_end(u), self.right_end(v)
                out = self._emit(out, lo, hi, pos, min_bp)
                continue
            e = self._ibd2_end_word(v)
            lo = WORD * u if u != self.w0 else self.lo
            if p.ibd2_start_refine and u != self.w0:
                m = int(self.ibs0[u - 1])
                if m:
                    lo = WORD * (u - 1) + _last_bit(m) + 1
            if e == self.w1:
                hi = self.hi
            else:
                m = int(self.ibs0[e])
                hi = WORD * e + (63 if m == 0 else self._pick(m, self.p.ibd2_right))
            out = self._emit(out, lo, hi, pos, min_bp)
        return out

    def _ibd2_end_word(self, v):
        p = self.p
        if v + p.ibd2_tail >= self.w1:
            return self.w1
        e = v + p.ibd2_ext
        if not p.ibd2_ext_last and e >= self.w1:
            e = self.w1 - 1 if self.w1 - 1 >= v else v
        return min(e, self.w1)

    # --- shared emit ---------------------------------------------------
    def _emit(self, out, lo, hi, pos, min_bp):
        p = self.p
        if p.clip_before_len:
            if out:
                lo = max(lo, out[-1][1] + 1)
            if lo <= hi and pos[hi] - pos[lo] >= min_bp:
                out.append((lo, hi))
        else:
            if lo <= hi and pos[hi] - pos[lo] >= min_bp:
                if out:
                    lo = max(lo, out[-1][1] + 1)
                if lo <= hi:
                    out.append((lo, hi))
        return out


# ---------------------------------------------------------------------------
# pair aggregation — the two rulers
# ---------------------------------------------------------------------------

def call_pair(ds, i, j, p=BASE, min_bp=SEGLEN):
    """Returns (ibd1_bp, ibd2_bp, longest_bp, max_ibd2_wordaligned)."""
    pos = ds.pos
    ibd1_bp = ibd2_bp = longest = 0
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        c2 = sc.ibd2(pos, min_bp)
        c1 = sc.ibd1(pos, min_bp)
        # What the IBD1 pass is cut by: the IBD2 calls that survived `--seglength`
        # ("kept", committed) or every IBD2 call the caller made, floor or no floor
        # ("all"). They are the same set at the default floor.
        cut = c2 if p.ibd1_cut == "kept" else sc.ibd2(pos, 0)
        for lo, hi in c2:
            ln = int(pos[hi] - pos[lo])
            ibd2_bp += ln
            longest = max(longest, ln)
        for lo, hi in c1:
            ln = int(pos[hi] - pos[lo])
            longest = max(longest, ln)
            if p.ibd1_sub == "pieces":
                ibd1_bp += sum(v for v in
                               (int(pos[b] - pos[a]) for a, b in _pieces((lo, hi), cut))
                               if v >= min_bp)
            else:
                ibd1_bp += ln - _overlap((lo, hi), cut, pos)
    return ibd1_bp, ibd2_bp, longest, max_ibd2(ds, i, j, p)


def _pieces(c, others):
    """`c` with the IBD2 calls cut out — `docs/research/18-ibd1-caller.md` §6.1.

    The cut excludes the IBD2 call's own end markers, so `[lo, hi]` cut by `[a, b]`
    leaves `[lo, a-1]` and `[b+1, hi]`; each piece then faces `--seglength` on its own.
    """
    lo, hi = c
    out, cur = [], lo
    for a, b in sorted(others):
        if b < lo or a > hi:
            continue
        if a > cur:
            out.append((cur, a - 1))
        cur = max(cur, b + 1)
    if cur <= hi:
        out.append((cur, hi))
    return out


def _overlap(c, others, pos):
    tot = 0
    for lo, hi in others:
        a, b = max(c[0], lo), min(c[1], hi)
        if a < b:
            tot += int(pos[b] - pos[a])
    return tot


IBS_IBD2_DIRTY = 5       # het-vs-hom mismatches that make a word break an --ibs run
IBS_IBD2_HETHET = 95     # HetHet markers one confirmation chunk must carry
IBS_IBD2_CHUNK_MIS = 5   # mismatches that close a chunk
IBS_IBD2_CHUNK_WORDS = 3  # words a chunk must span to be confirmable
IBS_IBD2_EXT_MIS = 1     # mismatches the measured interval may reach past the confirmation
IBS_IBD2_MIN_WORDS = 3   # words a call's measured interval must span


def _chunk_extend(n1, conf, b, ext_mis=IBS_IBD2_EXT_MIS):
    """How far past the confirmed end `conf` the measured interval reaches: on through
    the run's words until the mismatches picked up would exceed `ext_mis`."""
    cum, e = 0, conf
    for k in range(conf + 1, b + 1):
        cum += int(n1[k])
        if cum > ext_mis:
            break
        e = k
    return e


def ibd2_words(sc, restart="fit"):
    """`--ibs`'s own IBD2 caller — the mirror of `Scan::ibd2_words`.

    NOT the `.seg` caller: opposite homozygotes and missing calls are irrelevant here at
    any density, only het-vs-hom mismatches break a run, and the call runs straight
    through IBS0 words that `Scan::ibd2` would split on.  Returns word intervals
    `(lo, hi)` inclusive, in global word indices.

    The shape is the **quantised confirmation scan** of
    `docs/research/16-segment-extension.md`: a run is confirmed in chunks of five
    het-vs-hom mismatches, each needing 95 HetHet over at least three words, and the
    reported interval stops at the last confirmed chunk plus a one-mismatch overhang.
    Reproduces `MaxIBD2` 158/158 and `Pr_IBD2` 158/158 on the corpus.

    `restart` picks the rule for where a new segment begins after a chunk is refused:
    "fit"   — after the 4th mismatch's word when the refusing word holds only the 5th,
              else after the refusing word (what the constructed fixtures measure, and
              what the Rust engine does);
    "after" — always after the refusing word;
    "at"    — always the refusing word.
    The corpus cannot separate "fit" from "after"; see §7 of the write-up.
    """
    n = sc.n
    if n == 0:
        return []
    w0, w1 = sc.w0, sc.w1
    n1, nhh = sc.n1, sc.nhh
    clean = [int(n1[w0 + k]) < IBS_IBD2_DIRTY for k in range(n)]
    # A lone dirty word between two clean ones is absorbed; two in a row are not. Read
    # from `clean`, never from the running copy, so dirty words cannot chain in.
    ok = list(clean)
    for k in range(1, max(0, n - 1)):
        if not clean[k] and clean[k - 1] and clean[k + 1]:
            ok[k] = True

    raw = []
    for a, b in _runs(ok):
        lo_w, hi_w = w0 + a, w0 + b
        # The scan runs one word past the run's last clean word: that word is what makes
        # the mismatch counter fire, and its HetHet counts towards the chunk it closes.
        scan_last = min(hi_w + 1, w1)
        exempt = hi_w + 1 >= w1        # the run reaches the segment's own last two words
        u = lo_w
        mis = het = 0
        cstart = lo_w
        conf = None
        last_mis = None                # last word holding a mismatch, this chunk
        k = u
        while k <= scan_last:
            m, h = int(n1[k]), int(nhh[k])
            before = mis
            mis += m
            het += h
            if mis >= IBS_IBD2_CHUNK_MIS:
                good = (het >= IBS_IBD2_HETHET
                        and k - cstart + 1 >= IBS_IBD2_CHUNK_WORDS)
                # A chunk closed by the usable segment's own last word is exempt.
                if good or (exempt and k >= scan_last):
                    conf, mis, het, last_mis, cstart = k, 0, 0, None, k + 1
                else:
                    if conf is not None:
                        raw.append((u, _chunk_extend(n1, conf, hi_w)))
                    if restart == "at":
                        u = max(k, u + 1)
                    elif restart == "after":
                        u = k + 1
                    elif (before == IBS_IBD2_CHUNK_MIS - 1 and m == 1
                          and last_mis is not None):
                        u = last_mis + 1
                    else:
                        u = k + 1
                    mis = het = 0
                    conf = None
                    last_mis = None
                    cstart = u
                    k = u
                    continue
            if m:
                last_mis = k
            k += 1
        if exempt:                     # the HetHet test is waived at the segment's end
            raw.append((u, w1))
        elif conf is not None:
            raw.append((u, w1 if hi_w + 2 >= w1 else _chunk_extend(n1, conf, hi_w)))

    out = []
    for u, e in raw:
        if u > e:
            continue
        if out:
            u = max(u, out[-1][1] + 1)
        if u > e or e + 1 - u < IBS_IBD2_MIN_WORDS:
            continue
        out.append((u, e))
    return out


def max_ibd2(ds, i, j, p=BASE):
    """`--ibs`'s `MaxIBD2`: the longest `ibd2_words` call, measured word-aligned."""
    pos = ds.pos
    best = 0
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        for u, e in ibd2_words(sc):
            best = max(best, int(pos[WORD * e + 63] - pos[WORD * u]))
    return best


def pr_ibd2(ds, i, j, p=BASE):
    """`--ibs`'s `Pr_IBD2`: the word-aligned `ibd2_words` total over `D`.

    A second aggregate over the same calls as `MaxIBD2`, which grades only the longest
    member.  The 10 Mb rule gates the **pair**, not the call: if no single call reaches
    `LONG`, `Pr_IBD2` is 0 however much shorter material was called.
    """
    pos = ds.pos
    tot = 0
    best = 0
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        for u, e in ibd2_words(sc):
            ln = int(pos[WORD * e + 63] - pos[WORD * u])
            tot += ln
            best = max(best, ln)
    if best < LONG:
        return 0.0
    return tot / ds.denom


_PRT = None


def pr_targets():
    """Reference `Pr_IBD2` for every pair whose `MaxIBD2` is non-zero."""
    global _PRT
    if _PRT is not None:
        return _PRT
    import os
    out = []
    base = os.path.join(kd.ROOT, "tests", "parity", "work", "ibs")
    for name in kd.DATASETS:
        ds = kd.load(name)
        idx = {(f, i): k for k, (f, i) in enumerate(ds.fam)}
        for ext, wide in ((".ibs", False), (".ibs0", True)):
            path = os.path.join(base, "ibs_%s%s" % (name, ext))
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                head = fh.readline().split()
                if "Pr_IBD2" not in head:
                    continue
                c, cm = head.index("Pr_IBD2"), head.index("MaxIBD2")
                for line in fh:
                    f = line.split()
                    if f[cm] == "-9" or float(f[cm]) <= 0:
                        continue
                    if wide:
                        i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
                    else:
                        i, j = idx[(f[0], f[1])], idx[(f[0], f[2])]
                    out.append((name, min(i, j), max(i, j), f[c]))
    _PRT = out
    return out


def score_pr(p=BASE):
    ok = 0
    err = 0.0
    tg = pr_targets()
    for name, i, j, want in tg:
        ds = kd.load(name)
        g = pr_ibd2(ds, i, j, p)
        if "%.4f" % g == want:
            ok += 1
        err += g - float(want)
    return ok, len(tg), err / len(tg)


def max_ibd2_words(ds, i, j, p=BASE):
    """Same as `max_ibd2` but returning `(u, e, bp)` of the winning call."""
    pos = ds.pos
    best = (None, None, 0)
    for seg in ds.segs:
        sc = SegScan(ds, i, j, seg, p)
        if sc.n == 0:
            continue
        prev = None
        for lo, hi in sc.ibd2(pos, 0):
            e = min(hi // WORD + p.ibs_pad, sc.w1)
            u = max(lo // WORD, sc.w0)
            if prev is not None:
                u = max(u, prev + 1)
            if u > e:
                continue
            prev = e
            ln = int(pos[WORD * e + 63] - pos[WORD * u])
            if ln > best[2]:
                best = (u, e, ln)
    return best


# ---------------------------------------------------------------------------
# scoring
# ---------------------------------------------------------------------------

def score_seg(p=BASE, datasets=None, suffix="__ibdseg", min_bp=SEGLEN, per_ds=False):
    """`.seg` scorecard: exact rows on all four printed columns, plus the pair set."""
    rows = exact = both = ibd1ok = ibd2ok = extra = missing = 0
    err = 0.0
    worst = 0.0
    by = {}
    for name in (datasets or kd.DATASETS):
        ds = kd.load(name)
        ref = ds._read_seg(suffix) if suffix != "__ibdseg" else ds.ref
        d = ds.denom
        e = b = r = 0
        for i, j in ds.pairs():
            i1, i2, lg, _ = call_pair(ds, i, j, p, min_bp)
            got = lg >= LONG
            want = (i, j) in ref
            if not want:
                if got:
                    extra += 1
                continue
            if not got:
                missing += 1
                continue
            rows += 1
            r += 1
            a1, a2, ap, at = ref[(i, j)]
            g1, g2 = i1 / d, i2 / d
            # `.seg`'s own PropIBD rule, unless a retired parameterisation asked for the
            # `.kin` one. `inf_type` still reads the full-precision value, as the engine
            # does — this is a printing rule, not a decision.
            gp = seg_prop_ibd(g1, g2) if p.seg_prop == "printed" else g2 + g1 / 2
            ok1 = kd.fmt4(g1) == a1
            ok2 = kd.fmt4(g2) == a2
            ibd1ok += ok1
            ibd2ok += ok2
            if ok1 and ok2:
                both += 1
                b += 1
            if ok1 and ok2 and kd.fmt4(gp) == ap \
                    and kd.inf_type(g1, g2, g2 + g1 / 2) == at:
                exact += 1
                e += 1
            err += abs(gp - ap)
            worst = max(worst, abs(gp - ap))
        by[name] = (e, b, r)
    out = dict(rows=rows, exact=exact, both=both, ibd1=ibd1ok, ibd2=ibd2ok,
               extra=extra, missing=missing,
               mae=err / rows if rows else 0.0, worst=worst)
    if per_ds:
        out["by"] = by
    return out


def max_targets(nonzero=True):
    """Reference `MaxIBD2` values: [(dataset, i, j, bp)] over `.ibs` **and** `.ibs0`.

    `.ibs` is the within-family table (`FID ID1 ID2`) and `.ibs0` the between-family one
    (`FID1 ID1 FID2 ID2`); both carry the column and together they cover every pair the
    reference grades.
    """
    import os
    out = []
    base = os.path.join(kd.ROOT, "tests", "parity", "work", "ibs")
    for name in kd.DATASETS:
        ds = kd.load(name)
        idx = {(f, i): k for k, (f, i) in enumerate(ds.fam)}
        for ext, wide in ((".ibs", False), (".ibs0", True)):
            path = os.path.join(base, "ibs_%s%s" % (name, ext))
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                head = fh.readline().split()
                if "MaxIBD2" not in head:
                    continue
                c = head.index("MaxIBD2")
                for line in fh:
                    f = line.split()
                    v = float(f[c])
                    if nonzero and v <= 0:
                        continue
                    if wide:
                        i, j = idx[(f[0], f[1])], idx[(f[2], f[3])]
                    else:
                        i, j = idx[(f[0], f[1])], idx[(f[0], f[2])]
                    out.append((name, min(i, j), max(i, j), int(round(v))))
    return out


def score_max(p=BASE, targets=None):
    tg = targets if targets is not None else max_targets()
    ok = 0
    bad = []
    for name, i, j, t in tg:
        ds = kd.load(name)
        g = max_ibd2(ds, i, j, p)
        if g == t:
            ok += 1
        else:
            bad.append((name, i, j, t, g))
    return ok, len(tg), bad


def main():
    import sys
    tg = max_targets()
    print("MaxIBD2 targets:", len(tg))
    ok, n, bad = score_max(BASE, tg)
    print("MaxIBD2 exact: %d/%d" % (ok, n))
    s = score_seg(BASE, per_ds=True)
    print(".seg: exact %(exact)d  both %(both)d  ibd1 %(ibd1)d  ibd2 %(ibd2)d  "
          "of %(rows)d   extra %(extra)d missing %(missing)d  MAE %(mae).5f "
          "worst %(worst).4f" % s)
    for k, v in s["by"].items():
        print("   %-12s exact %3d  both %3d  of %3d" % (k, v[0], v[1], v[2]))
    if "-v" in sys.argv:
        for row in bad:
            print("   MISS %-12s %3d,%-3d want %10d got %10d  d=%+d"
                  % (row[0], row[1], row[2], row[3], row[4], row[4] - row[3]))


if __name__ == "__main__":
    main()
