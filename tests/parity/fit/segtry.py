"""Does the chunk geometry carry the .seg pass too?  No — see
`docs/research/16-segment-extension.md` §9. Four exact rows better and a much
better worst row, but the mean absolute PropIBD error nearly triples, which is
what a *partly* right rule looks like. NOT committed: `Scan::ibd2` still runs the
rule of `14-ibd2-geometry.md`."""
import sys; sys.path.insert(0,'/Users/wgu/Desktop/open-king/tests/parity/fit')
import numpy as np, engine as E, kingdata as kd

WORD = 64

def ibd2_chunk_seg(self, pos, min_bp):
    """.seg IBD2 with the chunk scan deciding the run's confirmed end."""
    p = self.p
    if self.n == 0:
        return []
    sl = slice(self.w0, self.w1 + 1)
    clean = (self.n0[sl] == 0) & (self.n1[sl] < 5)
    ok = clean.copy()
    n0 = self.n0[sl]
    for k in range(1, self.n - 1):
        if not clean[k] and clean[k-1] and clean[k+1] and n0[k] == 0:
            ok[k] = True
    out = []
    for a, b in E._runs(ok):
        lo_w, hi_w = self.w0 + a, self.w0 + b
        scan_last = min(hi_w + 1, self.w1)
        exempt = hi_w + 1 >= self.w1
        u, mis, het, conf, cstart = lo_w, 0, 0, None, lo_w
        k = u
        while k <= scan_last:
            m, h = int(self.n1[k]), int(self.nhh[k])
            mis += m; het += h
            if mis >= 5:
                if (het >= 95 and k - cstart + 1 >= 3) or (exempt and k >= scan_last):
                    conf, mis, het, cstart = k, 0, 0, k + 1
                else:
                    if conf is not None:
                        out = self._emit_words(out, u, E._chunk_extend(self.n1, conf, hi_w), pos, min_bp)
                    u, mis, het, conf, cstart = k + 1, 0, 0, None, k + 1
                    k = u
                    continue
            k += 1
        if exempt:
            out = self._emit_words(out, u, self.w1, pos, min_bp)
        elif conf is not None:
            out = self._emit_words(out, u, self.w1 if hi_w + 2 >= self.w1
                                   else E._chunk_extend(self.n1, conf, hi_w), pos, min_bp)
    return out

def _emit_words(self, out, u, e, pos, min_bp):
    if u > e or e + 1 - u < 3:
        return out
    if not self.informative(self.cum2, u, min(e, self.w1)):
        return out
    lo = WORD * u if u != self.w0 else self.lo
    if e == self.w1:
        hi = self.hi
    else:
        m = int(self.ibs0[e])
        hi = WORD * e + (63 if m == 0 else self._pick(m, self.p.ibd2_right))
    return self._emit(out, lo, hi, pos, min_bp)

E.SegScan._emit_words = _emit_words
print('committed .seg :', E.score_seg())
E.SegScan.ibd2 = ibd2_chunk_seg
print('chunk .seg     :', E.score_seg())
