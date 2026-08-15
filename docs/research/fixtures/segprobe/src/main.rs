//! Scratch scorer for `--build`'s `Join3/Join2`, and a raw IBD-interval dump.
//!
//! Not committed (lives under the gitignored fixtures work dir). Rebuilds the tool
//! `docs/research/fixtures/avfs.py` describes as "not preserved".
//!
//!     segprobe join <bed> R,N1,N2 [more triples...]
//!     segprobe pairs <bed> A,B [more pairs...]     # dump IBD intervals
//!     segprobe scan  <bed> R S1,S2,S3,...          # every pair from a sibship
//!     segprobe all   <bed>                         # every (R;N1,N2) with R 2nd to both

use std::collections::HashMap;

use king_cli::cli;
use king_cli::load;
use king_core::ibdseg::{self, Usable};
use king_io::Genotypes;

struct Ctx {
    variant: u32,
    pos: Vec<i64>,
    segs: Vec<Usable>,
    denom: i64,
    seglen: i64,
    g: Genotypes,
    idx: HashMap<String, usize>,
    iids: Vec<String>,
}

/// Half-open bp intervals of the union of the pair's called IBD1 and IBD2 segments.
///
/// `variant` selects which set the ratio is computed over; 0 is the one the module doc
/// carries (every IBD1 and IBD2 *call*, merged, endpoints at the boundary markers).
fn union_iv(c: &Ctx, i: usize, j: usize) -> Vec<(i64, i64)> {
    let v = c.variant;
    let mut out: Vec<(i64, i64)> = Vec::new();
    for &seg in &c.segs {
        if seg.words() == 0 {
            continue;
        }
        let scan = ibdseg::Scan::new(&c.g, i, j, seg);
        let merge = v != 2;
        let ibd2 = scan.ibd2(&c.pos, c.seglen, merge);
        let ibd1 = scan.ibd1(&c.pos, c.seglen, merge);
        for call in &ibd2 {
            out.push(bounds(c, call.lo, call.hi));
        }
        if v == 5 {
            // IBD2 only.
            continue;
        }
        for call in &ibd1 {
            if v == 1 {
                // Exactly the pieces `IBD1Seg` sums: the IBD2 calls cut out, each piece
                // facing the --seglength floor on its own.
                let mut cur = call.lo;
                for o in &ibd2 {
                    if o.hi < call.lo || o.lo > call.hi {
                        continue;
                    }
                    if o.lo > cur && c.pos[o.lo - 1] - c.pos[cur] >= c.seglen {
                        out.push(bounds(c, cur, o.lo - 1));
                    }
                    cur = cur.max(o.hi + 1);
                }
                if cur <= call.hi && c.pos[call.hi] - c.pos[cur] >= c.seglen {
                    out.push(bounds(c, cur, call.hi));
                }
            } else {
                out.push(bounds(c, call.lo, call.hi));
            }
        }
    }
    out.sort_unstable();
    // Merge the IBD1/IBD2 overlaps into one interval set.
    let mut merged: Vec<(i64, i64)> = Vec::new();
    for iv in out {
        match merged.last_mut() {
            Some(last) if iv.0 <= last.1 => last.1 = last.1.max(iv.1),
            _ => merged.push(iv),
        }
    }
    merged
}

/// One call's bp interval.  Variant 4 pushes each endpoint out to the midpoint of the
/// marker gap it sits on, which is the obvious "refined endpoint" convention.
fn bounds(c: &Ctx, lo: usize, hi: usize) -> (i64, i64) {
    if c.variant != 4 {
        return (c.pos[lo], c.pos[hi]);
    }
    let a = if lo == 0 { c.pos[lo] } else { (c.pos[lo] + c.pos[lo - 1]) / 2 };
    let b = if hi + 1 >= c.pos.len() { c.pos[hi] } else { (c.pos[hi] + c.pos[hi + 1]) / 2 };
    (a, b)
}

fn intersect(a: &[(i64, i64)], b: &[(i64, i64)]) -> Vec<(i64, i64)> {
    let (mut x, mut y, mut out) = (0, 0, Vec::new());
    while x < a.len() && y < b.len() {
        let lo = a[x].0.max(b[y].0);
        let hi = a[x].1.min(b[y].1);
        if lo < hi {
            out.push((lo, hi));
        }
        if a[x].1 < b[y].1 {
            x += 1;
        } else {
            y += 1;
        }
    }
    out
}

fn total(v: &[(i64, i64)]) -> i64 {
    v.iter().map(|(a, b)| b - a).sum()
}

fn ctx(bed: &str) -> Ctx {
    let args: Vec<String> = vec!["-b".into(), bed.into(), "--ibdseg".into()];
    let parsed = cli::parse(&args);
    let mut sink = Vec::new();
    let loaded = load::load(&parsed.options, &mut sink).ok().expect("load");
    let sexchr = i64::from(parsed.options.int(cli::Opt::Sexchr));
    let (mut chr, mut pos) = (Vec::new(), Vec::new());
    for &k in &loaded.fileset.kept {
        let v = &loaded.fileset.variants[k];
        if load::classify(&v.chrom, sexchr).is_autosomal() {
            chr.push(load::chromosome_code(&v.chrom, sexchr));
            pos.push(v.bp);
        }
    }
    let segs = ibdseg::usable_segments(&chr, &pos);
    let denom = ibdseg::denominator(&segs, &pos);
    // Positions restart per chromosome, so lay the chromosomes end to end before any
    // interval arithmetic: within-chromosome differences are untouched and no two
    // chromosomes can overlap.
    let mut gpos = Vec::with_capacity(pos.len());
    let (mut base, mut prev) = (0i64, None);
    for k in 0..pos.len() {
        if prev != Some(chr[k]) {
            if let Some(&last) = gpos.last() {
                base = last + 1_000_000_000;
            }
            prev = Some(chr[k]);
        }
        gpos.push(base + pos[k]);
    }
    let pos = gpos;
    let iids: Vec<String> = loaded
        .fileset
        .samples
        .iter()
        .map(|s| s.iid.clone())
        .collect();
    let idx = iids
        .iter()
        .enumerate()
        .map(|(n, s)| (s.clone(), n))
        .collect();
    Ctx {
        variant: std::env::var("VARIANT").ok().and_then(|v| v.parse().ok()).unwrap_or(0),
        pos,
        segs,
        denom,
        seglen: king_cli::analysis::ibdseg::seglength_bp(&parsed.options),
        g: loaded.fileset.genotypes,
        idx,
        iids,
    }
}

fn join(c: &Ctx, r: usize, n1: usize, n2: usize) -> (i64, i64) {
    let a = union_iv(c, r, n1);
    let b = union_iv(c, r, n2);
    let s = union_iv(c, n1, n2);
    let j2 = intersect(&a, &b);
    let j3 = intersect(&j2, &s);
    (total(&j2), total(&j3))
}

fn main() {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let mode = argv[0].clone();
    let c = ctx(&argv[1]);
    match mode.as_str() {
        "join" => {
            for t in &argv[2..] {
                let n: Vec<usize> = t.split(',').map(|x| c.idx[x]).collect();
                let (j2, j3) = join(&c, n[0], n[1], n[2]);
                println!(
                    "{t}\tJoin2={j2}\tJoin3={j3}\tratio={:.6}\tprinted={:.3}",
                    if j2 == 0 { f64::NAN } else { j3 as f64 / j2 as f64 },
                    if j2 == 0 { f64::NAN } else { j3 as f64 / j2 as f64 }
                );
            }
        }
        "pairs" => {
            for t in &argv[2..] {
                let n: Vec<usize> = t.split(',').map(|x| c.idx[x]).collect();
                let iv = union_iv(&c, n[0], n[1]);
                let s = ibdseg::pair_segments(&c.g, &c.pos, &c.segs, n[0], n[1], c.seglen);
                println!(
                    "{t}\tIBD1Seg={:.4}\tIBD2Seg={:.4}\tPropIBD={:.4}\tunion={}bp\tn={}",
                    s.ibd1_seg(c.denom),
                    s.ibd2_seg(c.denom),
                    s.prop_ibd(c.denom),
                    total(&iv),
                    iv.len()
                );
                for (a, b) in iv {
                    println!("    {a}\t{b}\t{}", b - a);
                }
            }
        }
        "scan" => {
            let r = c.idx[&argv[2]];
            let sibs: Vec<usize> = argv[3].split(',').map(|x| c.idx[x]).collect();
            for a in 0..sibs.len() {
                for b in 0..sibs.len() {
                    if a == b {
                        continue;
                    }
                    let (j2, j3) = join(&c, r, sibs[a], sibs[b]);
                    println!(
                        "{} ; {} {}\tJoin2={j2}\tJoin3={j3}\tratio={:.4}",
                        c.iids[r],
                        c.iids[sibs[a]],
                        c.iids[sibs[b]],
                        j3 as f64 / j2 as f64
                    );
                }
            }
        }
        _ => eprintln!("modes: join | pairs | scan"),
    }
}
