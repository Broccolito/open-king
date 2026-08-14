#!/usr/bin/env python3
"""Reconstruct <prefix>splitped.txt from a .fam and diff against a golden file."""
import sys


def key(s):
    """KING's ID comparator (BEHAVIOR.md Q6): digit runs by (len, bytes),
    non-digits ASCII upper-folded, non-digit before digit, prefix first."""
    out = []
    i = 0
    while i < len(s):
        if s[i].isdigit():
            j = i
            while j < len(s) and s[j].isdigit():
                j += 1
            out.append((1, len(s[i:j]), s[i:j]))
            i = j
        else:
            out.append((0, 0, s[i].upper()))
            i += 1
    return out


def build(famfile):
    fam = [l.split() for l in open(famfile) if l.strip()]
    order = []                      # FIDs in first-appearance order
    members = {}                    # fid -> list of dicts
    for f in fam:
        fid, iid, fa, mo, sex, phe = f[0], f[1], f[2], f[3], f[4], f[5]
        if fid not in members:
            members[fid] = []
            order.append(fid)
        members[fid].append(dict(fid=fid, iid=iid, fa=fa, mo=mo, sex=sex,
                                 phe=('0' if phe == '-9' else phe), dummy='0'))
    out = []
    counter = [0]
    allids = {f[1]: f for f in fam}
    phantoms = {}
    for fid in order:               # KING<n> numbering follows .fam family order
        ms = members[fid]
        if len(ms) < 2 and all(m['fa'] == '0' and m['mo'] == '0' for m in ms):
            continue                # a lone parentless founder is dropped
        have = {m['iid'] for m in ms}
        phantom = []
        for m in ms:
            for slot, sex in (('fa', '1'), ('mo', '2')):
                p = m[slot]
                if p != '0' and p not in have:
                    src = allids.get(p)
                    psex = src[4] if src else sex
                    phantom.append(dict(fid=fid, iid=p, fa='0', mo='0',
                                        sex=psex, phe='0', dummy='1'))
                    have.add(p)
            if (m['fa'] == '0') != (m['mo'] == '0'):
                slot, sex = ('fa', '1') if m['fa'] == '0' else ('mo', '2')
                counter[0] += 1
                name = f"KING{counter[0]}"
                phantom.append(dict(fid=fid, iid=name, fa='0', mo='0',
                                    sex=sex, phe='0', dummy='1'))
                m[slot] = name
                have.add(name)
        phantoms[fid] = phantom
    for fid in sorted(phantoms, key=key):
        ms, phantom = members[fid], phantoms[fid]
        people = phantom + ms
        by = {p['iid']: p for p in people}
        # connected components over parent-child edges
        parent = {p['iid']: p['iid'] for p in people}

        def find(x):
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[rb] = ra
        for p in people:
            for slot in ('fa', 'mo'):
                if p[slot] != '0' and p[slot] in by:
                    union(p[slot], p['iid'])
        depth = {}

        def d(x):
            if x in depth:
                return depth[x]
            p = by[x]
            ps = [p[s] for s in ('fa', 'mo') if p[s] != '0' and p[s] in by]
            depth[x] = 0 if not ps else 1 + max(d(q) for q in ps)
            return depth[x]
        people.sort(key=lambda p: (d(p['iid']), key(p['iid'])))
        adj = {p['iid']: [] for p in people}
        for p in people:
            for slot in ('fa', 'mo'):
                q = p[slot]
                if q != '0' and q in by:
                    adj[q].append(p['iid'])
                    adj[p['iid']].append(q)
        for a in adj:
            adj[a].sort(key=lambda x: [q['iid'] for q in people].index(x))
        comps = []
        done = set()
        for p in people:
            if p['iid'] in done:
                continue
            queue = [p['iid']]
            done.add(p['iid'])
            comp = []
            while queue:
                cur = queue.pop(0)
                comp.append(by[cur])
                for nb in adj[cur]:
                    if nb not in done:
                        done.add(nb)
                        queue.append(nb)
            comps.append(comp)
        if len(comps) == 1:
            comps = [people]
        for ci, comp in enumerate(comps):
            newfid = fid if len(comps) == 1 else f"{fid}_S{ci + 1}"
            for p in comp:
                out.append(f"{fid} {p['iid']} {newfid} {p['iid']} "
                           f"{p['fa']} {p['mo']} {p['sex']} {p['phe']} {p['dummy']}")
    return "\n".join(out) + ("\n" if out else "")


if __name__ == "__main__":
    got = build(sys.argv[1])
    if len(sys.argv) > 2:
        want = open(sys.argv[2]).read()
        print("MATCH" if got == want else "DIFFER")
        if got != want:
            g, w = got.splitlines(), want.splitlines()
            for i in range(max(len(g), len(w))):
                a = g[i] if i < len(g) else "<none>"
                b = w[i] if i < len(w) else "<none>"
                if a != b:
                    print(f"  line {i+1}:\n    ours   {a}\n    golden {b}")
    else:
        sys.stdout.write(got)
