import numpy as np, kingdata as kd, ibsmeasure as M
WORD=64; PC=np.bitwise_count; DIRTY=5
def segs(ds,i,j):
    ibs0, ibs1, _, _ = ds.masks(i,j)
    n0,n1 = PC(ibs0), PC(ibs1)
    pos=ds.pos; out=[]
    for _,lo,hi in ds.segs:
        w0,w1=-(-lo//WORD),(hi+1)//WORD-1
        if w1<w0: continue
        clean=[(n0[w]==0) and (n1[w]<DIRTY) for w in range(w0,w1+1)]
        ok=list(clean)
        for k in range(1,len(clean)-1):
            if not clean[k] and clean[k-1] and clean[k+1] and n0[w0+k]==0: ok[k]=True
        k=0
        while k<len(ok):
            if not ok[k]: k+=1; continue
            k0=k
            while k<len(ok) and ok[k]: k+=1
            u,v=w0+k0,w0+k-1
            e=w1 if v+2>=w1 else v+1
            out.append((w0,w1,u,v,e,int(pos[WORD*e+WORD-1]-pos[WORD*u])))
    return out
for name in ("nuclear","dups","monomorphic","threegen"):
    ds=kd.load(name); t=M.targets(ds)
    print(f"== {name}  D={ds.denom}")
    for (i,j),(tm,tp) in sorted(t.items()):
        if float(tm)<=0: continue
        s=segs(ds,i,j)
        tot=sum(x[5] for x in s); mx=max((x[5] for x in s),default=0)
        print(f"  {i},{j} target max={tm} pr={tp} | ours max={mx}.000 pr={tot/ds.denom:.4f} "
              f"({len(s)} segs) targetTotal={float(tp)*ds.denom:.0f} ourTotal={tot}")
        for x in s: print("       ", x)
