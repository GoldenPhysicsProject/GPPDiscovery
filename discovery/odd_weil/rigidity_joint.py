"""Joint perturbations: build J (near-null sensitivities). Directions in ker J escape the first-order
constraints; measure how far one can move along them, especially composite-only directions."""
import mpmath as mp, pickle, sys
from mp_weil import Gsym, vm
mp.mp.dps=50
out=pickle.load(open('mp_N14.pkl','rb')); N=14
def is_pp(n):
    for p in range(2,n+1):
        if n%p==0:
            while n%p==0: n//=p
            return n==1
for lam in [float(x) for x in sys.argv[1:]]:
    d=out[lam]; L=mp.log(lam); A=mp.matrix(d['A'])
    E,Q=mp.eigsy(A)
    ns=[n for n in range(2,int(lam**2)+1) if mp.log(n)<2*L]
    Es={}
    for n in ns:
        Gn=mp.matrix(N,N)
        for j in range(1,N+1):
            for k in range(j,N+1):
                Gn[j-1,k-1]=Gn[k-1,j-1]=Gsym(j,k,mp.log(n),L)
        Es[n]=-2*Gn/mp.sqrt(n)
    Ais=Q*mp.diag([1/mp.sqrt(e) for e in E])*Q.T
    def extent(dirv):
        M=sum((dirv[i]*Es[n] for i,n in enumerate(ns)),mp.zeros(N,N))
        M=Ais*M*Ais; ev,_=mp.eigsy((M+M.T)/2)
        hi=min([-1/e for e in ev if e<0],default=mp.inf); lo=max([-1/e for e in ev if e>0],default=-mp.inf)
        return lo,hi
    # Gram-like sensitivity: S_{n,m} = sum_i over ALL eigvecs weighted 1/E_i -> use generalized problem:
    # maximize |t| over unit composite-only directions: sweep composites jointly via principal directions
    comps=[n for n in ns if not is_pp(n)]; pps=[n for n in ns if is_pp(n)]
    print(f"\nlam={lam}: coefficients={len(ns)} (prime powers {pps}, composites {comps})")
    # best composite-only direction: the one minimizing the max generalized eigenvalue -> search via
    # projecting out near-null sensitivities: Jc v=0 where Jc rows = q_i^T E_n q_i for near-null q_i
    k=sum(1 for e in E if e/L<1e-8)
    import itertools
    J=mp.matrix(k,len(ns))
    for i in range(k):
        qi=Q[:,i]
        for c,n in enumerate(ns): J[i,c]=(qi.T*Es[n]*qi)[0]
    # nullspace of J restricted to composites, else full
    for label,idx in [('composite-only',[ns.index(n) for n in comps]),('all coefficients',list(range(len(ns))))]:
        if not idx: continue
        Jr=mp.matrix([[J[i,c] for c in idx] for i in range(k)]) if k>0 else None
        JtJ=Jr.T*Jr
        ew,EV=mp.eigsy(JtJ)
        scale=max(abs(x) for x in ew)
        nullcols=[r for r in range(len(idx)) if abs(ew[r])<scale*mp.mpf(10)**-40]
        rank=len(idx)-len(nullcols); free=len(nullcols)
        best=(mp.mpf(0),None)
        if free>0:
            for r in nullcols:
                v=[mp.mpf(0)]*len(ns)
                for a,c in enumerate(idx): v[c]=EV[a,r]
                lo,hi=extent(v); w=min(-lo,hi)
                if w>best[0]: best=(w,(lo,hi))
        # also axis directions for reference
        axis=max(min(-extent([1 if c==a else 0 for c in range(len(ns))])[0],extent([1 if c==a else 0 for c in range(len(ns))])[1]) for a in idx)
        print(f"  {label}: near-null constraints k={k}, rank {rank}, free directions {free}; "
              f"best two-sided extent along free dirs = {mp.nstr(best[0],3) if best[1] else 'none'}; best single axis = {mp.nstr(axis,3)}")
