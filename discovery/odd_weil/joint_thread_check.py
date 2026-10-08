"""Independent check of Meta Muse's 2026-10-08 thread-extension claims (feasible-set program).

Same objects as prime_blind.py: for each lambda, A_lam(y) = W_inf - sum_{n<lam^2, log n<2L} y_n G_n, b(n)=y_n sqrt(n)/2.
margin t = lambda_min(A)/L. Lambda (b = von Mangoldt) has margin ~ 0 (boundary point).

Tests (all numerical, binary64 + Clarabel; NOT interval certified):
 T0  single-lambda optimum, uncapped          (reproduces odd_weil correction table)
 T1  single-lambda optimum, cap |b|<=B
 T2  freeze b* (optimum at lambda0) on n<lam0^2, free new coefficients at lambda1: does the margin survive?
 T3  joint optimum over lambda-sets with shared coefficients, cap |b|<=B: margin, distance from Lambda
 T4  joint margin as the lambda-set grows (flat vs shrinking)
"""
import numpy as np, cvxpy as cp, json, sys
from dh_control import build, Gsym_matrix
def vm(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return np.log(p) if m==1 else 0.0
    return 0.0
def Nof(lam): return int(round(16+10*np.log2(lam/3.0)))
_cache={}
def data(lam):
    if lam in _cache: return _cache[lam]
    N=Nof(lam); L=np.log(lam)
    W=build('zeta',lam,N)['W']; W=(W+W.T)/2
    ns=[n for n in range(2,int(lam**2)+1) if np.log(n)<2*L]
    Gs={n:Gsym_matrix(N,L,np.log(n)) for n in ns}
    _cache[lam]=(N,L,W,ns,Gs); return _cache[lam]
def margin_of(lam,b):
    """exact lambda_min(A_lam)/L for coefficient dict b (n->b(n)); missing n -> 0."""
    N,L,W,ns,Gs=data(lam)
    A=W-sum((2*b.get(n,0.0)/np.sqrt(n))*Gs[n] for n in ns)
    return np.linalg.eigvalsh((A+A.T)/2).min()/L
def solve(lams,cap=None,frozen=None,solver=cp.CLARABEL):
    allns=sorted(set(n for lam in lams for n in data(lam)[3]))
    y={n:cp.Variable() for n in allns}; t=cp.Variable(); cons=[]
    for lam in lams:
        N,L,W,ns,Gs=data(lam)
        M=W-sum(y[n]*Gs[n] for n in ns)-t*L*np.eye(N)
        cons.append((M+M.T)/2>>0)
    for n in allns:
        if cap is not None: cons+= [y[n]*np.sqrt(n)/2<=cap, y[n]*np.sqrt(n)/2>=-cap]
        if frozen and n in frozen: cons.append(y[n]==2*frozen[n]/np.sqrt(n))
    prob=cp.Problem(cp.Maximize(t),cons); prob.solve(solver=solver)
    if prob.status not in ('optimal','optimal_inaccurate'): return None,None,prob.status
    b={n:float(y[n].value)*np.sqrt(n)/2 for n in allns}
    return float(t.value),b,prob.status
def dist(b,upto):
    return float(np.sqrt(sum((b[n]-vm(n))**2 for n in b if n<upto)))
out={}
print("T0 single-lambda uncapped")
for lam in [3.0,4.0,6.0,8.0]:
    t,b,st=solve([lam]); print(f"  lam={lam}: t*={t:.3e} ({st}); recheck={margin_of(lam,b):.3e}; Lambda margin={margin_of(lam,{n:vm(n) for n in data(lam)[3]}):.2e}; max|b|={max(abs(v) for v in b.values()):.3g}")
    out[f'T0_{lam}']=t
print("T1 single-lambda, cap |b|<=10")
for lam in [3.0,4.0,6.0,8.0]:
    t,b,st=solve([lam],cap=10); print(f"  lam={lam}: t*={t:.3e}; recheck={margin_of(lam,b):.3e}; dist from Lambda (all n)={dist(b,10**9):.3f}")
    out[f'T1_{lam}']=t
print("T2 freeze optimum from lam0, extend to lam1 (new coefficients free, cap 10 on new ones)")
for lam0,lam1,cap0 in [(3.0,4.0,None),(3.0,4.0,10),(4.0,5.0,10)]:
    t0,b0,_=solve([lam0],cap=cap0)
    frozen={n:b0[n] for n in data(lam0)[3]}
    t1,b1,st=solve([lam1],cap=10,frozen=frozen)
    print(f"  lam0={lam0}->lam1={lam1} (cap0={cap0}): start margin {t0:.3e}; extended max margin = {t1 if t1 is None else f'{t1:.3e}'} ({st})")
    out[f'T2_{lam0}_{lam1}_{cap0}']=t1
print("T3/T4 joint optimum over growing lambda-sets, cap |b|<=10")
sets=[[3.0,4.0,5.0],[3.0,4.0,5.0,6.0],[3.0,4.0,5.0,6.0,8.0],[3.0,4.0,5.0,6.0,8.0,10.0],[3.0,4.0,5.0,6.0,8.0,10.0,12.0]]
for S in sets:
    t,b,st=solve(S,cap=10)
    if t is None: print("  ",S,st); continue
    ms=[margin_of(l,b) for l in S]
    print(f"  set up to {S[-1]}: joint t*={t:.3e}; per-lambda recheck min={min(ms):.3e} max={max(ms):.3e}; dist from Lambda (n<9)={dist(b,9):.3f}; (all n)={dist(b,10**9):.3f}")
    out[f'T3_{S[-1]}']=dict(t=t,dist9=dist(b,9),distall=dist(b,10**9),recheck_min=min(ms))
json.dump(out,open('results/joint_thread_check.json','w'),indent=1)
