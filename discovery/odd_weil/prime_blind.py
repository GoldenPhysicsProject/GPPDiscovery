"""Lowest-energy prime-blind state. Input: ONLY the Gamma(s/2) archimedean form and the positions log n.
   minimize tr(Z W_inf)  s.t.  Z >= 0, L tr Z = 1, tr(Z G_n) = 0 for all 2 <= n < lambda^2.
Dual:  maximize t  s.t.  W_inf - sum y_n G_n >= t L I.  Then b(n) := y_n sqrt(n)/2 should be Lambda(n)."""
import numpy as np, cvxpy as cp, sys, pickle
from dh_control import build, Gsym_matrix
def vm(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return np.log(p) if m==1 else 0.0
    return 0.0
res={}
for lam,N in [(3.0,16),(4.0,20),(6.0,26),(8.0,30)]:
    L=np.log(lam)
    W=build('zeta',lam,N)['W']; W=(W+W.T)/2
    ns=[n for n in range(2,int(lam**2)+1) if np.log(n)<2*L]
    Gs=[Gsym_matrix(N,L,np.log(n)) for n in ns]
    # solve the dual (cleaner numerically): maximize t over y
    y=cp.Variable(len(ns)); t=cp.Variable()
    M=W-sum(y[i]*Gs[i] for i in range(len(ns)))-t*L*np.eye(N)
    prob=cp.Problem(cp.Maximize(t),[ (M+M.T)/2>>0 ])
    prob.solve(solver=cp.CLARABEL)
    Z=prob.constraints[0].dual_value
    b=y.value*np.sqrt(ns)/2
    ev,V=np.linalg.eigh(Z); rankZ=int(np.sum(ev>1e-6*ev.max()))
    print(f"\nlambda={lam} N={N}: min prime-blind energy = max_b lambda_min = {t.value: .3e}   rank of optimal state = {rankZ}")
    print("   n   recovered b(n)   Lambda(n)")
    for n,bn in zip(ns,b):
        print(f"  {n:3d}  {bn: .6f}      {vm(n):.6f}")
    res[lam]=dict(Z=Z,b=b,ns=ns,L=L,N=N,t=t.value)
pickle.dump(res,open('prime_blind.pkl','wb'))
