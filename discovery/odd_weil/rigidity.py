"""Rigidity test at 50 digits: perturb a single coefficient b(n) -> Lambda(n)+t and find
the interval of t keeping the odd semilocal Weil form PSD (zeta, N=14)."""
import mpmath as mp, pickle, sys
from mp_weil import Gsym, vm
mp.mp.dps=50
out=pickle.load(open('mp_N14.pkl','rb'))
N=14
for lam in [float(x) for x in sys.argv[1:]]:
    d=out[lam]; L=mp.log(lam)
    A=mp.matrix(d['A'])
    E,Q=mp.eigsy(A)
    # A^{-1/2}
    Ais=Q*mp.diag([1/mp.sqrt(e) for e in E])*Q.T
    nnull=sum(1 for e in E if e/L<1e-8)
    print(f"\nlam={lam}: eigenvalues/L below 1e-8: {nnull}   (min {mp.nstr(E[0]/L,4)})")
    print("  n   Lambda(n)   allowed t in [t_lo, t_hi]  (b(n)=Lambda(n)+t keeps A PSD)")
    n=2
    while mp.log(n)<2*L:
        Gn=mp.matrix(N,N)
        for j in range(1,N+1):
            for k in range(j,N+1):
                Gn[j-1,k-1]=Gn[k-1,j-1]=Gsym(j,k,mp.log(n),L)
        Edir=-2*Gn/mp.sqrt(n)          # dA/db(n)
        M=Ais*Edir*Ais
        ev,_=mp.eigsy((M+M.T)/2)
        # A + t Edir >= 0  <=>  1 + t*ev_i >= 0 for all i
        hi=min([-1/e for e in ev if e<0],default=mp.inf)
        lo=max([-1/e for e in ev if e>0],default=-mp.inf)
        print(f"  {n:2d}  {float(vm(n)) if vm(n) else 0:8.4f}   [{mp.nstr(lo,3):>10}, {mp.nstr(hi,3):>10}]")
        n+=1
