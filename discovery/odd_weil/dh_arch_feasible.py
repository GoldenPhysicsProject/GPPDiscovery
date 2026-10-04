"""Same archimedean place as DH (Gamma((s+1)/2), q=5). Compare DH's coefficients with the
Euler-product point b(n) = Lambda(n) Re chi(n), chi mod 5 with chi(2)=i (sqrt of L(chi)L(chibar))."""
import numpy as np
from dh_control import build, Gsym_matrix, dirichlet_coeffs
def vm(n):
    for p in range(2,n+1):
        if n%p==0:
            m=n
            while m%p==0: m//=p
            return np.log(p) if m==1 else 0.0
    return 0.0
rechi={0:0.0,1:1.0,2:0.0,3:0.0,4:-1.0}
for lam in [6.0,8.0,10.0,12.0]:
    L=np.log(lam); N=int(np.ceil(110*L/np.pi))
    d=build('dh',lam,N)
    W=d['W']; M=int(lam**2)
    A2=W.copy()
    for n in range(2,M+1):
        if np.log(n)>=2*L: continue
        b=vm(n)*rechi[n%5]
        if b: A2-=2*b/np.sqrt(n)*Gsym_matrix(N,L,np.log(n))
    e1=np.linalg.eigvalsh(d['A']/L)[0]; e2=np.linalg.eigvalsh(A2/L)[0]
    print(f"lam={lam:5.1f}  same Gamma((s+1)/2), q=5:  DH coefficients min eig = {e1: .3e}   Lambda*Re(chi) min eig = {e2: .3e}")
