import numpy as np, pickle
from dh_control import build, Gsym_matrix, dirichlet_coeffs
def is_pp(n):
    for p in range(2,n+1):
        if n%p==0:
            while n%p==0: n//=p
            return n==1
for lam in [8.0,10.0]:
    L=np.log(lam); N=int(np.ceil(110*L/np.pi))
    d=build('dh',lam,N); ev,V=np.linalg.eigh(d['A']); v=V[:,0]
    _,b=dirichlet_coeffs('dh',int(lam**2)+1)
    pp=comp=0.0; top=[]
    for n in range(2,int(lam**2)+1):
        if np.log(n)>=2*L or b[n]==0: continue
        t=-2*b[n]/np.sqrt(n)*v@Gsym_matrix(N,L,np.log(n))@v
        (pp:=pp+t) if is_pp(n) else (comp:=comp+t)
        top.append((t,n))
    top.sort()
    print(f"lam={lam}: A(v)={ev[0]:.4f} = arch {v@d['W']@v:.4f} + prime-power terms {pp:.4f} + composite terms {comp:.4f}")
    print("   most negative single terms (value, n):",[(round(t,3),n) for t,n in top[:6]])
