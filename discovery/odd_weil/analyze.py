import pickle, numpy as np, mpmath as mp
from scipy.linalg import eigh
from scipy.optimize import brentq
out=pickle.load(open('mp_N14.pkl','rb'))
gam=[14.134725,21.022040,25.010858,30.424876,32.935062,37.586178,40.918719,43.327073,48.005151,49.773832]
for lam,d in out.items():
    L=d['L']; N=len(d['v']); v=np.array(d['v'])
    a=np.arange(1,N+1)*np.pi/L; sg=(-1.0)**np.arange(1,N+1)
    R=lambda r: np.sum(v*sg*2*a/(r**2-a**2))
    F=lambda r: np.sin(r*L)*R(r)
    rs=np.linspace(0.05,55,200000); vals=np.array([F(r) for r in rs])
    z=[brentq(F,rs[i],rs[i+1]) for i in range(len(rs)-1) if vals[i]*vals[i+1]<0]
    grid=[m*np.pi/L for m in range(N+1,40)]
    nontriv=[x for x in z if min(abs(x-g) for g in grid)>1e-3 and min(abs(x-ak) for ak in a)>1e-6]
    print(f"lam={lam}: zeros of F(ir) (excluding trivial m*pi/L grid):",
          ' '.join(f"{x:.4f}" for x in nontriv[:8]))
    # Laplacian split margin
    Lap=np.array(d['Lap'],dtype=float); Mass=np.array(d['Mass'],dtype=float)
    mu=eigh(-Mass,Lap,eigvals_only=True)
    print(f"   max mu = {mu[-1]:.12f} (A>=0 iff mu<=1);  modes of -Mass top eig vector dominant k:",
          np.argmax(np.abs(eigh(-Mass)[1][:,-1]))+1)
print("zeta gammas:", gam[:8])
