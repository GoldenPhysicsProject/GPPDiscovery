import numpy as np, mpmath as mp, pickle, os
from odd_weil import build, S_basis
if not os.path.exists('zeros.pkl'):
    g=[float(mp.zetazero(n).imag) for n in range(1,1501)]
    pickle.dump(g,open('zeros.pkl','wb'))
g=np.array(pickle.load(open('zeros.pkl','rb')))
for lam in [2.0,3.0]:
    L=np.log(lam); N=6
    d=build(N,L)
    S=np.array([S_basis(k,g,L) for k in range(1,N+1)])
    Q=2*S@S.T
    a=np.arange(1,N+1)*np.pi/L
    # F_k(1/2)=int sin(a u) e^{u/2} du over [-L,L]
    from scipy.integrate import quad
    F=np.array([quad(lambda u:np.sin(ak*u)*np.exp(u/2),-L,L)[0] for ak in a])
    Az=Q+2*np.outer(F,F)
    print(lam, np.abs(Az-d['A']).max(), np.abs(d['A']).max())
