import pickle, numpy as np, mpmath as mp
from scipy.optimize import brentq
out=pickle.load(open('mp_N14.pkl','rb'))
g=pickle.load(open('zeros.pkl','rb'))[:8]
for lam in [2.0,2.5,3.0,4.0]:
    d=out[lam]; L=d['L']; v=np.array(d['v']); N=len(v)
    a=np.arange(1,N+1)*np.pi/L; sg=(-1.0)**np.arange(1,N+1)
    R=lambda r: np.sum(v*sg*2*a/(r**2-a**2))
    errs=[]
    for gg in g[:7]:
        try: errs.append(brentq(R,gg-0.05,gg+0.05,xtol=1e-14)-gg)
        except Exception: errs.append(np.nan)
    print(lam, ' '.join(f"{e:+.1e}" for e in errs))
