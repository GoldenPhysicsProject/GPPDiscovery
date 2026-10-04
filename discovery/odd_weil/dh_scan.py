import numpy as np, pickle, sys
from dh_control import build, S_basis
from scipy.optimize import brentq
dz=pickle.load(open('dh_zeros.pkl','rb')); online=np.array(dz['online'])
res={}
for lam in [3.0,4.0,6.0,8.0,10.0,12.0]:
    L=np.log(lam); N=int(np.ceil(110*L/np.pi))
    row=[]
    for kind in ['zeta','chi5real','dh']:
        d=build(kind,lam,N)
        ev,V=np.linalg.eigh(d['A']/L)
        row.append(ev[0])
        if kind=='dh': res[lam]=dict(ev=ev,V=V,L=L,N=N)
    print(f"lam={lam:5.1f} N={N:3d}  min eig: zeta={row[0]: .2e}  chi5={row[1]: .2e}  DH={row[2]: .3e}",flush=True)
pickle.dump(res,open('dh_scan.pkl','wb'))
