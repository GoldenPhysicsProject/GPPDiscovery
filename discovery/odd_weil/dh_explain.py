import numpy as np, pickle, mpmath as mp
dz=pickle.load(open('dh_zeros.pkl','rb')); online=np.array(dz['online']); off=dz['off'][0]
res=pickle.load(open('dh_scan.pkl','rb'))
def Fz(v,z,L):
    k=np.arange(1,len(v)+1); a=k*np.pi/L
    return np.sum(v*(-2*a*(-1.0)**k*np.sinh(z*L)/(z**2+a**2)))
for lam in [4.0,8.0,10.0]:
    r=res[lam]; L=r['L']; v=r['V'][:,0]; ev=r['ev'][0]
    rs=np.linspace(0.5,120,6000); spec=np.array([abs(Fz(v,1j*x,L))**2 for x in rs])
    peak=rs[np.argmax(spec)]
    on=sum(2*abs(Fz(v,1j*g,L))**2 for g in online)
    z=off-0.5  # delta + i gamma
    quart=sum((Fz(v,w,L)*Fz(v,-w,L)) for w in [z,np.conj(z)])*2  # rho,conj and mirrors
    print(f"lam={lam}: eig*L={ev*L: .4e}  spectral peak r={peak:.2f}  "
          f"zero-side: on-line(<100)={on.real: .4e}  off-line quartet={quart.real: .4e}  sum={on.real+quart.real: .4e}")
# validation: near-null DH vector at lam=4 vanishes at DH on-line zeros
r=res[4.0]; L=r['L']; v=r['V'][:,0]
from scipy.optimize import brentq
f=lambda x: (Fz(v,1j*x,L)).imag
errs=[]
for g in online[:8]:
    try: errs.append(brentq(f,g-0.1,g+0.1)-g)
    except Exception: errs.append(np.nan)
print("lam=4 DH near-null vector: zero errors vs DH on-line zeros:", ' '.join(f"{e:+.1e}" for e in errs))
