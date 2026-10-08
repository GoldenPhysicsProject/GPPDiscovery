"""Check ftheta.build()['Q'] against the zero side: sum over zeta zeros rho=1/2+i gamma of H(rho+theta)+H(rho-theta)."""
import numpy as np, mpmath as mp, os, pickle
from ftheta import build
if not os.path.exists('zeros2000.pkl'):
    pickle.dump([float(mp.zetazero(n).imag) for n in range(1, 2001)], open('zeros2000.pkl', 'wb'))
g = np.array(pickle.load(open('zeros2000.pkl', 'rb')))
def Fz(N, L, z):  # int_{-L}^{L} sin(a x) e^{z x} dx, complex z
    a = np.arange(1, N+1)*np.pi/L
    E = lambda zz: np.exp(zz*L)*(zz*np.sin(a*L)-a*np.cos(a*L)) - np.exp(-zz*L)*(-zz*np.sin(a*L)-a*np.cos(a*L))
    return E(z)/(z*z+a*a)
for th, lam in [(0.0, 3), (0.2, 3), (0.3, 4)]:
    L = np.log(lam); N = 8
    d = build(lam, N, th)
    Z = np.zeros((N, N))
    for gm in g:
        for z in [1j*gm+th, 1j*gm-th, -1j*gm+th, -1j*gm-th]:
            A, B = Fz(N, L, z), Fz(N, L, -z)
            Z += np.real(0.5*(np.outer(A, B)+np.outer(B, A)))
    print(f"theta={th} lam={lam}: max|Q_primes - Q_zeros| = {np.abs(d['Q']-Z).max():.2e}  (scale {np.abs(Z).max():.2e}),"
          f" min eig zero-side {np.linalg.eigvalsh(Z/L)[0]: .3e} vs prime-side {np.linalg.eigvalsh(d['Q']/L)[0]: .3e}")
