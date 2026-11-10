"""GW-EF-2 (NUMERICAL, run 2026-11-10): Phi(autocorr f)(s) = Psi f(s) Psi f(1-s) = Psi f(s) conj(Psi f(iota s)),
iota s = 1 - conj(s), for real, non-even f.  Grid quadrature (trapezoid, h = 1e-3 on [-14, 14])."""
import numpy as np
h = 1e-3
y = np.arange(-14, 14 + h/2, h)
f = (1 + y + 0.3*y**2) * np.exp(-y**2 - 0.4*y)
# autocorr g(x) = int f(y) f(y-x) dy on the lag grid x = k h
g = np.correlate(f, f, mode="full") * h          # g[k] = sum_y f(y) f(y - x_k)-type; orientation checked below
lags = (np.arange(len(g)) - (len(f) - 1)) * h
# np.correlate(f, f, 'full')[k] = sum_n f[n + lag] f[n]  with lag = k - (N-1): equals int f(y+x) f(y) dy = g(-x) ... fix:
# int f(y) f(y-x) dy = int f(z+x) f(z) dz  (z = y - x)  = correlate value at lag +x.  So g(x) is correlate at lag x. OK.
def Psi(s):  return np.sum(f * np.exp((s - 0.5) * y)) * h
def Phi(s):  return np.sum(g * np.exp((s - 0.5) * lags)) * h
for s in [0.5 + 14.1347j, 0.8 + 3.0j, 0.2 - 2.0j, 1.3 + 0.7j, 0.35 + 21.022j]:
    lhs = Phi(s); r1 = Psi(s) * Psi(1 - s); r2 = Psi(s) * np.conj(Psi(1 - np.conj(s)))
    print(f"s={s}: |Phi - Psi(s)Psi(1-s)|={abs(lhs-r1):.2e}  |Phi - Psi(s)conj Psi(iota s)|={abs(lhs-r2):.2e}  |Phi|={abs(lhs):.3e}")
