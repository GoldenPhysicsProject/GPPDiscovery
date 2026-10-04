"""Scale all prime weights: b = c * Lambda. Find the interval of c where the odd form stays PSD."""
import numpy as np
from dh_control import build
for lam,N in [(2.0,14),(3.0,16),(4.0,20),(6.0,26),(8.0,30),(10.0,34)]:
    d=build('zeta',lam,N); W=d['W']; Pr=d['Pr']; L=d['L']
    f=lambda c: np.linalg.eigvalsh((W-c*Pr+(W-c*Pr).T)/2/L)[0]
    cs=np.linspace(0,1.5,301); v=np.array([f(c) for c in cs])
    ok=cs[v>-1e-11]
    print(f"lambda={lam:4.1f}: PSD for c in [{ok.min():.3f}, {ok.max():.3f}]   "
          f"min eig at c=0.95: {f(0.95): .2e}, c=1: {f(1.0): .2e}, c=1.01: {f(1.01): .2e}, c=1.05: {f(1.05): .2e}")
