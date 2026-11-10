"""GW-EF-3 mechanism control (NUMERICAL, run 2026-11-10).  SYNTHETIC zero set: NOT a claim about zeta.

Z' = (first 300 zeros of zeta, on the line, +/- gamma)  U  (an invented off-line quadruple at height gamma_1 = 14.1347..:
rho0 = 1/2 + d + i g1, iota rho0 = 1/2 - d + i g1, and their complex conjugates), d = 0.1.
Q(f) = sum_{rho in Z'} F(rho) F(1 - rho),  F(s) = int f(u) e^{(s-1/2)u} du,  f real.

f_sigma(u) = u exp(-u^2/(2 sigma^2)) cos(g1 u)   (odd Gaussian packet at height g1: F(rho0) ~ -F(iota rho0))
g_sigma(u) =   exp(-u^2/(2 sigma^2)) cos(g1 u)   (even packet: the pairing has the other sign)
Genuine zeros contribute |F|^2 >= 0 (they are on the line); the invented pair contributes the pairedForm value
2 Re( F(rho0) conj F(iota rho0) ) (+ its conjugate pair), which the odd packet makes negative.
"""
import numpy as np
from mpmath import zetazero
g1 = float(zetazero(1).imag)
true_g = np.array([float(zetazero(j).imag) for j in range(1, 301)])
d = 0.1
syn = [0.5 + d + 1j*g1, 0.5 - d + 1j*g1, 0.5 + d - 1j*g1, 0.5 - d - 1j*g1]
h = 2e-3; sig_max = 12
u = np.arange(-9*sig_max, 9*sig_max + h/2, h)
def F_of(fu):
    return lambda s: np.sum(fu * np.exp((s - 0.5) * u)) * h
def Q_parts(fu):
    F = F_of(fu)
    q_true = sum(F(0.5 + 1j*t) * F(0.5 - 1j*t) for t in true_g) + sum(F(0.5 - 1j*t) * F(0.5 + 1j*t) for t in true_g)
    q_syn = sum(F(r) * F(1 - r) for r in syn)
    return q_true.real, q_syn.real, abs(q_true.imag) + abs(q_syn.imag)
print("sigma | packet | Q_true(on-line zeros) | Q_synthetic_pair | Q_total | max|Im|")
for sg in [1, 2, 3, 5, 8]:
    env = np.exp(-u**2 / (2*sg**2)) * np.cos(g1*u)
    for name, fu in [("odd", u*env), ("even", env)]:
        qt, qs, im = Q_parts(fu)
        print(f"{sg:5} | {name:4} | {qt:+.4e} | {qs:+.4e} | {qt+qs:+.4e} | {im:.1e}")
