"""Shifted-product control F_theta(s) = xi(s-theta) xi(s+theta) in the odd fixed-window Weil form,
with and without the radial Euler compensation Lambda(n) -> Lambda(n)(1 - 1/n).

Exact reduction (derived in 2026-10-08_compensation_vs_shifted_control.md):
  zero-sum of F_theta on h = f*f~  =  2 * E_zeta[ h(u) cosh(theta u) ],
where E_zeta is zeta's explicit-formula zero functional:
  E_zeta[k] = W_inf[k] - 2 sum_n Lambda(n) n^{-1/2} k(log n) + Pole[k],
  Pole[k]_ij = -F_i(theta+1/2)F_j(theta+1/2) - F_i(theta-1/2)F_j(theta-1/2),  F_i(c) = int f_i e^{cx}.
Equivalently 2 Re Q_zeta(e^{-theta x} f, e^{theta x} f): a cross term of zeta's zero form.
"""
import sys, numpy as np
sys.path.insert(0, '../odd_weil')
from numpy.polynomial.legendre import leggauss
from dh_control import Gsym_matrix, GAMMA_E

def vm(M):
    lam = np.zeros(M + 1)
    for p in range(2, M + 1):
        if all(p % q for q in range(2, int(p**.5) + 1)):
            q = p
            while q <= M:
                lam[q] = np.log(p); q *= p
    return lam

def Fvec(N, L, c):
    a = np.arange(1, N + 1) * np.pi / L
    # int_{-L}^{L} sin(a x) e^{c x} dx
    return np.array([(np.exp(c*L)*(c*np.sin(ak*L)-ak*np.cos(ak*L)) - np.exp(-c*L)*(-c*np.sin(ak*L)-ak*np.cos(ak*L)))/(c*c+ak*ak) for ak in a])

def build(lam, N, theta, comp=False, nq=1200):
    L = np.log(lam); G0 = L*np.eye(N)
    t, w = leggauss(nq); us = L*(t+1); ws = L*w
    W = (-np.log(np.pi) - GAMMA_E)*G0
    for u, wt in zip(us, ws):
        W += wt*2*(np.exp(-2*u)*G0 - np.exp(-0.5*u)*np.cosh(theta*u)*Gsym_matrix(N, L, u))/(-np.expm1(-2*u))
    W += -G0*np.log(-np.expm1(-4*L))
    M = int(np.floor(lam**2)); Lam = vm(M)
    Pr = np.zeros((N, N)); Corr = np.zeros((N, N))
    for n in range(2, M+1):
        if Lam[n] == 0 or np.log(n) >= 2*L: continue
        Gn = Gsym_matrix(N, L, np.log(n))*np.cosh(theta*np.log(n))
        Pr += 2*Lam[n]/np.sqrt(n)*Gn
        Corr += 2*Lam[n]/n**1.5*Gn
    Fp, Fm = Fvec(N, L, theta+.5), Fvec(N, L, theta-.5)
    Pole = -np.outer(Fp, Fp) - np.outer(Fm, Fm)
    Q = 2*(W - Pr + Pole)              # true zero-sum form of F_theta (theta=0: 2x zeta)
    Qc = Q + 2*Corr                    # compensated: Lambda(n) -> Lambda(n)(1-1/n)
    A_plain = 2*(W - Pr)               # pole-free 'A' convention used in odd_weil notes
    return dict(Q=Q, Qc=Qc, A=A_plain, Ac=A_plain+2*Corr, Corr=2*Corr, L=L)

if __name__ == '__main__':
    for th in [0.0, 0.1, 0.2, 0.3, 0.375]:
        for lam in [2, 3, 4, 6, 8]:
            L = np.log(lam); N = int(np.ceil(60*L/np.pi))
            d = build(lam, N, th)
            e = lambda X: np.linalg.eigvalsh(X/L)[0]
            print(f"theta={th:5.3f} lam={lam:2d} N={N:3d}  minQ={e(d['Q']): .3e}  minQcomp={e(d['Qc']): .3e}  "
                  f"minA={e(d['A']): .3e}  minAcomp={e(d['Ac']): .3e}  ||Corr||={np.linalg.norm(d['Corr']/L,2):.3f}", flush=True)
