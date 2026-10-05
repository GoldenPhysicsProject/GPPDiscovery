"""
Davenport-Heilbronn control for the odd semilocal Weil form.

DH: f(s) = sum a(n) n^{-s}, a(n) by n mod 5: 1->1, 2->kappa, 3->-kappa, 4->-1, 0->0,
kappa = (sqrt(10-2 sqrt5) - 2)/(sqrt5 - 1).
Completed Lambda(s) = (5/pi)^{s/2} Gamma((s+1)/2) f(s) = Lambda(1-s): the same
mirror symmetry D as an L-function, real coefficients, NO Euler product.

Zero-free side (semilocal, exact finite sum because g has compact support):
  A(f) = [log(q/pi) - gammaE] G(0)
         + int_0^inf 2[e^{-2u} G(0) - e^{-(1/2 + kappa_inf) u} G(u)]/(1-e^{-2u}) du
         - 2 sum_n b(n) n^{-1/2} G(log n),
  -f'/f = sum b(n) n^{-s}.
For zeta: q = 1, kappa_inf = 0 (Gamma(s/2)), b = Lambda, plus pole terms (removed by
working with A = Q - P).  For DH: q = 5, kappa_inf = 1 (Gamma((s+1)/2)), no poles.
"""
import numpy as np
from numpy.polynomial.legendre import leggauss

GAMMA_E = 0.5772156649015329


def dirichlet_coeffs(kind, M):
    a = np.zeros(M + 1)
    if kind == "zeta":
        a[1:] = 1.0
    elif kind == "dh":
        kap = (np.sqrt(10 - 2 * np.sqrt(5)) - 2) / (np.sqrt(5) - 1)
        tab = {0: 0.0, 1: 1.0, 2: kap, 3: -kap, 4: -1.0}
        for n in range(1, M + 1):
            a[n] = tab[n % 5]
    elif kind == "chi4":  # odd character mod 4 (Q(i)), Gamma((s+1)/2)
        tab = {0: 0.0, 1: 1.0, 2: 0.0, 3: -1.0}
        for n in range(1, M + 1):
            a[n] = tab[n % 4]
    elif kind == "chi5real":  # real character mod 5 (Legendre), has Euler product
        tab = {0: 0.0, 1: 1.0, 2: -1.0, 3: -1.0, 4: 1.0}
        for n in range(1, M + 1):
            a[n] = tab[n % 5]
    # Dirichlet inverse of a
    inv = np.zeros(M + 1)
    inv[1] = 1.0 / a[1]
    for n in range(2, M + 1):
        s = 0.0
        for d in range(1, n):
            if n % d == 0:
                s += inv[d] * a[n // d]
        inv[n] = -s / a[1]
    # b = (a * log) conv inv
    al = a * np.log(np.maximum(np.arange(M + 1), 1))
    b = np.zeros(M + 1)
    for n in range(1, M + 1):
        s = 0.0
        for d in range(1, n + 1):
            if n % d == 0:
                s += al[d] * inv[n // d]
        b[n] = s
    return a, b


def Gsym_matrix(N, L, u):
    """Gsym[j,k](u) for scalar u in [0, 2L]."""
    if u >= 2 * L:
        return np.zeros((N, N))
    a = np.arange(1, N + 1) * np.pi / L
    lo, hi = -L, L - u
    aj = a[:, None]
    ak = a[None, :]

    def I(c, d):
        with np.errstate(divide="ignore", invalid="ignore"):
            out = (np.sin(c + d * hi) - np.sin(c + d * lo)) / d
        out = np.where(np.abs(d) < 1e-12, np.cos(c) * (hi - lo), out)
        return out

    C = 0.5 * (I(aj * u, aj - ak) - I(aj * u, aj + ak))
    return 0.5 * (C + C.T)


def build(kind, lam, N, nq=1200):
    L = np.log(lam)
    q, kinf = {"zeta": (1.0, 0.0), "dh": (5.0, 1.0), "chi5real": (5.0, 0.0), "chi4": (4.0, 1.0)}[kind]
    G0 = L * np.eye(N)
    # Archimedean, u-space Gauss-Legendre on [0, 2L]
    t, w = leggauss(nq)
    us = L * (t + 1)  # [0, 2L]
    ws = L * w
    W = (np.log(q / np.pi) - GAMMA_E) * G0
    for u, wt in zip(us, ws):
        Gu = Gsym_matrix(N, L, u)
        W += wt * 2 * (np.exp(-2 * u) * G0 - np.exp(-(0.5 + kinf) * u) * Gu) / (-np.expm1(-2 * u))
    W += -G0 * np.log(-np.expm1(-4 * L))  # tail u > 2L where G = 0
    M = int(np.floor(lam**2))
    _, b = dirichlet_coeffs(kind, M + 1)
    Pr = np.zeros((N, N))
    Lap = np.zeros((N, N))
    Smass = 0.0
    for n in range(2, M + 1):
        if np.log(n) >= 2 * L or b[n] == 0:
            continue
        c = b[n] / np.sqrt(n)
        Gn = Gsym_matrix(N, L, np.log(n))
        Pr += 2 * c * Gn
        # passive split works for either sign: -2c<Tf,f> = |c|(||f -+ Tf||^2) - 2|c| ||f||^2
        Lap += abs(c) * (2 * G0 - 2 * np.sign(c) * Gn)
        Smass += abs(c)
    A = W - Pr
    return dict(A=A, W=W, Pr=Pr, Lap=Lap, Mass=W - 2 * Smass * G0, S=Smass, L=L, b=b)


def S_basis(k, r, L):
    a = k * np.pi / L
    return L * (np.sinc((r - a) * L / np.pi) - np.sinc((r + a) * L / np.pi))
