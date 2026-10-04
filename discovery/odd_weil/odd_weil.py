"""
Odd-sector semilocal Weil form, zero-free side, and its Laplacian split.

Conventions (standard explicit formula, Iwaniec-Kowalski style):
  f real odd, supported in [-L, L], L = log(lambda).
  F(z) = int f(u) e^{zu} du.
  G(u) = int f(u+b) f(b) db   (autocorrelation, even), h(r) = |F(ir)|^2.
  Q(f) = sum_rho F(rho-1/2) F(1/2-rho)          (zero side)
  P(f) = 2 F(1/2) F(-1/2)                        (pole side)
  A(f) = Q - P = -log(pi) G(0)
               + (1/2pi) int h(r) Re psi(1/4 + i r/2) dr
               - 2 sum_n Lambda(n) n^{-1/2} G(log n)

Basis: phi_k(u) = sin(k pi u / L), k = 1..N on [-L, L].
"""
import numpy as np
from scipy.special import psi
from numpy.polynomial.legendre import leggauss


def vonmangoldt(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0:
                m //= p
            return np.log(p) if m == 1 else 0.0
    return 0.0


def S_basis(k, r, L):
    """F_k(ir) = i * S_k(r), S_k(r) = int_{-L}^{L} sin(a u) sin(r u) du."""
    a = k * np.pi / L
    return L * (np.sinc((r - a) * L / np.pi) - np.sinc((r + a) * L / np.pi))


def arch_matrix(N, L, R=4000.0, dr=0.002):
    r = np.arange(0.0, R, dr) + dr / 2  # midpoint rule on [0, R]
    w = np.real(psi(0.25 + 0.5j * r))
    S = np.array([S_basis(k, r, L) for k in range(1, N + 1)])
    M = (S * w) @ S.T * dr * 2 / (2 * np.pi)  # even integrand, x2
    # tail r > R: S_j S_k ~ 4 a_j a_k (-1)^{j+k} sin^2(rL)/r^4, avg sin^2 = 1/2
    a = np.arange(1, N + 1) * np.pi / L
    sg = (-1.0) ** np.arange(1, N + 1)
    # int_R^inf log(r/2)/r^4 dr = (3 log(R/2) + 1)/(9 R^3)
    tail = (3 * np.log(R / 2) + 1) / (9 * R**3)
    M += np.outer(a * sg, a * sg) * 4 * 0.5 * tail * 2 / (2 * np.pi)
    G0 = L * np.eye(N)
    return M - np.log(np.pi) * G0


def corr_matrix(N, L, x, nq=400):
    """C_jk(x) = int phi_j(u+x) phi_k(u) du, overlap u in [-L, L-x], x >= 0."""
    if x >= 2 * L:
        return np.zeros((N, N))
    t, wt = leggauss(nq)
    lo, hi = -L, L - x
    u = 0.5 * (hi - lo) * t + 0.5 * (hi + lo)
    wu = 0.5 * (hi - lo) * wt
    a = np.arange(1, N + 1)[:, None] * np.pi / L
    Pj = np.sin(a * (u + x))
    Pk = np.sin(a * u)
    return (Pj * wu) @ Pk.T


def prime_terms(N, L):
    lam2 = np.exp(2 * L)
    terms = []
    for n in range(2, int(np.ceil(lam2)) + 1):
        if np.log(n) >= 2 * L:
            break
        Ln = vonmangoldt(n)
        if Ln == 0:
            continue
        C = corr_matrix(N, L, np.log(n))
        Csym = 0.5 * (C + C.T)
        terms.append((n, Ln / np.sqrt(n), Csym))
    return terms


def build(N, L):
    W = arch_matrix(N, L)
    terms = prime_terms(N, L)
    Pr = sum(2 * c * Cs for _, c, Cs in terms) if terms else np.zeros((N, N))
    A = W - Pr
    Gram = L * np.eye(N)
    S = sum(c for _, c, _ in terms)
    # Laplacian: sum c ||f - T f||^2 = sum c (2 Gram - 2 Csym)
    Lap = sum(c * (2 * Gram - 2 * Cs) for _, c, Cs in terms) if terms else 0 * Gram
    Mass = W - 2 * S * Gram
    return dict(A=A, W=W, Pr=Pr, Lap=Lap, Mass=Mass, S=S, terms=terms, Gram=Gram)


if __name__ == "__main__":
    import sys
    for lam in [1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0]:
        L = np.log(lam)
        N = 14
        d = build(N, L)
        A = d["A"]
        ev = np.linalg.eigvalsh(A / L)  # normalize by Gram
        evM = np.linalg.eigvalsh(d["Mass"] / L)
        evW = np.linalg.eigvalsh(d["W"] / L)
        print(f"lam={lam:4.1f} primes<{lam**2:5.1f} S={d['S']:.3f} "
              f"minA={ev[0]: .3e} minW={evW[0]: .3e} minMass={evM[0]: .3e} "
              f"resid={np.abs(A - (d['Lap'] + d['Mass'])).max():.1e}")
