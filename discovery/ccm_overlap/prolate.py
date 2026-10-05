"""Prolate eigenfunctions h_{n,lam} of PW_lam = -d/dx((lam^2-x^2) d/dx) + (2 pi lam x)^2 on [-lam, lam]
(Connes-Consani-Moscovici, arXiv:2511.22755, eq. 7.5), in an orthonormal Legendre basis.
x = lam s:  -d/ds((1-s^2) d/ds) + c^2 s^2,  c = 2 pi lam^2.  p_j = sqrt((2j+1)/2) P_j.
Even sector (indices j = 0,2,4,...) is tridiagonal.  h_n labelled like Hermite functions:
n even <-> (n/2)-th even eigenvalue (h_0: lowest even, h_4: third even)."""
import numpy as np
import mpmath as mp
from scipy.linalg import eigh_tridiagonal


def _a(j):
    return mp.mpf(j + 1) / mp.sqrt(mp.mpf((2 * j + 1) * (2 * j + 3))) if j >= 0 else mp.mpf(0)


def even_tridiag(c2, K):
    """diag and offdiag of the even-sector matrix, indices j=2k, k=0..K-1 (c2 = c^2, mp)."""
    d, e = [], []
    for k in range(K):
        j = 2 * k
        d.append(mp.mpf(j * (j + 1)) + c2 * (_a(j - 1) ** 2 + _a(j) ** 2))
        if k < K - 1:
            e.append(c2 * _a(j) * _a(j + 1))
    return d, e


def thomas_solve(d, e, sigma, b):
    n = len(d)
    dd = [d[i] - sigma for i in range(n)]
    cp = [mp.mpf(0)] * n
    bp = [mp.mpf(0)] * n
    piv = dd[0] if dd[0] != 0 else mp.mpf(10) ** (-mp.mp.dps + 5)
    cp[0] = e[0] / piv if n > 1 else 0
    bp[0] = b[0] / piv
    for i in range(1, n):
        den = dd[i] - e[i - 1] * cp[i - 1]
        if den == 0:
            den = mp.mpf(10) ** (-mp.mp.dps + 5)
        if i < n - 1:
            cp[i] = e[i] / den
        bp[i] = (b[i] - e[i - 1] * bp[i - 1]) / den
    x = [mp.mpf(0)] * n
    x[-1] = bp[-1]
    for i in range(n - 2, -1, -1):
        x[i] = bp[i] - cp[i] * x[i + 1]
    return x


def eigpair(lam, idx, dps, K=None):
    """index idx within the even sector. Returns (chi, coeffs d_k on p_{2k}) at precision dps."""
    mp.mp.dps = dps
    lam = mp.mpf(lam)
    c2 = (2 * mp.pi * lam ** 2) ** 2
    if K is None:
        K = int(max(120, 1.2 * float(2 * mp.pi * lam ** 2)) + 40 + dps // 2)
    d, e = even_tridiag(c2, K)
    df = np.array([float(x) for x in d]); ef = np.array([float(x) for x in e])
    w, v = eigh_tridiagonal(df, ef, select="i", select_range=(idx, idx))
    x = [mp.mpf(float(t)) for t in v[:, 0]]
    sigma = mp.mpf(float(w[0]))
    for it in range(12):
        y = thomas_solve(d, e, sigma + mp.mpf(10) ** (-(10 + 8 * it)) if it < 2 else sigma, x)
        nrm = mp.sqrt(sum(t * t for t in y))
        x = [t / nrm for t in y]
        Tx = [d[i] * x[i] + (e[i - 1] * x[i - 1] if i > 0 else 0) + (e[i] * x[i + 1] if i < K - 1 else 0) for i in range(K)]
        sigma = sum(a * b for a, b in zip(x, Tx))
        res = mp.sqrt(sum((Tx[i] - sigma * x[i]) ** 2 for i in range(K)))
        if res < mp.mpf(10) ** (-(dps - 8)):
            break
    # sign convention: coefficient of p_0 positive for n = 0, otherwise leave (overall signs cancel in h_lam)
    tail = max(abs(t) for t in x[-8:])
    return sigma, x, res, tail, K


def make_hlam(lam, dps):
    """Return (g, K, info): Legendre coefficients (on p_{2k}) of H = d4_0 h_0 - d0_0 h_4  (zero integral)."""
    chi0, d0, r0, t0, K0 = eigpair(lam, 0, dps)
    chi4, d4, r4, t4, K4 = eigpair(lam, 2, dps, K=K0)
    g = [d4[0] * a - d0[0] * b for a, b in zip(d0, d4)]
    info = dict(chi0=chi0, chi4=chi4, res0=r0, res4=r4, tail0=t0, tail4=t4, K=K0,
                I0=mp.sqrt(2) * mp.mpf(lam) * d0[0], I4=mp.sqrt(2) * mp.mpf(lam) * d4[0])
    return g, K0, info


def eval_H(g, lam, y):
    """H(y) for |y| <= lam via Legendre recurrence (even coefficients only); 0 outside."""
    lam = mp.mpf(lam)
    if abs(y) > lam:
        return mp.mpf(0)
    s = y / lam
    K = len(g)
    # explicit loop (clear version)
    Pjm1, Pj = mp.mpf(1), s       # P_0, P_1
    tot = g[0] * mp.sqrt(mp.mpf(1) / 2)
    for j in range(1, 2 * K - 1):
        Pjp1 = ((2 * j + 1) * s * Pj - j * Pjm1) / (j + 1)
        Pjm1, Pj = Pj, Pjp1       # now Pj = P_{j+1}
        jj = j + 1
        if jj % 2 == 0 and jj // 2 < K:
            tot += g[jj // 2] * mp.sqrt(mp.mpf(2 * jj + 1) / 2) * Pj
    return tot
