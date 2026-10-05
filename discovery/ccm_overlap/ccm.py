# Daniel/Claude chat output (ccm_zeta_spectral_triple.py), copied verbatim as ccm.py for the trial-vector overlap test.
"""Connes-Consani-Moscovici 'Zeta Spectral Triples' (arXiv:2511.22755), built from the paper's formulas.
Matrix of the truncated Weil form QW_lambda^N in the basis V_n, |n|<=N, L = 2 log(lambda):
   tau = W02 - WR - sum_p Wp          (eq. 3.10, 3.13, Lemma 4.1, eq. 4.3, eq. 3.15)
Smallest eigenvector xi (should be even), then zeros of
   xihat(z) = 2 L^{-1/2} sin(zL/2) sum_j xi_j/(z - 2 pi j/L)      (eq. 5.25)
approximate the zeta ordinates gamma_n.  Theorem 5.10: all zeros of xihat are real."""
import mpmath as mp, sys, time

def build(lam, N, dps):
    mp.mp.dps = dps
    lam = mp.mpf(lam); L = 2*mp.log(lam); pi = mp.pi
    idx = list(range(-N, N+1))
    rho = lambda x: mp.exp(x/2)/(mp.exp(x)-mp.exp(-x))
    # von Mangoldt terms for 1 < k <= e^L = lam^2
    K = int(mp.floor(lam**2 + mp.mpf('1e-30')))
    def vm(k):
        for p in range(2, k+1):
            if k % p == 0:
                m = k
                while m % p == 0: m //= p
                return mp.log(p) if m == 1 else mp.mpf(0)
        return mp.mpf(0)
    primes_terms = [(k, vm(k)) for k in range(2, K+1) if vm(k) != 0]
    def omega(n, m, y):                         # q(U_n,U_m)(y), Lemma 2.3, y in [0,L]
        if n != m: return (mp.sin(2*pi*m*y/L) - mp.sin(2*pi*n*y/L))/(pi*(n-m))
        return 2*(L-y)*mp.cos(2*pi*n*y/L)/L
    # archimedean pieces: alpha(n) for off-diagonal, diagonal integral directly from (3.15)
    alpha = {n: mp.quad(lambda x: mp.sin(2*pi*n*x/L)*rho(x), [0, L])/pi for n in range(0, N+1)}
    for n in range(1, N+1): alpha[-n] = -alpha[n]
    tail = mp.log((mp.exp(L)+1)/(mp.exp(L)-1))/2
    diagR = {}
    for n in range(0, N+1):
        integ = mp.quad(lambda y: (mp.exp(y/2)*2*(1-y/L)*mp.cos(2*pi*n*y/L) - 2)/(mp.exp(y)-mp.exp(-y)), [0, L])
        diagR[n] = (mp.log(4*pi) + mp.euler) + integ - 2*tail
        diagR[-n] = diagR[n]
    T = mp.matrix(2*N+1, 2*N+1)
    for a, n in enumerate(idx):
        for b, m in enumerate(idx):
            w02 = 32*L*mp.sinh(L/4)**2*(L**2 - 16*pi**2*m*n)/((L**2+16*pi**2*m**2)*(L**2+16*pi**2*n**2))
            wr = diagR[n] if n == m else (alpha[m]-alpha[n])/(n-m)
            wp = sum(c*mp.mpf(k)**(-mp.mpf(1)/2)*omega(n, m, mp.log(k)) for k, c in primes_terms)
            T[a, b] = w02 - wr - wp
    return T, L, idx, [k for k,_ in primes_terms]

def analyze(lam, N, dps, nzeros=8, verbose=True):
    t0 = time.time()
    T, L, idx, ks = build(lam, N, dps)
    E, Q = mp.eigsy(T)
    order = sorted(range(len(E)), key=lambda i: E[i])
    e0, e1 = E[order[0]], E[order[1]]
    xi = [Q[r, order[0]] for r in range(len(idx))]
    even = max(abs(xi[r] - xi[len(idx)-1-r]) for r in range(len(idx)))
    xihat = lambda z: 2*mp.sin(z*L/2)*sum(x/(z - 2*mp.pi*j/L) for x, j in zip(xi, idx))/mp.sqrt(L)
    gam = [mp.zetazero(k).imag for k in range(1, nzeros+1)]
    out = []
    for g in gam:
        try:
            z = mp.findroot(xihat, g)
            out.append((g, z, abs(z-g)))
        except Exception as ex:
            out.append((g, None, None))
    if verbose:
        print(f"lambda={lam}  N={N}  dps={dps}  L={mp.nstr(L,6)}  prime powers k<=lambda^2: {ks}  ({time.time()-t0:.0f}s)")
        print(f"  smallest eigenvalue eps_N = {mp.nstr(e0,5)}   next = {mp.nstr(e1,5)}   gap ratio = {mp.nstr(e1/e0,4)}")
        print(f"  ground state even (xi_j = xi_-j): max |xi_j - xi_-j| = {mp.nstr(even,3)}")
        for g, z, err in out:
            print(f"    gamma = {mp.nstr(g,15):>18}   zero of xihat = {mp.nstr(z,15) if z else '---':>18}   |error| = {mp.nstr(err,3) if err is not None else '---'}")
    return e0, e1, even, out

if __name__ == "__main__":
    lam, N, dps = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    analyze(lam, N, dps)
