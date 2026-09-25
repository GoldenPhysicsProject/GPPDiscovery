"""Claude (chat seat) - R_inf(lambda): fixed-lambda behaviour of the compensated-Weyl remainder.
R_N(z) = C(z) - G(z),  C = B/A + (L/pi^2) sin^2(pi z) A'/A,  G = b_W(z) - (eps+L) sin(2 pi z)/(2 pi).
Relative |R|/|G| at z in {0.5, 1.3+0.5i, 3.2+i}, every N (odd and even) up to Nmax.
Builder: CCM, Zeta Spectral Triples, arXiv:2511.22755 (ccm.py)."""
import mpmath as mp, json, sys, time, argparse
from ccm import build

def setup(lamstr, N, dps):
    mp.mp.dps = dps
    lam = mp.mpf(lamstr); L = 2*mp.log(lam); pi = mp.pi
    h   = lambda s: mp.sin(2*pi*s)/(2*pi)
    psi = lambda y, s: mp.sin(2*pi*(1-y/L)*s)/pi
    pL  = lambda s: (L/pi**2)*s*(mp.cosh(L/2)-mp.cos(2*pi*s))/(s**2+(L/(4*pi))**2)
    K = int(mp.floor(lam**2 + mp.mpf('1e-30')))
    def vm(k):
        for p in range(2, k+1):
            if k % p == 0:
                m = k
                while m % p == 0: m //= p
                return mp.log(p) if m == 1 else mp.mpf(0)
        return mp.mpf(0)
    PK = [(k, vm(k)) for k in range(2, K+1) if vm(k) != 0]
    bP  = lambda s: sum(c/mp.sqrt(k)*psi(mp.log(k), s) for k, c in PK)
    rho = lambda y: mp.exp(y/2)/(mp.exp(y)-mp.exp(-y))
    r   = lambda y: 1/(mp.exp(y)-mp.exp(-y))
    c0 = mp.log(4*pi)+mp.euler; tau = mp.log((mp.exp(L)+1)/(mp.exp(L)-1))/2
    def bR(s):
        I = mp.quad(lambda y: rho(y)*psi(y, s) - r(y)*2*h(s), [0, L])
        return I + (c0-2*tau)*h(s)
    bW = lambda s: pL(s) - bR(s) - bP(s)
    T, Lb, idx, ks = build(lamstr, N, dps)
    return L, bW, T, idx

ZS = [mp.mpf('0.5'), mp.mpc('1.3', '0.5'), mp.mpc('3.2', '1')]

def rel_R(lamstr, N, dps):
    L, bW, T, idx = setup(lamstr, N, dps); n = len(idx); pi = mp.pi
    E, Q = mp.eigsy(T); order = sorted(range(n), key=lambda i: E[i]); k = order[0]
    eps, e1 = E[k], E[order[1]]
    xi = [Q[r, k] for r in range(n)]; s = sum(xi); xi = [x/s for x in xi]
    beta = [idx[a]*T[a, idx.index(0)] for a in range(n)]
    A  = lambda z: sum(xi[j]/(idx[j]-z) for j in range(n))
    Ap = lambda z: sum(xi[j]/(idx[j]-z)**2 for j in range(n))
    B  = lambda z: sum(beta[j]*xi[j]/(idx[j]-z) for j in range(n))
    rL = lambda z: (L/pi**2)*mp.sin(pi*z)**2
    out = []
    for z in ZS:
        G = bW(z) - (eps+L)*mp.sin(2*pi*z)/(2*pi)
        C = B(z)/A(z) + rL(z)*Ap(z)/A(z)
        out.append(abs(C-G)/abs(G))
    even = max(abs(xi[r]-xi[n-1-r]) for r in range(n))
    return eps, e1, even, out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--lam", required=True); ap.add_argument("--Ns", default="6-30,32,34,36,38,40,44,48")
    ap.add_argument("--dps-base", type=int, default=60); ap.add_argument("--dps-slope", type=int, default=2)
    ap.add_argument("--budget-min", type=float, default=300); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    Ns = []
    for part in a.Ns.split(","):
        if "-" in part: lo, hi = map(int, part.split("-")); Ns += list(range(lo, hi+1))
        else: Ns.append(int(part))
    t0 = time.time()
    with open(a.out, "w") as f:
        for N in Ns:
            if (time.time()-t0)/60 > a.budget_min: print("budget reached before N =", N); break
            dps = a.dps_base + a.dps_slope*N; t1 = time.time()
            eps, e1, even, rel = rel_R(a.lam, N, dps)
            rec = {"lam": a.lam, "N": N, "parity": "odd" if N % 2 else "even", "dps": dps,
                   "eps": mp.nstr(eps, 12), "gap_ratio": mp.nstr(e1/eps, 6), "even_defect": mp.nstr(even, 3),
                   "rel_R": [mp.nstr(x, 6) for x in rel], "max_rel_R": mp.nstr(max(rel), 6),
                   "sec": round(time.time()-t1, 1)}
            f.write(json.dumps(rec)+"\n"); f.flush(); print(json.dumps(rec)); sys.stdout.flush()
