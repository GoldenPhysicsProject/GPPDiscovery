"""Claude (chat seat) - is the CCM ground state Riemann's Phi?
Profile f(x) = sum_j (-1)^j xi_j cos(2 pi j x / L) on [-L/2, L/2] (log coordinate), whose window Fourier transform is xihat.
Phi(u) = sum_n (2 pi^2 n^4 e^{9u/2} - 3 pi n^2 e^{5u/2}) exp(-pi n^2 e^{2u}), even; Xi(t) = int Phi(u) cos(tu) du.
Reports per (lambda, N): relative L2 distance of normalized profiles, moments <x^2>, <x^4>, <x^6> vs Phi, positivity,
and the Rayleigh quotient of Phi (projected on the N modes) against the ground energy eps."""
import mpmath as mp, json, argparse, time
from ccm import build
def Phi(u):
    u = abs(u); return sum((2*mp.pi**2*n**4*mp.exp(9*u/2)-3*mp.pi*n**2*mp.exp(5*u/2))*mp.exp(-mp.pi*n*n*mp.exp(2*u)) for n in range(1, 12))
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--lam", required=True); ap.add_argument("--Ns", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args(); t0 = time.time()
    with open(a.out, "w") as fo:
        for N in [int(x) for x in a.Ns.split(",")]:
            T, L, idx, ks = build(a.lam, N, 60+2*N); mp.mp.dps = 60+2*N
            E, Q = mp.eigsy(T); k = min(range(len(idx)), key=lambda i: E[i]); xi = [Q[r, k] for r in range(len(idx))]
            order = sorted(range(len(idx)), key=lambda i: E[i]); gap = E[order[1]]/E[k]
            f = lambda x: sum(c*(-1)**j*mp.cos(2*mp.pi*j*x/L) for c, j in zip(xi, idx))
            I = mp.quad(f, [-L/2, 0, L/2]); IP = mp.quad(Phi, [-L/2, 0, L/2])
            mom = lambda g, n, J: mp.quad(lambda x: x**n*g(x), [-L/2, 0, L/2])/J
            l2 = mp.sqrt(mp.quad(lambda x: (f(x)/I-Phi(x)/IP)**2, [-L/2, 0, L/2])/mp.quad(lambda x: (Phi(x)/IP)**2, [-L/2, 0, L/2]))
            c = mp.matrix([mp.quad(lambda x: Phi(x)*(-1)**j*mp.cos(2*mp.pi*j*x/L), [-L/2, 0, L/2])/L for j in idx])
            ray = (c.T*T*c)[0]/(c.T*c)[0]; ov = abs((c.T*mp.matrix(xi))[0])/mp.sqrt((c.T*c)[0])
            xs = [-L/2+L*i/120 for i in range(121)]; mn = min(f(x)/I for x in xs)
            rec = {"lam": a.lam, "N": N, "eps": mp.nstr(E[k], 8), "gap_ratio": mp.nstr(gap, 6), "relL2_f_vs_Phi": mp.nstr(l2, 6),
                   "moments_f": [mp.nstr(mom(f, n, I), 10) for n in (2, 4, 6)], "moments_Phi": [mp.nstr(mom(Phi, n, IP), 10) for n in (2, 4, 6)],
                   "min_profile": mp.nstr(mn, 4), "rayleigh_Phi": mp.nstr(ray, 8), "rayleigh_over_eps": mp.nstr(ray/E[k], 6),
                   "overlap_Phi_ground": mp.nstr(ov, 12), "sec": round(time.time()-t0)}
            fo.write(json.dumps(rec)+"\n"); fo.flush(); print(json.dumps(rec), flush=True); mp.mp.dps = 30
