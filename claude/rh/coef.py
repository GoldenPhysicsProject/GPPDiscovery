"""Claude (chat seat) - Polya-Jensen route: log-Taylor coefficients of the CCM xihat_lambda vs Xi at t = 0.
f_lambda = exp(-a_lambda t^2) xihat_lambda is Laguerre-Polya when a_lambda >= 0 (xihat real-rooted, CCM Thm 5.10).
Reports: a_lambda = sigma_1(Xi) - sigma_1(xihat) (must be >= 0), and relative differences of the t^{2k} log-coefficients, k = 1..6."""
import mpmath as mp, json, argparse, time
from ccm import build
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--lam", required=True); ap.add_argument("--Ns", required=True)
    ap.add_argument("--out", required=True); a = ap.parse_args(); K = 6; t0 = time.time()
    mp.mp.dps = 60
    xi = lambda s: s*(s-1)/2*mp.pi**(-s/2)*mp.gamma(s/2)*mp.zeta(s)
    cX = mp.taylor(lambda t: mp.log(mp.re(xi(mp.mpf(1)/2+1j*t))), 0, 2*K)
    tgt = [cX[2*k] for k in range(1, K+1)]
    Ns = []
    for part in a.Ns.split(","):
        if ":" in part: lo, hi, st = map(int, part.split(":")); Ns += list(range(lo, hi+1, st))
        else: Ns.append(int(part))
    with open(a.out, "w") as f:
        for N in Ns:
            T, L, idx, ks = build(a.lam, N, 60+2*N); mp.mp.dps = 60+2*N
            E, Q = mp.eigsy(T); k0 = min(range(len(idx)), key=lambda i: E[i]); x = [Q[r, k0] for r in range(len(idx))]
            lxh = lambda z: mp.log(2*mp.sin(z*L/2)*sum(c/(z-2*mp.pi*j/L) for c, j in zip(x, idx))/mp.sqrt(L))
            c = mp.taylor(lxh, mp.mpf('1e-40'), 2*K)
            rel = [(c[2*k]-tgt[k-1])/tgt[k-1] for k in range(1, K+1)]
            rec = {"lam": a.lam, "N": N, "eps": mp.nstr(E[k0], 8), "a_lambda": mp.nstr(c[2]-tgt[0], 10),
                   "rel_logcoef_t2_to_t12": [mp.nstr(r, 6) for r in rel], "sec": round(time.time()-t0)}
            f.write(json.dumps(rec)+"\n"); f.flush(); print(json.dumps(rec), flush=True); mp.mp.dps = 60
