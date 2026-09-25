"""Claude (chat seat) - radial Higgs stiffness and field-vs-current test.
Part A: d/dsigma log|xi(sigma+it)| = Re(xi'/xi) split into channels
        zeta: Re(zeta'/zeta), gamma: Re(psi(s/2))/2 - log(pi)/2, pole: Re(1/s + 1/(s-1)).
        Lagarias: RH iff total > 0 for all sigma > 1/2.  Per sigma: min total, min zeta channel,
        min relative margin total/(|zeta|+|gamma|+|pole|) and where it occurs.
Part B: connected field  F_x(t) = sum_{n<=x} n^{-s} - int_1^x u^{-s} du          (RH-blind: -> zeta(s)+1/(1-s))
        connected current J_x(t) = sum_{n<=x} Lambda(n) n^{-s} - int_1^x u^{-s} du (zeros enter via x^{rho-s}/(rho-s))
        at s = 1/2 + it, t in [1, 60]; sup norms and values at gamma_1 as x grows."""
import mpmath as mp, numpy as np, json, time, argparse

def part_a(sigmas, tmax, nt):
    mp.mp.dps = 25
    ts = [mp.mpf(tmax)*(k+0.5)/nt for k in range(nt)]
    res = []
    for sg in sigmas:
        sg = mp.mpf(sg); rows = []
        for t in ts:
            s = mp.mpc(sg, t)
            z = mp.zeta(s); zp = mp.zeta(s, derivative=1)
            ch_z = mp.re(zp/z); ch_g = mp.re(mp.digamma(s/2))/2 - mp.log(mp.pi)/2; ch_p = mp.re(1/s + 1/(s-1))
            rows.append((float(t), float(ch_z), float(ch_g), float(ch_p)))
        a = np.array(rows); tot = a[:,1]+a[:,2]+a[:,3]
        i_tot = int(np.argmin(tot)); rel = tot/(np.abs(a[:,1])+np.abs(a[:,2])+np.abs(a[:,3]))
        i_r = int(np.argmin(rel))
        res.append({"sigma": float(sg), "min_total": tot[i_tot], "t_at_min_total": a[i_tot,0],
                    "min_zeta_channel": float(a[:,1].min()), "frac_t_zeta_negative": float((a[:,1] < 0).mean()),
                    "min_relative_margin": float(rel[i_r]), "t_at_min_relative_margin": a[i_r,0],
                    "min_gamma_channel": float(a[:,2].min()), "min_pole_channel": float(a[:,3].min())})
        print(json.dumps(res[-1]), flush=True)
    return res

def part_b(xs, t_lo, t_hi, nt, chunk=4):
    X = max(xs)
    n = np.arange(1, X+1, dtype=np.float64)
    lam = np.zeros(X+1)                              # von Mangoldt sieve
    isp = np.ones(X+1, bool); isp[:2] = False
    for p in range(2, int(X**0.5)+1):
        if isp[p]: isp[p*p::p] = False
    for p in np.nonzero(isp)[0]:
        pk = p
        while pk <= X: lam[pk] = np.log(p); pk *= p
    lam = lam[1:]
    ts = np.linspace(t_lo, t_hi, nt)
    g1 = 14.134725141734693
    ts = np.sort(np.append(ts, g1))
    bounds = [0] + sorted(xs)
    F = np.zeros((len(xs), len(ts)), complex); J = np.zeros_like(F)
    logn = np.log(n); w = n**-0.5
    for c0 in range(0, len(ts), chunk):
        tt = ts[c0:c0+chunk]
        accF = np.zeros(len(tt), complex); accJ = np.zeros(len(tt), complex)
        for j in range(len(xs)):
            lo, hi = bounds[j], bounds[j+1]
            ph = np.exp(-1j*np.outer(tt, logn[lo:hi]))*w[lo:hi]
            accF += ph.sum(1); accJ += (ph*lam[lo:hi]).sum(1)
            x = float(bounds[j+1]); s = 0.5+1j*tt
            cond = (x**(1-s)-1)/(1-s)
            F[j, c0:c0+len(tt)] = accF - cond; J[j, c0:c0+len(tt)] = accJ - cond
    ig = int(np.argmin(abs(ts-g1)))
    mp.mp.dps = 20
    zref = [complex(mp.zeta(mp.mpc(0.5, t)) + 1/(1-mp.mpc(0.5, t))) for t in ts[::max(1, len(ts)//40)]]
    out = []
    for j, x in enumerate(xs):
        dF = np.abs(F[j, ::max(1, len(ts)//40)] - np.array(zref)).max()
        out.append({"x": x, "L": float(np.log(x)), "sup_field": float(np.abs(F[j]).max()),
                    "sup_current": float(np.abs(J[j]).max()), "t_sup_current": float(ts[np.argmax(np.abs(J[j]))]),
                    "field_at_gamma1": abs(F[j, ig]), "current_at_gamma1": abs(J[j, ig]),
                    "field_minus_limit_max": float(dF)})
        print(json.dumps(out[-1]), flush=True)
    return out

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--sigmas", default="0.5,0.51,0.55,0.6,0.75,1.0")
    ap.add_argument("--tmax", type=float, default=120); ap.add_argument("--nt", type=int, default=2400)
    ap.add_argument("--xs", default="1000,10000,100000,1000000,10000000")
    ap.add_argument("--bt", type=int, default=3000); ap.add_argument("--out", required=True)
    a = ap.parse_args(); t0 = time.time()
    A = part_a([float(x) for x in a.sigmas.split(",")], a.tmax, a.nt)
    B = part_b([int(x) for x in a.xs.split(",")], 1.0, 60.0, a.bt)
    json.dump({"part_a_radial_stiffness": A, "part_b_field_vs_current": B, "sec": round(time.time()-t0, 1)},
              open(a.out, "w"), indent=1)
