"""Claude (chat seat) - the 1D -> 2D ladder, numerically.
A  radial Wigner delay split: Re(xi'/xi) = gamma + pole + primitive-prime P'(s) + (m>=2 sector),
   P'(s) = sum_k mu(k) (zeta'/zeta)(ks)  [Moebius inversion of log zeta].  Does gamma+pole+primitive carry positivity alone?
B  Hecke L(s, lambda^{4J}) of Q(i), lambda(alpha) = alpha/|alpha|: spin-J zero spectra on Re s = 1/2 (celestial helicity ladder).
   Lambda(s) = pi^{-s} Gamma(s + 2J) L(s), L = (1/4) sum_{alpha != 0} lambda^{4J}(alpha) N(alpha)^{-s}; root number fitted and checked.
C  Picard S-matrix phi(s) = Lambda_K(s-1)/Lambda_K(s), K = Q(i): (i) lattice check sum_c Phi(c) N(c)^{-s} = zeta_K(s-1)/zeta_K(s);
   (ii) unitarity and Wigner delay on s = 1 + i nu; (iii) resonances = zeros of zeta and L(chi_{-4}) at depth 1/2 below the unitary line.
   Hand result: phi(s) = -s/(2-s) xi_K(2-s)/xi_K(s); delay = -2/(1+nu^2) [vacuum bound state] + 2 Re(xi_K'/xi_K)(1+i nu) [> 0 unconditionally]."""
import mpmath as mp, json, argparse, time, math

def mobius(n):
    r, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r

def partA(sigmas, tmax, nt, K=25):
    mp.mp.dps = 20; out = []
    mus = {k: mobius(k) for k in range(1, K+1)}
    for sg in sigmas:
        rows = []
        for i in range(nt):
            t = mp.mpf(tmax)*(i+0.5)/nt; s = mp.mpc(sg, t)
            zl = lambda w: mp.zeta(w, derivative=1)/mp.zeta(w)
            full = mp.re(zl(s))
            prim = mp.re(sum(mus[k]*zl(k*s) for k in range(1, K+1) if mus[k] != 0))
            g = mp.re(mp.digamma(s/2))/2 - mp.log(mp.pi)/2; pl = mp.re(1/s + 1/(s-1))
            rows.append((float(t), float(full), float(prim), float(g), float(pl)))
        tot = [r[1]+r[3]+r[4] for r in rows]; gpp = [r[2]+r[3]+r[4] for r in rows]; m2 = [r[1]-r[2] for r in rows]
        decisive = sum(1 for a, b in zip(tot, gpp) if a > 0 and b <= 0)
        i0 = min(range(nt), key=lambda i: gpp[i])
        out.append({"sigma": sg, "min_total": min(tot), "min_gamma_pole_primitive": gpp[i0], "t_at_min_gpp": rows[i0][0],
                    "frac_gpp_negative": sum(1 for x in gpp if x <= 0)/nt, "m2_sector_min": min(m2), "m2_sector_max": max(m2),
                    "points_where_m2_sector_is_decisive": decisive})
        print(json.dumps(out[-1]), flush=True)
    return out

def hecke_coeffs(J, Nmax):
    a = [mp.mpf(0)]*(Nmax+1); R = int(math.isqrt(Nmax))+1
    for x in range(-R, R+1):
        for y in range(-R, R+1):
            n = x*x+y*y
            if 0 < n <= Nmax:
                a[n] += mp.re((mp.mpc(x, y)/mp.sqrt(n))**(4*J))/4
    return a

def hecke_Lambda(J, a, s, eps):
    k = 2*J; tot = mp.mpf(0)
    for n in range(1, len(a)):
        if a[n] == 0: continue
        x = mp.pi*n
        tot += a[n]*((x)**(-s)*mp.gammainc(s+k, x) + eps*(x)**(-(1-s))*mp.gammainc(1-s+k, x))
    return tot

def partB(Js, tmax, dt, Nmax):
    mp.mp.dps = 25; out = []
    for J in Js:
        a = hecke_coeffs(J, Nmax); s0 = mp.mpf(4)
        direct = mp.pi**(-s0)*mp.gamma(s0+2*J)*sum(a[n]*mp.mpf(n)**(-s0) for n in range(1, Nmax+1))
        eps = min([1, -1], key=lambda e: abs(hecke_Lambda(J, a, s0, e)-direct))
        chk = abs(hecke_Lambda(J, a, s0, eps)-direct)/abs(direct)
        sfe = mp.mpc('0.3', '7.1'); fe = abs(hecke_Lambda(J, a, sfe, eps)-eps*hecke_Lambda(J, a, 1-sfe, eps))
        Z = lambda t: hecke_Lambda(J, a, mp.mpf(1)/2+1j*t, eps)
        ts = [dt*(i+1) for i in range(int(tmax/dt))]; vals = [Z(t) for t in ts]
        zeros = []; imag_max = max(abs(mp.im(v)) / (abs(v)+mp.mpf('1e-300')) for v in vals)
        for i in range(len(ts)-1):
            if mp.re(vals[i])*mp.re(vals[i+1]) < 0:
                zeros.append(float(mp.findroot(lambda t: mp.re(Z(t)), (ts[i], ts[i+1]), solver='anderson')))
        rec = {"J": J, "spin_4J": 4*J, "root_number": eps, "direct_check_rel": mp.nstr(chk, 3), "func_eq_residual": mp.nstr(fe, 3),
               "max_rel_imag_on_line": mp.nstr(imag_max, 3), "zeros": [round(z, 8) for z in zeros], "n_zeros_below_tmax": len(zeros)}
        out.append(rec); print(json.dumps(rec), flush=True)
    return out

def partC(nus, s_checks, Cmax):
    mp.mp.dps = 25; res = {}
    zK = lambda s: mp.zeta(s)*mp.dirichlet(s, [0, 1, 0, -1])
    # Phi(c) = |(Z[i]/c)^x| via prime-ideal factorisation of the norm
    def gauss_Phi(x, y):
        n = x*x+y*y; val = mp.mpf(n); m = n; p = 2
        primes = []
        while p*p <= m:
            if m % p == 0:
                primes.append(p)
                while m % p == 0: m //= p
            p += 1
        if m > 1: primes.append(m)
        for p in primes:
            if p == 2: val *= mp.mpf(1)/2
            elif p % 4 == 3: val *= (1-mp.mpf(1)/(p*p))
            else:
                a0 = next(a for a in range(1, p) if (a*a+1) % p == 0)
                d1 = (x + a0*y) % p == 0
                d2 = (x - a0*y) % p == 0
                if d1: val *= (1-mp.mpf(1)/p)
                if d2: val *= (1-mp.mpf(1)/p)
        return val
    checks = []
    for s in s_checks:
        s = mp.mpf(s); S = mp.mpf(0)
        for x in range(1, Cmax+1):          # one representative per unit class: x > 0, y >= 0
            for y in range(0, Cmax+1):
                if x*x+y*y <= Cmax*Cmax: S += gauss_Phi(x, y)*mp.mpf(x*x+y*y)**(-s)
        checks.append({"s": float(s), "lattice_sum": mp.nstr(S, 10), "zetaK(s-1)/zetaK(s)": mp.nstr(zK(s-1)/zK(s), 10)})
        print(json.dumps(checks[-1]), flush=True)
    res["lattice_check"] = checks
    LK = lambda s: mp.pi**(-s)*mp.gamma(s)*zK(s)
    xK = lambda s: s*(s-1)*LK(s)
    unit = []; delay = []; arith = []; split_err = []
    for nu in nus:
        s = mp.mpc(1, nu); S = LK(2-s)/LK(s); unit.append(float(abs(abs(S)-1)))
        d = 2*mp.re(mp.diff(LK, s)/LK(s)); delay.append(float(d))
        ar = 2*mp.re(mp.diff(xK, s)/xK(s)); arith.append(float(ar))
        split_err.append(float(abs(d - (ar - 2/(1+nu**2)))))
    res["unitarity_max_dev"] = max(unit); res["wigner_delay_min"] = min(delay); res["wigner_delay_at_nu0"] = delay[0]
    res["arithmetic_delay_min"] = min(arith); res["vacuum_arith_split_max_err"] = max(split_err)
    res["wigner_delay_samples"] = [(round(float(n), 2), round(d, 6), round(ar, 6)) for n, d, ar in list(zip(nus, delay, arith))[::max(1, len(nus)//25)]]
    Lz = []
    mp.mp.dps = 20
    Zchi = lambda t: mp.re((4/mp.pi)**((mp.mpf(1)/2+1j*t)/2)*mp.gamma((mp.mpf(3)/2+1j*t)/2)*mp.dirichlet(mp.mpf(1)/2+1j*t, [0, 1, 0, -1]))
    ts = [0.05*(i+1) for i in range(800)]; v = [Zchi(t) for t in ts]
    for i in range(len(ts)-1):
        if v[i]*v[i+1] < 0: Lz.append(round(float(mp.findroot(Zchi, (ts[i], ts[i+1]), solver='anderson')), 8))
    res["resonances"] = {"zeta_zeros_t<40": [round(float(mp.zetazero(k).imag), 8) for k in range(1, 7)],
                         "L_chi4_zeros_t<40_on_line": Lz, "resonance_depth_below_unitary_line": 0.5}
    print(json.dumps(res), flush=True)
    return res

if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--part", default="ABC"); ap.add_argument("--small", action="store_true")
    ap.add_argument("--out", required=True); a = ap.parse_args(); t0 = time.time(); R = {}
    if "A" in a.part: R["A"] = partA([0.55, 0.6, 0.75, 0.9], 12 if a.small else 100, 24 if a.small else 1000, 12 if a.small else 25)
    if "B" in a.part: R["B"] = partB([1, 2] if a.small else [1, 2, 3, 4, 5, 6], 12 if a.small else 30, 0.25 if a.small else 0.05, 80 if a.small else 300)
    if "C" in a.part: R["C"] = partC([mp.mpf(i)/10 for i in range(1, 60 if a.small else 400)], [3, 4], 12 if a.small else 150)
    R["sec"] = round(time.time()-t0); json.dump(R, open(a.out, "w"), indent=1, default=str)
