"""Fixed-window prime inequality (Daniel, RH_Fixed_Window_and_Real_Axis_Closure, 2026-09-23)
as a no-blowup test, run on zeta, a Euler-product control without pole (Legendre mod 5),
and the Davenport-Heilbronn function (no Euler product, off-line zeros).

S_l(x) = sum_n b(n) n^{-1/2} K(log(n/x)),  K(u) = (1-|u|/l)_+,  b = coefficients of -f'/f.
Explicit formula: S_l(x) = [pole term] - sum_rho Khat(rho-1/2) x^{rho-1/2} + small,
Khat(z) = l (sinh(l z/2)/(l z/2))^2.  Pole term for zeta: (8/l)(cosh(l/2)-1) sqrt(x).
Paper: RH <=> S_{log2}(x) - A sqrt(x) >= -1 for all x >= 2.  Evidence here, nothing proved.
"""
import json, math, os, sys, time
import numpy as np
try:
    from numba import njit
except ImportError:  # smoke tests only
    def njit(**k):
        return lambda f: f
import mpmath as mp

OUT = os.path.join(os.path.dirname(__file__), "results")
os.makedirs(OUT, exist_ok=True)
N = int(float(sys.argv[1])) if len(sys.argv) > 1 else int(2e7)
KAP = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)


@njit(cache=True)
def mangoldt(N):
    lam = np.zeros(N + 1)
    sieve = np.ones(N + 1, np.bool_)
    for p in range(2, N + 1):
        if sieve[p]:
            for k in range(p * p, N + 1, p):
                sieve[k] = False
            lp = math.log(p); q = p
            while q <= N:
                lam[q] = lp
                q *= p
    return lam


@njit(cache=True)
def logderiv_coeffs(N, tab):
    # f = sum a(n) n^-s, a periodic mod 5 via tab, a(1)=1.  -f' = f * (-f'/f):
    # sum_{d|n} b(d) a(n/d) = a(n) log n.
    acc = np.zeros(N + 1)
    for n in range(2, N + 1):
        acc[n] = tab[n % 5] * math.log(n)
    b = np.zeros(N + 1)
    for d in range(2, N + 1):
        b[d] = acc[d]
        bd = b[d]
        if bd != 0.0:
            m = 2
            while d * m <= N:
                acc[d * m] -= bd * tab[m % 5]
                m += 1
    return b


def window(b, ell, xs):
    n = np.arange(len(b), dtype=np.float64)
    w = np.zeros(len(b), np.longdouble)
    w[1:] = b[1:] / np.sqrt(n[1:])
    W = np.cumsum(w); WL = np.cumsum(w * np.log(np.maximum(n, 1)).astype(np.longdouble))
    lo = np.floor(xs * math.exp(-ell)).astype(np.int64)
    mid = np.floor(xs).astype(np.int64)
    hi = np.floor(xs * math.exp(ell)).astype(np.int64)
    lx = np.log(xs).astype(np.longdouble)
    p1 = (W[mid] - W[lo]) * (1 - lx / ell) + (WL[mid] - WL[lo]) / ell
    p2 = (W[hi] - W[mid]) * (1 + lx / ell) - (WL[hi] - WL[mid]) / ell
    return (p1 + p2).astype(np.float64)


def khat(z, ell):
    u = ell * z / 2
    return ell * (mp.sinh(u) / u) ** 2


def dh(s):
    s = mp.mpf(1) * s
    return 5 ** (-s) * (mp.zeta(s, mp.mpf(1) / 5) + KAP * mp.zeta(s, mp.mpf(2) / 5)
                        - KAP * mp.zeta(s, mp.mpf(3) / 5) - mp.zeta(s, mp.mpf(4) / 5))


def main():
    t0 = time.time(); res = {"N": N}
    mp.mp.dps = 30
    guesses = [0.808517 + 85.699348j, 0.650830 + 114.163343j, 0.574356 + 166.479306j, 0.724258 + 176.702461j]
    zeros = []
    for g in guesses:
        r = mp.findroot(dh, mp.mpc(g))
        zeros.append(complex(r))
    res["dh_offline_zeros"] = [[z.real, z.imag, abs(complex(dh(z)))] for z in zeros]

    lam = mangoldt(N)
    tab_dh = np.array([0.0, 1.0, KAP, -KAP, -1.0])
    tab_l5 = np.array([0.0, 1.0, -1.0, -1.0, 1.0])
    bdh = logderiv_coeffs(N, tab_dh)
    bl5 = logderiv_coeffs(N, tab_l5)
    lp = np.zeros(N + 1); lp[1:] = np.arange(1, N + 1)
    res["chi5_check_maxerr"] = float(np.max(np.abs(bl5[2:] - lam[2:] * tab_l5[np.arange(2, N + 1) % 5])))
    # growth of DH coefficients: sup |b(n)| / n^{sigma-1} test via max over dyadic blocks
    grow = []
    k = 2
    while 2 * k <= N:
        blk = np.abs(bdh[k:2 * k]); grow.append([k, float(blk.max()), float(np.mean(blk))])
        k *= 2
    res["dh_coeff_growth_dyadic"] = grow
    res["t_coeffs"] = time.time() - t0

    for ell in [math.log(2), 1.0, 2.0]:
        X = int(N / math.exp(ell)) - 2
        xs = np.arange(2, X + 1, dtype=np.float64)
        A = (8 / ell) * (math.cosh(ell / 2) - 1)
        out = {"ell": ell, "A": A, "Xmax": X}
        dz = window(lam, ell, xs) - A * np.sqrt(xs)
        i = int(np.argmin(dz)); out["zeta_min_deficit"] = [float(xs[i]), float(dz[i])]
        out["zeta_max"] = float(dz.max())
        for name, b in [("chi5", bl5), ("dh", bdh)]:
            s = window(b, ell, xs)
            blocks = []
            k = 2
            while 2 * k <= X:
                seg = s[k - 2:2 * k - 2]
                blocks.append([k, float(seg.min()), float(seg.max())])
                k *= 2
            out[name + "_dyadic_minmax"] = blocks
            if name == "zeta":
                pass
        # zeta dyadic too
        blocks = []; k = 2
        while 2 * k <= X:
            seg = dz[k - 2:2 * k - 2]; blocks.append([k, float(seg.min()), float(seg.max())]); k *= 2
        out["zeta_dyadic_minmax"] = blocks
        # prediction from the 4 off-line DH zeros (and conjugates): -2 Re sum Khat(rho-1/2) x^{rho-1/2}
        grid = np.unique(np.geomspace(1e3, X, 4000).astype(np.int64)).astype(np.float64)
        sdh = window(bdh, ell, grid)
        pred = np.zeros_like(grid)
        amps = []
        for z in zeros:
            kz = complex(khat(mp.mpc(z) - 0.5, ell)); amps.append([abs(kz), z.real - 0.5])
            pred += -2 * np.real(kz * np.exp((z - 0.5) * np.log(grid)))
        resid = sdh - pred
        out["dh_offline_amp_|Khat|_and_exponent"] = amps
        out["dh_pred_rms"] = float(np.sqrt(np.mean(pred ** 2)))
        out["dh_actual_rms"] = float(np.sqrt(np.mean(sdh ** 2)))
        out["dh_resid_rms"] = float(np.sqrt(np.mean(resid ** 2)))
        out["dh_corr_actual_pred"] = float(np.corrcoef(sdh, pred)[0, 1])
        # x at which first off-line term reaches amplitude 1
        a0, e0 = amps[0]
        out["dh_x_where_offline_term_reaches_1"] = float((1 / (2 * a0)) ** (1 / e0))
        res["ell=%.4f" % ell] = out
        print(json.dumps(out)[:1500], flush=True)
    res["t_total"] = time.time() - t0
    json.dump(res, open(os.path.join(OUT, "fixed_window.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
