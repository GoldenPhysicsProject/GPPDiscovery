#!/usr/bin/env python3
"""One (c, N, T, dps) point of the truncated-Weil ground-state scan.

Emits a single JSON line. Designed to be one cell of a GitHub Actions matrix
so the whole scan runs in parallel wall-clock time instead of serially.

lambda_min(c) is the smallest eigenvalue of the Connes-van Suijlekom
Proposition 4.1 Galerkin matrix Q(c) restricted to the even sector.
RH  <=>  Weil positivity  <=>  lambda_min >= 0 for every cutoff.

The question this scan exists to answer is NOT whether lambda_min > 0 (it is,
assuming RH) but HOW FAST it decays: is log lambda_min linear in log c, and if
so what is the constant? Two-point estimates gave a natural-log slope near -64;
-2*gamma_1 = -28.27 and -4*gamma_1 = -56.55. Six converged points settle it.
"""
import argparse, json, platform, sys, time
import mpmath as mp
import connes_cvs as cc

p = argparse.ArgumentParser()
p.add_argument("--c", type=int, required=True)
p.add_argument("--N", type=int, required=True)
p.add_argument("--T", type=int, default=300)
p.add_argument("--dps", type=int, default=90)
p.add_argument("--out", default="result.json")
a = p.parse_args()

t0 = time.time()
mp.mp.dps = a.dps
Q = cc.build_galerkin_matrix(a.c, N=a.N, T=a.T, dps=a.dps)
lam, vec = cc.compute_ground_state(Q)

rec = dict(
    c=a.c, N=a.N, T=a.T, dps=a.dps,
    dim=2 * a.N + 1,
    a=float(mp.log(a.c)),
    lam=mp.nstr(lam, 25),
    lam_abs_log10=float(mp.log(abs(lam), 10)) if lam != 0 else None,
    sign=int(mp.sign(lam)),
    seconds=round(time.time() - t0, 1),
    python=platform.python_version(),
    mpmath=mp.__version__,
    connes_cvs=getattr(cc, "__version__", "?"),
)
with open(a.out, "w") as f:
    json.dump(rec, f)
print(json.dumps(rec), flush=True)
