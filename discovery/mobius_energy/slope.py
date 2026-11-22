"""Measured slope of Q_c(log N) in log N, c = 1/zeta(1/2), by least squares on a log-spaced grid (NUMERICAL).
Also reports the amplitude of B = S - c and the slope for windows, to expose the almost-periodic wobble.
Usage: python3 -I slope.py [NMAX]"""
import sys, json, math
import numpy as np
from mpmath import mp, zeta
NMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10_000_000
mp.dps = 30
c0 = float(1 / zeta(mp.mpf(1) / 2))
mu = np.ones(NMAX + 1, dtype=np.int8); mu[0] = 0
isp = np.ones(NMAX + 1, dtype=bool); isp[:2] = False
for i in range(2, int(NMAX ** 0.5) + 1):
    if isp[i]: isp[i * i::i] = False
for p in np.flatnonzero(isp):
    mu[p::p] *= -1
    if p * p <= NMAX: mu[p * p::p * p] = 0
n = np.arange(NMAX + 1, dtype=np.float64)
a = np.zeros(NMAX + 1); a[1:] = mu[1:] / np.sqrt(n[1:])
S = np.cumsum(a); b = S - c0
dl = np.zeros(NMAX + 1); dl[1:-1] = np.log((n[1:-1] + 1) / n[1:-1])
Q = np.cumsum(dl * b * b)
grid = np.unique(np.logspace(3, math.log10(NMAX), 4000).astype(int))
x = np.log(grid.astype(float)); y = Q[grid - 1]
res = {}
for lo in (3, 4, 5):
    m = grid >= 10 ** lo
    A = np.vstack([x[m], np.ones(m.sum())]).T
    slope, icpt = np.linalg.lstsq(A, y[m], rcond=None)[0]
    resid = y[m] - (slope * x[m] + icpt)
    res[f"slope_from_1e{lo}"] = float(slope)
    res[f"resid_rms_from_1e{lo}"] = float(np.sqrt((resid ** 2).mean()))
res["max_abs_b_from_1e4"] = float(np.abs(b[10000:]).max())
res["rms_b_over_log_scale_1e4_to_NMAX"] = float(math.sqrt(np.trapezoid(b[10000:NMAX] ** 2 / n[10000:NMAX], n[10000:NMAX]) / (math.log(NMAX) - math.log(10000))))
res["c0"] = c0; res["NMAX"] = NMAX
json.dump(res, open("slope_result.json", "w"), indent=1)
print(json.dumps(res, indent=1))
