"""Task 3(i): collect shards -> kappa; partial sums; tail estimate.
kappa = 2 sum_{gamma>0} 1/(gamma^2 |zeta'(rho)|^2)  (mean square of B(u) = sum 2 Re[e^{i gamma u}/(i gamma zeta'(rho))]).
Tail beyond the last zero T: tail ~ 2 m / T with m = mean of 1/|zeta'|^2 per unit height, from the last windows.
The tail estimate assumes 1/|zeta'|^2 has a stable mean density in height (Gonek-Hejhal type conjecture, UNPROVEN).
Usage: python kappa_collect.py SHARD_DIR OUT_DIR"""
import sys, csv, json, glob, math, os
from pathlib import Path
import numpy as np

sd, out = Path(sys.argv[1]), Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
rows = []
for p in sorted(glob.glob(str(sd / '**' / 'kappa_*.csv'), recursive=True)):
    with open(p) as f:
        r = csv.reader(f); next(r)
        rows += [(int(a[0]), a[1], float(a[1]), float(a[2]), float(a[3]), float(a[5])) for a in r]
rows.sort()
N = np.array([r[0] for r in rows]); assert (N == np.arange(1, len(N) + 1)).all(), "gap in zero indices"
g = np.array([r[2] for r in rows]); zr = np.array([r[3] for r in rows]); zi = np.array([r[4] for r in rows])
term = np.array([r[5] for r in rows]); assert (np.diff(g) > 0).all()
n = len(N); T = g[-1]
with open(out / 'kappa_all.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'gamma', 'zp_re', 'zp_im', 'term'])
    for r in rows: w.writerow([r[0], r[1], r[3], r[4], r[5]])
cum = np.cumsum(term)
inv2 = 1.0 / (zr ** 2 + zi ** 2)
# local mean density m of 1/|zeta'|^2 per unit height over windows of the last 40% of zeros
wins = []
lo = int(0.6 * n); step = max(1, (n - lo) // 4)
for a in range(lo, n - step + 1, step):
    b = a + step
    wins.append(inv2[a:b].sum() / (g[b - 1] - g[a - 1]))
m_bar = float(np.mean(wins)); m_spread = float(max(wins) - min(wins))
tail = 2 * m_bar / T
tail_err = 2 * m_spread / T
kappa_partial = float(cum[-1])
cuts = [c for c in (10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000) if c <= n]
tab = [dict(zeros=c, gamma_max=float(g[c - 1]), partial=float(cum[c - 1])) for c in cuts]
summ = dict(zero_source=os.environ.get('ZERO_SOURCE', 'platt'), zeros=int(n), gamma_max=float(T), kappa_partial=kappa_partial, m_bar=m_bar, m_spread=m_spread,
            tail_estimate=tail, tail_err=tail_err, kappa=kappa_partial + tail, kappa_err=tail_err,
            partial_sums=tab, windows_m=wins)
json.dump(summ, open(out / 'kappa_summary.json', 'w'), indent=1)
with open(out / 'kappa.md', 'w') as f:
    f.write("# kappa = 2 sum_{gamma>0} 1/(gamma^2 |zeta'(rho)|^2)  (NUMERICAL)\n\n")
    src = os.environ.get('ZERO_SOURCE', 'platt')
    f.write(("Zeros from Platt's tables read through the unmodified `platt_zeros.py`; " if src == 'platt' else
             "Zeros from `mpmath.zetazero` (fallback: Platt's tables were not available to the runner); ") +
            "zeta'(rho) by mpmath at 32 digits "
            "(residual |zeta(rho)| < 1e-20 sampled every 500 zeros). Assumes RH only through the labelling of rho as "
            "1/2 + i gamma; the tail estimate additionally assumes a stable mean density of 1/|zeta'|^2 (UNPROVEN).\n\n")
    f.write("| zeros | gamma_max | partial sum |\n|---|---|---|\n")
    for t in tab: f.write(f"| {t['zeros']:,} | {t['gamma_max']:.1f} | {t['partial']:.10f} |\n")
    f.write(f"\nTail estimate beyond gamma = {T:.1f}: {tail:.3e} +- {tail_err:.1e} (m = {m_bar:.4f} per unit height, spread {m_spread:.4f}).\n\n")
    f.write(f"**kappa = {kappa_partial + tail:.8f} +- {tail_err:.1e}** (partial {kappa_partial:.8f} + tail).\n")
print(json.dumps({k: v for k, v in summ.items() if k not in ('windows_m', 'partial_sums')}, indent=1))
