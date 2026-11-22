"""Task 3(iii): B_T(u) = sum_{0<gamma<=T} 2 Re[ e^{i gamma u} / (i gamma zeta'(rho)) ] against the sieved
B(u) = S(floor(e^u)) - c on u in [10, 18], with 10^3, 10^4, 10^5 (and intermediate) zeros.
Needs the sieved S on the u-grid: run mu_S_grid first (grid file written by this script with --grid).
Usage: python reconstruct.py --grid GRID.txt            (writes the K list)
       python reconstruct.py KAPPA_ALL.csv S_GRID.csv OUT_DIR"""
import sys, csv, json, math
from pathlib import Path
import numpy as np

U0, U1, G = 10.0, 18.0, 20001
u = np.linspace(U0, U1, G)
K = np.floor(np.exp(u)).astype(np.int64)
if sys.argv[1] == '--grid':
    assert (np.diff(K) > 0).all()
    Path(sys.argv[2]).write_text('\n'.join(str(int(k)) for k in K) + '\n'); sys.exit(0)
from mpmath import mp, zeta, mpf
mp.dps = 30
c = float(1 / zeta(mpf(1) / 2))
kap, sg, out = sys.argv[1], sys.argv[2], Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
S = {int(r['K']): float(r['S']) for r in csv.DictReader(open(sg))}
B = np.array([S[int(k)] for k in K]) - c
rows = list(csv.DictReader(open(kap)))
g = np.array([float(r['gamma']) for r in rows]); zr = np.array([float(r['zp_re']) for r in rows]); zi = np.array([float(r['zp_im']) for r in rows])
coef = 2.0 / (1j * g * (zr + 1j * zi))
n = len(g)
cuts = [x for x in (10, 30, 100, 300, 1000, 3000, 10000, 30000, 100000, 300000, 1000000) if x <= n]
BT = np.zeros(G); res = []
rmsB = float(np.sqrt(np.mean(B ** 2)))
ci = 0; chunk = 256; idx = 0
saved = {}
while idx < n:
    nxt = min(n, idx + chunk)
    if ci < len(cuts): nxt = min(nxt, cuts[ci])
    ph = np.outer(g[idx:nxt], u)
    BT += (coef[idx:nxt, None] * (np.cos(ph) + 1j * np.sin(ph))).real.sum(axis=0)
    idx = nxt
    if ci < len(cuts) and idx == cuts[ci]:
        err = BT - B
        res.append(dict(zeros=idx, T=float(g[idx - 1]), rms_error=float(np.sqrt(np.mean(err ** 2))), rms_B=rmsB,
                        relative=float(np.sqrt(np.mean(err ** 2)) / rmsB),
                        corr=float(np.corrcoef(BT, B)[0, 1])))
        saved[idx] = BT.copy(); ci += 1
json.dump(res, open(out / 'reconstruct_summary.json', 'w'), indent=1)
# overlay data on a coarse subsample (every 10th grid point) for the three headline cuts
head = [x for x in (1000, 10000, 100000) if x in saved]
with open(out / 'reconstruct_overlay.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['u', 'K', 'B_sieve'] + [f'B_T_{h}zeros' for h in head])
    for i in range(0, G, 10):
        w.writerow([f"{u[i]:.5f}", int(K[i]), f"{B[i]:.6f}"] + [f"{saved[h][i]:.6f}" for h in head])
with open(out / 'reconstruct.md', 'w') as f:
    f.write("# Explicit-formula reconstruction of B(u) (NUMERICAL; assumes RH and simple zeros for the form of the series)\n\n")
    f.write(f"B(u) = S(floor(e^u)) - c, c = 1/zeta(1/2) = {c:.10f}. B_T(u) = sum_(0<gamma<=T) 2 Re[e^(i gamma u)/(i gamma zeta'(rho))] "
            f"(the pole of x^s/(s zeta(s+1/2)) at s = i gamma). u-grid: {G} equispaced points on [{U0:g}, {U1:g}] (N from {K[0]:,} to {K[-1]:,}). "
            f"rms of the sieved B on the grid: {rmsB:.5f}. The sieved B is a step function with jumps mu(k)/sqrt(k) at integers k, "
            "so a floor on the mismatch is expected from the frequencies beyond T.\n\n")
    f.write("| zeros | T = gamma_max | L2 mismatch (rms) | relative to rms B | correlation |\n|---|---|---|---|---|\n")
    for r in res:
        f.write(f"| {r['zeros']:,} | {r['T']:.1f} | {r['rms_error']:.5f} | {r['relative']:.4f} | {r['corr']:.5f} |\n")
print(json.dumps(res, indent=1))
