"""Task 2(i) post-processing.  Q_c(log N) = WSS - 2c WS + c^2 W (energy identity), from mu_scan.csv.
Usage: python direct_post.py mu_scan.csv OUT_DIR"""
import sys, csv, math, json
from pathlib import Path
import numpy as np
from mpmath import mp, zeta, mpf
mp.dps = 30
c = float(1 / zeta(mpf(1) / 2))
rows = list(csv.reader(open(sys.argv[1]))); hdr = rows[0]; rows = rows[1:]
out = Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
NF = 21
Ns, L, W, WS, WSS, maxb = [], [], [], [], [], []
for r in rows:
    v = [float(x) for x in r]
    Ns.append(int(r[0])); L.append(math.log(int(r[0]))); W.append(v[1])
    WS.append([v[2 + 2 * s] for s in range(NF)]); WSS.append([v[3 + 2 * s] for s in range(NF)]); maxb.append(v[-1])
L = np.array(L); W = np.array(W); WS = np.array(WS); WSS = np.array(WSS)
Qmu = WSS[:, 0] - 2 * c * WS[:, 0] + c * c * W
copt_mu = WS[:, 0] / W
Qmu_opt = WSS[:, 0] - WS[:, 0] ** 2 / W
with open(out / 'direct_Q_mu.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'L', 'Q', 'Q_over_L', 'Q_over_L2', 'max_abs_b_k_from_1e4', 'c_opt_logmean_of_S', 'Q_at_c_opt_over_L'])
    for i, N in enumerate(Ns):
        w.writerow([N, f"{L[i]:.6f}", f"{Qmu[i]:.8f}", f"{Qmu[i]/L[i]:.8f}", f"{Qmu[i]/L[i]**2:.8f}", f"{maxb[i]:.6f}",
                    f"{copt_mu[i]:.8f}", f"{Qmu_opt[i]/L[i]:.8f}"])
# random comparison: s = 1..20
Qr0 = WSS[:, 1:] - 0 + 0                      # c = 0
Qr0 = WSS[:, 1:]
Qropt = WSS[:, 1:] - WS[:, 1:] ** 2 / W[:, None]   # best constant c = log-mean of S (most generous to the random functions)
with open(out / 'direct_Q_random.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['N', 'L', 'mean_Q_over_L_c0', 'min_Q_over_L_c0', 'max_Q_over_L_c0', 'mean_Q_over_L2_c0',
                'mean_Q_over_L_copt', 'min_Q_over_L_copt', 'max_Q_over_L_copt', 'mean_Q_over_L2_copt'])
    for i, N in enumerate(Ns):
        a, b = Qr0[i] / L[i], Qropt[i] / L[i]
        w.writerow([N, f"{L[i]:.6f}", f"{a.mean():.6f}", f"{a.min():.6f}", f"{a.max():.6f}", f"{(Qr0[i]/L[i]**2).mean():.6f}",
                    f"{b.mean():.6f}", f"{b.min():.6f}", f"{b.max():.6f}", f"{(Qropt[i]/L[i]**2).mean():.6f}"])
# local slope of Q(L) between checkpoints, and a least-squares slope from N >= 1e4
sel = np.array(Ns) >= 10 ** 4
A = np.vstack([L[sel], np.ones(sel.sum())]).T
slope, icpt = np.linalg.lstsq(A, Qmu[sel], rcond=None)[0]
res = Qmu[sel] - (slope * L[sel] + icpt)
loc = [(Ns[i], Ns[i + 4], (Qmu[i + 4] - Qmu[i]) / (L[i + 4] - L[i])) for i in range(len(Ns) - 4)]
slopes_r = []
for s in range(20):
    q = Qropt[sel, s] if False else Qropt[:, s][sel]
    sl, _ = np.linalg.lstsq(A, q, rcond=None)[0]
    slopes_r.append(float(sl))
summ = dict(c=c, N_max=Ns[-1], Q_over_L_at_Nmax=float(Qmu[-1] / L[-1]), Q_over_L2_at_Nmax=float(Qmu[-1] / L[-1] ** 2),
            ls_slope_from_1e4=float(slope), ls_intercept=float(icpt), ls_resid_rms=float(np.sqrt((res ** 2).mean())),
            local_slopes_decade=loc, max_abs_b_overall=float(maxb[-1]),
            random_copt_ls_slope_mean=float(np.mean(slopes_r)), random_copt_ls_slope_min=float(min(slopes_r)),
            random_copt_ls_slope_max=float(max(slopes_r)))
json.dump(summ, open(out / 'direct_summary.json', 'w'), indent=1)
with open(out / 'direct.md', 'w') as f:
    f.write("# Direct Moebius side: Q_c(log N) (NUMERICAL)\n\n")
    f.write(f"c = 1/zeta(1/2) = {c:.12f}. b_k = S(k) - c, Q(log N) = sum_(k<N) log((k+1)/k) b_k^2 (the energy identity; "
            "Lean: `GppMobiusEcho.energy_identity`). Segmented sieve in C to N = 10^8; random comparison = 20 functions "
            "f(n) = prod_(p|n) eps_p on squarefree n (eps_p = +-1 hashed), 0 otherwise, run at c = 0 and at their own "
            "best constant c_opt (log-mean of S, the choice most favourable to them).\n\n")
    f.write("| N | L = log N | Q | Q/L | Q/L^2 | max abs b_k (k >= 1e4) |\n|---|---|---|---|---|---|\n")
    for i, N in enumerate(Ns):
        if i % 2 == 0 or i == len(Ns) - 1:
            f.write(f"| {N:,} | {L[i]:.3f} | {Qmu[i]:.5f} | {Qmu[i]/L[i]:.5f} | {Qmu[i]/L[i]**2:.5f} | {maxb[i]:.4f} |\n")
    f.write(f"\nLeast-squares fit Q = slope*L + intercept over N >= 1e4: slope = **{slope:.6f}**, intercept {icpt:.4f}, "
            f"residual rms {summ['ls_resid_rms']:.4f}. Q/L converges to the slope from below/above like intercept/L, so the "
            "slope, not Q/L at a finite N, is the quantity compared with kappa.\n\n")
    f.write("## Random comparison (20 functions)\n\n| N | mean Q/L (c=0) | mean Q/L^2 (c=0) | mean Q/L (c_opt) | min..max Q/L (c_opt) | mean Q/L^2 (c_opt) |\n|---|---|---|---|---|---|\n")
    rr = list(csv.reader(open(out / 'direct_Q_random.csv')))[1:]
    for i, r in enumerate(rr):
        if i % 2 == 0 or i == len(rr) - 1:
            f.write(f"| {int(r[0]):,} | {float(r[2]):.4f} | {float(r[5]):.5f} | {float(r[6]):.4f} | {float(r[7]):.4f}..{float(r[8]):.4f} | {float(r[9]):.5f} |\n")
    f.write(f"\nLeast-squares slopes of Q(L) for the random functions at c_opt (N >= 1e4): mean {summ['random_copt_ls_slope_mean']:.4f}, "
            f"range {summ['random_copt_ls_slope_min']:.4f}..{summ['random_copt_ls_slope_max']:.4f}.\n")
print(json.dumps({k: v for k, v in summ.items() if k != 'local_slopes_decade'}, indent=1))
