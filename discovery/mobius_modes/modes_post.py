"""Task 2(ii) post-processing + independent cross-check of the incremental F(y) by the direct mode form
sum_m J(m) T(y,m)^2 (loop over squarefree m, T by strided sums) and, for small y, the double sum definition.
Usage: python modes_post.py modes_F.csv OUT_DIR [Ycheck ...]"""
import sys, csv, math, json
from pathlib import Path
import numpy as np
inp, out = sys.argv[1], Path(sys.argv[2]); out.mkdir(parents=True, exist_ok=True)
checks = [int(x) for x in sys.argv[3:]] or [100, 1000, 10000, 100000]
rows = list(csv.DictReader(open(inp)))
byy = {int(r['y']): r for r in rows}
Y = max(checks)
mu = np.zeros(Y + 1, dtype=np.int64); mu[1] = 1
for i in range(1, Y + 1):
    if mu[i]: mu[2 * i::i] -= mu[i]
J = np.zeros(Y + 1)
for k in range(1, Y + 1):
    m = np.arange(k, Y + 1, k); J[m] += mu[m // k] * math.sqrt(k)
a = np.zeros(Y + 1); n = np.arange(1, Y + 1); a[1:] = mu[1:] / np.sqrt(n)
chk = []
for y in checks:
    F = 0.0
    for m in range(1, y + 1):
        if mu[m] == 0: continue
        T = a[m:y + 1:m].sum(); F += J[m] * T * T
    Fc = float(byy[y]['F']); chk.append(dict(y=y, F_mode_direct=F, F_incremental_C=Fc, rel_diff=abs(F - Fc) / Fc))
    if y <= 400:
        G = 0.0
        for d in range(1, y + 1):
            for e in range(1, y + 1):
                if mu[d] and mu[e]: G += mu[d] * mu[e] * math.sqrt(math.gcd(d, e)) / math.sqrt(d * e)
        chk[-1]['F_double_sum_definition'] = G
yy = sorted(byy)
sel = [y for y in yy if y in (10, 100, 400, 900, 1600) or y in (1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000, 10000000)]
with open(out / 'modes_F_table.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['y', 'F', 'F_over_sqrt_y', 'S', 'S_sq', 'F_over_S_sq'])
    for y in yy:
        r = byy[y]; w.writerow([y, r['F'], r['F_over_sqrt_y'], r['S'], r['S_sq'], r['F_over_S_sq']])
R = np.array([float(byy[y]['F_over_sqrt_y']) for y in yy]); Ys = np.array(yy)
big = Ys >= 1600
with open(out / 'modes.md', 'w') as f:
    f.write("# F(y) through the mode form (NUMERICAL evidence for an UNCONDITIONAL identity)\n\n")
    f.write("F(y) = sum_(d,e<=y) mu(d)mu(e) sqrt(gcd(d,e))/sqrt(de) = sum_m J(m) T(y,m)^2 (Lean: `GppMobiusModes.mode_decomposition`). "
            "Computed for every y <= Y by the incremental update T(y+1,m) = T(y,m) + a(y+1) for m | y+1 (O(Y log Y) total).\n\n")
    f.write("| y | F(y) | F/sqrt(y) | S(y) | S(y)^2 | F/S^2 |\n|---|---|---|---|---|---|\n")
    for y in sel:
        r = byy[y]; f.write(f"| {y:,} | {float(r['F']):.5f} | {float(r['F_over_sqrt_y']):.5f} | {float(r['S']):.5f} | {float(r['S_sq']):.5f} | {float(r['F_over_S_sq']):.1f} |\n")
    f.write(f"\nF/sqrt(y) over y in [1600, {yy[-1]:,}]: min {R[big].min():.5f}, max {R[big].max():.5f}, last {R[-1]:.5f}. "
            f"min J(m) over m <= {yy[-1]:,}: {float(rows[-1]['minJ_all']):.6f} (J >= 0 is the Lean theorem `J_nonneg`).\n\n")
    f.write("## Cross-checks\n\n| y | direct mode form | incremental (C) | rel diff | double-sum definition |\n|---|---|---|---|---|\n")
    for c in chk:
        f.write(f"| {c['y']:,} | {c['F_mode_direct']:.10f} | {c['F_incremental_C']:.10f} | {c['rel_diff']:.1e} | {c.get('F_double_sum_definition','')} |\n")
json.dump(dict(checks=chk, ratio_min=float(R[big].min()), ratio_max=float(R[big].max()), ratio_last=float(R[-1]), y_last=int(yy[-1])),
          open(out / 'modes_summary.json', 'w'), indent=1)
print(json.dumps(chk, indent=1)); print('F/sqrt(y) for y>=1600: min %.5f max %.5f last %.5f' % (R[big].min(), R[big].max(), R[-1]))
