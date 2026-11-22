"""Task 3(iv), pass 2: for the smallest-gap pairs, exact heights through the UNMODIFIED reader, then |zeta'| at both
ends and the pair's contribution/share of kappa.  Usage: python gaps_zeta.py PLATT_DIR GAPS_CSV KAPPA_SUMMARY_JSON OUT_CSV [NPAIRS=100]"""
import sys, csv, json
sys.path.insert(0, sys.argv[1])
import mpmath
from mpmath import mp, mpf, mpc
import platt_zeros as pz

rows = list(csv.DictReader(open(sys.argv[2])))
kap = json.load(open(sys.argv[3]))['kappa']
npairs = int(sys.argv[5]) if len(sys.argv) > 5 else 100
res = []
for r in rows[:npairs]:
    N = int(r['N_low'])
    z = list(pz.zeros_starting_at_N(N, 2)); assert z[0][0] == N and z[1][0] == N + 1
    mp.dps = 32
    g0, g1 = mpf(z[0][1]), mpf(z[1][1])
    d0 = abs(mpmath.zeta(mpc(0.5, g0), derivative=1)); d1 = abs(mpmath.zeta(mpc(0.5, g1), derivative=1))
    t0 = 2 / (g0 ** 2 * d0 ** 2); t1 = 2 / (g1 ** 2 * d1 ** 2)
    res.append([int(r['rank']), N, mpmath.nstr(g0, 30), mpmath.nstr(g1, 30), mpmath.nstr(g1 - g0, 12), r['normalized_gap'],
                mpmath.nstr(d0, 12), mpmath.nstr(d1, 12), mpmath.nstr(t0 + t1, 6), mpmath.nstr((t0 + t1) / kap, 6)])
    print(res[-1][:2], res[-1][5], res[-1][6:8], flush=True)
with open(sys.argv[4], 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['rank', 'N_low', 'gamma_low', 'gamma_high', 'gap', 'normalized_gap', 'abs_zp_low', 'abs_zp_high',
                'contribution_to_kappa', 'share_of_kappa'])
    w.writerows(res)
tot = sum(float(r[8]) for r in res)
print('total contribution of these', len(res), 'pairs to kappa:', tot, 'share', tot / kap)
