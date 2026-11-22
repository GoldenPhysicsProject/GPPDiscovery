"""Task 3(iv), pass 1: over ALL zeros in the tables, the 1000 smallest normalized gaps
g_n = (gamma_{n+1} - gamma_n) log(gamma_n / 2 pi) / (2 pi).
Gaps are built from Platt's INTEGER deltas (units 2**-101), so they are exact to float64 relative precision
independent of the height (a float64 gamma at t ~ 1e7 would not resolve them).  Block-crossing gaps use exact
Python ints.  Usage: python gaps.py ZETA_FOLDER OUT_CSV"""
import sys, csv, math, heapq
from pathlib import Path
import numpy as np
import plattio

zeta = Path(sys.argv[1]); out = sys.argv[2]
files = sorted((zeta / 'data').glob('zeros_*.dat'), key=lambda q: int(q.stem.split('_')[1]))
KEEP = 1000
cand = []          # (g_norm, N_low, gap)
prev_last_int = None; prev_last_f = None; total = 0
for p in files:
    for b, off, t0, t1, N0, N1 in plattio.blocks(p):
        n = N1 - N0
        if n == 0: continue
        arr = plattio.read_block(p, off, n)
        d = plattio.delta_float(arr)
        gam = t0 + np.cumsum(d)                       # float64 heights, only used inside log()
        gap = d.copy()
        t0i = plattio.t0_int(t0)
        first_int = t0i + plattio.entry_int(arr, 0)
        if prev_last_int is not None:
            gap[0] = float(first_int - prev_last_int) * 2.0 ** -101
            low = np.concatenate(([prev_last_f], gam[:-1]))
            Nlow = N0 + np.arange(n)                  # 1-based index of the lower zero of gap i is N0 + i
        else:
            gap = gap[1:]; low = gam[:-1]; Nlow = N0 + 1 + np.arange(n - 1)
        g = gap * np.log(low / (2 * np.pi)) / (2 * np.pi)
        k = min(KEEP, len(g))
        sel = np.argpartition(g, k - 1)[:k]
        cand += [(float(g[i]), int(Nlow[i]), float(gap[i])) for i in sel]
        cand = heapq.nsmallest(KEEP, cand)
        prev_last_int = t0i + plattio.delta_sum_int(arr); prev_last_f = float(gam[-1])
        total += len(g)
    print(p.name, 'gaps so far', total, 'current 1000th smallest', cand[-1][0] if len(cand) == KEEP else None, flush=True)
cand.sort()
with open(out, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['rank', 'N_low', 'gap', 'normalized_gap'])
    for r, (g, N, gap) in enumerate(cand, 1): w.writerow([r, N, repr(gap), repr(g)])
print('gaps scanned:', total, 'smallest normalized gap:', cand[0])
