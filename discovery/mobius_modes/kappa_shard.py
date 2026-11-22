"""Task 3(i): |zeta'(rho)| and zeta'(rho) at 30 digits for zeros N in [start, start+count), zeros read through the
UNMODIFIED platt_zeros.py reader.  Usage: python kappa_shard.py PLATT_DIR START COUNT OUT.csv
(PLATT_DIR holds platt_zeros.py with zeta_folder pointed at the download)."""
import sys, csv, time
sys.path.insert(0, sys.argv[1])
import mpmath
from mpmath import mp, mpf, mpc
import platt_zeros as pz

start, count, out = int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
zs = list(pz.zeros_starting_at_N(start, count))      # reader sets its own (high) precision
assert zs and zs[0][0] == start, (zs[:1], start)
mp.dps = 32
t0 = time.time()
with open(out, 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['N', 'gamma', 'zp_re', 'zp_im', 'abs_zp', 'term'])
    for k, (N, g) in enumerate(zs):
        g = mpf(g)
        rho = mpc(mpf('0.5'), g)
        d = mpmath.zeta(rho, derivative=1)
        term = 2 / (g * g * (d.real ** 2 + d.imag ** 2))
        w.writerow([N, mpmath.nstr(g, 32), mpmath.nstr(d.real, 30), mpmath.nstr(d.imag, 30),
                    mpmath.nstr(abs(d), 25), mpmath.nstr(term, 25)])
        if k == 20 and k < len(zs) - 1:
            print('rate', (time.time() - t0) / 21, 's/zero', flush=True)
        if k % 500 == 0:
            # residual check on a sample: |zeta(rho)| must be ~ 0 at a Platt zero
            r = abs(mpmath.zeta(rho))
            assert r < mpf(10) ** -20, (N, r)
print('done', len(zs), 'zeros in', round(time.time() - t0, 1), 's')
