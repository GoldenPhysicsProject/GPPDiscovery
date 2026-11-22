"""Parseval check for the Moebius energy slope (NUMERICAL).

If B(u) = sum_{n <= e^u} mu(n)/sqrt(n) - 1/zeta(1/2) is, under RH with simple zeros, the almost-periodic
function  B(u) ~ sum_{gamma>0} 2 Re[ e^{i gamma u} / (rho zeta'(rho)) ],  rho = 1/2 + i gamma,
then its mean square, hence the slope dQ/d(log N), is  m2 = 2 * sum_{gamma>0} 1/|rho zeta'(rho)|^2.
This script sums that over the first NZ zeros (mpmath) and compares it with the slope measured from the
sieve in scan.py. Agreement is a consistency check between the sieve and the explicit formula, not a
proof of anything; the truncation tail is estimated, not bounded.

Usage: python3 -I parseval.py [NZ]
"""
import sys, json, math
from mpmath import mp, zetazero, zeta, mpf, mpc

mp.dps = 20
NZ = int(sys.argv[1]) if len(sys.argv) > 1 else 300
partial = []
s = mpf(0)
for n in range(1, NZ + 1):
    rho = zetazero(n)
    d = zeta(rho, derivative=1)
    s += 2 / abs(rho * d) ** 2
    if n in (10, 25, 50, 100, 200, 300, 500, 1000) or n == NZ:
        partial.append((n, float(rho.imag), float(s)))
        print(f"N_zeros={n:5d} gamma_max={float(rho.imag):9.3f}  m2={float(s):.6f}", flush=True)
json.dump({"NZ": NZ, "partial_sums": partial}, open("parseval_result.json", "w"), indent=1)
