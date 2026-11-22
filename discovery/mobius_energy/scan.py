"""Numerical receipt for the Moebius energy Q_c(log N) = sum_{k<N} log((k+1)/k) * b_k^2,
b_k = sum_{n<=k} mu(n)/sqrt(n) - c   (GppMobiusEcho in GPPVerify; the energy of the Moebius echo packet).

NUMERICAL, not a proof. It records how Q_c(log N) and max_{k<=N}|b_k| grow for c = 1/zeta(1/2)
(the packet's natural constant), for c = 0, and for a wrong constant, out to N = NMAX. The point is to
see whether the growth looks subpower (consistent with RH; hypothesis (ii) of rh_of_storage_work_causal
with the trivial storage stor = 0) and how much c matters.

Usage: python3 -I scan.py [NMAX]      (default 3e7)
"""
import sys, json, math
import numpy as np
from mpmath import mp, zeta

NMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 30_000_000
mp.dps = 30
c0 = float(1 / zeta(mp.mpf(1) / 2))          # 1/zeta(1/2)

# mu(n) for n <= NMAX by a prime sieve over int8
mu = np.ones(NMAX + 1, dtype=np.int8)
mu[0] = 0
is_p = np.ones(NMAX + 1, dtype=bool); is_p[:2] = False
for i in range(2, int(NMAX ** 0.5) + 1):
    if is_p[i]:
        is_p[i * i::i] = False
for p in np.flatnonzero(is_p):
    mu[p::p] *= -1
    if p * p <= NMAX:
        mu[p * p::p * p] = 0

n = np.arange(NMAX + 1, dtype=np.float64)
a = np.zeros(NMAX + 1)
a[1:] = mu[1:] / np.sqrt(n[1:])
S = np.cumsum(a)                              # S[k] = sum_{n<=k} a(n)
dl = np.zeros(NMAX + 1)
dl[1:-1] = np.log((n[1:-1] + 1) / n[1:-1])    # log((k+1)/k), k = 1..NMAX-1

rows = []
checkpoints = sorted(set([10 ** e for e in range(2, 9) if 10 ** e <= NMAX] + [NMAX]))
for name, c in [("c=1/zeta(1/2)", c0), ("c=0", 0.0), ("c=-0.5", -0.5)]:
    b = S - c
    Q = np.cumsum(dl * b * b)                 # Q[N-1] = sum_{k<N}; Q at log N = Q[N-1]
    run_max = np.maximum.accumulate(np.abs(b[1:]))
    for N in checkpoints:
        rows.append({"c": name, "cval": c, "N": N, "S(N)": float(S[N]), "b_N": float(b[N]),
                     "max|b_k|<=N": float(run_max[N - 1]), "Q(log N)": float(Q[N - 1]),
                     "log N": math.log(N)})
out = {"c0": c0, "NMAX": NMAX, "rows": rows}
# local growth exponent of Q between consecutive checkpoints: d log Q / d log N
for name in ("c=1/zeta(1/2)", "c=0", "c=-0.5"):
    rs = [r for r in rows if r["c"] == name]
    for r1, r2 in zip(rs, rs[1:]):
        r2["exponent_vs_prev"] = (math.log(r2["Q(log N)"]) - math.log(r1["Q(log N)"])) / (
            math.log(r2["N"]) - math.log(r1["N"]))
json.dump(out, open("scan_result.json", "w"), indent=1)
print("c0 = 1/zeta(1/2) =", c0)
for r in rows:
    print(f'{r["c"]:>14} N={r["N"]:>9} S={r["S(N)"]:+.6f} b_N={r["b_N"]:+.6f} max|b|={r["max|b_k|<=N"]:.4f} '
          f'Q={r["Q(log N)"]:.5g} exp={r.get("exponent_vs_prev", float("nan")):.3f}')
