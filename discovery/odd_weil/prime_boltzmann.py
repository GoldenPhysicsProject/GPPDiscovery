"""Prime gas as a Boltzmann gas (Claude chat seat, 2026-10-05).

Part 1 (exact, discrete): fusion gas on {2..N}, collisions n + m <-> nm, detailed balance with pi(n) ∝ n^-beta
(energy log n is additive under fusion). Linearized collision (Dirichlet) form
    Q(f) = sum_{n,m>=2, nm<=N} pi(n) pi(m) (f(n) + f(m) - f(nm))^2  >= 0   (H-theorem).
Collision invariants = completely additive f = one free parameter per prime: prediction nullity = pi(N).
Their generalized Gibbs states pi(n) ∝ prod_p p^{-beta_p v_p(n)} are exactly the Euler-product (Euler-coherent) states.

Part 2 (continuum log-line, semilocal odd Weil form): collisions = jumps of the log-scale by log n with rate
c_n = b_n n^{-1/2}. Exact identity -2 c_n G_n = |c_n| (2 G0 - 2 sgn(c_n) G_n) - 2|c_n| G0 gives
    A = (W - 2 S G0)  +  Lap,     Lap = sum |c_n| (2 G0 - 2 sgn c_n G_n)  ⪰ 0   (collision / H-theorem form),
    'Mass' = W - 2 S G0 (archimedean continuum minus the scalar collision mass).
Collision strength needed:  c* = max_v  -<v,Mass v> / <v,Lap v>.  A ⪰ 0  <=>  c* <= 1.
c* < 1: the physical collision term over-pays (slack). c* > 1: collisions too weak, positivity fails.
Compared across zeta (Euler product, molecular chaos), chi5real (Euler product), DH (no Euler product).
"""
import json, sys, numpy as np, scipy.linalg as sl
from dh_control import build
out = {}
# Part 1
def part1(N, beta=1.5):
    n = np.arange(2, N + 1); idx = {v: i for i, v in enumerate(n)}
    pi = n.astype(float) ** (-beta)
    rows = []
    for a in n:
        for b in n:
            if a <= b and a * b <= N:
                r = np.zeros(len(n)); r[idx[a]] += 1; r[idx[b]] += 1; r[idx[a * b]] -= 1
                rows.append(np.sqrt(pi[idx[a]] * pi[idx[b]]) * r)
    R = np.array(rows); Q = R.T @ R
    ev = np.linalg.eigvalsh(Q)
    nullity = int(np.sum(ev < 1e-12 * ev[-1]))
    nprimes = sum(1 for v in n if all(v % q for q in range(2, int(v**0.5) + 1)))
    return dict(N=N, nullity=nullity, primes_le_N=nprimes, min_eig=float(ev[0]), first_nonzero=float(ev[nullity]))
out["part1"] = [part1(N) for N in (30, 60, 120)]
print(out["part1"], flush=True)
# Part 2
p2 = {}
for kind in ["zeta", "chi5real", "dh"]:
    for lam in [3.0, 4.0, 5.0, 6.0]:
        d = build(kind, lam, 20)
        A, Lap, Mass = [(X + X.T) / 2 for X in (d["A"], d["Lap"], d["Mass"])]
        lap_ev = np.linalg.eigvalsh(Lap)
        cs = sl.eigh(-Mass, Lap + 1e-14 * np.eye(len(Lap)), eigvals_only=True)
        p2[f"{kind}_{lam}"] = dict(c_star=float(cs[-1]), minA=float(np.linalg.eigvalsh(A)[0]),
                                   minMass=float(np.linalg.eigvalsh(Mass)[0]), minLap=float(lap_ev[0]), S=float(d["S"]))
        print(kind, lam, p2[f"{kind}_{lam}"], flush=True)
out["part2"] = p2
json.dump(out, open(sys.argv[1], "w"), indent=1)
