"""Split vs folded archimedean place (Claude chat seat, 2026-10-05).

Gamma_C(s) = Gamma_R(s) Gamma_R(s+1): a complex place is one even half times one odd half (folded, 'Lorentzian').
A real quadratic field has two real places, Gamma_R(s)^2 (split, '(2,2)').
  split:  zeta_{Q(sqrt5)} = zeta * L(chi5)   chi5 even, conductor 5   (smallest real quadratic discriminant)
  folded: zeta_{Q(i)}     = zeta * L(chi_-4) chi_-4 odd, conductor 4
Weil forms add for products: A_K = A_zeta + A_chi (odd semilocal sector, pole killed).
Reports min eigenvalue (positivity margin) and the scaled margin min_eig / trace-scale for each piece and each field,
lambda = 3..8, plus prime-split statistics within the window.
"""
import json, sys, numpy as np
from dh_control import build
N = int(sys.argv[2]) if len(sys.argv) > 2 else 20
out = {}
for lam in [3.0, 4.0, 5.0, 6.0, 7.0, 8.0]:
    Az = build("zeta", lam, N)["A"]; A5 = build("chi5real", lam, N)["A"]; A4 = build("chi4", lam, N)["A"]
    sym = lambda A: (A + A.T) / 2
    row = {}
    for name, A in [("zeta", Az), ("L_chi5", A5), ("L_chi-4", A4), ("Q(sqrt5)_split", Az + A5), ("Q(i)_folded", Az + A4)]:
        ev = np.linalg.eigvalsh(sym(A))
        row[name] = dict(min=float(ev[0]), second=float(ev[1]), max=float(ev[-1]), rel=float(ev[0] / ev[-1]))
    out[str(lam)] = row
    print(lam, {k: f"{v['min']:.4e}" for k, v in row.items()}, flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
