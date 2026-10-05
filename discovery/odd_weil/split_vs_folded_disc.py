"""Split vs folded across discriminants (Claude chat seat, 2026-10-05).
Real quadratic (two real places, Gamma_R^2, 'split') d = 5, 8, 12, 13, 17
Imaginary quadratic (one complex place, Gamma_C = Gamma_R(s)Gamma_R(s+1), 'folded') d = -3, -4, -7, -8, -11
zeta_K = zeta * L(chi_d); odd semilocal Weil form A_K = A_zeta + A_chi_d. Report min eigenvalue at lambda = 3, 3.5, 4, 4.5
(above the float64 floor ~1e-12 seen in the first run). Question: is the split margin systematically above the folded one
at matched conductor size?"""
import json, sys, numpy as np
from dh_control import build
N = 20
out = {}
for lam in [3.0, 3.5, 4.0, 4.5]:
    Az = build("zeta", lam, N)["A"]
    row = {}
    for d in [5, 8, 12, 13, 17, -3, -4, -7, -8, -11]:
        A = Az + build(f"kron:{d}", lam, N)["A"]; A = (A + A.T) / 2
        row[str(d)] = float(np.linalg.eigvalsh(A)[0])
    out[str(lam)] = row
    print(lam, row, flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
