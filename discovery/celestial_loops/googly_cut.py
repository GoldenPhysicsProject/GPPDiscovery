"""Is the cut's flat (stabilizer) direction Penrose's googly projection G? (Claude chat seat, 2026-10-05)

Split signature (2,2), real independent spinors: l = lam (x) lamt, lam = sqrt(w)(1,z), lamt = sqrt(w)(1,zb).
Cut: z6 = -1/zb5, zb6 = -1/z5, w5 = M/(2(1+u)), w6 = M u/(2(1+u)), u = z5 zb5.
Checks (sympy, exact):
 1. P^{AA'} = sum_i lam_i^A lamt_i^{A'} is (M/2) * I-type (CM frame) and invertible: the pair is a massive (two-twistor) object.
 2. P maps leg 5's primed spinor to leg 6's unprimed spinor: P.lamt5 ∝ lam6 (and P.lamt6 ∝ lam5).
 3. Penrose TN43 flat-space projectors at the point x^{AA'} = c P^{AA'} (incidence omega = i x pi):
    L(Z) = (i x pi, pi), G(Z) = (omega - i x pi, 0); verify L+G=1, L^2=L, G^2=G, LG=GL=0. The twistor through x with
    pi = lamt5 has omega ∝ lam6: L reads leg 5's primed half, its S_A part is leg 6's unprimed half.
 4. The flat direction z5->e^t z5, zb5->e^-t zb5 (with the induced z6, zb6) leaves P invariant exactly: it lies in the
    little group of P (SO(1,2) in split signature, SO(3) in Lorentzian). Report the 4x4 Lorentz generator it induces.
"""
import json, sys, os
import sympy as sp
out = {}
z5, zb5, M, t, c = sp.symbols('z5 zb5 M t c', positive=True)
u = z5*zb5
w5, w6 = M/(2*(1+u)), M*u/(2*(1+u))
z6, zb6 = -1/zb5, -1/z5
lam = lambda w, z: sp.Matrix([1, z]) * sp.sqrt(w)
l5, lt5, l6, lt6 = lam(w5, z5), lam(w5, zb5), lam(w6, z6), lam(w6, zb6)
Pm = sp.simplify(l5*lt5.T + l6*lt6.T)
out["P^{AA'} matrix"] = str(Pm)
out["det P"] = str(sp.simplify(Pm.det()))
eps = sp.Matrix([[0, 1], [-1, 0]])
v = sp.simplify(Pm*eps*lt5)  # contract primed index with epsilon
out["P.lamt5 vs lam6 (ratio components)"] = str(sp.simplify(sp.Matrix([v[0]/l6[0], v[1]/l6[1]])))
v2 = sp.simplify(Pm*eps*lt6)
out["P.lamt6 vs lam5 (ratio components)"] = str(sp.simplify(sp.Matrix([v2[0]/l5[0], v2[1]/l5[1]])))
# 3. projectors on 4-dim twistor space Z = (omega^A, pi_A'), x = c P
x = c*Pm*eps
Lp = sp.Matrix(sp.BlockMatrix([[sp.zeros(2), sp.I*x], [sp.zeros(2), sp.eye(2)]]))
Gp = sp.Matrix(sp.BlockMatrix([[sp.eye(2), -sp.I*x], [sp.zeros(2), sp.zeros(2)]]))
out["L+G=1"] = str(sp.simplify(Lp+Gp-sp.eye(4)) == sp.zeros(4))
out["L^2=L, G^2=G, LG=0, GL=0"] = str([sp.simplify(Lp*Lp-Lp) == sp.zeros(4), sp.simplify(Gp*Gp-Gp) == sp.zeros(4),
                                        sp.simplify(Lp*Gp) == sp.zeros(4), sp.simplify(Gp*Lp) == sp.zeros(4)])
Z = sp.Matrix([0, 0, lt5[0], lt5[1]])
LZ = sp.simplify(Lp*Z)
out["omega-part of L(Z) for pi=lamt5, ratio to lam6"] = str(sp.simplify(sp.Matrix([LZ[0]/l6[0], LZ[1]/l6[1]])))
# 4. flat direction preserves P
sub = {z5: sp.exp(t)*z5, zb5: sp.exp(-t)*zb5}
Pt = sp.simplify(Pm.subs(sub, simultaneous=True))
out["P invariant under flat direction"] = str(sp.simplify(Pt - Pm) == sp.zeros(2))
gen5 = sp.simplify(sp.diff(l5.subs(sub, simultaneous=True), t).subs(t, 0))
out["d/dt lam5 at t=0 (generator on unprimed spinor)"] = str(gen5)
p = sys.argv[1] if len(sys.argv) > 1 else "googly_cut.json"
os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
json.dump(out, open(p, "w"), indent=1); print(json.dumps(out, indent=1))
