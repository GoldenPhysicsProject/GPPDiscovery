"""Celestial shadow sewing: numerical checks of the exact pieces (Claude chat seat, 2026-10-04).

Recovered thread: one-loop integrands from the Delta_5 + Delta_6 = 2 shadow locus of a six-point
tree (bridge: CLAUDE_CELESTIAL_LOOPS.md). This script checks only the pieces that are exact:

 1. P(lam) = Gamma(1+i lam) Gamma(1-i lam) = pi lam / sinh(pi lam)       (not pi/sinh: old-chat error)
 2. P_hat(xi) = int P(lam) e^{-i lam xi} dlam = (pi/2) sech^2(xi/2)
 3. Mellin-Plancherel sewing: for f, g on (0,inf) with Mellin transforms f~(D) = int w^{D-1} f dw,
        int_0^inf w^k f(w) g(w) dw = int dlam/2pi  f~(c+i lam) g~(k+1-c-i lam)
    i.e. sewing with measure w^k dw pairs dimensions on D5 + D6 = k+1. Massless Lorentz-invariant
    phase space d^3p/(2E) ~ w dw d^2z has k = 1, so the pairing locus is exactly D5 + D6 = 2,
    the shadow locus. Checked for k = 0, 1, 2 with f = e^{-a w}, g = e^{-b w}.
 4. For k = 1, f = e^{-a w}, g = e^{-b w}: the spectral integrand is exactly P(lam) a^{-1-i lam} b^{-1+i lam},
    so the sewing weight P(lam) appears by itself, and the closed form is 1/(a+b)^2.

Nothing here is a proof of the loop-extraction claim; it checks the measure-level identities it uses.
"""
import json, sys
import mpmath as mp
mp.mp.dps = 30
out = {}

# 1
lams = [mp.mpf(x) for x in ("0.1", "0.5", "1", "2.5", "7")]
out["gamma_product_max_relerr"] = float(max(abs(mp.gamma(1+1j*l)*mp.gamma(1-1j*l) - mp.pi*l/mp.sinh(mp.pi*l))/(mp.pi*l/mp.sinh(mp.pi*l)) for l in lams))
out["old_chat_form_pi_over_sinh_relerr_at_lam1"] = float(abs(mp.gamma(1+1j)*mp.gamma(1-1j) - mp.pi/mp.sinh(mp.pi))/(mp.pi/mp.sinh(mp.pi)))

# 2
P = lambda l: mp.pi*l/mp.sinh(mp.pi*l) if l != 0 else mp.mpf(1)
errs = []
for xi in (0, 0.7, 2.0, 5.0):
    num = 2*mp.quad(lambda l: P(l)*mp.cos(l*xi), [0, mp.inf])
    errs.append(abs(num - mp.pi/2*mp.sech(xi/2)**2))
out["P_hat_max_abserr"] = float(max(errs))

# 3, 4
rows = []
for k in (0, 1, 2):
    for (a, b) in ((1, 2), (1, 1), (0.3, 4.0)):
        a, b = mp.mpf(a), mp.mpf(b)
        lhs = mp.gamma(k+1)/(a+b)**(k+1)
        c = mp.mpf(k+1)/2   # symmetric contour: Re D5 = Re D6 = (k+1)/2
        integrand = lambda l: (mp.gamma(c+1j*l)*a**(-(c+1j*l)) * mp.gamma(k+1-c-1j*l)*b**(-(k+1-c-1j*l))).real
        rhs = mp.quad(integrand, [-mp.inf, 0, mp.inf])/(2*mp.pi)
        rows.append({"k": k, "a": float(a), "b": float(b), "pair_locus_D5+D6": k+1,
                     "lhs": float(lhs), "rhs": float(rhs), "relerr": float(abs(lhs-rhs)/lhs)})
out["plancherel_sewing"] = rows
out["k1_closed_form_check"] = [float(abs((1/(r["a"]+r["b"])**2) - r["lhs"])) for r in rows if r["k"] == 1]

ok = out["gamma_product_max_relerr"] < 1e-20 and out["P_hat_max_abserr"] < 1e-12 and all(r["relerr"] < 1e-12 for r in rows)
out["all_pass"] = bool(ok)
json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "sewing.json", "w"), indent=1)
print(json.dumps(out, indent=1))
sys.exit(0 if ok else 1)
