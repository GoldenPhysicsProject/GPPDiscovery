"""Two-particle cut in (2,2) split signature (Claude chat seat, 2026-10-05).

Question: is the shadow-locus cut weight P(lam) = pi lam/sinh(pi lam) the Fisher (Kubo-Mori) kernel of the
Bisognano-Wichmann state, kappa_{2pi}(k) = pi k / sinh(pi k), with k the charge of a SINGLE chiral boost?
In Lorentzian signature the boost z -> e^th z moves z and zbar together, so the cut variable
x = log(w5/w6) shifts by 2 th (charge 2 lam, beta = pi). In (2,2) z, zbar are independent reals.

Checks (sympy exact + mpmath numeric):
 1. q(z,zb) = (1+z zb, z+zb, z-zb, 1-z zb), eta = diag(-1,1,-1,1): q.q = 0, q(z).q(w) = -2(z-w)(zb-wb)?? (computed, not assumed)
 2. Solve P = w5 q5 + w6 q6 for P = (M,0,0,0): get z6, zb6, w5, w6 (antipodal map in split signature).
 3. Jacobian of (w5, w6, z6, zb6) -> l5 + l6 and the reduced cut density rho(z5, zb5) with on-shell measure
    w dw dz dzb on each leg (cone measure; its SO(2,2) normalization constant is computed too).
 4. Show rho depends on (z5,zb5) only through u = z5 zb5 times the flat measure of the stabilizer direction
    (a-b, with a=log|z|, b=log|zb|), and that the reduced density in x = -log u is the logistic sech^2(x/2)/4.
 5. Charges: chiral dilation z->e^th z shifts x by th; diagonal (Lorentzian-type) dilation shifts by 2th;
    the anti-diagonal one (z->e^th z, zb->e^-th zb) is the stabilizer (x invariant).
 6. Fourier transform of the logistic in x is P(lam) => P = kappa_{2pi}(k) for the chiral charge k = lam.
"""
import json, sys, os
import sympy as sp
import mpmath as mp

out = {}
z, zb, w, M = sp.symbols('z zb w M', real=True)
eta = sp.diag(-1, 1, -1, 1)
def q(a, b): return sp.Matrix([1 + a*b, a + b, a - b, 1 - a*b])
dot = lambda A, B: sp.expand((A.T * eta * B)[0])
z2, zb2 = sp.symbols('z2 zb2', real=True)
out["q.q"] = str(sp.simplify(dot(q(z, zb), q(z, zb))))
out["q(z).q(w)"] = str(sp.factor(dot(q(z, zb), q(z2, zb2))))

# 2. solve momentum conservation
w5, w6, z5, zb5, z6, zb6 = sp.symbols('w5 w6 z5 zb5 z6 zb6', real=True)
P = sp.Matrix([M, 0, 0, 0])
eqs = list(P - w5*q(z5, zb5) - w6*q(z6, zb6))
sol = sp.solve(eqs, [w5, w6, z6, zb6], dict=True)
out["solutions"] = [{str(k): str(sp.simplify(v)) for k, v in s.items()} for s in sol]

# 3. Jacobian of the constraint map in the eliminated variables
F = (w5*q(z5, zb5) + w6*q(z6, zb6))
J = F.jacobian([w5, w6, z6, zb6])
good = [s for s in sol if sp.simplify(s[w5]) != 0][0]
detJ = sp.simplify(J.det().subs(good))
# cone measure on each leg: d^4 l delta(l^2) = c0 * w dw dz dzb ; get c0 from Jacobian of (w,z,zb,m2)->l
mm = sp.symbols('mm', real=True)
lmap = w*q(z, zb) + sp.Matrix([mm, 0, 0, -mm])  # transverse-ish deformation to compute delta(l^2) jacobian
Jl = lmap.jacobian([w, z, zb, mm]).subs(mm, 0)
l2 = dot(lmap, lmap)
dl2dm = sp.simplify(sp.diff(l2, mm).subs(mm, 0))
c0 = sp.simplify(sp.Abs(Jl.det()) / sp.Abs(dl2dm))
out["cone_measure_factor_c0 (d^4l delta(l^2) = c0 dw dz dzb)"] = str(c0)
rho = sp.simplify((c0.subs({w: w5, z: z5, zb: zb5}) * c0.subs({w: w6, z: z6, zb: zb6})).subs(good) / sp.Abs(detJ))
out["detJ"] = str(detJ)
out["rho(z5,zb5) up to (2pi) conventions"] = str(rho)

# 4. reduce: a = log z5, b = log zb5 (first quadrant), u = e^{a+b}; dz dzb = e^{a+b} da db
a, b, x = sp.symbols('a b x', real=True)
rho_ab = sp.simplify(rho.subs({z5: sp.exp(a), zb5: sp.exp(b)}) * sp.exp(a + b))
out["rho in (a,b) incl. Jacobian"] = str(rho_ab)
s_, d_ = sp.symbols('s d', real=True)  # s=a+b, d=a-b ; da db = ds dd /2
rho_sd = sp.simplify(rho_ab.subs({a: (s_ + d_)/2, b: (s_ - d_)/2}) / 2)
out["depends on d (stabilizer)?"] = str(sp.simplify(sp.diff(rho_sd, d_)) != 0)
rho_x = sp.simplify(rho_sd.subs(s_, -x))
logistic = sp.sech(x/2)**2 / 4
ratio = sp.simplify(sp.rewrite(rho_x / logistic, sp.exp) if hasattr(sp, 'rewrite') else (rho_x/logistic).rewrite(sp.exp))
out["rho_x / logistic"] = str(ratio)

# 5. charges
th = sp.symbols('th', real=True)
X = -sp.log(z5*zb5)
out["shift under chiral z->e^th z"] = str(sp.simplify(sp.expand_log(X.subs(z5, sp.exp(th)*z5) - X, force=True)))
out["shift under diagonal z,zb->e^th"] = str(sp.simplify(sp.expand_log(X.subs({z5: sp.exp(th)*z5, zb5: sp.exp(th)*zb5}) - X, force=True)))
out["shift under anti-diagonal"] = str(sp.simplify(sp.expand_log(X.subs({z5: sp.exp(th)*z5, zb5: sp.exp(-th)*zb5}) - X, force=True)))

# 6. Fourier of logistic = P(lam) = kappa_{2pi}(lam)
mp.mp.dps = 30
errs = []
for lam in (mp.mpf('0.3'), mp.mpf('1.7'), mp.mpf('4.1')):
    ft = mp.quad(lambda t: mp.cos(lam*t) * mp.sech(t/2)**2 / 4, [-mp.inf, 0, mp.inf])
    P = mp.pi*lam/mp.sinh(mp.pi*lam)
    kap = (2*mp.pi*lam/2)/mp.sinh(2*mp.pi*lam/2)
    errs.append(float(max(abs(ft - P), abs(P - kap))))
out["fourier(logistic)=P=kappa_2pi max err"] = max(errs)

p = sys.argv[1] if len(sys.argv) > 1 else "split.json"
os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
json.dump(out, open(p, "w"), indent=1)
print(json.dumps(out, indent=1))
