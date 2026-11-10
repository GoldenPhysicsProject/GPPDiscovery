"""E6 sanity check: the repo's archimedeanTerm vs the digamma line integral.

archimedeanTerm g = -(gamma + log pi) g(0) + 2 * int_0^inf (e^{-2u} g(0) - e^{-u/2} g(u)) / (1 - e^{-2u}) du
Weil/Bombieri form  = -g(0) log pi + (1/2pi) int h(r) Re psi(1/4 + i r/2) dr,  h(r) = int g(u) e^{iru} du.
Evidence only (not a proof): numeric agreement for even test functions.
"""
from mpmath import mp, mpf, quad, exp, euler, log, pi, digamma, re, inf, sqrt, cosh, cos
mp.dps = 30

def check(g, h, name):
    g0 = g(mpf(0))
    arch = -(euler + log(pi)) * g0 + 2 * quad(lambda u: (exp(-2*u) * g0 - exp(-u/2) * g(u)) / (1 - exp(-2*u)), [0, 1, 10, inf])
    weil = -g0 * log(pi) + quad(lambda r: h(r) * re(digamma(mpf(1)/4 + 1j*r/2)), [-inf, -10, 0, 10, inf]) / (2*pi)
    print(f"{name:28s} archimedeanTerm = {mp.nstr(arch, 20)}   Weil digamma form = {mp.nstr(weil, 20)}   diff = {mp.nstr(abs(arch-weil), 5)}")

for s in (mpf(1), mpf('0.6'), mpf('1.7')):
    check(lambda u, s=s: exp(-u*u/(2*s*s)),
          lambda r, s=s: s*sqrt(2*pi)*exp(-(s*r)**2/2), f"Gaussian sigma={s}")
# even, non-Gaussian: g(u) = 1/cosh(u)^2 ; h(r) = pi r / sinh(pi r / 2)  (h(0)=2)
def hsech2(r):
    return mpf(2) if r == 0 else pi*r/mp.sinh(pi*r/2)
check(lambda u: 1/cosh(u)**2, hsech2, "sech^2")
