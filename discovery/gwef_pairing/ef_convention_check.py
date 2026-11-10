"""GW-EF-0 / GW-EF-1 convention check (NUMERICAL, run 2026-11-10).

The tree's weilForm  g = poleTerm g + archimedeanTerm g - primeSide g
  poleTerm g        = int g(u) 2cosh(u/2) du
  archimedeanTerm g = -(gamma + log pi) g(0) + 2 int_0^inf (e^{-2u} g(0) - e^{-u/2} g(u))/(1 - e^{-2u}) du
  primeSide g       = sum_n Lambda(n)/sqrt(n) (g(log n) + g(-log n))
is compared with the zero side  sum_rho Phi g(rho),  Phi g(s) = int g(u) e^{(s-1/2)u} du,
for g = autocorr f, f = exp(-y^2/(2a^2))  (g(x) = a sqrt(pi) exp(-x^2/(4a^2))).
On the line, Phi g(1/2 + i t) = 2 pi a^2 exp(-a^2 t^2).  Zeros: mpmath zetazero (all simple, on the line
for the range used), counted as +gamma and -gamma.  This tests the identity only at on-line zeros.
"""
from mpmath import mp, mpf, zetazero, exp, sqrt, pi, log, quad, euler, cosh, inf, mangoldt
mp.dps = 30
def run(a, nz):
    a = mpf(a)
    g = lambda x: a*sqrt(pi)*exp(-x**2/(4*a**2))
    pole = quad(lambda u: g(u)*2*cosh(u/2), [-inf, 0, inf])
    arch_int = quad(lambda u: (exp(-2*u)*g(0)-exp(-u/2)*g(u))/(1-exp(-2*u)), [0, 1, 5, inf])
    arch = -(euler+log(pi))*g(0) + 2*arch_int
    prime = sum(mangoldt(n)/sqrt(n)*(g(log(n))+g(-log(n))) for n in range(2, 4000) if mangoldt(n) != 0)
    weil = pole + arch - prime
    Phi = lambda t: 2*pi*a**2*exp(-a**2*t**2)
    zs = sum(2*Phi(zetazero(j).imag) for j in range(1, nz+1))
    return pole, arch, prime, weil, zs
if __name__ == "__main__":
    for a, nz in [(0.2, 40), (0.1, 80), (0.05, 160)]:
        pole, arch, prime, weil, zs = run(a, nz)
        print(f"a={a}: pole={float(pole):.10f} arch={float(arch):.10f} prime={float(prime):.3e} "
              f"weil={float(weil):.12f} zeros={float(zs):.12f} diff={float(weil-zs):.2e}")
