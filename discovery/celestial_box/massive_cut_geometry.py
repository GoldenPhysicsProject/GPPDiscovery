"""Equal-mass two-particle cut geometry in celestial stereographic coordinates.

This is the correct first step for D-dimensional unitarity: a D-dimensional
massless cut momentum appears in its four-dimensional component as a massive
momentum with mass mu.  The angular cut geometry remains a round celestial
sphere; the only phase-space change is beta = sqrt(1-4 mu^2/M^2).

This script verifies the exact COM parametrization, antipodal pairing, and
normalization

    dPhi_2 = beta * d^2 z / [8 pi^2 (1+|z|^2)^2],
    int dPhi_2 = beta/(8 pi).

It intentionally does NOT identify the massive spectral transform with the
massless Mellin Gamma/Beta weight: massive conformal-primary wavefunctions are
H^3 harmonic transforms, not simple energy Mellin modes.
"""
import mpmath as mp
mp.mp.dps = 50


def beta(M, mu):
    return mp.sqrt(1 - 4 * mu * mu / (M * M))


def nvec(z):
    r2 = abs(z) ** 2
    return mp.matrix([
        (z + z.conjugate()).real / (1 + r2),
        (-1j * (z - z.conjugate())).real / (1 + r2),
        (1 - r2) / (1 + r2),
    ])


def antipode(z):
    return -1 / z.conjugate()


def minkowski_sq(p):
    return p[0] ** 2 - sum(p[i] ** 2 for i in range(1, 4))


def momentum(M, mu, z):
    E = M / 2
    k = E * beta(M, mu)
    n = nvec(z)
    return mp.matrix([E, k * n[0], k * n[1], k * n[2]])


for M, mu, z in [
    (mp.mpf('5.0'), mp.mpf('0.7'), mp.mpc('0.3', '0.4')),
    (mp.mpf('3.7'), mp.mpf('0.2'), mp.mpc('-0.8', '0.6')),
    (mp.mpf('8.1'), mp.mpf('1.3'), mp.mpc('1.1', '-0.2')),
]:
    p5 = momentum(M, mu, z)
    p6 = momentum(M, mu, antipode(z))
    assert abs(minkowski_sq(p5) - mu**2) < mp.mpf('1e-45')
    assert abs(minkowski_sq(p6) - mu**2) < mp.mpf('1e-45')
    total = p5 + p6
    assert abs(total[0] - M) < mp.mpf('1e-45')
    assert max(abs(total[i]) for i in range(1, 4)) < mp.mpf('1e-45')

# Round-sphere normalization in stereographic coordinates.
M = mp.mpf('5.0'); mu = mp.mpf('0.7')
b = beta(M, mu)
I = mp.quad(lambda r: b * (2 * mp.pi * r) /
            (8 * mp.pi**2 * (1 + r*r)**2), [0, mp.inf])
assert abs(I - b/(8*mp.pi)) < mp.mpf('1e-45')

# Massless limit of the geometric measure only.
for mu0 in [mp.mpf('1e-2'), mp.mpf('1e-4'), mp.mpf('1e-8')]:
    assert abs(beta(M, mu0) - 1) < 2 * mu0**2

print('PASS: p5^2 = p6^2 = mu^2 at generic celestial points')
print('PASS: antipodal z6=-1/conj(z5) gives p5+p6=(M,0,0,0)')
print('PASS: int dPhi2 = beta/(8*pi)')
print('PASS: geometric cut measure tends to the massless round-sphere measure')
