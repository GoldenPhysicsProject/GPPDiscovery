"""Test the 'two halves of the hologram' split of xi.

Riemann's proof splits the Mellin integral of psi(x)=sum exp(-pi n^2 x) at x=1
(the fixed point of the inversion x -> 1/x).  This gives exactly
    xi(s) = E(s) + E(1-s),
    E(s)  = 1/4 + (1/2) s(s-1) sum_n (pi n^2)^(-s/2) Gamma(s/2, pi n^2).
If |E(1-s)| < |E(s)| on Re s > 1/2 then xi has no zeros there (RH).
This script is a falsifier: it looks for points violating that inequality and
for zeros of E in the right half-plane.
"""
import mpmath as mp
mp.mp.dps = 30

def E(s, N=8):
    tot = mp.mpf(0)
    for n in range(1, N + 1):
        a = mp.pi * n * n
        tot += a ** (-s / 2) * mp.gammainc(s / 2, a)
    return mp.mpf(1) / 4 + s * (s - 1) / 2 * tot

def xi(s):
    return s * (s - 1) / 2 * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

# sanity: the split reproduces xi
for s in [mp.mpc(2, 0), mp.mpc(0.5, 14.134725), mp.mpc(0.7, 30), mp.mpc(0.3, 5)]:
    print("check", s, mp.nstr(abs(E(s) + E(1 - s) - xi(s)), 5), mp.nstr(abs(xi(s)), 5))

worst = []
for sig in [0.51, 0.55, 0.6, 0.75, 1.0, 1.25, 1.5]:
    bad = 0; mx = 0; mxt = None
    t = 0.0
    while t <= 120:
        s = mp.mpc(sig, t)
        r = abs(E(1 - s)) / abs(E(s))
        if r > mx:
            mx, mxt = r, t
        if r >= 1:
            bad += 1
        t += 0.25
    print(f"sigma={sig}: max |E(1-s)|/|E(s)| = {mp.nstr(mx, 8)} at t={mxt}; violations={bad}")
