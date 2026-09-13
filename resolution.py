#!/usr/bin/env python3
"""The archimedean truncation criterion for the CvS Galerkin scan.

WHAT THIS IS

`connes_cvs.build_galerkin_matrix(c, N, T, dps)` computes the archimedean
piece of each matrix entry as

    psi_arch(x) = (1/2pi^2) * INT_{-T}^{T} h_plus(tau) * Re S_hat(tau, x) d_tau

with the quadrature split at the integrand's kinks, tau = 0 and tau = +/- alpha_x,
where

    alpha_x = 2*pi*x / log(c)

and x runs over the basis indices -N .. N.  (See connes_cvs/operator.py,
`psi_arch`: the split list is literally {0, alpha_x, -alpha_x} filtered to
-T < s < T.)

So the feature that carries basis mode x's archimedean weight sits at spectral
height alpha_x.  The truncation at |tau| <= T keeps it only while alpha_x <= T.
For the outermost mode x = N that is

    2*pi*N / log(c)  <=  T          equivalently   N <= T*log(c)/(2*pi)

Past that point the top modes are computed with their archimedean contribution
truncated away, the Galerkin matrix is not the matrix anyone intended, and
lambda_1 can and does go negative.

WHY IT MATTERS FOR THE SCAN

1. T is NOT a free "make it big enough" parameter.  It is tied to N and c.
   Holding T fixed while raising N walks off the cliff, and it does so SOONER
   FOR SMALLER c, because log(c) is in the denominator.  That is why c=3 is the
   value that has always misbehaved.

2. A c-scan at fixed (N, T) does not have uniform numerical quality.  At
   N=36, T=300 the margin T*log(c)/(2*pi*N) runs from 1.46 at c=3 to 3.90 at
   c=19.  A slope fitted across that range is fitting a varying truncation
   quality as well as the mathematics.

   The fix: hold the MARGIN fixed, not T.  Choose T = kappa * 2*pi*N/log(c)
   for a fixed kappa (>= 2, say), so every point in the scan sits the same
   distance from the cliff.  `plan_T` below does that.

STATUS: this is a numerical-validity criterion read off the implementation and
checked against every point on record.  It is not a theorem about the operator.
"""
import math


def tau_max(c: float, N: int) -> float:
    """Spectral height of the outermost basis mode's archimedean kink."""
    return 2 * math.pi * N / math.log(c)


def n_max(c: float, T: float) -> float:
    """Largest N whose archimedean kink still lies inside [-T, T]."""
    return T * math.log(c) / (2 * math.pi)


def margin(c: float, N: int, T: float) -> float:
    """T / tau_max.  >= 1 is inside the window; the larger the safer."""
    return T / tau_max(c, N)


def within_window(c: float, N: int, T: float) -> bool:
    return margin(c, N, T) >= 1.0


def plan_T(c: float, N: int, kappa: float = 3.0) -> int:
    """T giving a fixed margin kappa, so a c-scan has uniform quality."""
    return int(math.ceil(kappa * tau_max(c, N)))


if __name__ == "__main__":
    print("Every (c, N, T) on record in this repo, against the criterion:\n")
    print(f"{'c':>3} {'N':>3} {'T':>4} {'tau_max':>9} {'margin':>7}  {'verdict':<8} observed lambda_1")
    recorded = [
        (3, 20, 100, "-0.332          (README, called 'marginal basis resolution')"),
        (3, 12, 60, "-0.2444         (2026-09-13)"),
        (3, 16, 60, "-0.4855         (2026-09-13)"),
        (3, 20, 60, "-0.4901         (2026-09-13)"),
        (3, 36, 120, "-0.4901         (2026-09-13)"),
        (3, 36, 160, "-0.4897         (2026-09-13)"),
        (3, 36, 190, "-0.3779         (2026-09-13)"),
        (3, 36, 210, "+4.656e-8       (2026-09-13)"),
        (3, 36, 300, "+4.951e-8"),
        (3, 44, 300, "+4.939e-8       (2026-09-13)"),
        (5, 36, 300, "+8.904e-18"),
        (5, 44, 300, "+8.616e-18      (2026-09-13)"),
        (7, 36, 300, "+8.093e-28"),
        (11, 36, 300, "+5.191e-46"),
        (13, 36, 300, "+5.602e-52"),
        (17, 36, 300, "+3.242e-60"),
        (19, 36, 300, "+2.661e-63"),
        (19, 44, 300, "+1.038e-69"),
        (19, 52, 300, "+2.949e-75"),
    ]
    for c, N, T, obs in recorded:
        m = margin(c, N, T)
        print(f"{c:>3} {N:>3} {T:>4} {tau_max(c,N):9.1f} {m:7.2f}  "
              f"{'PAST' if m < 1 else 'within':<8} {obs}")
    print("\nA uniform-margin c-scan at N=36, kappa=3:\n")
    print(f"{'c':>3} {'T':>6}")
    for c in (3, 5, 7, 11, 13, 17, 19, 23, 29):
        print(f"{c:>3} {plan_T(c, 36):>6}")
