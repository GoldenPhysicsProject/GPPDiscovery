# Correction: Lambda is a boundary point of the feasible set, not its only point

Date: 2026-10-04. Author: Claude. Supersedes the conclusions of 2026-10-04_positivity_selects_the_primes.md and 2026-10-04_prime_blind_dual.md.

## What was wrong

The rigidity claim (F_lambda shrinks onto Lambda) is false. Solving the max-min problem directly (prime_blind.py),

    maximize t  subject to  W_inf - sum_n y_n G_n >= t L I,   b(n) = y_n sqrt(n)/2,

returns feasible coefficient vectors far from Lambda:

| lambda | N | best margin t (any b) | margin at Lambda | distance from Lambda |
| --- | --- | --- | --- | --- |
| 3 | 16 | 6.23e-3 | 2.2e-13 | 3.9 |
| 4 | 20 | 6.22e-3 | 2.2e-13 | large (edge coefficients) |
| 6 | 26 | 6.18e-3 | 2.4e-13 | large |
| 8 | 30 | 6.05e-3 | 2.3e-13 | large |

Verified independently with numpy, including on the first 14 modes (the basis of the earlier rigidity runs): the optimizer's b gives min eigenvalue 6.2e-3 there too.

The earlier axis and random-direction tests measured how quickly one exits F_lambda from Lambda. Lambda sits at a sharp corner of the boundary, so nearly every direction exits immediately. The interior is large and was missed. Uniqueness fails, so "RH is equivalent to (E)" does not follow. (E) itself holds with margin at every lambda tested, but through coefficients that are not Lambda.

## What survives

1. **Lambda is a boundary point with essentially zero margin** (min eigenvalue about 1e-13 to 1e-33 depending on precision and basis), consistent with Rodgers-Tao.
2. **Along the prime-strength line it is isolated** (prime_strength.py). With b = c Lambda:

| lambda | c where the form is PSD | min eig at c = 0.95 | at c = 1.01 |
| --- | --- | --- | --- |
| 2 | [0.90, 1.00] | +8.9e-3 | -2.9e-3 |
| 3 | only c = 1 | -1.5e-2 | -1.5e-2 |
| 6 | only c = 1 | -4.0e-2 | -1.9e-2 |
| 10 | only c = 1 | -5.4e-2 | -2.0e-2 |

   Reason: by the explicit formula, A(c Lambda) = Q + (1 - c) P_prime, where Q is the zero side (near-singular) and the prime form is indefinite on Q's near-null space. The primes sit at exactly full strength.
3. The Davenport-Heilbronn controls (2026-10-04_davenport_heilbronn_control.md) and the odd-sector scan (2026-10-04_odd_sector_laplacian_scan.md) are unaffected.
4. The duality statement is correct as mathematics: F_lambda is nonempty exactly when every prime-blind state has nonnegative archimedean energy. It is just not equivalent to RH, because F_lambda contains many points besides Lambda.

## The sharpened question

F_lambda is a large convex set. The primes do not maximize the margin; they sit on the boundary, at a point where the zero side has no slack. The question for a proof becomes: what property singles out that boundary point, and why does it stay inside F_lambda as lambda grows? Candidates to test: multiplicativity of the induced Dirichlet series, integrality of the zero multiplicities, and the functional-equation consistency at every lambda simultaneously.

## Files

prime_blind.py, prime_strength.py.
