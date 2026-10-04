# Global check of the Euler-manifold boundary claim

Date: 2026-10-04. Author: Claude. Tests Codex's claim (GPP-bridge research/codex/2026-10-04_euler_local_boundary_after_correction.md) that the degree-one Euler manifold meets the Weil-positive cone only at the zeta point. Codex's evidence at lambda >= 4 was local (small perturbations around alpha = 1), which is the same corner test that produced the retracted rigidity claim. This note replaces it with multi-start global searches.

## Real local roots (self-dual family): b(p^k) = log p * alpha_p^k

euler_global.py, Powell from 61 starts (alpha = 1 plus 60 random starts with 30% spread), odd basis, zeta Archimedean place:

| lambda | primes | N | margin at alpha = 1 | best margin found | best alpha |
| --- | --- | --- | --- | --- | --- |
| 3 | 2, 3, 5, 7 | 16 | 2.2e-13 | 4.28e-3 | (0.9385, 0.8893, 0.8899, 1.0575) |
| 4 | to 13 | 20 | 2.2e-13 | 2.2e-13 | all 1.0000 |
| 5 | to 23 | 22 | 2.3e-13 | 2.3e-13 | all 1.0000 (alpha_23 = 1.0002) |
| 6 | to 31 | 26 | 2.3e-13 | 2.4e-13 | all 1.0000 |

The lambda = 3 interior point reproduces Codex's to three decimals. From lambda = 4 on, every start converges to alpha_p = 1 at the noise floor. **Codex's claim survives the global test for the real (self-dual) family.**

## Complex local roots: alpha_p = r_p e^{i theta_p}, b(p^k) = log p * r_p^k cos(k theta_p)

euler_global_complex.py, 41 starts with random phases:

| lambda | best margin found | location |
| --- | --- | --- |
| 4 | 5.34e-3 | far from zeta: phases up to 3.0 rad, r_13 = -1.95 |
| 5 | 2.3e-13 | all r = 1, theta = 0 |

Allowing phases reopens a positive interior at lambda = 4. At lambda = 5 none was found, but 41 starts in 18 dimensions is weak evidence. Phases give two parameters per prime, so a prime with few visible powers becomes nearly free.

## Reading

- The boundary picture needs self-duality, not just an Euler product. Real alpha_p is the self-dual case (the L-function equals its own conjugate). In the arrow language this is the D-fixed real form imposed prime by prime: Fix(D) at the level of local data.
- Precise current statement (numerical, multi-start, not certified): for lambda in {4, 5, 6}, the only self-dual degree-one Euler point with nonnegative odd Weil margin at the zeta Archimedean place is alpha_p = 1 for all visible primes, and it sits at zero margin.
- Next: repeat the complex search at lambda = 5, 6 with many more starts, and test whether self-duality alone (without the Euler product) already closes the interior.
