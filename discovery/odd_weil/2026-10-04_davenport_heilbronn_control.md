# Davenport-Heilbronn control: same mirror, no Euler product, negative energy

Date: 2026-10-04. Author: Claude. Status: double-precision numerics, cross-checked against the zero side. No proof.

## Question

Does the odd semilocal Weil form go negative when the mirror symmetry D is kept but the Euler product is removed?

## Objects

Three Dirichlet series, each with an exact functional equation (the same D):

| name | Euler product | Gamma factor, conductor | zeros |
| --- | --- | --- | --- |
| zeta | yes | Gamma(s/2), q = 1 | on line (known range) |
| L(s, chi_5), real character mod 5 | yes | Gamma(s/2), q = 5 | on line (known range) |
| Davenport-Heilbronn (DH) | no | Gamma((s+1)/2), q = 5 | 52 on-line zeros below 100; off-line zero at 0.808517 + 85.699348 i |

DH coefficients by n mod 5: 1, kappa, -kappa, -1, 0 with kappa = (sqrt(10 - 2 sqrt 5) - 2)/(sqrt 5 - 1). Its "von Mangoldt" coefficients b(n), from -f'/f = sum b(n) n^{-s}, live on composites (6, 12, 14, 18, 21, 26, 28, ...) and have the wrong size on prime powers (b(4) = -1.44, b(9) = -2.29 versus Lambda = 0.69, 1.10).

The form is the zero-free side, built from b(n) and the Gamma factor only:

    A(f) = [log(q/pi) - gamma_E] G(0)
           + int_0^inf 2[e^{-2u} G(0) - e^{-(1/2 + k) u} G(u)]/(1 - e^{-2u}) du
           - 2 sum_{n < lambda^2} b(n) n^{-1/2} G(log n),

with k = 0 for Gamma(s/2) and k = 1 for Gamma((s+1)/2). Basis: odd sines on [-log lambda, log lambda], frequencies up to about 110.

## Validation

- The zeta case reproduces the 50-digit value at lambda = 2 to nine digits (1.36841509e-5).
- The u-space Archimedean kernel for Gamma((s+1)/2) matches the r-space psi(3/4 + ir/2) integral to 1.7e-8.
- At lambda = 4 the DH near-null vector's Fourier transform vanishes at the DH on-line zeros (errors 1.7e-10, 1.7e-9, 3.3e-8, ...). The DH matrix is the true DH form.

## Results: smallest eigenvalue (double precision; noise floor about 2e-13)

| lambda | N | zeta | L(chi_5) | DH |
| --- | --- | --- | --- | --- |
| 3 | 39 | 2.1e-13 | 3.4e-5 | 4.8e-4 |
| 4 | 49 | 2.1e-13 | 3.8e-12 | 9.4e-11 |
| 6 | 63 | 1.9e-13 | 1.9e-13 | -3.2e-9 |
| 8 | 73 | 1.6e-13 | 1.9e-13 | -1.01 |
| 10 | 81 | 2.2e-13 | 2.0e-13 | -1.38 |
| 12 | 88 | 2.2e-13 | 2.2e-13 | -1.78 |

Both Euler-product functions stay at the positive noise floor. DH goes decisively negative once lambda >= 8.

## The negative direction is the off-line zero

At lambda = 8 the negative eigenvector's spectrum peaks at r = 84.6. Evaluating the zero side on it:

| lambda | eigenvalue x L | on-line zeros < 100 | off-line quartet | sum |
| --- | --- | --- | --- | --- |
| 8 | -2.104 | +0.116 | -2.255 | -2.139 |
| 10 | -3.184 | +0.143 | -3.336 | -3.193 |

The zero-free matrix, built without zeros, found the off-line zero at height 85.7 and reports its negative energy. The small gap is zeros above 100.

## Which places drive it

Energy of the negative vector, split by source:

| lambda | A(v) | Archimedean | prime-power n | composite n |
| --- | --- | --- | --- | --- |
| 8 | -2.10 | +8.72 | -4.41 | -6.42 |
| 10 | -3.18 | +9.70 | -3.71 | -9.17 |

The most negative single terms come from n = 14, 21, 6, 26 (composites, which have no Euler-product counterpart) and n = 9 (a prime power with the wrong weight). Dropping the composite terms would leave the Archimedean term ahead (8.72 - 4.41 > 0). That is an illustration of where the imbalance sits, not a valid operation on a real L-function.

## Reading

- The mirror symmetry D is identical across all three. Only DH has a negative mode. Symmetry alone does not place zeros on the line, now shown numerically on the same matrices.
- The failure is located: the "prime" side of DH carries terms at composite n and wrong prime-power weights. The Euler product, which forces b(n) = Lambda(n) chi(n), is what keeps the prime side balanced against the Archimedean place.
- In the arrow language: DH has the global reversal symmetry (Theorem 4.1 analog) but violates the place-by-place structure, and a D-odd observable (the off-line quartet) appears. That is Hypothesis H2 failing, detected by a zero-free computation.

## Files

dh_control.py (matrix builder, three kinds), dh_validate.py, dh_zeros.py, dh_scan.py, dh_explain.py, dh_breakdown.py.
