# Odd-sector semilocal Weil form: Laplacian split and zero detection

Date: 2026-10-04. Author: Claude. Status: numerical experiment at 50 digits. No proof; no interval certification.

## Setup

- Test space: real odd f on [-L, L], L = log(lambda), basis sin(k pi u / L), k = 1..14.
- Form: A(f) = Q(f) - P(f), the prime-plus-Archimedean side of the explicit formula, built with **no zero data**:

      A(f) = W_inf(f) - 2 sum_{n < lambda^2} Lambda(n) n^{-1/2} G(log n),
      G(u) = int f(u+b) f(b) db,
      W_inf(f) = (-log pi - gamma_E) G(0)
                 + int_0^inf 2[e^{-2u} G(0) - e^{-u/2} G(u)] / (1 - e^{-2u}) du.

  The u-space Archimedean kernel comes from psi(s) = -gamma + int_0^inf (e^{-t} - e^{-st})/(1 - e^{-t}) dt with t = 2u.
- Laplacian split: with c_n = Lambda(n)/sqrt(n) and S = sum c_n,

      A = Lap + Mass,   Lap(f) = sum c_n ||f - T_n f||^2 >= 0,   Mass(f) = W_inf(f) - 2 S ||f||^2.

## Validation

The zero-free matrix was compared entrywise with the zero side, 2 sum_{gamma>0} |F(i gamma)|^2 + 2 F(1/2)^2, using the first 1500 zeta zeros (validate.py). Maximum discrepancy 1.3e-7 at lambda = 2 and 5.2e-8 at lambda = 3, consistent with zero truncation. The normalization is correct.

## Results (mp.dps = 50, N = 14)

| lambda | prime powers used | min eig A | 2nd eig A | min eig Mass |
| --- | --- | --- | --- | --- |
| 2.0 | 2, 3 | 1.37e-5 | 5.2e-2 | -2.86 |
| 2.5 | 2, 3, 4, 5 | 3.18e-14 | 2.2e-9 | -5.29 |
| 3.0 | 2 to 8 | 3.80e-20 | 2.5e-15 | -7.45 |
| 4.0 | 2 to 13 | 5.21e-27 | 2.8e-22 | -11.32 |
| 5.0 | 2 to 23 | 4.62e-31 | 2.3e-26 | -15.89 |
| 6.0 | 2 to 32 | 1.93e-33 | 1.2e-28 | -19.83 |

Eigenvalues are normalized by the Gram factor L. The slowdown at lambda = 5, 6 is basis saturation at N = 14, not a change in the arithmetic.

1. **A is positive at every lambda tested**, with a margin that collapses super-exponentially. Consistent with Rodgers-Tao zero slack.
2. **The mass term is badly negative and the odd class does not tame it.** Its worst direction is the lowest mode k = 1, and it grows like -3.3 lambda, tracking -2S.
3. **The passive Laplacian compensates with zero slack.** The largest generalized eigenvalue mu of (-Mass, Lap) is 0.99999307 at lambda = 2 and equals 1 to 12 digits for lambda >= 2.5. Positivity is an exact balance between the positive network and the negative mass, attained on one direction.
4. **That direction detects the zeta zeros.** Let v be the near-null eigenvector. The zeros of its Fourier transform F_v(ir), excluding the trivial grid m pi/L, match the zeta ordinates:

| lambda | gamma1 err | gamma2 err | gamma3 err | gamma4 err | gamma5 err | gamma6 err |
| --- | --- | --- | --- | --- | --- | --- |
| 2.0 | 1.0e-4 | 1.6e-3 | 8.9e-3 | none | none | none |
| 2.5 | 4.5e-13 | 2.6e-11 | 4.2e-10 | 3.9e-8 | 2.8e-7 | 6.1e-6 |
| 3.0 | 3.6e-15 | 1.1e-14 | 7.1e-15 | 6.1e-12 | 1.3e-10 | 5.9e-8 |
| 4.0 | 0 | 3.6e-15 | 2.1e-14 | 4.3e-11 | 4.1e-9 | 1.1e-4 |

At lambda = 3, prime powers up to 8 and fourteen sine modes reproduce gamma1 to gamma3 to double precision. This is the Connes-Consani zeta-cycle phenomenon, here in the odd, pole-free sector.

## Reading

- Point 4 is forced by the zero-side expansion: on odd f, A = 2 sum |F(i gamma)|^2 + 2 F(1/2)^2, so a near-null vector must nearly vanish at the low zeros and at s = 1/2. It certifies that the matrices are right. It is not independent evidence for RH beyond the lambdas computed.
- Points 2 and 3 settle the question asked: no termwise argument of the form "positive network plus controlled mass" can work. The mass is not controlled; it is exactly cancelled, and only on the zero-detecting direction.
- Constructive consequence for the Schur/network route: the parent network cannot be generic. Its ground state must be the zero-detecting vector v. The concrete target becomes an identity exhibiting v (or its lambda -> infinity limit) as the harmonic/ground state of a positive operator built from primes and Gamma. The Connes-Consani-Moscovici prolate operator is the existing candidate; compare v against its eigenfunctions next.

## Files

odd_weil.py (double precision), mp_weil.py and run_mp.py (50-digit build), validate.py (zero-side check), analyze.py (Fourier zeros, Laplacian margin), zerr.py (zero errors).
