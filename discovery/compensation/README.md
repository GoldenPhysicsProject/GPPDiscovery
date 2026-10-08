# compensation: radial Euler compensation against the shifted-product control

Status: numerical. Nothing here is proved. Floating point, double precision, N modes (odd sine basis on [-L, L], L = log lambda), noise floor about 3e-13.

## What is computed

`ftheta.py` builds the odd fixed-window Weil form of F_theta(s) = xi(s - theta) xi(s + theta) (theta = 0 is zeta squared up to normalization, so the form is twice zeta's), from the prime side only:

    Q = 2 (W_inf - sum_n Lambda(n) n^{-1/2} cosh(theta log n) 2 G_n + Pole),

where `G_n` is the translation autocorrelation at log n, `W_inf` the archimedean form, and `Pole` the rank-two term from s = theta + 1/2 and theta - 1/2. Radial compensation replaces Lambda(n) by Lambda(n)(1 - 1/n), i.e. adds `Corr = 2 sum Lambda(n) n^{-3/2} cosh(theta log n) 2 G_n` (the `Qc` and `Ac` columns).

Exact reduction used (derivation: zero side of F_theta on h = f * f~ equals 2 E_zeta[h(u) cosh(theta u)], a cross term of zeta's zero form, 2 Re Q_zeta(e^{-theta x} f, e^{theta x} f)). Taken on trust from the zero side only after the check below.

## Validation (done)

`validate_zero_side.py` compares `Q` against the zero side summed over the first 2000 zeta zeros (N = 8):

| theta | lambda | max abs diff, scale | min eig, zero side | min eig, prime side |
| --- | --- | --- | --- | --- |
| 0.0 | 3 | 9.3e-8, 4.9 | 1.7e-16 | 3.8e-13 |
| 0.2 | 3 | 9.4e-8, 5.0 | -4.534e-2 | -4.534e-2 |
| 0.3 | 4 | 5.7e-8, 6.4 | -2.280e-1 | -2.280e-1 |

The agreement is limited by the 2000-zero truncation. Output: `validate_zero_side.out`.

## Scan

`ftheta_scan.out`: theta in {0, .1, .2, .3, .375}, lambda in {2, 3, 4, 6, 8}, N = ceil(60 L / pi). Columns: min eigenvalue of Q/L, of the compensated form, of the pole-free A form and its compensated version, and the operator norm of Corr/L.

## Readings, with status

1. Exact, checked numerically: the local identity J_p(s) - J_p(s+1) = log p sum_k (1 - p^{-k}) p^{-ks} (residual 1e-28 at 30 digits), and the global sum equals -d/ds log(zeta(s)/zeta(s+1)).
2. Numerical: the uncompensated F_theta form goes negative for every theta > 0 tested and is at the noise floor for theta = 0. The onset deepens with lambda.
3. Numerical, expected: the compensated form is not PSD even at theta = 0 (min eigenvalue about -0.6 to -0.95). Compensation removes the n^{-3/2} part of the prime side, which is not a positive form (`Corr` has norm about 0.7 to 2 in units of L, bounded above by 2 * (-zeta'/zeta)(3/2) = 3.01). So the positive-coefficient object R_zeta = zeta(s)/zeta(s+1) does not carry a positive Weil form by itself; the stable shifted remainder must be kept in any positivity argument.
4. Not established: that compensation separates zeta from F_theta. At theta = 0.375 the compensated minimum is lower than the plain one at lambda = 2, 3, higher at lambda = 4, and lower again at lambda = 6, 8 (-2.285 vs -2.238, -3.564 vs -3.503). The effect on F_theta is mixed in sign and small relative to the form itself, so on this family compensation does not act as a discriminator.

## Controls and caveats

Davenport-Heilbronn is not yet run through the compensated form (needs its own archimedean place and non-Euler coefficients; `../odd_weil/dh_control.py` has them). The theta grid and lambda grid are one sampling convention; the cheapest separating test is a theta very close to 0 and a lambda between the tested values.
