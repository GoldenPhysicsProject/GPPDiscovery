# kappa (zero side) against the Q(L)/L plateau (Moebius side) (NUMERICAL)

Two independent computations of the same constant, the mean square of B(u) per unit u: the sieve gives the slope of Q(L) = integral_0^L B^2; the zeros give kappa = 2 sum_(gamma>0) 1/(gamma^2 |zeta'(rho)|^2). Agreement is evidence for the almost-periodic (explicit-formula) picture and says nothing about proving RH.

| quantity | value | note |
|---|---|---|
| kappa (zeros, 100,000 of them, gamma <= 74921) | **0.029057** +- 4.8e-07 | partial 0.029055 + tail 2.63e-06 |
| Q(L) least-squares slope (sieve, N from 1e4 to 100,000,000) | **0.028949** | residual rms 0.0016 (almost-periodic wobble of Q around the line) |
| Q/L at N = 100,000,000 | 0.159073 | approaches the slope like intercept/L, intercept 2.3957 |
| Q/L^2 at N = 100,000,000 | 0.008636 | goes to 0: linear, not quadratic, growth |

Difference slope - kappa = -0.000109 (-0.37% of kappa).

## Local slopes of Q(L) over one decade of N (sieve)

| N from | N to | slope |
|---|---|---|
| 1,000 | 10,000 | 0.02990 |
| 1,778 | 17,783 | 0.02792 |
| 3,162 | 31,623 | 0.02905 |
| 5,623 | 56,234 | 0.02846 |
| 10,000 | 100,000 | 0.02980 |
| 17,783 | 177,828 | 0.03113 |
| 31,623 | 316,228 | 0.02845 |
| 56,234 | 562,341 | 0.02918 |
| 100,000 | 1,000,000 | 0.02819 |
| 177,828 | 1,778,279 | 0.02709 |
| 316,228 | 3,162,278 | 0.02897 |
| 562,341 | 5,623,413 | 0.02968 |
| 1,000,000 | 10,000,000 | 0.02854 |
| 1,778,279 | 17,782,794 | 0.03051 |
| 3,162,278 | 31,622,777 | 0.02810 |
| 5,623,413 | 56,234,133 | 0.02830 |
| 10,000,000 | 100,000,000 | 0.02943 |
