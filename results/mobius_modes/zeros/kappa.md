# kappa = 2 sum_{gamma>0} 1/(gamma^2 |zeta'(rho)|^2)  (NUMERICAL)

Zeros from `mpmath.zetazero` (fallback: Platt's tables were not available to the runner); zeta'(rho) by mpmath at 32 digits (residual |zeta(rho)| < 1e-20 sampled every 500 zeros). Assumes RH only through the labelling of rho as 1/2 + i gamma; the tail estimate additionally assumes a stable mean density of 1/|zeta'|^2 (UNPROVEN).

| zeros | gamma_max | partial sum |
|---|---|---|
| 10 | 49.8 | 0.0253351395 |
| 30 | 101.3 | 0.0272089502 |
| 100 | 236.5 | 0.0283107663 |
| 300 | 541.8 | 0.0287224837 |
| 1,000 | 1419.4 | 0.0289314726 |
| 3,000 | 3533.3 | 0.0290033199 |
| 10,000 | 9877.8 | 0.0290382846 |
| 30,000 | 25755.3 | 0.0290500474 |
| 100,000 | 74920.8 | 0.0290546901 |

Tail estimate beyond gamma = 74920.8: 2.630e-06 +- 4.8e-07 (m = 0.0985 per unit height, spread 0.0181).

**kappa = 0.02905732 +- 4.8e-07** (partial 0.02905469 + tail).
