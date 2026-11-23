# Explicit-formula reconstruction of B(u) (NUMERICAL; assumes RH and simple zeros for the form of the series)

B(u) = S(floor(e^u)) - c, c = 1/zeta(1/2) = -0.6847652361. B_T(u) = sum_(0<gamma<=T) 2 Re[e^(i gamma u)/(i gamma zeta'(rho))] (the pole of x^s/(s zeta(s+1/2)) at s = i gamma). u-grid: 20001 equispaced points on [10, 18] (N from 22,026 to 65,659,969). rms of the sieved B on the grid: 0.16978. The sieved B is a step function with jumps mu(k)/sqrt(k) at integers k, so a floor on the mismatch is expected from the frequencies beyond T.

| zeros | T = gamma_max | L2 mismatch (rms) | relative to rms B | correlation |
|---|---|---|---|---|
| 10 | 49.8 | 0.06114 | 0.3601 | 0.93290 |
| 30 | 101.3 | 0.04312 | 0.2540 | 0.96720 |
| 100 | 236.5 | 0.02707 | 0.1595 | 0.98720 |
| 300 | 541.8 | 0.01817 | 0.1070 | 0.99426 |
| 1,000 | 1419.4 | 0.01126 | 0.0663 | 0.99780 |
| 3,000 | 3533.3 | 0.00745 | 0.0439 | 0.99904 |
| 10,000 | 9877.8 | 0.00444 | 0.0261 | 0.99966 |
| 30,000 | 25755.3 | 0.00272 | 0.0161 | 0.99987 |
| 100,000 | 74920.8 | 0.00161 | 0.0095 | 0.99995 |
