# Smallest normalized zero gaps over all of Platt's tables (NUMERICAL)

Data: all 17 files, 63,150,760 zeros, heights 14 to 27,746,000 (see `inventory.md`; zero #1 and one zero per file agree with `mpmath.zetazero` to 25+ digits, all files contiguous). Normalized gap = (gamma_(n+1) - gamma_n) log(gamma_n / 2 pi) / (2 pi). Gaps are built from Platt's integer deltas (units 2^-101) and agree with differences of the exact reader heights to 3e-12 relative over the 100 pairs below; float64 heights at t ~ 1e7 could not resolve them.

The 1000 smallest range from 0.001235 to 0.025327 (full list: `smallest_gaps_1000.csv`).

| rank | N (lower zero) | gamma | normalized gap | abs zeta' (lower) | abs zeta' (upper) | pair contribution to kappa |
|---|---|---|---|---|---|---|
| 1 | 41,820,581 | 18882503.9013 | 0.001235 | 0.0190 | 0.0190 | 3.10e-11 |
| 2 | 24,360,732 | 11412127.1122 | 0.002021 | 0.0143 | 0.0143 | 1.50e-10 |
| 3 | 52,823,900 | 23482717.9702 | 0.004212 | 0.0713 | 0.0714 | 1.42e-12 |
| 4 | 42,513,405 | 19174199.8839 | 0.004531 | 0.0532 | 0.0531 | 3.85e-12 |
| 5 | 49,079,247 | 21924545.3943 | 0.004887 | 0.0412 | 0.0411 | 4.91e-12 |
| 6 | 8,546,950 | 4316832.4671 | 0.004970 | 0.0279 | 0.0278 | 2.77e-10 |
| 7 | 39,436,579 | 17876433.1678 | 0.004998 | 0.0765 | 0.0765 | 2.14e-12 |
| 8 | 1,115,578 | 663318.5083 | 0.005447 | 0.1751 | 0.1756 | 2.96e-10 |
| 9 | 33,349,976 | 15289808.6664 | 0.005496 | 0.1299 | 0.1295 | 1.02e-12 |
| 10 | 61,515,842 | 27074091.1555 | 0.005857 | 0.1092 | 0.1092 | 4.58e-13 |

Over the 100 smallest-gap pairs: |zeta'(rho)| ranges 0.0143 to 0.477 (small gaps go with small |zeta'|, as expected). Their total contribution to kappa = 2 sum 1/(gamma^2 |zeta'|^2) is 4.97e-09, a share of **1.71e-07** of kappa = 0.029057 (the zero-sum value from `kappa.md`; the share column in the CSV uses it). It is negligible because 1/gamma^2 at gamma ~ 1e7 is about 1e-14 even when |zeta'| is 0.01.
