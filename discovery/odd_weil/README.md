# odd_weil: the odd semilocal Weil form, zero-free

Thread started 2026-10-04. Nothing here is proved; every result is numerical evidence.

## What this thread computes

The Weil quadratic form restricted to real odd test functions supported in [-log lambda, log lambda], built from the zero-free side of the explicit formula only (Archimedean Gamma term plus the prime or coefficient terms at log n). Odd tests kill the two pole channels, so any negative eigenvalue is a genuine zero signal. Basis: sin(k pi u / L).

## Status by note (read in this order)

| note | status |
| --- | --- |
| 2026-10-04_odd_sector_laplacian_scan.md | stands. 50-digit scan; zeta-cycle zero detection; positive-network split is exact balance, not slack |
| 2026-10-04_davenport_heilbronn_control.md | stands. Same functional equation, no Euler product: the form goes negative exactly at DH's off-line zero (0.8085 + 85.6993i); composite "fake prime" terms drive it |
| 2026-10-04_positivity_selects_the_primes.md | **retracted**: rigidity was a corner artifact |
| 2026-10-04_prime_blind_dual.md | **retracted** conclusion (the duality statement itself is correct) |
| 2026-10-04_correction_lambda_is_a_boundary_point.md | stands. Feasible set is large; Lambda sits on its boundary at zero margin; along b = c Lambda only c = 1 works for lambda >= 3 |
| 2026-10-04_euler_manifold_global_check.md | stands, multi-start evidence. Real local roots (self-dual degree-one Euler data): only alpha_p = 1 is feasible for lambda = 4..6. Complex roots reopen an interior at lambda = 4 |

## Scripts

- Builders: odd_weil.py (double), mp_weil.py + run_mp.py (50-digit), dh_control.py (zeta, L(chi_5), Davenport-Heilbronn; u-space Archimedean kernel).
- Validation: validate.py (zero side vs zero-free side), dh_validate.py, dh_zeros.py.
- Experiments: analyze.py, zerr.py, dh_scan.py, dh_explain.py, dh_breakdown.py, dh_arch_feasible.py, rigidity.py, rigidity_joint.py, maxreach.py, locking.py, prime_blind.py, prime_strength.py, euler_global.py, euler_global_complex.py, euler_complex_hard.py.
- Data: zeros.pkl (1500 zeta zeros), dh_zeros.pkl, mp_N14.pkl (50-digit matrices), dh_scan.pkl, prime_blind.pkl.

Requires numpy, scipy, mpmath, cvxpy + clarabel.

## Related

Codex's parallel experiment: GPPDiscovery2 scripts/check_euler_local_boundary.py (real local roots, local perturbations); euler_global.py here is the global multi-start check of the same claim.
