# CCM trial-vector overlap (Poisson-prolate k_lambda)

Evidence only. Each row is one (lambda, N). eta-normalised (eta^T k = 1).

| lambda | N | dps | eps0 | first-zero err | kTk | M_C | kTk*M_C | eta.k raw | max overlap r>=1 | even-mode overlaps (>1e-12) |
|---|---|---|---|---|---|---|---|---|---|---|

## Low eigenvalues and clusters


## Convergence in N (first even overlap, kTk)


## Scale-free reading (eta-normalisation is ill-posed: eta.k_raw -> 0 as N grows)

kTk_sf = kTk*(eta.k_raw)^2/max|c|^2 ; cos_r <= |ov_r|*|eta.k_raw|/max|c| (upper bound on the cosine with mode r).

| lambda | N | kTk_sf | M_C | kTk_sf*M_C | cos(e2) bound | cos(e4) bound | lambda_even1 |
|---|---|---|---|---|---|---|---|

Verdict rule: GO if kTk*M_C stays bounded while the first even overlap decays no faster
than the eigenvalue (alpha <= 1/2); NO-GO if kTk*M_C grows with N; else INCONCLUSIVE.
A verdict only counts once the N-sequence has converged.
