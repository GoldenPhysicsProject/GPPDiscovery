# weil-decay

Search-and-discovery repo for the ground-state decay of the **truncated Weil
quadratic form**.

## The question

Connes–van Suijlekom (Prop 4.1) and Connes–Consani–Moscovici (Lemma 5.1) build,
for a prime cutoff `c` and band `N`, a finite `(2N+1) x (2N+1)` Galerkin matrix
`Q(c)`. CvS Theorem 6.1: the zeros of its characteristic function lie on the
critical line **for every finite c**. Criticality is a theorem. Connes (2026,
§6) poses **convergence** as the open question.

`lambda_min(c)` is the smallest eigenvalue of `Q(c)` on the even sector.
Weil positivity, hence RH, is `lambda_min >= 0` for every `c`.

We are not testing whether it is positive. We are measuring **how fast it
decays**:

> Is `log lambda_min` linear in `log c`, and if so, what is the constant?

Measured so far (unconverged, N=20, T=100): `3.67e-27` at c=7, `8.09e-18` at
c=5. Nine orders in one step. A two-point natural-log slope of about `-64`,
against `-2*gamma_1 = -28.27` and `-4*gamma_1 = -56.55`.

If the constant is arithmetic, this is a quantity that *sees the zeros*, which
is more than can be said for most reformulations of RH.

## Why N-convergence comes first

Going from dim 33 to dim 41 at c=7 moved `lambda_min` from `2.4e-25` to
`3.67e-27`, a factor of 65. So `N=20` is **not converged** and any slope fitted
to it is meaningless. The scan does `c=7` at `N = 20, 28, 36` before anything
else. If it is still moving at 36, raise N before trusting the c-scan.

## Running

`python point.py --c 7 --N 36 --T 300 --dps 90` for one point, or push / use
workflow_dispatch to run the whole matrix in parallel on Actions. Results land
in `results.jsonl` and `RESULTS.md`, committed back automatically.

## Caveat on small c

At `c=3, N=20, T=100` we saw `lambda_min = -0.332`. Almost certainly marginal
basis resolution (41 modes on an interval of length 2.2, archimedean quadrature
truncated at T=100); the reference implementation reports the same at c=23, 29
and attributes it to exactly that. The scan re-runs c=3 at N=36, T=300 to
settle it. If it survives, that is a much bigger deal than the slope.

## Status

Discovery only. Nothing here is proved. Anything that becomes a theorem goes to
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify) for Lean 4
formalization, and the first thing that should graduate is **certified
enclosures** for `lambda_min > 0` via interval arithmetic, not another
equivalence.
