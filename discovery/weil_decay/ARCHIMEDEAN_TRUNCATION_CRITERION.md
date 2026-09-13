# The archimedean truncation criterion

**2026-09-13. Numerical-validity finding. Nothing here is proved about the operator —
this is a statement about when `connes_cvs.build_galerkin_matrix` computes the matrix
it is supposed to compute.**

## The short version

`RESULTS.md` and this repo's README have carried, since the beginning, one point that did
not fit: at `c=3, N=20, T=100` the smallest even-sector eigenvalue came out **negative**,
`lambda_1 = -0.332`, and was set aside as "almost certainly marginal basis resolution,"
with no criterion attached. There is a criterion, it is one line, and it separates every
point on record.

`connes_cvs` computes the archimedean piece of each matrix entry as

```
psi_arch(x) = (1/2pi^2) * INT_{-T}^{T} h_plus(tau) * Re S_hat(tau, x) d_tau
```

splitting the quadrature at the integrand's kinks `tau = 0` and `tau = +/- alpha_x`, where

```
alpha_x = 2*pi*x / log(c)
```

and `x` runs over the basis indices `-N .. N`. (`connes_cvs/operator.py`, `psi_arch`: the
split list is literally `{0, alpha_x, -alpha_x}` filtered to `-T < s < T`.)

So the feature carrying basis mode `x`'s archimedean weight sits at spectral height
`alpha_x`, and the truncation `|tau| <= T` keeps it only while `alpha_x <= T`. For the
outermost mode `x = N`:

> **The Galerkin matrix is the intended matrix only while**
>
> ```
> margin := T * log(c) / (2*pi*N)  >=  1
> ```

Below 1, the top basis modes are built with their archimedean contribution truncated away.
The matrix is not the operator's matrix, and `lambda_1` goes negative.

## Every point on record, sorted by margin

| c | N | T | tau_max | margin | lambda_1 |
|---|---|---|---------|--------|----------|
| 3 | 20 | 60 | 114.4 | **0.52** | `-0.4901` |
| 3 | 16 | 60 | 91.5 | **0.66** | `-0.4855` |
| 3 | 36 | 120 | 205.9 | **0.58** | `-0.4901` |
| 3 | 36 | 160 | 205.9 | **0.78** | `-0.4897` |
| 3 | 12 | 60 | 68.6 | **0.87** | `-0.2444` |
| 3 | 20 | 100 | 114.4 | **0.87** | `-0.332` (the README's anomaly) |
| 3 | 36 | 300 | 205.9 | 1.46 | `+4.951e-8` |
| 5 | 36 | 300 | 140.5 | 2.13 | `+8.904e-18` |
| 7 | 36 | 300 | 116.2 | 2.58 | `+8.093e-28` |
| 19 | 52 | 300 | 111.0 | 2.70 | `+2.949e-75` |
| 11 | 36 | 300 | 94.3 | 3.18 | `+5.191e-46` |
| 19 | 44 | 300 | 93.9 | 3.20 | `+1.038e-69` |
| 13 | 36 | 300 | 88.2 | 3.40 | `+5.602e-52` |
| 17 | 36 | 300 | 79.8 | 3.76 | `+3.242e-60` |
| 19 | 36 | 300 | 76.8 | 3.91 | `+2.661e-63` |

Sorted by margin the table separates completely: **every point below 1 is negative, every
point at or above 1.46 is positive.** No exceptions.

## The controlled test

The table above mixes `c`, `N` and `T`, so on its own it only shows a correlation. The
decisive run holds `c = 3` and `N = 36` fixed and varies **only `T`**:

| c | N | T | margin | lambda_1 |
|---|---|---|--------|----------|
| 3 | 36 | 120 | 0.5825 | `-0.4901288932404565165938318` |
| 3 | 36 | 160 | 0.7767 | `-0.4897489413922827120029554` |
| 3 | 36 | 190 | **0.9228** | `-0.3779171010840920867089124` |
| 3 | 36 | 210 | **1.0200** | `+4.656080582364100789220654e-8` |
| 3 | 36 | 300 | 1.4571 | `+4.950881954831601328856756e-8` |

Same matrix size, same precision, same `c`. The sign of the ground eigenvalue is controlled
by `T` alone, and it flips across `margin = 1`.

The `T = 190` / `T = 210` pair is the sharpest form of the test this apparatus allows: a
**10% change in a single parameter**, straddling `margin = 1`, flips the sign and moves
`lambda_1` by seven orders of magnitude.

It is also not a continuous drift with `T`. Above the line the value has already settled —
`T=210` gives `+4.656e-8` against `T=300`'s `+4.951e-8`, a 6% difference, while `T=190` sits
seven decades away on the other side. Below the line `lambda_1` is pinned near `-0.49`, lifts
toward zero as the margin approaches 1, crosses, and then stops moving. That is a threshold,
not a trend.

## Magnitude tracks the margin too, not just the sign

Margins `0.87, 0.66, 0.52` at `c=3, T=60` give `lambda_1 = -0.2444, -0.4855, -0.4901`.
More telling: two points with *different* `(N, T)` but the *same* margin `0.87` —
`(N=20, T=100)` and `(N=12, T=60)` — give `-0.332` and `-0.2444`. The margin, not `N` or
`T` separately, looks like the controlling variable.

## Why this matters beyond bookkeeping

**1. `T` is not a free parameter.** It is tied to `N` and `c`. Holding `T` fixed while
raising `N` walks off the cliff, and it does so **sooner for smaller `c`**, because
`log(c)` is in the denominator. That is the entire reason `c=3` has always been the value
that misbehaves. It was never `c=3`; it was the margin.

**2. The c-scan does not have uniform numerical quality.** At `N=36, T=300` the margin runs
from **1.46** at `c=3` to **3.91** at `c=19` — the small-`c` end sits nearly three times
closer to the cliff. A decay constant fitted across that range is fitting a varying
truncation quality along with the mathematics.

This is a live candidate explanation for the curvature. Consecutive natural-log slopes of
`lambda_1` across the `N=36` c-scan are

```
c= 3 ->  5 :  -43.93
c= 5 ->  7 :  -68.72
c= 7 -> 11 :  -92.68
c=11 -> 13 :  -82.24
c=13 -> 17 :  -70.71
c=17 -> 19 :  -63.88
```

— they rise then fall, and a linear fit leaves **3.55 decades** of residual in `log10`.
There is no single decay constant in this data. The slope is steepest exactly where the
margin is changing fastest.

**The fix for scan design: hold the margin fixed, not `T`.** Set
`T = kappa * 2*pi*N / log(c)` for a fixed `kappa` (>= 2), so every point in a c-scan sits
the same distance from the cliff. `resolution.plan_T(c, N, kappa)` does this. At `N=36`,
`kappa=3` that means `T = 618` at `c=3` down to `T = 231` at `c=19` — i.e. the current
fixed `T=300` is simultaneously too small at the bottom of the scan and wasteful at the top.

## What is now recorded automatically

`point.py` writes `tau_max`, `margin` and `within_window` into every result record. A point
with `margin < 1` can no longer be mistaken for a measurement of anything.

## Two open questions this does NOT answer

**(a) `lambda_1` is not converging in `N` at `c=19`, and that is inside the valid region.**
All three points (`N=36, 44, 52`, margins 3.91, 3.20, 2.70) are well above 1, and yet

```
N=36 -> 44 : d log10(lambda_1) = -6.409     d log10(gap) = -6.050
N=44 -> 52 : d log10(lambda_1) = -5.547     d log10(gap) = -5.325
```

Both the ground eigenvalue and the gap are still falling geometrically, with decrement
ratios 0.866 and 0.880. If that continues, "`lambda_min(c) > 0` at finite `c`" is a
statement about the truncation rather than about the operator. Extrapolating the two
decrements geometrically puts the limit near `1e-110`, but two decrements is not enough to
call it; it is enough to say the quantity is **not** converged and no slope fitted to it
means anything yet.

*(Answered below: this is the prime-spacing criterion, not the archimedean one. `c=19` at
`N=36, 44` does not resolve the prime comb.)*

**(b) `lambda_1` and the gap are not independent.** `log10(lambda_1) - log10(gap)` across
the `N=36` c-scan runs `-4.44, -5.74, -6.36, -6.74, -6.79, -6.71, -6.67` — it **plateaus**
near `-6.7` for `c >= 11`. The Davis-Kahan strategy needs `r_c / delta_c`, and `delta_c` is
not free to be large while `lambda_1` is small; the two track each other with a
`c`-independent ratio at the large-`c` end. Worth understanding before that route is costed.

## Reproduce

```
python3 resolution.py                                   # the criterion and the table
python3 point.py --c 3 --N 36 --T 120 --dps 90          # margin 0.58 -> negative
python3 point.py --c 3 --N 36 --T 300 --dps 90          # margin 1.46 -> positive
```

Environment for the 2026-09-13 runs: Python 3.11.15, mpmath 1.4.1, connes-cvs 0.2.2.
(The earlier records in `results.jsonl` were produced with mpmath 1.3.0; the archimedean
criterion is a property of the quadrature layout, not of the mpmath version, and the
`c=3, N=36, T=300` point reproduces.)

---

# The prime-spacing criterion (a second, independent one)

Added 2026-09-13, same session. The criterion above is about `T`. This one is about `N`
alone, and it is the larger effect.

## Statement

The Weil functional's prime piece is supported at `t = +/- log(p^k)` for prime powers
`p^k <= c`. The Galerkin basis is Fourier modes `exp(2*pi*i*k*t/(2 log c))` for `|k| <= N`:
`2N+1` degrees of freedom on an interval of length `2 log c`, so its resolution in `t` is

```
delta_t = 2 log(c) / (2N + 1).
```

Two prime powers closer together than `delta_t` are not separated by the basis. Hence

> ```
> N  >=  log(c) / min_gap(c)  -  1/2,    min_gap(c) = min adjacent gap of {log p^k : p^k <= c}
> ```

`spacing.py` computes it.

| c | #p^k | tightest pair | min_gap | N_req | N=36 | N=52 |
|---|------|---------------|---------|-------|------|------|
| 3 | 2 | 2/3 | 0.4055 | 2.2 | ok | ok |
| 5 | 4 | 4/5 | 0.2231 | 6.7 | ok | ok |
| 7 | 5 | 4/5 | 0.2231 | 8.2 | ok | ok |
| 11 | 8 | 8/9 | 0.1178 | 19.9 | ok | ok |
| 13 | 9 | 8/9 | 0.1178 | 21.3 | ok | ok |
| 17 | 11 | **16/17** | 0.0606 | 46.2 | **NO** | ok |
| 19 | 12 | **16/17** | 0.0606 | 48.1 | **NO** | ok |
| 23 | 13 | 16/17 | 0.0606 | 51.2 | NO | ok |
| 29 | 16 | 16/17 | 0.0606 | 55.0 | NO | **NO** |
| 31 | 17 | 16/17 | 0.0606 | 56.1 | NO | **NO** |
| 37 | 19 | **31/32** | 0.0317 | 113.2 | NO | **NO** |

## The test

`N`-convergence of `lambda_1` at fixed `c`, `T=300`, going `N=36 -> 44`:

| c | N_req | headroom 36/N_req | change in lambda_1 |
|---|-------|-------------------|--------------------|
| 3 | 2.2 | 16.4 | **-0.23%** |
| 5 | 6.7 | 5.4 | **-3.24%** |
| 19 | 48.1 | **0.75** | **factor 3.9e-7** (6.4 decades) |

Converged, converged, catastrophically not converged — predicted in advance from prime-power
spacings alone, with no reference to the eigenvalue computation.

## Consequences

**1. The `N=36` c-scan mixes five converged points with two that are not.** `c = 3, 5, 7, 11,
13` are resolved; `c = 17, 19` are not, and they are the large-`c` end where any slope fit is
anchored. The seven-point scan in `RESULTS.md` should not be fitted as it stands.

**2. Even the resolved points do not give a constant slope.** Consecutive natural-log slopes
among `c = 3, 5, 7, 11, 13` are `-43.93, -68.72, -92.68, -82.24`. Still rising then falling.
So removing the unresolved points does not rescue the single-decay-constant picture; the
honest statement is that there are at most five usable values of `c` here and they do not
determine a constant.

**3. Headroom is not constant across the scan either** — `36/N_req` runs
`16.4, 5.4, 4.4, 1.81, 1.69` for `c = 3 .. 13`, and `c=5` already drifts 3.2% at headroom 5.4.
A c-scan should hold **headroom** fixed, exactly as it should hold the archimedean **margin**
fixed: vary `N` with `c` rather than pinning it, and set `T` from `N` and `c`.

**4. The binding constraint is arithmetic.** What forces `N` up is `16/17` — a power of two
landing next to a prime. It stays the binding pair from `c=17` all the way to `c=31`, and then
`31/32` takes over at `c=37` and `N_req` jumps from 56 to 113. The cost of this computation is
governed by how well powers of small primes approximate other primes. That is not a generic
quadrature nuisance; it is the arithmetic pushing back on the discretization, and it means the
scan gets abruptly harder at specific `c`, not smoothly harder.

## Scan design that follows

For a c-scan that measures the same thing at every point:

```
N(c) = ceil(eta * (log(c)/min_gap(c) - 1/2))      # fixed headroom eta, say 4
T(c) = ceil(kappa * 2*pi*N(c)/log(c))             # fixed archimedean margin kappa, say 3
```

Both helpers are in `spacing.py` and `resolution.py`. Note the two criteria pull in opposite
directions in `c`: `N_req` grows with `c` (more prime powers, tighter gaps) while the
archimedean threshold `T*log(c)/(2*pi*N)` gets *easier* with `c` (bigger `log c`). At fixed
headroom the required `T` is therefore roughly `kappa * 2*pi*eta/min_gap(c)` — set by the
prime-power spacing alone, independent of `log c`.

### The cost this implies

`eta = 4`, `kappa = 3`:

| c | N_req | N | T | dim |
|---|-------|---|---|-----|
| 3 | 2.2 | 9 | 155 | 19 |
| 5 | 6.7 | 27 | 317 | 55 |
| 7 | 8.2 | 33 | 320 | 67 |
| 11 | 19.9 | 80 | 629 | 161 |
| 13 | 21.3 | 86 | 633 | 173 |
| 17 | 46.2 | 185 | 1231 | 371 |
| 19 | 48.1 | 193 | 1236 | 387 |
| 23 | 51.2 | 205 | 1233 | 411 |
| 29 | 55.0 | 221 | 1238 | 443 |
| 31 | 56.1 | 225 | 1236 | 451 |

A uniform-quality scan needs `mp.eigsy` on 370-450 dimensional matrices at 90 dps, where the
present work runs at dim 73-121 and already costs 20-40 minutes a point. **The measurement the
thread actually wants is out of reach with this apparatus at `c >= 17`**, and the reason is
arithmetic: `16/17`. Lowering `eta` buys it back at the cost of accuracy, and `c=5` already
drifts 3.2% at `eta = 5.4`, so there is not much room. Worth costing a different eigenvalue
method (shift-invert on the even sector, or a sparse/structured factorization of `Q`) before
buying more `c`.
