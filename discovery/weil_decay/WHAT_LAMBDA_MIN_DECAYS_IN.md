# What lambda_min(c) actually decays in

**2026-09-13, same session as the resolution criteria. Read
`ARCHIMEDEAN_TRUNCATION_CRITERION.md` and `CRITERION_RUNS.md` first: everything below uses
points that satisfy both criteria, and several earlier conclusions in this repo do not.**

**Nothing here is proved. This is numerics on a Galerkin truncation.**

## The headline

The README asks: *is `log lambda_min` linear in `log c`, and if so what is the constant?*

The answer is no, and the frame is wrong. `lambda_min` decays **exponentially in `c`**, and
the rate is **insensitive to the arithmetic content**.

```
lambda_min(c)  ~  A * exp(-K c),     K ~ 11.4
```

## The four unit steps

All at `N=36`, margin ~3.0, headroom >= 4.8, so both resolution criteria hold at every point.

| step | d ln(lambda) | new prime powers |
|------|--------------|------------------|
| 3 -> 4 | **-11.093** | `4 = 2^2` |
| 4 -> 5 | **-11.266** | `5` |
| 5 -> 6 | **-12.042** | **none** |
| 6 -> 7 | **-11.159** | `7` |

Every unit step in `c` costs the same, and the step that adds **no arithmetic at all** costs
the most. (The `6 -> 7` step compares against the `c=7` point at margin 2.58 rather than 3.0;
the archimedean difference is 0.013 in `ln`, negligible against 11.16.)

## The method that made this visible: non-prime cutoffs

`c` is a **cutoff**, not required to be prime. This repo has only ever scanned
`c = 3, 5, 7, 11, 13, 17, 19`. On that set the candidate explanatory variables are nearly
collinear -- `c`, `psi(c)`, `pi(c)` and `log c` all move together -- so no amount of
in-sample fit quality can separate them.

| c | p^k <= c | psi(c) | log c |
|---|----------|--------|-------|
| 3 | 2, 3 | 1.792 | 1.099 |
| 4 | 2, 3, 4 | 2.485 | 1.386 |
| 5 | 2, 3, 4, 5 | 4.094 | 1.609 |
| **6** | **2, 3, 4, 5** | **4.094** | **1.792** |
| 7 | 2, 3, 4, 5, 7 | 6.040 | 1.946 |

`c=6` carries exactly the prime content of `c=5`. It is the single point that separates
arithmetic from geometry, and it had never been run.

## Out-of-sample model selection

Fit each two-parameter law on `c=3` and `c=5` **only** (which determines it exactly), then
predict a point never used in the fit.

| model | predict c=4 | predict c=6 | predict c=7 | predict c=11 |
|-------|-------------|-------------|-------------|--------------|
| **exp(-K c)** | **1.09x** | **2.1x** | **2.0x** | **11x** |
| power law `c^-s` | 4.5x | 52x | 4.2e3 | 1.5e13 |
| exp(-K psi) | 78x | 1.6e5 | 63x | 170x |
| exp(-K sqrt c) | - | - | 132x | 2.1e7 |
| exp(-K c log c) | - | - | 25x | 2.6e9 |
| exp(-K c^2) | - | - | 3.8e4 | 2.2e30 |

`exp(-K c)` predicts `c=4` to within a factor of **1.09**, out of sample, on a quantity that
has fallen 20 orders of magnitude by `c=7`.

## Why the power law looked plausible for so long

A pure exponential has **log-log slope `-K c`**, which grows with `c`. That is exactly the
"curvature" that made the power-law reading look broken. Using a single `K` fitted on
`c = 3,5,7,11` to predict the measured log-log slopes:

| pair | c_eff | predicted `-K c_eff` | measured | error | headroom |
|------|-------|----------------------|----------|-------|----------|
| 3 -> 5 | 3.92 | -42.80 | -43.93 | +2.6% | 5.36 |
| 5 -> 7 | 5.94 | -64.99 | -68.72 | +5.7% | 4.38 |
| 7 -> 11 | 8.85 | -96.75 | -92.68 | -4.2% | 1.81 |
| 11 -> 13 | 11.97 | -130.89 | -82.24 | -37% | 1.69 |
| 13 -> 17 | 14.91 | -163.02 | -70.71 | -57% | 0.78 |
| 17 -> 19 | 17.98 | -196.59 | -63.88 | -68% | 0.75 |

Every well-resolved pair predicted to within 6% by one constant. The failures begin exactly
where the headroom does and grow monotonically with how bad it gets.

## The consequence for the "-64 ~ -4 gamma_1" coincidence

The README records a two-point natural-log slope of about `-64`, set against
`-2 gamma_1 = -28.27` and `-4 gamma_1 = -56.55`, and reads the proximity as a sign that the
quantity "sees the zeros."

**That slope is `-K c` evaluated at `c ~ 6`.** It is not a constant. Measure at `c = 3,5`
instead and you get `-44`, which is within 4% of `-3 gamma_1`. Measure at `c = 7,11` and you
get `-93`, near `-6.6 gamma_1`. The apparent `gamma_1` multiple is whatever the measurement
location makes it.

## The shape of the spectrum: a ladder, and what counts its rungs

Dumping the full low spectrum rather than two eigenvalues (`spectrum.py`):

| c | even sector, log10 | rungs below the bulk |
|---|--------------------|----------------------|
| 3 | -7.29, -2.86, **-0.17**, 0.04, 0.08 | 2 |
| 5 | -16.96, -11.24, -6.27, -2.60, **-0.22**, -0.07 | 4 |
| 6 | -22.14, -15.88, -10.63, -6.12, -2.40, **-0.30**, -0.11 | 5 |
| 7 | -27.01, -20.59, -14.93, -10.11, -5.74, -2.42, **-0.28**, -0.13 | 6 |

The operator is not a continuum being scaled. It is an **O(1) bulk plus a ladder of
exponentially small eigenvalues**, with a clean 2.1-2.7 decade gap between the lowest rung
and the bulk edge, in both parity sectors. The counts are identical in even and odd.

| c | rungs (even / odd) | prime powers <= c | c - 1 |
|---|--------------------|-------------------|-------|
| 3 | 2 / 2 | 2 | **2** |
| 5 | 4 / 4 | 4 | **4** |
| **6** | **5 / 5** | 4 | **5** |
| 7 | 6 / 6 | 5 | **6** |

**The count is `c - 1`, the number of integers in `[2, c]`, not the number of prime powers.**
`Lambda(6) = 0` and the operator produces a rung at `c=6` regardless.

`c = 3` and `c = 5` cannot distinguish the two hypotheses, because every integer in `[2,5]`
is a prime power; `6` is the first that is not. The first reading of this table -- "one rung
per prime power" -- was made on exactly those two points and is wrong.

**Still open:** all four points have integer `c`, so `floor(c) - 1` and a smooth geometric
quantity that happens to pass through those values are not yet distinguishable. `c = 6.5`
discriminates (an integer count gives 5, the same as `c=6`).

## Thread 2: what the archimedean place is doing

A point below the archimedean threshold computes a different, well-defined operator (see
`CRITERION_RUNS.md`). Its ground state, at matched margins ~0.85:

| c | pi(c) | sub-threshold lambda_1 |
|---|-------|------------------------|
| 3 | 2 | -0.4901288932404565 |
| 4 | 2 | -0.8033547640309103 |
| 5 | 3 | -1.0268339297106557 |
| 6 | 3 | -1.2615716118640772 |
| 7 | 4 | -1.5848179322855760 |

**Linear in `c`:** `lambda_1 = -0.26476 c + 0.29046`, max residual **0.037** over a range of
1.09. Linear in `pi(c)` gives max residual 0.165, four and a half times worse.

Both degenerate pairs are broken directly. `pi(3) = pi(4) = 2` yet the values differ by
**0.313**; `pi(5) = pi(6) = 3` yet they differ by **0.235**. If the prime count drove this,
those pairs would coincide. They do not.

Set that against the correct operator over the same range, where `ln lambda_1 = -11.4 c + b`.
Cutting the top archimedean modes does not perturb the ground state; it **changes the
functional form of its `c`-dependence from exponential to linear**. The archimedean place is
not a correction to the exponential smallness. It is the source of it.

`c=11` was run for this and is **discarded**: headroom 1.81, so it fails the prime-comb
criterion and is compromised by a defect unrelated to the one under study.

## Mechanism: stated as a hypothesis, not a result

The archimedean kinks sit at `2 pi x / log c`, so growing `c` compresses the whole archimedean
structure toward the origin. Combined with the sub-threshold result -- delete the archimedean
modes and the exponential becomes linear -- the `c`-dependence may be **wholly archimedean**,
entering through `log c` in the kink positions, with the primes setting where the rungs are
but not how far down they go. **Untested.**

A tempting wrong explanation, recorded so nobody re-derives it: *"nested test-function spaces,
so lambda_min is monotone decreasing by the variational principle."* That argument applies in
the `N -> infinity` continuum. At fixed `N` the Galerkin trial space is always `2N+1`
dimensional, and changing `c` gives a **different** subspace of that dimension, not a larger
one. Neither is nested in the other. The monotonicity is observed, not explained.

## What this does to the thread's premise

The thread exists because `lambda_min(c)` was hoped to be "a quantity that sees the zeros."
Three independent observables now say it mostly sees the box:

1. The decay is exponential in the cutoff, at a rate insensitive to prime content -- the step
   that adds no arithmetic costs the most of the four.
2. The archimedean-truncated operator is linear in `c` and blind to `pi(c)` at both
   degenerate pairs.
3. The ladder rung count tracks integers in `[2,c]`, not prime powers.

That does not make the quantity uninteresting, and it does not touch Weil positivity itself.
It does mean the **decay constant** is the wrong thing to be measuring for arithmetic content,
and that a scan costing 20-40 minutes a point has been aimed at it.

## Reproduce

```
python3 spectra.py     # variable comparison on the published scan
python3 modelcmp.py    # out-of-sample model selection
python3 shape.py       # spectrum ratios and separate decay rates
python3 point.py    --c 6 --N 36 --T 379 --dps 90   # the decisive composite cutoff
python3 spectrum.py --c 7 --N 28 --T 420 --dps 90 --keep 10
```


---

# The mechanism: levels are born at the bulk and descend

`c = 6.5`, run at matched margin 4.48, settles the rung count and in doing so gives the
whole picture. **`floor(c) - 1` is also wrong.**

The even sector shows 6 rungs and the odd 5, against `floor(6.5) - 1 = 5`. The disagreement
is the result: the sixth even level sits at `-1.12`, only 0.98 decades above the bulk, where
at every integer `c` the lowest rung stands 2.1-2.7 decades clear. It is a level **in
transit**, and a fixed threshold catches it in one sector and not the other.

Rungs do not appear at integers. They emerge continuously.

| c | lvl1 | lvl2 | lvl3 | lvl4 | lvl5 | lvl6 | bulk |
|---|------|------|------|------|------|------|------|
| 6.0 | -22.14 | -15.88 | -10.63 | -6.12 | -2.40 | -- | -0.30 |
| 6.5 | -24.65 | -18.35 | -12.83 | -8.08 | -4.05 | **-1.12** | -0.14 |
| 7.0 | -27.01 | -20.59 | -14.93 | -10.11 | -5.74 | **-2.42** | -0.28 |

A level detaches from the bulk between `c=6` and `c=6.5`, sits at `-1.12` halfway, and by
`c=7` has reached `-2.42` -- exactly where `c=6`'s newest rung was. **One level born per unit
`c`.** Meanwhile every existing level descends:

| level | descent per unit c, in ln |
|-------|--------------------------|
| 1 (the ground state) | **-11.20** |
| 2 | -10.84 |
| 3 | -9.91 |
| 4 | -9.19 |
| 5 | -7.67 |

**This closes the circle.** `K = 11.20` measured as the descent rate of the bottom level is
the same `K = 11.4` measured from `lambda_min(c)` across five cutoffs. `lambda_min` is simply
the oldest rung -- born around `c ~ 2` and falling at a constant rate since -- so

```
ln lambda_min(c)  ~  -K (c - 2) + const.
```

The exponential law, the rung count growing by one per unit `c`, and the `floor(c)-1`
coincidence at integer `c` are one phenomenon seen three ways. None of it is arithmetic:
levels are born at a rate of one per unit cutoff and descend at a fixed rate, whether or not
a prime power is crossed.

## Sampling bias, recorded once rather than six times

Six claims were made and retracted in this session:

1. the sub-threshold magnitude tracks the archimedean margin -- refuted by a third point at
   the same margin;
2. the varying archimedean error tilts the slope -- it is 0.08% of the signal;
3. the rigid `lambda_1`/gap ratio is structure -- it is an under-resolution artifact;
4. the decay is a product over prime powers, `exp(-K psi)` -- refuted by `c=6`;
5. one ladder rung per prime power -- refuted by `c=7`;
6. `floor(c) - 1` rungs -- refuted by `c=6.5`.

These are not six independent mistakes. **Every one chose the more arithmetic of two
explanations that the sample could not separate**, and in each case the sample could not
separate them because it was drawn from a set on which the candidates are collinear: primes
only (1, 4, 5), integers only (6), or points already disqualified for another reason (3).

The operational rule that would have caught all six:

> Before claiming variable `X` explains `Y`, list what else is collinear with `X` **on the
> points you actually have**, and find the cheapest point that breaks the collinearity. If
> no such point exists yet, the claim is not available.

Every one of the six was killed by a single run costing 10-40 minutes.
