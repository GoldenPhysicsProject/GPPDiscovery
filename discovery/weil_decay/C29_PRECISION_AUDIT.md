# c=29 ground-state precision audit

Date: 2026-08-24  
Implementation: `connes-cvs==0.2.2`, `mpmath==1.3.0`, `python-flint==0.9.0`

## Question

The older `E(c)` scan reported

```text
c=29, N=28, T=300, dps=60:
lambda_1 = -7.3152019431628632973e-61
```

Because the claimed eigenvalue is smaller than `10^-60`, its sign sits beyond
the requested 60-decimal working precision.  This audit tests whether the
negative sign survives higher precision, matrix diagnostics, larger basis
sizes, and different Archimedean cutoffs.

The executable reproduction is [`precision_audit.py`](precision_audit.py).
It uses the package's exact matrix formulas and parallelizes only the
independent `psi(k), psi'(k)` evaluations.

## Precision sweep at fixed c=29, N=28, T=300

| dps | lambda_1 | sign | comment |
|---:|---:|:---:|---|
| 60 | `-7.315201943162863297331931e-61` | - | reproduced old result |
| 90 | `+1.59354525025261901834003210284e-62` | + | sign flips |
| 120 | `+1.5935452502526190183400321027707511976e-62` | + | stable |
| 150 | `+1.5935452502526190183400321027707511976e-62` | + | stable |

At 60 dps there are `-0.136` guard digits beyond the eigenvalue scale.  The
eigenpair residual is `1.295e-61`, comparable to the alleged eigenvalue itself,
and Cholesky fails.  At 120 dps there are `58.2` guard digits; the diagnostics
are:

| diagnostic | value |
|---|---:|
| lambda_1 | `1.59354525025261901834003210277e-62` |
| lambda_2 | `3.40591348870452672919958541817e-56` |
| spectral gap | `3.40591189515927647658056707814e-56` |
| eigenpair residual (2-norm) | `1.6384612196721549e-121` |
| relative residual | `1.4520936457994941e-122` |
| Rayleigh quotient minus lambda_1 | `-1.7437099537018207e-122` |
| symmetry defect | `0` |
| parity defect | `0` |
| Cholesky | succeeds |
| Cholesky reconstruction defect | `6.05092486695206e-123` |

This is decisive for the old sign: the negative value is a precision artifact,
not a stable negative mode of that finite matrix.

## N sweep at adequate precision, c=29, T=300

| N | dps | lambda_1 | sign |
|---:|---:|---:|:---:|
| 20 | 120 | `1.1514939759527984e-50` | + |
| 24 | 120 | `6.0402935827598580e-57` | + |
| 28 | 120 | `1.5935452502526190e-62` | + |
| 32 | 140 | `6.2314625171601147e-68` | + |
| 36 | 160 | `3.5685685204312120e-73` | + |

Every adequately resolved point is positive, but `lambda_1` still falls by
roughly five to six orders of magnitude for each `N += 4`.  Therefore the
Galerkin sequence is **not N-converged**, and these data do not justify a
cutoff-decay fit.

## Archimedean T sweep, c=29, N=28

| T | dps | lambda_1 | sign |
|---:|---:|---:|:---:|
| 200 | 120 | `1.3923821502104579e-62` | + |
| 300 | 120 | `1.5935452502526190e-62` | + |
| 400 | 120 | `1.6421073469104179e-62` | + |
| 600 | 120 | `1.8521380308819276e-62` | + |
| 800 | 140 | `1.8732931281302787e-62` | + |

The sign is stable, but the magnitude is not fully T-converged.  The `T=600`
to `T=800` change is about 1.1%.

## Classification

1. **Closed:** the published `c=29, N=28, T=300, dps=60` negative sign is a
   working-precision artifact.
2. **Not closed:** finite-matrix positivity has not been certified by interval
   arithmetic.  High-precision `eigsy` plus Cholesky is strong numerical
   evidence, not a proof.
3. **Not closed:** neither the N limit nor the Archimedean T limit is converged.
4. **Invalid for now:** fitting a c-decay exponent from these smallest
   eigenvalues.

## Next rigorous step

Use outward-rounded interval arithmetic or a certified LDL/Cholesky enclosure
for each finite even-sector matrix.  Precision must be adaptive: at minimum,
request substantially more than `-log10(abs(lambda_1))` digits and require the
residual/enclosure radius to be much smaller than `lambda_1`.  Only after N and
T convergence are quantitatively controlled should the c-scaling question be
reopened.

