# Regulated scalar box: exact cut -> dispersion -> dilogarithm derivation

Starting from the audited celestial cut paper, write the Euclidean invariants as

- `s = -S`, `u = -U`, with `S,U>0`,
- `m = mu^2 > 0`.

The exact cut entering the fixed-`u` dispersion relation is

`C(s',-U,m) = (1/(8*pi))*4/sqrt(d(d+4c))*atanh(sqrt(d/(d+4c)))`,

with `d = U s'`, `c = m(s'+m)`, and

`J(-S,-U) = 8*pi * int_0^inf C(s',-U,m)/(s'+S) ds'`.

Define

`r^2 = U s' / ((U+4m)s' + 4m^2)`.

Then

`s' = 4 m^2 r^2 / (U-(U+4m)r^2)`,

and the square-root/Jacobian combination collapses exactly:

`ds'/sqrt(U s' ((U+4m)s'+4m^2)) = 2 dr/(U-(U+4m)r^2)`.

After combining with `(s'+S)^-1`, all extra factors cancel and one gets

`J(-S,-U) = 8/(S U) int_0^R atanh(r)/(1-kappa^2 r^2) dr`,

where

`R^2 = U/(U+4m)`,

`kappa^2 = [S(U+4m)-4m^2]/(S U)`.

This is already the analytic explanation for the polylogarithmic box structure: `atanh r` is logarithmic and the remaining denominator has only two simple linear factors after partial fractions.

Put

`y=(1+r)/(1-r)`, `alpha=(1-kappa)/(1+kappa)`.

For

`F_q(y)=log(y) log(1+y/q) + Li_2(-y/q)`,

one finds

`int_0^R atanh(r)/(1-kappa^2 r^2) dr`

`= (1/(4 kappa))*[F_alpha(y)-F_(alpha^-1)(y)-F_alpha(1)+F_(alpha^-1)(1)]`,

so

`J(-S,-U) = 2/(kappa S U) * [F_alpha(y)-F_(alpha^-1)(y)-F_alpha(1)+F_(alpha^-1)(1)]`.

The companion script `regulated_box_dilog.py` verifies this against both the reduced dispersion integral and an independent direct Feynman-parameter evaluation. At the paper anchor `(S,U,m)=(3,2,0.5)` all routes give

`1.043736598981870404239036240933802454240557498...`

to ~50 digits.

## Immediate next steps

1. Derive the `m -> 0+` asymptotics in this parametrization. Note the exact endpoint identity
   
   `1-kappa^2 R^2 = 4 m^2/[S(U+4m)]`,
   
   which controls the double-log IR singularity.
2. Rewrite the finite part in standard scalar-box variables/cross-ratios and reduce the dilog arguments by standard identities.
3. Replace scalar tree factors by helicity-dependent Yang-Mills/gravity cut trees and derive numerator factors before performing dispersion.
4. Seek a Mehler-Fock representation of the cut kernel so that the principal-series block and the dispersion kernel can be composed analytically.
