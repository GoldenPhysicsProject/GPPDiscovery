# D-dimensional unitarity boundary for the celestial box program

## Correction

The mass parameter in the current internally regulated scalar-box paper is **not** automatically the same object as the `mu^2` appearing in D-dimensional generalized unitarity.

Current scalar-box construction:

- the two `s`-channel cut propagators are exactly massless;
- the two uncut propagators carry mass `m = mu^2`;
- therefore the celestial cut geometry is the massless two-particle sphere with antipodal directions.

By contrast, in `D=4-2 epsilon` unitarity one splits

\[
\ell_D=\ell_4+\ell_\perp,\qquad \mu^2=-\ell_\perp^2\ge0.
\]

A D-dimensional **massless** cut condition becomes

\[
\ell_D^2=0\quad\Longrightarrow\quad \ell_4^2=\mu^2.
\]

Thus the four-dimensional cut states are effectively massive. Rational all-plus gauge/gravity amplitudes are invisible to strictly four-dimensional cuts but are recovered by D-dimensional cuts; this is standard in the unitarity literature (e.g. Bern-Dixon-Perelstein-Rozowsky, hep-th/9809160; Glover-Williams, arXiv:0810.2964).

## Consequence for the program

The honest route now bifurcates:

1. **4D cut-constructible sector:** continue using the proven massless celestial two-particle cut, principal-series phase-space factor, conformal blocks and dispersion reconstruction.
2. **Rational / D-dimensional sector:** derive the equal-mass two-particle cut with four-dimensional mass `mu`, then identify its correct celestial harmonic basis. This may require the massive celestial basis rather than the massless `S^2` principal-series basis.

Do not claim that the present internally regulated scalar box already captures D-dimensional rational terms. It does not: its cut legs are massless by construction.

## Immediate technical target

Derive the exact equal-mass two-particle phase-space measure first. In the center-of-mass frame `P=(M,0)`, with cut mass `mu` and `M>2mu`, each cut momentum has

\[
E=M/2,\qquad |\mathbf p|=\frac{M}{2}\sqrt{1-\frac{4\mu^2}{M^2}},
\]

and the integrated two-body phase space is

\[
\int d\Pi_2=\frac{1}{8\pi}\sqrt{1-\frac{4\mu^2}{M^2}}.
\]

The angular direction is still an `S^2`, but energy is no longer encoded by a null celestial scaling `omega q(z)`. The next problem is therefore to determine the Mellin/harmonic transform appropriate to this massive shell and its `mu -> 0` reduction to the already-proved Gamma/Beta phase-space weight.
