# Exact bridge: cut weight, Wiener-Hopf Fourier window, conical prefactor

Let

`s = 1/2 + i t`, `lambda = 2 t`, `k = 2 pi t`.

The principal-series cut weight is

`P(2t)=|Gamma(1+2it)|^2=|Gamma(2s)|^2`.

The Fourier transform of `P(lambda)=pi lambda/sinh(pi lambda)` is

`P_hat(k)=pi/(2 cosh^2(k/2))`.

Its Wiener-Hopf Gamma factorization is

`P_hat_+(k)=(1/sqrt(2pi))*Gamma(1/2 - i k/(2pi))^2`,

`P_hat_-(k)=(1/sqrt(2pi))*Gamma(1/2 + i k/(2pi))^2`.

Hence on `k=2pi t`,

`P_hat(2pi t)=|Gamma(s)|^4/(2pi)`.

The conical-reduction prefactor of the scalar principal-series block satisfies

`|c(lambda)|^2 = 4 |Gamma(1+i lambda)|^2 / |Gamma(1/2+i lambda/2)|^4`.

At `lambda=2t`,

`|c(2t)|^2 = 4 P(2t)/|Gamma(s)|^4`.

Combining with the Fourier formula gives the exact identity

`P(2t) = (pi/2) |c(2t)|^2 P_hat(2pi t)`.

Using the already-proved prefactor modulus

`|c(2t)|^2 = (4t/pi)coth(pi t)`,

this is equivalently

`P(2t) = 2 t coth(pi t) P_hat(2pi t)`.

Direct hyperbolic reduction:

`2 t coth(pi t) * [pi/(2 cosh^2(pi t))] = 2 pi t/sinh(2 pi t) = P(2t)`.

## Interpretation boundary

This is an exact special-function identity. It shows that the phase-space Plancherel weight, its logistic/Fourier thermal window, and the principal-series block normalization are algebraically coupled. It does not by itself establish a Weil-criterion or RH statement.

## Why it matters

The earlier Wiener-Hopf program sought a manifest-square Archimedean functional from `P_hat_+ P_hat_-`. The bridge above shows that the same factors are already tied exactly to the conical block normalization and to the cut-derived phase-space density. A future Weil-square construction should therefore track this identity rather than treating the block and spectral sectors independently.
