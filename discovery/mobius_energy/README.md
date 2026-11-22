# Möbius energy — numerical receipts (NUMERICAL, not proved)

Everything here is evidence, never a theorem. It tests a picture, nothing more.

Setup: a(n)=μ(n)/√n, b_k = Σ_{n≤k} a(n) − c, c = 1/ζ(1/2) ≈ −0.684765, and
Q(log N) = Σ_{k<N} log((k+1)/k)·b_k² (the energy identity; Lean: `GppMobiusEcho.Q`).
In Lean (GPPVerify #659–#667) hypothesis (i) of `rh_of_storage_work` is proved; RH then follows from
subpower growth of Q along some N_j→∞ (hypothesis (ii)). These scripts look at how Q actually grows.

| file | what | result |
|---|---|---|
| `scan.py` | sieve μ to N=1e7, Q for c=1/ζ(½), 0, −0.5 | c=1/ζ(½): b_k stays bounded (max ≲0.45), Q linear in log N; wrong c grows faster |
| `slope.py` | least-squares slope of Q vs log N | **0.0289–0.0291** (N=1e7 and 6e7, fits from 1e3…1e5); residual rms ≈0.002 |
| `parseval.py` | mpmath, 2Σ_{γ>0} 1/|ρζ′(ρ)|² over first NZ zeros | 0.02531 (10), 0.02695 (25), 0.02776 (50), 0.02829 (100), 0.02858 (200), 0.02870 (300 zeros, γ≤541.8) |

Reading (mine, not verified): the measured slope matches the almost-periodic / explicit-formula prediction
B(u) ≈ Σ_{γ>0} 2Re[e^{iγu}/(ρζ′(ρ))] (assumes RH and simple zeros), whose mean square is
2Σ 1/|ρζ′(ρ)|². The partial sums are still rising toward the measured slope. So in the true world
(ii) looks like a mean-square/almost-periodicity statement about Möbius sums — RH-strength, and
this check says nothing about proving it. Reproduce: `python3 -I scan.py`, `slope.py [NMAX]`, `parseval.py [NZ]`.

## Correction (2026-11-22, same day)

The explicit-formula term above was written with denominator `rho zeta'(rho)`. The correct one is `(rho - 1/2) zeta'(rho) = i gamma zeta'(rho)`:
the pole of `x^s / (s zeta(s + 1/2))` sits at `s = i gamma` (not `s = rho`), so `B(u) = sum_{gamma>0} 2 Re[e^{i gamma u} / (i gamma zeta'(rho))]`
and the mean square is `kappa = 2 sum_{gamma>0} 1/(gamma^2 |zeta'(rho)|^2)`. The change is `~ 1/(4 gamma^2)` relative per zero (0.08% on the first
zero, 0.085% on a 40-zero partial sum) and does not affect the agreement quoted above at the digits shown, but the
`parseval.py` here uses the wrong denominator. The corrected computation, on Platt's tables read through the reference reader, is in
`discovery/mobius_modes/` and `results/mobius_modes/`.
