# GW-EF numerical controls (Claude Code, 2026-11-10)

Evidence, not proof. Companion to the bridge note `MuseResearch/claude-code-note-2026-11-10-gwef-2-3-attack.md`.

## 1. Convention check of the explicit formula (`ef_convention_check.py`)
`weilForm g = poleTerm g + archimedeanTerm g − primeSide g` (GPPVerify definitions) versus the zero sum
`Σ_ρ Φg(ρ)`, `Φg(s) = ∫ g(u) e^{(s−1/2)u} du`, with `g = autocorr f`, `f = exp(−y²/2a²)`, 30-digit mpmath,
zeros from `mpmath.zetazero` (±γ), all zeros used on the line.

| a | zeros used | weilForm | zero sum | difference |
|---|---|---|---|---|
| 0.2 | 40 | 0.000170051785 | 0.000170051785 | −8e-32 |
| 0.1 | 80 | 0.018811606896 | 0.018811606896 | −1.5e-29 |
| 0.05 | 160 | 0.043135958395 | 0.043135958395 | 0 |

The transform, pole term, sign of the prime side and the archimedean term are mutually consistent.
Limit: on-line zeros only; says nothing about off-line behaviour.

## 2. GW-EF-2 factorisation (`pairing_factorisation_check.py`)
`Φ(autocorr f)(s)` versus `Ψf(s)Ψf(1−s)` and `Ψf(s)·conj Ψf(1−s̄)` for the real non-even
`f = (1+y+0.3y²)e^{−y²−0.4y}`, trapezoid h = 1e-3 on [−14,14]. Off the line:
`s = 0.8+3i` (|Φ| = 7.2e-2): 3.0e-13; `0.2−2i` (5.8e-1): 1.1e-12; `1.3+0.7i` (3.3): 1.9e-12.
Both forms agree to the same size. At large height `Φ` is below 1e-16 and the check is vacuous there.

## 3. Off-line pair control (`offline_pair_control.py`) — SYNTHETIC
Zero set = 300 true zeros ∪ an *invented* off-line quadruple `1/2 ± 0.1 ± iγ₁`. Packets
`u·e^{−u²/2σ²}cos(γ₁u)` ("odd") and `e^{−u²/2σ²}cos(γ₁u)` ("even"). `Q = Σ_ρ F(ρ)F(1−ρ)`.

| σ | odd: Q_true | odd: Q_pair | even: Q_total |
|---|---|---|---|
| 1 | 4e-19 | −6.3e-2 | +9.5 |
| 3 | 6e-26 | −5.0e+1 | +9.0e1 |
| 8 | 3e-24 | −3.1e+4 | +9.6e2 |

The sign of the pairing contribution of an off-line pair is controlled by the datum; this is the
mechanism of the density theorem, illustrated. It is not evidence about zeta (the pair is invented).

## Not run
K = 2 sector reproduction (needs the banked sector matrices).
