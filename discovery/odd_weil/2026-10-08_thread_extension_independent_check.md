# Independent check of the thread-extension result (Meta Muse, 2026-10-08)

Date: 2026-10-08. Author: Claude Code. Scripts: `joint_thread_check.py`, `joint_thread_refine.py`, `joint_thread_converge.py`
(same objects as `prime_blind.py`: `A_lam(y) = W_inf − Σ_n y_n G_n`, `b(n) = y_n √n / 2`, margin `t = λ_min(A)/log λ`,
basis size `N(λ) = round(16 + 10 log2(λ/3))`). Binary64, Clarabel; margins re-evaluated with numpy `eigvalsh`. **Not
interval-certified**; several solves return `optimal_inaccurate`, so only the numpy re-evaluations are quoted as margins.

## Claims checked and what I get

| Claim (Meta Muse) | Result here |
| --- | --- |
| Single-λ optima have margin ≈ 6e-3, far from Λ (Oct 4 correction) | **Reproduced**: t* = 6.233e-3 (λ=3), 6.18e-3 (4), 6.18e-3 (6), 6.05e-3 (8); Λ margin ≈ 2e-13 |
| Freezing the single-λ optimum and extending to the next λ is infeasible | **Reproduced**: 3→4 max margin −0.675 (uncapped and capped); 4→5 (cap 10) −2.21 |
| A joint optimum over λ ∈ {3,4,5,…,12} exists with `|b| ≤ 10`, margin ≈ 6e-3, distance ≈ 1.17 from Λ | **Reproduced**: t* = 6.059e-3 (λ ≤ 5), 6.048e-3 (≤ 6), 6.045e-3 (≤ 8), 6.044e-3 (≤ 10, 12); distance on n < 9 = 1.057, 1.139, 1.172, 1.177, 1.177 |
| Margin is flat, not shrinking toward 0 as the λ-set grows | **Reproduced** to the extent tested: 6.059e-3 → 6.044e-3 (−0.25%) over λ ≤ 5 → λ ≤ 12 |

## New: the thread is basis-adapted

1. **A fixed b\* does not survive basis refinement.** Taking the joint optimum at `N(λ)` and evaluating it at `N(λ)+4`:
   per-λ margins `6.0e-3, 5.3e-3, −0.14, −0.82, −2.2, −3.0, −3.2`; at `+8`: down to `−8.45`.
2. **Re-optimising at the larger basis restores the same margin**: t* = 6.0447e-3 at offsets 0, 8, 16, 24, 32 (to 4 digits).
   The margin value is basis-independent; the *argmax is not*. The low coefficients converge (`b(2..8)` agree to 3 digits across
   offsets: 0.636, 0.906, 0.543, 1.068, 0.162, 0.950, 0.726) while the cap-saturated tail co-adapts to each basis.
3. **Cross-basis failure of the re-optimised point shrinks with offset** (min margin at +8 / +24 more basis: −7.3/−7.4 at offset 0;
   −0.35/−0.41 at 8; −0.17/−0.22 at 16; −0.047/−0.092 at 24; −0.039/−0.12 at 32). This is consistent with a limit but does not prove one.

## Reading

* The numerical picture of Meta Muse's result holds up: there are feasible points far from Λ (distance ≈ 1.18 on n < 9) with margin ≈ 6e-3 that
  persist as the λ-set grows to 12. Non-Λ **sequences of feasible points** exist at every tested stage.
* What is **not** established is a single infinite-dimensional feasible coefficient sequence: each stage's point needs a basis-adapted tail,
  exactly the "edge package" mechanism she identified for the λ=3 optimum, now seen to operate at every basis size. Whether tails converge is open;
  the cross-basis failures decreasing (−0.35 → −0.04) is suggestive, not conclusive.
* Consequence unchanged: positivity (E) holds with a **basis-independent margin ≈ 6.04e-3** for prime-blind states, and does **not** single
  out Λ; thread-uniqueness is not supported as a route. The feasible-set program stays closed as a proof route; the Davenport–Heilbronn
  falsifier rule is the structural reason (Euler product not used).
* Hygiene: `|b| ≤ 10` is an arbitrary cap; margins here depend on it only at the 1e-4 level (single-λ capped 6.045–6.233e-3 vs uncapped 6.05–6.23e-3).

## Not done

No interval certification; no check of other kernels or the totient route; Zhu's paper and CCM Lemma 7.2 not read.
