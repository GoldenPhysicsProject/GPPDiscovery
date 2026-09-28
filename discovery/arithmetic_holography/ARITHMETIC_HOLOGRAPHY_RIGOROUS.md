# Arithmetic holography — a rigorous formulation

Claude Code, 2026-09-28. The picture is Daniel's and Codex's (GPPDiscovery2 `codex/discovery-workbench`,
09-27 notes on arithmetic holography, Bohr–Hardy, BPY intertwiner, density-matrix Fredholm). This note
makes it precise, separates what is a theorem from what is open, and records a numerical falsifier.

Labels: **THEOREM** (proved, with source), **LEAN** (certified in GPPVerify), **OPEN** (RH-strength),
**KILLED** (falsified here). Nothing below proves RH.

---

## 0. The picture in one paragraph

The spectral variable `s` lives on the Riemann sphere `S² = Ĉ`. The critical line together with `∞`
is a great circle `L` — the equator — and it is the common boundary of two hemispheres
`H₊ = {Re s > 1/2}` and `H₋ = {Re s < 1/2}`. On the other side of a Mellin duality sits the
multiplicative number line `ℝ₊`, split at its fixed point `x = 1` into the two halves `x > 1` (where
the integers live) and `x < 1` (their reciprocals). The equator is the hologram of the point `x = 1`,
the hemispheres are the holograms of the two half-lines, and the mirror symmetries correspond
(`x ↦ 1/x` ↔ `s ↦ 1 − s̄`). Inside that, the integers are themselves the Fourier modes of a circle,
and the prime structure lives on a compact prime torus into which the critical line embeds densely:
a hologram inside a hologram. The functional equation says the hologram is **symmetric**. RH says it is
also **positive** (reflection-positive, i.e. unitary). Everything below makes these sentences exact.

---

## 1. The sphere, the equator, and the mirror

Let `τ(s) = 1 − s̄`. It is an antiholomorphic involution of `Ĉ` with fixed set exactly `L`, and it
exchanges `H₊` and `H₋`. So `Ĉ` is the Schottky double of `H₊` across `L`.

**THEOREM 1 (mirror symmetry).** `ξ(τ s) = conj ξ(s)`.
*Proof.* `ξ(1 − s) = ξ(s)` (functional equation) and `ξ(s̄) = conj ξ(s)` (reality). Compose. ∎

Consequences:
- `ξ` is real on the equator (`Ξ(t) = ξ(1/2 + it) ∈ ℝ`).
- The zero set `Z` is `τ`-invariant. A zero on `L` is `τ`-fixed (a *boundary mode*). An off-line zero
  `ρ` comes with its distinct mirror `τρ = 1 − ρ̄` (a *bulk mirror pair*).
- **RH ⟺ `Z ⊂ Fix(τ) = L`.**

**LEAN.** Zero-pairing under `ρ ↦ 1 − ρ̄`: `SpectralWeil.fePartner_zero_of_strip`; Mathlib has
`completedRiemannZeta_one_sub` and `riemannZeta_conj`.

---

## 2. First hologram: the multiplicative line and the sphere (Mellin)

Let `ℝ₊` act on half-densities `L²(ℝ₊, dx)` by unitary dilations `(U_a f)(x) = a^{1/2} f(ax)`.

**THEOREM 2 (Mellin–Plancherel).** `Mf(s) = ∫ f(x) x^{s−1} dx` is unitary (up to `2π`) from
`L²(ℝ₊, dx)` onto `L²(L)`. The unitary characters of the dilation group on half-densities are exactly
`x^{−s}` with `Re s = 1/2`. *The equator is the unitary dual of scaling.*

**THEOREM 3 (Paley–Wiener, Mellin form).** `supp f ⊂ [1, ∞)` iff `Mf` extends to the Hardy space
`H²(H₋)`; `supp f ⊂ (0, 1]` iff `Mf ∈ H²(H₊)`. Hence
`L²(ℝ₊) = L²(0,1] ⊕ L²[1,∞)` is carried onto `L²(L) = H²(H₊) ⊕ H²(H₋)`.
(Standard: the Laplace form of Paley–Wiener in `u = log x`. Which hemisphere matches which half
depends only on the sign convention in `x^{s−1}` versus `x^{−s}`.)

**THEOREM 4 (the mirror is the inversion).** `(Jf)(x) = x^{−1} f(1/x)` is a unitary involution of
`L²(ℝ₊, dx)`. It exchanges the two half-lines, fixes `x = 1`, and satisfies `M(Jf)(s) = Mf(1 − s)`.
With complex conjugation it realizes `τ`.

**Dictionary (exact):**

| boundary (multiplicative line) | bulk (spectral sphere) |
|---|---|
| dilation group `ℝ₊` on half-densities | its unitary dual = equator `L` |
| half-line `x > 1` / `x < 1` | hemisphere `H₋` / `H₊` (Hardy spaces) |
| fixed point `x = 1` | fixed circle `L` |
| inversion `x ↦ 1/x` | reflection `s ↦ 1 − s̄` |
| support in a half-line | analyticity in a hemisphere |

---

## 3. Second hologram: integers as modes, primes as a torus

**3a. The circle and theta (Riemann).** `ψ(x) = Σ_{n≥1} e^{−π n² x}` is the heat trace of a circle
whose Fourier modes are the integers (`ℤ` = Pontryagin dual of `S¹`).

**THEOREM 5 (boundary self-duality ⟹ bulk mirror symmetry).** Poisson summation for `ℤ ⊂ ℝ`, where
`ℤ` is its own dual lattice, gives `θ(1/x) = √x θ(x)` with `θ = 1 + 2ψ`. Splitting the Mellin integral at
`x = 1` and folding the lower half onto the upper half with this identity gives
`ξ(s) = ½ + ½ s(s−1) ∫₁^∞ ψ(x)(x^{s/2} + x^{(1−s)/2}) dx/x`, which is manifestly `τ`-symmetric.
*The bulk mirror symmetry of Theorem 1 is the hologram of the self-duality of the integer lattice.*
(Riemann 1859. The adelic version is Tate's thesis, where `ℚ ⊂ 𝔸` is a self-dual lattice.)

**3b. The prime torus and the dense embedding.** By unique factorization, `ℚ₊^×` is free abelian on
the primes. Its Pontryagin dual is the compact prime torus `K = ∏_p S¹`, and its characters are
indexed by `ℚ₊^×`; the positive integers form the positive cone.

**THEOREM 6 (hologram inside a hologram).** Dualizing the dense inclusion `ℚ₊^× ↪ ℝ₊` gives the dense
Kronecker flow `ℝ → K`, `t ↦ (p^{−it})_p`. The critical line is thus embedded densely in the prime
torus, and its unitary characters pull back to `n^{−it}`: the integers sit densely in the scale line,
and the equator sits densely in the torus. (Pontryagin duality; Kronecker's theorem, since the
`log p` are ℚ-linearly independent.)

**THEOREM 7 (Codex, 09-27; Bohr lift).** On the Hardy space `H²(K)`, `Z_s = Σ n^{−s} χ_n` and
`M_s = Σ μ(n) n^{−s} χ_n` satisfy `Z_s M_s = 1` in `L¹(K)` for all `Re s > 1/2`, with no zero input.
Scalarizing back to `ζ(s)` means evaluating at the identity of `K`, which is unbounded on `H²` for
`Re s ≤ 1` (Hedenmalm–Lindqvist–Seip). *The internal hologram is healthy up to the equator. The
obstruction is entirely in projecting the torus back to the line.*

**THEOREM 8 (Codex, 09-27; connected intertwiner).** `T_s M_s = A_∞(s)·1 + C_s M_s`, where `C_s` from
the prime-torus Hardy space into the BPY Gaussian field is bounded and intertwines prime shifts with
Gaussian decimations. The BPY field carries the functional equation, since `E[Q^{s/2}] = 2ξ(s)`. The
only singular piece is the rank-one coherent channel `A_∞(s)·1`.

---

## 4. Third level: the celestial sphere (structural parallel, not an identification)

The celestial Mellin transform in energy `ω` uses a spectral plane with principal series
`Re Δ = 1` and shadow `Δ ↦ 2 − Δ`. Under `Δ = 2s` these become `L` and `s ↦ 1 − s`
(GPPVerify: `CelestialRiemannAffineBridge`). Both planes are the Mellin duals of a dilation
subgroup: a boost along the scattering direction acts on `ω` by dilation. **Exact:** the arithmetic
sphere is the scale sector of celestial data. **Not established:** that `ξ` is a celestial correlator
of any physical theory. Keep two different spheres apart: the celestial sphere of *directions*
`z ∈ CP¹` is not the spectral `Δ`-sphere; the arithmetic hologram corresponds to the latter.
The nearest exact bridge is the BPY thermal Gaussian field, whose scale-Mellin transform is `ξ`.

---

## 5. RH as a holographic statement — four equivalent rigorous forms

All four are classical or proved in the programme. Each is equivalent to RH, so none is an
unconditional shortcut.

**(A) No bulk mirror pairs:** `Z ⊂ L` (§1).

**(B) Reflection positivity across the equator (Weil 1952).** Let `W` be the explicit-formula
distribution on `ℝ₊`. RH iff `W(f ⋆ f̃) ≥ 0` for all test functions `f`, where
`f̃(x) = conj f(1/x)` is the inversion (Theorem 4). The spectral side is
`Σ_ρ m_ρ F(ρ) conj F(τρ)`, with `F` the Mellin transform of `f`.
**LEAN (finite core):** `GppHolographicRP.reflectionForm_nonneg_iff`. For a finite `τ`-invariant set
with positive multiplicities, `Σ m_z v_z conj v_{τz} ≥ 0` for all `v` iff `τ` fixes every point.
`mirror_pair_indefinite` shows a mirror pair is a hyperbolic plane with one positive and one ghost
direction. *Reflection positivity on the boundary ⟺ spectral support on the equator.*

**(C) Positive bulk vacuum (Codex 09-27, density-matrix equivalence).** RH iff
`ξ(1/2 + z)/ξ(1/2) = det(I + R z² ρ)` for a trace-class `ρ ≥ 0`.
**LEAN (finite core):** `GppPositiveFredholmFactor.fredholm_zero_on_imaginary_axis`.

**(D) Tempered reconstruction (Codex 09-27 + von Koch).** RH iff the primitive boundary distribution
`e^{x/2}dx − Σ_p (log p) p^{−1/2} δ_{log p}` is tempered, iff
`Σ_{p≤x} (log p)/√p = 2√x + O(logᴺ x)`. In holographic terms, the bulk log-field on `H₊` is the
Paley–Wiener extension of tempered half-line boundary data.

---

## 6. What was killed today

**KILLED — "each half dominates its own side".** Riemann's canonical split at `x = 1` gives
`ξ(s) = E(s) + E(1 − s)` with
`E(s) = 1/4 + ½ s(s−1) Σ_n (πn²)^{−s/2} Γ(s/2, πn²)`.
The identity was checked to `10⁻³¹`. If `|E(1−s)| < |E(s)|` held on `H₊`, a Hermite–Biehler argument
would give RH. It fails badly: on `σ ∈ {0.51, …, 1.5}`, `t ∈ [0, 120]`, about 57% of sampled points
violate it. The maximum ratio is at `t ≈ 6` for every `σ`, reaching `2.146` at `σ = 1.5`, deep in the
zero-free region. Script: `riemann_split_hb.py`.

*Lesson.* The geometric midpoint `x = 1` is the right fixed point for the **symmetry**, but the
naive half is the wrong object for **positivity**. Some Hermite–Biehler split does exist if RH holds
(de Branges–Lagarias), but it must be dynamical, not the bare cut. That matches form (B): the
mechanism is a positive *pairing across* the equator, not domination within each hemisphere.

---

## 7. Constraints any closing mechanism must satisfy

1. **Euler product and self-duality together.** Beurling generalized primes have an Euler product
   and a PNT but zeros off the line (Diamond–Montgomery–Vorhauer). Davenport–Heilbronn functions have
   a functional equation but no Euler product, and zeros off the line. So the positivity in (B) must
   couple hologram 3b (primes, torus) to hologram 3a (lattice self-duality) in one step. Theorem 7
   alone does not: it holds for Beurling systems. Theorem 8 is the first object here that couples
   both.
2. **Exactness (zero margin).** De Bruijn–Newman: `Λ ≥ 0` (Rodgers–Tao), and RH ⟺ `Λ = 0`. If RH is
   true, it holds with no room to spare, so the mechanism must be an exact positivity (a norm square,
   a KMS or OS identity), not a coercive estimate with slack. Codex's finite CCM margins shrinking to
   `~10⁻¹⁷` are consistent with this.
3. **Positivity of the density is not enough.** GPPVerify `FisherZeroLogConcavityNoGo`: a positive,
   arbitrarily strongly log-concave density can have off-axis characteristic-function zeros. The
   positivity must be of the *two-point pairing across the mirror* (form (B)), not of a one-point
   density.
4. **No factorization through two bounded maps at half density.** This is Codex's critical-split
   no-go, and the Bohr/Hedenmalm–Lindqvist–Seip evaluation threshold. The projection from the torus
   to the line must be one connected, graded map.

---

## 8. The rigorous programme

The picture is now fully rigorous except at one joint:

> **OPEN (holographic reflection positivity).** Construct, without zero data, a Hilbert space `𝓗`
> with a unitary involution `Θ` implementing `x ↦ 1/x` (or `s ↦ 1 − s̄`), and a vector map
> `f ↦ a_f`, such that `W(f ⋆ f̃) = ⟨a_f, Θ a_f⟩_𝓗` with `Θ ≥ 0` on the relevant cone
> (Osterwalder–Schrader form). By §5(B) this is equivalent to RH.

Why this is the right target in the holographic language:
- In QFT, holographic and OS reconstruction need **reflection symmetry** plus **reflection
  positivity**. The first is Riemann's functional equation (Theorems 1 and 5, proved). The second is
  Weil's criterion (RH).
- The OS-reconstructed Hamiltonian of such an arithmetic boundary theory is automatically
  self-adjoint. Its spectrum is the equator, which is the Hilbert–Pólya operator obtained as a
  *consequence* rather than postulated.
- Physical mechanisms that give positivity *exactly*, as constraint 2 demands:
  1. **KMS midpoint positivity.** For a β-KMS state, `ω(a* σ_{iβ/2}(a)) ≥ 0`. At the Bost–Connes
     critical temperature β = 1, the midpoint `β/2 = 1/2` *is* the critical line. The half-density
     weight is exactly this midpoint (her heuristic 2.6; `KMSCriticalLimit` in Lean). Known gap: the
     Bost–Connes KMS₁ functional is not Weil's `W`; the Archimedean place must be sewn in.
  2. **Gaussian OS positivity.** A Gaussian measure is reflection-positive iff its covariance is. The
     BPY field is Gaussian, and Codex's two-copy BPY colligation already states RH as a contraction
     `‖C_ω‖ ≤ 1` (09-20), which is exactly an OS-positivity statement.
  3. **TFD doubling.** Thermofield-double cross-covariances are reflection-positive with respect to
     the sheet exchange (her 09-27 parity/TFD work, `TFDParitySewing` in Lean). The two TFD sheets are
     the two hemispheres.
- Constraint 1 decides among them. The winning construction must be a product over places (Euler)
  of local reflection-positive factors whose *global* pairing is fixed by the self-dual lattice
  (Poisson). The known failure mode to avoid is lesson 11.3: local prime positivity enters Weil's
  formula with a minus sign. A mechanism that makes each prime's contribution a TFD *cross*-term with
  respect to `Θ`, rather than a diagonal term, is what the three candidates above have in common.

**Next concrete steps:**
1. Write Weil's `W` explicitly as `⟨a_f, Θ a_f⟩` candidates in the TFD/BPY doubled field, and compute
   the defect `W(f⋆f̃) − ⟨a_f, Θ a_f⟩` on finite prime sets. This is an exact finite computation.
2. Apply constraint 1 to each candidate: run the same construction on a Beurling system and on a
   Davenport–Heilbronn function. A candidate that stays "positive" there is carrying no arithmetic.
3. Formalize in Lean each exact finite identity that survives, alongside
   `HolographicReflectionPositivity`, `PositiveFredholmFactor` and `TFDParitySewing`.

No RH proof is claimed. The holographic picture is made rigorous up to one open statement, and that
statement is RH-equivalent.
