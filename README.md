# GPPDiscovery

Numeric and exploratory research workbench for the Golden Physics Project's shadow
framework — the standing discovery repo behind
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify), which carries the
formalized, Lean 4-checked results.

## The framework, briefly

The central claim, developed at length in Daniel Toupin's *On the Nature of Nature*
(book manuscript, evolving; companion papers and preprints in Drive; running commentary
at [goldenphysics.org](https://goldenphysics.org)): shadow symmetry in celestial
holography — the conformal-dimension involution `Δ ↔ 2-Δ` — is time reversal `T`. This
identification traces to a single root, self-dual Haar measure on `(R⁺, ×)`
(`dμ(ω) = dμ(ω⁻¹)`), grown by three Cayley–Dickson doublings into the Grassmannian
`Gr(2,4)`. Under the identification `Δ = 2s`, shadow symmetry becomes the Riemann zeta
functional equation `s ↔ 1-s`, and celestial unitarity `Re(Δ) = 1` becomes the critical
line `Re(s) = 1/2` — the same fixed-point set of the same involution, read in two
languages. The same root structure is argued to reach the Birch–Swinnerton-Dyer
conjecture, the Yang–Mills mass gap, the Born rule and measurement as Haar projection,
and the Standard Model's gauge group and generation count from the division-algebra
tower.

This is a live research program, not a closed result, and its own primary source keeps
an explicit ledger of what's rigorously proven, what's argued but not yet peer-reviewed,
what's open, and what's conjectural — this repo inherits that discipline rather than
restating the book's claims as settled. **Nothing in this repo is proved.** Numeric
evidence, derivations worked by hand, and literature checks live here; the moment a
result is solid enough to state as an actual theorem, it gets formalized in Lean in
GPPVerify (no `sorry`, no axiom asserting the open claim, CI-verified — see that repo's
own `CLAUDE.md`). This repo's own `CLAUDE.md` documents the discovery → formalization
workflow and the branch-hygiene discipline that keeps work from getting orphaned between
the two.

## Active threads

- **`weil_decay/`** (root scripts + `discovery/weil_decay/`) — the truncated Weil
  quadratic form. Connes–van Suijlekom and Connes–Consani–Moscovici build, for a prime
  cutoff `c` and band `N`, a finite Galerkin matrix `Q(c)` whose zeros provably sit on
  the critical line for every finite `c`; convergence as `N → ∞` is the open question
  this thread measures. See "The Weil-decay question" below for the current numbers.
- **`discovery/celestial_box/`** — extending the audited scalar-box cut/dispersion
  construction toward pure Einstein gravity; tracking exactly where D-dimensional
  unitarity subtleties (rational terms invisible to 4D cuts) start to matter.
- **`discovery/wiener_hopf/`** — the exact bridge between the principal-series cut
  weight, the Wiener–Hopf Fourier window, and the conical prefactor at
  `s = 1/2 + it`.
- **`discovery/positive_reals_cft/`** — isolating the precise representation-theoretic
  structure behind treating `(R⁺, ×)` as a one-dimensional principal-series system.
- **`discovery/number_thermodynamics/`** — the canonical Gibbs distribution
  `P_β(n) = n⁻β/ζ(β)` on the positive integers and its thermodynamic reading.

## The Weil-decay question

`λ_min(c)`, the smallest even-sector eigenvalue of `Q(c)`, is non-negative for every
finite `c` (Weil positivity, hence RH, is this holding in the limit). We are not testing
whether it's positive — we're measuring **how fast it decays**:

> Is `log λ_min` linear in `log c`, and if so, what is the constant?

Measured so far (unconverged, N=20, T=100): `3.67e-27` at c=7, `8.09e-18` at c=5. A
two-point natural-log slope of about `-64`, against `-2γ₁ = -28.27` and
`-4γ₁ = -56.55`. If the constant is arithmetic, this is a quantity that *sees the
zeros*, which is more than can be said for most reformulations of RH.

**N-convergence comes first.** Going from dim 33 to dim 41 at c=7 moved `λ_min` from
`2.4e-25` to `3.67e-27`, a factor of 65 — N=20 is not converged and any slope fitted to
it is meaningless. The scan does c=7 at N = 20, 28, 36 before anything else.

**Caveat on small c.** At c=3 (N=20, T=100) we saw `λ_min = -0.332` — almost certainly
marginal basis resolution, per the reference implementation's own account at c=23, 29.
The c=29 case has since been resolved at higher precision: stably positive
(`≈1.59e-62`, 120-digit residual `1.64e-121`) at 90–150 digits, though N and T are still
unconverged there too, so this doesn't yet rescue a c-slope fit. See
[`discovery/weil_decay/C29_PRECISION_AUDIT.md`](discovery/weil_decay/C29_PRECISION_AUDIT.md).

## Running

`python point.py --c 7 --N 36 --T 300 --dps 90` for one point, or push / use
`workflow_dispatch` to run the whole matrix in parallel on Actions. Results land in
`results.jsonl` and `RESULTS.md`, committed back automatically.

## Status

Discovery only, across every thread above. See `CLAUDE.md` for the process rules
(branch hygiene, when and how a result graduates to GPPVerify).
