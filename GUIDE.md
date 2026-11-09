# Guide to GPPDiscovery

GPPDiscovery is the standing numeric workbench of the Golden Physics Project and Claude's seat repository. The framework is shadow holography: one involution, the celestial shadow map Δ ↦ 2 − Δ and the zeta functional equation s ↦ 1 − s, read as a single dictionary between arithmetic and celestial scattering, with the Riemann Hypothesis as the open reflection-positivity question inside it. The full write-up is `discovery/arithmetic_holography/ARITHMETIC_HOLOGRAPHY_RIGOROUS.md` in this repository.

This guide is the map: where things are, how results are recorded, and how
the repository relates to the rest of the Golden Physics Project. It
describes the layout as of 2026-11-09. The workbench moves; the map is
revised when the layout changes.

## Start here

- `README.md`: the repository's own introduction.
- `RESULTS.md`: the human-readable results ledger.
- `results.jsonl`: the raw result rows behind the ledger.
- `E_RESULTS.md` and `E_results.jsonl`: the E-series results ledger and its raw rows.
- `CLAUDE.md`: the seat's operating notes for this repository.
## Layout

- Root scripts: `spectrum.py`, `spacing.py`, `fit.py`, `efit.py`, `point.py`, `resolution.py`, and the Cayley–Fock bridge scripts `cayley_fock_bridge.py`, `cayley_fock_multichannel.py`.
- `discovery/`: the research threads, including `arithmetic_holography/`.
- `claude/`: the seat's working files.
- `lean_sweep/`: Lean-side sweep results.
- `results/`: saved outputs.

## Conventions

Results carry epistemic labels: THEOREM, PROVED LEMMA, CONDITIONAL,
NUMERICAL, OPEN, KILLED. Numerical work is evidence, not proof; what
hardens is promoted to
[GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify), the
kernel-checked Lean 4 library, through the build seat's audit. Candidate
positivity mechanisms are tested against the controls F_θ = ξ(s−θ)ξ(s+θ)
and the Davenport–Heilbronn function before they are taken seriously.
Working notes, heuristics, and internal strategy live in the private
GPP-bridge repository; this repository is public by design and is the
publication of the workbench's math, code, and results.

## The family

- [GPPVerify](https://github.com/GoldenPhysicsProject/GPPVerify): the formalization library.
- [GPPDiscovery](https://github.com/GoldenPhysicsProject/GPPDiscovery): Claude's workbench.
- [GPPDiscovery2](https://github.com/GoldenPhysicsProject/GPPDiscovery2): Codex's workbench.
- [GPPDiscovery3](https://github.com/GoldenPhysicsProject/GPPDiscovery3): Grok's (Astra's) workbench.
- [GPPDiscovery4](https://github.com/GoldenPhysicsProject/GPPDiscovery4): Turbo's workbench.
- [GPPDiscovery5](https://github.com/GoldenPhysicsProject/GPPDiscovery5): Grok's lab (guest seats).

License: Apache-2.0. Citation: see `CITATION.cff`. Author: Daniel Toupin,
Golden Physics Project.
