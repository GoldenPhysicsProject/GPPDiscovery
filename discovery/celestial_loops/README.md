# celestial_loops — Δ₅ + Δ₆ = 2 shadow sewing

The full record (claim, status table, what was wrong in old drafts, next steps) lives in
`GPP-bridge/CLAUDE_CELESTIAL_LOOPS.md`. This folder holds the code.

- `check_sewing.py`: exact pieces only. Γ(1+iλ)Γ(1−iλ) = πλ/sinh πλ; P̂(ξ) = (π/2)sech²(ξ/2);
  Mellin–Plancherel sewing with weight ω^k dω pairs dimensions on Δ₅+Δ₆ = k+1 (massless phase space: k = 1, so 2).
- Runs on `.github/workflows/claude-celestial.yml`; output `results/sewing.json`.

Nothing here is proved; the loop-extraction claim beyond the cut is argued/open (see the bridge file).
