# Number thermodynamics on the convergent half-plane

For real beta > 1, define the canonical distribution on positive integers

`P_beta(n) = n^(-beta) / zeta(beta)`.

This is normalized because `zeta(beta)=sum_{n>=1} n^(-beta)`.

Interpret the arithmetic energy as

`E_n = log n`, so that `exp(-beta E_n)=n^(-beta)`.

Then the exact canonical thermodynamic quantities are:

- Partition function: `Z(beta)=zeta(beta)`.
- Free energy: `F(beta)=-(1/beta) log zeta(beta)`.
- Mean arithmetic energy:

  `U(beta) = E_beta[log n] = - d/d beta log zeta(beta) = -zeta'(beta)/zeta(beta)`.

- Shannon/Gibbs entropy:

  `S_N(beta) = -sum_n P_beta(n) log P_beta(n)`

  `= log zeta(beta) + beta U(beta)`

  `= log zeta(beta) - beta zeta'(beta)/zeta(beta)`.

- Energy fluctuation / susceptibility:

  `Var_beta(log n) = d^2/d beta^2 log zeta(beta)`

  `= zeta''(beta)/zeta(beta) - (zeta'(beta)/zeta(beta))^2 >= 0`.

Thus `log zeta(beta)` is convex on beta>1.  This variance is also the Fisher information of the one-parameter Gibbs family with respect to beta, giving a canonical information-geometric metric on the real thermodynamic domain.

Unique factorization makes the same ensemble a product of independent geometric prime occupations. If

`n = prod_p p^(k_p)`, then

`P_beta(k_p=k) = (1-p^(-beta)) p^(-beta k)`

and

`E[k_p] = 1/(p^beta-1)`.

The total mean energy decomposes exactly as

`U(beta)=sum_p log(p)/(p^beta-1)=sum_{p,k>=1} log(p) p^(-k beta) = -zeta'/zeta(beta)`.

This is the precise sense in which the Euler product is a noninteracting bosonic prime gas in its convergence domain.

## Physics dictionary to test, not assume

- `log n`: additive position/energy coordinate of multiplicative number space.
- prime `p`: elementary bosonic mode with energy `log p`.
- `log zeta`: pressure / Massieu potential.
- `d^2 log zeta`: fluctuation susceptibility / information metric.
- beta=1 pole: limiting/Hagedorn-type singularity of the naive prime gas.
- analytic completion + functional equation: candidate interacting/renormalized continuation beyond the naive thermodynamic domain.

The last two lines are interpretations. All formulas above them are exact standard consequences of the Dirichlet series and canonical Gibbs calculus for beta>1.

## Hankel/Gram positivity check (supports `formalization_queue` item `22070a50`)

The higher-moment hierarchy `m_r(sigma) = (-1)^r F^(r)(sigma)` where
`F(sigma) = -zeta'/zeta(1+sigma) = U(1+sigma)` (mean arithmetic energy above) is
structurally a moment/Gram matrix: for any real vector `c_0..c_N`,

`sum_{i,j} c_i c_j m_{i+j}(sigma) = sum_n Lambda(n) n^{-(1+sigma)} (sum_i c_i (log n)^i)^2 >= 0`,

so the Hankel matrix `(m_{i+j})_{0<=i,j<=N}` is positive semidefinite for every
`sigma > 0`, strictly positive definite in practice since the polynomial vanishes on at
most finitely many `log n` values. `discovery/number_thermodynamics/hankel_check.py`
numerically confirms this: truncating the von Mangoldt sum at `n_max=200000` and testing
`sigma in {2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01}`, every leading principal minor of
the `7x7` Hankel matrix (`N=6`) is strictly positive at every tested `sigma`, including
deep into the `sigma -> 0+` regime where convergence is slowest. This is a numeric sanity
check ahead of the Lean formalization, not a substitute for it -- the actual proof route
is the sum-of-squares argument above, not determinant expansion (which is numerically
unstable and combinatorially blows up for larger `N`).
