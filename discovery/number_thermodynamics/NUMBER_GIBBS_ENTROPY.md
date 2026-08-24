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
