> **RETRACTED (2026-10-04, later the same day).** The rigidity claim below is false. Solving the dual problem directly
> (prime_blind.py) finds coefficient vectors far from Lambda (for example b(8) = 4.41 at lambda = 3) that keep the odd
> form positive with margin 6.2e-3, versus 2e-13 at Lambda; verified independently, including on the same 14-mode
> basis. The small reaches measured here reflect Lambda sitting on the boundary of F_lambda at a sharp corner, not
> F_lambda being small. Uniqueness (R) fails at every lambda tested, so the conclusion "RH is equivalent to (E)" does
> not follow. What survives: Lambda lies on the boundary of F_lambda with essentially zero margin, and the
> Davenport-Heilbronn controls. See 2026-10-04_correction_lambda_is_a_boundary_point.md.

# Positivity selects the primes: rigidity and an existence-plus-uniqueness program

Date: 2026-10-04. Author: Claude, for the Golden Physics Project. Status: numerical rigidity at 50 digits in a 14-mode truncation, plus a proof program. No proof of RH.

## 1. The rigidity experiment

Fix the zeta Archimedean place (Gamma(s/2), pole terms removed by working in the odd sector). Treat the "prime" coefficients b(n), n < lambda^2, as unknowns:

    A_lambda(b) = W_inf - sum_n b(n) (2/sqrt n) G_n,   G_n = autocorrelation evaluated at log n.

A_lambda(b) is affine in b. Define the feasible set

    F_lambda = { b : A_lambda(b) is positive semidefinite on odd tests supported in [-log lambda, log lambda] }.

Question: how large is F_lambda around b = Lambda?

### Single coefficients (rigidity.py)

Moving one b(n) = Lambda(n) + t, the allowed interval for t at lambda = 4 is:

| n | Lambda(n) | allowed t |
| --- | --- | --- |
| 2 | 0.6931 | [-9.2e-23, 3.1e-26] |
| 3 | 1.0986 | [-4.0e-23, 5.0e-25] |
| 6 | 0 | [-3.9e-23, 6.9e-20] |
| 10 | 0 | [-4.6e-17, 2.0e-14] |
| 12 | 0 | [-2.3e-13, 5.4e-11] |
| 15 | 0 | [-1.4e-5, 9.3e-4] |

Both ends are bounded, for prime powers and composites alike. Slack appears only at the edge n near lambda^2, where the shift log n barely overlaps the support.

### Joint directions (rigidity_joint.py)

Joint perturbations can cancel first-order constraints. Directions in the null space of the near-null sensitivity matrix were tested to second order:

| lambda | coefficients | near-null constraints | free first-order directions | largest two-sided movement found |
| --- | --- | --- | --- | --- |
| 3 | 7 | 3 | 4 | 7.6e-11 |
| 4 | 14 | 5 | 9 | 1.9e-16 |
| 5 | 23 | 6 | 17 | 8.8e-19 |
| 6 | 34 | 7 | 27 | 4.6e-21 |

The pinning tightens as lambda grows. The near-null constraints are the zeta-cycle eigenvalues; each tiny eigenvalue is a half-space that cuts through within its own size of Lambda.

Caveats: N = 14 modes (more modes add constraints, so this understates rigidity); directions tested are null-space eigen-directions, not an exhaustive search of F_lambda; numerics, not interval-certified.

## 2. Reading

Positivity does not merely tolerate the primes. In every direction tested, it selects them: the von Mangoldt weights sit at a near-corner of a convex set that shrinks as lambda grows. Composites are forbidden from both sides. This is the converse of the Davenport-Heilbronn control: there, composite coefficients produced a negative mode; here, any composite coefficient produces one.

In the arrow language: H2 (no D-odd observable) is not just consistent with the Euler product. It numerically forces the Euler-product weights.

## 3. The program this suggests

The sets F_lambda are convex (A is affine in b and the PSD cone is convex) and nested (tests on a smaller interval are tests on a larger one, and see fewer n). If they are also bounded, compactness gives a nonempty intersection whenever each is nonempty:

    F_inf = intersection over lambda of F_lambda.

RH is exactly the statement that Lambda belongs to F_inf. That splits the problem into three pieces with different characters:

1. **Uniqueness (numerically supported here).** F_inf contains at most one point. Target: prove the shrinking from the zeta-cycle eigenvalues. Grade D, but now with a concrete mechanism to analyze.
2. **Existence (the hard core).** Show F_lambda is nonempty for every lambda by a construction that does not assume RH. Equivalently: a positive measure (the "zeros") whose explicit-formula partner is the zeta Archimedean distribution minus a measure supported on {log n}. This is a Fourier-pair problem of crystalline type. Kurasov-Sarnak built positive crystalline measures from Lee-Yang stable polynomials, and Olevskii-Ulanovskii characterized unit-mass Fourier quasicrystals; check precisely how far their constructions reach once an absolutely continuous Archimedean part is allowed. This is where the Lee-Yang route from the v2 roadmap plugs in. Grade D.
3. **Identification (classical).** Show the unique point of F_inf is Lambda. Hamburger's theorem (1921): a Dirichlet series with zeta's functional equation and finitely many poles is a constant multiple of zeta. The remaining gap is integrality: one needs the positive zero measure to have integer masses so that a Dirichlet series can be rebuilt from it. Grade D, but a well-posed classical question.

If 1 and 3 hold, then the existence of any feasible point at all, by any construction, would force it to be Lambda, and RH would follow. Positivity would then be inherited from existence rather than proved by estimates. This is the "everything equals zero" picture in operational form: the arithmetic balance is unique, so if any balanced object exists, it is the primes.

## 4. Next steps

- Repeat rigidity_joint.py with N = 24 and 30 at 60 digits, and a convex-optimization search (maximize distance from Lambda inside F_lambda) instead of fixed directions.
- Test uniqueness against Davenport-Heilbronn's Archimedean place (Gamma((s+1)/2), q = 5): does F_inf for that place contain L(chi_5)-type points and exclude DH's coefficients?
- Literature check on crystalline measures with an absolutely continuous Fourier part.

## Files

rigidity.py (single-coefficient intervals), rigidity_joint.py (null-space directions); both use mp_N14.pkl and mp_weil.py.
