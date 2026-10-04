> **RETRACTED (2026-10-04, later the same day).** The rigidity claim below is false. Solving the dual problem directly
> (prime_blind.py) finds coefficient vectors far from Lambda (for example b(8) = 4.41 at lambda = 3) that keep the odd
> form positive with margin 6.2e-3, versus 2e-13 at Lambda; verified independently, including on the same 14-mode
> basis. The small reaches measured here reflect Lambda sitting on the boundary of F_lambda at a sharp corner, not
> F_lambda being small. Uniqueness (R) fails at every lambda tested, so the conclusion "RH is equivalent to (E)" does
> not follow. What survives: Lambda lies on the boundary of F_lambda with essentially zero margin, and the
> Davenport-Heilbronn controls. See 2026-10-04_correction_lambda_is_a_boundary_point.md.

# The prime-blind dual: RH as rigidity plus an archimedean statement

Date: 2026-10-04. Author: Claude, for the Golden Physics Project. Status: a reformulation with numerical support at 50 digits. One step is a theorem target (R), one is the remaining open core (E). No proof of RH.

## 1. Setup

Odd tests on [-L, L], L = log lambda. For any real coefficient vector b on {2 <= n < lambda^2}:

    A_lambda(b) = W_inf - sum_n b(n) W_n,   W_n = (2/sqrt n) G_n,

where G_n is the autocorrelation evaluated at log n and W_inf is the archimedean (Gamma) form. The feasible set

    F_lambda = { b : A_lambda(b) >= 0 }

is convex and nested in lambda. RH is equivalent to Lambda in F_lambda for all lambda.

## 2. Rigidity (R): F_lambda collapses onto Lambda

**Locking table** (locking.py): the largest one-sided move of a single coefficient that keeps A PSD.

| lambda | b(6) | b(8) | b(12) | b(15) | b(24) |
| --- | --- | --- | --- | --- | --- |
| 2.5 | 0.34 | absent | absent | absent | absent |
| 3 | 2.9e-11 | 2.7e-5 | absent | absent | absent |
| 4 | 6.9e-20 | 2.2e-17 | 5.4e-11 | 9.3e-4 | absent |
| 5 | 1.3e-24 | 5.2e-23 | 3.2e-19 | 4.4e-16 | 0.023 |
| 6 | 3.3e-27 | 5.0e-26 | 2.8e-23 | 4.1e-21 | 3.3e-14 |

Each coefficient has slack only while log n sits at the edge of the support, then locks. Random directions chosen to raise every near-null eigenvalue at once (maxreach.py) never reach beyond 3e-7.

**Truncation is conservative.** Positivity on a larger test space is a stronger condition, so the true F_lambda is contained in the 14-mode one. These numbers are upper bounds on the true slack, up to floating error.

**Theorem target (R).** For each n, sup over b in F_lambda of |b(n) - Lambda(n)| tends to 0 as lambda grows. Proposed mechanism: on the subspace K of tests whose transforms vanish at the first zeros, the explicit formula gives A_lambda(b) = Q + sum (Lambda - b)(n) W_n with Q tiny on K (unconditionally, from the tail of the zero sum). If no nonzero combination sum delta(n) W_n restricted to K is positive semidefinite, which by the theorem of alternatives means some positive definite Z on K is orthogonal to every W_n, then delta is forced below |Q on K| divided by the margin. That condition is a finite, checkable statement about zero geometry and prime shifts. Grade D.

## 3. Duality: existence becomes prime-free (E)

By the theorem of alternatives for linear matrix inequalities (subject to the usual regularity condition), F_lambda is empty exactly when there is a positive semidefinite Z (a mixed test state) with

    tr(Z W_n) = 0 for every 2 <= n < lambda^2,   tr(Z W_inf) < 0.

Since tr(Z W_n) = (2/sqrt n) G_Z(log n), the first condition says the state is blind to every log n. So:

**(E)** For every lambda, every prime-blind mixed state (a sum of autocorrelations of odd tests on [-log lambda, log lambda] that vanishes at every log n) has nonnegative archimedean energy.

(E) contains no primes, no von Mangoldt weights and no zeros: only the Gamma factor and the positions log 2, log 3, log 4, and so on. When lambda < sqrt 2, no log n fits in the support, every state is prime-blind, and (E) is archimedean Weil positivity on a short interval. That is the regime Connes and Consani proved (check their exact constant before citing). (E) is the natural extension of their theorem to every lambda.

## 4. The chain

- RH implies (E) directly: for prime-blind Z, tr(Z W_inf) = tr(Z A(Lambda)) >= 0.
- (E) implies F_lambda is nonempty for every lambda (duality).
- Nonempty, nested and, by (R), shrinking onto Lambda, the sets F_lambda force Lambda into every F_lambda. That is RH.

So, given (R): **RH is equivalent to (E).** The identification problem disappears because (R) centers the collapse at Lambda itself.

## 5. Consistency check at a second archimedean place

With the Davenport-Heilbronn archimedean place (Gamma((s+1)/2), q = 5) held fixed (dh_arch_feasible.py):

| lambda | DH coefficients, min eigenvalue | Lambda(n) Re chi(n), min eigenvalue |
| --- | --- | --- |
| 6 | -3.2e-9 | 3.5e-6 |
| 8 | -1.01 | 6.0e-10 |
| 10 | -1.38 | 1.7e-13 |
| 12 | -1.78 | 1.0e-13 |

The Euler-product point (the coefficients of the square root of L(chi)L(chi-bar)) is feasible; DH's coefficients are not. The feasible set at a given archimedean place picks out Euler-product data, as (R) predicts. Note the half-integer zero multiplicities of that point: (R)-type uniqueness is about the coefficient vector, not about integrality.

## 6. Status

| item | grade |
| --- | --- |
| F_lambda convex, nested; truncation conservative | A |
| Duality between F_lambda nonempty and (E) | A, given the standard LMI regularity condition |
| Locking and maxreach numerics | numerical, 50 digits, 14 modes |
| (R) rigidity theorem | D, mechanism stated |
| (E) prime-blind archimedean positivity for all lambda | D, the remaining core; proven for lambda < sqrt 2 |

## Files

locking.py, maxreach.py, rigidity.py, rigidity_joint.py, dh_arch_feasible.py.
