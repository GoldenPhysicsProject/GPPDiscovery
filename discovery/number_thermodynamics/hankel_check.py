"""
Numeric check of Hankel/Gram positivity for F(sigma) = -zeta'/zeta(1+sigma).

m_r(sigma) = (-1)^r F^(r)(sigma) = sum_{n>=2} Lambda(n) (log n)^r n^{-(1+sigma)}

formalization_queue item 22070a50 ("Zeta prime-response Hankel positivity") asks for
det (m_{i+j})_{0<=i,j<=N} > 0 for sigma > 0. This is the finite von Mangoldt Hankel/Gram
matrix -- exactly the "energy fluctuation / susceptibility" hierarchy from
NUMBER_GIBBS_ENTROPY.md (m_0=U-ish, m_1=Var, higher m_r are higher cumulant-type moments
of log n under the Gibbs weight n^{-(1+sigma)}/zeta(1+sigma) up to normalization).

Positivity is structurally guaranteed since (m_{i+j}) is literally a moment/Gram matrix:
for any real vector c_0..c_N, sum_{i,j} c_i c_j m_{i+j}
  = sum_n Lambda(n) n^{-(1+sigma)} (sum_i c_i (log n)^i)^2 >= 0,
strictly positive whenever the polynomial sum_i c_i x^i is not identically zero on the
support {log n : n=p^k}. This script checks the numeric determinant sign directly as a
sanity check ahead of/alongside the Lean formalization (which should prove positivity via
exactly this sum-of-squares argument, not by determinant expansion).
"""
import math

def von_mangoldt(n):
    if n < 2:
        return 0.0
    m = n
    for p in range(2, int(math.isqrt(n)) + 1):
        if m % p == 0:
            while m % p == 0:
                m //= p
            return math.log(p) if m == 1 else 0.0
    return math.log(n)  # n itself is prime


def moments(sigma, N, n_max):
    """m_r(sigma) for r = 0..2N, truncated sum over n=2..n_max."""
    m = [0.0] * (2 * N + 1)
    for n in range(2, n_max + 1):
        L = von_mangoldt(n)
        if L == 0.0:
            continue
        ln = math.log(n)
        w = L * n ** (-(1 + sigma))
        lnp = 1.0
        for r in range(2 * N + 1):
            m[r] += w * lnp
            lnp *= ln
    return m


def hankel_leading_minors(m, N):
    """Leading principal minors of the (N+1)x(N+1) Hankel matrix (m_{i+j})."""
    H = [[m[i + j] for j in range(N + 1)] for i in range(N + 1)]
    minors = []
    for k in range(1, N + 2):
        sub = [row[:k] for row in H[:k]]
        minors.append(det(sub))
    return minors


def det(A):
    n = len(A)
    A = [row[:] for row in A]
    d = 1.0
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(A[r][i]))
        if abs(A[piv][i]) < 1e-300:
            return 0.0
        if piv != i:
            A[i], A[piv] = A[piv], A[i]
            d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            for c in range(i, n):
                A[r][c] -= f * A[i][c]
    return d


if __name__ == "__main__":
    N_MAX = 6
    N_TRUNC = 200000
    print(f"n_max={N_TRUNC}, testing Hankel size N+1 up to {N_MAX+1}")
    for sigma in [2.0, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01]:
        m = moments(sigma, N_MAX, N_TRUNC)
        minors = hankel_leading_minors(m, N_MAX)
        signs = ["+" if x > 0 else ("0" if x == 0 else "-") for x in minors]
        print(f"sigma={sigma:<6} minors_sign={''.join(signs)} "
              f"minors={[f'{x:.3e}' for x in minors]}")
