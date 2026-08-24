#!/usr/bin/env python3
"""Precision and residual audit for the truncated Weil ground state.

The original ``E(c)`` scan requested ``dps=60`` and reported a negative
``lambda_1`` at ``c=29, N=28`` of size about ``10^-61``.  A sign at or below
the requested decimal precision is not meaningful without a precision sweep.

This script rebuilds the exact Connes--van Suijlekom Galerkin matrix using the
``connes-cvs==0.2.2`` implementation, parallelizing only the independent
``(psi(k), psi'(k))`` evaluations.  It then reports:

* the two lowest even-sector eigenvalues and their gap;
* the ground-state Rayleigh quotient and eigenpair residual;
* symmetry and parity-commutator defects of the assembled matrix;
* whether high-precision Cholesky succeeds, and its reconstruction defect.

These are numerical diagnostics, not interval certificates.  A stable positive
value plus a small residual rules out the old negative value as a stable mode,
but a proof of positivity still requires rigorous enclosures.
"""

from __future__ import annotations

import argparse
import json
import multiprocessing as multiprocessing_pool
import platform
import time

import mpmath as mp

from connes_cvs import __version__ as connes_cvs_version
from connes_cvs.operator import HAS_FLINT
from connes_cvs.sweep import _compute_psi_pair_worker, _init_worker

if HAS_FLINT:
    from flint import ctx as flint_ctx


def build_matrix_parallel(c: int, n: int, t_cutoff: int, dps: int, workers: int) -> mp.matrix:
    """Reproduce ``build_galerkin_matrix`` with parallel psi-pair evaluation."""
    mp.mp.dps = dps
    if HAS_FLINT:
        flint_ctx.prec = int(dps * 3.5)

    indices = list(range(-n, n + 1))
    with multiprocessing_pool.Pool(
        workers,
        initializer=_init_worker,
        initargs=(c, dps, t_cutoff),
    ) as pool:
        pairs = pool.map(_compute_psi_pair_worker, indices)

    psi = {k: mp.mpf(value) for k, value, _ in pairs}
    psi_deriv = {k: mp.mpf(value) for k, _, value in pairs}
    dim = 2 * n + 1
    matrix = mp.matrix(dim, dim)
    for i in range(dim):
        m_index = i - n
        for j in range(dim):
            n_index = j - n
            if m_index == n_index:
                matrix[i, j] = psi_deriv[n_index]
            else:
                matrix[i, j] = (psi[m_index] - psi[n_index]) / (m_index - n_index)

    for i in range(dim):
        for j in range(i + 1, dim):
            average = (matrix[i, j] + matrix[j, i]) / 2
            matrix[i, j] = average
            matrix[j, i] = average
    return matrix


def even_projection(matrix: mp.matrix) -> tuple[mp.matrix, mp.matrix]:
    """Return the even-sector matrix and its isometric embedding."""
    dim = matrix.rows
    n = (dim - 1) // 2
    embedding = mp.matrix(dim, n + 1)
    embedding[n, 0] = 1
    inv_sqrt_two = 1 / mp.sqrt(2)
    for k in range(1, n + 1):
        embedding[n + k, k] = inv_sqrt_two
        embedding[n - k, k] = inv_sqrt_two
    return embedding.T * matrix * embedding, embedding


def vector_norm(vector: mp.matrix) -> mp.mpf:
    return mp.sqrt(mp.fsum(abs(vector[i]) ** 2 for i in range(vector.rows)))


def matrix_frobenius_norm(matrix: mp.matrix) -> mp.mpf:
    return mp.sqrt(
        mp.fsum(abs(matrix[i, j]) ** 2 for i in range(matrix.rows) for j in range(matrix.cols))
    )


def max_abs_matrix(matrix: mp.matrix) -> mp.mpf:
    return max(abs(matrix[i, j]) for i in range(matrix.rows) for j in range(matrix.cols))


def as_decimal(value: mp.mpf, digits: int = 50) -> str:
    return mp.nstr(value, digits)


def audit(c: int, n: int, t_cutoff: int, dps: int, workers: int) -> dict[str, object]:
    started = time.time()
    matrix = build_matrix_parallel(c, n, t_cutoff, dps, workers)
    build_seconds = time.time() - started

    dim = matrix.rows
    symmetry_defect = max(
        abs(matrix[i, j] - matrix[j, i]) for i in range(dim) for j in range(dim)
    )
    parity_defect = max(
        abs(matrix[i, j] - matrix[dim - 1 - i, dim - 1 - j])
        for i in range(dim)
        for j in range(dim)
    )

    even_matrix, _ = even_projection(matrix)
    eig_started = time.time()
    eigenvalues, eigenvectors = mp.eigsy(even_matrix)
    eig_seconds = time.time() - eig_started
    lambda_one = eigenvalues[0]
    lambda_two = eigenvalues[1]
    ground_vector = eigenvectors[:, 0]
    residual_vector = even_matrix * ground_vector - lambda_one * ground_vector
    residual = vector_norm(residual_vector)
    vector_size = vector_norm(ground_vector)
    matrix_size = matrix_frobenius_norm(even_matrix)
    residual_scale = matrix_size * vector_size + abs(lambda_one) * vector_size
    relative_residual = residual / residual_scale if residual_scale != 0 else mp.mpf(0)
    rayleigh = (ground_vector.T * even_matrix * ground_vector)[0] / (
        ground_vector.T * ground_vector
    )[0]

    cholesky_ok = False
    cholesky_min_diagonal = None
    cholesky_defect = None
    try:
        lower = mp.cholesky(even_matrix)
        cholesky_ok = True
        cholesky_min_diagonal = min(lower[i, i] for i in range(lower.rows))
        cholesky_defect = max_abs_matrix(even_matrix - lower * lower.T)
    except ValueError:
        pass

    lambda_log10 = mp.log10(abs(lambda_one)) if lambda_one != 0 else mp.ninf
    return {
        "c": c,
        "N": n,
        "T": t_cutoff,
        "dps": dps,
        "workers": workers,
        "dim": dim,
        "lambda_1": as_decimal(lambda_one, min(dps, 100)),
        "lambda_2": as_decimal(lambda_two, min(dps, 100)),
        "gap": as_decimal(lambda_two - lambda_one, min(dps, 100)),
        "lambda_1_sign": int(mp.sign(lambda_one)),
        "lambda_1_log10_abs": float(lambda_log10),
        "digits_beyond_lambda_scale": float(dps + lambda_log10),
        "rayleigh_minus_lambda_1": as_decimal(rayleigh - lambda_one),
        "eigenpair_residual_2": as_decimal(residual),
        "eigenpair_relative_residual": as_decimal(relative_residual),
        "matrix_frobenius_norm": as_decimal(matrix_size),
        "symmetry_defect_max": as_decimal(symmetry_defect),
        "parity_defect_max": as_decimal(parity_defect),
        "cholesky_ok": cholesky_ok,
        "cholesky_min_diagonal": (
            as_decimal(cholesky_min_diagonal) if cholesky_min_diagonal is not None else None
        ),
        "cholesky_reconstruction_defect_max": (
            as_decimal(cholesky_defect) if cholesky_defect is not None else None
        ),
        "build_seconds": round(build_seconds, 3),
        "eigensolve_seconds": round(eig_seconds, 3),
        "total_seconds": round(time.time() - started, 3),
        "python": platform.python_version(),
        "mpmath": mp.__version__,
        "connes_cvs": connes_cvs_version,
        "python_flint": HAS_FLINT,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--c", type=int, default=29)
    parser.add_argument("--N", type=int, default=28)
    parser.add_argument("--T", type=int, default=300)
    parser.add_argument("--dps", type=int, default=120)
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    print(json.dumps(audit(args.c, args.N, args.T, args.dps, args.workers), indent=2))


if __name__ == "__main__":
    main()
