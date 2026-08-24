#!/usr/bin/env python3
"""Finite-prime CAR/Hodge experiment and Cayley--Dickson dimension check.

This is discovery code, not a proof.  It constructs the n-channel fermionic
creation/annihilation operators by the Jordan--Wigner representation, checks CAR,
and verifies numerically that for Q = sum_i z_i a_i^dagger,

    D = Q + Q^dagger,      D^2 = (sum_i |z_i|^2) I.

The Hilbert-space dimension is 2^n, exactly the same doubling sequence as the
vector-space dimensions of the Cayley--Dickson tower.  The script does NOT identify
the two multiplication laws.
"""

from __future__ import annotations

import numpy as np

I2 = np.eye(2, dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
CREATE = np.array([[0, 0], [1, 0]], dtype=complex)
ANNIHILATE = CREATE.conj().T


def kron_all(xs):
    out = np.array([[1.0 + 0.0j]])
    for x in xs:
        out = np.kron(out, x)
    return out


def jw_create(n: int, i: int) -> np.ndarray:
    return kron_all([Z if j < i else CREATE if j == i else I2 for j in range(n)])


def jw_annihilate(n: int, i: int) -> np.ndarray:
    return jw_create(n, i).conj().T


def check(n: int, z: np.ndarray, tol: float = 5e-12) -> dict:
    dim = 2**n
    creators = [jw_create(n, i) for i in range(n)]
    annih = [c.conj().T for c in creators]
    eye = np.eye(dim, dtype=complex)
    zero = np.zeros((dim, dim), dtype=complex)

    car_err = 0.0
    for i in range(n):
        for j in range(n):
            target = eye if i == j else zero
            car_err = max(car_err, np.linalg.norm(annih[i] @ creators[j] + creators[j] @ annih[i] - target))
            car_err = max(car_err, np.linalg.norm(creators[i] @ creators[j] + creators[j] @ creators[i]))
            car_err = max(car_err, np.linalg.norm(annih[i] @ annih[j] + annih[j] @ annih[i]))

    Q = sum((z[i] * creators[i] for i in range(n)), start=zero.copy())
    nilpotence_err = np.linalg.norm(Q @ Q)
    D = Q + Q.conj().T
    energy = float(np.sum(np.abs(z) ** 2))
    square_err = np.linalg.norm(D @ D - energy * eye)
    eig = np.linalg.eigvalsh(D)

    assert car_err < tol, (n, car_err)
    assert nilpotence_err < tol, (n, nilpotence_err)
    assert square_err < tol, (n, square_err)

    return {
        "n": n,
        "fock_dim": dim,
        "cayley_dickson_dim": 2**n,
        "car_error": car_err,
        "Q_squared_error": nilpotence_err,
        "D_squared_error": square_err,
        "energy": energy,
        "eig_min": float(eig[0]),
        "eig_max": float(eig[-1]),
    }


def main():
    rng = np.random.default_rng(20260823)
    for n in range(1, 7):
        z = rng.normal(size=n) + 1j * rng.normal(size=n)
        print(check(n, z))


if __name__ == "__main__":
    main()
