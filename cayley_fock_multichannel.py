#!/usr/bin/env python3
"""Finite-channel CAR/Koszul experiment for the Codex Cayley--Fock front.

Construct Jordan--Wigner creation/annihilation matrices for n fermionic channels and
verify, for reproducible complex holonomies z_i,

    {a_i, a_j^†} = delta_ij I,
    {a_i, a_j} = 0,
    Q = sum_i z_i a_i^†,
    Q^2 = 0,
    D = Q + Q^†,
    D^2 = (sum_i |z_i|^2) I.

The dimension is 2^n, the same doubling sequence as Cayley--Dickson vector-space
dimensions.  This script tests the CAR/Koszul operator identity; it does not identify the
Cayley--Dickson multiplication law with the exterior algebra.
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


def verify(n: int):
    cs = [jw_create(n, i) for i in range(n)]
    ans = [jw_annihilate(n, i) for i in range(n)]
    ident = np.eye(2**n, dtype=complex)
    zero = np.zeros_like(ident)

    car_err = 0.0
    for i in range(n):
        for j in range(n):
            target = ident if i == j else zero
            car_err = max(car_err, np.linalg.norm(ans[i] @ cs[j] + cs[j] @ ans[i] - target))
            car_err = max(car_err, np.linalg.norm(cs[i] @ cs[j] + cs[j] @ cs[i]))
            car_err = max(car_err, np.linalg.norm(ans[i] @ ans[j] + ans[j] @ ans[i]))

    z = np.array([complex(0.37 + 0.11 * (i + 1), -0.23 + 0.07 * i) for i in range(n)])
    Q = sum((z[i] * cs[i] for i in range(n)), start=np.zeros_like(ident))
    D = Q + Q.conj().T
    energy = float(np.sum(np.abs(z) ** 2))

    q2_err = np.linalg.norm(Q @ Q)
    d2_err = np.linalg.norm(D @ D - energy * ident)
    return n, 2**n, car_err, q2_err, d2_err, energy


if __name__ == "__main__":
    print("n dim CAR_err Q2_err D2_err energy")
    for n in range(1, 7):
        row = verify(n)
        print(f"{row[0]} {row[1]} {row[2]:.3e} {row[3]:.3e} {row[4]:.3e} {row[5]:.12g}")
