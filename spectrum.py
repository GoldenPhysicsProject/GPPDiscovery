#!/usr/bin/env python3
"""Dump the full even and odd spectra, not just lambda_1 and lambda_2.

The bottom of the CvS Galerkin spectrum is not a continuum. It is an O(1) bulk plus a
LADDER of exponentially small eigenvalues, separated from the bulk by 2+ decades, with
the same count in each parity sector. This script exposes that; point.py, which reports
only the lowest two even eigenvalues and the lowest odd one, cannot.

Note --c is a float: c is a CUTOFF and is not required to be prime, or an integer. Non-
integer and non-prime cutoffs are the only way to separate what depends on the prime
content {p^k <= c} from what depends on the support length 2 log c, because on the
primes those quantities are nearly collinear.
"""
import argparse, json, math, platform, time
import mpmath as mp
import connes_cvs as cc
import resolution, spacing


def sector_spectrum(Q, parity):
    DIM = Q.rows
    N = (DIM - 1) // 2
    inv_sqrt2 = 1 / mp.sqrt(2)
    if parity == "even":
        V = mp.matrix(DIM, N + 1)
        V[N, 0] = mp.mpf(1)
        for k in range(1, N + 1):
            V[N + k, k] = inv_sqrt2
            V[N - k, k] = inv_sqrt2
    else:
        V = mp.matrix(DIM, N)
        for k in range(1, N + 1):
            V[N + k, k - 1] = inv_sqrt2
            V[N - k, k - 1] = -inv_sqrt2
    eigs, _ = mp.eigsy(V.T * Q * V)
    return sorted([eigs[i] for i in range(len(eigs))])


p = argparse.ArgumentParser()
p.add_argument("--c", type=float, required=True)
p.add_argument("--N", type=int, required=True)
p.add_argument("--T", type=int, default=300)
p.add_argument("--dps", type=int, default=90)
p.add_argument("--keep", type=int, default=10)
p.add_argument("--out", default="spectrum.json")
a = p.parse_args()

t0 = time.time()
mp.mp.dps = a.dps
Q = cc.build_galerkin_matrix(a.c, N=a.N, T=a.T, dps=a.dps)
ev_e = sector_spectrum(Q, "even")
ev_o = sector_spectrum(Q, "odd")
k = a.keep
rec = dict(
    c=a.c, N=a.N, T=a.T, dps=a.dps, dim=2 * a.N + 1,
    margin=resolution.margin(a.c, a.N, a.T),
    headroom=a.N / spacing.n_required(int(a.c)),
    even=[mp.nstr(x, 22) for x in ev_e[:k]],
    odd=[mp.nstr(x, 22) for x in ev_o[:k]],
    even_log10=[float(mp.log(abs(x), 10)) if x != 0 else None for x in ev_e[:k]],
    odd_log10=[float(mp.log(abs(x), 10)) if x != 0 else None for x in ev_o[:k]],
    even_ratio_log10=[float(mp.log(abs(x / ev_e[0]), 10)) for x in ev_e[:k]],
    odd_ratio_log10=[float(mp.log(abs(x / ev_e[0]), 10)) for x in ev_o[:k]],
    seconds=round(time.time() - t0, 1),
    python=platform.python_version(), mpmath=mp.__version__,
)
json.dump(rec, open(a.out, "w"))
print(json.dumps(rec), flush=True)
