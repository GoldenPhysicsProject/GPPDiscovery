#!/usr/bin/env python3
"""Second resolution criterion: the basis must resolve the prime-power support.

The Weil functional's prime piece is supported on t = +/- log(p^k) for prime powers
p^k <= c. The Galerkin basis is Fourier modes e_k(t) = exp(2*pi*i*k*t/(2 log c)) for
|k| <= N, i.e. a uniform grid on the interval (-log c, log c) of length 2 log c with
2N+1 degrees of freedom. Its resolution in t is

    delta_t = 2 log(c) / (2N + 1).

Two prime powers closer together than delta_t are not separated by the basis. So the
basis resolves the prime comb only while

    delta_t <= min_gap(c) = min over adjacent prime powers of log(p^k / q^j)

equivalently

    N >= log(c)/min_gap(c) - 1/2.

This is INDEPENDENT of the archimedean truncation criterion in resolution.py. That one
is about T; this one is about N alone.
"""
import math


def _primes_up_to(n):
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for i in range(2, int(n ** 0.5) + 1):
        if sieve[i]:
            sieve[i*i::i] = [False] * len(sieve[i*i::i])
    return [i for i, v in enumerate(sieve) if v]


def prime_powers_up_to(c):
    out = []
    for p in _primes_up_to(int(c)):
        q = p
        while q <= c:
            out.append(q)
            q *= p
    return sorted(out)


def min_gap(c):
    pp = prime_powers_up_to(c)
    logs = [math.log(x) for x in pp]
    return min(logs[i + 1] - logs[i] for i in range(len(logs) - 1)), pp


def n_required(c):
    g, _ = min_gap(c)
    return math.log(c) / g - 0.5


def plan_N(c, eta=4.0):
    """N giving a fixed headroom eta over the prime-spacing requirement."""
    return int(math.ceil(eta * n_required(c)))


if __name__ == "__main__":
    print(f"{'c':>3} {'log c':>7} {'#p^k':>5} {'tightest pair':>16} {'min_gap':>8} {'N_req':>7} "
          f"{'N=36 ok':>8} {'N=52 ok':>8}")
    for c in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        g, pp = min_gap(c)
        logs = [math.log(x) for x in pp]
        i = min(range(len(logs) - 1), key=lambda j: logs[j + 1] - logs[j])
        nreq = n_required(c)
        print(f"{c:>3} {math.log(c):7.4f} {len(pp):>5} {str(pp[i])+'/'+str(pp[i+1]):>16} "
              f"{g:8.4f} {nreq:7.1f} {'yes' if 36>=nreq else 'NO':>8} {'yes' if 52>=nreq else 'NO':>8}")
