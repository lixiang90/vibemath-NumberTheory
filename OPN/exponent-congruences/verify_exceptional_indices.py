"""Exclude every small prime divisor of the exceptional even exponent plus one."""

import heapq
import json
import sys
from fractions import Fraction
from functools import cache
from math import isqrt
from pathlib import Path


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_one_exception import prime_sieve
from verify_exponent import require, verify_primes


# For each exceptional prime, every prime index from 5 up to this bound is excluded,
# with the upper endpoint omitted.
BOUNDS = {
    7: 263, 13: 607, 19: 263, 31: 263, 61: 181, 97: 263, 127: 193, 331: 193,
    3169: 131, 5419: 131, 1674289: 131, 3348577: 109, 934415109937: 131,
    3737657091169: 97,
}


@cache
def abundance(p, exponent):
    return Fraction(p**(exponent+1)-1, (p-1)*p**exponent)


def check_vertices(data):
    vertices = tuple(map(int, data["primes"]))
    require(vertices == tuple(sorted(set(vertices))) and vertices[0] >= 5, "invalid vertex list")
    small = [p for p in vertices if p <= 10**12]
    divisors = prime_sieve(isqrt(max(small)))
    for p in small:
        for r in divisors:
            if r*r > p:
                break
            require(p % r != 0, f"composite vertex {p}")
    verify_primes((p for p in vertices if p > 10**12), data)
    require(set(BOUNDS) <= set(vertices), "missing exceptional primes")
    return vertices


def index_roots(base, index, vertices):
    # A prime dividing base-1 divides the geometric sum exactly when it equals index.
    return {p for p in vertices if (p == index and base % p == 1)
            or ((p-1) % index == 0 and base % p != 1 and pow(base, index, p) == 1)}


def close(exception, index, euler, roots, edges):
    known = roots | {exception, 13}
    if euler:
        known.add(euler)
    pending = list(known)
    heapq.heapify(pending)
    ordinary = []
    value = abundance(3, 8)
    while pending:
        p = heapq.heappop(pending)
        if p == exception:
            value *= abundance(p, min(index-1, 8))
            targets = ()
        elif p == euler:
            value *= abundance(p, 1)
            targets = edges[p][1]
        else:
            value *= abundance(p, 2)
            targets = edges[p][0]
            ordinary.append(p)
        if value > 2:
            return ordinary
        for r in targets:
            if r not in known:
                known.add(r)
                heapq.heappush(pending, r)
    raise ValueError(f"unclosed case for E={exception}, index={index}, Euler={euler}")


def verify():
    data = json.loads((HERE / "exceptional-index-primes.json").read_text(encoding="ascii"))
    vertices = check_vertices(data)
    edges = {p: ([r for r in vertices if (p*p+p+1) % r == 0],
                 [r for r in vertices if (p+1) % r == 0]) for p in vertices}
    indices = prime_sieve(max(BOUNDS.values()))
    count = branches = 0
    for e, upper in BOUNDS.items():
        require(upper in indices, "nonprime published endpoint")
        for d in (p for p in indices if 5 <= p < upper):
            roots = index_roots(e, d, vertices)
            require(bool(roots), f"missing index root for {e}, {d}")
            ordinary = close(e, d, 0, roots, edges)
            choices = [p for p in ordinary if p % 4 == 1]
            for p in choices:
                close(e, d, p, roots, edges)
            count += 1
            branches += 1+len(choices)
    print(f"{count} prime-index exclusions, {branches} Euler branches, {len(vertices)} vertices.")
    print(f"{len(data['elliptic'])} elliptic steps and {len(data['lucas'])} Lucas certificates.")


if __name__ == "__main__":
    verify()
