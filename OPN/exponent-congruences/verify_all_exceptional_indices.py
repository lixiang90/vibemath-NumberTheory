"""Independent check using the whole forcing closure for each possible Euler prime."""

import json
import sys
from functools import cache
from math import isqrt
from pathlib import Path

from flint import fmpz

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from verify_exceptional_indices import BOUNDS


@cache
def terms(p, exponent):
    denominator = fmpz(p)**exponent
    return (denominator*p-1)//(p-1), denominator


def abundant(vertices, exception, index, euler):
    numerator, denominator = terms(3, 8)
    for p in sorted(vertices):
        exponent = 2
        if p == exception:
            exponent = min(index-1, 8)
        elif p == euler:
            exponent = 1
        n, d = terms(p, exponent)
        numerator *= n
        denominator *= d
        if numerator > 2*denominator:
            return True
    return False


def closure(roots, exception, euler, ordinary_edges, euler_edges):
    seen = set(roots)
    pending = list(roots)
    while pending:
        p = pending.pop()
        if p == exception:
            continue
        adjacent = euler_edges[p] if p == euler else ordinary_edges[p]
        for t in adjacent:
            if t not in seen:
                seen.add(t)
                pending.append(t)
    return seen


def verify():
    data = json.loads((HERE / "exceptional-index-primes.json").read_text(encoding="ascii"))
    vertices = tuple(map(int, data["primes"]))
    if len(vertices) != len(set(vertices)) or not all(p > 3 and fmpz(p).is_prime() for p in vertices):
        raise ValueError("invalid prime universe")
    ordinary = {p: [t for t in vertices if (p*p+p+1) % t == 0] for p in vertices}
    special = {p: [t for t in vertices if (p+1) % t == 0] for p in vertices}
    choices = [0]+[p for p in vertices if p % 4 == 1]
    assignments = 0
    for e, upper in BOUNDS.items():
        seeds = {}
        for d in range(5, upper, 2):
            if any(d % t == 0 for t in range(3, isqrt(d)+1, 2)):
                continue
            value = (e**d-1)//(e-1)
            seeds[d] = {p for p in vertices if value % p == 0}
        for p in choices:
            if p == e:
                continue
            base = closure({e, 13, p}-{0}, e, p, ordinary, special)
            if not abundant(base, e, 5, p):
                for d, roots in seeds.items():
                    implied = base | closure(roots, e, p, ordinary, special)
                    if not abundant(implied, e, d, p):
                        raise ValueError(f"open assignment for E={e}, index={d}, Euler={p}")
            assignments += len(seeds)
    print(f"{assignments} Euler assignments verified on {len(vertices)} prime vertices.")


if __name__ == "__main__":
    verify()
