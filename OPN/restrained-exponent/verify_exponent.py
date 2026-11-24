"""Replay the 899-index exclusion theorem using only the Python standard library."""

import argparse
import json
from fractions import Fraction
from math import gcd, isqrt, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def trial_prime(n):
    return n >= 2 and all(n % k for k in range(2, isqrt(n) + 1))


def add(p, q, a, n):
    """Affine addition. Every division must be valid at every prime dividing n."""
    if p is None:
        return q
    if q is None:
        return p
    x, y = p
    u, v = q
    if x == u:
        if (y + v) % n == 0:
            return None
        require(y == v, f"ambiguous point addition modulo {n}")
        slope = (3 * x * x + a) * pow(2 * y, -1, n) % n
    else:
        slope = (v - y) * pow(u - x, -1, n) % n
    new_x = (slope * slope - x - u) % n
    return new_x, (slope * (x - new_x) - y) % n


def multiply(k, point, a, n):
    result = None
    while k:
        if k % 2:
            result = add(result, point, a, n)
        k //= 2
        if k:
            point = add(point, point, a, n)
    return result


def verify_primes(nodes, proofs):
    elliptic = proofs["elliptic"]
    lucas = proofs["lucas"]
    require(not (elliptic.keys() & lucas.keys()), "duplicate primality proof")
    seen = set()

    def prove(n):
        if n in seen:
            return
        if n < 10_000:
            require(trial_prime(n), f"composite leaf {n}")
            return
        if str(n) in elliptic:
            trace, s, a, x, y = map(int, elliptic[str(n)])
            m = n + 1 - trace
            require(trace * trace < 4 * n and s > 0 and m % s == 0, f"bad curve order for {n}")
            q = m // s
            # This integer upper bound for (n**(1/4) + 1)**2 suffices for Hasse.
            require((isqrt(isqrt(n)) + 2) ** 2 < q < n, f"elliptic size bound for {n}")
            prove(q)
            x, y = x % n, y % n
            b = (y * y - x**3 - a * x) % n
            require(gcd(6 * (4 * a**3 + 27 * b**2), n) == 1, f"singular curve for {n}")
            point = multiply(s, (x, y), a, n)
            require(point is not None and multiply(q, point, a, n) is None, f"bad point for {n}")
        else:
            base, recorded = lucas[str(n)]
            factors = [(int(q), e) for q, e in recorded]
            require(len({q for q, _ in factors}) == len(factors), f"duplicate factor for {n}")
            require(all(2 <= q < n and e > 0 for q, e in factors), f"invalid factors for {n}")
            require(prod(q**e for q, e in factors) == n - 1, f"incomplete factorization for {n}")
            require(pow(base, n - 1, n) == 1, f"bad Lucas base for {n}")
            for q, _ in factors:
                prove(q)
                require(gcd(pow(base, (n - 1) // q, n) - 1, n) == 1, f"bad Lucas order for {n}")
        seen.add(n)

    for n in nodes:
        prove(n)
    require({str(n) for n in seen} == elliptic.keys() | lucas.keys(), "unused primality proof")


def verify(chains, proofs):
    # The inputs to the theorem are fixed here, independently of the search files.
    excluded_small = {5, 7, 11, 13}
    expected = {n for n in range(3, 900, 2) if trial_prime(n)} - {857}
    roots = {int(d): int(r) for d, r in chains["roots"].items()}
    extras = {int(d): int(r) for d, r in chains["extra_indices"].items()}
    require(set(roots) == expected and set(extras) == {857**2}, "incomplete index coverage")
    require(chains["unresolved"] == [857] and chains["inspected_through"] == 899,
            "wrong claimed range")
    require(chains["terminal_bound"] == max(excluded_small), "wrong base exclusions")
    steps = {int(r): row for r, row in chains["steps"].items()}
    verify_primes(steps, proofs)

    visiting = set()
    closed = set()

    def exclude(r):
        if r in excluded_small or r in closed:
            return
        require(r > 13 and r in steps and r not in visiting, f"open or cyclic node {r}")
        visiting.add(r)
        row = steps[r]
        branches = {"square", "euler"} if r % 12 == 1 else {"square"}
        require(set(row) == branches, f"missing or surplus role for {r}")
        for role in branches:
            child = int(row[role])
            required = r * r + r + 1 if role == "square" else r + 1
            require(child >= 5 and required % child == 0, f"bad {role} edge {r} -> {child}")
            exclude(child)
        visiting.remove(r)
        closed.add(r)

    for d, r in (roots | extras).items():
        require(r >= 5 and pow(3, d, r) == 1, f"bad initial divisor for index {d}")
        if d in extras:
            require(pow(3, 857, r) != 1, "the square-index witness has smaller order")
        exclude(r)
    require(closed == steps.keys(), "unused exclusion nodes")
    density = prod(Fraction(d - 1, d) for d in roots) * Fraction(857**2 - 1, 857**2)
    return {"prime_indices": len(roots), "additional_index": 857**2, "internal_nodes": len(closed),
            "elliptic_steps": len(proofs["elliptic"]), "lucas_proofs": len(proofs["lucas"]),
            "least_unexcluded_prime": 857, "density_floor_scaled_1e12": 10**12 * density // 1}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    directory = parser.parse_args().directory
    chains = json.loads((directory / "index-chains.json").read_text(encoding="utf-8"))
    proofs = json.loads((directory / "primality.json").read_text(encoding="utf-8"))
    print(json.dumps(verify(chains, proofs), indent=2))
