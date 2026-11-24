"""Exact enumeration of smooth integers avoiding specified divisors."""

from collections import defaultdict


def smooth_factorization(value, primes):
    if value < 2:
        raise ValueError("a forbidden index must exceed 1")
    factors = {}
    for p in primes:
        if value % p != 0:
            continue
        exponent = 0
        while value % p == 0:
            value //= p
            exponent += 1
        factors[p] = exponent
        if value == 1:
            break
    if value != 1:
        raise ValueError("an index has a prime factor outside the claimed window")
    return factors


def surviving_products(primes, obstructions):
    """Return (product, largest prime) pairs, including (1, 1)."""
    indices = set(obstructions)
    allowed = [p for p in primes if p not in indices]
    caps = {}
    mixed = defaultdict(list)
    for d in sorted(indices):
        if d in primes:
            continue
        factors = smooth_factorization(d, primes)
        if any(p in indices for p in factors):
            continue
        if len(factors) == 1:
            p, j = next(iter(factors.items()))
            caps[p] = min(caps.get(p, j - 1), j - 1)
        else:
            mixed[max(factors)].append(d)
    if set(caps) != set(allowed):
        raise ValueError("missing valuation cap for a surviving prime")

    def visit(start, product, largest):
        yield product, largest
        for i in range(start, len(allowed)):
            p = allowed[i]
            candidate = product
            for _ in range(caps[p]):
                candidate *= p
                if any(candidate % d == 0 for d in mixed[p]):
                    break
                yield from visit(i + 1, candidate, p)

    return visit(0, 1, 1)
