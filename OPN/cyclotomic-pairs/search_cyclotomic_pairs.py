"""Search Phi_3(x) = p**a*r, Phi_3(y) = p**b*r for prime y, integers x, a > b > 0."""

import argparse
import json
from math import isqrt


def primes_through(limit):
    sieve = bytearray(b"\x01") * (limit + 1)
    sieve[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit - p*p) // p + 1)
    return [p for p in range(2, limit + 1) if sieve[p]]


def root_classes(prime):
    if prime == 3:
        return (1,)
    if prime % 3 != 1:
        return ()
    for base in range(2, prime):
        root = pow(base, (prime - 1) // 3, prime)
        if root != 1:
            return root, prime - 1 - root
    raise ArithmeticError("cube roots of unity missing")


def factorizations(limit):
    primes = primes_through(limit)
    remainder = {y: y*y+y+1 for y in primes}
    factors = {y: [] for y in primes}
    for p in primes:
        for root in root_classes(p):
            for y in range(root, limit + 1, p):
                if y not in remainder:
                    continue
                value = remainder[y]
                exponent = 0
                while value % p == 0:
                    value //= p
                    exponent += 1
                if exponent == 0:
                    raise ArithmeticError("bad sieve residue")
                remainder[y] = value
                factors[y].append((p, exponent))
    for y, value in remainder.items():
        # Every residual prime exceeds limit, and Phi_3(y) < (limit + 1)**2.
        if value != 1:
            factors[y].append((value, 1))
    return factors


def assignments(factors):
    if len(factors) == 1:
        p, exponent = factors[0]
        if exponent >= 2:
            yield p, exponent - 1, p
    if len(factors) == 2:
        for (p, b), (r, exponent) in (factors, factors[::-1]):
            if exponent == 1:
                yield p, b, r


def search(y_limit, x_limit):
    if y_limit < 2 or x_limit < 1:
        raise ValueError("expected y_limit >= 2 and x_limit >= 1")
    bound = x_limit * x_limit + x_limit + 1
    squares_mod_63 = {v*v % 63 for v in range(63)}
    profiles = factorizations(y_limit)
    triples = tested = 0
    hits = []
    for y, factors in profiles.items():
        for p, b, r in assignments(factors):
            if p <= 3:
                # Phi_3 is odd and its 3-adic valuation is at most one.
                continue
            triples += 1
            a, value = b + 1, (y*y + y + 1) * p
            while value <= bound:
                tested += 1
                discriminant = 4 * value - 3
                if discriminant % 63 in squares_mod_63:
                    root = isqrt(discriminant)
                    if root * root == discriminant:
                        hits.append({"x": (root - 1) // 2, "y": y, "p": p, "r": r, "a": a, "b": b})
                a += 1
                value *= p
    return {"prime_y_at_most": y_limit, "integer_x_at_most": x_limit,
            "prime_y_checked": len(profiles), "assignments": triples,
            "discriminants_checked": tested, "integer_hits": hits}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--y-limit", type=int, default=10**6)
    parser.add_argument("--x-limit", type=int, default=10**100)
    args = parser.parse_args()
    print(json.dumps(search(args.y_limit, args.x_limit), indent=2))
