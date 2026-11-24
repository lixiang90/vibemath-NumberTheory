"""Verify index obstructions for odd perfect numbers with one exceptional even exponent."""

import json
import sys
from fractions import Fraction
from functools import cache
from math import isqrt, prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLDER = HERE.parent / "restrained-exponent"
sys.path.insert(0, str(OLDER))
from smooth_products import surviving_products
from verify_exponent import require, verify, verify_primes
from verify_smooth import read, verify_new_graph


def prime_sieve(limit):
    sieve = bytearray([1]) * (limit + 1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit - p*p) // p + 1)
    return [p for p in range(limit + 1) if sieve[p]]


def certified_universe():
    data = json.loads((HERE / "one-exception-primes.json").read_text(encoding="utf-8"))
    primes = tuple(map(int, data["primes"]))
    require(primes == tuple(sorted(set(primes))) and primes[0] >= 5, "invalid prime list")
    small = [p for p in primes if p <= 10**12]
    divisors = prime_sieve(isqrt(max(small)))
    for p in small:
        for t in divisors:
            if t*t > p:
                break
            require(p % t != 0, f"composite universe element {p}")
    verify_primes((p for p in primes if p > 10**12), data)
    return primes


def prove_abundance_cases(primes):
    normal = {}
    special = {}
    unit = {p: Fraction(p + 1, p) for p in primes}
    regular_gain = {}
    exceptional_gain = {}
    for p in primes:
        value = p*p + p + 1
        normal[p] = frozenset(q for q in primes if value % q == 0)
        special[p] = frozenset(q for q in primes if (p + 1) % q == 0)
        regular_gain[p] = Fraction(value, p*p) / unit[p]
        exceptional_gain[p] = Fraction(p**5 - 1, (p - 1) * p**4) / unit[p]

    @cache
    def base_factors(group):
        return prod(unit[p] for p in group)

    def prove(roots, exception_available):
        visits = 0

        def visit(known, pending, value, allow_exception, allow_euler):
            nonlocal visits
            visits += 1
            if value > 2:
                return
            require(bool(pending), f"insufficient forcing from {sorted(roots)}")
            p = min(pending)
            rest = pending - {p}
            extra = normal[p] - known
            visit(known | extra, rest | extra, value * regular_gain[p] * base_factors(extra),
                  allow_exception, allow_euler)
            if allow_exception:
                visit(known, rest, value * exceptional_gain[p], False, allow_euler)
            if allow_euler and p % 4 == 1:
                extra = special[p] - known
                visit(known | extra, rest | extra, value * base_factors(extra),
                      allow_exception, False)

        visit(roots, roots, Fraction(121, 81) * base_factors(roots), exception_available, True)
        return visits

    base_visits = sum(prove(frozenset([p]), False) for p in (7, 11, 13, 31))
    indices = [3*p for p in prime_sieve(856) if p >= 3]
    total = 0
    max_visits = (0, 0)
    for d in indices:
        roots = frozenset(p for p in primes if pow(3, d, p) == 1)
        visits = prove(roots, True)
        total += visits
        max_visits = max(max_visits, (visits, d))
    print(f"{len(primes)} certified primes, {base_visits} calls for the four initial primes")
    size, index = max_visits
    print(f"{len(indices)} indices of the form 3p, {total} calls, {size} in the largest tree (d={index})")


def transfer_old_graphs():
    base = read(OLDER, "index-chains.json")
    verify(base, read(OLDER, "primality.json"))
    for stem in ("smooth", "shared"):
        graph = read(OLDER, f"{stem}-chains.json")
        primes, forbidden = verify_new_graph(graph, read(OLDER, f"{stem}-primality.json"), base)
        base["steps"].update(graph["steps"])
        base["roots"].update(graph["roots"])
    all_roots = base["roots"] | base["extra_indices"]

    @cache
    def trace(p, euler_already_used):
        if p in (7, 11, 13):
            return
        if p == 5:
            require(euler_already_used, "5 could be the first Euler vertex")
            return
        if not euler_already_used:
            require(p % 12 != 5, f"possible omitted Euler branch at {p}")
        row = base["steps"][str(p)]
        trace(int(row["square"]), euler_already_used)
        if not euler_already_used and p % 4 == 1:
            trace(int(row["euler"]), True)

    for r in all_roots.values():
        trace(int(r), False)
    require(all(p in forbidden for p in prime_sieve(856)), "gap below 857")
    require(857**2 in forbidden, "missing square of 857")
    candidates = [n for n, _ in surviving_products(primes, forbidden) if n > 1]
    print(f"{len(all_roots)} transferred divisor obstructions")
    print(f"{len(candidates)} positive 2850-smooth candidates above 1, largest {max(candidates)}")
    require(max(candidates) == 2383**15, "smooth maximum changed")


if __name__ == "__main__":
    prove_abundance_cases(certified_universe())
    transfer_old_graphs()
