"""Rebuild the divisibility trees for two exceptional non-Euler exponents."""

import argparse
import json
import sys
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import ceil, isqrt
from pathlib import Path

from flint import fmpz


LIMIT = 853
SCALE = 10**15
FREE_OF_THREE = {5, 7, 13, 71, 103, 541}
WITH_THREE = {541}
REPUNITS = {d: (3**d - 1) // 2 for d in (71, 103, 541)}
SECOND_EXCEPTIONS = {7, 13, 19, 31, 61, 97, 127, 331, 3169, 5419,
                     1674289, 3348577, 934415109937, 3737657091169}


def small_primes(limit):
    return [p for p in range(5, limit + 1, 2) if all(p % k for k in range(2, isqrt(p) + 1))]


@cache
def lower(p, e):
    return SCALE * (p**(e + 1) - 1) // ((p - 1) * p**e)


def product(x, y):
    return x * y // SCALE


def required_checks():
    primes = set(small_primes(LIMIT))
    if not (FREE_OF_THREE | WITH_THREE) <= primes:
        raise ArithmeticError("nonprime endpoint in the tables")
    first = (primes - FREE_OF_THREE) | {25, 49, 169, 71**2, 541**2, 103**4}
    first.update(p*q for p, q in combinations(sorted(FREE_OF_THREE), 2))
    second = {27, 3*541**2} | {3*p for p in primes - WITH_THREE} | {9*p for p in WITH_THREE}
    return {1: first, 2: second}


def check_smooth_indices(checks):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "restrained-exponent"))
    from smooth_products import surviving_products

    primes = small_primes(LIMIT)
    first = {n for n, _ in surviving_products(primes, checks[1])}
    second = {3*n for n, _ in surviving_products([3, *primes], {d//3 for d in checks[2]})}
    if first != {1, 5, 7, 13, 71, 103, 103**2, 103**3, 541} or second != {3, 9, 1623}:
        raise ArithmeticError("unexpected list of smooth parts")
    # The initial 1 and the five exact small indices are excluded in the text.
    left = (first | second) - {1, 3, 5, 7, 9, 13}
    left = {n for n in left if n % 3 == 0 or n-1 >= 1150002 or all(n % d for d in REPUNITS)}
    if left != {1623}:
        raise ArithmeticError("unexpected smooth index after the support bound")
    print(f"remaining smooth index: {sorted(left)}")


def verify():
    path = Path(__file__).with_name("two-index-primes.json")
    vertices = tuple(map(int, json.loads(path.read_text(encoding="ascii"))))
    if vertices != tuple(sorted(set(vertices))) or min(vertices) < 5:
        raise ArithmeticError("invalid vertex list")
    if not all(fmpz(p).is_prime() for p in vertices):
        raise ArithmeticError("composite vertex")

    @cache
    def adjacent(p, euler):
        n = p + 1 if euler else p*p + p + 1
        return frozenset(t for t in vertices if n % t == 0)

    def exclude(index, budget, seeds=None, forbidden=frozenset(), fixed_exception=0):
        roots = (frozenset(p for p in vertices if pow(3, index, p) == 1)
                 if seeds is None else frozenset(seeds))
        if not roots:
            raise ArithmeticError(f"no initial factors for index {index}")
        nodes = 0

        def visit(known, pending, exceptions, has_euler, abundance):
            nonlocal nodes
            nodes += 1
            bound = abundance
            for p in sorted(pending):
                bound = product(bound, lower(p, 1))
            if bound > 2*SCALE:
                return
            if not pending:
                raise ArithmeticError(f"unresolved index {index} with budget {budget}")
            p = min(pending)
            rest = pending - {p}
            if p == fixed_exception:
                visit(known, rest, exceptions, has_euler, product(abundance, lower(p, 4)))
                return
            ordinary = adjacent(p, False) - known
            visit(known | ordinary, rest | ordinary, exceptions, has_euler,
                  product(abundance, lower(p, 2)))
            if exceptions < budget and p not in forbidden:
                visit(known, rest, exceptions + 1, has_euler, product(abundance, lower(p, 4)))
            if not has_euler and p % 4 == 1:
                forced = adjacent(p, True) - known
                visit(known | forced, rest | forced, exceptions, True,
                      product(abundance, lower(p, 1)))

        visit(roots, roots, 0, False, lower(3, min(index - 1, 18)))
        return nodes

    print(f"{len(vertices)} verified primes")
    checks = required_checks()
    for budget, indices in checks.items():
        total = sum(exclude(d, budget) for d in sorted(indices))
        print(f"outside exceptions <= {budget}: {len(indices)} divisors, {total} states")
    for index, r in REPUNITS.items():
        if r not in vertices:
            raise ArithmeticError("missing prime repunit")
        special = exclude(index, 1, forbidden={r})
        special += sum(exclude(index, 0, seeds=(p,), fixed_exception=r) for p in (7, 31))
        print(f"index {index} and the two induction bases: {special} states")
    if not SECOND_EXCEPTIONS <= set(vertices) or REPUNITS[541] in SECOND_EXCEPTIONS:
        raise ArithmeticError("invalid list of second exceptions")
    x = 3**541
    large = (x*x + x + 1) // 13
    if large not in vertices or 13*REPUNITS[541]*large != (3**1623-1)//2:
        raise ArithmeticError("incorrect initial factorization at index 1623")
    cases = exclude(1623, 2, forbidden={REPUNITS[541]})
    cases += exclude(1623, 2, forbidden=SECOND_EXCEPTIONS)
    print(f"both exception locations at index 1623: {cases} states")
    check_smooth_indices(checks)


def weighted_support_bound():
    import verify_ternary_support as old

    k, y, end = 1150000, 3*10**8, 6000
    count, total = old.PREFIX[y, 1]
    s = (Fraction(total, old.SCALE) + Fraction(k-count, y-1) + Fraction(1, old.P-1))
    correction = (1-Fraction(1, y))*Fraction(total-count, old.SCALE) - Fraction(count, y)
    h = old.log_interval(y)[1]*s - correction
    upper = s + Fraction(old.PREFIX[10**10, 5][1], old.SCALE) + Fraction(1, 10**6)
    for (a, i, j), b in zip(old.LOW, [r[0] for r in old.LOW[1:]]+[end], strict=True):
        u = 10**a
        density, error = old.sieve_constants(i, j)
        upper += (density*h*Fraction(u, u-1)*old.log_interval(Fraction(b, a))[1]
                  + Fraction(11*error*k, 5*(u-1)))
    if (old.P**7 <= 675**16 or end*old.log_interval(10)[0]/8 <= 16+Fraction(k,674)
            or correction <= 0 or upper >= old.log_interval(Fraction(16,15))[0]):
        raise ArithmeticError("weighted support bound failed")
    print(f"K={k}, log upper <= {ceil(upper*old.SCALE)}/{old.SCALE}")


def check_tree_transfer(path):
    from verify_role_tree import verify_tree
    from verify_ternary_support import HIGH, support_bound

    print(verify_tree(path))
    for line in path.read_text(encoding="ascii").splitlines():
        if "^" in line and int(line.split("^")[0]) in REPUNITS.values():
            raise ArithmeticError("exceptional prime is an internal node")
    if any(r <= 3*10**11 or r % 3 != 1 for r in REPUNITS.values()):
        raise ArithmeticError("exceptional prime in the small target range")
    weighted_support_bound()
    print(support_bound(10**9, 3*10**11, 6000000, HIGH, Fraction(64,51)))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tree", type=Path, help="also check the three prime repunits and the support estimates")
    args = parser.parse_args()
    verify()
    if args.tree:
        check_tree_transfer(args.tree)
