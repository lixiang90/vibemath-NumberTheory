"""Exclude v_3(N) = 2 when at most one non-Euler index is prime to 3."""

import sys
from collections import defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations
from math import comb, prod
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cyclotomic-pairs"))
from search_cyclotomic_pairs import factorizations


B = 1_000
X = 1_000_000
SCALE = 10**12
F = Fraction
FAR_TAIL = F(1, 200)
TOTAL_TAIL = F(3, 200)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


@cache
def upper(p):
    return F(p, p - 1)


@cache
def abundancy(p, exponent):
    return F(p**(exponent + 1) - 1, p**exponent * (p - 1))


def v_three(n):
    exponent = 0
    while n % 3 == 0:
        n //= 3
        exponent += 1
    return exponent


def verify_tail(profiles):
    far = F(comb(19, 4) - 1, X)
    far += 2 * (1 + F(2, X)) * comb(19, 3) * F(2, 5)**16
    require(far < FAR_TAIL, "far smooth tail")
    weights = defaultdict(int)
    for p, factors in profiles.items():
        if B < p and p % 6 == 5:
            label = max(t for t, exponent in factors)
            weights[label] += (SCALE + p - 2) // (p - 1)
    peak = max(weights, key=weights.get)
    require(F(weights[peak], SCALE) < F(13, 10_000), "middle interval")
    require(F(4, B) + 4*F(13, 10_000) + FAR_TAIL < TOTAL_TAIL, "total tail")
    print(f"Largest middle-interval weight is at {peak}: {weights[peak]}/{SCALE}")


def make_test(profiles):
    by_support = defaultdict(list)
    for p, factors in profiles.items():
        if p % 6 == 5:
            by_support[tuple(t for t, exponent in factors)].append(p)
    trial = [p for p in profiles if 5 <= p < 2_000]

    @cache
    def children(p, euler):
        if not euler and p in profiles:
            return frozenset(t for t, exponent in profiles[p] if t > 3)
        n = p + 1 if euler else p*p + p + 1
        return frozenset(t for t in trial if n % t == 0)

    @cache
    def small_preimages(targets):
        return frozenset(p for k in range(len(targets) + 1)
                         for sub in combinations(targets, k)
                         for p in by_support.get(sub, ()))

    def possible(seed, five, role13):
        def visit(pending, regular, exception, euler, used):
            known = pending | regular | {exception, euler}
            known = known - {0}
            if not five and 5 in known:
                return False
            free = 2 - used + (exception == 0) + (euler == 0)
            if sum(p % 3 == 1 for p in pending) > free:
                return False
            low = F(13, 9) * prod(abundancy(p, 2) for p in regular)
            low *= prod(abundancy(p, 1) for p in pending)
            if exception:
                low *= abundancy(exception, 4)
            if euler:
                low *= abundancy(euler, 1)
            if low > 2:
                return False

            if not pending:
                if free:
                    return True
                targets = tuple(sorted(p for p in known if p % 3 == 1))
                for p in regular:
                    value = p*p + p + 1
                    for t in (3, *targets):
                        while value % t == 0:
                            value //= t
                    if value != 1:
                        return False
                allowed = known | small_preimages(targets)
                if not five:
                    allowed = allowed - {5}
                return F(13, 9)*prod(map(upper, allowed)) >= 2*(1 - FAR_TAIL)

            p = min(pending)
            rest = pending - {p}
            role = role13 if p == 13 else None
            normal_copies = used + (p % 3 == 1)
            if role in (None, "R") and normal_copies <= 2:
                more = children(p, False) - known
                if visit(rest | more, regular | {p}, exception, euler, normal_copies):
                    return True
            if role in (None, "E") and not exception:
                if visit(rest, regular, p, euler, used):
                    return True
            euler_copies = used + v_three(p + 1)
            if role in (None, "P") and not euler and p % 4 == 1 and euler_copies <= 2:
                more = children(p, True) - known
                if visit(rest | more, regular, exception, p, euler_copies):
                    return True
            return False

        seeds = {seed, 13}
        if five:
            seeds.add(5)
        return visit(seeds, set(), 0, 0, 0)

    return possible


def main():
    profiles = factorizations(X)
    verify_tail(profiles)
    possible = make_test(profiles)
    small = [p for p in profiles if 5 <= p <= B]
    for five, role in ((True, None), (False, "R"), (False, "E"), (False, "P")):
        allowed = [p for p in small if possible(p, five, role)]
        bound = F(13, 9) * prod(map(upper, allowed)) / (1 - TOTAL_TAIL)
        require(bound < 2, f"open case for 5={five}, role of 13={role}")
        scaled = -(-10**6 * bound.numerator // bound.denominator)
        print(f"5={five}, role13={role}, small primes={allowed}, I(N) <= {scaled}/1000000")


if __name__ == "__main__":
    main()
