"""Check the two triangular-window identities and local obstructions by enumeration."""

from collections import Counter
from math import isqrt


def least_factors(limit):
    factors = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if factors[p] == p:
            for n in range(p * p, limit + 1, p):
                if factors[n] == n:
                    factors[n] = p
    return factors


def check(x, h):
    m = 4 * x
    factors = least_factors(m + 2 * h)
    chen = []
    twins = []
    for p in range(x + 1, 2 * x + 1):
        if factors[p] != p:
            continue
        n = p + 2
        s = factors[n]
        if s == n:
            twins.append(p)
        elif factors[n // s] != n // s or s**11 <= p**3:
            continue
        chen.append(p)
    chen_set = set(chen)
    assert all(p + 2 not in chen_set for p in chen)

    primes = [int(n >= 2 and s == n) for n, s in enumerate(factors)]
    prefix = [0]
    for flag in primes:
        prefix.append(prefix[-1] + flag)

    def count_after(y):
        return prefix[y + h + 1] - prefix[y + 1]

    for name, bases in (("Chen", chen), ("twin", twins)):
        base_set = set(bases)
        represented = Counter()
        for shift in range(1, 2 * h):
            n = m + shift
            represented[n] = sum(primes[n - p] for p in bases)
            if n % 2 or n % 3 == 2:
                assert represented[n] == 0
        triangular = sum((h - abs(n - m - h)) * r for n, r in represented.items())
        occupancy_sum = overlap = second_moment = 0
        # On (j-1,j) the first interval contains j,...,j+h-1 and
        # the reflected interval contains M-j+1,...,M-j+h.
        for j in range(x + 2 - h, 2 * x + 1):
            occupancy = sum(n in base_set for n in range(j, j + h))
            occupancy_sum += occupancy
            overlap += occupancy * count_after(m - j)
            second_moment += occupancy**2
        pair_moment = h * len(bases) + 2 * sum(
            (h - d) * sum(p + d in base_set for p in bases) for d in range(1, h)
        )
        assert occupancy_sum == h * len(bases)
        assert triangular == overlap and second_moment == pair_moment
        print(x, h, name, len(bases), triangular, second_moment, sum(r > 0 for r in represented.values()))


if __name__ == "__main__":
    for parameters in ((1000, 19), (10000, 31), (30000, 45)):
        check(*parameters)
