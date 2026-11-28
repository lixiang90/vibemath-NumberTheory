"""Certify the forbidden-prime table for 3-free odd perfect numbers."""

import sys
from fractions import Fraction as F
from functools import cache
from math import isqrt, prod
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cyclotomic-pairs"))
from search_cyclotomic_pairs import factorizations, primes_through

X = 10**7
Y = 10**9
SCALE = 10**15
ROUND = 10**6
# Forbidden primes, excluded total exceptional budget K (including Euler),
# certified upward bounds for U and D in millionths.
ROWS = (
    ((5,), 9, 1936425, 11516),
    ((7,), 6, 1942639, 891),
    ((11,), 6, 1958320, 1628),
    ((13, 19), 6, 1981963, 857),
    ((7, 11), 9, 1934684, 6257),
    ((7, 13), 8, 1964449, 2331),
    ((11, 13), 7, 1976318, 2286),
    ((5, 7), 13, 1852217, 43941),
    ((5, 11), 12, 1889507, 52857),
    ((5, 7, 11), 15, 1741695, 103283),
)


def require(ok, reason):
    if not ok:
        raise ArithmeticError(reason)


def cap(p):
    return F(p, p - 1)


def sixth_root_lower(t):
    lo, hi = 1000, 1000 * t
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**6 <= t * 1000**6:
            lo = mid
        else:
            hi = mid
    return F(lo, 1000)


def finite_sum(bases, x, y):
    """Integer upper sum for vectors with m <= y^2, in units 1/SCALE."""
    y2 = y * y

    def visit(i, m):
        if i == len(bases):
            return SCALE // max(x, isqrt(m)) + 1, 1
        total = count = 0
        while m <= y2:
            contribution, n = visit(i + 1, m)
            total += contribution
            count += n
            m *= bases[i]
        return total, count

    return visit(0, 1)


@cache
def tail(bases):
    finite, count = finite_sum(bases, X, Y)
    # Y^(-2/3) = 10^(-6).
    far = prod(z / (z - 1) for z in map(sixth_root_lower, bases)) / 10**6
    return F(len(bases), X) + F(X + 2, X) * (F(finite, SCALE) + far), count


def weights(limit):
    w, label = {}, {}
    for p, factors in factorizations(limit).items():
        if p <= 3:
            continue
        r = p if p % 3 == 1 else max(t for t, _ in factors)
        w[r] = w.get(r, F(1)) * cap(p)
        if any(p in forbidden for forbidden, *_ in ROWS):
            label[p] = r
    return w, label


def self_test():
    # Independent direct enumeration, including the vector (0,0).
    values = [7**a * 13**b for a in range(10) for b in range(10)
              if 7**a * 13**b <= 40**2]
    require(finite_sum((7, 13), 10, 40) ==
            (sum(SCALE // max(10, isqrt(m)) + 1 for m in values), len(values)),
            "smooth-product enumeration")
    # Compare all small cyclotomic profiles to direct trial division.
    for p, profile in factorizations(1000).items():
        m = p*p + p + 1
        expected = []
        for q in range(2, isqrt(m) + 1):
            power = 0
            while m % q == 0:
                power += 1
                m //= q
            if power:
                expected.append((q, power))
        if m > 1:
            expected.append((m, 1))
        require(list(profile) == expected, f"sieve profile for {p}")


def main():
    self_test()
    w, label = weights(X)
    # At most len(F) weights are decreased. The top K after decreasing
    # can therefore be computed from the first K+len(F) original weights.
    size = max(k + len(f) for f, k, *_ in ROWS)
    prefix = sorted(w, key=w.get, reverse=True)[:size]
    require(len(prefix) == size, "weight prefix")
    target_primes = [t for t in primes_through(200) if t % 3 == 1]
    print("F | K excluded | U < | D < | number of smooth products | top labels")
    for forbidden, k, u_num, d_num in ROWS:
        adjusted = {r: w[r] for r in prefix}
        for p in forbidden:
            if label[p] in adjusted:
                adjusted[label[p]] /= cap(p)
        order = sorted(adjusted, key=adjusted.get, reverse=True)[:k]
        upper = prod(adjusted[r] for r in order)
        bases = tuple(t for t in target_primes if t not in forbidden)[:k]
        require(len(bases) == k, "target-prime list exhausted")
        delta, n = tail(bases)
        require(upper < F(u_num, ROUND), "finite cap")
        require(delta < F(d_num, ROUND), "tail cap")
        require(u_num + 2*d_num < 2*ROUND, "no strict abundance separation")
        print(forbidden, k, F(u_num, ROUND), F(d_num, ROUND), n, order)
    print(f"Verified all {len(ROWS)} forbidden-prime bounds.")


if __name__ == "__main__":
    main()
