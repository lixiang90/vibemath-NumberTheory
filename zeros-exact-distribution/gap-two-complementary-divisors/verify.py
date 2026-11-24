"""Finite arithmetic checks for the complementary-divisor formula at gap two.

Run with Python 3.9 or later. Only the standard library is used.
"""

from array import array
from fractions import Fraction
from math import exp, fsum, isqrt, log, prod


def smallest_prime_factors(limit):
    spf = array("I", range(limit + 1))
    for q in range(2, isqrt(limit) + 1):
        if spf[q] == q:
            for n in range(q * q, limit + 1, q):
                if spf[n] == n:
                    spf[n] = q
    return spf


def squarefree_factors(n, spf):
    factors = []
    while n > 1:
        q = spf[n]
        n //= q
        if n % q == 0:
            return None
        factors.append(q)
    return factors


def signed_divisors(factors):
    divisors = [(1, 1)]
    for q in factors:
        divisors += [(d * q, -sign) for d, sign in divisors]
    return divisors


def exact_ratio(m, factors):
    """Return exp(2 T(m)) as an exact rational."""
    return prod(
        (Fraction(m, d * d) ** sign
         for d, sign in signed_divisors(factors) if d * d < m),
        start=Fraction(1),
    )


def check_pointwise(spf, limit):
    checked = 0
    for m in range(3, limit + 1, 2):
        factors = squarefree_factors(m, spf)
        if factors is None:
            continue
        value = exact_ratio(m, factors)
        r = len(factors)
        if r % 2 == 1:
            assert value == (m if r == 1 else 1), (m, value)
        if r == 2:
            assert value == factors[0] ** 2, (m, value)
        if r == 4:
            a, b, c, d = factors
            magnitude = max(1, min(a * a, Fraction(a * b * c, d)))
            assert value == Fraction(1, 1) / magnitude**2, (m, value)
        if r >= 3 and factors[-1] ** 2 > m:
            assert value == 1, (m, value)
        checked += 1
    return checked


def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def check_six_factors():
    cases = [
        ([5, 7, 11, 13, 17, 23], Fraction(13731656640625, 279841)),
        ([3, 5, 7, 11, 13, 14009], Fraction(196252081, 225450225)),
    ]
    for factors, expected in cases:
        m = prod(factors)
        assert all(is_prime(q) for q in factors) and is_prime(m - 2)
        assert exact_ratio(m, factors) == expected
        print(f"p={m - 2}  factors={factors}  exp(2T)={expected}")


def bump(u):
    if 1 < u < 2:
        return exp(4 - 1 / ((u - 1) * (2 - u)))
    return 0.0


def bump_integral(n):
    # Composite Simpson rule. The values at both endpoints are zero.
    return fsum((4 if j % 2 else 2) * bump(1 + j / n)
                for j in range(1, n)) / (3 * n)


def check_window(x, spf, integral):
    z = exp(log(x) ** (1 / 8))
    twin_count = 0
    twin_weights = []
    all_weights = []
    even_weights = {2: [], 4: [], 6: []}
    for p in range(x + 1, 2 * x):
        if spf[p] != p:
            continue
        m = p + 2
        weight = log(p) * bump(p / x) / integral
        if spf[m] == m:
            twin_count += 1
            twin_weights.append(weight * log(m))
        if spf[m] <= z:
            continue
        factors = squarefree_factors(m, spf)
        if factors is None:
            continue
        t = fsum(sign * log(m / (d * d)) / 2
                 for d, sign in signed_divisors(factors) if d * d < m)
        all_weights.append(weight * t)
        if len(factors) % 2 == 0:
            even_weights[min(len(factors), 6)].append(weight * t)
    mass = fsum(twin_weights) / x
    total = fsum(all_weights) / x
    even = [fsum(weights) / x for weights in even_weights.values()]
    error = abs(mass - 2 * (total - fsum(even)))
    assert error < 1e-12 and even[1] <= 1e-14, (x, error, even)
    print(f"{x:>8} {twin_count:>6} {mass:.7f} {total:.7f} "
          f"{even[0]:.7f} {even[1]:.7f} {even[2]:.3g} {error:.2g}")


if __name__ == "__main__":
    windows = [1_000, 10_000, 100_000, 1_000_000]
    spf = smallest_prime_factors(2 * max(windows) + 2)
    print("Odd squarefree integers checked:", check_pointwise(spf, 100_000))
    check_six_factors()
    integral = bump_integral(16384)
    assert abs(integral - bump_integral(8192)) < 1e-14
    print("       X  count       M/X     A_sf/X       E2/X       E4/X     E6+/X   residual")
    for x in windows:
        check_window(x, spf, integral)
