"""Exact local computations for X^3 - X - 1 and the level-23 weight-one form."""

from collections import Counter
from fractions import Fraction
from math import gcd, isqrt


def primes_up_to(bound):
    sieve = bytearray([1]) * (bound + 1)
    sieve[:2] = bytes(2)
    for p in range(2, isqrt(bound) + 1):
        if sieve[p]:
            sieve[p * p::p] = bytes((bound - p * p) // p + 1)
    return [p for p in range(2, bound + 1) if sieve[p]]


def eta_coefficients(bound):
    """Coefficients of q E(q) E(q^23), E(q) = product (1-q^n)."""
    terms = [(0, 1)]
    k = 1
    while k * (3 * k - 1) // 2 < bound:
        for exponent in (k * (3 * k - 1) // 2, k * (3 * k + 1) // 2):
            if exponent < bound:
                terms.append((exponent, (-1) ** k))
        k += 1
    terms.sort()
    result = [0] * (bound + 1)
    for j, second in terms:
        shift = 1 + 23 * j
        if shift > bound:
            break
        for i, first in terms:
            if shift + i > bound:
                break
            result[shift + i] += first * second
    return result


def multiply_mod(f, g, p):
    """Multiply triples in F_p[X]/(X^3-X-1)."""
    a, b, c = f
    d, e, h = g
    cubic = b * h + c * e
    quartic = c * h
    return ((a * d + cubic) % p, (a * e + b * d + cubic + quartic) % p,
            (a * h + b * e + c * d + quartic) % p)


def trace_at_prime(p):
    if p == 23:
        return 1
    # Quadratic reciprocity gives chi_{-23}(p) = (p/23), including p=2.
    if pow(p, 11, 23) != 1:
        return 0
    # For an even Frobenius the cubic is either split or irreducible.
    exponent = p
    x = (0, 1, 0)
    power, result = x, (1, 0, 0)
    while exponent:
        if exponent % 2:
            result = multiply_mod(result, power, p)
        power = multiply_mod(power, power, p)
        exponent //= 2
    return 2 if result == x else -1


def euler_coefficients(bound, primes):
    """Compute a(n) multiplicatively from every local Euler polynomial."""
    coeff = [0] * (bound + 1)
    coeff[1] = 1
    for p in primes:
        trace = trace_at_prime(p)
        det = 0 if p == 23 else (1 if pow(p, 11, 23) == 1 else -1)
        pk = p
        previous, current = 1, trace
        while pk <= bound:
            for m in range(1, bound // pk + 1):
                if m % p and coeff[m]:
                    coeff[m * pk] = coeff[m] * current
            pk *= p
            previous, current = current, trace * current - det * previous
    return coeff


def roots(p):
    return [r for r in range(p) if (r**3 - r - 1) % p == 0]


def main():
    bound = 100_000
    primes = primes_up_to(bound)
    eta = eta_coefficients(bound)
    euler = euler_coefficients(bound, primes)
    if eta != euler:
        n = next(i for i in range(bound + 1) if eta[i] != euler[i])
        raise ArithmeticError(f"Euler and eta coefficients disagree at {n}")
    for p in primes_up_to(1000):
        if eta[p] != len(roots(p)) - 1:
            raise ArithmeticError(f"Direct root count disagrees at {p}")
    counts = Counter(eta[p] for p in primes if p != 23)
    print(f"Euler and eta coefficients agree for 1 <= n <= {bound}")
    print("p <= 100000, unramified prime counts", dict(sorted(counts.items())))
    print("p, residue modulo 23, roots, trace")
    for p in (2, 5, 13, 23, 59):
        print(p, p % 23, roots(p), eta[p])
    for modulus in (1, 23, 69):
        rows = {}
        for p in primes:
            if p == 23 or gcd(p, modulus) != 1:
                continue
            row = rows.setdefault(p % modulus, [0, 0, 0])
            row[0] += 1
            row[1] += eta[p]
            row[2] += eta[p] ** 2
        total = sum(n for n, _, _ in rows.values())
        residual = sum(squares - Fraction(s * s, n) for n, s, squares in rows.values())
        print(f"m={modulus}, finite best mean square = {float(residual / total):.8f}")


if __name__ == "__main__":
    main()
