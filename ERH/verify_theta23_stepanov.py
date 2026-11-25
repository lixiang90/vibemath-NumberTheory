"""Exact checks for the regular CM theta packet and the two-weight content at -23.

Run with Python 3.9+ from this directory, e.g. python -B verify_theta23_stepanov.py.
All comparisons use integer arithmetic. The cubic factorization test is in cubic23.py.
"""

from collections import Counter
from functools import reduce
from itertools import combinations
from math import comb, gcd, isqrt

from cubic23 import primes_up_to, trace_at_prime


def norm_representatives(n):
    """Generators modulo sign, in the basis 1, omega=(1+sqrt(-23))/2."""
    for b in range(isqrt(4 * n // 23) + 1):
        auxiliary = 4 * n - 23 * b * b
        t = isqrt(auxiliary)
        if t * t != auxiliary or (t - b) % 2:
            continue
        yield (t - b) // 2, b
        if b and t:
            yield (-t - b) // 2, b


def multiply(x, y):
    a, b = x
    c, e = y
    return a * c - 6 * b * e, a * e + b * c + b * e


def power(x, n):
    result = (1, 0)
    while n:
        if n % 2:
            result = multiply(result, x)
        x = multiply(x, x)
        n //= 2
    return result


def coefficient(reps, degree):
    terms = [power(x, degree) for x in reps]
    if sum(b for _, b in terms) != 0:
        raise ArithmeticError("Harmonic coefficient is not rational")
    return sum(a for a, _ in terms)


def row(p, degree, reps, chi):
    def corrected(e):
        return coefficient(reps, e) + (1 - chi) * (1 - p**e)

    return corrected(degree)**3, corrected(3 * degree)


def square_real(x):
    """Square a+b*sqrt(69)."""
    a, b = x
    return a * a + 69 * b * b, 2 * a * b


def twice_traces(delta, m, principal):
    """Twice the three first traces, in Q(sqrt(69)), for split p."""
    a, b = power(delta, 2 * m)
    t = 2 * a + b
    if principal:
        return [(2 * t, 0)] * 3
    return [(2 * t, 0), (-t, -b), (-t, b)]


def energy(traces):
    result = [0, 0]
    for (a, b), (c, e) in combinations(traces, 2):
        x, y = square_real((a - c, b - e))
        result[0] += x
        result[1] += y
    return tuple(result)


def check_packet(p, trace, chi, principal):
    if chi == 1:
        delta = next(x for x in norm_representatives(p**3) if gcd(*x) == 1)
    for m in (1, 2):
        d = 6 * m
        first = twice_traces(delta, m, principal) if chi == 1 else [(0, 0)] * 3
        # The r-th row below is 2^r times the ordinary Frobenius-power trace.
        # Newton's recurrence uses the determinant, including its sign at inert p.
        previous, current = [(2, 0)] * 3, first
        determinant = 4 * chi * p**d
        for r in range(1, 7):
            expected = 18 * 4**r * p**(r * d) if trace == -1 and r % 3 else 0
            if energy(current) != (expected, 0):
                raise ArithmeticError(f"Frobenius trace energy at {p}, m={m}, r={r}")
            following = []
            for base, value, (a, b) in zip(first, current, previous):
                c, e = real_product(base, value)
                following.append((c - determinant * a, e - determinant * b))
            previous, current = current, following


def real_product(x, y):
    a, b = x
    c, e = y
    return a * c + 69 * b * e, a * e + b * c


def determinant_with_one(x, y):
    def sub(a, b):
        return a[0] - b[0], a[1] - b[1]

    return sub(real_product(sub(x[1], x[0]), sub(y[2], y[0])),
               real_product(sub(x[2], x[0]), sub(y[1], y[0])))


def check_labels():
    x = twice_traces((1, 1), 1, False)
    y = twice_traces((1, 1), 2, False)
    # Scaling columns by (2,4) or (2,2) scales the minors by 8 or 4.
    minors = [determinant_with_one(x, list(map(square_real, x))),
              determinant_with_one(x, y)]
    for (a, b), scale, expected in zip(minors, (8, 4), (135, 576)):
        if a != 0 or abs(b) != scale * expected:
            raise ArithmeticError("Hecke-label determinant at 2")
    return x, y


def linear_product(coefficients, a, b):
    result = [0] * (len(coefficients) + 1)
    for i, value in enumerate(coefficients):
        result[i] += a * value
        result[i + 1] += b * value
    return result


def jet(coefficients, point, order):
    return sum(comb(i, order) * c * point**(i - order)
               for i, c in enumerate(coefficients) if i >= order)


def valuation(n, p):
    if not n:
        raise ValueError("valuation requires a nonzero integer")
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    return exponent


def check_polynomial(rows):
    degree = 6
    multiplicities = {2: 2, 5: 1, 13: 1, 47: 1, 59: 2, 101: 1}
    coefficients = [1]
    expected_content = 1
    expected_root = 0
    for p, count in multiplicities.items():
        a, b = rows[p, degree]
        for _ in range(count):
            coefficients = linear_product(coefficients, a, b)
        trace = trace_at_prime(p)
        if trace == -1:
            expected_root += count
            expected_content *= p**(3 * degree * count)
        elif trace == 2:
            expected_content *= 3**count
    size = sum(multiplicities.values())
    content = reduce(gcd, coefficients)
    root = next(j for j in range(size + 1) if jet(coefficients, -1, j))
    if content != expected_content or root != expected_root:
        raise ArithmeticError("Product content or horizontal multiplicity")
    first = jet(coefficients, -1, root)
    for p, count in multiplicities.items():
        if trace_at_prime(p) == 2 and valuation(first, p) < degree * count:
            raise ArithmeticError(f"Missing principal-class divisibility at {p}")
    for s in range(1, size + 2):
        jet_gcd = reduce(gcd, (jet(coefficients, t, j) for t in range(s)
                               for j in range(size // s + 1)))
        quotient, remainder = divmod(jet_gcd, content)
        if remainder:
            raise ArithmeticError("Jets lost coefficient content")
        # A positive quotient with all prime factors < s is the required equality
        # at every prime >= s, including primes outside the example support.
        for p in primes_up_to(s - 1):
            while quotient % p == 0:
                quotient //= p
        if quotient != 1:
            raise ArithmeticError(f"Primitive jet content at {s} points")
        folded = []
        for order in (s, s + 1):
            for j in range(size // s + 1):
                remainder = [0] * order
                for i in range(j, size + 1):
                    remainder[(i - j) % order] += comb(i, j) * coefficients[i]
                folded.extend(remainder)
        if reduce(gcd, folded) != content:
            raise ArithmeticError(f"Cyclotomic Hasse jets at orders {s}, {s+1}")
    binary_content = 1
    for bit in range(max(multiplicities.values()).bit_length()):
        layer = [1]
        for p, count in multiplicities.items():
            if count & (1 << bit):
                layer = linear_product(layer, *rows[p, degree])
        binary_content *= reduce(gcd, layer)**(1 << bit)
    if binary_content != content:
        raise ArithmeticError("Binary layers did not recover weighted content")
    return size, root, content


def check_first_jet(prime_coefficients):
    depth = 11
    multiplicities = {2: 5, 3: 1, 13: 6, 47: 2, 59: 3, 101: 4, 167: 1}
    actual = expected = 1
    for bit in range(max(multiplicities.values()).bit_length()):
        support = sorted(p for p, count in multiplicities.items() if count & (1 << bit))
        block_size = isqrt(len(support) - 1) + 1
        for start in range(0, len(support), block_size):
            block = support[start:start + block_size]
            polynomial = [1]
            for p in block:
                if pow(p, 11, 23) != 1:
                    raise ArithmeticError("Non-split auxiliary support in first-jet example")
                norm = p**depth
                polynomial = linear_product(polynomial, norm + prime_coefficients[p], norm)
            all_jets = [jet(polynomial, -1, j) for j in range(len(block) + 1)]
            index = next(j for j, value in enumerate(all_jets) if value)
            if index != sum(trace_at_prime(p) == -1 for p in block):
                raise ArithmeticError("Order of first nonzero block jet")
            actual *= abs(all_jets[index])**(1 << bit)
    for p, count in multiplicities.items():
        base = p**depth if trace_at_prime(p) == -1 else abs(prime_coefficients[p])
        expected *= base**count
    if expected != actual or prime_coefficients[59] != -253942:
        raise ArithmeticError("Selected first-jet product")
    return actual.bit_length()


def main():
    bound = 1000
    counts = Counter()
    rows = {}
    prime_coefficients = {}
    for p in primes_up_to(bound):
        if p == 23:
            continue
        trace = trace_at_prime(p)
        counts[trace] += 1
        chi = -1 if pow(p, 11, 23) != 1 else 1
        reps = list(norm_representatives(p * p))
        prime_reps = list(norm_representatives(p))
        b = coefficient(prime_reps, 6)
        prime_coefficients[p] = b
        check_packet(p, trace, chi, bool(prime_reps))
        for depth in (1, 6, 17):
            norm = p**depth
            error = b + 1 - chi
            expected_content = norm if trace == -1 else 1
            if gcd(norm, error) != expected_content or (error == 0) != (trace == -1):
                raise ArithmeticError(f"Fixed-weight row at {p}, depth={depth}")
        for d in (2, 6, 12, 18):
            a, b = row(p, d, reps, chi)
            rows[p, d] = a, b
            expected = {-1: p**(3 * d), 0: 1, 2: 3}[trace]
            if gcd(a, b) != expected:
                raise ArithmeticError(f"Row content at {p}, d={d}")
            if (a == b) != (trace == -1):
                raise ArithmeticError(f"Horizontal root at {p}, d={d}")
            if trace == 2 and valuation(a - b, p) != d:
                raise ArithmeticError(f"Principal row difference at {p}, d={d}")
            height = max(abs(a), abs(b))
            if not p**(3 * d) <= 8 * height <= 216 * p**(3 * d):
                raise ArithmeticError(f"Archimedean height at {p}, d={d}")

    jet_bits = check_first_jet(prime_coefficients)
    x, y = check_labels()
    b36 = coefficient(list(norm_representatives(36)), 6)
    if b36 % 6 != 1:
        raise ArithmeticError("Expected nonmultiplicativity at norm 36")
    size, root, content = check_polynomial(rows)
    print(f"Unramified p <= {bound}: {sum(counts.values())}, "
          f"(split, transposition, 3-cycle) = {counts[2], counts[0], counts[-1]}")
    print("Trace energies checked at weights 7, 13 and Frobenius powers 1 through 6.")
    print("Row contents, horizontal roots and heights checked at d=2,6,12,18.")
    print(f"Twice the traces at p=2, expressed in (1, sqrt(69)): {x}, {y}")
    print(f"B_6(36) = {b36}, residue mod 6 = {b36 % 6}; B_6(59) = {prime_coefficients[59]}")
    print(f"Fixed-weight block jets reconstruct the {jet_bits}-bit absolute product.")
    print(f"Product degree = {size}, multiplicity at -1 = {root}, content = {content}")
    print("Integer-grid jets verified away from primes < s. Two cyclotomic grids give exact gcds")
    print("for s=1,...,degree+1, including small characteristics. Binary layers recombine exactly.")


if __name__ == "__main__":
    main()
