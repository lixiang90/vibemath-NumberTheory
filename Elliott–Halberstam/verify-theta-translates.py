"""Finite checks for translated cyclic character weights and the Stepanov ramp."""

import math
import random
from fractions import Fraction
from functools import cache

from flint import acb, acb_series, arb, ctx, dirichlet_char

from verify import primitive_inducer, units, zero_ordinates


def ball(q):
    return arb(q.numerator) / q.denominator


@cache
def root(order, k):
    return acb(arb(2 * (k % order)) / order).exp_pi_i()


def prime_powers(limit):
    for p in range(2, limit + 1):
        if any(p % d == 0 for d in range(2, math.isqrt(p) + 1)):
            continue
        n = p
        while n <= limit:
            yield n, arb(p).log()
            n *= p


def moebius(n):
    factors = 0
    for p in range(2, math.isqrt(n) + 1):
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            factors += 1
    return (-1) ** (factors + int(n > 1))


def cyclic_groups(q, residues):
    groups = {}
    for label in residues:
        powers = [1]
        j = label
        while j != 1:
            powers.append(j)
            j = j * label % q
        groups.setdefault(frozenset(powers), powers)
    return groups


def assembled(groups, measurements):
    result = acb(0)
    for outer in groups:
        for inner in groups:
            if inner <= outer:
                result += moebius(len(outer) // len(inner)) * measurements[inner]
    return result


def ramp(n, y):
    return max(Fraction(y, n) - 1, Fraction(0))


def measured_level(q, b, powers, values, perturbation):
    k = len(powers)
    if k == 1:
        return acb(0), arb(0)
    chi = dirichlet_char(q, powers[1])
    index, remainder = divmod(k * int(chi.chi_exponent(b)), int(chi.group().exponent()))
    assert remainder == 0
    nu = [Fraction(3 + int(u == index), 3 * k + 1) for u in range(k)]
    error = [Fraction(0)] * k
    error[(index + 1) % k] += perturbation / (20 * (3 * k + 1) ** 2)
    error[index] -= perturbation / (20 * (3 * k + 1) ** 2)
    assert sum(nu) == 1 and sum(error) == 0
    assert min(n + e for n, e in zip(nu, error)) >= 0
    level, noise_energy, signal_energy = acb(0), arb(0), arb(0)
    for j, label in enumerate(powers[1:], 1):
        coefficient = sum(ball(n + e) * root(k, -u * j)
                          for u, (n, e) in enumerate(zip(nu, error)))
        noise = coefficient - root(k, -index * j) / (3 * k + 1)
        noise_energy += abs(noise) ** 2
        signal_energy += abs(values[label]) ** 2
        level += (3 * k + 1) * coefficient * values[label]
    exact_energy = k * sum(e**2 for e in error)
    assert noise_energy.contains(ball(exact_energy))
    return level, (3 * k + 1) * (ball(exact_energy) * signal_energy).sqrt()


def check_modulus(q):
    y, h = 240, Fraction(1, 8)
    y2 = (1 + h) * y
    residues = units(q)
    g = len(residues)
    groups = cyclic_groups(q, residues)
    powers = list(prime_powers(int(y2)))
    chars = {l: dirichlet_char(q, l) for l in residues if l != 1}
    primitive = {l: dirichlet_char(*primitive_inducer(chi, residues)) for l, chi in chars.items()}
    values = {y: {}, y2: {}}
    extra = {}
    for l, chi in primitive.items():
        jet = acb_series([1, 1], 2).dirichlet_l(chi)
        linear = jet[1] / jet[0]
        for t in values:
            values[t][l] = ball(Fraction(t)) * linear + sum(
                weight * chi(n) * ball(ramp(n, t)) for n, weight in powers)
        extra[l] = sum(weight * chi(n) * ball((ramp(n, y2) - (1 + h) * ramp(n, y)) / h)
                       for n, weight in powers if math.gcd(q, n) != 1)

    max_error = max_ratio = 0.0
    for b in residues:
        exact = {}
        for t, z in values.items():
            direct = sum(chi(b).conjugate() * z[l] for l, chi in chars.items()) / g
            exact_levels, noisy_levels, bounds = {}, {}, {}
            for group, gen in groups.items():
                exact_levels[group], _ = measured_level(q, b, gen, z, Fraction(0))
                noisy_levels[group], bounds[group] = measured_level(
                    q, b, gen, z, Fraction(1 if t == y else -1))
            exact[t] = assembled(groups, exact_levels) / g
            assert (exact[t] - direct).contains(0)
            diff = abs(assembled(groups, noisy_levels) / g - direct)
            bound = sum(abs(moebius(len(outer) // len(inner))) * bounds[inner]
                        for outer in groups for inner in groups if inner <= outer) / g
            # Cauchy can be an equality for characters of order two.
            assert diff.upper() <= bound.lower() + arb("1e-24")
            max_ratio = max(max_ratio, float(diff.upper() / bound.lower()))

        recovered = (exact[y2] - ball(1 + h) * exact[y]) / ball(h)
        correction = sum(chars[l](b).conjugate() * v for l, v in extra.items()) / g
        arithmetic = arb(0)
        for n, weight in powers:
            ap = Fraction(int(n % q == b)) - Fraction(int(math.gcd(q, n) == 1), g)
            arithmetic += weight * ball(ap * (ramp(n, y2) - (1 + h) * ramp(n, y)) / h)
        discrepancy = recovered - correction - arithmetic
        assert discrepancy.contains(0)
        max_error = max(max_error, float(abs(discrepancy).upper()))
    print(f"q={q:2d}, residues={g:2d}, cyclic subgroups={len(groups):2d}, "
          f"ramp error <= {max_error:.2e}, max noise/bound={max_ratio:.6f}")


def check_invisible_data():
    p, m = 383, 16
    assert all(p % t for t in range(2, math.isqrt(p) + 1)) and p % 6 == 5
    g = p - 1
    signs, other = [1], [1]
    while len(signs) < m:
        signs, other = signs + other, signs + [-v for v in other]
    values = [0] * g
    for j, sign in enumerate(signs, start=g // 8):
        for u, factor in ((1, 1), (-1, 1), (3, -1), (-3, -1)):
            index = u * j % g
            assert values[index] == 0
            values[index] = factor * sign
    assert all(values[j] == values[-j % g] for j in range(g))
    for d in range(1, g + 1):
        if g % d == 0:
            assert d * sum(values[::d]) == 0
    assert sum(v**4 for v in values) == 4 * m
    print(f"real reciprocal-even model p={p}, m={m}, all cyclic sums zero, g*m*M4/c^4=4")


def check_zero_labels():
    p = 5
    roots = {l: zero_ordinates(p, l) for l in (2, 3, 4)}
    for l in (2, 3):
        gamma = min(x for x in roots[l] if x > 0)
        distance = min(abs(gamma - x) for t, values in roots.items() if t != l for x in values)
        assert distance > 0.1
        residues = {b: dirichlet_char(p, l)(b).conjugate() for b in units(p)}
        for j in units(p):
            inverse = sum(dirichlet_char(p, j)(b) * v for b, v in residues.items()) / (p - 1)
            assert inverse.contains(int(j == l))
        print(f"q=5, label={l}, first positive ordinate={gamma:.12f}, "
              f"separation from other nonprincipal labels > {math.floor(100*distance)/100}")


def check_sharp_profile():
    for k in (7, 17, 41):
        squares = {u * u % k for u in range(1, k)}
        sign = [0 if u == 0 else (1 if u in squares else -1) for u in range(k)]
        noise = [sum(s * root(k, -u * j) for u, s in enumerate(sign)) for j in range(1, k)]
        assert all((abs(value)**2).contains(k) for value in noise)
        total = sum(abs(value) for value in noise)
        assert (total - (k - 1) * arb(k).sqrt()).contains(0)
    print("sharp profile: quadratic Gauss magnitudes and loss verified at orders 7, 17, 41")


def check_balanced_selection():
    # A cyclic sample with two test functions. B counts marked adjacent pairs.
    # Integer arithmetic is exact after clearing the Bernoulli probabilities' denominator.
    n, k, fav = 1024, 4, 3
    d0 = n // 2
    rate8 = [4 if j % k == fav else 3 for j in range(n)]
    mean_counts8 = [sum(rate8[u::k]) for u in range(k)]
    mean_total8 = sum(mean_counts8)
    b0_64 = sum(p * rate8[(j + 1) % n] for j, p in enumerate(rate8))
    # Values of two mean-zero Fourier polynomials on the four quotient classes.
    tests = [(2, -3, 4, -3), (3, -2, -3, 2)]
    rand = random.Random(15299)
    for trial in range(1, 1001):
        mask = [rand.randrange(8) < p for p in rate8]
        counts = [sum(mask[u::k]) for u in range(k)]
        size = sum(counts)
        errors8 = [8 * v - a for v, a in zip(counts, mean_counts8)]
        bad = sum(flag and mask[(j + 1) % n] for j, flag in enumerate(mask))
        if (sum(e*e for e in errors8) > 64 * 8 * d0
                or (8 * size - mean_total8)**2 > 64 * 8 * d0
                or 64 * bad > 4 * b0_64):
            continue
        ratios = []
        for W in tests:
            w0 = Fraction(sum(a * v for a, v in zip(mean_counts8, W)), mean_total8)
            residual = sum(Fraction(e, 8) * (v - w0) for e, v in zip(errors8, W))
            variance = sum(Fraction(p * (8 - p), 64) * (W[j % k] - w0)**2
                           for j, p in enumerate(rate8))
            energy = Fraction(sum(v*v for v in W), k)
            assert variance <= 2 * d0 * energy
            assert residual == sum(c * v for c, v in zip(counts, W)) - size * w0
            ratios.append(residual**2 / variance)
        if max(ratios) <= 8:
            assert d0 // 2 <= size <= 2 * d0
            profile = [Fraction(c, size) - Fraction(a, mean_total8)
                       for c, a in zip(counts, mean_counts8)]
            assert sum(e**2 for e in profile) * d0 <= 128
            print(f"selection N={n}, k={k}, trial={trial}, size={size}, bad={bad}, "
                  f"E[bad]={float(Fraction(b0_64,64)):.1f}, "
                  f"max marked squared error/variance={float(max(ratios)):.6f}")
            return
    raise ArithmeticError("No simultaneous sample was found")


if __name__ == "__main__":
    ctx.prec = 120
    check_zero_labels()
    check_invisible_data()
    check_sharp_profile()
    check_balanced_selection()
    for modulus in (8, 15, 35, 45):
        check_modulus(modulus)
