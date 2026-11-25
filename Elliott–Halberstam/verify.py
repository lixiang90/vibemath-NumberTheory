"""Finite checks for the Elliott–Halberstam zero-sum reduction (python-flint 0.9)."""

import cmath
import math
from collections import defaultdict
from functools import cache

from flint import acb, arb, ctx, dirichlet_char


HEIGHT = 35
WIDTH = 0.30
ENDPOINT = 1024


def units(modulus):
    return [a for a in range(1, modulus + 1) if math.gcd(a, modulus) == 1]


def certify_right_half_plane(chi, height, left=0.5, right=2.0):
    """Enclose L(sigma+it) on the whole horizontal segment by subdivision."""
    middle = (left + right) / 2
    box = acb(arb(middle, (right - left) / 2), height)
    if chi.l_function(box).real > 0:
        return
    if right - left < 2**-20:
        raise ArithmeticError("Could not certify the horizontal argument")
    certify_right_half_plane(chi, height, left, middle)
    certify_right_half_plane(chi, height, middle, right)


def zero_count(chi):
    """Argument-principle count in the full strip, including both signs."""
    for height in (-HEIGHT, HEIGHT):
        certify_right_half_plane(chi, height)
    top = acb(0.5, HEIGHT)
    gamma_phase = ((top + chi.parity()) / 2).lgamma().imag
    delta = HEIGHT * (arb(chi.modulus()) / arb.pi()).log() + 2 * gamma_phase
    delta += chi.l_function(top).arg() - chi.l_function(top.conjugate()).arg()
    count = (delta / arb.pi()).unique_fmpz()
    if count is None:
        raise ArithmeticError("Winding number does not isolate an integer")
    return int(count)


def hardy_sign(chi, t):
    value = chi.hardy_z(t).real
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise ArithmeticError("Indeterminate sign at a sampling point")


@cache
def zero_ordinates(conductor, label):
    chi = dirichlet_char(conductor, label)
    expected = zero_count(chi)
    for subdivisions in (8, 16, 32, 64):
        brackets = []
        left = -HEIGHT
        sign = hardy_sign(chi, left)
        for index in range(1, 2 * HEIGHT * subdivisions + 1):
            right = -HEIGHT + index / subdivisions
            next_sign = hardy_sign(chi, right)
            if next_sign != sign:
                brackets.append((left, right, sign))
            left, sign = right, next_sign
        if len(brackets) == expected:
            break
    else:
        raise ArithmeticError("Sign changes do not account for the winding number")

    zeros = []
    for left, right, sign in brackets:
        for _ in range(38):
            middle = (left + right) / 2
            if hardy_sign(chi, middle) == sign:
                left = middle
            else:
                right = middle
        zeros.append((left + right) / 2)
    return zeros


def primitive_inducer(chi, residues):
    r = chi.conductor()
    modulus_exponent = chi.group().exponent()
    for label in units(r):
        candidate = dirichlet_char(r, label)
        conductor_exponent = candidate.group().exponent()
        if candidate.is_primitive() and all(
            chi.chi_exponent(a) * conductor_exponent
            == candidate.chi_exponent(a) * modulus_exponent for a in residues
        ):
            return r, label
    raise ArithmeticError("Missing primitive character")


def zero_terms(r, label):
    single_edge = 0j
    dyadic = 0j
    for gamma in zero_ordinates(r, label):
        rho = 0.5 + 1j * gamma
        weight = cmath.exp(1j * gamma * math.log(ENDPOINT) + WIDTH**2 * rho**2 / 2) / rho
        single_edge += weight
        dyadic += weight * (1 - 2 ** (-rho))
    return single_edge, -math.sqrt(ENDPOINT) * dyadic


def trivial_zero_terms(parity):
    """Residues at 0, -1, -2, -3 crossed on moving to Re(s)=-7/2."""
    value = math.log(2) if parity == 0 else 0
    for m in range(1, 4):
        if m % 2 == parity:
            value += (2**m - 1) * ENDPOINT ** (-m) * math.exp(WIDTH**2 * m**2 / 2) / m
    return value


def prime_power_weights():
    """Convolve the two log endpoints with the same Gaussian."""
    limit = math.ceil(ENDPOINT * math.exp(11 * WIDTH))
    prime = bytearray([1]) * (limit + 1)
    result = []
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        for multiple in range(2 * p, limit + 1, p):
            prime[multiple] = 0
        n = p
        while n <= limit:
            t = math.log(n / ENDPOINT) / (WIDTH * math.sqrt(2))
            weight = (math.erfc(t) - math.erfc(t + math.log(2) / (WIDTH * math.sqrt(2)))) / 2
            result.append((n, p, math.log(p) * weight))
            n *= p
    return result


def fourth_moment(modulus, residues, coefficients):
    phi = len(residues)
    direct = 0.0
    for a in residues:
        z = sum(complex(dirichlet_char(modulus, label)(a)).conjugate() * v
                for label, v in coefficients.items()) / phi
        direct += abs(z) ** 4 / phi
    convolution = defaultdict(complex)
    for first, v in coefficients.items():
        for second, w in coefficients.items():
            convolution[first * second % modulus] += v * w
    via_products = sum(abs(v) ** 2 for v in convolution.values()) / phi**4
    diagonal = (2 * sum(abs(v)**2 for v in coefficients.values())**2
                - sum(abs(v)**4 for v in coefficients.values())) / phi**4
    error = abs(direct - via_products)
    if error > 1e-12 * max(1.0, direct):
        raise ArithmeticError("Fourth moment normalization mismatch")
    return error, via_products - diagonal


def check_character_relations():
    """Conrey labels identify the character group with the units modulo q."""
    checked = 0
    for q in (8, 15, 16, 24, 35, 45):
        labels = units(q)
        for m, n in ((1, 0), (0, -2), (2, -4), (3, 6), (4, 4)):
            d = math.gcd(m, n)
            fibres = defaultdict(int)
            for eta in labels:
                for nu in labels:
                    fibres[pow(eta, m, q) * pow(nu, n, q) % q] += 1
            image = {pow(label, d, q) for label in labels}
            expected_size = len(labels) * sum(pow(label, d, q) == 1 for label in labels)
            if set(fibres) != image or any(size != expected_size for size in fibres.values()):
                raise ArithmeticError("Character relation fibre mismatch")
            checked += 1
    print(f"exact character-group fibre checks: {checked}")


def check_modulus(q, primes):
    residues = units(q)
    direct = {a: 0.0 for a in residues}
    for n, _, weight in primes:
        if n % q in direct:
            direct[n % q] += weight
    mean = sum(direct.values()) / len(residues)
    reconstructed = {a: 0j for a in residues}
    coefficients = {}
    for label in residues:
        chi = dirichlet_char(q, label)
        if chi.is_principal():
            continue
        r, primitive_label = primitive_inducer(chi, residues)
        primitive = dirichlet_char(r, primitive_label)
        z, dyadic = zero_terms(r, primitive_label)
        dyadic -= trivial_zero_terms(primitive.parity())
        for n, p, weight in primes:
            if q % p == 0 and r % p != 0:
                dyadic -= complex(primitive(n)) * weight
        for a in residues:
            reconstructed[a] += complex(chi(a)).conjugate() * dyadic / len(residues)
        # Small fixed cuts exercise both filters in the finite algebra check.
        if r > 3 and chi.order() > 2:
            coefficients[label] = z
    ap_error = max(abs(direct[a] - mean - reconstructed[a]) for a in residues)
    if ap_error > 2e-7:
        raise ArithmeticError(f"Smoothed explicit formula mismatch for q={q}: {ap_error}")
    moment_error, off_diagonal = fourth_moment(q, residues, coefficients)
    print(f"q={q:2d}  kept={len(coefficients)}  AP error={ap_error:.3e}  "
          f"fourth-moment error={moment_error:.3e}  off-diagonal={off_diagonal:+.8e}")


def main():
    ctx.prec = 90
    check_character_relations()
    primes = prime_power_weights()
    for q in (5, 10, 15):
        check_modulus(q, primes)
    total = 0
    for r in (3, 5, 15):
        counts = []
        for label in units(r):
            if dirichlet_char(r, label).is_primitive():
                count = len(zero_ordinates(r, label))
                counts.append(f"{label}={count}")
                total += count
        print(f"conductor {r}, counts by Conrey label: " + ", ".join(counts))
    print(f"certified zeros counted across characters: {total}")


if __name__ == "__main__":
    main()
