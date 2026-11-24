"""Recover the normalized Euler exponent modulo 3**k from a formal factorization."""

import argparse
import json
from math import gcd


def strip_three(n):
    if n <= 0:
        raise ValueError("a positive integer is required")
    valuation = 0
    while n % 3 == 0:
        valuation += 1
        n //= 3
    return valuation, n


def geometric_sum(base, count, modulus):
    """Sum of the first count powers, reduced without dividing by base - 1."""
    total, power = 0, 1
    for bit in bin(count)[2:]:
        total = total * (1 + power) % modulus
        power = power * power % modulus
        if bit == "1":
            total = (total + power) % modulus
            power = power * base % modulus
    return total


def normalized_block(base, height, modulus):
    """Return S_(3**height)(base)/3**height and base**(3**height) modulo 3**k."""
    if base % 3 != 1:
        raise ValueError("a normalized block needs a base equal to 1 modulo 3")
    unit = 1
    base %= modulus
    while height and base != 1:
        # Replacing base by base + modulus preserves this quotient modulo modulus.
        unit = unit * ((base * base + base + 1) // 3) % modulus
        base = pow(base, 3, modulus)
        height -= 1
    return unit, base


def euler_residue(p, a, components, k):
    """Use p for the Euler base, a for v_3(N), and (r, even exponent) pairs.

    The computation also applies to pairwise coprime formal bases. The OPN
    application assumes that the specified bases are prime. A negative
    valuation budget returns None.
    """
    if p <= 1 or p % 4 != 1 or p % 3 == 0 or a < 0 or a % 2 or k < 1:
        raise ValueError("invalid Euler base, 3-exponent, or precision")
    modulus = 3**k
    g, reduced_p_plus_one = strip_three(p + 1)
    h = a - g
    unit = reduced_p_plus_one * pow(p, -1, modulus) % modulus
    bases = [p]
    for r, exponent in components:
        if (r < 5 or gcd(r, 6) != 1 or exponent < 2 or exponent % 2
                or any(gcd(r, other) != 1 for other in bases)):
            raise ValueError("invalid ordinary component or overlapping bases")
        bases.append(r)
        residue = pow(r, -1, modulus)
        length = exponent + 1
        if r % 3 == 1:
            valuation, length = strip_three(length)
            h -= valuation
            block, residue = normalized_block(residue, valuation, modulus)
            unit = unit * block % modulus
        unit = unit * geometric_sum(residue, length, modulus) % modulus
    if h < 0:
        return None
    sigma_three = (pow(3, a + 1, modulus) - 1) * pow(2, -1, modulus)
    target = 2 * pow(unit * sigma_three, -1, modulus) % modulus
    scale, base = normalized_block(pow(p * p, -1, modulus), h, modulus)
    m, place = 0, 1
    for _ in range(k):
        next_modulus = 3 * place
        for digit in (0, 1, 2):
            candidate = m + digit * place
            value = scale * geometric_sum(base, candidate, next_modulus)
            if (value - target) % next_modulus == 0:
                m = candidate
                break
        else:
            raise ArithmeticError("the ternary lift failed")
        place = next_modulus
    return {"h": h, "m_residue": m, "modulus": modulus}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("p", type=int, help="Euler base")
    parser.add_argument("a", type=int, help="exponent of 3 (zero if absent)")
    parser.add_argument("k", type=int, help="ternary precision")
    parser.add_argument("--factor", nargs=2, action="append", type=int, default=[],
                        metavar=("BASE", "EVEN_EXPONENT"))
    args = parser.parse_args()
    print(json.dumps(euler_residue(args.p, args.a, args.factor, args.k), indent=2))
