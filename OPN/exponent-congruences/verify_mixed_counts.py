"""Exact upper bounds for the logarithmic abundancy in the mixed-index paper."""

from fractions import Fraction
from math import comb, prod


MIN_BASE = 10**8
MAX_DEGREE = 20
BLOCKS = 512
SCALE = 10**14
SMALL_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)


def group_bound(size, denominator, density, discrepancy):
    """Upper sum over positions 0,...,size-1, rounded up block by block."""
    step = (size + BLOCKS - 1) // BLOCKS
    units = 0
    for first in range(0, size, step):
        length = min(step, size - first)
        num = length * density.numerator
        den = denominator * density.numerator + max(0, first - discrepancy) * density.denominator
        units += (SCALE * num + den - 1) // den
    return Fraction(units, SCALE)


def upper_bound(ell, root, k):
    if root**(ell - 1) != MIN_BASE or ell not in (3, 5):
        raise ValueError("expected the exact root for ell = 3 or 5")
    primes = [p for p in SMALL_PRIMES if p != ell]
    density = Fraction(ell - 2, ell) * prod(Fraction(p - 1, p) for p in primes)
    discrepancy = (ell - 2) * 2**len(primes)
    finite = Fraction(k, MIN_BASE - 1)
    for j in range(1, MAX_DEGREE + 1):
        size = comb(k + j - 1, j)
        denominator = max(MIN_BASE - 1, root**j - 2)
        finite += group_bound(size, denominator, density, discrepancy)

    first = MAX_DEGREE + 1
    ratio = Fraction(k + first, (first + 1) * root)
    if ratio >= 1 or root**first < MIN_BASE:
        raise ArithmeticError("geometric tail estimate unavailable")
    tail = Fraction(MIN_BASE, MIN_BASE - 2) * Fraction(comb(k + first - 1, first), root**first)
    return finite + tail / (1 - ratio)


def main():
    # First five positive terms of 2 atanh(1/3) are a lower bound for log(2).
    log_two_lower = sum(Fraction(2, (2*i + 1) * 3**(2*i + 1)) for i in range(5))
    for ell, root, k, units in ((3, 10000, 17864, 693108), (5, 100, 230, 692821)):
        result = upper_bound(ell, root, k)
        bound = Fraction(units, 10**6)
        if not result < bound < log_two_lower:
            raise ArithmeticError(f"comparison failed for ell={ell}, k={k}")
        print(f"ell={ell}, k<={k}: upper bound < {units}/1000000 < log(2)")


if __name__ == "__main__":
    main()
