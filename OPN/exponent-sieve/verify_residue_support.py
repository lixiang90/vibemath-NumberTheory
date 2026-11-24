"""Verify the reciprocal sum and count in the seven permitted classes modulo 15."""

import json
from fractions import Fraction as F
from math import isqrt


LOW = 10**8
BREAK = 8_000_000_000
END = 455_000_000_000_000_000
SCALE = 10**16


def ceil_fraction(x):
    return -(-x.numerator // x.denominator)


def small_primes(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[:2] = b"\0\0"
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p::p] = bytes((limit - p*p)//p + 1)
    return [p for p in range(3, limit + 1, 2) if flags[p]]


def finite_interval(lo, hi, segment=2_000_000):
    """Count lo < p <= hi, p not 1 mod 15, and upper-bound sum 1/(p-1)."""
    if lo < 5 or hi < lo or segment < 1:
        raise ValueError("invalid range")
    primes = small_primes(isqrt(hi))
    count, weighted = 0, 0
    first = lo + 1 + (lo & 1)
    while first <= hi:
        size = min(segment, (hi - first)//2 + 1)
        last = first + 2*(size - 1)
        flags = bytearray(b"\x01") * size
        for p in primes:
            if p*p > last:
                break
            offset = (-first*((p + 1)//2)) % p
            if first + 2*offset == p:
                offset += p
            if offset < size:
                flags[offset::p] = bytes((size - 1 - offset)//p + 1)
        # first + 2*i == 1 mod 15 iff i == excluded mod 15.
        excluded = (8*(1 - first)) % 15
        flags[excluded::15] = bytes(len(range(excluded, size, 15)))
        left = 0
        while left < size:
            denominator = first + 2*left - 1
            right = min(size, left + max(1, denominator//8000))
            number = flags.count(1, left, right)
            count += number
            weighted += (SCALE*number + denominator - 1)//denominator
            left = right
        first = last + 2
    return count, weighted


def log_bounds(x):
    """Outward rational bounds, using log y = 2 atanh((y-1)/(y+1))."""
    x = F(x)
    if x < 1:
        raise ValueError("logarithm input below one")
    power = x.numerator.bit_length() - x.denominator.bit_length()
    if 2**power > x:
        power -= 1

    def series(y):
        z = F(y - 1, y + 1)
        lower = 2*sum(z**(2*j + 1)/(2*j + 1) for j in range(20))
        return lower, lower + 2*z**41/(41*(1 - z*z))

    lo2, hi2 = series(F(2))
    lo, hi = series(x / 2**power)
    precision = 10**12
    return (F(int((power*lo2 + lo)*precision), precision),
            F(ceil_fraction((power*hi2 + hi)*precision), precision))


def verify():
    count, weighted = finite_interval(LOW, BREAK)
    b0, b1 = log_bounds(BREAK)
    x0, x1 = log_bounds(END)
    ratio = log_bounds(x1/b0)[1]
    density = F(7, 8)
    error = F(7, 840)*(1/b0 - 1/x1 + 1/b0**2 + 1/x0**2)
    tail = F(BREAK, BREAK-1)*(density*ratio + error)
    upper_weight = F(weighted, SCALE) + tail + F(1, LOW)
    lower_count = (density*(END*(1/x1 + 1/x1**2) - BREAK*(1/b0 + 1/b0**2))
                   - F(7, 840)*(END/x0**2 + BREAK/b0**2))
    if upper_weight >= log_bounds(2)[0] or lower_count <= 10_030_000_000_000_000:
        raise ArithmeticError("the abundance or prime-count comparison failed")
    return {"lower_endpoint": LOW, "analytic_from": BREAK, "upper_endpoint": END,
            "small_admissible_primes": count, "finite_weight_scale": SCALE,
            "finite_weight_upper_numerator": weighted,
            "total_weight_upper_in_trillionths": ceil_fraction(upper_weight*10**12),
            "analytic_count_lower_floor": int(lower_count)}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
