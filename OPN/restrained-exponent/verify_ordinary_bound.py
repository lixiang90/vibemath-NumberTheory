"""Bound the ordinary abundancy by rational intervals and a finite sieve sum."""

import json
from fractions import Fraction
from itertools import combinations
from math import ceil, floor, prod
from pathlib import Path


SCALE = 10**18
PRIMES = (7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103)
PUBLISHED_PREFIXES = {2*10**11: (917142659, 62696136715648189),
                      2*10**12: (8391832364, 72401332985083531)}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def unit_log(x):
    """For 1 <= x <= 2 use the positive atanh series and a geometric tail."""
    x = Fraction(x)
    require(1 <= x <= 2, "logarithm argument outside unit interval")
    t = (x - 1) / (x + 1)
    power = t
    t2 = t*t
    lower = Fraction(0)
    for k in range(30):
        lower += 2*power/(2*k+1)
        power *= t2
    return lower, lower + 2*power/(61*(1-t2))


def log_interval(value):
    x = Fraction(value)
    require(x >= 1, "expected logarithm of a number at least 1")
    k = x.numerator.bit_length() - x.denominator.bit_length()
    if 2**k > x:
        k -= 1
    lo, hi = unit_log(x / 2**k)
    lo2, hi2 = unit_log(Fraction(2))
    return (Fraction(floor(SCALE*(lo+k*lo2)), SCALE),
            Fraction(ceil(SCALE*(hi+k*hi2)), SCALE))


def majorant():
    """Pointwise polynomial from Stormwhip's verify.py, checked on every mask."""
    terms = [(subset, (-1)**k) for k in range(3) for subset in combinations(PRIMES, k)]
    terms.extend((c, -1) for c in combinations(PRIMES, 3) if 7 in c and prod(c) <= 10000)
    require(len(terms) == 97 and max(prod(c) for c, _ in terms) <= 10000, "modulus window")
    masks = [(sum(1 << PRIMES.index(p) for p in c), s) for c, s in terms]
    for b in range(1 << len(PRIMES)):
        value = sum(sign for mask, sign in masks if b & mask == mask)
        require(value >= int(b == 0), "negative sieve majorant")
    density = sum(sign * prod(Fraction(2, p-1) for p in c) for c, sign in terms)
    # The empty monomial counts the two reduced classes modulo 3 on (X,Z].
    classes = sum(2**len(c) if c else 2 for c, _ in terms)
    require(density == Fraction(64959667, 169646400) and classes == 434, "majorant constants")
    return density, Fraction(classes, 840)


def bound_for_prefix(prefix, a, z):
    x = prefix["upper"]
    require(prefix["cutoff"] == 3000000 and prefix["scale"] == SCALE, "wrong finite sieve")
    require(8*10**9 <= x < z and a > 0, "invalid analytic interval")
    count, total = prefix["survivors"], prefix["scaled_reciprocal_upper"]
    require(count >= 0 and total >= count*((SCALE+x-2)//(x-1)), "invalid prefix sum")
    if x in PUBLISHED_PREFIXES:
        require((count, total) == PUBLISHED_PREFIXES[x], "published prefix record changed")
    density, error = majorant()
    lx, _ = log_interval(x)
    lz, uz = log_interval(z)
    _, ratio = log_interval(uz / lx)
    integration_error = 1/lx - 1/uz + 3/(2*lx**2) + 1/(2*lz**2)
    middle = Fraction(x, x-1) * (density*ratio + error*integration_error)
    high = Fraction(2*a+1, z-1)
    special = Fraction(1, prefix["cutoff"]-1)
    upper = Fraction(total, SCALE) + middle + high + special
    target, _ = log_interval(Fraction(4, 3))
    require(upper < target, "prefix and analytic tail do not exclude the claimed exponent")
    return {"assumed_a_at_most": a, "upper_sieve_endpoint": x, "last_split": z,
            "upper_log_product_scaled_1e18": ceil(upper*SCALE),
            "lower_log_four_thirds_scaled_1e18": floor(target*SCALE)}


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    directory = parser.parse_args().directory
    for x, a in [(2*10**11, 10**17), (2*10**12, 15*10**17)]:
        data = json.loads((directory / f"ordinary-prefix-{x}.json").read_text(encoding="utf-8"))
        require(data["upper"] == x, "wrong sieve endpoint")
        print(json.dumps(bound_for_prefix(data, a, 250*a), indent=2))
