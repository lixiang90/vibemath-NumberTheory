"""Exact constants for the explicit interval containing ordinary primes."""

import argparse
import json
from fractions import Fraction
from math import gcd, isqrt, prod
from pathlib import Path


def require(condition, message):
    if not condition:
        raise ValueError(message)


def interval_certificate(prefix):
    m = 5460
    units = {r for r in range(m) if gcd(r, m) == 1}
    forbidden = {r for r in units if r % 3 == 1 and gcd(r*r+r+1, 91) > 1
                 and (r % 12 != 1 or gcd(r+1, 455) > 1)}
    require(len(units) == 1152 and len(forbidden) == 172, "progression count")
    alpha = Fraction(2*len(forbidden), len(units))
    require(alpha == Fraction(43, 144), "sieve dimension")
    extra = [q for q in range(7, 10001, 6) if q % m not in forbidden
             and all(q % d for d in range(2, isqrt(q)+1))]
    require(prod(q-1 for q in extra) > 4*prod(q-3 for q in extra), "finite Euler product")
    require(Fraction(29, 10)**3 > 24 and 1/(1-Fraction(344, 840*20)) < Fraction(103, 100)
            and Fraction(29, 10)*Fraction(103, 100) < 3, "Mertens constant")
    require(1+alpha < Fraction(13, 10) and 25**13 < 66**10, "counting constant")
    require(Fraction(101, 100)*199/alpha < 675 and 3800**12 < 10**43, "tail constant")
    x = 2*10**12
    require(prefix["cutoff"] == 3000000 and prefix["upper"] == x and prefix["scale"] == 10**18,
            "wrong finite interval")
    total = prefix["scaled_reciprocal_upper"]
    require(total == 72401332985083531 and prefix["survivors"] == 8391832364, "finite sieve record")
    k = 70_000_000_000
    upper = Fraction(total, 10**18) + Fraction(k, x) + Fraction(675, 3800) + Fraction(1, 2999999)
    lower = sum(Fraction((-1)**(j+1), j*3**j) for j in range(1, 5))
    require(upper < lower, "logarithmic comparison")
    return {"excluded_classes": len(forbidden), "finite_extra_primes": len(extra),
            "ordinary_primes_strictly_more_than": k, "lower_endpoint_exclusive": x,
            "log_of_upper_endpoint": 10**12, "upper_log_if_at_most_k": str(upper),
            "lower_log_four_thirds": str(lower)}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("prefix", nargs="?", type=Path,
                     default=Path(__file__).resolve().parent.parent /
                     "restrained-exponent/ordinary-prefix-2000000000000.json")
    data = json.loads(cli.parse_args().prefix.read_text(encoding="utf-8"))
    print(json.dumps(interval_certificate(data), indent=2))
