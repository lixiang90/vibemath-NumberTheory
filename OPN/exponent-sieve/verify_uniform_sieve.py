"""Exact arithmetic for the explicit ternary cyclotomic sieve."""

from fractions import Fraction as F
from math import isqrt, prod


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify(prefix=None):
    split = [q for q in range(7, 101, 6) if all(q % d for d in range(2, isqrt(q)+1))]
    correction = prod(1 + F(5*q+3, (q-3)*(q+1)**2) for q in split)
    require(correction < F(123, 100) and F(123, 100)*F(64, 63) < F(5, 4),
            "exceptional-prime product")
    c_upper = F(3, 2)*F(22, 7)**2/6*prod(1-F(1, p*p) for p in (3, 7, 13))
    require(c_upper < F(15, 7), "convolution factor")
    character_sum = sum(F(1 if k % 3 == 1 else -1, k) for k in range(1, 75) if k % 3)
    # The last positive term at 73 must be followed by its negative term at 74.
    require(character_sum > F(3, 5), "ternary character threshold")
    e_lower = F(27, 10)
    require(e_lower**15 > 156 and 13**6*74 < e_lower**21, "sieve level")
    small_term = F(23**2, 4)/e_lower**15
    require(small_term < F(1, 100), "three small-prime terms")
    close_to_one = 1 + F(1, 10**6)
    require(1+1/(e_lower**30-1) < close_to_one, "partial summation endpoint")
    require(1+F(2, 3_000_000-3) < close_to_one, "Euler exception")
    main = (1+F(1, 144))*F(200, 7)
    require((5*main+F(1, 100))*close_to_one**2 < 144, "general tail constant")
    require((close_to_one*main+F(1, 100))*close_to_one**2 < 29, "restrained tail constant")
    require(F(144*30, 30-7) < 188 and F(3*30, 2) < 188, "all-endpoint constant")
    # e^(23/10) < 10 by a positive Taylor sum and a geometric remainder.
    term, partial = F(1), F(1)
    for k in range(1, 17):
        term *= F(23, 10*k)
        partial += term
    require(partial+term*F(23, 170)/(1-F(23, 180)) < 10, "log 10")
    # For log(3/2) and log(2) the atanh-series arguments are 1/5 and 1/3.
    log_three_halves = 2*(F(1, 5) + F(1, 5**3*3)/(1-F(1, 25)))
    log_two_lower = 2*(F(1, 3) + F(1, 3**3*3) + F(1, 3**5*5))
    b, z = F(552), F(828)
    prime_error = F(2, 840)*(1/b**2 + 1/z**2 + 1/b - 1/z)
    require(log_three_halves < F(406, 1000) and
            close_to_one*(F(406, 1000) + prime_error) < F(407, 1000), "prime interval")
    require(log_two_lower-F(144, 552-7)-F(407, 1000) > F(1, 50), "high exponent count")
    x = 2*10**12
    # Use the larger of the two published upper bounds for ordinary primes <= x.
    finite = F(73_802_690_054_260_520, 10**18)
    if prefix is not None:
        require((prefix["upper"], prefix["cutoff"], prefix["scale"]) == (x, 3_000_000, 10**18),
                "invalid prefix range")
        require(prefix["survivors"] == 8_391_832_364 and
                prefix["scaled_reciprocal_upper"] == 72_401_332_985_083_531 and
                F(prefix["scaled_reciprocal_upper"], prefix["scale"]) <= finite, "invalid prefix sum")
    intervals = []
    for exponent, count in [(75, 70_000_000_000), (100, 160_000_000_000)]:
        upper = finite + F(count, x) + F(29, exponent*F(23, 10)-7) + F(1, 2_999_999)
        require(upper < F(31, 108), "interval does not force the claimed count")
        intervals.append((exponent, count, upper))
    # Four alternating terms in log(1+1/3) end on a lower bound.
    require(sum((-1)**(k+1)*F(1, k*3**k) for k in range(1, 5)) == F(31, 108),
            "lower bound for log(4/3)")
    return intervals


if __name__ == "__main__":
    import argparse
    import json
    from pathlib import Path

    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prefix", type=Path, help="also check the independently sieved prefix record")
    args = cli.parse_args()
    prefix = json.loads(args.prefix.read_text(encoding="utf-8")) if args.prefix else None
    for exponent, count, upper in verify(prefix):
        print(f"(2*10^12, 10^{exponent}] requires more than {count:,} ordinary primes"
              f" (log upper {upper} < 31/108)")
