"""Check the unsmoothed splitting-prime sum for the discriminant-229 example.

Requires python-flint 0.9.0. Run with python verify_cm229_split_primes.py.
"""

from math import isqrt
from operator import gt, lt

from flint import arb, ctx, nmod_poly


LIMIT = 100_000_000


def odd_prime_sieve(limit):
    """Entry j describes the integer 2*j+1."""
    prime = bytearray(b"\x01") * ((limit + 1) // 2)
    prime[0] = 0
    for p in range(3, isqrt(limit) + 1, 2):
        if prime[p // 2]:
            start = p * p // 2
            count = 1 + (len(prime) - 1 - start) // p
            prime[start::p] = bytes(count)
    return prime


def principal_norms(limit):
    """Mark integers p = 1 mod 4 in the reduced norm sector."""
    represented = bytearray(limit // 4 + 1)
    for b in range(isqrt(225 * limit // 229) + 1):
        offset = 229 * b * b
        least_p = max(2, (offset + 224) // 225)
        lower = isqrt(offset + 4 * least_p - 1) + 1
        upper = isqrt(offset + 4 * limit)
        for residue in range(8):
            if (residue * residue - offset) % 16 != 4:
                continue
            first = lower + (residue - lower) % 8
            for m in range(first, upper + 1, 8):
                p = (m * m - offset) // 4
                represented[p // 4] = 1
    return represented


def splits_cubic(p):
    x = nmod_poly([0, 1], p)
    cubic = nmod_poly([1, -4, 0, 1], p)
    return x.pow_mod(p, cubic) == x


def checked_extreme(best, value, point, improves):
    old = best[0]
    if improves(value, old):
        return value, point
    if improves(old, value):
        return best
    raise ArithmeticError(f"extremum comparison unresolved at {point}")


def run():
    ctx.prec = 128
    prime = odd_prime_sieve(LIMIT)
    represented = principal_norms(LIMIT)
    theta = arb(0)
    count = candidates = 0
    maximum = minimum = (-arb(2).sqrt() / 12, 2)
    first_above = last_above = None
    last_crossing = None
    for p in range(5, LIMIT + 1, 4):
        if p == 229 or not prime[p // 2]:
            continue
        candidates += 1
        split = splits_cubic(p)
        if split != bool(represented[p // 4]):
            raise ArithmeticError(f"norm and Frobenius tests disagree at {p}")
        if not split:
            continue
        count += 1
        main = arb(p) / 12
        root = arb(p).sqrt()
        minimum = checked_extreme(minimum, (theta - main) / root, p, lt)
        theta += arb(p).log()
        maximum = checked_extreme(maximum, (theta - main) / root, p, gt)
        if theta > main:
            first_above = first_above or p
            last_above = p
            last_crossing = 12 * theta
        elif not theta < main:
            raise ArithmeticError(f"sign unresolved at {p}")
    endpoint = (theta - arb(LIMIT) / 12) / arb(LIMIT).sqrt()
    minimum = checked_extreme(minimum, endpoint, LIMIT, lt)
    if not (candidates == 2_880_503 and count == 479_605
            and (first_above, last_above) == (37, 7_207_729)
            and 7_207_871 < last_crossing < 7_207_872
            and maximum[1] == 873_157 and minimum[1] == 55_511_413
            and 4_619_508 < 10**7 * maximum[0] < 4_619_509
            and -21_541_849 < 10**7 * minimum[0] < -21_541_848):
        raise ArithmeticError("results differ from the published enclosures")
    print(f"Range [2, {LIMIT}]. Tested {candidates} primes, {count} split.")
    print(f"First and last above-main jumps: {first_above}, {last_above}")
    print(f"Final crossing 12*theta(last above-main jump): {last_crossing}")
    print(f"Maximum {maximum[0]} at {maximum[1]}")
    print(f"Infimum {minimum[0]} at the endpoint or left limit of {minimum[1]}")
    print(f"Normalized error at the upper endpoint: {endpoint}")


if __name__ == "__main__":
    run()
