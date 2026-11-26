"""Certify archimedean costs for three complete windows (python-flint 0.9.0).

Keep verify_cm229_binomial.py, verify_cm229_split_primes.py and cm229-zeros64.csv
in this directory. The zero isolation itself is checked by verify_cm229.py.
"""

import csv
from collections import defaultdict
from fractions import Fraction
from math import isqrt
from pathlib import Path

from flint import acb, arb, ctx

from verify_cm229_binomial import STARTS, local_data, trace_coefficient
from verify_cm229_split_primes import odd_prime_sieve


LOWER_MILLIONTHS = ((1_953_553, 11_287_457, 131_589_139),
                    (5_190_081, 19_506_001, 412_129_950),
                    (8_863_359, 27_461_267, 835_188_995))


def mixed_degrees(p, order):
    if p == 2:
        return 2, 2, 2
    if p == 229:
        return 2, 1, 1
    image = {1: (0, 1, 2), 2: (1, 0, 2), 3: (1, 2, 0)}[order]
    gaussian_inert = p % 4 == 3
    odd = pow(p, 114, 229) != 1
    # The minority sign is fixed at real embedding 0.
    a_number = sum((r_image == 0) != ((r == 0) ^ gaussian_inert)
                   for r, r_image in enumerate(image))
    splits = odd or a_number == (3 if gaussian_inert else 0)
    f_b = 2 if odd or gaussian_inert else 1
    return 1, f_b if splits else 3 * f_b, f_b


def add_carry(changes, q, start, stop, weight):
    """Add the jumps of weight * 1_{n mod q >= ceil(q/2)} on [start, stop)."""
    half = (q + 1) // 2
    if start % q >= half:
        changes[start] += weight
    after = start + 1
    for n in range(after + (half - after) % q, stop, q):
        changes[n] += weight
    for n in range(after + (-after) % q, stop, q):
        changes[n] -= weight


def binomial_order(n, q):
    """The difference of the two generalized factorial orders."""
    a, b, value = 2 * n, n, 0
    while a >= q:
        a //= q
        b //= q
        value += a - 2 * b
    return value


def carry_changes(stop, samples):
    signed = [arb(0) for _ in range(stop)]
    numerator = [arb(0) for _ in range(stop)]
    surcharge = [arb(0) for _ in range(stop)]
    exact = {n: Fraction(1) for n in samples}
    multiplier = {n: 1 for n in samples}
    sample_limit = 2 * max(samples)
    limit = 2 * (stop - 1)
    sieve = odd_prime_sieve(limit)
    primes = [2] + [p for p in range(3, limit + 1, 2) if sieve[p // 2]]
    for p in primes:
        e, f_m, f_b, order = local_data(p)
        if mixed_degrees(p, order) != (e, f_m, f_b):
            raise ArithmeticError(f"mixed CM test disagrees with residue degrees at {p}")
        if p <= sample_limit:
            for n in samples:
                v, rem = divmod(3 * binomial_order(n, p**f_m) - binomial_order(n, p**f_b), e)
                if rem:
                    raise ArithmeticError(f"nonrational sample valuation at n={n}, p={p}")
                exact[n] *= Fraction(p) ** v
                if p * p <= 2 * n:
                    loss = 2 * binomial_order(n, p) - v
                    if loss < 0:
                        raise ArithmeticError(f"nonintegral multiplier at n={n}, p={p}")
                    multiplier[n] *= p**loss
        logp = arb(p).log()
        events = defaultdict(int)
        cutoff = (p * p + 1) // 2
        q, j = p, 1
        while q <= limit:
            weight, remainder = divmod(3 * (j % f_m == 0) - (j % f_b == 0), e)
            if remainder or 2 * weight != trace_coefficient(p, j, order):
                raise ArithmeticError(f"rational valuation or Artin trace disagrees at {p}^{j}")
            if weight:
                add_carry(events, q, 0, stop, weight)
            if cutoff < stop and weight != 2:
                add_carry(surcharge, q, cutoff, stop, (2 - weight) * logp)
            q *= p
            j += 1
        previous = 0
        for n, change in sorted(events.items()):
            current = previous + change
            if change:
                signed[n] += change * logp
            gain = max(current, 0) - max(previous, 0)
            if gain:
                numerator[n] += gain * logp
            previous = current
    return signed, numerator, surcharge, exact, multiplier


def costs(windows, changes):
    signed, numerator, surcharge, exact, multiplier = changes
    weights = [y**3 for y, _ in windows]
    totals = [[arb(0) for _ in range(3)] for _ in windows]
    log_r, log_u, extra = arb(0), arb(0), arb(0)
    for n in range(1, len(signed)):
        log_r += signed[n]
        log_u += numerator[n]
        extra += surcharge[n]
        if n in exact:
            r = exact[n]
            direct = arb(r.numerator).log() - arb(r.denominator).log()
            if not log_r.overlaps(direct) or not log_u.overlaps(arb(r.numerator).log()):
                raise ArithmeticError(f"direct fraction and running product disagree at {n}")
            if not extra.overlaps(arb(multiplier[n]).log()):
                raise ArithmeticError(f"large-prime multiplier disagrees at {n}")
        for (left, right), w, total in zip(windows, weights, totals):
            if left <= n < right:
                k = isqrt(w * w // n**3)
                total[0] += k * log_r.max(0)
                total[1] += k * (log_r + extra).max(0)
                total[2] += k * log_u
    return [[2 * value / w for value in row] for row, w in zip(totals, weights)]


def check_residue():
    keys = ("lower_numerator", "upper_numerator", "denominator")
    with Path(__file__).with_name("cm229-zeros64.csv").open(encoding="utf-8", newline="") as stream:
        intervals = [tuple(int(row[key]) for key in keys) for row in csv.DictReader(stream)]
    lo, hi, den = min(intervals)
    if (lo, hi, den) != (13_347_347, 13_347_348, 2**24) or len(intervals) != 192:
        raise ArithmeticError("unexpected first zero or zero count")
    degree = 80
    if degree * hi >= 64 * den or any(d != den for _, _, d in intervals):
        raise ArithmeticError("harmonic tests outside the certified height or on a different grid")
    gap = min(max(k * lo - b, a - k * hi)
              for k in range(2, degree + 1) for a, b, _ in intervals)
    if gap != 13_542:
        raise ArithmeticError(f"unexpected separation for the first 80 harmonics: {gap}")
    height = arb(lo + hi, hi - lo) / (2 * den)
    rho = acb(arb(1) / 2, height)
    residue = -((rho * arb(2).log()).exp() - 2) * rho.zeta() / rho
    half_size = abs(residue) / 2
    if not arb("0.48564") < half_size < arb("0.48565"):
        raise ArithmeticError(f"first-pole residue outside its certified interval: {half_size}")
    print("Half modulus of the first residue:", half_size.str(13))
    lower_test = 4 * degree * half_size / ((degree + 1) * arb.pi())
    if not arb("0.61070") < lower_test < arb("0.61072"):
        raise ArithmeticError(f"Fejer bound outside its published interval: {lower_test}")
    print("Fejer lower bound:", lower_test.str(12), "harmonic gap in grid units:", gap)


def main():
    ctx.prec = 160
    check_residue()
    windows = [(y, isqrt(y**3 - 1) + 1) for y in STARTS]
    samples = (5, 15, 20, 37, 115, 229, 997, 1024)
    changes = carry_changes(max(b for _, b in windows), samples)
    examples = {5: Fraction(1, 3), 15: Fraction(4, 255), 20: Fraction(37**2, 5)}
    for n, r in examples.items():
        if changes[3][n] != r:
            raise ArithmeticError(f"sample generator disagrees at {n}")
        print(f"r({n}) = {r}")
    for (left, right), row, lower in zip(windows, costs(windows, changes), LOWER_MILLIONTHS):
        for value, lo in zip(row, lower):
            if not lo < value * 10**6 < lo + 1:
                raise ArithmeticError(f"cost outside its published enclosure at {left}")
        print(f"[{left}, {right})  min, large primes only, integral coordinates (per 2W):",
              *(value.str(20) for value in row))


if __name__ == "__main__":
    main()
