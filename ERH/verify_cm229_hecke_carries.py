"""Cubic Hecke carry sums, the prime-square drift and a three-frequency test.

Requires python-flint 0.9.0, cm229-zeros64.csv and the preceding CM229 verifiers.
"""

import csv
from itertools import product
from math import isqrt
from pathlib import Path

from flint import acb, arb, ctx, fmpz

from verify_cm229_binomial import STARTS, local_data, trace_coefficient
from verify_cm229_optimal_map import add_carry, binomial_order
from verify_cm229_split_primes import odd_prime_sieve


CUTOFFS = ((1, 2), (4, 7), (2, 3))
LOWER_MILLIONTHS = ((320_451, 304_625, 329_555, 692_335),
                    (1_074_672, 334_594, 1_203_290, 657_061),
                    (2_209_373, 315_753, 2_493_294, 621_041))


def hecke_trace(e, f_m, f_b, j):
    """Sum the degree-weighted cubic character over B-primes above p."""
    if j % f_b:
        return 0
    if f_m == f_b or (j // f_b) % 3 == 0:
        return 4 // e
    return -(2 // e)


def add_prime_above_cutoff(changes, p, a, b, weight):
    """Restrict the carry to (2n)^a < p^b, including its last downward jump."""
    stop = len(changes)
    power = p**b
    if power <= (2 * (stop - 1))**a:
        stop = int(fmpz(power - 1).root(a)) // 2 + 1
    add_carry(changes, p, 0, stop, weight)
    if stop < len(changes) and (stop - 1) % p >= (p + 1) // 2:
        changes[stop] -= weight


def changes_and_checks(stop, samples):
    tails = [[arb(0) for _ in range(stop)] for _ in CUTOFFS]
    full = [arb(0) for _ in range(stop)]
    checks = {n: [arb(0) for _ in range(4)] for n in samples}
    limit = 2 * (stop - 1)
    sieve = odd_prime_sieve(limit)
    primes = [2] + [p for p in range(3, limit + 1, 2) if sieve[p // 2]]
    for p in primes:
        e, f_m, f_b, order = local_data(p)
        logp = arb(p).log()
        # Chebotarev weight from the splitting type, independent of the Hecke trace.
        w = 3 * (e == 1 and f_m == 1) - (e == 1 and f_b == 1)
        if e == 1 and hecke_trace(e, f_m, f_b, 1) != 2 * w:
            raise ArithmeticError(f"incorrect degree-one splitting identity at {p}")
        if w:
            for changes, (a, b) in zip(tails, CUTOFFS):
                add_prime_above_cutoff(changes, p, a, b, w * logp)
        q, j = p, 1
        while q <= limit:
            trace = hecke_trace(e, f_m, f_b, j)
            if trace != trace_coefficient(p, j, order):
                raise ArithmeticError(f"Hecke and Artin traces differ at {p}^{j}")
            if trace:
                add_carry(full, q, 0, stop, (trace * logp) / 2)
            q *= p
            j += 1
        for n, values in checks.items():
            if p > 2 * n:
                continue
            carry = (2 * n) // p - 2 * (n // p)
            if w and carry:
                for i, (a, b) in enumerate(CUTOFFS):
                    if (2 * n)**a < p**b:
                        values[i] += w * logp
            v = 3 * binomial_order(n, p**f_m) - binomial_order(n, p**f_b)
            values[-1] += (v * logp) / e
    return tails + [full], checks


def window_masses(kappa):
    windows = [(y, isqrt(y**3 - 1) + 1) for y in STARTS]
    samples = {n for a, b in windows for n in (a, b - 1)}
    stop = max(b for _, b in windows)
    changes, checks = changes_and_checks(stop, samples)
    values = [arb(0) for _ in changes]
    masses = [[arb(0) for _ in range(4)] for _ in windows]
    budgets = [y**3 for y, _ in windows]
    for n in range(1, stop):
        for i, delta in enumerate(changes):
            values[i] += delta[n]
        if n in checks and not all(x.overlaps(y) for x, y in zip(values, checks[n])):
            raise ArithmeticError(f"finite product and carry jumps disagree at {n}")
        if n < windows[0][0]:
            continue
        raw, four_sevenths, two_thirds, log_r = values
        data = (raw.max(0), abs(raw - four_sevenths), two_thirds.max(0),
                abs(raw - log_r + kappa * arb(n).sqrt()))
        for (left, right), w, row in zip(windows, budgets, masses):
            if left <= n < right:
                k = isqrt(w * w // n**3)
                for i, value in enumerate(data):
                    row[i] += k * value
    return [(y, right, [x / w for x in row])
            for (y, right), w, row in zip(windows, budgets, masses)]


def three_frequencies(kappa):
    with Path(__file__).with_name("cm229-zeros64.csv").open(encoding="utf-8", newline="") as f:
        keys = ("lower_numerator", "upper_numerator", "denominator")
        rows = sorted(tuple(int(r[k]) for k in keys) for r in csv.DictReader(f))
    den = 2**24
    if len(rows) != 192 or any(d != den for _, _, d in rows):
        raise ArithmeticError("unexpected zero certificate grid")
    chosen = rows[:3]
    gap_zero, gap_poles, max_height = 64 * den, 64 * den, 0
    for coeffs in product(range(-2, 3), range(-1, 2), range(-1, 2)):
        if sum(map(abs, coeffs)) < 2:
            continue
        lo = sum(min(v * a, v * b) for v, (a, b, _) in zip(coeffs, chosen))
        hi = sum(max(v * a, v * b) for v, (a, b, _) in zip(coeffs, chosen))
        if hi < 0:
            lo, hi = -hi, -lo
        gap_zero = min(gap_zero, lo)
        gap_poles = min(gap_poles, min(max(lo - b, a - hi) for a, b, _ in rows))
        max_height = max(max_height, hi)
    if (gap_zero, gap_poles, max_height) != (6_222_089, 487_800, 107_284_867):
        raise ArithmeticError("a frequency or its distance to a pole differs from the certificate")
    sizes = []
    for a, b, _ in chosen:
        rho = acb(arb(1) / 2, arb(a + b, b - a) / (2 * den))
        sizes.append(abs((2 - (rho * arb(2).log()).exp()) * rho.zeta() / rho) / 2)
    bound = (arb(4) / 3 * sizes[0] + sizes[1] + sizes[2] - kappa) * 3 / 64
    if not arb("0.0160647") < bound < arb("0.0160652"):
        raise ArithmeticError(f"positive mass bound differs from the certificate: {bound}")
    print("positive mass lower bound", bound.str(12), "gap to a zero in grid units", gap_poles)


def main():
    ctx.prec = 160
    kappa = (arb(2).sqrt() - 2) * (arb(1) / 2).zeta()
    if not arb("0.85545586538") < kappa < arb("0.85545586540"):
        raise ArithmeticError(f"incorrect square-prime constant: {kappa}")
    print("square-prime drift", kappa.str(18))
    three_frequencies(kappa)
    for (left, right, row), bounds in zip(window_masses(kappa), LOWER_MILLIONTHS, strict=True):
        for value, lo in zip(row, bounds, strict=True):
            if not lo < value * 10**6 < lo + 1:
                raise ArithmeticError(f"window {left} differs from its published enclosure")
        print(f"[{left}, {right})  raw positive, middle absolute, top positive, correction error:",
              *(x.str(18) for x in row))


if __name__ == "__main__":
    main()
