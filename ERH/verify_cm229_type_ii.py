"""Check the cutoffs in the cubic Type II sum on the whole [4096,262144) window.

Requires python-flint 0.9.0 and the companion CM229 programs in the same directory.
"""

import csv
from array import array
from math import isqrt
from pathlib import Path

from flint import arb, ctx, fmpz

from verify_cm229 import coefficients
from verify_cm229_bilinear import CHI4, bilinear_prefix, local_traces, prefix_sums
from verify_cm229_split_primes import odd_prime_sieve


LOWER_MILLIONTHS = ((868_852, 745_264, 42_464, 1_792_361, 913_569),
                    (653_600, 276_567, 79_605, 911_352, 1_385_341))
PRIME_MILLIONTHS = ((606_624, 532_182, 1_129_625), (566_701, 390_309, 852_363))


def prime_data(a, limit):
    """Euler derivatives and, independently, the largest prime at every product."""
    sieve = odd_prime_sieve(limit)
    largest = array("I", [0]) * (limit + 1)
    prime_logs = {}
    powers = []
    cubes = []
    for p in range(2, limit + 1):
        if p != 2 and (p % 2 == 0 or not sieve[p // 2]):
            continue
        for t in range(p, limit + 1, p):
            largest[t] = p
        logp = arb(p).log()
        prime_logs[p] = logp
        for q, j, trace in local_traces(a, p, limit):
            if trace:
                pair = (q, trace * logp)
                powers.append(pair)
                if j >= 3:
                    cubes.append(pair)
    return powers, cubes, largest, prime_logs


def largest_prime_prefix(a, largest, prime_logs, block, limit):
    cutoff = isqrt(4 * block)
    values = [arb(0) for _ in range(limit + 1)]
    for q in range(4, limit + 1):
        p = largest[q]
        if p <= cutoff or q == p:
            continue
        h = q // p
        if h >= p or a[h] * a[p] != a[q]:
            raise ArithmeticError(f"multiplicativity or largest prime fails at {q} in {block}")
        values[q] = -(1 + CHI4[q % 4]) * a[q] * prime_logs[p] / 2
    prefix_sums(values)
    return values


def shared_prime_prefix(a, prime_logs, z, v, limit):
    """Compute the p | h part from the Euler polynomial in (5) of the prime-range note."""
    max_cofactor = limit // (v + 1)**2
    inverse = [0, 1] + [0] * (z - 1)
    for d in range(2, z + 1):
        inverse[d] = -sum(inverse[j] * a[d//j] for j in range(1, d) if d % j == 0)
    truncated_product = [0]
    for j in range(1, max_cofactor + 1):
        truncated_product.append(sum(inverse[d] * a[j//d] for d in range(1, min(z, j) + 1)
                                     if j % d == 0))
    values = [arb(0) for _ in range(limit + 1)]
    for p, logp in prime_logs.items():
        if p*p > limit:
            break
        if p <= v:
            continue
        determinant = a[p]**2 - a[p*p]
        for q, coefficient in ((p*p, -a[p]**2), (p**3, a[p] * determinant)):
            for j in range(1, limit // q + 1):
                t = j*q
                weight = (1 + CHI4[t % 4]) * coefficient * truncated_product[j]
                if weight:
                    values[t] += weight * logp / 2
    prefix_sums(values)
    return values


def check_window():
    y = 4096
    budget, stop = y**3, isqrt(y**3 - 1) + 1
    a = coefficients(2 * (stop - 1))["real"]
    powers, cubes, largest, logs = prime_data(a, len(a) - 1)
    primes = [(p, a[p] * logp) for p, logp in logs.items() if a[p]]
    m = int((arb(2*y).log()**2).ceil().unique_fmpz())
    r = (arb(2).log() / m).exp()
    masses = [[arb(0) for _ in range(5)] for _ in range(2)]
    prime_masses = [[arb(0) for _ in range(3)] for _ in range(2)]
    cutoffs = []
    block = y
    while block < stop:
        end = min(2 * block, stop)
        limit = 2 * (end - 1)
        sixth = int(fmpz(block).root(6))
        analytic_z = arb(block) ** (arb(1) / 6) / (arb(2*block).log()**3)
        original = max(1, int(analytic_z.floor().unique_fmpz()))
        balanced = int(fmpz(block).root(4))
        cutoffs.append((block, original, sixth, balanced))
        largest_model = largest_prime_prefix(a, largest, logs, block, limit)
        prime_model = bilinear_prefix(a, primes, 1, limit, ell_min=isqrt(4*block))
        if not all(x.overlaps(p) for x, p in zip(largest_model, prime_model, strict=True)):
            raise ArithmeticError(f"prime-only convolution and largest prime disagree in {block}")
        fourth_model = bilinear_prefix(a, powers, balanced, limit)
        for row, prime_row, z in zip(masses, prime_masses, (original, sixth), strict=True):
            v = isqrt(block) // z
            base = bilinear_prefix(a, powers, z, limit)
            strip = bilinear_prefix(a, [(q, b) for q, b in powers if q <= v], z, limit)
            high = bilinear_prefix(a, cubes, z, limit, ell_min=v)
            # Compute the surviving prime/square events independently of this subtraction.
            squares = []
            for p, logp in logs.items():
                if p*p > limit:
                    break
                if p*p > v:
                    squares.append((p*p, (2*a[p*p] - a[p]**2) * logp))
            square_sum = bilinear_prefix(a, squares, z, limit, ell_min=v)
            coprime = bilinear_prefix(a, primes, z, limit, ell_min=v, coprime=True)
            shared = shared_prime_prefix(a, logs, z, v, limit)
            partition = zip(base, strip, high, square_sum, shared, coprime, strict=True)
            if not all((b - s - h - sq - div).overlaps(c) for b, s, h, sq, div, c in partition):
                raise ArithmeticError(f"coprime-prime partition fails in {block} at cutoff {z}")
            extreme = bilinear_prefix(a, primes, z, limit, ell_min=int(fmpz(block**2).root(3)))
            for n in range(block, end):
                right = 2 * n
                left = int((right / r).floor().unique_fmpz())
                raw = base[right] - base[left]
                parts = (raw.max(0), abs(strip[right] - strip[left]), abs(high[right] - high[left]),
                         abs(raw - (fourth_model[right] - fourth_model[left])),
                         abs(raw - (largest_model[right] - largest_model[left])))
                k = isqrt(budget**2 // n**3)
                for j, value in enumerate(parts):
                    row[j] += k * value
                coprime_band = coprime[right] - coprime[left]
                extreme_band = extreme[right] - extreme[left]
                for j, value in enumerate((coprime_band.max(0), extreme_band.max(0),
                                           abs(raw - extreme_band))):
                    prime_row[j] += k * value
        block = end
    print("blocks, original z, sixth root and fourth root:", cutoffs)
    if m != 82 or any(z != 1 for _, z, _, _ in cutoffs):
        raise ArithmeticError("unexpected original window cutoffs")
    for name, row, lower in zip(("original", "sixth root"), masses, LOWER_MILLIONTHS, strict=True):
        row = [x / budget for x in row]
        if not all(lo < x * 10**6 < lo + 1 for x, lo in zip(row, lower, strict=True)):
            raise ArithmeticError(f"{name} window differs from the published intervals")
        print(name, "positive, strip, high powers, fourth-root difference, largest-prime difference:",
              *(x.str(18) for x in row))
    for name, row, lower in zip(("original", "sixth root"), prime_masses, PRIME_MILLIONTHS, strict=True):
        row = [x / budget for x in row]
        if not all(lo < x * 10**6 < lo + 1 for x, lo in zip(row, lower, strict=True)):
            raise ArithmeticError(f"{name} prime ranges differ from their published intervals")
        print(name, "coprime positive, prime above two-thirds positive, original-minus-latter absolute:",
              *(x.str(18) for x in row))


def frequency_check():
    with Path(__file__).with_name("cm229-zeros64.csv").open(encoding="utf-8", newline="") as f:
        rows = sorted((int(r["lower_numerator"]), int(r["upper_numerator"]), int(r["denominator"]))
                      for r in csv.DictReader(f))
    if len(rows) != 192 or any(d != 2**24 for _, _, d in rows):
        raise ArithmeticError("unexpected low-zero data")
    a, b, den = rows[0]
    gap = min(max(j*a - hi, lo - j*b) for j in range(2, 81) for lo, hi, _ in rows)
    if (a, b, gap) != (13_347_347, 13_347_348, 13_542) or 80*b >= 64*den:
        raise ArithmeticError("a harmonic intersects a zero or leaves the certified range")
    constant = 80 * arb(2).log() / (81 * arb.pi() * arb(2).sqrt())
    if not arb("0.15408684349") < constant < arb("0.15408684350"):
        raise ArithmeticError("the low-zero coefficient differs from its published interval")
    print("first-zero harmonic gap in grid units", gap, "positive-mass liminf bound", constant.str(18))


if __name__ == "__main__":
    ctx.prec = 160
    frequency_check()
    check_window()
