"""Check positive masses on complete power windows. Requires python-flint 0.9.0
and verify_cm229_split_primes.py in the same directory.
"""

from math import isqrt, lcm

from flint import arb, ctx, nmod_poly

from verify_cm229_split_primes import odd_prime_sieve, splits_cubic


STARTS = (1024, 4096, 10_000)
INDEX_BOUNDS = ((1955656, 1955657), (5191168, 5191169), (8864300, 8864301))


def local_data(p):
    """Return e, f_M, f_B and the cubic Frobenius order (0 at ramification)."""
    if p == 2:
        return 2, 2, 2, 0
    if p == 229:
        return 2, 1, 1, 0
    quadratic = 1 if pow(p, 114, 229) == 1 else 2
    gaussian = 1 if p % 4 == 1 else 2
    relative_order = 2 if quadratic == 2 else (1 if splits_cubic(p) else 3)
    # The trace calculation uses a separate factorization of f modulo p.
    _, factors = nmod_poly([1, -4, 0, 1], p).factor()
    order = max(poly.degree() for poly, _ in factors)
    return 1, lcm(relative_order, gaussian), lcm(quadratic, gaussian), order


def trace_coefficient(p, j, order):
    if p == 2:
        return 2 if j % 2 == 0 else 0
    if p == 229:
        return 2
    if j % order == 0:
        trace = 2
    else:
        trace = 0 if order == 2 else -1
    twist = 1 if p % 4 == 1 or j % 2 == 0 else -1
    return (1 + twist) * trace


def prime_data(limit):
    sieve = odd_prime_sieve(limit)
    primes = [2] + [p for p in range(3, limit + 1, 2) if sieve[p // 2]]
    increments = [arb(0) for _ in range(limit + 1)]
    prime_powers = []
    split_primes = []
    for p in primes:
        e, f_m, f_b, order = local_data(p)
        logp = arb(p).log()
        if e == 1 and f_m == 1:
            split_primes.append((p, logp))
        power, j = p, 1
        while power <= limit:
            degree_difference = (12 // e) * (j % f_m == 0) - (4 // e) * (j % f_b == 0)
            trace = trace_coefficient(p, j, order)
            if degree_difference != 2 * trace:
                raise ArithmeticError(f"prime-ideal and trace data disagree at {p}^{j}")
            if trace:
                coefficient = trace * logp
                prime_powers.append((power, degree_difference * logp / 2))
                for n in range(power, limit + 1, power):
                    increments[n] += coefficient
            power *= p
            j += 1
    return increments, prime_powers, split_primes


def index_at(n, prime_powers):
    """Direct log of the normalized quotient of binomial ideal norms."""
    return sum((weight * ((2 * n) // q - 2 * (n // q))
                for q, weight in prime_powers if q <= 2 * n), arb(0))


def binomial_masses(increments, prime_powers, windows):
    ends = {t for window in windows for t in window}
    masses = [arb(0) for _ in windows]
    value = arb(0)
    for n in range(1, max(ends) + 1):
        value += increments[2 * n - 1] + increments[2 * n] - 2 * increments[n]
        if n in ends and not value.overlaps(index_at(n, prime_powers)):
            raise ArithmeticError(f"direct factorial valuations disagree at {n}")
        for j, (left, right) in enumerate(windows):
            if left <= n < right:
                masses[j] += value.max(0) / (n * arb(n).sqrt())
    return masses


def positive_plateau(theta, lo, hi):
    """Integral of [theta/4 - t/48]_+ t^(-3/2) on one plateau."""
    if 12 * theta <= lo:
        return arb(0)
    cutoff = (12 * theta).max(lo).min(hi)
    a, b = arb(lo).sqrt(), cutoff.sqrt()
    return (b - a) * (theta / (2 * a * b) - arb(1) / 24)


def split_mass(lo, hi, split_primes):
    theta, total, start = arb(0), arb(0), lo
    for p, logp in split_primes:
        if p >= hi:
            break
        if p > lo:
            total += positive_plateau(theta, start, p)
            start = p
        theta += logp
    return total + positive_plateau(theta, start, hi)


def main():
    ctx.prec = 160
    windows = [(y, y * isqrt(y)) for y in STARTS]
    increments, powers, primes = prime_data(2 * max(b for _, b in windows))
    masses = binomial_masses(increments, powers, windows)
    for (a, b), mass, (low, high) in zip(windows, masses, INDEX_BOUNDS):
        fixed = split_mass(a, b, primes)
        if not low < mass * 10**6 < high:
            raise ArithmeticError(f"binomial mass outside the published interval at {a}")
        if a == 10_000:
            if not 855_272 < fixed * 10**8 < 855_273:
                raise ArithmeticError("positive prime mass outside the published interval")
        elif fixed != 0:
            raise ArithmeticError(f"the fixed prime integral is not exactly zero at {a}")
        print(f"[{a}, {b})  binomial {mass.str(22)}, fixed prime integral {fixed.str(22)}")


if __name__ == "__main__":
    main()
