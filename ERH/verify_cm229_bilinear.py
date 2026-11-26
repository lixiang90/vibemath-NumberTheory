"""Verify Vaughan's identity through 4096 and evaluate its bilinear remainder.

Requires python-flint 0.9.0 and the companion files verify_cm229.py,
verify_tetra3136.py, verify_cm229_binomial.py and verify_cm229_split_primes.py.
"""

from math import isqrt

from flint import arb, ctx, fmpz

from verify_cm229 import coefficients
from verify_cm229_binomial import local_data, trace_coefficient
from verify_cm229_split_primes import odd_prime_sieve


CHI4 = (0, 1, 0, -1)


def local_traces(a, p, limit):
    """Return b_+(p^j)/log(p) from the logarithmic derivative of the Euler series."""
    previous = []
    q, j = p, 1
    while q <= limit:
        trace = j * a[q] - sum(t * a[p**(j-i)] for i, t in enumerate(previous, 1))
        yield q, j, trace
        previous.append(trace)
        q *= p
        j += 1


def verify(limit=4096):
    a = coefficients(limit)["real"]
    divisors = [[] for _ in a]
    for d in range(1, limit + 1):
        for t in range(d, limit + 1, d):
            divisors[t].append(d)

    inverse = [0] * (limit + 1)
    inverse[1] = 1
    for t in range(2, limit + 1):
        inverse[t] = -sum(a[d] * inverse[t // d] for d in divisors[t] if d > 1)
    prime_traces = {}
    traces_both = {}
    logs = [{} for _ in a]
    for p in range(2, limit + 1):
        if len(divisors[p]) != 2:
            continue
        order = local_data(p)[3]
        for q, j, value in local_traces(a, p, limit):
            prime_traces[q] = value
            traces_both[q] = trace_coefficient(p, j, order)
            for t in range(q, limit + 1, q):
                logs[t][p] = j

    for z in (2, 4, 8):
        residual = [0] * (limit + 1)
        for t in range(1, limit + 1):
            residual[t] = int(t == 1) - sum(inverse[d] * a[t // d] for d in divisors[t] if d <= z)
        if any(residual[1:z+1]):
            raise ArithmeticError(f"residual has small support for cutoff {z}")
        tests = 0
        for t in range(2, limit + 1):
            for p in logs[t]:
                part_one = sum(inverse[d] * a[t // d] * logs[t // d].get(p, 0)
                               for d in divisors[t] if d <= z)
                small, bilinear = 0, 0
                for j in range(1, logs[t][p] + 1):
                    q = p**j
                    remaining = t // q
                    if q <= z:
                        small += prime_traces[q] * sum(inverse[d] * a[remaining // d]
                                                       for d in divisors[remaining] if d <= z)
                    elif remaining > z:
                        bilinear += residual[remaining] * prime_traces[q]
                target = prime_traces.get(t, 0)
                small_b = target if t <= z else 0
                total = part_one - small + small_b + bilinear
                if total != target or (1 + CHI4[t % 4]) * total != traces_both.get(t, 0):
                    raise ArithmeticError(f"coefficient of log({p}) fails at {t}, cutoff {z}")
                tests += 1
        print(f"cutoff {z}: {tests} prime-log coefficients verified for both factors up to {limit}")


def prefix_sums(values):
    for i in range(1, len(values)):
        values[i] += values[i-1]


def bilinear_prefix(a, powers, z, limit):
    """Build (9) at product indices with both variables strictly above z."""
    inverse = [0, 1] + [0] * (z - 1)
    for d in range(2, z + 1):
        inverse[d] = -sum(a[j] * inverse[d // j] for j in range(2, d + 1) if d % j == 0)
    max_h = limit // (z + 1)
    residual = [0] * (max_h + 1)
    residual[1] = 1
    for d in range(1, z + 1):
        if inverse[d]:
            for j in range(1, max_h // d + 1):
                residual[d*j] -= inverse[d] * a[j]
    events = [arb(0) for _ in range(limit + 1)]
    for q, b in powers:
        if z < q <= max_h:
            for h in range(z + 1, limit // q + 1):
                t = h*q
                weight = (1 + CHI4[t % 4]) * residual[h]
                if weight:
                    events[t] += weight * b / 2
    prefix_sums(events)
    return events


def finite_window():
    y = 4096
    stop = isqrt(y**3 - 1) + 1
    limit = 2 * (stop - 1)
    a = coefficients(limit)["real"]
    sieve = odd_prime_sieve(limit)
    theta = [arb(0) for _ in range(limit + 1)]
    powers = []
    for p in range(2, limit + 1):
        if p != 2 and (p % 2 == 0 or not sieve[p // 2]):
            continue
        e, f_m, f_b, _ = local_data(p)
        logp = arb(p).log()
        theta[p] = (3 * (e == 1 and f_m == 1) - (e == 1 and f_b == 1)) * logp
        for q, _, trace in local_traces(a, p, limit):
            if trace:
                powers.append((q, trace * logp))
    prefix_sums(theta)

    m = int((arb(2*y).log() ** 2).ceil().unique_fmpz())
    r = (arb(2).log() / m).exp()
    square_factor = 1 - 1 / r.sqrt()
    budget = y**3
    endpoint_mass, bilinear_mass, error = arb(0), arb(0), arb(0)
    block = y
    while block < stop:
        end = min(2 * block, stop)
        z = int(fmpz(block**4).root(25))
        bilinear = bilinear_prefix(a, powers, z, 2 * (end - 1))
        for n in range(block, end):
            left = int((2*n / r).floor().unique_fmpz())
            prime_sum = theta[2*n] - theta[left]
            bilinear_sum = bilinear[2*n] - bilinear[left]
            k = isqrt(budget**2 // n**3)
            endpoint_mass += k * prime_sum.max(0)
            bilinear_mass += k * bilinear_sum.max(0)
            error += k * abs(prime_sum + arb(2*n).sqrt() * square_factor - bilinear_sum)
        block = end
    row = [x / budget for x in (endpoint_mass, bilinear_mass, error)]
    lower = (49_468_097, 67_649_279, 166_678_923)
    if m != 82 or not all(lo < x * 10**8 < lo + 1 for x, lo in zip(row, lower, strict=True)):
        raise ArithmeticError("endpoint or bilinear sum differs from its published enclosure")
    print(f"[{y}, {stop}), m={m}: endpoint positive, bilinear positive, corrected absolute error:",
          *(x.str(18) for x in row))


if __name__ == "__main__":
    ctx.prec = 160
    verify()
    finite_window()
