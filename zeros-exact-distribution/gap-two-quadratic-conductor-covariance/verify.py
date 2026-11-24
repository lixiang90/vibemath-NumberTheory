"""Two disjoint zero windows for prime-conductor real quadratic L-functions.

Run with Python, NumPy, SciPy and python-flint. The zero search uses the
Mellin transform of the quadratic theta series. The prime-side comparison,
the larger conductor sample, and the exact sign enumeration are separate.
"""

from itertools import product
from math import ceil, exp, log, pi, sqrt

import numpy as np
from flint import acb, arb, ctx, dirichlet_char
from scipy.optimize import brentq
from scipy.special import loggamma
from scipy.stats import irwinhall


B, BOX, ORDER = 3.45, 0.75, 2
RADIUS = 6.2
HEIGHT = 13 * pi / log(sqrt(17 * 19))
TOP = HEIGHT + RADIUS + 0.43
TAGS = {17: (17, 19), 53: (5, 11, 53), 143: (5, 11, 13, 29)}


def primes_up_to(limit):
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[:2] = False
    for p in range(2, int(sqrt(limit)) + 1):
        if sieve[p]:
            sieve[p * p::p] = False
    return np.flatnonzero(sieve)


def log_kernel(u):
    law = irwinhall(2 * ORDER)
    return (law.cdf((u + B) / (2 * BOX) + ORDER)
            - law.cdf((u - B) / (2 * BOX) + ORDER))


def zero_kernel(t):
    return 2 * B * np.sinc(B * t / pi) * np.sinc(BOX * t / pi) ** (2 * ORDER)


def prime_coefficient(p, h):
    base = log(p)
    return sum(base * exp(-j * base / 2) * log_kernel(j * base)
               * np.cos(j * h * base)
               for j in range(1, int((B + 2 * ORDER * BOX) / base) + 1, 2))


def character_column(ds, p):
    # The conductors are 3 modulo 4, so reciprocity has a fixed sign.
    residues = np.arange(1, (p + 1) // 2, dtype=np.int64)
    symbols = -np.ones(p, dtype=np.int8)
    symbols[0] = 0
    symbols[residues**2 % p] = 1
    return symbols[ds % p] * (1 if p % 4 == 1 else -1)


class QuadraticTheta:
    """Completed primitive odd quadratic L from its theta transformation."""

    def __init__(self, d, step=0.04):
        self.d = d
        y_max = log(105 * d / pi)
        self.y = np.linspace(0, y_max, ceil(y_max / step) + 1)
        ns = np.arange(1, ceil(sqrt(110 * d / pi)) + 1, dtype=np.int64)
        symbols = np.array([pow(int(n), (d - 1) // 2, d) for n in ns])
        signs = np.where(symbols == d - 1, -1, symbols)
        theta = np.exp(-np.outer(np.exp(self.y), pi * ns**2 / d)) @ (ns * signs)
        self.weights = 2 * (self.y[1] - self.y[0]) * np.exp(0.75 * self.y) * theta
        self.weights[[0, -1]] *= 0.5

    def completed_real(self, t):
        return np.cos(np.asarray(t)[..., None] * self.y / 2) @ self.weights

    def l_value(self, s):
        completed = np.cosh((s - 0.5) * self.y / 2) @ self.weights
        return completed * np.exp(-((s + 1) / 2) * log(self.d / pi)
                                  - loggamma((s + 1) / 2))

    def winding_count(self, panels):
        # Follow arg L from Re(s)=2, where its principal argument is continuous.
        for _ in range(5):
            xs = 0.5 + 1.5 * np.linspace(1, 0, panels + 1)**2
            phases = np.unwrap([np.angle(self.l_value(x + 1j * TOP)) for x in xs])
            if np.max(np.abs(np.diff(phases))) < 1.25:
                break
            panels *= 2
        else:
            raise ArithmeticError("horizontal phase sampling needs refinement")
        theta = TOP / 2 * log(self.d / pi) + loggamma(0.75 + 0.5j * TOP).imag
        value = 2 * (theta + phases[-1]) / pi
        assert abs(value - round(value)) < 1e-6
        return round(value)

    def zeros(self):
        expected = self.winding_count(32)
        assert self.winding_count(64) == expected
        assert expected % 2 == 0 and self.completed_real(0) > 0
        step = 0.11
        for _ in range(6):
            mesh = np.linspace(0, TOP, ceil(TOP / step) + 1)
            signs = np.signbit(self.completed_real(mesh))
            crossings = np.flatnonzero(signs[1:] != signs[:-1])
            if 2 * len(crossings) == expected:
                break
            step /= 2
        else:
            raise ArithmeticError(f"unresolved crossings for conductor {self.d}")
        roots = np.array([brentq(self.completed_real, mesh[j], mesh[j + 1],
                                 xtol=2e-12) for j in crossings])
        return np.r_[-roots[::-1], roots]


def centered(x):
    return x - np.mean(x, axis=0)


def contrasts(left, right, mark):
    # Population covariance and third cumulant, conditioned on the one product.
    stats = []
    for sign in (1, -1):
        a, b = centered(left[mark == sign]), centered(right[mark == sign])
        stats.append(np.array([-np.mean(a * b), np.mean(a**3)]))
    return (stats[0] - stats[1]) / 2


def arithmetic_columns(ds, primes, coeffs):
    values = np.zeros((len(ds), 2))
    for p, coeff in zip(primes, coeffs):
        values += character_column(ds, int(p))[:, None] * coeff
    return centered(values)


def row_mark(ds, factors):
    return np.prod([character_column(ds, p) for p in factors], axis=0)


def predicted(factors, coeff):
    if len(factors) == 2:
        p, q = factors
        return np.array([-(coeff[p][0] * coeff[q][1] + coeff[q][0] * coeff[p][1]), 0])
    if len(factors) == 3:
        return np.array([0, 6 * np.prod([coeff[p][0] for p in factors])])
    return np.zeros(2)


def exact_sign_check(coeff):
    labels = sorted({p for factors in TAGS.values() for p in factors} | {7})
    signs = np.array(list(product((-1, 1), repeat=len(labels))))
    values = signs @ np.array([coeff[p] for p in labels])
    for n, factors in TAGS.items():
        assert np.prod(factors) == n * (n + 2)
        mark = np.prod(signs[:, [labels.index(p) for p in factors]], axis=1)
        assert np.allclose(contrasts(*values.T, mark), predicted(factors, coeff),
                           rtol=0, atol=2e-13)
    print(f"exact parity checks: {len(signs)} sign vectors, {len(TAGS)} rows")


def balance_multiplier(x):
    x = np.mod(x, 1.0)
    t = (np.minimum(x, 1 - x) - 1 / 11) * 11
    cutoff = np.asarray(t >= 1, dtype=float)
    middle = (t > 0) & (t < 1)
    left = np.exp(-1 / t[middle])
    right = np.exp(-1 / (1 - t[middle]))
    cutoff[middle] = left / (left + right)
    return np.divide(cutoff, 4 * x * (1 - x), out=np.zeros_like(cutoff), where=x != 0)


def calibrated_check():
    grid_size = 2**14
    grid = np.arange(grid_size) / grid_size
    harmonics = 512
    coeff = np.fft.rfft(balance_multiplier(grid)).real[:harmonics + 1] / grid_size
    coeff[1:] *= 2
    indices = np.arange(len(coeff))

    def evaluate(xs):
        return np.cos(2 * pi * np.outer(xs, indices)) @ coeff

    ratios = np.linspace(3 / 11, 0.5, 201)
    # Direct sum of the three permutations with a frequency in the last slot.
    calibrated = (2 * ratios * (1 - ratios)
                  * (evaluate(ratios) + evaluate(1 - ratios) + evaluate([1])[0]))
    ideal_error = np.max(np.abs(calibrated - 1))
    assert ideal_error < 1e-10
    print(f"factor-balance multiplier: 201 splits, max norm error {ideal_error:.3g}")

    # Check the difference between log(n) and log(n+2) in both orientations.
    start = 10**6
    primes = primes_up_to(2 * start + 2)
    prime_set = set(primes)
    small = primes[(primes >= start**(3/11)) & (primes <= sqrt(2 * start + 2))]
    for offset in (0, 2):
        checked = []
        for p in primes[(primes > start) & (primes <= 2 * start)]:
            n = int(p) - offset
            other = n + 2 - offset
            divisors = small[(small**2 <= other) & (other % small == 0)]
            if len(divisors) != 1:
                continue
            s = int(divisors[0])
            t = other // s
            if t == s or t not in prime_set or s <= n**(3/11):
                continue
            logs = np.log([p, s, t])
            relative = (2 * np.prod(logs) * np.sum(evaluate(logs / log(n)))
                        * sqrt(n / (n + 2)) / (log(n)**2 * log(n + 2)))
            checked.append(abs(relative - 1))
            if len(checked) == 40:
                break
        assert len(checked) == 40 and max(checked) < 1e-5
        print(f"prime endpoint n+{offset}: 40 semiprime rows, max relative error {max(checked):.3g}")


def short_majorant_check():
    # Each real ball contains an entire subinterval, including its endpoints.
    panels = 4096
    lower, upper = 2.0, 0.0
    for j in range(panels):
        theta = arb(3) / 11 + (arb(j) + arb(0.5, 0.5)) * (arb(5) / (22 * panels))
        multiplier = (arb(752) / 1000
                      + arb(465) / 3000 * (1 + 2 * (2 * theta).cos_pi())
                      + arb(72) / 3000 * (1 + 2 * (4 * theta).cos_pi()))
        weight = 6 * theta * (1 - theta) * multiplier
        assert weight > arb(100044) / 100000
        assert weight < arb(100372) / 100000
        lower = min(lower, float(weight.lower()))
        upper = max(upper, float(weight.upper()))
    print(f"three heights, {panels} Arb intervals: [{lower:.8f}, {upper:.8f}] (rounded)")


def main():
    ctx.prec = 90
    assert 2 * RADIUS < HEIGHT and B >= 2 * ORDER * BOX
    primes = primes_up_to(int(exp(B + 2 * ORDER * BOX)))
    primes = primes[primes > 3]
    coeffs = np.array([[prime_coefficient(p, h) for h in (0, HEIGHT)] for p in primes])
    coeff = dict(zip(primes, coeffs))
    exact_sign_check(coeff)
    calibrated_check()
    short_majorant_check()

    all_ds = primes_up_to(12_000_000)
    all_ds = all_ds[all_ds % 48 == 23]  # d=7 (mod 16), d=-1 (mod 3)
    ds = all_ds[(all_ds > 800) & (all_ds <= 4200)]
    zero_traces, count, arb_error = [], 0, 0.0
    for i, d in enumerate(map(int, ds)):
        theta = QuadraticTheta(d)
        roots = theta.zeros()
        count += len(roots)
        if i in (0, len(ds) // 2, len(ds) - 1):
            chi = dirichlet_char(d, d - 1)
            assert chi.order() == 2 and chi.parity() and chi.is_primitive()
            for s in (1.21 + 3.25j, 0.63 + 9.37j, 1.5 + 1j * TOP):
                arb_error = max(arb_error, abs(theta.l_value(s) - complex(chi.l_function(acb(s)))))
        zero_traces.append([sum(zero_kernel(roots[np.abs(roots - h) <= RADIUS] - h))
                            for h in (0, HEIGHT)])
    assert arb_error < 3e-6
    spectral = -0.5 * centered(np.array(zero_traces) - np.log(ds[:, None]))
    arithmetic = arithmetic_columns(ds, primes, coeffs)
    trace_errors = np.max(np.abs(spectral - arithmetic), axis=0)
    assert np.all(trace_errors < 0.003)
    mark = row_mark(ds, TAGS[17])
    measured = contrasts(*spectral.T, mark)[0]
    direct = contrasts(*arithmetic.T, mark)[0]
    assert abs(measured - direct) < 0.002
    print(f"d in (800,4200]: {len(ds)}, zeros |t| <= {TOP:.5f}: {count}, "
          f"Arb max |L difference|: {arb_error:.3g}")
    print(f"h={HEIGHT:.8f}, R={RADIUS}, max trace errors: {trace_errors}")
    print(f"row 17 finite covariance contrast: zeros {measured:.8f}, arithmetic {direct:.8f}")

    ds = all_ds[all_ds > sqrt(12_000_000)]
    values = arithmetic_columns(ds, primes, coeffs)
    print(f"larger arithmetic sample: {len(ds)} conductors")
    for n, factors in TAGS.items():
        mark = row_mark(ds, factors)
        print(f"row {n}: n+={np.sum(mark == 1)}, measured {contrasts(*values.T, mark)}, "
              f"independent reference {predicted(factors, coeff)}")


if __name__ == "__main__":
    main()
