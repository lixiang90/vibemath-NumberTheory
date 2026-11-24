"""Check four height windows with primitive quadratic twists at 3q, 12q and 24q.

Run from the field directory with numpy and scipy:
    python -B gap-two-separated-zero-windows/verify.py
Add --arb to check the Dirichlet logarithmic derivatives against python-flint.
"""

import argparse
import math
import runpy
from functools import lru_cache
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
SHARED = runpy.run_path(str(ROOT / "gap-two-single-prime-conductor" / "verify.py"))
PREVIOUS = runpy.run_path(str(ROOT / "gap-two-chen-zero-pairs" / "verify.py"))
factor_weights = PREVIOUS["factor_weights"]
exact_weight = PREVIOUS["exact_weight"]
prime_weights = SHARED["prime_weights"]
primitive_root = SHARED["primitive_root"]
base_window = SHARED["window"]
spectral_shifts = SHARED["spectral_shifts"]


@lru_cache(maxsize=None)
def sharpened_shifts(order):
    base_points, base_weights = spectral_shifts(order)
    j = np.arange(1 - order, order)
    weights = (order - np.abs(j)) / order**2
    points = (base_points[None, :] + 2 * j[:, None]).ravel()
    masses = (base_weights[None, :] * weights[:, None]).ravel()
    assert abs(masses.sum() - 1) < 5e-14
    return points, masses


def transform(s, order):
    points, masses = sharpened_shifts(order)
    return (2 * np.pi / (32 * SHARED["J"])) * (
        masses @ np.exp(-np.cos((s + 1j * points) / 32))
    )


def window(u, order):
    v = math.log(u)
    fejer = abs(np.exp(2j * np.arange(order) * v).mean()) ** 2
    return base_window(u, order) * fejer


def band_radius(q, center, order):
    height_log = 2 + math.log(24 * q * (3 + abs(center)))
    budget = math.log(24 * q) + math.log(height_log)
    return 10 * order - 2 + 32 * math.log(64 * budget)


def centers_for(q, order):
    unit = 320
    while True:
        centers = [unit * j for j in (1, 3, 5, 7)]
        radii = [band_radius(q, t, order) for t in centers]
        if centers[0] > radii[0] and all(
            t + r < u - s for t, r, u, s in zip(centers, radii, centers[1:], radii[1:])
        ):
            return centers
        unit *= 2


def omega(n, discriminant):
    if discriminant == 1:
        return np.ones_like(n)
    minus_four = np.take([0, 1, 0, -1], n % 4)
    if discriminant == -4:
        return minus_four
    eight = np.take([0, 1, 0, -1, 0, -1, 0, 1], n % 8)
    return eight if discriminant == 8 else eight * minus_four


def auxiliary(n, discriminant):
    return np.take([0, 1, -1], n % 3) * omega(n, discriminant)


class BandCalculation:
    """The character-dependent part of each individual zero trace."""

    right = 4.5

    def __init__(self, q, discriminant, order, scales):
        self.q = q
        self.discriminant = discriminant
        self.modulus = 3 * abs(discriminant)
        self.order = order
        self.scales = np.array(sorted(scales), dtype=float)
        self.log_scales = np.log(self.scales)
        self.residues = primitive_root(q)
        self.indices = {int(a): i for i, a in enumerate(self.residues)}
        self.odd = np.arange(1, q - 1, 2)
        # Uniform absolute tail of L and L' on Re(s)=4.5. Use whole periods modulo q.
        self.cut = q * math.ceil(300 * max(self.scales) ** (9 / 7) / q)
        ns = np.arange(1, self.cut + 1).reshape((-1, q))
        self.logs = np.log(ns)
        self.coefficients = auxiliary(ns, discriminant) * np.exp(
            -self.right * self.logs
        )

    def derivatives(self, amplitudes):
        sums = amplitudes.sum(axis=0)[self.residues - 1]
        slopes = (amplitudes * self.logs).sum(axis=0)[self.residues - 1]
        return -np.fft.ifft(slopes)[1::2] / np.fft.ifft(sums)[1::2]

    def check_arb(self, center):
        from flint import acb, acb_series, ctx, dirichlet_char

        ctx.prec = 110
        q, m = self.q, self.modulus
        conductor = m * q
        g = int(self.residues[1])
        # This lift evaluates the q-component alone.
        unit = g + q * (((1 - g) * pow(q, -1, m)) % m)
        s = acb(self.right, center)
        slopes = {}
        for index in range(1, conductor):
            if math.gcd(index, conductor) != 1:
                continue
            chi = dirichlet_char(conductor, index)
            if not chi.is_primitive() or chi.parity() != int(self.discriminant < 0):
                continue
            j = round(np.angle(complex(chi(unit))) * (q - 1) / (2 * np.pi)) % (q - 1)
            if j % 2 == 0:
                continue
            assert j not in slopes
            values = acb_series([s, 1], prec=2).dirichlet_l(chi)
            slopes[j] = complex(values[1] / values[0])
        expected = np.array([slopes[int(j)] for j in self.odd])
        computed = self.derivatives(
            self.coefficients * np.exp(-1j * center * self.logs)
        )
        assert np.max(np.abs(expected - computed)) < 2e-10

    def contour(self, center):
        radius = band_radius(self.q, center, self.order)
        # Nested trapezoidal grids on the two infinite vertical lines, truncated
        # where the halo bound in the note applies. Gamma terms cancel on centering.
        panels = 2 * math.ceil(radius / 0.16)
        step = 2 * radius / panels
        stride = np.exp(-1j * step * self.logs)
        base = self.coefficients * np.exp(-1j * center * self.logs)
        carrier = np.exp(1j * center * self.log_scales)
        right, left = self.right, 1 - self.right
        fine = np.zeros((len(self.scales), len(self.odd)), dtype=complex)
        coarse = np.zeros_like(fine)
        for k in range(panels + 1):
            v = -radius + k * step
            if k % 96 == 0:
                amplitudes = base * np.exp(-1j * v * self.logs)
            deriv = self.derivatives(amplitudes)
            deriv -= deriv.mean()
            right_factor = (
                carrier
                * np.exp((right + 1j * v) * self.log_scales)
                * transform(right + 1j * v, self.order)
            )
            left_factor = (
                carrier
                * np.exp((left + 1j * v) * self.log_scales)
                * transform(left + 1j * v, self.order)
            )
            integrand = np.outer(right_factor, deriv) + np.outer(
                left_factor, deriv.conj()
            )
            if k in (0, panels):
                integrand *= 0.5
            fine += integrand
            if k % 2 == 0:
                coarse += integrand
            amplitudes *= stride
        return fine * (step / (2 * np.pi)), coarse * (step / np.pi)

    def arithmetic(self, center, rates):
        traces = []
        for y in self.scales:
            ns = np.flatnonzero(rates[: math.floor(32 * y) + 1])
            direct = np.zeros(self.q, dtype=complex)
            dual = np.zeros_like(direct)
            for n in ns:
                weight = auxiliary(n, self.discriminant) * rates[n]
                if not weight or n % self.q == 0:
                    continue
                ratio = n / y
                direct[n % self.q] += (
                    weight
                    * window(ratio, self.order)
                    * np.exp(-1j * center * math.log(ratio))
                )
                ratio = 1 / (n * y)
                dual[n % self.q] += (
                    weight
                    / n
                    * window(ratio, self.order)
                    * np.exp(-1j * center * math.log(ratio))
                )
            first = (self.q - 1) * np.fft.ifft(direct[self.residues])[1::2]
            second = np.fft.fft(dual[self.residues])[1::2]
            raw = -first - second
            traces.append(raw - raw.mean())
        return np.array(traces)

    def differences(self, traces, a, y):
        row = traces[int(np.searchsorted(self.scales, y))]
        phase = np.exp(-2j * np.pi * self.odd * self.indices[a % self.q] / (self.q - 1))
        return -0.5 * (phase[:, None] - phase[None, :]) * (row[:, None] - row[None, :])

    def projection(self, traces, a, y):
        row = traces[int(np.searchsorted(self.scales, y))]
        phase = np.exp(-2j * np.pi * self.odd * self.indices[a % self.q] / (self.q - 1))
        return -(phase @ row) / len(row)


def distinct_product(left, right):
    """Average over ordered label pairs on disjoint conjugacy orbits."""
    count = left.shape[0]
    conjugate = np.arange(count - 1, -1, -1)
    right = right.conj()
    total = right.sum()
    rows, cols = right.sum(axis=1), right.sum(axis=0)
    result = 0j
    for i in range(count):
        for j in range(i + 1, count):
            blocked = list({i, j, int(conjugate[i]), int(conjugate[j])})
            allowed = (
                total
                - rows[blocked].sum()
                - cols[blocked].sum()
                + right[np.ix_(blocked, blocked)].sum()
            )
            result += 2 * left[i, j] * allowed
    return result / count**4


def statistics(family, left, right, weights):
    q, x = family.q, family.q + 2
    diff = family.differences

    def pair(a, y, b, z):
        left_pair, right_pair = diff(left, a, y), diff(right, b, z)
        value = distinct_product(left_pair, right_pair)
        if q == 17:
            count = len(left_pair)
            orbit = [frozenset((i, count - 1 - i)) for i in range(count)]
            direct = (
                sum(
                    left_pair[i, j] * right_pair[r, s].conjugate()
                    for i, j, r, s in np.ndindex(count, count, count, count)
                    if (orbit[i] | orbit[j]).isdisjoint(orbit[r] | orbit[s])
                )
                / count**4
            )
            assert abs(value - direct) < 1e-10
        return value

    main = pair(2, x, 2, x)
    correction = omega(x, family.discriminant)
    for d, weight in weights.items():
        main += correction * weight * pair(d, d, 2 * pow(d, -1, q), x / d)
    calibration = pair(5, 5, 5, 5) / math.log(5) ** 2
    return main / math.log(x) ** 2, calibration


def mixed_statistics(first, second, weights):
    left, left_traces = first
    right, right_traces = second
    x = left.q + 2

    def sign(n):
        return omega(n, left.discriminant) * omega(n, right.discriminant)

    def pair(a, y, b, z):
        return (
            left.projection(left_traces, a, y)
            * right.projection(right_traces, b, z).conjugate()
        )

    main = sign(x) * pair(2, x, 2, x)
    for d, weight in weights.items():
        main += (
            omega(x, right.discriminant)
            * sign(d)
            * weight
            * pair(d, d, 2 * pow(d, -1, left.q), x / d)
        )
    return main / math.log(x) ** 2, sign(5) * pair(5, 5, 5, 5) / math.log(5) ** 2


def check_factor_intervals():
    """Check the factor-restricted gate with exact integer endpoints."""
    intervals = (1000, 10000, 100000)
    limit = 2 * max(intervals) + 2
    primes, rates = prime_weights(limit)
    small_primes = [int(p) for p in np.flatnonzero(primes[: math.isqrt(limit) + 1])]
    # The illustrative interval is [3/11 + 1/132, 1/2 - 1/132].
    powers = [d**132 for d in range(math.isqrt(limit) + 1)]
    print("X       twins     balanced semiprimes     retained semiprimes")
    for start in intervals:
        twins = balanced = retained = 0
        for qi in np.flatnonzero(primes[start + 1 : 2 * start + 1]) + start + 1:
            q = int(qi)
            if q % 6 != 5:
                continue
            x = q + 2
            lower, upper = x**37, x**65
            is_twin = bool(primes[x])
            in_band = False
            if is_twin:
                twins += 1
            else:
                r = next(p for p in small_primes if x % p == 0)
                if primes[x // r] and r**4 >= x:
                    balanced += 1
                    in_band = lower <= powers[r] <= upper
                    retained += int(in_band)
            weights = {
                d: b
                for d, b in factor_weights(x).items()
                if lower <= powers[d] <= upper
            }
            g = exact_weight(x, rates, weights) / math.log(x) ** 2
            assert (g > 5 / 12) == (is_twin or in_band)
        print(f"{start:<7} {twins:5}     {balanced:19}     {retained:19}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", nargs="+", type=int, default=[17, 23, 29, 53, 173])
    parser.add_argument("--arb", action="store_true")
    args = parser.parse_args()
    primes, rates = prime_weights(32 * (max(args.q) + 2) + 2)
    interval_values = {}
    print(
        "q    x      g       Re mixed Z12     two-family four     four-family four     max |trace diff|",
        flush=True,
    )
    for q in args.q:
        if q < 17 or not primes[q] or q % 6 != 5:
            raise ValueError("expected a prime q >= 17 congruent to 5 mod 6")
        x = q + 2
        order = math.ceil(4 * math.log(math.log(q)) / math.log(4))
        weights = factor_weights(x)
        scales = {5, x, *weights, *(x / d for d in weights)}
        centers = centers_for(q, order)
        analytic, arithmetic, errors = {}, {}, []
        for disc, slots in ((1, (0, 1)), (-4, (2, 3)), (8, (1,)), (-8, (3,))):
            family = BandCalculation(q, disc, order, scales)
            for idx in slots:
                center = centers[idx]
                if args.arb:
                    family.check_arb(center)
                fine, coarse = family.contour(center)
                other = family.arithmetic(center, rates)
                errors.append(np.max(np.abs(fine - other)))
                assert np.max(np.abs(fine - coarse)) < 2e-5, (
                    q,
                    disc,
                    center,
                    "quadrature",
                )
                analytic[disc, idx] = family, fine
                arithmetic[disc, idx] = family, other

        def results(inputs):
            first, _ = statistics(
                inputs[1, 0][0], inputs[1, 0][1], inputs[1, 1][1], weights
            )
            second, cal = statistics(
                inputs[-4, 2][0], inputs[-4, 2][1], inputs[-4, 3][1], weights
            )
            mixed_first, _ = mixed_statistics(inputs[1, 0], inputs[8, 1], weights)
            mixed_second, mixed_cal = mixed_statistics(
                inputs[-4, 2], inputs[-8, 3], weights
            )
            return np.array([first, second, cal, mixed_first, mixed_second, mixed_cal])

        measured, independent = results(analytic), results(arithmetic)
        assert max(errors) < 3e-5, (q, errors)
        assert np.max(np.abs(measured - independent)) < 3e-6
        first, second, cal, mixed_first, mixed_second, mixed_cal = measured
        quartic = (first * (2 * second - cal)).real
        mixed_quartic = (mixed_first * (2 * mixed_second - mixed_cal)).real
        g = exact_weight(x, rates, weights) / math.log(x) ** 2
        if q in (17, 23, 29):
            interval_values[q] = mixed_quartic if mixed_first.real > 5 / 12 else 0.0
        print(
            f"{q:<4} {x:<4} {g:7.4f}     {mixed_first.real:12.7f}   {quartic:14.7f}     {mixed_quartic:14.7f}"
            f"       {max(errors):.2e}",
            flush=True,
        )
    if len(interval_values) == 3:
        estimate = sum(interval_values.values())
        exact = sum(bool(primes[q] and primes[q + 2]) for q in range(17, 33))
        assert round(estimate) == exact == 2
        print(
            f"(16,32]: gated four-family sum = {estimate:.8f}, integer count = {exact}"
        )
    check_factor_intervals()
    print("All four bands, all three conductors, and both signed products checked.")


if __name__ == "__main__":
    main()
