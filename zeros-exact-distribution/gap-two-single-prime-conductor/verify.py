"""Check the single-conductor zero statistic in the accompanying research note.

Run: python verify.py
Optional independent logarithmic-derivative check: python verify.py --arb
Dependencies: numpy, scipy. The optional check also uses python-flint.
"""

import argparse
import math
from functools import lru_cache

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.integrate import quad
from scipy.interpolate import BSpline
from scipy.special import bernoulli, k0


BAND = 32
J = 2 * k0(1.0)
SAMPLES = (127, 139, 157, 241)


def cutoff(q, window_order):
    return 8 * window_order + BAND * math.log(32 * math.log(q))


@lru_cache(maxsize=None)
def spectral_shifts(window_order):
    if window_order == 0:
        return np.array([0.0]), np.array([1.0])
    count = 2 * window_order
    # The sum of count uniforms on [0,1] has this normalized B-spline density.
    density = BSpline.basis_element(np.arange(count + 1), extrapolate=False)
    nodes, weights = leggauss(16)
    points = (np.arange(count)[:, None] + 0.5 * (nodes + 1)).ravel()
    masses = np.tile(0.5 * weights, count) * density(points)
    assert abs(masses.sum() - 1) < 3e-14
    return 8 * (points - window_order), masses


def transform(s, window_order):
    shifts, masses = spectral_shifts(window_order)
    base = np.exp(-np.cos((s + 1j * shifts) / BAND))
    return (2 * np.pi / (BAND * J)) * (masses @ base)


def window(u, window_order):
    log_u = math.log(u)
    frequency = BAND * abs(log_u)
    integral, _ = quad(lambda t: math.exp(-math.cosh(t)), 0, 5.5,
                       weight="cos", wvar=frequency, epsabs=3e-15, epsrel=3e-14)
    return (2 * integral / J) * np.sinc(4 * log_u / np.pi)**(2 * window_order)


def arithmetic_side(q, weights, window_order, twisted, periods=36):
    x = q - 2
    direct = -weights[2] * window(2 / x, window_order)
    for k in range(1, periods + 1):
        for sign in (-1, 1):
            n = k * q + 2 * sign
            if weights[n]:
                # The auxiliary factor is minus the quadratic character modulo 3.
                multiplier = (0, -1, 1)[n % 3] if twisted else 1
                direct -= sign * multiplier * weights[n] * window(n / x, window_order)
    # The proof bounds the entire reciprocal sum by 3*x**(-12).
    tail = 6 * (math.log(q) / (11 * periods**11) + 1 / (10 * periods**10))
    return direct, tail + 3 * x**(-12)


def prime_weights(limit):
    prime = np.ones(limit + 1, dtype=bool)
    prime[:2] = False
    for n in range(2, math.isqrt(limit) + 1):
        if prime[n]:
            prime[n*n::n] = False
    weights = np.zeros(limit + 1)
    for p in np.flatnonzero(prime):
        power = int(p)
        while power <= limit:
            weights[power] = math.log(p)
            power *= int(p)
    return prime, weights


def primitive_root(q):
    # For the small demonstration primes, enumerate the cyclic subgroup.
    for g in range(2, q):
        a, residues = 1, []
        while not residues or a != 1:
            residues.append(a)
            a = a * g % q
        if len(residues) == q - 1:
            return np.array(residues)
    raise ValueError("expected a prime modulus")


class CharacterFamily:
    """The odd family of conductor q or its primitive even twist of conductor 3q."""

    order = 10
    coefficients = [float(b) / math.factorial(2*k)
                    for k, b in enumerate(bernoulli(20)[::2]) if k]

    def __init__(self, q, window_order, twisted):
        residues = primitive_root(q)
        self.q = q
        self.window_order = window_order
        self.twisted = twisted
        self.conductor = (3 if twisted else 1) * q
        self.height = cutoff(self.conductor, window_order)
        self.log_conductor = math.log(self.conductor)
        self.log_x = math.log(q - 2)
        if twisted:
            # Input primes q are 1 mod 3. Lift each residue to 1 and 2 mod 3.
            plus = residues + q * ((1 - residues) % 3)
            minus = residues + q * ((2 - residues) % 3)
            self.u = np.concatenate((plus, minus)) / self.conductor
        else:
            self.u = residues / q
        odd = np.arange(1, q - 1, 2)
        two = int(np.flatnonzero(residues == 2)[0])
        self.phase = (2 / (q - 1)) * np.exp(-2j * np.pi * odd * two / (q - 1))
        if twisted:
            self.phase *= -1
        self.terms = 36 + math.ceil(self.height / 2)
        self.logs = np.log(np.arange(self.terms)[:, None] + self.u)
        self.base = self.terms + self.u
        self.base_log = np.log(self.base)

    def logarithmic_derivatives(self, s):
        finite = np.exp(-s * self.logs)
        zeta = np.sum(finite, axis=0)
        slope = -np.sum(finite * self.logs, axis=0)
        tail = np.exp(-s * self.base_log)
        integral = self.base * tail / (s - 1)
        zeta += integral + tail / 2
        slope -= integral * (self.base_log + 1 / (s - 1)) + tail * self.base_log / 2

        product, derivative = 1, 0
        for j in range(2 * self.order - 1):
            derivative = derivative * (s + j) + product
            product *= s + j
            if j % 2 == 0:
                factor = self.coefficients[j // 2] * tail / self.base**(j + 1)
                zeta += factor * product
                slope += factor * (derivative - product * self.base_log)
        if self.twisted:
            split = self.q - 1
            zeta = zeta[:split] - zeta[split:]
            slope = slope[:split] - slope[split:]
        values = np.fft.ifft(zeta)[1::2]
        derivatives = np.fft.ifft(slope)[1::2]
        return derivatives / values - self.log_conductor

    def factor(self, s):
        return np.exp(self.log_x * s) * transform(s, self.window_order)

    def vertical(self, t):
        s = 1.5 + 1j * t
        derivatives = self.logarithmic_derivatives(s)
        right = self.factor(s) * (self.phase @ derivatives)
        # Reflection leaves the gamma and conductor term times sum(phase) = 0.
        left = -self.factor(-0.5 + 1j*t) * (self.phase @ derivatives.conj())
        return (right - left).real

    def contour(self, height, order):
        vertical = integrate_panels(self.vertical, 0, height, order, width=2)

        def horizontal(sigma):
            s = sigma + 1j * height
            return (self.factor(s) * (self.phase @ self.logarithmic_derivatives(s))).imag

        top = integrate_panels(horizontal, -0.5, 1.5, order, width=0.125)
        return (vertical - top) / np.pi

    def check_arb(self):
        from flint import acb, acb_series, ctx, dirichlet_char

        ctx.prec = 120
        for s in (1.5 + 4.25j, -0.5 + 4.25j, 0.6 + 1j * self.height):
            weighted = 0j
            count = 0
            for i in range(1, self.conductor):
                if math.gcd(i, self.conductor) != 1:
                    continue
                chi = dirichlet_char(self.conductor, i)
                if chi.is_primitive() and chi.parity() == int(not self.twisted):
                    series = acb_series([acb(s.real, s.imag), 1], prec=2).dirichlet_l(chi)
                    weighted += complex(chi(2)).conjugate() * complex(series[1] / series[0])
                    count += 1
            assert count == (self.q - 1) // 2
            weighted *= 2 / (self.q - 1)
            computed = self.phase @ self.logarithmic_derivatives(s)
            assert abs(weighted - computed) < 2e-8 * max(1, abs(weighted)), (s, weighted, computed)


@lru_cache(maxsize=None)
def quadrature(order):
    return leggauss(order)


def integrate_panels(func, start, stop, order, width):
    nodes, weights = quadrature(order)
    count = math.ceil((stop - start) / width)
    step = (stop - start) / count
    values = []
    for panel in range(count):
        midpoint = start + (panel + 0.5) * step
        values.extend(w * func(midpoint + 0.5 * step * t) for t, w in zip(nodes, weights))
    return math.fsum(values) * step / 2


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arb", action="store_true", help="also compare L'/L to Arb")
    args = parser.parse_args()
    primes, weights = prime_weights(36 * max(SAMPLES) + 2)
    print("q     cond   N    T          arithmetic      contour(12)    contour(16)    trace/log(x)")
    for q in SAMPLES:
        x = q - 2
        assert primes[q] and q % 6 == 1
        local_order = math.ceil(3 * math.log(math.log(q)) / math.log(4))
        for window_order, twisted in ((0, False), (local_order, False), (local_order, True)):
            family = CharacterFamily(q, window_order, twisted)
            if args.arb:
                family.check_arb()
            direct, arithmetic_tail = arithmetic_side(q, weights, window_order, twisted)
            assert arithmetic_tail < 1e-14
            first, second = (family.contour(family.height, order) for order in (12, 16))
            print(f"{q:<5} {family.conductor:<6} {window_order:<3} {family.height:8.4f}  {direct:14.9f} "
                  f"{first:14.9f} {second:14.9f} {second/math.log(x):13.8f}", flush=True)
            assert abs(direct - second) < 2e-8 and abs(first - second) < 2e-8
            assert bool(primes[x]) == (second > 2 * math.log(x) / 3)
            nearest = weights[x]
            if not twisted:
                nearest -= weights[q+2] * window((q+2) / x, window_order)
            error_bound = 4.0**(-window_order) * (0.002 * math.log(q) + 0.004)
            assert abs(direct - nearest) < error_bound + arithmetic_tail + 2e-12
            if twisted:
                log_x = math.log(x)
                candidates = [0] + [log_x / j for j in range(1, x.bit_length(), 2) if 5**j <= x]
                recovered = min(candidates, key=lambda value: abs(value-second))
                assert abs(recovered - weights[x]) < 1e-12
    print("All arithmetic, contour, and threshold checks passed.")


if __name__ == "__main__":
    main()
