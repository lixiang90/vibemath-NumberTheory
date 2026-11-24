"""Reproduce the distinct-character zero-pair computation in the research note.

Run with numpy and scipy: python -B gap-two-chen-zero-pairs/verify.py
The adjacent gap-two-single-prime-conductor/verify.py supplies the Mellin
window and the Euler--Maclaurin/FFT evaluation of Dirichlet L'/L.
"""

import argparse
import math
import runpy
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.special import digamma


SHARED = runpy.run_path(str(Path(__file__).resolve().parents[1]
                            / "gap-two-single-prime-conductor" / "verify.py"))
CharacterFamily = SHARED["CharacterFamily"]
window = SHARED["window"]
transform = SHARED["transform"]
prime_weights = SHARED["prime_weights"]
primitive_root = SHARED["primitive_root"]


def factor_weights(x):
    result = {}
    for d in range(5, math.isqrt(x) + 1):
        if math.gcd(d, 6) != 1 or d**4 < x:
            continue
        orientations = 1 if d*d == x else 2
        result[d] = orientations * math.log(x)**2 / (4 * math.log(d) * math.log(x/d))
    return result


def exact_weight(x, rates, weights):
    return rates[x]**2 + sum(value * rates[d] * rates[x//d]
                             for d, value in weights.items() if x % d == 0)


def integrate(func, start, stop, order, panel_width):
    nodes, weights = leggauss(order)
    panels = math.ceil((stop-start) / panel_width)
    step = (stop-start) / panels
    total = 0
    for j in range(panels):
        mid = start + (j + 0.5) * step
        for node, weight in zip(nodes, weights):
            total = total + weight * func(mid + 0.5 * step * node)
    return total * step / 2


class TraceCalculation:
    def __init__(self, q):
        self.q = q
        self.x = q + 2
        self.order = math.ceil(4 * math.log(math.log(q)) / math.log(4))
        self.family = CharacterFamily(q, self.order, twisted=True)
        self.weights = factor_weights(self.x)
        self.scales = np.array(sorted({self.x, *self.weights, *(self.x/d for d in self.weights)}))
        self.indices = {y: i for i, y in enumerate(self.scales)}
        self.log_scales = np.log(self.scales)
        self.residues = primitive_root(q)
        self.indices_mod_q = {int(n): i for i, n in enumerate(self.residues)}
        self.odd = np.arange(1, q-1, 2)

    def factor(self, s):
        return np.exp(s * self.log_scales) * transform(s, self.order)

    def contour(self, gauss_order):
        height = self.family.height

        def vertical(t):
            right, left = 1.5 + 1j*t, -0.5 + 1j*t
            right_slope = self.family.logarithmic_derivatives(right)
            archimedean = (math.log(self.family.conductor / np.pi)
                           + 0.5 * (digamma(left/2) + digamma((1-left)/2)))
            left_slope = -right_slope.conj() - archimedean
            return np.outer(self.factor(right), right_slope) - np.outer(self.factor(left), left_slope)

        def upper(sigma):
            s = sigma + 1j*height
            return np.outer(self.factor(s), self.family.logarithmic_derivatives(s))

        sides = integrate(vertical, 0, height, gauss_order, 1.5)
        top = integrate(upper, -0.5, 1.5, gauss_order, 0.125)
        # Reflection in the real axis reverses the order of the character list.
        all_residues = ((sides + sides[:, ::-1].conj()) / (2*np.pi)
                        + (top[:, ::-1].conj() - top) / (2j*np.pi))
        # Even primitive L-functions have a simple trivial zero at zero.
        return all_residues - transform(0, self.order)

    def arithmetic_traces(self, rates):
        # On Re(s)=1/2 the gamma integral includes the correction at s=0.
        def gamma_integrand(t):
            s = 0.5 + 1j*t
            kernel = math.log(self.family.conductor/np.pi) + digamma(s/2).real
            return kernel * self.factor(s).real / np.pi

        common = integrate(gamma_integrand, 0, self.family.height, 16, 1.5)
        result = []
        for i, y in enumerate(self.scales):
            ns = np.flatnonzero(rates[:math.floor(32*y)+1])
            direct = np.zeros(self.q)
            dual = np.zeros(self.q)
            for n in ns:
                eta = (0, 1, -1)[n % 3]
                if eta == 0 or n % self.q == 0:
                    continue
                direct[n % self.q] += eta * rates[n] * window(n/y, self.order)
                dual[n % self.q] += (eta * rates[n] / n) * window(1/(n*y), self.order)
            right = (self.q-1) * np.fft.ifft(direct[self.residues])[1::2]
            reflected = (self.q-1) * np.fft.ifft(dual[self.residues])[1::2].conj()
            result.append(common[i] - right - reflected)
        return np.array(result)

    def terms(self, traces, a, y):
        angle = -2j*np.pi*self.odd*self.indices_mod_q[a % self.q] / (self.q-1)
        return (-2 / (self.q-1)) * np.exp(angle) * traces[self.indices[y]]

    def pair_sums(self, traces):
        p = self.terms(traces, 2, self.x)
        matrix = np.outer(p, p)
        for d, weight in self.weights.items():
            u = self.terms(traces, d, d)
            v = self.terms(traces, 2 * pow(d, -1, self.q), self.x/d)
            matrix += weight * np.outer(u, v)
        same = self.odd[:, None] == self.odd[None, :]
        conjugate = self.odd[:, None] + self.odd[None, :] == self.q-1
        mask = same | conjugate
        values = np.array([matrix.sum(), matrix[mask].sum(), matrix[~mask].sum()])
        assert np.max(np.abs(values.imag)) < 2e-8, values
        return values.real


def check_arithmetic_classes(limit):
    primes, rates = prime_weights(limit+2)
    counts = dict(prime=0, balanced=0, other=0)
    for q_raw in np.flatnonzero(primes):
        q = int(q_raw)
        if q < 17 or q > limit or q % 6 != 5:
            continue
        x = q+2
        weights = factor_weights(x)
        value = exact_weight(x, rates, weights) / math.log(x)**2
        if primes[x]:
            kind, expected = "prime", 1.0
        elif any(x % d == 0 and primes[d] and primes[x//d] for d in weights):
            kind, expected = "balanced", 0.5
        else:
            assert -1e-12 < value <= 13/36 + 1e-12, (q, value)
            counts["other"] += 1
            continue
        counts[kind] += 1
        assert abs(value - expected) < 1e-12, (q, kind, value)
    print(f"Arithmetic classification through {limit}: {counts}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--q", nargs="+", type=int, default=[29, 47, 53, 173, 383])
    args = parser.parse_args()
    primes, rates = prime_weights(32*(max(args.q)+2)+2)
    check_arithmetic_classes(100_000)
    print("q    x     T       G/L^2       full/L^2    removed/L^2   cross/L^2   max trace error")
    for q in args.q:
        if q < 17 or not primes[q] or q % 6 != 5:
            raise ValueError("q must be a prime at least 17 congruent to 5 modulo 6")
        calc = TraceCalculation(q)
        arithmetic = calc.arithmetic_traces(rates)
        coarse = calc.contour(12)
        fine = calc.contour(16)
        err = np.max(np.abs(fine - arithmetic))
        assert err < 8e-7 and np.max(np.abs(fine-coarse)) < 8e-7, (q, err)
        full, removed, cross = calc.pair_sums(fine)
        direct = calc.pair_sums(arithmetic)
        assert np.max(np.abs(direct - [full, removed, cross])) < 5e-7
        target = exact_weight(calc.x, rates, calc.weights)
        assert abs(target-full) < 5e-7, (q, target, full)
        scale = math.log(calc.x)**2
        print(f"{q:<4} {calc.x:<5} {calc.family.height:6.2f}  {target/scale:11.8f}"
              f" {full/scale:12.8f} {removed/scale:13.8f} {cross/scale:12.8f} {err:13.2e}",
              flush=True)
    print("Individual traces, full sums, and distinct-character sums agree.")


if __name__ == "__main__":
    main()
