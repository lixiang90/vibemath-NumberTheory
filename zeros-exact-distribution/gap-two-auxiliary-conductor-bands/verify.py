"""Check the fixed-conductor formula on actual Dirichlet zeros.

Run with Python, numpy, scipy and python-flint. The small-conductor window
has a known taper, which is divided out before the prime indicator is formed.
No prime or prime-power data enter the zero search or the reconstruction.
"""

import math

import numpy as np
from flint import acb, ctx, dirichlet_char
from scipy.optimize import brentq
from scipy.special import loggamma
from scipy.stats import irwinhall


ORDER = 2
BOX = 0.16
PLATEAU_CENTER = 1.55
RADIUS = 48.0
HEIGHT = 51.37


def time_weight(t):
    return (2 * PLATEAU_CENTER * np.sinc(PLATEAU_CENTER * t / np.pi)
            * np.sinc(BOX * t / np.pi) ** (2 * ORDER))


def log_weight(u):
    law = irwinhall(2 * ORDER)
    right = (u + PLATEAU_CENTER) / (2 * BOX) + ORDER
    left = (u - PLATEAU_CENTER) / (2 * BOX) + ORDER
    return law.cdf(right) - law.cdf(left)


def moebius(n):
    value = 1
    p = 2
    while p * p <= n:
        if n % p == 0:
            n //= p
            if n % p == 0:
                return 0
            value = -value
        p += 1
    return -value if n > 1 else value


def roots_by_integer(limit):
    roots = {n: [] for n in range(2, limit + 1)}
    for r in range(2, limit + 1):
        power, exponent = r, 1
        while power <= limit:
            roots[power].append((r, exponent))
            power *= r
            exponent += 1
    return roots


class OddFamily:
    def __init__(self, modulus):
        self.q = modulus
        characters = [dirichlet_char(modulus, a) for a in range(1, modulus)]
        self.chars = [c for c in characters if c.parity()]
        assert all(c.is_primitive() for c in self.chars)
        residues = np.arange(1, modulus)
        self.arguments = [acb(int(a)) / modulus for a in residues]
        self.table = np.array([[complex(c(int(a))) for a in residues]
                               for c in self.chars])
        gauss = self.table @ np.exp(2j * np.pi * residues / modulus)
        self.phase = np.exp(-0.5j * np.angle(gauss / (1j * math.sqrt(modulus))))

    def l_values(self, s):
        # One Hurwitz evaluation per residue serves every character.
        values = np.array([complex(acb(s).zeta(a)) for a in self.arguments])
        return self.q ** (-s) * np.einsum("ij,j->i", self.table, values)

    def hardy_values(self, t):
        theta = (0.5 * t * math.log(self.q / np.pi)
                 + loggamma(0.75 + 0.5j * t).imag)
        z = self.phase * np.exp(1j * theta) * self.l_values(0.5 + 1j * t)
        assert np.max(np.abs(z.imag)) < 1e-9
        return z.real

    def winding_counts(self, steps):
        # Continue arg L to the critical line along the two horizontal edges.
        # On Re(s)=2 its continuous argument is the principal argument.
        def edge(t):
            for panels in (steps, 2 * steps, 4 * steps, 8 * steps):
                xs = 0.5 + 1.5 * np.linspace(1, 0, panels + 1) ** 2
                values = np.array([self.l_values(x + 1j * t) for x in xs])
                phases = np.unwrap(np.angle(values), axis=0)
                if np.max(np.abs(np.diff(phases, axis=0))) < 1.25:
                    return phases[-1]
            raise ArithmeticError("horizontal argument sampling did not stabilize")

        gamma_phase = (HEIGHT * math.log(self.q / np.pi)
                       + 2 * loggamma(0.75 + 0.5j * HEIGHT).imag)
        counts = (gamma_phase + edge(HEIGHT) - edge(-HEIGHT)) / np.pi
        rounded = np.rint(counts)
        assert np.max(np.abs(counts - rounded)) < 1e-8
        return rounded.astype(int)

    def locate_zeros(self):
        totals = self.winding_counts(24)
        assert np.array_equal(totals, self.winding_counts(48))
        mesh = np.linspace(-HEIGHT, HEIGHT, math.ceil(2 * HEIGHT / 0.12) + 1)
        grid = np.array([self.hardy_values(t) for t in mesh])
        result = []
        for i, chi in enumerate(self.chars):
            def hardy(t):
                return complex(chi.hardy_z(t)).real

            xs, ys = mesh, grid[:, i]
            for _ in range(5):
                crossings = np.flatnonzero(np.signbit(ys[:-1]) != np.signbit(ys[1:]))
                if len(crossings) == totals[i]:
                    break
                # Refine only characters for which the independent count differs.
                mid = (xs[:-1] + xs[1:]) / 2
                xs = np.sort(np.r_[xs, mid])
                ys = np.array([hardy(t) for t in xs])
            assert len(crossings) == totals[i], (self.q, chi.number(), totals[i])
            roots = [brentq(hardy, xs[j], xs[j + 1], xtol=2e-12) for j in crossings]
            result.append(np.array(roots))
        print(f"v={self.q}, odd characters={len(result)}, symmetric zero count={sum(totals)}")
        return result

    def traces(self, zeros, h):
        assert abs(h) + RADIUS < HEIGHT
        return np.array([time_weight(t[np.abs(t - h) <= RADIUS] - h).sum()
                         for t in zeros])

    def prime_indicators(self, traces, h, limit):
        # All factors are defined from integers and from the recovered Fourier modes.
        assert limit * math.floor(math.exp(PLATEAU_CENTER + 2 * ORDER * BOX)) < self.q - 1
        rates = {}
        for r in range(2, limit + 1):
            projection = -np.mean(self.table[:, r - 1].conj() * traces)
            rates[r] = (math.sqrt(r) * np.exp(1j * h * math.log(r))
                        * projection / log_weight(math.log(r)))
        return {n: sum(moebius(j) * rates[r] for r, j in roots) / math.log(n)
                for n, roots in roots_by_integer(limit).items()}


def arithmetic_traces(family, h):
    # Independent arithmetic check. A character mean removes the gamma term.
    powers = {2: 2, 3: 3, 4: 2, 5: 5, 7: 7, 8: 2}
    assert math.floor(math.exp(PLATEAU_CENTER + 2 * ORDER * BOX)) == 8
    return -2 * sum((family.table[:, n - 1] * np.exp(-1j * h * math.log(n))).real
                    * math.log(p) * log_weight(math.log(n)) / math.sqrt(n)
                    for n, p in powers.items())


def check_arcs(u, max_q, width):
    # Evaluate arcs by two independent methods: exact integration of all
    # Fourier modes, and direct quadrature of the prime exponential sum.
    assert 2 * max_q**2 * width < u
    prime = np.ones(2 * u + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, math.isqrt(2 * u) + 1):
        if prime[p]:
            prime[p * p::p] = False
    primes = np.flatnonzero(prime & (np.arange(2 * u + 1) > u))
    logs = np.log(primes)
    indicator = np.zeros(u)
    indicator[primes - u - 1] = logs
    spectrum = np.fft.rfft(indicator, n=2 * u)
    correlations = np.fft.irfft(np.abs(spectrum) ** 2, n=2 * u)
    lags = np.arange(1 - u, u)

    def integrated(shift):
        indices = lags + shift
        ramanujan = np.zeros(len(indices), dtype=int)
        for d in range(1, max_q + 1):
            coefficient = d * sum(moebius(j) for j in range(1, max_q // d + 1))
            ramanujan[indices % d == 0] += coefficient
        kernel = (2 * width / u) * np.sinc(2 * width * indices / u)
        return np.dot(correlations[np.abs(lags)], kernel * ramanujan)

    direct = np.zeros(2)
    fractions = {}
    for q in range(1, max_q + 1):
        numerators = np.array([a for a in range(q) if math.gcd(a, q) == 1], dtype=int)
        fractions[q] = (numerators, np.exp(4j * np.pi * numerators / q), primes % q)
    nodes, weights = np.polynomial.legendre.leggauss(36)
    for node, weight in zip(nodes, weights):
        beta = node * width / u
        masses = logs * np.exp(2j * np.pi * beta * primes)
        for q, (numerators, phases, residues) in fractions.items():
            grouped = (np.bincount(residues, weights=masses.real, minlength=q)
                       + 1j * np.bincount(residues, weights=masses.imag, minlength=q))
            energy = np.abs((q * np.fft.ifft(grouped))[numerators]) ** 2
            direct[0] += weight * (width / u) * energy.sum()
            direct[1] += (weight * (width / u)
                          * (energy @ phases * np.exp(4j * np.pi * beta)).real)
    mass, shifted = integrated(0), integrated(2)
    assert max(abs(direct - np.array([mass, shifted]))) < 1e-8 * u
    arithmetic = indicator[:-2] @ indicator[2:]
    assert abs(correlations[2] - arithmetic) < 1e-8 * u
    print(f"  U={u}, K={max_q}, W={width}, major energy share={mass / correlations[0]:.6f}, "
          f"major/U={shifted/u:.6f}, minor/U={(arithmetic-shifted)/u:.6f}, total/U={arithmetic/u:.6f}")


def main():
    ctx.prec = 75
    predictions = []
    for q, h, limit in ((43, 1.25, 5), (71, -1.75, 7)):
        family = OddFamily(q)
        zeros = family.locate_zeros()
        traces = family.traces(zeros, h)
        expected = arithmetic_traces(family, h)
        error = np.max(np.abs(traces - traces.mean() - expected))
        assert error < 0.001
        predicted = family.prime_indicators(traces, h, limit)
        prime_error = max(abs(value - (n in {2, 3, 5, 7}))
                          for n, value in predicted.items())
        assert prime_error < 0.01
        print(f"  centered trace error={error:.3g}, prime indicator error={prime_error:.3g}")
        predictions.append(predicted)
        if q == 71:
            same_family = sum((predicted[n] * predicted[n + 2].conjugate()).real
                              for n in range(2, 6))
            assert round(same_family) == 2
            print(f"  one conductor, one height, gap-two total={same_family:.8f}")
            # All products m*n of supported prime powers are below q-1.
            assert 8**2 < q - 1
            h2 = 0.85
            other = family.traces(zeros, h2)
            observed = np.mean((traces - traces.mean()) * (other - other.mean()))
            exact = np.mean(expected * arithmetic_traces(family, h2))
            diagonal = 2 * sum((math.log(p) * log_weight(math.log(n))) ** 2 / n
                               * math.cos((h - h2) * math.log(n))
                               for n, p in ((2, 2), (3, 3), (4, 2), (5, 5), (7, 7), (8, 2)))
            assert abs(exact - diagonal) < 1e-12 and abs(observed - diagonal) < 0.001
            print(f"  raw two-height covariance={observed:.9f}, diagonal={diagonal:.9f}")
    left, right = predictions
    values = {n: (left[n] * right[n + 2].conjugate()).real for n in range(2, 6)}
    for n, value in values.items():
        assert round(value) == (n in (3, 5))
    assert round(sum(values.values())) == 2
    print("  gap-two rows:", " ".join(f"{n}:{v:.8f}" for n, v in values.items()))
    print(f"  total={sum(values.values()):.8f}, nearest integer=2")
    for u, max_q in ((8192, 16), (32768, 32), (131072, 64)):
        check_arcs(u, max_q, 4)


if __name__ == "__main__":
    main()
