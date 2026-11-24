"""Finite checks for the Chen factor-ratio bounds and conditioned zero traces.

Run from the checkout with Python, NumPy, SciPy and python-flint.
The physical-row experiment computes prime-side masses. The separate
quadratic-conductor experiment locates zeros with the theta integral.
"""

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import log, pi
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import digamma
from scipy.stats import irwinhall


ROOT = Path(__file__).resolve().parents[1]
A = Fraction(3, 11)
D = 11


def load_article(slug):
    spec = spec_from_file_location(slug.replace("-", "_"), ROOT / slug / "verify.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def weight(t):
    return 2 * t * (1 - t)


def lower_mass(mass, cap):
    end = cap * (Fraction(1, 2) - A)
    u = min(end, mass)
    return (weight(A) * u + (1-2*A) * u**2 / cap
            - 2 * u**3 / (3 * cap**2) + (mass-u) / 2)


def upper_mass(mass, cap):
    return mass / 2 - 2 * mass**3 / (3 * cap**2)


def primitives(t):
    return t**2 - Fraction(2, 3) * t**3


def exact_checks():
    partial = Fraction(3, 4) * Fraction(15, 16) * Fraction(35, 36)
    assert partial == Fraction(175, 256) and 16 * partial < D
    assert weight(A) == Fraction(48, 121)
    # Integrate bounded piecewise-constant densities exactly, including both
    # extremizers (packing the full density at either end of the interval).
    half = Fraction(1, 2)
    step = (half-A) / 6
    integrals = [primitives(A + (k+1)*step) - primitives(A + k*step) for k in range(6)]
    hits = [0, 0]
    for levels in product((0, Fraction(D, 2), D), repeat=6):
        mass = sum(levels) * step
        statistic = sum(density * area for density, area in zip(levels, integrals))
        low, high = lower_mass(mass, D), upper_mass(mass, D)
        assert low <= statistic <= high
        hits[0] += low == statistic
        hits[1] += high == statistic
    assert min(hits) >= 7
    print(f"exact density tests: {3**6} profiles; lower/upper equality {hits}")


def physical_windows():
    previous = load_article("gap-two-coherent-zero-heights")
    sizes = (1024, 8192, 32768)
    limit = int(max(sizes)**1.5 / 2) + 2
    factors = previous.smallest_factors(limit)
    witness = previous.test_witness(limit)
    beta = float(weight(A))
    print(f"finite four-sign witness: {witness}")
    for y in sizes:
        twins, semis = previous.admitted_rows(y, factors, witness)
        length = digamma(int(y**1.5/2) + 1) - digamma(y)

        def row_weights(ns):
            return np.log(ns) * np.log(ns + 2) / (ns * length)

        n, p, s, t = semis.T
        assert np.all((s < t) & (s*t + p == 2*n+2))
        ratios = np.log(s) / np.log(s*t)
        weights = row_weights(n)
        twin_mass = float(row_weights(twins).sum())
        semi_mass = float(weights.sum())
        measured = float(weights @ weight(ratios))
        total = twin_mass + semi_mass
        new = max(0, total - brentq(lambda x: float(lower_mass(x, D)) - measured, 0, 10))
        old = max(0, total - measured / beta)
        # These finite comparisons check the arithmetic reference only.
        assert 0 < old < new < twin_mass
        assert float(lower_mass(semi_mass, D)) < measured < float(upper_mass(semi_mass, D))
        print(f"Y={y}: twins={len(twins)}, semis={len(n)}, J={twin_mass:.7f}, "
              f"b={semi_mass:.7f}, M0={measured:.7f}, old={old:.7f}, new={new:.7f}")


def zero_traces():
    previous = load_article("gap-two-quadratic-conductor-covariance")
    band, box = 1.55, 0.30
    support = band + 4 * box
    tags = (5, 7, 11, 13)
    assert log(13) < support < log(17) and band > 4 * box
    law = irwinhall(4)

    def kernel(u):
        return law.cdf((u + band) / (2*box) + 2) - law.cdf((u - band) / (2*box) + 2)

    # Higher odd powers of each free prime are outside the arithmetic support.
    target = -4 * sum(log(p) / np.sqrt(p) * kernel(log(p)) for p in tags)
    ds = previous.primes_up_to(24_000)
    ds = ds[(ds > 800) & (ds % 48 == 23)]
    signs = np.stack([previous.character_column(ds, p) for p in tags], axis=1)
    groups, count = [], 0
    for sign in (1, -1):
        chosen = ds[np.all(signs == sign, axis=1)]
        assert len(chosen) >= 4
        stats = []
        for d in chosen:
            roots = previous.QuadraticTheta(int(d)).zeros()
            count += len(roots)
            trace = np.sum(2*band*np.sinc(band*roots/pi)*np.sinc(box*roots/pi)**4)
            stats.append(trace - log(d))
        groups.append(np.array(stats))
    measured = groups[0].mean() - groups[1].mean()
    spread = max(np.ptp(group) for group in groups)
    assert abs(measured-target) < 0.003 and spread < 0.006
    print(f"actual zero groups: {len(groups[0])}/{len(groups[1])}, {count} zeros, "
          f"height {previous.TOP:.5f}, max within-group range {spread:.8f}")
    print(f"mean trace difference: zeros {measured:.8f}, prime side {target:.8f}")


if __name__ == "__main__":
    exact_checks()
    physical_windows()
    zero_traces()
