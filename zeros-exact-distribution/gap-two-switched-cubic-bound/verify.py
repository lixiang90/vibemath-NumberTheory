"""Numerical checks for the switched sieve and angular cubic zero statistic.

Run with Python, NumPy, SciPy and python-flint in the repository checkout.
The larger integer samples use the arithmetic side. The small conductor
sample computes both the zeros and the prime side independently.
"""

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd, lcm, log, log1p, pi, prod
from pathlib import Path

import numpy as np
from flint import arb, ctx
from scipy.special import digamma
from scipy.stats import irwinhall


A = 3 / 11
BINS = ((A, 0.35), (0.35, 0.43), (0.43, 0.5))
ROOT = Path(__file__).resolve().parents[1]


def article(slug):
    spec = spec_from_file_location(slug.replace("-", "_"), ROOT / slug / "verify.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def exact_sieve(primes):
    radius = 31
    small = [int(p) for p in primes if 2 < p <= radius]
    factors = {d: [p for p in small if d % p == 0] for d in range(1, radius+1, 2)}
    factors = {d: ps for d, ps in factors.items() if prod(ps) == d}
    weights = {d: Fraction(1, prod(p-2 for p in ps)) for d, ps in factors.items()}
    total = sum(weights.values())
    lambdas = {}
    for d, ps in factors.items():
        partial = sum(h for e, h in weights.items() if e*d <= radius and gcd(e, d) == 1)
        lambdas[d] = (-1)**len(ps) * prod(Fraction(p-1, p-2) for p in ps) * partial / total
    assert lambdas[1] == 1 and max(map(abs, lambdas.values())) == 1
    quadratic = sum(ld * le / prod(p-1 for p in small if lcm(d, e) % p == 0)
                    for d, ld in lambdas.items() for e, le in lambdas.items())
    assert quadratic == 1 / total
    prime_set = set(map(int, primes))
    products = [int(s*t) for s in primes if 37 <= s < 65 for t in primes if 149 <= t < 239]
    for shift in (-2, 2):
        count = sum((v+shift) in prime_set for v in products)
        bound = sum(sum(w for d, w in lambdas.items() if (v+shift) % d == 0)**2
                    for v in products)
        assert bound >= count > 0
        print(f"finite sieve shift {shift:+}: {count} prime endpoints, upper sum {float(bound):.6f}")
    print(f"Selberg quadratic form at R={radius}: {quadratic}")


def jackson_multipliers(kmax):
    order = kmax // 2 + 1
    triangle = 1 - np.abs(np.arange(1-order, order, dtype=float)) / order
    coeff = np.convolve(triangle, triangle)
    center = 2*order - 2
    result = coeff[center:] / coeff[center]
    assert len(result) == kmax + 1 and result[0] == 1
    return result


def arc_coefficients(x, interval, modes, multipliers):
    """Real cosine coefficients of the even four-arc mask after convolution."""
    lo, hi = interval
    positive = modes[1:]
    phase = 2*pi*x[:, None]*positive
    harmonics = (np.sin(phase*hi) - np.sin(phase*lo) + np.sin(phase*(1-lo))
                 - np.sin(phase*(1-hi))) / (pi*positive)
    return np.column_stack((2*x*(hi-lo), harmonics)) * multipliers


def finite_cubic(semis, y, kmax):
    """Prime-power permutation sum for the same fixed three zero centers per mode."""
    ll = log(y)
    length = digamma(int(y**1.5/2)+1) - digamma(y)
    xi = ll**-3
    modes = np.arange(kmax+1)
    heights = 2*pi*modes/(8*ll)
    heights[0] = 2*xi
    kernel = jackson_multipliers(kmax)
    output = np.zeros(len(BINS))
    for start in range(0, len(semis), 64):
        block = semis[start:start+64]
        ns = block[:, 0].astype(float)
        ps = block[:, 1:].astype(float)
        logs = np.log(ps)

        def coeff(h):
            z = np.exp(1j * logs[..., None]*np.asarray(h))
            r = z*z/ps[..., None]
            return (logs/np.sqrt(ps))[..., None] * (z*(1-r**5)/(1-r)).real

        left, mid = coeff([0])[:, :, 0], coeff([xi])[:, :, 0]
        moving = coeff(heights)
        cubic = np.zeros((len(ns), len(modes)))
        for j in range(3):
            i, l = (j+1) % 3, (j+2) % 3
            cubic += (left[:, i]*mid[:, l] + left[:, l]*mid[:, i])[:, None] * moving[:, j]
        cubic /= length * np.log(ns)[:, None]
        for i, interval in enumerate(BINS):
            output[i] += np.sum(arc_coefficients(np.log(ns)/(8*ll), interval, modes, kernel)*cubic)
    return output


def physical_check(previous, c2):
    sizes = (8192, 32768)
    limit = int(max(sizes)**1.5/2)+2
    factors = previous.smallest_factors(limit)
    witness = previous.test_witness(limit)
    print("finite witness:",witness)
    for y in sizes:
        twins, semis = previous.admitted_rows(y, factors, witness)
        n, _, s, t = semis.T
        theta = np.log(s)/np.log(s*t)
        length = digamma(int(y**1.5/2)+1)-digamma(y)
        w = np.log(n)*np.log(n+2)/(n*length)
        twin_mass = np.sum(np.log(twins)*np.log(twins+2)/(twins*length))
        masses = np.array([np.sum((w*2*theta*(1-theta))[(lo <= theta) & (theta < hi)])
                           for lo, hi in BINS])
        estimates = [finite_cubic(semis, y, k) for k in (128, 1024, 2048)]
        assert abs(np.sum(masses)-np.sum(w*2*theta*(1-theta))) < 1e-12
        assert max(abs(masses-estimates[-1])) < 0.002
        def bin_bound(interval, mass):
            lo, hi = interval
            last = min(hi, lo + max(0, mass)/(2*c2))
            return c2*(log(last/(1-last))-log(lo/(1-lo)))
        total = twin_mass + sum(w)
        lower_whole = total - bin_bound((A, 0.5), sum(masses))
        lower_bins = total - sum(bin_bound(interval, m) for interval, m in zip(BINS, masses))
        assert 0 < lower_whole < lower_bins < twin_mass
        print(f"Y={y}, triples={len(semis)}, arithmetic bin masses: {masses}")
        for k, estimate in zip((128, 1024, 2048), estimates):
            print(f"  shared modes K={k}: {estimate}")
        print(f"  whole/bin lower proxies and actual J: {lower_whole:.7f}, "
              f"{lower_bins:.7f}, {twin_mass:.7f}")


def actual_zero_check(previous):
    ds = previous.primes_up_to(12_000)
    ds = ds[(ds>800) & (ds % 48 == 23)]
    y, kmax = 17, 16
    modes = np.arange(kmax+1)
    xi = 0.09
    probe = 2*pi*modes/(8*log(y))
    probe[0] = 2*xi
    heights = np.r_[0, xi, probe]
    assert max(heights) + previous.RADIUS < previous.TOP
    zero_rows, nroots = [], 0
    for d in ds:
        roots = previous.QuadraticTheta(int(d)).zeros()
        nroots += len(roots)
        offsets = roots[:, None] - heights
        values = previous.zero_kernel(offsets) * (abs(offsets) <= previous.RADIUS)
        zero_rows.append(np.sum(values, axis=0)-log(d))
    zeros = -0.5*previous.centered(np.array(zero_rows))
    limit = int(np.exp(previous.B+4*previous.BOX))
    ps = previous.primes_up_to(limit)
    ps = ps[ps>3]
    logs = np.log(ps)
    coefs = np.zeros((len(ps),len(heights)))
    for odd in range(1,int(log(limit)/logs[0])+1,2):
        u = odd*logs
        g = (irwinhall.cdf((u+previous.B)/(2*previous.BOX)+2,4)
             - irwinhall.cdf((u-previous.B)/(2*previous.BOX)+2,4))
        coefs += (logs*np.exp(-u/2)*g)[:,None]*np.cos(u[:,None]*heights)
    signs = np.column_stack([previous.character_column(ds,int(p)) for p in ps])
    arithmetic = previous.centered(signs @ coefs)
    differences = []
    for n, tags in ((53,(5,11,53)), (83,(5,17,83)), (17,(17,19))):
        mark = previous.row_mark(ds,tags)
        weights = np.vstack([arc_coefficients(np.array([log(n)/(8*log(y))]), interval,
                            modes,jackson_multipliers(kmax))[0] for interval in BINS])
        def reads(data):
            moments = []
            for symbol in (1,-1):
                group = previous.centered(data[mark==symbol])
                moments.append(np.mean(group[:,0,None]*group[:,1,None]*group[:,2:],axis=0))
            return weights @ (moments[0]-moments[1])/(2*log(n))
        differences.extend(abs(reads(zeros)-reads(arithmetic)))
    assert max(differences)<0.001
    print(f"actual zero test: {len(ds)} conductors, {nroots} zeros, "
          f"9 same-family cubic comparisons; max difference {max(differences):.8f}")


def main():
    previous = article("gap-two-quadratic-conductor-covariance")
    primes = previous.primes_up_to(1_000_000)
    exact_sieve(primes)
    odd = primes[primes>2].astype(float)
    c2 = np.exp(np.sum(np.log1p(-1/(odd-1)**2)))
    d = 1/112
    base = c2*(1-2*A)
    improved = base - c2*(2*d-log1p(2*d))
    # Enclose the infinite product independently of the floating estimates.
    # The omitted sum of (p-1)^(-2) is bounded by 1/(limit-1).
    ctx.prec = 100
    product_hi = arb(1)
    for p in primes[1:]:
        product_hi *= 1-arb(1)/(int(p)-1)**2
    product_lo = product_hi*(1-arb(1)/(1_000_000-1))
    exact_factor = arb(5)/11-arb(1)/56+(arb(57)/56).log()
    assert product_lo*exact_factor > arb(29996)/100000
    assert product_hi*exact_factor < arb(29997)/100000
    print("Arb and the product tail place the improved ceiling in (0.29996, 0.29997)")
    print(f"truncated-product constants: C2={c2:.8f}, whole caps {base:.8f}/{improved:.8f}, "
          f"b cap {c2*log(8/3):.8f}")
    np.set_printoptions(precision=7)
    physical_check(article("gap-two-coherent-zero-heights"), c2)
    actual_zero_check(previous)


if __name__ == "__main__":
    main()
