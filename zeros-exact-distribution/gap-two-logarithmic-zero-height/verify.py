"""Arithmetic identities and a numerical Riemann-zero trace for the paper.

Requires Python 3.10+, NumPy and python-flint. Pass --contour to also integrate
the full Dirichlet families in the X=4 and X=5 windows. Finite rational checks use the
standard library. Fourier sums and quadrature accumulation use binary64.
"""

from collections import defaultdict
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import ceil, factorial, gcd, isqrt, lcm, log, pi, prod
import sys

import numpy as np
from flint import acb, acb_series, ctx, dirichlet_char


def factors(n):
    result = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            result[d] = result.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        result[n] = 1
    return result


def signed_divisors(primes):
    values = [(1, 1)]
    for p in primes:
        values += [(d * p, -sign) for d, sign in values]
    return values


def mobius(n):
    f = factors(n)
    return 0 if any(v > 1 for v in f.values()) else (-1) ** len(f)


def mangoldt(n):
    f = factors(n)
    return log(next(iter(f))) if len(f) == 1 else 0.0


def odd_part(n):
    while n % 2 == 0:
        n //= 2
    return n


def is_prime(n):
    return n > 1 and factors(n) == {n: 1}


def check_square_divisors():
    checked = 0
    for m in range(3, 20001, 2):
        f = factors(m)
        coefficients = defaultdict(int)
        for c, sign_c in signed_divisors(p for p, k in f.items() if k >= 2):
            for b, sign_b in signed_divisors(p for p in f if c % p):
                for v, sign_v in signed_divisors(factors(c)):
                    assert lcm(c * c, b * v) == c * c * b
                    coefficients[b * v] += sign_c * sign_b * sign_v
        sf = int(all(k == 1 for k in f.values()))
        assert all(coefficients[d] == sf * sign for d, sign in signed_divisors(f))
        checked += 1
    print("Odd integers in square-divisor check:", checked)


@lru_cache(None)
def model_ramp(x):
    # Integrate the positive part over 2+3 uniform boxes, by finite differences.
    widths = [Fraction(1, 2)] * 2 + [Fraction(1, 3)] * 3
    center = x + sum(widths)
    power = len(widths) + 1
    numerator = sum(
        (-1) ** j * max(center - 2 * sum(subset), 0) ** power
        for j in range(power) for subset in combinations(widths, j)
    )
    return numerator / (factorial(power) * prod(2 * w for w in widths))


def check_reflection():
    count = 0
    for r in (2, 3, 4, 5):
        for lengths in combinations(range(2, 12), r):
            mid = Fraction(sum(lengths), 2)
            value = sum((-1) ** j * model_ramp(mid - sum(subset))
                        for j in range(r + 1) for subset in combinations(lengths, j))
            if r == 2:
                assert 0 <= value <= min(lengths)
            elif r == 4:
                assert value <= 0
            else:
                assert value == 0
            count += 1
    print("Rational length configurations checked:", count)


G = 1 / 64
A, B = log(5 / 4), log(8 / 5)


def sinhc(z):
    return np.sinc(z / (1j * pi))


class Window:
    def __init__(self, x, modes=4096):
        self.x = x
        self.orders = (max(2, ceil(200 * log(log(x)))), max(2, ceil(20 * log(x))))
        self.norm = self.transform(1, normalize=False).real.item()
        self.frequencies = pi * np.arange(1, modes + 1)
        self.coefficients = self.transform(1j * self.frequencies)

    def transform(self, s, normalize=True):
        s = np.asarray(s, dtype=complex)
        interval = np.divide(np.exp(A*s) * np.expm1((B-A)*s), s,
                             out=np.full_like(s, B-A), where=(s != 0))
        result = interval * sinhc(G*s)**2
        for order in self.orders:
            result *= sinhc(G*s/order)**order
        return result / self.norm if normalize else result

    def weight(self, n):
        y = log(n / self.x)
        if y <= A - 4*G or y >= B + 4*G:
            return 0.0
        if A + 4*G <= y <= B - 4*G:
            return 1 / self.norm
        fourier = np.sum(self.coefficients * np.exp(-1j * self.frequencies * y))
        return (self.transform(0).real / 2 + fourier.real).item()


def check_arithmetic(x):
    window = Window(x)
    dmax = 2*x + 2
    weights = {n: window.weight(n) for n in range(x, 2*x + 1)}
    exceptional = {a**k for a in range(2, isqrt(dmax)+1)
                   for k in range(2, dmax.bit_length()) if a**k <= dmax}
    twins, powers, odd_sum = 0.0, 0.0, 0.0
    for n, w in weights.items():
        term = mangoldt(n) * mangoldt(odd_part(n + 2)) * w
        odd_sum += term
        if n in exceptional or n + 2 in exceptional:
            powers += term
        if is_prime(n) and is_prime(n + 2):
            twins += log(n) * log(n + 2) * w
    ap_sum = sum(-mobius(d)*log(d)*sum(mangoldt(n)*w for n, w in weights.items()
                                        if (n+2) % d == 0)
                 for d in range(3, dmax+1, 2))
    assert abs(ap_sum - odd_sum) < 1e-10 and abs(odd_sum - powers - twins) < 1e-10

    euler_by_power = [0.0, 0.0]
    euler_by_character = 0j
    for d in range(3, dmax+1, 2):
        mu = mobius(d)
        if not mu:
            continue
        primes = factors(d)
        phi = prod(p - 1 for p in primes)
        for p in primes:
            e = d // p
            n = p
            while n <= 2*x:
                if n >= x and weights[n]:
                    w = weights[n] * log(p)
                    if (n+2) % e == 0:
                        euler_by_power[int(n != p)] += -mu * log(d) * w / (p-1)
                    lift = n + e * ((1-n) * pow(e, -1, p) % p)
                    for j in range(1, d):
                        if gcd(j, d) == 1:
                            chi = dirichlet_char(d, j)
                            if chi.conductor() % p:
                                phase = complex(chi((-2) % d)).conjugate()
                                euler_by_character += -mu*log(d)/phi * phase*complex(chi(lift % d))*w
                n *= p
    assert abs(euler_by_character - sum(euler_by_power)) < 1e-10
    single = sum(log(p)**2 * w / (p-1) for p, w in weights.items() if is_prime(p))
    assert abs(single - euler_by_power[0]) < 1e-10
    print(f"X={x:3}  M={twins:.8f}  P={powers:.8f}  "
          f"I2+={euler_by_power[1]:.8f}  I1={single:.8f}")
    return twins, powers, sum(euler_by_power)


def check_zeta():
    window = Window(4)
    high_modes = Window(4, modes=8192)
    arithmetic = sum(mangoldt(n) * window.weight(n) for n in range(4, 9))
    assert max(abs(window.weight(n)-high_modes.weight(n)) for n in range(4, 9)) < 1e-12
    trivial = sum(float((4.0**(-2*k) * window.transform(-2*k)).real) for k in range(1, 20))
    # The first zero beyond 3000 closes the last truncation interval.
    roots = []
    start = 1
    while not roots or roots[-1].imag <= 3000:
        block = acb.zeta_zeros(start, 256)
        roots.extend(complex(z) for z in block)
        start += len(block)
    zeros = np.asarray(roots)
    y = zeros.imag
    terms = 2 * (np.exp(log(4)*zeros) * window.transform(zeros)).real
    for height in (64, 256, 1024, 3000):
        mask = y <= height
        spectral = 4 - np.sum(terms[mask]) - trivial
        print(f"zeta T={height:4}  count={np.count_nonzero(mask):4}  "
              f"trace={spectral:.11f}  arithmetic={arithmetic:.11f}  "
              f"diff={abs(spectral-arithmetic):.3g}")
    assert abs(spectral-arithmetic) < 1e-8


def check_full_contour(x):
    # For both verification windows the nonzero-weight moduli are prime.
    window = Window(x)
    truth, powers, local = check_arithmetic(x)
    single = sum(log(p)**2 * window.weight(p)/(p-1)
                 for p in range(x, 2*x+1) if is_prime(p))
    principal = 0.0
    chars = []
    for q in range(3, 2*x+3, 2):
        if not mobius(q):
            continue
        assert is_prime(q)
        weight = log(q) / (q-1)
        principal += weight
        for j in range(2, q):
            chi = dirichlet_char(q, j)
            assert chi.is_primitive()
            chars.append((weight*complex(chi((-2) % q)).conjugate(), chi))

    def integrand(s):
        series = acb_series([acb(s.real, s.imag), 1], prec=2)
        zeta = series.zeta()
        logarithmic = principal * complex(zeta[1]/zeta[0])
        for weight, chi in chars:
            value = series.dirichlet_l(chi)
            logarithmic += weight * complex(value[1]/value[0])
        return logarithmic * x**s * window.transform(s)

    def quadrature(fn, lo, hi, step, nodes, weights):
        panels = ceil((hi-lo)/step)
        width = (hi-lo)/panels
        return width/2 * sum(w*fn(lo + width*(i+(node+1)/2))
                             for i in range(panels) for node, w in zip(nodes, weights))

    arch = log(3)/2 * sum((x**(-float(k))*window.transform(-k)).real
                          for k in range(1, 30))
    for height in (32.25, 96.25, 256.25):
        values = []
        for order in (8, 12):
            nodes, weights = np.polynomial.legendre.leggauss(order)
            vertical = quadrature(lambda t: (integrand(1.5+1j*t)-
                                              integrand(-0.5+1j*t)).real,
                                  0, height, 2, nodes, weights)
            top = quadrature(lambda sigma: integrand(sigma+1j*height).imag,
                             -0.5, 1.5, 1/8, nodes, weights)
            integral = (vertical-top)/pi
            values.append(-integral-arch-powers-local)
        detector = ((log(3)/2 - principal)*window.transform(0).real
                    - integral - powers - (local-single))
        print(f"full family X={x} T={height}  reconstruction={values[-1]:.8f}  "
              f"arithmetic={truth:.8f}  detector={detector:.8f}  "
              f"quadrature change={abs(values[0]-values[1]):.2g}")
        assert abs(values[0]-values[1]) < 2e-4
    assert abs(values[-1]-truth) < 0.05


if __name__ == "__main__":
    ctx.prec = 100
    check_square_divisors()
    check_reflection()
    for x in (4, 5, 20, 40, 100):
        check_arithmetic(x)
    check_zeta()
    if "--contour" in sys.argv:
        for x in (4, 5):
            check_full_contour(x)
