"""Interval certificate for Artin representation 2.108.3t2.b.a.

Run with Python 3.12, python-flint 0.9.0 and SymPy 1.14.0:
    python verify.py
"""

from fractions import Fraction
from functools import lru_cache
from math import isqrt

from flint import acb, arb, ctx
from sympy import Matrix, Rational, binomial


LEFT_ENDPOINTS = [
    47338829, 76170593, 103070721, 120905884, 145045445, 163818832,
    178901731, 203238085, 214940952, 233679149, 251575522, 267024897,
    283239096, 296142421, 320405682, 324286195, 344396641, 359724200,
    376676067, 384649618, 404670159, 416223128, 433679140, 442515649,
    459340549, 473686619, 487720354, 498596666, 511940655, 531535614,
    538534604, 550873446, 570069772, 577154072, 595206536, 604175246,
    617052665, 632808995, 644705244, 656797506, 666499321, 683817929,
    695164458, 706518251, 717295088, 731920727, 746277060, 754247694,
    768467298, 779854665, 792401817, 805810862, 818500030, 824670454,
]


@lru_cache(None)
def coefficient(n):
    if n == 1:
        return 1
    if n % 2 == 0 or n % 3 == 0:
        return 0
    p = next((j for j in range(5, isqrt(n) + 1, 2) if n % j == 0), n)
    m, power = n, 0
    while m % p == 0:
        m //= p
        power += 1
    if p % 3 == 2:
        trace, det = 0, -1
    else:
        trace = 2 if pow(2, (p - 1) // 3, p) == 1 else -1
        det = 1
    previous, current = 1, trace
    for _ in range(1, power):
        previous, current = current, trace * current - det * previous
    return current * coefficient(m)


def as_ball(r):
    return arb(int(r.numerator)) / int(r.denominator)


def evaluator(precision, cutoff):
    ctx.prec = precision
    pi = arb.pi()
    C = arb(108).sqrt() / (2 * pi)
    terms = []
    for n in range(1, cutoff + 1):
        if coefficient(n):
            x = arb(n) / C
            terms.append((x, x.log(), coefficient(n)))
    tail = 2 * (C + C**2) * (-(cutoff + 1) / C).exp() / (1 - (-1 / C).exp())

    @lru_cache(None)
    def completed(sigma, ordinate):
        s = acb(arb(sigma.numerator) / sigma.denominator,
                arb(ordinate.numerator) / ordinate.denominator)
        total = acb(0)
        for x, log_x, value in terms:
            first = (-log_x * s).exp() * acb(x).gamma_upper(s)
            if sigma == Fraction(1, 2):
                second = first.conjugate()
            else:
                second = (-log_x * (1 - s)).exp() * acb(x).gamma_upper(1 - s)
            total += value * (first + second)
        return total + acb(arb(0, tail), arb(0, tail))
    return pi, C, completed


def certify_count():
    pi, C, completed = evaluator(384, 288)
    half = Fraction(1, 2)
    height = 50
    assert completed(half, height).real > 0
    degree = 35
    mesh = [Rational(j, degree) for j in range(degree + 1)]
    basis = Matrix([[binomial(degree, k) * x**k * (1 - x)**(degree - k)
                     for k in range(degree + 1)] for x in mesh])
    inverse = basis.inv()
    rows = [[as_ball(v) for v in inverse.row(k)] for k in range(degree + 1)]
    error = arb(1200) / 16**36
    least = arb(1)
    for block in range(24):
        values = []
        for j in range(degree + 1):
            sigma = half + Fraction(block, 16) + Fraction(j, 16 * degree)
            rotation = acb(0, -pi * as_ball(sigma - half) / 2).exp()
            values.append((rotation * completed(sigma, height)).real)
        for row in rows:
            coefficient_ball = sum(w * val for w, val in zip(row, values))
            assert coefficient_ball > error
            least = min(least, coefficient_ball.lower())

    at_right = completed(2, height)
    s_right = acb(2, height)
    raw = at_right * (-s_right * C.log()).exp() * s_right.rgamma()
    rotated = acb(0, -3 * pi / 4).exp() * at_right
    assert raw.real > 0 and rotated.real > 0
    argument = (height * C.log() + s_right.lgamma().imag + raw.arg()
                - 3 * pi / 4 - rotated.arg())
    count = argument / pi
    assert count.contains(54) and count > 53.99 and count < 54.01

    denominator = 2**24
    assert len(LEFT_ENDPOINTS) == 54 and 0 < LEFT_ENDPOINTS[0] and LEFT_ENDPOINTS[-1] + 1 < height * denominator
    assert all(a + 1 < b for a, b in zip(LEFT_ENDPOINTS, LEFT_ENDPOINTS[1:]))
    for a in LEFT_ENDPOINTS:
        lo = completed(half, Fraction(a, denominator)).real
        hi = completed(half, Fraction(a + 1, denominator)).real
        assert lo * hi < 0
    print("top lower bound before interpolation error", least)
    print("interpolation error", error)
    print("winding count", count)
    print("disjoint sign-change intervals", len(LEFT_ENDPOINTS))


def certify_velocities():
    _, _, completed = evaluator(640, 512)
    half = Fraction(1, 2)
    step = Fraction(1, 2**88)
    step_ball = as_ball(step)
    width = Fraction(1, 2**175)
    rounding = 2**200

    def central_value(x):
        return completed(half, x).real

    def newton_guess(x):
        left, value, right = [central_value(u) for u in (x-step, x, x+step)]
        change = value * (2 * step_ball) / (right-left)
        integer = ((as_ball(x)-change).mid() * rounding).floor().unique_fmpz()
        assert integer is not None
        return Fraction(int(integer), rounding)

    velocities = []
    for seed in LEFT_ENDPOINTS:
        lo, hi = Fraction(seed, 2**24), Fraction(seed+1, 2**24)
        x = (lo+hi)/2
        for _ in range(4):
            x = newton_guess(x)
        assert lo < x-width < x+width < hi
        assert central_value(x-width) * central_value(x+width) < 0
        left, value, right = [central_value(u) for u in (x-step, x, x+step)]
        error1 = 224 * step_ball**2 / (1-step_ball**2) + 448 * as_ball(width)
        error2 = 448 * step_ball**2 / (1-step_ball**2) + 1344 * as_ball(width)
        first = (right-left)/(2*step_ball) + arb(0, error1)
        second = (right-2*value+left)/step_ball**2 + arb(0, error2)
        assert first > 0 or first < 0
        velocity = second / first
        assert velocity > 0 or velocity < 0
        velocities.append(velocity)
        print(len(velocities), float(x), velocity.str(14))

    gaps = [right-left for left, right in zip(velocities, velocities[1:])]
    assert len(velocities) == 54 and len(gaps) == 53
    assert sum(v < 0 for v in velocities) == 50
    assert sum(v > 0 for v in velocities) == 4
    assert sum(g < 0 for g in gaps) == 29 and sum(g > 0 for g in gaps) == 24
    assert min(abs(g).lower() for g in gaps) > arb("0.00969")
    print("first and third gap derivatives", gaps[0].str(14), gaps[2].str(14))


if __name__ == "__main__":
    certify_count()
    certify_velocities()
