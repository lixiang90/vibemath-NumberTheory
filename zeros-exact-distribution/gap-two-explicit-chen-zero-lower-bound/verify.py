"""Reproduce the bounds, oriented sign statistics and finite zero comparisons.

Run from the published tree with Python 3, python-flint, NumPy and SciPy.
Only the analytic constants are certified using interval arithmetic.
"""

from fractions import Fraction as Q
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import log, sqrt
from pathlib import Path

import numpy as np
from flint import acb, arb, ctx
from scipy.integrate import quad
from scipy.special import digamma


ROOT = Path(__file__).resolve().parents[1]
A, D = Q(3, 11), Q(2, 7)
ANCHOR, END = Q(561, 2000), Q(117, 400)
GLOBAL_INPUTS = [(Q("0.296"), Q("0.0766")), (Q("0.3000"), Q("0.0129")),
                 (Q("0.3008"), Q("0.00025")), (Q("0.300816"), Q("0.000002"))]


def article(slug):
    spec = spec_from_file_location(slug.replace("-", "_"), ROOT / slug / "verify.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ball(value):
    value = Q(value)
    return arb(value.numerator) / value.denominator


def marked_level(u):
    return min((1 + u)/2, (2 - 3*u)/2)


def check_factorizations():
    # Exact arithmetic at the limiting levels, including both sides of the switch.
    for i in range(161):
        u = Q(1, 7) + Q(i, 1120)
        theta = marked_level(u)
        residual = theta-u
        for eta in [Q(k, 240) for k in range(81)] + [residual]:
            if eta <= residual:
                q1, q2, q3 = eta, u, residual-eta
            else:
                assert eta >= u
                q1, q2, q3 = eta, theta-eta, 0
            assert min(q1, q2, q3) >= 0 and q1+q2+q3 == theta
            assert q1 <= eta and 2*eta+q2+2*q3 <= 1 and 2*eta+5*q2+2*q3 <= 2
    print("Marked-prime modulus allocations passed at 161 rational parameters.")


@lru_cache(maxsize=None)
def triple_cost(c):
    """Integrate in the smaller prime first, then certify the remaining integral."""
    lo, hi = Q(1, 3), 1-2*c
    cuts = sorted({lo, hi} | {t for t in ((1-c)/2, Q(3, 8), Q(2, 5)) if lo < t < hi})
    total = arb(0)
    for left, right in zip(cuts, cuts[1:]):
        mid = (left+right)/2

        def integrand(w, analytic):
            small = 1-2*w if mid < (1-c)/2 else ball(c)
            ratio = (1-w-small)/small
            if mid < Q(3, 8):
                level = (2+3*w)/5
            elif mid < Q(2, 5):
                level = 1-w
            else:
                level = ball(Q(3, 5))
            return ratio.log(analytic=analytic)/(w*(1-w)*level)

        value = acb.integral(integrand, ball(left), ball(right))
        assert 0 in value.imag and value.real.rad() < arb("1e-27")
        total += value.real
    return total


def full_coefficient(c, r):
    assert A <= c < D and c <= r <= Q(1, 3)
    base = ball(Q(7, 4))*arb(3).log() - 2*arb(2).log() - ball(3*c/(2-5*c)).log()
    strip = ball(Q(5, 3))*ball(r*(1-c)/(c*(1-r))).log()
    return base-triple_cost(c)-strip


def check_integral_reduction():
    """Independently integrate the ordered two-dimensional prime-factor region."""
    for c in (A, Q(7, 25), ANCHOR):
        def inner(u):
            right = (1-u)/2
            cuts = [v for v in (1-u-3/8, 1-u-2/5) if u < v < right]
            def integrand(v):
                w = 1-u-v
                level = max(3/5, min(1-w, (2+3*w)/5))
                return 1/(u*v*w*level)
            return quad(integrand, u, right, points=cuts)[0]
        direct = quad(inner, float(c), 1/3, points=[3/10, 5/16])[0]
        assert abs(direct-float(triple_cost(c))) < 1e-13


def f(r):
    return 2*r*(1-r)


def mixed_bound(c2, include_global):
    grid = {A + (END-A)*i/128 for i in range(129)}
    if include_global:
        grid.update(r for r, _ in GLOBAL_INPUTS)
    result, last = arb(0), Q(0)
    for r in sorted(grid):
        positive = c2*full_coefficient(min(r, ANCHOR), r)/4
        positive = positive.max(arb(0))
        other = positive
        if include_global:
            for theta, lower in GLOBAL_INPUTS:
                if r <= theta:
                    other = other.max(c2*ball(lower)/8)
        result += ball(f(r)-last)*(positive+other)
        last = f(r)
    return result


def constant_checks(primes):
    ctx.prec = 144
    partial = arb(1)
    for p in primes[primes > 2]:
        partial *= 1 - 1/arb(int(p)-1)**2
    # The product is taken up to 10^6. Every omitted factor is bounded via integers.
    c2 = partial.union(partial*(1 - arb(1)/999999))
    check_integral_reduction()
    targets = [A, Q(11, 40), Q(7, 25), D, Q(36, 125), Q(29, 100), END]
    floors = ["0.026434", "0.023652", "0.017144", "0.009396", "0.006323",
              "0.003645", "0.000314"]
    for r, floor in zip(targets, floors):
        value = c2*full_coefficient(min(r, ANCHOR), r)/4
        assert value > arb(floor) and value < arb(floor)+arb("0.000001")
        print(f"r={str(r):>7}, one orientation > {floor}")
    external_floors = ["0.0063210", "0.0010645", "0.0000206300", "0.00000016504"]
    for (theta, lower), floor in zip(GLOBAL_INPUTS, external_floors):
        value = c2*ball(lower)/8
        assert value > arb(floor)
        print(f"external global input, r={float(theta):.6f}, C2*b/8 = {value}")
    mixed = mixed_bound(c2, False)
    extended = mixed_bound(c2, True)
    assert mixed > arb("0.0214420") and mixed < arb("0.0214422")
    assert extended > arb("0.0214754") and extended < arb("0.0214756")
    assert 2*f(END) == Q("0.827775") and 2*f(GLOBAL_INPUTS[-1][0]) == Q("0.841302936576")
    print("Mixed lower constants, base / with the four external inputs:", mixed, extended)
    # J=0, and each B has density 0.9 on [1/3,2/5]. The simple bound uses
    # H(c) <= H(A), lambda <= 5/8 and u(1-u) <= 1/4 throughout.
    model_mass = Q(9, 10)*(Q(2, 5)-Q(1, 3))
    rough_at_a = full_coefficient(A, A)+triple_cost(A)
    assert model_mass == Q(3, 50) and c2*rough_at_a/4 < arb("0.048") < ball(model_mass)
    assert c2*ball(GLOBAL_INPUTS[0][1])/8 < ball(model_mass)
    assert c2/(4*ball(Q(5, 8))*ball(Q(1, 4))) > 1
    print("Zero-twin model fits both directional tails and density caps.")


def exact_sign_check():
    # Unrelated background primes are present in every linear observation.
    signs = np.array(list(product((-1, 1), repeat=7)), dtype=int)
    obs = signs @ np.array([[3, 5, 7, 11, 13, 17, 19],
                            [2, 7, 3, 13, 5, 19, 11],
                            [5, 3, 11, 2, 19, 13, 7]]).T
    for count in range(1, 6):
        mark = np.prod(signs[:, :count], axis=1)
        groups = [obs[mark == sign] for sign in (1, -1)]
        d1 = Q(int(groups[0][:, 0].sum()-groups[1][:, 0].sum()), len(signs))
        def covariance(group):
            x, y = group[:, 1], group[:, 2]
            return Q(int(x@y), len(group)) - Q(int(x.sum())*int(y.sum()), len(group)**2)
        d2 = (covariance(groups[0])-covariance(groups[1]))/2
        assert d1 == (3 if count == 1 else 0)
        assert d2 == (2*3+7*5 if count == 2 else 0)
    print("128 exact sign vectors: first and second endpoint contrasts checked.")


def finite_windows(physical):
    sizes = (8192, 32768)
    limit = int(max(sizes)**1.5/2)+2
    factors = physical.smallest_factors(limit)
    witness = physical.test_witness(limit)
    print("Finite four-sign witness:", witness)
    for y in sizes:
        twins, semis = physical.admitted_rows(y, factors, witness)
        n, p, s, t = semis.T
        norm = digamma(int(y**1.5/2)+1)-digamma(y)

        def mass(rows):
            return float(np.sum(np.log(rows)*np.log(rows+2)/rows)/norm)

        twin_mass = mass(twins)
        raw = np.log(s)/np.log(s*t)
        for r in (A, END):
            keep = np.maximum(float(A), raw) >= float(r)
            plus = mass(n[keep & (p == n)])
            minus = mass(n[keep & (p == n+2)])
            assert abs(plus+minus-mass(n[keep])) < 1e-13
            print(f"Y={y}, r={float(r):.6f}, twins={twin_mass:.7f}, "
                  f"forward/backward = {twin_mass+plus:.7f}, {twin_mass+minus:.7f}")


def finite_zeros(trace):
    conductors = trace.primes_up_to(12000)
    conductors = conductors[(conductors>800) & (conductors%48 == 23)]
    heights = np.array([0, 0.09, 0.18])
    assert max(heights)+trace.RADIUS < trace.TOP
    rows, root_count = [], 0
    for d in conductors:
        roots = trace.QuadraticTheta(int(d)).zeros()
        root_count += len(roots)
        diffs = roots[:, None]-heights
        kern = trace.zero_kernel(diffs)*(abs(diffs) <= trace.RADIUS)
        rows.append(np.sum(kern, axis=0)-log(d))
    zero_side = -0.5*trace.centered(np.array(rows))
    primes = trace.primes_up_to(int(np.exp(trace.B+4*trace.BOX)))
    primes = primes[primes>3]
    characters = np.column_stack([trace.character_column(conductors, int(p)) for p in primes])
    coefs = np.vstack([trace.prime_coefficient(int(p), heights) for p in primes])
    prime_side = trace.centered(characters @ coefs)

    def directed_read(data, prime_end, other):
        def groups(endpoint):
            # Euler's criterion uses the prime conductor even for a composite endpoint.
            residues = [pow(int(endpoint), (int(d)-1)//2, int(d)) for d in conductors]
            marks = np.array([1 if a == 1 else -1 if a else 0 for a in residues])
            return [data[marks == sign] for sign in (1, -1)]
        def mean_contrast(group, i):
            return (group[0][:, i].mean()-group[1][:, i].mean())/2
        prime = sqrt(prime_end)*mean_contrast(groups(prime_end), 0)/log(prime_end)
        group = groups(other)
        first = sqrt(other)*mean_contrast(group, 1)/log(other)
        covs = [np.mean(np.prod(trace.centered(g[:, 1:3]), axis=1)) for g in group]
        second = sqrt(other)*(covs[0]-covs[1])/(2*log(other)**2)
        return np.array([prime*first, prime*second])

    errors = []
    for n in (65, 71, 77, 83, 89):
        for p, other in ((n, n+2), (n+2, n)):
            errors.append(abs(directed_read(zero_side, p, other)
                              - directed_read(prime_side, p, other)))
    maximum = np.max(errors, axis=0)
    assert np.shape(errors) == (10, 2) and max(maximum) < 0.001
    print(f"{len(conductors)} conductors / {root_count} numerical zeros / "
          f"20 products, max zero-prime discrepancies (twins / semiprimes) = "
          f"{maximum[0]:.9g} / {maximum[1]:.9g}")


if __name__ == "__main__":
    trace = article("gap-two-quadratic-conductor-covariance")
    check_factorizations()
    constant_checks(trace.primes_up_to(1_000_000))
    exact_sign_check()
    finite_windows(article("gap-two-coherent-zero-heights"))
    finite_zeros(trace)
