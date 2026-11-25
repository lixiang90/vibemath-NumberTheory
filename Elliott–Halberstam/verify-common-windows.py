"""Common-window selection, a signed four-character sum, and its norm bound."""
import cmath
import math
from fractions import Fraction
from functools import cache
from itertools import combinations

from flint import acb, arb, ctx, dirichlet_char

from verify import isolate_ordinates


def common_window(upper, ratio, starts):
    """Return [y,2y] avoiding the half-open intervals [a,ratio*a)."""
    starts = set(starts)
    while any(upper / (2 * ratio) < a <= upper for a in starts):
        lower = upper / (2 * ratio)
        starts = {a for a in starts if a <= lower or a > upper}
        upper = lower
    return upper / 2


def check_windows():
    upper = Fraction(72)
    starts = tuple(map(Fraction, (0, 1, 2, 3, 4, 6, 8, 12, 18, 24, 36, 72, 73)))
    checked = 0
    for ratio in (Fraction(3, 2), Fraction(2), Fraction(3)):
        for size in range(len(starts) + 1):
            for selected in combinations(starts, size):
                y = common_window(upper, ratio, selected)
                assert upper / (2 * (2 * ratio)**size) <= y and 2 * y <= upper
                assert all(2 * y < a or ratio * a <= y for a in selected)
                checked += 1
    print(f"exact interval checks: {checked}")


def horizontal_argument(chi, height):
    """Continue arg L from 2+it to 1/2+it with certified open half-planes."""
    @cache
    def endpoint(sigma):
        return chi.l_function(acb(sigma, height))

    def change(left, right):
        mid = (left + right) / 2
        image = chi.l_function(acb(arb(mid, (right - left) / 2), height))
        if image.real > 0 or image.real < 0 or image.imag > 0 or image.imag < 0:
            return (endpoint(left) / endpoint(right)).arg()
        if right - left < 2**-22:
            raise ArithmeticError("horizontal segment contains an unresolved value")
        return change(mid, right) + change(left, mid)

    return endpoint(2).arg() + change(0.5, 2)


def strip_count(chi, height):
    s = acb(0.5, height)
    phase = height * (arb(chi.modulus()) / arb.pi()).log()
    phase += 2 * ((s + chi.parity()) / 2).lgamma().imag
    phase += horizontal_argument(chi, height) - horizontal_argument(chi, -height)
    integer = (phase / arb.pi()).unique_fmpz()
    if integer is None:
        raise ArithmeticError("winding number failed to isolate an integer")
    return int(integer)


def signed_quartets(z, cutoff):
    """Use additive exponents of cyclic character labels modulo a prime."""
    g = len(z)
    all_terms = kept = 0j
    forbidden_powers = {(m, n) for m in range(-cutoff, cutoff + 1)
                        for n in range(-cutoff, cutoff + 1) if (m, n) != (0, 0)}
    for eta in range(g):
        for nu in range(g):
            terms = sum(z[j] * z[(j + eta + nu) % g] *
                        (z[(j + eta) % g] * z[(j + nu) % g]).conjugate()
                        for j in range(g))
            all_terms += terms
            short = any(g // math.gcd(g, m * eta + n * nu) <= cutoff
                        for m, n in forbidden_powers)
            if not short:
                kept += terms
    return all_terms / g**4, kept / g**4


def check_zero_sums():
    p, generator, height, width = 47, 5, 12, 0.30
    g = p - 1
    labels = [pow(generator, j, p) for j in range(g)]
    assert len(set(labels)) == g
    roots = [[] for _ in range(g)]
    for j in range(1, g // 2):
        chi = dirichlet_char(p, labels[j])
        assert chi.conductor() == p and chi.order() > 2
        roots[j] = isolate_ordinates(chi, height, strip_count(chi, height))
        roots[-j] = [-t for t in roots[j]]
    print(f"q={p}, {g-2} nonreal characters, certified zero count={sum(map(len, roots))}")
    for y in (256, 512, 1024):
        z = []
        for ordinates in roots:
            value = 0j
            for t in ordinates:
                rho = 0.5 + 1j * t
                value += cmath.exp(1j * t * math.log(y) + width**2 * rho**2 / 2) / rho
            z.append(value)
        transformed = []
        for a in range(1, p):
            transformed.append(sum(complex(dirichlet_char(p, label)(a)).conjugate() * v
                                   for label, v in zip(labels, z)) / g)
        m2 = sum(abs(v)**2 for v in transformed) / g
        m4 = sum(abs(v)**4 for v in transformed) / g
        expanded, generic = signed_quartets(z, 2)
        err = abs(expanded - m4)
        assert err < 1e-12 and abs(generic.imag) < 1e-12
        assert abs(m2 - sum(abs(v)**2 for v in z) / g**2) < 1e-12
        cap = max(abs(v) for v in transformed)
        assert m4 <= cap**2 * m2 + 1e-14
        print(f"y={y:4d}  M4={m4:.9e}  H={generic.real:+.9e}  "
              f"cap^2*M2={cap**2*m2:.9e}  expansion error={err:.2e}")


def rudin_shapiro(length):
    p, q = [1], [1]
    while len(p) < length:
        p, q = p + q, p + [-v for v in q]
    assert len(p) == length
    return p


def check_fourier_model():
    p, m = 12289, 64
    assert all(p % d for d in range(2, math.isqrt(p) + 1))
    g = p - 1
    signs = rudin_shapiro(m)
    c = 1 / 8
    removed = []
    energy = 0.0
    max_coefficient = 0.0
    for j in range(g):
        w = cmath.exp(2j * math.pi * j / g)
        value = 0j
        for sign in reversed(signs):
            value = sign + w * value
        value *= c / math.sqrt(m)
        max_coefficient = max(max_coefficient, abs(value))
        if g // math.gcd(j, g) <= 4:
            removed.append((j, value))
        else:
            energy += abs(value)**2 / g**2
    fourth = second = peak = 0.0
    for s in range(g):
        v = c * signs[s] / math.sqrt(m) if s < m else 0.0
        v -= sum(z * cmath.exp(-2j * math.pi * j * s / g) for j, z in removed) / g
        peak = max(peak, abs(v))
        fourth += abs(v)**4 / g
        second += abs(v)**2 / g
    assert abs(second - energy) < 1e-14
    assert max_coefficient <= math.sqrt(2) * c + 1e-13
    assert peak >= c / (2 * math.sqrt(m))
    assert fourth >= c**4 / (16 * g * m)
    print(f"model p={p}, m={m}: g*m*M4/c^4={g*m*fourth/c**4:.8f}, "
          f"sqrt(m)*peak/c={math.sqrt(m)*peak/c:.8f}")


def check_full_order_projection():
    p, m, c = 1009, 64, 1 / 8
    assert all(p % d for d in range(2, math.isqrt(p) + 1))
    g = p - 1
    assert g == 2**4 * 3**2 * 7
    divisor_bound = 2**3
    frequencies = [k for k in range(g) if math.gcd(k, g) == 1]
    density = len(frequencies) / g
    roots = [cmath.exp(2j * math.pi * k / g) for k in range(g)]
    signs = rudin_shapiro(m)
    spectrum = [c / math.sqrt(m) * sum(s * roots[(k * j) % g]
                                        for j, s in enumerate(signs)) for k in range(g)]

    def retained_energy(shift):
        return sum(abs(spectrum[(k + shift) % g])**2 for k in frequencies) / g

    shift = max(range(g), key=retained_energy)
    assert retained_energy(shift) >= density * c**2 - 1e-14
    projected = [sum(spectrum[(k + shift) % g] * roots[(-k * j) % g] for k in frequencies) / g
                 for j in range(g)]
    real_parts = ([v.real for v in projected], [v.imag for v in projected])
    values = max(real_parts, key=lambda seq: sum(v**2 for v in seq[:m]))
    fourth = sum(v**4 for v in values) / g
    assert max(map(abs, values)) <= divisor_bound * c / math.sqrt(m) + 1e-14
    assert fourth >= density**4 * c**4 / (4 * g * m) - 1e-18

    excluded = 0.0
    for k in range(g):
        coefficient = sum(v * roots[(k * j) % g] for j, v in enumerate(values))
        assert abs(coefficient) <= math.sqrt(2) * c + 1e-12
        if math.gcd(k, g) != 1:
            excluded = max(excluded, abs(coefficient))
    assert excluded < 1e-12
    print(f"full order p={p}, order={g}, frequencies={len(frequencies)}, shift={shift}, "
          f"g*m*M4/c^4={g*m*fourth/c**4:.8f}, excluded spectrum={excluded:.2e}")


if __name__ == "__main__":
    ctx.prec = 96
    check_windows()
    check_zero_sums()
    check_fourier_model()
    check_full_order_projection()
