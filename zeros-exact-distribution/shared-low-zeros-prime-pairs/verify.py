"""Certify the finite example Q=53, X=13 with python-flint and scipy."""

from math import ceil, isqrt

from flint import acb, arb, ctx, dirichlet_char
from scipy.optimize import brentq


Q, X, ORDER, HEIGHT = 53, 13, 17, 88


def z_sign(chi, t):
    value = chi.hardy_z(t).real
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise ArithmeticError(f"unresolved Hardy Z sign at {t}")


def positive_roots(chi, step):
    """Each returned ball contains a sign-changing zero in its own grid cell."""
    zeros = []
    lo, left_sign = 0.0, z_sign(chi, 0)
    for j in range(1, ceil(HEIGHT / step) + 1):
        hi = min(j * step, HEIGHT)
        right_sign = z_sign(chi, hi)
        if left_sign != right_sign:
            root = brentq(lambda t: float(chi.hardy_z(t).real.mid()), lo, hi, xtol=1e-12)
            width = 2.0**-32
            assert lo < root - width < root + width < hi
            assert z_sign(chi, root - width) != z_sign(chi, root + width)
            zeros.append(arb(root, width))
        lo, left_sign = hi, right_sign
    return zeros


def horizontal_argument(chi):
    """Continue arg L from 2+iT to 1/2+iT with certified half-plane images."""
    def value(sigma):
        return chi.l_function(acb(sigma, HEIGHT))

    def increment(start, end, first, last, depth):
        middle = (start + end) / 2
        enclosure = value(arb(middle, (start - end) / 2))
        if (enclosure / first).real > 0:
            return (last / first).arg()
        if depth >= 25:
            raise ArithmeticError("could not enclose a horizontal segment")
        center = value(middle)
        return (increment(start, middle, first, center, depth + 1)
                + increment(middle, end, center, last, depth + 1))

    first = value(2)
    return first.arg() + increment(2.0, 0.5, first, value(0.5), 0)


def main():
    ctx.prec = 88
    chars = {i: dirichlet_char(Q, i) for i in range(1, Q) if dirichlet_char(Q, i).parity()}
    assert len(chars) == (Q - 1) // 2 and all(c.is_primitive() for c in chars.values())
    pi = arb.pi()
    theta = (acb(arb(3)/4, arb(HEIGHT)/2).lgamma().imag
             + arb(HEIGHT)/2 * (Q/pi).log())
    arguments = {i: horizontal_argument(chi) for i, chi in chars.items()}
    roots = {i: positive_roots(chi, 0.25) for i, chi in chars.items()}
    for i in chars:
        inverse = pow(i, -1, Q)
        if i > inverse:
            continue
        count = ((2*theta + arguments[i] + arguments[inverse])/pi).unique_fmpz()
        if count is None:
            raise ArithmeticError("argument principle did not isolate an integer")
        for step in (0.125, 0.0625, 0.03125):
            if len(roots[i]) + len(roots[inverse]) == count:
                break
            roots[i] = positive_roots(chars[i], step)
            roots[inverse] = positive_roots(chars[inverse], step)
        assert len(roots[i]) + len(roots[inverse]) == count

    length = arb(X).log()
    halfwidth = (arb(Q-X-1)/X).log()/2
    width = length + 2*halfwidth
    assert HEIGHT >= arb(1).exp()*ORDER/halfwidth
    tail = (40 * arb(X).sqrt() * width * (1+HEIGHT)
            * arb(4*Q*(2+HEIGHT)).log() * arb(-ORDER).exp())
    assert tail < arb(2).log()/8

    def transform(t):
        return (width * acb(0, length*t/2).exp()
                * (width*t/2).sinc() * (halfwidth*t/ORDER).sinc()**ORDER)

    positive = {i: sum((transform(t) for t in row), acb(0)) for i, row in roots.items()}
    spectra = {i: positive[i] + positive[pow(i, -1, Q)].conjugate() for i in chars}
    flags = set()
    max_deviation = arb(0)
    for a in range(2, X+1):
        total = sum((chi(a).conjugate() * spectra[i] for i, chi in chars.items()), acb(0))
        signal = -2 * arb(a).sqrt() * total / (Q-1)
        assert signal.imag.contains(0)
        prime = all(a % p for p in range(2, isqrt(a)+1))
        base = next(p for p in range(2, a+1) if a % p == 0)
        cofactor = a
        while cofactor % base == 0:
            cofactor //= base
        expected = arb(base).log() if cofactor == 1 else arb(0)
        error = abs(signal.real - expected).upper()
        assert error < tail
        max_deviation = max(max_deviation, error)
        threshold = 3 * arb(a).log()/4
        if signal.real > threshold:
            flags.add(a)
        else:
            assert signal.real < threshold
        assert prime == (a in flags)
    goldbach = [sum(n-p in flags for p in flags if 2*p <= n) for n in range(4, X+1, 2)]
    twin = sum(p+2 in flags for p in flags)
    assert goldbach == [1, 1, 1, 2, 1] and twin == 3
    print("Odd characters:", len(chars), "labelled zeros:", 2*sum(map(len, roots.values())))
    print("Tail bound:", tail, "maximum certified deviation:", max_deviation)
    print("Primes:", sorted(flags), "unordered Goldbach counts:", goldbach, "twin pairs:", twin)


if __name__ == "__main__":
    main()
