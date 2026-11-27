"""Compute exponential sums for the two CM fields from rational Euler coefficients."""

from math import gcd, isqrt

from flint import arb, ctx, fmpz

from verify_cm229 import coefficients as cubic_coefficients
from verify_cm7168 import local_data


POINTS = (16, 64, 256, 1024)
BAD = 2*7*229


def product(a, b):
    n = len(a)-1
    out = [0]*(n+1)
    for j in range(1, n+1):
        if a[j]:
            for k in range(1, n//j+1):
                out[j*k] += a[j]*b[k]
    return out


def inverse(a):
    if a[1] != 1:
        raise ValueError("the inverse needs constant Dirichlet coefficient one")
    n = len(a)-1
    out = [0, 1]+[0]*(n-1)
    for j in range(1, n+1):
        if out[j]:
            for k in range(2, n//j+1):
                out[j*k] -= out[j]*a[k]
    return out


def abelian_coefficients(limit, discriminants):
    units = [int(gcd(j, BAD) == 1) for j in range(limit+1)]
    out = [0, 1]+[0]*(limit-1)
    for disc in discriminants:
        character = [int(fmpz(disc).jacobi(j)) if units[j] else 0 for j in range(limit+1)]
        out = product(out, character)
    return out


def power_convolution(a, b, m):
    """Coefficients of A(s) B(ms), including the complete last admissible power."""
    n = len(a)-1
    out = [0]*(n+1)
    for h in range(1, len(b)):
        power = h**m
        if power > n:
            break
        if b[h]:
            for j in range(1, n//power+1):
                out[power*j] += b[h]*a[j]
    return out


def filtered(a, denominator, m):
    out = power_convolution(a, inverse(denominator), m)
    if power_convolution(out, denominator, m) != a:
        raise ArithmeticError("failed to reconstruct every numerator coefficient")
    return out


def coefficients(limit):
    first = cubic_coefficients(limit)
    old = product(first["real"], first["cm_twist"])
    new, _, _ = local_data(limit)
    for p in (2, 7, 229):
        for j in range(p, limit+1, p):
            old[j] = new[j] = 0
    bound = isqrt(limit)
    rational_denominators = {
        ("229", 3): abelian_coefficients(bound, (1, -4, 229, -916)),
        ("7168", 2): abelian_coefficients(bound, (8, 28)),
        ("7168", 4): abelian_coefficients(bound, (1, 56)),
    }
    prefixes = {"229": old[:bound+1], "7168": new[:bound+1]}
    blocks = {}
    for label, a in (("229", old), ("7168", new)):
        for m in (2, 3, 4, 7):
            denominator = rational_denominators.get((label, m), prefixes[label])
            blocks[label, m] = filtered(a, denominator, m)
    return blocks


def main():
    ctx.prec = 192
    blocks = coefficients(96*max(POINTS))
    for x in POINTS:
        top = 96*x
        step = (-arb(1)/x).exp()
        weight = arb(1)
        sums = {name: arb(0) for name in blocks}
        for j in range(1, top+1):
            weight *= step
            for name, coeff in blocks.items():
                if coeff[j]:
                    sums[name] += coeff[j]*weight
        # After n=96X the majorant n^3 exp(-n/X) is decreasing.
        tail = x*arb(-96).exp()*(top**3+3*x*top**2+6*x*x*top+6*x**3)
        if not tail < arb("3e-24"):
            raise ArithmeticError("tail bound is too large")
        print("X", x, "tail", tail.str(6))
        for name, total in sums.items():
            value = total+arb(0, tail)
            print(name, value.str(16))


if __name__ == "__main__":
    main()
