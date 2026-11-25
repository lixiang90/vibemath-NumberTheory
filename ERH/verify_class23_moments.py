"""Exact class-character, ideal parity, and quadratic-twist checks for discriminant -23."""

from fractions import Fraction
from math import isqrt

from cubic23 import eta_coefficients, primes_up_to, trace_at_prime


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def convolve(a, b):
    bound = len(a) - 1
    out = [0] * (bound + 1)
    for d in range(1, bound + 1):
        if a[d]:
            for m in range(1, bound // d + 1):
                out[d * m] += a[d] * b[m]
    return out


def ideal_sequences(bound, primes):
    # Local factors of the signed and unsigned squarefree-ideal sums.
    signed, unsigned = [0] * (bound + 1), [0] * (bound + 1)
    signed[1] = unsigned[1] = 1
    for p in primes:
        trace = trace_at_prime(p)
        if p == 23:
            factors = ((p, -1, 1),)
        elif trace == 0:
            factors = ((p * p, -1, 1),)
        else:
            factors = ((p, -trace, trace), (p * p, 1, 1))
        for power, signed_local, unsigned_local in factors:
            for m in range(1, bound // power + 1):
                if m % p:
                    signed[m * power] = signed[m] * signed_local
                    unsigned[m * power] = unsigned[m] * unsigned_local
    return signed, unsigned


def check_coefficients(bound):
    primes = primes_up_to(bound)
    form = eta_coefficients(bound)
    signed, unsigned = ideal_sequences(bound, primes)
    unit = [0] * (bound + 1)
    unit[1] = 1
    require(convolve(signed, form) == unit, "Signed squarefree series times L")

    square_form = [0] * (bound + 1)
    all_squares = [0] * (bound + 1)
    for j in range(1, isqrt(bound) + 1):
        square_form[j*j] = form[j]
        all_squares[j*j] = 1
    require(convolve(unsigned, square_form) == form, "Unsigned series times L(2s)")

    # Independent ideal counts from (1 * chi_{-23})(n).
    legendre = {a*a % 23 for a in range(1, 23)}
    chi = [0 if n % 23 == 0 else (1 if n % 23 in legendre else -1)
           for n in range(bound + 1)]
    ideals = convolve([0] + [1] * bound, chi)
    ramified = unit.copy()
    ramified[23] = 1
    left = convolve(convolve([x*x for x in form], all_squares), ramified)
    require(left == convolve(ideals, form), "Rankin square identity including p=23")
    require(all((b - a) % 2 == 0 for a, b in zip(signed, unsigned)), "Parity integrality")

    # Apply four separate additive valuations in place of log; all checks are exact.
    for ell in (2, 5, 23, 59):
        valuation = [0] * (bound + 1)
        for m in range(ell, bound + 1, ell):
            valuation[m] = 1 + valuation[m // ell]
        derivative = [value * order for value, order in zip(form, valuation)]
        target = [0] * (bound + 1)
        trace = trace_at_prime(ell)
        det = chi[ell]
        previous, current, power = 2, trace, ell
        while power <= bound:
            target[power] = current
            previous, current = current, trace * current - det * previous
            power *= ell
        require(convolve(signed, derivative) == target, f"Logarithmic derivative at {ell}")
        even = [(b + a) // 2 for a, b in zip(signed, unsigned)]
        parity = convolve(even, derivative)
        blind = convolve(unsigned, derivative)
        require(all(2*x-y == z for x, y, z in zip(parity, blind, target)),
                f"Parity decomposition at {ell}")
    print("Exact identities through norm", bound)


def finite_twists():
    # Orient every nonprincipal split prime so its difference is +i*sqrt(3).
    split = [p for p in primes_up_to(101) if p != 23 and pow(p, 11, 23) == 1]
    nonprincipal = [p for p in split if trace_at_prime(p) == -1]
    diagonal = 3 * len(nonprincipal)
    require(diagonal == 3 * (len(split) - sum(trace_at_prime(p) == 2 for p in split)),
            "Class-character diagonal")
    primes = primes_up_to(1_000_000)
    for lower, upper in ((1_000, 10_000), (10_000, 100_000), (100_000, 1_000_000)):
        moduli = [q for q in primes if lower < q <= upper and q % 8 == 3]
        energy = cross = 0
        for q in moduli:
            values = [1 if pow(p, (q-1)//2, q) == 1 else -1 for p in nonprincipal]
            energy += 3 * sum(values)**2
            cross += 6 * sum(a*b for j, a in enumerate(values) for b in values[:j])
        mean = Fraction(energy, len(moduli))
        require(mean - diagonal == Fraction(cross, len(moduli)), "Off-diagonal expansion")
        print(lower, upper, len(moduli), len(split), diagonal, mean, float(mean))


if __name__ == "__main__":
    check_coefficients(32768)
    finite_twists()
