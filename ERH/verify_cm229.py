"""Zero and height certificates for F=Q[x]/(x^3-4x+1) and E=F(i).

Requires python-flint 0.9.0 and verify_tetra3136.py in the same directory.
Run --regenerate to rebuild the dyadic interval table from the two Euler products.
"""

import argparse
import csv
from math import isqrt
from pathlib import Path

from flint import arb, fmpz_poly, nmod_poly

from verify_tetra3136 import (
    CUTOFF, HALF, HEADER, PI, Completed, GammaFactors,
    half_winding, locate, verify_signs,
)


FILE = Path(__file__).with_name("cm229-zeros64.csv")
FACTORS = {"real": GammaFactors(229, (0, 0)), "cm_twist": GammaFactors(3664, (1, 1))}


def coefficients(limit):
    """Local series of zeta_F/zeta and its primitive twist by chi_{-4}."""
    is_prime = bytearray([1]) * (limit + 1)
    is_prime[0:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if is_prime[p]:
            is_prime[p*p::p] = bytes((limit - p*p) // p + 1)

    real = [1] * (limit + 1)
    real[0] = 0
    for p in range(2, limit + 1):
        if not is_prime[p]:
            continue
        poly = nmod_poly([1, -4, 0, 1], p)
        degrees = tuple(sorted(f.degree() for f, _ in poly.factor()[1]))
        if p == 229 and degrees != (1, 1):
            raise ArithmeticError("Ramified local factor at 229")
        e, power = 1, p
        while power <= limit:
            if p == 229:
                value = 1
            elif degrees == (1, 1, 1):
                value = e + 1
            elif degrees == (1, 2):
                value = int(e % 2 == 0)
            elif degrees == (3,):
                value = (1, -1, 0)[e % 3]
            else:
                raise ArithmeticError(f"Unexpected factorization at {p}: {degrees}")
            for k in range(1, limit // power + 1):
                if k % p:
                    real[k * power] *= value
            power *= p
            e += 1
    chi4 = (0, 1, 0, -1)
    return {"real": real, "cm_twist": [a * chi4[n % 4] for n, a in enumerate(real)]}


def completed_on_real_axis(function, s):
    mag = (function.rule.frequency * (s - HALF)).exp()
    return ((function.poly(mag) + function.poly(1 / mag)).real
            - function.coeff[0] + function.rest)


def log_derivative_at_one(function):
    """Certified logarithmic derivative of either Artin factor at one."""
    h = arb(2) ** -64
    slope = (completed_on_real_axis(function, 1 + h) - completed_on_real_axis(function, 1 - h)) / (2*h)
    cauchy = function.factors.strip_bound(2, arb(1)) * h*h / (1 - h*h)
    log_derivative = (slope + arb(0, cauchy)) / completed_on_real_axis(function, arb(1))
    gamma_term = sum(((arb(1) + shift) / 2).digamma() / 2 for shift in function.factors.shifts)
    return log_derivative - function.factors.log_scale - gamma_term


def height_difference(function):
    """Stable Faltings height of a mixed type minus that of an induced type."""
    log_derivative = log_derivative_at_one(function)
    qtr = arb(1) / 4
    beta_log_derivative = (arb.const_euler() + (2 * PI).log()
                           - 2 * (qtr.lgamma() - (3 * qtr).lgamma()))
    difference = (2 * log_derivative - 4 * beta_log_derivative + arb(229).log()) / 3
    return log_derivative, beta_log_derivative, difference


def unit_errors():
    roots = fmpz_poly([1, -4, 0, 1]).complex_roots()
    if len(roots) != 3 or any(mult != 1 or not z.imag.is_zero() for z, mult in roots):
        raise ArithmeticError("Unit calculation requires three distinct real roots")
    values = sorted(z.real for z, _ in roots)
    logs = [(abs(x).log(), abs(x - 2).log()) for x in values]
    determinant = logs[0][0] * logs[1][1] - logs[1][0] * logs[0][1]
    bounds = [(abs(s) + abs(t)) / 2 for s, t in logs]
    if not determinant > 0 or not all(b < arb("1.297") for b in bounds):
        raise ArithmeticError("Unit lattice or coordinate bounds")
    if not sum(bounds) < arb("3.342"):
        raise ArithmeticError("Sum of the unit coordinate bounds")
    return bounds


def read_rows(path=FILE):
    blocks = {name: [] for name in FACTORS}
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle)
        if tuple(next(reader)) != ("factor", *HEADER):
            raise ValueError("Unexpected interval header")
        for name, lo, hi, den in reader:
            blocks[name].append((int(lo), int(hi), int(den)))
    return blocks


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true")
    regenerate = parser.parse_args().regenerate
    blocks = {} if regenerate else read_rows()
    a = coefficients(CUTOFF)
    print("Unit coordinate bounds:", [b.str(9) for b in unit_errors()])
    for name, factors in FACTORS.items():
        function = Completed(a[name], factors, radius=512)
        if regenerate:
            blocks[name] = locate(function)
        rows = blocks[name]
        verify_signs(function, rows)
        half, margin = half_winding(function)
        if (2 * half).unique_fmpz() != 2 * len(rows):
            raise ArithmeticError(f"{name}: winding {2*half} vs {len(rows)} intervals")
        print(name, "positive zeros", len(rows), "half winding", half.str(27),
              "top margin", margin.str(8), "first and last", rows[0], rows[-1])
        if name == "cm_twist":
            derivative, beta, height = height_difference(function)
            if not arb("1.1323555935686836618893") < height < arb("1.1323555935686836618894"):
                raise ArithmeticError(f"Height difference outside the reported interval: {height}")
            print("Odd log derivative at 1:", derivative.str(26))
            print("Dirichlet beta log derivative at 1:", beta.str(26))
            print("Mixed minus induced height:", height.str(26))
            induced = 3 * (beta / 2 + arb(4).log() / 4
                           - arb.const_euler() / 2 - (2 * PI).log())
            print("Induced and mixed heights:", induced.str(26), (induced + height).str(26))
    intervals = sorted((lo, hi, name) for name, rows in blocks.items() for lo, hi, _ in rows)
    adjacent = list(zip(intervals, intervals[1:]))
    if any(left[1] >= right[0] for left, right in adjacent):
        raise ArithmeticError("Two isolating intervals overlap")
    cross_gap = min(right[0] - left[1] for left, right in adjacent if left[2] != right[2])
    print("Smallest cross-factor gap between interval endpoints (units of 2^-24):", cross_gap)
    if regenerate:
        with FILE.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(("factor", *HEADER))
            for name, rows in blocks.items():
                writer.writerows((name, *row) for row in rows)


if __name__ == "__main__":
    main()
