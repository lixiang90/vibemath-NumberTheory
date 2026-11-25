"""Verify every zero of the level-23 Artin factor up to height 1024 (python-flint 0.9.0)."""

import argparse
import csv
import sys
from math import factorial

from flint import acb, acb_poly, arb, ctx


ctx.prec = 288
HEIGHT = 1024
MESH = 1024
LAST = 9 * MESH
HALF = arb(1) / 2
PI = arb.pi()
C = arb(23).sqrt() / (2 * PI)
ALPHA = PI / 2 - arb(1) / 16
PHASE = acb(0, ALPHA / 2).exp()
STRIP = arb(1) / 32
NORM = 2**29
FREQUENCY = 2 * PI * MESH
ALIAS = 2 * NORM * (-(FREQUENCY - HEIGHT) * STRIP).exp() / (1 - (-FREQUENCY * STRIP).exp())
edge = STRIP * arb(9).exp()
TAIL = 8 * (edge + 1) * (-edge).exp() / STRIP**4
ERROR = arb(0, ALIAS + TAIL)


def make_polynomial():
    ray = acb(0, 1) * acb(0, ALPHA).exp() / arb(23).sqrt()
    coefficients = []
    for j in range(LAST + 1):
        tau = ray * (arb(j) / MESH).exp()
        value = tau.modular_eta() * (23 * tau).modular_eta()
        coefficients.append(value / 2 if j == 0 else value)
    return acb_poly(coefficients)


def top_value(top_poly, sigma):
    # The height phase is already in the coefficients. The argument is real.
    left = PHASE * top_poly((arb(sigma) / MESH).exp())
    if sigma == HALF:
        return acb(2 * left.real / MESH + ERROR)
    dual = PHASE * top_poly((arb(1 - sigma) / MESH).exp())
    return (left + dual.conjugate()) / MESH + acb(ERROR, ERROR)


def line_intervals(poly, grid):
    # Chirp convolution evaluates the same trapezoid at every t=m/grid.
    length = HEIGHT * grid
    chirps = [acb(0, arb(k * k) / (2 * MESH * grid)).exp()
              for k in range(max(length, LAST) + 1)]
    coeff = [poly[j] * (arb(j) / (2 * MESH)).exp() * chirps[j]
             for j in range(LAST + 1)]
    kernel = [chirps[abs(k)].conjugate() for k in range(-LAST, length + 1)]
    convolution = acb_poly(coeff) * acb_poly(kernel)
    intervals = []
    previous = None
    for m in range(length + 1):
        value = 2 * (PHASE * chirps[m] * convolution[m + LAST]).real / MESH + ERROR
        if value > 0:
            sign = 1
        elif value < 0:
            sign = -1
        else:
            raise ArithmeticError(f"Unresolved sign at {m}/{grid}")
        if previous is not None and previous != sign:
            intervals.append((m - 1, m, grid))
        previous = sign
    return intervals


def winding_count(poly):
    mesh, degree = 256, 40
    segments = 3 * mesh // 2
    top_poly = acb_poly([poly[j] * acb(0, arb(HEIGHT * j) / MESH).exp()
                         for j in range(LAST + 1)])
    values = [top_value(top_poly, HALF + arb(j) / mesh) for j in range(segments + 1)]
    rows = [values]
    for _ in range(degree):
        row = rows[-1]
        rows.append([b - a for a, b in zip(row, row[1:])])
    remainder = arb(NORM * factorial(degree)) / mesh**(degree + 1)
    change, least = arb(0), arb(1)
    for j in range(segments):
        forward = j + degree <= segments
        anchor = j if forward else j + 1
        value = values[anchor]
        real_coordinate = abs(value.real) > abs(value.imag)
        margin = abs(value.real if real_coordinate else value.imag) - remainder
        for k in range(1, degree + 1):
            delta = rows[k][anchor if forward else anchor - k]
            margin -= abs(delta.real if real_coordinate else delta.imag) / k
        if not margin > 0:
            raise ArithmeticError(f"Unresolved top edge segment {j}: {margin}")
        least = min(least, margin.lower())
        change += (values[j + 1] / values[j]).arg()
    s = acb(2, HEIGHT)
    unrotate = (acb(0, ALPHA) * (s - HALF)).exp()
    l_value = values[-1] * unrotate * (-s * C.log()).exp() * s.rgamma()
    if not l_value.real > 0:
        raise ArithmeticError("Unresolved Euler-product argument")
    angle = HEIGHT * C.log() + s.lgamma().imag + l_value.arg() - 3 * ALPHA / 2 - change
    return angle / PI, least


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", action="store_true", help="emit the certified intervals as CSV")
    args = parser.parse_args()
    if not (ALPHA + STRIP).cos() / C > STRIP:
        raise ArithmeticError("Invalid decay rate at the edge of the integration strip")
    strip_norm = 4 / ((1 - (-STRIP).exp())**2 * STRIP**3)
    if not strip_norm < NORM:
        raise ArithmeticError("Integration strip norm exceeds the assumed bound")
    poly = make_polynomial()
    winding, margin = winding_count(poly)
    intervals = line_intervals(poly, 64)
    number = len(intervals)
    if not (winding.contains(number) and winding > number - HALF and winding < number + HALF):
        raise ArithmeticError(f"{number} line intervals vs. strip winding {winding}")
    if args.csv:
        table = csv.writer(sys.stdout, lineterminator="\n")
        table.writerow(("lower_numerator", "upper_numerator", "denominator"))
        table.writerows(intervals)
    else:
        print("Quadrature alias bound", ALIAS.str(6))
        print("Quadrature truncation bound", TAIL.str(6))
        print("Top edge coordinate margin", margin.str(6))
        print("Winding enclosure", winding.str(65))
        print("Certified positive-height intervals", number)
        print("First and last", intervals[0], intervals[-1])


if __name__ == "__main__":
    main()
