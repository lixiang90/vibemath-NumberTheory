"""Certify the critical-strip zeros of the level-23 Artin factor (python-flint 0.9.0)."""

import argparse
import csv
import sys
from math import factorial

from flint import acb, acb_poly, arb, ctx


ctx.prec = 288
HALF = arb(1) / 2
PI = arb.pi()
C = arb(23).sqrt() / (2 * PI)


class Quadrature:
    """Parameters and rigorous error bounds for the rotated Mellin integral."""

    def __init__(self, height):
        cutoff, norm_power = {1024: (9, 29), 8192: (12, 44)}[height]
        self.height = self.mesh = height
        self.last = cutoff * self.mesh
        self.alpha = PI / 2 - arb(64) / height
        self.phase = acb(0, self.alpha / 2).exp()
        strip = arb(32) / height
        self.norm = 2**norm_power
        if not (self.alpha + strip).cos() / C > strip:
            raise ArithmeticError("Invalid decay rate at the edge of the integration strip")
        strip_norm = 4 / ((1 - (-strip).exp())**2 * strip**3)
        if not strip_norm < self.norm:
            raise ArithmeticError("Integration strip norm exceeds the assumed bound")
        frequency = 2 * PI * self.mesh
        self.alias = (2 * self.norm * (-(frequency - height) * strip).exp()
                      / (1 - (-frequency * strip).exp()))
        edge = strip * arb(cutoff).exp()
        self.tail = 8 * (edge + 1) * (-edge).exp() / strip**4
        self.error = arb(0, self.alias + self.tail)


DEFAULT = Quadrature(1024)


def make_polynomial(rule=DEFAULT):
    ray = acb(0, 1) * acb(0, rule.alpha).exp() / arb(23).sqrt()
    coefficients = []
    for j in range(rule.last + 1):
        tau = ray * (arb(j) / rule.mesh).exp()
        value = tau.modular_eta() * (23 * tau).modular_eta()
        coefficients.append(value / 2 if j == 0 else value)
    return acb_poly(coefficients)


def top_value(top_poly, sigma, rule=DEFAULT):
    # The height phase is already in the coefficients. The argument is real.
    left = rule.phase * top_poly((arb(sigma) / rule.mesh).exp())
    if sigma == HALF:
        return acb(2 * left.real / rule.mesh + rule.error)
    dual = rule.phase * top_poly((arb(1 - sigma) / rule.mesh).exp())
    return (left + dual.conjugate()) / rule.mesh + acb(rule.error, rule.error)


def line_intervals(poly, grid, rule=DEFAULT):
    # Chirp convolution evaluates the same trapezoid at every t=m/grid.
    length = rule.height * grid
    chirps = [acb(0, arb(k * k) / (2 * rule.mesh * grid)).exp()
              for k in range(max(length, rule.last) + 1)]
    coeff = [poly[j] * (arb(j) / (2 * rule.mesh)).exp() * chirps[j]
             for j in range(rule.last + 1)]
    kernel = [chirps[abs(k)].conjugate() for k in range(-rule.last, length + 1)]
    convolution = acb_poly(coeff) * acb_poly(kernel)
    intervals = []
    previous = None
    for m in range(length + 1):
        value = 2 * (rule.phase * chirps[m] * convolution[m + rule.last]).real / rule.mesh + rule.error
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


def winding_count(poly, rule=DEFAULT):
    mesh, degree = 256, 40
    segments = 3 * mesh // 2
    top_poly = acb_poly([poly[j] * acb(0, arb(rule.height * j) / rule.mesh).exp()
                         for j in range(rule.last + 1)])
    values = [top_value(top_poly, HALF + arb(j) / mesh, rule) for j in range(segments + 1)]
    rows = [values]
    for _ in range(degree):
        row = rows[-1]
        rows.append([b - a for a, b in zip(row, row[1:])])
    remainder = arb(rule.norm * factorial(degree)) / mesh**(degree + 1)
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
    s = acb(2, rule.height)
    unrotate = (acb(0, rule.alpha) * (s - HALF)).exp()
    l_value = values[-1] * unrotate * (-s * C.log()).exp() * s.rgamma()
    if not l_value.real > 0:
        raise ArithmeticError("Unresolved Euler-product argument")
    angle = rule.height * C.log() + s.lgamma().imag + l_value.arg() - 3 * rule.alpha / 2 - change
    return angle / PI, least


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", action="store_true", help="emit the certified intervals as CSV")
    parser.add_argument("--height", type=int, choices=(1024, 8192), default=1024,
                        help="upper height of the certified rectangle (default: 1024)")
    args = parser.parse_args()
    rule = Quadrature(args.height)
    poly = make_polynomial(rule)
    winding, margin = winding_count(poly, rule)
    intervals = line_intervals(poly, 64, rule)
    number = len(intervals)
    if not (winding.contains(number) and winding > number - HALF and winding < number + HALF):
        raise ArithmeticError(f"{number} line intervals vs. strip winding {winding}")
    if args.csv:
        table = csv.writer(sys.stdout, lineterminator="\n")
        table.writerow(("lower_numerator", "upper_numerator", "denominator"))
        table.writerows(intervals)
    else:
        print("Quadrature alias bound", rule.alias.str(6))
        print("Quadrature truncation bound", rule.tail.str(6))
        print("Top edge coordinate margin", margin.str(6))
        print("Winding enclosure", winding.str(65))
        print("Certified positive-height intervals", number)
        print("First and last", intervals[0], intervals[-1])


if __name__ == "__main__":
    main()
