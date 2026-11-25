"""Verify the critical-strip zero count up to height 64 using python-flint 0.9.0."""

import argparse
import csv
import sys
from math import factorial
from pathlib import Path

from flint import acb, arb, ctx

from cubic23 import eta_coefficients


HEIGHT = 64
DENOMINATOR = 2**24
DATA = Path(__file__).with_name("zeros23.csv")
HEADER = ("lower_numerator", "upper_numerator", "denominator")
ctx.prec = 352
PI = arb.pi()
C = arb(23).sqrt() / (2 * PI)
HALF = arb(1) / 2
CUTOFF = 128
TERMS = []
for n, coefficient in enumerate(eta_coefficients(CUTOFF)):
    if coefficient:
        x = arb(n) / C
        TERMS.append((acb(x), x.log(), coefficient))
TAIL = 2 * (C + C**2) * (-(CUTOFF + 1) / C).exp() / (1 - (-1 / C).exp())
ERROR = arb(0, TAIL)


def completed(sigma, t):
    s = acb(sigma, t)
    total = acb(0)
    on_line = sigma == HALF
    for x, log_x, coefficient in TERMS:
        first = (-s * log_x).exp() * x.gamma_upper(s)
        if on_line:
            total += 2 * coefficient * first.real
        else:
            dual = ((s - 1) * log_x).exp() * x.gamma_upper(1 - s)
            total += coefficient * (first + dual)
    return total + acb(ERROR, 0 if on_line else ERROR)


def line_sign(numerator):
    value = completed(HALF, arb(numerator) / DENOMINATOR).real
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise ArithmeticError(f"Indeterminate sign at {numerator}/{DENOMINATOR}")


def locate():
    """Recreate the data using exact dyadic endpoints and certified signs."""
    table = csv.writer(sys.stdout, lineterminator="\n")
    table.writerow(HEADER)
    step = DENOMINATOR // 8
    last_sign = line_sign(0)
    for end in range(step, HEIGHT * DENOMINATOR + 1, step):
        new_sign = line_sign(end)
        if new_sign != last_sign:
            lo, hi = end - step, end
            while hi - lo > 1:
                mid = (lo + hi) // 2
                if line_sign(mid) == last_sign:
                    lo = mid
                else:
                    hi = mid
            table.writerow((lo, hi, DENOMINATOR))
        last_sign = new_sign


def critical_line_count():
    with DATA.open(encoding="utf-8-sig", newline="") as handle:
        rows = csv.reader(handle)
        if tuple(next(rows)) != HEADER:
            raise ValueError("Unexpected interval table header")
        previous = 0
        count = 0
        for row in rows:
            lo, hi, denominator = map(int, row)
            if not (denominator == DENOMINATOR and previous < lo < hi < HEIGHT * denominator
                    and hi == lo + 1):
                raise ValueError(f"Invalid interval {row}")
            if line_sign(lo) == line_sign(hi):
                raise ArithmeticError(f"No certified sign change in {row}")
            previous = hi
            count += 1
    return count


def top_argument():
    """Lift the argument from 1/2+64i to 2+64i, certifying every mesh segment."""
    mesh, degree = 256, 40
    segments = 3 * mesh // 2
    values = [completed(HALF + arb(j) / mesh, HEIGHT) for j in range(segments + 1)]
    differences = [values]
    for _ in range(degree):
        row = differences[-1]
        differences.append([b - a for a, b in zip(row, row[1:])])
    remainder = arb(12 * factorial(degree)) / mesh**(degree + 1)
    change = arb(0)
    least_margin = arb(1)
    for j in range(segments):
        forward = j + degree <= segments
        anchor = j if forward else j + 1
        value = values[anchor]
        use_real = abs(value.real) > abs(value.imag)
        margin = abs(value.real if use_real else value.imag) - remainder
        for k in range(1, degree + 1):
            delta = differences[k][anchor if forward else anchor - k]
            margin -= abs(delta.real if use_real else delta.imag) / k
        if not margin > 0:
            raise ArithmeticError(f"Unresolved top edge segment {j}: {margin}")
        least_margin = min(least_margin, margin.lower())
        change += (values[j + 1] / values[j]).arg()
    return change, values[-1], least_margin


def ramp_coefficient():
    # On the unit disk about 1 the folded Mellin integral has modulus <= 12.
    h = arb(2) ** -40
    slope = (completed(1 + h, 0).real - completed(1 - h, 0).real) / (2 * h)
    derivative = slope + arb(0, 12 * h**2 / (1 - h**2))
    value = C.log() - arb.const_euler() - derivative / completed(arb(1), 0).real
    if not (arb("-1.017064725656269247184") < value < arb("-1.017064725656269247183")):
        raise ArithmeticError(f"Riesz main coefficient outside its certificate: {value}")
    return value


def verify():
    change, right, margin = top_argument()
    s = acb(2, HEIGHT)
    l_value = right * (-s * C.log()).exp() * s.rgamma()
    if not l_value.real > 0:
        raise ArithmeticError("Ambiguous Euler-product argument at the right endpoint")
    right_argument = HEIGHT * C.log() + s.lgamma().imag + l_value.arg()
    winding = (right_argument - change) / PI
    count = critical_line_count()
    if not (winding.contains(count) and winding > count - HALF and winding < count + HALF):
        raise ArithmeticError(f"Winding enclosure {winding} vs. {count} sign changes")
    print("Riesz main coefficient", ramp_coefficient().str(25))
    print("Mellin tail bound", TAIL.str(6))
    print("Least certified top-edge coordinate margin", margin.str(6))
    print("Winding enclosure", winding.str(65))
    print(f"{count} disjoint line intervals, each of width 1/{DENOMINATOR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locate", action="store_true", help="print regenerated intervals as CSV")
    if parser.parse_args().locate:
        locate()
    else:
        verify()
