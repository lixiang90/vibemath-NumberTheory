"""Verify the critical-strip zero count for the Artin factor 3.3136.4t4.a.a.

Requires python-flint 0.9.0. Run with --regenerate to reconstruct the CSV.
All sign and winding decisions use enclosing Arb balls.
"""

import argparse
import csv
from math import factorial, isqrt
from pathlib import Path

from flint import acb, acb_poly, arb, ctx, nmod_poly


HEIGHT = 64
DEN = 2**24
FILE = Path(__file__).with_name("tetra3136-zeros64.csv")
HEADER = ("lower_numerator", "upper_numerator", "denominator")
PERIOD = 256
SAMPLES = 2048
MODES = 480
CUTOFF = 32768
RIGHT = 48
EXTRA = 22
RADIUS = 384
ctx.prec = 896
PI = arb.pi()
HALF = arb(1) / 2
STEP = arb(PERIOD) / SAMPLES
FREQ = 2 * PI / PERIOD
DECAY = arb(3) / 2


def artin_coefficients(limit):
    """Expand the Euler product using the distinct residue degrees of prime ideals."""
    smallest = list(range(limit + 1))
    for p in range(2, isqrt(limit) + 1):
        if smallest[p] == p:
            for m in range(p * p, limit + 1, p):
                if smallest[m] == m:
                    smallest[m] = p

    # The power basis is integral, also at 2 and 7. Factor multiplicities are
    # ramification indices and do not repeat a prime ideal in a zeta factor.
    prime_powers = {}
    for p in range(2, limit + 1):
        if smallest[p] != p:
            continue
        f = nmod_poly([2, 0, 2, -2, 1], p)
        degrees = sorted(g.degree() for g, _ in f.factor()[1])
        if p == 2:
            allowed = ([1],)
        elif p == 7:
            allowed = ([1, 1],)
        elif p % 7 in (1, 6):
            allowed = ([1, 1, 1, 1], [2, 2])
        else:
            allowed = ([1, 3],)
        if degrees not in allowed:
            raise ArithmeticError(f"Unexpected residue degrees at {p}: {degrees}")
        power = p
        max_exp = 0
        while power <= limit:
            max_exp += 1
            power *= p
        # zeta_K / zeta has local numerator 1-X.
        series = [1, -1] + [0] * (max_exp - 1)
        for degree in degrees:
            for e in range(degree, max_exp + 1):
                series[e] += series[e - degree]
        prime_powers[p] = series

    coeff = [0, 1] + [0] * (limit - 1)
    for n in range(2, limit + 1):
        p = smallest[n]
        m, e = n, 0
        while m % p == 0:
            m //= p
            e += 1
        coeff[n] = coeff[m] * prime_powers[p][e]
    return coeff


def gamma_factor(s):
    return ((arb(56).log() - 3 * PI.log() / 2) * s).exp() / PI * (
        (s / 2).gamma() * ((s + 1) / 2).gamma() ** 2
    )


def strip_bound(c, angle):
    """C(c,angle) from contour rotation in the three gamma integrals."""
    c = arb(c)
    return (c.zeta() ** 3 * gamma_factor(c)
            / angle.cos() ** ((3 * c + 2) / 2))


def point_bound(c):
    return 2 * strip_bound(c, arb(1)) / (DECAY * PERIOD)


def boundary_polynomial():
    """Values of the truncated Dirichlet series on RIGHT+i*j*STEP."""
    coeff = artin_coefficients(CUTOFF)
    last = RADIUS * SAMPLES // PERIOD
    values = [acb(1) for _ in range(last + 1)]
    block = 64
    for n, a_n in enumerate(coeff):
        if n < 2 or not a_n:
            continue
        angle = -arb(n).log() * STEP
        unit = acb(0, angle).exp()
        large_step = acb(0, block * angle).exp()
        fine = [acb(1)]
        for _ in range(1, block):
            fine.append(fine[-1] * unit)
        weight = acb(arb(a_n) / (arb(n) ** RIGHT))
        for start in range(0, last + 1, block):
            for k in range(min(block, last + 1 - start)):
                values[start + k] += weight * fine[k]
            weight *= large_step
    return values


def fourier_coefficients():
    boundary = boundary_polynomial()
    raw_tail = arb(2).zeta() ** 3 * arb(CUTOFF) ** (2 - RIGHT)
    box = arb(0, raw_tail)
    bins = [acb(0) for _ in range(SAMPLES)]
    for j, partial in enumerate(boundary):
        s = acb(RIGHT, j * STEP)
        value = gamma_factor(s) * (partial + acb(box, box))
        bins[j % SAMPLES] += value
        if j:
            bins[-j % SAMPLES] += value.conjugate()
    sampled = acb.dft(bins)
    omitted = (2 * strip_bound(RIGHT, arb(1)) * (-DECAY * (RADIUS + STEP)).exp()
               / (SAMPLES * (1 - (-DECAY * STEP).exp())))
    aliases = (2 * point_bound(RIGHT + EXTRA)
               * (-EXTRA * FREQ * (SAMPLES - MODES)).exp()
               / (1 - (-EXTRA * FREQ * SAMPLES).exp()))
    coeff_error = arb(0, omitted + aliases)
    a = RIGHT - HALF
    coefficients = []
    for m in range(MODES + 1):
        coefficients.append((sampled[m].real / SAMPLES + coeff_error)
                            * (-a * FREQ * m).exp())

    rate = (a - 3 * HALF) * FREQ
    high_modes = 2 * point_bound(RIGHT) * (-rate * (MODES + 1)).exp() / (1 - (-rate).exp())
    images = (2 * strip_bound(2, arb(1)) * (-DECAY * (PERIOD - HEIGHT)).exp()
              / (1 - (-DECAY * PERIOD).exp()))
    return coefficients, arb(0, high_modes + images)


class Completed:
    """Fourier enclosures valid for |Re z| <= 64, |Im z| <= 3/2."""
    def __init__(self):
        self.coeff, self.rest = fourier_coefficients()
        self.poly = acb_poly(self.coeff)

    def line(self, numerator):
        phase = acb(2 * arb(numerator) / (DEN * PERIOD)).exp_pi_i()
        return 2 * self.poly(phase).real - self.coeff[0] + self.rest

    def top(self, offset):
        # z=64-i*offset, and exp(i*FREQ*64)=i exactly.
        mag = (FREQ * offset).exp()
        return (self.poly(acb(0, mag)) + self.poly(acb(0, -1 / mag))
                - self.coeff[0] + acb(self.rest, self.rest))

    def grid(self, divisions):
        length = PERIOD * divisions
        data = [acb(0) for _ in range(length)]
        data[0] = acb(self.coeff[0])
        for m in range(1, MODES + 1):
            data[m] = data[-m] = acb(self.coeff[m])
        points = acb.dft(data)
        return [z.real + self.rest for z in points[:HEIGHT * divisions + 1]]


def sign(value):
    if value > 0:
        return 1
    if value < 0:
        return -1
    raise ArithmeticError(f"Unresolved critical-line sign: {value}")


def locate(completed):
    divisions = 64
    signs = [sign(v) for v in completed.grid(divisions)]
    rows = []
    for j in range(1, len(signs)):
        if signs[j - 1] == signs[j]:
            continue
        lo, hi = (j - 1) * DEN // divisions, j * DEN // divisions
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if sign(completed.line(mid)) == signs[j - 1]:
                lo = mid
            else:
                hi = mid
        rows.append((lo, hi, DEN))
    return rows


def verify_signs(completed, rows):
    previous = 0
    for lo, hi, denominator in rows:
        if not (denominator == DEN and previous < lo < hi < HEIGHT * DEN and hi - lo == 1):
            raise ArithmeticError(f"Invalid isolating interval: {(lo, hi, denominator)}")
        if sign(completed.line(lo)) == sign(completed.line(hi)):
            raise ArithmeticError(f"Equal signs at {lo}, {hi}")
        previous = hi


def half_winding(completed):
    mesh, degree = 64, 20
    segments = 3 * mesh // 2
    values = [completed.top(arb(j) / mesh) for j in range(segments + 1)]
    differences = [values]
    for _ in range(degree):
        row = differences[-1]
        differences.append([right - left for left, right in zip(row, row[1:])])

    # Unit disks about [1/2,2]+64i lie between -2 and 3 in real part.
    angle = arb(3) / 2
    disk_bound = strip_bound(3, angle) * (-3 * angle * (HEIGHT - 1) / 2).exp()
    remainder = disk_bound * factorial(degree) / mesh ** (degree + 1)
    change = arb(0)
    margins = []
    for j in range(segments):
        forward = (j + degree <= segments)
        anchor = j if forward else j + 1
        real_margin, imag_margin = abs(values[anchor].real), abs(values[anchor].imag)
        for k in range(1, degree + 1):
            delta = differences[k][anchor if forward else anchor - k]
            real_margin -= abs(delta.real) / k
            imag_margin -= abs(delta.imag) / k
        margin = real_margin.max(imag_margin) - remainder
        if not margin > 0:
            raise ArithmeticError(f"Unresolved top segment {j}: {margin}")
        margins.append(margin.lower())
        change += (values[j + 1] / values[j]).arg()

    s = acb(2, HEIGHT)
    L_right = values[-1] / gamma_factor(s)
    if not (3 * arb(2).zeta().log() < PI / 2 and L_right.real > 0):
        raise ArithmeticError("Euler-product phase is ambiguous")
    gamma_phase = (HEIGHT * (arb(56).log() - 3 * PI.log() / 2)
                   + (s / 2).lgamma().imag + 2 * ((s + 1) / 2).lgamma().imag)
    return (gamma_phase + L_right.arg() - change) / PI, min(margins)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true", help="reconstruct the interval file")
    args = parser.parse_args()
    completed = Completed()
    if args.regenerate:
        rows = locate(completed)
    else:
        with FILE.open(encoding="utf-8", newline="") as handle:
            reader = csv.reader(handle)
            if tuple(next(reader)) != HEADER:
                raise ValueError("Unexpected interval file header")
            rows = [tuple(map(int, row)) for row in reader]
    verify_signs(completed, rows)
    half_count, margin = half_winding(completed)
    if (2 * half_count).unique_fmpz() != 2 * len(rows):
        raise ArithmeticError(f"Full winding {2 * half_count} vs. {len(rows)} positive intervals")
    if args.regenerate:
        with FILE.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, lineterminator="\n")
            writer.writerow(HEADER)
            writer.writerows(rows)
    print("Half of full rectangle count:", half_count.str(28))
    print("Least top-edge coordinate margin:", margin.str(8))
    print(f"{len(rows)} positive simple zeros; first and last dyadic intervals:")
    print(rows[0], rows[-1])


if __name__ == "__main__":
    main()
