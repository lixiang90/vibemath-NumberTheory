"""Certify the two Artin factors of the discriminant-229 cubic and its CM twist.

Run with python-flint 0.9.0. Pass --regenerate to reconstruct the interval CSV.
"""

import argparse
import csv
from pathlib import Path

from flint import acb, arb, ctx, dirichlet_char

from verify_cm229 import coefficients, log_derivative_at_one, read_rows
from verify_tetra3136 import (
    HEADER, Completed, FourierRule, GammaFactors, half_winding, locate, verify_signs,
)


FILE = Path(__file__).with_name("cm229-zeros256.csv")
RULE = FourierRule(height=256, samples=8192, modes=2000, right=176, extra=90,
                   denominator=2**14)
PARAMETERS = {"real": (229, 0), "cm_twist": (3664, 1)}


def check_old(new, old):
    low = [row for row in new if row[1] <= 64 * row[2]]
    if len(low) != len(old):
        raise ArithmeticError("old and new counts differ through height 64")
    for (lo, hi, den), (a, b, d) in zip(low, old, strict=True):
        if not (lo*d <= a*den and b*den <= hi*d):
            raise ArithmeticError("old interval is not inside the new interval")


def tail_mass(completed, rows, old):
    cutoff = 64 * RULE.denominator
    merged = list(old) + [row for row in rows if row[0] >= cutoff]
    if len(merged) != len(rows):
        raise ArithmeticError("an isolating interval crosses the height-64 cut")
    heights = [arb(a).union(arb(b))/d for a, b, d in merged]
    low = sum(2/(arb(1)/4+t*t) for t in heights)
    f = completed.factors
    gamma_slope = f.log_scale + sum(((arb(1)+shift)/2).digamma()/2 for shift in f.shifts)
    total = 2*(log_derivative_at_one(completed) + gamma_slope)
    if not total > low:
        raise ArithmeticError("partial zero sum exceeds the total mass")
    return low, total-low


def rational_masses():
    pi, euler = arb.pi(), arb.const_euler()
    q = arb(1)/4
    beta_slope = euler+(2*pi).log()-2*(q.lgamma()-(3*q).lgamma())
    return 2+euler-(4*pi).log(), (4/pi).log()-euler+2*beta_slope


class CompletedBeta:
    def __init__(self):
        self.rule = RULE
        self.factors = GammaFactors(4, (1,))
        self.chi = dirichlet_char(4, 3)

    def top(self, offset):
        s = acb(arb(1)/2+offset, self.rule.height)
        return self.factors.value(s)*s.dirichlet_l(self.chi)


def rational_counts():
    # FLINT certifies individual indexed zeta zeros, including their index.
    index = 1
    while acb.zeta_zero(2*index).imag < RULE.height:
        index *= 2
    lo, hi = index, 2*index
    while hi-lo > 1:
        mid = (lo+hi)//2
        if acb.zeta_zero(mid).imag < RULE.height:
            lo = mid
        else:
            hi = mid
    if not acb.zeta_zero(lo).imag < RULE.height < acb.zeta_zero(hi).imag:
        raise ArithmeticError("height touches an indexed zeta zero interval")
    beta_half, _ = half_winding(CompletedBeta(), mesh=128, degree=32)
    integer = beta_half.unique_fmpz()
    if integer is None or (2*beta_half).unique_fmpz() != 2*integer:
        raise ArithmeticError("beta rectangle winding is ambiguous")
    return lo, int(integer)


def check_boundary_poles():
    rho = acb.zeta_zero(1)
    if not (rho.real == arb(1)/2 and 14 < rho.imag < 15):
        raise ArithmeticError("first zeta zero is outside the asserted interval")
    chi = dirichlet_char(4, 3)
    beta = []
    for t in (arb(6), arb(61)/10):
        s = acb(arb(1)/2, t)
        u = (s+1)/2
        value = (u*(4/arb.pi()).log()).exp()*u.gamma()*s.dirichlet_l(chi)
        if not value.imag.contains(0):
            raise ArithmeticError("completed beta is not real on the critical line")
        beta.append(value.real)
    if not beta[0] > 0 > beta[1]:
        raise ArithmeticError("beta zero not bracketed by opposite signs")
    print("rational zeros for the boundary poles", rho.imag.str(15), "and beta in (6,6.1)")


def check_at_three():
    # At the two primes of norm 3, psi is omega and its conjugate.
    # tau contributes (-1)^n on the ideal of norm 3^n.
    for m, n, expected in ((1, 2, (-2, 1)), (2, 1, (2, -1))):
        actual = []
        for direction in (1, -1):
            counts = [0, 0, 0]
            for a in range(m+1):
                for b in range(n+1):
                    if a+b <= 2 and m+n-a-b <= 2:
                        power = 2*a-m + direction*(2*b-n)
                        counts[power % 3] += (-1)**n
            if counts[1] != counts[2]:
                raise ArithmeticError("coefficient is not invariant under conjugation")
            actual.append(counts[0]-counts[1])
        if tuple(actual) != expected:
            raise ArithmeticError("incorrect local coefficient at 3")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true")
    args = parser.parse_args()
    ctx.prec = 2688
    coeff = coefficients(12288)
    old = read_rows()
    blocks = {} if args.regenerate else read_rows(FILE)
    masses = {}
    with ctx.workprec(192):
        check_boundary_poles()
        check_at_three()
        nz, nd = rational_counts()
        print("rational zeta and beta positive zeros through 256", nz, nd, flush=True)
    for name, (conductor, shift) in PARAMETERS.items():
        print(name, "constructing the completed function", flush=True)
        factors = GammaFactors(conductor, (shift, shift))
        function = Completed(coeff[name], factors, radius=1792, rule=RULE)
        half, margin = half_winding(function, mesh=128, degree=32)
        print(name, "half winding", half.str(28), "top margin", margin.str(8), flush=True)
        if args.regenerate:
            blocks[name] = locate(function)
        rows = blocks[name]
        verify_signs(function, rows)
        if (2*half).unique_fmpz() != 2*len(rows):
            raise ArithmeticError(f"{name}: full winding {2*half} vs. {len(rows)} intervals")
        check_old(rows, old[name])
        masses[name] = tail_mass(function, rows, old[name])
        print(name, "low and remaining mass", *(v.str(9) for v in masses[name]), flush=True)
        print(name, "positive zeros", len(rows), "last", rows[-1], flush=True)

    combined = sorted((lo, hi, name) for name, rows in blocks.items() for lo, hi, _ in rows)
    cross = [b[0]-a[1] for a, b in zip(combined, combined[1:]) if a[2] != b[2]]
    if min(cross) <= 0:
        raise ArithmeticError("overlapping intervals for the two factors")
    print("least separation between different factors", min(cross), "/", RULE.denominator)
    print("field counts with multiplicity", "F", nz+len(blocks["real"]),
          "E", nz+nd+sum(map(len, blocks.values())))
    mz, md = rational_masses()
    ratio = (1+(arb(3)/(4*RULE.height))**2).sqrt()
    endpoints = {"F": (".866731", ".018486"), "E": ("4.70862", ".040372")}
    for key, rational, names in (("F", mz, ("real",)), ("E", mz+md, tuple(PARAMETERS))):
        critical = rational+sum(masses[name][0] for name in names)
        remainder = ratio*sum(masses[name][1] for name in names)
        a, c = map(arb, endpoints[key])
        if not (critical < a and remainder < c):
            raise ArithmeticError(f"{key}: zero masses exceed the reported constants")
        print(key, "absolute formula coefficients", critical.str(10), remainder.str(10))
    upper_b = arb(".87333263744590506611657458485")
    for d, log_max, bound in ((".87333", 10**6, ".92"), (".85", 100, ".84")):
        d = arb(d)
        envelope = arb("4.70862")*arb(1024)**(arb(1)/2-d)
        envelope += arb(".040372")*((upper_b-d)*log_max).exp()
        if not envelope < arb(bound):
            raise ArithmeticError("finite-range improvement fails at its endpoint bounds")
        print("field zero sum relative to target power", d.str(7), envelope.str(10))
    if args.regenerate:
        with FILE.open("w", encoding="utf-8", newline="") as destination:
            writer = csv.writer(destination, lineterminator="\n")
            writer.writerow(("factor", *HEADER))
            for name, rows in blocks.items():
                writer.writerows((name, *row) for row in rows)


if __name__ == "__main__":
    main()
