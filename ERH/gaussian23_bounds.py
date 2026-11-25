"""Check the ball-arithmetic constants and the overlapping level-23 zero tables."""

import csv
from pathlib import Path

from flint import arb, ctx


ctx.prec = 192
HERE = Path(__file__).resolve().parent
T = 1024


def intervals(name):
    with (HERE / name).open(encoding="utf-8-sig", newline="") as handle:
        return [
            (int(row["lower_numerator"]), int(row["upper_numerator"]), int(row["denominator"]))
            for row in csv.DictReader(handle)
        ]


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def exp(x):
    return arb(x).exp()


def main():
    old, new = intervals("zeros23.csv"), intervals("zeros23-1024.csv")
    low = [entry for entry in new if entry[1] <= 64 * entry[2]]
    require(len(old) == len(low) == 59 and len(new) == 1845, "Zero counts disagree")
    for (lo, hi, d), (a, b, q) in zip(old, low):
        require(a * d <= lo * q and hi * q <= b * d, "Old and new zero intervals disagree")
    lo, hi, d = old[0]
    require(lo > 5 * d and old[1][0] > 7 * old[1][2], "First ordinate bounds failed")

    eighth = arb(1) / 8
    half = arb(1) / 2
    quarter = arb(1) / 4
    for v, cutoff in [(arb(1), 30), (quarter, 60)]:
        left_line = exp(v / 8) * (50 * (2 * arb.pi() / v).sqrt() + 8 / v) / (2 * arb.pi())
        require(left_line < cutoff, "Left vertical integral")

    low_tail = 5 * 53 * exp(eighth - arb(49) / 2)
    high_tail = 5 * (4 + T**2) * exp(half - 50)
    first_error = low_tail + high_tail + 30 * exp(-30)
    amplitude_lo = 2 * exp(eighth - (arb(hi) / d)**2 / 2)
    amplitude_hi = 2 * exp(eighth - (arb(lo) / d)**2 / 2)
    require(first_error < arb("6.89e-9"), "One-zero error")

    f_error = exp(-15) + amplitude_hi + first_error
    k_error = 12 * exp(eighth) + 30 * exp(-30)
    powers_two = 15 * exp(quarter)
    powers_three = 24 * exp(arb(1) / 9)
    require(f_error < arb("1e-5") and k_error < 14, "Fixed-variance logarithmic derivatives")
    require(powers_two < 20 and powers_three < 27, "Prime power constants")
    ramified = 3 * arb(23).log() * exp(-15)
    absolute_error = (k_error + 2 * arb("1e-5") + 6 * powers_two + ramified) / 6
    require(absolute_error < arb("21.527"), "Split-prime constant")

    square_error = 400 * exp(-45) + 5 * (4 + T**2) * exp(eighth - arb(T**2) / 8)
    bias_error = (first_error + (1 + arb(23).log()) * exp(-90)
                  + 2 * powers_three * exp(-30) + square_error)
    require(bias_error < arb("7e-9"), "Bias error")
    center = -exp(eighth) / 2
    require(center - amplitude_hi - bias_error > arb("-0.566580"), "Lower bias bound")
    require(center + amplitude_hi + bias_error < arb("-0.566568"), "Upper bias bound")

    variable_tail = 5 * (4 + T**2) * exp(-arb(99) / 2 - eighth)
    variable_l = 5 * 29 * exp(-arb(25) / 2) + variable_tail + exp(-15) + 30 * exp(-30)
    require(variable_l < arb("0.001"), "Variable-variance Artin bound")
    variable_split = ((12 + 30 * exp(-30) + 2 * arb("0.001") + 4 * exp(-15)) / 6
                      + 15 * (half + (2 / arb.pi()).sqrt()))
    threshold = 1 / (T**2 - 3 * quarter)
    require(variable_split < 22, "All-scale split-prime constant")
    require(threshold < arb("9.537e-7"), "Admissible variance slope")

    print("Leading amplitude bounds", amplitude_lo.str(13), amplitude_hi.str(13))
    print("Single-zero approximation error <", first_error.upper().str(8))
    print("Fixed-variance split-prime constant <", absolute_error.upper().str(8))
    print("Secondary term approximation error <", bias_error.upper().str(8))
    print("All-scale split-prime constant <", variable_split.upper().str(8))
    print("Upper bound for the minimum variance slope", threshold.str(12))
    print("log Y endpoints at b=1 and b=31/32", T**2 - 100, arb(16) * (T**2 - 100) / 15)
    b = (777 - arb(1081).sqrt()) / 852
    print("log Y endpoint at Liu's B_2", ((T**2 - 100) / (2 * b - 1)).str(16))
    mixed = arb("0.8733326374459052")
    mixed_end = (T**2 - 100) / (2 * mixed - 1)
    require(mixed_end > 1404211, "Fixed-variance endpoint with mixed moments")
    print("log Y endpoint at the rounded mixed-moment cap", mixed_end.str(16))
    for real_cap, rounded in [(arb(31) / 32, "8.941e-7"), (b, "7.122e-7"),
                              (mixed, "7.121e-7")]:
        slope = (2 * real_cap - 1) / (T**2 - real_cap**2 + quarter)
        require(slope < arb(rounded), "Slope with external zero-free half-plane")
        print("Variance slope at real-part cap", real_cap.str(12), slope.str(12))


if __name__ == "__main__":
    main()
