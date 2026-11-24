"""Rational certificates for residue-class and exponent counts in odd perfect numbers."""

import json
from fractions import Fraction as F
from functools import cache
from math import comb, gcd

from verify_residue_support import BREAK, LOW, SCALE, finite_interval, log_bounds


EPS = F(1, 10**7)
H = 5040
PARAMETERS = {
    "ternary": (6, (5,), 2),
    "quintic": (10, (3, 7, 9), 4),
    "middle": (30, (7, 13, 19), 8),
    "eleven": (30, (11,), 8),
    "both": (30, (17, 23, 29), 8),
}
K3, K5 = F(16, 5), F(81, 20)
D3, D5 = F(19, 2), F(12)


def require(test, reason):
    if not test:
        raise ArithmeticError(reason)


def round_up(value, scale):
    return F(-((-value.numerator*scale)//value.denominator), scale)


def exp_upper(x):
    """Use a positive Taylor polynomial and majorize the remaining ratios."""
    x = F(x)
    require(0 <= x <= 2, "exponential argument outside the certified range")
    term, partial = F(1), F(1)
    for j in range(1, 19):
        term *= x/j
        partial += term
    return round_up(partial + term*x/19/(1-x/20), 10**12)


@cache
def early_sum(modulus, classes):
    scale = 10**12
    by_residue = [0]*modulus
    for j in range(1, H+1):
        by_residue[j % modulus] += (scale+j-1)//j
    largest = max(sum(by_residue[j] for j in range(modulus)
                      if (a+j*r) % modulus in classes)
                  for a in range(modulus) for r in range(modulus) if gcd(r,modulus) == 1)
    return F(largest, scale)


def root_constant(name, log_endpoint):
    m, classes, phi = PARAMETERS[name]
    initial_log, _ = log_bounds(F(H, m))
    require(1 < initial_log < log_endpoint, "incorrect Brun–Titchmarsh cutoff")
    ratio_log = log_bounds(F(log_endpoint)/initial_log)[1]
    tail = F(2*len(classes), phi)*(1+EPS)**2*(F(1, log_endpoint)+ratio_log)
    return round_up(early_sum(m, classes) + tail, 10**9)


def check_large_sieve():
    for c, K in [(F(7, 50)*F(1149, 1000), K3), (F(6, 55)*F(1165, 1000), K5)]:
        base = (1 + F(15, 576) + F(15, 10**28))/(2*c)
        require((base+F(1, 10**4))*(1+EPS)**2 < K, "one-class tail coefficient")
    require(F(70**2, 24*2**35) < F(1, 10**4), "small sieve primes")
    require(70 > 28*log_bounds(10)[1], "range of the sieve denominator")
    require(F(10**8, 10**8-4) < 1+EPS, "exceptional prime loss")
    require((1+F(2, 10**8))*log_bounds(2)[1] < F(7, 10), "total reciprocal weight")
    require(exp_upper(2*(1+EPS)*F(7, 10)) < F(41, 10), "unconditional cubic loss")


def one_exponent(ell, denominator, y):
    ternary = ell == 3
    K, d, classes = (K3, D3, 4) if ternary else (K5, D5, 6)
    C = root_constant("ternary" if ternary else "quintic", y)
    fraction = F(1, denominator)
    initial = (1 + (1 if ternary else 4) + (ell-1)*C)*fraction
    tail = classes*K*exp_upper((ell-1)*(1+EPS)*fraction)/(y-d)
    upper = round_up(initial+tail, 10**12)
    require(upper < log_bounds(2)[0], "exception count contradiction")
    return {"ell": ell, "denominator": denominator, "log_cut": y,
            "upper_in_trillionths": int(upper*10**12)}


def residue_bound(power, denominator, y, z):
    log_h = power*log_bounds(10)[0]
    require(70 <= y <= z and max(8*(y+F(1,100)), 4*(z+F(1,100))) < log_h,
            "all forced residues must lie in the finite interval")
    C = root_constant("ternary", y)
    M = root_constant("middle", y)
    Mz = root_constant("middle", z)
    V = root_constant("quintic", z)
    # A + B/5 <= 1/denominator. B counts exponents missed by both divisibilities.
    weight5 = F(5,denominator) + K3*F(41,10)/(log_h-D3)
    loss5 = exp_upper(4*(1+EPS)*weight5)
    weight3 = max(F(5), 4*(1+Mz))/denominator + 3*K5*loss5/(z-D5)
    loss3 = exp_upper(2*(1+EPS)*weight3)
    coeff_a = 8*C*(1+M) + 8*(1+C) + 4*(1+V)
    coeff_b = 2 + 2*C + 4*(1+V)
    total = max(coeff_a,5*coeff_b)/denominator + 4*K3*loss3/(y-D3) + 6*K5*loss5/(z-D5)
    upper = round_up(total, 10**12)
    require(upper < log_bounds(2)[0], "residue count contradiction")
    return {"decimal_cutoff_power": power, "denominator": denominator, "y": y, "z": z,
            "upper_in_trillionths": int(upper*10**12)}


def mixed_complement(prefix):
    """Bound the opposite three target classes via the smaller factor of Phi_5(q)."""
    size, cutoff, y = 760, 10**47, 20_000
    count, numerator = prefix
    require((count, numerator) == (135_759_787, 800_384_533_846_518),
            "three-class finite interval mismatch")
    b, _ = log_bounds(BREAK)
    x0, x1 = log_bounds(cutoff)
    error = F(3,840)*(1/b**2 + 1/x0**2 + 1/b - 1/x1)
    small = F(numerator,SCALE) + F(BREAK,BREAK-1)*(F(3,8)*log_bounds(x1/b)[1]+error)
    C, U, V = [root_constant(name,y) for name in ("ternary","eleven","both")]
    loss3 = exp_upper(2*(1+EPS)*F(size,LOW))
    small5 = F(1,LOW) + 2*(1+U)*F(size,LOW)
    loss5 = exp_upper(4*(1+EPS)*(small5+K3*loss3/(y-D3)))
    multiset_weight = sum(F(4*comb(size+j-1,j), max(cutoff,100**j-3))
                          for j in range(1,100))
    # Above j=99, replace the denominator by 100**j/2 and use a geometric tail.
    ratio = F(size+100, 101*100)
    multiset_weight += F(8*comb(size+99,100),100**100)/(1-ratio)
    upper = round_up(small + multiset_weight + F(size+6,LOW) + 2*(1+C)*F(size,LOW)
                     + 4*V*small5 + 4*K3*loss3/(y-D3) + 6*K5*loss5/(y-D5), 10**12)
    require(upper < log_bounds(2)[0], "three-class support comparison")
    return {"excluded_total_support": size, "decimal_cutoff_power": 47, "log_cut": y,
            "small_count": count, "finite_weight_numerator": numerator,
            "upper_in_trillionths": int(upper*10**12)}


def verify(prefix):
    check_large_sieve()
    return {"exponent_counts": [one_exponent(3,22,170), one_exponent(5,70,320)],
            "residue_counts": [residue_bound(400,700,115,230), residue_bound(3500,580,130,2000)],
            "complementary_classes": mixed_complement(prefix)}


if __name__ == "__main__":
    prefix = finite_interval(LOW, BREAK, extra_excluded=(4,7,11,13))
    print(json.dumps(verify(prefix), indent=2))
