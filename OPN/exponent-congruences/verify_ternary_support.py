"""Verify prime-support and abundancy bounds with exact interval arithmetic."""

import argparse
import json
import subprocess
import sys
from fractions import Fraction as F
from functools import cache
from math import ceil, factorial, floor, isqrt, prod
from pathlib import Path

OLDER = Path(__file__).resolve().parents[1] / 'restrained-exponent'
sys.path.insert(0, str(OLDER))
from verify_ordinary_bound import PUBLISHED_PREFIXES, SCALE, log_interval, majorant, require


P = 3_000_000
# (endpoint, residue mod 6) -> (number of candidates, sum of ceiling weights).
PREFIX = {
    (300_000_000, 1): (909_158, 14_719_406_968_152_615),
    (10**10, 5): (26_199_184, 24_422_990_910_195_945),
    (300_000_000_000, 1): (676_766_694, 32_152_938_648_615_242),
    (10**11, 5): (235_925_513, 29_883_045_391_089_345),
}
# Each triple specifies the lower endpoint 10^d, the number of primes fully
# excluded, and the total number of primes in the Bonferroni polynomial.
LOW = ((10, 0, 2), (13, 2, 12), (20, 12, 100), (40, 50, 1227),
       (100, 160, 1227), (400, 800, 1227), (2000, 1227, 1227))
HIGH = ((11, 0, 0), (14, 1, 6), (20, 9, 50), (40, 40, 800),
        (100, 160, 1227), (400, 800, 1227), (2000, 1227, 1227))
PRIMES = [p for p in range(5, 9974) if all(p % d for d in range(2, isqrt(p)+1))]


def pair_sum(values):
    return (sum(values, F(0))**2 - sum(v*v for v in values)) / 2


@cache
def sieve_constants(i, j):
    require(0 <= i <= j <= len(PRIMES) == 1227 and PRIMES[-1] < P, 'invalid sieve cut')
    counts = [3 if p % 3 == 1 else 1 for p in PRIMES[:j]]
    weights = [F(n, p) for n, p in zip(counts, PRIMES[:j])]
    tail = weights[i:]
    density = F(1, 6) * prod(1-w for w in weights[:i]) * (1-sum(tail)+pair_sum(tail))
    labels = counts[i:]
    pairs = (sum(labels)**2-sum(v*v for v in labels)) // 2
    error = prod(1+n for n in counts[:i]) * (1+sum(labels)+pairs)
    require(0 < density <= 1 and error >= 1, 'invalid sieve constants')
    return F(ceil(density*SCALE), SCALE), error


def support_bound(k, y, end, bands, required):
    first = 10**bands[0][0]
    n, weight = PREFIX[y, 1]
    s = F(weight, SCALE) + F(k-n, y-1) + F(1, P-1)
    h = log_interval(y)[1]*s
    partial = s + F(PREFIX[first, 5][1], SCALE)
    contributions = []
    for (a, i, j), b in zip(bands, [row[0] for row in bands[1:]]+[end], strict=True):
        u = 10**a
        require(10**10 <= u and a < b, 'invalid interval')
        rho, error = sieve_constants(i, j)
        part = rho*h*F(u, u-1)*log_interval(F(b, a))[1] + F(11*error*k, 5*(u-1))
        contributions.append(ceil(part*SCALE))
    require(P**7 > 675**16 and end*log_interval(10)[0]/8 > 16+F(k, 674), 'smooth tail')
    require(sum(F(16**i, factorial(i)) for i in range(24)) > 2_000_000, 'exponential tail')
    upper = partial + F(sum(contributions), SCALE) + F(1, 10**6)
    lower = log_interval(required)[0]
    require(upper < lower, 'insufficient separation')
    return dict(cap=k, log_upper_scaled=ceil(SCALE*upper),
                log_target_lower_scaled=floor(SCALE*lower), bands_scaled=contributions)


def factor_counts():
    x = 2*10**11
    count, w = PUBLISHED_PREFIXES[x]
    small = F(w, SCALE) + F(1_280_000_000-count, x)
    require(small < log_interval(F(16, 15))[0], 'universal factor count')
    x, z = 2*10**12, 5*10**17
    density, err = majorant()
    lx = log_interval(x)[0]
    lz, uz = log_interval(z)
    middle = F(x, x-1)*(density*log_interval(uz/lx)[1]
                         + err*(1/lx-1/uz+3/(2*lx**2)+1/(2*lz**2)))
    big = F(PUBLISHED_PREFIXES[x][1], SCALE) + middle + F(4*10**15, z)
    require(big < log_interval(F(64, 51))[0], 'factor count with Euler prime at least 17')
    return dict(universal_log_upper_scaled=ceil(small*SCALE),
                nonfive_log_upper_scaled=ceil(big*SCALE))


def recompute(binary, threads, *, reference=False):
    require(threads >= 1, 'expected positive thread count')
    for (upper, residue), (count, total) in PREFIX.items():
        if reference:
            cmd = [str(binary.resolve()), str(upper), str(residue), str(threads)]
            values = tuple(map(int, subprocess.check_output(cmd, text=True).split()))
            require(values == (upper, residue, count, total), 'reference sieve mismatch')
            continue
        cmd = [str(binary.resolve()), str(P), str(upper), str(threads), str(residue)]
        row = json.loads(subprocess.check_output(cmd, text=True))
        require(row == dict(cutoff=P, upper=upper, residue_mod_6=residue, scale=SCALE,
                            survivors=count, scaled_reciprocal_upper=total), 'prefix mismatch')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sieve', type=Path, help='compiled sieve_ordinary_prefix executable')
    parser.add_argument('--reference-sieve', type=Path, help='compiled sieve_support_classes executable')
    parser.add_argument('--threads', type=int, default=4)
    args = parser.parse_args()
    if args.sieve:
        recompute(args.sieve, args.threads)
    if args.reference_sieve:
        recompute(args.reference_sieve, args.threads, reference=True)
    print(json.dumps([support_bound(10**6, 3*10**8, 6000, LOW, F(16, 15)),
                      support_bound(10**9, 3*10**11, 6_000_000, HIGH, F(64, 51)),
                      factor_counts()], indent=2))
