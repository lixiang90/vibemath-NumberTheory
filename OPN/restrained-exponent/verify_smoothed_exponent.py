"""Certify the smoothed character sieve with integers and rational intervals."""

import argparse
import csv
import json
from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations
from math import ceil, floor, isqrt, prod
from pathlib import Path

from verify_ordinary_bound import PUBLISHED_PREFIXES, SCALE, log_interval, require
from verify_refined_exponent import BASELINE, DISJOINT, chain_totals, check_hash, trial_prime
from verify_weighted_prefix import B, SUFFIX, X, table_suffix


ETA, K, U = F(1, 40), 8, F(57, 2)
DELTA = K*ETA
HEIGHT = 10**8
CONDUCTOR_LIMIT = 400_000
TINY = F(1, 10**12)
FOUR_BILLION = (21_900_925, 44_001_737_926_517_942)
EXTENSION_HASH = "e059582abddaeaa5354e89c6216435615cfe06e9871514e00d2bbbcb3e0dd37f"


def sieve_tree(primes, negative_nodes):
    """Positive nodes omit selected negative children; all positive children are compulsory."""
    allowed = set(negative_nodes)
    require(len(allowed) == len(negative_nodes), "duplicate negative sieve node")
    terms = {}

    def expand(factors, end, sign):
        d = prod(factors)
        require(d <= CONDUCTOR_LIMIT, "sieve exceeds Platt's verified conductor")
        terms[d] = (sign, factors)
        for i, p in enumerate(primes[:end]):
            if sign < 0 or d*p in allowed:
                expand(factors+(p,), i, -sign)

    expand((), len(primes), 1)
    require(allowed == {d for d, (s, _) in terms.items() if s < 0},
            "unreachable or incorrect negative sieve node")
    return terms


def fourier_weights(terms):
    """Combine the characters induced by the same primitive character exactly."""
    coefficients = defaultdict(F)
    norms = {}
    original_norm = F(0)
    for sign, factors in terms.values():
        original_norm += F(4, 3)**len(factors)
        for length in range(len(factors)+1):
            for kept in combinations(factors, length):
                f = prod(kept)
                coefficients[f] += sign*prod(F(2, p-1) for p in factors if p not in kept)
                if f not in norms:
                    norms[f] = prod((F(4, 3)-F(2, p-1) for p in kept), start=F(1))
    weights = {f: abs(v)*norms[f] for f, v in coefficients.items() if v}
    require(0 < coefficients[1] <= 1, "invalid main term")
    return coefficients[1], weights, original_norm


def check_zero_kernel():
    """Monotone interval enclosures below 80 and a power-law comparison beyond."""
    for j in range(80*64):
        left, right = F(j, 64), F(j+1, 64)
        cutoff = (1/(ETA*left))**K if left > 1/ETA else 1
        lower = 3*F(3, 2)/(F(9, 4)+right**2) + 2*U/(U**2+right**2)
        require(cutoff**2 <= lower**2*(F(1, 4)+left**2), "low zero kernel comparison")
    require(80 >= U and F(80, 2**K) <= U and 80*ETA == 2, "high kernel comparison")


@cache
def log_upper(n):
    # Replace the argument by a rational upper endpoint with 11 significant binary bits.
    step = 1 << max(0, n.bit_length()-11)
    return rounded_log(((n+step-1)//step)*step)


@cache
def rounded_log(n):
    return log_interval(n)[1]


@cache
def primitive_error(f):
    require(1 <= f <= CONDUCTOR_LIMIT, "unverified conductor")
    b1 = log_interval_b()-DELTA
    # exp(DELTA) <= 1+DELTA+DELTA^2/(2*(1-DELTA/3)).
    exponential = 1+DELTA+DELTA**2/(2*(1-DELTA/3))
    low_kernel = 3*(log_upper(2*f)/2+F(7, 2))+2*(log_upper(29*f)/2+F(7, 2))
    high_kernel = log_upper(HEIGHT+f)/2+F(7, 2)
    bound = 2/b1*(exponential/isqrt(B)*low_kernel + 3*high_kernel*(F(f, HEIGHT)/ETA)**K)
    return F(ceil((bound+TINY)*SCALE), SCALE)


@cache
def log_interval_b():
    return log_interval(B)[0]


def smooth_mass(z):
    b1 = log_interval_b()-DELTA
    zprime = log_interval(z)[1] + K*ETA**2/2
    require(z > B and b1 > 30, "invalid smoothing interval")
    return log_interval(zprime/b1)[1] - (1-F(B, z))/zprime


def support_estimate(cap, z, finite_count, finite_weight, density, error, tail_power=None):
    """Assume at most cap ordinary primes globally or, if specified, in (B,10^tail_power]."""
    adjustment = finite_count if tail_power is None else 0
    require(cap > adjustment and (tail_power is None or 10**tail_power >= z),
            "invalid support cap or endpoints")
    estimate = F(finite_weight, SCALE) + F(cap-adjustment, z-1) + F(1, 2_999_999)
    estimate += F(B, B-1)*F(z, z-1)*(density*smooth_mass(z)+error)
    if tail_power is not None:
        log_y = tail_power*log_interval(10)[0]
        require(log_y >= 30, "outside the ternary tail theorem")
        estimate += 29/(log_y-7)
    return estimate


@cache
def prime64(n):
    """Deterministic twelve-base test in the range of the extension table."""
    if not 2 <= n < 2**64:
        return False
    bases = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    if any(n % p == 0 for p in bases):
        return n in bases
    d = n-1
    while d % 2 == 0:
        d //= 2
    for p in bases:
        value = pow(p, d, n)
        if value == 1 or value == n-1:
            continue
        for _ in range(((n-1)//d).bit_length()-2):
            value = value*value % n
            if value == n-1:
                break
        else:
            return False
    return True


def extension_totals(path):
    check_hash(path, EXTENSION_HASH)
    with path.open(encoding="ascii", newline="") as stream:
        reader = csv.reader(stream, delimiter="\t")
        groups = dict(zip(next(reader), next(reader), strict=True))
        interval = dict(zip(next(reader), next(reader), strict=True))
        rows = {r["category"]: r for r in csv.DictReader(stream, delimiter="\t")}
    require([int(interval[k]) for k in ("lower", "upper", "second_limit")] ==
            [10**9, 4*10**9, 100_000], "extension interval")
    totals = {k: (int(r["count"]), int(r["scaled_sum"])) for k, r in rows.items()}
    require(totals["baseline"] == tuple(a-b for a, b in zip(FOUR_BILLION, BASELINE))
            and sum(int(interval[k]) for k in ("prime_quotients", "composite_quotients"))
            == totals["baseline"][0], "extension and independent prime sieve disagree")
    require(totals["unresolved"] == (0, 0), "reported unresolved count differs")
    removed = weight = 0
    for k in DISJOINT:
        r = rows[k]
        s, t, h, j, v, u = (int(r[field]) for field in
                             ("r", "ell", "excluded", "special", "other", "other_excluded"))
        require(10**9 < s <= 4*10**9 and prime64(s), "invalid source")
        for target, excluded in [(t, h)]+([(v, u)] if k == "two_specials" else []):
            require(target > 3_000_000 and prime64(target) and (s*s+s+1) % target == 0,
                    "invalid forced prime")
            require(3 < excluded <= 100_000 and prime64(excluded)
                    and (target*target+target+1) % excluded == 0, "invalid ordinary exclusion")
        if k in ("two_specials", "special_conditional"):
            require(t % 12 == 1 and (k != "two_specials" or (v != t and v % 12 == 1)),
                    "invalid forced Euler primes")
        else:
            require((t < s) == (k == "backward"), "wrong direction")
            require((j == 0 and t % 4 == 3) or
                    (3 < j <= 100_000 and prime64(j) and (t+1) % j == 0),
                    "invalid Euler exclusion")
        removed += totals[k][0]
        weight += totals[k][1]
    restored = int(groups["max_group_scaled_sum"])
    require(0 < int(groups["special_groups"]) <= totals["special_conditional"][0] and
            0 < restored <= totals["special_conditional"][1], "invalid conditional correction")
    require(removed <= totals["baseline"][0] and weight <= totals["baseline"][1]
            and weight > restored, "excessive finite removals")
    return removed, weight-restored


def verify(directory, sieve=None, replay_indices=False):
    def read(name):
        return json.loads((directory / name).read_text(encoding="utf-8"))

    record = read("smoothed-sieve.json")
    require(record["prime_limit"] == 800, "incorrect prime list")
    primes = [p for p in range(7, 801, 3) if trial_prime(p)]
    terms = sieve_tree(primes, record["negative_nodes"])
    density, weights, original_norm = fourier_weights(terms)
    require((len(terms), max(terms), len(record["negative_nodes"])) == (1551, 398209, 637),
            "unexpected sieve tree")
    check_zero_kernel()
    error = sum(w*primitive_error(f) for f, w in weights.items()) + original_norm*TINY
    prefix = read(f"ordinary-prefix-{X}.json")
    require((prefix["upper"], prefix["cutoff"], prefix["scale"]) == (X, 3_000_000, SCALE)
            and (prefix["survivors"], prefix["scaled_reciprocal_upper"]) == PUBLISHED_PREFIXES[X],
            "incorrect finite prefix")
    suffix = table_suffix(sieve) if sieve else SUFFIX
    level, removed, weight = chain_totals(directory / "chain-3000000.tsv")
    require(level == 3_000_000, "incorrect first chain table")
    h = prefix["survivors"]+suffix[0]-removed
    v = prefix["scaled_reciprocal_upper"]+suffix[1]-weight
    require((h, v) == (188_812_885_018, 52_797_916_403_925_437), "first finite correction")
    additional, saved = extension_totals(directory / "chain-1b-4b.tsv")
    target = floor(log_interval(F(4, 3))[0]*SCALE)
    bounds = []

    def compare(a, count, total, interval_power=None):
        z, cap = 410*a, 2*a+1
        value = ceil(SCALE*support_estimate(cap, z, count, total, density, error, interval_power))
        require(value < target, "abundancy does not contradict the proposed cap")
        bounds.append({"support_cap_excluded": cap, "z": z, "interval_power": interval_power,
                       "finite_weight_scaled": total, "upper_scaled_1e18": value})

    compare(10**26, h, v)
    h -= additional
    v -= saved
    require((h, v) == (188_798_941_483, 46_312_154_605_862_714), "second finite correction")
    a = 38*10**25
    compare(a, h, v)
    count = 2*a
    log_size = a*log_interval(3)[0] + 2*count*(log_interval(2*count)[0]-1)
    digits_exponent = 408*10**26
    require(log_size > digits_exponent*log_interval(10)[1], "size bound from distinct odd primes")
    for size, power in [(19*10**21, 250), (28*10**24, 1000), (16*10**25, 3000)]:
        compare(size, h, v, power)
    previous = list(map(int, read("shared-candidates.json")["values"]))
    require(previous == sorted([2383**j for j in range(7, 16)]
                               + [2579*2383**j for j in range(6, 11)]), "previous candidate set")
    if replay_indices:
        from verify_shared import replay
        require(replay(directory, sieve)["after_weighted_bound"] == len(previous),
                "previous chain certificates")
    survivors = [m for m in previous if m > a+1]
    require(survivors == sorted([2383**j for j in range(8, 16)]
                                + [2579*2383**j for j in range(7, 11)]), "remaining families")
    return {"terms": len(terms), "largest_modulus": max(terms), "density": str(density),
            "error_upper_scaled_1e18": ceil(SCALE*error), "extension_removed_count": additional,
            "extension_removed_weight_scaled": saved, "finite_superset_count_lower": h,
            "finite_weight_upper_scaled": v, "target_lower_scaled_1e18": target, "bounds": bounds,
            "excluded_exponent_at_most": a, "log10_N_exceeds": digits_exponent,
            "eliminated_smooth_indices": [str(m) for m in previous if m not in survivors],
            "remaining_smooth_indices": len(survivors), "least_remaining_index": str(survivors[0])}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    cli.add_argument("--sieve", type=Path, help="also check all rows of the original sieve.tsv")
    cli.add_argument("--replay-indices", action="store_true", help="replay all earlier index chains")
    args = cli.parse_args()
    print(json.dumps(verify(args.directory, args.sieve, args.replay_indices), indent=2))
