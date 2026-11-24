"""Verify the 116-term sieve certificate and the refined restrained exponent bound."""

import argparse
import csv
import hashlib
import json
from fractions import Fraction as F
from functools import cache
from math import ceil, floor, isqrt, prod
from pathlib import Path

from verify_ordinary_bound import PUBLISHED_PREFIXES, SCALE, log_interval, require
from verify_weighted_prefix import B, SUFFIX, X, table_suffix


PRIMES = (7, 13, 19, 31, 37, 43, 61, 67, 73, 79, 97, 103, 109)
BMOR_HASH = "ce5d67887f3e95a4b76b5cfcf9037187c0399357d6848a2f24e65cadf03f443c"
CHAIN_LEVELS = {
    10_000: (213 * 10**20, "3f28abb1bcce622907a563f8f8ca11d0acf75e09fd2bdebe709096cb04370ae3"),
    100_000: (267 * 10**20, "d9547b1b1ce8358973b3f133cdd647cd42e9d8b55f2ddc6082b554d6cc52e9b3"),
    3_000_000: (328 * 10**20, "fa3af87533b095c81964a7c08d9d78217c4812cd1a990c37303b5213bd72f1c2"),
}
DISJOINT = ("backward", "forward_ordinary", "forward_special", "two_specials", "special_conditional")
BASELINE = (5_819_954, 36_514_781_544_466_812)


def check_hash(path, digest):
    data = path.read_bytes().replace(b"\r\n", b"\n")
    require(hashlib.sha256(data).hexdigest() == digest, f"source digest mismatch: {path.name}")


@cache
def trial_prime(n):
    return n == 2 or (n > 2 and n % 2 == 1 and
                      all(n % p for p in range(3, isqrt(n)+1, 2)))


def check_witness(row, level):
    name = row["category"]
    source, target, ordinary, special, other, other_ordinary = (
        int(row[k]) for k in ("r", "ell", "excluded", "special", "other", "other_excluded"))
    require(3_000_000 < source <= 10**9 and trial_prime(source), "invalid source prime")

    def forced_prime(t, h):
        require(t > 3_000_000 and trial_prime(t) and (source**2+source+1) % t == 0,
                "first link of example")
        require(3 < h <= level and trial_prime(h) and (t*t+t+1) % h == 0,
                "ordinary exclusion in example")

    forced_prime(target, ordinary)
    if name in ("two_specials", "special_conditional"):
        require(target % 12 == 1, "invalid forced Euler prime")
        if name == "two_specials":
            forced_prime(other, other_ordinary)
            require(other != target and other % 12 == 1, "two distinct Euler primes needed")
    else:
        require((name == "backward") == (target < source), "incorrect direction of witness")
        if special:
            require(3 < special <= level and trial_prime(special) and (target+1) % special == 0,
                    "Euler exclusion in example")
        else:
            require(target % 4 == 3, "Euler congruence in example")


def chain_totals(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        group = dict(zip(next(reader), next(reader), strict=True))
        bounds = dict(zip(next(reader), next(reader), strict=True))
        records = list(csv.DictReader(handle, delimiter="\t"))
    rows = {row["category"]: (int(row["count"]), int(row["scaled_sum"])) for row in records}
    level = int(bounds["second_limit"])
    require(level in CHAIN_LEVELS, "unsupported second-step cutoff")
    check_hash(path, CHAIN_LEVELS[level][1])
    require((int(bounds["lower"]), int(bounds["upper"])) == (3_000_000, 10**9), "chain interval")
    require(rows["baseline"] == BASELINE and rows["unresolved"] == (0, 0), "chain baseline")
    require(BASELINE[0] == int(bounds["prime_quotients"]) + int(bounds["composite_quotients"]),
            "quotient counts")
    for row in records:
        if row["category"] in DISJOINT:
            check_witness(row, level)
    max_group = int(group["max_group_scaled_sum"])
    require(0 < int(group["special_groups"]) <= rows["special_conditional"][0]
            and 0 < max_group <= rows["special_conditional"][1], "invalid conditional groups")
    count = sum(rows[name][0] for name in DISJOINT)
    weight = sum(rows[name][1] for name in DISJOINT)
    require(0 < count <= BASELINE[0] and max_group < weight <= BASELINE[1], "finite removals")
    return level, count, weight-max_group


def check_bmor(terms, path):
    check_hash(path, BMOR_HASH)
    with path.open(encoding="ascii", newline="") as handle:
        wanted = {3 if d == 1 else d: F(c) for d, _, c in terms}
        found = {int(row["#q"]): F(row["c_pi"]) for row in csv.DictReader(handle, delimiter="\t")
                 if int(row["#q"]) in wanted}
    require(found == wanted, "arithmetic-progression constant mismatch")


def majorant(terms):
    require(len(terms) == 116 and len({d for d, _, _ in terms}) == len(terms),
            "expected distinct sieve monomials")
    values = [0] * (1 << len(PRIMES))
    density, error = F(0), F(0)
    for d, sign, constant in terms:
        subset = [p for p in PRIMES if d % p == 0]
        require(1 <= d <= 100_000 and prod(subset) == d, "invalid squarefree modulus")
        require(sign in (-1, 1) and F(constant) > 0, "invalid coefficient")
        mask = sum(1 << PRIMES.index(p) for p in subset)
        values[mask] = sign
        roots = 2**len(subset)
        density += sign * F(roots, prod(p-1 for p in subset))
        error += (roots if subset else 2) * F(constant)
    # Subset sums evaluate the polynomial at all Boolean incidence vectors.
    for i in range(len(PRIMES)):
        bit = 1 << i
        for mask in range(len(values)):
            if mask & bit:
                values[mask] += values[mask ^ bit]
    require(values[0] == 1 and min(values) >= 0, "invalid pointwise sieve majorant")
    require((density, error) == (F(2_492_150_407, 6_616_209_600), F(2_759_763, 5_000_000)),
            "published density or error constant changed")
    return density, error


def refined_bound(prefix, terms, chain, suffix=SUFFIX):
    require((prefix["cutoff"], prefix["upper"], prefix["scale"]) == (3_000_000, X, SCALE),
            "local sieve range")
    local = prefix["survivors"], prefix["scaled_reciprocal_upper"]
    require(local == PUBLISHED_PREFIXES[X] and suffix == SUFFIX, "finite prefix mismatch")
    level, excluded_count, excluded_weight = chain
    a = CHAIN_LEVELS[level][0]
    count = local[0] + suffix[0] - excluded_count
    weight = local[1] + suffix[1] - excluded_weight
    require(0 < count < 2*a+1 and weight > 0, "invalid finite corrections")
    density, error = majorant(terms)
    z = 290*a
    t = F(B, B-1)
    b, _ = log_interval(B)
    _, u = log_interval(z)
    _, ratio = log_interval(u/b)
    li_lower = z/u - B/b + (z-B)/u**2
    main = density*(t*ratio-li_lower/(z-1))
    remainder = error*(t/b**2+t*t*(1/b-1/u))
    upper = F(weight, SCALE) + F(2*a+1-count, z-1) + main + remainder + F(1, 2_999_999)
    target, _ = log_interval(F(4, 3))
    require(upper < target, "bound fails the abundancy comparison")
    return {"second_step_limit": level, "excluded_a_at_most": a, "split": z, "certificate_terms": len(terms),
            "boolean_inputs": 1 << len(PRIMES), "density": str(density), "error": str(error),
            "finite_count_lower_for_candidate_superset": count, "finite_weight_upper_scaled": weight,
            "upper_scaled_1e18": ceil(upper*SCALE), "target_lower_scaled_1e18": floor(target*SCALE)}


if __name__ == "__main__":
    directory = Path(__file__).resolve().parent
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prefix", type=Path, default=directory / f"ordinary-prefix-{X}.json")
    cli.add_argument("--certificate", type=Path, default=directory / "bmor-majorant.json")
    cli.add_argument("--sieve", type=Path, help="check the original 4943-row sieve.tsv")
    cli.add_argument("--chain", type=Path, default=directory / "chain-3000000.tsv",
                     help="two-step sieve table, including the original author's chain.tsv")
    cli.add_argument("--bmor", type=Path, help="compare each c_pi with the authors' full table")
    args = cli.parse_args()
    prefix = json.loads(args.prefix.read_text(encoding="utf-8"))
    terms = json.loads(args.certificate.read_text(encoding="ascii"))
    if args.bmor:
        check_bmor(terms, args.bmor)
    print(json.dumps(refined_bound(prefix, terms, chain_totals(args.chain),
                                   table_suffix(args.sieve) if args.sieve else SUFFIX), indent=2))
