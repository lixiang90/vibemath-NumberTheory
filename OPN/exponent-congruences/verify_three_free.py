"""Reproduce the five-exponent exclusion using the Python standard library."""

import sys
from collections import defaultdict
from fractions import Fraction
from itertools import combinations, count, product
from math import isqrt, prod
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cyclotomic-pairs"))
from search_cyclotomic_pairs import factorizations


X = 10**6
Y = 10**16
F = Fraction


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def ceiling_ratio(p):
    return F(p, p - 1)


def sixth_root_lower(p):
    """Round the positive sixth root down to a multiple of 1/1000."""
    low, high = 1000, p * 1000
    while high - low > 1:
        mid = (low + high) // 2
        if mid**6 <= p * 1000**6:
            low = mid
        else:
            high = mid
    return F(low, 1000)


def coarse_tail(targets, outside_count):
    # X^(-2/3) = 1/10000. Each missing exceptional base exceeds X.
    roots = [sixth_root_lower(t) for t in targets]
    series = prod(z / (z - 1) for z in roots)
    return F(outside_count, X) + F(X + 2, X * 10000) * series


def small_groups(profiles):
    groups = defaultdict(list)
    for p, factors in profiles.items():
        if p > 3:
            label = p if p % 3 == 1 else max(r for r, exponent in factors)
            groups[label].append(p)
    return groups


def possible_labels(groups):
    weights = {r: prod(ceiling_ratio(p) for p in ps) for r, ps in groups.items()}
    ordered = sorted(weights, key=weights.get, reverse=True)
    four = prod(weights[r] for r in ordered[:4])
    tail = coarse_tail((7, 13, 19, 31, 37), 5)
    require(tail < F(7, 500), "initial tail bound")
    threshold = 2 * (1 - F(7, 500))
    require(four < threshold, "four labels can reach the threshold")
    pool = sorted(r for r in ordered if weights[r] * four >= threshold)
    eligible = [rs for rs in combinations(pool, 5)
                if prod(weights[r] for r in rs) >= threshold]
    print("Four largest weights, product below 1.909:", ordered[:4], four < F(1909, 1000))
    print("Labels after the four-weight comparison:", pool)
    print("Eligible sets of five labels:", eligible)
    return eligible


def possible_exceptions(label_sets, profiles, groups):
    branches = set()
    for labels in label_sets:
        for size in range(1, len(labels) + 1):
            for targets in combinations(labels, size):
                if not any(t % 12 == 1 for t in targets):
                    continue
                target_set = set(targets)
                regular = [p for r in targets for p in groups[r] if p % 3 == 2
                           and all(t in target_set for t, exponent in profiles[p])]
                known = prod(ceiling_ratio(p) for p in (*targets, *regular))
                options = [[p for p in groups[r] if p % 3 == 2] + [0]
                           for r in labels if r not in target_set]
                # 0 stands for an unknown exceptional prime greater than X.
                for other in product(*options):
                    upper = known * prod(ceiling_ratio(p) for p in other if p)
                    tail = coarse_tail(targets, other.count(0))
                    if upper >= 2 * (1 - tail):
                        branches.add((targets, tuple(sorted(other))))
    return sorted(branches)


def smooth_preimages(targets):
    """List every integer q = 5 mod 6, q <= Y, with Phi_3(q) supported on targets."""
    bound = Y * Y + Y + 1
    hits = []

    def visit(i, value):
        if i == len(targets):
            disc = 4 * value - 3
            root = isqrt(disc)
            q = (root - 1) // 2
            if root * root == disc and q % 6 == 5:
                hits.append(q)
            return
        while value <= bound:
            visit(i + 1, value)
            value *= targets[i]

    visit(0, 1)
    return sorted(hits)


def exponent_intervals(p, role):
    upper = ceiling_ratio(p)
    if role == "regular":
        yield F(1), F(1)  # A smooth preimage need not divide N.
    for e in count(1 if role == "Euler" else 2, 4 if role == "Euler" else 2):
        if ((e + 1) % 3 == 0) != (role == "regular"):
            continue
        gap = F(1, (p - 1) * p**e)
        value = upper - gap
        if gap < F(1, 10**10):
            yield value, upper
            return
        yield value, value


def verify_intervals(targets, other, regular, euler):
    bases = [(p, "Euler" if p == euler else "exception") for p in (*targets, *other) if p]
    bases.extend((p, "regular") for p in regular)
    choices = [list(exponent_intervals(p, role)) for p, role in sorted(bases)]
    # At Y use s=1/4 and 7^(1/4)>8/5. Include any unlocated exceptions.
    tail = F(140, isqrt(Y)) + F(other.count(0), X)
    lower_cut, upper_cut = 2 - F(1, 100000), 2 + F(1, 100000)
    suffix = [F(1)] * (len(choices) + 1)
    for i in reversed(range(len(choices))):
        suffix[i] = suffix[i + 1] * max(hi for lo, hi in choices[i])

    def visit(i, low, high):
        if low > upper_cut or high * suffix[i] < lower_cut * (1 - tail):
            return
        require(i < len(choices), f"unresolved interval for {targets}, {other}, {euler}")
        for lo, hi in choices[i]:
            visit(i + 1, low * lo, high * hi)

    visit(0, F(1), F(1))
    return prod(map(len, choices))


def main():
    require(F(8, 5)**4 < 7 and F(Y + 2, Y) * F(8, 3)**5 < 140, "fine tail bound")
    profiles = factorizations(X)
    groups = small_groups(profiles)
    labels = possible_labels(groups)
    branches = possible_exceptions(labels, profiles, groups)
    require(bool(branches), "unexpected empty intermediate enumeration")
    print("T, other exceptions (0 if above X), regular preimages, interval combinations")
    for targets, other in branches:
        regular = smooth_preimages(targets)
        for euler in targets:
            if euler % 12 == 1:
                boxes = verify_intervals(targets, other, regular, euler)
                print(targets, other, regular, boxes, "Euler prime", euler)
    print(f"Verified {len(branches)} final branches. All exclude 2.")


if __name__ == "__main__":
    main()
