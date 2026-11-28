"""Exact role exclusions for an odd perfect number without 3 and with b=5."""

import sys
from pathlib import Path
from itertools import combinations
from math import prod

sys.path.insert(0, str(Path(__file__).resolve().parent))
from verify_three_free import F, X, ceiling_ratio, coarse_tail, factorizations, require, small_groups
from search_cyclotomic_pairs import primes_through

P = 200
K = 6
BASES = (5, 7, 11)
LABELS = {7, 19, 31}
MIDDLE_CAP = F(503, 500)
# Fixed T, forbidden T, forced exceptions outside T, number of rows,
# and strict caps for U (in thousandths) and D (in millionths).
ROWS = (
    ((7,), (), 2, 988, 1976, 6156),
    ((7, 19, 31), (13,), 1, 376, 1940, 10167),
    ((7, 31), (13,), 1, 376, 1940, 10167),
    ((7, 19), (13, 31), 1, 4688, 1973, 8932),
)


def certificate():
    profiles = factorizations(X)
    groups = small_groups(profiles)
    # The permanent labels come from actual bases, regardless of role.
    expected = {5: 31, 7: 7, 11: 19}
    for p in BASES:
        label = p if p % 3 == 1 else max(r for r, _ in profiles[p])
        require(label == expected[p], "base label")
    require(prod(r**e for r, e in profiles[11]) == 133, "Phi3(11)")
    small = sorted(p for p in profiles if 5 <= p <= P)
    candidates = sorted({p for p in small if p % 3 == 1} |
                        {r for p in small if p % 3 == 2 for r, _ in profiles[p]})
    require(len(candidates) == 35, "number of explicit target primes")
    optional_bases = [p for p in small if p % 3 == 2 and p not in BASES]
    weights = {r: prod(ceiling_ratio(p) for p in ps if p > P)
               for r, ps in groups.items()}
    require(all(w < MIDDLE_CAP for w in weights.values()), "middle label cap")
    small_targets = [p for p in primes_through(P) if p % 3 == 1]

    for fixed, banned, forced, count_expected, u_num, d_num in ROWS:
        require(set(fixed) <= set(candidates) and not set(fixed) & set(banned),
                "target partition")
        slots = K - forced - len(fixed)
        available = [p for p in candidates if p not in fixed and p not in banned]
        bound, cases = F(0), 0
        for j in range(slots + 1):
            for chosen in combinations(available, j):
                targets = set(fixed) | set(chosen)
                if not (targets & {13, 19}):
                    continue
                labels = LABELS | targets
                require(len(labels) <= K, "number of labels")
                mids = prod(weights.get(r, 1) for r in labels)
                mids *= MIDDLE_CAP ** (K - len(labels))
                regular, remaining = [], []
                for p in optional_bases:
                    dst = regular if all(r in targets for r, _ in profiles[p]) else remaining
                    dst.append(p)
                stable = prod(ceiling_ratio(p) for p in
                              [*BASES, *regular, *(p for p in targets if p != 7 and p <= P)])
                for unknown in range(slots - j + 1):
                    if not unknown and not any(p % 12 == 1 for p in targets):
                        continue
                    # Small regular bases cost no extra E slot in this
                    # relaxation; only the remaining small bases need one.
                    extra = remaining[:slots - j - unknown]
                    upper = stable * mids * prod(map(ceiling_ratio, extra))
                    bound = max(bound, upper)
                    cases += 1

        minima = [p for p in small_targets if p not in banned][:K-forced]
        require(len(minima) == K-forced, "enough minimum targets")
        delta = coarse_tail(minima, K)
        upper_cap, tail_cap = F(u_num, 1000), F(d_num, 10**6)
        require(cases == count_expected and 0 < bound < upper_cap, "finite cap")
        require(delta < tail_cap and upper_cap + 2 * tail_cap < 2, "complete cap")
        print(f"T starts with {fixed}, T excludes {banned}, "
              f"{forced} fixed exceptions outside T: {cases} cases, "
              f"U < {upper_cap}, D < {tail_cap}")


if __name__ == "__main__":
    certificate()
