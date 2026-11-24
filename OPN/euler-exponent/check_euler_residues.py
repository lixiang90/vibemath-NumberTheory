"""Compare ternary lifting with the original divisor-sum congruence."""

import json
from itertools import combinations
from math import prod

from euler_residues import euler_residue


def valuation(n):
    h = 0
    while n % 3 == 0:
        h += 1
        n //= 3
    return h


def check():
    cases = comparisons = impossible = 0
    profiles = [[], [(7, 2)], [(19, 4), (11, 6)], [(7, 8), (11, 2)]]
    for p in (5, 13, 17, 37, 41, 61, 73, 89):
        for factors in profiles:
            for a in (0, 2, 4, 6):
                rest = 3**a * prod(r**f for r, f in factors)
                sigma_rest = (3**(a + 1) - 1) // 2
                for r, f in factors:
                    sigma_rest *= sum(r**j for j in range(f + 1))
                h = a - valuation(sigma_rest * (p + 1))
                for k in (1, 2, 3, 4):
                    cases += 1
                    predicted = euler_residue(p, a, factors, k)
                    if h < 0:
                        if predicted is not None:
                            raise ValueError("negative valuation budget accepted")
                        impossible += 1
                        continue
                    modulus = 3**k
                    full_modulus = 3**(a + k)
                    solutions = []
                    for m0 in range(1, modulus):
                        if m0 % 3 == 0:
                            continue
                        # One odd representative in every unit class gives e = 1 mod 4.
                        m = m0 if m0 % 2 else m0 + modulus
                        e = 2 * 3**h * m - 1
                        # Exact division after modular powering evaluates sigma(p**e).
                        v = pow(p, e + 1, (p - 1) * full_modulus)
                        sigma_p = (v - 1) // (p - 1)
                        rhs = 2 * rest * pow(p, e, full_modulus)
                        if (sigma_rest * sigma_p - rhs) % full_modulus == 0:
                            solutions.append(m0)
                        comparisons += 1
                    if len(solutions) != 1 or predicted != {
                        "h": h, "m_residue": solutions[0], "modulus": modulus
                    }:
                        raise ValueError(f"disagreement for {(p, a, factors, k)}")
    factors = [(7, 2), (11, 2), (13, 2)]
    p = 22021
    if p != 19**2 * 61 or (p + 1) * prod(r*r+r+1 for r in [3, 7, 11, 13]) != (
        2 * p * prod(r*r for r in [3, 7, 11, 13])
    ):
        raise ArithmeticError("bad formal control")
    for k in range(1, 13):
        if euler_residue(p, 2, factors, k) != {"h": 0, "m_residue": 1, "modulus": 3**k}:
            raise ArithmeticError("Descartes control failed")
    mod9_cases = 0
    # One prime from each unit class mod 9. Vary their presence independently.
    ordinary = (19, 11, 31, 5, 7, 17)
    for p in (37, 13, 61):
        t = (p - 1) // 3
        for size in range(len(ordinary) + 1):
            for support in combinations(ordinary, size):
                weight = sum({2: 1, 5: 2, 7: 2}.get(r % 9, 0) for r in support)
                for a in (4, 6):
                    h = a - sum(r % 3 == 1 for r in support)
                    wanted = (1 - 3 * (1 + t + weight + (h > 0) * (t + t*t))) % 9
                    actual = euler_residue(p, a, [(r, 2) for r in support], 2)
                    if actual["m_residue"] != wanted:
                        raise ValueError("closed mod-9 formula disagrees with lifting")
                    mod9_cases += 1
    # Cover h=0 with four ordinary primes congruent to 1 modulo 3.
    for p in (37, 13, 61):
        support = (7, 19, 31, 43)
        weight = sum({2: 1, 5: 2, 7: 2}.get(r % 9, 0) for r in support)
        wanted = (1 - 3 * (1 + (p - 1) // 3 + weight)) % 9
        if euler_residue(p, 4, [(r, 2) for r in support], 2)["m_residue"] != wanted:
            raise ValueError("closed mod-9 formula failed for h=0")
        mod9_cases += 1
    return {"profiles_and_precisions": cases, "negative_valuation_budgets": impossible,
            "original_congruences_checked": comparisons, "descartes_precisions": 12,
            "closed_mod9_checks": mod9_cases}


if __name__ == "__main__":
    print(json.dumps(check(), indent=2))
