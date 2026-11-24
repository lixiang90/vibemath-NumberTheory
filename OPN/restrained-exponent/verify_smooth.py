"""Verify forbidden divisors of a+1, then enumerate every 2850-smooth survivor."""

import argparse
import json
from pathlib import Path

from smooth_products import smooth_factorization, surviving_products
from verify_exponent import require, trial_prime, verify, verify_primes
from verify_ordinary_bound import bound_for_prefix


def read(directory, name):
    return json.loads((directory / name).read_text(encoding="utf-8"))


def verify_new_graph(chains, proofs, old):
    roots = {int(d): int(r) for d, r in chains["roots"].items()}
    steps = {int(r): row for r, row in chains["steps"].items()}
    primes = [p for p in range(2, 2851) if trial_prime(p)]
    terminal = {5, 7, 11, 13} | set(map(int, old["steps"]))
    old_indices = {2, *map(int, old["roots"]), *map(int, old["extra_indices"])}
    require(not (steps.keys() & terminal), "a previously excluded prime was redeclared")
    verify_primes(steps, proofs)
    done, visiting = set(), set()

    def exclude(r):
        if r in terminal or r in done:
            return
        require(r in steps and r not in visiting, f"missing or cyclic node {r}")
        visiting.add(r)
        row = steps[r]
        roles = {"square", "euler"} if r % 12 == 1 else {"square"}
        require(r > 13 and row.keys() == roles, f"invalid roles for {r}")
        for role, child in row.items():
            t = int(child)
            value = r*r+r+1 if role == "square" else r+1
            require(t >= 5 and value % t == 0, f"invalid {role} edge for {r}")
            exclude(t)
        visiting.remove(r)
        done.add(r)

    indices = set(old_indices)
    for d in sorted(roots):
        require(d > 1 and not any(d % n == 0 for n in indices), "redundant or invalid index")
        r = roots[d]
        factors = smooth_factorization(d, primes)
        require(r >= 5 and pow(3, d, r) == 1, f"bad power for index {d}")
        require(all(pow(3, d//p, r) != 1 for p in factors), f"order smaller than {d}")
        exclude(r)
        indices.add(d)
    require(done == steps.keys(), "unused exclusion nodes")
    return primes, indices


def replay(directory):
    old = read(directory, "index-chains.json")
    verify(old, read(directory, "primality.json"))
    chains = read(directory, "smooth-chains.json")
    primes, indices = verify_new_graph(chains, read(directory, "smooth-primality.json"), old)
    candidates = sorted(surviving_products(primes, indices))
    bound = 15*10**17
    prefix = read(directory, "ordinary-prefix-2000000000000.json")
    require(prefix["upper"] == 2*10**12, "wrong sieve endpoint")
    arithmetic = bound_for_prefix(prefix, bound, 250*bound)
    large = [(m, p) for m, p in candidates if m > bound+1]
    recorded = read(directory, "smooth-candidates.json")
    require([str(m) for m, _ in large] == recorded, "candidate list mismatch")
    maximum_857 = max(m for m, _ in candidates if m % 857 == 0)
    require(len(candidates) == 1305 and candidates[-1][0] == 2383**27, "wrong smooth range")
    require(maximum_857 == 857*2383**2*2591 and maximum_857 <= bound, "857 branch still open")
    require(len(large) == 173 and min(p for _, p in large) == 2383, "largest prime bound")
    return {"new_forbidden_indices": len(chains["roots"]), "new_exclusion_nodes": len(chains["steps"]),
            "prime_support": [p for p in primes if p not in indices],
            "total_smooth_survivors": len(candidates), "maximum": str(candidates[-1][0]),
            "maximum_multiple_of_857": maximum_857,
            "survivors_after_size_bound": len(large), "least_possible_largest_prime": min(p for _, p in large),
            "candidates_with_largest_prime_below_2383": [str(m) for m, p in large if p < 2383],
            "survivors_above_stormwhip_threshold": sum(m > 114*10**20 + 1 for m, _ in candidates),
            "size_bound": arithmetic}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    print(json.dumps(replay(parser.parse_args().directory), indent=2))
