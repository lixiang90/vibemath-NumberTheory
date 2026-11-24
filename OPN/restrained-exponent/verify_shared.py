"""Replay the three chain certificates and independently re-enumerate the smooth indices."""

import argparse
import json
from pathlib import Path

from smooth_products import surviving_products
from verify_ordinary_bound import require
from verify_smooth import read, replay as replay_previous, verify_new_graph
from verify_weighted_prefix import A, SUFFIX, table_suffix, weighted_bound


def replay(directory, table=None):
    previous = replay_previous(directory)
    base = read(directory, "index-chains.json")
    middle = read(directory, "smooth-chains.json")
    base["steps"].update(middle["steps"])
    base["roots"].update(middle["roots"])
    graph = read(directory, "shared-chains.json")
    proofs = read(directory, "shared-primality.json")
    primes, forbidden = verify_new_graph(graph, proofs, base)
    fresh = {int(d) for d in graph["roots"]}
    candidates = sorted(m for m, _ in surviving_products(primes, forbidden) if m > 15*10**17+1)
    old_list = list(map(int, read(directory, "smooth-candidates.json")))
    require(previous["survivors_after_size_bound"] == len(old_list), "wrong initial candidate count")
    require(candidates == [m for m in old_list if all(m % d for d in fresh)],
            "enumeration and filtering disagree")
    weighted = weighted_bound(read(directory, "ordinary-prefix-2000000000000.json"),
                              table_suffix(table) if table is not None else SUFFIX)
    final = [m for m in candidates if m > A+1]
    record = read(directory, "shared-candidates.json")
    require([str(m) for m in final] == record["values"], "final list mismatch")
    expected = sorted(c*2383**j for c, low, high in record["families"] for j in range(low, high+1))
    require(final == expected, "factorized families do not give the exact list")
    return {"fresh_forbidden_indices": len(fresh), "new_exclusion_nodes": len(graph["steps"]),
            "elliptic_steps": len(proofs["elliptic"]), "lucas_certificates": len(proofs["lucas"]),
            "previous_candidates": len(old_list), "after_new_chains": len(candidates),
            "after_weighted_bound": len(final), "families": record["families"],
            "least": str(min(final)), "greatest": str(max(final)), "weighted_estimate": weighted}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("directory", nargs="?", type=Path, default=Path(__file__).resolve().parent)
    cli.add_argument("--table", type=Path, help="source sieve.tsv; verify its interval totals too")
    args = cli.parse_args()
    print(json.dumps(replay(args.directory, args.table), indent=2))
