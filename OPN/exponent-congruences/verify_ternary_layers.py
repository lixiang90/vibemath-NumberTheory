"""Replay the two-step count bound and apply the prime-power layer inequality."""

import argparse
import json
import sys
from fractions import Fraction
from math import ceil
from pathlib import Path


OLDER = Path(__file__).resolve().parent.parent / "restrained-exponent"
sys.path.insert(0, str(OLDER))
from verify_ordinary_bound import SCALE, require
from verify_refined_exponent import chain_totals, check_bmor, refined_bound
from verify_weighted_prefix import SUFFIX, X, table_suffix


def verify(sieve=None, bmor=None):
    prefix = json.loads((OLDER / f"ordinary-prefix-{X}.json").read_text(encoding="utf-8"))
    terms = json.loads((OLDER / "bmor-majorant.json").read_text(encoding="ascii"))
    if bmor:
        check_bmor(terms, bmor)
    chain = chain_totals(OLDER / "chain-3000000.tsv")
    base = refined_bound(prefix, terms, chain, table_suffix(sieve) if sieve else SUFFIX)
    a = base["excluded_a_at_most"]
    limit = 65_625 * 10**18
    require(a == 328*10**20 and base["split"] == 290*a, "incorrect reference row")
    require(limit-2*a == 25*10**18, "count correction")
    upper = Fraction(base["upper_scaled_1e18"], SCALE)+Fraction(limit-2*a, base["split"]-1)
    target = Fraction(base["target_lower_scaled_1e18"], SCALE)
    require(upper < target, "count bound does not close")
    exponent_bound, rem = divmod(2*limit, 7)
    require(rem == 0 and exponent_bound == 1875*10**19, "exponents 2 and 8")

    record = json.loads((OLDER / "shared-candidates.json").read_text(encoding="utf-8"))
    values = sorted([2383**j for j in range(7, 16)]+[2579*2383**j for j in range(6, 11)])
    require(list(map(int, record["values"])) == values and
            261*10**18 < exponent_bound < min(values)-1, "smooth index range")
    print(json.dumps({"excluded_M_at_most": limit+1, "upper_scaled_1e18": ceil(upper*SCALE),
                      "target_lower_scaled_1e18": base["target_lower_scaled_1e18"],
                      "excluded_v3_at_most": exponent_bound, "smooth_candidates": len(values)},
                     indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sieve", type=Path, help="Stormwhip's original sieve.tsv")
    parser.add_argument("--bmor", type=Path, help="BMOR c_all_rounded.txt")
    args = parser.parse_args()
    verify(args.sieve, args.bmor)
