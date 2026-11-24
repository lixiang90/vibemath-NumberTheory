"""Recompute the weighted exponent bound, optionally reading the original TSV."""

import argparse
import csv
import hashlib
import json
from fractions import Fraction
from math import ceil, floor
from pathlib import Path

from verify_ordinary_bound import PUBLISHED_PREFIXES, SCALE, log_interval, majorant, require


X = 2*10**12
B = 50_290_000_000_000
A = 261*10**18
ROWS = 4829
SUFFIX = (180_426_232_393, 12_351_431_182_971_148)
TABLE_SHA256 = "68d2f87fabcc38920d891040829f3f42d72e387b50066d71ede230a483d9fb8f"


def table_suffix(path):
    digest = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    require(digest == TABLE_SHA256, "source table digest mismatch")
    endpoint, count, weight, rows = 3_000_000, 0, 0, 0
    with path.open(encoding="utf-8", newline="") as file:
        for row in csv.DictReader(file, delimiter="\t"):
            lo, hi, limit, n, total = (int(row[k]) for k in
                ("lower_exclusive", "upper_inclusive", "residue_limit", "survivors",
                 "scaled_reciprocal_upper"))
            require(lo == endpoint and lo < hi and not lo < X < hi, "noncontiguous table")
            require(0 < limit <= 3_000_000 and n >= 0 and total >= n*((SCALE+hi-2)//(hi-1)),
                    "invalid sieve row")
            endpoint = hi
            if lo >= X:
                rows += 1
                count += n
                weight += total
    require(endpoint == B and rows == ROWS and (count, weight) == SUFFIX, "suffix mismatch")
    return count, weight


def weighted_bound(prefix, suffix=SUFFIX, a=A):
    require(prefix["upper"] == X and prefix["cutoff"] == 3_000_000 and prefix["scale"] == SCALE,
            "incorrect local sieve interval")
    local = prefix["survivors"], prefix["scaled_reciprocal_upper"]
    require(local == PUBLISHED_PREFIXES[X] and suffix == SUFFIX, "sieve totals mismatch")
    count, total = (u+v for u, v in zip(local, suffix))
    z = 262*a
    require(z > B and 2*a+1 >= count, "invalid weighted interval")
    density, error = majorant()
    t = Fraction(B, B-1)
    lb, _ = log_interval(B)
    _, uz = log_interval(z)
    _, ratio = log_interval(uz/lb)
    li_lower = z/uz-B/lb+(z-B)/uz**2
    main = density*(t*ratio-li_lower/(z-1))
    remainder = error*(t/lb**2+t*t*(1/lb-1/uz))
    upper = Fraction(total, SCALE) + Fraction(2*a+1-count, z-1) + main + remainder
    upper += Fraction(1, 2_999_999)
    target, _ = log_interval(Fraction(4, 3))
    require(upper < target, "weighted estimate does not imply the assumed exponent bound")
    return {"excluded_a_at_most": a, "local_endpoint": X, "table_endpoint": B, "split": z,
            "combined_count": count, "combined_weight": total,
            "upper_scaled_1e18": ceil(upper*SCALE), "target_lower_scaled_1e18": floor(target*SCALE)}


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--prefix", type=Path, default=Path(__file__).resolve().with_name(
        "ordinary-prefix-2000000000000.json"))
    cli.add_argument("--table", type=Path, help="also verify all rows of the published sieve.tsv")
    args = cli.parse_args()
    prefix = json.loads(args.prefix.read_text(encoding="utf-8"))
    print(json.dumps(weighted_bound(prefix, table_suffix(args.table) if args.table else SUFFIX), indent=2))
