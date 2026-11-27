"""Certify the complete-splitting prime error of the quartic CM example on [2,10^8]."""

from math import gcd
from operator import gt, lt

from flint import arb, ctx, fmpz, nmod_poly

from verify_cm229_split_primes import LIMIT, checked_extreme, odd_prime_sieve


def run():
    ctx.prec = 128
    primes = odd_prime_sieve(LIMIT)
    residues = [j for j in range(1, 56, 2) if gcd(j, 56) == 1
                and fmpz(8).jacobi(j) == fmpz(28).jacobi(j) == 1]
    theta, split, candidates = arb(0), 0, 0
    maximum = minimum = (-arb(2).sqrt()/8, 2)
    first_above = last_above = crossing = None
    for start in range(0, LIMIT+1, 56):
        for residue in residues:
            p = start+residue
            if p > LIMIT or not primes[p//2]:
                continue
            candidates += 1
            root = fmpz(2).sqrtmod(p)
            positive = (-3-root).jacobi(p) == 1
            x = nmod_poly([0, 1], p)
            poly = nmod_poly([7, 0, 6, 0, 1], p)
            if positive != (x.pow_mod(p, poly) == x):
                raise ArithmeticError(f"two splitting tests disagree at {p}")
            if not positive:
                continue
            split += 1
            bound, rootp = arb(p)/8, arb(p).sqrt()
            minimum = checked_extreme(minimum, (theta-bound)/rootp, p, lt)
            theta += arb(p).log()
            maximum = checked_extreme(maximum, (theta-bound)/rootp, p, gt)
            if theta > bound:
                first_above = first_above or p
                last_above, crossing = p, 8*theta
            elif not theta < bound:
                raise ArithmeticError(f"unresolved sign at {p}")
    end = (theta-arb(LIMIT)/8)/arb(LIMIT).sqrt()
    minimum = checked_extreme(minimum, end, LIMIT, lt)
    if (candidates, split, first_above, last_above, maximum[1], minimum[1]) != (
            1439448, 719724, 47, 28806233, 9791, 91009):
        raise ArithmeticError("counts or extremum locations differ from the published values")
    if not (28806382 < crossing < 28806383 and arb(".37236742") < maximum[0] < arb(".37236743")
            and arb("-2.14203574") < minimum[0] < arb("-2.14203573")
            and arb("-.89838775") < end < arb("-.89838774")):
        raise ArithmeticError("the normalized errors lie outside the reported enclosures")
    print("split candidates and completely split primes", candidates, split)
    print("first and last above-main jumps", first_above, last_above, "last crossing", crossing)
    print("maximum", maximum[0].str(15), "at", maximum[1])
    print("infimum", minimum[0].str(15), "at or just before", minimum[1])
    print("normalized endpoint", end.str(15))


if __name__ == "__main__":
    run()
