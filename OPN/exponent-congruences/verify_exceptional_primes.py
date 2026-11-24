"""Classify the sole exceptional even exponent outside 3 using twenty primes."""

from fractions import Fraction
from math import isqrt, prod


PRIMES = (5, 7, 13, 19, 31, 61, 79, 97, 127, 181, 331, 3169, 5233, 5419, 42073,
          1674289, 3348577, 26042690887, 934415109937, 3737657091169)
OUTSIDE = 0
CHAIN = {13, 61, 97, 3169, 3348577, 3737657091169}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def abundance(q, exponent):
    return Fraction(q**(exponent+1)-1, (q-1)*q**exponent)


def closure(exception, euler):
    known = {13, exception, euler} - {OUTSIDE}
    pending = list(known)
    while pending:
        q = pending.pop()
        if q == exception:
            continue
        divisor_sum = q+1 if q == euler else q*q+q+1
        for t in PRIMES:
            if t not in known and divisor_sum % t == 0:
                known.add(t)
                pending.append(t)
    value = abundance(3, 6)
    for q in known:
        exponent = 2
        if q == exception:
            exponent = 4
        elif q == euler:
            exponent = 1
        value *= abundance(q, exponent)
    return value


def verify():
    require(len(PRIMES) == len(set(PRIMES)) == 20, 'incorrect number of vertices')
    for p in PRIMES:
        require(p > 3 and p % 2 == 1 and all(p % d for d in range(3, isqrt(p)+1, 2)),
                f'nonprime vertex {p}')
    euler_choices = (OUTSIDE,) + tuple(p for p in PRIMES if p % 4 == 1)
    special_rows = {5: {13, 61}, 13: {7, 19, 127, 5419}, 61: {31, 331}, 97: {7},
                    3348577: {1674289, 934415109937},
                    181: set(), 3169: set(), 1674289: set(), 3737657091169: set()}
    survivors = 0
    rejected = []
    possible = set()
    for p in euler_choices:
        row = set()
        for q in (OUTSIDE,) + PRIMES:
            if q == p and q != OUTSIDE:
                continue
            value = closure(q, p)
            if value > 2:
                rejected.append(value)
            else:
                row.add(q)
        require(row == special_rows.get(p, CHAIN), f'unexpected row for Euler prime {p}')
        possible.update(row)
        survivors += len(row)
    require((len(rejected), survivors) == (226, 35), 'incorrect case counts')
    require(min(rejected) == Fraction(195540893746, 97706208015), 'minimum excess changed')
    require(len(possible) == 14 and all(q % 6 == 1 for q in possible), 'exception list')
    require(abundance(3, 6)*abundance(5, 2)*abundance(13, 1)*abundance(31, 1) > 2,
            '5 cannot have an even exponent')
    require(abundance(3, 6)*abundance(5, 5)*abundance(13, 2) > 2, 'Euler exponent at 5')
    require(prod(abundance(q, e) for q, e in [(3, 6), (7, 2), (13, 5), (19, 2), (127, 2)]) > 2,
            'Euler exponent at 13')
    print(f'{len(rejected)} abundant and {survivors} surviving assignments; minimum rejected value {min(rejected)}')
    print('possible exceptional primes outside 3:', *sorted(possible))


if __name__ == '__main__':
    verify()
