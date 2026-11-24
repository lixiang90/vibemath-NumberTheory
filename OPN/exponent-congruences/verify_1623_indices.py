"""Reconstruct the exponent-index obstructions at the two forced repunit factors."""

import heapq
import json
from functools import cache
from math import prod
from pathlib import Path

from flint import fmpz

from verify_exceptional_indices import BOUNDS
from verify_exceptional_primes import PRIMES, abundance
from verify_two_indices import SCALE, SECOND_EXCEPTIONS, lower, product, small_primes

HERE = Path(__file__).resolve().parent
X = 3**541
R = (X-1)//2
L = (X*X+X+1)//13


def require(value, message):
    if not value:
        raise ArithmeticError(message)


def divides_sum(p, base, length):
    if base % p == 1:
        return length % p == 0
    return pow(base, length, p) == 1


def verify():
    require(fmpz(R).is_prime() and fmpz(L).is_prime() and 13*R*L == (3**1623-1)//2,
            'initial factorization')
    require(R % 3 == L % 3 == 1 and R != L, 'repunit residues')
    former = json.loads((HERE/'exceptional-index-primes.json').read_text(encoding='ascii'))
    require(R not in set(map(int, former['primes'])) | set(PRIMES), 'exception inside a transferred proof')
    require(set(BOUNDS) == SECOND_EXCEPTIONS, 'second exception table')
    vertices = tuple(map(int, json.loads((HERE/'1623-index-primes.json').read_text(encoding='ascii'))))
    require(vertices == tuple(sorted(set(vertices))) and min(vertices) >= 5, 'vertex order')
    require(R not in vertices and L not in vertices and set(BOUNDS) <= set(vertices), 'unexpected vertex')
    require(all(fmpz(p).is_prime() for p in vertices), 'composite index vertex')

    @cache
    def adjacent(p, euler):
        value = p+1 if euler else p*p+p+1
        return tuple(q for q in vertices if value % q == 0)

    def close(e, pi, seeds):
        known = set(seeds) | {13, e}
        if pi:
            known.add(pi)
        pending = list(known)
        heapq.heapify(pending)
        value = lower(3, 18)
        ordinary = []
        while pending:
            p = heapq.heappop(pending)
            if p == e:
                value = product(value, lower(p, 4))
                targets = ()
            elif p == pi:
                value = product(value, lower(p, 1))
                targets = adjacent(p, True)
            else:
                ordinary.append(p)
                value = product(value, lower(p, 2))
                targets = adjacent(p, False)
            if value > 2*SCALE:
                return ordinary
            for q in targets:
                if q not in known:
                    known.add(q)
                    heapq.heappush(pending, q)
        raise ArithmeticError(f'unclosed exponent index for E={e}, Euler={pi}, roots={sorted(seeds)}')

    euler_roots = [p for p in vertices if (L+1) % p == 0]
    for e in sorted(BOUNDS):
        close(e, 0, euler_roots)
    print('L excluded as Euler prime for all 14 second exceptions.')

    for name, base, bound, factor in [('R', R, 41, 1), ('L', L, 59, 3)]:
        count = 0
        indices = small_primes(bound-1)
        for prime in indices:
            d = factor*prime
            seeds = [p for p in vertices if divides_sum(p, base, d)]
            require(bool(seeds), f'no roots for {name}, {prime}')
            for e in sorted(BOUNDS):
                ordinary = close(e, 0, seeds)
                count += 1
                for pi in ordinary:
                    if pi % 4 == 1:
                        close(e, pi, seeds)
                        count += 1
        print(f'{name}: {len(indices)} prime indices, {count} Euler assignments.')
    require(fmpz(7822361).is_prime() and R % 7822361 != 1 and pow(R, 7, 7822361) == 1,
            'seventh-order witness')
    surviving = (19, 127, 5419)
    require(min(p**(BOUNDS[p]-1) for p in surviving) == 19**262, 'small exceptional power')
    for terms in ([(3, 6), (7, 8), (13, 1), (19, 2), (127, 2)],
                  [(3, 18), (7, 8), (13, 5), (19, 16)]):
        require(prod(abundance(p, j) for p, j in terms) > 2, 'small-exponent alternatives')
    lower_size = 3**1622 * 7**2 * 13 * 19**262 * R**40 * L**2
    require(10**12454 < lower_size < 10**12455, 'size at the exact index')
    require(1623*857-1 == 1390910 and 3**1390910 > 10**12454, 'larger 3-adic index')
    print('N > 10^12454; both small-exponent alternatives checked.')


if __name__ == '__main__':
    verify()
