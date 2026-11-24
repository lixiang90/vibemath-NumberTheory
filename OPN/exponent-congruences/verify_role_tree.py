"""Replay the Ochem–Zelinsky divisibility trees with the Euler prime unrestricted mod 3."""

import argparse
import hashlib
import json
from functools import cache
from math import isqrt
from pathlib import Path

from flint import fmpz


LIMIT = 3_000_000
DIGEST = '57b69c549df9fab23ecfda1de66d0794edc1d4a2fbaa33bb882f631351a49eb7'


def require(test, message):
    if not test:
        raise ValueError(message)


def prime_flags(limit):
    flags = bytearray([1]) * (limit + 1)
    flags[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if flags[p]:
            flags[p*p::p] = bytes((limit - p*p) // p + 1)
    return flags


def verify_tree(path):
    data = path.read_bytes().replace(b'\r\n', b'\n')
    require(hashlib.sha256(data).hexdigest() == DIGEST, 'unexpected tree.txt digest')
    flags = prime_flags(LIMIT)

    @cache
    def larger_prime(q):
        return bool(fmpz(q).is_prime())

    def is_prime(q):
        return q >= 2 and (bool(flags[q]) if q <= LIMIT else larger_prime(q))

    blocks = data.decode('ascii').strip().split('\n\n')
    require(blocks[0].startswith('Loaded '), 'missing header')
    roots = [q for q in range(11, LIMIT + 1) if flags[q]]
    branches = 0

    for root, block in zip(roots, blocks[1:], strict=True):
        edges = {}
        nodes = set()
        for line in block.splitlines():
            base, tail = line.strip().split('^')
            q = int(base)
            kind, labels = tail.split(' => ')
            require(q >= root and is_prime(q) and kind in ('2', 'e'), 'invalid node')
            nodes.add(q)
            value = q*q + q + 1 if kind == '2' else q + 1
            for number in labels.split(' F ')[0].split():
                target = int(number)
                require(target >= 2 and value % target == 0, 'incorrect divisor')
                if target >= 5:
                    edges.setdefault((q, kind), set()).add(target)
        require(root in nodes, 'missing root')

        # None covers an Euler prime outside the block. A root 5 mod 12
        # is being tested only as a non-Euler prime. Else every role is open.
        choices = [None] + sorted(q for q in nodes if q % 4 == 1
                                  and not (q == root and q % 12 == 5))
        for euler in choices:
            branches += 1
            seen = {root}
            pending = [root]
            closed = False
            while pending and not closed:
                q = pending.pop()
                kind = 'e' if q == euler else '2'
                for target in edges.get((q, kind), ()):
                    if target < root and flags[target]:
                        if target % 12 != 5 or euler is not None:
                            closed = True
                            break
                    if target in nodes and target not in seen:
                        seen.add(target)
                        pending.append(target)
            require(closed, f'open branch for root {root}, Euler choice {euler}')
    return dict(limit=LIMIT, root_primes=len(roots), euler_choices=branches,
                larger_primes=larger_prime.cache_info().currsize)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tree', type=Path, help='tree.txt from the Ochem–Zelinsky data page')
    print(json.dumps(verify_tree(parser.parse_args().tree), indent=2))
