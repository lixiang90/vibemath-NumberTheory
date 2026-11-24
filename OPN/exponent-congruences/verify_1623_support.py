"""Test the Euler prime and the second exception when v_3(N) = 1622."""

import gzip
import heapq
import json
import sys
from collections import defaultdict, deque
from fractions import Fraction as F
from functools import cache
from math import ceil, prod
from pathlib import Path

from flint import fmpz

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'cyclotomic-pairs'))
from search_cyclotomic_pairs import factorizations, primes_through
from verify_exceptional_primes import CHAIN, PRIMES
from verify_ternary_support import log_interval
from verify_two_indices import SCALE, lower, product

X = 3_000_000
K = 1625
SMALL = 997
R = (3**541-1)//2
SPECIAL = {5, 13, 61, 97, 181, 3169, 1674289, 3348577, 3737657091169}
FIXED = ((13, 5), (61, 5), (7, 13), (7, 97), (31, 61), (331, 61),
         (1674289, 3348577), (934415109937, 3348577))


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


@cache
def finite_primes():
    raw = json.loads((HERE/'two-index-primes.json').read_text(encoding='ascii'))
    vertices = tuple(sorted({int(p) for p in raw if int(p) < 10**17}))
    require(R not in vertices and set(PRIMES) <= set(vertices), 'initial vertex set')
    require(all(fmpz(p).is_prime() for p in vertices), 'composite initial vertex')
    return vertices


@cache
def finite_edges(p, euler):
    value = p+1 if euler else p*p+p+1
    return tuple(q for q in finite_primes() if value % q == 0)


def abundance_test(e, pi, seeds=()):
    """With pi=0, record each ordinary pivot used before the first strict excess."""
    known = {13, e, *seeds}
    if pi:
        known.add(pi)
    pending = list(known)
    heapq.heapify(pending)
    ordinary = set()
    value = lower(3, 18)
    while pending:
        p = heapq.heappop(pending)
        if p == e:
            value = product(value, lower(p, 16))
            targets = ()
        elif p == pi:
            value = product(value, lower(p, 1))
            targets = finite_edges(p, True)
        else:
            ordinary.add(p)
            value = product(value, lower(p, 2))
            targets = finite_edges(p, False)
        if value > 2*SCALE:
            return True, ordinary, known
        for q in targets:
            if q not in known:
                known.add(q)
                heapq.heappush(pending, q)
    return False, ordinary, known


def read_extra():
    with gzip.open(HERE/'1623-support-edges.txt.gz', 'rt', encoding='ascii') as src:
        edges = {p: targets for p, *targets in (map(int, line.split()) for line in src)}
    require(len(edges) == 30892 and sum(map(len, edges.values())) == 32071, 'edge table size')
    return edges


class Graph:
    def __init__(self, extra):
        profiles = factorizations(X)
        self.primes = [p for p in profiles if p >= 5]
        self.reverse = defaultdict(list)
        high = set()
        for p in self.primes:
            for q, _ in profiles[p]:
                if q > 3:
                    self.reverse[q].append(p)
                if q > X:
                    high.add(q)
        del profiles
        vertices = set(extra)
        for p, targets in extra.items():
            require(p > X and targets and p != R, 'invalid additional source')
            for q in targets:
                require(q > SMALL and (p*p+p+1) % q == 0, 'invalid additional edge')
                vertices.add(q)
                self.reverse[q].append(p)
        require(all(fmpz(p).is_prime() for p in vertices), 'composite additional vertex')
        high.update(p for p in vertices if p > X)
        small_primes = primes_through(SMALL)
        primorial = fmpz(prod(p for p in small_primes if p % 3 == 1))
        for p in high:
            if p == R:
                continue
            common = primorial.gcd(p*p+p+1)
            for q, _ in common.factor():
                self.reverse[int(q)].append(p)


def fixed_case(graph, e, pi):
    closed, _, core = abundance_test(e, pi)
    require(not closed, 'fixed case closed at its roots')
    bad = {p for p in graph.primes if p <= SMALL and p not in core and abundance_test(e, pi, (p,))[0]}
    queue = list(bad)
    while queue:
        for p in graph.reverse.get(queue.pop(), ()):
            if p not in (e, pi, R) and p not in bad:
                bad.add(p)
                queue.append(p)
    require(not (core & bad), 'fixed case closed by divisibility')
    weights = [0, 0]
    for p in graph.primes:
        if p not in core and p not in bad:
            weights[int(p % 3 == 2)] += weight(p)
    return interval_test(core, weights, pi, bad, False)


def free_euler_case(graph, e):
    possible = {}
    for p in graph.primes:
        if p > SMALL:
            break
        closed, ordinary, _ = abundance_test(e, 0, (p,))
        if closed:
            choices = {r for r in ordinary if r % 4 == 1 and r not in SPECIAL and r != e}
            possible[p] = frozenset(r for r in choices if not abundance_test(e, r, (p,))[0])
    quick = [(q, choices) for q, choices in possible.items() if q < 200]

    @cache
    def can_be_euler(p):
        return (p % 4 == 1 and p not in SPECIAL and p not in (e, R)
                and all((p+1) % q or p in choices for q, choices in quick))

    queue = deque(possible)
    # Missing entries allow every Euler prime. Finite entries can only shrink.
    while queue:
        q = queue.popleft()
        for p in graph.reverse.get(q, ()):
            if p in (e, R):
                continue
            old = possible.get(p)
            if old == frozenset():
                continue
            new = possible[q] | {p} if can_be_euler(p) else possible[q]
            if old is not None:
                new = old & new
            if new != old:
                possible[p] = new
                queue.append(p)

    core = {13, e}
    pending = list(core)
    while pending:
        p = pending.pop()
        if p == e or (p % 4 == 1 and p not in SPECIAL):
            continue
        for q in PRIMES:
            if q not in core and (p*p+p+1) % q == 0:
                core.add(q)
                pending.append(q)
    base = [0, 0]
    extra = [defaultdict(int), defaultdict(int)]
    for p in graph.primes:
        if p in core:
            continue
        i = int(p % 3 == 2)
        if p not in possible:
            base[i] += weight(p)
        else:
            for r in possible[p]:
                # The actual Euler component is already charged to the known product.
                if r != p:
                    extra[i][r] += weight(p)
    weights = [v+max(d.values(), default=0) for v, d in zip(base, extra, strict=True)]
    bad = {p for p in graph.primes if p <= SMALL and possible.get(p) == frozenset()}
    return interval_test(core, weights, 0, bad, True)


def weight(p):
    return (SCALE+p-2)//(p-1)


@cache
def log_ratio(u, v):
    return log_interval(F(v-1, u-1))[1]


def interval_test(core, weights, pi, bad, variable):
    targets = {p for p in core if p % 3 == 1 and p <= 10000}
    count = len(targets) + int(variable)
    w = (F(weights[0]+sum(weight(p) for p in core if p % 3 == 1 and p <= X and p not in targets), SCALE)
         + F(K, X-1))
    v = F(weights[1], SCALE)
    fixed = F(3, 2)*prod(F(p, p-1) for p in core)
    if pi == 5:
        fixed *= F(24, 25)
    if variable:
        fixed *= F(17, 16)
    require(fixed < 2, 'fixed abundance')
    required = log_interval(2/fixed)[0]

    wheels = []
    for r in primes_through(199):
        if r < 5:
            continue
        if r % 3 == 2 or r in targets:
            wheels.append((r, 1))
        elif r in bad:
            wheels.append((r, 3))
    wheels.sort(key=lambda item: F(item[1], item[0]*(2 if item[1] == 3 else 1)), reverse=True)
    density, error = F(1), 1
    prefixes = [(density, error)]
    for r, c in wheels:
        density *= 1-F(c, r)
        error *= 1+c
        prefixes.append((density, error))

    j = 1
    while j**3 < K*K*w:
        j += 1
    d = 1
    while (d+1)**3 <= 1+1/w:
        d += 1
    require(d > 1, 'large-target threshold')
    exponent = ceil(F(j*d, d-1))
    results = []
    for h in (12, 14, 16, 18, 20):
        points = (X, 10**8, 10**12, 10**20, 10**(3*h))
        medium = sum(min(rho*w*log_ratio(u, z)/3 + F(2*err*K, u-1) for rho, err in prefixes)
                     for u, z in zip(points, points[1:]))
        far = 2*F(9, 4)**count*F(11, 4)**exponent/10**h
        results.append(w+v+F(3*(5**count-1), X-1)+medium+far)
    upper = min(results)
    require(upper < required, f'insufficient separation: {float(upper)} against {float(required)}')
    return tuple(ceil(SCALE*f) for f in (w, v, upper)), int(SCALE*required)


def verify():
    graph = Graph(read_extra())
    for e, pi in FIXED:
        print(e, pi, *fixed_case(graph, e, pi), flush=True)
    for e in sorted(CHAIN):
        print(e, 'other', *free_euler_case(graph, e), flush=True)


if __name__ == '__main__':
    verify()
