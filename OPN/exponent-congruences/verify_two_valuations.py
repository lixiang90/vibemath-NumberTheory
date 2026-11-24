"""Check small v_3(N) with at most two exceptional non-Euler exponents."""

import argparse
import sys
from bisect import bisect_right
from collections import Counter, defaultdict
from functools import cache
from heapq import nlargest
from math import isqrt
from pathlib import Path

from flint import fmpz

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "cyclotomic-pairs"))
from search_cyclotomic_pairs import factorizations, primes_through


X = 10**6
Y = 10**13
S = 10**18  # fixed-point denominator; each rounding direction is explicit


def ceil_div(n, d):
    return (n + d - 1) // d


def upper_product(x, y):
    return ceil_div(x * y, S)


def lower_product(x, y):
    return x * y // S


def phi(p):
    return p * p + p + 1


def sigma(p, e):
    return (p**(e + 1) - 1) // (p - 1)


def strip(n, primes):
    for p in primes:
        while n % p == 0:
            n //= p
    return n


def valuation(n, p):
    k = 0
    while n % p == 0:
        n //= p
        k += 1
    return k


@cache
def abundance(p, e):
    n, d = S * sigma(p, e), p**e
    return n // d, ceil_div(n, d)


@cache
def ceiling(p):
    return ceil_div(S * p, p - 1)


@cache
def weight(p):
    return ceil_div(S, p - 1)  # log(p/(p-1)) <= 1/(p-1)


def floor_root(n, k):
    if k == 1:
        return n
    if k == 2:
        return isqrt(n)
    low, high = 1, 1 << ceil_div(n.bit_length(), k)
    while high - low > 1:
        mid = (low + high) // 2
        if mid**k <= n:
            low = mid
        else:
            high = mid
    return low


@cache
def series_factor(p, k):
    scale = 10**6
    root = floor_root(p * scale**k, k)
    return ceil_div(S * root, root - scale)


@cache
def rankin_constant(k):
    scale = 10**6
    root = floor_root(X**(k - 2) * scale**k, k)
    return ceil_div(S * (X + 2) * scale, X * root)


@cache
def rankin_product(targets, k):
    value = rankin_constant(k)
    for t in targets:
        value = upper_product(value, series_factor(t, k))
    return value


@cache
def far_tail(targets, free, cut):
    best = ceil_div(S * 3 * (5**(len(targets) + free) - 1), X)
    for k in (3, 4, 5, 6, 7, 8, 9):
        value = rankin_product(targets, k)
        if free:
            factor = series_factor(max(7, cut + 1), k)
            for _ in range(free):
                value = upper_product(value, factor)
        best = min(best, value)
    return best


@cache
def exp_upper(d):
    # For 0 <= x <= 1, terms starting at degree 6 sum to <= (7/6)x^6/6!.
    if not 0 <= d <= S:
        raise ArithmeticError("exponential bound outside its interval")
    total = term = S
    for k in range(1, 7):
        term = ceil_div(term * d, S * k)
        total += term if k < 6 else ceil_div(7 * term, 6)
    return total


def can_reach_two(high, tail):
    return high >= 2*S or tail >= S or high * exp_upper(tail) >= 2*S*S


class PrimeData:
    def __init__(self, max_cut, max_slots):
        self.profiles = factorizations(X)
        self.small = tuple(p for p in self.profiles if 3 < p <= max_cut)
        self.trial = tuple(p for p in self.profiles if p % 3 == 1 and p <= 30000)
        self.qs = tuple(p for p in self.profiles if p % 6 == 5)
        self.labels = {}
        groups = defaultdict(list)
        totals = Counter()
        for q in self.qs:
            label = self.profiles[q][-1][0]
            self.labels[q] = label
            groups[label].append(q)
            totals[label] += weight(q)
        self.suffix = {}
        for t, qs in groups.items():
            sums = [0] * (len(qs) + 1)
            for i in reversed(range(len(qs))):
                sums[i] = sums[i + 1] + weight(qs[i])
            self.suffix[t] = (qs, sums)
        self.best = {}
        for cut in (3, *self.small):
            if cut in self.labels:
                totals[self.labels[cut]] -= weight(cut)
            self.best[cut] = nlargest(max_slots, ((w, t) for t, w in totals.items() if t > cut))

    @cache
    def factors(self, p, euler):
        if not euler and p in self.profiles:
            return tuple(t for t, _ in self.profiles[p] if t > 3)
        n = p + 1 if euler else phi(p)
        if n.bit_length() > 220:
            return tuple(t for t in self.small if n % t == 0)
        # Partial factorization suffices. Certify primality and divisibility of every edge used.
        factors = []
        for t, _ in fmpz(n).factor_smooth(32, 0):
            t = int(t)
            if t > 3 and fmpz(t).is_prime():
                if n % t:
                    raise ArithmeticError("invalid edge")
                factors.append(t)
        return tuple(factors)

    @cache
    def extended_factors(self, p):
        if p in self.profiles:
            return tuple(t for t, _ in self.profiles[p] if t > 3)
        n = phi(p)
        factors = set(self.factors(p, False))
        if p == 468214602455154100807:
            factors.add(102404822560567927)
        residual = strip(n, (3, *factors))
        if residual > 1:
            if residual.bit_length() <= 112:
                factors.update(int(t) for t, _ in fmpz(residual).factor())
            else:
                factors.update(t for t in self.trial if residual % t == 0)
        if any(n % t or not fmpz(t).is_prime() for t in factors):
            raise ArithmeticError("invalid supplementary edge")
        return tuple(sorted(factors))

    @cache
    def tail(self, cut, targets, free):
        finite = 0
        for t in targets:
            if t in self.suffix:
                qs, sums = self.suffix[t]
                finite += sums[bisect_right(qs, cut)]
        eligible = [w for w, t in self.best[cut] if t not in targets]
        finite += sum(eligible[:free])
        return finite + ceil_div(S * free, cut) + far_tail(targets, free, cut)


@cache
def smooth_preimages(targets, cut, regular_target):
    bound = phi(Y)
    hits = []

    def visit(i, n):
        if i == len(targets):
            root = isqrt(4*n - 3)
            q = (root - 1) // 2
            if root*root == 4*n - 3 and cut < q <= Y and q % 6 == (1 if regular_target else 5):
                hits.append(q)
            return
        while n <= bound:
            visit(i + 1, n)
            n *= targets[i]

    visit(0, 3 if regular_target else 1)
    return tuple(sorted(hits))


def known_support_tail(targets, others, cut):
    finite = sum(weight(q) for q in smooth_preimages(targets, cut, False) if q not in others)
    return finite + ceil_div(S * 3 * (5**len(targets) - 1), Y)


@cache
def target_gap(targets, cut, free):
    if free == 1:
        if any(q not in targets and fmpz(q).is_prime() for q in smooth_preimages(targets, cut, True)):
            return 0
        return Y
    if free != 2:
        return 0
    bound = 10**7
    support = (3, *targets)
    for t in primes_through(bound):
        if t <= cut or t % 3 != 1 or t in targets:
            continue
        residual = strip(phi(t), support)
        if residual == 1:
            return 0
        k = 1
        while t**k <= residual:
            u = floor_root(residual, k)
            if u**k == residual and u >= t and strip(phi(u), (*support, t)) == 1:
                return 0
            k += 1
    return bound


def mixed_tail(targets, free, bound):
    root = isqrt(bound)
    value = S
    for _ in range(free):
        value = upper_product(value, series_factor(bound, 2))
    value -= S
    for t in targets:
        value = upper_product(value, series_factor(t, 2))
    return ceil_div(value * root, root - 2) + ceil_div(S * free, bound)


class ValuationCheck:
    def __init__(self, data, a, cut):
        self.data = data
        self.a = a
        self.cut = cut
        self.small = tuple(p for p in data.small if p <= cut)
        self.budget = 2 - int((a + 1) % 3 != 0)
        self.initial = sigma(3, a)
        self.roots = frozenset(int(p) for p, _ in fmpz(self.initial).factor() if p > 3)
        if any(self.initial % p or not fmpz(p).is_prime() for p in self.roots):
            raise ArithmeticError("invalid initial divisors")
        self.counts = Counter()
        self.unresolved = set()

    @cache
    def finite_tail(self, targets, free, others):
        if not free:
            return known_support_tail(targets, others, self.cut)
        known = set(targets)
        groups = Counter()
        finite = 0
        for q in self.data.qs:
            if q <= self.cut or q in others:
                continue
            new = [t for t, _ in self.data.profiles[q] if t not in known]
            if len(new) > free or any(t <= self.cut for t in new):
                continue
            if new:
                groups[max(new)] += weight(q)
            else:
                finite += weight(q)
        return (finite + sum(nlargest(free, groups.values())) + ceil_div(S * free, self.cut)
                + far_tail(targets, free, self.cut))

    def possible_targets(self, seeds, targets, free):
        extra = set(seeds) - targets
        pending = list(extra)
        while pending:
            t = pending.pop()
            if t <= self.cut or len(extra) > free:
                return False
            for p in self.data.extended_factors(t):
                if p not in targets and p not in extra:
                    extra.add(p)
                    pending.append(p)
        return len(extra) < free or all(strip(phi(t), (3, *targets, *extra)) == 1 for t in extra)

    @cache
    def residue_tail(self, targets, free, others):
        # The Euler prime and both exceptions have been assigned. Every new t is ordinary.
        bound = 100000
        known = set(targets)
        eligible = {t for t in self.data.profiles if self.cut < t <= bound and t % 3 == 1
                    and t not in known and self.possible_targets((t,), known, free)}
        groups = Counter()
        for q in self.data.qs:
            if q <= self.cut or q in others:
                continue
            new = {t for t, _ in self.data.profiles[q]} - known
            if not new or len(new) > free or any(t <= bound and t not in eligible for t in new):
                continue
            if q <= 150000 and not self.possible_targets(new, known, free):
                continue
            groups[max(new)] += weight(q)
        weights = [weight(t) for t in eligible] + [ceil_div(S, bound)] * free
        reciprocal = sum(nlargest(free, weights))
        if 8**17 * (X-1) <= 3**17 * (Y-1):
            raise ArithmeticError("logarithmic range exceeds 17")
        # Two classes modulo 6t for each new t, followed by the smooth-support count.
        return (known_support_tail(targets, others, self.cut) + sum(nlargest(free, groups.values()))
                + ceil_div(20 * reciprocal, 3) + ceil_div(2 * S * free, X-1)
                + ceil_div(S * 3 * (5**(len(targets) + free) - 1), Y))

    def exponent_options(self, p, role, absent, incoming):
        first, step = {"Euler": (1, 4), "ordinary": (2, 6), "exception": (4, 2)}[role]
        options = []
        e = first
        while True:
            if e >= incoming and (role != "exception" or (e + 1) % 3):
                lo, hi = abundance(p, e)
                if (p - 1) * p**e > Y:
                    options.append((lo, ceiling(p), 0))
                    return options
                n = sigma(p, e)
                if all(n % t for t in absent):
                    extra = valuation(n, 3)
                    if role == "ordinary" and p % 3 == 1:
                        extra -= 1
                    if role == "Euler":
                        extra -= valuation(p + 1, 3)
                    options.append((lo, hi, extra))
            e += step

    def terminal(self, regular, exceptions, euler, v, free, targets):
        self.counts["terminal"] += 1
        known = regular | exceptions | ({euler} if euler else set())
        absent = tuple(t for t in self.small if t not in known)
        mandatory = [self.initial, *(phi(q) for q in regular)]
        if euler:
            mandatory.append(euler + 1)
        options = []
        for p in sorted(known):
            incoming = sum(valuation(n, p) for n in mandatory)
            role = "Euler" if p == euler else "ordinary" if p in regular else "exception"
            rows = self.exponent_options(p, role, absent, incoming)
            if not rows:
                return
            options.append(rows)
        others = tuple(sorted(q for q in regular if q > self.cut and q % 3 == 2))
        tail = self.finite_tail(targets, free, others)
        suffix = [S] * (len(options) + 1)
        for i in reversed(range(len(options))):
            suffix[i] = upper_product(suffix[i + 1], max(hi for lo, hi, j in options[i]))

        def visit(i, lo, hi, extra):
            if lo > 2*S or extra + v > self.a or not can_reach_two(upper_product(hi, suffix[i]), tail):
                return False
            if i == len(options):
                return True
            return any(visit(i + 1, lower_product(lo, l), upper_product(hi, h), extra + j)
                       for l, h, j in options[i])

        if not visit(0, *abundance(3, self.a), 0):
            return
        if self.a == 2 and len(exceptions) == self.budget and euler and free in (1, 2):
            bound = target_gap(targets, self.cut, free)
            if bound:
                self.counts[f"gap {free}"] += 1
                tail = known_support_tail(targets, others, self.cut) + mixed_tail(targets, free, bound)
                if not visit(0, *abundance(3, self.a), 0):
                    return
        if self.a == 8 and len(exceptions) == self.budget and euler:
            self.counts["residue"] += 1
            tail = self.residue_tail(targets, free, others)
            if not visit(0, *abundance(3, self.a), 0):
                return
        self.unresolved.add((tuple(sorted(exceptions)), euler))

    def run(self):
        def visit(pending, regular, exceptions, euler, v, ix, lo, hi):
            used = regular | exceptions | ({euler} if euler else set())
            free = self.a - v + self.budget - len(exceptions) + int(not euler)
            if sum(p % 3 == 1 for p in pending) > free:
                return
            minimum = lo
            for p in pending:
                minimum = lower_product(minimum, abundance(p, 1)[0])
            if minimum > 2*S:
                return
            cut = self.small[ix - 1] if ix else 3
            if pending:
                p = min(pending)
                rest = pending - {p}

                def children(is_euler):
                    new = set(self.data.factors(p, is_euler)) - used - pending
                    return None if any(t <= cut for t in new) else rest | new

                inc = int(p % 3 == 1)
                high = upper_product(hi, ceiling(p))
                if v + inc <= self.a:
                    child = children(False)
                    if child is not None:
                        visit(child, regular | {p}, exceptions, euler, v + inc, ix,
                              lower_product(lo, abundance(p, 2)[0]), high)
                if len(exceptions) < self.budget:
                    visit(rest, regular, exceptions | {p}, euler, v, ix,
                          lower_product(lo, abundance(p, 4)[0]), high)
                copies = valuation(p + 1, 3)
                if not euler and p % 4 == 1 and v + copies <= self.a:
                    child = children(True)
                    if child is not None:
                        visit(child, regular, exceptions, p, v + copies, ix,
                              lower_product(lo, abundance(p, 1)[0]), high)
                return
            targets = tuple(sorted(p for p in used if p % 3 == 1))
            if not free and any(strip(phi(q), (3, *targets)) != 1 for q in regular):
                return
            if not can_reach_two(hi, self.data.tail(cut, targets, free)):
                return
            if ix == len(self.small):
                self.terminal(regular, exceptions, euler, v, free, targets)
                return
            q = self.small[ix]
            if q not in used:
                visit({q}, regular, exceptions, euler, v, ix + 1, lo, hi)
            visit(set(), regular, exceptions, euler, v, ix + 1, lo, hi)

        visit(self.roots, set(), set(), 0, 0, 0, *abundance(3, self.a))


CUTOFFS = {2: 997, 4: 499, 6: 499, 8: 997, 12: 997}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("valuations", type=int, nargs="*", metavar="a")
    args = parser.parse_args()
    selected = args.valuations or tuple(CUTOFFS)
    if any(a not in CUTOFFS for a in selected):
        parser.error(f"supported valuations are {tuple(CUTOFFS)}")
    data = PrimeData(max(CUTOFFS[a] for a in selected), max(selected) + 3)
    for a in selected:
        check = ValuationCheck(data, a, CUTOFFS[a])
        check.run()
        print(f"a={a}, B={check.cut}, {dict(check.counts)}, remaining={sorted(check.unresolved)}", flush=True)
        if check.unresolved:
            raise ArithmeticError(f"unresolved valuation {a}")


if __name__ == "__main__":
    main()
