"""Check the cyclotomic traces in the common-translate note (python-flint 0.9)."""

import math
from functools import cache

from flint import acb, arb, ctx, dirichlet_char, fmpz


@cache
def units(r):
    return tuple(a for a in range(1, r + 1) if math.gcd(a, r) == 1)


def divisors(n):
    return (d for d in range(1, n + 1) if n % d == 0)


@cache
def root(r, exponent):
    return acb(arb(2 * (exponent % r)) / r).exp_pi_i()


def check_zero(value):
    assert value.contains(0), value
    radius = float(abs(value).upper())
    assert radius < 1e-40, radius
    return radius


def prime_powers(y):
    for p in range(2, y + 1):
        if fmpz(p).is_prime():
            n = p
            while n <= y:
                yield n, arb(p).log()
                n *= p


def field_degree(r, chi):
    """Count the image of Gal(Q(zeta_lcm)/Q) in the two cyclotomic factors."""
    two_part = r & -r
    order = int(chi.order())
    image = {int(chi.chi_exponent(a)) for a in units(r)
             if (a - 1) % two_part == 0 or (a + 1) % two_part == 0}
    c = len(image)
    assert order % c == 0
    e = max(2, two_part.bit_length() - 1)
    two_level = 2 ** (e + 2)
    pairs = {(min(a % two_level, (-a) % two_level), int(chi.chi_exponent(a % r)))
             for a in units(math.lcm(two_level, r))}
    assert len(pairs) == 2**e * c
    return order // c


@cache
def trace_blocks(r):
    """Form each generator projection from relative traces of the root of unity."""
    group = units(r)
    g = len(group)
    blocks = {}
    for label in group:
        chi = dirichlet_char(r, label)
        if not chi.is_primitive() or r == 1:
            continue
        o = int(chi.order())
        generators = frozenset(pow(label, k, r) for k in units(o))
        if generators in blocks:
            continue
        ex = int(chi.group().exponent())
        code = {a: int(chi.chi_exponent(a)) for a in group}
        coordinates = {a: acb(0) for a in group}
        for d in divisors(o):
            mu = int(fmpz(o // d).moebius_mu())
            if not mu:
                continue
            kernel = [h for h in group if (code[h] * (o // d)) % ex == 0]
            for a in group:
                coordinates[a] += mu * sum(root(r, a * h % r) for h in kernel) / len(kernel)
        # The normalized trace on the embedding space contributes 1/g.
        # The two vector scalings contribute g/r.
        mixed = {t: sum(coordinates[a].conjugate() * coordinates[a * t % r]
                        for a in group) / r for t in group}
        for t in group:
            order_t = ex // math.gcd(ex, code[t])
            exact = arb(int(fmpz(order_t).moebius_mu()) * len(generators)) / (g * len(units(order_t)))
            check_zero(mixed[t] - exact)
            fourier = sum(dirichlet_char(r, h)(t) for h in generators) / g
            check_zero(mixed[t] - fourier)
        blocks[generators] = (field_degree(r, chi), mixed)
    return blocks


def check_counts(q, endpoints):
    group = units(q)
    blocks = [(r, mixed) for r in divisors(q) if r > 1
              for _, mixed in trace_blocks(r).values()]
    worst = 0.0
    powers = list(prime_powers(max(endpoints)))
    for y in endpoints:
        counts = {a: arb(0) for a in group}
        for n, weight in powers:
            if n <= y and math.gcd(n, q) == 1:
                counts[n % q] += weight
        mean = sum(counts.values()) / len(group)
        for b in group:
            measured = acb(0)
            for r, mixed in blocks:
                inv = pow(b, -1, r)
                measured += len(units(r)) * sum(counts[a] * mixed[a * inv % r] for a in group)
            difference = measured / len(group) - (counts[b] - mean)
            worst = max(worst, check_zero(difference))
    return len(blocks), worst


def check_level_obstruction():
    blocks = trace_blocks(17)
    generators, (_, kernel) = next((s, v) for s, v in blocks.items() if len(s) == 8)
    group = units(17)
    for t in group:
        check_zero(kernel[t] - arb(int(t == 1) - int(t == 16)) / 2)
    for label in group:
        chi = dirichlet_char(17, label)
        coefficient = sum(chi(t).conjugate() * kernel[t] for t in group) / 16
        energy = sum(abs(kernel[t] - coefficient * chi(t))**2 for t in group) / 16
        check_zero(energy - arb(8 - int(label in generators)) / 256)


def check_split_weil(r, label, w):
    """Use the trace fibers of F_(w^2) and a cyclic subgroup of F_w^*."""
    g = len(units(r))
    chi = dirichlet_char(r, label)
    o = int(chi.order())
    assert fmpz(w).is_prime() and (w - 1) % o == 0 and pow(-131, (w - 1) // 2, w) == 1
    alpha = next(a for a in range(2, w) if pow(a, o, w) == 1
                 and all(pow(a, d, w) != 1 for d in divisors(o) if d < o))
    subgroup = [pow(alpha, j, w) for j in range(o)]
    remaining = set(range(1, w))
    representatives = set()
    while remaining:
        t = min(remaining)
        representatives.add(t)
        remaining.difference_update(h * t % w for h in subgroup)
    # In a quadratic power basis with trace-zero generator, Tr(s+t*sqrt(d))=2s.
    sizes = [0] * w
    for s in range(w):
        for _ in range(w):
            sizes[2 * s % w] += 1
    assert sizes == [w] * w
    cells = [{h * t % w for t in representatives} for h in subgroup]
    assert sum(len(cell) for cell in cells) == len(set.union(*cells)) == w - 1
    assert all(sum(sizes[t] for t in cell) * o == w * (w - 1) for cell in cells)
    # The projected vector has values sqrt(w / (g*(w-1))) * coeff[t].
    coeff = [0] * w
    for j, cell in enumerate(cells):
        h_order = o // math.gcd(o, j)
        value = int(fmpz(h_order).moebius_mu()) * len(units(o)) // len(units(h_order))
        for t in cell:
            coeff[t] = value
    ex = int(chi.group().exponent())
    gen = frozenset(pow(label, i, r) for i in units(o))
    kernel = trace_blocks(r)[gen][1]
    for a in units(r):
        j, remainder = divmod(o * int(chi.chi_exponent(a)), ex)
        assert remainder == 0
        local_inner = sum(coeff[t] * coeff[subgroup[j] * t % w] for t in range(w))
        check_zero(kernel[a] - arb(local_inner) / (g * (w - 1)))
    return o, len(representatives)


def check_fourth_moment(q, cutoff, lo, hi):
    """Check the high-conductor fourth moment with a genuine interval of prime powers."""
    residues = units(q)
    g = len(residues)
    coeff = {}
    powers = [(n, wt) for n, wt in prime_powers(hi) if lo < n <= hi]
    for r in divisors(q):
        if r <= cutoff:
            continue
        for label in units(r):
            primitive = dirichlet_char(r, label)
            if not primitive.is_primitive():
                continue
            label_q = next(l for l in residues if all(
                int(dirichlet_char(q, l).chi_exponent(a)) * int(primitive.group().exponent())
                == int(primitive.chi_exponent(a)) * int(dirichlet_char(q, l).group().exponent())
                for a in residues))
            coeff[label_q] = sum(wt * primitive(n) for n, wt in powers)
    left = sum(abs(sum(dirichlet_char(q, l)(a).conjugate() * z for l, z in coeff.items()) / g)**4
               for a in residues)
    by_product = {a: acb(0) for a in residues}
    for l, z in coeff.items():
        for k, v in coeff.items():
            by_product[l * k % q] += z * v
    right = sum(abs(z)**2 for z in by_product.values()) / g**3
    check_zero(left - right)


def check_splitting():
    """Check the finite congruences and cover, independent of asymptotic prime windows."""
    upper, u, small = 105, 131, (3, 5, 7)
    base = u * math.prod(small)
    candidates = []
    for k in range(1, 20000):
        v = base * k - 1
        if v % 16 == 7 and fmpz(v).is_prime():
            candidates.append(v)
        if len(candidates) == 48:
            break
    assert len(candidates) == 48
    primes = [p for p in range(3, upper + 1) if fmpz(p).is_prime()]
    cover = {}
    for v in candidates:
        assert u % 8 == 3 and all((v + 1) % p == 0 for p in (*small, u))
        cover[v] = {q for q in range(1, upper + 1)
                    if all(pow(-v, (p - 1) // 2, p) == 1 for p in primes if q % p == 0)}
    missing = set(range(1, upper + 1))
    chosen = []
    while missing:
        v = max(candidates, key=lambda p: len(missing & cover[p]))
        assert missing & cover[v]
        chosen.append(v)
        missing -= cover[v]
    return chosen


if __name__ == "__main__":
    ctx.prec = 192
    moduli = (5, 8, 12, 15, 16, 17, 25, 35, 45, 56, 63, 65)
    total = 0
    maximum = 0.0
    for modulus in moduli:
        blocks, error = check_counts(modulus, (217, 512))
        total += 2 * len(units(modulus))
        maximum = max(maximum, error)
        print(f"q={modulus:2}  cyclic blocks={blocks:2}  enclosing error <= {error:.3g}")
    check_level_obstruction()
    for r, label, w in ((17, 3, 113), (17, 4, 113), (13, 2, 61), (25, 2, 41)):
        order, size = check_split_weil(r, label, w)
        print(f"r={r}, w={w}: {order} disjoint translates, {size} trace residues per cell")
    for q in (15, 28, 35, 45):
        check_fourth_moment(q, 5, 73, 217)
    assert any(d > 1 for d, _ in trace_blocks(16).values())
    print(f"{total} centered counts; maximum enclosing error {maximum:.3g}")
    print("q=17: best squared error for one level character is 7/256")
    print("v covering all q <= 105 with u=131:", check_splitting())
