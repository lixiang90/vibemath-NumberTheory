"""Exact cofactor identities and Arb box averages for the CM229 coprime-prime sum.

Run from any directory with python-flint 0.9.0 and the companion ERH files.
"""

from bisect import bisect_right
from math import gcd, isqrt

from flint import arb, ctx, fmpz

from verify_cm229 import coefficients
from verify_cm229_bilinear import CHI4
from verify_cm229_split_primes import odd_prime_sieve


Y = 4096
WINDOW_MILLIONTHS = ((532_182, 339_071, 655_849), (390_309, 261_359, 480_923))


def prime_list(limit):
    odd = odd_prime_sieve(limit)
    return [2] + [p for p in range(3, limit + 1, 2) if odd[p // 2]]


def inverse_coefficients(a, limit):
    u = [0] * (limit + 1)
    u[1] = 1
    for d in range(1, limit + 1):
        if u[d]:
            for j in range(2, limit // d + 1):
                u[d*j] -= u[d] * a[j]
    return u


def cofactor_tail(a, u, z, limit):
    c = [0] * (limit + 1)
    for d in range(z + 1, limit + 1):
        if u[d]:
            for j in range(1, limit // d + 1):
                c[d*j] += u[d] * a[j]
    return c


def check_cofactors(a, u, primes, limit):
    prime_set = {p for p in primes if p <= limit}
    for p in primes:
        if p > limit:
            break
        for d in range(1, limit // p + 1):
            if d % p and u[d*p] != -u[d] * a[p]:
                raise ArithmeticError(f"inverse coefficient at d={d}, p={p}")
    for z in (1, 2, 4, 7, 16, 32, 64):
        c = cofactor_tail(a, u, z, limit)
        small = [int(h == 1) for h in range(limit + 1)]
        for d in range(1, z + 1):
            for j in range(1, limit // d + 1):
                small[d*j] -= u[d] * a[j]
        if small != c:
            raise ArithmeticError(f"two truncations disagree at cutoff {z}")
        smooth = [1] * (limit + 1)
        for p in primes:
            if p > z:
                break
            power = p
            while power <= limit:
                for h in range(power, limit + 1, power):
                    smooth[h] *= p
                power *= p
        for h in range(2, limit + 1):
            s = smooth[h]
            if c[h] != a[h//s] * (c[s] - int(s == 1)):
                raise ArithmeticError(f"smooth-part formula at h={h}, z={z}")
            if z < h <= 2*z and c[h] != u[h]:
                raise ArithmeticError(f"first cofactor interval at h={h}, z={z}")
            if z < h <= z*z and s < h and h not in prime_set and c[h]:
                raise ArithmeticError(f"mixed composite at h={h}, z={z}")
    print(f"both convolution tails, smooth parts and coprime inverse identities through {limit}")


def quadratic_character(n):
    residue = pow(n, 114, 229)
    return -1 if residue == 228 else residue


def check_norm_support(a, u, primes, limit):
    possible = [False] * (limit + 1)
    partial = [0] * (limit + 1)
    for d in range(1, limit + 1):
        if u[d]:
            for j in range(2, limit // d + 1):
                h = d*j
                partial[h] -= u[d]*a[j]
                if partial[h]:
                    possible[h] = True
    norm = [True] * (limit + 1)
    for p in primes:
        if p > limit:
            break
        if quadratic_character(p) == -1:
            odd_power = p
            while odd_power <= limit:
                for j in range(1, limit // odd_power + 1):
                    if j % p:
                        norm[odd_power*j] = False
                odd_power *= p*p
    if possible[2:] != norm[2:]:
        raise ArithmeticError("union of cofactor supports differs from ideal norms")
    if limit == 8192 and sum(possible) != 2047:
        raise ArithmeticError("cofactor support count differs from the reported count")
    print("possible cofactors (2 through limit):", sum(possible))

    cut = 300_000
    if primes[-1] < cut:
        raise ArithmeticError("prime list does not cover the density product")
    value = -sum(quadratic_character(j) * (2*(arb.pi()*j/229).sin()).log()
                 for j in range(1, 229)) / arb(229).sqrt()
    product = arb(1)
    for p in primes:
        if p > cut:
            break
        if quadratic_character(p) == -1:
            product *= arb(p*p)/(p*p - 1)
    lower = (value * 229*product/(228*arb.pi())).sqrt()
    upper = lower * (arb(cut + 1)/cut).sqrt()
    if not arb("0.6886417") < lower < upper < arb("0.6886429"):
        raise ArithmeticError("density constant differs from the published enclosure")
    print("density constant enclosed between", lower.str(18), "and", upper.str(18))


def floor_ball(x):
    value = x.floor().unique_fmpz()
    if value is None:
        raise ArithmeticError(f"unresolved integer endpoint {x}")
    return int(value)


def prime_events(block, z, limit, a, u, primes):
    lower = int(fmpz(block**2).root(3))
    if (lower + 1)**2 <= limit:
        raise ArithmeticError("large prime is not guaranteed to be unique on this block")
    c = cofactor_tail(a, u, z, limit // (lower + 1))
    events = []
    for p in primes[bisect_right(primes, lower):bisect_right(primes, limit//(z + 1))]:
        if not a[p]:
            continue
        logp = arb(p).log()
        for h in range(z + 1, limit // p + 1):
            q = h*p
            weight = (1 + CHI4[q % 4]) * c[h] * a[p]
            if weight:
                if gcd(h, p) != 1:
                    raise ArithmeticError(f"noncoprime event at {q}")
                events.append((q, weight * logp / 2, h, p))
    events.sort(key=lambda event: event[0])
    positions = [event[0] for event in events]
    if any(x >= y for x, y in zip(positions, positions[1:])):
        raise ArithmeticError("two large-prime representations for one product")
    sums, scaled = [arb(0)], [arb(0)]
    for q, weight, _, _ in events:
        sums.append(sums[-1] + weight)
        scaled.append(scaled[-1] + weight / arb(q).sqrt())
    return events, positions, sums, scaled


def sector_prefixes(events, block, end, z, primes):
    small_primes = primes[:bisect_right(primes, z)]
    prime_h = set(primes[:bisect_right(primes, z*z)])
    sectors = [([], [arb(0)]), ([], [arb(0)])]
    for q, weight, h, p in events:
        if q > 2*(end - 1) or p*z*z <= 4*block:
            continue
        if h in prime_h:
            positions, sums = sectors[0]
        else:
            t = h
            for ell in small_primes:
                while t % ell == 0:
                    t //= ell
            if t != 1 or h >= z*z:
                raise ArithmeticError(f"unexpected surviving composite cofactor {h}")
            positions, sums = sectors[1]
        positions.append(q)
        sums.append(sums[-1] + weight)
    return sectors


def check_window(a, u, primes):
    y = Y
    budget, stop = y**3, isqrt(y**3 - 1) + 1
    m = int((arb(2*y).log()**2).ceil().unique_fmpz())
    r = (arb(2).log() / m).exp()
    rows = [[arb(0) for _ in range(3)] for _ in range(2)]
    sector_masses = [arb(0) for _ in range(3)]
    block = y
    while block < stop:
        end = min(2*block, stop)
        K = arb(block) ** (arb(7)/15)
        width = 1/K
        growth = width.exp()
        limit = floor_ball(2*(end - 1)*growth)
        original = max(1, floor_ball(arb(block)**(arb(1)/6) / arb(2*block).log()**3))
        sixth = int(fmpz(block).root(6))
        for variant, (row, z) in enumerate(zip(rows, (original, sixth), strict=True)):
            events, positions, sums, scaled = prime_events(block, z, limit, a, u, primes)
            sectors = sector_prefixes(events, block, end, z, primes) if variant == 1 else []

            def index(x):
                return bisect_right(positions, floor_ball(x))

            def primitive(x):
                hi, lo = index(2*x), index(2*x/r)
                return (2*arb(2).sqrt() * (scaled[hi] - scaled[lo]/r.sqrt())
                        - 2 * (sums[hi] - sums[lo]) / x.sqrt())

            for n in range(block, end):
                left = floor_ball(2*n/r)
                raw = sums[bisect_right(positions, 2*n)] - sums[bisect_right(positions, left)]
                x = arb(n)
                box = K*x.sqrt()*(primitive(x*growth) - primitive(x/growth))/2
                if n in (block, (block + end)//2, end - 1):
                    direct = arb(0)
                    for q, weight, _, _ in events[index(2*x/(r*growth)):index(2*x*growth)]:
                        start = (arb(q)/(2*n)).log().max(-width)
                        finish = (r*q/(2*n)).log().min(width)
                        if not start < finish:
                            raise ArithmeticError(f"empty box overlap at q={q}, n={n}")
                        direct += weight*K*((-start/2).exp() - (-finish/2).exp())
                    if not box.overlaps(direct):
                        raise ArithmeticError(f"box primitive disagrees with direct sum at n={n}")
                weight = isqrt(budget**2 // n**3)
                for j, value in enumerate((raw.max(0), box.max(0), abs(raw - box))):
                    row[j] += weight*value
                if sectors:
                    prime, smooth = (s[bisect_right(v, 2*n)] - s[bisect_right(v, left)]
                                     for v, s in sectors)
                    for j, part in enumerate((prime, smooth, prime + smooth)):
                        sector_masses[j] += weight*part.max(0)
        block = end
    for row, name, bounds in zip(rows, ("original", "sixth root"), WINDOW_MILLIONTHS, strict=True):
        values = [x/budget for x in row]
        if not all(lo < x*10**6 < lo + 1 for x, lo in zip(values, bounds, strict=True)):
            raise ArithmeticError(f"whole-window intervals differ for {name}")
        print(name, "sharp positive, box positive, absolute difference:",
              *(x.str(20) for x in values))
    sector_values = [x/budget for x in sector_masses]
    gap = sector_values[0] + sector_values[1] - sector_values[2]
    bounds = (177_860, 184_547, 266_649, 95_758)
    if not all(lo < x*10**6 < lo + 1 for x, lo in zip((*sector_values, gap), bounds, strict=True)):
        raise ArithmeticError("two cofactor classes differ from their published intervals")
    print("sixth-root cutoff, prime positive, smooth positive, joint positive, cancellation:",
          *(x.str(20) for x in (*sector_values, gap)))


if __name__ == "__main__":
    ctx.prec = 160
    event_limit = 3*isqrt(Y**3)
    primes = prime_list(event_limit // 2)
    a = coefficients(event_limit // 2)["real"]
    u = inverse_coefficients(a, 8192)
    check_cofactors(a, u, primes, 8192)
    check_norm_support(a, u, primes, 8192)
    check_window(a, u, primes)
