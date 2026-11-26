"""Certified signed energies for (10a) and cofactor values at two isolated zeros.

Requires python-flint 0.9.0 and the companion CM229 files. All interval ends,
quadrature nodes, Euler coefficients and error bounds are computed in the script.
"""

from math import comb, factorial

from flint import acb, arb, ctx, fmpz

from verify_cm229 import coefficients, read_rows
from verify_cm229_bilinear import CHI4
from verify_cm229_coprime_smoothing import cofactor_tail, floor_ball, inverse_coefficients, prime_list


ENERGY_LOWER_BILLIONTHS = {
    (False, "log"): (5179261, 3736268, 8022517),
    (False, "full"): (51100428, 25969693, 89883597),
    (True, "log"): (5084165, 3938653, 3327013, 3375333),
    (True, "full"): (28500458, 25665250, 30764813, 28801352),
}


def original_cutoff(T):
    return max(1, floor_ball(arb(T)**(arb(1)/6) / arb(2*T).log()**3))


def dyadic_blocks(top):
    d = 1
    while d <= top:
        yield d
        d *= 2


def rectangles(T, z, p0, r):
    result = []
    for d in dyadic_blocks(8*T//(p0+1)):
        h0 = max(d, z+1)
        if h0 >= 2*d:
            continue
        for p in dyadic_blocks(8*T//(z+1)):
            lo = max(p, p0+1)
            if lo >= 2*p or h0*lo >= 8*T or (2*d-1)*(2*p-1) <= T/r:
                continue
            result.append((d, p))
    return result


def signed_energy(T, z, short):
    m = int((arb(2*T).log()**2).ceil().unique_fmpz())
    v = arb(2).log()/m
    K = arb(T)**(arb(7)/15)
    if not v*(1+K) < 1:
        raise ArithmeticError("the weight of (10a) is not constant on the full band")
    p0 = int(fmpz(T*T).root(3))
    rects = rectangles(T, z, p0, v.exp())
    last_h = max(2*d-1 for d, _ in rects)
    last_p = max(2*p-1 for _, p in rects)
    a = coefficients(max(last_p, last_h))["real"]
    primes = prime_list(last_p)
    prime_set = set(primes)
    c = cofactor_tail(a, inverse_coefficients(a, last_h), z, last_h)
    h_parts, p_parts = {}, {}
    # Bands are the original dyadic rectangles. The short test additionally clips h.
    for d in sorted({d for d, _ in rects}):
        parts = ([], []) if short else ([],)
        for h in range(max(d, z+1), min(2*d, z*z+1) if short else 2*d):
            if not c[h]:
                continue
            group = int(short and h not in prime_set)
            if group and any(h % q == 0 for q in primes if z < q <= h):
                raise ArithmeticError("surviving composite is not smooth")
            parts[group].append((h, c[h]/arb(h).sqrt(), CHI4[h % 4], arb(h).log()))
        h_parts[d] = parts
    for p in sorted({p for _, p in rects}):
        p_parts[p] = [(q, a[q]*arb(q).log()/arb(q).sqrt(), CHI4[q % 4], arb(q).log())
                      for q in primes if max(p, p0+1) <= q < 2*p and a[q]]

    # For any union of the cells this bounds both its coefficient l1 norm and
    # the square root of the corresponding disjoint-cell energy majorant.
    mass = arb(0)
    q_min, q_max = None, 0
    for d, p in rects:
        for part in h_parts[d]:
            if not part or not p_parts[p]:
                continue
            mass += sum(abs(x[1]) for x in part) * sum(abs(x[1]) for x in p_parts[p])
            lo, hi = part[0][0]*p_parts[p][0][0], part[-1][0]*p_parts[p][-1][0]
            q_min = lo if q_min is None else min(q_min, lo)
            q_max = max(q_max, hi)
    omega = (arb(q_max)/q_min).log()

    def pair(part, t):
        plus, minus = acb(0), acb(0)
        for _, coef, chi, logn in part:
            value = coef * acb(0, -t*logn).exp()
            plus += value
            minus += chi*value
        return plus, minus

    def cells(t):
        hd = {d: [pair(part, t) for part in parts] for d, parts in h_parts.items()}
        bp = {p: pair(part, t) for p, part in p_parts.items()}
        values = [[acb(0), acb(0)] for _ in range(2 if short else 1)]
        for d, p in rects:
            for g, hp in enumerate(hd[d]):
                for branch in range(2):
                    values[g][branch] += hp[branch]*bp[p][branch]/2
        return values

    # Check the signed half sum independently before integrating.
    for t in (arb(0), arb(3)/7):
        direct = acb(0)
        for d, p in rects:
            for part in h_parts[d]:
                for h, ch, _, _ in part:
                    for q, aq, _, _ in p_parts[p]:
                        proj = 1 + CHI4[(h*q) % 4]
                        if proj:
                            direct += proj*ch*aq*acb(0, -t*arb(h*q).log()).exp()/2
        if not direct.overlaps(sum((sum(row) for row in cells(t)), acb(0))):
            raise ArithmeticError("direct mod-four projection differs from the factorization")

    def energy_groups(t):
        vals = cells(t)
        if not short:
            plus, minus = vals[0]
            return [abs(plus)**2, abs(minus)**2, abs(plus+minus)**2]
        pplus, pminus = vals[0]
        splus, sminus = vals[1]
        return [sum(abs(y)**2 for row in vals for y in row),
                abs(pplus+pminus)**2 + abs(splus+sminus)**2,
                abs(pplus+splus)**2 + abs(pminus+sminus)**2,
                abs(pplus+pminus+splus+sminus)**2]

    # Integer samples after modulation. Check the cardinal-series identity
    # independently of the energy quadrature, with its absolute tail bound.
    def modulated(t):
        value = sum((sum(row) for row in cells(t)), acb(0))
        return acb(0, t*arb(2*T).log()).exp()*value

    radius, power = 160, 14
    b, test_t = arb(1)/4, arb(3)/7
    samples = [modulated(arb(j)) for j in range(radius+1)]
    last = floor_ball(arb(2*T).log())
    peak = max(range(last+1), key=lambda j: float(abs(samples[j]).mid()))
    if not all(abs(samples[j]) < abs(samples[peak]) for j in range(last+1) if j != peak):
        raise ArithmeticError("maximum on the logarithmic sampling grid is not separated")
    peak_index, peak_lo = (0, 2791384673) if short else (9, 5966869266)
    if last != 9 or peak != peak_index or not peak_lo < abs(samples[peak])*10**9 < peak_lo+1:
        raise ArithmeticError("the sampled maximum differs from the published value")
    interp = acb(0)
    for j in range(-radius, radius+1):
        delta = test_t-j
        freq = arb.pi()*delta
        base = b*delta/power
        kernel = freq.sin()/freq*(base.sin()/base)**power
        interp += kernel*(samples[j] if j >= 0 else samples[-j].conjugate())
    tail = 2*mass/(arb.pi()*power)*(power/(b*(radius-1)))**power
    if not abs(interp-modulated(test_t)) < tail < arb("2.1e-6" if short else "6.5e-6"):
        raise ArithmeticError("finite integer sampling misses its certified tail bound")
    print("short" if short else "original", "integer grid 0 to", last, "peak at", peak,
          "modulus", abs(samples[peak]).str(16), "sampling tail bound", tail.str(4))

    order = 20
    rule = [arb.legendre_p_root(order, j, weight=True) for j in range(order)]
    for name, cutoff in (("log", arb(2*T).log()), ("full", K)):
        panels = int((cutoff/2).ceil().unique_fmpz())
        half_width = cutoff/(2*panels)
        sums = [arb(0) for _ in range(4 if short else 3)]
        for j in range(panels):
            center = (2*j+1)*half_width
            for node, weight in rule:
                for i, value in enumerate(energy_groups(center+half_width*node)):
                    sums[i] += weight*value
        # Taylor remainder of degree 2*order-1; integration and the positive rule
        # each cost one uniform remainder. Integrands are even in t.
        err = 4*cutoff*v*v*mass**2*(omega*half_width)**(2*order)/factorial(2*order)
        sums = [2*v*v*half_width*s + arb(0, err) for s in sums]
        if not err < arb("8e-25"):
            raise ArithmeticError("the quadrature error exceeds the published bound")
        for value, lo in zip(sums, ENERGY_LOWER_BILLIONTHS[short, name], strict=True):
            if not lo < 10**9*value < lo+1:
                raise ArithmeticError("the energy is outside the published interval")
        print("short" if short else "original", name, "rectangles", len(rects),
              "cutoff", cutoff.str(12), "quadrature bound", err.str(4))
        print("energies", *(s.str(16) for s in sums))


def zero_modes():
    # A concrete original z, and a specified finite prefactor for H.
    T = 2**150
    logT = arb(2*T).log()
    z = original_cutoff(T)
    H = floor_ball(z*logT**2/1024)
    if not z < H or 2*H > z*z:
        raise ArithmeticError("cofactor ranges overlap the mixed-composite region")
    a = coefficients(2*H)["real"]
    c = cofactor_tail(a, inverse_coefficients(a, 2*H), z, 2*H)
    prime_set = set(prime_list(2*H))
    m = 12*int(logT.ceil().unique_fmpz()) + 20
    if (z, H, m) != (29, 310, 1280):
        raise ArithmeticError("unexpected cutoff in the finite zero test")
    # Exact finite CDF for the sum of m continuous independent uniforms.
    with ctx.workprec(3072):
        binomials = [(-1)**j*comb(m, j) for j in range(m//2+1)]
        fact = fmpz(factorial(m))
        weights = {}
        for h in range(z+1, 2*H):
            if not c[h]:
                continue
            if h <= H:
                weight = arb(1)
            else:
                w = m*(arb(h)/H).log()/arb(2).log()
                small = w.min(m-w)
                last = floor_ball(small)
                tail = sum(binomials[j]*(small-j)**m for j in range(last+1)) / fact
                weight = (1-tail if w < arb(m)/2 else tail).max(0).min(1)
            weights[h] = weight
    print("zero tests: T=2^150, z, H, m =", z, H, m)
    for name, intervals in read_rows().items():
        lo, hi, den = intervals[0]
        gamma = arb(lo).union(arb(hi))/den
        s = acb(arb(1)/2, gamma)
        prime, smooth = acb(0), acb(0)
        for h, weight in weights.items():
            coef = c[h] if name == "real" else CHI4[h % 4]*c[h]
            value = coef*weight*(-s*arb(h).log()).exp()
            if h in prime_set:
                prime += value
            else:
                smooth += value
        centers = (("0.384940", "0.270366", "0.001758", "0.050398", "0.38670", "0.320763")
                   if name == "real" else
                   ("1.051617", "0.231197", "0.064777", "-0.071290", "1.116394", "0.159908"))
        total = prime+smooth
        coords = (prime.real, prime.imag, smooth.real, smooth.imag, total.real, total.imag)
        if not all(abs(value-arb(center)) < arb("1e-5")
                   for value, center in zip(coords, centers, strict=True)):
            raise ArithmeticError("zero-mode evaluation is outside the published enclosure")
        print(name, "first zero", gamma.str(10), "prime", prime.str(9), "smooth", smooth.str(9),
              "sum", total.str(9))


if __name__ == "__main__":
    ctx.prec = 160
    signed_energy(4096, original_cutoff(4096), False)
    signed_energy(4096, 4, True)
    zero_modes()
