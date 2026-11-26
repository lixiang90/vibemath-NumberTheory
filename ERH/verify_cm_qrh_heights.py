"""Check the two finite-prime height formulas for Q(i) and the CM229 sextic.

Requires python-flint 0.9.0 and the CM229 companion files in this directory.
All parameters, finite Euler factors, spectral sums and bounds use intervals.
"""

from flint import acb, arb, ctx, nmod_poly

from verify_cm229 import FACTORS, coefficients, log_derivative_at_one, read_rows
from verify_cm229_split_primes import odd_prime_sieve
from verify_tetra3136 import CUTOFF, Completed


POINTS = (256, 4096, 65536)
ERROR_FLOORS = {
    False: ((440208, 681058), (-296446, -297794), (15160, -9353)),
    True: ((16262349, 17202168), (-20238363, -17016608), (3245038, 3836463)),
}


def prime_sums(cubic):
    top = max(POINTS)
    primes = odd_prime_sieve(top)
    sums = {x: [arb(0), arb(0)] for x in POINTS}
    for p in range(2, top+1):
        if p != 2 and (p % 2 == 0 or not primes[p//2]):
            continue
        degrees = [1]
        if cubic:
            factors = nmod_poly([1, -4, 0, 1], p).factor()[1]
            degrees = [g.degree() for g, _ in factors]
        chi = (0, 1, 0, -1)[p % 4]
        degrees_cm = []
        for f in degrees:
            if not chi:
                degrees_cm.append(f)
            elif chi**f == 1:
                degrees_cm.extend([f, f])
            else:
                degrees_cm.append(2*f)
        j, power, logp = 1, p, arb(p).log()
        while power <= top:
            base = sum(f for f in degrees if j % f == 0)
            full = sum(f for f in degrees_cm if j % f == 0)
            relative = base*chi**j
            if full != base+relative:
                raise ArithmeticError("local relative factor disagrees with residue degrees")
            for x, pair in sums.items():
                if power <= x:
                    weight = logp*(arb(1)/power-arb(1)/x)
                    pair[0] += relative*weight
                    pair[1] += full*weight
            j += 1
            power *= p
    return sums


def spectral_check(x, actual, positive_heights, mass):
    weights = [2/(arb(1)/4+t*t) for t in positive_heights]
    partial = -sum(w*(t*arb(x).log()).cos()
                   for w, t in zip(weights, positive_heights, strict=True))/arb(x).sqrt()
    remaining = mass-sum(weights)
    if not remaining>0:
        raise ArithmeticError("isolated zeros exceed the total real spectral mass")
    tail = arb(x)**(-arb(1)/8)*4/arb(7).sqrt()*remaining
    if not abs(actual-partial)<tail:
        raise ArithmeticError("explicit formula differs from the independently bounded zero sum")


def check_field(cubic, ell_even, ell_odd):
    d, df = (3, 229) if cubic else (1, 1)
    de = 4**d*df**2
    lc = arb(de//df).log()
    two_pi = (2*arb.pi()).log()
    euler = arb.const_euler()
    qtr = arb(1)/4
    ell_beta = euler+two_pi-2*(qtr.lgamma()-(3*qtr).lgamma())
    gamma_f = euler+ell_even
    h0 = lc/4-d*euler/2-d*two_pi
    height = h0+(ell_beta+ell_odd)/2
    const0 = arb(de).log()-2*d*(euler+two_pi)
    # Sum of Re(1/(rho(1-rho))) for zeta_F is twice the logarithmic
    # derivative at 1 of the completed entire function with its poles removed.
    mass = 2+euler-(4*arb.pi()).log()
    ordinates = [z.imag for z in acb.zeta_zeros(1, 48)]
    if cubic:
        mass += 2*(ell_even+FACTORS["real"].log_scale+(arb(1)/2).digamma())
        ordinates.extend(arb(lo).union(arb(hi))/den for lo, hi, den in read_rows()["real"])
    print("degree", 2*d, "gamma_F", gamma_f.str(18), "average height", height.str(18))
    for (x, (s_eta, s_full)), floors in zip(prime_sums(cubic).items(), ERROR_FLOORS[cubic],
                                             strict=True):
        q, logx = arb(1)/x, arb(x).log()
        odd_arch = q*q.atanh()+(-q*q).log1p()/2
        cm_arch = q+(1-q)*(-q).log1p()
        rel = (h0-s_eta/2-q*(lc/4+d*two_pi/2)-d*odd_arch/2)/(1+q)
        ek = (logx-1-s_full+q*((d-1)*(1+logx)-const0)-d*cm_arch)/(1+q)
        full = h0+(ek-gamma_f)/2
        bound = (arb(30)/7)*(arb(de).log()+3)*q**(arb(1)/8)/(1+q)
        if not abs(full-height)<bound:
            raise ArithmeticError("height outside the explicit all-field bound")
        spectral_check(x, 2*(1+q)*(full-rel), ordinates, mass)
        for error, floor in zip((rel-height, full-height), floors, strict=True):
            if not floor < 10**9*error < floor+1:
                raise ArithmeticError("finite-prime height error outside the reported interval")
        print(x, "relative error", (rel-height).str(16),
              "full error", (full-height).str(16))


def main():
    # The existing Fourier certificates use 896 bits and control the
    # differential quotient by a Cauchy remainder on a fixed disk.
    ctx.prec = 896
    coeff = coefficients(CUTOFF)
    at_one = {name: log_derivative_at_one(Completed(coeff[name], factors, radius=512))
              for name, factors in FACTORS.items()}
    with ctx.workprec(160):
        check_field(False, arb(0), arb(0))
        check_field(True, at_one["real"], at_one["cm_twist"])


if __name__ == "__main__":
    main()
