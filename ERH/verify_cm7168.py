"""Certify the quartic CM field of discriminant 7168, its reflex and their closure.

Requires python-flint 0.9.0. Use --regenerate to reconstruct the zero table.
"""

import argparse
from dataclasses import replace
from math import gcd, prod
from pathlib import Path

from flint import acb_series, arb, ctx, dirichlet_char, fmpz, fmpz_poly, nmod_poly

from verify_cm229 import log_derivative_at_one, read_rows
from verify_cm229_256 import (
    RATIONAL_DEN, RULE, CompletedRealDirichlet, indexed_zeta_count,
    indexed_zeta_intervals, separation, write_rows,
)
from verify_cm229_split_primes import odd_prime_sieve
from verify_tetra3136 import Completed, GammaFactors, half_winding, locate, verify_signs


FILE = Path(__file__).with_name("cm7168-zeros256.csv")
ZERO_RULE = replace(RULE, denominator=RATIONAL_DEN)
CHARACTERS = {8: 5, 28: 27, 56: 27}
NAMES = ("zeta", "chi8", "chi28", "chi56", "dihedral")
POINTS = (256, 4096, 65536)
UPPER_B = ".87333263744590506611657458485"


def denominator(degrees):
    return prod(fmpz_poly([1]+[0]*(degree-1)+[-1]) for degree in degrees)


def local_factor(p):
    """Return trace, determinant, and the residue degrees in the base and quartic."""
    quartic = [f.degree() for f, _ in nmod_poly([7, 0, 6, 0, 1], p).factor()[1]]
    if p == 2:
        trace, det, base = 0, 0, [1]
    elif p == 7:
        trace, det, base = 1, 0, [1, 1]
    else:
        two, seven = fmpz(2).jacobi(p), fmpz(7).jacobi(p)
        base = [1, 1] if two == 1 else [2]
        trace, det = 0, int(two*seven)
        if two == seven == 1:
            root = fmpz(2).sqrtmod(p)
            trace = int(2*(-3-root).jacobi(p))
            if trace != 2*(-3+root).jacobi(p):
                raise ArithmeticError("choice of square root changes the Frobenius class")
    if denominator(quartic) != denominator(base)*fmpz_poly([1, -trace, det]):
        raise ArithmeticError(f"incompatible Euler factors at {p}")
    return trace, det, base, quartic


def local_data(limit):
    """Build a(n), explicit-formula sums, and central splitting counts independently."""
    top = max(limit, *POINTS)
    sieve = odd_prime_sieve(top)
    coeff = [0]+[1]*limit
    sums = {x: [arb(0), arb(0)] for x in POINTS}
    counts = {x: [0, 0, arb(0), arb(0)] for x in POINTS}
    for p in range(2, top+1):
        if p != 2 and (p % 2 == 0 or not sieve[p//2]):
            continue
        trace, det, base, quartic = local_factor(p)
        logp = arb(p).log()
        if p not in (2, 7) and abs(trace) == 2:
            side = int(trace < 0)
            for x, record in counts.items():
                if p <= x:
                    record[side] += 1
                    record[2+side] += logp
        power, j = p, 1
        previous, current = 0, 1
        while power <= top:
            previous, current = current, trace*current-det*previous
            if power <= limit:
                for m in range(1, limit//power+1):
                    if m % p:
                        coeff[m*power] *= current
            base_log = sum(f for f in base if j % f == 0)
            quartic_log = sum(f for f in quartic if j % f == 0)
            for x, pair in sums.items():
                if power <= x:
                    weight = logp*(arb(1)/power-arb(1)/x)
                    pair[0] += (quartic_log-base_log)*weight
                    pair[1] += quartic_log*weight
            power *= p
            j += 1
    return coeff, sums, counts


def certify(function, rows, regenerate, divisions):
    half, margin = half_winding(function, mesh=128, degree=32)
    if regenerate:
        rows = locate(function, divisions=divisions)
    verify_signs(function, rows)
    if (2*half).unique_fmpz() != 2*len(rows):
        raise ArithmeticError("the isolated line zeros do not exhaust the rectangle")
    print("count", len(rows), "top margin", margin.str(8), "first and last",
          rows[0], rows[-1], flush=True)
    return rows


def rational_certificates(blocks, regenerate):
    nz = indexed_zeta_count(ZERO_RULE.height)
    zeta_rows = indexed_zeta_intervals(nz)
    if not regenerate and blocks["zeta"] != zeta_rows:
        raise ArithmeticError("indexed zeta zeros differ from the table")
    blocks["zeta"] = zeta_rows
    slopes = {}
    for conductor, index in CHARACTERS.items():
        chi = dirichlet_char(conductor, index)
        if chi.parity() != 0 or any(chi(a) != fmpz(conductor).jacobi(a)
                                    for a in range(1, conductor, 2) if gcd(a, conductor) == 1):
            raise ArithmeticError("wrong primitive quadratic character")
        name = "chi"+str(conductor)
        print(name, flush=True)
        blocks[name] = certify(CompletedRealDirichlet(chi), blocks.get(name), regenerate, 16)
        jet = acb_series([1, 1], prec=2).dirichlet_l(chi)
        if not (jet[0].real > 0 and jet[0].imag.contains(0) and jet[1].imag.contains(0)):
            raise ArithmeticError("invalid Dirichlet L-jet at one")
        slopes[name] = jet[1].real/jet[0].real
    return slopes


def low_zeros(rows):
    return [arb(a).union(arb(b))/den for a, b, den in rows]


def spectral_constants(blocks, slopes, ell):
    pi, euler = arb.pi(), arb.const_euler()
    masses = {"zeta": 2+euler-(4*pi).log()}
    for conductor in CHARACTERS:
        name = "chi"+str(conductor)
        masses[name] = 2*slopes[name]+(conductor/pi).log()+(arb(1)/2).digamma()
    total = 2*ell+arb(896).log()-2*pi.log()+2*arb(1).digamma()
    low = {name: sum(2/(arb(1)/4+t*t) for t in low_zeros(rows))
           for name, rows in blocks.items()}
    masses["dihedral"] = total
    if any(not masses[name] > low[name] for name in blocks):
        raise ArithmeticError("partial mass exceeds total mass")
    ratio = (1+(arb(3)/(4*ZERO_RULE.height))**2).sqrt()
    tail = ratio*(total-low["dihedral"])
    if not 0 < tail < arb(".020155"):
        raise ArithmeticError("dihedral tail exceeds the stated bound")
    for base, bound in (("chi8", "2.492753"), ("chi28", "2.830066")):
        critical = masses["zeta"]+masses[base]+low["dihedral"]
        if not critical < arb(bound):
            raise ArithmeticError("critical-line coefficient exceeds the stated bound")
        print(base, "critical and tail coefficients", critical.str(13), tail.str(13))
    target = arb(".87333")
    scale_bound = (arb("2.492753")*arb(1024)**(arb(1)/2-target)
                   +arb(".020155")*(1350000*(arb(UPPER_B)-target)).exp())
    if not scale_bound < arb(".897"):
        raise ArithmeticError("finite scale-range bound failed")
    print("finite scale-range coefficient", scale_bound.str(15))
    print("relative critical coefficient", low["dihedral"].str(13),
          "base masses", *(masses[name].str(13) for name in ("zeta", "chi8", "chi28")))
    return masses, low, tail


def check_geometric_data(sums, counts, blocks, slopes, ell, masses, low, tail):
    euler, two_pi = arb.const_euler(), (2*arb.pi()).log()
    gamma_f = euler+slopes["chi8"]
    h0 = arb(896).log()/4-euler-2*two_pi
    height = h0+ell/2
    if not arb("-2.834176474691678") < height < arb("-2.834176474691677"):
        raise ArithmeticError("CM height differs from its published enclosure")
    c_k = arb(7168).log()-4*(euler+two_pi)
    print("Artin log derivative", ell.str(20), "gamma_F", gamma_f.str(20),
          "CM height", height.str(20))
    ordinates = {name: low_zeros(blocks[name]) for name in ("zeta", "chi8", "dihedral")}
    for x, (s_rel, s_full) in sums.items():
        q, logx = arb(1)/x, arb(x).log()
        arch_rel = q*q.atanh()+(-q*q).log1p()/2
        arch_full = q+(1-q)*(-q).log1p()
        rel = (h0-s_rel/2-q*(arb(896).log()/4+two_pi)-arch_rel)/(1+q)
        ek = (logx-1-s_full+q*(1+logx-c_k)-2*arch_full)/(1+q)
        full = h0+(ek-gamma_f)/2
        partial = {name: -sum(2*(t*logx).cos()/(arb(1)/4+t*t) for t in ts)*q.sqrt()
                   for name, ts in ordinates.items()}
        tail_rel = tail*arb(x)**(arb(UPPER_B)-1)
        tail_rat = (masses["zeta"]+masses["chi8"]-low["zeta"]-low["chi8"])*q.sqrt()
        if not (abs(2*(1+q)*(rel-height)-partial["dihedral"]) < tail_rel
                and abs(2*(1+q)*(full-height)-sum(partial.values())) < tail_rel+tail_rat):
            raise ArithmeticError("prime sums and zero sums disagree")
        plus, minus, theta_plus, theta_minus = counts[x]
        print(x, "height errors", (rel-height).str(16), (full-height).str(16),
              "central counts", plus, minus, "centered theta",
              (theta_plus-arb(x)/8).str(12), (theta_minus-arb(x)/8).str(12))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true")
    regenerate = parser.parse_args().regenerate
    blocks = {} if regenerate else read_rows(FILE, NAMES)
    with ctx.workprec(192):
        slopes = rational_certificates(blocks, regenerate)
        coeff, sums, counts = local_data(12288)
    with ctx.workprec(2688):
        function = Completed(coeff, GammaFactors(896, (1, 1)), radius=1792, rule=ZERO_RULE)
        print("dihedral", flush=True)
        blocks["dihedral"] = certify(function, blocks.get("dihedral"), regenerate, 64)
        ell = log_derivative_at_one(function)
    with ctx.workprec(192):
        gap = separation(blocks)
        print("minimum gap among five factors", gap, "/", RATIONAL_DEN)
        sizes = {name: len(rows) for name, rows in blocks.items()}
        if tuple(sizes[name] for name in NAMES) != (111, 195, 246, 274, 498) or gap != 1111:
            raise ArithmeticError("zero counts or separation differ from the published table")
        print("positive zeros", sizes, "quartic", sizes["zeta"]+sizes["chi8"]+sizes["dihedral"],
              "reflex", sizes["zeta"]+sizes["chi28"]+sizes["dihedral"],
              "closure with multiplicity", sum(sizes.values())+sizes["dihedral"])
        masses, low, tail = spectral_constants(blocks, slopes, ell)
        check_geometric_data(sums, counts, blocks, slopes, ell, masses, low, tail)
    if regenerate:
        write_rows(FILE, blocks)


if __name__ == "__main__":
    main()
