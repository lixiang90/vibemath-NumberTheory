"""Integer and rational certificates for the 3-or-5 least-prime bound."""

import json
from fractions import Fraction as F
from math import factorial, isqrt, prod

from quintic_character import harmonic_blocks


def require(condition, message):
    if not condition:
        raise ValueError(message)


def scaled_lower(values, count, scale=10**12):
    return F(sum((scale*values[n % 5])//n for n in range(1, count+1)), scale)


def minus_exp_lower(x, degree=15):
    """Odd Taylor truncation is a lower bound on exp(-x) for x >= 0."""
    require(x >= 0 and degree % 2 == 1, "exponential domain")
    return sum((-x)**k/factorial(k) for k in range(degree+1))


def quintic_denominator():
    scale = 10**14
    minima = harmonic_blocks(10_000_000, [1, 40], scale)
    require(minima[0][1] >= scale*16//100 and minima[1][1] >= scale*3//10,
            "finite harmonic bounds for the quintic characters")
    m = 10_000
    psi = scaled_lower((0, 1, -1, -1, 1), m) - F(2, m)
    real = scaled_lower((0, 1, 0, 0, -1), m) - F(3, 2*m)
    imag = scaled_lower((0, 0, 1, -1, 0), m) - F(3, 2*m)
    require(min(psi, real, imag) > 0 and psi*(real**2+imag**2) > F(339, 1000),
            "character L-values")
    u = 215
    require(u**3 < 10_000_000 and
            F(339, 1000) - F(8, u) - F(33, 4*u**2) - F(9, 2*u**3) > F(3, 10),
            "infinite harmonic bound")
    primes = [p for p in range(2, 101) if all(p % i for i in range(2, isqrt(p)+1))]
    c = F(5, 4)
    for p in primes:
        if p % 5 in (2, 3):
            c /= 1-F(1, p**4)
        if p % 5 == 4:
            c /= (1-F(1, p**2))**2
    log_tail = F(2, 3*100**3) + F(3, 100*9)
    require(c/(1-log_tail) < F(11, 8), "finite and infinite convolution factors")
    require(F(3, 20)/F(11, 8) == F(6, 55), "quintic quadratic coefficient")
    split = [q for q in range(11, 1001, 10) if all(q % p for p in range(2, isqrt(q)+1))]
    correction = prod(F(q-1, q-5) / (1+F(1, q))**4 for q in split)
    require(F(25, 2)/(1-F(5, 11))+2 < 25 and
            correction/(1-F(25, 100*99)) < F(5, 4), "uniform exception product")
    small_term = F(21**2, 12)/F(27, 10)**15
    main = 20*(1+F(1, 144))*F(4*55, 6)
    require((main+small_term)*(1+F(1, 10**6))**2 < 740 and
            13**30*40**14 < F(27, 10)**135, "uniform quintic reciprocal tail")
    require(F(740*30, 21) < 1058 and 3*(188+1058+2) == 3744, "all-rank constants")
    return minima


def large_sieve_constants():
    e_lower = sum(F(1, factorial(j)) for j in range(7))
    primes = [q for q in range(2, 101) if all(q % p for p in range(2, isqrt(q)+1))]

    def log_upper(q):
        # Certify log(q) < j/100 by a strict integer-power comparison.
        j = 1
        while e_lower**j <= q**100:
            j += 1
        return F(j, 100)

    for ell, gain, shift in [(3, F(1149, 1000), F(77, 100)), (5, F(1165, 1000), 1)]:
        split = [q for q in primes if q % ell == 1]
        mass = prod((1-F(1, q))**ell/(1-F(ell, q)) for q in split)
        mean = sum(F(ell*(ell-1), (q-ell)*(q-1))*log_upper(q) for q in split)
        require(mass > gain and mean < shift, "finite sieve-prime gain")
    # 2*log(24.1) + 2*a_l + 2*mean_l, using the upper bounds just proved.
    require(F(241, 10)**600*74**100 < e_lower**2343 and
            F(241, 10)**3000*40**1400 < e_lower**14730, "logarithmic shifts")
    tiny = F(1, 10**6)
    require(F(10**14, 24) - F(10**15, 241) >= 1, "floor of sieve level")
    require(F(64**2, 24)/e_lower**32 < F(1, 10_000), "small primes before sieving")
    require((1+F(1, 10**28-1))**2 < 1+tiny, "partial summation at X")
    require(F(10**28, 10**28-5) < 1+tiny, "loss at a large prime")
    require((1+tiny)*F(6, 840)*(F(2, 73**2)+F(1, 73))+3*tiny < F(1, 8000),
            "progression and multiplicity error")
    z0 = F(1, 1000)
    constants = [(3, F(7, 50)*F(1149, 1000), F(321, 100), F(319, 100)),
                 (5, F(6, 55)*F(1165, 1000), F(407, 100), F(404, 100))]
    for ell, coefficient, target, bounded_target in constants:
        main = F(1, 8)*(1+F(15, 24**2)+F(15, 10**28))*4/coefficient
        log_loss = (1+tiny)*(ell-1)*z0+F(1, 8000)
        require((main+F(1, 10_000))*(1+tiny)/(1-log_loss) < target, "integrated sieve constant")
        # t>=10^28. (log t)^2/sqrt(t) is decreasing there and 28*log(10)<65.
        small = F((10**8+1)*65**2, 24*10**14)
        require((main*(1+tiny)+small)*(1+tiny) < bounded_target, "bounded-exponent sieve")
    # Verify log 2 < 0.699, log 2 > 0.693 and 2.302 < log 10 < 65/28.
    log_two_lower = 2*sum(F(1, 3**(2*j+1)*(2*j+1)) for j in range(3))
    log_two_upper = log_two_lower + F(2, 7*3**7)/(1-F(1, 9))
    require(F(693, 1000) < log_two_lower and (1+tiny)*log_two_upper < F(7, 10),
            "total reciprocal mass")
    y = F(1151, 500)
    require(sum(y**j/factorial(j) for j in range(17)) +
            y**17/factorial(17)/(1-y/18) < 10, "logarithmic threshold")
    require(sum(F(65, 28)**j/factorial(j) for j in range(9)) > 10, "upper log 10")
    bound = (1+tiny)*(F(1, 2500)+F(2, 840)*(F(2, 73**2)+F(1, 73)))
    require(bound < F(1, 2000), "short prime interval")
    b = 28*y
    require(4*F(319, 100)/(b-F(935, 100)) + 6*F(404, 100)/(b-F(1182, 100)) +
            F(1, 2000) < F(693, 1000), "bounded-exponent contradiction")


def mass_iterations():
    z = F(1, 1000)
    # Both entries of a row use only the upper bounds in the preceding row.
    table = [(700, 700), (79, 98), (64, 76), (61, 74)]
    for (previous3, previous5), (next3, next5) in zip(table, table[1:]):
        u, v = F(previous3, 1000)+z, min(F(7, 10), F(3*previous5, 1000))+z
        upper3 = F(321, 100)/(F(736, 10)-F(19, 2))*(F(5, 2)-F(3, 2)*minus_exp_lower(2*v/3))
        upper5 = F(407, 100)/(F(736, 10)-12)*(F(3, 2)-F(1, 2)*minus_exp_lower(4*u))
        require(upper3 < F(next3, 1000) and upper5 < F(next5, 1000), "upper iteration")
    require(F(4*table[-1][0]+6*table[-1][1], 1000)+z < F(693, 1000), "abundancy contradiction")
    return table


def verify():
    minima = quintic_denominator()
    large_sieve_constants()
    return {"harmonic_scale": 10**14, "finite_harmonic_minima": minima,
            "per_class_upper_pairs_scaled_1000": mass_iterations(), "mixed_least_prime_below": 10**32,
            "non_euler_exponent_cap": 10**8, "capped_least_prime_below": 10**28}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
