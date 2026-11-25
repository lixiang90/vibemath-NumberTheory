"""Checks for the factor-dependent sieve levels and the cubic inverse statistic.

Run with NumPy, SciPy and python-flint from the published directory tree.
The integer samples and the finite-conductor zero comparison use different data.
"""

from fractions import Fraction
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from math import exp, factorial, log, pi
from pathlib import Path

import numpy as np
from flint import arb, arb_series, ctx, fmpz
from scipy.integrate import quad
from scipy.special import digamma


ROOT = Path(__file__).resolve().parents[1]
A = Fraction(3, 11)
BREAKS = (A, Fraction(1, 3), Fraction(3, 8), Fraction(2, 5), Fraction(1, 2))
OMEGA = 22 * pi / 39


def article(slug):
    spec = spec_from_file_location(slug.replace("-", "_"), ROOT / slug / "verify.py")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def level(u):
    return max(Fraction(3, 5), min((2 + 3*u) / 5, 1-u))


def sieve_constants(primes):
    for lo, hi in zip(BREAKS, BREAKS[1:]):
        for step in range(101):
            u = lo + (hi-lo)*step/100
            two_factor = min((2+3*u)/5, 1-u)
            d2 = two_factor-u
            assert d2 >= 0 and 2*u+d2 <= 1 and 2*u+5*d2 <= 2
            assert 2*u+d2 == 1 or 2*u+5*d2 == 2
            assert level(u) >= Fraction(3, 5)
    assert level(BREAKS[2]) == Fraction(5, 8)
    ctx.prec = 180
    partial = arb(1)
    for p in primes[1:]:
        partial *= 1-arb(1)/(int(p)-1)**2
    lower = partial*(1-arb(1)/999_999)

    def log_ratio(a, b):
        return (arb(a)/b).log()

    cubic = arb(53)/198 + arb(8)/3*log_ratio(25, 24)
    two_factor = arb(5)/3*log_ratio(275, 248) + log_ratio(5, 4)
    unweighted = (arb(5)/3*log_ratio(4, 3) + arb(5)/2*log_ratio(9, 8)
                  + log_ratio(16, 15) - arb(3)/2*log_ratio(25, 24)
                  + log_ratio(10, 9) + arb(1)/15 + arb(5)/3*log_ratio(3, 2))/2
    assert lower*cubic > arb(24857)/100_000 and partial*cubic < arb(24858)/100_000
    assert lower*unweighted > arb(53639)/100_000 and partial*unweighted < arb(53640)/100_000
    assert lower*two_factor > arb(26101)/100_000 and partial*two_factor < arb(26102)/100_000
    for integrand, expected in ((lambda u: 1/float(level(u)), cubic),
                                (lambda u: 1/(2*float(level(u))*u*(1-u)), unweighted)):
        integral = sum(quad(integrand, float(lo), float(hi))[0]
                       for lo, hi in zip(BREAKS, BREAKS[1:]))
        assert abs(integral-float(expected)) < 1e-13
    print("Arb with the infinite product tail: cubic ceiling in (0.24857, 0.24858), "
          "semiprime ceiling in (0.53639, 0.53640)")
    print(f"left-band ceilings with the truncated product: "
          f"{float(partial)*10/99:.8f} cubic, "
          f"{float(partial)*5/6*log(4/3):.8f} semiprime")
    c2 = float(partial)
    for mass in (0.02, 0.10):
        def threshold(density):
            odds = float(A/(1-A))*exp(mass/density)
            return 2*density*(odds/(1+odds)-float(A))
        old, new = threshold(c2), threshold(5*c2/6)
        assert 0 < old < new
        print(f"semiprime mass {mass:.2f}, old/new threshold {old:.8f}/{new:.8f}")


def inverse_series(degree, two_over_pi):
    """Taylor coefficients of 1/arccos(y), for ball or float arithmetic."""
    zero = two_over_pi*0
    arcsine = [zero for _ in range(degree+1)]
    values = [two_over_pi]
    if degree >= 1:
        arcsine[1] = two_over_pi
    for j in range(3, degree+1, 2):
        arcsine[j] = arcsine[j-2]*(j-2)**2/((j-1)*j)
    for n in range(1, degree+1):
        values.append(sum((arcsine[j]*values[n-j] for j in range(1,n+1,2)), zero))
    return values


def cosine_weights(degree):
    series = np.array(inverse_series(degree, 2/pi))
    weights = np.polynomial.chebyshev.poly2cheb(series)
    assert min(weights) >= 0 and abs(sum(weights)-sum(series)) < 1e-12
    return weights


def power_corrected_weights(y, cutoff, degree, offset):
    """Taylor series of the odd-power inverse, followed by cosine conversion."""
    ctx.prec = max(180, 2*degree+80)
    ctx.cap = degree+1
    ll = arb(y).log()
    omega = arb.pi()/(arb(3)/11+arb(3)*cutoff/2)
    variable = arb_series([0, 1], prec=degree+1)
    ratio = variable.acos()/omega
    inverse = arb_series([1], prec=degree+1)
    for j in range(3, cutoff+1, 2):
        mu = int(fmpz(j).moebius_mu())
        if mu:
            inverse += mu*(-arb(j-1)*ll*ratio/2).exp()/j
    series = inverse/(ratio*(arb(offset)*ll*ratio).cos())
    assert all(series[j].rad() < arb("1e-12") for j in range(degree+1))
    return np.polynomial.chebyshev.poly2cheb([float(series[j]) for j in range(degree+1)])


def flattened_counts(y, twins, semis):
    """Normalized permutation sums for each order on the arithmetic side."""
    roots = np.unique(np.r_[twins, twins+2, semis[:,1:].ravel()])
    ll, cutoff, degree, width = log(y), 5, 1100, 0.0001
    omega = pi/(float(A)+1.5*cutoff)
    assert omega/ll > 18*width
    logs = np.log(roots)
    responses = []
    for offset in (0, 4*width, 8*width):
        coeffs = power_corrected_weights(y,cutoff,degree,offset)
        values = np.zeros(len(roots))
        for j in range(1,26,2):
            angle = j*logs
            response = np.polynomial.chebyshev.chebval(np.cos(omega*angle/ll),coeffs)
            values += np.exp(-(j-1)*logs/2)*np.cos(offset*angle)*response
        responses.append(values*logs/ll)
    for k, labels in ((2,np.column_stack((twins,twins+2))), (3,semis[:,1:])):
        positions = np.searchsorted(roots, labels)
        estimates = np.zeros(len(labels))
        for order in permutations(range(k)):
            estimates += np.prod([responses[i][positions[:,j]] for i,j in enumerate(order)],axis=0)
        estimates /= factorial(k)
        approx = sum(estimates)
        assert round(approx)==len(labels)
        print(f"  order {k} joint sum before rounding={approx:.7f}, "
              f"integer target={len(labels)}, max row error={max(abs(estimates-1)):.8f}")


def exact_label_check():
    """Include a free background sign and check both centered conditional moments."""
    labels = (5, 7, 11, 13, 17, 19)
    signs = tuple(product((-1,1),repeat=len(labels)))
    coefficients = [[Fraction(p+i, p*(i+1)) for p in labels] for i in (1,2,3)]
    observations = [[sum(a*s for a,s in zip(row, sample)) for row in coefficients]
                    for sample in signs]
    for count in (2,3,4,5):
        for k in (2,3):
            moments = []
            for sign in (1,-1):
                group = [row for sample,row in zip(signs,observations)
                         if np.prod(sample[:count])==sign]
                means = [sum(row[i] for row in group)/len(group) for i in range(k)]
                centered = [np.prod([row[i]-means[i] for i in range(k)]) for row in group]
                moments.append(sum(centered)/len(group))
            expected = sum(np.prod([coefficients[i][j] for i,j in enumerate(order)])
                           for order in permutations(range(k))) if count==k else 0
            assert (moments[0]-moments[1])/2 == expected
    print("64 exact sign configurations, two moment orders, four label counts verified")


def kernel_check():
    rho = (2*arb.pi()/13).cos()
    degrees = (8, 16, 32, 64, 96)
    coefficients = inverse_series(max(degrees), 2/arb.pi())
    for k in degrees:
        tail = 1/rho.acos() - sum(a*rho**j for j, a in enumerate(coefficients[:k+1]))
        assert tail > 0
        if k == 64:
            assert tail < arb(1)/5000
        weights = cosine_weights(k)
        angles = np.linspace(2*pi/13, 11*pi/13, 1001)
        values = np.polynomial.chebyshev.chebval(np.cos(angles), weights)
        assert abs(max(abs(1/angles-values))-float(tail)) < 1e-12
        print(f"K={k:2}, interval maximum kernel error={float(tail):.10g}, "
              f"cosine coefficient sum={sum(weights):.6f}")


def finite_arithmetic(rows, y, degree, all_powers):
    """Sum the six permutations for each of the genuine three-factor rows."""
    ll = log(y)
    inv = cosine_weights(degree)
    xi = ll**-3
    modes = OMEGA*np.arange(degree+1)/ll
    modes[0] = 2*xi
    length = digamma(int(y**1.5/2)+1)-digamma(y)
    total = 0.0
    for start in range(0, len(rows), 128):
        block = rows[start:start+128]
        ns = block[:,0].astype(float)
        ps = block[:,1:].astype(float)
        logs = np.log(ps)

        def odd_coeff(heights):
            z = np.exp(1j*logs[...,None]*heights)
            if all_powers:
                z = z/(1-z*z/ps[...,None])
            return (logs/np.sqrt(ps))[...,None]*z.real

        left = odd_coeff(np.array([0]))[:,:,0]
        middle = odd_coeff(np.array([xi]))[:,:,0]
        right = odd_coeff(modes)
        cubic = np.zeros((len(block), len(modes)))
        for j in range(3):
            p, q = (j+1)%3, (j+2)%3
            cubic += (left[:,p]*middle[:,q]+left[:,q]*middle[:,p])[:,None]*right[:,j,:]
        read = (OMEGA*np.log(ns)/(2*ll))*(cubic@inv) - cubic[:,0]/6
        total += np.sum(read/np.log(ns))
    return total/length


def integer_check(previous):
    sizes = (8192, 32768)
    limit = int(max(sizes)**1.5/2)+2
    factors = previous.smallest_factors(limit)
    witness = previous.test_witness(limit)
    print("finite integer witness:", witness)
    for y in sizes:
        twins, semis = previous.admitted_rows(y, factors, witness)
        ns = semis[:,0].astype(float)
        length = digamma(int(y**1.5/2)+1)-digamma(y)
        mass = sum(np.log(ns)*np.log(ns+2)/ns)/length
        values = [finite_arithmetic(semis, y, k, False) for k in (8, 24, 64)]
        odd_read = finite_arithmetic(semis, y, 64, True)
        assert abs(values[-1]-mass) < 0.0002
        assert abs(values[-1]-mass) < abs(values[0]-mass)
        assert abs(odd_read-mass) < 0.01
        print(f"Y={y}: {len(semis)} rows, true semiprime mass={mass:.8f}, "
              f"first-power reads={[round(float(v),8) for v in values]}, all odd={odd_read:.8f}")
        flattened_counts(y,twins,semis)


def zero_check(previous):
    ds = previous.primes_up_to(12000)
    ds = ds[(ds>800) & (ds%48 == 23)]
    y, degrees, xi = 32, (4, 12, 24), 0.09
    moving = OMEGA*np.arange(max(degrees)+1)/log(y)
    moving[0] = 2*xi
    simple_heights = np.r_[0, xi, moving]
    cutoff, degrees_full = 3, (24, 64)
    omega = pi/(float(A)+1.5*cutoff)
    centers = omega*np.arange(max(degrees_full)+1)/log(y)
    offsets = (0,0.008,0.016)
    shifted = [np.abs(centers+sign*d) for d in offsets for sign in (-1,1)]
    heights = np.r_[simple_heights, np.concatenate(shifted)]
    assert max(heights)+previous.RADIUS < previous.TOP
    zero_rows = []
    nroots = 0
    for d in ds:
        roots = previous.QuadraticTheta(int(d)).zeros()
        nroots += len(roots)
        displacements = roots[:,None]-heights
        kernel = previous.zero_kernel(displacements)*(abs(displacements)<=previous.RADIUS)
        zero_rows.append(np.sum(kernel,axis=0)-log(d))
    zero_side = -0.5*previous.centered(np.array(zero_rows))
    ps = previous.primes_up_to(int(np.exp(previous.B+4*previous.BOX)))
    ps = ps[ps>3]
    signs = np.column_stack([previous.character_column(ds,int(p)) for p in ps])
    coefs = np.vstack([previous.prime_coefficient(int(p),heights) for p in ps])
    prime_side = previous.centered(signs@coefs)
    errors, full_errors = [], []
    linear_forms = []
    for k in degrees_full:
        matrix = np.zeros((len(heights),len(offsets)))
        for i, offset in enumerate(offsets):
            c = power_corrected_weights(y,cutoff,k,offset)
            for sign_index in (0,1):
                start = len(simple_heights)+(2*i+sign_index)*len(centers)
                matrix[start:start+k+1,i] += c/2
        linear_forms.append((zero_side@matrix,prime_side@matrix))
    for n in (65, 71, 77, 83, 89):
        # Direct Euler criterion. No factorization labels enter either zero statistic.
        marks = np.array([1 if pow(n*(n+2),(int(d)-1)//2,int(d))==1 else -1 for d in ds])
        def cubic(data):
            means = []
            for sign in (1, -1):
                group = previous.centered(data[marks==sign])
                means.append(np.mean((group[:,0]*group[:,1])[:,None]*group[:,2:],axis=0))
            return (means[0]-means[1])/(2*log(n))
        differences = cubic(zero_side)-cubic(prime_side)
        for k in degrees:
            weights = OMEGA*log(n)/(2*log(y))*cosine_weights(k)
            weights[0] -= 1/6
            errors.append(abs(weights@differences[:k+1]))
        for zm,pm in linear_forms:
            for rank in (2,3):
                def contrast(data):
                    moments = []
                    for sign in (1,-1):
                        group = previous.centered(data[marks==sign,:rank])
                        moments.append(np.mean(np.prod(group,axis=1)))
                    return (moments[0]-moments[1])/2
                scale = np.sqrt(n*(n+2))/(factorial(rank)*log(y)**rank)
                full_errors.append(scale*abs(contrast(zm)-contrast(pm)))
    assert len(errors)==15 and max(errors)<0.001
    assert len(full_errors)==20 and max(full_errors)<0.005, (len(full_errors),max(full_errors))
    print(f"{len(ds)} conductors, {nroots} numerical zeros, "
          f"15 inverse statistics: maximum zero/prime difference={max(errors):.9f}")
    print(f"20 power-corrected covariance/cubic comparisons: max difference={max(full_errors):.9f}")


if __name__ == "__main__":
    previous = article("gap-two-quadratic-conductor-covariance")
    sieve_constants(previous.primes_up_to(1_000_000))
    exact_label_check()
    kernel_check()
    integer_check(article("gap-two-coherent-zero-heights"))
    zero_check(previous)
