"""Check the common-height formulas, the squareclass variance and actual zeros.

Run with Python, NumPy, SciPy and python-flint from the repository checkout.
The zero experiment reuses the theta computation in the preceding article.
"""

from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import comb, log, pi, sqrt
from pathlib import Path

import numpy as np
from flint import fmpz
from scipy.special import digamma
from scipy.stats import irwinhall


def gate_check():
    def gates(factor_rows):
        labels = sorted({p for row in factor_rows for endpoint in row for p in endpoint})
        values = np.array(list(product(range(4), repeat=len(labels))), dtype=np.int8)
        result = []
        for row in factor_rows:
            colors = [np.bitwise_xor.reduce(values[:, [labels.index(p) for p in end]], axis=1)
                      for end in row]
            result.append((colors[0] == 1) & (colors[1] == 1))
        return np.array(result, dtype=float)

    # The two semiprimes have a common factor 5.
    indicators = gates([((53,), (5, 11)), ((83,), (5, 17))])
    assert np.array_equal(indicators.mean(axis=1), np.full(2, 1 / 16))
    share = np.mean(indicators[0] * indicators[1])
    assert share == 1 / 256
    coefficients = np.array([2, -3])
    assert np.var(coefficients @ indicators) == 15 / 256 * (coefficients @ coefficients)
    # Across scales, the second left endpoint equals the first product.
    crossing = gates([((5,), (7,)), ((5, 7), (37,))])
    assert not np.any(crossing[0] * crossing[1])
    print(f"two-row probability: shared factor {share}, reused product 0")


def previous_zero_code():
    path = (Path(__file__).resolve().parents[1]
            / "gap-two-quadratic-conductor-covariance" / "verify.py")
    spec = spec_from_file_location("previous_quadratic", path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def zero_check():
    prev = previous_zero_code()
    ds = prev.primes_up_to(4200)
    ds = ds[(ds > 800) & (ds % 48 == 23)]
    radius, band, box = 3.45, 4.4, 1.05
    hs = (0.0, 7.01, prev.HEIGHT)
    assert min(np.diff(hs)) > 2 * radius and hs[-1] + radius < prev.TOP
    data = []
    count = 0
    for d in ds:
        roots = prev.QuadraticTheta(int(d)).zeros()
        count += len(roots)
        def trace(h):
            t = roots[np.abs(roots - h) <= radius] - h
            return np.sum(2*band*np.sinc(band*t/pi)*np.sinc(box*t/pi)**4)
        data.append([trace(h) for h in hs])
    zeros = -0.5 * prev.centered(np.array(data) - np.log(ds[:, None]))
    limit = int(np.exp(band + 4*box))
    assert np.min(ds)**2 > limit
    primes = prev.primes_up_to(limit)
    primes = primes[primes > 3]
    logs = np.log(primes)
    signs = np.column_stack([prev.character_column(ds, int(p)) for p in primes])
    coefs = np.zeros((len(primes), len(hs)))
    for m in range(1, int(log(limit)/logs[0]) + 1, 2):
        u = m * logs
        g = irwinhall.cdf((u+band)/(2*box) + 2, 4) - irwinhall.cdf((u-band)/(2*box) + 2, 4)
        coefs += (logs*np.exp(-u/2)*g)[:, None] * np.cos(np.outer(u, hs))
    sums = prev.centered(signs @ coefs)

    def measured(points, signs, n):
        by_sign = []
        for value in (1, -1):
            x, middle, y = prev.centered(points[signs == value]).T
            by_sign.append([np.mean(x*y), np.mean(x*middle*y) / log(n)])
        diff = (np.array(by_sign[0]) - by_sign[1]) / 2
        return np.cos(hs[-1] * log(n)) * diff * [1, 2]

    for n in (17, 53):
        mark = prev.row_mark(ds, prev.TAGS[n])
        actual = measured(zeros, mark, n)
        arithmetic = measured(sums, mark, n)
        assert np.max(np.abs(actual - arithmetic)) < 0.006
        print(f"actual zeros row {n}: {actual}, same-family arithmetic {arithmetic}")
    print(f"found {count} zeros across {len(ds)} conductors; three disjoint windows")


def smallest_factors(limit):
    factors = np.zeros(limit + 1, dtype=np.int32)
    for p in range(2, int(sqrt(limit)) + 1):
        if factors[p] == 0:
            tail = factors[p*p::p]
            tail[tail == 0] = p
    return factors


def test_witness(max_endpoint):
    # Both primes exceed every endpoint. The original congruences hold for 3,5,7.
    q = max_endpoint + 1
    q += (3 - q) % 8
    while not fmpz(q).is_prime():
        q += 8
    product_modulus = 3 * 5 * 7 * q
    d = (-1 + product_modulus * ((8 * pow(product_modulus, -1, 16)) % 16))
    while d <= q or not fmpz(d).is_prime():
        d += 16 * product_modulus
    assert q % 8 == 3 and d % 16 == 7 and (d + 1) % product_modulus == 0
    return q, d


def admitted_rows(y, factors, witness):
    end = int(y**1.5 / 2)
    left = np.arange(2*y + 1 + (5 - (2*y+1)) % 6, end + 1, 6, dtype=np.int64)
    right = left + 2
    twin = (factors[left] == 0) & (factors[right] == 0)
    twins = left[twin]
    semi_rows = []
    for prime_end, other in ((left, right), (right, left)):
        s = factors[other]
        possible = (factors[prime_end] == 0) & (s > left**(3/11))
        n, p, s, t = left[possible], prime_end[possible], s[possible], other[possible] // s[possible]
        good = (t != s) & (factors[t] == 0)
        semi_rows.append(np.array([n[good], p[good], s[good], t[good]]).T)
    semis = np.vstack(semi_rows)
    if witness is not None:
        q, d = map(fmpz, witness)
        def keep(ns):
            return np.fromiter((fmpz(int(n)).jacobi(q) == fmpz(int(n+2)).jacobi(q) == 1
                                and fmpz(int(n)).jacobi(d) == fmpz(int(n+2)).jacobi(d) == -1
                                for n in ns), dtype=bool, count=len(ns))
        twins = twins[keep(twins)]
        semis = semis[keep(semis[:, 0])]
    return twins, semis


def actual_gate_sample(factors):
    twins, semis = admitted_rows(128, factors, None)
    ns = np.r_[twins, semis[:, 0]]
    primes = [q for q in range(1027, 4097, 8) if fmpz(q).is_prime()]
    upper = 60_000_000
    moments = []
    samples = 0
    for q in primes:
        modulus = 3*q
        residue = -1 + modulus * ((8*pow(modulus, -1, 16)) % 16)
        step = 16*modulus
        d = residue + max(0, (int(sqrt(upper)) - residue)//step + 1)*step
        candidates = [v for v in range(d, upper+1, step) if fmpz(v).is_prime()]
        passed_q = [int(n) for n in ns if fmpz(int(n)).jacobi(q) == fmpz(int(n+2)).jacobi(q) == 1]
        counts = np.array([sum(fmpz(n).jacobi(v) == fmpz(n+2).jacobi(v) == -1
                               for n in passed_q) for v in candidates])
        assert len(counts)
        moments.append([np.mean(counts), np.mean(counts**2)])
        samples += len(counts)
    first, second = np.mean(moments, axis=0)
    print(f"actual two-stage screen: n={len(ns)}, q={len(primes)}, pairs={samples}, "
          f"mean/variance {first:.6f}/{second-first**2:.6f}, "
          f"independent {len(ns)/16:.6f}/{15*len(ns)/256:.6f}")


def coefficient(primes, logs, h):
    # Physical odd powers through the ninth, for an independent local arithmetic test.
    e = np.exp(1j * h * logs)
    ratio = e**2 / primes
    return logs / np.sqrt(primes) * (e * (1 - ratio**5) / (1 - ratio)).real


def height_check(y, factors, witness):
    twins, semis = admitted_rows(y, factors, witness)
    assert len(twins) and len(semis)
    ll = log(y)
    length = digamma(int(y**1.5/2) + 1) - digamma(y)
    tl = np.log(twins)
    tps = np.array([twins, twins + 2], dtype=float)
    tlogs = np.log(tps)
    n = semis[:, 0].astype(float)
    sl = np.log(n)
    sps = semis[:, 1:].T.astype(float)
    slogs = np.log(sps)
    theta = slogs[1] / (slogs[1] + slogs[2])
    jt = tl*np.log(twins+2) / (twins*length)
    sw = sl*np.log(n+2) / (n*length)
    bm = 2*theta*(1-theta)*sw
    layers = {1.0: [], 0.4: [], 0.0: []}
    max_shape_error = np.zeros(2)
    norm_t, norm_s = coefficient(tps, tlogs, 0), coefficient(sps, slogs, 0)
    middle = coefficient(sps, slogs, 1/ll**3)
    last = 96
    near_v, near_b = [], []
    for k in range(1, last + 1):
        h = 2*pi*k / (8*ll)
        at, ass = coefficient(tps, tlogs, h), coefficient(sps, slogs, h)
        pair = (norm_t[0]*at[1] + at[0]*norm_t[1]) / length
        triple = sum(ass[j] * (norm_s[(j+1)%3]*middle[(j+2)%3]
                               + middle[(j+1)%3]*norm_s[(j+2)%3]) for j in range(3)) / (sl*length)
        for lam, entries in layers.items():
            v = np.cos(lam*h*tl) @ pair
            b = (2*np.cos(lam*h*sl)) @ triple
            entries.append([v,b])
            # Each prime contributes its sum and difference of logarithmic phases.
            v_model = jt @ (np.cos(h*(lam-1)*tl) + np.cos(h*(lam+1)*tl))
            b_model = np.zeros_like(n)
            for x in (sl, slogs[1], slogs[2]):
                b_model += np.cos(h*(lam*sl-x)) + np.cos(h*(lam*sl+x))
            errors = np.abs([v-v_model, b-bm@b_model])
            max_shape_error = np.maximum(max_shape_error, errors)
        near_lambda = 1 + 2/last
        near_v.append(np.cos(near_lambda*h*tl) @ pair)
        near_b.append(2*np.cos(near_lambda*h*sl) @ triple)
    center = np.array([sum(jt), sum(bm)])
    observed = np.mean(layers[1.0][last//2:], axis=0)
    unmatched = np.mean(layers[0.4][last//2:], axis=0)
    profile = np.array([jt @ np.sinc(2*2*tl/(8*ll)), bm @ np.sinc(2*2*sl/(8*ll))])
    near = np.array([np.mean(near_v), np.mean(near_b)])
    assert np.all(max_shape_error < np.array([8/(2*y), 22/(2*y)**(3/11)]) * center)
    print(f"Y={y}, gate={witness is not None}, rows={len(twins)}/{len(n)}, "
          f"centers {center}, high matched {observed}, high unmatched {unmatched}")
    print(f"detuned averages {near}, sinc prediction {profile}, "
          f"maximum finite prime-power/endpoint errors {max_shape_error}")
    energies = [np.sum((np.array(layers[lam])-target)**2, axis=0)
                for lam, target in ((0.0, 0), (0.4, 0), (1.0, center))]
    print(f"96-mode squared energies, lambda=0,0.4,1: {np.array(energies).ravel()}")
    m, k0 = 16, 64
    offsets = np.arange(-m, m+1)
    weights = np.array([comb(2*m, int(m+i))/4**m for i in offsets])
    assert sum(weights) == 1
    test_angles = np.linspace(-0.49, 0.49, 201)
    kernel = np.cos(2*pi*np.outer(test_angles, k0+offsets)) @ weights
    assert np.max(abs(kernel - np.cos(2*pi*k0*test_angles)*np.cos(pi*test_angles)**(2*m))) < 2e-14
    block = weights @ np.array(layers[1.0])[k0+offsets-1]
    low = [sum(np.array(layers[lam])[:, 1] < center[1]/2) for lam in (1.0, 0.4)]
    print(f"matched binomial block {block}, cubic below half its reference at {low} heights")


def main():
    np.set_printoptions(precision=7, suppress=True)
    gate_check()
    zero_check()
    sizes = [1024, 8192, 32768]
    limit = int(max(sizes)**1.5/2) + 2
    factors = smallest_factors(limit)
    actual_gate_sample(factors)
    witness = test_witness(limit)
    print(f"finite witness Q,d0={witness}")
    for y in sizes:
        height_check(y, factors, witness)
    height_check(sizes[1], factors, None)


if __name__ == "__main__":
    main()
