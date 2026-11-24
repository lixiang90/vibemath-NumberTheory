"""Check the auxiliary-prime switch for gap two.

Run from the field directory with numpy, scipy and python-flint:
    python -B gap-two-chen-auxiliary-switch/verify.py
Add --arb for an independent check of the contour's logarithmic derivatives.
"""

import argparse
import itertools
import math
import runpy
from pathlib import Path

import numpy as np
from flint import fmpz


ROOT = Path(__file__).resolve().parents[1]
OLD = runpy.run_path(str(ROOT / "gap-two-separated-zero-windows" / "verify.py"))
factor_weights = OLD["factor_weights"]
prime_weights = OLD["prime_weights"]
BandCalculation = OLD["BandCalculation"]
omega = OLD["omega"]
SIGNS = ((-1, -1), (-1, 1), (1, -1), (1, 1))


def band_weights(x):
    # The small factor has logarithmic exponent between 3/10 and 9/20.
    return {
        d: b for d, b in factor_weights(x).items() if d**10 >= x**3 and d**20 <= x**9
    }


def legendre(n, modulus):
    return int(fmpz(int(n)).jacobi(modulus))


def auxiliary_families():
    # Numerical samples in congruence classes. The sizes are not intended to
    # implement the asymptotic power-window and cohomology hypotheses.
    for slot in range(48):
        q = 1_000_003 + 8 * (2500 * slot)
        while not fmpz(q).is_prime():
            q += 8
        ds = []
        k = 0
        while len(ds) < 24:
            d = (8 + 16 * k) * (3 * 5 * 7 * q) - 1
            if fmpz(d).is_prime():
                assert d % 16 == 7 and (d + 1) % (105 * q) == 0
                ds.append(d)
            k += 1
        assert q % 8 == 3
        yield q, ds


def exhaustive_switch(triples):
    labels = sorted({n for row in triples for n in row})
    assignments = [
        dict(zip(labels, v)) for v in itertools.product((-1, 1), repeat=len(labels))
    ]
    count = len(assignments) ** 2
    left = np.zeros((count, len(SIGNS)), dtype=np.int64)
    right = np.zeros_like(left)
    for i, (cq, cd) in enumerate(itertools.product(assignments, repeat=2)):
        for p, r, s in triples:
            if cq[p] == 1 and cq[r] == cq[s] and cd[r] != cd[s]:
                cell = SIGNS.index((cq[r], cd[r]))
                side = left if cd[p] == -1 else right
                side[i, cell] += 1
    diff, active = left - right, left + right
    assert np.all(64 * left.sum(axis=0) == len(triples) * count)
    assert np.all(diff.sum(axis=0) == 0)
    # Each actual factor weight is 1/2. Test the diagonal and all mixed terms.
    gram = diff.T @ diff
    assert np.array_equal(
        128 * gram, 4 * len(triples) * count * np.eye(len(SIGNS), dtype=np.int64)
    )
    assert np.array_equal(
        (diff**4).sum(axis=0), (3 * active**2 - 2 * active).sum(axis=0)
    )
    if len(triples) == 3:
        assert np.any(np.count_nonzero(diff, axis=1) >= 2)
    print(
        f"All four sign cells, their 1/128 diagonal and zero cross moments checked on {count} exact assignments."
    )


def incidence_check():
    # A Page exception links the squarefree labels u and ker(u*e) by an involution.
    labels = [5, 7, 11, 35, 55, 77, 385]
    coefs = {u: j - 3 for j, u in enumerate(labels)}
    for exceptional in labels:
        total = 0
        for u, c in coefs.items():
            partner = u * exceptional // math.gcd(u, exceptional) ** 2
            assert partner * exceptional // math.gcd(partner, exceptional) ** 2 == u
            total += abs(c * coefs.get(partner, 0))
        assert total <= sum(c * c for c in coefs.values())


def arithmetic_rows(start, primes, rates):
    ps = [
        int(n)
        for n in np.flatnonzero(primes[start + 1 : 2 * start + 1]) + start + 1
        if n % 6 == 5
    ]
    endpoints, factors = [], []
    for p in ps:
        x = p + 2
        if rates[x]:
            endpoints.append((p, x, (rates[x] / math.log(x)) ** 2, primes[x]))
        for d, b in band_weights(x).items():
            if x % d or not rates[d] or not rates[x // d]:
                continue
            value = b * rates[d] * rates[x // d] / math.log(x) ** 2
            factors.append((p, x, d, value, primes[d] and primes[x // d]))
    return endpoints, factors


def sample_block(start, families, primes, rates):
    ends, rows = arithmetic_rows(start, primes, rates)
    twins = sum(int(e[3]) for e in ends)
    semis = sum(int(r[4]) for r in rows)
    # One character evaluation per distinct integer, then use array indexing.
    vals = sorted(
        {n for p, x, _, _ in ends for n in (p, x)}
        | {n for p, x, d, _, _ in rows for n in (p, x, d)}
    )
    index = {v: i for i, v in enumerate(vals)}
    ep, ex = ([index[e[j]] for e in ends] for j in (0, 1))
    fp, fx, fd = ([index[r[j]] for r in rows] for j in (0, 1, 2))
    endpoint_weights = np.array([e[2] for e in ends])
    factor_weights_array = np.array([r[3] for r in rows])
    is_semiprime = np.array([r[4] for r in rows])
    means, raw_moments = [], []
    standardized = []
    fixed_left, fixed_switch = [], []
    chosen_cell = SIGNS.index((-1, 1))
    for q, ds in families:
        sq = np.array([legendre(n, q) for n in vals], dtype=np.int8)
        eligible_end = (sq[ep] == 1) & (sq[ex] == 1)
        eligible_factor = (sq[fp] == 1) & (sq[fx] == 1)
        local_means, local_moments = [], []
        for d_aux in ds:
            sd = np.array([legendre(n, d_aux) for n in vals], dtype=np.int8)
            good_end = eligible_end & (sd[ep] == -1) & (sd[ex] == -1)
            main = float(endpoint_weights[good_end].sum())
            valid = eligible_factor & (sd[fx] == -1)
            cell = 2 * (sq[fd] == 1).astype(int) + (sd[fd] == 1)

            def sum_cells(mask):
                return np.bincount(
                    cell[mask], weights=factor_weights_array[mask], minlength=4
                )

            minus = sum_cells(valid & (sd[fp] == -1))
            plus = sum_cells(valid & (sd[fp] == 1))
            local_means.append(main + minus.mean())
            local_moments.append(np.mean((minus - plus) ** 2))
            fixed_left.append(main + minus[chosen_cell])
            fixed_switch.append(minus[chosen_cell] - plus[chosen_cell])
            chosen = valid & is_semiprime & (cell == chosen_cell)
            switched = -0.5 * int(sd[fp][chosen].sum())
            standardized.append(switched / math.sqrt(semis / 128))
        means.append(np.mean(local_means))
        raw_moments.append(np.mean(local_moments))
    # Here row weights are 1. Restore mu_U explicitly in the mean/second-moment formula.
    # It cancels when reporting the inferred counts 128*M2 and 16*(M1-M2).
    mean, moment = np.mean(means), np.mean(raw_moments)
    print(
        f"{start:7} {twins:6} {semis:8} {128 * moment:12.3f} {16 * (mean - moment):12.3f}"
    )
    samples = np.array(standardized)
    print(
        f"  fixed (-,+) semiprime cell: E W^2 = {np.mean(samples**2):.3f},"
        f" E W^4 = {np.mean(samples**4):.3f}, fraction at zero = {np.mean(samples == 0):.3f}"
    )
    # In unit row weights the factors of mu cancel from the calibration in (14).
    calibrated = np.asarray(fixed_switch) / math.sqrt(np.mean(fixed_left))
    variance_ratio = semis / (8 * twins + semis)
    model_mass = (
        math.erf(0.25 / math.sqrt(2 * variance_ratio)) if variance_ratio else 1.0
    )
    print(
        f"  fixed (-,+), full factor sum: fraction |V| <= 1/4 ="
        f" {np.mean(np.abs(calibrated) <= 0.25):.3f}, Gaussian comparison = {model_mass:.3f}"
    )


def contour_check(families, arb):
    # Pick one actual pair with opposite last endpoint signs, a common smaller
    # factor and all the other endpoint conditions in force for p=53 and 83.
    q_aux, d_aux = next(
        (q, d)
        for q, ds in families
        if all(legendre(n, q) == 1 for n in (53, 55, 83, 85))
        for d in ds
        if [legendre(n, d) for n in (53, 55, 83, 85)]
        in ([-1, -1, 1, -1], [1, -1, -1, -1])
    )
    print(f"Contours at 320 and 960, auxiliary Q = {q_aux}, D = {d_aux}")
    primes, rates = prime_weights(32 * 85 + 2)
    for p in (23, 29, 53, 83):
        x = p + 2
        order = math.ceil(4 * math.log(math.log(p)) / math.log(4))
        weights = band_weights(x)
        scales = {x, *weights, *(x / d for d in weights)}
        traces = []
        for disc, height in ((1, 320), (8, 960)):
            fam = BandCalculation(p, disc, order, scales)
            if arb:
                fam.check_arb(height)
            val, coarse = fam.contour(height)
            assert np.max(np.abs(val - coarse)) < 2e-5
            assert np.max(np.abs(val - fam.arithmetic(height, rates))) < 3e-5
            traces.append((fam, val))
        (one, first), (eight, second) = traces

        def pair(a, y, b, z):
            return (
                one.projection(first, a, y) * eight.projection(second, b, z).conjugate()
            )

        main = omega(x, 8) * pair(2, x, 2, x) / math.log(x) ** 2
        observed = np.zeros(4, dtype=complex)
        expected = np.zeros(4)
        for d, b in weights.items():
            i = SIGNS.index((legendre(d, q_aux), legendre(d, d_aux)))
            observed[i] += (
                omega(x, 8)
                * omega(d, 8)
                * b
                * pair(d, d, 2 * pow(d, -1, p), x / d)
                / math.log(x) ** 2
            )
            if x % d == 0:
                expected[i] += b * rates[d] * rates[x // d] / math.log(x) ** 2
        discrepancy = max(
            abs(main - (rates[x] / math.log(x)) ** 2),
            np.max(np.abs(observed - expected)),
        )
        assert discrepancy < 1e-5
        print(
            f"p={p:2} endpoint={main.real:.7f} factor cells={observed.real.round(7)} max error={discrepancy:.2e}"
        )
    assert primes[53] and primes[83]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arb", action="store_true")
    args = parser.parse_args()
    families = list(auxiliary_families())
    shared_factor = ((53, 5, 11), (83, 5, 17))
    exhaustive_switch(shared_factor)
    exhaustive_switch((*shared_factor, (89, 7, 13)))
    incidence_check()
    primes, rates = prime_weights(200_002)
    print("      U  twins    semis       128*M2   16*(M1-M2)")
    for start in (1000, 10000, 100000):
        sample_block(start, families, primes, rates)
    contour_check(families, args.arb)


if __name__ == "__main__":
    main()
