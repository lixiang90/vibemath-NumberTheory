"""Replicate the fixed-gap multiscale quadratic-variation checks.

From the zeros-exact-distribution directory run
    python -B gap-two-window-quadratic-variation/verify.py --arb
The contour evaluator and auxiliary-prime families are shared with the
preceding articles. Dependencies are numpy, scipy and python-flint.
"""

import argparse
import itertools
import math
import runpy
from fractions import Fraction
from pathlib import Path

import numpy as np
from flint import fmpz
from scipy.special import gammainc


ROOT = Path(__file__).resolve().parents[1]
PRIOR = runpy.run_path(str(ROOT / "gap-two-chen-auxiliary-switch" / "verify.py"))
SIGNS = PRIOR["SIGNS"]


def weight(p):
    return np.log(p) * np.log(p + 2) / p


def product_variance(x):
    cross = np.column_stack(
        [x[:, i] * x[:, j] for i, j in itertools.combinations(range(4), 2)]
    )
    result = np.mean(cross**2, axis=1) - np.prod(x, axis=1) / 3
    assert np.all(result >= -1e-14)
    return result


def check_reused_prime():
    # The earlier endpoint 53 is the larger factor of the later neighbour 1219.
    rows = ((53, 5, 11), (1217, 23, 53))
    labels = sorted({v for row in rows for v in row})
    assert all(fmpz(v).is_prime() for v in labels)
    assert all(p + 2 == r * s for p, r, s in rows)
    signs = [
        dict(zip(labels, v)) for v in itertools.product((-1, 1), repeat=len(labels))
    ]
    doubled = np.zeros((len(signs) ** 2, 2, 4), dtype=np.int64)
    for i, (q, d) in enumerate(itertools.product(signs, repeat=2)):
        for j, (p, r, s) in enumerate(rows):
            if q[p] == 1 and q[r] == q[s] and d[r] != d[s]:
                cell = SIGNS.index((q[r], d[r]))
                doubled[i, j, cell] = -(j + 2) * d[p]
    n = len(doubled)
    for cell, (sigma_q, _) in enumerate(SIGNS):
        early, late = doubled[:, 0, cell], doubled[:, 1, cell]
        assert 32 * int(early @ early) == n * 4
        assert 32 * int(late @ late) == n * 9
        assert int(early @ late) == int((early**2) @ late) == 0
        # Relative to the product of expectations, the exact square covariance
        # effect is 0 or 2 for the mixed moment, depending on the shared Q sign.
        mixed = int((early**2) @ (late**2))
        assert 1024 * mixed == n * 4 * 9 * (1 + sigma_q)
    print(f"Cross-scale endpoint/factor collision checked on {n} exact assignments.")


def check_square_statistics():
    # In the Gaussian comparison H/b is the mean of four squared normals,
    # hence Gamma(shape=2, scale=1/2). Check the constants in rational arithmetic.
    moments = [Fraction(math.factorial(k + 1), 2**k) for k in range(5)]
    assert moments[1] == 1 and moments[2] - moments[1] ** 2 == Fraction(1, 2)
    assert moments[2] / 3 == Fraction(1, 2)
    assert (moments[4] - moments[2] ** 2) / 9 == Fraction(7, 12)
    for cells in (4, 64, 1024):
        # Three disjoint deterministic windows with constant b and zero twin
        # mass. This model uses the exact gamma CDF, with no random simulation.
        one = gammainc(2 * cells, 2 * cells)
        maximum = 1 - (1 - one) ** 3
        print(
            f"Gaussian squares, {cells:4} cells: one sign {one:.6f}, best of three {maximum:.6f}"
        )
    assert abs(maximum - 7 / 8) < 0.005


def check_common_endpoint():
    # The twin 59,61 contributes to all four coordinates. The neighbour of 53
    # contributes to a single factor mark. Give the two rows weights 4 and 9.
    labels = (59, 61, 53, 5, 11)
    signs = [
        dict(zip(labels, v)) for v in itertools.product((-1, 1), repeat=len(labels))
    ]
    doubled = np.zeros((len(signs) ** 2, 4), dtype=np.int64)
    for i, (q, d) in enumerate(itertools.product(signs, repeat=2)):
        if q[59] == q[61] == 1 and d[61] == -1:
            doubled[i] -= 4 * d[59]
        if q[53] == 1 and q[5] == q[11] and d[5] != d[11]:
            doubled[i, SIGNS.index((q[5], d[5]))] -= 3 * d[53]
    n = len(doubled)
    assert np.all(doubled.sum(axis=0) == 0)
    gram = doubled.T @ doubled
    assert np.array_equal(
        128 * gram,
        n * (256 * np.ones((4, 4), dtype=np.int64) + 36 * np.eye(4, dtype=np.int64)),
    )
    print(
        f"Common twin contribution to all four covariances checked on {n} exact assignments."
    )
    product_variance(doubled / 2)
    for d, b in ((0, 1), (1, 0), (2, 3)):

        def wick(ids):
            if not ids:
                return 1
            return sum(
                (d + b * (ids[0] == ids[i])) * wick(ids[1:i] + ids[i + 1 :])
                for i in range(1, len(ids))
            )

        second = wick((0, 0, 3, 3))
        variance = (d + b) ** 2 + d**2
        assert second - d * d == variance
        assert second - Fraction(wick((0, 1, 2, 3)), 3) == variance


def cell_ids(p, parts, first, dyads):
    boundaries = first * 2 ** np.arange(dyads + 1)
    parent = np.searchsorted(boundaries[1:], p, side="left")
    lo = boundaries[parent]
    assert np.all((lo < p) & (p <= 2 * lo))
    return parent * parts + (parts * (p - lo) - 1) // lo


def real_auxiliaries():
    # Y=4096 gives exactly four full dyadic blocks in (2Y,Y^(3/2)/2].
    first, dyads = 8192, 4
    last = first * 2**dyads
    primes, rates = PRIOR["prime_weights"](last + 2)
    ends, factors = [], []
    for i in range(dyads):
        te, bf = PRIOR["arithmetic_rows"](first * 2**i, primes, rates)
        ends.extend(te)
        factors.extend(bf)
    ep = np.array([row[0] for row in ends])
    fp = np.array([row[0] for row in factors])
    fd = np.array([row[2] for row in factors])
    tw = np.array([row[3] for row in ends])
    sem = np.array([row[4] for row in factors])
    ew = weight(ep) * np.array([row[2] for row in ends])
    gf = np.array([row[3] for row in factors])
    fw = weight(fp)
    assert ep[tw].size == sum(
        primes[p + 2] for p in range(first + 1, last + 1) if primes[p]
    )
    print(
        f"Intervals ({first},{last}]: {sum(tw)} twins, {sum(sem)} internal semiprimes."
    )
    print(
        f"Physical masses T/16 = {weight(ep[tw]).sum() / 16:.6f}, B/128 = {fw[sem].sum() / 128:.6f}"
    )
    levels = (1, 4, 16, 64)
    bins = {k: cell_ids(fp, k, first, dyads) for k in levels}
    end_bins = {k: cell_ids(ep, k, first, dyads) for k in levels}
    moments = {k: [] for k in levels}
    quartets = {k: [] for k in levels}
    candidates = []
    paired = []
    final_dyad = ep > last // 2
    end_coefficients = ew / np.sqrt(weight(ep))
    labels, inverse = np.unique(
        np.concatenate((ep, ep + 2, fp, fp + 2, fd)), return_inverse=True
    )
    e, ex, f, fx, d = np.split(inverse, np.cumsum([len(ep), len(ep), len(fp), len(fp)]))
    for q, family in PRIOR["auxiliary_families"]():
        sq = np.array([PRIOR["legendre"](v, q) for v in labels])
        eg, fg = (sq[e] == 1) & (sq[ex] == 1), (sq[f] == 1) & (sq[fx] == 1)
        for aux in family:
            sd = np.array([PRIOR["legendre"](v, aux) for v in labels])
            em = eg & (sd[e] == -1) & (sd[ex] == -1)
            fm = fg & (sd[fx] == -1)
            cells = 2 * (sq[d] == 1).astype(int) + (sd[d] == 1)
            candidates.append(ew[em].sum() + (fw * gf)[fm & (sd[f] == -1)].sum() / 4)
            amplitudes = -np.sqrt(fw) * sd[f] * gf * fm
            end_amplitudes = -end_coefficients * sd[e] * eg * (sd[ex] == -1)
            for k in levels:
                z = np.bincount(
                    4 * bins[k] + cells, weights=amplitudes, minlength=4 * k * dyads
                ).reshape(-1, 4)
                h = np.mean(z * z, axis=1)
                moments[k].append((h.sum(), np.dot(h, h) / 3))
                common = np.bincount(
                    end_bins[k], weights=end_amplitudes, minlength=k * dyads
                )
                x = z + common[:, None]
                quartets[k].append(
                    (np.dot(x[:, 0], x[:, 3]), product_variance(x).sum())
                )
                if k == 1:
                    paired.append(x[-1, [0, 3]])
    print(" K        E H       sd(H)     model s     rms s_hat       E(A-H)")
    for k in levels:
        bs = np.bincount(bins[k][sem], weights=fw[sem] / 128, minlength=k * dyads)
        h, vhat = np.asarray(moments[k]).T
        print(
            f"{k:2d} {h.mean():11.6f} {h.std():11.6f} {np.linalg.norm(bs) / math.sqrt(2):11.6f}"
            f" {math.sqrt(vhat.mean()):13.6f} {np.mean(candidates) - h.mean():13.6f}"
        )
    twin_variance = weight(ep[tw & final_dyad]).sum() / 8
    semi_variance = fw[sem & (fp > last // 2)].sum() / 128
    rho = twin_variance / (twin_variance + semi_variance)
    first_mark, second_mark = np.array(paired).T
    print(
        f"Last dyad, two different marks: covariance correlation {np.corrcoef(first_mark, second_mark)[0, 1]:.6f} vs {rho:.6f};"
        f" same-sign frequency {np.mean(first_mark * second_mark > 0):.6f} vs {0.5 + math.asin(rho) / math.pi:.6f}"
    )
    print(" K     E product       sd(product)     model s      rms feasible s")
    for k in levels:
        ds = np.bincount(
            end_bins[k][tw], weights=weight(ep[tw]) / 8, minlength=k * dyads
        )
        bs = np.bincount(bins[k][sem], weights=fw[sem] / 128, minlength=k * dyads)
        products, estimates = np.asarray(quartets[k]).T
        scale = math.sqrt(np.sum((ds + bs) ** 2 + ds**2))
        print(
            f"{k:2d} {products.mean():13.6f} {products.std():17.6f} {scale:11.6f} {math.sqrt(estimates.mean()):17.6f}"
        )


def contour_mesh(arb, heights):
    q, d = 1_000_003, 47_880_143_639
    primes, rates = PRIOR["prime_weights"](32 * 91 + 2)
    sampled = [p for p in range(51, 91) if primes[p] and p % 6 == 5]
    sums = np.zeros((2, 2, 4))
    common = np.zeros((2, 2))
    left = np.zeros(2)
    maximum_leak = 0
    for p in sampled:
        x, logx = p + 2, math.log(p + 2)
        band = PRIOR["band_weights"](x)
        scales = {x, *band, *(x / f for f in band)}
        obs = []
        for disc, h in zip((1, 8), heights):
            family = PRIOR["BandCalculation"](p, disc, 2, scales)
            if arb:
                family.check_arb(h)
            trace, coarse = family.contour(h)
            assert np.max(np.abs(trace - coarse)) < 3e-5
            assert np.max(np.abs(trace - family.arithmetic(h, rates))) < 3e-5
            obs.append((family, trace))
        (one, u), (eight, v) = obs

        def product(c, y, b, z):
            return (
                one.projection(u, c, y) * eight.projection(v, b, z).conjugate()
            ).real

        om = PRIOR["omega"]
        endpoint = np.array(
            (om(x, 8) * product(2, x, 2, x) / logx**2, (rates[x] / logx) ** 2)
        )
        factors = np.zeros((2, 4))
        leak = 0
        for f, value in band.items():
            k = SIGNS.index((PRIOR["legendre"](f, q), PRIOR["legendre"](f, d)))
            actual = (
                value
                * om(x, 8)
                * om(f, 8)
                * product(f, f, 2 * pow(f, -1, p) % p, x / f)
                / logx**2
            )
            exact = value * rates[f] * rates[x // f] / logx**2 if x % f == 0 else 0
            factors[:, k] += (actual, exact)
            leak += abs(actual - exact)
        maximum_leak = max(maximum_leak, leak)
        gate = (
            PRIOR["legendre"](p, q) == PRIOR["legendre"](x, q) == 1
            and PRIOR["legendre"](x, d) == -1
        )
        sign = PRIOR["legendre"](p, d)
        sums[:, int(p > 70), :] -= math.sqrt(weight(p)) * gate * sign * factors
        common[:, int(p > 70)] -= math.sqrt(weight(p)) * gate * sign * endpoint
        left += weight(p) * gate * (sign == -1) * (endpoint + factors.mean(axis=1))
    # Parseval on each mark's two physical cells, computed from the contour.
    h = np.sum(sums * sums, axis=(1, 2)) / 4
    freq = np.fft.fft(sums[0], axis=0)
    assert np.isclose(h[0], np.sum(abs(freq) ** 2) / 8)
    assert abs(h[0] - h[1]) < 3e-5 and abs(left[0] - left[1]) < 3e-5
    print(f"Heights {heights}, two physical cells, five prime endpoints.")
    print(f"  max sum of absolute divisor errors = {maximum_leak:.3e}")
    print(
        f"  H(contour) = {h[0]:.9f}, H(arithmetic) = {h[1]:.9f}, absolute error = {abs(h[0] - h[1]):.3e}"
    )
    print(
        f"  A-H(contour) = {left[0] - h[0]:.9f}, A-H(arithmetic) = {left[1] - h[1]:.9f}"
    )
    return sums[0], left[0] - h[0], sums + common[:, :, None]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arb", action="store_true")
    args = parser.parse_args()
    check_reused_prime()
    check_square_statistics()
    check_common_endpoint()
    real_auxiliaries()
    heights = (320, 960, 1920, 2880)
    for p in (53, 59, 71, 83, 89):
        radii = [PRIOR["OLD"]["band_radius"](p, h, 2) for h in heights]
        assert all(
            heights[i] + radii[i] < heights[i + 1] - radii[i + 1] for i in range(3)
        )
    first, r_first, x_first = contour_mesh(args.arb, heights[:2])
    second, r_second, x_second = contour_mesh(args.arb, heights[2:])
    shift_error = np.max(abs(first - second))
    assert shift_error < 3e-5 and abs(r_first - r_second) < 3e-5
    assert np.max(abs(x_first - x_second)) < 3e-5
    print(
        f"Four disjoint height windows: max change of Z = {shift_error:.3e}, change of A-H = {abs(r_first - r_second):.3e}"
    )
    print(f"  max change of common-endpoint X = {np.max(abs(x_first - x_second)):.3e}")
    paired_products = np.sum(x_first[:, :, 0] * x_second[:, :, 3], axis=1)
    assert abs(paired_products[0] - paired_products[1]) < 3e-5
    print(
        f"  four-window product = {paired_products[0]:.9f}, arithmetic = {paired_products[1]:.9f}, absolute error = {abs(paired_products[0] - paired_products[1]):.3e}"
    )


if __name__ == "__main__":
    main()
