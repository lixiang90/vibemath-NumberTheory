"""Translated quadratic characters on one prime modulus.

Run python verify.py for all checks, or python verify.py --arithmetic-only
for the translation and arc checks. The zero search reuses the sibling
gap-two-auxiliary-conductor-bands/verify.py and requires python-flint.
"""

import argparse
import math
import runpy
from pathlib import Path

import numpy as np
from scipy.stats import irwinhall


Q = 107
D = 449399
ALPHA = math.sqrt(2) - 1


def primes_through(limit):
    prime = np.ones(limit + 1, dtype=bool)
    prime[:2] = False
    for p in range(2, math.isqrt(limit) + 1):
        if prime[p]:
            prime[p * p::p] = False
    return prime


class Translation:
    def __init__(self, modulus):
        self.modulus = modulus
        residues = np.arange(1, (modulus + 1) // 2, dtype=np.int64) ** 2 % modulus
        self.symbol = -np.ones(modulus)
        self.symbol[0] = 0
        self.symbol[residues] = 1
        assert len(np.unique(residues)) == (modulus - 1) // 2
        self.spectrum = np.fft.rfft(self.symbol)

    def sums(self, labels, weights):
        # Cyclic convolution with the reversed array gives lambda(t+j).
        coefficients = np.zeros(self.modulus)
        coefficients[-labels % self.modulus] = weights
        return np.fft.irfft(np.fft.rfft(coefficients) * self.spectrum, n=self.modulus)


def probes(group, labels, left, right, weights, alpha=ALPHA):
    phase = 2 * np.pi * (6 * labels - 1) * alpha
    gram = 0.0
    independent_log = 0.0
    group_log = 0.0
    characteristic = []
    for axis in (np.cos(phase), np.sin(phase)):
        direct = np.sqrt(weights) * axis
        l = direct * left
        r = direct * right
        gram += (group.sums(labels, l) @ group.sums(labels, r)
                 + l.sum() * r.sum()) / group.modulus
        s = np.sign(axis) * np.arctan(np.sqrt(np.tanh(weights * axis**2)))
        l, r = s * left, s * right
        minus = np.cos(group.sums(labels, l - r))
        plus = np.cos(group.sums(labels, l + r))
        characteristic.extend((minus, plus))
        assert min(minus.mean(), plus.mean()) > 0.01
        group_log += 0.5 * (np.log(minus.mean()) - np.log(plus.mean()))
        independent_log += 0.5 * np.sum(np.log(np.cos(l - r)) - np.log(np.cos(l + r)))
    exact = weights @ (left * right)
    assert abs(gram - exact) < 1e-11
    return gram, independent_log, group_log, characteristic


def major_coefficient(values, max_q, halfwidth):
    size = len(values)
    assert 2 * max_q**2 * halfwidth < 1
    freqs = np.fft.rfft(values, n=2 * size)
    autocorr = np.fft.irfft(abs(freqs) ** 2, n=2 * size)
    lags = np.arange(1 - size, size)
    harmonics = lags + 2
    ramanujan = np.zeros(len(harmonics))
    for q in range(1, max_q + 1):
        numerators = np.array([a for a in range(q) if math.gcd(a, q) == 1])
        phases = np.exp(2j * np.pi * np.outer(np.arange(q), numerators) / q)
        ramanujan += np.rint(phases.sum(axis=1).real)[harmonics % q]
    return np.dot(autocorr[np.abs(lags)],
                  2 * halfwidth * np.sinc(2 * halfwidth * harmonics) * ramanujan)


def arithmetic(group, u, report_samples=False):
    prime = primes_through(2 * u + 2)
    labels = np.arange((u + 1) // 6, (2 * u + 1) // 6 + 1)
    ns = 6 * labels - 1
    use = (ns > u) & (ns + 2 <= 2 * u)
    for p in (2, 3, 5, 7):
        use &= (ns % p != 0) & ((ns + 2) % p != 0)
    labels = labels[use]
    ns = 6 * labels - 1
    left, right = prime[ns], prime[ns + 2]
    w = 0.2 * math.log(u) ** 2 / u
    gram, independent, complete, means = probes(group, labels, left, right,
                                                np.full(len(labels), w))
    count = int(np.dot(left.astype(int), right.astype(int)))
    assert abs(independent / w - count) < 1e-10
    assert round(complete / w) == count
    print(f"U={u}, cells={len(labels)}, left={left.sum()}, right={right.sum()}, twins={count}, "
          f"Gram/w={gram/w:.8f}, translated log/w={complete/w:.8f}")
    if report_samples:
        phase = 2 * np.pi * ns * ALPHA
        rows = [np.sign(c) * np.arctan(np.sqrt(np.tanh(w * c**2))) * left
                for c in (np.cos(phase), np.sin(phase))]
        moments = []
        for coeff in rows:
            observed = group.sums(labels, coeff)
            var = np.mean(observed**2)
            assert abs(var - (coeff @ coeff - coeff.sum()**2 / group.modulus)) < 1e-11
            cumulant = np.mean(observed**4) - 3 * var**2
            moments.append((var, cumulant, -2 * sum(coeff**4)))
        var, cumulant, independent_fourth = max(moments, key=lambda triple: triple[0])
        print(f"  larger-variance left projection: variance={var:.9f}, fourth cumulant={cumulant:.9f}, "
              f"independent-sign fourth={independent_fourth:.9f}")
        stride = 2 * u + 3
        for copies in (4, 16, 32):
            assert copies * stride + int(labels[-1]) < group.modulus
            rms = []
            for values in means:
                sampled = sum(np.roll(values, -b * stride) for b in range(copies)) / copies
                rms.append(np.sqrt(np.mean((sampled - values.mean())**2)))
            print(f"  copies={copies}, max characteristic RMS over origins={max(rms):.6f}")
        # A genuine interval with both marginals present and no twin row.
        cuts = np.r_[-1, np.flatnonzero(left & right), len(labels)]
        runs = [slice(a + 1, b) for a, b in zip(cuts, cuts[1:])
                if left[a + 1:b].any() and right[a + 1:b].any()]
        span = max(runs, key=lambda section: section.stop - section.start)
        sub = labels[span]
        _, product, measured, _ = probes(group, sub, left[span], right[span], np.full(len(sub), w))
        assert abs(product) < 1e-13 and round(measured / w) == 0
        print(f"  no-twin starts [{6*sub[0]-1},{6*sub[-1]-1}], "
              f"left={left[span].sum()}, right={right[span].sum()}, translated log/w={measured/w:.8f}")
        def has_chen_factors(p):
            x = p + 2
            return any(x % r == 0 and prime[r] and prime[x // r] and r**11 > p**3
                       for r in range(2, math.isqrt(x) + 1))
        semi = np.array([has_chen_factors(int(n)) for n in ns])
        assert not np.any(semi & right)
        observed = []
        for indicator in (right | semi, semi):
            _, product, measured, _ = probes(group, labels, left, indicator, np.full(len(labels), w))
            expected = np.count_nonzero(left & indicator)
            assert abs(product / w - expected) < 1e-9 and round(measured / w) == expected
            observed.append(measured / w)
        print(f"  Chen log/w={observed[0]:.8f}, semiprime log/w={observed[1]:.8f}")

    theta = np.log(np.arange(u + 1, 2 * u + 1)) * prime[u + 1:2 * u + 1]
    # At most the two boundary cells are incomplete. For this check we
    # only weight complete cells, since the tested gap has both ends inside.
    t = 1234
    signs = group.symbol[(t + labels) % group.modulus]
    assert np.all(signs != 0)
    masked = theta.copy()
    masked[ns - u - 1] *= 1 + signs
    masked[ns + 2 - u - 1] *= 1 - signs
    pair_mass = theta[:-2] @ theta[2:]
    assert masked[:-2] @ masked[2:] == 0
    max_q = max(1, math.isqrt(u) // 5)
    halfwidth = 2 / u
    for q in range(1, max_q + 1):
        assert abs(ALPHA * q - round(ALPHA * q)) / q > halfwidth
    full_major = major_coefficient(theta, max_q, halfwidth)
    masked_major = major_coefficient(masked, max_q, halfwidth)
    print(f"  raw major/U={full_major/u:.6f}, raw minor/U={(pair_mass-full_major)/u:.6f}, "
          f"opposite major/U={masked_major/u:.6f}, opposite minor/U={-masked_major/u:.6f}")
    angles = 2 * np.pi * ALPHA * np.arange(u + 1, 2 * u + 1)
    values = theta * np.exp(1j * angles)
    a = values[ns - u - 1]
    b = values[ns + 2 - u - 1]
    ls = group.sums(labels, a.real) + 1j * group.sums(labels, a.imag)
    rs = group.sums(labels, b.real) + 1j * group.sums(labels, b.imag)
    contrast = np.mean(abs(values.sum() + ls + rs)**2 - abs(values.sum() + ls - rs)**2)
    correction = (abs(a.sum() + b.sum())**2 - abs(a.sum() - b.sum())**2) / group.modulus
    recovered = (contrast + correction) / (4 * math.cos(4 * np.pi * ALPHA))
    assert abs(recovered - pair_mass) < 1e-9 * u
    print(f"  single-frequency energy contrast/U={recovered/u:.9f}, weighted twins/U={pair_mass/u:.9f}")


def zeros_test(group):
    from flint import ctx

    previous = Path(__file__).resolve().parents[1] / "gap-two-auxiliary-conductor-bands" / "verify.py"
    source = runpy.run_path(str(previous))
    ctx.prec = 75
    family = source["OddFamily"](Q)
    zeros = family.locate_zeros()
    center, radius, box, half = 1.25, 48.0, 0.14, 2.2
    traces = []
    for row in zeros:
        ts = row[np.abs(row - center) <= radius] - center
        transform = 2 * half * np.sinc(half * ts / np.pi) * np.sinc(box * ts / np.pi) ** 4
        traces.append(transform.sum())
    law = irwinhall(4)
    roots = source["roots_by_integer"](13)
    labels = np.array([1, 2])
    ns = 6 * labels - 1
    indicators = {}
    for n in (*ns, *(ns + 2)):
        total = 0j
        for r, j in roots[int(n)]:
            support = range(2, math.floor(math.exp(half + 4 * box)) + 1)
            assert [m for m in support if m % Q in (r, Q - r)] == [r]
            assert all(m * r % Q not in (1, Q - 1) for m in support)
            taper = (law.cdf((math.log(r) + half) / (2 * box) + 2)
                     - law.cdf((math.log(r) - half) / (2 * box) + 2))
            projected = -np.mean(family.table[:, r - 1].conj() * traces)
            total += (source["moebius"](j) * math.sqrt(r) * np.exp(1j * center * math.log(r))
                      * projected / taper)
        indicators[n] = (total / math.log(n)).real
    left = np.array([indicators[n] for n in ns])
    right = np.array([indicators[n + 2] for n in ns])
    w = 0.2
    gram, independent, complete, _ = probes(group, labels, left, right, np.full(2, w))
    assert round(gram / w) == round(independent / w) == round(complete / w) == 2
    print(f"actual zeros: max prime-indicator error={max(abs(x-1) for x in indicators.values()):.6g}, "
          f"Gram/w={gram/w:.8f}, translated log/w={complete/w:.8f}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--arithmetic-only", action="store_true")
    args = parser.parse_args()
    assert Q % 8 == 3 and D % 16 == 7 and (D + 1) % (3 * 5 * 7 * Q) == 0
    assert primes_through(D)[[Q, D]].all()
    group = Translation(D)
    for u in (512, 2048, 8192):
        arithmetic(group, u, report_samples=(u == 2048))
    if not args.arithmetic_only:
        zeros_test(group)


if __name__ == "__main__":
    main()
