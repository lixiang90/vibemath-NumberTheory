---
title: "从有限个 ζ 零点认证相邻平方之间的素数"
date: "2026-11-24"
abstract: "在 ζ-RH 下给出正核谱尾的显式质量及纯素数计数下界。用前 2048 个正零点和区间算术认证 255 个平方窗口，并与独立整数筛逐窗核对。"
---

# 从有限个 ζ 零点认证相邻平方之间的素数

假设 Riemann ζ 函数满足 RH。用正频率上的前 $K$ 个零点，可以给开区间 $(n^2,(n+1)^2)$ 内的真素数个数算出严格下界。方法是把截断后全部未知相位的最坏贡献，换成可从极短尺度精确求得的正核总质量。取 $K=2048$，在 $2\le n\le256$ 的 255 个区间中，这个下界全为正。其中 $n=256$ 时，单靠谱下界就给出至少 15 个素数，随后整数筛数得 53 个。这个有限实验检验从零点到算术的误差控制。

## 1. 显式公式与正核尾部

记非平凡零点为 $1/2+i\gamma$，两个符号的纵坐标都按重数计。令 $\gamma_E$ 为欧拉常数，$\Lambda$ 为 von Mangoldt 函数，并写
$$
A(y)=\sum_{m\le y}\Lambda(m)(y/m-1),\quad
V(y)=\sum_{j\ge1}\frac{y^{-2j}}{2j(2j+1)},\quad
H(u)=\sum_{\gamma}\frac{e^{i\gamma u}}{1/4+\gamma^2}.
$$
在 $y=e^u\ge1$ 有
$$
H(u)=\frac{A(y)-y\log y+(1+\gamma_E)y-\log(2\pi)+V(y)}{\sqrt y},
\qquad C_0=H(0)=2+\gamma_E-\log(4\pi).
$$
将 $-(\zeta'/\zeta)(s)y^s/[s(s-1)]$ 在 $\Re s>1$ 反演，向左移线。$s=1$ 的二阶极点、$s=0$ 和负偶数的留数分别为 $y(\log y-1-\gamma_E)$、$\log(2\pi)$ 和 $-V(y)$。非平凡零点给出上式的绝对收敛级数，$y=1$ 的等式由连续性得到。

置 $x=n(n+1)$，$u=\log x$，$d=\log((n+1)/n)$。对 $e^{-u/2}A(e^u)$ 施算子 $D_dF(u)=F(u+d)-2\cosh(d/2)F(u)+F(u-d)$。单个整数在窗口外对该差分的贡献为零，在窗口内留下权
$$
w_n(m)=
\begin{cases}
1/n-n/m,&n^2<m\le x,\\
(n+1)/m-1/(n+1),&x<m<(n+1)^2.
\end{cases}
$$
于是令 $k_d(t)=2(\cosh(d/2)-\cos(dt))/(1/4+t^2)\ge0$，得到
$$
W_n=\sum_{n^2<m<(n+1)^2}w_n(m)\Lambda(m)
 =d-\sum_\gamma k_d(\gamma)e^{i\gamma u}-R_n .
$$

这里主项 $2\sqrt{x}\,d\sinh(d/2)=d$。令 $\alpha_j=2j+1/2$，则
$$
0\le R_n=\sum_{j\ge1}
 \frac{2(\cosh(\alpha_jd)-\cosh(d/2))x^{-\alpha_j}}{2j(2j+1)}
 \le\frac1{3n^5(1-n^{-4})}.
$$

最后一步用 $xe^{-d}=n^2$ 逐项估计几何级数。

核的总质量恰为 $2(\cosh(d/2)C_0-H(d))$。由于 $e^d<2$，计算 $H(d)$ 时 $A(e^d)=0$，而
$$
V(y)=1+\tfrac12\big((y-1)\log(1-1/y)-(y+1)\log(1+1/y)\big),\quad y>1.
$$

若 $\gamma_1,\ldots,\gamma_K$ 是所用正纵坐标，令
$$
Z_K=2\sum_{j=1}^K k_d(\gamma_j)\cos(\gamma_j u),\qquad
B_K=2(\cosh(d/2)C_0-H(d))-2\sum_{j=1}^K k_d(\gamma_j).
$$

则实数 $\sum_\gamma k_d(\gamma)e^{i\gamma u}$ 与 $Z_K$ 的差绝对值至多为 $B_K$。这个上界只利用高零点的剩余总质量。

## 2. 从算术权扣除高次幂

连续平方之间没有偶次幂。对其中每个整数底数的奇次幂 $b^j$，$b\ge2,\ j\ge3$，按所有表示累加 $w_n(b^j)\log b$，得到 $P_n$。这覆盖了每个可能的非素数 von Mangoldt 项。于是纯素数加权和至少为
$$
L_n=d-Z_K-B_K-\frac1{3n^5(1-n^{-4})}-P_n.
$$
权在 $m=x$ 达到最大值 $1/x$。因此区间内有至少 $\max(0,\lceil xL_n/(2\log(n+1))\rceil)$ 个素数。若要沿此路线覆盖所有充分大的 $n$，一个充分条件是在这些尺度对完整相位和证明 $\sum_\gamma k_d(\gamma)e^{i\gamma u}\le(1-\eta)d$，其中 $\eta>0$ 固定。每个指数 $j\ge3$ 的可能底数在该窗口内至多有常数个，故高次幂代价为 $O(\log^2 n/n^2)=o(d)$。这里的逐平方相位估计仍有待证明。

## 3. 区间算术复算

下表的前两列数值分别向下、向上截取。整体最小的 $L_n/d$ 下界超过 $0.402793$。最后一列由独立的整数筛计算，并在估计谱余项之后才使用。

| $n$ | $L_n/d$ 下界 | $B_K/d$ 上界 | 证出的素数个数下界 | 筛得的实际个数 |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0.827637 | 0.004377 | 1 | 2 |
| 14 | 0.402793 | 0.025339 | 2 | 4 |
| 32 | 0.994620 | 0.056300 | 5 | 9 |
| 128 | 0.674993 | 0.231703 | 9 | 24 |
| 256 | 0.637900 | 0.431550 | 15 | 53 |

运行下列完整代码使用 Python 3.12 与 python-flint 0.9.0。默认采用 128 比特球算术，第 2048 个正纵坐标约为 $2565.975270478284$。程序同时检验每个窗口的 $W_n$ 双侧包围和纯素数个数下界。[FLINT 的 ζ 零点隔离](https://flintlib.org/doc/acb_dirichlet.html)

```python
"""Certify primes between squares from zeta ordinates. Requires python-flint 0.9.0.

Run: python verify_square_windows.py
The analytic identity uses RH. All bounds use Arb enclosures. The final sieve
checks the finite experiment independently and is not an input to the bound.
"""

import argparse
from math import isqrt

from flint import acb, arb, ctx


def positive_int(value):
    n = int(value)
    if n < 1:
        raise argparse.ArgumentTypeError("expected a positive integer")
    return n


def square_weight(n, m):
    """Hyperbolic difference of the ramp on n^2 < m < (n+1)^2."""
    if m <= n * (n + 1):
        return arb(m - n * n) / (n * m)
    return arb((n + 1) ** 2 - m) / ((n + 1) * m)


def odd_powers(last):
    """All integer-base odd powers, including duplicates, within each window."""
    groups = [[] for _ in range(last + 1)]
    stop = (last + 1) ** 2
    base = 2
    while base**3 < stop:
        m = base**3
        while m < stop:
            n = isqrt(m)
            if n * n < m:
                groups[n].append((m, base))
            m *= base * base
        base += 1
    return groups


def ramp_at_short_shift(x, euler, log_two_pi):
    """H(log x) for 1 < x < 2, where the prime sum is empty."""
    trivial = 1 + ((x - 1) * (1 - 1 / x).log() - (x + 1) * (1 + 1 / x).log()) / 2
    return (x * (1 + euler - x.log()) - log_two_pi + trivial) / x.sqrt()


def spectral_bounds(n, heights, mass, euler, log_two_pi, powers):
    x = arb(n + 1) / n
    d = x.log()
    center = n * (n + 1)
    u = arb(center).log()
    ch = (d / 2).cosh()
    full_mass = 2 * (ch * mass - ramp_at_short_shift(x, euler, log_two_pi))
    used_mass, phase_sum = arb(0), arb(0)
    for gamma, factor in heights:
        weight = factor * (ch - (gamma * d).cos())
        used_mass += weight
        phase_sum += weight * (gamma * u).cos()
    tail = full_mass - used_mass
    if not tail > 0:
        raise ArithmeticError("failed to resolve the remaining positive spectral mass")
    trivial_bound = 1 / (arb(3) * n**5 * (1 - arb(1) / n**4))
    lower = d - phase_sum - tail - trivial_bound
    upper = d - phase_sum + tail
    proper_powers = sum((square_weight(n, m) * arb(b).log() for m, b in powers), arb(0))
    prime_lower = lower - proper_powers
    largest_weight = 2 * arb(n + 1).log() / center
    min_count = max(0, int((prime_lower / largest_weight).lower().ceil().unique_fmpz()))
    return lower, upper, prime_lower / d, tail / d, min_count


def sieve_check_data(last):
    """Independent prime counts and weighted sums, obtained only from integers."""
    limit = (last + 1) ** 2 - 1
    prime = bytearray([1]) * (limit + 1)
    prime[:2] = bytes(2)
    for p in range(2, isqrt(limit) + 1):
        if prime[p]:
            prime[p * p :: p] = bytes((limit - p * p) // p + 1)
    weights = [arb(0) for _ in range(last + 1)]
    counts = [0] * (last + 1)
    for p in range(2, limit + 1):
        if not prime[p]:
            continue
        counts[isqrt(p)] += 1
        log_p = arb(p).log()
        m = p
        while m <= limit:
            n = isqrt(m)
            if m != n * n:
                weights[n] += log_p * square_weight(n, m)
            m *= p
    return weights, counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--last", type=positive_int, default=256)
    parser.add_argument("--zeros", type=positive_int, default=2048)
    args = parser.parse_args()
    if args.last < 2:
        parser.error("--last must be at least 2 (n=1 is immediate)")
    ctx.prec = 128
    zeros = acb.zeta_zeros(1, args.zeros)
    if len(zeros) != args.zeros:
        raise ArithmeticError("unexpected number of zero enclosures")
    half = arb(1) / 2
    heights = []
    previous = arb(0)
    for z in zeros:
        if z.real != half or not z.imag > previous:
            raise ArithmeticError("failed to isolate ordered critical-line zeros")
        previous = z.imag
        heights.append((z.imag, 4 / (z.imag**2 + half**2)))
    print(f"zeros={len(zeros)}; last ordinate={previous.str(16)}; precision={ctx.prec} bits")
    euler = arb.const_euler()
    mass = 2 + euler - (4 * arb.pi()).log()
    log_two_pi = (2 * arb.pi()).log()
    powers = odd_powers(args.last)
    results = {}
    for n in range(2, args.last + 1):
        results[n] = spectral_bounds(n, heights, mass, euler, log_two_pi, powers[n])

    weights, counts = sieve_check_data(args.last)
    good = 0
    for n, (lower, upper, ratio, _, min_count) in results.items():
        if not lower < weights[n] or not weights[n] < upper or min_count > counts[n]:
            raise ArithmeticError(f"independent sieve disagrees at n={n}")
        good += int(ratio > 0)
    # Midpoints choose the row to display; certified minima are taken on lower endpoints.
    worst = min(results, key=lambda n: float(results[n][2]))
    min_ratio = min(row[2].lower() for row in results.values())
    print(f"spectrally certified {good}/{len(results)} intervals; sieve and enclosures agree on all")
    print(f"min certified ratio >= {min_ratio.str(14)}; smallest midpoint at n={worst}")
    print("n   spectral prime lower/d     omitted mass/d          certified count   actual count")
    for n in sorted({2, 8, 32, 128, args.last, worst} & results.keys()):
        _, _, ratio, tail_ratio, min_count = results[n]
        print(n, ratio.str(12), tail_ratio.str(12), min_count, counts[n], sep="  ")


if __name__ == "__main__":
    main()
```
