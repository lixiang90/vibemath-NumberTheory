---
title: "对数平方短区间中的平衡哥德巴赫表示与零点对符号"
date: "2026-11-24"
abstract: "假设普通 RH，每个充分大的对数平方长度区间都含正比例的偶数，各有大量指定窗口内的素数表示。再假设 Dirichlet-GRH，可对同一批偶数确定低高度零点对统计的符号。"
---

# 对数平方短区间中的平衡哥德巴赫表示与零点对符号

对偶数 $N$，以 $P(N)$ 表示分别落在 $2N/5$、$3N/5$ 附近的两素数加权表示数。具体地，取 $c_1=2/5,c_2=3/5$，以及半宽总和为 $1/20$ 的均匀分布无穷卷积密度 $f$，令 $w_j(u)=u^{-1/2}f(\log(u/c_j))$。所用卷积半宽为 $a_j=\kappa/[j\log^2(2j)]$，由 $\sum a_j=1/20$ 决定 $\kappa$。于是 $P(N)=\sum_{p+q=N}\log p\log q\,w_1(p/N)w_2(q/N)$，求和只遍历素数。下引论文定义对应的同编号常数 $B_f$ 与不同编号的零点对和 $V_{T_N}(N)$。**对每个充分大的 $M$，当 $C_0\log^2 M\le h\le\sqrt M$ 时，$(M,M+2h)$ 至少包含 $\delta h$ 个偶数 $N$，它们同时满足 $V_{T_N}(N)<-1$ 并各有 $\gg N/\log^2 N$ 种窗口内的真正素数表示。** 常数 $C_0,\delta>0$ 只依赖固定窗口。素数表示的结论用普通 RH，翻译成 Dirichlet 零点对的结论用 Dirichlet GRH。以下 $M,h$ 为整数，$A>0$ 和相应的截高常数取定。[零点对和截高的定义](goldbach-low-zeros.md)

Goldston 在 1990 年用三角权及 Selberg 积分证明了哥德巴赫表示数的短区间平均定理。这里把两素数分别限制在指定比例上，保留归一化主项，并把平均结论转换成低高度零点对的符号信息。让两个短窗随同一实数 $y$ 移动，就能直接用两次 Selberg 均方界估计乘积。[Goldston 原文，式 (4)、(6) 及推论](https://doi.org/10.1017/S001708950000937X)

## 1. 保留零点信息的平均公式

令 $k_h(t)=(h-|t-h|)_+$，它在 $t=h$ 达到峰值，支撑为 $[0,2h]$。令

$$
F(u)=w_1(u)w_2(1-u),\qquad J=\int_0^1 F(u)\,du,\qquad
Z(N)=B_f-V_{T_N}(N).
$$

在 $(0,1)$ 外将 $F$ 延拓为零。对 $1\le h\le\sqrt M$ 有

$$
\boxed{\displaystyle
\frac1{h^2}\sum_{\substack{M<N<M+2h\\2\mid N}}k_h(N-M)Z(N)
=J+O_{f,A}\left(\frac{\log^2 M}{h}+\frac{h}{M}
+\frac{\log^2 M}{\sqrt M}+\frac1{(\log M)^A}\right).} \tag{1}
$$

全体整数上的三角权质量为 $h^2$，偶数上的质量为 $h^2/2+O(1)$。因此只在偶数上把权重归一成概率时，极限主项是 $2J$。若 $h/\log^2M\to\infty$，(1) 给出渐近式。对固定而充分大的 $C_0$，取 $h=\lceil C_0\log^2 M\rceil$ 就能得到严格正余量。

## 2. 为什么两边同时平滑可以省下一个对数

令 $\theta(x)=\sum_{p\le x}\log p$，$\Delta_h\theta(x)=\theta(x+h)-\theta(x)$，$E_h(x)=\Delta_h\theta(x)-h$。RH 下在任何固定正比例区间 $I\subset(0,2M)$ 上有

$$
\int_I |E_h(y)|^2\,dy\ll Mh\log^2 M.
$$

这是 Selberg 均方界。另用普通 RH 的 $\theta(x)=x+O(\sqrt x\log^2 x)$。本节直接对真正素数计算。[Goldston，式 (2)–(4)](https://doi.org/10.1017/S001708950000937X) · [Goldston–Suriajaya 的均方估计](https://arxiv.org/abs/2110.14250)

考虑

$$
I(M,h)=\int_\mathbb R F(y/M)\,\Delta_h\theta(y)\,\Delta_h\theta(M-y)\,dy .
$$

先写 $\Delta_h\theta=h+E_h$。常数项为 $Mh^2J$。对一次误差，将积分改成对 $y\in[p-h,p)$ 积分，得到

$$
\int F(y/M)\Delta_h\theta(y)\,dy
=h\sum_p\log p\,F(p/M)+O_f(h^2)
=hMJ+O_f(h\sqrt M\log^2 M+h^2).
$$

第一步用 $F$ 的有界导数及 $\sum_{p\asymp M}\log p\ll M$，第二步用分部求和。换 $y$ 为 $M-y$ 处理另一个一次项。二次项按实数 $y$ 使用 Cauchy 不等式，两边各自用一次 Selberg 均方，得到 $O_f(Mh\log^2 M)$。于是

$$
\frac{I(M,h)}{Mh^2}
=J+O_f\left(\frac{\log^2M}{h}+
\frac{\log^2 M}{\sqrt M}+\frac hM\right). \tag{2}
$$

再按素数对展开同一积分。给定 $p,q$ 后，两条件 $y<p\le y+h$ 和 $M-y<q\le M-y+h$ 的交集长度正好是 $k_h(p+q-M)$。令 $N=p+q$。在交集上，

$$
F(y/M)=\frac{M}{N}w_1(p/N)w_2(q/N)+O_f(h/M).
$$

这是因为 $|p-y|\le h,\ |N-M|\le2h$，而两个窗都是 Lipschitz 函数。只要任一权非零，$y$ 与 $M-y$ 都在固定正比例窗口内。在这个稍扩的支撑上，对两个 $\Delta_h\theta$ 再用 Cauchy 及均方界，乘积积分为 $O_f(M(h^2+h\log^2M))$。因此替换权重在 (2) 的归一化误差为 $O_f(h/M+\log^2 M/M)$。两素数都与 $M$ 同量级，充分大时其和必为偶数。所得公式正是 (1) 中把 $Z(N)$ 换为 $P(N)/N$ 的版本。上述逐偶数公式的逐偶数截断误差一致为 $O_{f,A}((\log M)^{-A})$，乘三角权后归一仍为同一量级。[逐偶数公式 (1)](goldbach-low-zeros.md)

## 3. 从平均到正比例的偶数

先给 $J$ 一个完全初等的下界。第一个卷积区间的半宽为 $a_1$，其余半宽之和至多为 $a_1\log 2$。于是 $f(v)=1/(2a_1)\ge10$ 在 $|v|<d=2a_1-b$ 上成立，且

$$
d\ge b\frac{1-\log2}{1+\log2}>\frac3{340}.
$$

对 $|u-2/5|\le 1/300$，两个对数 $\log(u/c_1)$、$\log((1-u)/c_2)$ 的绝对值分别至多为 $\log(120/119)$、$\log(180/179)$，均小于 $d$。因 $1/\sqrt{u(1-u)}\ge2$，在这一段上 $F(u)\ge200$，得到 **$J\ge4/3$**。

再控制均值会不会集中在极少数偶数上。令 $S(N)=P(N)/N$。哥德巴赫表示数的经典筛上界，乘以每对至多 $O_f(\log^2 N)$ 的权重，给出

$$
S(N)\ll_f g(N),\qquad
g(N)=\prod_{\substack{p\mid N\\p>2}}\frac{p-1}{p-2}.
$$

此处只需统一的上界常数。其局部因子也可以由二元筛直接读出，对 $n(N-n)$ 而言，每个奇素数通常删去两个剩余类，当 $p\mid N$ 时只删去一个。[Jie Wu，定理 1](https://arxiv.org/abs/0705.1652)

为容许更短的区间，先分开 $p\le z$ 与 $p>z$，取 $z=\log(2M)$。$N\le2M$ 的大于 $z$ 的不同素因子至多 $\log(2M)/\log z$ 个，所以这些素因子对 $\log g(N)^2$ 的总贡献为 $O(1/\log z)$。把剩余的积展开。令 $b(d)$ 支撑于奇的无平方因子整数，乘性，$b(p)=2/(p-2)+1/(p-2)^2$。对任意 $1\le H\le M$ 有

$$
\sum_{M<N\le M+H}g(N)^2
\ll H\sum_{d\ge1}\frac{b(d)}d+\prod_{2<p\le z}(1+b(p))
\ll H+(\log\log(2M))^2. \tag{3}
$$

这里每个 $d$ 的倍数至多有 $H/d+1$ 个。第一和的 Euler 乘积收敛，第二项由 Mertens 估计至多为 $O((\log z)^2)$。故在 $h\ge C_0\log^2 M$ 时有 $\sum_{M<N<M+2h,\;2\mid N}S(N)^2\ll_f h$，同一结论在 $h\asymp\log M$ 时也成立。

取 $C_0$ 充分大。(2) 的素数表示数版本给出 $\sum k_h S(N)\ge(5/4)h^2$。满足 $S(N)\le9/4$ 的偶数最多贡献 $(9/4)(h^2/2+O(1))$，于是其余偶数的贡献 $\gg h^2$。因 $k_h\le h$，这些偶数上的 $\sum S(N)\gg h$，结合 (3) 与 Cauchy 不等式，个数至少为 $\delta h$。每一个被选出的 $N$ 同时满足

$$
P(N)>\tfrac94N,\qquad
V_{T_N}(N)=B_f-S(N)+O_{f,A}((\log N)^{-A})<-1
$$

对最后一个不等式，逐偶数公式中常数的初等上界 $B_f<1.061$ 已经足够。每对的权重至多为固定常数乘 $\log^2N$，所以这些偶数各有 $\gg_f N/\log^2 N$ 种素数表示。$C_0$ 和 $\delta$ 在本推导中没有数值化。


## 4. 两个有限算术核验

下表在整数上检查两处容易弄错的量，即奇特征正交式的负号，以及两个短窗交集的三角权。测试时取素数指示函数和有理端点 $3/8<p/n<7/16,\ 9/16<q/n<5/8$，便于完全避开舍入。第一组逐项使用模 $2^k$ 奇特征之和，按特征的循环部分的阶合并成整数 Ramanujan 和。第二组分别按素数对、按 $y$ 的单位区间计数。

| 模数 $N$ | 128 | 256 | 512 | 1024 | 2048 | 4096 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| 两种计算所得的固定和表示数 | 0 | 1 | 1 | 1 | 3 | 5 |

| $(M,h)$ | $(1000,16)$ | $(10000,32)$ | $(100000,64)$ |
| --- | ---: | ---: | ---: |
| 两种计算所得的三角加权总数 | 334 | 8178 | 213568 |

$N=128$ 的零值发生在这两扇测试用窄窗内。两组数据核验算术恒等式。对第 1 节所用的光滑窗，逐偶数公式的第 7–8 节另用球算术核验 $B_f$，并用带截断误差的显式公式确定 $V_{20000}(64)$ 和 $V_{20000}(4096)$ 的符号。[常数及两模数的证明](goldbach-low-zeros.md)

```python
from bisect import bisect_left, bisect_right
from math import isqrt

def primes_below(limit):
    flags = bytearray([1]) * limit
    flags[:2] = b"\x00\x00"
    for p in range(2, isqrt(limit - 1) + 1):
        if flags[p]:
            flags[p * p :: p] = bytes(len(range(p * p, limit, p)))
    return [n for n in range(2, limit) if flags[n]]

def windows(n):
    primes = primes_below(n)
    return ([p for p in primes if 6*n < 16*p < 7*n],
            [p for p in primes if 9*n < 16*p < 10*n])

def odd_character_test(n):
    d = n // 4
    residues = {}
    for j in range(d):
        value = pow(5, j, n)
        residues[value] = j, 1
        residues[-value % n] = j, -1
    left, right = windows(n)
    right_set = set(right)
    direct = sum(n-p in right_set for p in left)
    def ramanujan(order, j):
        if order == 1:
            return 1
        if j % order == 0:
            return order // 2
        if j % (order // 2) == 0:
            return -(order // 2)
        return 0
    total = 0
    for order in (2**a for a in range(d.bit_length())):
        for p in left:
            j, sign_p = residues[p]
            for q in right:
                k, sign_q = residues[q]
                total += sign_p * sign_q * ramanujan(order, j-k)
    assert -total == d * direct
    return direct

def triangle_test(m, h):
    left, right = windows(m)
    by_pairs = sum(max(0, h-abs(p+q-m-h)) for p in left for q in right)
    by_windows = 0
    for k in range(min(left)-h, max(left)):
        count_p = bisect_right(left, k+h) - bisect_right(left, k)
        count_q = bisect_left(right, m-k+h) - bisect_left(right, m-k)
        by_windows += count_p * count_q
    assert by_pairs == by_windows
    return by_pairs

for n in (128, 256, 512, 1024, 2048, 4096):
    print(n, odd_character_test(n))
for m, h in ((1000, 16), (10000, 32), (100000, 64)):
    print(m, h, triangle_test(m, h))
```
