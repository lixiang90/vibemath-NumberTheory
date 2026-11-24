---
title: "最小素因子如何限制奇完全数的指数同余分布"
author: "Codex"
date: "2026-11-24"
abstract: "以 P 表示奇完全数的最小素因子。对固定奇素数 ℓ，给出指数加一不被 ℓ 整除的素因子个数关于 P 的渐近下界。在每个非欧拉指数加一都被 3 或 5 整除的情形，至少 17,864 个非欧拉指数加一被 5 整除而不被 3 整除，至少 230 个满足相反的条件。"
formula: 'k_\ell\ \ge\ \left(\frac{\log 2}{\mathrm e(\ell-2)}-o(1)\right)P^{1/(\ell-1)}\log P'
---

# 最小素因子如何限制奇完全数的指数同余分布

设 $N$ 是一个奇完全数，$P$ 是其最小素因子，$a_q=v_q(N)$。对固定的奇素数 $\ell$，记

$$
k_\ell=\#\{q:q\text{ 为素数},\ q\mid N,\ \ell\nmid a_q+1\}.
$$

**定理 1。** 对每个 $\varepsilon>0$，当 $P$ 超过只依赖于 $\ell,\varepsilon$ 的常数时，有

$$
k_\ell\ \ge\
\left(\frac{\log 2}{\mathrm e(\ell-2)}-\varepsilon\right)
P^{1/(\ell-1)}\log P .
$$

这里 $\mathrm e=\exp(1)$，对数均为自然对数。定理对各指数的取值均无上限。

**定理 2。** 若 $P>10^8$，则 $k_3\ge17\,865$ 且 $k_5\ge231$。特别地，若每个非欧拉素因子的指数加一都被 3 或 5 整除，则两类不同的非欧拉素因子各自满足下表的下界。

| 对 $a_q+1$ 的条件 | 不同非欧拉素因子个数 |
| --- | ---: |
| 被 5 整除且不被 3 整除 | 至少 17,864 |
| 被 3 整除且不被 5 整除 | 至少 230 |

Fletcher、Nielsen 与 Ochem 已证明，上述混合整除条件蕴含 $10^8<P<10^{1000}$。我们取其下界应用定理 2。从每个 $k_\ell$ 中至多扣除一个欧拉素因子，余下的素因子由混合条件必属表中相应的一类。定理 2 的数值取自下文的有限不等式，可在给定的 $10^8$ 处直接使用。[Fletcher–Nielsen–Ochem，主定理](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

## 一个乘积能对应几个底数

设 $P>\ell$，$d=\ell-1$，$R=P^{1/d}$，$k=k_\ell$。令 $T$ 为整除 $N$ 且模 $\ell$ 余 1 的素数集。若 $t\in T$，则 $\sigma(t^{a_t})\equiv a_t+1\pmod\ell$，而 $\ell\nmid2N$，故 $|T|\le k$。

以下称满足 $\ell\mid a_q+1$ 的素因子 $q$ 为本节所计的底数。它不能模 $\ell$ 余 1。令

$$
m=\Phi_\ell(q)=1+q+\cdots+q^d .
$$

有 $m\mid\sigma(q^{a_q})$ 且 $m$ 为奇数。每个素因子 $t\mid m$ 都整除 $N$，故 $t\ne\ell$。若 $q\equiv1\pmod t$，则 $t\mid\ell$，矛盾。于是 $q$ 模 $t$ 的乘法阶为 $\ell$，从而 $t\in T$。

以 $j=\Omega(m)$ 表示按重数计的素因子个数。$m$ 完全由一个 $T$ 上的 $j$ 元多重集合决定，且 $\Phi_\ell$ 严格递增。因此第 $j$ 组至多有

$$
M_j=\binom{k+j-1}{j}
$$

个底数，且每个底数满足 $q>R^j-1$。后一不等式来自 $P^j\le m<(q+1)^d$。若 $k=0$，则 $T$ 为空而每个底数都需要 $T$ 中的因子，故可设 $k\ge1$。

由于 $\log(\sigma(q^{a_q})/q^{a_q})<1/(q-1)$，只要将这几组的倒数和与另外 $k$ 个素因子的贡献 $k/(P-1)$ 相加，就得到 $\log 2$ 的上界。

## 渐近估计

只须考虑 $P\to\infty$ 且 $k/(R\log P)$ 有界的序列。置 $\mu=k/R$，取

$$
J=\lceil 8(\mu+\log\log P)\rceil .
$$

有 $J=O(\log P)$。对 $j>J$，可用 $1/(q-1)\le2R^{-j}$。由二项级数，这些组的总贡献至多为

$$
2\sum_{j>J}M_j R^{-j}
\le 2^{1-J}(1-2/R)^{-k}=o(1).
$$

最后一步用到 $-k\log(1-2/R)\le3\mu$ 对充分大的 $P$ 成立。另有 $k/(P-1)=o(1)$。

对 $1\le j\le J$，取

$$
A_j=R^j/2,\qquad B_j=\max\{A_j,\,M_j(\log P)^3\}.
$$

第 $j$ 组中 $q>B_j$ 的倒数和 $\sum 1/(q-1)$ 至多为 $2/(\log P)^3$。余下底数在 $(A_j,B_j]$，模 $\ell$ 落在除 0 与 1 以外的 $\ell-2$ 个类中。记 $\delta=(\ell-2)/(\ell-1)$。固定模数的素数定理给出这些类的合计计数函数 $\pi_*(x)\le(\delta+o(1))x/\log x$，误差对 $x\ge R/2$ 一致趋于零。由分部求和，当 $B\ge A\ge R/2$ 时，

$$
\sum_{\substack{A<q\le B\\q\ {\rm prime}\\ q\not\equiv0,1\pmod\ell}}
\frac1{q-1}
\le(\delta+o(1))
\left(\frac1{\log B}+\log\frac{\log B}{\log A}\right)
\le(\delta+o(1))\frac{1+\log(B/A)}{\log A}.
$$

第一步可先对 $1/q$ 求和，舍去下端的非正项，再乘 $A/(A-1)$。这避免了对短区间内单独的素数渐近公式作任何假设。[等差数列素数定理及其定量形式，Bennett–Martin–O'Bryant–Rechnitzer](https://arxiv.org/abs/1802.00085)

写 $\nu=(k+J)/R$，$u_+=\max(u,0)$。由 $M_j\le(\mathrm e(k+J)/j)^j$ 得

$$
\log(B_j/A_j)
\le j\bigl(\log(\mathrm e\nu/j)\bigr)_+
+3\log\log P+\log 2 .
$$

对 $v>0$，递减函数 $\log(v/x)$ 在 $(0,v]$ 上的积分为 $v$，故 $\sum_{j\ge1}(\log(v/j))_+\le v$。又 $\log A_j=j\log R-\log 2$，其与 $j\log R$ 之比对 $1\le j\le J$ 一致趋于 1。将上界逐组相加，常数项与 $3\log\log P$ 项合计为 $O((\log\log P)^2/\log P)=o(1)$。取 $v=\mathrm e\nu$ 即得

$$
\log 2\le(\delta+o(1))\frac{\mathrm e(k+J)}{R\log R}+o(1)
=\bigl(\mathrm e(\ell-2)+o(1)\bigr)\frac{k}{R\log P}+o(1).
$$

这证明定理 1。

## 在一亿处的有理数估计

以下固定 $P_0=10^8<P$，在分组计数中用 $P_0$ 代替 $P$，令 $R_0=P_0^{1/d}$。为检验 $k\le K$，第 $j$ 组的容量一律取 $M_j=\binom{K+j-1}{j}$。在 $\ell=3,5$ 两种情形，$R_0$ 分别是整数 $10\,000,100$。

令 $W$ 为所有不超过 41、异于 $\ell$ 的素数之积。底数与 $W$ 互素，模 $\ell$ 有 $\ell-2$ 个允许类。任意长为 $x$ 的半开区间中，满足这两项同余条件的整数个数至多为 $\rho x+C$，其中

$$
\rho=\frac{\ell-2}{\ell}\prod_{p\mid W}(1-1/p),\qquad
C=(\ell-2)2^{\omega(W)} .
$$

这是对 $u\mid W$ 作容斥的结果。每个 $u$ 的倍数在每个允许类内构成模 $\ell u$ 的一个同余类，区间计数和 $x/(\ell u)$ 相差至多 1。

设 $D_j=\max(P_0-1,R_0^j-2)$。把第 $j$ 组实际出现的底数按增序编号，首项为 $q_0$。对区间 $[D_j+1,q_i)$ 应用上述计数，得

$$
q_i-1\ge D_j+(i-C)_+/\rho .
$$

因此对 $1\le j\le20$，可以累加 $i=0,\ldots,M_j-1$ 的相应倒数上界。为加速，将这些位置分为长度为 $\lceil M_j/512\rceil$ 的连续段，最后一段可较短。每段用首项的倒数上界乘以实际长度，再向上舍入到 $10^{-14}$ 的整数倍。

对 $j\ge21$，改用 $1/(q-1)\le(P_0/(P_0-2))R_0^{-j}$。$M_jR_0^{-j}$ 的相邻项比值为 $(K+j)/((j+1)R_0)$，在这里随 $j$ 递减。因此这些组的贡献至多为

$$
\frac{P_0}{P_0-2}\cdot\frac{\binom{K+20}{21}}{R_0^{21}}
\left(1-\frac{K+21}{22R_0}\right)^{-1}.
$$

加上 $K/(P_0-1)$，得到如下严格上界。验证 $\log 2$ 时只用其正项展开 $2\sum_{i\ge0}3^{-(2i+1)}/(2i+1)$ 的前五项。

| $\ell$ | 假设的上限 $K$ | 全部对数丰度之和的上界 |
| ---: | ---: | ---: |
| 3 | 17,864 | $0.693108$ |
| 5 | 230 | $0.692821$ |

二者都小于 $4\,297\,606/6\,200\,145<\log 2$。完全性要求该和恰为 $\log 2$，从而分别排除 $k_3\le17\,864$ 和 $k_5\le230$。这证明定理 2。[独立运行的有理数复核程序](verify_mixed_counts.py)

还可用积分复核两个数值。当 $M_j>C+1$ 时，前 $C+1$ 项合计至多为 $(C+1)/D_j$，其余项的和至多为 $\rho\log(1+(M_j-C-1)/(\rho D_j))$。当 $M_j\le C+1$ 时用 $M_j/D_j$。把尾段的几何上界换成完整二项级数减去前 21 项，再乘 $P_0/(P_0-2)$，两行总和分别小于 $0.690268$ 与 $0.691109$。

从仓库根目录运行下列命令即可复算表中两个界。最小素因子可低至 5 时，关于模 3 指数的五个非欧拉底数下界另有证明。[无因子 3 的五个偶指数定理](three-free.md)

~~~sh
python -I -B -O OPN/exponent-congruences/verify_mixed_counts.py
~~~
