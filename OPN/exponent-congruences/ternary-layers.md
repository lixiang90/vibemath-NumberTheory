---
title: "允许任意多个八次分量时，奇完全数中 3 的指数超过 1.875 × 10²²"
author: "Codex"
date: "2026-11-24"
abstract: "若奇完全数至多有一个非欧拉指数不在 {2,8} 中，则这个指数必在 3 上，且超过 1.875 × 10²²。将每个分量中由三的幂给出的分圆因子逐层计数，可证至多一层的全部素因子属于 {3,π}，其中 π 是欧拉素数。结合公开的二步筛，另给出不同素因子数的下界。"
formula: '\#\{q\mid N:q\ne\pi,\ e_q\notin\{2,8\}\}\le1\quad\Longrightarrow\quad v_3(N)>1.875\cdot10^{22}'
---

# 允许任意多个八次分量时，奇完全数中 3 的指数超过 1.875 × 10²²

设 $N$ 为奇完全数，$\sigma(N)=2N$，写 $e_q=v_q(N)$，以 $\pi$ 表示唯一具有奇指数的欧拉素因子。分量指精确整除 $N$ 的素数幂 $q^{e_q}$。

**定理 1。** 若至多一个非欧拉指数不属于 $\{2,8\}$，则恰好有一个这样的指数，它是 $a=v_3(N)$。此时 $\pi\equiv1\bmod12$，并且

$$
a>1.875\cdot10^{22}.
$$

若 $a+1$ 的素因子均不超过 2850，则还必须有

$$
a+1\in\{2383^j:7\le j\le15\}\ \cup\
\{2579\cdot2383^j:6\le j\le10\}.
$$

两族此前在其余非欧拉分量全为平方的情形得到。定理 1 允许任意多个八次分量。[平方分量的结构定理](https://arxiv.org/html/2607.19746v1#S4) · [两族指数的整除证书](../restrained-exponent/shared-indices.md)

## 例外所在的素数与欧拉素数

更一般地，假设至多一个非欧拉素因子 $q$ 的 $e_q+1$ 不是 3 的正整数次幂。McDaniel 的定理要求至少一个非欧拉指数加一与 3 互素，故这个素因子恰好有一个，记为 $E$。单例外定理给出 $3\mid N$。若 $E\ne3$，则该定理还给出 $3\parallel a+1$ 和 $a\ge2570$。但此时 $a+1$ 必须是 3 的幂，只能等于 3，矛盾。故 $E=3$。[McDaniel 定理及另一证明](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf) · [单例外定理](one-exception.md)

若 $\pi\equiv5\bmod12$，则 $(3/\pi)=-1$。因为 $a+1$ 为奇数，$\pi\mid\sigma(3^a)$ 将使 3 模 $\pi$ 的乘法阶为奇数，从而 3 是平方剩余，矛盾。对其余每个非欧拉素因子写 $e_q+1=3^k$，它的约数和分解为

$$
\sigma(q^{e_q})=\prod_{j=0}^{k-1}\Phi(q^{3^j}),\qquad \Phi(x)=x^2+x+1 .
$$

每一项的非 3 素因子都模 3 余 1，因为 $x$ 模该素数的阶为 3。这样 $\pi$ 不能来自这些分量，也不能来自自身，与 $\pi\mid\sigma(N)$ 矛盾。欧拉条件 $\pi\equiv1\bmod4$ 遂给出 $\pi\equiv1\bmod12$。已验证的小素数排除又给出 $\pi$ 及所有 $q\notin\{3,\pi\}$ 均大于 $P=3\,000\,000$。[适用于任意 $e_q\equiv2\bmod6$ 的小素数排除](one-exception-support.md)

## 分圆因子逐层计数

本节只要求 $3\mid N$ 且 $3\mid e_q+1$ 对 $R=\{q:q\mid N,\ q\text{ 为素数},\ q\notin\{3,\pi\}\}$ 中所有 $q$ 成立。令 $R_1=\{q\in R:q\equiv1\bmod3\}$、$M=|R|$、$k_q=v_3(e_q+1)$。设 $h=v_3(e_\pi+1)$、$W=(e_\pi+1)/3^h$。以 $L$ 表示 $(\pi+1)/2$ 中模 3 余 1 的素因子数，按重数计。则有

$$
\sum_{q\in R} k_q+h\tau(W)+L\ \le\ 1+\sum_{t\in R_1}e_t . \tag{1}
$$

这里 $\tau$ 表示正约数个数。为证 (1)，对每个 $q\in R$ 取 $k_q$ 层 $\Phi(q^{3^j})$，$0\le j<k_q$。它们的乘积整除 $\sigma(q^{e_q})$，且所有底数 $x=q^{3^j}$ 两两不同。每层的非 3 素因子必属于 $R_1\cup\{\pi\}$。若一层没有 $R_1$ 素因子，则 $\pi\equiv1\bmod3$，并且它的值只能是 $\pi,3\pi,3\pi^2$。事实上，$v_3(\Phi(x))$ 为 0 或 1。对正奇数 $x>1$，Bennett–Levin 命题 1 应用于 $(x^3-1)/(x-1)=y^u$，排除 $u\ge2$。另有结论表明，$\Phi(x)=3y^u$ 对 $y>1,\ u\ge3$ 无解。前一结论允许 $x$ 为奇合数，这一点在 $j>0$ 时必需。[Bennett–Levin，命题 1 及四个已知解](https://arxiv.org/pdf/1312.4037) · [Bugeaud–Mihăilescu，定理 2](https://doi.org/10.7146/math.scand.a-15038) · [两条结论在平方分量中的用法](https://arxiv.org/html/2607.19746v1#S4)

这三个值在所有层中合计至多出现一次。相同的值由严格单调性排除。对不同值，证明如下。

- 若 $\Phi(y)=\pi$、$\Phi(x)=3\pi$，写 $x=3s+1$，即 $3s(s+1)=y(y+1)$。$y$ 是异于 3 的奇素数的幂，故 $y$ 必整除 $s,s+1$ 之一。于是 $s\ge y-1$，而 $3(y-1)y>y(y+1)$。
- 若 $\Phi(x)=3\pi^2$，记 $t=(2x+1)/3$，则 $(2\pi-1)(2\pi+1)=3t^2$。两因子为互素奇数，且 $3\mid2\pi+1$，故 $2\pi-1=v^2,\ 2\pi+1=3u^2$。若还有 $\Phi(y)=\pi$，便有 $(v-y-1)(v+y+1)=y^2$。两因子为互素正奇数，因 $y$ 为素数幂，只能是 $1,y^2$，从而 $2(y+1)=y^2-1$、$y=3$，与底数的选择矛盾。
- 若同时有 $\Phi(x)=3\pi^2$ 与 $\Phi(y)=3\pi$，再取 $w=(2y+1)/3$，则 $w^2=2u^2-1,\ v^2=3u^2-2$。于是 $1,u^2,w^2,v^2$ 成非恒定的四平方等差数列，矛盾。[四平方等差数列定理](https://arxiv.org/abs/0712.3850) · [平方分量的三个例外值](../ordinary-components/counting.md)

因此至少 $\sum k_q-1$ 层各向 $R_1$ 贡献一份素因子重数。再取 $d\mid e_\pi+1$ 且 $3\mid d$，共有 $h\tau(W)$ 个。由 Bang–Zsigmondy 定理，每个 $d$ 给出一个 $\pi$ 模其乘法阶恰为 $d$ 的素数。它属于 $R_1$，不同 $d$ 不会给出同一素数。这里 $\pi$ 为奇素数且 $d\ge3$，不会遇到定理的例外。$(\pi+1)/2$ 另贡献 $L$ 份，其阶均为 2。汇总各来源，素数 $t$ 的总重数为 $v_t(\sigma(N))=e_t$，即得 (1)。[Zsigmondy 原定理](https://doi.org/10.1007/BF01692444)

若 $\pi\equiv1\bmod3$，对 $\sigma(N)=2N$ 比较 3 的赋值，有 $a=h+\sum_{t\in R_1}k_t$。现在回到定理 1，令 $b_1,b_2$ 分别表示 $R$ 中模 3 余 1、余 2 的八次分量数。式 (1) 给出

$$
M+b_2+h\tau(W)+L \le 1+2(a-h)+3b_1,\qquad
2b_1\le a-h . \tag{2}
$$

特别地 $M\le1+7a/2$。保留 $b_1,b_2$ 的式 (2) 还说明，如果 $2a$ 小于下节不同素因子数的下界，必须有相应数目的八次分量来自余 1 的素数。

## 不同素因子数的下界

**定理 2。** 假设 $3\mid N$，每个 $q\in R$ 都有 $3\mid e_q+1$，且 $\pi>P$。取 $D=6.5625\cdot10^{22}$，则

$$
M>D+1 . \tag{3}
$$

由小素数排除，每个 $q\in R$ 都大于 $P$，$\Phi(q)$ 的非 3 素因子也都大于 $P$。若这样的因子 $t$ 取非欧拉身份，则 $\Phi(t)\mid2N$。若取欧拉身份，则 $t\equiv1\bmod4$ 且 $t+1\mid2N$。这给出二步筛所用的全部整除条件。具体地，若 $t$ 的两种身份都强迫小素数，或两个不同的 $t$ 都必须取欧拉身份，即可删除 $q$。若只强迫某一欧拉素数，则按该素数把候选分组，至多保留一组。间隔上界所用的 $\log I(q^{e_q})<1/(q-1)$ 对任意正指数都成立，其中 $I(m)=\sigma(m)/m$。[小素数排除](one-exception-support.md) · [Stormwhip 的原始二步筛](https://github.com/stormwhip/OddPerfectNumberBounds/blob/dc6e8353293d8922a8ade8ffa625f4aeb2c72ae6/chain.cpp) · [截点提高到 300 万的二步筛](../restrained-exponent/refined-lower-bound.md)

为核对 (3)，反设 $M\le D+1$。沿用 300 万截点的公开证书，令 $A=328\cdot10^{20}$、$B=50\,290\,000\,000\,000$、$Z=290A$，以 $w(t)=1/(t-1)-1/(Z-1)$ 为权。所有 $q$ 的倒数权重 $\sum_{q\in R}1/(q-1)$ 至多为

$
\frac{M}{Z-1}+\sum_{\substack{q\in R\\ q\le Z}} w(q).
$

在 $(P,B]$ 上，二步筛给出的上界是 $V/10^{18}-H/(Z-1)$，其中 $V=52\,797\,916\,403\,925\,437$ 且 $H=188\,812\,885\,018$。在 $(B,Z]$ 上，116 项筛多项式只用各个小素数能否整除 $\Phi(q)$。欧拉分量的对数丰度至多为 $1/(P-1)$。因此改变 $M$ 的上界只会改变式中的第一项。[300 万截点数据](../restrained-exponent/chain-3000000.tsv) · [116 项筛多项式](../restrained-exponent/bmor-majorant.json) · [等差数列素数计数的显式估计](https://arxiv.org/abs/1802.00085)

当第一项的分子取 $2A+1$ 时，原证书给出总上界 $0.287678023229577646$。换成 $D+1$，保留相同的 $B,Z$、有限筛和尾段估计，便有

$
\log I(N/3^a)\le
0.287678023229577646+\frac{25\cdot10^{18}}{Z-1}
\le0.287680651488618858 .
$

完全性与 $I(3^a)<3/2$ 却要求左端大于 $\log(4/3)>0.287682072451780927$。这就证明 (3)。验证器重算公开的筛系数、代表见证和有理数区间，再加上这个分子修正。[原数值证书](../restrained-exponent/verify_refined_exponent.py) · [本次复算](verify_ternary_layers.py)

## 指数下界和光滑情形

将 (2) 与 (3) 合并，得 $7a>2D$，即 $a>1.875\cdot10^{22}$。若已知 $a+1$ 为 2850-光滑数，已有 903 个整除障碍都适用于例外在 3 的情形。按这些障碍枚举，再要求 $a>261\cdot10^{18}$，只余定理 1 所列 14 个整数。本次所得下界强于这一截点，而 14 个数中最小的 $2383^7$ 仍大于新下界加一。因此它们给出当前子类的必要候选。[903 个障碍的适用条件](one-exception.md) · [14 个整数的完整枚举](../restrained-exponent/shared-indices.md)

指数加一为更高的 3 次幂时也有定量结论。设除一个非欧拉素因子外，其余的指数加一均为 3 的正整数次幂，且仅对 $t\in R_1$ 要求 $k_t\le K$。由第一节 $E=3,\ \pi\equiv1\bmod12$，故 (3) 适用。从 (1) 减去 $\sum_{q\in R}(k_q-1)$，并用 $(3^k-k)/k$ 对 $k\ge1$ 递增，可得

$$
M+\sum_{q\in R\setminus R_1}(k_q-1)+h\tau(W)+L
\le 1+\sum_{t\in R_1}(3^{k_t}-k_t)
\le 1+\frac{3^K-K}{K}(a-h).
$$

于是 $(3^K-K)(a-h)>KD$。模 3 余 2 的非欧拉素因子可具有任意高的这类指数。

在仓库根目录运行以下命令可复算数值界。第二条另核对原始分段表和等差数列误差常数表。

~~~sh
python -I -B -O OPN/exponent-congruences/verify_ternary_layers.py
python -I -B -O OPN/exponent-congruences/verify_ternary_layers.py --sieve sieve.tsv --bmor c_all_rounded.txt
~~~
