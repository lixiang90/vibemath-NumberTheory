---
title: "3 或 5 整除各非欧拉指数加一时，最小素因子的两个上界"
author: "Codex"
date: "2026-11-24"
abstract: "设奇完全数的每个非欧拉素因子的指数加一都被 3 或 5 整除。将 Fletcher、Nielsen 与 Ochem 所得的最小素因子上界从 10¹⁰⁰⁰ 降至 10³²。若每个非欧拉指数还不超过一亿，上界为 10²⁸。对任意奇完全数，若最小素因子超过 10³²，则指数加一与 15 互素的素因子的倒数权重之和大于 1/1000。"
formula: '10^8<P<10^{32}'
---

# 3 或 5 整除各非欧拉指数加一时，最小素因子的两个上界

设 $N$ 是奇完全数，$P$ 是其最小素因子，$e_q=v_q(N)$。Fletcher、Nielsen 与 Ochem 证明，当每个非欧拉素因子 $q$ 都满足 $3\mid e_q+1$ 或 $5\mid e_q+1$ 时，$10^8<P<10^{1000}$。本文将其上界改进为

$$
10^8<P<10^{32}.
$$

这个结论允许 $e_q$ 任意大。若额外假定每个非欧拉指数不超过 $10^8$，则 $P<10^{28}$，欧拉指数仍可任意大。特别地，若 $N$ 不被任何素数的六次幂整除，其非欧拉指数均为 2 或 4，从而 $10^8<P<10^{28}$。下界 $10^8$ 沿用原论文。[Fletcher–Nielsen–Ochem 主定理](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

我们证明如下更一般的定量结论。对任意奇完全数，令 $E=\{q\mid N:q\text{ 为素数},\ \gcd(e_q+1,15)=1\}$。若 $P>10^{32}$，则

$$
\sum_{q\in E}\frac1{q-1}>\frac1{1000}.
$$

因此 $\#E>(P-1)/1000$。令 $Y=(P-1)^{2501/2500}$，还有 $\#\{q\in E:q>Y\}>(Y-1)/2000$。这些素因子中至多有一个是欧拉素因子，其余各指数均至少为 6。在文首的假设下，$E$ 至多含欧拉素因子，倒数权重至多为 $1/(P-1)$，故推出所述最小素因子上界。

## 两类指数各自的筛

先设 $P>10^{32}$，置 $X=P-1,\ b=\log X>73.6$。对 $\ell=3,5$，令 $R_\ell=\{q\mid N:q\text{ 为素数},\ \ell\mid e_q+1\}$，定义

$$
u_\ell=\max_{a\bmod15}\ \sum_{\substack{q\in R_\ell\\q\equiv a\bmod15}}\frac1{q-1},
\qquad z=\sum_{q\in E}\frac1{q-1}.
$$

只须排除 $z\le z_0=1/1000$。全部素因子的倒数权重 $W$ 满足

$$
\log2<W<7/10 .
$$

下界来自 $\log(\sigma(q^{e_q})/q^{e_q})<1/(q-1)$。上界来自 $\prod_{q\mid N}(1+1/q)\le2$，以及 $1/(q-1)\le((q+1)/(q-1))\log(1+1/q)$。本节所有误差估计都可取 $\epsilon=10^{-6}$。

若 $s\in R_\ell$，则 $\Phi_\ell(s)\mid2N$。因 $\ell<P$，还必须有 $s\not\equiv1\bmod\ell$。所以 $R_3$ 在模 15 的既约类中只占 $2,8,11,14$ 四类，$R_5$ 只占 $2,4,7,8,13,14$ 六类。特别地，

$$
W\le 4u_3+6u_5+z.
$$

对于 $t\ge X$，令 $v=\log t,\ Q=\lfloor\sqrt t/24\rfloor$。固定一个允许类 $a$，写 $s=a+15j$，相应 $j$ 属于一个长度至多 $t/15+1$ 的整数区间。对 $q\nmid15$，仿射替换保持禁用类数。对 $q\equiv1\bmod\ell$ 且 $q\nmid N$ 禁用 0 和 $\ell-1$ 个非平凡单位根，其余只禁用 0。

设 $G_{0,\ell}(Q)$ 为对所有 $q\equiv1\bmod\ell$ 都禁用单位根时的完整分母，含 $q=3,5$ 处各一个零类。利用小素数的高次幂所得的分母引理给出

$$
G_{0,\ell}(Q)\ge c'_\ell(\log Q-a'_\ell)_+^2,\qquad
\begin{aligned}
 c'_3&=(7/50)(1.149), &a'_3&=(\log74)/6+0.77,\\
 c'_5&=(6/55)(1.165), &a'_5&=(7/15)\log40+1 .
\end{aligned}
$$

在我们的范围内 $Q>74$。对 $q\mid N,\ q\equiv1\bmod\ell$ 放宽筛后，分母的总损失因子至多为

$$
M_\ell(Q)=\prod_{\substack{q\mid N,\ q\le Q\\q\equiv1\bmod\ell}}
       \frac{q-1}{q-\ell}.
$$

从分母删除 3 和 5 的零类再损失至多 $(3/2)(5/4)=15/8$。这些比较均来自平方自由分母在单个素数处为 $A+g(q)B$，且 $0\le B\le A$。先计入至多 $Q$ 个小素数，再对单个类应用大筛。因此 $A_{\ell,a}(t)=\#\{s\in R_\ell:s\le t,\ s\equiv a\bmod15\}$ 满足

$$
A_{\ell,a}(t)\le Q+
 \frac18\left(1+\frac{15}{576}+\frac{15}{t}\right)
 \frac{4}{c'_\ell}\frac{tM_\ell(Q)}{(v-d_\ell)^2},
\qquad d_3=9.35,\quad d_5=11.82.
$$

这里 $Q\ge\sqrt t/24.1$，并且 $2\log(24.1)+2a'_\ell<d_\ell$。这些分母和计数不等式在稍后所用的 $X\ge10^{28}$ 时仍有效。[三次分母的证明](ternary-tail-constants.md) · [小素数的加强式与五次分母](quintic-denominator.md) · [大筛不等式](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

## 为什么要分出模 15 余 1 的素数

若 $q\mid N$ 且 $q\equiv1\bmod\ell$，则 $\sigma(q^{e_q})\equiv e_q+1\bmod\ell$。因 $\ell\nmid2N$，有 $q\notin R_\ell$。特别地，模 15 余 1 的素因子都在 $E$ 中。它们在 $\log M_\ell$ 中的贡献至多为 $(1+\epsilon)(\ell-1)z_0$。

对大于 $X$ 的其余模 $\ell$ 余 1 的素数，其模 15 剩余类列在下表中。记 $m$ 为 $\{3,5\}$ 中异于 $\ell$ 的另一项。

| $\ell$ | 排除模 15 余 1 后的可能类 | 类数 $k_\ell$ | $\alpha_\ell=(\ell-1)k_\ell/8$ |
| --- | --- | ---: | ---: |
| 3 | $4,7,13$ | 3 | $3/4$ |
| 5 | $11$ | 1 | $1/2$ |

这些类中实际整除 $N$ 的素数属于 $R_m\cup E$。它们在 $\log M_\ell$ 中的贡献不超过 $(1+\epsilon)m_\ell$，其中

$$
m_3=2\big(\min(0.7,3u_5)+z_0\big),\qquad m_5=4(u_3+z_0).
$$

这里三个模 3 余 1 的类至多贡献 $3u_5$ 的 $R_5$ 权重，而模 5 余 1 且异于 1 的单个模 15 类至多贡献 $u_3$ 的 $R_3$ 权重。也可以将这些类中 $(X,Q]$ 的全部素数计入，得到随 $Q$ 变化的第二个上界。

具体地，对模 15 的一个既约类 $a$，Bennett、Martin、O'Bryant 与 Rechnitzer 的定理 1.3 给出 $|\pi(u;15,a)-\operatorname{Li}(u)/8|\le u/(840\log^2 u)$，对 $u\ge X$ 有效。在 $Q\ge X$ 时，以 $1/u$ 分部求和，得

$$
\sum_{\substack{X<q\le Q\\q\equiv a\bmod15}}\frac1q
\le \frac18\log\frac{\log Q}{b}+
\frac1{840}\left(\frac2{b^2}+\frac1b\right).
$$

结合 $\log((q-1)/(q-\ell))\le(1+\epsilon)(\ell-1)/q$，得到下面的界。在 $Q<X$ 时空乘积为 1。

$$
M_\ell(Q)\le
\exp\!\left((1+\epsilon)(\ell-1)z_0+\frac1{8000}\right)
\min\!\left(e^{m_\ell},\
 \max\left(1,\left(\frac{\log Q}{b}\right)^{\alpha_\ell}\right)\right).
$$

为核对误差，记 $J=\max(0,\alpha_\ell\log(\log Q/b))$。除去模 15 余 1 的素数后，两个上界中的较小者至多为 $(1+\epsilon)\min(m_\ell,J)+(1+\epsilon)(6/840)(2/b^2+1/b)$。因 $m_\ell<3$，右式小于 $\min(m_\ell,J)+1/8000$。这一估计同时利用了小 $Q$ 时的素数分布和大 $Q$ 时实际素因子总量的限制。[等差数列素数的显式误差](https://arxiv.org/abs/1802.00085)

## 倒数积分与三轮比较

对每个允许类的计数界分部求和，该类的权重至多为 $\int_X^\infty A_{\ell,a}(t)/(t-1)^2\,dt$，且 $(X/(X-1))^2<1+\epsilon$。置 $u=(v-d_\ell)/(b-d_\ell)$。因 $v\ge b$，有 $\log Q/b\le v/(2b)\le u/2$。对 $0<\alpha<1$ 和 $s\ge0$，

$$
\begin{aligned}
F_\alpha(s)
&=\int_1^\infty
 \min\big(e^s,\max(1,(u/2)^\alpha)\big)\frac{du}{u^2}\\
&=\frac12+\frac{1-\alpha\exp(-s(1-\alpha)/\alpha)}{2(1-\alpha)} .
\end{aligned}
$$

分点为 $2$ 和 $2\exp(s/\alpha)$。又 $Q(v-d_\ell)^2/t<1/10000$，可把计数界中的首项并入系数。取上界常数 $K_3=3.21,\ K_5=4.07$，有

$$
u_\ell\le\frac{K_\ell}{b-d_\ell}F_{\alpha_\ell}(m_\ell).
$$

具体地，这两个 $K_\ell$ 均严格大于

$$
\left[
\frac18\left(1+\frac{15}{576}+\frac{15}{10^{28}}\right)
\frac4{c'_\ell}+\frac1{10000}
\right](1+\epsilon)
\exp\left((1+\epsilon)(\ell-1)z_0+\frac1{8000}\right).
$$

函数 $F_\alpha$ 单调递增。以 $b>73.6$ 和 $(u_3,u_5)<(0.700,0.700)$ 开始，每一步仅使用上一行的界。为方便，可将 $d_3,d_5$ 分别放宽为 $9.5,12$。代入 $F_{3/4}(s)=5/2-(3/2)e^{-s/3}$ 和 $F_{1/2}(s)=3/2-(1/2)e^{-s}$，得到

| 轮次 | $u_3$ 上界 | $u_5$ 上界 |
| ---: | ---: | ---: |
| 1 | 0.079 | 0.098 |
| 2 | 0.064 | 0.076 |
| 3 | 0.061 | 0.074 |

于是 $\log2<W\le4(0.061)+6(0.074)+0.001=0.689<\log2$，矛盾。验证器用有理数计算指数函数的泰勒下界，并对表中的小数取向外界。[常数与迭代复核程序](verify_mixed_prime_bound.py)

## 远处的例外素因子

取 $Y=X^{2501/2500}$。对模 3 的两个既约类使用同一素数定理，再对 $1/t$ 分部求和。把分母换成 $q-1$ 后，有

$$
\sum_{X<q\le Y,\ q\ {\rm prime}}\frac1{q-1}
\le(1+\epsilon)\left(\log(2501/2500)+\frac2{840}
\left(\frac2{b^2}+\frac1b\right)\right)<\frac1{2000}.
$$

而 $E$ 至少提供 $1/1000$ 的权重，故在 $Y$ 以外还须提供大于 $1/2000$ 的权重。每项小于 $1/(Y-1)$，得到文首的个数下界。

## 非欧拉指数有界时

现在只假设 $P>10^{28}$，并假设每个非欧拉指数至多为 $H=10^8$。对模 $\ell$ 余 1 的每个筛素数 $q\le Q$，若 $q\mid N$ 且非欧拉，则 $v_q(2N)\le H$。又 $\prod_{s\in R_\ell}\Phi_\ell(s)\mid2N$，所以至多有 $H$ 个不同的 $s\in R_\ell$ 落入模 $q$ 的非平凡 $\ell$ 次单位根。把这些例外在单个模 15 类中全部先计入，最多增加 $H\pi(Q)$ 项。其余底数可使用完整单位根筛，欧拉素因子若属于筛素数，仅损失一个因子 $(q-1)/(q-\ell)<1+\epsilon$。

因此对单个类的计数界，可把前面的 $Q$ 改为 $(H+1)Q$，同时把 $M_\ell(Q)$ 改为 $1+\epsilon$。因 $v^2e^{-v/2}$ 在 $v>4$ 时递减，且 $28\log10<65$，对 $t\ge X$ 有

$$
\frac{(H+1)Q(v-d_\ell)^2}{t}
\le \frac{(10^8+1)\,65^2}{24\cdot10^{14}}.
$$

同样分部求和，此时积分中不再有随 $Q$ 变化的损失因子。以 $3.19,4.04$ 分别给两个系数取上界。若 $z\le1/2000$，则

$$
W\le\frac{4(3.19)}{28(2.302)-9.35}
   +\frac{6(4.04)}{28(2.302)-11.82}+\frac1{2000}
<0.693<\log2 ,
$$

其中用了 $\log10>2.302$。因此，若 $P>10^{28}$ 且所有非欧拉指数均不超过一亿，则 $E$ 的倒数权重超过 $1/2000$。在混合指数假设下至多有一个例外，所以 $P<10^{28}$。

## 按大小排列的素因子

混合条件还给出每个固定次序的显式上界。对 $Y\ge3$，$R_3\cup R_5$ 中大于 $Y$ 的素因子的倒数权重至多为 $(188+1058)/\log Y$。在一般的部分丰度间距引理中取联合尾和常数 1246，至多一个欧拉例外，用本文的 $B_1=10^{32}$ 起步。对 $j\ge2$ 置 $t=j-1$ 和

$$
B_j=\exp\!\left(3744\,(4t B_{j-1}^{\,t})^{2^t}\right).
$$

若 $N$ 至少有 $j$ 个不同素因子，第 $j$ 小者至多为 $B_j$。非欧拉指数均不超过一亿时还可取 $B_1=10^{28}$。[部分丰度间距与递推的证明](prime-tails.md) · [三次尾界](ternary-tail-constants.md) · [五次尾界](quintic-denominator.md)

在本目录运行以下命令即可重算五次卷积前一千万项、所有显式常数和三轮倒数和比较。仅需 Python 标准库。

~~~sh
python -B -O verify_mixed_prime_bound.py
~~~
