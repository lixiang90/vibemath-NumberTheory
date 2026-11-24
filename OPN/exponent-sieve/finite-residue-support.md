---
title: "从一个剩余类的素因子个数得到奇完全数的有限性"
author: "Codex"
date: "2026-11-24"
abstract: "若每个非欧拉指数加一都被 3 或 5 整除，则限制模 15 余 11 的素因子个数就能给整个奇完全数一个显式上界。对任意固定素数 ℓ≥5，也给出只依赖于 ℓ 在 N 中的指数及指数加一不被 ℓ 整除的素因子个数的上界。证明使用分圆值、Evertse 的 Thue–Mahler 定理及 Nielsen 的上界。"
formula: '\#\{q\mid N:q\equiv11\bmod15\}\le b\quad\Longrightarrow\quad N<2^{4^{K(b)}}'
---

# 从一个剩余类的素因子个数得到奇完全数的有限性

设 $N$ 是奇完全数，$e_q=v_q(N)$，本文集合中的 $q\mid N$ 均指素因子。先假设每个非欧拉指数加一都被 3 或 5 整除。令

$$
b=\#\{q\mid N:q\equiv11\pmod{15}\}.
$$

**定理 1。** 对每个给定的 $b_0\ge0$，满足 $b\le b_0$ 的上述奇完全数只有有限多个。下文给出只依赖于 $b_0$ 的显式上界。素因子取值和所有指数都可变化。同样，限制模 15 余 $4,7,13$ 的素因子总数，也能推出有效有限性。

**定理 2。** 对一般奇完全数，固定素数 $\ell\ge5$，置

$$
k_\ell=\#\{q\mid N:\ell\nmid e_q+1\},\qquad
a_\ell=v_\ell(N),\qquad h_\ell=a_\ell+k_\ell .
$$

则 $\omega(N)\le h_\ell+B_\ell(h_\ell)$，并有 $N<2^{4^{h_\ell+B_\ell(h_\ell)}}$，其中 $B_\ell$ 是下节的显式函数。因此固定 $\ell$ 且限制 $h_\ell$ 即得有效有限性。$\ell=3$ 时已有采用三次 Thue 方程的更小上界。[三次情形](../exponent-congruences/cubic-support.md) · [Nielsen，奇完全数关于不同素因子数的上界](https://mathdept.byu.edu/~pace/upper_bound_web.pdf)

## 有限个素数支持的分圆值

记 $F_\ell(X,Y)=\sum_{j=0}^{\ell-1}X^{\ell-1-j}Y^j$，$\Phi_\ell(x)=F_\ell(x,1)$。给定任意 $s$ 个素数组成的集合 $S$，使 $\Phi_\ell(x)$ 的全部素因子属于 $S$ 的正整数 $x>1$ 至多有 $B_\ell(s)$ 个。可以取

$$
B_\ell(0)=0,\qquad
B_\ell(s)=\left\lfloor\frac{(5\cdot10^6(\ell-1))^{s+1}-3}{2}\right\rfloor
\quad(s\ge1).
$$

事实上，$F_\ell$ 在 $\mathbb Q$ 上不可约，次数 $\ell-1\ge4$。这也可对 $\Phi_\ell(X+1)$ 用 Eisenstein 判别法直接验证。Evertse 定理 1 对 $\mathbb Q$、集合 $S$ 中的有限素位及实无限位给出至多 $(5\cdot10^6(\ell-1))^{s+1}$ 个射影解。其条件是坐标属于 $\mathbb Z[S^{-1}]$，且二元型取值为该环的单位。因 $F_\ell$ 对称，每个计入的 $x$ 至少给出斜率 $x,1/x$ 两个不同解。斜率 $0,\infty,-1$ 还各给出一个值为 1 的解。不同正整数 $x>1$ 的这两种斜率互不相同，故得到所需上界。空集合无法支持大于 1 的分圆值。[Evertse，The number of solutions of the Thue–Mahler equation，定理 1](https://pub.math.leidenuniv.nl/~evertsejh/95-Thue-Mahler.pdf)

当 $5\le\ell\le23$ 为素数时，还可以取 $B_\ell(1)=1$。Bugeaud 与 Mihăilescu 的定理 1 表明，正底数的 Nagell–Ljunggren 方程在已列明的解以外，长度的最小素因子至少为 29。因此 $\Phi_\ell(x)=r^m$、$m\ge2$ 在本范围内只可能为 $\Phi_5(3)=11^2$。当 $m=1$ 时，严格单调性给出至多一个 $x$。特殊值 $r=11,\ell=5$ 无法再给出 $m=1$，因为 $\Phi_5(2)=31$。[Bugeaud–Mihăilescu，定理 1](https://irma.math.unistra.fr/~bugeaud/travaux/PredaYann1.pdf)

Yamada 此前已将高次分圆值的解数估计用于所有非欧拉指数相等的形状，并对 $\Phi_\ell(x)=p^m r$ 证明了远小于一般 Thue–Mahler 界的解数上限。这里保留任意多个不同的指数，用整个素数集支持的分圆值来控制 $\omega(N)$。[Yamada，定理 1.1 与推论 1.2](https://www.iam.fmph.uniba.sk/amuc/ojs/index.php/amuc/article/download/1321/873/)

## 固定一个素数阶

记 $R_\ell=\{q\mid N:\ell\mid e_q+1\}$，$T_\ell=\{r\mid N:r\equiv1\bmod\ell\}$，$t=|T_\ell|$。对每个 $q\in R_\ell$ 有 $\Phi_\ell(q)\mid\sigma(q^{e_q})\mid2N$。此分圆值为奇数。若 $r\mid\Phi_\ell(q)$ 且 $r\ne\ell$，则 $q$ 模 $r$ 的乘法阶等于 $\ell$，故 $r\in T_\ell$。另一方面，$\ell\mid\Phi_\ell(q)$ 当且仅当 $q\equiv1\bmod\ell$。

因此每个 $R_\ell\cap T_\ell$ 中的素数都向 $\sigma(N)$ 贡献一个因子 $\ell$，总数至多为 $a_\ell$。$T_\ell$ 中余下素数均计入 $k_\ell$，故 $t\le h_\ell$。对 $R_\ell\setminus T_\ell$ 中的每个 $q$，$\Phi_\ell(q)$ 完全由 $T_\ell$ 支持，故这一部分至多有 $B_\ell(t)$ 个。相加得到

$$
\omega(N)\le k_\ell+\min(a_\ell,t)+B_\ell(t).
$$

代入 $t\le h_\ell$ 并用 Nielsen 的 $N<2^{4^{\omega(N)}}$，即得定理 2。这个较细的不等式还表明，只限制 $k_\ell$ 和 $t$ 也足够，即使 $a_\ell$ 本身可任意大。[Nielsen 的上界](https://mathdept.byu.edu/~pace/upper_bound_web.pdf)

## 由模 15 余 11 控制两种分圆值

回到定理 1。Fletcher、Nielsen 与 Ochem 已证明该混合指数条件蕴含最小素因子 $P>10^8$。因此 $3,5\nmid N$，$T_3\cap R_3$ 和 $T_5\cap R_5$ 都为空。任意素因子若模 15 余 1，则既不在 $R_3$ 中也不在 $R_5$ 中，只能是唯一的欧拉素因子。所以 $|T_5|\le b+1$。又因为 5 不整除 $N$，所有 $q\in R_5$ 的 $\Phi_5(q)$ 全由 $T_5$ 支持，故 $|R_5|\le u=B_5(b+1)$。[Fletcher–Nielsen–Ochem 主定理](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf)

每个 $T_3$ 中的非欧拉素因子都属于 $R_5$，因此 $|T_3|\le u+1$。对 $q\in R_3$ 有 $q\equiv5\bmod6$，且 $q^2+q+1$ 的全部素因子属于 $T_3$。已证的三次分圆值界给出 $|R_3|\le3(5^{u+1}-1)$。计入可能在两类以外的欧拉素因子，得到显式上界

$$
K(b)=1+u+3(5^{u+1}-1),\qquad
\omega(N)\le K(b),\qquad N<2^{4^{K(b)}} .
$$

函数随 $b$ 单调，故定理 1 的首个结论成立。反过来，若 $c$ 是模 15 余 $4,7,13$ 的素因子总数，则 $|T_3|\le c+1$，$|R_3|\le v=3(5^{c+1}-1)$，$|T_5|\le v+1$，从而 $\omega(N)\le1+v+B_5(v+1)$。两种上界分别使用一侧的素因子个数。[三次分圆值界，定理 1](../exponent-congruences/cubic-support.md)

例如，文献中的下界以及对两类指数的计数给出 $|R_5|\ge17\,864$。若 $|T_5|\le1$，则单素数支持的改进界给出 $|R_5|\le1$。故至少有两个不同的素因子模 5 余 1，其中至少一个模 15 余 11。在只允许非欧拉指数 2 和 4 的情形，这个单一剩余类中有超过 $10^{15}$ 个素因子。[两类指数的计数](../exponent-congruences/least-prime.md) · [七个剩余类的计数与指数重数](mixed-factor-count.md)
