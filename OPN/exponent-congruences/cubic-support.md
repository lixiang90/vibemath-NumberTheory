---
title: "第三分圆多项式的光滑值与奇完全数的有限性"
author: "Codex"
date: "2026-11-24"
abstract: "给定 s 个模 3 余 1 的素数，证明至多 3(5^s−1) 个正整数 q 模 6 余 5 且 q²+q+1 只含这些素因子。用素数能否表示为 x²+27y² 细化此界。由此得到一族奇完全数的有效有限性，并用 77 个三次 Thue 方程完整求出四个素数集对应的全部底数。另给出两个分圆值之商为立方时的同余限制。"
formula: '\#\{q\equiv5\bmod6:\operatorname{rad}(q^2+q+1)\mid\prod_{p\in S}p,\ q>0\}\le3(5^{|S|}-1)'
---

# 第三分圆多项式的光滑值与奇完全数的有限性

设 $\Phi(q)=q^2+q+1$，$S$ 是由 $s$ 个模 3 余 1 的素数组成的集合。若一个正整数的素因子都属于 $S$，称它在 $S$ 上光滑。记

$$
b=\#\{p\in S:p\text{ 不能表示为 }u^2+27v^2,\ u,v\in\mathbb Z\},
\qquad B(S)=5^s+2(-1)^b5^{s-b}-3 .
$$

**定理 1。** 对所有 $q>0$、$q\equiv5\pmod6$，使 $\Phi(q)$ 在 $S$ 上光滑的底数至多有 $B(S)$ 个。特别地，上界 $3(5^s-1)$ 只依赖于素数的个数。底数可以是合数。

**定理 2。** 如果 $N$ 是奇完全数，记 $a=v_3(N)$，$k$ 为指数加一不被 3 整除的不同素因子数，$h=a+k$。那么

$$
\omega(N)\le h+\mathbf1_{\,3\mid N,\ 3\mid a+1}+3(5^h-1).
$$

因此，对任意给定的整数 $H\ge0$，满足 $v_3(N)+k\le H$ 的奇完全数只有有限多个。一个只用 $H$ 表示的上界是 $N<2^{4^{H+1+3(5^H-1)}}$。此处使用 Nielsen 关于不同素因子数的上界。下文还给出保留各素因子指数的形式。[Nielsen，An upper bound for odd perfect numbers](https://mathdept.byu.edu/~pace/upper_bound_web.pdf)

分圆多项式值只含少数给定素因子时，相应的丢番图方程可用于限制奇完全数。Yamada 对较高素数阶的 $\Phi_\ell(x)=p^m r$ 给出过解数上界。这里研究 $\Phi_3$，利用艾森斯坦整数的分解和模 2、模 3 的限制，把解数归约到二元三次型表示 1 的问题。[Yamada，An exponential Diophantine equation related to odd perfect numbers](https://www.iam.fmph.uniba.sk/amuc/ojs/index.php/amuc/article/download/1321/873/) · [Bennett，定理 1.4](https://personal.math.ubc.ca/~bennett/B-Cubic.pdf)

## 从艾森斯坦整数得到三次型

令 $\zeta^2+\zeta+1=0$，$R=\mathbb Z[\zeta]$。范数为 $\operatorname{N}(x+y\zeta)=x^2-xy+y^2$。该环关于范数是欧几里得环，单位为 $\pm1,\pm\zeta,\pm\zeta^2$。每个 $p\in S$ 分裂为 $\pi_p\bar\pi_p$。可把这两个因子都规范为模 $3R$ 余 1，因为六个单位恰好代表 $(R/3R)^\times$ 的六个元素，而共轭保持余类 1。

若 $\Phi(q)$ 在 $S$ 上光滑，注意到 $\Phi(q)=\operatorname{N}(q-\zeta)$。在 $q-\zeta$ 的分解中，$\pi_p,\bar\pi_p$ 至多出现其一，同时出现会使有理素数 $p$ 整除系数 $-1$。各指数除以 3 取余，便得

$$
q-\zeta=\zeta^2 d w^3,\qquad
d=\prod_{p\in S}d_p,\qquad
d_p\in\{1,\pi_p,\bar\pi_p,\pi_p^2,\bar\pi_p^2\},\qquad w\in R .
$$

这里的单位由 $q\equiv2\pmod3$ 决定。模 $3R$ 有 $q-\zeta\equiv\zeta^2$，$d\equiv1$，且每个可逆元素的立方为 $\pm1$，故前置单位只能是 $\pm\zeta^2$，负号可并入 $w^3$。再取模 $2R$，所得商环是四元域。$q$ 为奇数，故 $q-\zeta\equiv\zeta^2$，且 $w^3\equiv1$，这要求 $d\equiv1\pmod{2R}$。

写 $\gamma=\zeta^2d=g+h_0\zeta$、$w=x+y\zeta$。令

$$
U=x^3-3xy^2+y^3,\quad V=3xy(x-y),\quad
F_\gamma(x,y)=h_0 U+(g-h_0)V .
$$

由于 $w^3=U+V\zeta$，比较两项系数可得 $F_\gamma(x,y)=-1$ 与 $q=gU-h_0V$。三次型的判别式为 $81(g^2-gh_0+h_0^2)^2>0$。又 $h_0\equiv-1\pmod3$，故 $F_\gamma(x,1)$ 为具有三个不同根的三次多项式。Bennett 的定理 1.4 用于 $-F_\gamma$，给出至多十个整数对。$w,\zeta w,\zeta^2 w$ 三者不同，对应同一 $q$，不同 $q$ 的这三元组互不相交。所以每种允许的 $d$ 至多给出三个底数。[Bennett 原定理](https://personal.math.ubc.ca/~bennett/B-Cubic.pdf)

其中 $d=1$ 的情形可以删去。此时 $F_\gamma=-U$，Bennett 在定理 1.4 后的表中列出判别式 81 的 $U(x,y)=1$ 恰有六解。$w=1$ 与 $w=2-\zeta$ 的两个三元组已经给出六解，对应 $q=-1,-19$。两者均不计入正底数。

## 模 2 条件留下几种选择

设 $\pi_p=c+e\zeta\equiv1\pmod{3R}$。它模 $2R$ 余 1 当且仅当 $p$ 能表示为 $u^2+27v^2$。事实上，此时 $6\mid e$，有 $p=(c-e/2)^2+27(e/6)^2$。反过来，从 $p=u^2+27v^2$ 取因子 $(u+3v)+6v\zeta$，选取正负号使它模 $3R$ 余 1。它的常数项为奇数，因而模 $2R$ 也余 1。此性质在取共轭时保持。

当 $p$ 可如此表示时，五种 $d_p$ 均模 $2R$ 余 1。其余的 $p$ 在四元域乘法群的三个元素上分别给出 $1,2,2$ 种选择，其中第一项对应单位元。设只用 $j$ 个后一类素数时乘积为单位元的选法有 $n_j$ 种。则 $n_0=1$ 且 $n_{j+1}=2\cdot5^j-n_j$，因为已有积为 1 时有一种选择，其余各有两种选择。因此

$$
n_b=(5^b+2(-1)^b)/3.
$$

总共有 $5^{s-b}n_b$ 种允许的 $d$。删去 $d=1$ 后乘以 3，便得定理 1。

这一计算也约束方程 $\Phi(q)=p^j r$，其中 $p,r>3$ 为不同素数，$q$ 为奇素数，$j\ge1$。此时 $q\equiv5\pmod6$。若只有 $p$ 可表示为 $u^2+27v^2$，则无解。若只有 $r$ 可如此表示，则 $3\mid j$。若两者均不可如此表示，则必有 $3\nmid j$。这些断言来自 $d\equiv1\pmod{2R}$，因为 $r$ 的指数固定为 1。它们也可同时用于共用 $p,r$ 的两个分圆值。[Ochem–Zelinsky，问题 10](https://arxiv.org/html/2607.19746v1#S5) · [两个指数之差为偶数时的限制](../cyclotomic-pairs/even-gaps.md)

## 约数和中的重数与有限性

记 $T=\{t:t\mid N,\ t\text{ 为素数},\ t\equiv1\pmod3\}$。设 $Q$ 为满足 $q\equiv5\pmod6,\ 3\mid v_q(N)+1$ 的素因子集合。对 $q\in Q$ 有 $\Phi(q)\mid\sigma(q^{v_q(N)})$，且 $\Phi(q)$ 的每个素因子都属于 $T$，因为相应的乘法阶为 3。对任何 $S\subseteq T$ 都有

$$
|Q|\le B(S)+\sum_{t\in T\setminus S}v_t(N).
$$

确实，$\Phi(q)$ 在 $S$ 上光滑的 $q$ 至多有 $B(S)$ 个。其余每个 $q$ 至少贡献一个 $T\setminus S$ 中的素因子。汇总所有来源时，每个 $t$ 的总重数至多为 $v_t(\sigma(N))=v_t(N)$。

取 $S=T$。每个 $t\in T$ 若指数加一被 3 整除，便至少向 $\sigma(N)$ 贡献一个 3，所以这样的 $t$ 至多有 $a=v_3(N)$ 个。其余 $t$ 包含在 $k$ 个已计的素因子中，故 $|T|\le h$。计入这 $k$ 个素因子、指数加一被 3 整除的 $t\in T$ 以及 $Q$ 后，只可能剩下 3 自身。3 只有在 $3\mid N$ 且 $3\mid a+1$ 时才需另计。代入 $B(T)\le3(5^{|T|}-1)$ 即得定理 2。

例如，当每个非欧拉指数加一都被 3 或 5 整除时，已有 $3\nmid N$，并有至少 230 个非欧拉指数加一被 3 整除而不被 5 整除。这些素因子全属 $Q$。任取两个 $T$ 中素数，取作 $S$ 后有 $B(S)\le72$，故剩余 $T$ 中各素因子的指数之和至少为 158。若 $|T|\le3$，必有 $|T|=3$ 且三者都可表示为 $u^2+27v^2$，因为 $s=3,\ b=1,2,3$ 时 $B(S)$ 分别是 $72,132,120$。[混合指数的计数下界](least-prime.md)

## 两个分圆值之商为立方

**命题 3。** 设 $x,y$ 为同奇偶的正整数。若 $\Phi(y)\mid x-y$，且 $\Phi(x)/\Phi(y)$ 为整数的立方，则 $x=y$。

写 $x-y=m\Phi(y)$。由于 $\Phi(y)=(y-\zeta)(y-\zeta^2)$，商 $\alpha=(x-\zeta)/(y-\zeta)=1+m(y-\zeta^2)$ 属于 $R$，其范数为立方。没有有理素数整除 $\alpha$，否则也将整除 $x-\zeta$ 的系数 $-1$。因此惰性素数不能出现，每对分裂素因子最多出现其一，而 3 上方的分歧素因子至多出现一次。由范数为立方，各素因子实际出现的指数都是 3 的倍数，故 $\alpha$ 是单位乘立方。$x,y$ 同奇偶，$x-\zeta$ 与 $y-\zeta$ 在四元域中相同且非零。因而 $\alpha\equiv1\pmod{2R}$，前置单位只能为 $\pm1$。把符号并入立方，写 $\alpha=(u+v\zeta)^3$。

比较 $(y-\zeta)(u+v\zeta)^3=x-\zeta$ 的 $\zeta$ 系数，得到

$$
u^3-3(y+1)u^2v+3yuv^2+v^3=1.
$$

取 $n=3y+1\ge4$、$X=v$、$Y=-u$，左端为 $X^3-(n-1)X^2Y-(n+2)XY^2-Y^3$。Mignotte 对 Thomas 猜想的证明给出此方程右端为 1 时的全部解 $(1,0),(0,-1),(-1,1)$。对应的 $u+v\zeta$ 均满足立方为 1，故 $\alpha=1$，命题得证。[Mignotte，Verification of a conjecture of E. Thomas](https://doi.org/10.1006/jnth.1993.1043)

固定不同的素数 $p,r>3$，考虑使 $\Phi(q)=p^j r$、$j\ge1$ 的素数 $q$。当 $p,r$ 都能表示为 $u^2+27v^2$ 时，这样的 $q$ 至多十二个。当 $p$ 不能如此表示时，至多四个。若只有 $p$ 可如此表示，上文已排除所有解。

为证这些计数，把每个 $q$ 按 $j\bmod3$ 以及在 $p,r$ 上方各取了哪一个共轭因子来分组。一组内若有 $x>y$，记其指数为 $a>b$。则 $y-\zeta$ 整除 $x-\zeta$，故也整除有理整数 $x-y$。每个分裂素数在有理整数中的两个共轭因子的重数相等，所以 $\Phi(y)\mid x-y$。又因 $3\mid a-b$，两分圆值之商为立方。命题 3 给出矛盾。因此每组至多有一个底数，总计至多 $3\cdot2\cdot2=12$ 个。若只有 $r$ 可表示，则 $3\mid j$，只余四组。若两者均不可表示，则 $j$ 只能模 3 余 1 或 2。对每一余类，选定 $p$ 上方的因子后，$d\equiv1\pmod{2R}$ 从 $r$ 上方的两个共轭因子中唯一选定一个，故仍只余四组。

在两者均不可表示的情形，若有一对 $x>y$ 且 $3\mid a-b$，则两个素数上方的选择必须同时改取共轭。此时 $y-\zeta^2$ 整除 $x-\zeta$，并整除两者之和 $x+y+1$，故 $\Phi(y)\mid x+y+1$。商模 6 余 5，因为 $x,y\equiv5\pmod6$ 而 $\Phi(y)\equiv1\pmod6$。于是 $x\ge5\Phi(y)-y-1$。由 $y\ge5$ 有 $(y+1)/\Phi(y)\le6/31$，因此

$$
p^{a-2b}=\frac{r\Phi(x)}{\Phi(y)^2}>
r\left(5-\frac6{31}\right)^2>23r .
$$

## 四个指定素数集的全部正底数

对固定 $S$，上述三次型给出完整求解算法。还可在解方程前使用一个模 $p$ 检验。若 $d_p\ne1$，在 $p$ 上方选取整除 $\gamma$ 的素理想。由

$$
F_\gamma(x,y)=\frac{\gamma w^3-\bar\gamma\,\bar w^3}{\zeta-\bar\zeta}
$$

可知 $F_\gamma$ 模该理想等于 $h_0\bar w^3$，比例系数由 $x=1,y=0$ 确定。$h_0$ 模 $p$ 非零，且 $-1$ 是立方，故须有 $h_0^{(p-1)/3}\equiv1\pmod p$。对通过检验的三次型求解后，还须将 $q^2+q+1$ 中 $S$ 内的素因子除尽，只保留余因子为 1 的正 $q\equiv5\pmod6$。这一步允许 $w$ 起初取任意艾森斯坦整数。

令 $A=\{7,13,19,31\}$。下表包括每个集合对应的全部正底数，对 $q$ 的大小没有限制。

| $S$ | 模 2 后的选法 | 实际求解的三次型 | 全部 $q$ |
| --- | ---: | ---: | --- |
| $A$ | 205 | 9 | 5, 11, 191, 653 |
| $A\cup\{67\}$ | 1,045 | 24 | 5, 11, 29, 191, 653 |
| $A\cup\{79\}$ | 1,045 | 19 | 5, 11, 23, 191, 653 |
| $A\cup\{307\}$ | 1,025 | 25 | 5, 11, 17, 191, 653 |

计算采用 PARI 2.15.4 的 `thue(thueinit(poly,1),-1)`。参数 1 要求无条件认证。程序保留 $d=1$ 作对照，对每个返回解重算两个系数与原分圆值。使用 PARI/GP 可以用 `gp -q -f smooth_cubics.gp` 重算四行。[复算程序](smooth_cubics.gp) · [PARI 对 thue 与 thueinit 的说明](https://pari.math.u-bordeaux.fr/dochtml/html/Polynomials_and_power_series.html#thue)

这四行也给出五个偶指数定理中最后一步的另一证明。此前得到的五个分支里，$T=A$ 时另有一个已定素因子 17 或 23，故 $\omega(N)\le4+1+4=9$。另三行只有五个 $T$ 中素数和至多五个表内底数，故 $\omega(N)\le10$。这些分支均有 $3\nmid N$，此时 Nielsen 已证明 $\omega(N)\ge12$。[五个分支的推导](three-free.md) · [Nielsen 的十二个不同素因子定理](https://arxiv.org/abs/math/0602485)
