---
title: "混合指数的奇完全数至少有 172414 个素因子模 15 余 11"
author: "Codex"
date: "2026-11-24"
abstract: "若每个非欧拉指数加一都被 3 或 5 整除，则在 10³⁵⁰⁰ 以内至少有 172414 个不同素因子模 15 余 11，在 10⁴⁰⁰ 以内至少有 142857 个。对于任意最小素因子 P 大于一亿的奇完全数，分别计数指数加一不能被 3 整除和不能被 5 整除的素因子，得到下界 (P−1)/22 和 (P−1)/70。证明结合两层分圆整除、Brun–Titchmarsh 不等式与大筛。"
formula: 'b(10^{3500})+\tfrac15 k>(P-1)/580'
---

# 混合指数的奇完全数至少有 172414 个素因子模 15 余 11

设 $N$ 是奇完全数，$P$ 是其最小素因子，$e_q=v_q(N)$。以下涉及 $q\mid N$ 的集合只计不同素数。记

$$
b(H)=\#\{q\mid N:q\le H,\ q\equiv11\bmod15\},\qquad
E=\{q\mid N:\gcd(e_q+1,15)=1\},\qquad k=|E|.
$$

**定理 1。** 对任何 $P>10^8$ 的奇完全数，有

$$
b(10^{400})+\frac{k}{5}>\frac{P-1}{700},\qquad
b(10^{3500})+\frac{k}{5}>\frac{P-1}{580}.
$$

若每个非欧拉指数加一都被 3 或 5 整除，则 $k\le1$。Fletcher、Nielsen 与 Ochem 已证明在此条件下 $P>10^8$。故模 15 余 11 的素因子，在 $10^{400}$ 以内至少有 **142857** 个，在 $10^{3500}$ 以内至少有 **172414** 个。这里允许所有指数任意大。[最小素因子的已知下界](https://www.lirmm.fr/~ochem/opn/OPNS_Adam_Pace.pdf) · [限制这个剩余类时的有效有限性](finite-residue-support.md)

**定理 2。** 对任意 $P>10^8$ 的奇完全数，令 $k_\ell=\#\{q\mid N:\ell\nmid e_q+1\}$。则

$$
k_3>(P-1)/22,\qquad k_5>(P-1)/70.
$$

特别地，在上述混合指数条件下，至少有 4,545,454 个非欧拉素因子满足 $5\mid e_q+1,\ 3\nmid e_q+1$，至少有 1,428,571 个满足 $3\mid e_q+1,\ 5\nmid e_q+1$。这里从 $k_\ell$ 中各扣除了至多一个欧拉素因子。[此前用乘积计数得到的 17864 与 230](../exponent-congruences/least-prime.md)

## 一个素因子能支持多大的倒数和

置 $X=P-1$。记 $W(S)=\sum_{q\in S}1/(q-1)$。有 $\log2<W(\{q:q\mid N\})<0.7$。下界来自 $\log(\sigma(q^{e_q})/q^{e_q})<1/(q-1)$。上界来自 $\prod_{q\mid N}(1+1/q)\le2$，并且 $1/(q-1)\le(1+2/X)\log(1+1/q)$。

固定一个素数 $r\ge P$ 及其一个非零剩余类 $a$。对下表的一行 $(m,C)$，记 $G_C(y)/(r-1)$ 为满足 $r<q\le\exp(y)$、$q\equiv a\bmod r$ 且 $q\bmod m\in C$ 的素数的倒数权重上界。若计入 $q<r$，每个根另外至多贡献 $1/X$。设 $h=5040,\ \epsilon=10^{-7}$，可以取

$$
\begin{aligned}
J_C&=\max_{\substack{u,v\bmod m\\(v,m)=1}}
 \sum_{\substack{1\le j\le h\\u+jv\bmod m\in C}}\frac1j,\\
G_C(y)&=J_C+\frac{2|C|}{\varphi(m)}(1+\epsilon)^2
\left(\frac1y+\log\frac{y}{\log(h/m)}\right).
\end{aligned}
$$

| 记号 | $m$ | $C$ | 所计的素数 |
| --- | ---: | --- | --- |
| $G_3$ | 6 | $5$ | 指数加一被 3 整除的可能底数 |
| $G_5$ | 10 | $3,7,9$ | 指数加一被 5 整除的可能底数 |
| $G_M$ | 30 | $7,13,19$ | 模 3 余 1 且指数加一被 5 整除的可能底数 |

在 $r<q\le hr$ 写 $q=a_0+jr$，其中 $1\le a_0<r$。由 $q-1\ge jr$ 得到 $J_C/(r-1)$。对更大的 $q$，中国剩余定理给出模 $mr$ 的 $|C|$ 个既约类。Brun–Titchmarsh 不等式将每个类的计数 $\pi(t;mr,c)$ 控制在 $2t/(\varphi(m)(r-1)\log(t/(mr)))$ 以下。对 $1/(t-1)$ 分部求和，舍去下端的非正项，置 $v=\log(\exp(y)/(mr))$。余下的上界是上式第二项除以 $r-1$，其中 $y$ 暂换成 $v$。因为 $v\ge\log(h/m)>1$ 时 $1/v+\log(v/\log(h/m))$ 随 $v$ 增加，且 $v\le y$，公式成立。若 $\exp(y)<hr$，第一项已经足够。[Montgomery–Vaughan，The large sieve](https://www.cambridge.org/core/journals/mathematika/article/abs/large-sieve/2B0269C6D6662EED351C6C06AA8D62ED)

## 单个剩余类的尾和

置 $R_\ell=\{q\mid N:\ell\mid e_q+1\}$，$T_\ell=\{r\mid N:r\equiv1\bmod\ell\}$。对 $q\in R_\ell$，有 $\Phi_\ell(q)\mid 2N$。由于 $3,5\nmid N$，$q\not\equiv1\bmod\ell$，且 $\Phi_\ell(q)$ 的所有素因子均属于 $T_\ell$。这说明 $R_3,R_5$ 在模 15 的既约类中分别只占四类、六类。

记 $M_\ell=\prod_{r\in T_\ell}(r-1)/(r-\ell)$。以下使用的常数为

$$
K_3=3.20,\quad K_5=4.05,\qquad d_3=9.5,\quad d_5=12.
$$

对任一模 15 的类和 $u\ge70$，该类中大于 $\exp(u)$ 的 $R_\ell$ 素因子的倒数权重至多为 $K_\ell M_\ell/(u-d_\ell)$。为核对适用范围，在已证的单类计数式中取 $Q=\lfloor\sqrt t/24\rfloor$，将部分乘积 $M_\ell(Q)$ 增大到 $M_\ell$。当 $t\ge\exp(70)>10^{28}$ 时，

$$
\#\{q\in R_\ell:q\le t,\ q\equiv a\bmod15\}
\le Q+\frac{1+15/576+15/t}{2c'_\ell}
       \frac{tM_\ell}{(\log t-d_\ell)^2},
$$

其中 $c'_3=(7/50)(1.149),\ c'_5=(6/55)(1.165)$。该计数式来自对 $t$ 的大筛，仅需 $3,5\nmid N$。对首项有 $Q(\log t)^2/t<10^{-4}$。将计数式除以 $(t-1)^2$ 从 $\exp(u)$ 积分到无穷，即得尾和。又有

$$
\log M_\ell\le(\ell-1)(1+\epsilon)W(T_\ell),\qquad M_3<4.1.
$$

后一式用了 $W(T_3)<0.7$。[单类计数式及分母常数](mixed-least-prime.md) · [五次分母证明](quintic-denominator.md)

## 单个指数的例外数

由于 $T_\ell\cap R_\ell$ 为空，有 $|T_\ell|\le k_\ell$。写 $f=k_\ell/X$。对 $q\in R_5,\ q\le\exp(y)$，选取某个 $r\mid\Phi_5(q)$，则 $q$ 落在模 $r$ 的四个根之一。加上 $r<q$ 与 $r>q$ 的贡献，得 $W(R_5\cap[1,\exp(y)])\le4(1+G_5(y))f$。

三次情形可以省一个常数。$\Phi_3(q)=q^2+q+1<(q+1)^2$，若是合数，必有素因子 $r<q$，因为 $r=q$ 不可能。若是素数，则这个值在 $T_3$ 内，且不同的 $q$ 给出不同的值。因此相应的有限和至多为 $(1+2G_3(y))f$。

最后加上 $k_\ell$ 个例外的权重 $f$ 和 $R_\ell$ 的各类尾和。以 $f$ 的假设上界代入，得到

| $\ell$ | $f$ 的上界 | $y$ | 全部权重上界，向上取整 |
| ---: | ---: | ---: | ---: |
| 3 | $1/22$ | 170 | $0.670261$ |
| 5 | $1/70$ | 320 | $0.690424$ |

两行都小于 $\log2$，故定理 2 成立。

## 从三次来源追到模 15 余 11

固定表中的 $H,y,z$，记 $A=b(H)/X,\ B=k/X$。如果 $t\in T_3\setminus E$，则 $t\in R_5$，且 $\Phi_5(t)\equiv2\bmod3$。它的素因子均模 5 余 1，因而至少有一个素因子 $r\equiv11\bmod15$。对 $t\le\exp(v),\ v=y$ 或 $z$，有 $r< (t+1)^4<H$。所以

$$
W(T_3\cap[1,\exp(v)])\le B+4(1+G_M(v))A.
$$

设 $q\in R_3,\ q\le\exp(y)$。若 $\Phi_3(q)$ 为合数，沿其一个小于 $q$ 的素因子 $t\in T_3$ 使用 $2G_3(y)/(t-1)$。若 $\Phi_3(q)=t$ 是素数且 $t\in E$，则这些 $q$ 至多有 $k$ 个。其余素数值 $t$ 属于 $T_3\setminus E$，于是存在 $r\equiv11\bmod15$ 整除

$$
F(q)=\Phi_5(\Phi_3(q)).
$$

$F$ 是首一八次多项式，模 $r$ 至多有八个根。又 $r\le F(q)<(q+1)^8<H$。所以截点以下的三次来源总权重至多为

$$
2G_3(y)\bigl[B+4(1+G_M(y))A\bigr]+B+8(1+G_3(y))A .
$$

对 $q\in R_5,\ q\le\exp(z)$，任取 $r\mid\Phi_5(q)$，则 $r<H$，并且 $r$ 模 15 余 11 或 1。余 1 的素因子全部属于 $E$。因此五次来源的有限权重至多为 $4(1+G_5(z))(A+B)$。

尾和中的 $M_3,M_5$ 还涉及 $H$ 以外的素因子。$T_5$ 中大于 $H$ 且不在 $E$ 中的素因子必在 $R_3$ 的单一类 $11\bmod15$ 内。类似地，$T_3\setminus E$ 在 $R_5$ 中只占三个类。因此依次有

$$
\begin{aligned}
W(T_5)&\le A+B+\frac{3.20(4.1)}{\log H-9.5},\\
W(T_3)&\le B+4(1+G_M(z))A+\frac{3(4.05)M_5}{z-12}.
\end{aligned}
$$

这里每个例外均可由 $B$ 一次计入。合并两个有限权重、四个三次尾和、六个五次尾和以及全部 $E$ 的权重，得到

$$
\begin{aligned}
W(\{q:q\mid N\}) &\le c_A A+c_B B+
\frac{4K_3M_3}{y-d_3}+\frac{6K_5M_5}{z-d_5},\\
c_A &=8G_3(y)(1+G_M(y))+8(1+G_3(y))+4(1+G_5(z)),\\
c_B &=2+2G_3(y)+4(1+G_5(z)).
\end{aligned}
$$

若 $A+B/5\le1/D$，则任一非负线性式 $uA+vB$ 至多为 $\max(u,5v)/D$。依次代入两个损失因子与最后的不等式，得到如下矛盾。

| $H$ | $D$ | $y$ | $z$ | 全部权重上界，向上取整 |
| --- | ---: | ---: | ---: | ---: |
| $10^{400}$ | 700 | 115 | 230 | $0.685605$ |
| $10^{3500}$ | 580 | 130 | 2000 | $0.689827$ |

两行都有 $\max(8(y+0.01),4(z+0.01))<\log H$，保证上述所选素因子确实落在各自的有限区间。验证器使用有理数，逐剩余类计算 $J_C$，对对数和指数函数的余项向外取界。[复算程序](verify_residue_incidence.py) · [相反方向三个剩余类的下界](complementary-residues.md)

在本目录运行下列命令即可复算两张常数表及相反方向的有限素数筛。

~~~sh
python -B -O verify_residue_incidence.py
~~~
